## R/chain-diagnostics.R
## --------------------------------------------------------------------------
## Standalone MCMC chain diagnostics: effective sample size (Geyer initial
## monotone sequence, multi-chain autocovariance) + split-Rhat + Monte-Carlo
## standard error, callable on an arbitrary draws object. Provided so callers
## (and the paper family -- P2/P6 hand-roll a Geyer ESS + Rhat stack across
## many scripts) do not each re-implement it. Formulas follow Vehtari, Gelman,
## Simpson, Carpenter & Buerkner (2021) / BDA3 ch.11.
## --------------------------------------------------------------------------


#' Coerce a draws object to an iterations x chains x parameters array.
#' Accepts: a numeric vector (1 chain, 1 param); a matrix (iterations x
#' parameters, 1 chain); a 3-D array already in [iter, chain, param] layout; or
#' a list of per-chain matrices (each iterations x parameters, identical
#' colnames). Parameter names are taken from the column names where present.
#' @noRd
.cd_as_array <- function(draws) {
  if (is.array(draws) && length(dim(draws)) == 3L) return(draws)
  if (is.list(draws) && !is.data.frame(draws)) {
    mats <- lapply(draws, function(m) as.matrix(m))
    ni <- nrow(mats[[1L]]); np <- ncol(mats[[1L]])
    if (!all(vapply(mats, function(m) nrow(m) == ni && ncol(m) == np, TRUE)))
      stop("chain_diagnostics: all chains must have identical dimensions.")
    A <- array(0, dim = c(ni, length(mats), np),
               dimnames = list(NULL, NULL, colnames(mats[[1L]])))
    for (c in seq_along(mats)) A[, c, ] <- mats[[c]]
    return(A)
  }
  m <- as.matrix(draws)
  if (is.null(dim(m)) || length(dim(m)) == 1L) m <- matrix(m, ncol = 1L)
  array(m, dim = c(nrow(m), 1L, ncol(m)),
        dimnames = list(NULL, NULL, colnames(m)))
}

#' Autocovariance (lags 0..n-1) of a single chain via FFT (biased 1/n estimator,
#' the convention used by the Geyer/Vehtari ESS).
#' @noRd
.cd_autocov <- function(x) {
  n <- length(x)
  xc <- x - mean(x)
  nfft <- stats::nextn(2L * n)
  f <- stats::fft(c(xc, rep(0, nfft - n)))
  ac <- Re(stats::fft(f * Conj(f), inverse = TRUE)) / nfft
  ac[seq_len(n)] / n
}

#' MCMC diagnostics: effective sample size, split-Rhat, MCSE
#'
#' Computes, per parameter, the effective sample size (ESS) via Geyer's initial
#' monotone positive sequence on the multi-chain-combined autocorrelation, the
#' split-Rhat potential-scale-reduction factor, the posterior mean/sd, and the
#' Monte-Carlo standard error (\code{sd / sqrt(ESS)}).
#'
#' @param draws One of: a numeric vector (one chain, one parameter); a matrix
#'   (iterations x parameters, one chain); a 3-D array \code{[iterations,
#'   chains, parameters]}; or a list of per-chain matrices (each iterations x
#'   parameters, identical column names). Column names become parameter names.
#' @param split Logical (default \code{TRUE}); if \code{TRUE} each chain is
#'   split in half before the Rhat/ESS computation (split-Rhat, which also
#'   detects within-chain non-stationarity a plain Rhat misses).
#' @return A \code{data.frame}, one row per parameter, with columns
#'   \code{param}, \code{mean}, \code{sd}, \code{ess}, \code{rhat}, \code{mcse}.
#' @examples
#' set.seed(1)
#' draws <- matrix(rnorm(4000), ncol = 2, dimnames = list(NULL, c("a", "b")))
#' chain_diagnostics(draws)
#' @seealso \code{sbc_uniformity_test} (SBC rank-uniformity diagnostics)
#' @export
chain_diagnostics <- function(draws, split = TRUE) {
  A <- .cd_as_array(draws)
  ni <- dim(A)[1L]; nc <- dim(A)[2L]; np <- dim(A)[3L]
  pnames <- dimnames(A)[[3L]]
  if (is.null(pnames)) pnames <- paste0("p", seq_len(np))
  if (ni < 4L) stop("chain_diagnostics: need at least 4 iterations per chain.")

  ## Split each chain in half (split-Rhat): 2*nc half-chains of length m.
  if (split && ni >= 4L) {
    m <- ni %/% 2L
    B <- array(0, dim = c(m, 2L * nc, np))
    for (c in seq_len(nc)) {
      B[, 2L * c - 1L, ] <- A[seq_len(m), c, ]
      B[, 2L * c,       ] <- A[m + seq_len(m), c, ]
    }
    A <- B; ni <- m; nc <- 2L * nc
  }

  out <- data.frame(param = pnames, mean = NA_real_, sd = NA_real_,
                    ess = NA_real_, rhat = NA_real_, mcse = NA_real_,
                    stringsAsFactors = FALSE)

  for (p in seq_len(np)) {
    X <- matrix(A[, , p], nrow = ni, ncol = nc)      # iterations x chains
    chain_means <- colMeans(X)
    chain_vars  <- apply(X, 2L, stats::var)
    W <- mean(chain_vars)                             # within-chain variance
    grand_mean <- mean(chain_means)

    if (nc > 1L) {
      B <- ni * stats::var(chain_means)               # between-chain variance
      var_plus <- ((ni - 1) / ni) * W + B / ni
      rhat <- if (W > 0) sqrt(var_plus / W) else NA_real_
    } else {
      var_plus <- W
      rhat <- NA_real_                                # Rhat undefined for 1 chain
    }

    ## Multi-chain-combined autocorrelation: rho_t = 1 - (W - mean_m gamma_{t,m})/var_plus.
    if (W <= 0 || var_plus <= 0) {
      ess <- as.numeric(ni * nc)                      # constant -> treat as iid
    } else {
      acov <- vapply(seq_len(nc), function(c) .cd_autocov(X[, c]), numeric(ni))
      mean_acov <- rowMeans(acov)                     # length ni, lags 0..ni-1
      rho <- 1 - (W - mean_acov) / var_plus           # rho[1] = lag 0 = 1
      ## Geyer initial monotone positive sequence on paired sums.
      max_pairs <- (ni - 1L) %/% 2L
      tau <- 1                                        # = rho_0 (lag-0 sum term)
      if (max_pairs >= 1L) {
        P <- vapply(seq_len(max_pairs), function(k) rho[2L * k] + rho[2L * k + 1L], 0)
        ## truncate at the first non-positive pair sum
        first_neg <- which(P <= 0)
        kmax <- if (length(first_neg)) first_neg[1L] - 1L else length(P)
        if (kmax >= 1L) {
          Pk <- P[seq_len(kmax)]
          Pk <- cummin(Pk)                            # Geyer monotone: non-increasing pair sums
          tau <- tau + 2 * sum(Pk)                    # tau = 1 + 2*sum_{t>=1} rho_t
        }
      }
      ess <- if (tau > 0) (ni * nc) / tau else as.numeric(ni * nc)
      ess <- min(ess, ni * nc)                        # ESS cannot exceed the draw count
    }

    out$mean[p] <- grand_mean
    out$sd[p]   <- sqrt(var_plus)
    out$ess[p]  <- ess
    out$rhat[p] <- rhat
    out$mcse[p] <- if (is.finite(ess) && ess > 0) sqrt(var_plus) / sqrt(ess) else NA_real_
  }
  out
}
