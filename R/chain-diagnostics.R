## R/chain-diagnostics.R
## --------------------------------------------------------------------------
## Standalone MCMC chain diagnostics: rank-normalised split-Rhat, bulk/tail
## effective sample size and Monte-Carlo standard error, callable on an
## arbitrary draws object. Provided so callers (and the paper family, which
## hand-roll a Geyer ESS + Rhat stack across many scripts) do not each
## re-implement it. Formulas follow Vehtari, Gelman, Simpson, Carpenter &
## Buerkner (2021); the estimators themselves are the package's single set in
## R/diag-helpers.R (.d5_*), shared with D5 and the estimation runners.
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

#' MCMC diagnostics: effective sample size, split-Rhat, MCSE
#'
#' Computes, per parameter, the rank-normalised split-\eqn{\hat R}, the bulk
#' and tail effective sample sizes of Vehtari, Gelman, Simpson, Carpenter &
#' Buerkner (2021), the ESS of the posterior mean, the posterior mean/sd, and
#' the Monte-Carlo standard error of the mean (\code{sd / sqrt(ess)}).
#'
#' Every ESS is the multi-chain estimator of Vehtari et al. (2021, eqs 10-11):
#' the between-chain-aware autocorrelation
#' \eqn{\hat\rho_t = 1 - (W - \bar\gamma_t) / \widehat{var}^+}, truncated by
#' Geyer's initial positive sequence on the pairs
#' \eqn{(\hat\rho_0 + \hat\rho_1), (\hat\rho_2 + \hat\rho_3), \ldots} and made
#' monotone by Geyer's initial monotone sequence. It is NOT capped at the draw
#' count: antithetic chains (negative lag-1 autocorrelation, common for NUTS)
#' legitimately have ESS > N, bounded as in Stan by \eqn{N \log_{10} N}.
#'
#' @param draws One of: a numeric vector (one chain, one parameter); a matrix
#'   (iterations x parameters, one chain); a 3-D array \code{[iterations,
#'   chains, parameters]}; or a list of per-chain matrices (each iterations x
#'   parameters, identical column names). Column names become parameter names.
#' @param split Logical (default \code{TRUE}); if \code{TRUE} each chain is
#'   split in half before the Rhat/ESS computation (split-Rhat, which also
#'   detects within-chain non-stationarity a plain Rhat misses).
#' @return A \code{data.frame}, one row per parameter, with columns
#'   \code{param}, \code{mean}, \code{sd}, \code{ess} (ESS of the mean: the
#'   estimator above on the raw draws, the one \code{mcse} uses),
#'   \code{ess_bulk} (the same on rank-normalised draws), \code{ess_tail}
#'   (minimum ESS of the 5\% and 95\% quantile indicators), \code{rhat}
#'   (maximum of the rank-normalised split-\eqn{\hat R} of the draws and of
#'   the draws folded about their median; \code{NA} for a single unsplit
#'   chain) and \code{mcse}. Constant draws give \code{NA} ESS and R-hat.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#'   Buerkner, P.-C. (2021). Rank-normalization, folding, and localization: an
#'   improved \eqn{\hat R} for assessing convergence of MCMC. \emph{Bayesian
#'   Analysis} 16(2), 667-718.
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

  ## One estimator set for the whole package (R/diag-helpers.R, verified there
  ## against the `posterior` package). This function used to carry its own
  ## copy, which paired autocorrelations as (rho1+rho2), (rho3+rho4) instead of
  ## Geyer's (rho0+rho1), (rho2+rho3), clamped ESS at the draw count (wrong for
  ## antithetic chains) and had neither rank normalisation nor the folded
  ## R-hat.
  sp <- if (split) .d5_split else identity

  out <- data.frame(param = pnames, mean = NA_real_, sd = NA_real_,
                    ess = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
                    rhat = NA_real_, mcse = NA_real_,
                    stringsAsFactors = FALSE)

  for (p in seq_len(np)) {
    X  <- matrix(A[, , p], nrow = ni, ncol = nc)      # iterations x chains
    Xs <- sp(X)
    n  <- nrow(Xs)
    W  <- mean(apply(Xs, 2L, stats::var))             # within-chain variance
    var_plus <- ((n - 1) / n) * W +
      (if (ncol(Xs) > 1L) stats::var(colMeans(Xs)) else 0)
    out$mean[p] <- mean(X)
    out$sd[p]   <- sqrt(var_plus)

    ## Constant (or non-finite) draws carry no autocorrelation information.
    if (.d5_degenerate(Xs)) next

    ess_mean <- .d5_ess_basic(Xs)
    out$ess[p]      <- ess_mean
    out$ess_bulk[p] <- .d5_ess_basic(.d5_zscale(Xs))
    out$ess_tail[p] <- min(vapply(c(0.05, 0.95), function(q) {
      qv <- stats::quantile(X, q, names = FALSE)
      .d5_ess_basic(sp(X <= qv) + 0)
    }, numeric(1)))
    if (ncol(Xs) > 1L)
      out$rhat[p] <- max(.d5_rhat_basic(.d5_zscale(Xs)),
                         .d5_rhat_basic(.d5_zscale(
                           sp(abs(X - stats::median(X))))))
    out$mcse[p] <- if (is.finite(ess_mean) && ess_mean > 0)
      sqrt(var_plus) / sqrt(ess_mean) else NA_real_
  }
  out
}
