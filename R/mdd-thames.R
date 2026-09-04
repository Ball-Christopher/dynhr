## R/mdd-thames.R
## --------------------------------------------------------------------------
## THAMES: Truncated Harmonic Mean Estimator for the marginal likelihood.
##
## Reference: Metodiev, M., Perrot-Dockès, M., Ouadah, S., Irons, N. J., &
## Raftery, A. E. (2023). Easily Computed Marginal Likelihoods from Posterior
## Simulation Using the THAMES Estimator.
## https://doi.org/10.48550/arXiv.2301.08847
## --------------------------------------------------------------------------

#' THAMES marginal-likelihood estimator
#'
#' Estimates \eqn{\log Z = \log \int L(\theta) p(\theta) \, d\theta} from
#' posterior draws and their unnormalized log-posterior values, using the
#' Truncated Harmonic Mean Estimator (THAMES) of Metodiev et al. (2023).
#'
#' @details
#' **Estimator.** Let \eqn{u(\theta_i) = L(\theta_i) p(\theta_i)} be the
#' unnormalized posterior evaluated at draw \eqn{\theta_i}.  By the
#' reciprocal-importance identity, for the uniform density \eqn{f} on an
#' ellipsoidal region \eqn{A}:
#' \deqn{\frac{1}{Z} = \frac{1}{N_A \cdot \text{Vol}(A)} \sum_{i : \theta_i \in A} \frac{1}{u(\theta_i)}}
#' where \eqn{N_A} is the total number of draws (not just those in \eqn{A}).
#' The region \eqn{A} is the Mahalanobis ellipsoid centred at the posterior
#' mean \eqn{\hat\mu} with shape \eqn{\hat\Sigma} and radius \eqn{r}:
#' \deqn{A = \{ \theta : (\theta - \hat\mu)^\top \hat\Sigma^{-1} (\theta - \hat\mu) \le r^2 \}}
#' with \eqn{\text{Vol}(A) = V_d \, r^d \, \sqrt{\det \hat\Sigma}}, where
#' \eqn{V_d = \pi^{d/2} / \Gamma(d/2 + 1)}.
#'
#' On the log scale (for numerical stability), using log-sum-exp:
#' \deqn{\widehat{\log Z} = \log N_A + \log \text{Vol}(A) - \text{logsumexp}_{i \in A}\{-\log u_i\}}
#'
#' **Radius default.** The default sets \eqn{r^2} equal to the empirical
#' median of the squared Mahalanobis distances \eqn{(\theta_i - \hat\mu)^\top
#' \hat\Sigma^{-1} (\theta_i - \hat\mu)}, so approximately half the draws
#' fall inside \eqn{A}.  Alternatively, specify \code{radius} or
#' \code{quantile} directly.
#'
#' **Bias control (data splitting).** When \code{split = TRUE} (default), the
#' first half of draws is used to estimate \eqn{(\hat\mu, \hat\Sigma)} and
#' define \eqn{A}; the second half supplies the harmonic-mean sum.  Both
#' orderings are averaged (mean on the log-reciprocal scale, before
#' inverting), which reduces estimator bias at the cost of a factor-of-two
#' effective sample size per half.
#'
#' **SE note.** The reported standard error treats the in-region draws as
#' approximately i.i.d. (delta method on the Monte Carlo variance of
#' \eqn{1/u_i}). For autocorrelated MCMC chains the true SE is larger; treat
#' \code{se} as a lower bound in that case.
#'
#' @param draws           Numeric matrix of posterior draws, \eqn{N \times d}
#'   (rows = draws, columns = parameters).
#' @param log_post_values Length-\eqn{N} numeric vector of unnormalized
#'   log-posterior values \eqn{\log u(\theta_i) = \log L_i + \log p_i} at
#'   each draw.  In dynhr these are the \code{post_logpost} or
#'   \code{logpost_trace} fields returned by the samplers.
#' @param radius          Optional positive scalar: the ellipsoid radius
#'   \eqn{r}.  If \code{NULL} (default), determined by \code{quantile}.
#' @param quantile        Quantile of the squared Mahalanobis distances used
#'   to set \eqn{r^2} when \code{radius = NULL}.  Default 0.5 (median).
#' @param split           Logical.  If \code{TRUE} (default), use data
#'   splitting to reduce bias: build the ellipsoid on the first half of draws
#'   and evaluate the estimator on the second, then average both orderings.
#' @param se_method       Character: \code{"iid"} (default) or \code{"batch"}.
#'   \itemize{
#'     \item \code{"iid"}: delta-method SE assuming the in-region draws are
#'       approximately i.i.d. (lower bound for autocorrelated MCMC chains).
#'     \item \code{"batch"}: batch-means SE that captures autocorrelation.
#'       With \eqn{N_{\text{eval}}} in-region contributions
#'       \eqn{v_i = \mathbf{1}\{i \in A\} \exp(-\log u_i) / (N_{\text{eval}} \cdot \text{Vol}(A))}
#'       in chain order, partition into \eqn{b = \lfloor\sqrt{N_{\text{eval}}}\rfloor}
#'       consecutive batches of size \eqn{m = \lfloor N_{\text{eval}}/b\rfloor};
#'       \eqn{\mathrm{SE}(1/Z) = N_{\text{eval}} \sqrt{\widehat{\mathrm{Var}}(\bar{M})/b}}
#'       where \eqn{\bar{M}} is the vector of batch means of \eqn{v_i}.
#'       Honest for autocorrelated chains; use this choice for MCMC output.
#'   }
#'   The default \code{"iid"} keeps all existing tests bit-identical.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{log_mdd}{Scalar estimate of \eqn{\log Z}.}
#'     \item{se}{Monte Carlo standard error of \code{log_mdd} (lower bound
#'       under autocorrelation; see Details).}
#'     \item{n_in_region}{Number of draws falling inside the ellipsoid (for
#'       the full-data estimator, or the second-half for \code{split = TRUE}).}
#'     \item{weight_ess}{Effective sample size of the in-region importance
#'       weights \eqn{1/u_i} (Kong-Liu-Wong). The THAMES reciprocal sum is
#'       dominated by the largest \eqn{1/u_i}; a small \code{weight_ess} means a
#'       few low-density draws dominate.}
#'     \item{reliability}{\code{weight_ess / n_in_region} in \eqn{[0,1]}. Near 1
#'       = well-behaved weights; near 0 = the estimate is dominated by a handful
#'       of draws and should be distrusted (the heavy-weight failure mode that
#'       PSIS Pareto-\eqn{\hat{k}} flags).}
#'     \item{reliable}{\code{TRUE} when \code{reliability > 0.1}. Inspect
#'       \code{reliability} directly for borderline cases.}
#'     \item{radius}{The ellipsoid radius \eqn{r} used.}
#'     \item{d}{Dimension \eqn{d}.}
#'     \item{n_used}{Number of draws after dropping non-finite
#'       \code{log_post_values}.}
#'   }
#'
#' @references
#' Metodiev, M., Perrot-Dockès, M., Ouadah, S., Irons, N. J., & Raftery,
#' A. E. (2023). Easily Computed Marginal Likelihoods from Posterior
#' Simulation Using the THAMES Estimator.
#' \doi{10.48550/arXiv.2301.08847}
#'
#' @examples
#' set.seed(42)
#' d <- 2
#' S <- matrix(c(1, 0.5, 0.5, 2), nrow = d)
#' m <- c(1, -1)
#' N <- 5000
#' draws <- MASS::mvrnorm(N, mu = m, Sigma = S)
#' # unnormalized log-posterior = log of unnormalized Gaussian kernel
#' log_u <- -0.5 * mahalanobis(draws, center = m, cov = S)
#' res <- thames_mdd(draws, log_u)
#' # analytic: logZ = (d/2)*log(2*pi) + 0.5*log(det(S))
#' logZ_true <- (d / 2) * log(2 * pi) + 0.5 * log(det(S))
#' cat(sprintf("Estimated logZ = %.4f  (truth = %.4f, se = %.4f)\n",
#'             res$log_mdd, logZ_true, res$se))
#'
#' @export
thames_mdd <- function(draws,
                       log_post_values,
                       radius    = NULL,
                       quantile  = 0.5,
                       split     = TRUE,
                       se_method = c("iid", "batch")) {

  se_method <- match.arg(se_method)

  ## ---- input validation -----------------------------------------------
  if (!is.matrix(draws))
    draws <- as.matrix(draws)
  log_post_values <- as.numeric(log_post_values)

  N_raw <- nrow(draws)
  d     <- ncol(draws)

  if (length(log_post_values) != N_raw)
    stop("length(log_post_values) must equal nrow(draws)")
  if (d < 1L)
    stop("draws must have at least one column")

  ## Drop non-finite log-posterior values
  ok <- is.finite(log_post_values)
  if (!all(ok)) {
    n_drop <- sum(!ok)
    message(sprintf("thames_mdd: dropping %d draw(s) with non-finite log_post_values",
                    n_drop))
    draws            <- draws[ok, , drop = FALSE]
    log_post_values  <- log_post_values[ok]
  }
  N <- nrow(draws)
  if (N <= 2L * d)
    stop(sprintf(
      "thames_mdd requires N > 2d; have N = %d, d = %d (after dropping non-finite values)",
      N, d))

  ## ---- log-volume helper ----------------------------------------------
  ## log Vol(ellipsoid) = log(V_d) + d*log(r) + 0.5*log(det(Sigma_hat))
  ## log(V_d) = (d/2)*log(pi) - lgamma(d/2 + 1)
  .log_vol <- function(r, log_det_sigma) {
    (d / 2) * log(pi) - lgamma(d / 2 + 1) + d * log(r) + 0.5 * log_det_sigma
  }

  ## ---- core THAMES estimate on a given split --------------------------
  ## build_idx: indices used to build (mu_hat, Sigma_hat)
  ## eval_idx:  indices used in the harmonic-mean sum
  .thames_half <- function(build_idx, eval_idx) {
    draws_b <- draws[build_idx, , drop = FALSE]
    draws_e <- draws[eval_idx,  , drop = FALSE]
    lp_e    <- log_post_values[eval_idx]

    ## Posterior mean and covariance from the build half
    mu_hat    <- colMeans(draws_b)
    Sigma_hat <- cov(draws_b)  # (N_b-1) denominator — fine for large N_b

    ## Guard against singular covariance
    ev <- eigen(Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev <= 0)) {
      ridge <- max(abs(ev)) * .Machine$double.eps^0.5 * d
      warning(sprintf(
        "thames_mdd: posterior covariance singular (min eigenvalue %.2e); adding ridge %.2e",
        min(ev), ridge))
      diag(Sigma_hat) <- diag(Sigma_hat) + ridge
    }

    ## Cholesky for fast Mahalanobis distances
    L <- tryCatch(
      chol(Sigma_hat),
      error = function(e) stop("thames_mdd: covariance matrix not positive definite after ridge correction")
    )
    log_det_sigma <- 2 * sum(log(diag(L)))

    ## Squared Mahalanobis distances for the EVAL half
    ## maha2_i = (x_i - mu)^T Sigma^{-1} (x_i - mu)
    ## = || L^{-T} (x_i - mu) ||^2 via backsolve
    ## maha2 = v' Sigma^{-1} v = || U^{-T} v ||^2 with Sigma = U'U (R's chol).
    ## Solve U' w = v (transpose = TRUE) to get w = U^{-T} v; using the plain
    ## solve U w = v would compute ||U^{-1} v||^2 = v'(UU')^{-1}v, a DIFFERENT
    ## (mis-rotated) ellipsoid. Same volume (det unchanged) so the estimator
    ## stays unbiased, but the wrong orientation inflates variance badly on
    ## strongly correlated posteriors (the DSGE case).
    centered_e <- sweep(draws_e, 2L, mu_hat, "-")
    Z_e        <- t(backsolve(L, t(centered_e), transpose = TRUE))   # N_e x d
    maha2_e    <- rowSums(Z_e^2)

    ## Determine radius from build half
    if (!is.null(radius)) {
      r <- radius
    } else {
      ## median squared Maha distance on the build half
      centered_b <- sweep(draws_b, 2L, mu_hat, "-")
      Z_b        <- t(backsolve(L, t(centered_b), transpose = TRUE))
      maha2_b    <- rowSums(Z_b^2)
      r <- sqrt(stats::quantile(maha2_b, probs = quantile, names = FALSE))
    }

    if (r <= 0)
      stop("thames_mdd: computed radius is <= 0; check your draws or specify radius manually")

    ## Which eval draws fall inside the ellipsoid?
    in_A <- maha2_e <= r^2
    n_in <- sum(in_A)
    if (n_in == 0L)
      stop("thames_mdd: no draws fall inside the ellipsoid; try a larger radius or quantile")

    N_e        <- length(eval_idx)
    log_vol    <- .log_vol(r, log_det_sigma)
    lp_in      <- lp_e[in_A]

    ## log(1/Z)^ = -log(N_e) - log(Vol) + logsumexp(-lp_in)
    ## => logZ^ = log(N_e) + log(Vol) - logsumexp(-lp_in)
    neg_lp_in  <- -lp_in
    lse        <- .logsumexp(neg_lp_in)
    log_recip  <- -log(N_e) - log_vol + lse

    ## SE on log(1/Z)^, two methods:
    ##
    ## "iid": delta-method on in-region draws (i.i.d. assumption — lower bound
    ##   under autocorrelation).
    ## "batch": batch-means SE in chain order captures autocorrelation.
    ##   v_i = 1{i in A} * exp(-lp_e[i]) / (N_e * Vol(A))
    ##   in chain order over ALL N_e eval draws (not just n_in).
    ##   (1/Z)^ = sum(v_i) = N_e * mean(v_i); so
    ##   SE(1/Z)^ = N_e * sqrt(Var(batch mean of v_i) / b)
    ##   SE(logZ^) = SE(1/Z)^ / (1/Z)^ = SE(1/Z)^ * Vol(A) / sum(exp(-lp_in))
    ##   where b = floor(sqrt(N_e)) batches of m = floor(N_e/b) draws each.
    ##
    ## Both methods work on a relative scale to avoid overflow.
    if (n_in > 1L) {
      if (identical(se_method, "iid")) {
        ## iid: delta-method on in-region draws (original code)
        rel     <- neg_lp_in - max(neg_lp_in)   # log(w_i / w_max)
        w_rel   <- exp(rel)
        mean_w  <- mean(w_rel)
        var_w   <- stats::var(w_rel)
        se_half <- sqrt(var_w / (n_in * mean_w^2))   # SE of mean(1/u) / w_max
        ## SE of log(mean(1/u)) via delta method = SE(mean(1/u)) / mean(1/u)
        se_logz <- se_half / mean_w
      } else {
        ## batch: batch-means on chain-ordered v_i in the eval half.
        ## v_i = 1{i in A} * exp(-lp_e[i]) in original scale (not /N_e*Vol,
        ## which is a constant that cancels in the ratio below).
        ## Work on a relative scale: anchor at max(-lp_in).
        log_anchor <- max(neg_lp_in)   # max over in-region values
        ## Length-N_e vector: exp(-lp_e[i] - log_anchor) for i in A; 0 otherwise
        v_full <- numeric(N_e)
        v_full[in_A] <- exp(neg_lp_in - log_anchor)

        b <- max(2L, floor(sqrt(N_e)))
        m <- floor(N_e / b)
        ## Drop the trailing (N_e - b*m) draws (standard batch-means convention)
        v_trim    <- v_full[seq_len(b * m)]
        batch_mat <- matrix(v_trim, nrow = m, ncol = b)    # m x b
        batch_means <- colMeans(batch_mat)                 # length b

        ## Var(mean v_i) ~ var(batch_means) / b (batch-means estimator)
        ## SE((1/Z)^ on relative scale) = N_e * sqrt(var(batch_means)/b)
        var_bm   <- stats::var(batch_means)                # corrected var of b values
        se_recip_rel <- sqrt(var_bm / b)                   # SE of mean(v_full)

        ## Convert to SE on logZ via delta method:
        ## (1/Z)^ (relative) = N_e * mean(v_full) = mean_w (relative scale)
        mean_w   <- mean(v_full)   # mean over all N_e draws (many zeros)
        if (mean_w <= 0) {
          se_logz <- NA_real_
        } else {
          se_logz <- se_recip_rel / mean_w
        }
      }
    } else {
      se_logz <- NA_real_
    }

    ## Reliability: effective sample size of the in-region importance weights
    ## w_i = 1/u_i (Kong-Liu-Wong 1994). The THAMES reciprocal sum is dominated
    ## by the largest 1/u_i (the lowest-density in-region draws); a small
    ## weight_ess / n_in means a few draws dominate and the estimate is
    ## unreliable -- the same heavy-weight failure mode that PSIS Pareto-khat
    ## flags, captured here by the robust ESS workhorse (no GPD fit needed).
    if (n_in > 1L) {
      w_rel      <- exp(neg_lp_in - max(neg_lp_in))   # 1/u_i up to a constant
      weight_ess <- (sum(w_rel))^2 / sum(w_rel^2)
    } else {
      weight_ess <- as.numeric(n_in)
    }

    list(log_recip  = log_recip,
         log_mdd    = -log_recip,
         se         = se_logz,
         n_in       = n_in,
         weight_ess = weight_ess,
         r          = r,
         N_e        = N_e)
  }

  ## ---- run estimator --------------------------------------------------
  if (split) {
    half1 <- seq_len(N %/% 2L)
    half2 <- seq(N %/% 2L + 1L, N)

    res_a <- .thames_half(build_idx = half1, eval_idx = half2)
    res_b <- .thames_half(build_idx = half2, eval_idx = half1)

    ## Average the two log(1/Z)^ estimates (then invert)
    log_recip_avg <- 0.5 * (res_a$log_recip + res_b$log_recip)

    ## Combined SE: propagate as sqrt(se_a^2 + se_b^2) / 2
    se_comb <- if (!is.na(res_a$se) && !is.na(res_b$se))
      0.5 * sqrt(res_a$se^2 + res_b$se^2)
    else
      NA_real_

    w_ess <- min(res_a$weight_ess, res_b$weight_ess)   # conservative
    list(
      log_mdd     = -log_recip_avg,
      se          = se_comb,
      n_in_region = res_b$n_in,        # second-half eval (conventional report)
      weight_ess  = w_ess,
      reliability = w_ess / res_b$n_in,
      reliable    = (w_ess / res_b$n_in) > 0.1,
      radius      = res_b$r,
      d           = d,
      n_used      = N
    )
  } else {
    res <- .thames_half(build_idx = seq_len(N), eval_idx = seq_len(N))
    list(
      log_mdd     = res$log_mdd,
      se          = res$se,
      n_in_region = res$n_in,
      weight_ess  = res$weight_ess,
      reliability = res$weight_ess / res$n_in,
      reliable    = (res$weight_ess / res$n_in) > 0.1,
      radius      = res$r,
      d           = d,
      n_used      = N
    )
  }
}


#' THAMES MDD from a dynhr_chains object
#'
#' Convenience wrapper around \code{\link{thames_mdd}} that extracts posterior
#' draws and their theta-space unnormalized log-posterior values from a
#' \code{dynhr_chains} object (as returned by \code{\link{run_full_estimation}}
#' or the public samplers \code{\link{mcmc}}, \code{\link{nuts}},
#' \code{\link{smc}}).
#'
#' @details
#' **Draw extraction.** The combined post-warmup draw matrix is taken from
#' \code{chains$chain} (rows = draws in chain order, columns = parameters).
#' For multi-chain results (when \code{chains$chain_list} is present), the
#' per-chain \code{post_logpost} vectors are concatenated in chain order to
#' give the matching log-posterior vector, preserving the within-chain
#' autocorrelation structure that \code{se_method = "batch"} exploits.
#' For single-chain results, \code{chains$post_logpost} is used directly.
#'
#' **Multi-chain handling.** Chains are pooled (concatenated in chain 1, 2,
#' ... order) into a single draw matrix before calling \code{thames_mdd}.
#' Pooling is valid for the point estimate; the batch-means SE will capture
#' within-chain autocorrelation but may slightly understate the cross-chain
#' variance component.  For a per-chain breakdown, call \code{thames_mdd}
#' directly on each \code{chain_list[[k]]$chain} and
#' \code{chain_list[[k]]$post_logpost}.
#'
#' **Never errors.** Any failure (too few draws, singular covariance, missing
#' log-posterior) returns \code{list(log_mdd = NA_real_, se = NA_real_, ...)}
#' with a \code{message()}, so it cannot break or slow estimation.
#'
#' @param chains  A \code{dynhr_chains} object.
#' @param ...     Additional arguments forwarded to \code{\link{thames_mdd}}
#'   (e.g. \code{radius}, \code{quantile}, \code{split}, \code{se_method}).
#'   Default: \code{se_method = "batch"}.
#'
#' @return A named list with elements \code{log_mdd}, \code{se},
#'   \code{n_in_region}, \code{radius}, \code{d}, \code{n_used}.
#'   All fields are \code{NA} when the computation fails.
#'
#' @seealso \code{\link{thames_mdd}}
#'
#' @examples
#' \dontrun{
#' result <- run_full_estimation(...)
#' thames_mdd_from_chains(result$chains)
#' }
#' @export
thames_mdd_from_chains <- function(chains, ...) {
  ## Safe return value for any failure path.
  .na_result <- list(log_mdd     = NA_real_,
                     se          = NA_real_,
                     n_in_region = NA_integer_,
                     radius      = NA_real_,
                     d           = NA_integer_,
                     n_used      = NA_integer_)

  tryCatch({
    ## ---- Extract draws matrix -------------------------------------------
    if (!inherits(chains, "dynhr_chains") && !is.list(chains))
      stop("'chains' must be a dynhr_chains object or list")

    draws <- chains$chain
    if (!is.matrix(draws) || nrow(draws) == 0L || ncol(draws) == 0L) {
      message("thames_mdd_from_chains: chains$chain is absent or empty; returning NA")
      return(.na_result)
    }

    ## ---- Extract log-posterior vector in chain order --------------------
    ## Multi-chain: chain_list[[k]]$post_logpost, concatenated in order.
    ## Single-chain: chains$post_logpost directly.
    cl <- chains$chain_list
    if (!is.null(cl) && length(cl) >= 1L) {
      ## Multi-chain: collect per-chain post_logpost values.
      ## Concatenate in chain order (same order as combined chain matrix).
      lp_parts <- lapply(cl, function(ch) ch$post_logpost)
      ok_parts <- !sapply(lp_parts, is.null)
      if (!any(ok_parts)) {
        message("thames_mdd_from_chains: no post_logpost in chain_list; returning NA")
        return(.na_result)
      }
      ## Concatenate only the non-NULL parts, but only use chains that also
      ## contributed to chains$chain (same number of rows).
      lp_vec <- unlist(lp_parts[ok_parts])
    } else {
      ## Single-chain path
      lp_vec <- chains$post_logpost
      if (is.null(lp_vec)) {
        message("thames_mdd_from_chains: chains$post_logpost not found; returning NA")
        return(.na_result)
      }
    }
    lp_vec <- as.numeric(lp_vec)

    ## ---- Dimension check ------------------------------------------------
    n_draws <- nrow(draws)
    if (length(lp_vec) != n_draws) {
      message(sprintf(
        paste0("thames_mdd_from_chains: length(logpost) = %d != nrow(draws) = %d",
               " (possible chain-list mismatch); returning NA"),
        length(lp_vec), n_draws))
      return(.na_result)
    }

    ## ---- Call thames_mdd with batch SE by default -----------------------
    dots <- list(...)
    if (is.null(dots$se_method)) dots$se_method <- "batch"

    do.call(thames_mdd, c(list(draws = draws, log_post_values = lp_vec), dots))

  }, error = function(e) {
    message("thames_mdd_from_chains: ", conditionMessage(e), "; returning NA")
    .na_result
  })
}
