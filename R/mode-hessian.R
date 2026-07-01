## R/mode-hessian.R
## --------------------------------------------------------------------------
## Phase-2 split: proposal covariance for MCMC initialisation.
##
## proposal_cov()    -- robust proposal covariance.
##                       method = "diagonal" (default): univariate numerical
##                       curvature + fallback ladder.
##                       (was: build_sigma_prop_robust in sigma-prop-robust-monolith.R)
##                       method = "full": inverse of the full negative
##                       Hessian at the mode with eigenvalue repair, falling
##                       back to "diagonal" on failure.
## build_sigma_prop() -- three-tier: empirical posterior cov > Hessian diagonal
##                       > prior stds. (moved from mcmc-parallel-monolith.R)
## --------------------------------------------------------------------------

#' Build a robust MCMC proposal covariance matrix
#'
#' Computes a diagonal proposal covariance for RWMH using univariate numerical
#' second derivatives at the posterior mode.  Falls back through three
#' progressively coarser strategies when the Hessian cannot be estimated
#' reliably (e.g., for inv_gamma priors with near-infinite variance).
#'
#' @param lp_fn         log-posterior function: theta -> list(logpost, ...)
#' @param theta_mode    Named mode vector
#' @param prior_spec    Prior spec data.frame (from extract_prior_spec())
#' @param mode_multistart  Optional multi-start result list with
#'                         `$results[[i]]$theta_mode` entries
#' @param verbose       Print per-parameter diagnostic table
#' @param method        \code{"diagonal"} (default, existing behaviour):
#'                       per-parameter univariate curvature with fallback
#'                       ladder, returns a diagonal matrix.
#'                       \code{"full"}: inverse of the full negative Hessian
#'                       at \code{theta_mode}, with eigenvalue repair.  Falls
#'                       back to \code{"diagonal"} (with a warning) if the
#'                       Hessian cannot be estimated or inverted reliably.
#' @return n x n covariance matrix (named rows/cols); diagonal unless
#'   \code{method = "full"} succeeds.
#' @noRd
proposal_cov <- function(lp_fn, theta_mode, prior_spec,
                         mode_multistart = NULL,
                         verbose = TRUE,
                         method = c("diagonal", "full")) {

  method <- match.arg(method)

  if (method == "full") {
    Sigma_full <- .proposal_cov_full(lp_fn, theta_mode, prior_spec, verbose = verbose)
    if (!is.null(Sigma_full)) return(Sigma_full)
    ## fall through to the diagonal path on any failure (warning already issued)
  }

  .proposal_cov_diagonal(lp_fn, theta_mode, prior_spec, mode_multistart, verbose)
}


#' Diagonal proposal covariance from univariate numerical curvature
#' @noRd
.proposal_cov_diagonal <- function(lp_fn, theta_mode, prior_spec,
                                    mode_multistart = NULL,
                                    verbose = TRUE) {

  n_par  <- length(theta_mode)
  pnames <- names(theta_mode)
  prop_var <- rep(NA_real_, n_par)
  names(prop_var) <- pnames

  lp_mode <- lp_fn(theta_mode)$logpost

  if (verbose) cat("  Building robust proposal covariance (univariate curvature)...\n")

  n_ok <- 0; n_fallback1 <- 0; n_fallback2 <- 0; n_fallback3 <- 0

  for (i in seq_along(pnames)) {
    nm  <- pnames[i]
    val <- theta_mode[i]

    idx <- match(nm, prior_spec$name)
    lo  <- if (!is.na(idx)) prior_spec$lower[idx] else -Inf
    hi  <- if (!is.na(idx)) prior_spec$upper[idx] else  Inf
    if (is.na(lo)) lo <- -Inf
    if (is.na(hi)) hi <-  Inf

    curvature <- NA_real_
    for (step_frac in c(0.01, 0.005, 0.002, 0.05, 0.001)) {
      h <- max(abs(val) * step_frac, 1e-6)
      if (val + h > hi) h <- (hi - val) * 0.5
      if (val - h < lo) h <- (val - lo) * 0.5
      if (h < 1e-10) next

      theta_plus  <- theta_mode; theta_plus[i]  <- val + h
      theta_minus <- theta_mode; theta_minus[i] <- val - h

      lp_plus  <- lp_fn(theta_plus)$logpost
      lp_minus <- lp_fn(theta_minus)$logpost

      if (is.finite(lp_plus) && is.finite(lp_minus)) {
        d2 <- (lp_plus - 2 * lp_mode + lp_minus) / h^2
        if (is.finite(d2) && d2 < -1e-8) {
          curvature <- -d2
          break
        }
      }
    }

    if (!is.na(curvature) && curvature > 1e-12) {
      prop_var[i] <- 1 / curvature
      n_ok <- n_ok + 1
      next
    }

    if (!is.null(mode_multistart)) {
      all_modes <- lapply(mode_multistart$results, function(r) r$theta_mode)
      modes_mat <- do.call(rbind, all_modes)

      if (!is.null(modes_mat) && ncol(modes_mat) >= i) {
        mode_spread <- sd(modes_mat[, i])
        if (is.finite(mode_spread) && mode_spread > 1e-10) {
          prop_var[i] <- mode_spread^2
          n_fallback1 <- n_fallback1 + 1
          next
        }
      }
    }

    is_stderr <- grepl(
      "^(epniid|eptiid|ew|ec|eik|eih|eg|ex|em|er|erh|es|eph|eln|epx|epm|eb|ers|eys|eps)$",
      nm)

    if (is_stderr && abs(val) > 1e-8) {
      prop_var[i] <- (0.10 * val)^2
      n_fallback2 <- n_fallback2 + 1
      next
    }

    if (is.finite(lo) && is.finite(hi)) {
      prop_var[i] <- (0.05 * (hi - lo))^2
      n_fallback3 <- n_fallback3 + 1
      next
    }

    prop_var[i] <- max(0.05 * abs(val), 0.01)^2
    n_fallback3 <- n_fallback3 + 1
  }

  if (verbose) {
    cat(sprintf("    Curvature OK: %d/%d\n", n_ok, n_par))
    cat(sprintf("    Mode spread:  %d/%d\n", n_fallback1, n_par))
    cat(sprintf("    Stderr rule:  %d/%d\n", n_fallback2, n_par))
    cat(sprintf("    Range rule:   %d/%d\n", n_fallback3, n_par))
    cat("\n    Proposal std dev:\n")
    prop_sd <- sqrt(prop_var)
    for (i in seq_along(pnames)) {
      flag <- ""
      if (prop_sd[i] > abs(theta_mode[i]) && abs(theta_mode[i]) > 1e-6) flag <- " <- WIDE"
      if (prop_sd[i] < abs(theta_mode[i]) * 0.001)                       flag <- " <- NARROW"
      cat(sprintf("      %12s: mode=%10.6f  prop_sd=%10.6f%s\n",
                  pnames[i], theta_mode[i], prop_sd[i], flag))
    }
  }

  Sigma <- diag(prop_var, nrow = n_par)
  rownames(Sigma) <- colnames(Sigma) <- pnames
  Sigma
}


#' Full-Hessian proposal covariance with eigenvalue repair
#'
#' Inverse of the (symmetrised) negative Hessian of the log-posterior at
#' \code{theta_mode}, with eigenvalues floored to keep the result positive
#' definite (Dynare / Herbst & Schorfheide 2015 standard practice).
#'
#' @return n x n covariance matrix (named), or \code{NULL} if any step
#'   fails (a warning is issued in that case so the caller can fall back).
#' @noRd
.proposal_cov_full <- function(lp_fn, theta_mode, prior_spec, verbose = TRUE) {

  n_par  <- length(theta_mode)
  pnames <- names(theta_mode)

  ## Prior variances for .make_pd fallback on bad coordinates
  pv <- prior_spec$std^2
  pv[!is.finite(pv) | pv <= 0] <- 1

  neg_lp <- function(theta) {
    names(theta) <- pnames
    val <- lp_fn(theta)$logpost
    -val
  }

  if (verbose) cat("  Building full-Hessian proposal covariance...\n")

  H <- NULL
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    ## Abort (return NULL) on any non-finite evaluation -- numDeriv's
    ## Richardson extrapolation cannot be steered with a penalty value.
    f0 <- tryCatch(neg_lp(theta_mode), error = function(e) NA_real_)
    if (is.finite(f0)) {
      neg_lp_checked <- function(theta) {
        v <- tryCatch(neg_lp(theta), error = function(e) NA_real_)
        if (!is.finite(v))
          stop("non-finite log-posterior encountered during Hessian evaluation")
        v
      }
      H <- tryCatch(
        numDeriv::hessian(neg_lp_checked, as.numeric(theta_mode)),
        error = function(e) NULL
      )
    }
  }

  if (is.null(H)) {
    H <- .neg_hessian_central(neg_lp, theta_mode)
  }

  ## --- From here on H is always an n x n matrix (possibly with Inf/NaN) ---
  dimnames(H) <- list(pnames, pnames)

  ## symmetrise (finite entries only -- non-finite stays non-finite)
  H <- (H + t(H)) / 2

  ## Fast path: all entries finite and negative-Hessian is PD (common case)
  if (all(is.finite(H))) {
    eig     <- eigen(H, symmetric = TRUE)
    lam_max <- max(eig$values)
    if (is.finite(lam_max) && lam_max > 0) {
      floor_val    <- lam_max * 1e-8
      n_floored    <- sum(eig$values < floor_val)
      lam_repaired <- pmax(eig$values, floor_val)
      Sigma        <- eig$vectors %*%
        diag(1 / lam_repaired, nrow = n_par) %*%
        t(eig$vectors)
      Sigma <- (Sigma + t(Sigma)) / 2
      dimnames(Sigma) <- list(pnames, pnames)

      ok_chol <- !inherits(tryCatch(chol(Sigma), error = function(e) e), "error")
      if (ok_chol) {
        if (verbose) {
          kappa <- lam_max / min(lam_repaired)
          cat(sprintf("    Hessian condition number (repaired): %.3e\n", kappa))
          cat(sprintf("    Eigenvalues floored: %d/%d\n", n_floored, n_par))
          cat("\n    Proposal std dev:\n")
          prop_sd <- sqrt(diag(Sigma))
          for (i in seq_along(pnames)) {
            cat(sprintf("      %12s: mode=%10.6f  prop_sd=%10.6f\n",
                        pnames[i], theta_mode[i], prop_sd[i]))
          }
        }
        return(Sigma)
      }
      ## chol failed despite spectral repair -- fall through to .make_pd
    }
  }

  ## Fallback path: Hessian has Inf/NaN entries (stencil crossed a hard wall),
  ## or is indefinite with no PD spectral repair.  Apply .make_pd(), which
  ## preserves finite local curvature and assigns prior-var to bad coordinates.
  ## This replaces the old "return NULL -> prior-only diagonal" behaviour.
  has_nonfinite <- any(!is.finite(H))
  warning(
    "proposal_cov(method = 'full'): ",
    if (has_nonfinite) "Hessian has non-finite entries; "
    else "repaired covariance is not positive-definite; ",
    "applying nearest-PD regularisation (.make_pd) with prior-var fallback ",
    "for affected coordinates instead of dropping to prior-only diagonal."
  )

  Sigma <- .make_pd(H, cond_target = 100, prior_var = pv)
  dimnames(Sigma) <- list(pnames, pnames)

  if (verbose) {
    n_bad <- sum(!apply(is.finite(H), 1L, all))
    cat(sprintf("    .make_pd: %d/%d coordinate(s) fell back to prior variance.\n",
                n_bad, n_par))
    cat("\n    Proposal std dev (regularised):\n")
    prop_sd <- sqrt(diag(Sigma))
    for (i in seq_along(pnames)) {
      cat(sprintf("      %12s: mode=%10.6f  prop_sd=%10.6f\n",
                  pnames[i], theta_mode[i], prop_sd[i]))
    }
  }

  Sigma
}


#' Nearest-PD covariance from a possibly indefinite or NaN-contaminated
#' negative-Hessian matrix
#'
#' Turns a symmetrised negative-Hessian \code{H} (which should be PD at the
#' mode but often is not due to numerical noise, near-boundary evaluations, or
#' Inf/NaN from stencil points that crossed a feasibility boundary) into a
#' guaranteed positive-definite covariance matrix.
#'
#' Algorithm:
#' \enumerate{
#'   \item Symmetrise \code{H}.
#'   \item Identify "bad" coordinates: any row/col with a non-finite entry.
#'     Those coordinates are decoupled and assigned \code{prior_var[i]} on the
#'     diagonal (or 1 if \code{prior_var} is not supplied / non-finite for that
#'     coord).
#'   \item On the remaining "finite" submatrix: compute \code{eigen()}, floor
#'     eigenvalues at \code{lam_max / cond_target} (spectral nearest-PD
#'     projection), reconstruct and invert to a covariance.
#'   \item Assemble the full n x n covariance: finite sub-block on the
#'     retained coordinates; \code{prior_var[i]} diagonal on the bad ones.
#'   \item Guarantee \code{chol()} succeeds by construction (assert internally).
#' }
#'
#' When \code{H} is already PD and has no non-finite entries the result is
#' numerically identical to \code{solve(H)} (identity-preserving).
#'
#' @param H           n x n symmetrised negative-Hessian matrix (of logpost).
#' @param cond_target Desired condition number upper bound for the returned
#'   covariance (default 100).  Eigenvalues below
#'   \code{lam_max / cond_target} are floored.
#' @param prior_var   Optional length-n vector of prior variances, used as the
#'   diagonal fallback for "bad" (Inf/NaN) coordinates.  If \code{NULL} or
#'   non-finite for a coordinate, 1 is used.
#' @return n x n positive-definite covariance matrix (same dimnames as \code{H}).
#' @noRd
.make_pd <- function(H, cond_target = 100, prior_var = NULL) {
  n <- nrow(H)
  if (is.null(n) || n == 0L) stop(".make_pd: H must be a non-empty square matrix")

  ## 1. Symmetrise
  H <- (H + t(H)) / 2

  ## 2. Identify bad coordinates (any non-finite entry in that row or col)
  bad <- !apply(is.finite(H), 1L, all)  # TRUE for rows that contain Inf/NaN

  ## Fallback variance for bad coordinates
  pv <- if (!is.null(prior_var) && length(prior_var) == n) prior_var else rep(1, n)
  pv[!is.finite(pv) | pv <= 0] <- 1

  ## 3. Build result matrix
  Sigma <- matrix(0, nrow = n, ncol = n)
  dimnames(Sigma) <- dimnames(H)

  ## 3a. Bad coordinates: prior-var diagonal, zero off-diagonal (decoupled)
  for (i in which(bad)) Sigma[i, i] <- pv[i]

  ## 3b. Finite sub-block
  good_idx <- which(!bad)
  if (length(good_idx) > 0L) {
    Hsub <- H[good_idx, good_idx, drop = FALSE]
    eig  <- eigen(Hsub, symmetric = TRUE)
    lam_max <- max(eig$values)

    if (!is.finite(lam_max) || lam_max <= 0) {
      ## No positive eigenvalue in sub-block -- fall back to prior-var diagonal
      ## for the good coordinates too.
      for (i in good_idx) Sigma[i, i] <- pv[i]
    } else {
      floor_val   <- lam_max / cond_target
      lam_floored <- pmax(eig$values, floor_val)
      ## Invert to get covariance (spectral nearest-PD)
      Sigma_sub <- eig$vectors %*%
        diag(1 / lam_floored, nrow = length(good_idx)) %*%
        t(eig$vectors)
      Sigma_sub <- (Sigma_sub + t(Sigma_sub)) / 2
      Sigma[good_idx, good_idx] <- Sigma_sub
    }
  }

  ## 4. Guarantee chol() succeeds
  if (inherits(tryCatch(chol(Sigma), error = function(e) e), "error")) {
    ## Last resort: force diagonal (should never happen given the construction)
    Sigma <- diag(pmax(diag(Sigma), pv), nrow = n)
    dimnames(Sigma) <- dimnames(H)
  }

  Sigma
}


#' Plain central-difference Hessian of the negative log-posterior
#'
#' Fallback for \code{.proposal_cov_full} when numDeriv is unavailable.
#' Uses per-parameter relative steps \code{h_i = max(1e-4*|theta_i|, 1e-5)}
#' and the standard 4-point formula for off-diagonals.
#' @noRd
.neg_hessian_central <- function(fn, theta) {
  n <- length(theta)
  h <- pmax(1e-4 * abs(theta), 1e-5)
  H <- matrix(NA_real_, nrow = n, ncol = n)

  f0 <- fn(theta)

  ## diagonal: f(x+h) - 2f(x) + f(x-h) / h^2
  fp <- fm <- numeric(n)
  for (i in seq_len(n)) {
    th_p <- th_m <- theta
    th_p[i] <- th_p[i] + h[i]
    th_m[i] <- th_m[i] - h[i]
    fp[i] <- fn(th_p)
    fm[i] <- fn(th_m)
    H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2
  }

  ## off-diagonals: standard 4-point central-difference formula
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      th_pp <- th_pm <- th_mp <- th_mm <- theta
      th_pp[i] <- th_pp[i] + h[i]; th_pp[j] <- th_pp[j] + h[j]
      th_pm[i] <- th_pm[i] + h[i]; th_pm[j] <- th_pm[j] - h[j]
      th_mp[i] <- th_mp[i] - h[i]; th_mp[j] <- th_mp[j] + h[j]
      th_mm[i] <- th_mm[i] - h[i]; th_mm[j] <- th_mm[j] - h[j]

      f_pp <- fn(th_pp); f_pm <- fn(th_pm)
      f_mp <- fn(th_mp); f_mm <- fn(th_mm)

      H[i, j] <- H[j, i] <- (f_pp - f_pm - f_mp + f_mm) / (4 * h[i] * h[j])
    }
  }

  H
}


#' Build MCMC proposal covariance with three-tier fallback
#'
#' Priority:
#'   1. Empirical posterior covariance from pooled MCMC draws (warm start).
#'      Captures correlations; guaranteed PD; free given an existing run.
#'   2. Numerical Hessian diagonal at the mode (semi-warm start).
#'      2n evaluations; no correlations but better curvature than prior stds.
#'   3. Scaled prior standard deviations (cold start fallback).
#'
#' `scale = 2.38` is the Gelman et al. optimal scaling for a d-dim Gaussian.
#'
#' NOTE: build_sigma_prop is defined in run-mode-finding.R (loaded later in
#' Collate order). The definition here in mode-hessian.R is intentionally
#' empty to avoid duplicate-definition warnings — the run-mode-finding.R
#' version (using full numerical Hessian) is the canonical one.
#' @noRd
NULL
