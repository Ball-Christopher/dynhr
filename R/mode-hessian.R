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
#' @param grad_fn       Optional theta-space gradient of the log-posterior
#'                       (the mode stage's \code{make_posterior_grad()}
#'                       closure). With \code{method = "full"} the Hessian is
#'                       then the central difference of the gradient
#'                       (\code{grad_hessian()}, 2n calls) instead of
#'                       numDeriv on the log-posterior; ignored by
#'                       \code{"diagonal"}.
#' @return n x n covariance matrix (named rows/cols); diagonal unless
#'   \code{method = "full"} succeeds.
#' @noRd
proposal_cov <- function(lp_fn, theta_mode, prior_spec,
                         mode_multistart = NULL,
                         verbose = TRUE,
                         method = c("diagonal", "full"),
                         grad_fn = NULL) {

  method <- match.arg(method)

  if (method == "full") {
    Sigma_full <- .proposal_cov_full(lp_fn, theta_mode, prior_spec, verbose = verbose,
                                     grad_fn = grad_fn)
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

  if (verbose) .dynhr_cat("  Building robust proposal covariance (univariate curvature)...\n")

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
    .dynhr_cat(sprintf("    Curvature OK: %d/%d\n", n_ok, n_par))
    .dynhr_cat(sprintf("    Mode spread:  %d/%d\n", n_fallback1, n_par))
    .dynhr_cat(sprintf("    Stderr rule:  %d/%d\n", n_fallback2, n_par))
    .dynhr_cat(sprintf("    Range rule:   %d/%d\n", n_fallback3, n_par))
    .dynhr_cat("\n    Proposal std dev:\n")
    prop_sd <- sqrt(prop_var)
    for (i in seq_along(pnames)) {
      flag <- ""
      if (prop_sd[i] > abs(theta_mode[i]) && abs(theta_mode[i]) > 1e-6) flag <- " <- WIDE"
      if (prop_sd[i] < abs(theta_mode[i]) * 0.001)                       flag <- " <- NARROW"
      .dynhr_cat(sprintf("      %12s: mode=%10.6f  prop_sd=%10.6f%s\n",
                  pnames[i], theta_mode[i], prop_sd[i], flag))
    }
  }

  Sigma <- diag(prop_var, nrow = n_par)
  rownames(Sigma) <- colnames(Sigma) <- pnames
  Sigma
}


#' numDeriv Hessian of the negative log-posterior for .proposal_cov_full
#'
#' NULL on any non-finite evaluation -- numDeriv's Richardson extrapolation
#' cannot be steered with a penalty value. A PSKF posterior is differenced on
#' the pruning selection made at \code{theta_mode}, which is evaluated first
#' and so records it (.pskf_freeze_open): a second difference across a
#' selection switch would otherwise return jump / h^2. Functions that run no
#' PSKF filter are unaffected.
#' @noRd
.proposal_numderiv_hessian <- function(neg_lp, theta_mode) {
  fr <- .pskf_freeze_open()
  on.exit(.pskf_freeze_close(fr), add = TRUE)
  neg_lp <- .pskf_freeze_wrap(neg_lp)
  f0 <- tryCatch(neg_lp(theta_mode),
                 error = function(e) .dynhr_reraise_bug(e, NA_real_))
  if (!is.finite(f0)) return(NULL)
  neg_lp_checked <- function(theta) {
    v <- tryCatch(neg_lp(theta),
                  error = function(e) .dynhr_reraise_bug(e, NA_real_))
    if (!is.finite(v))
      stop("non-finite log-posterior encountered during Hessian evaluation")
    v
  }
  tryCatch(
    numDeriv::hessian(neg_lp_checked, as.numeric(theta_mode)),
    error = function(e) .dynhr_reraise_bug(e, NULL)
  )
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
.proposal_cov_full <- function(lp_fn, theta_mode, prior_spec, verbose = TRUE,
                               grad_fn = NULL) {

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

  if (verbose) .dynhr_cat("  Building full-Hessian proposal covariance...\n")

  H <- NULL
  bnd_eta <- NULL
  ## Exact gradient available: central differences of the gradient,
  ## 2n calls, instead of numDeriv's Richardson stencil on the log-posterior.
  ## The bound-active coordinates (same reach as the bound path below) are
  ## differenced one-sided inward from their own gradient component and
  ## decoupled. A non-finite interior entry leaves H NULL: the numDeriv path
  ## below then runs exactly as without a gradient.
  if (.step6_grad_usable(grad_fn)) {
    h_c    <- pmax(1e-4 * abs(theta_mode), 1e-5)
    bnds   <- .step6_bounds(prior_spec, pnames)
    at_bnd <- .step6_bound_active(theta_mode, bnds$lower, bnds$upper, h_c)
    coords <- which(!at_bnd)
    if (verbose)
      .dynhr_cat("    (central differences of the analytic gradient)\n")
    gh <- grad_hessian(grad_fn, theta_mode, coords = coords, at_bnd = at_bnd,
                       bnds = bnds)
    if (all(is.finite(gh$hess))) {
      if (any(at_bnd)) {
        bh <- .step6_bound_hessian(gh$hess, lp_fn, theta_mode, at_bnd, h_c, bnds,
                                   prior_spec, verbose = verbose, os = gh$os)
        H <- -bh$hess
        bnd_eta <- bh
      } else {
        H <- -gh$hess
      }
    } else if (verbose) {
      .dynhr_cat("    (gradient stencil non-finite; numDeriv on the log-posterior)\n")
    }
  }

  if (is.null(H) && requireNamespace("numDeriv", quietly = TRUE)) {
    H <- .proposal_numderiv_hessian(neg_lp, theta_mode)
  }

  if (is.null(H)) {
    ## Mode within the central stencil of a prior bound (a KKT-active
    ## constrained mode): the plain central stencil would step outside the
    ## support and return non-finite entries. Use the bound-aware Step-6
    ## Hessian instead (interior block central, bound coordinates one-sided
    ## inward and decoupled -- see .step6_bound_hessian()).
    h_c    <- pmax(1e-4 * abs(theta_mode), 1e-5)
    bnds   <- .step6_bounds(prior_spec, pnames)
    at_bnd <- .step6_bound_active(theta_mode, bnds$lower, bnds$upper, h_c)
    if (any(at_bnd)) {
      coords <- which(!at_bnd)
      H <- matrix(0, n_par, n_par)
      if (length(coords)) {
        th_sub_fn <- function(x) { th <- theta_mode; th[coords] <- x; neg_lp(th) }
        H[coords, coords] <- .neg_hessian_central(th_sub_fn, theta_mode[coords])
      }
      ## .step6_bound_hessian works on the LOG-posterior Hessian; H here is
      ## the NEGATIVE log-posterior Hessian.
      bh <- .step6_bound_hessian(-H, lp_fn, theta_mode, at_bnd, h_c, bnds,
                                 prior_spec, verbose = verbose)
      H <- -bh$hess
      bnd_eta <- bh
    } else {
      H <- .neg_hessian_central(neg_lp, theta_mode)
    }
  }

  ## Bound-active mode: carry the eta-space proposal covariance (the space
  ## the transform_params samplers use) as attr "Sigma_eta" -- unscaled,
  ## like this method's theta-space result. NULL (no attribute) otherwise.
  .attach_eta <- function(Sigma) {
    if (!is.null(bnd_eta))
      attr(Sigma, "Sigma_eta") <- .step6_eta_cov(Sigma, bnd_eta, theta_mode,
                                                 prior_spec, scale = 1)
    Sigma
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
          .dynhr_cat(sprintf("    Hessian condition number (repaired): %.3e\n", kappa))
          .dynhr_cat(sprintf("    Eigenvalues floored: %d/%d\n", n_floored, n_par))
          .dynhr_cat("\n    Proposal std dev:\n")
          prop_sd <- sqrt(diag(Sigma))
          for (i in seq_along(pnames)) {
            .dynhr_cat(sprintf("      %12s: mode=%10.6f  prop_sd=%10.6f\n",
                        pnames[i], theta_mode[i], prop_sd[i]))
          }
        }
        return(.attach_eta(Sigma))
      }
      ## chol failed despite spectral repair -- fall through to .make_pd
    }
  }

  ## Fallback path: Hessian has Inf/NaN entries (stencil crossed a hard wall),
  ## or is indefinite with no PD spectral repair.  Apply .make_pd(), which
  ## preserves finite local curvature and assigns prior-var to bad coordinates.
  ## This replaces the old "return NULL -> prior-only diagonal" behaviour.
  has_nonfinite <- any(!is.finite(H))
  .dynhr_warn(
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
    .dynhr_cat(sprintf("    .make_pd: %d/%d coordinate(s) fell back to prior variance.\n",
                n_bad, n_par))
    .dynhr_cat("\n    Proposal std dev (regularised):\n")
    prop_sd <- sqrt(diag(Sigma))
    for (i in seq_along(pnames)) {
      .dynhr_cat(sprintf("      %12s: mode=%10.6f  prop_sd=%10.6f\n",
                  pnames[i], theta_mode[i], prop_sd[i]))
    }
  }

  .attach_eta(Sigma)
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


#' Sanity-check an FD Hessian's condition number against a reference Hessian
#'
#' \code{numDeriv::hessian}'s DEFAULT step overshoots near-unit-root /
#' determinacy boundaries and inflates the reported posterior-Hessian
#' condition number by 5-7 orders of magnitude relative to an analytic or
#' BFGS-curvature reference (examples: nk_small FD 1.1e11 vs
#' analytic 3.6e5; sw2007-stress FD 2e16 vs analytic 2.4e8). This helper
#' compares the two condition numbers and warns when the FD one is
#' implausibly larger, so the trap is caught instead of silently poisoning a
#' proposal covariance or a Laplace marginal-likelihood calculation.
#'
#' Both matrices are symmetrised (\code{(H + t(H)) / 2}) before their
#' condition numbers are computed, since finite-difference Hessians are only
#' numerically (not exactly) symmetric. The condition number is the ratio of
#' the largest to the smallest eigenvalue magnitude, \code{|lambda_max| /
#' |lambda_min|}; a zero or non-finite denominator returns \code{Inf}
#' (documented, not an error -- a singular/degenerate Hessian is itself
#' informative and should not abort the check).
#'
#' @param H_fd   Finite-difference Hessian (e.g. from \code{numDeriv::hessian}
#'   or \code{num_hessian}/\code{num_hessian_mirai}); square numeric matrix.
#' @param H_ref  Reference Hessian to compare against -- typically an
#'   analytic Hessian (e.g. \code{posterior_hessian}) or a BFGS/quasi-Newton
#'   curvature estimate (e.g. \code{mode_res$H_bfgs}, appropriately inverted/
#'   signed); square numeric matrix, same dimension as \code{H_fd}.
#' @param ratio_tol  Warn when \code{kappa_fd / kappa_ref > ratio_tol}
#'   (default 1e3; the paper's examples show ratios of 1e5-1e8, so 1e3 is a
#'   conservative trigger well below the pathological cases while still
#'   quiet on ordinary numerical noise).
#' @param label_fd,label_ref  Labels used in the warning message to identify
#'   which matrix is which (default \code{"FD"} / \code{"reference"}).
#'
#' @return Invisibly, a list with elements \code{kappa_fd}, \code{kappa_ref},
#'   \code{ratio} (\code{kappa_fd / kappa_ref}), and \code{flagged} (logical,
#'   \code{TRUE} when \code{ratio > ratio_tol}).
#'
#' @export
check_hessian_conditioning <- function(H_fd, H_ref, ratio_tol = 1e3,
                                        label_fd = "FD",
                                        label_ref = "reference") {
  .check_square_numeric <- function(H, argname) {
    if (!is.matrix(H) || !is.numeric(H))
      stop(sprintf("check_hessian_conditioning: '%s' must be a numeric matrix.", argname))
    if (nrow(H) != ncol(H))
      stop(sprintf("check_hessian_conditioning: '%s' must be square.", argname))
  }
  .check_square_numeric(H_fd, "H_fd")
  .check_square_numeric(H_ref, "H_ref")

  if (nrow(H_fd) != nrow(H_ref))
    stop("check_hessian_conditioning: 'H_fd' and 'H_ref' must have the same dimension.")

  if (any(!is.finite(H_fd)))
    stop("check_hessian_conditioning: 'H_fd' contains non-finite entries.")
  if (any(!is.finite(H_ref)))
    stop("check_hessian_conditioning: 'H_ref' contains non-finite entries.")

  .kappa <- function(H) {
    Hs  <- (H + t(H)) / 2
    ev  <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
    mag <- abs(ev)
    lam_max <- max(mag)
    lam_min <- min(mag)
    if (!is.finite(lam_max) || !is.finite(lam_min) || lam_min <= 0) return(Inf)
    lam_max / lam_min
  }

  kappa_fd  <- .kappa(H_fd)
  kappa_ref <- .kappa(H_ref)

  ratio <- if (is.finite(kappa_fd) && is.finite(kappa_ref) && kappa_ref > 0) {
    kappa_fd / kappa_ref
  } else if (is.finite(kappa_ref) && !is.finite(kappa_fd)) {
    Inf
  } else {
    NA_real_
  }

  flagged <- isTRUE(is.finite(ratio) && ratio > ratio_tol) ||
    isTRUE(!is.finite(kappa_fd) && is.finite(kappa_ref))

  if (flagged) {
    .dynhr_warn(sprintf(
      paste0(
        "check_hessian_conditioning: %s condition number (%.3e) vastly exceeds ",
        "the %s condition number (%.3e), ratio %.3e > ratio_tol %.3e. This usually ",
        "means the finite-difference step overshot a near-unit-root or determinacy ",
        "boundary, not that the posterior is genuinely that ill-conditioned. ",
        "Recommend: (1) use the analytic posterior_hessian() instead of numDeriv's ",
        "default-step FD Hessian; (2) if FD is unavoidable, use a smaller, ",
        "feasibility-aware step size; default numDeriv steps were found to inflate ",
        "kappa by 5-7 orders of magnitude on near-boundary DSGE models ",
        "(nk_small: FD 1.1e11 vs analytic 3.6e5; sw2007-stress: 2e16 vs 2.4e8)."
      ),
      label_fd, kappa_fd, label_ref, kappa_ref, ratio, ratio_tol
    ))
  }

  invisible(list(
    kappa_fd  = kappa_fd,
    kappa_ref = kappa_ref,
    ratio     = ratio,
    flagged   = flagged
  ))
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


## --------------------------------------------------------------------------
## Step-6 proposal covariance at a bound-constrained (KKT) mode
## --------------------------------------------------------------------------
##
## A mode can sit ON a prior bound (a KKT-valid constrained mode: the
## log-posterior still increases toward the bound, the support stops it). The
## central FD stencil of the Step-6 Hessian then steps outside the support,
## the prior returns -Inf, and the whole row/column of that coordinate became
## non-finite -> .make_pd assigned PRIOR variances (NZSIM at the Dynare mc5
## mode: five parameters 1e-8..1e-13 from their bounds).
##
## What the proposal SHOULD be there depends on the sampler space:
##
##  * theta-space samplers (transform_params = FALSE): the Laplace
##    approximation is valid only on the FACE of the box. The interior
##    coordinates get the inverse curvature of the log-posterior restricted
##    to that face (bound coordinates held at the mode); each bound
##    coordinate gets the inverse of its ONE-SIDED (inward) curvature,
##    decoupled from the rest. That is a local step scale, not a posterior
##    variance: the marginal there is a truncated, gradient-dominated
##    (non-Gaussian) density, and the RWMH proposals that cross the bound are
##    rejected by the prior. A non-concave inward direction falls back to
##    the prior variance for that coordinate.
##
##  * eta-space samplers (transform_params = TRUE, the default): bounds are
##    at infinity, so an FD Hessian in eta never crosses them -- but the
##    theta-mode maps to eta ~ log(distance to bound) -> -Inf, where the
##    Jacobian-adjusted eta target has slope ~1 and curvature ~0 (the
##    delta-method conversion D^-1 Sigma_theta D^-1 then divides by
##    (dtheta/deta)^2 ~ dist^2, e.g. 1e-20). The eta target's own mode along
##    that coordinate is at distance delta* from the bound solving
##    c delta^2 - g delta - 1 = 0 (g, c: inward slope and minus curvature of
##    the log-posterior), where its curvature is -(1 + c delta*^2) in
##    [-2, -1] for a log-type transform: an O(1) eta proposal variance. The
##    bound coordinates therefore get the inverse negative curvature of the
##    Jacobian-adjusted eta log-posterior AT ITS CONDITIONAL MODE (1-D search,
##    other coordinates at the mode), decoupled; the interior block keeps
##    the delta-method conversion samplers already apply.
##
## At an interior mode (no stencil crosses a bound) none of this runs and
## Step 6 is bit-identical to before.

#' Prior-support box of the Step-6 parameters
#'
#' The support the eta-space transform uses (\code{.resolve_param_support}:
#' explicit prior-spec lower/upper, else the distribution's support).
#' Parameters without a prior-spec row are unbounded.
#' @noRd
.step6_bounds <- function(prior_spec, pnames) {
  lo <- stats::setNames(rep(-Inf, length(pnames)), pnames)
  hi <- stats::setNames(rep( Inf, length(pnames)), pnames)
  if (!is.data.frame(prior_spec) || is.null(prior_spec$name))
    return(list(lower = lo, upper = hi))
  for (k in seq_along(pnames)) {
    idx <- match(pnames[k], prior_spec$name)
    if (is.na(idx)) next
    sup <- .resolve_param_support(prior_spec[idx, , drop = FALSE])
    lo[k] <- sup$a
    hi[k] <- sup$b
  }
  list(lower = lo, upper = hi)
}

#' Coordinates whose FD stencil leaves the prior-support box
#'
#' @param reach Per-coordinate maximal stencil displacement (num_hessian's
#'   diagonal reaches \code{2 h_i}; .neg_hessian_central's \code{h_i}).
#' @return Named logical: \code{TRUE} where \code{theta_i -/+ reach_i} falls
#'   strictly outside \code{[lower_i, upper_i]} (the bound is active at the
#'   FD scale).
#' @noRd
.step6_bound_active <- function(theta, lower, upper, reach) {
  out <- as.logical((theta - reach < lower) | (theta + reach > upper))
  out[is.na(out)] <- FALSE
  names(out) <- names(theta)
  out
}

#' One-sided (inward) curvature and slope of a scalar function
#'
#' Steps away from the NEARER bound. With \code{s} the inward direction,
#' \code{f1 = f(theta + s h e_i)}, \code{f2 = f(theta + 2 s h e_i)}:
#' \code{curv = (f2 - 2 f1 + f0) / h^2} (first-order accurate) and
#' \code{slope = (4 f1 - f2 - 3 f0) / (2 h)}, the derivative along \code{s}
#' (second-order accurate). The step is shrunk only if \code{2 h} would
#' reach the FAR bound. A PSKF log-posterior is differenced on the pruning
#' selection made at \code{theta} (evaluated first, so it records; see
#' .pskf_freeze_open); other functions are unaffected.
#' @noRd
.step6_one_sided <- function(f, theta, i, h, lower, upper) {
  iw  <- .step6_inward(theta, i, h, lower, upper)
  s   <- iw$dir; h <- iw$h
  fr  <- .pskf_freeze_open()
  on.exit(.pskf_freeze_close(fr), add = TRUE)
  f   <- .pskf_freeze_wrap(f)
  f0  <- f(theta)
  th1 <- theta; th1[i] <- theta[[i]] + s * h
  th2 <- theta; th2[i] <- theta[[i]] + 2 * s * h
  f1  <- f(th1)
  f2  <- f(th2)
  list(curv  = (f2 - 2 * f1 + f0) / h^2,
       slope = (4 * f1 - f2 - 3 * f0) / (2 * h),
       dir   = s, h = h,
       bound = iw$bound)
}

#' Inward direction and step of a one-sided stencil at a bound coordinate
#'
#' Steps away from the NEARER bound; the step \code{h} is shrunk only if the
#' stencil's reach \code{2 h} would reach the FAR bound. Shared by the
#' log-posterior (\code{.step6_one_sided}) and gradient
#' (\code{.step6_grad_plan}) one-sided stencils.
#' @return list(dir = +1 / -1, h, bound = the nearer bound).
#' @noRd
.step6_inward <- function(theta, i, h, lower, upper) {
  room_up <- upper[[i]] - theta[[i]]
  room_dn <- theta[[i]] - lower[[i]]
  s    <- if (room_up >= room_dn) 1 else -1
  room <- max(room_up, room_dn)
  if (!(2 * h < room)) h <- room / 2.5
  list(dir = s, h = h, bound = if (s < 0) upper[[i]] else lower[[i]])
}

#' Bound-aware repair of a Step-6 log-posterior Hessian
#'
#' Decouples the bound-active coordinates (rows/columns zeroed) and sets
#' their diagonal to the one-sided inward curvature (or, with
#' \code{keep_diag = TRUE}, keeps a supplied finite analytic diagonal); a
#' non-concave or non-finite inward curvature falls back to
#' \code{-1 / prior_var}. See the block comment above for the rationale.
#'
#' @param hess LOG-posterior Hessian (n x n) whose interior block is valid.
#' @param h_vec Per-coordinate FD step.
#' @param os Optional precomputed one-sided results for the bound-active
#'   coordinates (a list named like \code{theta_mode}, as
#'   \code{.step6_grad_finish()} returns them from the gradient). NULL (the
#'   default): computed here from the log-posterior (\code{.step6_one_sided}).
#' @return list(hess, os (per-coordinate one-sided results), at_bnd, bnds, f)
#' @noRd
.step6_bound_hessian <- function(hess, lp_fn, theta_mode, at_bnd, h_vec, bnds,
                                 prior_spec, keep_diag = FALSE,
                                 verbose = TRUE, os = NULL) {
  pnames <- names(theta_mode)
  f <- function(th) {
    names(th) <- pnames
    r <- lp_fn(th)
    if (is.list(r)) r$logpost else r
  }
  pv <- if (is.data.frame(prior_spec) && !is.null(prior_spec$std))
    prior_spec$std[match(pnames, prior_spec$name)]^2 else rep(NA_real_, length(pnames))
  pv[!is.finite(pv) | pv <= 0] <- 1

  B <- which(at_bnd)
  d_in <- diag(hess)[B]
  hess[B, ] <- 0
  hess[, B] <- 0
  os_in <- os
  os <- stats::setNames(vector("list", length(pnames)), pnames)
  for (k in seq_along(B)) {
    b <- B[k]
    o <- if (!is.null(os_in) && !is.null(os_in[[b]])) os_in[[b]] else
      .step6_one_sided(f, theta_mode, b, h_vec[[b]], bnds$lower, bnds$upper)
    os[[b]] <- o
    d <- if (isTRUE(keep_diag) && is.finite(d_in[k])) d_in[k] else o$curv
    if (!(is.finite(d) && d < 0)) d <- -1 / pv[b]
    hess[b, b] <- d
  }

  if (verbose) {
    .dynhr_cat(sprintf(
      "  %d parameter(s) at a prior bound (KKT-active at the FD scale): one-sided inward curvature, decoupled.\n",
      length(B)))
    for (b in B)
      .dynhr_cat(sprintf("      %12s: mode=%.10g  bound=%.10g  distance=%.3e\n",
                         pnames[b], theta_mode[[b]], os[[b]]$bound,
                         abs(theta_mode[[b]] - os[[b]]$bound)))
  }
  list(hess = hess, os = os, at_bnd = at_bnd, bnds = bnds, f = f)
}

#' Eta-space variance of one bound-active coordinate
#'
#' Inverse negative curvature of the Jacobian-adjusted eta log-posterior at
#' its conditional mode along coordinate \code{i} (other coordinates held at
#' \code{theta_mode}). The mode is bracketed around the closed-form
#' \code{delta* = (g + sqrt(g^2 + 4c)) / (2c)} (distance from the bound;
#' \code{g}, \code{c}: inward slope and minus curvature from
#' \code{.step6_one_sided}) and located by \code{optimize()}; the curvature
#' is a central difference with an eta step of 1e-3 (the eta scale there is
#' O(1)). Returns 1 (the asymptotic unit curvature) if the curvature is not
#' negative and finite.
#' @noRd
.step6_eta_bound_var <- function(f, theta_mode, i, tr, os, bnds) {
  eta0 <- tr$to_unconstrained(theta_mode)
  fe <- function(e) {
    eta <- eta0; eta[i] <- e
    th <- theta_mode; th[i] <- tr$to_constrained(eta)[[i]]
    v <- f(th) + tr$log_jacobian(eta)
    if (is.finite(v)) v else -.Machine$double.xmax
  }
  width <- bnds$upper[[i]] - bnds$lower[[i]]
  c_in  <- -os$curv
  g_in  <- os$slope
  dstar <- if (is.finite(c_in) && c_in > 0 && is.finite(g_in)) {
    (g_in + sqrt(g_in^2 + 4 * c_in)) / (2 * c_in)
  } else if (is.finite(g_in) && g_in < 0) {
    -1 / g_in
  } else {
    0.1 * min(width, 1)
  }
  dstar <- min(dstar, 0.5 * width)
  eta_at <- function(d) {
    th <- theta_mode; th[i] <- os$bound + os$dir * d
    tr$to_unconstrained(th)[[i]]
  }
  lo_d <- dstar / 30
  hi_d <- min(dstar * 30, 0.99 * width)
  for (pass in 1:2) {
    iv  <- sort(c(eta_at(lo_d), eta_at(hi_d)))
    opt <- stats::optimize(fe, iv, maximum = TRUE, tol = 1e-8)
    if (min(abs(opt$maximum - iv)) > 1e-4 * diff(iv)) break
    lo_d <- lo_d / 30
    hi_d <- min(hi_d * 30, 0.99 * width)
  }
  he <- 1e-3
  e  <- opt$maximum
  ## PSKF: the three points share the selection made at e (evaluated first)
  fr <- .pskf_freeze_open()
  on.exit(.pskf_freeze_close(fr), add = TRUE)
  fz <- .pskf_freeze_wrap(fe)
  f_e  <- fz(e)
  curv <- (fz(e + he) - 2 * f_e + fz(e - he)) / he^2
  if (is.finite(curv) && curv < 0) -1 / curv else 1
}

#' Eta-space proposal covariance at a bound-active mode
#'
#' Interior block: the delta-method conversion samplers already apply
#' (\code{.cov_theta_to_eta} at \code{theta_mode}); bound-active
#' coordinates: \code{scale * .step6_eta_bound_var()}, decoupled.
#' @param Sigma Theta-space proposal covariance (block-diagonal across the
#'   interior / bound split).
#' @param bh Result of \code{.step6_bound_hessian()}.
#' @param scale The factor the theta-space result carries (2.38^2 / n for
#'   build_sigma_prop, 1 for proposal_cov(method = "full")).
#' @noRd
.step6_eta_cov <- function(Sigma, bh, theta_mode, prior_spec, scale) {
  pnames <- names(theta_mode)
  tr <- build_param_transform(prior_spec, pnames)
  Sig_eta <- .cov_theta_to_eta(Sigma, tr, theta_mode)
  for (b in which(bh$at_bnd)) {
    v <- .step6_eta_bound_var(bh$f, theta_mode, b, tr, bh$os[[b]], bh$bnds)
    Sig_eta[b, ] <- 0
    Sig_eta[, b] <- 0
    Sig_eta[b, b] <- scale * v
  }
  dimnames(Sig_eta) <- list(pnames, pnames)
  Sig_eta
}


## --------------------------------------------------------------------------
## Step-6 Hessian from the exact gradient
## --------------------------------------------------------------------------
##
## num_hessian() differences the LOG-POSTERIOR over every (i, j) pair:
## 2 n (n + 1) evaluations (NZSIM, n = 68: 9,384; SW2007, n = 36: 2,664).
## When the mode stage has an analytic gradient (make_posterior_grad(), the
## Step-5 grad_fn), the Hessian is instead the symmetrised central difference
## of that gradient,
##     H[, j] = (g(theta + h_j e_j) - g(theta - h_j e_j)) / (2 h_j),
##     H <- (H + t(H)) / 2,
## 2 n gradient calls (SW2007 under load_all: 1.2 s vs 11.7 s). It is also
## the MORE accurate of the two: the lp stencil loses ~eps |f| / h^2 to
## roundoff (|f| ~ 900..9000), the gradient stencil only eps |g| / h.
##
## Step h_j = 1e-5 * max(1, |theta_j|). Step sweep on SW2007 at the benchmark
## mode (adjoint_solution gradient; error = relative Frobenius distance from
## a Richardson extrapolation of the h = 1e-5 / 3e-6 pair):
##     h = 1e-3: 9.4e-4   1e-4: 9.4e-6   3e-5: 8.4e-7   1e-5: 9.4e-8
##     h = 1e-6: 4.2e-8   1e-7: 2.8e-7   1e-8: 4.6e-6
##     num_hessian (lp stencil, h = 1e-4): 5.7e-4
## i.e. truncation O(h^2) down to h ~ 1e-6, roundoff below; 1e-5 sits a
## decade above the roundoff floor (margin for a noisier gradient).
##
## Not used for a "hybrid" gradient (grad_method attribute): its likelihood
## score is itself a forward difference of the log-posterior, so 2 n gradient
## calls cost ~ as many log-posterior evaluations as num_hessian and the
## nested difference has the worse roundoff. A gradient stencil with a
## non-finite entry falls back to num_hessian.
##
## Bound-active coordinates (above): unchanged selection (the lp
## stencil's reach decides, so the same parameters are decoupled); their
## inward curvature is the one-sided difference of their own gradient
## component, phi''(0) = s (4 g1 - g2 - 3 g0) / (2 h) with g_k the component
## at theta + k s h e_b, and the inward slope s g0 -- both exact-gradient
## counterparts of .step6_one_sided's. The interior coordinates are
## differenced with the bound ones held at the mode, as before.

.step6_grad_h <- 1e-5

#' Can this gradient closure supply the Step-6 Hessian?
#' @noRd
.step6_grad_usable <- function(grad_fn) {
  is.function(grad_fn) && !identical(attr(grad_fn, "grad_method"), "hybrid")
}

#' A gradient evaluation as a plain numeric vector in theta's order
#'
#' NA (length n) when the closure returns anything but n values.
#' @noRd
.step6_grad_eval <- function(grad_fn, theta) {
  v  <- grad_fn(theta)
  nm <- names(v)
  if (!is.null(nm) && !is.null(names(theta)) && all(names(theta) %in% nm))
    v <- v[names(theta)]
  v <- as.numeric(v)
  if (length(v) != length(theta)) rep(NA_real_, length(theta)) else v
}

#' The gradient evaluation points of the Step-6 gradient stencil
#'
#' @param coords Interior coordinates (central stencil).
#' @param at_bnd Logical, bound-active coordinates (one-sided stencil).
#' @param h_g Per-coordinate step.
#' @return list(points = list of theta vectors, coords, B, h_g, inward =
#'   per-bound-coordinate \code{.step6_inward()} results, i0 = index of the
#'   unperturbed point (0 when no coordinate is bound-active)).
#' @noRd
.step6_grad_plan <- function(theta, coords, at_bnd, h_g, bnds) {
  pts <- vector("list", 0L)
  for (j in coords) {
    tp <- tm <- theta
    tp[j] <- theta[[j]] + h_g[[j]]
    tm[j] <- theta[[j]] - h_g[[j]]
    pts <- c(pts, list(tp, tm))
  }
  B  <- which(at_bnd)
  iw <- list()
  i0 <- 0L
  if (length(B)) {
    pts <- c(pts, list(theta))
    i0  <- length(pts)
    for (b in B) {
      w <- .step6_inward(theta, b, h_g[[b]], bnds$lower, bnds$upper)
      iw[[as.character(b)]] <- w
      t1 <- t2 <- theta
      t1[b] <- theta[[b]] + w$dir * w$h
      t2[b] <- theta[[b]] + 2 * w$dir * w$h
      pts <- c(pts, list(t1, t2))
    }
  }
  list(points = pts, coords = coords, B = B, h_g = h_g, inward = iw, i0 = i0,
       n = length(theta), pnames = names(theta))
}

#' Assemble the Step-6 Hessian from gradients at the plan's points
#'
#' @param G List of gradient vectors, one per \code{plan$points} entry.
#' @return list(hess = LOG-posterior Hessian, n x n: the interior block
#'   symmetrised, bound rows/columns 0; os = one-sided results for the
#'   bound-active coordinates, the \code{.step6_bound_hessian()} input).
#' @noRd
.step6_grad_finish <- function(plan, G) {
  n  <- plan$n
  cs <- plan$coords
  H  <- matrix(0, n, n)
  if (length(cs)) {
    cols <- matrix(NA_real_, n, length(cs))
    for (a in seq_along(cs)) {
      j <- cs[a]
      cols[, a] <- (G[[2L * a - 1L]] - G[[2L * a]]) / (2 * plan$h_g[[j]])
    }
    Bk <- cols[cs, , drop = FALSE]
    H[cs, cs] <- (Bk + t(Bk)) / 2
  }
  os <- NULL
  if (length(plan$B)) {
    os <- stats::setNames(vector("list", n), plan$pnames)
    g0 <- G[[plan$i0]]
    for (k in seq_along(plan$B)) {
      b  <- plan$B[k]
      w  <- plan$inward[[as.character(b)]]
      s  <- w$dir
      g1 <- G[[plan$i0 + 2L * k - 1L]][[b]]
      g2 <- G[[plan$i0 + 2L * k]][[b]]
      os[[b]] <- list(curv  = s * (4 * g1 - g2 - 3 * g0[[b]]) / (2 * w$h),
                      slope = s * g0[[b]],
                      dir   = s, h = w$h, bound = w$bound)
    }
  }
  list(hess = H, os = os)
}

#' Step-6 Hessian by central differences of an exact gradient (serial)
#'
#' See the block comment above. \code{coords}, as in \code{num_hessian}, are
#' the coordinates differenced; the others' rows/columns are 0.
#' @param grad_fn Theta-space gradient of the log-posterior
#'   (\code{make_posterior_grad()}).
#' @param h Base relative step (\code{h * max(1, |theta_j|)}).
#' @return list(hess, os) as \code{.step6_grad_finish()}.
#' @noRd
grad_hessian <- function(grad_fn, theta, h = .step6_grad_h,
                         coords = seq_along(theta),
                         at_bnd = rep(FALSE, length(theta)), bnds = NULL) {
  h_g  <- h * pmax(1, abs(theta))
  plan <- .step6_grad_plan(theta, coords, at_bnd, h_g, bnds)
  G    <- lapply(plan$points, function(th) .step6_grad_eval(grad_fn, th))
  .step6_grad_finish(plan, G)
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
