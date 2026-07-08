## R/diag-robust-inference.R
## --------------------------------------------------------------------------
## Identification-robust inference via Andrews & Mikusheva (2015) LM2 score
## test.
##
## Under weak or partial identification the Wald/Hessian interval is
## spuriously tight. The LM2 score test inverts a statistic that is
## asymptotically chi^2(k) WITHOUT assuming identification, because the
## per-period Kalman score {s_t} forms a martingale-difference sequence
## under H0 regardless of the strength of identification.
##
## Statistic:
##   LM2(theta0) = S_T(theta0)' Vhat(theta0)^{-1} S_T(theta0)
## where
##   S_T(theta0)  = sum_t s_t(theta0)  (full-sample likelihood score)
##   Vhat(theta0) = sum_t s_t s_t'     (OPG "meat" -- stays non-degenerate
##                                       even when the Hessian is rank-def.)
##
## Confidence set: { theta0 : LM2(theta0) <= qchisq(level, k) }
##
## Reference:
##   Andrews, I., & Mikusheva, A. (2015). Maximum likelihood inference in
##   weakly identified DSGE models. Quantitative Economics, 6, 123-152.
##   DOI: 10.3982/QE331
## --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Internal: Moore-Penrose pseudo-inverse with effective-rank estimation.
# ---------------------------------------------------------------------------
.robust_pinv <- function(M) {
  sv  <- svd(M)
  tol <- max(dim(M)) * max(sv$d) * .Machine$double.eps * 100
  dInv <- ifelse(sv$d > tol, 1 / sv$d, 0)
  rank <- sum(sv$d > tol)
  Minv <- sv$v %*% (dInv * t(sv$u))
  list(inv = Minv, rank = rank, d = sv$d)
}


# ---------------------------------------------------------------------------
# Core: LM2 statistic at a single theta0.
#
# Arguments:
#   theta0            Named numeric vector (the hypothesis to test).
#   S_T               Full-sample likelihood score (NOT posterior; prior
#                     contribution must already be stripped by the caller).
#   loglik_contrib_fn Function theta -> length-T numeric of per-period
#                     log-likelihood contributions (LIKELIHOOD only, no prior).
#   eps               Step scale for .d35_score_matrix FD.
#
# Returns: list(stat, df, Vhat, rank_Vhat, S_T_used)
#
# The per-period score matrix is built from loglik_contrib_fn (finite
# differences via .d35_score_matrix), then Vhat = crossprod(S_mat). If
# Vhat is singular a pseudo-inverse is used and df = rank(Vhat) < k.
# ---------------------------------------------------------------------------
.lm2_statistic <- function(theta0, S_T, loglik_contrib_fn, eps = 1e-4) {
  k <- length(S_T)

  # Per-period score matrix (T x k) by central FD of per-period log-lik.
  # .d35_score_matrix is in R/diag-deep-d35-softness.R (loaded before us).
  S_mat <- tryCatch(
    .d35_score_matrix(loglik_contrib_fn, theta0, eps = eps),
    error = function(e) NULL
  )

  if (is.null(S_mat) || !any(is.finite(S_mat))) {
    warning("LM2: per-period score matrix could not be computed; returning NA.")
    return(list(stat = NA_real_, df = k,
                Vhat = NULL, rank_Vhat = NA_integer_, S_T_used = S_T))
  }

  # Check for all-NA columns (params where FD failed).
  ok_cols <- which(apply(S_mat, 2, function(col) all(is.finite(col))))
  if (length(ok_cols) < k) {
    warning(sprintf("LM2: %d/%d score columns finite; using finite subset.",
                    length(ok_cols), k))
  }
  if (length(ok_cols) == 0) {
    warning("LM2: no finite score columns; returning NA.")
    return(list(stat = NA_real_, df = k,
                Vhat = NULL, rank_Vhat = 0L, S_T_used = S_T))
  }

  S_sub  <- S_mat[, ok_cols, drop = FALSE]
  S_T_sub <- S_T[ok_cols]
  Vhat   <- crossprod(S_sub)          # k_ok x k_ok OPG meat

  # Robust pseudo-inverse of Vhat.
  piv <- tryCatch(
    .robust_pinv(Vhat),
    error = function(e) {
      warning(sprintf("LM2: pseudo-inverse failed (%s); returning NA.", e$message))
      NULL
    }
  )
  if (is.null(piv)) {
    return(list(stat = NA_real_, df = length(ok_cols),
                Vhat = Vhat, rank_Vhat = NA_integer_, S_T_used = S_T))
  }

  stat <- as.numeric(S_T_sub %*% piv$inv %*% S_T_sub)
  if (!is.finite(stat)) stat <- NA_real_

  list(stat = stat, df = piv$rank, Vhat = Vhat,
       rank_Vhat = piv$rank, S_T_used = S_T)
}


# ---------------------------------------------------------------------------
#' Identification-robust confidence sets (Andrews-Mikusheva LM2)
#'
#' Inverts the LM2 score test of Andrews & Mikusheva (2015) to form
#' marginal confidence intervals that are valid under weak or partial
#' identification.  The intervals are computed by scanning a 1-D grid over
#' each requested parameter.  By default (\code{profile = TRUE}) the other
#' parameters are PROFILED OUT at each grid point (inner minimisation of the
#' statistic) so the result is the projection of the joint robust set onto the
#' parameter's axis -- this is what correctly returns an UNBOUNDED interval for
#' a jointly-flat direction (e.g. when only a sum of parameters is identified).
#' With \code{profile = FALSE} the other parameters are held at the mode (a
#' cheaper CONDITIONAL slice that is spuriously bounded for jointly-flat
#' directions, like Wald).
#'
#' The statistic \deqn{LM2(\theta_0) = S_T(\theta_0)' \hat{V}(\theta_0)^{-1}
#' S_T(\theta_0)} is asymptotically \eqn{\chi^2(k)} under \eqn{H_0: \theta =
#' \theta_0} regardless of identification strength, because the per-period
#' Kalman score is a martingale-difference sequence under H0.
#'
#' \strong{Score consistency requirement}: \code{S_T} (the gradient) must
#' equal \code{colSums} of the per-period score matrix.  When
#' \code{loglik_contrib_fn} is supplied together with the analytic
#' \code{grad_fn}, \code{robust_confidence_set} verifies this internally and
#' warns if the max absolute discrepancy exceeds \code{score_tol}.
#'
#' @param theta_mode Named numeric vector -- the parameter mode (point estimate).
#' @param loglik_contrib_fn Function \code{theta -> } length-T numeric of
#'   per-period \strong{log-likelihood} contributions (NOT posterior; no prior
#'   term).  This is the same interface as for D35.
#' @param grad_fn Optional function \code{theta -> } length-k numeric,
#'   the analytic \strong{likelihood} gradient \eqn{S_T(\theta)}.  If
#'   \code{NULL}, \eqn{S_T} is approximated by \code{colSums} of the
#'   per-period score matrix (consistent but slightly less efficient).
#' @param params Character vector of parameter names to form CIs for.
#'   Defaults to \code{names(theta_mode)}.
#' @param level Confidence level (default 0.95).
#' @param n_grid Number of grid points per parameter (default 50).
#' @param grid_width_sd Multiplier: grid spans mode +/- \code{grid_width_sd}
#'   nominal Hessian SEs.  Default 5.  Ignored when \code{grid_bounds} is
#'   supplied.
#' @param grid_bounds Named list of \code{c(lo, hi)} per parameter.  If
#'   supplied, overrides \code{grid_width_sd} for the named parameters.
#' @param hess_se Named numeric vector of nominal (Hessian-based) SEs.
#'   Used only to set grid bounds when \code{grid_bounds} is NULL.  If
#'   NULL and \code{grid_bounds} is also NULL, a default grid width of 0.2
#'   (absolute) is used.
#' @param eps Step scale for the per-period score matrix FD (default 1e-4).
#' @param score_tol Tolerance for the score-consistency check (default 1e-4).
#' @param profile Logical (default \code{TRUE}). When \code{TRUE}, profile the
#'   other parameters out at each grid point (projection of the joint robust
#'   set -- the honest weak-ID-robust marginal, correctly UNBOUNDED for a
#'   jointly-flat direction). When \code{FALSE}, hold the others at the mode (a
#'   cheaper conditional slice, spuriously bounded under joint flatness).
#'   Profiling costs an inner optimisation per grid point, so on large models
#'   restrict \code{params} to the weakly-identified subset.
#'   A warning is issued if \code{max(|colSums(S_mat) - S_T|)} exceeds this.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{$ci}}{data.frame(param, lo, hi, bounded, wald_lo, wald_hi,
#'     wald_width, robust_width) -- the robust marginal interval.}
#'   \item{\code{$grid}}{list, one element per param: data.frame(value, lm2,
#'     in_ci) for plotting the LM2 profile.}
#'   \item{\code{$critical_value}}{The chi^2 critical value used.}
#'   \item{\code{$lm2_at_mode}}{LM2 at the mode (should be near 0).}
#'   \item{\code{$score_consistency}}{Max abs discrepancy between analytic S_T
#'     and colSums(S_mat) at the mode.}
#' }
#'
#' @references Andrews, I., & Mikusheva, A. (2015). Maximum likelihood
#'   inference in weakly identified DSGE models. Quantitative Economics,
#'   6, 123-152.
#' @export
robust_confidence_set <- function(theta_mode,
                                  loglik_contrib_fn,
                                  grad_fn       = NULL,
                                  params        = NULL,
                                  level         = 0.95,
                                  n_grid        = 50L,
                                  grid_width_sd = 5,
                                  grid_bounds   = NULL,
                                  hess_se       = NULL,
                                  eps           = 1e-4,
                                  score_tol     = 1e-4,
                                  profile       = TRUE) {

  if (is.null(params)) params <- names(theta_mode)
  if (is.null(params)) stop("robust_confidence_set: theta_mode must be named or params must be supplied.")

  ## Fail-loud validation of the per-period loglik closure (pathological-DSGE
  ## paper gap #1a: a NULL field access -- e.g. `$ll_contrib` instead of
  ## kalman_filter()'s actual `$loglik_contrib` -- yields a length-0 vector
  ## that previously flowed through silently and produced a degenerate
  ## [-Inf, Inf] confidence set).
  ll0 <- tryCatch(loglik_contrib_fn(theta_mode), error = function(e)
    stop("robust_confidence_set: loglik_contrib_fn(theta_mode) errored: ",
         conditionMessage(e), call. = FALSE))
  if (!is.numeric(ll0) || length(ll0) < 2L)
    stop("robust_confidence_set: loglik_contrib_fn(theta) must return the ",
         "length-T numeric vector of PER-PERIOD log-likelihood contributions; ",
         "got ", if (is.null(ll0)) "NULL" else paste0(class(ll0)[1L],
         " of length ", length(ll0)), ". (Common footgun: the ",
         "kalman_filter(return_ll_contrib = TRUE) field is `loglik_contrib` -- ",
         "accessing a misspelled field like `ll_contrib` returns NULL.)",
         call. = FALSE)
  if (!any(is.finite(ll0)))
    stop("robust_confidence_set: loglik_contrib_fn(theta_mode) returned no ",
         "finite values.", call. = FALSE)
  if (stats::sd(ll0[is.finite(ll0)]) == 0)
    stop("robust_confidence_set: loglik_contrib_fn(theta_mode) returned an ",
         "all-equal vector -- this is not a per-period contribution vector ",
         "(a recycled scalar total?).", call. = FALSE)

  k    <- length(theta_mode)
  crit <- stats::qchisq(level, df = k)   # chi^2(k) critical value

  # ------------------------------------------------------------------
  # Step 1: build the per-period score matrix and full-sample score
  # at the mode, then check score consistency.
  # ------------------------------------------------------------------
  S_mat_mode <- tryCatch(
    .d35_score_matrix(loglik_contrib_fn, theta_mode, eps = eps),
    error = function(e) NULL
  )
  if (is.null(S_mat_mode)) stop("robust_confidence_set: per-period score matrix failed at mode.")

  S_T_fd <- colSums(S_mat_mode)   # S_T by summing per-period contributions

  if (!is.null(grad_fn)) {
    S_T_analytic <- tryCatch(grad_fn(theta_mode), error = function(e) NULL)
    if (!is.null(S_T_analytic) && length(S_T_analytic) == k) {
      consistency_gap <- max(abs(S_T_analytic - S_T_fd))
      if (consistency_gap > score_tol) {
        warning(sprintf(
          "robust_confidence_set: score inconsistency max|S_T_analytic - colSums(S_mat)| = %.2e (> tol %.2e). LM2 computed with colSums. Reconcile prior handling.",
          consistency_gap, score_tol))
      }
      # Use the analytic gradient as the authoritative S_T
      S_T_mode <- S_T_analytic
    } else {
      consistency_gap <- NA_real_
      S_T_mode <- S_T_fd
    }
  } else {
    consistency_gap <- NA_real_
    S_T_mode <- S_T_fd
  }

  # ------------------------------------------------------------------
  # Step 2: LM2 at the mode (should be near 0 -- sanity check).
  # ------------------------------------------------------------------
  lm2_at_mode_res <- .lm2_statistic(theta_mode, S_T_mode, loglik_contrib_fn, eps = eps)
  lm2_at_mode     <- lm2_at_mode_res$stat

  # ------------------------------------------------------------------
  # Step 3: scan 1-D grid for each requested parameter.
  # ------------------------------------------------------------------
  ci_rows   <- vector("list", length(params))
  grid_list <- vector("list", length(params))
  names(ci_rows) <- params
  names(grid_list) <- params

  for (j in seq_along(params)) {
    nm <- params[j]
    if (!(nm %in% names(theta_mode))) {
      warning(sprintf("robust_confidence_set: parameter '%s' not in theta_mode; skipping.", nm))
      next
    }
    idx <- which(names(theta_mode) == nm)
    mode_j <- theta_mode[idx]

    # Determine grid bounds for this parameter.
    if (!is.null(grid_bounds) && nm %in% names(grid_bounds)) {
      lo_j <- grid_bounds[[nm]][1]
      hi_j <- grid_bounds[[nm]][2]
    } else if (!is.null(hess_se) && nm %in% names(hess_se)) {
      se_j <- hess_se[[nm]]
      lo_j <- mode_j - grid_width_sd * se_j
      hi_j <- mode_j + grid_width_sd * se_j
    } else {
      # Fallback: +/- 0.3 around the mode (absolute).
      lo_j <- mode_j - 0.3
      hi_j <- mode_j + 0.3
    }
    grid_j <- seq(lo_j, hi_j, length.out = n_grid)

    # LM2 at an arbitrary point (full-sample score + OPG meat).
    lm2_at <- function(theta0) {
      if (!is.null(grad_fn)) {
        S_T_j <- tryCatch(grad_fn(theta0), error = function(e) NULL)
      } else {
        S_T_j <- NULL
      }
      if (is.null(S_T_j) || length(S_T_j) != k) {
        S_mat_j <- tryCatch(
          .d35_score_matrix(loglik_contrib_fn, theta0, eps = eps),
          error = function(e) NULL
        )
        S_T_j <- if (!is.null(S_mat_j)) colSums(S_mat_j) else rep(NA_real_, k)
      }
      res <- .lm2_statistic(theta0, S_T_j, loglik_contrib_fn, eps = eps)
      if (is.null(res) || !is.finite(res$stat)) NA_real_ else res$stat
    }

    # Scan the 1-D grid for parameter `nm`.
    #  * profile = FALSE: the SLICE -- hold all other params at the mode. This is
    #    a CONDITIONAL interval; for a jointly-flat direction (e.g. only psi1+psi2
    #    is identified) the slice is spuriously bounded, exactly like Wald.
    #  * profile = TRUE (default): the PROJECTION of the joint robust set onto
    #    the `nm` axis -- at each grid value minimise LM2 over the OTHER params.
    #    For psi1+psi2 this finds psi2 = c - psi1 (keeping the identified sum at
    #    its mode), so LM2 ~ 0 everywhere and the robust set is (correctly)
    #    unbounded. This is the honest weak-ID-robust marginal set. Cost: an
    #    inner optimisation per grid point (expensive on large models -- restrict
    #    `params` to the weakly-identified subset, e.g. those flagged by D38).
    free_idx <- setdiff(seq_len(k), idx)
    do_profile <- isTRUE(profile) && length(free_idx) >= 1L
    stats_j <- vapply(grid_j, function(v0) {
      theta0      <- theta_mode
      theta0[idx] <- v0
      if (!do_profile) return(lm2_at(theta0))
      obj <- function(fv) {
        th <- theta0; th[free_idx] <- fv
        # suppressWarnings: the inner search may probe invalid regions (e.g. a
        # negative variance), which produce NaN likelihoods/scores. Those are
        # penalised to 1e10 below, so the optimiser steers away -- the warnings
        # are expected noise, not a failure.
        v  <- suppressWarnings(lm2_at(th))
        if (!is.finite(v)) 1e10 else v
      }
      start <- theta_mode[free_idx]
      mn <- if (length(free_idx) == 1L) {
        rng <- 10 * (hi_j - lo_j)
        op <- tryCatch(stats::optimize(obj, lower = start - rng, upper = start + rng),
                       error = function(e) NULL)
        if (is.null(op)) lm2_at(theta0) else op$objective
      } else {
        op <- tryCatch(stats::optim(start, obj, method = "Nelder-Mead",
                                    control = list(maxit = 200L)),
                       error = function(e) NULL)
        if (is.null(op)) lm2_at(theta0) else op$value
      }
      mn
    }, numeric(1))

    in_ci_j <- !is.na(stats_j) & stats_j <= crit
    vals_in  <- grid_j[in_ci_j]

    lo_ci <- if (length(vals_in) == 0) NA_real_ else min(vals_in)
    hi_ci <- if (length(vals_in) == 0) NA_real_ else max(vals_in)

    # Flag whether the CI appears bounded (not hitting the grid edges).
    bounded <- !is.na(lo_ci) && lo_ci > (lo_j + 0.05 * (hi_j - lo_j)) &&
               !is.na(hi_ci) && hi_ci < (hi_j - 0.05 * (hi_j - lo_j))

    # Wald CI for comparison.
    if (!is.null(hess_se) && nm %in% names(hess_se)) {
      se_wald <- hess_se[[nm]]
      z <- stats::qnorm((1 + level) / 2)
      wald_lo <- mode_j - z * se_wald
      wald_hi <- mode_j + z * se_wald
    } else {
      wald_lo <- NA_real_; wald_hi <- NA_real_
    }
    wald_width   <- if (!is.na(wald_lo)) wald_hi - wald_lo else NA_real_
    robust_width <- if (!is.na(lo_ci) && !is.na(hi_ci)) hi_ci - lo_ci else Inf

    ci_rows[[j]] <- data.frame(
      param        = nm,
      lo           = if (is.na(lo_ci)) -Inf else lo_ci,
      hi           = if (is.na(hi_ci))  Inf else hi_ci,
      bounded      = bounded,
      wald_lo      = wald_lo,
      wald_hi      = wald_hi,
      wald_width   = wald_width,
      robust_width = robust_width,
      stringsAsFactors = FALSE
    )

    grid_list[[nm]] <- data.frame(
      value = grid_j,
      lm2   = stats_j,
      in_ci = in_ci_j,
      stringsAsFactors = FALSE
    )
  }

  ci_df <- do.call(rbind, Filter(Negate(is.null), ci_rows))
  rownames(ci_df) <- NULL

  structure(
    list(
      ci                 = ci_df,
      grid               = grid_list,
      critical_value     = crit,
      df                 = k,
      level              = level,
      lm2_at_mode        = lm2_at_mode,
      score_consistency  = consistency_gap
    ),
    class = c("dynhr_robust_ci", "list")
  )
}


# ---------------------------------------------------------------------------
#' @export
print.dynhr_robust_ci <- function(x, ...) {
  cat(sprintf("Andrews-Mikusheva (2015) LM2 robust confidence set (level = %.2f)\n",
              x$level))
  cat(sprintf("  chi^2(%d) critical value: %.3f\n", x$df, x$critical_value))
  if (!is.null(x$lm2_at_mode) && !is.na(x$lm2_at_mode))
    cat(sprintf("  LM2 at mode: %.4f (should be near 0)\n", x$lm2_at_mode))
  if (!is.null(x$score_consistency) && !is.na(x$score_consistency))
    cat(sprintf("  Score consistency check: max|S_T - colSums(S_mat)| = %.2e\n",
                x$score_consistency))
  cat("\nMarginal robust intervals (1-D grid scan, nuisance params profiled):\n")
  ci <- x$ci
  for (i in seq_len(nrow(ci))) {
    lo <- ci$lo[i]; hi <- ci$hi[i]
    lo_s <- if (!is.finite(lo)) "-Inf" else sprintf("%.4f", lo)
    hi_s <- if (!is.finite(hi))  "Inf" else sprintf("%.4f", hi)
    bnd  <- if (!is.na(ci$bounded[i]) && ci$bounded[i]) "bounded" else "UNBOUNDED (or at grid edge)"
    wald <- if (!is.na(ci$wald_lo[i]))
      sprintf("  Wald [%.4f, %.4f]", ci$wald_lo[i], ci$wald_hi[i]) else ""
    cat(sprintf("  %-15s [%s, %s]  %s%s\n",
                ci$param[i], lo_s, hi_s, bnd, wald))
  }
  invisible(x)
}
