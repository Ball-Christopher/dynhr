## R/kalman-filter-student.R
## --------------------------------------------------------------------------
## Student-t innovation likelihood for DSGE Kalman filters.
##
## Motivation: fat-tailed observation/innovation distributions are a common
## extension in DSGE estimation (e.g. Chib, Ramamurthy & Shephard 2010;
## Canova & Ferroni 2011; 2025 JES fat-tail-DSGE survey). This file adds
## a multivariate Student-t log-likelihood layer on top of the standard
## Gaussian Kalman recursions: the state mean / covariance recursions are
## unchanged (Gaussian KF), only the per-period log-likelihood contribution
## is replaced by the log-density of a multivariate-t.
##
## Implemented estimator:
##   kalman_filter_student_t(Y, dr, model, params, obs_vars,
##                           student_df, me_variance = 0,
##                           lik_init = "auto")
##
## Per-period log-likelihood (k_t = number of non-NA observables at t):
##
##   The Gaussian filter produces innovation v_t ~ N(0, F_t) where F_t is the
##   k_t × k_t innovation covariance.  Under the Student-t observation model:
##
##   1. Scale the t-distribution so that its covariance equals F_t:
##        Sigma_t = c_nu * F_t,   c_nu = (nu - 2) / nu
##      (variance matching). This REQUIRES nu > 2: the t covariance is
##      nu/(nu-2) * Sigma_t and is infinite for nu <= 2, so no scale matches
##      F_t there and c_nu -> 0 as nu -> 2+. The filter therefore aborts with
##      class `dynhr_error_student_df` for nu <= 2 instead of switching to an
##      arbitrary scale (the former c_nu = 1 fallback made the likelihood jump
##      discontinuously at nu = 2).
##
##   2. Log-density of multivariate-t with location 0, scale Sigma_t, df nu:
##        loglik_t = lgamma((nu + k_t) / 2)
##                 - lgamma(nu / 2)
##                 - (k_t / 2) * log(nu * pi)
##                 - 0.5 * log|Sigma_t|
##                 - (nu + k_t) / 2 * log(1 + (1/nu) * v_t' Sigma_t^{-1} v_t)
##
##   Substituting Sigma_t = c_nu * F_t:
##        log|Sigma_t| = k_t * log(c_nu) + log|F_t|
##        v_t' Sigma_t^{-1} v_t = (1/c_nu) * v_t' F_t^{-1} v_t
##
##   The Gaussian limit (nu -> Inf) recovers the standard KF loglik exactly.
##   The lgamma difference is evaluated as lgamma(k/2) - lbeta(nu/2, k/2)
##   (R's lbeta is cancellation-free for a large argument) and log(c_nu) as
##   log1p(-2/nu), so a large nu (e.g. 1e9) reproduces kalman_filter() to
##   ~1e-9 per period instead of losing ~1e-6 per period to lgamma(5e8)
##   cancellation.
##
## UPGRADE NOTE: this version keeps the GAUSSIAN Kalman recursions and only
## changes the per-period log-likelihood.  The natural upgrade is a robust
## filter (Masreliez-Martin or scale-mixture reweighting) that adjusts the
## Kalman gain by the per-period t-weight (nu + k) / (nu + Q_t) where
## Q_t = v_t' F_t^{-1} v_t.  That is the exact filter under the
## conditionally-Gaussian scale-mixture representation.  The present
## version is the defensible first step: a consistent likelihood under the
## t-observation model, commonly used in the DSGE literature.
##
## The per-period covariance update is the shared .kf_step_core() of
## R/kf-step.R (brief 23 D1): me_variance is TRUE i.i.d. measurement error, so
## it enters F_t AND the Joseph term K me K' -- matching kalman_filter()
## (brief 23 A2: the former hand-rolled update dropped K me K').
##
## Missing observations: the same NA-reduction as kalman_filter() is used.
## --------------------------------------------------------------------------


#' Validate a Student-t likelihood request (shared by the filter and by
#' make_log_posterior(), which must fail at closure build rather than map
#' every draw's abort to -Inf).
#' @noRd
.student_t_validate <- function(nu, model, me_extra = NULL, shock_scale = NULL,
                                known_shocks = NULL, caller) {
  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 2)
    .dynhr_abort(caller, ": student_df must be a finite ",
                 "scalar > 2 (got ", format(nu), "). The t scale is matched ",
                 "to the innovation covariance, Sigma_t = (nu - 2)/nu * F_t, ",
                 "which needs a finite t variance (nu > 2).",
                 class = "dynhr_error_student_df")

  ## Unsupported per-period inputs: abort instead of silently ignoring them.
  .refuse_obs_trends(model, caller)
  unsupported <- c(
    me_extra     = !is.null(me_extra),
    shock_scale  = !is.null(shock_scale),
    known_shocks = !is.null(known_shocks),
    filter_tunes = .student_t_has_filter_tunes(model),
    heteroskedastic_shocks = .student_t_has_het_shocks(model))
  if (any(unsupported))
    .dynhr_abort(caller, ": ",
                 paste(names(unsupported)[unsupported], collapse = ", "),
                 " not supported by the Student-t likelihood (it would be ",
                 "silently ignored). Use likelihood = \"gaussian\".",
                 class = "dynhr_error_student_t_unsupported")
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## kalman_filter_student_t()
## ---------------------------------------------------------------------------
##
## Runs the standard Gaussian Kalman filter recursions for state prediction
## and update (mean and covariance), but evaluates the per-period
## log-likelihood contribution under a multivariate Student-t density.
##
## This is NOT a robust filter: the Kalman gain is unmodified (it is optimal
## under Gaussianity).  This is the standard first-step approximation used
## widely in the DSGE fat-tail literature.
##
## Parameters:
##   Y           n_obs x T observation matrix (NAs allowed).
##   dr          decision rule (output of solve_perturbation).
##   model       compiled model object.
##   params      named numeric vector of parameter values.
##   obs_vars    character vector of observed variable names.
##   student_df  degrees of freedom nu; must be a finite scalar > 2 (the
##               variance-matching scale needs a finite t variance; nu >= 5
##               also gives a finite kurtosis). A large nu (1e9) is
##               numerically the Gaussian filter.
##   me_variance scalar measurement-error variance (default 0), true i.i.d.
##               noise as in kalman_filter().
##   lik_init    P0 initialization: "auto", "stationary", "kappa"
##               ("diffuse" is not supported here -- use the Gaussian
##                kalman_filter() for the exact diffuse phase, then switch
##                to Student-t for the post-diffuse tail if needed).
##   me_extra, shock_scale, known_shocks
##               NOT supported: the per-period ME / heteroskedastic-shock /
##               known-shock paths of kalman_filter() are not implemented
##               here. Supplying any of them non-NULL -- or a model carrying
##               filter_tunes or heteroskedastic_shocks, which the Gaussian
##               path turns into me_extra / shock_scale -- aborts with class
##               `dynhr_error_student_t_unsupported` rather than being
##               silently ignored.
##
## Returns a list with:
##   loglik     scalar total log-likelihood
##   n_obs      number of observables
##   n_T        number of time periods
##   method     "student_t"
##   student_df the nu used
##   lik_init   initialization used
##
#' @noRd
kalman_filter_student_t <- function(Y, dr, model, params, obs_vars,
                                    student_df,
                                    me_variance  = 0,
                                    lik_init     = "auto",
                                    me_extra     = NULL,
                                    shock_scale  = NULL,
                                    known_shocks = NULL) {
  nu <- student_df
  .student_t_validate(nu, model, me_extra, shock_scale, known_shocks,
                      caller = "kalman_filter_student_t")

  ## -- Extract state-space matrices (same as kalman_filter()) ---------------
  state_idx <- dr$state_idx
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)
  n_obs     <- length(obs_vars)

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("Observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d   <- dr$ys[obs_vars]

  Sigma_e <- .get_shock_cov(model, exo, params)
  QQ      <- tcrossprod(RR %*% Sigma_e, RR)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## Precompute Y - d
  Y_minus_d <- Y - d

  ## -- Scale factor for the t-distribution scale matrix --------------------
  ## Sigma_t = c_nu * F_t with c_nu = (nu - 2) / nu  (nu > 2, checked above).
  ## log(c_nu) via log1p for large-nu accuracy.
  log_c_nu <- log1p(-2 / nu)
  c_nu     <- (nu - 2) / nu

  ## -- Initialization -------------------------------------------------------
  if (lik_init == "auto") {
    tt_evals <- eigen(TT, only.values = TRUE)$values
    if (any(Mod(tt_evals) > 1 - 1e-6)) {
      P0_try <- tryCatch(solve_lyapunov(TT, QQ), error = function(e) NULL)
      ok_stat <- .kf_stationary_P0_ok(P0_try)   # relative rule (W77)
      lik_init <- if (ok_stat) "stationary" else "kappa"
      P0 <- if (ok_stat) P0_try else .build_P0(TT, QQ)
    } else {
      lik_init <- "stationary"
      P0 <- solve_lyapunov(TT, QQ)
    }
  } else if (lik_init == "stationary") {
    P0 <- solve_lyapunov(TT, QQ)
    if (anyNA(P0))
      stop("kalman_filter_student_t: lik_init = \"stationary\" failed ",
           "(unit-root TT). Use lik_init = \"auto\" or \"kappa\".",
           call. = FALSE)
  } else if (lik_init == "kappa") {
    P0 <- .build_P0(TT, QQ)
  } else {
    stop("kalman_filter_student_t: lik_init = \"diffuse\" is not supported; ",
         "use \"auto\", \"stationary\", or \"kappa\".", call. = FALSE)
  }

  ## -- Gaussian Kalman recursions + Student-t log-likelihood ----------------
  s      <- numeric(n_state)
  P      <- P0
  loglik <- 0

  for (t in seq_len(n_T)) {
    v <- Y_minus_d[, t] - as.numeric(ZZ %*% s)

    ## -- Handle missing observations (same NA-reduction as kalman_filter) --
    obs_ok <- which(!is.na(v))
    k_t    <- length(obs_ok)
    if (k_t == 0L) {
      ## Fully missing: pure prediction step, no likelihood contribution
      s <- drop(TT %*% s)
      P <- tcrossprod(TT %*% P, TT) + QQ
      P <- (P + t(P)) * 0.5
      next
    }

    ## Shared measurement/prediction update (R/kf-step.R): F_t, its log-det
    ## and Mahalanobis form, and the Joseph P-update including K me K'.
    st <- .kf_step_core(s, P, v[obs_ok], TT,
                        ZZ[obs_ok, , drop = FALSE], RR,
                        DD[obs_ok, , drop = FALSE], Sigma_e,
                        me_vec = rep(me_variance, k_t))
    if (is.null(st)) {
      loglik <- -Inf; break
    }

    ## -- Student-t scale matrix Sigma_t = c_nu * F_t ----------------------
    ## log|Sigma_t| = k_t*log(c_nu) + log|F_t|
    log_det_Sigma <- k_t * log_c_nu + st$logdet_F

    ## Mahalanobis^2 under Sigma_t: v' Sigma_t^{-1} v = (1/c_nu) * v' F_t^{-1} v
    Q_t <- st$maha / c_nu

    ## Multivariate-t log-density. lgamma((nu+k)/2) - lgamma(nu/2) is
    ## lgamma(k/2) - lbeta(nu/2, k/2), free of the lgamma(nu/2) cancellation.
    ll_t <- lgamma(k_t / 2) - lbeta(nu / 2, k_t / 2) -
      (k_t / 2) * log(nu * pi) -
      0.5 * log_det_Sigma -
      (nu + k_t) / 2 * log1p(Q_t / nu)

    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) {
      loglik <- -Inf; break
    }
    loglik <- loglik + ll_t

    ## -- Gaussian KF update step (state mean + covariance unchanged) -------
    s <- st$s
    P <- st$P
  }

  list(loglik     = loglik,
       n_obs      = n_obs,
       n_T        = n_T,
       method     = "student_t",
       student_df = nu,
       lik_init   = lik_init)
}


## TRUE iff the model carries filter_tunes rows (a bare data.frame of tunes or
## a spec list whose rows live in $tunes -- same test as the tpf guard in
## make_log_posterior; NB a data.frame IS a list, so test is.data.frame first).
.student_t_has_filter_tunes <- function(model) {
  ft <- model$filter_tunes
  if (is.null(ft)) return(FALSE)
  n <- if (is.data.frame(ft)) nrow(ft) else NROW(ft$tunes)
  isTRUE(n > 0L)
}

## TRUE iff the model's heteroskedastic_shocks spec has scale rows (the
## condition under which .build_shock_scale_matrix() returns a non-NULL
## shock_scale on the Gaussian path; parse_mod() attaches an EMPTY spec).
.student_t_has_het_shocks <- function(model) {
  spec <- model$heteroskedastic_shocks
  if (is.null(spec) || is.null(spec$scales)) return(FALSE)
  isTRUE(NROW(spec$scales) > 0L)
}
