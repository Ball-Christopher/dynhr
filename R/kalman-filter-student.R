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
##   1. Scale the t-distribution so that its covariance equals F_t when nu > 2:
##        Sigma_t = (nu - 2) / nu * F_t
##      (for nu <= 2 this is a location-scale t with infinite variance, still
##       a valid density and valid likelihood).
##
##   2. Log-density of multivariate-t with location 0, scale Sigma_t, df nu:
##        loglik_t = lgamma((nu + k_t) / 2)
##                 - lgamma(nu / 2)
##                 - (k_t / 2) * log(nu * pi)
##                 - 0.5 * log|Sigma_t|
##                 - (nu + k_t) / 2 * log(1 + (1/nu) * v_t' Sigma_t^{-1} v_t)
##
##   Substituting Sigma_t = c_nu * F_t with c_nu = (nu - 2) / nu:
##        log|Sigma_t| = k_t * log(c_nu) + log|F_t|
##        v_t' Sigma_t^{-1} v_t = (1/c_nu) * v_t' F_t^{-1} v_t
##
##   For nu > 2 this is exact; for nu <= 2 the scale convention still defines a
##   valid distribution (heavier-tailed, infinite variance).  The Gaussian limit
##   (nu -> Inf) recovers the standard KF loglik exactly (the lgamma ratio and
##   log(1 + Q/nu) both converge to their Gaussian counterparts).
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
## Missing observations: the same NA-reduction as kalman_filter() is used.
## --------------------------------------------------------------------------


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
##   student_df  degrees of freedom nu (must be > 0; nu >= 3 ensures finite
##               kurtosis; nu = 1e7 is numerically indistinguishable from
##               Gaussian).
##   me_variance scalar measurement-error variance (default 0).
##   lik_init    P0 initialization: "auto", "stationary", "kappa"
##               ("diffuse" is not supported here -- use the Gaussian
##                kalman_filter() for the exact diffuse phase, then switch
##                to Student-t for the post-diffuse tail if needed).
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
                                    me_variance = 0,
                                    lik_init    = "auto") {
  nu <- student_df
  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 0)
    stop("kalman_filter_student_t: student_df must be a positive finite scalar.",
         call. = FALSE)

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
  HH      <- tcrossprod(DD %*% Sigma_e, DD)
  SS      <- RR %*% Sigma_e %*% t(DD)
  me_diag <- me_variance * diag(n_obs)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## Precompute Y - d
  Y_minus_d <- Y - d

  ## -- Scale factor for the t-distribution scale matrix --------------------
  ## Sigma_t = c_nu * F_t  with  c_nu = (nu - 2) / nu  (for nu > 2).
  ## For nu <= 2 we use c_nu = 1 (equal scale to F_t): this is a
  ## convention choice for near-degenerate nu.  The density is still valid.
  c_nu <- if (nu > 2) (nu - 2) / nu else 1

  ## -- Initialization -------------------------------------------------------
  if (lik_init == "auto") {
    tt_evals <- eigen(TT, only.values = TRUE)$values
    if (any(Mod(tt_evals) > 1 - 1e-6)) {
      P0_try <- tryCatch(solve_lyapunov(TT, QQ), error = function(e) NULL)
      ok_stat <- !is.null(P0_try) && all(is.finite(P0_try)) &&
        min(Re(eigen((P0_try + t(P0_try)) / 2, symmetric = TRUE,
                     only.values = TRUE)$values)) > -1e-8
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
  tZZ    <- t(ZZ)

  for (t in seq_len(n_T)) {
    v <- Y_minus_d[, t] - as.numeric(ZZ %*% s)

    ## -- Handle missing observations (same NA-reduction as kalman_filter) --
    if (anyNA(v)) {
      obs_ok  <- which(!is.na(v))
      k_t     <- length(obs_ok)
      if (k_t == 0L) {
        ## Fully missing: pure prediction step, no likelihood contribution
        s <- drop(TT %*% s)
        P <- tcrossprod(TT %*% P, TT) + QQ
        P <- (P + t(P)) * 0.5
        next
      }
      ZZ_t  <- ZZ[obs_ok, , drop = FALSE]
      DD_t  <- DD[obs_ok, , drop = FALSE]
      HH_t  <- tcrossprod(DD_t %*% Sigma_e, DD_t)
      SS_t  <- RR %*% Sigma_e %*% t(DD_t)
      me_t  <- me_variance * diag(k_t)
      v_t   <- v[obs_ok]
      Ft    <- ZZ_t %*% P %*% t(ZZ_t) + HH_t + me_t
      Ft    <- (Ft + t(Ft)) * 0.5
    } else {
      k_t   <- n_obs
      ZZ_t  <- ZZ; DD_t <- DD; SS_t <- SS
      v_t   <- v
      Ft    <- ZZ_t %*% P %*% t(ZZ_t) + HH + me_diag
      Ft    <- (Ft + t(Ft)) * 0.5
    }

    ## -- Innovation covariance Cholesky (needed for KF update + t-density) --
    Fc <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) {
      loglik <- -Inf; break
    }
    Fi <- chol2inv(Fc)

    ## -- Student-t scale matrix Sigma_t = c_nu * F_t ----------------------
    ## log|Sigma_t| = k_t*log(c_nu) + log|F_t|
    log_det_Ft    <- 2 * sum(log(diag(Fc)))
    log_det_Sigma <- k_t * log(c_nu) + log_det_Ft

    ## Mahalanobis^2 under Sigma_t: v' Sigma_t^{-1} v = (1/c_nu) * v' F_t^{-1} v
    Q_t <- drop(crossprod(v_t, Fi %*% v_t)) / c_nu

    ## Multivariate-t log-density
    ll_t <- lgamma((nu + k_t) / 2) - lgamma(nu / 2) -
      (k_t / 2) * log(nu * pi) -
      0.5 * log_det_Sigma -
      (nu + k_t) / 2 * log1p(Q_t / nu)

    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) {
      loglik <- -Inf; break
    }
    loglik <- loglik + ll_t

    ## -- Gaussian KF update step (state mean + covariance unchanged) -------
    K <- (TT %*% P %*% t(ZZ_t) + SS_t) %*% Fi
    s <- drop(TT %*% s) + drop(K %*% v_t)

    TmKZ <- TT - K %*% ZZ_t
    RmKD <- RR - K %*% DD_t
    P <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD)
    P <- (P + t(P)) * 0.5
  }

  list(loglik     = loglik,
       n_obs      = n_obs,
       n_T        = n_T,
       method     = "student_t",
       student_df = nu,
       lik_init   = lik_init)
}
