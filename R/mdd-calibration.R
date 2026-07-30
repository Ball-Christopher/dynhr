## R/mdd-calibration.R
## --------------------------------------------------------------------------
## E3: marginal-likelihood (model-evidence) calibration harness.
##
## SBC validates POSTERIOR SHAPE; it says nothing about whether an evidence
## estimator (THAMES, SMC, Laplace) recovers the right NORMALIZING CONSTANT.
## This harness ties each estimator to a ground truth that is either
## closed-form (rung a) or a fine quadrature of the exact pipeline likelihood
## (rung b), so a regression in any estimator shows up as a truth-vs-estimate
## gap instead of silently drifting.
##
## Two built-in cases:
##   "conjugate_regression" (rung a) -- linear-Gaussian regression with a
##     conjugate Gaussian prior; log-evidence is closed form. Hand-rolled
##     closures (no dynhr model object) -- a pure sanity check on the
##     estimator FORMULAS.
##   "local_level"          (rung b) -- the local-level (random-walk-plus-
##     noise) linear-Gaussian state-space model, built via the real
##     parse_mod -> compile_model -> solve_steady -> solve_perturbation ->
##     run_mode_finding pipeline (mirrors tests/testthat/test-kalman-diffuse.R
##     .local_level_model()). One free parameter (the observation-noise
##     stderr) with a proper inv-gamma prior; ground truth is a quadrature
##     of prior x exact diffuse-KF likelihood over a fine grid.
## --------------------------------------------------------------------------


## ---- shared numeric helper -------------------------------------------

## log( sum(exp(x)) ), stable.
.mdd_logsumexp <- function(x) {
  mx <- max(x)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}


## ==========================================================================
## Case builder: conjugate_regression (rung a)
## ==========================================================================

## Builds ONE conjugate-regression subproblem for a single prior scale.
## Wrapped in its own function (not inlined in a for-loop body) so every
## closure captures a FRESH execution environment per prior_scale value --
## closures created directly inside a `for` loop share the loop's single
## environment and would all silently see the LAST iteration's data.
.mdd_conjugate_regression_one <- function(X, y, sigma2, prior_scale,
                                          n_draws, n_particles) {
  n <- nrow(X); d <- ncol(X)
  b0     <- rep(0, d)
  V0     <- diag(rep(prior_scale, d))
  V0_inv <- solve(V0)

  ## ---- analytic evidence: y ~ N(X b0, sigma2*I + X V0 X') --------------
  Sigma          <- sigma2 * diag(n) + X %*% V0 %*% t(X)
  resid          <- y - as.numeric(X %*% b0)
  chol_Sigma     <- chol(Sigma)
  log_det_Sigma  <- 2 * sum(log(diag(chol_Sigma)))
  quad           <- sum(backsolve(chol_Sigma, resid, transpose = TRUE)^2)
  logZ_true      <- -0.5 * n * log(2 * pi) - 0.5 * log_det_Sigma - 0.5 * quad

  ## ---- analytic conjugate posterior -------------------------------------
  XtX       <- crossprod(X)
  post_prec <- XtX / sigma2 + V0_inv
  post_cov  <- solve(post_prec)
  post_mean <- as.numeric(post_cov %*% (crossprod(X, y) / sigma2 + V0_inv %*% b0))
  post_chol <- chol(post_cov)

  log_lik <- function(beta) {
    r <- y - as.numeric(X %*% beta)
    -0.5 * n * log(2 * pi * sigma2) - 0.5 * sum(r^2) / sigma2
  }
  log_prior_fn <- function(beta) {
    r <- beta - b0
    -0.5 * d * log(2 * pi) - 0.5 * log(det(V0)) - 0.5 * as.numeric(t(r) %*% V0_inv %*% r)
  }
  log_post_fn <- function(beta) log_lik(beta) + log_prior_fn(beta)
  log_post_fn_smc <- function(beta) {
    ll <- log_lik(beta); lpr <- log_prior_fn(beta)
    list(loglik = ll, logprior = lpr, logpost = ll + lpr)
  }
  prior_sampler <- function() as.numeric(rnorm(d) %*% chol(V0)) + b0

  ## Minimal `dynhr_mode_result`-shaped list so the EXPORTED laplace_log_marglik()
  ## is exercised directly (not a hand-rolled copy of its formula). The
  ## posterior here is exactly Gaussian, so H = -post_prec exactly and the
  ## Laplace approximation should recover logZ_true to ~machine precision.
  mode_result <- list(
    theta_mode    = post_mean,
    hessian_exact = -post_prec,
    mode          = list(logpost = log_post_fn(post_mean))
  )

  list(
    truth            = logZ_true,
    d                = d,
    log_post_fn      = log_post_fn,
    log_post_fn_smc  = log_post_fn_smc,
    post_mean        = post_mean,
    post_cov         = post_cov,
    post_chol        = post_chol,
    prior_sampler    = prior_sampler,
    mode_result      = mode_result,
    n_draws          = n_draws,
    n_particles      = n_particles
  )
}

#' @param seed        RNG seed for the (fixed) simulated data set.
#' @param n           Number of observations.
#' @param d           Number of regressors (including intercept).
#' @param sigma2      Known residual variance.
#' @param beta_true   Length-`d` true coefficient vector used to simulate `y`.
#' @param prior_scale Numeric vector of conjugate prior variances (`V0 =
#'   diag(rep(prior_scale, d))`); one subproblem is built per entry, letting
#'   `mdd_calibration()` track a multi-value prior-scale sensitivity case in
#'   a single call.
#' @param n_draws     Posterior draws for THAMES (i.i.d. exact Gaussian draws).
#' @param n_particles SMC particle count.
#' @noRd
.mdd_build_conjugate_regression <- function(seed = 20260710L, verbose = FALSE,
                                            n = 40L, d = 3L, sigma2 = 1.3,
                                            beta_true = c(1.5, -0.8, 0.4),
                                            prior_scale = 4,
                                            n_draws = 5000L,
                                            n_particles = 2000L,
                                            ...) {
  stopifnot(length(beta_true) == d)
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (d - 1)), n, d - 1))
  y <- as.numeric(X %*% beta_true + rnorm(n, sd = sqrt(sigma2)))

  subproblems <- list()
  for (ps in prior_scale) {
    key <- sprintf("prior_scale=%g", ps)
    subproblems[[key]] <- .mdd_conjugate_regression_one(
      X, y, sigma2, ps, n_draws = n_draws, n_particles = n_particles)
  }
  list(subproblems = subproblems, axis = "prior_scale")
}


## ==========================================================================
## Case builder: local_level (rung b)
## ==========================================================================

## Builds ONE local-level subproblem for a single (T_obs, prior_sd) pair.
## Own function for the same closure-capture reason as
## .mdd_conjugate_regression_one() above.
.mdd_local_level_one <- function(sigma_eta, sigma_eps_true, T_obs,
                                 prior_mean, prior_sd, prior_lower, prior_upper,
                                 grid_n, n_iter, n_particles,
                                 n_mcmc_draws, n_mcmc_warmup, seed) {
  mod_txt <- sprintf("
  var mu y;
  varexo eta eps;
  parameters sigma_eta sigma_eps;
  sigma_eta = %.10f;
  sigma_eps = %.10f;
  model;
    mu = mu(-1) + eta;
    y = mu + eps;
  end;
  initval;
    mu = 0; y = 0;
  end;
  shocks;
    var eta; stderr sigma_eta;
    var eps; stderr sigma_eps;
  end;
  estimated_params;
    sigma_eps, inv_gamma_pdf, %.10f, %.10f, %.10f, %.10f;
  end;
  ", sigma_eta, sigma_eps_true, prior_mean, prior_sd, prior_lower, prior_upper)

  m  <- parse_mod(mod_txt, verbose = FALSE)
  cm <- compile_model(m, verbose = FALSE)
  ss <- solve_steady(cm, m$param_values, endo_names = m$var_names,
                     exo_names = m$varexo_names, verbose = FALSE)
  if (!isTRUE(ss$converged))
    stop("mdd_calibration('local_level'): steady state failed to converge")
  dr <- solve_perturbation(m, cm, ss$values, m$param_values, verbose = FALSE)
  if (!isTRUE(dr$bk_satisfied))
    stop("mdd_calibration('local_level'): Blanchard-Kahn not satisfied")

  set.seed(seed)
  mu_path <- cumsum(rnorm(T_obs, sd = sigma_eta))
  y <- mu_path + rnorm(T_obs, sd = sigma_eps_true)
  Y <- data.frame(y = y)

  slv <- list(model = m, compiled = cm, dr = dr)
  class(slv) <- "dynhr_solved"

  ## ---- ground truth: quadrature of prior x exact diffuse-KF likelihood ----
  ## Deliberately built from the SAME primitives as the closed-form oracle in
  ## test-kalman-diffuse.R (kalman_filter(..., lik_init = "diffuse")) and the
  ## SAME per-parameter prior density used inside make_log_posterior()
  ## (log_prior_density()), so "truth" is independent of run_mode_finding()'s
  ## internal wiring while still using the real KF/prior code paths.
  grid <- seq(prior_lower, prior_upper, length.out = grid_n)
  pp   <- m$param_values
  loglik_grid <- vapply(grid, function(v) {
    pp["sigma_eps"] <- v
    kalman_filter(matrix(y, nrow = 1), dr, m, pp, obs_vars = "y",
                  lik_init = "diffuse", me_variance = 0)$loglik
  }, numeric(1))
  logprior_grid <- vapply(grid, function(v)
    log_prior_density(v, "inv_gamma_pdf", prior_mean, prior_sd,
                      prior_lower, prior_upper),
    numeric(1))
  dxg  <- diff(grid)[1]
  logZ_quad <- .mdd_logsumexp(loglik_grid + logprior_grid) + log(dxg)

  ## ---- real pipeline: mode + (numerical, since the model is diffuse-init
  ## and posterior_hessian's analytic exact-Hessian path explicitly excludes
  ## lik_init = "diffuse") Hessian, fed into the SAME exported
  ## laplace_log_marglik() used elsewhere in the package. ----------------
  mf <- run_mode_finding(slv, Y, obs_vars = "y", n_iter = n_iter,
                         use_exact_hessian = TRUE, verbose = FALSE)
  if (is.null(mf$hessian_exact)) {
    lp_scalar <- function(th) {
      r <- mf$log_post_fn(th)
      if (is.list(r)) r$logpost else r
    }
    mf$hessian_exact <- if (requireNamespace("numDeriv", quietly = TRUE)) {
      numDeriv::hessian(lp_scalar, mf$theta_mode)
    } else {
      stats::optimHess(mf$theta_mode, lp_scalar)
    }
  }

  list(
    truth          = logZ_quad,
    T_obs          = T_obs,
    prior_sd       = prior_sd,
    m              = m,
    dr             = dr,
    Y              = Y,
    y              = y,
    slv            = slv,
    mf             = mf,
    prior_spec_df  = prior_spec(m),
    n_particles    = n_particles,
    n_mcmc_draws   = n_mcmc_draws,
    n_mcmc_warmup  = n_mcmc_warmup
  )
}

#' @param seed            RNG seed for data simulation (shared across the
#'   T_obs/prior_sd ladder so subproblems differ only in the axis varied).
#' @param sigma_eta       Known state-innovation stderr.
#' @param sigma_eps_true  True observation-noise stderr used to simulate data.
#' @param T_obs           Numeric vector of sample sizes. When of length > 1
#'   (with `prior_sd` left at its scalar default) this builds a T-ladder used
#'   to show the Laplace residual shrinking with T.
#' @param prior_mean,prior_sd,prior_lower,prior_upper Inv-gamma prior
#'   hyperparameters for `sigma_eps` (`prior_sd` may be a vector for a
#'   prior-scale sensitivity ladder; `T_obs`/`prior_sd` are recycled against
#'   each other, i.e. a PARALLEL, not Cartesian, ladder).
#' @param grid_n          Quadrature grid points for the ground truth.
#' @param n_iter          Mode-finding iteration budget.
#' @param n_particles     SMC particle count.
#' @param n_mcmc_draws,n_mcmc_warmup RWMH draws/warmup for THAMES input.
#' @noRd
.mdd_build_local_level <- function(seed = 123L, verbose = FALSE,
                                   sigma_eta = 0.5, sigma_eps_true = 0.3,
                                   T_obs = 40L,
                                   prior_mean = 0.3, prior_sd = 0.15,
                                   prior_lower = 0.01, prior_upper = 3,
                                   grid_n = 800L, n_iter = 2000L,
                                   n_particles = 1000L,
                                   n_mcmc_draws = 3000L, n_mcmc_warmup = 1000L,
                                   ...) {
  n_sub      <- max(length(T_obs), length(prior_sd))
  T_obs_v    <- rep_len(T_obs, n_sub)
  prior_sd_v <- rep_len(prior_sd, n_sub)

  subproblems <- list()
  for (i in seq_len(n_sub)) {
    key <- sprintf("T=%d_priorsd=%g", T_obs_v[i], prior_sd_v[i])
    subproblems[[key]] <- .mdd_local_level_one(
      sigma_eta = sigma_eta, sigma_eps_true = sigma_eps_true,
      T_obs = T_obs_v[i], prior_mean = prior_mean, prior_sd = prior_sd_v[i],
      prior_lower = prior_lower, prior_upper = prior_upper,
      grid_n = grid_n, n_iter = n_iter, n_particles = n_particles,
      n_mcmc_draws = n_mcmc_draws, n_mcmc_warmup = n_mcmc_warmup, seed = seed)
  }
  list(subproblems = subproblems, axis = if (length(T_obs) > 1) "T_obs" else "prior_sd")
}


## ==========================================================================
## Case builder: dsge_var (rung a, exact closed form; validates d15)
## ==========================================================================

## Exact conjugate matrix-variate Normal-Inverse-Wishart (MNIW) log marginal
## likelihood for a dummy-observation ("Minnesota"/DSGE-prior) BVAR, i.e. the
## Del Negro & Schorfheide (2004) DSGE-VAR evidence at a given tightness
## lambda:
##
##   ln p(Y|lambda) = -(n*T/2) ln(pi) + [ln Gamma_n(v1/2) - ln Gamma_n(v0/2)]
##     + (n/2) ln|XtX0| - (n/2) ln|XtX0 + XtX1|
##     + (v0/2) ln|S0|  - (v1/2) ln|S1|
##
## where (XtX0, XtY0, YtY0, T0) are the prior/dummy sufficient statistics
## implied by the DSGE model's population autocovariances (scaled by the
## dummy-observation count T0 = round(lambda * n_obs)), (XtX1, XtY1, YtY1,
## T1) are the REAL data's OLS sufficient statistics, v0 = T0 - k + n + 1,
## v1 = v0 + T1, S0/S1 are the prior/posterior residual sums of squares, and
## Gamma_n is the multivariate gamma function (the pi^{n(n-1)/4} prefactor
## is identical in Gamma_n(v1/2) and Gamma_n(v0/2) and cancels).
##
## Deliberately re-derived/re-coded from scratch here (own lgamma-ratio sum,
## own solve()-based residual sums of squares, no ridge term) rather than
## calling d15_dsge_var()'s internal .d15_lmvgamma_ratio()/matrix pipeline,
## so this is an INDEPENDENT check on the package formula, not a copy of it.
## The formula itself is additionally validated against brute-force 2D
## numerical quadrature (n=1, k=1 reduction) in test-mdd-calibration.R.
.mdd_dsge_var_exact_log_ml <- function(lam, n_obs, k_coef, T_eff,
                                       G, rhs, Gamma_0, XtX1, XtY1, YtY1) {
  T0 <- max(1L, round(lam * n_obs))
  XtX0 <- matrix(0, k_coef, k_coef)
  XtX0[1L, 1L]   <- T0
  XtX0[-1L, -1L] <- T0 * G
  XtY0 <- matrix(0, k_coef, n_obs)
  XtY0[-1L, ] <- T0 * rhs
  YtY0 <- T0 * Gamma_0

  v0 <- T0 - k_coef + n_obs + 1
  v1 <- (T0 + T_eff) - k_coef + n_obs + 1

  XtX_post <- XtX0 + XtX1
  XtY_post <- XtY0 + XtY1
  YtY_post <- YtY0 + YtY1

  S0 <- YtY0 - t(XtY0) %*% solve(XtX0, XtY0)
  A_post <- solve(XtX_post, XtY_post)
  S1 <- YtY_post - t(XtY_post) %*% A_post

  i <- seq_len(n_obs)
  lmvgamma_ratio <- sum(lgamma((v1 - i + 1) / 2) - lgamma((v0 - i + 1) / 2))

  lmvgamma_ratio +
    0.5 * n_obs * log(det(XtX0)) - 0.5 * n_obs * log(det(XtX_post)) +
    0.5 * v0 * log(det(S0))      - 0.5 * v1 * log(det(S1)) -
    0.5 * T_eff * n_obs * log(pi)
}

## Builds one subproblem PER lambda grid point of the DSGE-VAR diagnostic
## d15_dsge_var(): a small stationary bivariate "DSGE" that is itself a
## VAR(1) (so the model-implied prior moments coincide with the DGP) is
## solved through the real parse_mod -> compile_model -> solve_steady ->
## solve_perturbation pipeline, data are simulated from it, and
## d15_dsge_var() is run ONCE over the full lambda grid. d15's own grid
## search naturally yields a multi-value ladder (>= 2 distinct lambda with
## distinct truths, the CLAUDE.md multi-value convention), so unlike the
## other two builders no extra per-value looping is needed -- the grid IS
## the ladder.
.mdd_build_dsge_var <- function(seed = 2026L, verbose = FALSE,
                                rho_x = 0.6, rho_z = 0.5, phi = 0.3,
                                T_obs = 80L, var_lag = 1L,
                                lambda_grid = c(0.75, 1, 2, 5, 10, 20, 50, 100),
                                ...) {
  mod_txt <- sprintf("
  var x z;
  varexo ex ez;
  parameters rho_x rho_z phi sx sz;
  rho_x = %.10f; rho_z = %.10f; phi = %.10f; sx = 1; sz = 1;
  model;
    x = rho_x*x(-1) + ex;
    z = rho_z*z(-1) + phi*x(-1) + ez;
  end;
  initval; x = 0; z = 0; end;
  shocks;
    var ex; stderr sx;
    var ez; stderr sz;
  end;
  ", rho_x, rho_z, phi)

  m  <- parse_mod(mod_txt, verbose = FALSE)
  cm <- compile_model(m, verbose = FALSE)
  ss <- solve_steady(cm, m$param_values, endo_names = m$var_names,
                     exo_names = m$varexo_names, verbose = FALSE)
  if (!isTRUE(ss$converged))
    stop("mdd_calibration('dsge_var'): steady state failed to converge")
  dr <- solve_perturbation(m, cm, ss$values, m$param_values, verbose = FALSE)
  if (!isTRUE(dr$bk_satisfied))
    stop("mdd_calibration('dsge_var'): Blanchard-Kahn not satisfied")

  set.seed(seed)
  Sigma_e   <- diag(2)
  state_idx <- dr$state_idx
  T_mat <- dr$ghx[state_idx, , drop = FALSE]
  R_mat <- dr$ghu[state_idx, , drop = FALSE]
  n_state <- nrow(T_mat)
  obs_idx <- match(c("x", "z"), dr$endo_names)
  Z_mat <- dr$ghx[obs_idx, , drop = FALSE]
  D_mat <- dr$ghu[obs_idx, , drop = FALSE]

  s      <- matrix(0, T_obs + 1L, n_state)
  shocks <- matrix(rnorm(T_obs * 2L), T_obs, 2L)
  Y      <- matrix(NA_real_, T_obs, 2L)
  for (t in seq_len(T_obs)) {
    s[t + 1L, ] <- as.numeric(T_mat %*% s[t, ] + R_mat %*% shocks[t, ])
    Y[t, ]      <- as.numeric(Z_mat %*% s[t + 1L, ] + D_mat %*% shocks[t, ])
  }
  colnames(Y) <- c("x", "z")

  res <- suppressMessages(d15_dsge_var(
    dr = dr, data = Y, obs_names = c("x", "z"), sigma_e = Sigma_e,
    model = m, var_lag = var_lag, lambda_grid = lambda_grid))
  r <- res$result
  if (is.null(r))
    stop("mdd_calibration('dsge_var'): d15_dsge_var() returned no result (",
         res$summary, ")")

  n_obs <- r$n_obs; k_coef <- r$k_coef; T_eff <- r$T_eff
  G <- r$G_prior; rhs <- r$rhs_prior; Gamma_0 <- r$Gamma0_prior
  XtX1 <- r$XtX_data; XtY1 <- r$XtY_data; YtY1 <- r$YtY_data

  keep <- is.finite(r$log_ml_values)
  subproblems <- list()
  for (i in which(keep)) {
    lam   <- r$lambda_grid[i]
    truth <- .mdd_dsge_var_exact_log_ml(
      lam = lam, n_obs = n_obs, k_coef = k_coef, T_eff = T_eff,
      G = G, rhs = rhs, Gamma_0 = Gamma_0,
      XtX1 = XtX1, XtY1 = XtY1, YtY1 = YtY1)
    key <- sprintf("lambda=%g", lam)
    subproblems[[key]] <- list(truth = truth, d15_log_ml = r$log_ml_values[i])
  }

  list(subproblems = subproblems, axis = "lambda")
}


## ==========================================================================
## Estimator runners
## ==========================================================================

.mdd_run_estimator <- function(case, est, sp, seed) {
  switch(
    paste(case, est, sep = "::"),

    "conjugate_regression::thames" = {
      set.seed(seed)
      Z     <- matrix(rnorm(sp$n_draws * sp$d), sp$n_draws, sp$d)
      draws <- sweep(Z %*% sp$post_chol, 2, sp$post_mean, "+")
      lp    <- apply(draws, 1, sp$log_post_fn)
      res   <- thames_mdd(draws, lp, split = TRUE)
      list(estimate = res$log_mdd, se = res$se)
    },

    "conjugate_regression::smc" = {
      set.seed(seed)
      res <- dynhr_smc(log_post_fn = sp$log_post_fn_smc,
                       prior_sampler = sp$prior_sampler,
                       n_particles = sp$n_particles, seed_base = seed,
                       verbose = FALSE)
      list(estimate = res$log_marginal_lik, se = NA_real_)
    },

    "conjugate_regression::laplace" = {
      val <- laplace_log_marglik(sp$mode_result)
      list(estimate = val, se = NA_real_)
    },

    "local_level::thames" = {
      set.seed(seed)
      ch  <- mcmc(sp$mf$log_post_fn, sp$mf$theta_mode, sp$mf$Sigma_prop,
                  n_draws = sp$n_mcmc_draws, n_warmup = sp$n_mcmc_warmup,
                  verbose = FALSE)
      res <- thames_mdd_from_chains(ch)
      list(estimate = res$log_mdd, se = res$se)
    },

    "local_level::smc" = {
      set.seed(seed)
      res <- smc(sp$mf$log_post_fn, sp$prior_spec_df,
                n_particles = sp$n_particles, seed_base = seed, verbose = FALSE)
      list(estimate = res$log_marginal_lik, se = NA_real_)
    },

    "local_level::laplace" = {
      val <- laplace_log_marglik(sp$mf)
      list(estimate = val, se = NA_real_)
    },

    "dsge_var::d15" = {
      ## d15_dsge_var() is deterministic given (dr, data) -- the "estimate"
      ## was already computed once in .mdd_build_dsge_var() (the package's
      ## own lambda-grid evidence at this subproblem's lambda).
      list(estimate = sp$d15_log_ml, se = NA_real_)
    },

    stop("mdd_calibration: unknown case/estimator combination '", case, "::", est, "'")
  )
}


## ==========================================================================
## Public entry point
## ==========================================================================

#' Ground-truth calibration harness for marginal-likelihood estimators
#'
#' SBC (\code{sbc_matrix_result()} and friends) validates POSTERIOR SHAPE: it
#' checks that credible intervals have nominal coverage. It says nothing
#' about whether an evidence/marginal-likelihood estimator
#' (\code{\link{thames_mdd}}, \code{\link{smc}}, \code{\link{laplace_log_marglik}})
#' recovers the right NORMALIZING CONSTANT, because SBC ranks are invariant
#' to \eqn{\log Z}. \code{mdd_calibration()} is a dev/diagnostic tool that
#' closes that gap: given a case with a KNOWN ground-truth log-evidence, it
#' runs the requested estimator(s) and reports estimate vs. truth.
#'
#' @details
#' Two built-in cases:
#' \describe{
#'   \item{\code{"conjugate_regression"}}{Linear-Gaussian regression with a
#'     conjugate Gaussian prior on the coefficients (known residual
#'     variance). Log-evidence is closed form. Uses hand-rolled log-density
#'     closures (no dynhr model object) -- a pure check on the estimator
#'     FORMULAS/implementations, decoupled from the estimation pipeline.
#'     Supports a `prior_scale` VECTOR: one subproblem is built per entry
#'     (distinct \eqn{V_0} => distinct truth), letting a single call track a
#'     multi-value prior-scale sensitivity case.}
#'   \item{\code{"local_level"}}{The local-level (random-walk-plus-noise)
#'     state-space model, built through the real
#'     \code{\link{parse_mod}} -> \code{\link{compile_model}} ->
#'     \code{\link{solve_steady}} -> \code{\link{solve_perturbation}} ->
#'     \code{\link{run_mode_finding}} pipeline (mirrors the fixture in
#'     \code{tests/testthat/test-kalman-diffuse.R}). One free parameter (the
#'     observation-noise stderr `sigma_eps`) with a proper inverse-gamma
#'     prior; the model has a unit root, so \code{lik_init} auto-resolves to
#'     \code{"diffuse"} exact initialization. Ground truth is a quadrature
#'     (default 800-point grid) of prior x exact diffuse-KF likelihood, built
#'     from the SAME \code{\link{kalman_filter}}/prior-density primitives the
#'     pipeline uses, independent of \code{run_mode_finding()}'s internal
#'     wiring. This is the rung that exercises the real estimation
#'     machinery rather than hand-rolled closures. Supports a `T_obs` VECTOR
#'     (sample-size ladder, for showing the Laplace residual shrink with T)
#'     and/or a `prior_sd` vector (prior-scale sensitivity ladder); the two
#'     are recycled against each other (a parallel, not Cartesian, ladder).}
#'   \item{\code{"dsge_var"}}{The DSGE-VAR tightness diagnostic
#'     (\code{d15_dsge_var()}), Del Negro & Schorfheide (2004): a small
#'     stationary bivariate "DSGE" that is itself a VAR(1) is solved through
#'     the real \code{parse_mod} -> \code{compile_model} -> \code{solve_steady}
#'     -> \code{solve_perturbation} pipeline, data are simulated from it, and
#'     \code{d15_dsge_var()} is run once over a \code{lambda_grid}. Ground
#'     truth is the exact closed-form conjugate matrix-variate
#'     Normal-Inverse-Wishart (MNIW) evidence, re-derived/re-coded
#'     independently of \code{d15_dsge_var()}'s internal formula (own
#'     multivariate-gamma-ratio sum, own \code{solve()}-based residual sums
#'     of squares). One subproblem is built PER lambda grid point -- d15's
#'     own grid search yields a multi-value lambda ladder with distinct
#'     truths for free.}
#' }
#'
#' Estimators (\code{estimators}, any of):
#' \describe{
#'   \item{\code{"thames"}}{\code{\link{thames_mdd}}
#'     (\code{conjugate_regression}: applied to i.i.d. exact-Gaussian
#'     posterior draws) or \code{\link{thames_mdd_from_chains}}
#'     (\code{local_level}: applied to \code{\link{mcmc}} RWMH draws around
#'     the mode). Reports a self-reported Monte Carlo \code{se}.}
#'   \item{\code{"smc"}}{The internal data-tempered SMC evidence estimator
#'     (\code{conjugate_regression}: raw \code{dynhr_smc()}; \code{local_level}:
#'     the exported \code{\link{smc}()} wrapper). No closed-form \code{se}.}
#'   \item{\code{"laplace"}}{\code{\link{laplace_log_marglik}}, the EXPORTED
#'     function (not a re-derivation), applied to a \code{dynhr_mode_result}-
#'     shaped list. \code{conjugate_regression} feeds it the exact analytic
#'     Hessian (posterior is exactly Gaussian there, so the error should be
#'     ~machine epsilon -- a sanity check on the formula). \code{local_level}
#'     feeds it a numerical Hessian (numDeriv, or \code{stats::optimHess} as
#'     a fallback) of the real posterior, because \code{run_mode_finding}'s
#'     analytic exact-Hessian path explicitly excludes
#'     \code{lik_init = "diffuse"} models -- deterministic, so \code{n_reps}
#'     is forced to 1 regardless of the requested value.}
#'   \item{\code{"d15"}}{\code{dsge_var} only: the package's own
#'     \code{d15_dsge_var()} evidence at the subproblem's lambda -- deterministic
#'     given (dr, data), so \code{n_reps} is forced to 1.}
#' }
#'
#' For stochastic estimators (\code{"thames"}, \code{"smc"}), setting
#' \code{n_reps > 1} repeats the estimator with distinct seeds (holding the
#' simulated data fixed) and the returned \code{summary} adds a bias t-stat
#' (\code{mean(estimate - truth) / (sd(estimate) / sqrt(n_reps))}) and, for
#' THAMES, the empirical coverage of its self-reported 95\% interval
#' (\code{estimate +/- 1.96*se}).
#'
#' @param case       \code{"conjugate_regression"}, \code{"local_level"}, or
#'   \code{"dsge_var"}.
#' @param estimators Character vector, any of \code{"thames"}, \code{"smc"},
#'   \code{"laplace"}, \code{"d15"}. Default \code{c("thames", "laplace")}
#'   (the fast, near-deterministic pair); pass \code{"smc"} explicitly to
#'   include it, or \code{"d15"} for the \code{dsge_var} case.
#' @param n_reps     Number of replications for stochastic estimators
#'   (default 1). Ignored (forced to 1) for \code{"laplace"}/\code{"d15"}.
#' @param seed       Base RNG seed; replication \code{r} of estimator
#'   \code{estimators[j]} uses a deterministic derived seed so different
#'   estimators/reps do not share a stream.
#' @param verbose    Passed through to the case builder (currently unused
#'   there; reserved).
#' @param ...        Case-specific arguments forwarded to the builder --
#'   see Details. E.g. \code{prior_scale = c(4, 1)} for
#'   \code{conjugate_regression}, \code{T_obs = c(20L, 80L, 320L)} for
#'   \code{local_level}, or \code{lambda_grid = c(1, 5, 20)} for
#'   \code{dsge_var}.
#'
#' @return A list of class \code{dynhr_mdd_calibration} with elements:
#'   \describe{
#'     \item{results}{Tidy per-replication data.frame: \code{case},
#'       \code{subproblem}, \code{estimator}, \code{rep}, \code{truth},
#'       \code{estimate}, \code{se}, \code{error} (\code{estimate - truth}),
#'       \code{abs_error}.}
#'     \item{summary}{One row per (\code{subproblem}, \code{estimator}):
#'       \code{n_reps}, \code{mean_estimate}, \code{mean_error},
#'       \code{rel_error}, \code{sd_estimate}, \code{bias_t_stat},
#'       \code{mean_se}, \code{coverage_95} (THAMES only, else \code{NA}).}
#'   }
#'
#' @seealso \code{\link{thames_mdd}}, \code{\link{smc}},
#'   \code{\link{laplace_log_marglik}}
#' @export
mdd_calibration <- function(case = c("conjugate_regression", "local_level", "dsge_var"),
                            estimators = c("thames", "laplace"),
                            n_reps = 1L,
                            seed = 1L,
                            verbose = FALSE,
                            ...) {
  case       <- match.arg(case)
  estimators <- match.arg(estimators, choices = c("thames", "smc", "laplace", "d15"),
                          several.ok = TRUE)
  if (!is.numeric(n_reps) || length(n_reps) != 1L || n_reps < 1L)
    stop("mdd_calibration: n_reps must be a positive integer scalar.")
  n_reps <- as.integer(n_reps)

  problem <- switch(case,
    conjugate_regression = .mdd_build_conjugate_regression(seed = seed, verbose = verbose, ...),
    local_level           = .mdd_build_local_level(seed = seed, verbose = verbose, ...),
    dsge_var              = .mdd_build_dsge_var(seed = seed, verbose = verbose, ...)
  )

  rows <- list()
  for (sp_name in names(problem$subproblems)) {
    sp <- problem$subproblems[[sp_name]]
    for (est in estimators) {
      reps <- if (est %in% c("laplace", "d15")) 1L else n_reps
      for (r in seq_len(reps)) {
        ## Distinct, deterministic seed per (subproblem, estimator, rep).
        rep_seed <- seed + 1000L * which(sp_name == names(problem$subproblems)) +
          100L * match(est, c("thames", "smc", "laplace", "d15")) + r
        res <- .mdd_run_estimator(case, est, sp, seed = rep_seed)
        rows[[length(rows) + 1L]] <- data.frame(
          case        = case,
          subproblem  = sp_name,
          estimator   = est,
          rep         = r,
          truth       = sp$truth,
          estimate    = res$estimate,
          se          = res$se %||% NA_real_,
          error       = res$estimate - sp$truth,
          abs_error   = abs(res$estimate - sp$truth),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  results_df <- do.call(rbind, rows)
  rownames(results_df) <- NULL

  summary_rows <- list()
  for (sp_name in unique(results_df$subproblem)) {
    for (est in unique(results_df$estimator)) {
      sub <- results_df[results_df$subproblem == sp_name & results_df$estimator == est, ]
      if (nrow(sub) == 0L) next
      truth        <- sub$truth[1]
      n_r          <- nrow(sub)
      mean_est     <- mean(sub$estimate)
      mean_err     <- mean_est - truth
      sd_est       <- if (n_r > 1L) sd(sub$estimate) else NA_real_
      bias_t       <- if (n_r > 1L && is.finite(sd_est) && sd_est > 0)
        mean_err / (sd_est / sqrt(n_r)) else NA_real_
      mean_se      <- if (all(is.finite(sub$se))) mean(sub$se) else NA_real_
      coverage_95  <- if (identical(est, "thames") && all(is.finite(sub$se)))
        mean(sub$abs_error <= 1.96 * sub$se) else NA_real_

      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        case          = case,
        subproblem    = sp_name,
        estimator     = est,
        n_reps        = n_r,
        truth         = truth,
        mean_estimate = mean_est,
        mean_error    = mean_err,
        rel_error     = mean_err / max(abs(truth), 1e-8),
        sd_estimate   = sd_est,
        bias_t_stat   = bias_t,
        mean_se       = mean_se,
        coverage_95   = coverage_95,
        stringsAsFactors = FALSE
      )
    }
  }
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  structure(
    list(case = case, axis = problem$axis, results = results_df,
        summary = summary_df, problem = problem),
    class = "dynhr_mdd_calibration"
  )
}
