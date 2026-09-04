## R/global-likelihood.R
## --------------------------------------------------------------------------
## Bootstrap particle-filter likelihood on a GLOBAL (projection) solution, and
## its log-posterior factory.
##
## Before E1-B `solve_global()` had ZERO callers outside its own tests: the
## projection solver could produce a policy but no estimator could consume
## one, so every likelihood in the package was tied to a perturbation
## solution.  `make_log_posterior_global_pf()` closes that gap.  Per draw:
##
##   theta -> params -> steady state -> solve_global()  (projection policy)
##                                   -> solve_perturbation() (order 1, used
##                                      ONLY for the stationary particle
##                                      initialisation P0)
##         -> bootstrap particle filter on the NONLINEAR policy
##
## STATE AND TIMING.  A particle carries the state LAG s_{t-1} in LEVELS --
## exactly the vector `predict.GlobalSolution()` takes.  One period is
##
##   eps_t     ~ N(0, Sigma_e)
##   y_t        = policy(feed(s_{t-1}, eps_t))
##   w_t        = N(obs_t ; y_t[obs_vars], me_variance * I)
##   s_t        = y_t[state_names]
##
## where `feed()` is solve_global()'s AR(1) shock injection (see
## `compute_next_lag()` there): the projection policy is a function of the
## state lag alone, with eps_t = 0, and a realised shock enters through
## s_lag[nm] + eps[k]/rho on the AR(1) states.  Sharing the solver's own
## `ar1_rho`/`ar1_shock_idx` rather than re-deriving them is deliberate.
##
## This is the same state-space timing the Kalman filter uses
## (obs_t = d + ZZ x_{t-1} + DD u_t), so the two likelihoods are directly
## comparable at a near-linear calibration -- which is the certification
## oracle in test-global-likelihood.R.
##
## It is ALSO the timing `simulate.GlobalSolution()` uses, but only since
## F1-C: until then the simulator applied eps_t to the lag entering period
## t + 1, so its output was this recursion shifted one period, with a
## deterministic first row.  See the timing audit in the header of
## R/global-sbc.R and the regression in test-global-sbc.R.  The SBC harness
## reuses `.gpf_feed_lag()` below for its DGP so that the generative law and
## the filtered law cannot drift apart.
##
## me_variance > 0 IS REQUIRED.  A bootstrap filter weights particles by the
## measurement density; with a singular observation block every weight is 0 or
## Inf and the estimate is meaningless.  Same restriction (and same reason) as
## make_log_posterior_tpf().
##
## UNBIASEDNESS.  Multinomial-free systematic resampling
## (`.smc_systematic_resample`, shared with the SMC sampler and the TPF -- this
## file deliberately writes NO third resampler) keeps the marginal-likelihood
## estimate unbiased, so RWMH over the returned closure is a valid PMMH.  As
## with every particle likelihood, a NON-NULL `seed` freezes the cloud and
## breaks that pseudo-marginal validity; it is there for reproducible one-off
## evaluations and the tests.
## --------------------------------------------------------------------------


## Apply solve_global()'s AR(1) "feed" to a WHOLE particle matrix at once.
##
## `lag` is N x n_state (levels, columns ordered as g$state_names) and `eps`
## is N x n_exo.  Returns the N x n_state matrix to hand to
## predict.GlobalSolution().  Scalar-for-scalar identical to
## `g$compute_next_lag()`; vectorised because the filter applies it N times
## per period.
#' @noRd
.gpf_feed_lag <- function(g, lag, eps) {
  rho <- g$ar1_rho
  psi <- g$ar1_psi
  si  <- g$ar1_shock_idx
  for (j in seq_along(g$state_names)) {
    k <- si[[j]]
    if (!is.na(k)) lag[, j] <- lag[, j] + psi[[j]] * eps[, k] / rho[[j]]
  }
  lag
}


## Bootstrap particle-filter log-likelihood on a GlobalSolution.
##
## @param Y        n_obs x T observation matrix (NA allowed: a period's
##   missing dimensions are dropped from the weight, exactly as the Kalman
##   filter drops them from its innovation).
## @param g        GlobalSolution.
## @param obs_idx  Positions of obs_vars in g$all_endo_names.
## @param s0       Length-n_state initial state-lag MEAN (levels).
## @param L0       n_state x n_state factor with tcrossprod(L0) = P0.
## @param Le       n_exo x n_exo factor with tcrossprod(Le) = Sigma_e.
## @param me_var   Scalar measurement-error variance (> 0).
## @param N        Particle count.
## @param ess_frac Resample when ESS < ess_frac * N.
## @return Scalar log-likelihood estimate, or -Inf.
#' @noRd
.global_pf_loglik <- function(Y, g, obs_idx, s0, L0, Le, me_var, N,
                              ess_frac = 0.5) {
  n_state <- length(g$state_names)
  n_exo   <- length(g$shock_names)
  n_obs   <- nrow(Y)
  n_T     <- ncol(Y)

  ## Particles: state LAG in levels, N x n_state.
  lag <- matrix(s0, N, n_state, byrow = TRUE) +
    matrix(stats::rnorm(N * n_state), N, n_state) %*% t(L0)
  colnames(lag) <- g$state_names

  ll     <- 0
  log_N  <- log(N)
  const1 <- -0.5 * log(2 * pi * me_var)
  ## NORMALISED log weights, carried across periods so that the ESS-triggered
  ## resampling stays unbiased: the increment is log sum_i W_i w_i(t), not
  ## log mean_i w_i(t).  Between two resamplings the cloud is NOT equally
  ## weighted, and using the equal-weight increment there is the classic
  ## adaptive-resampling bias.
  logW <- rep(-log_N, N)

  for (t in seq_len(n_T)) {
    eps  <- matrix(stats::rnorm(N * n_exo), N, n_exo) %*% t(Le)
    feed <- .gpf_feed_lag(g, lag, eps)
    y_t  <- tryCatch(predict(g, feed), error = function(e) NULL)
    if (is.null(y_t) || !all(is.finite(y_t))) return(-Inf)

    yh <- y_t[, obs_idx, drop = FALSE]
    yv <- Y[, t]
    ok <- which(!is.na(yv))
    if (length(ok)) {
      dev <- sweep(yh[, ok, drop = FALSE], 2L, yv[ok], "-")
      lw  <- length(ok) * const1 - 0.5 * rowSums(dev * dev) / me_var
    } else {
      lw <- numeric(N)          # fully-missing period contributes nothing
    }
    if (!any(is.finite(lw))) return(-Inf)

    lz <- logW + lw
    mx <- max(lz)
    if (!is.finite(mx)) return(-Inf)
    su <- sum(exp(lz - mx))
    if (!is.finite(su) || su <= 0) return(-Inf)
    ll_t <- mx + log(su)
    ll   <- ll + ll_t
    logW <- lz - ll_t                  # renormalised in log space

    ## The time-t state (= the lag entering t+1) is the policy output; the
    ## particle cloud advances BEFORE any resampling, so a resampled cloud is
    ## already the period-t filtered state lag.
    lag <- y_t[, g$state_names, drop = FALSE]

    w   <- exp(logW)
    ess <- 1 / sum(w * w)
    if (ess < ess_frac * N) {
      idx  <- .smc_systematic_resample(w, N)
      lag  <- lag[idx, , drop = FALSE]
      logW <- rep(-log_N, N)
    }
  }
  ll
}


#' Log-posterior with a global-solution bootstrap particle filter
#'
#' Solves the model GLOBALLY (Chebyshev projection, \code{\link{solve_global}})
#' at every parameter draw and evaluates the likelihood with a bootstrap
#' particle filter run on the resulting NONLINEAR policy: the state transition
#' is the policy plus the drawn shocks and the measurement density is
#' \eqn{N(obs_t; y_t[obs\_vars], \sigma^2_{me} I)}.  Unlike every other
#' likelihood in the package this one never linearises the model, so it is the
#' path to take when the nonlinearity IS the object of interest.
#'
#' The returned estimate is unbiased for the marginal likelihood, so
#' \code{\link{pmmh}()} -- random-walk Metropolis driven by an unbiased
#' likelihood estimate -- is valid over this closure, provided
#' \code{seed = NULL} (the default).
#'
#' @param model      Parsed model from \code{\link{parse_mod}}.
#' @param data       Observation matrix (\code{n_obs x T} or \code{T x n_obs};
#'   re-oriented automatically, as in \code{\link{kalman_filter}}).
#' @param prior_spec Prior specification from \code{\link{extract_prior_spec}}.
#' @param obs_vars   Character vector of observed variable names.
#' @param compiled   \code{dynhr_compiled} from \code{\link{compile_model}}.
#' @param me_variance Measurement-error variance, a positive scalar.
#'   \strong{Required}: a bootstrap filter cannot weight particles without it.
#' @param n_particles Number of particles (default 1000).
#' @param poly_degree,n_quad,n_nodes,state_domain,solve_tol,solve_max_iter
#'   Passed to \code{\link{solve_global}} at every draw
#'   (\code{solve_tol}/\code{solve_max_iter} are its \code{tol}/\code{max_iter}).
#' @param ess_frac   Resample threshold as a fraction of \code{n_particles}
#'   (default 0.5); see Details in \code{\link{make_log_posterior_tpf}} for the
#'   same convention.
#' @param seed       Integer seed, or \code{NULL} (default) for fresh
#'   randomness at every evaluation.  A non-\code{NULL} seed is applied
#'   through \code{.with_local_seed()}, so the caller's global RNG stream is
#'   left exactly as it was found.
#' @param power      Power-posterior tempering exponent zeta;
#'   \code{NULL} (default) resolves the \code{power_posterior} option once at
#'   factory time, as in \code{\link{make_log_posterior}}.
#' @param system_priors Optional system-prior specification.
#' @return A \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}.
#' @seealso \code{\link{solve_global}}, \code{\link{make_log_posterior}},
#'   \code{\link{euler_errors}}
#' @export
make_log_posterior_global_pf <- function(model, data, prior_spec, obs_vars,
                                         compiled,
                                         me_variance = 0,
                                         n_particles = 1000L,
                                         poly_degree = 3L,
                                         n_quad      = 5L,
                                         n_nodes     = 7L,
                                         state_domain = NULL,
                                         solve_tol      = 1e-7,
                                         solve_max_iter = 200L,
                                         ess_frac    = 0.5,
                                         seed        = NULL,
                                         power       = NULL,
                                         system_priors = NULL) {
  force(data); force(prior_spec); force(obs_vars); force(me_variance)
  force(n_particles); force(seed); force(state_domain)
  force(poly_degree); force(n_quad); force(n_nodes)
  force(solve_tol); force(solve_max_iter); force(ess_frac)

  power <- .dynhr_opt("power_posterior", power, default = 1)
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0)
    stop("make_log_posterior_global_pf: `power` must be a finite scalar ",
         "in (0, 1].", call. = FALSE)
  power <- as.numeric(power)

  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance <= 0)
    stop("make_log_posterior_global_pf: `me_variance` must be a positive ",
         "scalar. A bootstrap particle filter weights particles by the ",
         "measurement density; with me_variance = 0 every weight is 0 or Inf ",
         "and the likelihood estimate is meaningless.", call. = FALSE)

  if (!is.numeric(n_particles) || length(n_particles) != 1L ||
      n_particles < 2)
    stop("make_log_posterior_global_pf: `n_particles` must be >= 2.",
         call. = FALSE)
  n_particles <- as.integer(n_particles)

  if (!is.numeric(ess_frac) || length(ess_frac) != 1L ||
      !is.finite(ess_frac) || ess_frac < 0 || ess_frac > 1)
    stop("make_log_posterior_global_pf: `ess_frac` must be in [0, 1].",
         call. = FALSE)

  if (!is.null(seed))
    warning("make_log_posterior_global_pf: a non-NULL `seed` freezes the ",
            "particle cloud across evaluations (marginal-likelihood variance ",
            "= 0). Fine for a one-off value, but it BREAKS the ",
            "pseudo-marginal validity of rwmh()/pmmh() -- use seed = NULL ",
            "for MCMC.", call. = FALSE)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  n_obs <- length(obs_vars)
  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  if (nrow(data) != n_obs) data <- t(data)
  Y <- data

  ## ---- FACTORY-TIME model-class validation --------------------------------
  ## solve_global() represents only a narrow class of models (every shock must
  ## sit in its own LINEAR AR(1) state process; see .global_shock_pairing()).
  ## Discovering that per draw would be useless twice over: the user learns
  ## nothing (the closure would simply return -Inf forever), and before the
  ## E1-B follow-up it did not even do that -- the name-based `eps_<state>`
  ## pairing silently dropped every unmatched shock and the filter returned a
  ## FINITE likelihood for a deterministic policy (-425 where the Kalman
  ## filter said 151 on the `y = beta*x + 0.2*x(+1) + u` oracle). Validate
  ## once, HERE, so an unsupported model is a build-time error naming the
  ## offending shock and equation.
  .gpf_validate_class(model, compiled, obs_vars)

  eval_one <- .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    solve_fn = function(model, compiled, sys_cache, ss, params, theta) {
      ## Order-1 solve first: it is the BK/stationarity gate (a projection
      ## solve on an explosive calibration wastes seconds before failing) AND
      ## supplies the stationary covariance used to initialise the particles.
      base <- .posterior_solve1(model, compiled, sys_cache, ss, params,
                                "spectral")
      if (is.null(base)) return(NULL)
      g <- tryCatch(
        solve_global(compiled, ss, params,
                     poly_degree  = poly_degree,
                     n_quad       = n_quad,
                     n_nodes      = n_nodes,
                     state_domain = state_domain,
                     tol          = solve_tol,
                     max_iter     = solve_max_iter,
                     verbose      = FALSE),
        error = function(e) NULL)
      if (is.null(g) || !isTRUE(g$converged)) return(NULL)
      list(dr = base$dr, sys = base$sys, global = g)
    },
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      g  <- sol$global
      dr <- sol$dr
      sn <- g$state_names

      obs_idx <- match(obs_vars, g$all_endo_names)
      if (anyNA(obs_idx)) return(NULL)

      Sigma_e <- .get_shock_cov(model, g$shock_names, params)
      Le <- tryCatch(t(chol(Sigma_e)), error = function(e)
        .tpf_psd_sqrt(Sigma_e))
      if (is.null(Le) || !all(is.finite(Le))) return(NULL)

      ## Stationary initialisation of the state LAG from the first-order
      ## solution -- the same P0 kalman_filter(lik_init = "stationary") uses,
      ## so the two likelihoods start from the same marginal.  dr's state
      ## ordering need not match the projection's, hence the explicit
      ## reordering rather than a positional assumption.
      TT <- dr$ghx[dr$state_idx, , drop = FALSE]
      RR <- dr$ghu[dr$state_idx, , drop = FALSE]
      P0 <- tryCatch(kf_stationary_init(TT, RR, Sigma_e),
                     error = function(e) NULL)
      if (is.null(P0) || !all(is.finite(P0))) return(NULL)
      pos <- match(sn, dr$state_vars)
      if (anyNA(pos)) return(NULL)
      P0 <- P0[pos, pos, drop = FALSE]
      L0 <- tryCatch(t(chol(P0)), error = function(e) .tpf_psd_sqrt(P0))
      if (is.null(L0) || !all(is.finite(L0))) return(NULL)

      s0 <- as.numeric(g$ss_vals[sn])
      if (anyNA(s0)) return(NULL)

      ll <- tryCatch(
        .global_pf_loglik(Y, g, obs_idx, s0, L0, Le, me_variance,
                          n_particles, ess_frac),
        error = function(e) -Inf)
      if (!is.finite(ll)) return(NULL)
      list(loglik = ll)
    },
    power          = power,
    needs_me_floor = FALSE,
    system_prior   = system_priors,
    pass_dots      = TRUE)

  function(theta) .with_local_seed(seed, eval_one(theta))
}


## Factory-time model-class check for make_log_posterior_global_pf().
##
## Runs the SAME .global_shock_pairing() contract solve_global() enforces per
## draw, once, at the baseline calibration, so an unsupported model errors at
## build time with the offending shock/equation named.  Also checks that the
## observables exist, because a bad obs_vars would otherwise surface as a
## silent per-draw rejection.
##
## If the baseline steady state does not converge we still validate: the
## pairing is read from the Jacobian's SPARSITY and its linear coefficients,
## which the initval-based guess resolves just as well.  Only a Jacobian that
## cannot be evaluated at all defers the check to the per-draw path (which
## rejects with -Inf, never a finite value).
#' @noRd
.gpf_validate_class <- function(model, compiled, obs_vars) {
  endo <- compiled$dynamic$endo_names
  miss <- setdiff(obs_vars, endo)
  if (length(miss))
    stop("make_log_posterior_global_pf: observable(s) not in the model: ",
         paste(miss, collapse = ", "), ".", call. = FALSE)

  params <- model$param_values
  ss     <- tryCatch(
    solve_steady_state(model, compiled, params, verbose = FALSE),
    error = function(e) NULL)
  ss_vals <- if (!is.null(ss) && isTRUE(ss$converged)) ss$ss else
    stats::setNames(rep(0, length(endo)), endo)

  state_names <- .global_state_names(compiled, endo)
  if (length(state_names) == 0L)
    stop("make_log_posterior_global_pf: the model has no state variables, so ",
         "there is no projection policy to filter with.", call. = FALSE)

  .global_shock_pairing(compiled, params, ss_vals, state_names,
                        context = "make_log_posterior_global_pf")
  invisible(TRUE)
}
