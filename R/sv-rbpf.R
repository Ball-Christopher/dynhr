## R/sv-rbpf.R
## --------------------------------------------------------------------------
## Rao-Blackwellized particle filter (RB-PF) for measurement-side stochastic
## volatility on the shocks of a linear DSGE, and its log-posterior factory.
##
## Design (verified in the scoping spikes; see the SV brief):
## Conditional on a volatility path {h_t}, an SV-on-shocks linear DSGE is
## exactly linear-Gaussian and its log-likelihood is a Kalman filter with
## shock_scale = exp(h_t / 2).  So we particle-filter ONLY the low-dimensional
## log-variance states h_t (dim = number of SV shocks, typically 3-7) and
## integrate the DSGE states analytically with kf_step() per particle.  Because
## the particle dimension is n_sv (not n_state), weights do not degenerate and
## no me_variance>0 crutch is required (unlike the TPF/OBC particle filters).
##
## The per-period marginal log-likelihood increment is a bootstrap PF estimate:
## propose h_t from its AR(1) prior, weight by the exact conditional Kalman
## likelihood exp(kf_step$ll), and take log-mean-exp.  Summed over t this is an
## UNBIASED estimator of the marginal likelihood -> valid inside pmmh() (the
## pseudo-marginal argument).  It is non-differentiable through resampling, so
## the estimator is gradient-free (PMMH/RWMH), as scoped.
##
## v1 scope: STATIONARY models only (shock_scale + diffuse init is a hard stop
## in kalman_filter; the same restriction applies here) and order-1 (linear)
## solutions.
##
## DEGENERATE PARTICLES PAST THE VOLATILITY OVERFLOW POINT (decided 2026-09-04,
## user; do not re-litigate without new evidence).  At extreme sigma_eta the
## filter returns a FINITE but meaningless value (-1.6e+132 nats at
## sigma_eta = 3000).  These are deliberately NOT converted to -Inf:
##
##   * Gating on rcond(Ft) -- the convention the other compiled KF kernels use
##     -- was MEASURED and is harmful twice over.  It moves the loglik by
##     18-48 nats at sigma_eta = 100, which the SBC's inv_gamma(0.35, 4) prior
##     genuinely reaches (max of 120k draws: 99.87); and it BIASES rather than
##     sanitises, because a near-singular Ft has a tiny determinant so
##     -0.5 log|F| is large POSITIVE -- the singular particles are
##     systematically the ones holding the cloud maximum, and gating them made
##     the sigma_eta = 300 estimate WORSE (-1974 -> -9.8e62).
##   * -Inf is strictly worse than finite garbage downstream.  The MH accept
##     step guards with `is.finite(log_alpha)` (R/sampler-rwmh.R:230), so a
##     chain that ever lands on -Inf can NEVER leave it (+Inf is not finite),
##     whereas from -1e132 a sane proposal gives a finite +1e132 and is
##     accepted immediately.
##   * The genuinely infeasible end already self-handles: once every particle
##     fails, the existing fail_period path returns -Inf with a warning.
##
## Full options table, the measured sweep and the falsifiers that WOULD
## justify revisiting are in memory/sv-rbpf-degenerate-particle-decision.md.
## Deferred there: a degenerate-particle COUNTER (visibility only, changes no
## number) -- worth adding next time src/sv_rbpf.cpp is open.
## --------------------------------------------------------------------------


#' TRUE when the compiled RB-PF sweep (sv_rbpf_loglik_cpp) is available.
#' @noRd
.HAS_RCPP_SV_RBPF <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("sv_rbpf_loglik_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}


#' Systematic resampling indices
#'
#' @param w Normalised weights (sum to 1), length N.
#' @return Integer vector of length N of resampled particle indices.
#' @noRd
.sv_systematic_resample <- function(w) {
  N <- length(w)
  positions <- (stats::runif(1) + 0:(N - 1L)) / N
  cumw <- cumsum(w)
  cumw[N] <- 1                                # guard against fp drift
  idx <- integer(N)
  i <- 1L
  for (j in seq_len(N)) {
    while (cumw[i] < positions[j]) i <- i + 1L
    idx[j] <- i
  }
  idx
}


#' RB-PF marginal log-likelihood for SV-on-shocks (internal)
#'
#' @param Y          n_obs x T observation matrix.
#' @param TT,ZZ,RR,DD State-space matrices (see \code{\link{kf_step}}).
#' @param Sigma_e    Baseline shock covariance, n_exo x n_exo.
#' @param d          Observation intercept, length n_obs (may be NULL).
#' @param P0         Baseline stationary prediction covariance (from
#'   \code{\link{kf_stationary_init}}).
#' @param sv_idx     Integer positions of the SV shocks in the n_exo shock
#'   vector (dr$exo_names order).
#' @param hyper      n_sv x 3 numeric matrix, columns (mu, rho, sigma_eta), rows
#'   aligned with \code{sv_idx}.
#' @param n_exo      Number of shocks.
#' @param n_particles Number of volatility particles.
#' @param me_diag    Optional length-n_obs baseline ME variances (NULL = none).
#'   TRUE i.i.d. observation noise: it enters F AND the Joseph state-covariance
#'   term (F4-A), matching \code{\link{kf_step}} and \code{\link{kalman_filter}}.
#' @return Scalar log-likelihood estimate, or \code{-Inf} if every particle is
#'   infeasible at some period.
#' @noRd
.sv_rbpf_loglik <- function(Y, TT, ZZ, RR, DD, Sigma_e, d, P0,
                            sv_idx, hyper, n_exo, n_particles,
                            me_diag = NULL) {
  ## Compiled full-sweep kernel (same RNG draw order as the R reference below,
  ## so a given set.seed produces the same particle cloud; parity-tested).
  if (.HAS_RCPP_SV_RBPF()) {
    res <- sv_rbpf_loglik_cpp(
      Y, TT, ZZ, RR, DD, Sigma_e,
      if (is.null(d)) numeric(nrow(Y)) else as.numeric(d),
      P0, as.integer(sv_idx),
      as.numeric(hyper[, "mu"]), as.numeric(hyper[, "rho"]),
      as.numeric(hyper[, "sigma_eta"]),
      as.integer(n_particles),
      if (is.null(me_diag)) numeric(0) else as.numeric(me_diag))
    ## The kernel REPORTS a fully infeasible cloud instead of warning from
    ## C++: Rcpp::warning() under options(warn = 2) becomes an R error, and
    ## R errors longjmp -- straight past ~RNGScope, leaving .Random.seed
    ## stale (and Armadillo buffers leaked) on exactly the runs a user asked
    ## to be strict. Warning here costs nothing and is longjmp-safe.
    if (isTRUE(as.integer(res$fail_period) > 0L))
      warning(sprintf(
        paste0("sv_rbpf_loglik: all %d particles failed at period %d ",
               "(non-PD forecast covariance or non-finite likelihood for ",
               "every particle) -- returning -Inf. This may indicate ",
               "linear-algebra failure rather than a genuine ",
               "zero-likelihood region."),
        as.integer(n_particles), as.integer(res$fail_period)),
        call. = FALSE)
    return(res$loglik)
  }

  N     <- as.integer(n_particles)
  n_sv  <- length(sv_idx)
  n_T   <- ncol(Y)
  n_s   <- nrow(TT)

  mu    <- hyper[, "mu"]
  rho   <- hyper[, "rho"]
  seta  <- hyper[, "sigma_eta"]

  ## Particle state: h (n_sv x N), s (n_s x N), P (list of N covariances).
  ## Init h from the stationary AR(1) marginal N(mu, seta^2/(1-rho^2)); (s,P)
  ## from the baseline stationary Kalman init (shared across particles).
  sd0 <- seta / sqrt(1 - rho^2)
  h_p <- matrix(stats::rnorm(n_sv * N, mean = rep(mu, N), sd = rep(sd0, N)),
                nrow = n_sv, ncol = N)
  s_p <- matrix(0, nrow = n_s, ncol = N)
  P_p <- rep(list(P0), N)

  scale_full <- rep(1, n_exo)                 # non-SV shocks stay at scale 1
  loglik <- 0

  for (t in seq_len(n_T)) {
    ## Propagate the log-variance (bootstrap proposal = AR(1) prior).
    h_p <- mu + rho * (h_p - mu) +
      matrix(stats::rnorm(n_sv * N, sd = rep(seta, N)), nrow = n_sv, ncol = N)
    scl <- exp(h_p / 2)                        # n_sv x N

    y_t   <- Y[, t]
    ll_t  <- numeric(N)
    s_new <- matrix(0, nrow = n_s, ncol = N)
    P_new <- vector("list", N)

    for (i in seq_len(N)) {
      scale_full[sv_idx] <- scl[, i]
      step <- kf_step(s_p[, i], P_p[[i]], y_t, TT, ZZ, RR, DD, Sigma_e,
                      scale = scale_full, d = d, me_diag = me_diag)
      if (is.null(step)) {                     # non-PD forecast cov -> weight 0
        ll_t[i] <- -Inf
        s_new[, i] <- s_p[, i]; P_new[[i]] <- P_p[[i]]
      } else {
        ll_t[i] <- step$ll
        s_new[, i] <- step$s; P_new[[i]] <- step$P
      }
    }

    ## Log-mean-exp period increment; bail if the whole cloud is infeasible.
    m <- max(ll_t)
    if (!is.finite(m)) return(-Inf)
    w_un <- exp(ll_t - m)
    loglik <- loglik + m + log(mean(w_un))

    ## Systematic resample proportional to weights.
    w   <- w_un / sum(w_un)
    idx <- .sv_systematic_resample(w)
    h_p <- h_p[, idx, drop = FALSE]
    s_p <- s_new[, idx, drop = FALSE]
    P_p <- P_new[idx]
  }

  loglik
}


## ---- Factory: make_log_posterior_sv_rbpf --------------------------------

#' Construct an SV-on-shocks RB-PF log-posterior function
#'
#' Returns a \code{function(theta)} that evaluates the log-posterior of a linear
#' DSGE with latent stochastic volatility on its shocks, using a
#' Rao-Blackwellized particle filter (see \code{\link{stochastic_volatility}}
#' for the model and \code{\link{kf_step}} for the analytic inner loop).
#'
#' The volatility hyperparameters (\code{mu}, \code{rho}, \code{sigma_eta} per
#' SV shock) are ordinary model parameters — declare them in
#' \code{estimated_params} to estimate them; \code{theta} flows through
#' \code{.apply_theta_to_params} exactly as for any other parameter.
#'
#' Because the marginal likelihood is a particle estimate, the closure is
#' intended for \code{\link{pmmh}()} (pseudo-marginal RWMH): use the default
#' \code{seed = NULL} so each evaluation draws a fresh, independent volatility
#' cloud — reusing one cloud across proposals would break the PMMH invariance.
#'
#' @param model       dynhr_mod from \code{\link{parse_mod}}.
#' @param data        Observation matrix (n_obs x T); columns are time periods.
#' @param prior_spec  Prior specification from \code{extract_prior_spec}.
#' @param obs_vars    Character vector of observed variable names.
#' @param compiled    dynhr_compiled from \code{\link{compile_model}}.
#' @param stochastic_volatility Optional call-level \code{sv_spec} (from
#'   \code{\link{stochastic_volatility}()}) overriding the model's
#'   \code{stochastic_volatility} block; \code{NULL} uses the model's block,
#'   \code{FALSE} disables it (which is an error here — an SV filter needs at
#'   least one SV shock).
#' @param n_particles Number of volatility particles (default 1000).
#' @param me_variance Optional baseline measurement-error variance (default 0;
#'   the RB-PF does not need it, but it is available for stochastically-singular
#'   observation blocks). It is TRUE i.i.d. observation noise -- the same law
#'   \code{\link{kalman_filter}} implements -- entering both the forecast
#'   covariance and the Joseph state-covariance update (F4-A).
#' @param seed        Integer RNG seed (default \code{NULL} = fresh randomness
#'   each call; required for valid PMMH).
#' @param power       Power-posterior tempering exponent zeta:
#'   \code{logpost = zeta * loglik + logprior}. \code{NULL} (default) resolves
#'   the global \code{power_posterior} option ONCE at factory time (falling
#'   back to 1), matching \code{make_log_posterior()} and
#'   \code{make_log_posterior_tpf()}.
#' @return A \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}.
#' @seealso \code{\link{stochastic_volatility}}, \code{\link{kf_step}},
#'   \code{\link{pmmh}}
#' @export
make_log_posterior_sv_rbpf <- function(model, data, prior_spec, obs_vars,
                                        compiled,
                                        stochastic_volatility = NULL,
                                        n_particles = 1000L,
                                        me_variance = 0,
                                        seed = NULL,
                                        power = NULL) {
  force(data); force(prior_spec); force(obs_vars)
  force(n_particles); force(me_variance); force(seed); force(power)

  ## Power-posterior exponent: explicit arg > global option > 1, resolved ONCE
  ## here (same contract as make_log_posterior() / make_log_posterior_tpf()).
  power <- .dynhr_opt("power_posterior", power, default = 1)
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power < 0)
    stop("make_log_posterior_sv_rbpf: 'power' must be a finite non-negative ",
         "scalar.", call. = FALSE)
  power <- as.numeric(power)

  model <- .resolve_stochastic_volatility(model, stochastic_volatility)
  sv_spec <- model$stochastic_volatility
  if (is.null(sv_spec) || is.null(sv_spec$sv) || nrow(sv_spec$sv) == 0L)
    stop("make_log_posterior_sv_rbpf: no stochastic_volatility entries. ",
         "Declare a stochastic_volatility block or pass one via the ",
         "'stochastic_volatility' argument.", call. = FALSE)

  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance < 0)
    stop("make_log_posterior_sv_rbpf: 'me_variance' must be a non-negative ",
         "finite scalar.", call. = FALSE)

  ## Hard stop on missing observations: kf_step has no missing-data handling,
  ## so an NA would silently poison every particle's weight and collapse the
  ## whole marginal likelihood to -Inf. Fail loud instead of silently wrong.
  ## (The exact Gaussian kalman_filter DOES skip NA dimensions; SV-on-shocks
  ## missing-data support is a documented v1 limitation.)
  if (anyNA(data))
    stop("make_log_posterior_sv_rbpf: the observation matrix contains NA. ",
         "The SV-on-shocks RB-PF does not support missing observations in v1 ",
         "(the per-particle Kalman step has no missing-data recursion). ",
         "Drop/interpolate the missing periods, or use the exact Gaussian ",
         "likelihood (which skips NA dimensions).", call. = FALSE)

  ## PMMH validity: a fixed seed makes the particle filter DETERMINISTIC (the
  ## same volatility cloud every call, marginal-likelihood SD = 0), which breaks
  ## the pseudo-marginal invariance of pmmh(). Warn; the default seed = NULL is
  ## correct. (Mirrors the TPF seed convention.)
  if (!is.null(seed))
    warning("make_log_posterior_sv_rbpf: a non-NULL 'seed' makes the RB-PF ",
            "deterministic across evaluations (marginal-likelihood variance = 0). ",
            "This is fine for a one-off likelihood value but BREAKS pmmh() ",
            "pseudo-marginal validity -- use seed = NULL (the default) for MCMC.",
            call. = FALSE)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  if (is.null(dim(data))) data <- matrix(data, nrow = length(obs_vars))
  n_obs   <- length(obs_vars)
  me_diag <- if (me_variance > 0) rep(me_variance, n_obs) else NULL

  ## Inner evaluator over the shared closure builder (R/posterior-closure.R).
  ## Kept separate from the returned closure because .with_local_seed() forces
  ## its `expr` as a promise, and a promise containing a top-level `return()`
  ## cannot be forced from another frame. `pass_dots = TRUE` returns the raw
  ## `function(theta, ...)` inner; the wrapper below restores the public
  ## `function(theta)` signature.
  eval_one <- .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      dr <- sol$dr
      ## ---- State-space matrices (identical extraction to kalman_filter) ---
      endo <- dr$endo_names; exo <- dr$exo_names
      state_idx <- dr$state_idx
      obs_idx   <- match(obs_vars, endo)
      if (any(is.na(obs_idx))) return(NULL)

      TT <- dr$ghx[state_idx, , drop = FALSE]
      RR <- dr$ghu[state_idx, , drop = FALSE]
      ZZ <- dr$ghx[obs_idx,   , drop = FALSE]
      DD <- dr$ghu[obs_idx,   , drop = FALSE]
      d  <- dr$ys[obs_vars]
      Sigma_e <- .get_shock_cov(model, exo, params)

      ## ---- Resolve SV spec against dr shock order + params ----------------
      sv_idx <- .sv_shock_index(sv_spec, exo)
      hyper  <- .sv_resolve_hyperparams(sv_spec, params)
      if (is.character(hyper)) return(NULL)   # domain-infeasible hyperparams

      P0 <- tryCatch(kf_stationary_init(TT, RR, Sigma_e),
                     error = function(e) NULL)
      if (is.null(P0) || !all(is.finite(P0))) return(NULL)

      Y <- if (nrow(data) != n_obs) t(data) else data

      loglik <- tryCatch(
        .sv_rbpf_loglik(Y, TT, ZZ, RR, DD, Sigma_e, d, P0,
                        sv_idx, hyper, length(exo), n_particles,
                        me_diag = me_diag),
        error = function(e) -Inf)
      if (!is.finite(loglik)) return(NULL)
      list(loglik = loglik)
    },
    power          = power,
    needs_me_floor = FALSE,
    pass_dots      = TRUE)

  function(theta) {
    ## Seeding is LOCAL: a non-NULL seed makes THIS evaluation reproducible
    ## without resetting the caller's global .Random.seed (A1) -- an outer
    ## sampler that draws around the likelihood keeps its own stream.
    .with_local_seed(seed, eval_one(theta))
  }
}
