## R/sv-rbpf.R
## --------------------------------------------------------------------------
## Rao-Blackwellized particle filter (RB-PF) for measurement-side stochastic
## volatility on the shocks of a linear DSGE, and its log-posterior factory.
##
## Design (verified in the scoping spikes; see the SV derivation):
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
#'   term, matching \code{\link{kf_step}} and \code{\link{kalman_filter}}.
#' @param n_burn     Pre-sample periods (default 0). Each particle's h and P
#'   are advanced \code{n_burn} unobserved periods
#'   (\eqn{P \leftarrow TT P TT' + RR \Sigma_e(h) RR'}) before period 1, so
#'   that P is the covariance of \eqn{s_0} conditional on the particle's own
#'   volatility history. The factory passes \code{\link{.sv_burn_in_length}}
#'   with \code{P0 = .sv_stationary_P0()}; 0 reproduces the shared-P0 start
#'   (used by the filter-mechanics tests and the grid/Chan oracles).
#' @return Scalar log-likelihood estimate, or \code{-Inf} if every particle is
#'   infeasible at some period.
#' @noRd
.sv_rbpf_loglik <- function(Y, TT, ZZ, RR, DD, Sigma_e, d, P0,
                            sv_idx, hyper, n_exo, n_particles,
                            me_diag = NULL, n_burn = 0L) {
  n_burn <- as.integer(n_burn)
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
      if (is.null(me_diag)) numeric(0) else as.numeric(me_diag),
      n_burn)
    ## The kernel REPORTS a fully infeasible cloud instead of warning from
    ## C++: Rcpp::warning() under options(warn = 2) becomes an R error, and
    ## R errors longjmp -- straight past ~RNGScope, leaving .Random.seed
    ## stale (and Armadillo buffers leaked) on exactly the runs a user asked
    ## to be strict. Warning here costs nothing and is longjmp-safe.
    if (isTRUE(as.integer(res$fail_period) > 0L))
      .dynhr_warn(sprintf(
        paste0("sv_rbpf_loglik: all %d particles failed at period %d ",
               "(non-PD forecast covariance or non-finite likelihood for ",
               "every particle) -- returning -Inf. This may indicate ",
               "linear-algebra failure rather than a genuine ",
               "zero-likelihood region."),
        as.integer(n_particles), as.integer(res$fail_period)),
        call. = FALSE)
    return(.sv_rbpf_finish(res$loglik, res$n_degenerate, n_particles,
                           ncol(Y)))
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

  ## Pre-sample burn-in: per-particle P becomes Var(s_0 | own h history).
  ## Same draw order as src/sv_rbpf.cpp.
  for (b in seq_len(n_burn)) {
    h_p <- mu + rho * (h_p - mu) +
      matrix(stats::rnorm(n_sv * N, sd = rep(seta, N)), nrow = n_sv, ncol = N)
    for (i in seq_len(N)) {
      scale_full[sv_idx] <- exp(h_p[, i] / 2)
      Pb <- TT %*% P_p[[i]] %*% t(TT) +
        RR %*% (Sigma_e * outer(scale_full, scale_full)) %*% t(RR)
      P_p[[i]] <- (Pb + t(Pb)) * 0.5
    }
  }
  scale_full[] <- 1

  loglik <- 0
  n_degenerate <- 0

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
        n_degenerate <- n_degenerate + 1
        ll_t[i] <- -Inf
        s_new[, i] <- s_p[, i]; P_new[[i]] <- P_p[[i]]
      } else {
        ll_t[i] <- step$ll
        s_new[, i] <- step$s; P_new[[i]] <- step$P
      }
    }

    ## Log-mean-exp period increment; bail if the whole cloud is infeasible.
    m <- max(ll_t)
    if (!is.finite(m))
      return(.sv_rbpf_finish(-Inf, n_degenerate, n_particles, n_T))
    w_un <- exp(ll_t - m)
    loglik <- loglik + m + log(mean(w_un))

    ## Systematic resample proportional to weights.
    w   <- w_un / sum(w_un)
    idx <- .sv_systematic_resample(w)
    h_p <- h_p[, idx, drop = FALSE]
    s_p <- s_new[, idx, drop = FALSE]
    P_p <- P_new[idx]
  }

  .sv_rbpf_finish(loglik, n_degenerate, n_particles, n_T)
}


#' Unconditional stationary state covariance under SV-on-shocks (internal)
#'
#' The RB-PF starts every particle from the same Gaussian
#' \eqn{s_0 \sim N(0, P_0)}. Under the SV model
#' \eqn{s_t = TT s_{t-1} + RR\, D_t e_t}, \eqn{D_t = diag(\exp(h_t/2))}
#' (non-SV shocks at scale 1), with \eqn{e_t} and the \eqn{h} paths mutually
#' independent and serially independent shocks,
#' \eqn{Var(s_0) = \sum_k TT^k RR\, E[D \Sigma_e D]\, RR' TT^{k\prime}}: the
#' Lyapunov solution with \eqn{\Sigma_e} replaced by
#' \eqn{\Sigma_e \circ W}, \eqn{W_{ii} = E[e^{h_i}] = e^{\mu_i + v_i/2}},
#' \eqn{W_{ij} = E[e^{h_i/2}] E[e^{h_j/2}]} (\eqn{i \ne j}, independent
#' AR(1) log-variances), \eqn{E[e^{h_i/2}] = e^{\mu_i/2 + v_i/8}},
#' \eqn{v_i = \sigma_{\eta,i}^2 / (1 - \rho_i^2)}. Using the baseline
#' \eqn{\Sigma_e} (the pre-0.9.4.77 behaviour) is this with \eqn{W = 1}, i.e.
#' it pretended \eqn{\mu = 0} and \eqn{\sigma_\eta = 0}.
#'
#' @param TT,RR,Sigma_e State transition, shock impact, baseline shock cov.
#' @param sv_idx Positions of the SV shocks in the shock vector.
#' @param hyper n_sv x 3 matrix with columns mu, rho, sigma_eta.
#' @param n_exo Number of shocks.
#' @return n_s x n_s covariance.
#' @noRd
.sv_stationary_P0 <- function(TT, RR, Sigma_e, sv_idx, hyper, n_exo) {
  v  <- hyper[, "sigma_eta"]^2 / (1 - hyper[, "rho"]^2)
  m1 <- rep(1, n_exo); m2 <- rep(1, n_exo)
  m1[sv_idx] <- exp(hyper[, "mu"] / 2 + v / 8)   # E[exp(h/2)]
  m2[sv_idx] <- exp(hyper[, "mu"] + v / 2)       # E[exp(h)]
  W <- outer(m1, m1)
  diag(W) <- m2
  kf_stationary_init(TT, RR, Sigma_e * W)
}


#' Pre-sample burn-in length for the RB-PF (internal)
#'
#' Starting the burn-in at the unconditional \code{.sv_stationary_P0()}, the
#' error left in each particle's \eqn{Var(s_0 \mid h)} after \eqn{B} steps is
#' \eqn{TT^B \Delta TT^{B\prime}}: a fraction \eqn{r^{2B}} of the
#' UNCONDITIONAL variance, with \eqn{r} the spectral radius of \code{TT}.
#' That variance is an average of \eqn{\exp(h)}, which sits far above a
#' typical particle's own when the log-volatility has a large stationary
#' variance \eqn{v = \sigma_\eta^2 / (1 - \rho^2)} (mean scale
#' \eqn{\exp(\mu + v/2)} against \eqn{\exp(\mu - v/2)} for a quiet particle:
#' a factor \eqn{\exp(v)}). \eqn{B} is therefore chosen so that
#' \eqn{r^{2B} \exp(v_{max}) \le} \code{tol}, with \eqn{v_{max}} the largest
#' \eqn{v} across the SV shocks. (A rule that ignored \eqn{v} left the
#' log-likelihood 0.2 nats high at \eqn{v = 18}.) \eqn{B} is at least 1 (with
#' \eqn{TT = 0} one step is exact: \eqn{Var(s_0 \mid h_0) = RR \Sigma_e(h_0)
#' RR'}) and at most \code{max_burn}; hitting the cap is announced once
#' (classed message \code{dynhr_sv_rbpf_burn_capped}) with the residual.
#' @param TT State transition.
#' @param tol Relative residual target.
#' @param max_burn Cap on the burn-in length.
#' @param hyper SV hyper-parameter matrix (columns \code{rho},
#'   \code{sigma_eta}); \code{NULL} takes \eqn{v_{max} = 0}.
#' @return Integer burn-in length.
#' @noRd
.sv_burn_in_length <- function(TT, tol = 1e-4, max_burn = 1000L, hyper = NULL) {
  if (length(TT) == 0L) return(0L)
  r <- max(Mod(eigen(TT, only.values = TRUE)$values))
  if (!is.finite(r) || r < 1e-12) return(1L)
  if (r >= 1) return(as.integer(max_burn))   # non-stationary: P0 already fails
  v_max <- 0
  if (!is.null(hyper) && nrow(hyper) > 0L) {
    v <- hyper[, "sigma_eta"]^2 / (1 - hyper[, "rho"]^2)
    v_max <- max(v[is.finite(v)], 0)
    if (any(!is.finite(v))) return(as.integer(max_burn))
  }
  ## In logs: exp(-v_max) underflows for v_max beyond ~700.
  log_resid <- function(B) 2 * B * log(r) + v_max     # log of r^(2B) exp(v_max)
  B_need <- max(1, ceiling((log(tol) - v_max) / (2 * log(r))))
  if (B_need > max_burn) {
    .dynhr_inform(
      sprintf(paste0("sv_rbpf: state persistence %.6f with log-volatility ",
                     "variance %.3g needs a %.0f-period pre-sample burn-in for ",
                     "the initial-state law; capped at %d (residual relative ",
                     "error 1e%+.0f in Var(s_0 | h))."),
              r, v_max, B_need, as.integer(max_burn),
              log_resid(max_burn) / log(10)),
      once = TRUE, key = "sv_rbpf_burn_capped",
      class = "dynhr_sv_rbpf_burn_capped")
    return(as.integer(max_burn))
  }
  as.integer(B_need)
}


#' Attach and report the degenerate particle-step count (internal)
#'
#' A degenerate particle-step is one that took a weight-0 path (non-finite,
#' non-PD or non-invertible forecast covariance, or a non-finite
#' log-likelihood). The count changes no number: the log-likelihood is
#' returned exactly as computed, with the count as attribute
#' \code{n_degenerate}. A positive count is announced once per run through a
#' classed message (\code{dynhr_sv_rbpf_degenerate}) so a PMMH/SMC loop is not
#' flooded.
#'
#' @param loglik Scalar log-likelihood estimate.
#' @param n_degenerate Number of weight-0 particle-steps.
#' @param n_particles,n_T Particle count and number of periods (for the
#'   message's denominator).
#' @return \code{loglik} with attribute \code{n_degenerate}.
#' @noRd
.sv_rbpf_finish <- function(loglik, n_degenerate, n_particles, n_T) {
  n_degenerate <- as.numeric(n_degenerate)
  if (isTRUE(n_degenerate > 0))
    .dynhr_inform(
      sprintf(paste0("sv_rbpf_loglik: %.0f of %.0f particle-steps (%.2g%%) ",
                     "degenerated (non-finite or non-PD forecast covariance) ",
                     "and received weight zero; the likelihood values are ",
                     "kept finite where the cloud survives but may not be ",
                     "meaningful there. Check the volatility ",
                     "hyperparameters (sigma_eta in particular)."),
              n_degenerate, as.numeric(n_particles) * as.numeric(n_T),
              100 * n_degenerate /
                (as.numeric(n_particles) * as.numeric(n_T))),
      once = TRUE, key = "sv_rbpf_degenerate",
      class = "dynhr_sv_rbpf_degenerate")
  attr(loglik, "n_degenerate") <- n_degenerate
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
#' @param data        Observation matrix, \code{T x n_obs} (rows are periods,
#'   columns the \code{obs_vars}); oriented by the rule documented in
#'   \code{\link{make_log_posterior_tpf}} (names decide; an unnamed
#'   \code{n_obs x T} or square sample is refused).
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
#'   covariance and the Joseph state-covariance update.
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

  ## A per-observable vector: refused with a classed error; an all-equal
  ## vector is the scalar.
  if (length(me_variance) > 1L)
    me_variance <- .kf_me_variance(me_variance, obs_vars,
                                   "make_log_posterior_sv_rbpf",
                                   allow_vector = FALSE)
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance < 0)
    stop("make_log_posterior_sv_rbpf: 'me_variance' must be a non-negative ",
         "finite scalar.", call. = FALSE)

  ## Particle count: a positive whole number. 0 / negative / NA used to reach
  ## the kernel and come back as a silent -Inf at every theta (and 2.5 was
  ## truncated, "10" coerced).
  n_particles <- .check_particle_count(n_particles, "make_log_posterior_sv_rbpf",
                                       "n_particles")

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
    .dynhr_warn("make_log_posterior_sv_rbpf: a non-NULL 'seed' makes the RB-PF ",
            "deterministic across evaluations (marginal-likelihood variance = 0). ",
            "This is fine for a one-off likelihood value but BREAKS pmmh() ",
            "pseudo-marginal validity -- use seed = NULL (the default) for MCMC.",
            call. = FALSE)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  ## `data` is T x n_obs (rows are periods), oriented by the package-wide
  ## rule (R/time-orientation.R); the RB-PF runs on its transpose
  ## (n_obs x T, one column per period).
  data <- t(.orient_time_matrix(data, obs_vars, "data",
                                "make_log_posterior_sv_rbpf"))
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

      ## Initial state law: the SV-unconditional stationary covariance
      ## (Sigma_e weighted by E[exp(h/2) exp(h'/2)], see .sv_stationary_P0),
      ## then a per-particle pre-sample burn-in so each particle starts from
      ## Var(s_0 | its own volatility history). Pre-0.9.4.77 every particle
      ## started from the BASELINE-Sigma_e Lyapunov P0 -- as if mu = 0 and
      ## sigma_eta = 0 -- a loglik bias of ~0.2-0.8 nats that moved with mu.
      P0 <- tryCatch(.sv_stationary_P0(TT, RR, Sigma_e, sv_idx, hyper,
                                       length(exo)),
                     error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(P0) || !all(is.finite(P0))) return(NULL)
      n_burn <- .sv_burn_in_length(TT, hyper = hyper)

      Y <- data                               # n_obs x T (oriented at build)

      loglik <- tryCatch(
        .sv_rbpf_loglik(Y, TT, ZZ, RR, DD, Sigma_e, d, P0,
                        sv_idx, hyper, length(exo), n_particles,
                        me_diag = me_diag, n_burn = n_burn),
        error = function(e) .dynhr_reraise_bug(e, -Inf))
      loglik <- as.numeric(loglik)            # drop the diagnostic attribute
      if (!is.finite(loglik)) return(NULL)
      list(loglik = loglik)
    },
    power          = power,
    needs_me_floor = FALSE,
    pass_dots      = TRUE)

  function(theta) {
    ## Seeding is LOCAL: a non-NULL seed makes THIS evaluation reproducible
    ## without resetting the caller's global .Random.seed -- an outer
    ## sampler that draws around the likelihood keeps its own stream.
    .with_local_seed(seed, eval_one(theta))
  }
}
