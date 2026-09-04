## R/obc-ppf.R
## --------------------------------------------------------------------------
## OBC Bootstrap Piecewise Particle Filter (PPF).
##
## Implements the bootstrap variant of the PPF for OCC-bin piecewise-linear
## models. Each particle independently resolves its own OBC regime at each
## period using the shared per-draw regime_cache.
##
## References:
##   Aruoba, Cuba-Borda, Higa-Flores, Schorfheide & Villalvazo (2021, RED)
##   Dynare 7.0 PPF documentation
##
## Provides:
##   ppf_likelihood()              -- per-evaluation bootstrap PF loglik
##   make_log_posterior_obc_ppf()  -- closure factory (mirrors PKF version)
## --------------------------------------------------------------------------


## ============================================================================
## Per-period bootstrap PPF step (single period, vectorised over particles)
## ============================================================================

#' Bootstrap PPF: process one time period
#'
#' Each particle draws eps_t^i ~ N(0, Sigma_e) and resolves its OBC regime.
#' Particle weight = Gaussian p(y_t | s_{t-1}^i, eps_t^i, regime^i).
#' Accumulates log_lik_contrib BEFORE resampling; resamples at end of period.
#'
#' @param particles n_state x N matrix of particles (s_{t-1}^i)
#' @param y_t       length-n_obs observation vector (may contain NAs)
#' @param L_e       n_exo x n_exo lower-triangular Cholesky of Sigma_e
#' @param Sigma_e   n_exo x n_exo shock covariance
#' @param dr_slack  Slack-regime DecisionRules
#' @param regime_cache R environment of per-regime policies
#' @param sys       System matrices (for lazy regime building)
#' @param specs     OBC spec list
#' @param obs_idx   Integer vector of observable indices
#' @param d_obs     length-n_obs observable steady-state mean
#' @param me_variance Scalar measurement error variance (must be > 0)
#' @return list(particles = n_state x N updated, log_lik_contrib = scalar)
#' @noRd
.ppf_run_period <- function(particles, y_t, L_e, Sigma_e,
                             dr_slack, regime_cache, sys, specs,
                             obs_idx, d_obs, me_variance) {

  N       <- ncol(particles)
  n_state <- nrow(particles)
  n_exo   <- ncol(L_e)

  ## -- Step 1: draw shocks for all particles --------------------------------
  shocks <- L_e %*% matrix(rnorm(n_exo * N), nrow = n_exo)   # n_exo x N

  ## -- Step 2: resolve regime and compute observation likelihoods -----------
  ## For each particle: check OBC regime from (s_{t-1}^i, eps^i)
  log_w <- numeric(N)

  for (i in seq_len(N)) {
    s_i   <- particles[, i]
    eps_i <- shocks[, i]

    ## Regime check: use pkf_check_binding with backward state = particle state
    ## (no backward smoothing in PPF; we use the particle's own s_{t-1})
    bind_flags <- pkf_check_binding(s_i, eps_i, specs, dr_slack)
    regime_idx <- obc_regime_idx(bind_flags)

    ## Ensure this regime is in the cache
    if (!exists(as.character(regime_idx), envir = regime_cache, inherits = FALSE))
      obc_ensure_policy(regime_idx, regime_cache, sys, dr_slack, specs, obs_idx)

    pol <- get(as.character(regime_idx), envir = regime_cache, inherits = FALSE)

    ## Observation prediction: y_pred = ZZ * s_{t-1} + DD * eps + d + c_obs
    y_pred <- drop(pol$ZZ %*% s_i) + drop(pol$DD %*% eps_i) + d_obs + pol$c_obs

    ## Bootstrap PPF weight: p(y_t | s_{t-1}^i, eps^i, regime^i)
    ## In the bootstrap PF, eps^i is drawn from the prior N(0, Sigma_e).
    ## Given eps^i and s_{t-1}^i, y_t is predicted exactly up to measurement
    ## error:  y_t = y_pred + noise_t,  noise_t ~ N(0, me_variance * I)
    ## So the weight is purely:
    ##   p(y_t | s_{t-1}^i, eps^i, regime^i) = N(y_t; y_pred, me_variance * I)
    F_i   <- me_variance * diag(length(y_pred))
    innov <- y_t - y_pred

    ## Handle missing observations
    obs_ok <- which(!is.na(innov))
    if (length(obs_ok) == 0L) {
      log_w[i] <- 0   # no information; weight = 1
    } else {
      innov_ok  <- innov[obs_ok]
      n_ok      <- length(obs_ok)
      ## F is diagonal: me_variance * I, so log_det = n_ok * log(me_variance)
      ## and quad form = sum(innov^2) / me_variance
      log_w[i]  <- -0.5 * n_ok * log(2 * pi) -
                   0.5 * n_ok * log(me_variance) -
                   0.5 * sum(innov_ok^2) / me_variance
    }

    ## Propagate state (TPF lesson: propagate AFTER weighting)
    ## s_t^i = TT * s_{t-1}^i + RR * eps^i + c_state
    shocks[, i] <- eps_i
  }

  ## -- Step 3: accumulate log_lik_contrib BEFORE resampling -----------------
  ## log p(y_t | Y_{1:t-1}) = log(mean exp(log_w))
  log_lik_contrib <- .smc_log_sum_exp(log_w) - log(N)

  ## -- Step 4: normalize and systematic resample ----------------------------
  log_w_c <- log_w - max(log_w)
  w_norm  <- exp(log_w_c)
  w_norm  <- w_norm / sum(w_norm)

  idx       <- .smc_systematic_resample(w_norm, N)
  shocks    <- shocks[, idx, drop = FALSE]

  ## -- Step 5: propagate resampled particles --------------------------------
  new_particles <- matrix(0, n_state, N)
  for (i in seq_len(N)) {
    orig_i <- idx[i]
    s_i    <- particles[, orig_i]
    eps_i  <- shocks[, i]

    ## Re-resolve regime after resampling (uses same eps, same s)
    bind_flags <- pkf_check_binding(s_i, eps_i, specs, dr_slack)
    regime_idx <- obc_regime_idx(bind_flags)
    if (!exists(as.character(regime_idx), envir = regime_cache, inherits = FALSE))
      obc_ensure_policy(regime_idx, regime_cache, sys, dr_slack, specs, obs_idx)
    pol <- get(as.character(regime_idx), envir = regime_cache, inherits = FALSE)

    new_particles[, i] <- drop(pol$TT %*% s_i) + drop(pol$RR %*% eps_i) + pol$c_state
  }

  list(particles = new_particles, log_lik_contrib = log_lik_contrib)
}


## ============================================================================
## Per-period COPF step (conditionally-optimal proposal)
## ============================================================================

#' COPF: process one time period using the conditionally-optimal Gaussian proposal
#'
#' For each particle, draws eps_t^i from the posterior N(mu_r^i, Omega_r) given
#' (y_t, s_{t-1}^i, regime r).  Weight = marginal N(v_r^i ; 0, F_r).
#' Fallback to bootstrap draw + bootstrap weight when verified regime != guessed.
#'
#' Accumulates log_lik_contrib BEFORE resampling (same as bootstrap PPF).
#'
#' @param particles     n_state x N matrix of particles (s_{t-1}^i)
#' @param y_t           length-n_obs observation vector (may contain NAs)
#' @param L_e           n_exo x n_exo lower-triangular Cholesky of Sigma_e
#' @param Sigma_e       n_exo x n_exo shock covariance
#' @param Sigma_e_inv   n_exo x n_exo inverse of Sigma_e (pre-computed once)
#' @param dr_slack      Slack-regime DecisionRules
#' @param regime_cache  R environment of per-regime policies (with COPF fields)
#' @param sys           System matrices (for lazy regime building)
#' @param specs         OBC spec list
#' @param obs_idx       Integer vector of observable indices
#' @param d_obs         length-n_obs observable steady-state mean
#' @param me_variance   Scalar measurement error variance (must be > 0)
#' @param copf_args     list(Sigma_e, Sigma_e_inv, me_variance) for cache builds
#' @param prev_regimes  Integer vector length N: verified regime of each
#'   particle's resampled ancestor from the previous period.  NULL or a
#'   zero-vector triggers all-slack guesses (t=1 behaviour).
#' @return list(particles = n_state x N updated, log_lik_contrib = scalar,
#'             n_fallback = integer count of bootstrap fallbacks,
#'             verified_regimes = integer vector length N of verified regimes
#'               after propagation, in resampled order — pass as prev_regimes
#'               to the next period call)
#' @noRd
.copf_run_period <- function(particles, y_t, L_e, Sigma_e, Sigma_e_inv,
                              dr_slack, regime_cache, sys, specs,
                              obs_idx, d_obs, me_variance, copf_args,
                              prev_regimes = NULL,
                              U_copf = NULL) { # CPM: list(z_copf = n_exo x N, z_fallback = n_exo x N)
                                               #      NULL -> draw fresh from RNG

  N       <- ncol(particles)
  n_state <- nrow(particles)
  n_exo   <- ncol(L_e)

  ## Pre-draw standard normals for COPF (consumed in mu + L_Omega * u)
  ## and fallback draws (used only on fallback particles).
  ## CPM: when U_copf is supplied, use its pre-drawn matrices (deterministic U structure).
  ## Both z_copf and z_fallback are always drawn/supplied; only the per-particle
  ## SELECTION (copf vs fallback) is data-dependent. The U structure itself is fixed.
  z_copf     <- if (!is.null(U_copf)) U_copf$z_copf     else matrix(rnorm(n_exo * N), nrow = n_exo)
  z_fallback <- if (!is.null(U_copf)) U_copf$z_fallback  else matrix(rnorm(n_exo * N), nrow = n_exo)

  ## Ancestor guesses: use prev_regimes if supplied and non-trivial,
  ## otherwise fall back to all-slack (regime 0) for all particles.
  ## This implements variance reduction: a particle whose ancestor was in the
  ## binding regime proposes from that regime's COPF distribution, not the
  ## slack regime's, slashing the mismatch rate during binding stretches.
  use_ancestor <- !is.null(prev_regimes) && length(prev_regimes) == N &&
                    any(prev_regimes != 0L)
  ## Pre-compute unique guessed regimes so we can batch-ensure their cache
  ## entries before the per-particle loop (avoids repeated env lookups in the
  ## common case where only 1-2 distinct regimes are present in the cloud).
  r_guesses <- if (use_ancestor) as.integer(prev_regimes) else rep(0L, N)
  for (rg in unique(r_guesses)) {
    if (!exists(as.character(rg), envir = regime_cache, inherits = FALSE))
      obc_ensure_policy(rg, regime_cache, sys, dr_slack, specs, obs_idx, copf_args)
  }

  log_w    <- numeric(N)
  shocks   <- matrix(0, n_exo, N)
  n_fallback <- 0L

  for (i in seq_len(N)) {
    s_i <- particles[, i]

    ## -- COPF path: guess regime from ancestor's verified regime
    ## (or all-slack at t=1 / when prev_regimes is NULL)
    r_guess <- r_guesses[i]

    ## Cache entry guaranteed to exist from the batch-ensure above
    pol_g <- get(as.character(r_guess), envir = regime_cache, inherits = FALSE)

    ## Innovation for guessed regime
    obs_offset_g <- d_obs + pol_g$c_obs
    v_g  <- y_t - drop(pol_g$ZZ %*% s_i) - obs_offset_g

    ## Handle NAs: use only observed components for COPF draw
    obs_ok <- which(!is.na(v_g))
    n_ok   <- length(obs_ok)

    if (n_ok == 0L) {
      ## All missing: draw from prior, weight 1
      eps_i    <- drop(L_e %*% z_copf[, i])
      shocks[, i] <- eps_i
      ## Propagate with guessed regime (regime 0)
      log_w[i] <- 0
    } else {
      ## Build COPF quantities restricted to observed components
      DD_ok   <- pol_g$DD[obs_ok, , drop = FALSE]
      v_ok    <- v_g[obs_ok]

      ## Omega_r_inv and Omega_r for partial obs case
      ## Use stored Omega if all obs present; recompute if partial obs
      if (n_ok == length(v_g) && !is.null(pol_g$Omega) && !is.null(pol_g$L_Omega)) {
        Omega_g   <- pol_g$Omega
        L_Omega_g <- pol_g$L_Omega
      } else {
        ## Partial obs: recompute Omega for the observed subset
        Omega_inv_g <- crossprod(DD_ok) / me_variance + Sigma_e_inv
        ch_Oi_g     <- tryCatch(chol(Omega_inv_g), error = function(e2) NULL)
        if (!is.null(ch_Oi_g)) {
          Omega_g <- chol2inv(ch_Oi_g)
          ch_O_g  <- tryCatch(chol(Omega_g), error = function(e2) NULL)
        } else {
          ch_O_g  <- NULL
        }
        if (is.null(ch_O_g)) {
          ## Degenerate Omega: fall back to a prior draw. The choice to fall
          ## back is deterministic given (s_i, y_t), so the bootstrap weight
          ## is a valid importance weight -- but the measurement density must
          ## use the regime actually implied by (s_i, eps_i), not the guess.
          eps_i <- drop(L_e %*% z_fallback[, i])
          shocks[, i] <- eps_i
          r_fb_flags <- pkf_check_binding(s_i, eps_i, specs, dr_slack)
          r_fb       <- obc_regime_idx(r_fb_flags)
          if (!exists(as.character(r_fb), envir = regime_cache, inherits = FALSE))
            obc_ensure_policy(r_fb, regime_cache, sys, dr_slack, specs, obs_idx, copf_args)
          pol_fb  <- get(as.character(r_fb), envir = regime_cache, inherits = FALSE)
          y_pred  <- drop(pol_fb$ZZ %*% s_i) + drop(pol_fb$DD %*% eps_i) +
                       d_obs + pol_fb$c_obs
          innov   <- y_t - y_pred
          innov_ok2 <- innov[which(!is.na(innov))]
          n_ok2 <- length(innov_ok2)
          log_w[i] <- if (n_ok2 == 0L) 0 else
            -0.5 * n_ok2 * log(2 * pi) - 0.5 * n_ok2 * log(me_variance) -
            0.5 * sum(innov_ok2^2) / me_variance
          n_fallback <- n_fallback + 1L
          next
        }
        L_Omega_g <- t(ch_O_g)
      }

      ## COPF draw: eps_i = mu_r^i + L_Omega_r * z
      mu_g  <- drop(Omega_g %*% (t(DD_ok) %*% v_ok / me_variance))
      eps_i <- mu_g + drop(L_Omega_g %*% z_copf[, i])

      ## -- Verify regime --
      r_verify_flags <- pkf_check_binding(s_i, eps_i, specs, dr_slack)
      r_verify       <- obc_regime_idx(r_verify_flags)

      if (r_verify == r_guess) {
        ## -- COPF accepted: weight = N(v_g ; 0, F_r) -----------------------
        ## Use stored F_inv and log_det_F if all obs present; otherwise recompute
        if (n_ok == length(v_g) && !is.null(pol_g$F_inv)) {
          F_inv_g     <- pol_g$F_inv
          log_det_F_g <- pol_g$log_det_F
        } else {
          F_g         <- DD_ok %*% Sigma_e %*% t(DD_ok) + me_variance * diag(n_ok)
          ch_Fg       <- tryCatch(chol(F_g), error = function(e2) NULL)
          if (is.null(ch_Fg)) {
            ## F Cholesky failed: weight -Inf
            shocks[, i] <- eps_i
            log_w[i]    <- -Inf
            next
          }
          F_inv_g     <- chol2inv(ch_Fg)
          log_det_F_g <- 2 * sum(log(diag(ch_Fg)))
        }

        log_w[i] <- -0.5 * n_ok * log(2 * pi) -
                    0.5 * log_det_F_g -
                    0.5 * drop(v_ok %*% (F_inv_g %*% v_ok))
        shocks[, i] <- eps_i

      } else {
        ## -- Regime mismatch: KEEP the COPF draw and use the generally-valid
        ## importance weight  p(y_t | s_i, eps_i, r_verify) p(eps_i) / q(eps_i).
        ## Redrawing from the prior here would make the effective proposal a
        ## draw-dependent mixture whose density the weights ignore, biasing
        ## the likelihood estimator by O(mismatch rate) per period.
        n_fallback <- n_fallback + 1L

        if (!exists(as.character(r_verify), envir = regime_cache, inherits = FALSE))
          obc_ensure_policy(r_verify, regime_cache, sys, dr_slack, specs, obs_idx, copf_args)
        pol_v <- get(as.character(r_verify), envir = regime_cache, inherits = FALSE)

        shocks[, i] <- eps_i

        ## log p(y_t | s_i, eps_i, r_verify): measurement density under the
        ## verified regime's policy
        obs_offset_v <- d_obs + pol_v$c_obs
        y_pred_v <- drop(pol_v$ZZ %*% s_i) + drop(pol_v$DD %*% eps_i) + obs_offset_v
        innov_v  <- y_t - y_pred_v
        obs_ok_v <- which(!is.na(innov_v))
        n_ok_v   <- length(obs_ok_v)
        log_p_y  <- if (n_ok_v == 0L) 0 else {
          innov_ok_v <- innov_v[obs_ok_v]
          -0.5 * n_ok_v * log(2 * pi) - 0.5 * n_ok_v * log(me_variance) -
            0.5 * sum(innov_ok_v^2) / me_variance
        }

        ## log p(eps_i) under the prior N(0, Sigma_e)
        u_p <- forwardsolve(L_e, eps_i)
        log_p_eps <- -sum(log(diag(L_e))) - 0.5 * sum(u_p^2)

        ## log q(eps_i) under the proposal N(mu_g, Omega_g) actually drawn from
        u_q <- forwardsolve(L_Omega_g, eps_i - mu_g)
        log_q_eps <- -sum(log(diag(L_Omega_g))) - 0.5 * sum(u_q^2)

        ## (the -0.5 * n_exo * log(2*pi) normalisers cancel in p/q)
        log_w[i] <- log_p_y + log_p_eps - log_q_eps
      }
    }
  }

  ## -- Accumulate log_lik_contrib BEFORE resampling -------------------------
  log_lik_contrib <- .smc_log_sum_exp(log_w) - log(N)

  ## -- Normalize and systematic resample ------------------------------------
  log_w_c <- log_w - max(log_w)
  w_norm  <- exp(log_w_c)
  w_norm  <- w_norm / sum(w_norm)

  idx    <- .smc_systematic_resample(w_norm, N)
  shocks <- shocks[, idx, drop = FALSE]

  ## -- Propagate resampled particles ----------------------------------------
  ## Also track each particle's verified regime so the next period can use it
  ## as its proposal guess (ancestor-regime variance reduction).
  new_particles    <- matrix(0, n_state, N)
  verified_regimes <- integer(N)
  for (i in seq_len(N)) {
    orig_i <- idx[i]
    s_i    <- particles[, orig_i]
    eps_i  <- shocks[, i]

    ## Re-resolve regime after resampling
    bind_flags <- pkf_check_binding(s_i, eps_i, specs, dr_slack)
    regime_idx <- obc_regime_idx(bind_flags)
    if (!exists(as.character(regime_idx), envir = regime_cache, inherits = FALSE))
      obc_ensure_policy(regime_idx, regime_cache, sys, dr_slack, specs, obs_idx, copf_args)
    pol <- get(as.character(regime_idx), envir = regime_cache, inherits = FALSE)

    new_particles[, i]    <- drop(pol$TT %*% s_i) + drop(pol$RR %*% eps_i) + pol$c_state
    verified_regimes[i]   <- regime_idx
  }

  list(particles = new_particles, log_lik_contrib = log_lik_contrib,
       n_fallback = n_fallback, verified_regimes = verified_regimes,
       z_copf_used = z_copf, z_fallback_used = z_fallback)  # CPM: record used draws
}


## ============================================================================
## Main PPF likelihood function
## ============================================================================

#' Bootstrap Piecewise Particle Filter log-likelihood for OBC models
#'
#' Evaluates the marginal log-likelihood log p(Y | theta) using a bootstrap
#' particle filter that treats the OBC regime sequence r_{1:T} as a latent
#' variable.
#'
#' At each period t, N particles carry s_{t-1}^i; each particle:
#'   1. Draws eps_t^i ~ N(0, Sigma_e) (prior proposal)
#'   2. Resolves its OBC regime from (s_{t-1}^i, eps_t^i) via the SHARED
#'      regime_cache (one cache per likelihood evaluation)
#'   3. Computes weight w_t^i = p(y_t | s_{t-1}^i, eps_t^i, regime_t^i)
#'   4. Accumulates log_lik += log(mean w_t^i) BEFORE resampling
#'   5. Systematically resamples and propagates the state
#'
#' @param Y             n_obs x T observation matrix
#' @param dr_slack      Slack-regime DecisionRules
#' @param regime_cache  FRESH R environment for this call (modified in-place)
#' @param sys           System matrices (for lazy regime building)
#' @param model         dynhr_mod
#' @param params        Named numeric parameter vector
#' @param obs_vars      Character vector of observed variable names
#' @param specs         OBC spec list from obc_parse_tags
#' @param obs_idx       Integer vector of observable indices
#' @param N             Number of particles (default 1000)
#' @param me_variance   Measurement error variance (must be > 0; default 1e-4)
#' @param return_particles Logical (default FALSE).  When TRUE, append the
#'   terminal \code{n_state x N} particle cloud as \code{$particles} in the
#'   return list.  Has no effect on the loglik path (byte-identical when FALSE).
#' @param seed          Integer RNG seed for reproducibility (NULL = no seed)
#' @return List with:
#'   $loglik       scalar log-likelihood
#'   $n_obs        integer
#'   $n_T          integer
#'   $particles    n_state x N terminal particle cloud (only when return_particles = TRUE)
#' @noRd
ppf_likelihood <- function(Y, dr_slack, regime_cache, sys,
                            model, params, obs_vars, specs,
                            obs_idx          = NULL,
                            N                = 1000L,
                            me_variance      = 1e-4,
                            proposal         = c("bootstrap", "copf"),
                            regime_guess     = c("ancestor", "slack"),
                            U_copf_list      = NULL,  # CPM: list of per-period U_copf structures
                            return_particles  = FALSE,
                            seed             = NULL) {

  proposal     <- match.arg(proposal)
  regime_guess <- match.arg(regime_guess)
  if (!is.null(seed)) {
    ## Local seed: deterministic in theta but leaves the CALLER's RNG stream
    ## untouched (an outer sampler must not replay the same proposals; see
    ## .with_local_seed in tpf-likelihood.R and NEWS 0.9.2.0003).
    .ge_ <- globalenv()
    .had_ <- exists(".Random.seed", envir = .ge_, inherits = FALSE)
    .old_ <- if (.had_) get(".Random.seed", envir = .ge_, inherits = FALSE) else NULL
    on.exit({
      if (.had_) assign(".Random.seed", .old_, envir = .ge_)
      else if (exists(".Random.seed", envir = .ge_, inherits = FALSE))
        rm(list = ".Random.seed", envir = .ge_)
    }, add = TRUE)
    set.seed(seed)
  }

  endo    <- dr_slack$endo_names
  exo     <- dr_slack$exo_names
  n_state <- length(dr_slack$state_idx)
  n_exo   <- length(exo)
  n_obs   <- length(obs_vars)

  if (is.null(obs_idx)) obs_idx <- match(obs_vars, endo)

  ## Data orientation: n_obs x T
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## Shock covariance + Cholesky factor
  Sigma_e <- .get_shock_cov(model, exo, params)
  L_e     <- t(chol(Sigma_e))   # lower-triangular L s.t. L %*% t(L) = Sigma_e

  ## COPF pre-computation: Sigma_e_inv + copf_args for cache population
  copf_args   <- NULL
  Sigma_e_inv <- NULL
  if (proposal == "copf") {
    ch_Se       <- tryCatch(chol(Sigma_e), error = function(e2) NULL)
    if (is.null(ch_Se))
      return(list(loglik = -Inf, n_obs = n_obs, n_T = n_T))
    Sigma_e_inv <- chol2inv(ch_Se)
    copf_args   <- list(Sigma_e = Sigma_e, Sigma_e_inv = Sigma_e_inv,
                        me_variance = me_variance)
    ## Ensure the already-built regime-0 entry gets COPF fields added
    ## (it was inserted before proposal was known; update it now)
    if (exists("0", envir = regime_cache, inherits = FALSE)) {
      pol0 <- get("0", envir = regime_cache, inherits = FALSE)
      if (is.null(pol0$Omega)) {
        ## Remove and re-build with COPF args
        rm(list = "0", envir = regime_cache)
        obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx, copf_args)
      }
    }
  }

  ## Observable steady-state means
  d_obs <- dr_slack$ys[obs_vars]

  ## ---- Initialise particles from slack-policy stationary distribution -----
  pol_s <- get("0", envir = regime_cache, inherits = FALSE)
  QQ_s  <- tcrossprod(pol_s$RR %*% Sigma_e, pol_s$RR)
  P_0   <- tryCatch(
    solve_lyapunov(pol_s$TT, QQ_s),
    error = function(e2) diag(1e-6, n_state)
  )
  if (any(!is.finite(P_0))) P_0 <- diag(1e-6, n_state)

  L_P <- tryCatch(t(chol(P_0)), error = function(e2) diag(1e-3, n_state))
  particles <- L_P %*% matrix(rnorm(n_state * N), nrow = n_state)

  ## ---- Main filter loop ---------------------------------------------------
  loglik           <- 0
  total_fallback   <- 0L
  prev_regimes     <- NULL   # NULL => all-slack guess at t=1 (COPF only)
  U_copf_realized  <- vector("list", n_T)  # CPM: collect used z_copf/z_fallback per period

  for (t in seq_len(n_T)) {
    y_t <- Y[, t]

    if (all(is.na(y_t))) next

    if (proposal == "bootstrap") {
      res <- .ppf_run_period(
        particles    = particles,
        y_t          = y_t,
        L_e          = L_e,
        Sigma_e      = Sigma_e,
        dr_slack     = dr_slack,
        regime_cache = regime_cache,
        sys          = sys,
        specs        = specs,
        obs_idx      = obs_idx,
        d_obs        = d_obs,
        me_variance  = me_variance
      )
    } else {
      ## Pass prev_regimes only when using ancestor guessing.
      ## regime_guess="slack" always passes NULL so r_guess=0 for all particles.
      ## CPM: thread per-period U_copf if supplied.
      U_copf_t <- if (!is.null(U_copf_list)) U_copf_list[[t]] else NULL
      res <- .copf_run_period(
        particles    = particles,
        y_t          = y_t,
        L_e          = L_e,
        Sigma_e      = Sigma_e,
        Sigma_e_inv  = Sigma_e_inv,
        dr_slack     = dr_slack,
        regime_cache = regime_cache,
        sys          = sys,
        specs        = specs,
        obs_idx      = obs_idx,
        d_obs        = d_obs,
        me_variance  = me_variance,
        copf_args    = copf_args,
        prev_regimes = if (regime_guess == "ancestor") prev_regimes else NULL,
        U_copf       = U_copf_t
      )
      total_fallback <- total_fallback + res$n_fallback
      ## Store verified regimes for the next period's guess
      prev_regimes <- res$verified_regimes
      ## Record used z_copf/z_fallback matrices (CPM: deterministic U structure)
      U_copf_realized[[t]] <- list(z_copf = res$z_copf_used, z_fallback = res$z_fallback_used)
    }

    loglik    <- loglik + res$log_lik_contrib
    particles <- res$particles

    if (!is.finite(loglik))
      return(list(loglik = -Inf, n_obs = n_obs, n_T = n_T))
  }

  out <- list(loglik = loglik, n_obs = n_obs, n_T = n_T)
  if (proposal == "copf") {
    out$n_fallback        <- total_fallback
    out$U_copf_realized   <- U_copf_realized  # CPM: used z_copf/z_fallback per period
  }
  ## return_particles = TRUE: append terminal n_state x N cloud (for ctx terminal-state
  ## dispatch; default FALSE preserves byte-identical loglik-only path).
  if (isTRUE(return_particles)) out$particles <- particles
  out
}


## ============================================================================
## Log-posterior factory
## ============================================================================

#' Create a PPF-based log-posterior evaluator for OBC models
#'
#' Bootstrap particle filter variant of \code{make_log_posterior_obc_pkf}.
#' Each evaluation builds a FRESH regime_cache (theta-dependent) and runs
#' the bootstrap PPF with N particles.
#'
#' me_variance > 0 is required (hard stop): the bootstrap weights degenerate
#' when me_variance = 0 with n_obs < n_exo.
#'
#' @param model       dynhr_mod
#' @param data        Observation matrix (T x n_obs or n_obs x T)
#' @param prior_spec  Prior specification data.frame
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled
#' @param specs       OBC spec list (from obc_parse_tags), or NULL to parse
#' @param me_variance Measurement error variance (must be > 0; default 1e-4)
#' @param N           Number of particles per evaluation (default 1000)
#' @param proposal    "bootstrap" (prior proposal) or "copf" (conditionally
#'   optimal proposal with regime verification and general-ratio fallback)
#' @param regime_guess For proposal = "copf": "ancestor" (default) guesses
#'   each particle's regime from its resampled ancestor's verified regime
#'   (all-slack at t = 1); "slack" always guesses all-slack. Pure variance
#'   reduction — weights are valid for any guess.
#' @param seed        Integer RNG seed (NULL = not fixed; each call differs)
#' @param power       Power-posterior (generalised-Bayes) tempering exponent
#'   \eqn{\zeta}: \code{$logpost} becomes
#'   \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)} while \code{$loglik}
#'   keeps the RAW (untempered) particle-filter estimate -- so PMMH's
#'   unbiasedness argument and any marginal-likelihood use of \code{$loglik}
#'   are unaffected. \code{NULL} (default) resolves the \code{power_posterior}
#'   package option ONCE, at factory time.
#' @return Function(theta) -> list(logpost, loglik, logprior)
#' @export
make_log_posterior_obc_ppf <- function(model, data, prior_spec, obs_vars,
                                        compiled, specs = NULL,
                                        me_variance  = 1e-4,
                                        N            = 1000L,
                                        proposal     = c("bootstrap", "copf"),
                                        regime_guess = c("ancestor", "slack"),
                                        seed         = NULL,
                                        power        = NULL) {

  ## Force promises (closure-capture safety; mirrors PKF and TPF factories)
  force(prior_spec); force(me_variance); force(N); force(seed)
  ## Resolve zeta ONCE here, not per draw (see .resolve_power_posterior).
  power <- .resolve_power_posterior(power, "make_log_posterior_obc_ppf")
  proposal     <- match.arg(proposal)
  regime_guess <- match.arg(regime_guess)

  ## Hard stop: me_variance = 0 degenerates weights for both bootstrap and COPF
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance <= 0) {
    stop(
      "make_log_posterior_obc_ppf: 'me_variance' must be a finite positive scalar.\n",
      "Bootstrap PF weights p(y_t | ...) degenerate at me_variance = 0 when ",
      "n_obs < n_exo. Recommended: me_variance >= 1e-6."
    )
  }

  if (is.null(specs)) specs <- obc_parse_tags(model)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence)) {
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  }
  endo    <- model$var_names
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("obs_vars contains names not found in model$var_names: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  ## Inner evaluator over the shared closure builder (R/posterior-closure.R):
  ## cold steady-state solve, always-eigen() stationarity guard, no system
  ## priors -- the OBC bootstrap/COPF particle filter is the only branch-
  ## specific part. `pass_dots = TRUE` returns the raw `function(theta, ...)`;
  ## the seeded wrapper below restores the public `function(theta)` signature.
  eval_one <- .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    stationarity = "eigen",
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      ## FRESH regime_cache per draw (theta-dependent matrices)
      regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
      obc_ensure_policy(0L, regime_cache, sol$sys, sol$dr, specs, obs_idx)

      pf <- ppf_likelihood(
        Y, sol$dr, regime_cache, sol$sys,
        model, params, obs_vars, specs,
        obs_idx      = obs_idx,
        N            = N,
        me_variance  = me_variance,
        proposal     = proposal,
        regime_guess = regime_guess,
        seed         = NULL   # the closure below already set the local seed
      )
      if (is.null(pf) || !is.finite(pf$loglik)) return(NULL)
      list(loglik = pf$loglik)
    },
    power          = power,
    warm_start     = FALSE,
    needs_me_floor = FALSE,
    pass_dots      = TRUE)

  ## ---- Closure: evaluated at each parameter draw -------------------------
  function(theta) {
    ## Fixed seed per theta for reproducibility (same seed = same loglik).
    ## LOCAL: deterministic in theta but leaves the CALLER's RNG stream
    ## untouched (an outer sampler must not replay the same proposals; see
    ## .with_local_seed in tpf-likelihood.R and NEWS 0.9.2.0003).
    .with_local_seed(seed, eval_one(theta))
  }
}
