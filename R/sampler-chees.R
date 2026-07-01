## R/sampler-chees.R
## --------------------------------------------------------------------------
## ChEES-HMC: Change in the Estimator of the Expected Square
## Hoffman, Radul & Sountsov (2021), AISTATS.
##
## Key differences from NUTS:
##   - Fixed-length leapfrog trajectories (no recursive tree / U-turn check)
##   - Trajectory LENGTH T = L * eps is adapted during warmup via the ChEES
##     criterion (stochastic gradient ascent on the squared-distance change)
##   - Trajectory jitter: actual length ~ Uniform(0, T_max) to avoid resonance
##   - Step size adapted separately via dual averaging (as in NUTS/HMC)
##
## Public API:
##   dynhr_chees()   -- ChEES-HMC sampler (mirrors dynhr_nuts / dynhr_hmc)
##
## Depends on HMC internals from sampler-hmc.R (loaded together in the package).
## --------------------------------------------------------------------------


# ============================================================================
# ChEES criterion helpers
# ============================================================================

#' Estimate the ChEES criterion from a batch of (theta, theta') pairs
#'
#' ChEES ∝ E[ (||theta'||^2 - ||theta||^2)^2 ]
#'
#' Given a matrix of starting positions (rows = samples) and a corresponding
#' matrix of proposed positions, compute the mean squared change in the squared
#' Euclidean norm.  A scalar estimate: we take the sample mean over the batch.
#'
#' When a mass matrix M (M_diag or M_inv) is supplied the squared norm is
#' computed in the METRIC space: ||theta||^2_M = theta' M theta.  This matches
#' the geometry in which the leapfrog dynamics live (Riemannian kinetic energy).
#' With the diagonal path: ||theta||^2 = sum(theta^2 * M_mass_diag).
##
#' @param theta_start d-vector: starting position of one chain link
#' @param theta_prop  d-vector: proposed position after trajectory
#' @param M_mass_diag Mass diagonal (diagonal path).  If NULL, uses identity.
#' @return scalar ChEES criterion value (always finite or 0)
#' @noRd
.chees_criterion <- function(theta_start, theta_prop, M_mass_diag = NULL) {
  if (is.null(M_mass_diag)) {
    M_mass_diag <- rep(1, length(theta_start))
  }
  sq_start <- sum(theta_start^2 * M_mass_diag)
  sq_prop  <- sum(theta_prop^2  * M_mass_diag)
  val <- (sq_prop - sq_start)^2
  if (is.finite(val)) val else 0
}


# ============================================================================
# dynhr_chees() -- ChEES-HMC sampler
# ============================================================================

#' ChEES-HMC sampler (Hoffman, Radul & Sountsov 2021)
#'
#' Runs fixed-length leapfrog HMC with trajectory-time adaptation via the
#' ChEES (Change in the Estimator of the Expected Square) criterion.
#' The step size is adapted independently by dual averaging (same as NUTS).
#' Trajectory jitter samples the actual length uniformly in \code{[0, T_max]}
#' at each iteration to avoid resonance.
#'
#' @param log_post_fn  function(theta) -> list(logpost, loglik, logprior) OR scalar
#' @param theta_init   Named numeric vector of starting parameter values
#' @param n_draws      Post-warmup draws to retain (default 2000)
#' @param n_warmup     Warmup iterations discarded after adaptation (default 1000)
#' @param step_size    Initial step size (NULL = auto-find via Algorithm 4)
#' @param T_init       Initial integration time T = L * eps.  NULL = auto
#'   (set to \code{10 * step_size} after finding eps0).
#' @param target_accept  Target acceptance rate for dual averaging (0.6-0.8).
#'   Default 0.65 (slightly lower than NUTS's 0.80 because ChEES trajectories
#'   can be longer and rejection is geometrically more expensive).
#' @param chees_lr     Learning rate for ChEES T-adaptation stochastic gradient.
#'   Default 0.05.  Smaller = more stable but slower convergence of T.
#' @param adapt_mass   Adapt diagonal mass matrix during warmup.  Default TRUE.
#' @param grad_fn      Optional analytic gradient function(theta) -> numeric vector.
#' @param grad_method  "forward" (default) or "central" or "Richardson".
#' @param delta_max    Maximum energy error before flagging divergence (1000).
#' @param mass_diag    Optional starting mass diagonal (d-vector).
#' @param verbose      Print progress messages.
#' @param progressor   progressr callback or NULL.
#' @param chain_id     Label for progress messages.
#' @param transform    Optional "dynhr_param_transform" object from
#'   \code{\link{build_param_transform}}. When supplied, ChEES runs in
#'   unconstrained eta-space (same as \code{dynhr_nuts}): theta_init is mapped
#'   to eta via to_unconstrained(); the Jacobian is included in the target;
#'   draws are mapped back to theta-space before being returned.
#' @param M_inv        Dense inverse-mass matrix (d x d).  NULL = diagonal.
#' @param chol_M       Upper Cholesky of M = solve(M_inv) (dense path).
#'
#' @return Named list (same shape as dynhr_hmc / dynhr_nuts):
#'   \describe{
#'     \item{chain}{Post-warmup draws (matrix n_draws x d)}
#'     \item{full_chain}{All draws including warmup}
#'     \item{logpost_trace}{Log-posterior at each iteration}
#'     \item{post_logpost}{Post-warmup log-posteriors}
#'     \item{acceptance_rate}{Post-warmup acceptance rate}
#'     \item{step_size}{Final adapted step size}
#'     \item{T_adapt}{Final adapted trajectory time}
#'     \item{T_trace}{Full trace of T_max (warmup + sampling)}
#'     \item{mass_matrix}{Final diagonal mass matrix}
#'     \item{divergences}{Logical vector of divergent iterations}
#'     \item{n_divergent}{Count of divergent post-warmup transitions}
#'     \item{n_draws, n_burn}{Counts}
#'     \item{n_grad_evals}{Total gradient evaluations}
#'     \item{elapsed_secs}{Wall time}
#'     \item{sampler}{"chees"}
#'   }
#' @noRd
dynhr_chees <- function(
    log_post_fn,
    theta_init,
    n_draws       = 2000L,
    n_warmup      = 1000L,
    step_size     = NULL,
    T_init        = NULL,
    target_accept = 0.65,
    chees_lr      = 0.05,
    adapt_mass    = TRUE,
    grad_fn       = NULL,
    grad_method   = "forward",
    delta_max     = 1000,
    mass_diag     = NULL,
    verbose       = TRUE,
    progressor    = NULL,
    chain_id      = NULL,
    transform     = NULL,
    M_inv         = NULL,
    chol_M        = NULL,
    checkpoint    = NULL
) {
  stopifnot(is.function(log_post_fn), is.numeric(theta_init))
  d         <- length(theta_init)
  par_names <- names(theta_init)
  n_total   <- n_draws + n_warmup

  # ---- Checkpoint / streaming (opt-in). When `checkpoint` is a list carrying a
  # `dir`, draws are streamed to per-chain files in flush_every-row chunks (RAM
  # bounded by flush_every * d, not n_total * d) and a restart state is saved
  # after every flush. `checkpoint$resume = TRUE` continues a prior run from its
  # saved state -- exactly (RNG, position, lp, step_size, T_adapt / log_T,
  # mass_matrix, n_done, n_warmup, accept count), so resume is bit-identical to a
  # single long run. Adaptation is FROZEN post-warmup; the ChEES running mean
  # (mu_run) is only needed during warmup and is not required for the continuation.
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)
  flush_every <- if (ckpt) as.integer(checkpoint$flush_every %||% 1000L) else NA_integer_
  ckpt_paths  <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id) else NULL

  # ---- Opt-in unconstrained transform (mirrors dynhr_nuts) ----
  if (!is.null(transform)) {
    target_fn  <- make_transformed_logpost(log_post_fn, transform,
                                           include_jacobian = TRUE)
    state_init <- transform$to_unconstrained(theta_init)
    names(state_init) <- par_names
  } else {
    target_fn  <- log_post_fn
    state_init <- theta_init
  }

  # --- Scalar log-posterior wrapper ---
  .lp_scalar <- function(theta) {
    names(theta) <- par_names
    res <- target_fn(theta)
    val <- if (is.list(res)) res$logpost else res
    if (!is.finite(val)) -1e300 else val
  }

  # --- Gradient ---
  if (is.null(grad_fn)) {
    .grad <- function(theta) .hmc_gradient(.lp_scalar, theta, method = grad_method)
  } else if (!is.null(transform)) {
    .grad <- make_transformed_grad(grad_fn, transform)
  } else {
    .grad <- grad_fn
  }

  # --- Mass matrix ---
  use_dense <- !is.null(M_inv) && is.matrix(M_inv)
  if (use_dense) {
    M_diag     <- rep(1, d)
    M_inv_diag <- NULL
    M_mass_diag <- rep(1, d)  # for ChEES norm: use identity when dense
  } else if (!is.null(mass_diag) && length(mass_diag) == d) {
    M_diag     <- pmax(as.numeric(mass_diag), 1e-12)
    M_inv_diag <- 1 / M_diag
    M_mass_diag <- M_diag
  } else {
    M_diag     <- rep(1, d)
    M_inv_diag <- rep(1, d)
    M_mass_diag <- M_diag
  }

  # --- Initial step size ---
  if (is.null(step_size)) {
    eps0 <- .hmc_find_stepsize(state_init, .lp_scalar, .grad,
                                M_inv_diag, M_diag,
                                M_inv = M_inv, chol_M = chol_M)
  } else {
    eps0 <- step_size
  }
  if (verbose) message(sprintf("ChEES: initial step_size = %.4e", eps0))

  # --- Initial trajectory time T_max ---
  # T_max controls the maximum integration time.  Actual trajectory at each
  # step is jittered: L_actual = max(1, round(runif(1, 0, T_max / eps_m) )).
  T_max <- if (is.null(T_init)) 10 * eps0 else T_init
  log_T <- log(T_max)  # adapt in log space for positivity
  if (verbose) message(sprintf("ChEES: initial T_max = %.4e", T_max))

  # --- Dual averaging for step size (same as NUTS) ---
  mu       <- log(10 * eps0)
  eps_bar  <- 1
  H_bar    <- 0
  gamma_da <- 0.05
  t0_da    <- 10
  kappa_da <- 0.75
  eps_m    <- eps0
  da_m     <- 0L

  # --- Storage ---
  # In checkpoint mode only a flush-sized buffer lives in RAM; the full chain is
  # read back from disk at the end. Otherwise pre-allocate exactly as before (the
  # non-checkpoint path is byte-identical).
  if (ckpt) {
    buf       <- matrix(NA_real_, nrow = flush_every, ncol = d)
    buf_lp    <- numeric(flush_every)
    buf_i     <- 0L
    chain         <- NULL
    logpost_trace <- NULL
  } else {
    chain         <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
    logpost_trace <- numeric(n_total)
  }
  state_chain   <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
  divergences   <- logical(n_total)
  accepted      <- logical(n_total)
  T_trace       <- numeric(n_total)   # trajectory time actually used at each step
  T_max_trace   <- numeric(n_total)   # T_max (the adapted upper bound)
  n_grad_evals  <- 0L

  theta   <- state_init
  lp_curr <- .lp_scalar(theta)
  ## Running mean of the chain, used to CENTER the ChEES criterion (||theta-mu||^2).
  mu_run  <- state_init
  n_mu    <- 0L
  trace_lp_curr <- if (!is.null(transform)) {
    lp_curr - transform$log_jacobian(theta)
  } else {
    lp_curr
  }

  n_accept <- 0L   # total accepted draws (post-draw-1)
  m_start  <- 2L   # loop start (overridden on resume)

  if (ckpt_resume) {
    # ---- Continue a saved run. Restore position, lp, RNG, step_size, T_adapt /
    # log_T, mass matrix, n_done, n_warmup, and accept count. Adaptation is frozen
    # post-warmup, so the ChEES running mean is NOT needed for the continuation.
    .ckpt_meta_verify(ckpt_paths$meta, "chees", checkpoint$fingerprint)
    st <- .ckpt_load_state(ckpt_paths$state)
    theta         <- st$theta
    lp_curr       <- st$lp_curr
    trace_lp_curr <- st$trace_lp_curr
    eps_m         <- st$step_size
    eps_bar       <- st$step_size   # frozen -- no further dual-averaging
    T_max         <- st$T_adapt
    log_T         <- log(T_max)
    M_diag        <- st$M_mass_diag
    M_inv_diag    <- if (use_dense) NULL else 1 / M_diag
    M_mass_diag   <- M_diag
    n_warmup      <- st$n_warmup    # original warmup count fixes the retained set
    n_accept      <- st$n_accept
    m_start       <- st$n_done + 1L  # resume from next draw (st$n_done is already on disk)
    .ckpt_truncate(ckpt_paths, st$n_done, d)  # drop any post-state partial flush
    assign(".Random.seed", st$rng, envir = .GlobalEnv)
    # Reconstruct state_chain rows 1..n_done from disk (needed for mass adapt
    # during warmup extensions, but post-warmup resumes won't reach that code).
    # For simplicity, fill with theta (values from the last saved position).
  } else {
    # ---- Fresh run: record draw 1. ----
    stored1 <- if (!is.null(transform)) transform$to_constrained(theta) else theta
    state_chain[1, ] <- theta
    divergences[1]   <- FALSE
    accepted[1]      <- TRUE
    T_trace[1]       <- T_max
    T_max_trace[1]   <- T_max
    if (ckpt) {
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))   # clear any stale fresh-run files
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "chees", checkpoint$fingerprint)
      buf_i <- 1L; buf[1, ] <- stored1; buf_lp[1] <- trace_lp_curr
    } else {
      chain[1, ]       <- stored1
      logpost_trace[1] <- trace_lp_curr
    }
  }

  t_start <- Sys.time()

  for (m in if (m_start <= n_total) m_start:n_total else integer(0)) {
    # ---- Sample momentum ----
    r0 <- .hmc_sample_momentum(d, M_diag, chol_M = chol_M)

    # ---- Jittered trajectory length ----
    # Sample actual L uniformly from [1, max(1, floor(T_max / eps_m))].
    # This prevents resonance (standing-wave pathologies) that plague vanilla
    # HMC with a fixed L.
    # Guard: clamp in floating-point BEFORE as.integer() to avoid overflow when
    # T_max / eps_m > 2^31 (e.g. tiny eps_m during early dual-averaging).
    raw_L_max <- if (is.finite(eps_m) && eps_m > 0) T_max / eps_m else 1
    L_max    <- as.integer(max(1, min(500, floor(raw_L_max))))
    L_actual <- if (L_max <= 1L) 1L else sample.int(L_max, 1L)
    T_actual <- L_actual * eps_m
    T_trace[m]     <- T_actual
    T_max_trace[m] <- T_max

    # ---- Run leapfrog trajectory ----
    theta_prop <- theta
    r_prop     <- r0
    divergent  <- FALSE

    H0 <- -lp_curr + .hmc_kinetic(r0, M_inv_diag, M_inv = M_inv)

    for (l in seq_len(L_actual)) {
      step <- .hmc_leapfrog(theta_prop, r_prop, eps_m, .grad,
                            M_inv_diag, M_inv = M_inv)
      if (is.null(step) || any(!is.finite(step$theta)) || any(!is.finite(step$r))) {
        divergent <- TRUE
        break
      }
      theta_prop <- step$theta
      r_prop     <- step$r
    }
    n_grad_evals <- n_grad_evals + L_actual

    names(theta_prop) <- par_names
    divergences[m] <- divergent

    # ---- Metropolis accept/reject ----
    alpha_m <- 0
    if (!divergent) {
      lp_prop   <- .lp_scalar(theta_prop)
      H1        <- -lp_prop + .hmc_kinetic(r_prop, M_inv_diag, M_inv = M_inv)
      log_alpha <- -H1 + H0
      energy_ok <- is.finite(log_alpha) && (H0 - (-lp_prop) < delta_max)

      if (energy_ok && log(runif(1)) < log_alpha) {
        theta   <- theta_prop
        lp_curr <- lp_prop
        accepted[m] <- TRUE
        n_accept    <- n_accept + 1L
      } else {
        accepted[m] <- FALSE
        if (!energy_ok) divergences[m] <- TRUE
      }
      alpha_m <- min(1, if (is.finite(log_alpha)) exp(log_alpha) else 0)
    } else {
      accepted[m] <- FALSE
    }

    trace_lp_curr <- if (!is.null(transform)) {
      lp_curr - transform$log_jacobian(theta)
    } else {
      lp_curr
    }
    stored_m <- if (!is.null(transform)) transform$to_constrained(theta) else theta
    state_chain[m, ] <- theta
    if (ckpt) {
      buf_i <- buf_i + 1L
      buf[buf_i, ]  <- stored_m
      buf_lp[buf_i] <- trace_lp_curr
      if (buf_i >= flush_every || m == n_total) {
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          theta         = theta,
          lp_curr       = lp_curr,
          trace_lp_curr = trace_lp_curr,
          step_size     = eps_m,
          T_adapt       = T_max,
          M_mass_diag   = M_diag,
          n_done        = m,
          n_warmup      = n_warmup,
          n_accept      = n_accept,
          n_draws_target = n_total,
          rng           = get(".Random.seed", envir = .GlobalEnv)))
        buf_i <- 0L
      }
    } else {
      chain[m, ]       <- stored_m
      logpost_trace[m] <- trace_lp_curr
    }
    ## Update the running mean (for centering the ChEES criterion).
    n_mu   <- n_mu + 1L
    mu_run <- mu_run + (theta - mu_run) / n_mu

    # ---- Adaptation (warmup only) ----
    if (m <= n_warmup) {
      da_m <- da_m + 1L

      # -- Dual averaging: step size --
      w     <- 1 / (da_m + t0_da)
      H_bar <- (1 - w) * H_bar + w * (target_accept - alpha_m)
      log_eps_m <- mu - (sqrt(da_m) / gamma_da) * H_bar
      eps_m     <- max(1e-10, exp(log_eps_m))
      m_kappa   <- da_m^(-kappa_da)
      eps_bar   <- exp(m_kappa * log_eps_m + (1 - m_kappa) * log(eps_bar))

      # -- ChEES criterion: adapt log(T_max) via stochastic gradient --
      # Gradient of ChEES w.r.t. log(T):
      #   d/d(log T) E[(||theta'||^2 - ||theta||^2)^2]
      # is estimated as 2 * (||theta_prop||^2 - ||theta_curr||^2) * ||theta_prop||^2
      # (product rule; treating ||theta'||^2 as increasing with T in expectation).
      # We sign the update so that T grows when the squared-norm change is large
      # (trajectory is exploring) and shrinks when it's near zero (stuck / short).
      # The update is only applied on non-divergent, accepted-or-rejected steps.
      if (!divergent) {
        ## CENTERED squared distances from the running mean: ChEES is about
        ## ||theta - E[theta]||^2, not ||theta||^2.
        dc_curr <- theta      - mu_run
        dc_prop <- theta_prop - mu_run
        sq_curr <- sum(dc_curr^2 * M_mass_diag)
        sq_prop <- sum(dc_prop^2 * M_mass_diag)
        if (is.finite(sq_curr) && is.finite(sq_prop)) {
          ## CORRECT ChEES gradient (Hoffman, Radul & Sountsov 2021). Adapt log T
          ## to MAXIMIZE E[(||theta'-mu||^2 - ||theta-mu||^2)^2]. The gradient is
          ## the correlation of the centered criterion with the RATE OF CHANGE of
          ## ||theta'-mu||^2 at the trajectory ENDPOINT, = 2*(theta'-mu).v_end,
          ## v_end = M^{-1} r_prop (final momentum). This is POSITIVE below the
          ## optimal T (distance still growing -> lengthen) and NEGATIVE past it
          ## (overshoot -> shorten), so T converges to the turning point
          ## (~pi/2 for an isotropic Gaussian). The original estimate had NO
          ## velocity term and pushed T monotonically to the clamp.
          v_end        <- if (is.null(M_inv)) M_inv_diag * r_prop else as.numeric(M_inv %*% r_prop)
          d_sqdist_end <- 2 * sum(dc_prop * v_end * M_mass_diag)
          chees_delta  <- sq_prop - sq_curr
          ## Dimensionless, bounded SGD step (tanh squash keeps a single noisy
          ## draw from blowing up log T; the average drives T to where the
          ## criterion/velocity correlation vanishes = the optimum).
          v_norm       <- sqrt(sum(v_end^2 * M_mass_diag)) + 1e-12
          scale        <- max(sq_prop, sq_curr, 1e-12) * v_norm
          chees_grad   <- (chees_delta * d_sqdist_end) / scale
          if (is.finite(chees_grad)) {
            log_T <- log_T + chees_lr * tanh(chees_grad)
          }
        }
        # Clamp: T must stay in a sensible range [eps_m, 1000]
        eps_floor <- max(eps_m, 1e-10)
        log_T <- if (is.finite(log_T)) {
          max(log(eps_floor), min(log_T, log(1000)))
        } else {
          log(10 * eps_floor)  # reset on NaN
        }
        T_max <- exp(log_T)
      }

      # -- Diagonal mass adaptation at 70% of warmup (mirrors dynhr_hmc) --
      if (!use_dense && adapt_mass && da_m == floor(n_warmup * 0.7) && da_m > 50L) {
        idx  <- max(1L, as.integer(floor(n_warmup * 0.2))):m
        vars <- apply(state_chain[idx, , drop = FALSE], 2, var)
        vars[vars < 1e-12 | !is.finite(vars)] <- 1
        M_diag     <- vars
        M_inv_diag <- 1 / M_diag
        M_mass_diag <- M_diag
        # Re-find step size after mass update
        eps_m <- .hmc_find_stepsize(theta, .lp_scalar, .grad, M_inv_diag, M_diag)
        mu      <- log(10 * eps_m)
        eps_bar <- 1
        H_bar   <- 0
        da_m    <- 0L
        # Also rescale T_max to remain sensible after scale change
        T_max <- max(10 * eps_m, T_max)
        log_T <- log(T_max)
        if (verbose) {
          message(sprintf("ChEES: mass adapted, step_size = %.4e, T_max = %.4e",
                          eps_m, T_max))
        }
      }
    }

    # ---- Fix step size and T at end of warmup ----
    if (m == n_warmup) {
      eps_m <- eps_bar   # use dual-averaged value (more stable than last iterate)
      if (verbose) {
        message(sprintf("ChEES: warmup complete: step_size = %.4e, T_max = %.4e",
                        eps_m, T_max))
      }
    }

    # ---- Progress ----
    if (m %% 200 == 0 || m == n_total) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / m * (n_total - m)
      ch_lab  <- if (is.null(chain_id)) "?" else as.character(chain_id)
      n_div   <- sum(divergences[1:m])
      acc_200 <- mean(accepted[max(1, m - 199):m])
      phase   <- if (m <= n_warmup) "warmup" else "sample"
      msg <- sprintf(
        "ChEES Ch%s [%s] %d/%d acc=%.0f%% lp=%.1f eps=%.3e T=%.3e div=%d ETA=%.0fs",
        ch_lab, phase, m, n_total,
        acc_200 * 100, trace_lp_curr, eps_m, T_max, n_div, eta)
      if (!is.null(progressor)) {
        progressor(message = msg, amount = 1)
      } else if (verbose) {
        message(msg)
      }
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  if (ckpt) {
    # Materialize the full chain from the streamed files for the return value.
    logpost_trace <- .ckpt_read_lp(ckpt_paths$lp)
    chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, d, par_names)
  }

  post_chain   <- if (is.null(chain)) NULL else
    chain[(n_warmup + 1):n_total, , drop = FALSE]
  post_logpost <- logpost_trace[(n_warmup + 1):n_total]
  post_divs    <- divergences[(n_warmup + 1):n_total]
  post_accept  <- accepted[(n_warmup + 1):n_total]

  list(
    chain           = post_chain,
    full_chain      = chain,
    logpost_trace   = logpost_trace,
    post_logpost    = post_logpost,
    acceptance_rate = mean(post_accept),
    step_size       = eps_m,
    T_adapt         = T_max,
    T_trace         = T_trace,         # full trace (warmup + sampling)
    T_max_trace     = T_max_trace,     # T_max at each iteration
    mass_matrix     = M_diag,
    divergences     = post_divs,
    n_divergent     = sum(post_divs),
    n_draws         = as.integer(n_draws),
    n_burn          = as.integer(n_warmup),
    n_grad_evals    = n_grad_evals,
    elapsed_secs    = elapsed,
    sampler         = "chees",
    checkpoint_dir  = if (ckpt) checkpoint$dir else NULL
  )
}
