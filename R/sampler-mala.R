## R/sampler-mala.R
## --------------------------------------------------------------------------
## Stage 1 (manifold-MCMC roadmap): Metropolis-Adjusted Langevin Algorithm
## preconditioned by a constant metric G (so G^{-1} ≈ posterior covariance).
##
## Stage 2 extension: position-dependent metric via metric_fn(theta).
##
## Proposal (preconditioned MALA):
##   mu(theta) = theta + (eps^2 / 2) * G_inv %*% grad(theta)
##   theta'    ~ N(mu(theta), eps^2 * G_inv)
##   i.e.  theta' = mu(theta) + eps * chol_Ginv' %*% rnorm(d)
##         where chol_Ginv = chol(G_inv)  (upper triangular)
##
## MH ratio (CONSTANT metric):
##   For constant G the Gaussian normalising constant is the same in both
##   directions (same covariance eps^2 G^{-1}), so only the quadratic forms
##   remain:
##     log q(a -> b) = -1/(2 eps^2) * (b - mu(a))' G (b - mu(a))
##     log_alpha     = (lp_prop - lp_curr) + (log_q_back - log_q_fwd)
##
## MH ratio (POSITION-DEPENDENT metric, simplified MMALA, Stage 2):
##   When metric_fn is supplied, G(theta) varies per step.  The proposal
##   covariances at theta and theta' DIFFER, so the Gaussian normalising
##   constants do NOT cancel.  The correct acceptance log-ratio is:
##
##     mu(theta)  = theta  + (eps^2/2) * G_inv(theta)  %*% grad(theta)
##     mu(theta') = theta' + (eps^2/2) * G_inv(theta') %*% grad(theta')
##
##     log q(theta -> theta') =
##       +0.5 * logdet(G(theta))
##       -0.5/eps^2 * (theta' - mu(theta))' G(theta) (theta' - mu(theta))
##       [+ 0.5*d*log(2*pi) + d*log(eps)  -- these CANCEL between directions]
##
##     log q(theta' -> theta) =
##       +0.5 * logdet(G(theta'))
##       -0.5/eps^2 * (theta  - mu(theta'))' G(theta') (theta  - mu(theta'))
##
##     log_alpha = (lp_prop - lp_curr)
##               + 0.5*(logdet(G(theta')) - logdet(G(theta)))   <- logdet terms
##               + (quadratic_back - quadratic_fwd)
##
##   The 0.5*logdet(G(.)) terms come from the normalising constant of
##   N(mu, eps^2 G^{-1}): log(2pi)^{d/2} |eps^2 G^{-1}|^{-1/2}
##   = const - 0.5 * log|eps^2 G^{-1}|
##   = const + 0.5 * logdet(G) - 0.5*d*log(eps^2).
##   The -0.5*d*log(eps^2) term CANCELS between directions; 0.5*logdet(G(.))
##   does NOT (different G at theta vs theta').
##   Omitting these terms silently biases the stationary distribution.
##
##   The dG/dtheta Christoffel drift terms are dropped (simplified MMALA):
##   dynhr has no third-order derivatives, and the MH correction ensures
##   correctness regardless (at the cost of some asymptotic efficiency).
##
## Step-size adaptation:
##   Dual averaging (Hoffman & Gelman 2014) targeting 0.574 acceptance
##   (the MALA optimum per Roberts & Rosenthal 1998 / Atchadé & Rosenthal
##   2003 for a d-dimensional target).  The dual-averaging code pattern is
##   identical to R/sampler-nuts.R.
##
## Interface mirrors dynhr_hmc / dynhr_nuts:
##   - same transform (eta-space) opt-in path
##   - same guard for non-finite grad / logpost (treat as auto-reject)
##   - same output list shape
## --------------------------------------------------------------------------


# ============================================================================
# Internal helpers
# ============================================================================

#' Log proposal density for one MALA step (constant metric, up to additive const)
#'
#' Computes log q(from -> to) = -1/(2 eps^2) * (to - mu_from)' G (to - mu_from),
#' omitting the normalising constant (same in both directions for constant G).
#'
#' @param from  Current position (numeric vector, length d)
#' @param to    Proposed position (numeric vector, length d)
#' @param mu_from  Proposal mean at `from`:  from + (eps^2/2) G_inv %*% grad(from)
#' @param G     Metric matrix (d x d SPD)
#' @param eps   Step size
#' @return Scalar log proposal density (up to additive constant)
#' @noRd
.mala_log_q <- function(from, to, mu_from, G, eps) {
  diff <- to - mu_from
  -0.5 / (eps^2) * sum(diff * as.numeric(G %*% diff))
}


#' Log proposal density for smMALA with POSITION-DEPENDENT metric
#'
#' Includes the normalising-constant term +0.5*logdet(G(from)) which does NOT
#' cancel when G varies between theta and theta'.  Use this when metric_fn is
#' position-dependent.
#'
#' log q(from -> to) = +0.5*logdet(G(from))
#'                     - (1/(2 eps^2)) * (to - mu_from)' G(from) (to - mu_from)
#'                     + [const in d, eps that cancels between forward/backward]
#'
#' @param to       Target position.
#' @param mu_from  Proposal mean at `from`.
#' @param G_from   Metric matrix at `from` (d x d SPD).
#' @param logdet_from  log|G(from)|.
#' @param eps      Step size.
#' @return Scalar including the 0.5*logdet(G(from)) normalising term.
#' @noRd
.mala_log_q_pd <- function(to, mu_from, G_from, logdet_from, eps) {
  diff <- to - mu_from
  quadratic <- sum(diff * as.numeric(G_from %*% diff))
  0.5 * logdet_from - 0.5 / (eps^2) * quadratic
}


#' Compute MALA proposal mean at position theta
#'
#' mu(theta) = theta + (eps^2 / 2) * G_inv %*% grad(theta)
#'
#' Returns NULL if the gradient is non-finite (off-support).
#'
#' @param theta  Current position
#' @param grad_fn  Gradient function (theta) -> numeric vector
#' @param G_inv  Inverse metric matrix (d x d SPD)
#' @param eps  Step size
#' @return Numeric vector (proposal mean) or NULL if non-finite gradient
#' @noRd
.mala_proposal_mean <- function(theta, grad_fn, G_inv, eps) {
  g <- grad_fn(theta)
  if (any(!is.finite(g))) return(NULL)
  theta + (eps^2 / 2) * as.numeric(G_inv %*% g)
}


# ============================================================================
# dynhr_mala() -- MALA sampler with dual-averaging step-size adaptation
# ============================================================================

#' MALA (Metropolis-Adjusted Langevin Algorithm) with constant or position-dependent metric
#'
#' Preconditioned MALA with either a fixed (Stage 1) or position-dependent
#' (Stage 2, simplified MMALA) SPD metric.  One gradient evaluation per
#' proposal step; step size adapted via dual averaging.
#'
#' @param log_post_fn  function(theta) -> list(logpost, loglik, logprior) or
#'   scalar.
#' @param theta_init   Named numeric starting vector.
#' @param n_draws      Post-warmup draws to retain.
#' @param n_warmup     Warmup iterations (adaptation and burn-in; discarded).
#' @param grad_fn      Optional analytic gradient function(theta) -> numeric
#'   vector.  If NULL, uses forward finite differences.
#' @param G            Metric matrix (\code{d x d} SPD; default identity).
#'   \code{G_inv} is derived as \code{solve(G)} if not supplied.
#'   Ignored when \code{metric_fn} is non-NULL.
#' @param G_inv        Inverse metric (\code{d x d} SPD; default identity).
#'   Ignored when \code{metric_fn} is non-NULL.
#' @param metric_fn    Optional \code{function(theta) -> list(G, G_inv, L, logdet)}.
#'   When supplied, the metric is evaluated at each proposed position (simplified
#'   MMALA, Stage 2).  The MH ratio includes the logdet normalising-constant
#'   correction (see header comment).  Fallback: when \code{metric_fn(theta)}
#'   fails (e.g. unit-root boundary), the constant \code{G}/\code{G_inv} is used
#'   instead.  Christoffel drift terms are dropped (no third derivatives).
#' @param eps          Initial step size (NULL = auto-find from a single
#'   grad-based heuristic).
#' @param adapt_step   Adapt step size via dual averaging during warmup
#'   (default TRUE).
#' @param target_accept  Target acceptance rate (default 0.574, the MALA
#'   optimum for a d-dimensional Gaussian).
#' @param transform    Optional \code{"dynhr_param_transform"} object from
#'   \code{build_param_transform}.  When non-NULL, MALA runs in unconstrained
#'   eta-space (same convention as \code{dynhr_nuts}).
#' @param verbose      Print progress messages.
#' @param progressor   progressr callback or NULL.
#' @param chain_id     Label for progress messages.
#' @param checkpoint   Optional list for streaming/restartable checkpointing.
#'   Fields: \code{dir} (directory path), \code{flush_every} (draws per flush,
#'   default 1000), \code{resume} (logical, default FALSE), \code{fingerprint}
#'   (config fingerprint for resume verification), \code{write_meta} (logical,
#'   default TRUE; set FALSE on parallel paths where the orchestrator writes
#'   meta.rds), \code{return_chain} (logical, default TRUE; set FALSE to keep
#'   draws on disk only).  When NULL (default), checkpointing is disabled and
#'   the non-checkpoint path is bit-identical to the pre-checkpoint code.
#'
#' @return Named list with fields compatible with \code{dynhr_hmc}:
#'   \code{chain}, \code{full_chain}, \code{logpost_trace},
#'   \code{post_logpost}, \code{acceptance_rate}, \code{step_size},
#'   \code{n_draws}, \code{n_burn}, \code{elapsed_secs}, \code{sampler},
#'   \code{checkpoint_dir}.
#' @noRd
dynhr_mala <- function(
    log_post_fn,
    theta_init,
    n_draws      = 2000L,
    n_warmup     = 1000L,
    grad_fn      = NULL,
    G            = NULL,
    G_inv        = NULL,
    metric_fn    = NULL,
    eps          = NULL,
    adapt_step   = TRUE,
    target_accept = 0.574,
    transform    = NULL,
    verbose      = TRUE,
    progressor   = NULL,
    chain_id     = NULL,
    checkpoint   = NULL
) {
  stopifnot(is.function(log_post_fn), is.numeric(theta_init))
  d         <- length(theta_init)
  par_names <- names(theta_init)
  n_total   <- n_draws + n_warmup

  # ---- Checkpoint / streaming (opt-in). When `checkpoint` is a list carrying
  # a `dir`, draws are streamed to per-chain files in flush_every-row chunks
  # (RAM bounded by flush_every * d, not n_total * d) and a restart state is
  # saved after every flush.  `checkpoint$resume = TRUE` continues a prior run
  # from its saved state -- exactly (RNG, position, lp, eps, n_done, n_accept,
  # n_warmup, met_curr for constant-metric runs), so a resume run is bit-
  # identical to the tail of a single long run.
  # The non-checkpoint path is unchanged (branch on `ckpt`).
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)
  flush_every <- if (ckpt) as.integer(checkpoint$flush_every %||% 1000L) else NA_integer_
  ckpt_paths  <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id) else NULL

  # ---- Metric setup -------------------------------------------------------
  # Position-dependent path (metric_fn): metric_fn(theta) -> list(G, G_inv, L, logdet).
  # Constant path: use supplied G/G_inv, default to identity.
  use_metric_fn <- !is.null(metric_fn) && is.function(metric_fn)

  if (is.null(G) && is.null(G_inv)) {
    G     <- diag(d)
    G_inv <- diag(d)
  } else if (is.null(G_inv)) {
    G_inv <- tryCatch(solve(G), error = function(e) {
      stop("dynhr_mala: supplied G is not invertible: ", conditionMessage(e), call. = FALSE)
    })
  } else if (is.null(G)) {
    G <- tryCatch(solve(G_inv), error = function(e) {
      stop("dynhr_mala: supplied G_inv is not invertible: ", conditionMessage(e), call. = FALSE)
    })
  }

  # Constant Cholesky (used when metric_fn is NULL, or as fallback)
  chol_Ginv <- tryCatch(chol(G_inv), error = function(e) {
    stop("dynhr_mala: G_inv is not positive definite (Cholesky failed): ",
         conditionMessage(e), call. = FALSE)
  })
  logdet_G_const <- sum(log(diag(chol(G))))   # log|G| for constant path

  # ---- Transform (eta-space) opt-in  ------------------------------------
  if (!is.null(transform)) {
    target_fn  <- make_transformed_logpost(log_post_fn, transform,
                                           include_jacobian = TRUE)
    state_init <- transform$to_unconstrained(theta_init)
    names(state_init) <- par_names
  } else {
    target_fn  <- log_post_fn
    state_init <- theta_init
  }

  # ---- Scalar log-posterior wrapper ------------------------------------
  .lp_scalar <- function(theta) {
    names(theta) <- par_names
    res <- target_fn(theta)
    val <- if (is.list(res)) res$logpost else res
    if (!is.finite(val)) -1e300 else val
  }

  # ---- Gradient function -----------------------------------------------
  if (is.null(grad_fn)) {
    .grad <- function(theta) .hmc_gradient(.lp_scalar, theta, method = "forward")
  } else if (!is.null(transform)) {
    .grad <- make_transformed_grad(grad_fn, transform)
  } else {
    .grad <- grad_fn
  }

  # ---- Helper: resolve metric at a position --------------------------------
  # Returns list(G, G_inv, chol_Ginv, logdet_G).
  # Falls back to the constant metric if metric_fn fails (unit-root etc.).
  .metric_at <- function(theta) {
    if (!use_metric_fn) {
      return(list(G = G, G_inv = G_inv, chol_Ginv = chol_Ginv,
                  logdet_G = logdet_G_const))
    }
    m <- tryCatch(metric_fn(theta), error = function(e) NULL)
    if (is.null(m) || !is.list(m) || is.null(m$G) || is.null(m$G_inv)) {
      # Fallback to constant metric
      return(list(G = G, G_inv = G_inv, chol_Ginv = chol_Ginv,
                  logdet_G = logdet_G_const))
    }
    # chol_Ginv: if L is chol(G) (upper triangular), then
    #   chol(G_inv) = solve(t(L))  (lower-triangular inverse of L transposed)
    #   But for sampling we need chol(G_inv) such that t(chol_Ginv)' t(chol_Ginv) = G_inv.
    # From L = chol(G) upper: G = t(L) L  -> G_inv = L^{-1} (L^{-1})^T
    # So chol(G_inv) = (t(L))^{-1} (upper triangular), i.e. solve(t(m$L)) transposed.
    # However, the sampling formula theta' = mu + eps * t(chol_Ginv) %*% z requires
    # Cov = eps^2 * t(chol_Ginv) %*% chol_Ginv = eps^2 * G_inv.
    # With chol_Ginv = chol(G_inv) (upper, so chol_Ginv^T chol_Ginv = G_inv) this works.
    # We compute chol(G_inv) directly for robustness.
    cGi <- tryCatch(chol(m$G_inv), error = function(e) NULL)
    if (is.null(cGi)) {
      return(list(G = G, G_inv = G_inv, chol_Ginv = chol_Ginv,
                  logdet_G = logdet_G_const))
    }
    logdet <- if (!is.null(m$logdet)) m$logdet else 2 * sum(log(diag(m$L)))
    list(G = m$G, G_inv = m$G_inv, chol_Ginv = cGi, logdet_G = logdet)
  }

  # ---- Initial step size -----------------------------------------------
  if (is.null(eps)) {
    eps <- .mala_find_stepsize(state_init, .lp_scalar, .grad, G_inv, G, chol_Ginv)
    if (verbose) message(sprintf("MALA: initial step_size = %.4e", eps))
  }

  # ---- Dual averaging parameters (identical to sampler-nuts.R:430) ------
  mu_da    <- log(10 * eps)
  eps_bar  <- 1
  H_bar    <- 0
  gamma_da <- 0.05
  t0_da    <- 10
  kappa_da <- 0.75
  da_m     <- 0L

  # ---- Storage -----------------------------------------------------------
  # In checkpoint mode only a flush-sized buffer lives in RAM; the full chain
  # is read back from disk at the end for the return value.  Otherwise,
  # pre-allocate the full chain (non-checkpoint path is bit-identical).
  if (ckpt) {
    buf    <- matrix(NA_real_, nrow = flush_every, ncol = d)
    buf_lp <- numeric(flush_every)
    buf_i  <- 0L
    chain         <- NULL
    logpost_trace <- NULL
  } else {
    chain         <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
    logpost_trace <- numeric(n_total)
  }
  accepted <- logical(n_total)

  theta   <- state_init
  lp_curr <- .lp_scalar(theta)
  trace_lp_curr <- if (!is.null(transform)) {
    lp_curr - transform$log_jacobian(theta)
  } else {
    lp_curr
  }
  n_accept  <- 0L
  t_start   <- Sys.time()
  m_start   <- 2L

  # Cache metric at current theta (avoid re-evaluation when rejected)
  met_curr <- .metric_at(theta)

  if (ckpt_resume) {
    # ---- Continue a saved run.  Refuse mismatched config, restore all
    # state needed for exact continuation: position, lp, RNG, draw count,
    # accept count, n_warmup, and the frozen post-warmup step size (eps).
    # For constant-metric runs the metric is already constructed from the
    # call args; for metric_fn runs it is recomputed each step so only eps
    # and the position matter.
    .ckpt_meta_verify(ckpt_paths$meta, "mala", checkpoint$fingerprint)
    st            <- .ckpt_load_state(ckpt_paths$state)
    theta         <- st$theta
    lp_curr       <- st$lp_curr
    trace_lp_curr <- st$trace_lp_curr
    eps           <- st$eps
    eps_bar       <- st$eps_bar
    H_bar         <- st$H_bar
    da_m          <- st$da_m
    n_accept      <- st$n_accept
    n_warmup      <- st$n_warmup    # original warmup fixes which rows are retained
    m_start       <- st$n_done + 1L
    met_curr      <- if (!is.null(st$met_curr)) st$met_curr else .metric_at(theta)
    .ckpt_truncate(ckpt_paths, st$n_done, d)
    assign(".Random.seed", st$rng, envir = .GlobalEnv)
  } else {
    # ---- Fresh run: record draw 1 (to the streaming buffer or in-RAM chain).
    stored1 <- if (!is.null(transform)) transform$to_constrained(theta) else theta
    if (ckpt) {
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))   # clear any stale files
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "mala", checkpoint$fingerprint)
      buf_i <- 1L; buf[1, ] <- stored1; buf_lp[1] <- trace_lp_curr
    } else {
      chain[1, ]       <- stored1
      logpost_trace[1] <- trace_lp_curr
    }
    accepted[1] <- TRUE
  }

  for (m in m_start:n_total) {
    # ---- Metric at current theta ----------------------------------------
    G_curr        <- met_curr$G
    G_inv_curr    <- met_curr$G_inv
    chol_Ginv_curr <- met_curr$chol_Ginv
    logdet_curr   <- met_curr$logdet_G

    # ---- Proposal mean at current theta --------------------------------
    mu_curr <- .mala_proposal_mean(theta, .grad, G_inv_curr, eps)

    if (is.null(mu_curr)) {
      # Off-support: non-finite gradient -- stay put (auto-reject)
      accepted[m] <- FALSE
      alpha_m     <- 0

    } else {
      # ---- Draw proposal  theta' = mu_curr + eps * t(chol_Ginv) %*% z --
      z          <- rnorm(d)
      theta_prop <- mu_curr + eps * as.numeric(t(chol_Ginv_curr) %*% z)
      names(theta_prop) <- par_names

      # ---- Log-posterior at proposal ------------------------------------
      lp_prop <- .lp_scalar(theta_prop)

      if (!is.finite(lp_prop)) {
        # Off-support proposal: auto-reject
        accepted[m] <- FALSE
        alpha_m     <- 0

      } else {
        # ---- Metric at proposal -----------------------------------------
        met_prop    <- .metric_at(theta_prop)
        G_prop      <- met_prop$G
        G_inv_prop  <- met_prop$G_inv
        logdet_prop <- met_prop$logdet_G

        # ---- Proposal mean at theta_prop (for reverse density) ---------
        mu_prop <- .mala_proposal_mean(theta_prop, .grad, G_inv_prop, eps)

        if (is.null(mu_prop)) {
          # Non-finite gradient at proposal: auto-reject
          accepted[m] <- FALSE
          alpha_m     <- 0

        } else {
          # ---- MH ratio ------------------------------------------------
          if (use_metric_fn) {
            # POSITION-DEPENDENT metric (smMALA, Stage 2):
            # log q(theta -> theta') = +0.5*logdet(G(theta))
            #                        - 1/(2 eps^2) (theta'-mu_curr)' G(theta) (theta'-mu_curr)
            # log q(theta' -> theta) = +0.5*logdet(G(theta'))
            #                        - 1/(2 eps^2) (theta-mu_prop)' G(theta') (theta-mu_prop)
            # log_alpha = (lp_prop - lp_curr)
            #           + (log_q_back - log_q_fwd)
            # The +0.5*logdet terms survive because G(theta) != G(theta').
            log_q_fwd  <- .mala_log_q_pd(theta_prop, mu_curr,  G_curr, logdet_curr, eps)
            log_q_back <- .mala_log_q_pd(theta,      mu_prop,  G_prop, logdet_prop, eps)
          } else {
            # CONSTANT metric (Stage 1): logdet cancels, use simple quadratic.
            log_q_fwd  <- .mala_log_q(theta, theta_prop, mu_curr, G_curr, eps)
            log_q_back <- .mala_log_q(theta_prop, theta, mu_prop, G_prop, eps)
          }

          log_alpha <- (lp_prop - lp_curr) + (log_q_back - log_q_fwd)

          # Clamp to guard -Inf / NaN
          if (!is.finite(log_alpha)) log_alpha <- -1e300

          alpha_m <- min(1, exp(log_alpha))

          if (log(runif(1)) < log_alpha) {
            theta         <- theta_prop
            lp_curr       <- lp_prop
            trace_lp_curr <- if (!is.null(transform)) {
              lp_curr - transform$log_jacobian(theta)
            } else {
              lp_curr
            }
            met_curr  <- met_prop   # cache for next iteration
            n_accept  <- n_accept + 1L
            accepted[m] <- TRUE
          } else {
            accepted[m] <- FALSE
            # met_curr stays the same (theta didn't change)
          }
        }
      }
    }

    # ---- Store draw m (streaming buffer or in-RAM chain) ---------------
    stored_m <- if (!is.null(transform)) transform$to_constrained(theta) else theta
    if (ckpt) {
      buf_i <- buf_i + 1L
      buf[buf_i, ]  <- stored_m
      buf_lp[buf_i] <- trace_lp_curr
      if (buf_i >= flush_every || m == n_total) {
        # Flush buffer to disk, then persist restart state (written atomically
        # AFTER draws so n_done never exceeds what is on disk).
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          theta         = theta,
          lp_curr       = lp_curr,
          trace_lp_curr = trace_lp_curr,
          eps           = eps,
          eps_bar       = eps_bar,
          H_bar         = H_bar,
          da_m          = da_m,
          n_accept      = n_accept,
          n_warmup      = n_warmup,
          n_done        = m,
          n_total_target = n_total,
          # Save the cached metric only for constant-metric runs; for
          # metric_fn runs it is recomputed from theta each step so we
          # save NULL and call .metric_at(theta) on resume.
          met_curr      = if (!use_metric_fn) met_curr else NULL,
          rng           = get(".Random.seed", envir = .GlobalEnv)
        ))
        buf_i <- 0L
      }
    } else {
      chain[m, ]       <- stored_m
      logpost_trace[m] <- trace_lp_curr
    }

    # ---- Dual averaging (step-size adaptation during warmup) -----------
    # Identical pattern to sampler-nuts.R:621-631.
    if (adapt_step && m <= n_warmup) {
      da_m <- da_m + 1L
      w    <- 1 / (da_m + t0_da)
      H_bar    <- (1 - w) * H_bar + w * (target_accept - alpha_m)
      log_eps  <- mu_da - (sqrt(da_m) / gamma_da) * H_bar
      eps      <- exp(log_eps)
      m_kappa  <- da_m^(-kappa_da)
      eps_bar  <- exp(m_kappa * log_eps + (1 - m_kappa) * log(eps_bar))
    }

    # Fix step size at end of warmup (use dual-averaged value)
    if (adapt_step && m == n_warmup) {
      eps <- eps_bar
      if (verbose) message(sprintf("MALA: warmup complete, final step_size = %.4e", eps))
    }

    # ---- Progress -------------------------------------------------------
    if (m %% 200 == 0 || m == n_total) {
      rate    <- n_accept / (m - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / m * (n_total - m)
      ch_lab  <- if (is.null(chain_id)) "?" else as.character(chain_id)
      phase   <- if (m <= n_warmup) "warmup" else "sample"
      msg     <- sprintf("MALA Ch%s [%s] %d/%d accept=%.0f%% lp=%.1f eps=%.3e ETA=%.0fs",
                         ch_lab, phase, m, n_total, rate * 100,
                         trace_lp_curr, eps, eta)
      if (!is.null(progressor)) {
        progressor(message = msg, amount = 1)
      } else if (verbose) {
        message(msg)
      }
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  if (ckpt) {
    # Materialise the full chain from the streamed files for the return value.
    # checkpoint$return_chain = FALSE skips this for very long runs (draws
    # remain on disk; read them with .ckpt_read_draws() or resume to extend).
    logpost_trace <- .ckpt_read_lp(ckpt_paths$lp)
    chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, d, par_names)
  }

  post_chain   <- if (is.null(chain)) NULL else
    chain[(n_warmup + 1):n_total, , drop = FALSE]
  post_logpost <- logpost_trace[(n_warmup + 1):n_total]

  list(
    chain           = post_chain,
    full_chain      = chain,
    logpost_trace   = logpost_trace,
    post_logpost    = post_logpost,
    acceptance_rate = n_accept / (n_total - 1),
    step_size       = eps,
    n_draws         = as.integer(n_draws),
    n_burn          = as.integer(n_warmup),
    elapsed_secs    = elapsed,
    sampler         = "mala",
    checkpoint_dir  = if (ckpt) checkpoint$dir else NULL
  )
}


# ============================================================================
# Step-size heuristic for MALA
# ============================================================================

#' Find a reasonable MALA initial step size
#'
#' Uses a single-step Langevin proposal from the starting position and doubles
#' or halves until the MH acceptance probability is approximately 50%.
#' Falls back to 0.01 if the gradient is not finite at theta_init.
#'
#' @noRd
.mala_find_stepsize <- function(theta, lp_fn, grad_fn, G_inv, G, chol_Ginv,
                                 target = 0.5, max_iter = 100L) {
  d   <- length(theta)
  eps <- 0.1  # initial trial

  # RNG-neutral: the trial proposals below draw rnorm(), but dynhr_mala()'s main
  # chain must continue from the SAME RNG position whether or not this search
  # ran -- otherwise the sampler is non-reproducible under a fixed seed (and not
  # checkpoint-streamable). Snapshot .Random.seed and restore it on exit. (This
  # replaces a `set.seed(NULL)` that reseeded from system entropy, which silently
  # made eps = NULL runs non-deterministic.)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    .rng0 <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", .rng0, envir = .GlobalEnv), add = TRUE)
  }

  # Compute gradient once (reused across iterations)
  mu0 <- .mala_proposal_mean(theta, grad_fn, G_inv, eps)
  if (is.null(mu0)) return(0.01)  # non-finite gradient at init

  lp0 <- lp_fn(theta)
  if (!is.finite(lp0)) return(0.01)

  # Single step with current eps
  .try_alpha <- function(eps_try) {
    g <- grad_fn(theta)
    if (any(!is.finite(g))) return(NA_real_)
    mu_try    <- theta + (eps_try^2 / 2) * as.numeric(G_inv %*% g)
    z         <- rnorm(d)
    prop      <- mu_try + eps_try * as.numeric(t(chol_Ginv) %*% z)
    lp_prop   <- lp_fn(prop)
    if (!is.finite(lp_prop)) return(0)
    mu_prop <- .mala_proposal_mean(prop, grad_fn, G_inv, eps_try)
    if (is.null(mu_prop)) return(0)
    log_q_fwd  <- .mala_log_q(theta, prop, mu_try, G, eps_try)
    log_q_back <- .mala_log_q(prop, theta, mu_prop, G, eps_try)
    log_alpha  <- (lp_prop - lp0) + (log_q_back - log_q_fwd)
    if (!is.finite(log_alpha)) return(0)
    min(1, exp(log_alpha))
  }

  alpha <- .try_alpha(eps)
  if (is.na(alpha)) return(0.01)

  # Determine direction: if alpha > target, grow eps; else shrink
  a <- if (alpha > target) 1 else -1

  for (k in seq_len(max_iter)) {
    eps_try <- eps * (2^a)
    alpha_try <- .try_alpha(eps_try)
    if (is.na(alpha_try)) break
    if (a * (alpha_try - target) <= 0) break  # crossed the target
    eps   <- eps_try
    alpha <- alpha_try
  }

  max(eps, 1e-10)
}


# ============================================================================
# .run_mala_batch() -- single-chain MALA, mirroring .run_hmc_batch()
# ============================================================================

#' Single-chain MALA batch, to be called from run_posterior_estimation
#'
#' @noRd
.run_mala_batch <- function(log_post_fn, theta_init,
                             n_draws = 2000L, n_warmup = 1000L,
                             G = NULL, G_inv = NULL,
                             metric_fn = NULL,
                             grad_fn = NULL,
                             transform = NULL,
                             verbose = TRUE, ...) {
  res <- dynhr_mala(
    log_post_fn = log_post_fn,
    theta_init  = theta_init,
    n_draws     = n_draws,
    n_warmup    = n_warmup,
    G           = G,
    G_inv       = G_inv,
    metric_fn   = metric_fn,
    grad_fn     = grad_fn,
    transform   = transform,
    verbose     = verbose,
    ...
  )
  list(
    chains = list(res),
    chain_stats = data.frame(
      chain         = 1L,
      accept_rate   = res$acceptance_rate,
      final_logpost = tail(res$post_logpost[is.finite(res$post_logpost)], 1L),
      stringsAsFactors = FALSE
    )
  )
}
