## R/sampler-rwmh.R
## --------------------------------------------------------------------------
## Phase-2 split from estimation-monolith.R.
##
## Random Walk Metropolis-Hastings (RWMH) sampler.
## Called by estimate_model() in estimate-monolith.R and by run_mcmc_parallel().
## --------------------------------------------------------------------------

#' Random Walk Metropolis-Hastings sampler
#'
#' @param log_post_fn  Function(theta) -> list(logpost, loglik, logprior)
#' @param theta0       Named initial parameter vector
#' @param Sigma_prop   Proposal covariance matrix (n_par x n_par). When
#'   `transform` is non-NULL this is interpreted as an ETA-SPACE (not
#'   theta-space) covariance -- see `transform` below.
#' @param n_draws      Total draws (including burn-in)
#' @param n_burn       Burn-in draws to discard from the returned chain
#' @param scale        Initial scaling factor for the proposal
#' @param target_rate  Target acceptance rate for adaptive scaling
#' @param adapt_every  Adapt scale every N draws (during burn-in)
#' @param verbose      Print progress messages
#' @param progressor   Optional progressr function (or NULL)
#' @param chain_id     Chain label for progress messages (or NULL)
#' @param transform    Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), the chain
#'   is run in UNCONSTRAINED eta-space:
#'   \itemize{
#'     \item `eta0 = transform$to_unconstrained(theta0)` is the starting
#'       state; proposals/adaptation operate on eta.
#'     \item the target is `lp_eta = make_transformed_logpost(log_post_fn,
#'       transform, include_jacobian = TRUE)`, i.e. the theta-space
#'       log-posterior PLUS the change-of-variables Jacobian, so eta-space
#'       draws are correctly distributed (Metropolis-Hastings ratios cancel
#'       the proposal density, so a symmetric eta-space random walk needs no
#'       extra correction beyond the Jacobian already folded into `lp_eta`).
#'     \item `Sigma_prop` (and its adaptive `scale`) are interpreted as an
#'       ETA-SPACE covariance -- the caller is responsible for converting a
#'       theta-space covariance via the delta method (see
#'       \code{run_posterior_estimation}'s `transform_params` path).
#'     \item each STORED chain row is `to_constrained(eta)` (theta-space, as
#'       always); `logpost_trace` stores the THETA-SPACE log-posterior (the
#'       Jacobian term is subtracted back out via one extra cheap call to
#'       `transform$log_jacobian()`), so traces remain directly comparable
#'       with the untransformed sampler's output.
#'   }
#'   When NULL (default), behaviour is bit-identical to before.
#' @param adapt_cov    Opt-in (default \code{FALSE}): Haario et al. (2001)
#'   adaptive proposal covariance. When \code{TRUE}, a history of
#'   SAMPLING-SPACE states (eta-space when `transform` is set, theta-space
#'   otherwise) is kept for the burn-in draws. Every `adapt_every` draws
#'   during burn-in, once at least `max(200, 10 * n_par)` states have been
#'   observed, the empirical covariance of the history so far is used to
#'   recompute the proposal Cholesky factor:
#'   \code{Sigma_emp = cov(state_hist) + 1e-8 * diag(n_par)} (Haario's
#'   epsilon regularisation, which keeps `Sigma_emp` non-degenerate even when
#'   a coordinate has not yet moved) and `L <- .robust_chol(Sigma_emp,
#'   n_par)`. The existing scalar `scale` adaptation continues to run on top
#'   of the adapted shape -- it tunes the overall magnitude, so Haario's
#'   classic fixed `2.38^2 / d` factor is subsumed by `scale` and is not
#'   applied separately. Adaptation stops at the end of burn-in (the
#'   covariance and `L` are frozen for the retained draws), so detailed
#'   balance holds for the post-burn-in chain and no diminishing-adaptation
#'   condition is required. When \code{FALSE} (default), no state history is
#'   allocated and no extra random numbers are drawn; behaviour is
#'   bit-identical to before.
#' @param n_blocks     Opt-in (default \code{1L}): randomized parameter
#'   blocking (Chib & Ramamurthy 2010; implementation follows Herbst &
#'   Schorfheide 2015, ch. 4). When `n_blocks > 1`, each iteration draws a
#'   random permutation of `1:n_par` (one `sample()` call) and splits it into
#'   `n_blocks` contiguous chunks (as equal in size as possible). For each
#'   block in turn, only those coordinates are proposed
#'   (`theta_prop[block] = theta_curr[block] + scale * L_b %*% z`, with `L_b
#'   = .robust_chol(Sigma_curr[block, block], length(block))`, recomputed
#'   per block per iteration -- cheap relative to a posterior evaluation),
#'   the FULL posterior is evaluated at the block-updated point, and the
#'   point is accepted/rejected; the state evolves WITHIN the iteration (a
#'   later block proposes from the post-earlier-block state). Each draw thus
#'   costs `n_blocks` posterior evaluations. `accepted[i]` records whether
#'   ANY block moved; the scalar `scale` adaptation instead uses the MEAN
#'   per-block acceptance rate (the any-block rate saturates near 1 as
#'   `n_blocks` grows, which would otherwise make `scale` explode). The mean
#'   block acceptance rate is returned as `block_accept_rate` (`NA` when
#'   `n_blocks == 1`). Blocks index whatever space the sampler is operating
#'   in (eta-space when `transform` is set) and, with `adapt_cov = TRUE`, use
#'   submatrices of the adapted `Sigma_curr`. When \code{1L} (default),
#'   behaviour is bit-identical to before.
#' @return List: chain, full_chain, logpost_trace, post_logpost,
#'         acceptance_rate, scale, n_draws, n_burn, block_accept_rate
#' @noRd
rwmh <- function(log_post_fn, theta0, Sigma_prop,
                 n_draws = 10000L, n_burn = 5000L,
                 scale = 0.5, target_rate = 0.25,
                 adapt_every = 100L, verbose = TRUE,
                 progressor = NULL, chain_id = NULL,
                 transform = NULL,
                 adapt_cov = FALSE, n_blocks = 1L,
                 checkpoint = NULL) {

  n_par     <- length(theta0)
  par_names <- names(theta0)

  if (!is.logical(adapt_cov) || length(adapt_cov) != 1L || is.na(adapt_cov))
    stop("adapt_cov must be a single logical value.")
  if (!is.numeric(n_blocks) || length(n_blocks) != 1L || is.na(n_blocks) ||
      n_blocks != as.integer(n_blocks) || n_blocks < 1L || n_blocks > n_par)
    stop("n_blocks must be a single integer in [1, n_par].")
  n_blocks <- as.integer(n_blocks)

  Sigma_curr <- Sigma_prop
  L <- .robust_chol(Sigma_curr, n_par)

  # ---- Checkpoint / streaming (opt-in). When `checkpoint` is a list carrying a
  # `dir`, draws are streamed to per-chain files in flush_every-row chunks (RAM
  # bounded by flush_every * n_par, not n_draws * n_par) and a restart state is
  # saved after every flush. `checkpoint$resume = TRUE` continues a prior run
  # from its saved state -- exactly (RNG, position, scale, proposal covariance),
  # so resume(N1)+resume(N2) is bit-identical to a single run of N1+N2.
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)
  flush_every <- if (ckpt) as.integer(checkpoint$flush_every %||% 1000L) else NA_integer_
  ckpt_paths  <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id) else NULL

  # ---- Transformed (opt-in) target: operate on eta = to_unconstrained(theta)
  if (!is.null(transform)) {
    lp_target  <- make_transformed_logpost(log_post_fn, transform, include_jacobian = TRUE)
    state0     <- transform$to_unconstrained(theta0)
    names(state0) <- par_names
  } else {
    lp_target  <- log_post_fn
    state0     <- theta0
  }

  # In checkpoint mode only a flush-sized buffer lives in RAM; the full chain is
  # read back from disk at the end for the return value. Otherwise pre-allocate
  # the full chain exactly as before (non-checkpoint path is byte-identical).
  if (ckpt) {
    buf    <- matrix(NA_real_, nrow = flush_every, ncol = n_par)
    buf_lp <- numeric(flush_every)
    buf_i  <- 0L
    chain  <- NULL
    logpost_trace <- NULL
  } else {
    chain         <- matrix(NA_real_, nrow = n_draws, ncol = n_par)
    colnames(chain) <- par_names
    logpost_trace <- numeric(n_draws)
  }
  accepted      <- logical(n_draws)

  # ---- Haario adaptive covariance (opt-in): history of SAMPLING-SPACE
  # states (eta-space when `transform` is set), preallocated to n_burn rows
  # since adaptation only runs during burn-in. When FALSE, no history is
  # allocated and no extra RNG draws occur anywhere below.
  state_hist <- if (adapt_cov) matrix(NA_real_, nrow = n_burn, ncol = n_par) else NULL
  adapt_min_n <- max(200L, 10L * n_par)

  # ---- Randomized blocking (opt-in): per-block acceptance bookkeeping.
  # `*_window` counters drive the scale adaptation (reset every adapt_every
  # draws); `*_total` counters accumulate over the whole run for the
  # reported `block_accept_rate`.
  block_accept_window <- 0L
  block_total_window  <- 0L
  block_accept_total  <- 0L
  block_total_total   <- 0L

  state_curr  <- state0
  result_curr <- lp_target(state_curr)
  lp_curr     <- result_curr$logpost
  if (!is.finite(lp_curr))
    stop("Initial parameter vector has -Inf log-posterior. Check starting values.")

  # Theta-space logpost for the trace (subtract the Jacobian back out when
  # transformed; identical to lp_curr otherwise).
  trace_lp_curr <- if (!is.null(transform)) {
    lp_curr - transform$log_jacobian(state_curr)
  } else {
    lp_curr
  }

  n_accept <- 0L
  i_start  <- 1L
  t_start  <- Sys.time()

  if (ckpt_resume) {
    # ---- Continue a saved run. First refuse a mismatched configuration (same
    # parameters are forced), then restore position, scale, proposal covariance,
    # accept count, burn-in and draw count, and -- crucially -- the RNG state,
    # so the continuation draws the SAME random stream a single long run would.
    .ckpt_meta_verify(ckpt_paths$meta, "rwmh", checkpoint$fingerprint)
    st <- .ckpt_load_state(ckpt_paths$state)
    state_curr    <- st$state_curr
    lp_curr       <- st$lp_curr
    trace_lp_curr <- st$trace_lp_curr
    scale         <- st$scale
    Sigma_curr    <- st$Sigma_curr
    L             <- .robust_chol(Sigma_curr, n_par)
    n_accept      <- st$n_accept
    n_burn        <- st$n_burn          # original burn-in fixes the retained set
    i_start       <- st$n_done
    .ckpt_truncate(ckpt_paths, st$n_done, n_par)  # drop any post-state partial flush
    assign(".Random.seed", st$rng, envir = .GlobalEnv)
  } else {
    # ---- Fresh run: record draw 1 (to the streaming buffer or the in-RAM chain).
    stored1 <- if (!is.null(transform)) transform$to_constrained(state_curr) else state_curr
    if (ckpt) {
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))    # clear any stale fresh-run files
      # The shared meta.rds is written once by the orchestrator on the parallel
      # path (write_meta = FALSE there) to avoid a multi-daemon write race.
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "rwmh", checkpoint$fingerprint)
      buf_i <- 1L; buf[1, ] <- stored1; buf_lp[1] <- trace_lp_curr
    } else {
      chain[1, ]       <- stored1
      logpost_trace[1] <- trace_lp_curr
    }
    accepted[1] <- TRUE
    if (adapt_cov) state_hist[1, ] <- state_curr
  }

  # Fresh: i = 2..n_draws (draw 1 already recorded). Resume: i = n_done+1..n_draws
  # (seq_len(0) when no additional draws are requested -> loop is skipped).
  for (i in i_start + seq_len(max(0L, n_draws - i_start))) {

    if (n_blocks == 1L) {
      ## ---- Unblocked path: identical to the pre-blocking code -----------
      z          <- rnorm(n_par)
      state_prop <- state_curr + scale * as.numeric(L %*% z)
      names(state_prop) <- par_names

      lp_prop   <- lp_target(state_prop)$logpost
      log_alpha <- lp_prop - lp_curr
      if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        state_curr <- state_prop; lp_curr <- lp_prop
        trace_lp_curr <- if (!is.null(transform)) {
          lp_curr - transform$log_jacobian(state_curr)
        } else {
          lp_curr
        }
        n_accept <- n_accept + 1L; accepted[i] <- TRUE
      }
    } else {
      ## ---- Randomized blocking (Chib & Ramamurthy 2010 / Herbst &
      ## Schorfheide 2015 ch. 4): one random permutation of 1:n_par per
      ## iteration, split into n_blocks contiguous (as-equal-as-possible)
      ## chunks. Each block proposes a move of only its coordinates against
      ## the FULL posterior; the state evolves within the iteration.
      perm   <- sample.int(n_par)
      ## Contiguous, as-equal-as-possible split of 1:n_par into n_blocks
      ## chunks (the first `n_par %% n_blocks` chunks get one extra element).
      base_size <- n_par %/% n_blocks
      n_extra   <- n_par %% n_blocks
      block_sizes <- rep(base_size, n_blocks) + c(rep(1L, n_extra), rep(0L, n_blocks - n_extra))
      block_grp   <- rep.int(seq_len(n_blocks), block_sizes)
      blocks      <- split(perm, block_grp)

      any_accept <- FALSE
      for (idx in blocks) {
        nb  <- length(idx)
        L_b <- .robust_chol(Sigma_curr[idx, idx, drop = FALSE], nb)
        z   <- rnorm(nb)

        state_prop <- state_curr
        state_prop[idx] <- state_curr[idx] + scale * as.numeric(L_b %*% z)
        names(state_prop) <- par_names

        lp_prop   <- lp_target(state_prop)$logpost
        log_alpha <- lp_prop - lp_curr

        block_total_window <- block_total_window + 1L
        block_total_total  <- block_total_total + 1L
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          state_curr <- state_prop; lp_curr <- lp_prop
          trace_lp_curr <- if (!is.null(transform)) {
            lp_curr - transform$log_jacobian(state_curr)
          } else {
            lp_curr
          }
          any_accept <- TRUE
          block_accept_window <- block_accept_window + 1L
          block_accept_total  <- block_accept_total + 1L
        }
      }

      if (any_accept) { n_accept <- n_accept + 1L; accepted[i] <- TRUE }
    }

    stored_i <- if (!is.null(transform)) transform$to_constrained(state_curr) else state_curr
    if (ckpt) {
      buf_i <- buf_i + 1L
      buf[buf_i, ]  <- stored_i
      buf_lp[buf_i] <- trace_lp_curr
      if (buf_i >= flush_every || i == n_draws) {
        # Flush the buffer to disk and persist the restart state (state written
        # atomically AFTER the draws, so n_done never exceeds what is on disk).
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          state_curr = state_curr, lp_curr = lp_curr, trace_lp_curr = trace_lp_curr,
          scale = scale, Sigma_curr = Sigma_curr, n_done = i, n_burn = n_burn,
          n_accept = n_accept, n_draws_target = n_draws,
          rng = get(".Random.seed", envir = .GlobalEnv)))
        buf_i <- 0L
      }
    } else {
      chain[i, ]       <- stored_i
      logpost_trace[i] <- trace_lp_curr
    }

    if (adapt_cov && i <= n_burn) state_hist[i, ] <- state_curr

    # ---- Haario adaptive covariance update (burn-in only; frozen
    # thereafter so the retained draws satisfy detailed balance under a
    # fixed proposal -- no diminishing-adaptation condition is needed).
    if (adapt_cov && i %% adapt_every == 0 && i <= n_burn && i >= adapt_min_n) {
      Sigma_emp <- stats::cov(state_hist[1:i, , drop = FALSE]) +
        1e-8 * diag(n_par)
      Sigma_curr <- Sigma_emp
      L <- .robust_chol(Sigma_curr, n_par)
    }

    # ---- Scalar scale adaptation. With blocking, the any-block acceptance
    # rate saturates near 1 as n_blocks grows (at least one of n_blocks
    # proposals is very likely to be accepted), which would drive `scale`
    # to explode; use the mean PER-BLOCK acceptance rate instead.
    if (i %% adapt_every == 0 && i <= n_burn) {
      if (n_blocks == 1L) {
        recent_rate <- mean(accepted[max(1, i - adapt_every + 1):i])
      } else {
        recent_rate <- if (block_total_window > 0L)
          block_accept_window / block_total_window else 0
        block_accept_window <- 0L
        block_total_window  <- 0L
      }
      if (recent_rate > target_rate + 0.05) scale <- scale * 1.1
      else if (recent_rate < target_rate - 0.05) scale <- scale * 0.9
    }

    if (i %% 1000 == 0 || i == n_draws) {
      rate    <- n_accept / (i - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / i * (n_draws - i)
      label   <- if (is.null(chain_id)) "?" else as.character(chain_id)
      msg <- sprintf("Ch%s %d/%d accept=%.0f%% lp=%.1f scale=%.3f ETA=%.0fs",
                     label, i, n_draws, rate * 100, trace_lp_curr, scale, eta)
      if (!is.null(progressor)) progressor(message = msg, amount = 1)
      else if (verbose) cat("  ", msg, "\n")
    }
  }

  if (ckpt) {
    # Materialize the full chain from the streamed files for the return value.
    # checkpoint$return_chain = FALSE skips this for very long runs (the draws
    # remain on disk; read them with .ckpt_read_draws() or resume to extend).
    logpost_trace <- .ckpt_read_lp(ckpt_paths$lp)
    chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, n_par, par_names)
  }
  post_chain   <- if (is.null(chain)) NULL else
    chain[(n_burn + 1):n_draws, , drop = FALSE]
  post_logpost <- logpost_trace[(n_burn + 1):n_draws]

  # ---- Mean per-block acceptance rate over the WHOLE run (for reporting);
  # NA when n_blocks == 1.
  block_accept_rate <- if (n_blocks == 1L) {
    NA_real_
  } else {
    block_accept_total / block_total_total
  }

  elapsed_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  list(
    chain             = post_chain,
    full_chain        = chain,
    logpost_trace     = logpost_trace,
    post_logpost      = post_logpost,
    acceptance_rate   = n_accept / (n_draws - 1),
    scale             = scale,
    n_draws           = n_draws,
    n_burn            = n_burn,
    block_accept_rate = block_accept_rate,
    checkpoint_dir    = if (ckpt) checkpoint$dir else NULL,
    elapsed_secs      = elapsed_secs,
    sampler           = "rwmh"
  )
}


## ---------------------------------------------------------------------------
## Correlated Pseudo-Marginal (CPM) RWMH sampler
##
## Implements Algorithm 1 of Deligiannidis, Doucet & Pitt (2018,
## J. Econometrics 206:799-829). At each step, the latent auxiliary variable
## U (the n_e x N standard-normal draws fed to the particle filter per period)
## is updated via an AR(1) (Crank-Nicolson) step jointly with theta, so that
## consecutive loglik evaluations are highly correlated. This dramatically
## reduces PMCMC loglik variance without changing the number of particles.
##
## Correctness: U is part of the JOINT state (theta, U). On acceptance BOTH
## theta and U move; on rejection BOTH stay. This ensures the chain targets
## the correct joint distribution whose theta-marginal is the posterior.
## ---------------------------------------------------------------------------

#' Correlated Pseudo-Marginal RWMH sampler (serial path only)
#'
#' @param log_post_fn  Function(theta, U_list = NULL) -> list(logpost, U_list).
#'   Must be the inner function returned by \code{make_log_posterior_tpf}.
#' @param theta0       Named initial parameter vector.
#' @param Sigma_prop   Proposal covariance (n_par x n_par).
#' @param n_draws      Total draws (including burn-in).
#' @param n_burn       Burn-in draws to discard.
#' @param scale        Initial proposal scale factor.
#' @param target_rate  Target acceptance rate for adaptive scaling.
#' @param adapt_every  Adapt scale every N draws (during burn-in).
#' @param rho_u        AR(1) correlation for U: U' = rho_u * U + sqrt(1 - rho_u^2) * Z.
#'   Range (0, 1). rho_u = 0 is independent (standard PMCMC); rho_u -> 1
#'   gives maximum correlation. Default 0.99.
#' @param verbose      Print progress messages.
#' @param progressor   Optional progressr function.
#' @param chain_id     Chain label for progress messages.
#' @return List: chain, full_chain, logpost_trace, post_logpost,
#'         acceptance_rate, scale, n_draws, n_burn, block_accept_rate.
#' @noRd
rwmh_cpm <- function(log_post_fn, theta0, Sigma_prop,
                     n_draws = 10000L, n_burn = 5000L,
                     scale = 0.5, target_rate = 0.25,
                     adapt_every = 100L,
                     rho_u = 0.99,
                     verbose = TRUE,
                     progressor = NULL, chain_id = NULL) {

  ## Input validation
  if (!is.numeric(rho_u) || length(rho_u) != 1L || !is.finite(rho_u) ||
      rho_u <= 0 || rho_u >= 1)
    stop("rwmh_cpm: rho_u must be a single numeric value in (0, 1).",
         " Got: ", rho_u, call. = FALSE)

  n_par     <- length(theta0)
  par_names <- names(theta0)

  L <- .robust_chol(Sigma_prop, n_par)

  chain         <- matrix(NA_real_, nrow = n_draws, ncol = n_par)
  colnames(chain) <- par_names
  logpost_trace <- numeric(n_draws)
  accepted      <- logical(n_draws)

  ## --- Prime U_curr and lp_curr from theta0 --------------------------------
  res0   <- log_post_fn(theta0, U_list = NULL)
  lp_curr <- res0$logpost
  if (!is.finite(lp_curr))
    stop("rwmh_cpm: Initial parameter vector has -Inf log-posterior. ",
         "Check starting values.", call. = FALSE)
  U_curr <- res0$U_list

  if (is.null(U_curr))
    stop("rwmh_cpm: log_post_fn did not return U_list. ",
         "Use make_log_posterior_tpf (TPF likelihood only).", call. = FALSE)

  ## Infer dimensions from U_curr.
  ## U_curr layout (3T+1 entries for Tier 10, 2T+1 for legacy Tier 9):
  ##   [[1]]:            init normals (n_s x N matrix)
  ##   [[2]].[[T+1]]:    per-period shock normals (n_e x N matrix each)
  ##   [[T+2]].[[2T+1]]: per-period phi=1 resampling z scalars (length-1 numeric each)
  ##   [[2T+2]].[[3T+1]]: per-period mid-stage z K-vectors (length-K numeric each)
  ## Note: legacy callers may supply only T+1 or 2T+1 entries; those are
  ## handled gracefully (missing slots treated as uncorrelated / fresh draws).
  n_U   <- length(U_curr)  # 3T+1 (or 2T+1 or T+1 for legacy)
  N     <- ncol(U_curr[[1L]])
  rho_c <- sqrt(1 - rho_u^2)  # complementary coefficient

  ## Classify each U slot by type for the AR(1) update:
  ##   "matrix" — init normals or shock normals (nrow x N)
  ##   "scalar" — phi=1 resampling z (length-1 numeric)
  ##   "kvec"   — mid-stage z K-vector (length-K numeric, K > 1)
  U_type <- vapply(U_curr, function(u) {
    if (is.matrix(u))                              "matrix"
    else if (is.numeric(u) && length(u) == 1L)    "scalar"
    else                                           "kvec"    # K-vector
  }, character(1L))

  ## The U_list=NULL priming call returns NA in the scalar/kvec z slots (the
  ## filter drew its uniforms internally and did not record a z). NA would
  ## propagate through the AR(1) update, silently disabling sorted-CPM resampling
  ## for the whole chain — draw fresh z ~ N(0,1) now for each NA-filled slot.
  U_curr <- lapply(U_curr, function(u) {
    if (!is.matrix(u) && is.numeric(u) && any(is.na(u)))
      rnorm(length(u))
    else u
  })

  ## Pre-record per-slot dimensions for the AR(1) update (only for matrix slots)
  U_dims <- lapply(U_curr, function(U_t) {
    if (is.matrix(U_t)) c(nrow(U_t), N) else c(length(U_t), 1L)
  })

  state_curr <- theta0
  chain[1, ]       <- state_curr
  logpost_trace[1] <- lp_curr
  accepted[1]      <- TRUE
  n_accept         <- 0L
  t_start          <- Sys.time()

  for (i in 2:n_draws) {

    ## -- Propose theta ---------------------------------------------------------
    z          <- rnorm(n_par)
    state_prop <- state_curr + scale * as.numeric(L %*% z)
    names(state_prop) <- par_names

    ## -- Propose U via AR(1) (pCN / autoregressive) update --------------------
    ## Matrix slots (init normals + shock normals): standard AR(1) on each element.
    ## Scalar slots (phi=1 resampling z): AR(1) on the scalar z ~ N(0,1) stored there;
    ##   the uniform used for resampling is u = pnorm(z), so correlating z gives
    ##   correlated resampling uniforms (Deligiannidis et al. sorted CPM).
    ## K-vector slots (mid-stage resampling z's): AR(1) on each element independently;
    ##   same N(0,1) marginal is preserved by the AR(1) update.
    U_prop <- mapply(function(U_t, dims, utype) {
      if (identical(utype, "matrix")) {
        nr <- dims[1L]
        rho_u * U_t + rho_c * matrix(rnorm(nr * N), nr, N)
      } else {
        ## Scalar z OR K-vector z: AR(1) update preserving N(0,1) marginal
        ## (length(U_t) == 1 for scalar, length(U_t) == K for kvec)
        rho_u * U_t + rho_c * rnorm(length(U_t))
      }
    }, U_curr, U_dims, U_type, SIMPLIFY = FALSE)

    ## -- Evaluate loglik at (theta_prop, U_prop) --------------------------------
    res_prop <- log_post_fn(state_prop, U_list = U_prop)
    lp_prop  <- res_prop$logpost

    ## -- Joint MH accept/reject (theta AND U move together) --------------------
    log_alpha <- lp_prop - lp_curr
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      state_curr <- state_prop
      lp_curr    <- lp_prop
      U_curr     <- res_prop$U_list  # CRITICAL: accept U_prop (not U_curr!)
      n_accept   <- n_accept + 1L
      accepted[i] <- TRUE
    }
    ## On rejection: state_curr, lp_curr, U_curr all stay unchanged.
    ## U_prop is discarded. This is the correct CPM stationary distribution.

    chain[i, ]       <- state_curr
    logpost_trace[i] <- lp_curr

    ## -- Adaptive scale (during burn-in) ---------------------------------------
    if (i %% adapt_every == 0L && i <= n_burn) {
      recent_rate <- mean(accepted[max(1L, i - adapt_every + 1L):i])
      if (recent_rate > target_rate + 0.05) scale <- scale * 1.1
      else if (recent_rate < target_rate - 0.05) scale <- scale * 0.9
    }

    if (i %% 1000 == 0 || i == n_draws) {
      rate    <- n_accept / (i - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / i * (n_draws - i)
      label   <- if (is.null(chain_id)) "?" else as.character(chain_id)
      msg <- sprintf("CPM Ch%s %d/%d accept=%.0f%% lp=%.1f scale=%.3f rho_u=%.3f ETA=%.0fs",
                     label, i, n_draws, rate * 100, lp_curr, scale, rho_u, eta)
      if (!is.null(progressor)) progressor(message = msg, amount = 1)
      else if (verbose) cat("  ", msg, "\n")
    }
  }

  post_chain   <- chain[(n_burn + 1):n_draws, , drop = FALSE]
  post_logpost <- logpost_trace[(n_burn + 1):n_draws]

  elapsed_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  list(
    chain             = post_chain,
    full_chain        = chain,
    logpost_trace     = logpost_trace,
    post_logpost      = post_logpost,
    acceptance_rate   = n_accept / (n_draws - 1),
    scale             = scale,
    n_draws           = n_draws,
    n_burn            = n_burn,
    block_accept_rate = NA_real_,
    elapsed_secs      = elapsed_secs,
    sampler           = "rwmh"
  )
}
