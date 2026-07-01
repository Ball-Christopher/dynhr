## R/sampler-dime.R
## --------------------------------------------------------------------------
## DIME MCMC: Differential-Independence Mixture Ensemble sampler.
##
## Reference: Boehl (2022/2024) "DIME MCMC: A Swiss Army Knife for Bayesian
## Inference". Reference implementation: gboehl/dime_sampler (Python, MIT).
##
## Algorithm overview (per-iteration, per-walker i):
##   With prob aimh_prob: draw AIMH independence proposal from a running
##     multivariate Student-t fitted to the ensemble (Hastings correction).
##   Otherwise: DE move  q = x_i + gamma*(x_a - x_b) + noise
##     where a, b != i are random ensemble indices (ter Braak 2006).
##   Metropolis accept/reject, including Hastings term for AIMH.
##   After all walkers are updated: update running mean/cov with exponential
##   decay rho (logspace weighted average, as in the Python reference).
##
## Key parameter defaults (match reference exactly):
##   gamma    = 2.38 / sqrt(2 * n_par)
##   sigma    = 1e-5  (DE noise std)
##   aimh_prob = 0.1
##   rho      = 0.999 (exponential memory decay)
##   df       = 10    (Student-t df for AIMH)
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## Internal helpers
## --------------------------------------------------------------------------

#' Log-density of multivariate Student-t(df, mu, Sigma)
#'
#' @param x     numeric vector (evaluation point)
#' @param mu    numeric vector (location)
#' @param chol_S lower-triangular Cholesky factor L s.t. L %*% t(L) = Sigma
#' @param df    degrees of freedom
#' @return scalar log-density
#' @noRd
.dime_mvt_logpdf <- function(x, mu, chol_S, df) {
  d <- length(x)
  ## Standardised residual via forward-solve (no extra package needed)
  z <- forwardsolve(chol_S, x - mu)
  maha <- sum(z^2)
  ## log |det(L)| = sum(log(diag(L)))
  log_det <- sum(log(diag(chol_S)))
  ## log-Gamma ratio + normalising constant for d-dim Student-t
  lgamma((df + d) / 2) - lgamma(df / 2) -
    (d / 2) * log(df * pi) -
    log_det -
    ((df + d) / 2) * log1p(maha / df)
}


#' Sample from multivariate Student-t(df, mu, L %*% t(L))
#'
#' Uses the normal-chi2 representation for reproducibility across platforms
#' (matching the Python dime_sampler approach).
#'
#' @param mu    location vector (length d)
#' @param chol_S lower-triangular Cholesky factor
#' @param df    degrees of freedom
#' @return numeric vector of length d
#' @noRd
.dime_mvt_sample <- function(mu, chol_S, df) {
  d <- length(mu)
  z <- rnorm(d)
  ## chi-squared scaling: x = z / sqrt(chi2/df) ~ Student-t_df(0, I)
  chi2 <- rchisq(1, df = df)
  mu + as.numeric(chol_S %*% z) / sqrt(chi2 / df)
}


#' Robust Cholesky for DIME (reuses dynhr's .robust_chol if available)
#' @noRd
.dime_chol <- function(S) {
  ## try to use the package helper; fall back to plain chol with jitter
  ch <- tryCatch(
    .robust_chol(S, nrow(S)),
    error = function(e) NULL
  )
  if (!is.null(ch)) return(ch)
  ## manual jitter fallback
  d   <- nrow(S)
  eps <- 1e-8
  for (k in 1:10) {
    cc <- tryCatch(t(chol(S + eps * diag(d))), error = function(e) NULL)
    if (!is.null(cc)) return(cc)
    eps <- eps * 10
  }
  stop("DIME: Cholesky decomposition failed even with heavy jitter.")
}


#' logaddexp in R (log(exp(a) + exp(b)), numerically stable)
#' @noRd
.dime_logaddexp <- function(a, b) {
  mx <- max(a, b)
  if (!is.finite(mx)) return(-Inf)
  mx + log(exp(a - mx) + exp(b - mx))
}


## --------------------------------------------------------------------------
## Core serial sampler
## --------------------------------------------------------------------------

#' DIME MCMC ensemble sampler (serial)
#'
#' Differential-Independence Mixture Ensemble MCMC (Boehl 2022/2024).
#' Gradient-free, swarm-based sampler robust to multimodality.
#'
#' @param log_post_fn  Function(theta) -> list(logpost, loglik, logprior)
#' @param prior_spec   Prior specification (named list or data.frame). Used to
#'   draw initial ensemble from the prior. Also accepted: a prior_sampler
#'   function via \code{prior_sampler} argument.
#' @param n_chain      Number of ensemble walkers. Recommended: >= 5 * n_par.
#'   Minimum: 3. Default: NULL (set to max(5 * n_par, 20)).
#' @param n_iter       Post-warmup iterations (per-walker). Default 1000.
#' @param n_burn       Burn-in (warm-up) iterations. Default 500.
#' @param aimh_prob    Probability of the AIMH independence move. Default 0.1.
#' @param sigma        DE noise standard deviation. Default 1e-5.
#' @param rho          Exponential memory decay for running statistics.
#'   Default 0.999.
#' @param df           Student-t df for the AIMH proposal. Default 10.
#' @param prior_sampler Optional: function() -> named numeric vector from the
#'   prior. Constructed from prior_spec if NULL.
#' @param verbose      Print progress every 100 iterations. Default TRUE.
#' @param progressor   Optional progressr callback (or NULL).
#' @param checkpoint   Optional list for memory-streamed, restartable
#'   checkpointing. Fields: \code{dir} (string, checkpoint directory),
#'   \code{flush_every} (integer, approximate row count between flushes;
#'   rounded up to a whole iteration boundary so each flush is exactly
#'   \code{flush_iters * n_chain} rows), \code{resume} (logical),
#'   \code{fingerprint} (from \code{.ckpt_fingerprint}),
#'   \code{write_meta} (logical, default TRUE),
#'   \code{return_chain} (logical, default TRUE). When NULL (default), the
#'   non-checkpoint path runs byte-identically to before.
#' @return List: chain, full_chain, logpost_trace, post_logpost,
#'   acceptance_rate, n_draws, n_burn, n_walkers, n_iter, elapsed_secs,
#'   n_eval, sampler = "dime", checkpoint_dir
#' @noRd
run_dime <- function(log_post_fn,
                     prior_spec,
                     n_chain      = NULL,
                     n_iter       = 1000L,
                     n_burn       = 500L,
                     aimh_prob    = 0.1,
                     sigma        = 1e-5,
                     rho          = 0.999,
                     df           = 10,
                     prior_sampler = NULL,
                     verbose      = TRUE,
                     progressor   = NULL,
                     checkpoint   = NULL) {

  stopifnot(is.function(log_post_fn))
  n_iter <- as.integer(n_iter)
  n_burn <- as.integer(n_burn)

  ## ---- Checkpoint / streaming setup (opt-in) ------------------------------
  ## When checkpoint is a list, draws are streamed to disk in flush-sized chunks
  ## and the restart state is persisted after every flush. The flush boundary is
  ## aligned to WHOLE ITERATIONS so that each flush is exactly flush_iters *
  ## n_chain rows (we can only compute n_chain after probing the prior).
  ## The `chain_store` 3-D array is replaced by a flat (flush_iters * n_chain)
  ## row buffer, bounding RAM to flush_window * n_chain * n_par rather than
  ## total_iter * n_chain * n_par.
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)

  ## ---- Build prior sampler ------------------------------------------------
  if (is.null(prior_sampler)) {
    prior_sampler <- .smc_make_prior_sampler(prior_spec)
  }

  ## ---- Probe dimension and parameter names --------------------------------
  theta0    <- prior_sampler()
  n_par     <- length(theta0)
  par_names <- names(theta0)

  if (is.null(n_chain)) n_chain <- max(5L * n_par, 20L)
  n_chain <- as.integer(n_chain)
  if (n_chain < 3L) stop("run_dime: n_chain must be >= 3.")

  gamma <- 2.38 / sqrt(2 * n_par)

  ## ---- Iteration-aligned flush boundary -----------------------------------
  ## flush_every (rows) / n_chain rounded up gives flush_iters (whole iters).
  ## The actual flush row count is flush_iters * n_chain >= flush_every.
  flush_iters <- if (ckpt) {
    fe <- as.integer(checkpoint$flush_every %||% 1000L)
    max(1L, as.integer(ceiling(fe / n_chain)))
  } else NA_integer_

  ckpt_paths <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id = NULL) else NULL

  ## ---- Initialise ensemble from prior (retry until all finite) ------------
  ## On RESUME this is overwritten from the saved state below; we still probe
  ## the prior once to validate n_par before loading the state.
  if (!ckpt_resume) {
    if (verbose) message(sprintf("DIME: initialising %d walkers x %d params...",
                                 n_chain, n_par))
  }
  ensemble  <- matrix(NA_real_, nrow = n_chain, ncol = n_par)
  lp_vec    <- rep(-Inf, n_chain)
  max_tries <- 200L
  if (!ckpt_resume) {
    for (i in seq_len(n_chain)) {
      tries <- 0L
      repeat {
        th  <- prior_sampler()
        lp  <- tryCatch(log_post_fn(th)$logpost, error = function(e) -Inf)
        tries <- tries + 1L
        if (is.finite(lp) || tries >= max_tries) break
      }
      if (!is.finite(lp))
        stop(sprintf("DIME: walker %d could not find finite log-posterior in %d tries.",
                     i, max_tries))
      ensemble[i, ] <- th
      lp_vec[i]     <- lp
    }
  }
  colnames(ensemble) <- par_names

  ## ---- Initialise running AIMH proposal statistics and scalar bookkeeping ---
  ## For FRESH runs (both ckpt and non-ckpt): derive stats from the just-drawn
  ## ensemble. For RESUME: these are overwritten from the saved state below, so
  ## we must NOT call .dime_chol on the uninitialized ensemble here.
  if (!ckpt_resume) {
    prop_mean  <- colMeans(ensemble)
    prop_cov   <- stats::cov(ensemble) + diag(1e-8, n_par)
    chol_S     <- .dime_chol(prop_cov)
    ## cumlweight tracks the log-sum of (rho-decayed) ensemble weights.
    ## Initialised with log(n_chain) (one batch) then decay by rho.
    cumlweight <- log(n_chain) + log(rho)
    n_accept     <- 0L
    n_eval       <- n_chain   # initial evaluations above
    t_start_iter <- 1L
  } else {
    ## Placeholders: overwritten immediately in the resume block below.
    prop_mean  <- NULL; prop_cov <- NULL; chol_S <- NULL; cumlweight <- NULL
    n_accept   <- 0L;  n_eval   <- 0L;   t_start_iter <- 1L
  }
  t_start    <- Sys.time()
  total_iter <- n_burn + n_iter

  ## ---- Fresh-vs-resume branch --------------------------------------------
  if (ckpt) {
    if (ckpt_resume) {
      ## Verify the saved meta matches current configuration, then restore state.
      .ckpt_meta_verify(ckpt_paths$meta, "dime", checkpoint$fingerprint)
      st <- .ckpt_load_state(ckpt_paths$state)
      ## Restore all sampler state
      ensemble     <- st$ensemble
      lp_vec       <- st$lp_vec
      prop_mean    <- st$prop_mean
      prop_cov     <- st$prop_cov
      chol_S       <- .dime_chol(prop_cov)   # recompute from saved cov
      cumlweight   <- st$cumlweight
      n_accept     <- st$n_accept
      n_burn       <- st$n_burn              # original burn-in fixes retained set
      n_eval       <- st$n_eval
      t_start_iter <- st$n_done + 1L        # resume from the next iteration
      ## n_done is in ITERATIONS; truncate draws/lp files to n_done * n_chain rows
      .ckpt_truncate(ckpt_paths, st$n_done * n_chain, n_par)
      assign(".Random.seed", st$rng, envir = .GlobalEnv)
    } else {
      ## Fresh checkpoint run: clear stale files and write meta once.
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "dime", checkpoint$fingerprint)
    }
  }

  ## ---- Storage ------------------------------------------------------------
  if (ckpt) {
    ## Flush-sized buffer (flat rows: flush_iters * n_chain x n_par).
    ## Each complete iteration appends n_chain rows; flushed every flush_iters iters.
    buf_rows <- flush_iters * n_chain
    buf      <- matrix(NA_real_, nrow = buf_rows, ncol = n_par)
    buf_lp   <- numeric(buf_rows)
    buf_i    <- 0L   # rows filled so far (reset after each flush)
    chain_store <- NULL
    lp_store    <- NULL
  } else {
    ## Non-checkpoint: pre-allocate the full 3-D array exactly as before.
    chain_store <- array(NA_real_, dim = c(total_iter, n_chain, n_par))
    lp_store    <- matrix(NA_real_, nrow = total_iter, ncol = n_chain)
  }

  ## ---- Main sampling loop -------------------------------------------------
  for (t in t_start_iter:total_iter) {

    ## Proposals for all walkers this iteration
    proposals  <- matrix(NA_real_, nrow = n_chain, ncol = n_par)
    is_aimh    <- logical(n_chain)
    factors    <- numeric(n_chain)   # Hastings correction (0 for DE)

    for (i in seq_len(n_chain)) {
      if (runif(1) < aimh_prob) {
        ## --- AIMH independence move ---
        q_prop <- .dime_mvt_sample(prop_mean, chol_S, df)
        ## Hastings correction: log q(x_curr) - log q(x_prop)
        lq_curr <- .dime_mvt_logpdf(ensemble[i, ], prop_mean, chol_S, df)
        lq_prop <- .dime_mvt_logpdf(q_prop,        prop_mean, chol_S, df)
        factors[i]  <- lq_curr - lq_prop
        is_aimh[i]  <- TRUE
        proposals[i, ] <- q_prop
      } else {
        ## --- Differential evolution move ---
        ## draw two distinct indices a, b != i
        others <- seq_len(n_chain)[-i]
        ab     <- sample(others, 2L)
        noise  <- rnorm(n_par, 0, sigma)
        proposals[i, ] <- ensemble[i, ] + gamma * (ensemble[ab[1], ] - ensemble[ab[2], ]) + noise
        ## factors[i] stays 0 (symmetric proposal)
      }
    }

    ## Evaluate log-posteriors for all proposals
    lp_prop <- vapply(seq_len(n_chain), function(i) {
      th <- proposals[i, ]
      names(th) <- par_names
      r <- tryCatch(log_post_fn(th)$logpost, error = function(e) -Inf)
      if (is.na(r) || !is.finite(r)) -Inf else r
    }, numeric(1))
    n_eval <- n_eval + n_chain

    ## Metropolis accept/reject per walker
    for (i in seq_len(n_chain)) {
      log_alpha <- lp_prop[i] - lp_vec[i] + factors[i]
      if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        ensemble[i, ] <- proposals[i, ]
        lp_vec[i]     <- lp_prop[i]
        if (t > n_burn) n_accept <- n_accept + 1L
      }
    }

    ## Update running mean/cov (in log-space, matching Python reference).
    ## Must happen BEFORE the flush/state-save so the persisted `prop_mean`,
    ## `prop_cov`, and `cumlweight` reflect the post-iteration-t values.
    ## (On resume, iteration t+1 will use these stats for AIMH proposals,
    ## so they must match what a single non-checkpoint run would have at t+1.)
    nmean  <- colMeans(ensemble)
    ncov   <- stats::cov(ensemble) + diag(1e-8, n_par)
    lweight <- log(n_chain)
    newcumlweight <- .dime_logaddexp(cumlweight, lweight)
    w_old  <- exp(cumlweight - newcumlweight)
    w_new  <- exp(lweight    - newcumlweight)
    prop_cov  <- w_old * prop_cov  + w_new * ncov
    prop_mean <- w_old * prop_mean + w_new * nmean
    cumlweight <- newcumlweight + log(rho)
    ## Re-Cholesky every iteration (cheap for typical DSGE sizes)
    chol_S <- tryCatch(.dime_chol(prop_cov), error = function(e) chol_S)

    ## Store current ensemble state
    if (ckpt) {
      ## Append the n_chain rows for this iteration to the flat buffer.
      ## Rows within the buffer are ordered: walker 1..n_chain for iter t,
      ## then walker 1..n_chain for iter t+1, etc. -- matching the non-ckpt
      ## full_chain layout (iteration varies slowly).
      row_start <- buf_i + 1L
      row_end   <- buf_i + n_chain
      buf[row_start:row_end, ] <- ensemble
      buf_lp[row_start:row_end]  <- lp_vec
      buf_i <- buf_i + n_chain

      ## Flush when the buffer is full OR this is the very last iteration.
      ## State is saved AFTER draws are on disk (n_done = t in ITERATIONS).
      ## The running stats (prop_mean/prop_cov/cumlweight) are already updated
      ## above, so the saved state is fully consistent with iteration t complete.
      if (buf_i >= buf_rows || t == total_iter) {
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          ensemble   = ensemble,
          lp_vec     = lp_vec,
          prop_mean  = prop_mean,
          prop_cov   = prop_cov,
          cumlweight = cumlweight,
          n_accept   = n_accept,
          n_burn     = n_burn,
          n_eval     = n_eval,
          n_done     = t,             # completed iterations (not rows)
          n_iter_target = n_iter,
          rng        = get(".Random.seed", envir = .GlobalEnv)
        ))
        buf_i <- 0L
      }
    } else {
      for (i in seq_len(n_chain)) chain_store[t, i, ] <- ensemble[i, ]
      lp_store[t, ] <- lp_vec
    }

    ## Progress reporting
    if (verbose && (t %% 100 == 0 || t == total_iter)) {
      phase <- if (t <= n_burn) "burn" else "post"
      rate  <- if (t > n_burn) n_accept / max(1L, (t - n_burn) * n_chain) else NA
      msg   <- sprintf("DIME [%s] iter %d/%d  lp=%.1f  accept=%.0f%%",
                       phase, t, total_iter,
                       mean(lp_vec[is.finite(lp_vec)]),
                       if (is.na(rate)) NA else rate * 100)
      if (!is.null(progressor)) progressor(message = msg, amount = n_chain)
      else if (verbose) message(msg)
    }
  }

  ## ---- Assemble output ----------------------------------------------------
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  if (ckpt) {
    ## Materialize the full chain + lp from the streamed files.
    ## checkpoint$return_chain = FALSE skips the draw read-back (draws stay on disk).
    lp_full <- .ckpt_read_lp(ckpt_paths$lp)  # length total_iter * n_chain
    full_chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, n_par, par_names)
  } else {
    ## Non-checkpoint: assemble from the in-memory 3-D array exactly as before.
    post_idx <- seq(n_burn + 1L, total_iter)
    ## full_chain: all iterations x all walkers, row-major (iteration varies slowly)
    full_chain <- matrix(NA_real_, nrow = total_iter * n_chain, ncol = n_par)
    colnames(full_chain) <- par_names
    for (t in seq_len(total_iter)) {
      rows <- ((t - 1L) * n_chain + 1L):(t * n_chain)
      full_chain[rows, ] <- chain_store[t, , ]
    }
    lp_full <- as.numeric(t(lp_store))   # length total_iter * n_chain
  }

  ## Post-burn slice: rows corresponding to iterations (n_burn+1):total_iter.
  ## In the flat layout each iteration occupies n_chain consecutive rows.
  post_row_start <- n_burn * n_chain + 1L
  post_row_end   <- total_iter * n_chain
  post_chain   <- if (is.null(full_chain)) NULL else
    full_chain[post_row_start:post_row_end, , drop = FALSE]
  post_logpost <- lp_full[post_row_start:post_row_end]

  accept_rate <- n_accept / max(1L, n_iter * n_chain)

  list(
    chain           = post_chain,
    full_chain      = full_chain,
    logpost_trace   = lp_full,
    post_logpost    = post_logpost,
    acceptance_rate = accept_rate,
    n_draws         = if (is.null(post_chain)) n_iter * n_chain else nrow(post_chain),
    n_burn          = n_burn,
    n_walkers       = n_chain,
    n_iter          = n_iter,
    elapsed_secs    = elapsed,
    n_eval          = n_eval,
    sampler         = "dime",
    checkpoint_dir  = if (ckpt) checkpoint$dir else NULL
  )
}


## --------------------------------------------------------------------------
## Parallel variant: farm ensemble log-posterior evaluations to mirai pool
## --------------------------------------------------------------------------

#' Run DIME MCMC on a pre-provisioned mirai daemon pool.
#'
#' Like run_dime() but farms the per-iteration ensemble log-posterior
#' evaluations to the pool via mirai_map (one task per walker per iteration).
#' The accept/reject and running-statistics update run on the host, so only
#' log_post_fn calls are parallelised -- the same pattern as run_smc_mirai().
#'
#' The pool must already be provisioned (via .mirai_pool_init or
#' .mirai_pool_closure) and torn down by the caller.
#'
#' @param log_post_fn log-posterior closure (used only to evaluate proposals;
#'   daemons use their .worker_lp).
#' @param prior_sampler function() -> named numeric from prior.
#' @param n_chain,n_iter,n_burn,aimh_prob,sigma,rho,df,par_names same as run_dime().
#' @param seed_base  base seed; iteration t uses seed_base + t * n_chain + i.
#' @return same structure as run_dime().
#' @noRd
.run_dime_mirai_inner <- function(
    log_post_fn,
    prior_sampler,
    n_chain,
    n_iter,
    n_burn,
    aimh_prob   = 0.1,
    sigma       = 1e-5,
    rho         = 0.999,
    df          = 10,
    par_names   = NULL,
    seed_base   = 1L,
    verbose     = TRUE,
    checkpoint  = NULL
) {
  n_par  <- length(prior_sampler())
  if (is.null(par_names)) par_names <- names(prior_sampler())
  gamma  <- 2.38 / sqrt(2 * n_par)

  ## Seed the ORCHESTRATOR RNG from seed_base. The DIME ensemble MCMC (proposal
  ## generation, accept/reject, AIMH adaptation) runs host-side; seed_base
  ## previously only seeded the per-walker daemon lp-evals, so a run was
  ## reproducible only if the CALLER had set.seed() -- the seed_base argument did
  ## not actually control it. Seed it here so seed_base fully determines a fresh
  ## run (a resume restores the saved RNG below, overriding this).
  set.seed(seed_base)

  ## ---- Checkpoint / streaming setup (opt-in) ------------------------------
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)

  ## Iteration-aligned flush boundary (same logic as run_dime()).
  flush_iters <- if (ckpt) {
    fe <- as.integer(checkpoint$flush_every %||% 1000L)
    max(1L, as.integer(ceiling(fe / n_chain)))
  } else NA_integer_

  ckpt_paths <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id = NULL) else NULL

  ## ---- Initialise ensemble (serial, prior draws) --------------------------
  if (!ckpt_resume) {
    if (verbose)
      cat(sprintf("  DIME: initialising %d walkers (serial)...\n", n_chain))
  }
  ensemble <- matrix(NA_real_, nrow = n_chain, ncol = n_par)
  lp_vec   <- rep(-Inf, n_chain)
  if (!ckpt_resume) {
    for (i in seq_len(n_chain)) {
      tries <- 0L
      repeat {
        th  <- prior_sampler()
        lp  <- tryCatch(log_post_fn(th)$logpost, error = function(e) -Inf)
        tries <- tries + 1L
        if (is.finite(lp) || tries >= 200L) break
      }
      if (!is.finite(lp))
        stop(sprintf("DIME: walker %d: no finite log-posterior in 200 prior draws.", i))
      ensemble[i, ] <- th
      lp_vec[i]     <- lp
    }
  }
  colnames(ensemble) <- par_names

  ## ---- Running AIMH statistics -------------------------------------------
  ## For FRESH runs: derive from initial ensemble.
  ## For RESUME: overwritten below from saved state.
  if (!ckpt_resume) {
    prop_mean  <- colMeans(ensemble)
    prop_cov   <- stats::cov(ensemble) + diag(1e-8, n_par)
    chol_S     <- .dime_chol(prop_cov)
    cumlweight <- log(n_chain) + log(rho)
    n_accept   <- 0L
    n_eval     <- n_chain
    t_start_iter <- 1L
  } else {
    prop_mean  <- NULL; prop_cov <- NULL; chol_S <- NULL; cumlweight <- NULL
    n_accept   <- 0L;  n_eval   <- 0L;   t_start_iter <- 1L
  }

  ## ---- Fresh-vs-resume branch --------------------------------------------
  if (ckpt) {
    if (ckpt_resume) {
      .ckpt_meta_verify(ckpt_paths$meta, "dime", checkpoint$fingerprint)
      st <- .ckpt_load_state(ckpt_paths$state)
      ensemble     <- st$ensemble
      lp_vec       <- st$lp_vec
      prop_mean    <- st$prop_mean
      prop_cov     <- st$prop_cov
      chol_S       <- .dime_chol(prop_cov)
      cumlweight   <- st$cumlweight
      n_accept     <- st$n_accept
      n_burn       <- st$n_burn
      n_eval       <- st$n_eval
      t_start_iter <- st$n_done + 1L
      .ckpt_truncate(ckpt_paths, st$n_done * n_chain, n_par)
      assign(".Random.seed", st$rng, envir = .GlobalEnv)
    } else {
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "dime", checkpoint$fingerprint)
    }
  }

  ## ---- Storage ------------------------------------------------------------
  total_iter <- n_burn + n_iter
  if (ckpt) {
    buf_rows <- flush_iters * n_chain
    buf      <- matrix(NA_real_, nrow = buf_rows, ncol = n_par)
    buf_lp   <- numeric(buf_rows)
    buf_i    <- 0L
    chain_store <- NULL
    lp_store    <- NULL
  } else {
    chain_store <- array(NA_real_, dim = c(total_iter, n_chain, n_par))
    lp_store    <- matrix(NA_real_, nrow = total_iter, ncol = n_chain)
  }

  t_start <- Sys.time()

  ## ---- Main loop ----------------------------------------------------------
  for (t in t_start_iter:total_iter) {
    ## Build proposals on host (deterministic given RNG state)
    proposals <- matrix(NA_real_, nrow = n_chain, ncol = n_par)
    factors   <- numeric(n_chain)
    for (i in seq_len(n_chain)) {
      if (runif(1) < aimh_prob) {
        q_prop <- .dime_mvt_sample(prop_mean, chol_S, df)
        lq_curr <- .dime_mvt_logpdf(ensemble[i, ], prop_mean, chol_S, df)
        lq_prop <- .dime_mvt_logpdf(q_prop,        prop_mean, chol_S, df)
        factors[i]     <- lq_curr - lq_prop
        proposals[i, ] <- q_prop
      } else {
        others <- seq_len(n_chain)[-i]
        ab     <- sample(others, 2L)
        noise  <- rnorm(n_par, 0, sigma)
        proposals[i, ] <- ensemble[i, ] +
          gamma * (ensemble[ab[1], ] - ensemble[ab[2], ]) + noise
      }
    }

    ## Evaluate all proposals in parallel on the pool
    iter_seed <- seed_base + t * n_chain
    ## proposals / par_names / iter_seed are free variables -> pass via `...`
    eval_task <- function(i) {
      lpf <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      set.seed(iter_seed + i)
      th <- proposals[i, ]; names(th) <- par_names
      lp <- tryCatch(lpf(th)$logpost, error = function(e) -Inf)
      if (is.na(lp) || !is.finite(lp)) -Inf else lp
    }
    ## Sever env via a data-bound child of the dynhr namespace (see the SMC
    ## eval_task) so the frame isn't serialised per task.
    .ev_env <- new.env(parent = asNamespace("dynhr"))
    list2env(list(proposals = proposals, par_names = par_names,
                  iter_seed = iter_seed), envir = .ev_env)
    environment(eval_task) <- .ev_env
    ev <- mirai::mirai_map(seq_len(n_chain), eval_task)[]
    lp_prop <- vapply(ev, function(x) if (is.numeric(x)) x else -Inf, numeric(1))
    n_eval <- n_eval + n_chain

    ## Accept/reject on host
    for (i in seq_len(n_chain)) {
      log_alpha <- lp_prop[i] - lp_vec[i] + factors[i]
      if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        ensemble[i, ] <- proposals[i, ]
        lp_vec[i]     <- lp_prop[i]
        if (t > n_burn) n_accept <- n_accept + 1L
      }
    }

    ## Update running statistics.
    ## Must happen BEFORE the flush/state-save so the persisted prop_mean,
    ## prop_cov, and cumlweight reflect post-iteration-t values (matching
    ## run_dime()'s ordering exactly).
    nmean    <- colMeans(ensemble)
    ncov     <- stats::cov(ensemble) + diag(1e-8, n_par)
    lweight  <- log(n_chain)
    newcuml  <- .dime_logaddexp(cumlweight, lweight)
    w_old    <- exp(cumlweight - newcuml)
    w_new    <- exp(lweight    - newcuml)
    prop_cov   <- w_old * prop_cov  + w_new * ncov
    prop_mean  <- w_old * prop_mean + w_new * nmean
    cumlweight <- newcuml + log(rho)
    chol_S     <- tryCatch(.dime_chol(prop_cov), error = function(e) chol_S)

    ## Store current ensemble state
    if (ckpt) {
      row_start <- buf_i + 1L
      row_end   <- buf_i + n_chain
      buf[row_start:row_end, ] <- ensemble
      buf_lp[row_start:row_end]  <- lp_vec
      buf_i <- buf_i + n_chain

      if (buf_i >= buf_rows || t == total_iter) {
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          ensemble      = ensemble,
          lp_vec        = lp_vec,
          prop_mean     = prop_mean,
          prop_cov      = prop_cov,
          cumlweight    = cumlweight,
          n_accept      = n_accept,
          n_burn        = n_burn,
          n_eval        = n_eval,
          n_done        = t,
          n_iter_target = n_iter,
          rng           = get(".Random.seed", envir = .GlobalEnv)
        ))
        buf_i <- 0L
      }
    } else {
      for (i in seq_len(n_chain)) chain_store[t, i, ] <- ensemble[i, ]
      lp_store[t, ] <- lp_vec
    }

    if (verbose && (t %% 100 == 0 || t == total_iter)) {
      phase <- if (t <= n_burn) "burn" else "post"
      rate  <- if (t > n_burn) n_accept / max(1L, (t - n_burn) * n_chain) else NA
      cat(sprintf("  DIME [%s] iter %d/%d  lp=%.1f  accept=%.0f%%\n",
                  phase, t, total_iter,
                  mean(lp_vec[is.finite(lp_vec)]),
                  if (is.na(rate)) NA else rate * 100))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  if (ckpt) {
    lp_full <- .ckpt_read_lp(ckpt_paths$lp)
    full_chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, n_par, par_names)
  } else {
    post_idx  <- seq(n_burn + 1L, total_iter)
    full_chain <- matrix(NA_real_, nrow = total_iter * n_chain, ncol = n_par)
    colnames(full_chain) <- par_names
    for (tt in seq_len(total_iter)) {
      rows <- ((tt - 1L) * n_chain + 1L):(tt * n_chain)
      full_chain[rows, ] <- chain_store[tt, , ]
    }
    lp_full <- as.numeric(t(lp_store))
  }

  post_row_start <- n_burn * n_chain + 1L
  post_row_end   <- total_iter * n_chain
  post_chain <- if (is.null(full_chain)) NULL else
    full_chain[post_row_start:post_row_end, , drop = FALSE]
  lp_post <- lp_full[post_row_start:post_row_end]

  list(
    chain           = post_chain,
    full_chain      = full_chain,
    logpost_trace   = lp_full,
    post_logpost    = lp_post,
    acceptance_rate = n_accept / max(1L, n_iter * n_chain),
    n_draws         = if (is.null(post_chain)) n_iter * n_chain else nrow(post_chain),
    n_burn          = n_burn,
    n_walkers       = n_chain,
    n_iter          = n_iter,
    elapsed_secs    = elapsed,
    n_eval          = n_eval,
    sampler         = "dime",
    checkpoint_dir  = if (ckpt) checkpoint$dir else NULL
  )
}


#' Run DIME MCMC in parallel on a mirai daemon pool (compile-per-daemon).
#'
#' Provisions a mirai pool via .mirai_pool_init (standard Gaussian model) or
#' .mirai_pool_closure (OBC/PKF or pre-built closure), then calls
#' .run_dime_mirai_inner() which farms per-iteration ensemble evaluations to
#' the daemons via mirai_map.
#'
#' @param parsed_model parsed dynare model (NULL for closure-only path).
#' @param Y observation matrix (NULL for closure-only path).
#' @param prior_spec prior specification data.frame.
#' @param obs_names observed variable names (NULL for closure-only path).
#' @param n_chain number of ensemble walkers.
#' @param n_iter  post-burn iterations.
#' @param n_burn  burn-in iterations.
#' @param aimh_prob,sigma,rho,df sampler hyperparameters.
#' @param seed_base base RNG seed.
#' @param n_cores worker count (NULL = auto, capped at n_chain).
#' @param me_variance measurement-error variance.
#' @param me_extra n_obs x T matrix of per-period extra ME variances.
#' @param log_post_fn optional pre-built log-posterior closure.
#' @param verbose print progress.
#' @return same structure as run_dime().
#' @noRd
run_dime_mirai <- function(
    parsed_model = NULL, Y = NULL, prior_spec, obs_names = NULL,
    n_chain      = NULL,
    n_iter       = 1000L,
    n_burn       = 500L,
    aimh_prob    = 0.1,
    sigma        = 1e-5,
    rho          = 0.999,
    df           = 10,
    seed_base    = 1L,
    n_cores      = NULL,
    me_variance  = 0,
    me_extra     = NULL,
    shock_scale  = NULL,
    system_priors = NULL,
    lik_init     = "auto",
    tpf_options  = list(),
    gradient_policy = "auto",
    log_post_fn  = NULL,
    ctx          = NULL,
    verbose      = TRUE,
    checkpoint   = NULL
) {
  ## Unpack ctx fields when provided (ctx wins over individual args).
  if (!is.null(ctx) && inherits(ctx, "dynhr_estimation_context")) {
    me_variance     <- ctx$me_variance
    me_extra        <- ctx$me_extra
    shock_scale     <- ctx$shock_scale
    system_priors   <- ctx$system_priors
    lik_init        <- ctx$lik_init        %||% "auto"
    tpf_options     <- ctx$tpf_options     %||% list()
    gradient_policy <- ctx$gradient_policy %||% "auto"
  }
  ## Probe n_par to set n_chain default
  prior_sampler <- .smc_make_prior_sampler(prior_spec)
  theta_probe   <- prior_sampler()
  n_par         <- length(theta_probe)
  if (is.null(n_chain)) n_chain <- max(5L * n_par, 20L)
  n_chain <- as.integer(n_chain)

  ## Cap daemons at n_chain (can't usefully have more daemons than walkers)
  n_cores <- .mirai_n_cores(n_cores, n_chain)
  if (verbose)
    cat(sprintf("  Parallel DIME (mirai): %d walkers on %d daemons\n",
                n_chain, n_cores))

  t_init <- proc.time()
  sh <- NULL
  if (!is.null(log_post_fn)) {
    .mirai_pool_closure(n_cores, log_post_fn)
  } else {
    sh <- .mirai_pool_init(n_cores, parsed_model, Y, prior_spec, obs_names,
                           me_variance, me_extra = me_extra,
                           shock_scale = shock_scale,
                           system_priors = system_priors,
                           lik_init = lik_init,
                           tpf_options = tpf_options)
  }
  on.exit({ mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }, add = TRUE)
  if (verbose)
    cat(sprintf("  Daemon init: %.1f sec\n",
                (proc.time() - t_init)[["elapsed"]]))

  ## Need a concrete log_post_fn for the host-side initialisation.
  ## If not supplied, build a minimal closure that calls .worker_lp via
  ## a direct eval in the main session's environment.
  if (is.null(log_post_fn)) {
    ## .mirai_pool_init compiled the model; we also need a host-side closure.
    ## Build it the same way run_smc_mirai does (if not OBC, rebuild on host).
    if (!is.null(parsed_model) && !is.null(Y)) {
      .mk_lp  <- utils::getFromNamespace("make_log_posterior", "dynhr")
      .cmpl   <- utils::getFromNamespace("compile_model", "dynhr")
      .hcm    <- .cmpl(parsed_model, verbose = FALSE)
      ## tpf_options entries are spread as individual `...` args (matching
      ## .mirai_pool_init); make_log_posterior has no `tpf_options` parameter.
      log_post_fn <- do.call(.mk_lp,
        c(list(parsed_model, Y, prior_spec, obs_names,
               .hcm, me_variance = me_variance,
               me_extra = me_extra,
               lik_init = lik_init),
          tpf_options))
    } else {
      stop("run_dime_mirai: must supply either log_post_fn or parsed_model + Y.")
    }
  }

  res <- .run_dime_mirai_inner(
    log_post_fn   = log_post_fn,
    prior_sampler = prior_sampler,
    n_chain       = n_chain,
    n_iter        = n_iter,
    n_burn        = n_burn,
    aimh_prob     = aimh_prob,
    sigma         = sigma,
    rho           = rho,
    df            = df,
    par_names     = names(theta_probe),
    seed_base     = seed_base,
    verbose       = verbose,
    checkpoint    = checkpoint
  )
  res
}
