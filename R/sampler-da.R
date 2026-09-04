## R/sampler-da.R
## --------------------------------------------------------------------------
## Delayed-acceptance (two-stage) Metropolis-Hastings.
##
## Christen, J. A. & Fox, C. (2005), "Markov chain Monte Carlo using an
## approximation", J. Comput. Graph. Statist. 14(4), 795-810.
##
## The package owns, for the SAME model, a cheap approximate posterior
## (order-1 Kalman, a pruned lower order, the HANK interpolant, ...) and an
## expensive exact one (order-2/3 pruned, TPF, SV-RBPF, the HANK
## likelihoods). Delayed acceptance uses the cheap one only as a SCREEN:
##
##   stage 1  propose theta' ~ q(theta, .) (symmetric random walk);
##            accept with probability min(1, cheap(theta') / cheap(theta)).
##            A stage-1 rejection costs ZERO expensive evaluations.
##   stage 2  having promoted theta', evaluate the expensive posterior and
##            accept with probability
##              min(1, [pi(theta')/pi(theta)] * [cheap(theta)/cheap(theta')]).
##
## EXACTNESS. Write q1(theta, theta') = q(theta, theta') a1(theta, theta') for
## the stage-1 kernel. Because a1 is a Metropolis ratio for the cheap target
## and q is symmetric, q1 satisfies detailed balance w.r.t. `cheap`:
##   cheap(theta) q1(theta, theta') = cheap(theta') q1(theta', theta).
## The stage-2 acceptance is therefore exactly the Metropolis-Hastings ratio
## for target pi under proposal kernel q1,
##   a2 = min(1, [pi(theta') q1(theta', theta)] / [pi(theta) q1(theta, theta')])
##      = min(1, [pi(theta') cheap(theta)]   / [pi(theta) cheap(theta')]),
## so the composed chain is reversible w.r.t. pi. The screen may be as biased
## as it likes -- it changes only the EFFICIENCY of the chain, never its
## stationary distribution. (`cheap` must be strictly positive wherever pi is,
## which is why an infinite/NA screen value is treated as a stage-1 rejection
## rather than silently promoted.)
##
## PSEUDO-MARGINAL USE. When the expensive target is an unbiased likelihood
## ESTIMATOR (a particle filter), exactness additionally requires that the
## INCUMBENT's stored estimate is never refreshed: the auxiliary randomness
## of the accepted evaluation is part of the chain state. `rwmh_da()` stores
## both the incumbent's expensive AND cheap values and re-uses them; the
## deliberately-wrong `recompute_incumbent = TRUE` variant exists only so the
## test suite can demonstrate that its exactness oracles have power.
## --------------------------------------------------------------------------


## Coerce a target's return value to a scalar log-density. Accepts either a
## bare numeric (handy for cheap screens) or the package's
## list(logpost, loglik, logprior) contract. Anything non-finite or absent
## becomes -Inf, i.e. "reject", never a silent NA that would poison the
## acceptance comparison.
.da_logpost_value <- function(x) {
  v <- if (is.list(x)) x$logpost else x
  if (is.null(v) || length(v) < 1L) return(-Inf)
  v <- suppressWarnings(as.numeric(v)[1L])
  if (is.na(v)) -Inf else v
}


## Evaluate `screen_fn` twice at theta0 with the caller's RNG state saved and
## restored, and refuse a screen whose value moves. A stochastic screen breaks
## the stage-1 detailed-balance argument above (q1 would no longer be a fixed
## kernel), and the failure mode is a silently WRONG stationary distribution,
## not a diagnostic. Returns the (deterministic) screen value at theta0.
.da_check_screen <- function(screen_fn, theta0) {
  ge  <- globalenv()
  had <- exists(".Random.seed", envir = ge, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = ge, inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = ge)
    else if (exists(".Random.seed", envir = ge, inherits = FALSE))
      rm(list = ".Random.seed", envir = ge)
  }, add = TRUE)

  a <- .da_logpost_value(screen_fn(theta0))
  b <- .da_logpost_value(screen_fn(theta0))
  if (!identical(a, b))
    stop("rwmh_da: `screen_fn` is NOT deterministic in theta -- two calls at ",
         "theta0 returned ", format(a), " and ", format(b), ". Delayed ",
         "acceptance requires a deterministic stage-1 screen (the stage-1 ",
         "kernel must be a fixed q1); a noisy screen changes the chain's ",
         "stationary distribution without any diagnostic showing it. If the ",
         "screen is a particle filter, seed it (seed = <int>) -- the ",
         "EXPENSIVE target is the one that must stay unseeded.",
         call. = FALSE)
  if (!is.finite(a))
    stop("rwmh_da: `screen_fn` is not finite at theta0 (", format(a), "). ",
         "The screen must be positive wherever the target is.", call. = FALSE)
  a
}


## Walk a closure's enclosing environments looking for a bound, non-NULL
## `seed`, stopping at the first namespace / global environment. Used to give
## a fast, explicit error for a fixed-seed particle-filter closure before
## paying for the K-evaluation preflight.
.da_closure_seed <- function(fn) {
  if (!is.function(fn)) return(NULL)
  e <- environment(fn)
  for (depth in seq_len(20L)) {
    if (is.null(e) || identical(e, emptyenv())) return(NULL)
    if (nzchar(environmentName(e))) return(NULL)   # namespace / global / base
    if (exists("seed", envir = e, inherits = FALSE)) {
      s <- tryCatch(get("seed", envir = e, inherits = FALSE),
                    error = function(err) NULL)
      return(if (is.null(s)) NULL else s)
    }
    e <- parent.env(e)
  }
  NULL
}


## Pseudo-marginal preflight for the EXPENSIVE target: (a) refuse a closure
## that carries a fixed seed, (b) run the package's particle-MCMC variance
## preflight, which itself stops when all K evaluations come out bit-identical
## (the fixed-seed blind spot closed by the 2026-09 RNG-hygiene fix). The
## preflight calls set.seed() internally, so the caller's stream is saved and
## restored around it -- a preflight must not move the chain's random numbers.
.da_pf_preflight <- function(log_post_fn, theta0, K = 20L, verbose = FALSE) {
  s <- .da_closure_seed(log_post_fn)
  if (!is.null(s) && (is.numeric(s) || is.integer(s)))
    stop("rwmh_da: the EXPENSIVE log-posterior closure was built with a ",
         "fixed seed (seed = ", format(s[1L]), "), so its particle filter is ",
         "deterministic in theta. Pseudo-marginal delayed acceptance needs a ",
         "STOCHASTIC expensive target whose randomness is fresh at every ",
         "evaluation: rebuild it with seed = NULL (e.g. ",
         "make_log_posterior_tpf(..., seed = NULL)).", call. = FALSE)

  ge  <- globalenv()
  had <- exists(".Random.seed", envir = ge, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = ge, inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = ge)
    else if (exists(".Random.seed", envir = ge, inherits = FALSE))
      rm(list = ".Random.seed", envir = ge)
  }, add = TRUE)

  .tpf_pmcmc_preflight(log_post_fn, theta0, K = K, verbose = verbose)
}


## ---- Chain-state packs (R/chain-state.R) ---------------------------------
## DA chains save/resume through the PUBLIC pack API rather than the internal
## streaming checkpoint used by rwmh(): one `dynhr_chain_state` pack holds the
## position, the expensive log-posterior, the RNG stream and the full adapter
## state; a sibling list-.rds holds the draws and is grown with
## mcmc_chain_extend(). No new field is needed in mcmc_chain_state() -- the
## DA-specific pieces (cheap incumbent value, proposal covariance, per-draw
## acceptance history, expensive-evaluation counters) all live inside the
## free-form `adapt_state`, and the parameter-name fingerprint in `meta`.

.da_ckpt_paths <- function(dir, chain_id) {
  tag <- if (is.null(chain_id)) "1" else as.character(chain_id)
  list(state  = file.path(dir, paste0("da-chain-", tag, ".state.rds")),
       draws  = file.path(dir, paste0("da-chain-", tag, ".draws.rds")))
}


#' Delayed-acceptance (two-stage) Random Walk Metropolis-Hastings
#'
#' @param log_post_fn EXPENSIVE target. Function(theta) -> list(logpost, ...)
#'   or a bare numeric log-density. This is the distribution the chain
#'   targets exactly.
#' @param screen_fn   CHEAP screen. Function(theta) -> list(logpost, ...) or a
#'   bare numeric log-density. Must be DETERMINISTIC in theta (checked). Its
#'   bias is irrelevant to correctness.
#' @param theta0      Named initial parameter vector.
#' @param Sigma_prop  Proposal covariance (n_par x n_par).
#' @param n_draws     Total draws (including burn-in).
#' @param n_burn      Burn-in draws to discard from the returned chain.
#' @param scale       Initial proposal scale factor.
#' @param target_rate Target OVERALL acceptance rate for adaptive scaling.
#' @param adapt_every Adapt the scale every N draws (during burn-in only).
#' @param adapt_cov   Opt-in Haario et al. (2001) adaptive proposal covariance
#'   over the burn-in states, exactly as in \code{rwmh()}; frozen at the end
#'   of burn-in so the retained draws satisfy detailed balance.
#' @param pseudo_marginal Logical (default \code{FALSE}). Set \code{TRUE} when
#'   \code{log_post_fn} is a particle-filter closure: the sampler then refuses
#'   a fixed-seed closure and runs \code{.tpf_pmcmc_preflight()} at
#'   \code{theta0} (result returned in \code{$pf_preflight}).
#' @param preflight_K Evaluations used by that preflight (default 20).
#' @param check_screen Logical (default \code{TRUE}): verify that
#'   \code{screen_fn} is deterministic at \code{theta0}.
#' @param recompute_incumbent Logical (default \code{FALSE}). DELIBERATELY
#'   WRONG diagnostic variant: re-evaluate the incumbent's expensive target at
#'   every stage-2 test instead of re-using the stored value. With a
#'   deterministic expensive target this is merely wasteful; with a particle
#'   filter it destroys pseudo-marginal exactness. Exists so the test suite
#'   can show its exactness oracles have power. Never use it for inference.
#' @param checkpoint  Optional \code{list(dir=, flush_every=, resume=)}. Draws
#'   and a \code{\link{mcmc_chain_state}} pack are written to \code{dir} every
#'   \code{flush_every} draws; \code{resume = TRUE} continues that run
#'   bit-identically.
#' @param verbose     Print progress messages.
#' @param progressor  Optional progressr function (or NULL).
#' @param chain_id    Chain label for progress messages / checkpoint files.
#' @return List with the same contract as \code{rwmh()} -- \code{chain},
#'   \code{full_chain}, \code{logpost_trace}, \code{post_logpost},
#'   \code{acceptance_rate}, \code{scale}, \code{n_draws}, \code{n_burn},
#'   \code{block_accept_rate}, \code{elapsed_secs}, \code{sampler} -- plus
#'   \code{screen_rate} (share of proposals killed at stage 1),
#'   \code{n_expensive} (expensive evaluations performed, including the one
#'   at \code{theta0}), \code{n_promoted} (proposals that reached stage 2)
#'   and \code{pf_preflight}.
#' @noRd
rwmh_da <- function(log_post_fn, screen_fn, theta0, Sigma_prop,
                    n_draws = 10000L, n_burn = 5000L,
                    scale = 0.5, target_rate = 0.25,
                    adapt_every = 100L, adapt_cov = FALSE,
                    pseudo_marginal = FALSE, preflight_K = 20L,
                    check_screen = TRUE,
                    recompute_incumbent = FALSE,
                    checkpoint = NULL,
                    verbose = TRUE, progressor = NULL, chain_id = NULL) {

  if (!is.function(log_post_fn))
    stop("rwmh_da: `log_post_fn` must be a function.", call. = FALSE)
  if (!is.function(screen_fn))
    stop("rwmh_da: `screen_fn` must be a function (the cheap stage-1 screen).",
         call. = FALSE)
  if (!is.numeric(theta0) || length(theta0) < 1L)
    stop("rwmh_da: `theta0` must be a non-empty numeric vector.", call. = FALSE)
  if (!is.logical(adapt_cov) || length(adapt_cov) != 1L || is.na(adapt_cov))
    stop("rwmh_da: `adapt_cov` must be a single logical value.", call. = FALSE)
  n_draws <- as.integer(n_draws); n_burn <- as.integer(n_burn)
  if (n_draws < 2L)
    stop("rwmh_da: `n_draws` must be at least 2.", call. = FALSE)
  if (n_burn < 0L || n_burn >= n_draws)
    stop("rwmh_da: `n_burn` must satisfy 0 <= n_burn < n_draws.", call. = FALSE)

  n_par     <- length(theta0)
  par_names <- names(theta0)

  Sigma_curr <- Sigma_prop
  L <- .robust_chol(Sigma_curr, n_par)

  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)
  flush_every <- if (ckpt) as.integer(checkpoint$flush_every %||% 1000L) else NA_integer_
  ckpt_paths  <- if (ckpt) .da_ckpt_paths(checkpoint$dir, chain_id) else NULL

  chain           <- matrix(NA_real_, nrow = n_draws, ncol = n_par)
  colnames(chain) <- par_names
  logpost_trace   <- numeric(n_draws)
  accepted        <- logical(n_draws)

  state_hist  <- if (adapt_cov) matrix(NA_real_, nrow = n_burn, ncol = n_par) else NULL
  adapt_min_n <- max(200L, 10L * n_par)

  ## ---- Bookkeeping ------------------------------------------------------
  n_accept    <- 0L   # proposals accepted (both stages passed)
  n_screened  <- 0L   # proposals killed at stage 1
  n_promoted  <- 0L   # proposals that reached stage 2
  n_expensive <- 0L   # expensive evaluations performed
  pf_pre      <- NULL
  i_start     <- 1L
  t_start     <- Sys.time()

  state_curr <- theta0

  if (!ckpt_resume) {
    ## ---- Fresh run: prime both incumbents ONCE ---------------------------
    if (isTRUE(check_screen)) {
      lp_cheap_curr <- .da_check_screen(screen_fn, theta0)
    } else {
      lp_cheap_curr <- .da_logpost_value(screen_fn(theta0))
    }
    if (isTRUE(pseudo_marginal))
      pf_pre <- .da_pf_preflight(log_post_fn, theta0, K = preflight_K,
                                 verbose = verbose)

    lp_exp_curr <- .da_logpost_value(log_post_fn(theta0))
    n_expensive <- n_expensive + 1L
    if (!is.finite(lp_exp_curr))
      stop("rwmh_da: initial parameter vector has -Inf expensive ",
           "log-posterior. Check starting values.", call. = FALSE)

    chain[1, ]       <- state_curr
    logpost_trace[1] <- lp_exp_curr
    accepted[1]      <- TRUE
    if (adapt_cov && n_burn >= 1L) state_hist[1, ] <- state_curr

    if (ckpt) {
      if (!dir.exists(checkpoint$dir)) dir.create(checkpoint$dir, recursive = TRUE)
      unlink(c(ckpt_paths$state, ckpt_paths$draws))
    }
  } else {
    ## ---- Resume: restore position, adapter state and the RNG stream ------
    if (!file.exists(ckpt_paths$state) || !file.exists(ckpt_paths$draws))
      stop("rwmh_da: resume = TRUE but no saved DA chain was found in '",
           checkpoint$dir, "'.", call. = FALSE)
    ## allow_mid_adaptation: this sampler DOES resume the full adapter state
    ## (scale, Sigma, acceptance window, Haario history) exactly, which is the
    ## precondition mcmc_chain_restore() asks about.
    st <- mcmc_chain_restore(ckpt_paths$state, allow_mid_adaptation = TRUE,
                             restore_rng = TRUE)
    if (!identical(st$meta$par_names, par_names))
      stop("rwmh_da: the saved DA chain was written for parameters (",
           paste(st$meta$par_names, collapse = ", "),
           ") but theta0 carries (", paste(par_names, collapse = ", "),
           "). Refusing to resume.", call. = FALSE)
    ad <- st$adapt_state
    state_curr    <- st$position
    lp_exp_curr   <- st$lp
    lp_cheap_curr <- ad$lp_cheap
    scale         <- st$scales
    Sigma_curr    <- ad$Sigma_curr
    L             <- .robust_chol(Sigma_curr, n_par)
    n_accept      <- ad$n_accept
    n_screened    <- ad$n_screened
    n_promoted    <- ad$n_promoted
    n_expensive   <- ad$n_expensive
    n_burn        <- ad$n_burn        # the ORIGINAL burn-in fixes the split
    pf_pre        <- ad$pf_preflight
    i_start       <- st$sweep

    if (i_start > n_draws)
      stop("rwmh_da: the saved chain already has ", i_start, " draws; ",
           "`n_draws` = ", n_draws, " would shorten it.", call. = FALSE)

    obj <- readRDS(ckpt_paths$draws)
    old_z  <- obj$z[seq_len(i_start), , drop = FALSE]   # drop any post-state rows
    old_lp <- obj$lp[seq_len(i_start)]
    chain[seq_len(i_start), ] <- old_z
    logpost_trace[seq_len(i_start)] <- old_lp
    accepted[seq_len(i_start)] <- ad$accepted[seq_len(i_start)]
    if (adapt_cov) {
      if (is.null(ad$state_hist))
        stop("rwmh_da: adapt_cov = TRUE but the saved state carries no ",
             "Haario history. Resume with the same adapt_cov setting.",
             call. = FALSE)
      keep <- seq_len(min(i_start, n_burn))
      if (length(keep)) state_hist[keep, ] <- ad$state_hist[keep, , drop = FALSE]
    }
  }

  ## ---- Streaming helper: append the newly drawn rows and re-stamp the pack.
  ## Draws are written BEFORE the state, so `sweep` never exceeds what is on
  ## disk; the resume path truncates any extra rows.
  flush_from <- if (ckpt_resume) i_start + 1L else 1L
  .da_flush <- function(i, from) {
    rows <- seq.int(from, i)
    if (!length(rows)) return(invisible(NULL))
    ## A flush must be RNG-NEUTRAL: mcmc_chain_save() checksums through a
    ## tempfile(), and any random-number consumption here would make a
    ## checkpointed run diverge from an uncheckpointed one.
    ge  <- globalenv()
    old <- get(".Random.seed", envir = ge, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = ge), add = TRUE)
    if (!file.exists(ckpt_paths$draws)) {
      saveRDS(list(z = chain[rows, , drop = FALSE], lp = logpost_trace[rows],
                   sweep = length(rows), done = FALSE),
              ckpt_paths$draws)
    } else {
      mcmc_chain_extend(ckpt_paths$draws, chain[rows, , drop = FALSE],
                        new_lp = logpost_trace[rows], done = FALSE)
    }
    st <- mcmc_chain_state(
      position = state_curr, lp = lp_exp_curr, sweep = i,
      scales = scale,
      adapt_state = list(
        lp_cheap = lp_cheap_curr, Sigma_curr = Sigma_curr,
        n_accept = n_accept, n_screened = n_screened,
        n_promoted = n_promoted, n_expensive = n_expensive,
        n_burn = n_burn, accepted = accepted[seq_len(i)],
        state_hist = state_hist, pf_preflight = pf_pre),
      ## Adaptation is still live during burn-in; the flag advertises that
      ## honestly and the resume path above restores the adapter exactly.
      adapt_frozen = (i > n_burn),
      meta = list(sampler = "rwmh_da", par_names = par_names))
    mcmc_chain_save(st, ckpt_paths$state)
    invisible(NULL)
  }

  ## ---- Main loop --------------------------------------------------------
  for (i in i_start + seq_len(max(0L, n_draws - i_start))) {

    z          <- rnorm(n_par)
    state_prop <- state_curr + scale * as.numeric(L %*% z)
    names(state_prop) <- par_names

    ## -- Stage 1: the CHEAP screen. A stage-1 rejection costs no expensive
    ##    evaluation, which is the whole point of the algorithm.
    lp_cheap_prop <- .da_logpost_value(screen_fn(state_prop))
    log_alpha1    <- lp_cheap_prop - lp_cheap_curr
    promoted      <- is.finite(log_alpha1) && log(runif(1)) < log_alpha1

    if (!promoted) {
      n_screened <- n_screened + 1L
    } else {
      n_promoted <- n_promoted + 1L

      ## -- Stage 2: the EXPENSIVE target, corrected for the screen ratio.
      ##    The incumbent's expensive value is the STORED one: with a particle
      ##    filter, refreshing it would break pseudo-marginal exactness.
      lp_exp_prop <- .da_logpost_value(log_post_fn(state_prop))
      n_expensive <- n_expensive + 1L

      lp_exp_ref <- if (isTRUE(recompute_incumbent)) {
        ## DELIBERATELY WRONG variant (diagnostic only) -- see the argument doc.
        n_expensive <- n_expensive + 1L
        .da_logpost_value(log_post_fn(state_curr))
      } else {
        lp_exp_curr
      }

      log_alpha2 <- (lp_exp_prop - lp_exp_ref) -
                    (lp_cheap_prop - lp_cheap_curr)
      if (is.finite(log_alpha2) && log(runif(1)) < log_alpha2) {
        state_curr    <- state_prop
        lp_exp_curr   <- lp_exp_prop
        lp_cheap_curr <- lp_cheap_prop
        n_accept      <- n_accept + 1L
        accepted[i]   <- TRUE
      }
    }

    chain[i, ]       <- state_curr
    logpost_trace[i] <- lp_exp_curr

    if (adapt_cov && i <= n_burn) state_hist[i, ] <- state_curr

    ## ---- Haario adaptive covariance (burn-in only, then frozen) ---------
    if (adapt_cov && i %% adapt_every == 0 && i <= n_burn && i >= adapt_min_n) {
      Sigma_curr <- stats::cov(state_hist[1:i, , drop = FALSE]) +
        1e-8 * diag(n_par)
      L <- .robust_chol(Sigma_curr, n_par)
    }

    ## ---- Scalar scale adaptation on the OVERALL acceptance rate, frozen
    ## after burn-in exactly as rwmh() does.
    if (i %% adapt_every == 0 && i <= n_burn) {
      recent_rate <- mean(accepted[max(1, i - adapt_every + 1):i])
      if (recent_rate > target_rate + 0.05) scale <- scale * 1.1
      else if (recent_rate < target_rate - 0.05) scale <- scale * 0.9
    }

    if (ckpt && (i %% flush_every == 0 || i == n_draws)) {
      .da_flush(i, flush_from)
      flush_from <- i + 1L
    }

    if (i %% 1000 == 0 || i == n_draws) {
      rate    <- n_accept / (i - 1)
      scr     <- n_screened / (i - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / i * (n_draws - i)
      label   <- if (is.null(chain_id)) "?" else as.character(chain_id)
      msg <- sprintf(
        "DA Ch%s %d/%d accept=%.0f%% screen=%.0f%% nexp=%d lp=%.1f scale=%.3f ETA=%.0fs",
        label, i, n_draws, rate * 100, scr * 100, n_expensive,
        lp_exp_curr, scale, eta)
      if (!is.null(progressor)) progressor(message = msg, amount = 1)
      else if (verbose) cat("  ", msg, "\n")
    }
  }

  post_chain   <- chain[(n_burn + 1):n_draws, , drop = FALSE]
  post_logpost <- logpost_trace[(n_burn + 1):n_draws]
  n_prop       <- n_draws - 1L
  elapsed_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  list(
    chain             = post_chain,
    full_chain        = chain,
    logpost_trace     = logpost_trace,
    post_logpost      = post_logpost,
    acceptance_rate   = n_accept / n_prop,
    screen_rate       = n_screened / n_prop,
    n_expensive       = n_expensive,
    n_promoted        = n_promoted,
    pf_preflight      = pf_pre,
    scale             = scale,
    n_draws           = n_draws,
    n_burn            = n_burn,
    block_accept_rate = NA_real_,
    checkpoint_dir    = if (ckpt) checkpoint$dir else NULL,
    elapsed_secs      = elapsed_secs,
    sampler           = "rwmh_da"
  )
}
