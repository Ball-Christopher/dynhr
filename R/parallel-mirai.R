## R/parallel-mirai.R
## --------------------------------------------------------------------------
## mirai + mori parallel backend for the estimation layer.
##
## A single persistent mirai daemon pool drives all three parallel consumers:
##   run_mcmc_mirai()  -- parallel RWMH chains
##   run_mode_mirai()  -- multi-start mode-finding
##   (SMC)             -- sampler-smc.R dispatches here for particle eval/mutation
##
## Two wins over the previous PSOCK (run_mcmc_parallel) and future
## (run_mode_parallel / future_lapply in sampler-smc.R) backends:
##
##   * Persistent daemons compile the model ONCE PER DAEMON via everywhere(),
##     instead of re-sourcing the package and recompiling the model on every
##     cluster init (PSOCK) or every task (future). The compiled derivative
##     closures are theta-independent and reused across all draws of all chains
##     on that daemon.
##   * The observation matrix Y is placed in OS shared memory once (mori::share)
##     and mapped zero-copy on each daemon (mori::map_shared), instead of being
##     serialised into every worker.
##
## The compiled model is a bundle of CLOSURES, which mori cannot share zero-copy
## (shared memory is for atomic vectors/matrices). So we ship the lightweight
## PARSED model through everywhere() and compile it on each daemon; only the
## atomic Y goes through mori. For SMC, where only the log-posterior CLOSURE is
## available (not the parsed model), everywhere() ships that closure once.
##
## Requires R 4.5.2 with mirai (>= 2.7.0) and mori (>= 0.2.0); see the
## mirai-mori-migration project note for the validated API surface.
## --------------------------------------------------------------------------


#' Resolve the daemon (worker) count.
#'
#' Mirrors the old PSOCK/future heuristic so worker counts are unchanged across
#' backends: detected logical cores, minus 2 above 4 cores, then capped at the
#' number of tasks.
#'
#' @param n_cores explicit worker count, or NULL to auto-resolve.
#' @param n_tasks number of tasks (chains / starts); caps the worker count.
#' @return a positive integer.
#' @noRd
.mirai_n_cores <- function(n_cores = NULL, n_tasks = NULL) {
  max_cores <- parallel::detectCores(logical = TRUE)
  if (is.null(n_cores)) {
    n_cores <- if (max_cores <= 4L) max_cores else max_cores - 2L
    n_cores <- max(1L, n_cores)
  }
  if (!is.null(n_tasks)) n_cores <- min(n_cores, n_tasks)
  max(1L, as.integer(n_cores))
}


## --------------------------------------------------------------------------
## Live progress: a nanonext back-channel + a single cli progress bar.
##
## mori shared memory is read-only (ALTREP / copy-on-write) and mirai_map's
## built-in .progress only reports whole-task completion, so daemons stream
## progress to the host over a nanonext push/pull socket on an ephemeral
## loopback port. The host renders ONE cli progress bar -- the only cursor
## model that behaves in RStudio (carriage-return only, no cursor-up), Windows
## conhost, real VT terminals, and piped logs alike: a fixed-width overall bar
## plus a width-aware per-chain token.
##
## MCMC streams the current draw every report interval (fine-grained bar).
## Mode-finding's optimizers (cmaes::cma_es, dfoptim::nmkb) are opaque blocking
## calls with no per-iteration hook, so each chain sends a single "done" ping
## and the bar advances task-by-task.
## --------------------------------------------------------------------------

#' Host: open a nanonext pull socket on an ephemeral loopback port.
#' @return list(sock, url); pass `url` to daemons, keep `sock` to drain.
#' @noRd
.progress_listener <- function() {
  sock <- nanonext::socket("pull")
  nanonext::listen(sock, url = "tcp://127.0.0.1:0", autostart = TRUE)
  list(sock = sock, url = nanonext::opt(sock$listener[[1L]], "url"))
}

#' A small fixed-width progress bar string.
#' @noRd
.progress_bar_str <- function(frac, w = 20L) {
  utf8 <- cli::is_utf8_output()
  on_ch  <- if (utf8) "\u2588" else "#"
  off_ch <- if (utf8) "\u00b7" else "-"
  k <- max(0L, min(w, as.integer(round(frac * w))))
  paste0(strrep(on_ch, k), strrep(off_ch, w - k))
}

#' Fit per-chain tokens to the available width: keep as many as fit, then "+k".
#' @noRd
.progress_fit_tokens <- function(toks, avail) {
  n <- length(toks)
  if (sum(nchar(toks)) + 2L * (n - 1L) <= avail)
    return(paste(toks, collapse = "  "))
  keep <- 0L; used <- 0L
  for (k in seq_along(toks)) {
    add <- nchar(toks[k]) + if (k > 1L) 2L else 0L
    if (used + add > avail - 6L) break
    used <- used + add; keep <- k
  }
  paste0(paste(toks[seq_len(max(1L, keep))], collapse = "  "), "  +", n - keep)
}

#' Render loop for parallel MCMC: one cli bar over total draws, per-chain bars.
#' Returns the collected results (handle[]).
#' @noRd
.render_mcmc_progress <- function(handle, prog_sock, n_chains, total_draws,
                                  render_hz = 8L) {
  draws <- integer(n_chains)                 # latest draw seen per chain
  total <- n_chains * total_draws
  env   <- environment()
  OVERALL_W <- 20L
  overall_bar <- .progress_bar_str(0, OVERALL_W)
  thread_str  <- ""

  refresh <- function(final = FALSE) {
    if (final) draws[] <- total_draws
    overall_bar <<- .progress_bar_str(sum(draws) / total, OVERALL_W)
    fr   <- draws / total_draws
    bars <- vapply(fr, function(f) .progress_bar_str(f, 6L), character(1))
    toks <- if (n_chains <= 8L)
      sprintf("C%d %s%3.0f%%", seq_len(n_chains), bars, 100 * fr)
    else sprintf("C%d %3.0f%%", seq_len(n_chains), 100 * fr)
    avail <- max(10L, cli::console_width() - (OVERALL_W + 34L))
    if (sum(nchar(toks)) + 2L * (n_chains - 1L) > avail)
      toks <- sprintf("C%d %3.0f%%", seq_len(n_chains), 100 * fr)
    thread_str <<- .progress_fit_tokens(toks, avail)
  }
  drain <- function() repeat {
    msg <- nanonext::recv(prog_sock, mode = "serial", block = FALSE)
    if (nanonext::is_error_value(msg)) break
    draws[msg[1L]] <<- msg[2L]
  }

  refresh()
  bar_id <- cli::cli_progress_bar(
    format = paste0("{cli::pb_spin} MCMC {overall_bar} {cli::pb_percent}  ",
                    "| {thread_str} | ETA {cli::pb_eta}"),
    format_done = paste0("{cli::col_green(cli::symbol$tick)} MCMC: ",
                         "{n_chains} chains x {total_draws} draws in {cli::pb_elapsed}."),
    total = total, clear = FALSE, .envir = env)

  delay <- 1 / render_hz
  while (any(mirai::unresolved(handle))) {
    drain(); refresh()
    cli::cli_progress_update(set = min(sum(draws), total - 1L),
                             id = bar_id, .envir = env)
    Sys.sleep(delay)
  }
  drain(); refresh(final = TRUE)
  cli::cli_progress_update(set = total, id = bar_id, .envir = env)
  try(cli::cli_progress_done(id = bar_id, .envir = env), silent = TRUE)
  handle[]
}

#' Render loop for multi-start mode-finding: per-chain progress from the
#' optimizer's evaluation count (fraction in [0,1]); a chain pings 1.0 when it
#' finishes (possibly before its full budget). cli bar total is scaled by 1000
#' so the overall bar advances smoothly with the mean chain fraction.
#' @noRd
.render_mode_progress <- function(handle, prog_sock, n_chains, render_hz = 6L) {
  frac <- numeric(n_chains)                    # latest fraction per chain [0,1]
  env  <- environment()
  SCALE <- 1000L
  OVERALL_W <- 20L
  overall_bar <- .progress_bar_str(0, OVERALL_W)
  thread_str  <- ""

  refresh <- function() {
    ## A chain pings 1.0 only on SUCCESS; one that errors freezes below 1, so a
    ## failed run honestly shows < 100% (never force-completed).
    overall_bar <<- .progress_bar_str(mean(frac), OVERALL_W)
    bars <- vapply(frac, function(f) .progress_bar_str(f, 6L), character(1))
    toks <- if (n_chains <= 8L)
      sprintf("C%d %s%3.0f%%", seq_len(n_chains), bars, 100 * frac)
    else sprintf("C%d %3.0f%%", seq_len(n_chains), 100 * frac)
    avail <- max(10L, cli::console_width() - (OVERALL_W + 34L))
    if (sum(nchar(toks)) + 2L * (n_chains - 1L) > avail)
      toks <- sprintf("C%d %3.0f%%", seq_len(n_chains), 100 * frac)
    thread_str <<- .progress_fit_tokens(toks, avail)
  }
  drain <- function() repeat {
    msg <- nanonext::recv(prog_sock, mode = "serial", block = FALSE)
    if (nanonext::is_error_value(msg)) break
    frac[msg[1L]] <<- max(frac[msg[1L]], msg[2L])   # monotone
  }

  refresh()
  bar_id <- cli::cli_progress_bar(
    format = paste0("{cli::pb_spin} Mode {overall_bar} {cli::pb_percent}  ",
                    "| {thread_str} | {cli::pb_elapsed}"),
    format_done = paste0("{cli::col_green(cli::symbol$tick)} Mode-finding: ",
                         "{sum(frac >= 0.999)}/{n_chains} chains succeeded in {cli::pb_elapsed}."),
    total = SCALE, clear = FALSE, .envir = env)

  delay <- 1 / render_hz
  while (any(mirai::unresolved(handle))) {
    drain(); refresh()
    cli::cli_progress_update(set = min(as.integer(mean(frac) * SCALE), SCALE - 1L),
                             id = bar_id, .envir = env)
    Sys.sleep(delay)
  }
  drain(); refresh()
  cli::cli_progress_update(set = as.integer(mean(frac) * SCALE),
                           id = bar_id, .envir = env)
  try(cli::cli_progress_done(id = bar_id, .envir = env), silent = TRUE)
  handle[]
}


#' Provision a daemon pool with the log-posterior built once per daemon.
#'
#' Starts \code{n_cores} mirai daemons, shares \code{Y} into OS shared memory,
#' and (on each daemon) loads dynhr, maps the shared \code{Y} zero-copy, compiles
#' the model once, and builds the log-posterior closure as \code{.worker_lp} in
#' the daemon's global environment. The caller is responsible for
#' \code{mirai::daemons(NULL)} teardown (typically via \code{on.exit}); the
#' returned share handle must be kept alive until the daemons have finished
#' reading.
#'
#' @param n_cores number of daemons.
#' @param parsed_model parsed (uncompiled) dynare model.
#' @param Y observation matrix.
#' @param prior_spec prior spec data.frame.
#' @param obs_names observed variable names.
#' @param me_variance measurement-error variance forwarded to make_log_posterior.
#' @param me_extra n_obs x T matrix of per-period extra ME variances (filter_tunes).
#' @param shock_scale n_exo x T matrix of per-period shock std scale factors.
#' Pin single-threaded BLAS for daemons spawned by the next mirai::daemons()
#'
#' Daemons inherit the parent's environment at launch, so we cap the BLAS thread
#' count in the parent here and restore it (via the returned closure) once the
#' pool is up. This matters most on macOS, whose default BLAS (Accelerate/vecLib)
#' uses Grand Central Dispatch and IGNORES OPENBLAS_NUM_THREADS: without
#' VECLIB_MAXIMUM_THREADS=1, each of n daemons spawns its own vecLib threads on
#' BLAS-heavy solves and oversubscribes the cores. On NZSIM's near-unit-root
#' mode-Hessian this was a 42x slowdown (1856s -> 44s on 4 cores); each daemon
#' already owns one core, so single-threaded BLAS per daemon is what we want.
#' @noRd
.mirai_blas_pin_vars <- c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS",
                          "OMP_NUM_THREADS", "MKL_NUM_THREADS")
.mirai_pin_blas_threads <- function() {
  old <- Sys.getenv(.mirai_blas_pin_vars, unset = NA_character_)
  do.call(Sys.setenv,
          stats::setNames(as.list(rep("1", length(.mirai_blas_pin_vars))),
                          .mirai_blas_pin_vars))
  function() for (v in .mirai_blas_pin_vars) {
    if (is.na(old[[v]])) Sys.unsetenv(v)
    else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
  }
}

#' @return the mori share handle for \code{Y} (keep alive until daemons finish).
#' @noRd
.mirai_pool_init <- function(n_cores, parsed_model, Y, prior_spec, obs_names,
                             me_variance = 0, me_extra = NULL,
                             shock_scale = NULL,
                             system_priors = NULL,
                             lik_init = "auto",
                             tpf_options = list(),
                             ctx = NULL) {
  ## Unpack ctx fields when provided (ctx wins over individual args).
  if (!is.null(ctx) && inherits(ctx, "dynhr_estimation_context")) {
    me_variance   <- ctx$me_variance
    me_extra      <- ctx$me_extra
    shock_scale   <- ctx$shock_scale
    system_priors <- ctx$system_priors
    lik_init      <- ctx$lik_init    %||% "auto"
    tpf_options   <- ctx$tpf_options %||% list()
  }
  .restore_blas <- .mirai_pin_blas_threads()
  on.exit(.restore_blas(), add = TRUE)
  mirai::daemons(n_cores)

  ## share() places Y in OS shared memory under an auto-generated UUID name
  ## (so repeated pools, mode -> mcmc, never collide on a stale segment) and
  ## returns a shared object that must be kept alive while daemons read it;
  ## shared_name() yields the key the daemons map by.
  sh  <- mori::share(Y)
  key <- mori::shared_name(sh)

  mirai::everywhere(
    {
      suppressMessages(library(dynhr))
      ## make_log_posterior is internal (non-exported); the `:::` operator does
      ## not resolve it reliably in the installed build, so reach it (and the
      ## compiler) via getFromNamespace, which does.
      .mk_lp <- utils::getFromNamespace("make_log_posterior", "dynhr")
      .cmpl  <- utils::getFromNamespace("compile_model", "dynhr")
      .worker_Y  <- mori::map_shared(key)
      .worker_cm <- .cmpl(parsed_model, verbose = FALSE)
      ## NOTE: `<<-`, not `<-`. A plain assignment inside everywhere() lands in
      ## the expression's local frame, which is reachable by bare-name lexical
      ## lookup but NOT by get0(envir = globalenv()). The mirai_map task bodies
      ## fetch this via get0(globalenv()), so the binding must reach the daemon's
      ## global environment -- `<<-` does that.
      ##
      ## tpf_options is a LIST of extra arguments forwarded to the TPF likelihood
      ## path via `...`. We use do.call() so the list items are spread as
      ## individual named arguments, not passed as a single `tpf_options=` arg.
      .worker_lp <<- do.call(.mk_lp,
        c(list(parsed_model, .worker_Y, prior_spec, obs_names,
               .worker_cm, me_variance = me_variance,
               me_extra = me_extra,
               shock_scale = shock_scale,
               system_priors = system_priors,
               lik_init = lik_init),
          tpf_options))
      ## Stash the ingredients make_posterior_grad() needs (parsed model,
      ## compiled model, mapped Y) as daemon globals too, so a chain task can
      ## build the analytic-gradient closure on demand without re-shipping or
      ## recompiling. Cheap: .worker_model is the lightweight parsed model
      ## (already a free variable here) and .worker_cm/.worker_Y are already
      ## built/mapped above for .worker_lp.
      .worker_model <<- parsed_model
      .worker_cm    <<- .worker_cm
      .worker_Y     <<- .worker_Y
    },
    .args = list(parsed_model = parsed_model, key = key,
                 prior_spec = prior_spec, obs_names = obs_names,
                 me_variance = me_variance, me_extra = me_extra,
                 shock_scale = shock_scale,
                 system_priors = system_priors,
                 lik_init = lik_init,
                 tpf_options = tpf_options)
  )[]
  ## ^ COLLECT (block) the everywhere() init. It recompiles the model on every
  ## daemon (seconds of work); leaving it uncollected returns while daemons are
  ## still initialising, and the in-flight init tasks' NNG contexts dangle if a
  ## later daemons() teardown races them -- the mechanism behind the intermittent
  ## 16-daemon dispatch hang (host pthread spins at 100%, daemons sit idle). The
  ## `[]` adds no wall-time (work queued after init already waits for it) and
  ## guarantees every daemon is ready before any task is dispatched.

  sh
}


#' Re-bind .worker_lp on an EXISTING daemon pool with a different lik_init.
#'
#' The expensive part of \code{.mirai_pool_init} is the per-daemon model
#' recompile (\code{.worker_cm}) and the shared-memory map of Y (\code{.worker_Y}).
#' The log-posterior closure \code{.worker_lp} -- which carries \code{lik_init} --
#' is cheap to rebuild from those. When successive mode-finding stages need
#' different filters (seed/Hessian at interior points want "stationary"; the
#' mode-search wants the robust "auto"), this lets a SINGLE pool serve all of
#' them: compile the model once, re-bind \code{.worker_lp} between stages, tear
#' the pool down once. Reuses the \code{.worker_model}/\code{.worker_cm}/
#' \code{.worker_Y} daemon globals set by \code{.mirai_pool_init}; does NOT
#' recompile or re-share. The caller must have a live pool (\code{.mirai_pool_init}
#' already run) and owns teardown.
#'
#' @return invisibly NULL.
#' @noRd
.mirai_rebind_worker_lp <- function(prior_spec, obs_names, me_variance = 0,
                                    me_extra = NULL, shock_scale = NULL,
                                    system_priors = NULL, lik_init = "auto",
                                    tpf_options = list()) {
  mirai::everywhere(
    {
      .mk_lp <- utils::getFromNamespace("make_log_posterior", "dynhr")
      ## Reuse the already-compiled model + mapped Y from the daemon globals;
      ## only the lik_init (and thus the filter setup) changes. `<<-` for the
      ## same global-env reason as in .mirai_pool_init.
      .worker_lp <<- do.call(.mk_lp,
        c(list(.worker_model, .worker_Y, prior_spec, obs_names,
               .worker_cm, me_variance = me_variance,
               me_extra = me_extra, shock_scale = shock_scale,
               system_priors = system_priors, lik_init = lik_init),
          tpf_options))
    },
    .args = list(prior_spec = prior_spec, obs_names = obs_names,
                 me_variance = me_variance, me_extra = me_extra,
                 shock_scale = shock_scale, system_priors = system_priors,
                 lik_init = lik_init, tpf_options = tpf_options)
  )[]
  invisible(NULL)
}


#' Run RWMH chains in parallel on a mirai daemon pool.
#'
#' Drop-in replacement for \code{run_mcmc_parallel()} (PSOCK). Each chain's task
#' reproduces the SEQUENTIAL path exactly -- \code{set.seed(seed_base + chain)},
#' compute the dispersed start, then draw with \code{rwmh()} continuing from that
#' RNG state -- so a mirai run is bit-for-bit identical to the sequential chains
#' under matched seeds. No \code{dynhr_files} sourcing: daemons load the
#' installed package and compile once via \code{\link{.mirai_pool_init}}.
#'
#' @param parsed_model parsed dynare model.
#' @param Y observation matrix.
#' @param prior_spec prior spec data.frame.
#' @param obs_names observed variable names.
#' @param theta_mode mode vector (chain 1 start; centre for chains 2..N).
#' @param Sigma_prop proposal covariance.
#' @param n_chains number of chains.
#' @param n_draws post-burn draws per chain.
#' @param n_burn burn-in draws per chain.
#' @param mh_scale initial RWMH scale.
#' @param target_accept target acceptance rate.
#' @param adapt_every adapt scale every N draws.
#' @param seed_base base RNG seed (chain k uses seed_base + k).
#' @param n_cores worker count (NULL = auto, capped at n_chains).
#' @param me_variance measurement-error variance.
#' @param me_extra n_obs x T matrix of per-period extra ME variances (filter_tunes).
#' @param log_post_fn optional pre-built log-posterior closure. When supplied,
#'   the pool ships this closure once via \code{\link{.mirai_pool_closure}}
#'   instead of recompiling a standard Gaussian posterior per daemon via
#'   \code{\link{.mirai_pool_init}} -- the path for OBC/PKF and cumulant
#'   models, where \code{parsed_model}/\code{Y}/\code{obs_names} may be NULL.
#' @param transform Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), each
#'   chain's \code{rwmh()} call receives \code{transform = transform} and
#'   \code{Sigma_prop} is interpreted as an ETA-SPACE covariance (the
#'   delta-method conversion is the caller's responsibility, as on the
#'   serial path -- see \code{run_posterior_estimation}'s
#'   \code{Sigma_prop_eta}). Dispersed starting points for chains 2..N are
#'   drawn in eta-space (\code{transform$to_unconstrained(theta_mode) +
#'   0.5 * chol(Sigma_prop) \%*\% z}, no boundary clamping -- eta has no
#'   boundary) and mapped back to theta-space via
#'   \code{transform$to_constrained()} for \code{rwmh()}'s \code{theta0}
#'   argument (which is always theta-space; \code{rwmh()} transforms
#'   internally). \code{transform} is a plain S3 list of closures and
#'   serializes to the daemons like any other free variable. Chain and
#'   \code{logpost_trace} are returned in theta-space, as on the serial path.
#' @param adapt_cov Opt-in (default \code{FALSE}); forwarded to each chain's
#'   \code{rwmh()} call -- Haario et al. (2001) adaptive proposal covariance
#'   (see \code{rwmh}'s `adapt_cov`).
#' @param n_blocks Opt-in (default \code{1L}); forwarded to each chain's
#'   \code{rwmh()} call -- randomized parameter blocking (see
#'   \code{rwmh}'s `n_blocks`).
#' @return list(chains, chain_stats, wall_time, n_cores) -- matches
#'   run_mcmc_parallel().
#' @noRd
run_mcmc_mirai <- function(
    parsed_model = NULL, Y = NULL, prior_spec, obs_names = NULL,
    theta_mode, Sigma_prop,
    n_chains      = 4L,
    n_draws       = 100000L,
    n_burn        = 50000L,
    mh_scale      = 1.50,
    target_accept = 0.25,
    adapt_every   = 200L,
    seed_base     = 42L,
    n_cores       = NULL,
    me_variance   = 0,
    me_extra      = NULL,
    shock_scale   = NULL,
    log_post_fn   = NULL,
    transform     = NULL,
    adapt_cov     = FALSE,
    n_blocks      = 1L,
    system_priors = NULL,
    lik_init      = "auto",
    tpf_options   = list(),
    gradient_policy = "auto",
    ctx           = NULL,
    checkpoint    = NULL,
    progress      = interactive()
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
  n_cores <- .mirai_n_cores(n_cores, n_chains)
  cat(sprintf("  Parallel MCMC (mirai): %d chains on %d daemons\n",
              n_chains, n_cores))
  cat(sprintf("  Per chain: %gk draws + %gk burn-in\n",
              n_draws / 1000, n_burn / 1000))

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
  cat(sprintf("  Daemon init: %.1f sec (load + compile + lp_fn)\n",
              (proc.time() - t_init)[["elapsed"]]))

  ## Live progress back-channel (host pull socket); daemons dial `prog_url`.
  prog_url <- NULL
  if (isTRUE(progress)) {
    pl <- .progress_listener(); prog_sock <- pl$sock; prog_url <- pl$url
    on.exit(try(nanonext::reap(prog_sock), silent = TRUE), add = TRUE)
  }

  total_draws <- n_draws + n_burn
  t_global <- Sys.time()

  ## One task per chain. The body mirrors the sequential .run_mcmc() loop body
  ## for exact RNG parity. .worker_lp is the daemon-local log-posterior.
  ##
  ## Per-chain seeds for TPF closures (Landmine 2):
  ## When log_post_fn is a TPF closure built with seed = NULL, the R RNG state
  ## set below drives each chain's particle paths independently -- chains are
  ## NOT identical.
  ## When log_post_fn is a TPF closure built with a non-NULL seed (frozen), the
  ## closure calls set.seed(seed) at each theta evaluation, making all chains'
  ## particle paths deterministic and identical.  For PMCMC independence, build
  ## the closure with seed = NULL (the preflight always does this).
  ## Full correlated pseudo-marginal (CPM, Deligiannidis et al.) with per-step
  ## correlated seeds is deferred -- see NEWS.md.
  chain_task <- function(ch) {
    log_post_fn <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    run_rwmh    <- utils::getFromNamespace("rwmh", "dynhr")
    ## mirai daemons default to the L'Ecuyer-CMRG generator; the sequential
    ## reference (and the rest of dynhr) uses Mersenne-Twister. Match it before
    ## seeding so a mirai chain reproduces the sequential chain bit-for-bit.
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed_base + ch)
    if (ch == 1L) {
      th0 <- theta_mode
    } else if (!is.null(transform)) {
      ## Disperse in eta-space (Sigma_prop here is Sigma_prop_eta) and map
      ## back to theta-space; eta has no boundary, so no clamping is needed.
      eta_mode <- transform$to_unconstrained(theta_mode)
      L   <- t(chol(Sigma_prop))
      z   <- rnorm(length(theta_mode))
      eta0 <- eta_mode + 0.5 * as.numeric(L %*% z)
      names(eta0) <- names(theta_mode)
      th0 <- transform$to_constrained(eta0)
    } else {
      L  <- t(chol(Sigma_prop))
      z  <- rnorm(length(theta_mode))
      th0 <- theta_mode + 0.5 * as.numeric(L %*% z)
      names(th0) <- names(theta_mode)
      for (i in seq_along(th0)) {
        th0[i] <- max(th0[i], prior_spec$lower[i] + 1e-6)
        th0[i] <- min(th0[i], prior_spec$upper[i] - 1e-6)
      }
    }
    lp0 <- log_post_fn(th0)$logpost
    if (!is.finite(lp0)) th0 <- theta_mode

    ## Build a progress reporter that streams the current draw to the host.
    ## rwmh calls progressor(message = "Ch<id> <i>/<n> ...") at each report
    ## interval; we parse <i> out and push c(ch, i). Side-effect only, so RNG
    ## parity with the sequential path is untouched. The push socket is opened
    ## once and reused; block=200 lets an early ping wait briefly for the pipe.
    progr <- NULL
    if (!is.null(prog_url)) {
      .psock <- nanonext::socket("push")
      nanonext::dial(.psock, url = prog_url, autostart = TRUE)
      on.exit(try(nanonext::reap(.psock), silent = TRUE), add = TRUE)
      progr <- function(message = NULL, amount = 1) {
        mm <- regmatches(message, regexec("([0-9]+)/([0-9]+)", message))[[1L]]
        if (length(mm) == 3L)
          nanonext::send(.psock, c(ch, as.integer(mm[2L])),
                         mode = "serial", block = 200)
      }
    }

    t0  <- proc.time()[["elapsed"]]
    res <- run_rwmh(log_post_fn, th0, Sigma_prop,
                    n_draws     = total_draws,
                    n_burn      = n_burn,
                    scale       = mh_scale,
                    target_rate = target_accept,
                    adapt_every = adapt_every,
                    verbose     = FALSE,
                    progressor  = progr,
                    chain_id    = ch,
                    transform   = transform,
                    adapt_cov   = adapt_cov,
                    n_blocks    = n_blocks,
                    checkpoint  = checkpoint)
    list(chain_id = ch, result = res,
         elapsed_min = (proc.time()[["elapsed"]] - t0) / 60)
  }

  ## Checkpoint: write the shared meta.rds ONCE here (the daemons set
  ## write_meta = FALSE to avoid a multi-writer race) and ship a per-daemon
  ## config. Each chain streams to its own chain_<id>.* files, so parallel
  ## writes never collide (parallel-safe by construction).
  ckpt_daemon <- NULL
  if (!is.null(checkpoint)) {
    if (!isTRUE(checkpoint$resume))
      .ckpt_meta_write(.ckpt_paths(checkpoint$dir)$meta, "rwmh", checkpoint$fingerprint)
    ckpt_daemon <- c(checkpoint, list(write_meta = FALSE))
  }

  ## chain_task references these as FREE variables (not formals), so they go
  ## through `...` -- mirai_map injects `...` objects into the task's evaluation
  ## environment, whereas `.args` would be matched as call arguments and error
  ## with "unused argument". `transform` (when non-NULL) is a plain S3 list of
  ## closures over local data and serializes like any other free variable.
  m_handle <- mirai::mirai_map(
    seq_len(n_chains), chain_task,
    seed_base = seed_base, theta_mode = theta_mode,
    Sigma_prop = Sigma_prop, prior_spec = prior_spec,
    total_draws = total_draws, n_burn = n_burn,
    mh_scale = mh_scale, target_accept = target_accept,
    adapt_every = adapt_every, prog_url = prog_url,
    transform = transform,
    adapt_cov = adapt_cov, n_blocks = n_blocks,
    checkpoint = ckpt_daemon
  )

  ## Collect: with a live bar, the render loop drains the socket while polling;
  ## otherwise just block on the results.
  raw <- if (isTRUE(progress))
    .render_mcmc_progress(m_handle, prog_sock, n_chains, total_draws)
  else m_handle[]

  wall_min <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))
  cat(sprintf("  All chains complete. Wall time: %.1f min\n", wall_min))

  chains <- vector("list", n_chains)
  chain_stats <- data.frame(
    chain = integer(), accept_rate = numeric(),
    final_logpost = numeric(), final_scale = numeric(),
    elapsed_min = numeric(), stringsAsFactors = FALSE
  )
  for (r in raw) {
    if (inherits(r, "miraiError") || inherits(r, "errorValue")) {
      cat(sprintf("  A chain FAILED: %s\n", as.character(r)))
      next
    }
    ch <- r$chain_id
    chains[[ch]] <- r$result
    chain_stats <- rbind(chain_stats, data.frame(
      chain         = ch,
      accept_rate   = r$result$acceptance_rate,
      final_logpost = tail(r$result$post_logpost, 1),
      final_scale   = r$result$scale,
      elapsed_min   = r$elapsed_min,
      stringsAsFactors = FALSE
    ))
  }

  list(chains = chains, chain_stats = chain_stats,
       wall_time = wall_min, n_cores = n_cores)
}


#' Run NUTS chains in parallel on a mirai daemon pool.
#'
#' K independent NUTS chains, one mirai task each, share a single daemon pool
#' that either compiles the (Gaussian) posterior once per daemon via
#' \code{\link{.mirai_pool_init}}, or -- when \code{log_post_fn} is supplied --
#' ships an already-built log-posterior closure once via
#' \code{\link{.mirai_pool_closure}} (the path for OBC/PKF and cumulant
#' models). Chain 1 starts at \code{theta_mode}; chains 2..N at dispersed
#' starts (\code{0.5 * chol(Sigma_prop) \%*\% z}, clamped to the prior box).
#' Every chain is preconditioned with a diagonal mass matrix
#' \code{1/diag(Sigma_prop)}, so \code{Sigma_prop} should approximate the
#' posterior covariance (e.g. the inverse-Hessian \code{V_mode}).
#'
#' @param parsed_model parsed dynare model.
#' @param Y observation matrix.
#' @param prior_spec prior spec data.frame.
#' @param obs_names observed variable names.
#' @param theta_mode mode vector (chain 1 start; centre for chains 2..N).
#' @param Sigma_prop posterior-covariance estimate (mass + dispersion source).
#' @param n_chains number of chains.
#' @param n_draws post-warmup draws per chain.
#' @param n_warmup warmup iterations per chain (step-size + mass adaptation).
#' @param max_treedepth maximum NUTS tree depth.
#' @param target_accept dual-averaging target acceptance rate.
#' @param seed_base base RNG seed (chain k uses seed_base + k).
#' @param n_cores worker count (NULL = auto, capped at n_chains).
#' @param me_variance measurement-error variance.
#' @param me_extra n_obs x T matrix of per-period extra ME variances (filter_tunes).
#' @param log_post_fn optional pre-built log-posterior closure. When supplied,
#'   the pool ships this closure once via \code{\link{.mirai_pool_closure}}
#'   instead of recompiling a standard Gaussian posterior per daemon via
#'   \code{\link{.mirai_pool_init}} -- the path for OBC/PKF and cumulant
#'   models, where \code{parsed_model}/\code{Y}/\code{obs_names} may be NULL.
#' @param transform Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), each
#'   chain's \code{dynhr_nuts()} call receives \code{transform = transform}
#'   (which expects a THETA-space \code{theta_init} and transforms
#'   internally -- see \code{dynhr_nuts}'s `transform` doc). \code{Sigma_prop}
#'   (and hence \code{mass_diag = 1/diag(Sigma_prop)}) is interpreted as an
#'   ETA-SPACE quantity -- the delta-method conversion from the theta-space
#'   proposal covariance is the caller's responsibility, as on the serial
#'   path (see \code{run_posterior_estimation}'s \code{Sigma_prop_eta} /
#'   eta-space \code{nuts_mass}). Dispersed starting points for chains 2..N
#'   are drawn in eta-space (\code{transform$to_unconstrained(theta_mode) +
#'   0.5 * chol(Sigma_prop) \%*\% z}, no boundary clamping) and mapped back to
#'   theta-space via \code{transform$to_constrained()} before being passed as
#'   \code{theta_init}. \code{transform} is a plain S3 list of closures over
#'   local data and serializes to the daemons like any other free variable.
#'   Chains and traces are returned in theta-space, as on the serial path.
#' @param analytic_grad Opt-in (default \code{FALSE}). When \code{TRUE} and
#'   the daemon pool was provisioned via \code{\link{.mirai_pool_init}} (i.e.
#'   \code{log_post_fn} is \code{NULL}, the standard-Gaussian path with a
#'   compiled model), each chain task builds an analytic/implicit gradient
#'   closure via \code{\link{make_posterior_grad}} from the daemon globals
#'   bound by \code{.mirai_pool_init} (\code{.worker_model}, \code{.worker_Y},
#'   \code{.worker_cm}), once per chain task, and passes it as \code{grad_fn}
#'   to \code{\link{dynhr_nuts}}. When \code{transform} is non-NULL, the
#'   gradient is composed with \code{\link{make_transformed_grad}} (chain-rule
#'   + log-Jacobian term), matching the serial path. Has no effect (silently
#'   ignored) when \code{log_post_fn} is supplied -- OBC/PKF and cumulant
#'   models have no analytic gradient.
#'
#'   Per-chain rebuild, not per-daemon caching: building the gradient closure
#'   does one model solve (and, for \code{grad_method = "implicit"}, one
#'   \code{solution_derivatives()} factorization) -- comparable to or cheaper
#'   than a single \code{dynhr_nuts} iteration -- and with typical
#'   \code{n_chains <= n_cores} each daemon serves one chain, so the rebuild
#'   cost is negligible relative to \code{n_draws + n_warmup} NUTS
#'   iterations. \code{.mirai_pool_init}'s \code{everywhere()} block runs once
#'   at pool provisioning, before chain assignment is known, so there is no
#'   natural place to build a per-daemon gradient cache there; a "first call
#'   on this daemon populates a global" cache inside the chain task itself
#'   would only help when \code{n_chains > n_cores} (multiple chains per
#'   daemon) and is left as a future optimisation if profiling shows it
#'   matters.
#' @param grad_method \code{"hybrid"} (default) or \code{"implicit"}; forwarded
#'   to \code{\link{make_posterior_grad}} when \code{analytic_grad = TRUE}.
#' @param progress show the live progress bar.
#' @return list(chains, chain_stats, wall_time, n_cores) -- matches
#'   run_mcmc_mirai (chain_stats also carries n_divergent, mean_treedepth).
#' @noRd
run_nuts_mirai <- function(
    parsed_model = NULL, Y = NULL, prior_spec, obs_names = NULL,
    theta_mode, Sigma_prop,
    n_chains      = 4L,
    n_draws       = 1000L,
    n_warmup      = 1000L,
    max_treedepth = 8L,
    target_accept = 0.80,
    seed_base     = 42L,
    n_cores       = NULL,
    me_variance   = 0,
    me_extra      = NULL,
    shock_scale   = NULL,
    log_post_fn   = NULL,
    transform     = NULL,
    analytic_grad = FALSE,
    grad_method   = c("hybrid", "implicit", "adjoint"),
    system_priors = NULL,
    lik_init      = "auto",
    tpf_options   = list(),
    gradient_policy = "auto",
    ctx           = NULL,
    checkpoint    = NULL,
    progress      = interactive()
) {
  grad_method <- match.arg(grad_method)
  ## Unpack ctx fields when provided (ctx wins over individual args).
  likelihood  <- "gaussian"
  freq_band   <- c(0, pi)
  if (!is.null(ctx) && inherits(ctx, "dynhr_estimation_context")) {
    me_variance     <- ctx$me_variance
    me_extra        <- ctx$me_extra
    shock_scale     <- ctx$shock_scale
    system_priors   <- ctx$system_priors
    likelihood      <- ctx$likelihood      %||% "gaussian"
    freq_band       <- ctx$freq_band       %||% c(0, pi)
    lik_init        <- ctx$lik_init        %||% "auto"
    tpf_options     <- ctx$tpf_options     %||% list()
    gradient_policy <- ctx$gradient_policy %||% "auto"
  }
  ## Guard centralised via .ctx_allows_analytic_gradient (me_extra/shock_scale
  ## are supported by the tv-aware tangent/adjoint path since Tier 7 item 3).
  if (isTRUE(analytic_grad)) {
    .tmp_ctx_nuts <- estimation_context(
      me_extra    = me_extra,
      shock_scale = shock_scale,
      likelihood  = likelihood
    )
    if (!.ctx_allows_analytic_gradient(.tmp_ctx_nuts)) {
      warning("run_nuts_mirai: analytic_grad ignored (likelihood/me_extra/shock_scale ",
              "combination is not supported by the analytic gradient path). Using the ",
              "numerical gradient.", call. = FALSE)
      analytic_grad <- FALSE
    }
    rm(.tmp_ctx_nuts)
  }
  n_cores <- .mirai_n_cores(n_cores, n_chains)
  cat(sprintf("  Parallel NUTS (mirai): %d chains on %d daemons\n",
              n_chains, n_cores))
  cat(sprintf("  Per chain: %d draws + %d warmup (max_treedepth=%d)\n",
              n_draws, n_warmup, max_treedepth))

  ## NOTE: when `transform` is non-NULL, Sigma_prop is the ETA-SPACE
  ## proposal covariance (already delta-method converted by the caller), so
  ## mass_diag here is correctly an eta-space mass matrix without further
  ## conversion -- mirroring the serial path's `nuts_mass / d_vec_nuts^2`.
  mass_diag <- 1 / pmax(diag(Sigma_prop), 1e-12)

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
  cat(sprintf("  Daemon init: %.1f sec (load + compile + lp_fn)\n",
              (proc.time() - t_init)[["elapsed"]]))

  prog_url <- NULL
  if (isTRUE(progress)) {
    pl <- .progress_listener(); prog_sock <- pl$sock; prog_url <- pl$url
    on.exit(try(nanonext::reap(prog_sock), silent = TRUE), add = TRUE)
  }

  total_draws <- n_draws + n_warmup
  t_global <- Sys.time()

  ## One task per chain. Mirrors run_mcmc_mirai for RNG/seed/dispersion parity.
  chain_task <- function(ch) {
    log_post_fn <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    run_nuts    <- utils::getFromNamespace("dynhr_nuts", "dynhr")
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed_base + ch)
    if (ch == 1L) {
      th0 <- theta_mode
    } else if (!is.null(transform)) {
      ## Disperse in eta-space (Sigma_prop here is Sigma_prop_eta) and map
      ## back to theta-space; eta has no boundary, so no clamping is needed.
      ## dynhr_nuts(..., transform = transform) expects a THETA-space
      ## theta_init and transforms internally.
      eta_mode <- transform$to_unconstrained(theta_mode)
      L   <- t(chol(Sigma_prop))
      z   <- rnorm(length(theta_mode))
      eta0 <- eta_mode + 0.5 * as.numeric(L %*% z)
      names(eta0) <- names(theta_mode)
      th0 <- transform$to_constrained(eta0)
    } else {
      L  <- t(chol(Sigma_prop))
      z  <- rnorm(length(theta_mode))
      th0 <- theta_mode + 0.5 * as.numeric(L %*% z)
      names(th0) <- names(theta_mode)
      for (i in seq_along(th0)) {
        th0[i] <- max(th0[i], prior_spec$lower[i] + 1e-6)
        th0[i] <- min(th0[i], prior_spec$upper[i] - 1e-6)
      }
    }
    lp0 <- log_post_fn(th0)$logpost
    if (!is.finite(lp0)) th0 <- theta_mode

    ## Analytic/implicit gradient (opt-in). Built ONCE per chain task from the
    ## daemon globals bound by .mirai_pool_init (.worker_model, .worker_Y,
    ## .worker_cm); see the `analytic_grad` doc above for the per-chain-vs-
    ## per-daemon tradeoff. Silently skipped (grad_fn stays NULL -> dynhr_nuts
    ## falls back to its numerical gradient) when the pool was provisioned via
    ## .mirai_pool_closure (log_post_fn supplied; no .worker_model/.worker_cm).
    grad_fn <- NULL
    if (isTRUE(analytic_grad)) {
      .worker_model <- get0(".worker_model", envir = globalenv(), inherits = FALSE)
      .worker_cm    <- get0(".worker_cm",    envir = globalenv(), inherits = FALSE)
      .worker_Y     <- get0(".worker_Y",     envir = globalenv(), inherits = FALSE)
      if (!is.null(.worker_model) && !is.null(.worker_cm) && !is.null(.worker_Y)) {
        .mk_grad <- utils::getFromNamespace("make_posterior_grad", "dynhr")
        base_grad_fn <- .mk_grad(.worker_model, .worker_Y, prior_spec, obs_names,
                                 .worker_cm, me_variance = me_variance,
                                 me_extra    = me_extra,
                                 shock_scale = shock_scale,
                                 grad_method = grad_method,
                                 likelihood  = likelihood,
                                 freq_band   = freq_band)
        grad_fn <- if (!is.null(transform)) {
          .mk_tgrad <- utils::getFromNamespace("make_transformed_grad", "dynhr")
          .mk_tgrad(base_grad_fn, transform)
        } else base_grad_fn
      }
    }

    progr <- NULL
    if (!is.null(prog_url)) {
      .psock <- nanonext::socket("push")
      nanonext::dial(.psock, url = prog_url, autostart = TRUE)
      on.exit(try(nanonext::reap(.psock), silent = TRUE), add = TRUE)
      progr <- function(message = NULL, amount = 1) {
        mm <- regmatches(message, regexec("([0-9]+)/([0-9]+)", message))[[1L]]
        if (length(mm) == 3L)
          nanonext::send(.psock, c(ch, as.integer(mm[2L])),
                         mode = "serial", block = 200)
      }
    }

    t0  <- proc.time()[["elapsed"]]
    res <- run_nuts(log_post_fn, th0,
                    n_draws       = n_draws,
                    n_warmup      = n_warmup,
                    mass_diag     = mass_diag,
                    max_treedepth = max_treedepth,
                    target_accept = target_accept,
                    verbose       = FALSE,
                    progressor    = progr,
                    chain_id      = ch,
                    transform     = transform,
                    grad_fn       = grad_fn,
                    checkpoint    = checkpoint)
    list(chain_id = ch, result = res,
         elapsed_min = (proc.time()[["elapsed"]] - t0) / 60)
  }

  ## Checkpoint: write the shared meta.rds ONCE here (daemons set write_meta =
  ## FALSE to avoid a multi-writer race); each chain streams to its own
  ## chain_<id>.* files (parallel-safe).
  ckpt_daemon <- NULL
  if (!is.null(checkpoint)) {
    if (!isTRUE(checkpoint$resume))
      .ckpt_meta_write(.ckpt_paths(checkpoint$dir)$meta, "nuts", checkpoint$fingerprint)
    ckpt_daemon <- c(checkpoint, list(write_meta = FALSE))
  }

  ## `transform` (when non-NULL) is a plain S3 list of closures over local
  ## data and serializes to the daemons like any other free variable here.
  m_handle <- mirai::mirai_map(
    seq_len(n_chains), chain_task,
    seed_base = seed_base, theta_mode = theta_mode, Sigma_prop = Sigma_prop,
    prior_spec = prior_spec, obs_names = obs_names, mass_diag = mass_diag,
    n_draws = n_draws, n_warmup = n_warmup, max_treedepth = max_treedepth,
    target_accept = target_accept, prog_url = prog_url, transform = transform,
    analytic_grad = analytic_grad, grad_method = grad_method,
    me_variance = me_variance, me_extra = me_extra, shock_scale = shock_scale,
    checkpoint = ckpt_daemon
  )

  raw <- if (isTRUE(progress))
    .render_mcmc_progress(m_handle, prog_sock, n_chains, total_draws)
  else m_handle[]

  wall_min <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))
  cat(sprintf("  All chains complete. Wall time: %.1f min\n", wall_min))

  chains <- vector("list", n_chains)
  chain_stats <- data.frame(
    chain = integer(), accept_rate = numeric(), final_logpost = numeric(),
    n_divergent = integer(), mean_treedepth = numeric(),
    elapsed_min = numeric(), stringsAsFactors = FALSE)
  for (r in raw) {
    if (inherits(r, "miraiError") || inherits(r, "errorValue")) {
      cat(sprintf("  A NUTS chain FAILED: %s\n", as.character(r)))
      next
    }
    ch <- r$chain_id
    chains[[ch]] <- r$result
    chain_stats <- rbind(chain_stats, data.frame(
      chain          = ch,
      accept_rate    = r$result$acceptance_rate %||% NA_real_,
      final_logpost  = tail(r$result$post_logpost, 1),
      n_divergent    = r$result$n_divergent %||% NA_integer_,
      mean_treedepth = r$result$mean_treedepth %||% NA_real_,
      elapsed_min    = r$elapsed_min, stringsAsFactors = FALSE))
  }

  list(chains = chains, chain_stats = chain_stats,
       wall_time = wall_min, n_cores = n_cores)
}


#' Multi-start parallel mode-finding on a mirai daemon pool.
#'
#' Drop-in replacement for \code{run_mode_parallel()} (future). Chain 1 starts at
#' \code{theta_init}; chains 2..N start at perturbations scaled by
#' \code{perturb_scale} times the prior std. Each task seeds with
#' \code{seed_base + chain*1000} then runs \code{.run_mode_finding()} over the
#' daemon-local log-posterior. With the explicit per-chain seed the results are
#' deterministic and reproducible (and match the old future backend, which also
#' set this seed explicitly inside each task).
#'
#' @param parsed_model parsed dynare model.
#' @param Y observation matrix.
#' @param prior_spec prior spec data.frame.
#' @param obs_names observed variable names.
#' @param theta_init starting vector (chain 1 = exact).
#' @param n_chains number of chains (NULL = n_cores).
#' @param nm_maxit mode-finding iteration budget per chain.
#' @param method optimizer name forwarded to .run_mode_finding().
#' @param perturb_scale perturbation as a fraction of prior std.
#' @param seed_base base RNG seed.
#' @param n_cores worker count (NULL = auto).
#' @param me_variance measurement-error variance.
#' @param me_extra n_obs x T matrix of per-period extra ME variances (filter_tunes).
#' @param log_post_fn optional pre-built log-posterior closure. When supplied,
#'   the pool ships this closure once via \code{\link{.mirai_pool_closure}}
#'   instead of recompiling a standard Gaussian posterior per daemon via
#'   \code{\link{.mirai_pool_init}} -- the path for OBC/PKF and cumulant
#'   models, where \code{parsed_model}/\code{Y}/\code{obs_names} may be NULL.
#' @param transform Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), each
#'   chain's \code{.run_mode_finding()} call receives \code{transform =
#'   transform} (Jacobian-free eta-space objective -- see
#'   \code{.run_mode_finding}'s `transform` doc; \code{theta_mode} is
#'   returned in theta-space as always). Dispersed starting points for
#'   chains 2..N are drawn as UNCONSTRAINED jitter in eta-space
#'   (\code{transform$to_unconstrained(theta_init) + rnorm(.,
#'   perturb_scale * prior_std_eta)}, no bound clamping -- eta has no
#'   boundary) and mapped back to theta-space via
#'   \code{transform$to_constrained()} before being handed to
#'   \code{.run_mode_finding()} as \code{theta0} (theta-space; it
#'   transforms internally). \code{transform} is a plain S3 list of
#'   closures over local data and serializes to the daemons like any other
#'   free variable.
#' @return list(chains, starts, logposts, best, best_chain, wall_time, n_cores)
#'   -- matches run_mode_parallel().
#' @noRd
run_mode_mirai <- function(
    parsed_model = NULL, Y = NULL, prior_spec, obs_names = NULL,
    theta_init,
    n_chains       = NULL,
    nm_maxit       = 5000L,
    method         = "cmaes_nmkb",
    perturb_scale  = 0.5,
    seed_base      = 42L,
    n_cores        = NULL,
    me_variance    = 0,
    me_extra       = NULL,
    shock_scale    = NULL,
    log_post_fn    = NULL,
    transform      = NULL,
    system_priors  = NULL,
    analytic_grad  = TRUE,
    grad_method    = "hybrid",
    likelihood     = "gaussian",
    freq_band      = NULL,
    newrat_H0_seed = NULL,
    lik_init       = "auto",
    pool_ready     = FALSE,
    progress       = interactive()
) {
  ## Analytic gradient needs the recompiled model on each daemon (.worker_model
  ## /.worker_cm), which only exist on the .mirai_pool_init path. A bare
  ## log_post_fn closure (.mirai_pool_closure) cannot build it, so disable.
  if (!is.null(log_post_fn)) analytic_grad <- FALSE
  n_cores <- .mirai_n_cores(n_cores, n_chains)
  if (is.null(n_chains)) n_chains <- n_cores
  n_cores <- min(n_cores, n_chains)

  par_names <- names(theta_init)
  prior_sds <- setNames(prior_spec$std,   prior_spec$name)
  prior_lo  <- setNames(prior_spec$lower, prior_spec$name)
  prior_hi  <- setNames(prior_spec$upper, prior_spec$name)

  ## Per-parameter dispersion scale. By default the jitter is a fraction of the
  ## PRIOR std -- but on a well-identified model the posterior is far tighter
  ## than the prior, so prior-std jitter throws restarts clean out of the mode's
  ## basin into BK-violating / non-PD regions where newrat converges to garbage
  ## (or fails on the first eval), wasting the whole multi-start. When a seed
  ## Hessian is available we instead disperse by the POSTERIOR std,
  ## sqrt(diag(H0^-1)), capped at the prior std so we never exceed prior support.
  ## This keeps restarts inside the high-probability region where newrat reliably
  ## reconverges -- turning the multi-start into a genuine robustness check
  ## (do nearby starts return to the same mode?) rather than a basin lottery.
  jitter_sds <- prior_sds
  if (!is.null(newrat_H0_seed)) {
    ## The seed is the neg-logpost Hessian, often NOT PD at theta_init (the prior
    ## mean is far from the mode, cond ~1e11), so a raw solve() gives negative /
    ## non-finite variances. .make_pd floors the eigenvalues to a target cond and
    ## falls back to prior_var on bad coordinates, returning a PD covariance whose
    ## diagonal is the (regularised) posterior variance per parameter.
    post_sd <- tryCatch({
      Sig <- .make_pd(newrat_H0_seed, cond_target = 100,
                      prior_var = (prior_sds[par_names])^2)
      v <- diag(Sig)
      if (all(is.finite(v)) && all(v > 0)) setNames(sqrt(v), par_names) else NULL
    }, error = function(e) NULL)
    if (!is.null(post_sd))
      jitter_sds <- pmin(post_sd, prior_sds[par_names], na.rm = TRUE)
  }

  ## Dispersed starts (identical logic to run_mode_parallel for parity).
  ## When `transform` is non-NULL, jitter is applied in ETA-SPACE (no bound
  ## clamping -- eta has no boundary) and mapped back to theta-space.
  starts <- vector("list", n_chains)
  starts[[1L]] <- theta_init
  eta_init <- if (!is.null(transform)) transform$to_unconstrained(theta_init) else NULL
  ## jitter_sds is a THETA-space sd. In the eta branch we jitter eta, so convert
  ## to an eta-space sd via the local Jacobian: deta ~= dtheta / (dtheta/deta).
  ## (Near a bound dtheta/deta -> 0, so the eta sd grows -- correct: a large eta
  ## step there maps to a small, in-bounds theta step.)
  eta_jitter_sds <- jitter_sds
  if (!is.null(transform)) {
    J <- tryCatch(abs(as.numeric(transform$dtheta_deta(eta_init))),
                  error = function(e) NULL)
    if (!is.null(J) && length(J) == length(par_names)) {
      names(J) <- par_names
      J[!is.finite(J) | J <= 0] <- 1
      eta_jitter_sds <- jitter_sds[par_names] / J
    }
  }
  for (ch in seq_len(n_chains)[-1L]) {
    set.seed(seed_base + ch)
    if (!is.null(transform)) {
      eta <- eta_init
      for (nm in par_names) {
        ps <- eta_jitter_sds[nm]
        if (!is.na(ps) && ps > 0) {
          eta[nm] <- eta[nm] + rnorm(1L, 0, perturb_scale * ps)
        }
      }
      starts[[ch]] <- transform$to_constrained(eta)
    } else {
      th <- theta_init
      for (nm in par_names) {
        ps <- jitter_sds[nm]
        if (!is.na(ps) && ps > 0) {
          th[nm] <- th[nm] + rnorm(1L, 0, perturb_scale * ps)
          lo <- prior_lo[nm]; hi <- prior_hi[nm]
          if (!is.na(lo) && is.finite(lo)) th[nm] <- max(th[nm], lo + 1e-6)
          if (!is.na(hi) && is.finite(hi)) th[nm] <- min(th[nm], hi - 1e-6)
        }
      }
      starts[[ch]] <- th
    }
  }

  ## Per-chain optimiser assignment.
  ##
  ## A SINGLE `method` (the common case) => pure MULTI-START of that method:
  ## every chain runs the same optimiser, chain 1 from theta_init (the exact
  ## serial replica) and chains 2..N from dispersed starts, keep-best. For the
  ## seeded trust-region methods (newrat / cmaes_newrat) this is the robust
  ## workhorse on near-unit-root DSGE models: the shared H0 seed + analytic
  ## gradient converge fast (~40 s/chain) and the dispersed restarts guard
  ## against a single bad basin -- at a per-chain cost no larger than the serial
  ## run it replaces.
  ##
  ## We deliberately do NOT fold the global / population searchers (cmaes*, jade,
  ## nmkb, nelder, combined) into the single-method default. On near-unit-root
  ## models their line search has no feasible-region backtracking, so they (a)
  ## stall at BK-violating points and return garbage logposts orders of magnitude
  ## below newrat, and (b) cost 30-50x the wall time of a newrat chain (a single
  ## cmaes_jade chain ran 36 min on NZSIM vs 40 s for newrat) -- so as an
  ## automatic default they only inflate wall time without ever winning the
  ## keep-best. They remain available as an OPT-IN diverse portfolio: pass a
  ## VECTOR `method` (e.g. c("newrat","newrat","cmaes_newrat","cmaes")) and the
  ## entry for each chain is honoured verbatim (recycled to n_chains).
  chain_methods <- if (length(method) > 1L) {
    rep_len(method, n_chains)
  } else {
    rep_len(method, n_chains)                       # pure multi-start
  }

  ## Per-method iteration budget. Pure `newrat` (the csminwel trust-region
  ## method) has a gradient-norm convergence test and stops early -- a big
  ## nm_maxit (e.g. 10000 "run to convergence") costs only the ~200 iters it
  ## actually needs. The global / population searchers (cmaes*, jade, nmkb,
  ## nelder, combined) have NO such early stop: they burn every iteration, and a
  ## cmaes generation is nm_maxit * lambda objective evals -- at nm_maxit=10000
  ## that is ~140k evals (hours), which would dominate the whole portfolio's wall
  ## time. Cap them at `global_maxit` so each global chain costs no more than a
  ## converged newrat chain; chain 1 (pure newrat) keeps the full budget so it
  ## stays the exact serial replica.
  global_maxit <- min(nm_maxit, 500L)
  configs <- lapply(seq_len(n_chains), function(ch) {
    meth <- chain_methods[[ch]]
    list(
      chain_id = ch, theta0 = starts[[ch]],
      seed = seed_base + ch * 1000L,
      nm_maxit = if (identical(meth, "newrat")) nm_maxit else global_maxit,
      method = meth
    )
  })

  cat(sprintf("  Parallel mode (mirai): %d chains on %d daemons | methods: %s | maxit: newrat=%d global=%d\n",
              n_chains, n_cores,
              paste(sprintf("C%d=%s", seq_len(n_chains), chain_methods), collapse = " "),
              nm_maxit, global_maxit))

  t_init <- proc.time()
  sh <- NULL
  if (isTRUE(pool_ready)) {
    ## Pool already provisioned by the CALLER (model compiled on each daemon);
    ## we only re-bind .worker_lp to this stage's lik_init and leave teardown to
    ## the caller. Saves a full per-daemon model recompile vs a fresh pool.
    .mirai_rebind_worker_lp(prior_spec, obs_names, me_variance = me_variance,
                            me_extra = me_extra, shock_scale = shock_scale,
                            system_priors = system_priors, lik_init = lik_init)
    cat(sprintf("  Daemon pool reused (lik_init = %s): %.1f sec\n",
                lik_init, (proc.time() - t_init)[["elapsed"]]))
  } else {
    if (!is.null(log_post_fn)) {
      .mirai_pool_closure(n_cores, log_post_fn)
    } else {
      sh <- .mirai_pool_init(n_cores, parsed_model, Y, prior_spec, obs_names,
                             me_variance, me_extra = me_extra,
                             shock_scale = shock_scale,
                             system_priors = system_priors, lik_init = lik_init)
    }
    on.exit({ mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }, add = TRUE)
    cat(sprintf("  Daemon init: %.1f sec\n", (proc.time() - t_init)[["elapsed"]]))
  }

  ## Live progress back-channel (host pull socket); daemons dial `prog_url`.
  prog_url <- NULL
  if (isTRUE(progress)) {
    pl <- .progress_listener(); prog_sock <- pl$sock; prog_url <- pl$url
    on.exit(try(nanonext::reap(prog_sock), silent = TRUE), add = TRUE)
  }

  t_global <- Sys.time()
  mode_task <- function(cfg) {
    log_post_fn <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    ## Match the main session's generator (daemons default to L'Ecuyer-CMRG)
    ## for reproducible, sequential-equivalent mode-finding.
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(cfg$seed)
    .psock <- NULL
    if (!is.null(prog_url)) {
      .psock <- nanonext::socket("push")
      nanonext::dial(.psock, url = prog_url, autostart = TRUE)
      on.exit(try(nanonext::reap(.psock), silent = TRUE), add = TRUE)
    }

    ## The optimizers (cmaes::cma_es, dfoptim::nmkb) are opaque blocking calls
    ## with no iteration hook, so we estimate progress from the OBJECTIVE
    ## EVALUATION count against the known evaluation budget. cmaes uses
    ## stopeval = max_iter * lambda; nmkb uses maxfeval = max_iter directly.
    ## We wrap the log-posterior to count calls (side-effect only -> no effect
    ## on the optimisation or RNG) and stream a fraction every 100 evals.
    n_dim  <- length(cfg$theta0)
    lambda <- 4L + floor(3L * log(n_dim))
    exp_evals <- switch(cfg$method,
      "cmaes"        = cfg$nm_maxit * lambda,
      "nmkb"         = cfg$nm_maxit,
      "nelder"       = cfg$nm_maxit,
      "cmaes_nmkb"   = round(0.8 * cfg$nm_maxit) * lambda + round(0.2 * cfg$nm_maxit),
      "combined"     = round(0.8 * cfg$nm_maxit) * lambda + round(0.2 * cfg$nm_maxit),
      "cmaes_jade"   = round(0.8 * cfg$nm_maxit) * lambda + round(0.2 * cfg$nm_maxit),
      ## csminwel: ~ (n_dim + 1) fcn evals per FD gradient + line search per iter
      "newrat"       = cfg$nm_maxit * (n_dim + 2L),
      "cmaes_newrat" = round(0.7 * cfg$nm_maxit) * lambda +
                       round(0.3 * cfg$nm_maxit) * (n_dim + 2L),
      cfg$nm_maxit * lambda)
    exp_evals <- max(1L, as.integer(exp_evals))
    .ec <- 0L
    lp_counted <- function(theta) {
      .ec <<- .ec + 1L
      if (!is.null(.psock) && .ec %% 100L == 0L)
        nanonext::send(.psock, c(cfg$chain_id, min(0.99, .ec / exp_evals)),
                       mode = "serial", block = 200)
      log_post_fn(theta)
    }

    ## Analytic gradient for the csminwel trust-region methods (newrat /
    ## cmaes_newrat). Built per chain from the daemon globals bound by
    ## .mirai_pool_init (.worker_model, .worker_cm, .worker_Y); without it
    ## csminwel falls back to a 68-eval finite-difference gradient that is far
    ## too slow to converge in parallel (the serial newrat reaches the mode
    ## precisely BECAUSE it has the analytic gradient). Composed with the param
    ## transform (chain rule) when transforming. NULL for non-gradient methods or
    ## when the pool was provisioned from a bare log_post_fn closure.
    grad_fn <- NULL
    if (isTRUE(analytic_grad) &&
        cfg$method %in% c("newrat", "cmaes_newrat", "lbfgsb", "combined")) {
      .wm <- get0(".worker_model", envir = globalenv(), inherits = FALSE)
      .wcm <- get0(".worker_cm",   envir = globalenv(), inherits = FALSE)
      .wY  <- get0(".worker_Y",    envir = globalenv(), inherits = FALSE)
      if (!is.null(.wm) && !is.null(.wcm) && !is.null(.wY)) {
        grad_fn <- tryCatch({
          .mkg <- utils::getFromNamespace("make_posterior_grad", "dynhr")
          ## THETA-space gradient, ALWAYS. .run_mode_finding (the consumer for
          ## these methods) expects grad_fn(theta) = d logpost / d theta and
          ## applies the eta chain rule ITSELF when `transform` is non-NULL
          ## (mode-orchestrate.R: -dtheta_deta(eta) * grad_fn(theta)). Wrapping
          ## with make_transformed_grad here would chain-rule TWICE (and feed an
          ## eta-expecting closure a theta argument) -> the gradient points the
          ## wrong way and newrat never leaves theta_init (logpost frozen at the
          ## prior-mean value). This was the real cause of the eta-space daemon
          ## "stuck" bug -- NOT a host/daemon or stale-install issue. The host
          ## full-API path (run-mode-finding.R) likewise passes the untransformed
          ## make_posterior_grad result; mirror it exactly.
          .mkg(.wm, .wY, prior_spec, obs_names, .wcm,
               me_variance = me_variance, me_extra = me_extra,
               shock_scale = shock_scale, grad_method = grad_method,
               likelihood = likelihood, freq_band = freq_band)
        }, error = function(e) NULL)
      }
    }

    ## Shared H0 seed (a precomputed neg-logpost Hessian matrix from the host):
    ## seed csminwel's initial curvature for the trust-region methods. Chain 1
    ## (theta0 == theta_init) gets the EXACT serial seed -> bit-identical to the
    ## serial newrat; dispersed chains reuse the same matrix as an approximate
    ## H0 (csminwel BFGS-refines it). NULL for global / gradient-free chains.
    hessian_fn <- NULL
    if (!is.null(newrat_H0_seed) &&
        cfg$method %in% c("newrat", "cmaes_newrat"))
      hessian_fn <- function(theta) newrat_H0_seed

    t0  <- proc.time()[["elapsed"]]
    run_mode <- utils::getFromNamespace(".run_mode_finding", "dynhr")
    ## tryCatch so a single bad start can't take down the whole chain.
    res <- tryCatch(
      run_mode(lp_counted, cfg$theta0, prior_spec,
               nm_maxit = cfg$nm_maxit, method = cfg$method,
               transform = transform, grad_fn = grad_fn,
               hessian_fn = hessian_fn, verbose = FALSE),
      error = function(e)
        list(theta_mode = cfg$theta0, logpost = -Inf, convergence = 1L,
             iterations = 0L, error = conditionMessage(e)))
    if (is.null(res) || is.null(res$theta_mode) ||
        !is.finite(res$logpost %||% NA_real_)) {
      res <- list(theta_mode = cfg$theta0, logpost = -Inf,
                  convergence = 1L, iterations = 0L,
                  error = res$error %||% "mode finding returned invalid result")
    }
    ## Completion ping (frac = 1) only on a genuine success, so a failed chain
    ## honestly shows < 100% rather than a full bar.
    ok <- is.finite(res$logpost) && is.null(res$error)
    if (!is.null(.psock) && ok)
      nanonext::send(.psock, c(cfg$chain_id, 1), mode = "serial", block = 500)
    list(chain_id = cfg$chain_id, theta_start = cfg$theta0, result = res,
         method = cfg$method, n_eval = .ec, nm_maxit = cfg$nm_maxit,
         elapsed_min = (proc.time()[["elapsed"]] - t0) / 60)
  }

  ## prior_spec / prog_url / transform are FREE variables in mode_task, so
  ## pass via `...` (see the note in run_mcmc_mirai on why `.args` would
  ## error here). `transform` (when non-NULL) is a plain S3 list of closures
  ## over local data and serializes to the daemons like any other free
  ## variable.
  ## Sever mode_task's env so mirai_map does not serialise run_mode_mirai's frame
  ## (parsed_model, Y, ...) with every chain. Bind the handful of light inputs the
  ## task needs into a fresh env whose parent is the dynhr namespace (for %||%
  ## etc.). NOTE: binding the DATA is essential -- a bare asNamespace() parent
  ## would resolve names like `prior_spec` / `transform` to the dynhr/base
  ## FUNCTIONS of the same name instead of the data ("closure not subsettable").
  ## .worker_* are fetched from the daemon globalenv via get0() inside the task.
  .mt_env <- new.env(parent = asNamespace("dynhr"))
  list2env(list(prior_spec = prior_spec, prog_url = prog_url,
                transform = transform, analytic_grad = analytic_grad,
                grad_method = grad_method, obs_names = obs_names,
                me_variance = me_variance, me_extra = me_extra,
                shock_scale = shock_scale, likelihood = likelihood,
                freq_band = freq_band, newrat_H0_seed = newrat_H0_seed),
           envir = .mt_env)
  environment(mode_task) <- .mt_env
  m_handle <- mirai::mirai_map(configs, mode_task)
  raw <- if (isTRUE(progress))
    .render_mode_progress(m_handle, prog_sock, n_chains)
  else m_handle[]

  wall_min <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))
  cat(sprintf("  All chains complete. Wall time: %.1f min\n", wall_min))

  ## Surface failed chains instead of silently collapsing them to -Inf (a
  ## swallowed error here propagates as an empty theta_mode -> 0x0 Sigma_prop
  ## downstream). Mirrors run_mcmc_mirai's per-chain failure reporting.
  .failed <- function(r) inherits(r, "miraiError") || inherits(r, "errorValue")
  for (i in seq_along(raw)) {
    r <- raw[[i]]
    if (.failed(r))
      cat(sprintf("  Mode chain %d FAILED: %s\n", i, as.character(r)))
    else if (!is.null(r$result$error))
      cat(sprintf("  Mode chain %d recovered from error: %s\n", i, r$result$error))
  }

  logposts <- vapply(raw, function(r) {
    if (.failed(r)) return(-Inf)
    lp <- r$result$logpost
    if (is.null(lp) || !is.finite(lp)) -Inf else as.numeric(lp)
  }, numeric(1))
  best_idx <- which.max(logposts)

  best_res <- if (.failed(raw[[best_idx]])) NULL else raw[[best_idx]]$result
  if (is.null(best_res) || is.null(best_res$theta_mode)) {
    warning("run_mode_mirai: all chains failed or returned no mode; ",
            "falling back to theta_init. See the per-chain errors above.",
            call. = FALSE)
    best_res <- list(theta_mode = theta_init, logpost = logposts[best_idx],
                     convergence = 1L, iterations = 0L,
                     error = "all mode chains failed")
  }
  cat(sprintf("  Best: chain %d  logpost=%.4f\n", best_idx, logposts[best_idx]))

  ## Per-chain diagnostics: method, wall time, objective-eval count, iterations,
  ## final quality. n_eval counts log-posterior (objective) calls only; for the
  ## gradient methods the analytic-gradient evaluations are additional.
  diagnostics <- do.call(rbind, lapply(seq_along(raw), function(i) {
    r <- raw[[i]]
    if (.failed(r))
      return(data.frame(chain = i, method = NA_character_, elapsed_min = NA_real_,
                        n_eval = NA_integer_, iterations = NA_integer_,
                        logpost = -Inf, converged = FALSE, stringsAsFactors = FALSE))
    data.frame(chain = i, method = r$method %||% NA_character_,
               elapsed_min = r$elapsed_min %||% NA_real_,
               n_eval = r$n_eval %||% NA_integer_,
               iterations = r$result$iterations %||% NA_integer_,
               logpost = logposts[i],
               converged = isTRUE((r$result$convergence %||% 1L) == 0L),
               stringsAsFactors = FALSE)
  }))

  list(
    chains     = lapply(raw, function(r)
                   if (.failed(r)) list(error = as.character(r)) else r$result),
    diagnostics = diagnostics,
    starts     = lapply(raw, function(r)
                   if (.failed(r)) NA else r$theta_start),
    logposts   = logposts,
    best       = best_res,
    best_chain = best_idx,
    wall_time  = wall_min,
    n_cores    = n_cores
  )
}


#' Provision a daemon pool that ships a log-posterior CLOSURE.
#'
#' For SMC, only the log-posterior closure (and prior sampler) are available, not
#' the parsed model, so we cannot compile per daemon. Instead we ship the closure
#' once via everywhere() (persistent daemons -> serialised once, not per task) as
#' \code{.worker_lp}, and optionally the prior sampler as \code{.worker_ps}. The
#' caller tears the pool down via \code{mirai::daemons(NULL)}.
#'
#' @param n_cores number of daemons.
#' @param log_post_fn the log-posterior closure.
#' @param prior_sampler optional prior sampler closure.
#' @return invisibly NULL.
#' @noRd
.mirai_pool_closure <- function(n_cores, log_post_fn, prior_sampler = NULL) {
  .restore_blas <- .mirai_pin_blas_threads()
  on.exit(.restore_blas(), add = TRUE)
  mirai::daemons(n_cores)
  mirai::everywhere(
    {
      suppressMessages(library(dynhr))
      ## `<<-` so the bindings reach the daemon globalenv (see the note in
      ## .mirai_pool_init); tasks retrieve them via get0(envir = globalenv()).
      .worker_lp <<- log_post_fn
      .worker_ps <<- prior_sampler
    },
    .args = list(log_post_fn = log_post_fn, prior_sampler = prior_sampler)
  )[]  # collect: block until every daemon has the closure (see .mirai_pool_init)
  invisible(NULL)
}


#' Start a bare daemon pool for SMC particle work.
#'
#' SMC's per-stage task closures (\code{.eval_particle}, \code{.mutate_one})
#' already carry everything they need (including the log-posterior) by lexical
#' capture, so the pool only has to exist with dynhr loaded on each daemon.
#' Loading dynhr once per daemon means the package namespace -- which the captured
#' closures reach into -- is present without re-loading per task. Caller tears
#' down via \code{mirai::daemons(NULL)}.
#'
#' @param n_cores number of daemons.
#' @return invisibly NULL.
#' @noRd
.smc_pool_setup <- function(n_cores) {
  mirai::daemons(n_cores)
  mirai::everywhere({ suppressMessages(library(dynhr)) })[]  # collect: see .mirai_pool_init
  invisible(NULL)
}


#' Parallel map for SMC, with reproducible per-task seeding.
#'
#' Replaces \code{future.apply::future_lapply(x, fn, future.seed = TRUE)} on a
#' mirai pool. Each task seeds deterministically with \code{seed_base + index}
#' before calling \code{fn}, so a run is reproducible. (This is a DIFFERENT RNG
#' stream than the old future L'Ecuyer one, so SMC results are reproducible but
#' not bit-identical to the previous future backend -- expected and documented.)
#' The closure \code{fn} is shipped once for the whole map, not per task.
#'
#' @param x vector mapped over (particle indices).
#' @param fn task closure of one index.
#' @param seed_base base RNG seed.
#' @return a list of results, one per element of \code{x}.
#' @noRd
.smc_pmap <- function(x, fn, seed_base = 0L) {
  task <- function(.i, .fn, .sb) {
    ## Match the main session's RNG (daemons default to L'Ecuyer-CMRG).
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(.sb + .i)
    .fn(.i)
  }
  ## Sever the env so mirai_map does not serialise .smc_pmap's frame (which
  ## holds `fn`) with every task on top of the once-per-task `.args` copy.
  environment(task) <- asNamespace("dynhr")
  mirai::mirai_map(x, task, .args = list(.fn = fn, .sb = seed_base))[]
}
