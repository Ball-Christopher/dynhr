## R/parallel-cluster.R
## --------------------------------------------------------------------------
## Phase-2 split from mcmc-parallel-monolith.R.
##
## run_mcmc_parallel()   -- parallel RWMH via PSOCK cluster + dynhr_files sourcing
## run_mode_parallel()   -- parallel multi-start mode-finding via future
## --------------------------------------------------------------------------


#' Replay the host's dynhr option state on every worker of a PSOCK cluster
#'
#' `.dynhr_opts` is a namespace-private ENVIRONMENT (and the `dynhr.*` switches
#' are base options), so a fresh PSOCK worker that merely `library(dynhr)`s
#' starts with an EMPTY option store: `power_posterior`, `debug_kf_errors`,
#' `dynhr.use_rcpp`, the hank backends and everything else set in the host
#' session were silently invisible to the workers, which then evaluated a
#' DIFFERENT posterior from the serial path. Only `me_variance` was ever shipped
#' (explicitly, as an argument).
#'
#' Call this ONCE per cluster, at init, AFTER `library(dynhr)` and after any
#' `dynhr_files` are sourced (a sourced copy of options.R installs its own empty
#' store, which would otherwise clobber the replayed values) -- never per task,
#' so the cost is O(n_workers).
#'
#' The `.dynhr_daemon_apply` FUNCTION OBJECT is shipped alongside the snapshot
#' rather than being called as `dynhr:::.dynhr_daemon_apply` on the worker: a
#' closure whose environment is a namespace serialises as a REFERENCE to that
#' namespace, so the body travels from the host while `.dynhr_opts` resolves in
#' the WORKER's own dynhr -- which is what must be written, and which also works
#' when the worker's installed dynhr predates the helper.
#'
#' @param cl A PSOCK cluster from `parallel::makeCluster()`.
#' @return (invisibly) the snapshot that was shipped.
#' @noRd
.dynhr_cluster_ship_options <- function(cl) {
  state    <- .dynhr_daemon_state()
  apply_fn <- .dynhr_daemon_apply
  parallel::clusterCall(cl, function(state, apply_fn) {
    apply_fn(state)
    ## A sourced dev copy of options.R installs its OWN `.dynhr_opts` store in
    ## globalenv, which shadows the namespace one for globalenv-resolved
    ## functions; keep both in sync.
    if (exists(".dynhr_daemon_apply", envir = globalenv(), inherits = FALSE))
      get(".dynhr_daemon_apply", envir = globalenv())(state)
    NULL
  }, state = state, apply_fn = apply_fn)
  invisible(state)
}


#' Run RWMH chains in parallel via a PSOCK cluster
#'
#' Workers re-source the package files and recompile the model (closures
#' cannot be serialised across PSOCK workers).  Pass all R source files the
#' workers need via `dynhr_files`.
#'
#' @param parsed_model   dynhr_mod (from parse_mod())
#' @param Y              Observation matrix (T ?-- n_obs)
#' @param prior_spec     Prior spec data.frame (from extract_prior_spec())
#' @param obs_names      Character vector of observed variable names
#' @param theta_mode     Named mode vector
#' @param Sigma_prop     Proposal covariance (n?--n)
#' @param n_chains       Number of chains
#' @param n_draws        Post-burn draws per chain
#' @param n_burn         Burn-in draws per chain
#' @param mh_scale       Initial RWMH scale
#' @param target_accept  Target acceptance rate for adaptive scaling
#' @param adapt_every    Adapt scale every N draws (during burn-in)
#' @param seed_base      Base RNG seed
#' @param n_cores        Worker count (NULL = n_chains or ncores-1)
#' @param dynhr_files    Character vector of R files to source on workers
#' @return list(chains, chain_stats, wall_time, n_cores)
#' @noRd
run_mcmc_parallel <- function(
    parsed_model, Y, prior_spec, obs_names,
    theta_mode, Sigma_prop,
    n_chains      = 4L,
    n_draws       = 100000L,
    n_burn        = 50000L,
    mh_scale      = 1.50,
    target_accept = 0.25,
    adapt_every   = 200L,
    seed_base     = 42L,
    n_cores       = NULL,
    dynhr_files   = NULL,
    me_variance   = NULL
) {
  me_variance <- .dynhr_opt("me_variance", me_variance, default = 0)
  if (is.null(dynhr_files) || length(dynhr_files) == 0)
    stop("dynhr_files must list all R source files workers need.")

  max_cores <- parallel::detectCores(logical = TRUE)
  if (is.null(n_cores)) {
    n_cores <- if (max_cores <= 4L) max_cores else max_cores - 2L
    n_cores <- max(1L, n_cores)
  }
  n_cores <- min(n_cores, n_chains)

  cat(sprintf("  Parallel MCMC: %d chains on %d workers (of %d logical processors detected)\n",
              n_chains, n_cores, max_cores))
  cat(sprintf("  Per chain: %dk draws + %dk burn-in\n",
              n_draws / 1000, n_burn / 1000))

  n_par <- length(theta_mode)
  L <- t(chol(Sigma_prop))

  starts <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    set.seed(seed_base + ch)
    if (ch == 1) {
      starts[[ch]] <- theta_mode
    } else {
      z  <- rnorm(n_par)
      th <- theta_mode + 0.5 * as.numeric(L %*% z)
      names(th) <- names(theta_mode)
      for (i in seq_len(n_par)) {
        th[i] <- max(th[i], prior_spec$lower[i] + 1e-6)
        th[i] <- min(th[i], prior_spec$upper[i] - 1e-6)
      }
      starts[[ch]] <- th
    }
  }

  total_draws <- n_draws + n_burn
  configs <- lapply(seq_len(n_chains), function(ch) list(
    chain_id      = ch,
    theta0        = starts[[ch]],
    seed          = seed_base + ch * 1000L,
    n_draws       = n_draws,
    n_burn        = n_burn,
    total_draws   = total_draws,
    scale         = mh_scale,
    target_accept = target_accept,
    adapt_every   = adapt_every
  ))

  t_init <- proc.time()
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  wd <- getwd()
  parallel::clusterCall(cl, setwd, wd)

  # Load the installed package first (provides compiled Rcpp functions),
  # then source dev R files to override with any updated pure-R functions.
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(dynhr)))
  for (f in dynhr_files)
    if (file.exists(f)) parallel::clusterCall(cl, source, f, local = FALSE)

  parallel::clusterExport(cl, c("parsed_model", "Y", "prior_spec", "obs_names",
                                "Sigma_prop"),
                          envir = environment())

  ## Ship the HOST's dynhr option state to every worker (see
  ## .dynhr_cluster_ship_options): without it the workers see an EMPTY option
  ## store and evaluate a DIFFERENT posterior from the serial path. AFTER the
  ## dynhr_files sourcing above, once per worker.
  .dynhr_cluster_ship_options(cl)

  parallel::clusterExport(cl, "me_variance", envir = environment())
  parallel::clusterEvalQ(cl, {
    .worker_cm <- compile_model(parsed_model, verbose = FALSE)
    .worker_lp <- make_log_posterior(parsed_model, Y, prior_spec,
                                     obs_names, .worker_cm,
                                     me_variance = me_variance)
  })

  dt_init <- (proc.time() - t_init)[["elapsed"]]
  cat(sprintf("  Cluster init: %.1f sec (source + compile + lp_fn)\n", dt_init))
  cat("  Chains running...\n")
  t_global <- Sys.time()

  raw <- parallel::parLapply(cl, configs, function(cfg) {
    set.seed(cfg$seed)
    t0 <- proc.time()[["elapsed"]]
    res <- rwmh(
      log_post_fn = .worker_lp,
      theta0      = cfg$theta0,
      Sigma_prop  = Sigma_prop,
      n_draws     = cfg$total_draws,
      n_burn      = cfg$n_burn,
      scale       = cfg$scale,
      target_rate = cfg$target_accept,
      adapt_every = cfg$adapt_every,
      verbose     = FALSE,
      progressor  = NULL,
      chain_id    = cfg$chain_id
    )
    if (is.null(res) || is.character(res) ||
        (!is.list(res$error) && !is.null(res$error))) {
      err_msg <- if (is.character(res)) res else "rwmh returned invalid result"
      res <- list(error = err_msg)
    }
    elapsed <- (proc.time()[["elapsed"]] - t0) / 60
    list(chain_id = cfg$chain_id, result = res, elapsed_min = elapsed)
  })

  wall_min <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))
  cat(sprintf("  All chains complete. Wall time: %.1f min\n", wall_min))

  chains <- vector("list", n_chains)
  chain_stats <- data.frame(
    chain = integer(), accept_rate = numeric(),
    final_logpost = numeric(), final_scale = numeric(),
    elapsed_min = numeric(), stringsAsFactors = FALSE
  )

  for (r in raw) {
    ch <- r$chain_id
    if (!is.null(r$result$error)) {
      cat(sprintf("  Chain %d FAILED: %s\n", ch, r$result$error))
      chain_stats <- rbind(chain_stats, data.frame(
        chain = ch, accept_rate = NA, final_logpost = NA,
        final_scale = NA, elapsed_min = r$elapsed_min,
        stringsAsFactors = FALSE
      ))
      next
    }
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

  for (i in seq_len(nrow(chain_stats))) {
    s <- chain_stats[i, ]
    if (is.na(s$accept_rate)) next
    cat(sprintf("  Chain %d: accept=%.1f%%, scale=%.3f, logpost=%.2f, %.1f min\n",
                s$chain, s$accept_rate * 100, s$final_scale,
                s$final_logpost, s$elapsed_min))
  }

  seq_min <- sum(chain_stats$elapsed_min, na.rm = TRUE)
  cat(sprintf("\n  Speedup: %.1fx (%.1f min wall vs %.1f min sequential)\n",
              seq_min / max(wall_min, 0.01), wall_min, seq_min))

  list(chains = chains, chain_stats = chain_stats,
       wall_time = wall_min, n_cores = n_cores)
}


#' Multi-start parallel mode-finding via future
#'
#' Runs `n_chains` independent optimisations; chain 1 starts at the exact
#' `theta_init`, chains 2..N start at perturbations scaled by `perturb_scale`
#' times the prior standard deviation.  Returns the best result.
#'
#' Requires the `future` and `future.apply` packages.
#'
#' @param parsed_model   dynhr_mod
#' @param Y              Observation matrix
#' @param prior_spec     Prior spec data.frame
#' @param obs_names      Observed variable names
#' @param theta_init     Named starting vector
#' @param n_chains       Number of chains (NULL = n_cores)
#' @param nm_maxit       Mode-finding iteration budget per chain
#' @param method         Optimizer name (passed to .run_mode_finding())
#' @param perturb_scale  Perturbation as fraction of prior std
#' @param seed_base      Base RNG seed
#' @param n_cores        Worker count (NULL = ncores - 1, capped at 8)
#' @param dynhr_files    R files to source on workers
#' @return list(chains, starts, logposts, best, best_chain, wall_time, n_cores)
#' @noRd
run_mode_parallel <- function(
    parsed_model, Y, prior_spec, obs_names,
    theta_init,
    n_chains       = NULL,
    nm_maxit       = NULL,
    method         = NULL,
    perturb_scale  = NULL,
    seed_base      = NULL,
    n_cores        = NULL,
    dynhr_files    = NULL,
    me_variance    = NULL
) {
  me_variance   <- .dynhr_opt("me_variance",   me_variance,   default = 0)
  nm_maxit      <- .dynhr_opt("nm_maxit",      nm_maxit,      default = 5000L)
  method        <- .dynhr_opt("mode_method",   method,        default = "cmaes_nmkb")
  perturb_scale <- .dynhr_opt("perturb_scale", perturb_scale, default = 0.5)
  seed_base     <- .dynhr_opt("seed_base",     seed_base,     default = 42L)
  if (is.null(dynhr_files) || length(dynhr_files) == 0)
    stop("dynhr_files must list all R source files workers need.")

  for (.p in c("future", "future.apply")) {
    if (!requireNamespace(.p, quietly = TRUE))
      stop(sprintf("Package '%s' required. Install with install.packages('%s')", .p, .p),
           call. = FALSE)
  }

  max_cores <- parallel::detectCores(logical = TRUE)
  if (is.null(n_cores)) {
    n_cores <- if (max_cores <= 4L) max_cores else max_cores - 2L
    n_cores <- max(1L, n_cores)
  }
  if (is.null(n_chains)) n_chains <- n_cores
  n_cores <- min(n_cores, n_chains)

  n_par     <- length(theta_init)
  par_names <- names(theta_init)
  prior_sds <- setNames(prior_spec$std,   prior_spec$name)
  prior_lo  <- setNames(prior_spec$lower, prior_spec$name)
  prior_hi  <- setNames(prior_spec$upper, prior_spec$name)

  starts <- vector("list", n_chains)
  starts[[1L]] <- theta_init

  for (ch in seq(2L, n_chains)) {
    set.seed(seed_base + ch)
    th <- theta_init
    for (nm in par_names) {
      ps <- prior_sds[nm]
      if (!is.na(ps) && ps > 0) {
        th[nm] <- th[nm] + rnorm(1L, 0, perturb_scale * ps)
        lo <- prior_lo[nm]; hi <- prior_hi[nm]
        if (!is.na(lo) && is.finite(lo)) th[nm] <- max(th[nm], lo + 1e-6)
        if (!is.na(hi) && is.finite(hi)) th[nm] <- min(th[nm], hi - 1e-6)
      }
    }
    starts[[ch]] <- th
  }

  configs <- lapply(seq_len(n_chains), function(ch) list(
    chain_id = ch,
    theta0   = starts[[ch]],
    seed     = seed_base + ch * 1000L,
    nm_maxit = nm_maxit,
    method   = method
  ))

  future::plan(future::multisession, workers = n_cores)
  on.exit(future::plan(future::sequential), add = TRUE)

  cat(sprintf("  Parallel mode: %d chains on %d cores | %s | maxit=%d/chain\n",
              n_chains, n_cores, method, nm_maxit))
  cat(sprintf("  Perturbation scale: %.2f x prior std  (chain 1 = exact calibration)\n",
              perturb_scale))

  t_global <- Sys.time()

  ## Same option-shipping problem as the PSOCK path above (see
  ## .dynhr_cluster_ship_options): a multisession worker starts with an empty
  ## `.dynhr_opts` and none of the host's `dynhr.*` base options. Ship the
  ## snapshot AND the apply function object as future globals; the function's
  ## namespace environment re-resolves in the worker's own dynhr.
  .dynhr_worker_opt_state <- .dynhr_daemon_state()
  .dynhr_worker_opt_apply <- .dynhr_daemon_apply

  raw <- future.apply::future_lapply(configs, function(cfg) {
    suppressPackageStartupMessages(library(dynhr))
    for (f in dynhr_files) if (file.exists(f)) source(f, local = FALSE)
    .dynhr_worker_opt_apply(.dynhr_worker_opt_state)
    if (exists(".dynhr_daemon_apply", envir = globalenv(), inherits = FALSE))
      get(".dynhr_daemon_apply", envir = globalenv())(.dynhr_worker_opt_state)
    .cm  <- compile_model(parsed_model, verbose = FALSE)
    .lp  <- make_log_posterior(parsed_model, Y, prior_spec, obs_names, .cm,
                               me_variance = me_variance)
    set.seed(cfg$seed)
    t0  <- proc.time()["elapsed"]
    res <- .run_mode_finding(.lp, cfg$theta0, prior_spec,
                              nm_maxit = cfg$nm_maxit,
                              method   = cfg$method,
                              verbose  = FALSE)
    if (is.null(res) || is.null(res$theta_mode) ||
        !is.finite(res$logpost %||% NA_real_)) {
      res <- list(theta_mode = cfg$theta0, logpost = -Inf,
                  convergence = 1L, iterations = 0L,
                  error = "mode finding returned invalid result")
    }
    elapsed <- (proc.time()["elapsed"] - t0) / 60
    list(chain_id    = cfg$chain_id,
         theta_start = cfg$theta0,
         result      = res,
         elapsed_min = elapsed)
  }, future.seed = TRUE)

  wall_min <- as.numeric(difftime(Sys.time(), t_global, units = "mins"))
  cat(sprintf("  All chains complete. Wall time: %.1f min\n", wall_min))

  logposts <- vapply(raw, function(r) {
    lp <- r$result$logpost
    if (is.null(lp) || !is.finite(lp)) -Inf else as.numeric(lp)
  }, numeric(1))

  best_idx <- which.max(logposts)

  cat(sprintf("\n  %-7s  %12s  %8s  %s\n", "chain", "logpost", "min", "start_delta"))
  for (r in raw) {
    lp_str  <- if (is.finite(logposts[r$chain_id])) sprintf("%12.4f", logposts[r$chain_id])
               else sprintf("%12s", "-Inf")
    err_str <- if (!is.null(r$result$error)) paste("[ERR]", r$result$error) else ""
    delta   <- r$theta_start - theta_init
    d_str   <- if (r$chain_id == 1L) "(calibration)"
               else sprintf("max|delta|=%.3f (%s)", max(abs(delta)),
                            names(which.max(abs(delta))))
    best_mk <- if (r$chain_id == best_idx) " <- best" else ""
    cat(sprintf("  ch %-4d  %s  %8.1f  %s%s%s\n",
                r$chain_id, lp_str, r$elapsed_min, d_str, best_mk, err_str))
  }

  cat(sprintf("\n  Best: chain %d  logpost=%.4f\n",
              best_idx, logposts[best_idx]))

  list(
    chains     = lapply(raw, function(r) r$result),
    starts     = lapply(raw, function(r) r$theta_start),
    logposts   = logposts,
    best       = raw[[best_idx]]$result,
    best_chain = best_idx,
    wall_time  = wall_min,
    n_cores    = n_cores
  )
}
