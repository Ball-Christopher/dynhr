## R/estimate-runner.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## .build_log_posterior()  -- thin wrapper / legacy fallback for make_log_posterior
## .run_mcmc()             -- sequential or parallel MCMC dispatcher
## .compute_convergence()  -- R-hat and ESS from a list of rwmh results
## estimate_model()        -- end-to-end DSGE estimation orchestrator
## .posterior_predictive() -- posterior predictive SD check
## --------------------------------------------------------------------------


## Thin compatibility wrapper around make_log_posterior().
## Prefer calling make_log_posterior() directly; this exists only so that
## legacy scripts that call .build_log_posterior() keep working.
.build_log_posterior <- function(model, data, prior_spec, obs_vars, compiled) {
  if (exists("make_log_posterior", mode = "function")) {
    lp <- make_log_posterior(model, data, prior_spec, obs_vars, compiled)
    return(function(theta) lp(theta))
  }
  if (exists("log_posterior_old", mode = "function")) {
    return(function(theta)
      log_posterior_old(theta, model, data, prior_spec, obs_vars, compiled))
  }
  stop("No log-posterior function found. Source dynhr_estimation.R first.")
}


## Sequential or parallel MCMC dispatcher.
## Delegates to run_mcmc_parallel() when config$parallel = TRUE and that
## function is available; otherwise runs rwmh() chains sequentially.
.run_mcmc <- function(log_post_fn, theta_mode, prior_spec, Sigma_prop,
                      config, parsed_model = NULL, Y = NULL,
                      obs_names = NULL, compiled = NULL) {
  n_chains <- config$n_chains
  n_draws  <- config$n_draws
  n_burn   <- config$n_burn

  if (config$parallel) {
    backend <- config$parallel_backend %||% "mirai"

    if (backend == "mirai") {
      cat("  Running parallel MCMC (mirai backend)...\n")
      return(run_mcmc_mirai(
        parsed_model = parsed_model, Y = Y,
        prior_spec = prior_spec, obs_names = obs_names,
        theta_mode = theta_mode, Sigma_prop = Sigma_prop,
        n_chains = n_chains, n_draws = n_draws, n_burn = n_burn,
        mh_scale      = config$mh_scale,
        target_accept = config$target_accept,
        adapt_every   = config$adapt_every,
        seed_base     = config$seed_base,
        n_cores       = config$n_cores
      ))
    }

    if (backend == "psock" && exists("run_mcmc_parallel", mode = "function")) {
      cat("  Running parallel MCMC (psock backend)...\n")
      return(run_mcmc_parallel(
        parsed_model = parsed_model, Y = Y,
        prior_spec = prior_spec, obs_names = obs_names,
        theta_mode = theta_mode, Sigma_prop = Sigma_prop,
        n_chains = n_chains, n_draws = n_draws, n_burn = n_burn,
        mh_scale      = config$mh_scale,
        target_accept = config$target_accept,
        adapt_every   = config$adapt_every,
        seed_base     = config$seed_base,
        dynhr_files   = config$dynhr_files
      ))
    }
  }

  cat(sprintf("  Running sequential MCMC: %d chains x %dk draws...\n",
              n_chains, n_draws / 1000))

  L <- t(chol(Sigma_prop))

  chains <- vector("list", n_chains)
  chain_stats <- data.frame(
    chain = integer(), accept_rate = numeric(),
    final_logpost = numeric(), final_scale = numeric(),
    elapsed_min = numeric(), stringsAsFactors = FALSE
  )

  for (ch in seq_len(n_chains)) {
    set.seed(config$seed_base + ch)

    if (ch == 1) {
      th0 <- theta_mode
    } else {
      z   <- rnorm(length(theta_mode))
      th0 <- theta_mode + 0.5 * as.numeric(L %*% z)
      names(th0) <- names(theta_mode)
      for (i in seq_along(th0)) {
        th0[i] <- max(th0[i], prior_spec$lower[i] + 1e-6)
        th0[i] <- min(th0[i], prior_spec$upper[i] - 1e-6)
      }
    }

    lp0 <- log_post_fn(th0)$logpost
    if (!is.finite(lp0)) th0 <- theta_mode

    cat(sprintf("  Chain %d/%d starting...\n", ch, n_chains))
    t0 <- Sys.time()
    chains[[ch]] <- rwmh(log_post_fn, th0, Sigma_prop,
                         n_draws     = n_draws + n_burn,
                         n_burn      = n_burn,
                         scale       = config$mh_scale,
                         target_rate = config$target_accept,
                         adapt_every = config$adapt_every,
                         verbose     = TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    if (!is.null(chains[[ch]])) {
      chain_stats <- rbind(chain_stats, data.frame(
        chain         = ch,
        accept_rate   = chains[[ch]]$acceptance_rate,
        final_logpost = tail(chains[[ch]]$post_logpost, 1),
        final_scale   = chains[[ch]]$scale,
        elapsed_min   = elapsed,
        stringsAsFactors = FALSE
      ))
      cat(sprintf("  Chain %d done: accept=%.1f%%, %.1f min\n",
                  ch, chains[[ch]]$acceptance_rate * 100, elapsed))
    }
  }

  list(chains = chains, chain_stats = chain_stats)
}


#' End-to-end DSGE estimation pipeline
#'
#' Runs the full pipeline: parse -> prior extraction -> compile -> steady state
#' -> mode-finding -> MCMC -> convergence diagnostics -> summary.
#'
#' @param config Named list with at least:
#'   \item{model_file}{Path to .mod file}
#'   \item{data_file}{Path to CSV data file}
#'   \item{obs_names}{Character vector of observed variable names}
#'   Optional keys: n_draws, n_burn, n_chains, mh_scale, target_accept,
#'   adapt_every, seed_base, nm_maxit, parallel, dynhr_files, output_dir,
#'   output_prefix, verbose, likelihood ("gaussian", "cumulant" or
#'   "whittle"), freq_band (radians, Whittle only).
#' @return Invisible list: config, parsed_model, compiled, prior_spec, mode,
#'   stoch_cal, stoch_mode, mcmc, convergence, data.
#' @noRd
estimate_model <- function(config) {

  defaults <- list(
    data_col_map  = NULL,
    n_draws       = 0L,
    n_burn        = 50000L,
    n_chains      = 4L,
    mh_scale      = 1.50,
    target_accept = 0.25,
    adapt_every   = 200L,
    seed_base     = 42L,
    nm_maxit      = 3500L,
    lbfgsb_maxit  = 200L,
    likelihood    = "gaussian",
    freq_band     = c(0, pi),
    parallel      = FALSE,
    parallel_backend = "mirai",
    n_cores       = NULL,
    dynhr_files   = NULL,
    output_dir    = ".",
    output_prefix = "est",
    verbose       = FALSE
  )
  for (nm in names(defaults))
    if (is.null(config[[nm]])) config[[nm]] <- defaults[[nm]]

  stopifnot(!is.null(config$model_file), !is.null(config$obs_names))

  vb     <- config$verbose
  prefix <- file.path(config$output_dir, config$output_prefix)

  cat("\n================================================================\n")
  cat("  estimate_model: end-to-end DSGE estimation\n")
  cat("================================================================\n\n")

  ## Step 1: Parse
  if (vb) cat("-- Step 1: Parse model --\n")
  parsed_model <- parse_mod(config$model_file)
  if (vb) cat(sprintf("  %d endo, %d exo, %d params, %d eqs\n",
                      length(parsed_model$var_names),
                      length(parsed_model$varexo_names),
                      length(parsed_model$param_values),
                      length(parsed_model$equations)))

  ## Step 2: Extract priors
  if (vb) cat("\n-- Step 2: Extract prior specification --\n")
  prior_spec <- extract_prior_spec(parsed_model, verbose = vb)

  ## Step 3: Compile
  if (vb) cat("-- Step 3: Compile model --\n")
  compiled <- compile_model(parsed_model, verbose = FALSE)

  ## Step 4: Solve at calibration
  if (vb) cat("-- Step 4: Solve at calibration --\n")
  stoch_cal <- stoch_simul(parsed_model, verbose = FALSE)
  bk_ok <- if (!is.null(stoch_cal$dr$bk_satisfied)) stoch_cal$dr$bk_satisfied
            else if (!is.null(stoch_cal$bk_satisfied)) stoch_cal$bk_satisfied
            else NA
  if (vb) cat(sprintf("  BK satisfied: %s\n", bk_ok))
  if (isFALSE(bk_ok)) stop("BK conditions not satisfied at calibration.")

  ## Step 5: Load data
  if (vb) cat("\n-- Step 5: Load data --\n")
  obs_names <- config$obs_names
  data_raw  <- read.csv(config$data_file)
  if (!is.null(config$data_col_map)) {
    for (mod_nm in names(config$data_col_map)) {
      dat_nm <- config$data_col_map[mod_nm]
      if (dat_nm %in% names(data_raw) && !(mod_nm %in% names(data_raw)))
        names(data_raw)[names(data_raw) == dat_nm] <- mod_nm
    }
  }
  stopifnot(all(obs_names %in% names(data_raw)))
  Y <- as.matrix(data_raw[, obs_names])
  if (vb) cat(sprintf("  Y: %d x %d\n", nrow(Y), ncol(Y)))

  ## Apply call-level filter_tunes override (config$filter_tunes).
  parsed_model <- .resolve_filter_tunes(parsed_model, config$filter_tunes)

  ## Expand observables for filter_tunes (no-op when no tunes are present).
  tunes_exp  <- .expand_observables_for_tunes(parsed_model, obs_names, Y)
  obs_names  <- tunes_exp$obs_vars
  Y          <- tunes_exp$Y
  me_extra   <- tunes_exp$me_extra

  ## Build shock_scale matrix from heteroskedastic_shocks block (Bug fix: was missing).
  shock_scale_mat <- .build_shock_scale_matrix(
    parsed_model, parsed_model$varexo_names, nrow(Y))

  ## Step 6: Build log-posterior
  if (vb) cat("\n-- Step 6: Build log-posterior --\n")
  log_post_fn <- make_log_posterior(parsed_model, Y, prior_spec,
                                    obs_names, compiled,
                                    me_extra      = me_extra,
                                    shock_scale   = shock_scale_mat,
                                    likelihood    = config$likelihood,
                                    freq_band     = config$freq_band,
                                    system_priors = config$system_priors %||% NULL)

  ## Step 7: Test at calibration
  if (vb) cat("\n-- Step 7: Test KF at calibration --\n")
  theta_cal <- setNames(numeric(nrow(prior_spec)), prior_spec$name)
  for (i in seq_len(nrow(prior_spec))) {
    nm <- prior_spec$name[i]
    theta_cal[i] <- if (nm %in% names(parsed_model$param_values))
                      as.numeric(parsed_model$param_values[[nm]])
                    else prior_spec$mean[i]
  }
  lp_cal <- log_post_fn(theta_cal)
  if (vb) cat(sprintf("  Log-posterior at calibration: %.4f\n", lp_cal$logpost))

  ## Step 9: Mode-finding
  if (vb) cat("\n-- Step 9: Mode-finding --\n")
  mode_res <- .run_mode_finding(log_post_fn, theta_cal, prior_spec,
                                nm_maxit     = config$nm_maxit,
                                lbfgsb_maxit = config$lbfgsb_maxit,
                                verbose      = vb)
  if (is.null(mode_res)) stop("Mode-finding failed.")

  theta_mode <- mode_res$theta_mode
  lp_mode    <- log_prior(theta_mode, prior_spec)
  ll_mode    <- mode_res$logpost - lp_mode

  mode_out <- list(
    theta_mode = theta_mode, logpost = mode_res$logpost,
    loglik = ll_mode, logprior = lp_mode,
    hessian = mode_res$hessian, V_mode = mode_res$V_mode,
    se_mode = mode_res$se_mode, prior_spec = prior_spec,
    obs_names = obs_names, method = "NM+L-BFGS-B"
  )
  saveRDS(mode_out, paste0(prefix, "_mode.rds"))
  if (vb) cat(sprintf("  Mode saved: %s_mode.rds\n", prefix))

  if (vb) {
    cat("\n  Mode estimates:\n")
    cat(sprintf("  %-14s %10s %10s %10s %10s\n",
                "Parameter", "Prior", "Mode", "SE", "Mode/Prior"))
    for (i in seq_along(theta_mode)) {
      se_i <- if (!is.null(mode_res$se_mode)) mode_res$se_mode[i] else NA
      cat(sprintf("  %-14s %10.4f %10.4f %10.4f %10.3f\n",
                  names(theta_mode)[i], prior_spec$mean[i],
                  theta_mode[i], se_i,
                  theta_mode[i] / prior_spec$mean[i]))
    }
    cat(sprintf("\n  loglik=%.2f  logprior=%.2f  logpost=%.2f\n",
                ll_mode, lp_mode, mode_res$logpost))
  }

  ## Step 10: Post-mode stoch_simul
  if (vb) cat("\n-- Step 10: Post-mode diagnostics --\n")
  orig_pv <- parsed_model$param_values
  for (nm in names(theta_mode))
    if (nm %in% names(parsed_model$param_values))
      parsed_model$param_values[[nm]] <- theta_mode[nm]
  stoch_mode <- stoch_simul(parsed_model, verbose = FALSE)
  parsed_model$param_values <- orig_pv

  ## Step 11-14: MCMC (optional)
  mcmc_res <- NULL
  conv_res <- NULL
  if (config$n_draws > 0) {
    if (vb) cat(sprintf("\n-- Step 11: MCMC (%d chains x %dk) --\n",
                        config$n_chains, config$n_draws / 1000))

    n_par      <- length(theta_mode)
    ## See .sampler_proposal_cov(): the V_mode branch below is unreachable on
    ## the .run_mode_finding() path, which is what made this proposal come
    ## from the PRIOR rather than the posterior Hessian.
    Sigma_prop <- if (!is.null(mode_res$V_mode)) {
                    S <- mode_res$V_mode * (2.38^2 / n_par)
                    dimnames(S) <- list(names(theta_mode), names(theta_mode))
                    S
                  } else {
                    .sampler_proposal_cov(log_post_fn, theta_mode, prior_spec,
                                          verbose = vb)
                  }

    mcmc_raw <- .run_mcmc(log_post_fn, theta_mode, prior_spec, Sigma_prop,
                          config, parsed_model, Y, obs_names, compiled)

    conv_res <- .compute_convergence(mcmc_raw$chains)
    combined <- conv_res$combined

    if (!is.null(conv_res$rhat) && vb) {
      cat(sprintf("  R-hat < 1.10: %d / %d\n",
                  sum(conv_res$rhat < 1.10), length(conv_res$rhat)))
      cat(sprintf("  Worst: %s = %.3f\n",
                  names(which.max(conv_res$rhat)), max(conv_res$rhat)))
      cat(sprintf("  Median ESS: %.0f\n", median(conv_res$ess)))
    }

    pp_sds <- .posterior_predictive(combined, parsed_model, obs_names, n_pp = 100)

    mcmc_res <- list(
      chains      = mcmc_raw$chains,
      combined    = combined,
      chain_stats = mcmc_raw$chain_stats,
      rhat        = conv_res$rhat,
      ess         = conv_res$ess,
      pp_sds      = pp_sds,
      config      = config,
      theta_mode  = theta_mode,
      V_mode      = mode_res$V_mode,
      prior_spec  = prior_spec,
      obs_names   = obs_names
    )
    saveRDS(mcmc_res, paste0(prefix, "_mcmc.rds"))
    if (vb) cat(sprintf("  MCMC saved: %s_mcmc.rds\n", prefix))
  } else {
    if (vb) cat("\n  MCMC skipped (n_draws = 0)\n")
  }

  ## Final summary
  cat("\n================================================================\n")
  cat("  ESTIMATION COMPLETE\n")
  cat("================================================================\n\n")
  cat(sprintf("  Model:     %s\n", config$model_file))
  cat(sprintf("  Data:      %s (%d x %d)\n", config$data_file, nrow(Y), ncol(Y)))
  cat(sprintf("  Params:    %d estimated\n", nrow(prior_spec)))
  cat(sprintf("  Mode:      logpost=%.2f\n", mode_res$logpost))
  if (!is.null(mcmc_res))
    cat(sprintf("  MCMC:      %d chains, accept=%s\n",
                config$n_chains,
                paste(sprintf("%.0f%%",
                              mcmc_res$chain_stats$accept_rate * 100),
                      collapse = ", ")))

  invisible(list(
    config       = config,
    parsed_model = parsed_model,
    compiled     = compiled,
    prior_spec   = prior_spec,
    mode         = mode_out,
    stoch_cal    = stoch_cal,
    stoch_mode   = stoch_mode,
    mcmc         = mcmc_res,
    convergence  = conv_res,
    data         = Y
  ))
}


## Posterior predictive SD check: simulates n_pp posterior draws and
## computes stoch_simul std_dev for each, for comparison with the data SDs.
.posterior_predictive <- function(combined, parsed_model, obs_names, n_pp = 100) {
  pp_sds <- matrix(NA_real_, n_pp, length(obs_names))
  colnames(pp_sds) <- obs_names
  orig_pv <- parsed_model$param_values
  idx     <- sample(nrow(combined), n_pp)

  for (k in seq_len(n_pp)) {
    theta_k <- combined[idx[k], ]
    for (nm in names(theta_k))
      if (nm %in% names(parsed_model$param_values))
        parsed_model$param_values[[nm]] <- theta_k[nm]

    r_k <- stoch_simul(parsed_model, verbose = FALSE)
    if (!is.null(r_k) && !is.null(r_k$moments$std_dev))
      for (v in obs_names)
        if (v %in% names(r_k$moments$std_dev))
          pp_sds[k, v] <- r_k$moments$std_dev[v]

    if (k %% 25 == 0) cat(sprintf("    %d/%d\n", k, n_pp))
  }
  parsed_model$param_values <- orig_pv
  pp_sds
}
