## R/run-full-estimation.R
## --------------------------------------------------------------------------
## run_full_estimation() -- Phase 4 top-level orchestrator.
##
## Runs the full pipeline in one call:
##   parse -> compile -> prior -> mode -> sampler -> diagnostics -> report
##
## Designed so that a new model can be estimated in <30 lines of user code:
##
##   result <- run_full_estimation(
##     mod_file  = "mymodel.mod",
##     data      = Y,           # matrix or path to CSV
##     obs_vars  = c("y", "c"),
##     output_dir = "output/",
##     sampler   = "rwmh",
##     n_draws   = 100000L
##   )
##   print(result)
##   write_report(result, "output/report.md")
## --------------------------------------------------------------------------


# ============================================================================
# S3 class: dynhr_estimation_result
# ============================================================================

#' Print an estimation result
#'
#' @param x    A \code{dynhr_estimation_result} object
#' @param ...  Ignored
#' @return \code{x} invisibly
#' @export
print.dynhr_estimation_result <- function(x, ...) {
  cat("\n<dynhr_estimation_result>\n")
  cat(sprintf("  Model     : %s\n", x$meta$mod_file))
  cat(sprintf("  Data      : %d x %d  (obs: %s)\n",
              nrow(x$data), ncol(x$data),
              paste(colnames(x$data), collapse = ", ")))
  cat(sprintf("  Parameters: %d estimated\n", nrow(x$prior_spec)))
  if (isTRUE(x$meta$use_obc))
    cat(sprintf("  OBC       : %d constraints (PKF)\n",
                length(x$obc_specs %||% list())))
  cat(sprintf("  Mode      : logpost = %.3f\n", x$mode$logpost))
  if (!is.null(x$chains)) {
    n_draws <- nrow(x$chains$chain)
    sampler <- toupper(x$chains$sampler %||% "?")
    cat(sprintf("  MCMC      : [%s] %d draws\n", sampler, n_draws))
    if (!is.null(x$convergence)) {
      rh <- x$convergence$rhat
      if (!is.null(rh))
        cat(sprintf("  R-hat     : max = %.3f  (>1.05: %d / %d)\n",
                    max(rh), sum(rh > 1.05), length(rh)))
      es <- x$convergence$ess
      if (!is.null(es))
        cat(sprintf("  ESS       : median = %.0f  min = %.0f\n",
                    median(es), min(es)))
    }
    if (!is.null(x$chains$log_marginal_lik))
      cat(sprintf("  log p(Y|M): %.3f\n", x$chains$log_marginal_lik))
    ## THAMES cross-check (guarded: cannot break print)
    thames_line <- tryCatch({
      tr <- thames_mdd_from_chains(x$chains)
      if (is.finite(tr$log_mdd))
        sprintf("  log p(Y|M) [THAMES]: %.3f +/- %.3f\n", tr$log_mdd, tr$se)
      else
        NULL
    }, error = function(e) NULL)
    if (!is.null(thames_line)) cat(thames_line)
  } else {
    cat("  MCMC      : skipped (n_draws = 0)\n")
  }
  if (!is.null(x$diagnostics))
    cat(sprintf("  Diagnostics: %d run\n", length(x$diagnostics)))
  invisible(x)
}


# ============================================================================
# Helper: Compute smoother and historical decomposition at mode
# ============================================================================

.compute_smoother_at_mode <- function(theta_mode, model, compiled, data, obs_vars,
                                       obc_specs, me_variance, max_inner) {
  # Extract parameters (incl. estimated shock stds via .apply_theta_to_params,
  # so the at-mode smoother/report uses the estimated stds, not frozen ones).
  params <- .apply_theta_to_params(model, theta_mode)

  # Solve steady state at mode
  ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
  if (is.null(ss_result) || !ss_result$converged)
    return(NULL)

  # Cache system structure and extract matrices at mode
  sys_cache <- cache_system_structure(compiled)
  ## Re-derive SSM-computed params for a consistent linearization point
  ## (no-op for non-SSM-parameter models; Tier 13 #1).
  params <- ss_result$params %||% params
  sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)

  # Solve perturbation (slack regime)
  dr_slack <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
  if (is.null(dr_slack) || !dr_slack$bk_satisfied)
    return(NULL)

  # Set up regime cache and seed with slack policy
  regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
  obs_idx <- match(obs_vars, model$var_names)
  obc_ensure_policy(0L, regime_cache, sys, dr_slack, obc_specs, obs_idx)

  # Prepare data
  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  # Run PKF with storage for smoother
  kf_res <- kalman_filter_obc_pkf(
    Y, dr_slack, regime_cache, sys,
    model, params, obs_vars, obc_specs,
    obs_idx         = obs_idx,
    regime_path_init = NULL,
    me_variance     = me_variance,
    max_inner       = max_inner,
    return_filtered = FALSE,
    return_shocks   = FALSE,
    return_store    = TRUE
  )

  if (is.null(kf_res) || !is.finite(kf_res$loglik))
    return(NULL)

  # Run fixed-interval smoother
  smoother_res <- pkf_smoother_obc(kf_res$kf_store)

  # Compute historical decomposition
  hd_res <- historical_decomposition_obc(
    smoother_res$smoothed_shocks,
    kf_res$regime_path,
    regime_cache,
    dr_slack
  )

  # Package results
  list(
    smoothed_states  = smoother_res$smoothed_states,
    smoothed_shocks  = smoother_res$smoothed_shocks,
    hd               = hd_res,
    regime_path      = kf_res$regime_path,
    loglik_mode      = kf_res$loglik
  )
}


# ============================================================================
# Main orchestrator
# ============================================================================

#' Run a full DSGE estimation pipeline
#'
#' Parses a \code{.mod} file, compiles the model, extracts priors, locates the
#' posterior mode, draws from the posterior with the chosen sampler, and
#' optionally runs the diagnostic battery.  Results are saved to
#' \code{output_dir} and returned as a structured list.
#'
#' @param mod_file   Path to a Dynare \code{.mod} file.  Alternatively, pass
#'   an already-parsed model via \code{model}.
#' @param data       Observation matrix (\eqn{T \times n_{\text{obs}}}), column
#'   names matching \code{obs_vars}.  May also be a path to a CSV file.
#' @param obs_vars   Character vector: names of the observable variables in the
#'   model that correspond to columns of \code{data}.
#' @param output_dir Directory for saved outputs (\code{_mode.rds},
#'   \code{_mcmc.rds}).  Created if it does not exist.  Defaults to the
#'   current working directory.
#' @param output_prefix Prefix prepended to saved file names (default
#'   \code{"dynhr_est"}).
#' @param sampler    Which sampler to use: \code{"rwmh"} (default),
#'   \code{"smc"}, \code{"nuts"}, or \code{"dime"}.
#' @param n_draws    Post-warmup draws to retain (default \code{10000L}).
#'   Set to \code{0} to run mode-finding only.
#' @param n_warmup   Warmup / burn-in draws, discarded (default \code{5000L}).
#'   Ignored for \code{smc}.
#' @param n_chains   Number of independent MCMC chains (default \code{4L}).
#'   Ignored for \code{smc} and \code{dime} (use \code{n_particles} /
#'   \code{n_walkers}).
#' @param n_particles Particle count for SMC (default \code{2000L}).  Ignored
#'   for \code{rwmh}, \code{nuts}, and \code{dime}.
#' @param n_walkers  Number of ensemble walkers for DIME (default \code{NULL}
#'   = max(5 * n_par, 20)). Ignored for all other samplers.
#' @param parallel Run multi-chain RWMH/NUTS/SMC in parallel on a persistent
#'   \code{mirai} daemon pool (default \code{FALSE}). Standard Gaussian
#'   Kalman-likelihood models recompile the posterior once per daemon via
#'   \code{.mirai_pool_init()}; OBC/PKF models instead ship the already-built
#'   log-posterior closure once via \code{.mirai_pool_closure()}.
#' @param parallel_backend Parallel backend for RWMH chains: \code{"mirai"}
#'   (default).  Reserved for future backends.
#' @param n_cores Worker (daemon) count when \code{parallel = TRUE}
#'   (\code{NULL} = auto-detect, capped at \code{n_chains}).
#' @param analytic_grad For \code{sampler = "nuts"} on a standard Gaussian
#'   model, use the exact analytic gradient (\code{\link{make_posterior_grad}}:
#'   a compiled Kalman score for the shock-std parameters plus a numerical
#'   gradient for the rest) instead of a fully numerical gradient. Default
#'   \code{FALSE}. Applies to the serial (single-chain) NUTS path; the parallel
#'   multi-chain path uses the numerical gradient.
#' @param n_mode_iter Maximum optimizer iterations for mode-finding (default
#'   \code{10000L}).
#' @param mode_method Optimizer sequence for mode-finding (default
#'   \code{"newrat"}, the csminwel quasi-Newton optimizer; most robust and
#'   cheapest on near-unit-root models).  See \code{\link{find_mode}} for the
#'   full list of options (e.g. \code{"cmaes_newrat"}, \code{"combined"}).
#' @param mode_n_starts Number of dispersed starting points for parallel
#'   multi-start mode-finding (Step 6) when \code{parallel = TRUE}
#'   (\code{NULL} = one per daemon). Ignored when \code{parallel = FALSE}.
#' @param me_variance Measurement-error variance added to the observation
#'   noise diagonal (default \code{0}).  Use a small positive value for
#'   stochastically singular models.
#' @param seed       Random seed for reproducibility (default \code{42L}).
#' @param run_diag   Run the diagnostic battery after sampling?  (default
#'   \code{FALSE}; set \code{TRUE} for a full run).
#' @param verbose    Print progress messages (default \code{TRUE}).
#' @param model      Optional: a pre-parsed \code{dynhr_mod} object (from
#'   \code{\link{parse_mod}}).  If supplied, \code{mod_file} is ignored.
#' @param compiled   Optional: a pre-compiled model (from
#'   \code{\link{compile_model}}).  Skips compilation if supplied.
#' @param data_col_map Optional named character vector mapping model observable
#'   names to CSV column names when they differ.
#' @param dates      Optional: date vector of length \code{nrow(data)}, passed
#'   to diagnostics for time-axis labelling.
#' @param obc        OBC mode: \code{NULL} (default) auto-detects from MCP tags
#'   in the model, \code{TRUE} forces OBC, \code{FALSE} forces standard
#'   Kalman filter even if tags are present.
#' @param obc_max_inner Maximum inner PKF iterations per period (default
#'   \code{10L}).  Ignored unless OBC is active.
#' @param obc_ppf_reweight If \code{TRUE} and OBC is active, re-weights the
#'   PKF-sampled chains by PPF importance weights at the posterior draws
#'   (default \code{FALSE}).  Ignored unless OBC is active.
#' @param compute_smoother If \code{TRUE} and OBC is active, runs the
#'   fixed-interval PKF smoother and historical decomposition at the posterior
#'   mode after sampling.  Results are stored in \code{$smoother} (default
#'   \code{FALSE}).
#' @param run_ramsey If \code{TRUE}, runs the Ramsey policy workflow at the
#'   posterior mode (or calibrated parameters when \code{n_draws = 0}).
#' @param ramsey_order Perturbation order used by Ramsey workflow (1 or 2).
#' @param ramsey_n_periods Simulation length for welfare evaluation.
#' @param ramsey_burn_in Burn-in for Ramsey welfare simulation.
#' @param likelihood  Likelihood type: \code{"gaussian"} (Kalman filter, default),
#'   \code{"cumulant"}, \code{"whittle"}, or \code{"tpf"} (Tempered Particle
#'   Filter).
#' @param lik_init    Kalman filter \code{P0} initialisation
#'   (default \code{"auto"}).  Ignored for non-Gaussian likelihoods.
#' @param freq_band   Numeric(2) \code{c(lo, hi)} in radians; Whittle band
#'   restriction (default \code{c(0, pi)}, ignored for other likelihoods).
#' @param system_priors  Optional \code{system_prior_spec} object
#'   (from \code{\link{system_prior_spec}}) providing system-prior
#'   log-density contributions, or \code{NULL}.
#' @param heteroskedastic_shocks  Call-level heteroskedastic-shocks override.
#'   \code{NULL} (default) uses the \code{heteroskedastic_shocks} block from
#'   the \code{.mod} file if present.  Mutually exclusive with \code{plan}.
#' @param tpf_options  Named list of options forwarded to
#'   \code{make_log_posterior_tpf} when \code{likelihood = "tpf"}.
#'   Typical keys: \code{n_particles}, \code{ess_target}, \code{n_mh},
#'   \code{seed}.
#' @param plan  A \code{\link{dynhr_plan}} object bundling
#'   \code{filter_tunes} and \code{heteroskedastic_shocks} overrides into a
#'   single argument.  Mutually exclusive with \code{filter_tunes} and
#'   \code{heteroskedastic_shocks}.
#' @param filter_tunes Call-level tune override.  \code{NULL} (default) uses
#'   the \code{filter_tunes} block from the \code{.mod} file if present.  Pass
#'   a \code{filter_tunes_spec} object (built with
#'   \code{\link{filter_tunes}()}) to override the mod-file block entirely.
#'   Pass \code{FALSE} to ignore the mod-file block.
#' @param ...        Additional arguments forwarded to the chosen sampler.
#'
#' @return A \code{dynhr_estimation_result} list with elements:
#'   \describe{
#'     \item{\code{model}}{Parsed \code{dynhr_mod}}
#'     \item{\code{compiled}}{Compiled model}
#'     \item{\code{prior_spec}}{Prior specification data.frame}
#'     \item{\code{data}}{Observation matrix}
#'     \item{\code{mode}}{Mode-finding result list}
#'     \item{\code{chains}}{A \code{\link{dynhr_chains}} object, or \code{NULL} if \code{n_draws = 0}}
#'     \item{\code{convergence}}{List with \code{rhat}, \code{ess}, \code{combined} (multi-chain only)}
#'     \item{\code{diagnostics}}{Named list of \code{dynhr_diagnostic} objects, or \code{NULL}}
#'     \item{\code{obc_specs}}{OBC constraint specifications list (when OBC active), else \code{NULL}}
#'     \item{\code{smoother}}{List with \code{smoothed_states}, \code{smoothed_shocks},
#'       \code{hd} (historical decomposition), \code{regime_path} at posterior mode.
#'       \code{NULL} unless \code{compute_smoother = TRUE} and OBC is active.}
#'     \item{\code{ramsey}}{A \code{dynhr_ramsey_result}, or \code{NULL} if
#'       \code{run_ramsey = FALSE}.}
#'     \item{\code{meta}}{Run metadata: mod_file, sampler, seed, timestamps, use_obc}
#'   }
#'
#' @seealso \code{\link{parse_mod}}, \code{\link{find_mode}}, \code{\link{mcmc}},
#'   \code{\link{smc}}, \code{\link{nuts}}, \code{run_all_diagnostics},
#'   \code{write_llm_report}
#'
#' @examples
#' \dontrun{
#' mod  <- system.file("extdata", "models", "rbc.mod", package = "dynhr")
#' data <- matrix(rnorm(100), nrow = 100, ncol = 1,
#'                dimnames = list(NULL, "y_obs"))
#' result <- run_full_estimation(
#'   mod_file  = mod,
#'   data      = data,
#'   obs_vars  = "y_obs",
#'   sampler   = "rwmh",
#'   n_draws   = 1000L,
#'   n_warmup  = 500L,
#'   seed      = 1L
#' )
#' print(result)
#' }
#' @export
run_full_estimation <- function(
    mod_file      = NULL,
    data          = NULL,
    obs_vars      = NULL,
    output_dir    = ".",
    output_prefix = "dynhr_est",
    sampler       = c("rwmh", "smc", "nuts", "dime"),
    n_draws       = 10000L,
    n_warmup      = 5000L,
    n_chains      = 4L,
    n_particles   = 2000L,
    n_walkers     = NULL,
    parallel      = FALSE,
    parallel_backend = "mirai",
    n_cores       = NULL,
    analytic_grad = FALSE,
    n_mode_iter   = 10000L,
    mode_method   = "newrat",
    mode_n_starts = NULL,
    me_variance   = 0,
    likelihood    = c("gaussian", "cumulant", "whittle", "tpf"),
    lik_init      = "auto",
    freq_band     = c(0, pi),
    system_priors = NULL,
    seed          = 42L,
    run_diag      = FALSE,
    verbose       = TRUE,
    model         = NULL,
    compiled      = NULL,
    data_col_map  = NULL,
    dates         = NULL,
    obc           = NULL,
    obc_max_inner = 10L,
    compute_smoother = FALSE,
    obc_ppf_reweight = FALSE,
    run_ramsey    = FALSE,
    ramsey_order  = 1L,
    ramsey_n_periods = 400L,
    ramsey_burn_in = 100L,
    filter_tunes  = NULL,
    heteroskedastic_shocks = NULL,
    tpf_options   = list(),
    plan          = NULL,
    ...) {
  likelihood <- match.arg(likelihood)

  sampler <- match.arg(sampler)
  t_start <- Sys.time()

  .vcat <- function(...) if (verbose) cat(...)

  .vcat("\n================================================================\n")
  .vcat("  run_full_estimation\n")
  .vcat("================================================================\n\n")

  # -------------------------------------------------------------------
  # Step 1: Parse
  # -------------------------------------------------------------------
  .vcat("-- Step 1: Parse model --\n")
  if (is.null(model)) {
    if (is.null(mod_file))
      stop("Provide either 'mod_file' or 'model'.")
    model <- parse_mod(mod_file, verbose = verbose)
  } else {
    mod_file <- model$mod_file %||% "(pre-parsed)"
  }
  .vcat(sprintf("  %d endo, %d exo, %d params, %d eqs\n",
                length(model$var_names),
                length(model$varexo_names),
                length(model$param_values),
                length(model$equations)))

  # -------------------------------------------------------------------
  # Step 2: Extract priors
  # -------------------------------------------------------------------
  .vcat("-- Step 2: Extract prior specification --\n")
  priors <- extract_prior_spec(model, verbose = verbose)
  .vcat(sprintf("  %d parameters to estimate\n", nrow(priors)))

  # -------------------------------------------------------------------
  # Step 3: Compile
  # -------------------------------------------------------------------
  .vcat("-- Step 3: Compile model --\n")
  if (is.null(compiled))
    compiled <- compile_model(model, verbose = FALSE)

  # -------------------------------------------------------------------
  # Step 4: Validate data
  # -------------------------------------------------------------------
  .vcat("-- Step 4: Load and validate data --\n")
  if (is.null(data))
    stop("Provide a data matrix or CSV path via 'data'.")
  if (is.character(data)) {
    data_raw <- read.csv(data)
    if (!is.null(data_col_map)) {
      for (mod_nm in names(data_col_map)) {
        dat_nm <- data_col_map[[mod_nm]]
        if (dat_nm %in% names(data_raw) && !(mod_nm %in% names(data_raw)))
          names(data_raw)[names(data_raw) == dat_nm] <- mod_nm
      }
    }
    if (!all(obs_vars %in% names(data_raw)))
      stop(sprintf("obs_vars not found in data: %s",
                   paste(setdiff(obs_vars, names(data_raw)), collapse = ", ")))
    data <- as.matrix(data_raw[, obs_vars])
  }
  if (is.null(colnames(data))) colnames(data) <- obs_vars
  # Validate that obs_vars are present in data columns
  if (!is.null(obs_vars) && !all(obs_vars %in% colnames(data)))
    stop(sprintf("obs_vars not found in data columns: %s",
                 paste(setdiff(obs_vars, colnames(data)), collapse = ", ")))
  .vcat(sprintf("  Data: %d x %d  (obs: %s)\n",
                nrow(data), ncol(data),
                paste(colnames(data), collapse = ", ")))

  ## Apply unified plan= if supplied.  Error if both plan and the individual
  ## args are non-NULL (ambiguous).
  ## Tier 8 item 10: plan adaptation now routes through estimation_context()
  ## so that the compiled specs live in ctx$plan for provenance.
  if (!is.null(plan)) {
    if (!inherits(plan, "dynhr_plan"))
      stop("run_full_estimation: 'plan' must be a dynhr_plan object.", call. = FALSE)
    if (!is.null(filter_tunes))
      stop("run_full_estimation: supply either 'plan' or 'filter_tunes', not both.",
           call. = FALSE)
    if (!is.null(heteroskedastic_shocks))
      stop("run_full_estimation: supply either 'plan' or 'heteroskedastic_shocks', not both.",
           call. = FALSE)
    ## Compile plan inside estimation_context (TPF check + spec resolution).
    .plan_ctx <- estimation_context(plan = plan,
                                    sample_start = model$sample_start,
                                    likelihood   = likelihood)
    filter_tunes           <- attr(.plan_ctx$plan, ".filter_tunes_spec")
    heteroskedastic_shocks <- attr(.plan_ctx$plan, ".shock_scale_spec")
    rm(.plan_ctx)
  }

  ## Apply call-level filter_tunes override (NULL = use mod-file block as-is;
  ## FALSE = disable; filter_tunes_spec = override).
  model <- .resolve_filter_tunes(model, filter_tunes)

  ## Apply call-level heteroskedastic_shocks override.
  model <- .resolve_heteroskedastic_shocks(model, heteroskedastic_shocks)

  ## Expand observables for filter_tunes (no-op when no tunes are present).
  tunes_exp <- .expand_observables_for_tunes(model, obs_vars, data)
  obs_vars  <- tunes_exp$obs_vars
  data      <- tunes_exp$Y
  me_extra  <- tunes_exp$me_extra

  ## Build shock_scale matrix from heteroskedastic_shocks block (NULL when unused).
  ## Row order aligns to dr$exo_names (Landmine 9); we use model$varexo_names
  ## here as a proxy and validate at KF call via the exo arg.
  shock_scale_mat <- .build_shock_scale_matrix(
    model, model$varexo_names, nrow(data))

  # -------------------------------------------------------------------
  # Step 5: Posterior
  # -------------------------------------------------------------------
  .vcat("-- Step 5: Build log-posterior --\n")

  # Auto-detect OBC: any equation carries an MCP tag of the form  mcp = '...' ?
  use_obc <- if (isTRUE(obc)) {
    TRUE
  } else if (isFALSE(obc)) {
    FALSE
  } else {
    any(vapply(model$equations,
               function(e) {
                 tg <- e$tag
                 !is.null(tg) && !is.na(tg) &&
                   grepl("^\\s*mcp\\s*=", tg, perl = TRUE)
               },
               logical(1L)))
  }

  obc_specs <- NULL
  if (use_obc) {
    .vcat("  OBC model detected -- using Pfeiffer-Ratto PKF log-posterior\n")
    obc_specs   <- obc_parse_tags(model)
    log_post_fn <- make_log_posterior_obc_pkf(
      model, data, priors, obs_vars, compiled,
      specs       = obc_specs,
      me_variance = me_variance,
      max_inner   = obc_max_inner
    )
  } else {
    log_post_fn <- make_log_posterior(model, data, priors, obs_vars, compiled,
                                      me_variance   = me_variance,
                                      likelihood    = likelihood,
                                      lik_init      = lik_init,
                                      me_extra      = me_extra,
                                      shock_scale   = shock_scale_mat,
                                      freq_band     = freq_band,
                                      system_priors = system_priors,
                                      ...)
  }

  ## Estimation context: single carrier for the per-estimation options.
  ## Also gates the mirai pool path: daemons recompile the STANDARD GAUSSIAN
  ## full-band posterior, so any other likelihood (whittle/cumulant/tpf) or a
  ## restricted freq_band must ship the already-built log_post_fn closure.
  ## tpf_options: thread the caller-supplied list into est_ctx so it reaches
  ## the mirai pool init path (previously omitted -- gap fixed here).
  ## When use_obc = TRUE, stamp the OBC filter type that was actually used so
  ## that conditional_forecast() dispatches to the correct terminal-state path
  ## (Tier 10 item 5; Tier 15 §C: ppf/copf stamping).
  ## run_full_estimation always calls make_log_posterior_obc_pkf; ctx$likelihood
  ## is "pkf".  If the caller pre-built a PPF/COPF posterior and passes
  ## obc_filter = "ppf"/"copf" via ..., that overrides the "pkf" default.
  .obc_filter_type <- if (use_obc) {
    ft <- list(...)$obc_filter %||% "pkf"
    match.arg(as.character(ft), c("pkf", "ppf", "copf"))
  } else likelihood
  est_ctx <- estimation_context(
    me_variance   = me_variance,
    likelihood    = .obc_filter_type,
    lik_init      = lik_init,
    me_extra      = me_extra,
    shock_scale   = shock_scale_mat,
    freq_band     = freq_band,
    system_priors = system_priors,
    tpf_options   = tpf_options,
    obc_specs     = obc_specs
    ## plan provenance: plan= not passed here because me_extra/shock_scale are
    ## already resolved matrices at this point; $plan is set directly below.
  )
  ## Store plan provenance in est_ctx (Tier 8 item 10).
  if (!is.null(plan)) est_ctx$plan <- plan
  pool_gaussian <- !use_obc && identical(likelihood, "gaussian") &&
                   isTRUE(all.equal(freq_band, c(0, pi)))

  # -------------------------------------------------------------------
  # Step 5.5: TPF PMCMC preflight (when likelihood = "tpf")
  # Run BEFORE mode-finding so the user sees the variance estimate early.
  # The preflight function must use seed = NULL to vary RNG each call.
  # (Landmine 1: a non-NULL seed would make the PF deterministic, giving SD=0.)
  # -------------------------------------------------------------------
  tpf_preflight_result <- NULL
  if (!use_obc && identical(likelihood, "tpf")) {
    pf_K      <- tpf_options$pmcmc_preflight_K    %||% 30L
    pf_skip   <- isTRUE(tpf_options$pmcmc_preflight_skip)
    pf_n_part <- tpf_options$n_particles %||% 1000L

    if (!pf_skip && pf_K > 0L) {
      .vcat(sprintf("-- Step 5.5: TPF PMCMC preflight (K = %d) --\n", pf_K))
      ## Build a seed=NULL version of the TPF closure for variance measurement.
      ## CRITICAL (Landmine 1): seed MUST be NULL so each evaluation draws
      ## different RNG streams; non-NULL seed makes the PF deterministic (SD=0).
      pf_tpf_args <- modifyList(tpf_options, list(seed = NULL))
      pf_log_post_fn <- do.call(
        make_log_posterior_tpf,
        c(list(model       = model,
               data        = if (ncol(data) == length(obs_vars)) t(data) else data,
               prior_spec  = priors,
               obs_vars    = obs_vars,
               compiled    = compiled,
               me_variance = me_variance),
          pf_tpf_args)
      )
      theta_pf <- setNames(priors$mean, priors$name)
      tpf_preflight_result <- .tpf_pmcmc_preflight(
        pf_log_post_fn, theta_pf, K = pf_K, verbose = verbose)

      if (!is.na(tpf_preflight_result$sd) && tpf_preflight_result$sd > 1) {
        n_needed <- ceiling(pf_n_part * tpf_preflight_result$n_needed_factor)
        warning(sprintf(
          "TPF loglik SD at prior mean = %.2f > 1 (Dynare threshold).\n",
          tpf_preflight_result$sd),
          "PMCMC acceptance will be dominated by loglik noise.\n",
          sprintf("Current n_particles = %d; to achieve SD < 1, raise to ~%d.\n",
                  pf_n_part, n_needed),
          call. = FALSE)
        tpf_preflight_result$n_needed <- n_needed
      } else if (!is.na(tpf_preflight_result$sd)) {
        .vcat(sprintf("  TPF loglik SD = %.3f (< 1 threshold, OK)\n",
                      tpf_preflight_result$sd))
      }
    }
  }

  # -------------------------------------------------------------------
  # Step 6: Mode-finding
  # -------------------------------------------------------------------
  .vcat(sprintf("-- Step 6: Mode-finding (method = '%s') --\n", mode_method))
  theta_init <- setNames(priors$mean, priors$name)
  set.seed(seed)

  ## Use .ctx_is_standard_gaussian for the mode-finding parallel-path guard.
  ## me_extra/shock_scale are already built above; freq_band from the new arg.
  .rfe_mode_ctx <- estimation_context(
    me_variance = me_variance,
    likelihood  = likelihood,
    me_extra    = me_extra,
    shock_scale = shock_scale_mat,
    freq_band   = freq_band
  )
  par_standard_mode <- .ctx_is_standard_gaussian(.rfe_mode_ctx, use_obc = use_obc)
  rm(.rfe_mode_ctx)
  use_par_mode <- isTRUE(parallel) && requireNamespace("mirai", quietly = TRUE)

  if (use_par_mode) {
    .vcat(sprintf("  Parallel multi-start mode-finding (mirai)%s\n",
                  if (!par_standard_mode) " [closure-shipped]" else ""))
    par_mode <- run_mode_mirai(
      parsed_model = if (par_standard_mode) model else NULL,
      Y            = if (par_standard_mode) data  else NULL,
      prior_spec   = priors,
      obs_names    = if (par_standard_mode) obs_vars else NULL,
      theta_init   = theta_init,
      n_chains     = mode_n_starts,
      nm_maxit     = n_mode_iter,
      method       = mode_method,
      n_cores      = n_cores,
      me_variance  = me_variance,
      me_extra     = me_extra,
      shock_scale  = shock_scale_mat,
      log_post_fn  = if (par_standard_mode) NULL else log_post_fn,
      progress     = verbose
    )
    mode_res <- par_mode$best
    mode_res$multistart <- par_mode
  } else {
    mode_res <- .run_mode_finding(log_post_fn, theta_init, priors,
                                  nm_maxit = n_mode_iter,
                                  method   = mode_method,
                                  verbose  = verbose)
  }
  if (is.null(mode_res) || !is.finite(mode_res$logpost))
    stop("Mode-finding failed or returned non-finite log-posterior.")
  theta_mode <- mode_res$theta_mode
  .vcat(sprintf("  logpost at mode: %.4f\n", mode_res$logpost))

  # Save mode
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, output_prefix)
  mode_path <- paste0(prefix, "_mode.rds")
  saveRDS(mode_res, mode_path)
  .vcat(sprintf("  Mode saved -> %s\n", mode_path))

  # -------------------------------------------------------------------
  # Step 7: Sampling
  # -------------------------------------------------------------------
  chains    <- NULL
  conv_res  <- NULL

  if (n_draws > 0L) {
    .vcat(sprintf("-- Step 7: Sampling [%s] --\n", toupper(sampler)))
    set.seed(seed + 1L)

    chains <- switch(sampler,

      rwmh = {
        n_par      <- length(theta_mode)
        opt_scale  <- 2.38^2 / n_par
        Sigma_prop <- if (!is.null(mode_res$V_mode))
                        mode_res$V_mode * opt_scale
                      else
                        diag(priors$std^2, nrow = n_par) * opt_scale
        rownames(Sigma_prop) <- colnames(Sigma_prop) <- names(theta_mode)

        if (n_chains == 1L) {
          res <- rwmh(log_post_fn, theta_mode, Sigma_prop,
                      n_draws     = n_draws + n_warmup,
                      n_burn      = n_warmup, ...)
          new_dynhr_chains(res, "rwmh")
        } else {
          # Parallel multi-chain RWMH via the mirai daemon pool. Standard
          # Gaussian models recompile the posterior per daemon; OBC/PKF
          # models ship the already-built log_post_fn closure once instead.
          use_par <- isTRUE(parallel) &&
                     identical(parallel_backend, "mirai") &&
                     requireNamespace("mirai", quietly = TRUE)
          if (use_par) {
            .vcat(sprintf("  Parallel RWMH (mirai): %d chains\n", n_chains))
            par_res <- run_mcmc_mirai(
              parsed_model = if (pool_gaussian) model else NULL,
              Y            = if (pool_gaussian) data else NULL,
              prior_spec   = priors, obs_names = obs_vars,
              theta_mode   = theta_mode, Sigma_prop = Sigma_prop,
              n_chains     = n_chains,
              n_draws      = n_draws, n_burn = n_warmup,
              seed_base    = seed, n_cores = n_cores,
              me_variance  = me_variance,
              me_extra     = me_extra,
              shock_scale  = shock_scale_mat,
              ctx          = est_ctx,
              log_post_fn  = if (pool_gaussian) NULL else log_post_fn,
              progress     = verbose)
            chain_list  <- par_res$chains
            chain_stats <- par_res$chain_stats
          } else {
            chain_list  <- vector("list", n_chains)
            chain_stats <- data.frame(
              chain = integer(), accept_rate = numeric(),
              final_logpost = numeric(), stringsAsFactors = FALSE
            )
            L <- .robust_chol(Sigma_prop, n_par)
            for (ch in seq_len(n_chains)) {
              set.seed(seed + ch)
              if (ch == 1L) {
                th0 <- theta_mode
              } else {
                z   <- rnorm(n_par)
                th0 <- theta_mode + 0.4 * as.numeric(L %*% z)
                names(th0) <- names(theta_mode)
                for (i in seq_along(th0)) {
                  th0[i] <- max(th0[i], priors$lower[i] + 1e-8)
                  th0[i] <- min(th0[i], priors$upper[i] - 1e-8)
                }
                if (!is.finite(log_post_fn(th0)$logpost)) th0 <- theta_mode
              }
              .vcat(sprintf("  Chain %d/%d...\n", ch, n_chains))
              chain_list[[ch]] <- rwmh(log_post_fn, th0, Sigma_prop,
                                       n_draws = n_draws + n_warmup,
                                       n_burn  = n_warmup, ...)
              if (!is.null(chain_list[[ch]]))
                chain_stats <- rbind(chain_stats, data.frame(
                  chain         = ch,
                  accept_rate   = chain_list[[ch]]$acceptance_rate,
                  final_logpost = tail(chain_list[[ch]]$post_logpost, 1),
                  stringsAsFactors = FALSE
                ))
            }
          }
          # Compute convergence
          conv_res <- .compute_convergence(chain_list)
          # Package the combined chain as a dynhr_chains object
          combined_res <- list(
            chain           = conv_res$combined,
            acceptance_rate = mean(chain_stats$accept_rate, na.rm = TRUE),
            sampler         = "rwmh",
            n_chains        = n_chains,
            chain_list      = chain_list,
            chain_stats     = chain_stats
          )
          new_dynhr_chains(combined_res, "rwmh")
        }
      },

      smc = {
        use_par <- isTRUE(parallel) &&
                   identical(parallel_backend, "mirai") &&
                   requireNamespace("mirai", quietly = TRUE)
        res <- if (use_par)
          run_smc_mirai(parsed_model = if (pool_gaussian) model else NULL,
                        Y = if (pool_gaussian) data else NULL,
                        prior_spec = priors, obs_names = obs_vars,
                        n_particles = n_particles, seed_base = seed,
                        n_cores = n_cores, me_variance = me_variance,
                        me_extra = me_extra,
                        shock_scale = shock_scale_mat,
                        ctx = est_ctx,
                        log_post_fn = if (pool_gaussian) NULL else log_post_fn,
                        verbose = verbose, ...)
        else
          dynhr_smc(log_post_fn, prior_spec = priors,
                    n_particles = n_particles, ...)
        ## Resample the weighted particle cloud to an equally-weighted draw
        ## matrix so all downstream consumers (diagnostics, Bayesian IRF,
        ## smoother) receive a standard draw matrix.  The original particles
        ## and weights remain in $particles and $smc_weights.
        if (.is_smc_weighted(res)) {
          res$chain <- as_posterior_draws(res, seed = seed)
          res$n_draws <- nrow(res$chain)
          if (verbose)
            .vcat(sprintf("  SMC: resampled %d particles to %d equally-weighted draws\n",
                          res$n_particles, res$n_draws))
        }
        new_dynhr_chains(res, "smc")
      },

      nuts = {
        # Diagonal mass preconditioning from the inverse-Hessian (V_mode) is
        # essential here: without it NUTS runs on an identity metric and a DSGE
        # posterior whose marginal variances span orders of magnitude forces a
        # tiny step size. (The single-chain path previously omitted this.)
        n_par <- length(theta_mode)
        V <- if (!is.null(mode_res$V_mode)) mode_res$V_mode
             else diag(priors$std^2, nrow = n_par)
        rownames(V) <- colnames(V) <- names(theta_mode)
        nuts_mass <- 1 / pmax(diag(V), 1e-12)

        # The analytic gradient path (tangent/adjoint KF) has no me_extra or
        # shock_scale support: the gradients would be of a different likelihood.
        # Fall back to the numerical gradient. Guard via .ctx_allows_analytic_gradient.
        {
          if (isTRUE(analytic_grad) && !.ctx_allows_analytic_gradient(est_ctx)) {
            warning("analytic_grad ignored: the analytic gradient path does ",
                    "not support per-period me_extra (filter_tunes) or ",
                    "shock_scale (heteroskedastic_shocks). Using the numerical gradient.",
                    call. = FALSE)
            analytic_grad <- FALSE
          }
        }

        use_par <- isTRUE(parallel) && n_chains > 1L &&
                   identical(parallel_backend, "mirai") &&
                   requireNamespace("mirai", quietly = TRUE)
        if (use_par) {
          # Analytic/implicit gradient on the parallel path requires a
          # standard Gaussian model: each daemon needs .worker_model/
          # .worker_cm/.worker_Y from .mirai_pool_init, which only runs when
          # parsed_model/Y are supplied (!use_obc). OBC/PKF ships a pre-built
          # log_post_fn closure instead (.mirai_pool_closure) and has no
          # analytic gradient -- mirrors the serial branch's
          # `analytic_grad && !use_obc` gate below.
          par_analytic_grad <- isTRUE(analytic_grad) && !use_obc
          par_grad_method <- .dynhr_opt("grad_method", default = "hybrid")
          if (isTRUE(analytic_grad) && use_obc)
            .vcat("  [NUTS] analytic_grad requires a standard Gaussian model (compiled per daemon); ",
                  "this OBC/PKF run uses the numerical gradient.\n")
          else if (par_analytic_grad)
            .vcat(sprintf("  [NUTS] parallel chains will build analytic gradients (%s) per daemon...\n",
                          par_grad_method))
          par_res <- run_nuts_mirai(
            parsed_model = if (pool_gaussian) model else NULL,
            Y = if (pool_gaussian) data else NULL,
            prior_spec = priors, obs_names = obs_vars,
            theta_mode = theta_mode, Sigma_prop = V,
            n_chains = n_chains, n_draws = n_draws, n_warmup = n_warmup,
            seed_base = seed, n_cores = n_cores,
            me_variance = me_variance,
            me_extra = me_extra,
            shock_scale = shock_scale_mat,
            ctx = est_ctx,
            log_post_fn = if (pool_gaussian) NULL else log_post_fn,
            analytic_grad = par_analytic_grad,
            grad_method = par_grad_method,
            progress = verbose)
          chain_list <- par_res$chains
          conv_res <- .compute_convergence(chain_list)
          new_dynhr_chains(list(
            chain           = conv_res$combined,
            acceptance_rate = mean(par_res$chain_stats$accept_rate, na.rm = TRUE),
            sampler         = "nuts", n_chains = n_chains,
            n_divergent     = sum(par_res$chain_stats$n_divergent, na.rm = TRUE),
            mean_treedepth  = mean(par_res$chain_stats$mean_treedepth, na.rm = TRUE),
            chain_list      = chain_list,
            chain_stats     = par_res$chain_stats), "nuts")
        } else {
          # Exact analytic gradient (C++ Kalman score for shock-std params +
          # numerical for the rest, or full implicit-differentiation gradient
          # when grad_method = "implicit") when requested on a standard
          # Gaussian model.
          nuts_grad <- if (isTRUE(analytic_grad) && !use_obc) {
            grad_method <- .dynhr_opt("grad_method", default = "hybrid")
            .vcat(sprintf("  [NUTS] building analytic gradient (%s)...\n", grad_method))
            make_posterior_grad(model, data, priors, obs_vars, compiled,
                                me_variance = me_variance,
                                me_extra    = me_extra,
                                shock_scale = shock_scale_mat,
                                grad_method = grad_method,
                                likelihood  = likelihood,
                                freq_band   = freq_band)
          } else NULL
          res <- dynhr_nuts(log_post_fn, theta_mode,
                            n_draws  = n_draws,
                            n_warmup = n_warmup,
                            mass_diag = nuts_mass, grad_fn = nuts_grad, ...)
          new_dynhr_chains(res, "nuts")
        }
      },

      dime = {
        use_par <- isTRUE(parallel) &&
                   identical(parallel_backend, "mirai") &&
                   requireNamespace("mirai", quietly = TRUE)
        res <- if (use_par)
          run_dime_mirai(
            parsed_model = if (pool_gaussian) model else NULL,
            Y            = if (pool_gaussian) data else NULL,
            prior_spec   = priors, obs_names = obs_vars,
            n_chain      = n_walkers,
            n_iter       = n_draws, n_burn = n_warmup,
            seed_base    = seed, n_cores = n_cores,
            me_variance  = me_variance,
            me_extra     = me_extra,
            shock_scale  = shock_scale_mat,
            ctx          = est_ctx,
            log_post_fn  = if (pool_gaussian) NULL else log_post_fn,
            verbose      = verbose)
        else
          run_dime(log_post_fn, prior_spec = priors,
                   n_chain = n_walkers,
                   n_iter  = n_draws, n_burn = n_warmup,
                   verbose = verbose, ...)
        new_dynhr_chains(res, "dime")
      }
    )

    # Compute convergence for single-chain samplers if not done yet
    if (is.null(conv_res) && sampler == "rwmh" && n_chains == 1L) {
      conv_res <- list(rhat = NULL, ess = NULL, combined = chains$chain)
    } else if (is.null(conv_res)) {
      conv_res <- list(rhat = NULL, ess = NULL, combined = chains$chain)
    }

    # Save chains
    chains_path <- paste0(prefix, "_chains.rds")
    saveRDS(chains, chains_path)
    .vcat(sprintf("  Chains saved -> %s\n", chains_path))

    # Log convergence
    if (!is.null(conv_res$rhat)) {
      .vcat(sprintf("  R-hat: max=%.3f  (>1.05: %d/%d)\n",
                    max(conv_res$rhat),
                    sum(conv_res$rhat > 1.05), length(conv_res$rhat)))
      .vcat(sprintf("  ESS:   median=%.0f  min=%.0f\n",
                    median(conv_res$ess), min(conv_res$ess)))
    }
  } else {
    .vcat("-- Step 7: Sampling skipped (n_draws = 0) --\n")
  }

  # -------------------------------------------------------------------
  # Step 8.25: Ramsey policy workflow (optional)
  # -------------------------------------------------------------------
  ramsey_res <- NULL
  if (isTRUE(run_ramsey)) {
    .vcat("-- Step 8.25: Ramsey policy workflow --\n")
    params_mode <- model$param_values
    for (nm in names(theta_mode)) {
      if (nm %in% names(params_mode)) params_mode[[nm]] <- theta_mode[[nm]]
    }
    ramsey_res <- ramsey_policy(
      model = model,
      compiled = compiled,
      params = params_mode,
      order = ramsey_order,
      n_periods = ramsey_n_periods,
      burn_in = ramsey_burn_in,
      verbose = FALSE
    )
  }

  # -------------------------------------------------------------------
  # Step 8: Diagnostics (optional)
  # -------------------------------------------------------------------
  diag_results <- NULL
  if (run_diag && !is.null(chains)) {
    .vcat("-- Step 8: Diagnostics --\n")
    dr_mode <- {
      orig_pv <- model$param_values
      for (nm in names(theta_mode))
        if (nm %in% names(model$param_values))
          model$param_values[[nm]] <- theta_mode[nm]
      dr <- solve_perturbation(compiled, solve_steady(compiled))
      model$param_values <- orig_pv
      dr
    }

    diag_results <- run_all_diagnostics(
      model       = model,
      compiled    = compiled,
      draws       = if (!is.null(chains)) chains$chain else NULL,
      chains_list = if (!is.null(chains)) list(chains$chain) else NULL,
      data        = data,
      dates       = dates,
      ramsey_result = ramsey_res
    )
    if (!is.null(diag_results))
      .vcat(sprintf("  %d diagnostics run\n", length(diag_results)))
  }

  # -------------------------------------------------------------------
  # Step 8.5: OBC smoother + historical decomposition at mode (optional)
  # -------------------------------------------------------------------
  smoother_res <- NULL
  if (compute_smoother && use_obc) {
    .vcat("-- Step 8.5: PKF smoother + historical decomposition at mode --\n")
    smoother_res <- .compute_smoother_at_mode(
      theta_mode, model, compiled, data, obs_vars,
      obc_specs, me_variance, obc_max_inner
    )
    if (!is.null(smoother_res))
      .vcat(sprintf("  Smoothed states: %d x %d, Smoothed shocks: %d x %d, HD: %d shocks+constraint\n",
                    nrow(smoother_res$smoothed_states), ncol(smoother_res$smoothed_states),
                    nrow(smoother_res$smoothed_shocks), ncol(smoother_res$smoothed_shocks),
                    length(smoother_res$hd$contributions)))
  }

  # -------------------------------------------------------------------
  # Step 8.6: PPF importance re-weighting (PKF adequacy check, opt-in)
  # -------------------------------------------------------------------
  ppf_reweight_res <- NULL
  if (use_obc && !is.null(chains) && isTRUE(obc_ppf_reweight)) {
    .vcat("-- Step 8.6: PPF importance re-weighting (PKF adequacy check) --\n")
    ppf_reweight_res <- tryCatch(
      ppf_reweight_posterior(
        chains, model, compiled, data, priors, obs_vars,
        specs       = obc_specs,
        n_thin      = 10L,
        n_particles = 1000L,
        me_variance = max(me_variance, 1e-6),
        seed        = seed
      ),
      error = function(e2) {
        warning("Step 8.6 ppf_reweight_posterior failed: ", conditionMessage(e2),
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(ppf_reweight_res)) {
      if (ppf_reweight_res$ess_fraction < 0.5)
        warning(sprintf(
          "PKF may be inadequate for this dataset (ESS/n = %.2f < 0.50). ",
          ppf_reweight_res$ess_fraction),
          "Consider re-estimating with make_log_posterior_obc_ppf ",
          "(proposal='bootstrap' or 'copf').",
          call. = FALSE)
      .vcat(sprintf("  PPF reweight: %s\n", ppf_reweight_res$verdict))
    }
  }

  # -------------------------------------------------------------------
  # Wrap up
  # -------------------------------------------------------------------
  t_end   <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

  .vcat("\n================================================================\n")
  .vcat(sprintf("  DONE in %.1f min\n", elapsed))
  .vcat("================================================================\n\n")

  result <- list(
    model           = model,
    compiled        = compiled,
    prior_spec      = priors,
    data            = data,
    mode            = mode_res,
    chains          = chains,
    convergence     = conv_res,
    diagnostics     = diag_results,
    obc_specs       = obc_specs,
    smoother        = smoother_res,
    ppf_reweight    = ppf_reweight_res,
    ramsey          = ramsey_res,
    tpf_preflight   = tpf_preflight_result,
    ctx             = est_ctx,
    meta            = list(
      mod_file    = mod_file,
      obs_vars    = obs_vars,
      sampler     = sampler,
      n_draws     = n_draws,
      n_warmup    = n_warmup,
      n_chains    = n_chains,
      seed        = seed,
      use_obc     = use_obc,
      run_ramsey  = isTRUE(run_ramsey),
      started_at  = format(t_start, "%Y-%m-%d %H:%M:%S"),
      elapsed_min = elapsed,
      dynhr_version = utils::packageVersion("dynhr")
    )
  )
  class(result) <- c("dynhr_estimation_result", "list")
  invisible(result)
}
