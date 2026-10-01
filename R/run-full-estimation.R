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
  if (is.null(x$mode))
    cat("  Mode      : not run (prior-initialised sampler)\n")
  else
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
    ## THAMES cross-check (guarded: a numerical failure cannot break print;
    ## a programming error is re-raised)
    thames_line <- tryCatch({
      tr <- thames_mdd_from_chains(x$chains)
      if (is.finite(tr$log_mdd))
        sprintf("  log p(Y|M) [THAMES]: %.3f +/- %.3f\n", tr$log_mdd, tr$se)
      else
        NULL
    }, error = function(e) .dynhr_reraise_bug(e, NULL))
    if (!is.null(thames_line)) cat(thames_line)
  } else {
    cat("  MCMC      : skipped (n_draws = 0)\n")
  }
  if (!is.null(x$diagnostics))
    cat(sprintf("  Diagnostics: %d run\n", length(x$diagnostics)))
  il <- .est_integrity_lines(x$provenance$integrity)
  if (length(il)) cat(il, sep = "\n")
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
  ## (no-op for non-SSM-parameter models).
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
#' @section One code path:
#' The arguments become a \code{\link{dynhr_estimation_spec}}
#' (\code{\link{as_estimation_spec}}), run by \code{\link{run_estimation}}:
#' the mode stage is exactly \code{\link{run_mode_finding}}'s and the sampler
#' stage exactly \code{\link{run_posterior_estimation}}'s, so
#' \code{run_full_estimation(...)} returns the mode and the draws of
#' \code{run_mode_finding()} followed by \code{run_posterior_estimation(seed =
#' seed)} with the same settings. In particular it honours the
#' result-changing options (\code{transform_params}, \code{rwmh_adapt_cov},
#' \code{rwmh_n_blocks}, \code{proposal_cov_method},
#' \code{use_exact_hessian}, ...; see \code{\link{dynhr_set_options}}), which
#' its own mode and sampler code ignored before dynhr 0.9.3.105. The
#' particle-likelihood variance preflight is evaluated at the mode.
#'
#' @section Arguments by stage:
#' The argument list is long because one call spans the whole pipeline.  Every
#' argument has a default; a minimal call needs only \code{mod_file},
#' \code{data} and \code{obs_vars}.  Grouped by the stage they act on:
#'
#' \describe{
#'   \item{\strong{Model and data}}{\code{mod_file}, \code{data},
#'     \code{obs_vars}, \code{data_col_map}, \code{dates}, \code{model},
#'     \code{compiled}}
#'   \item{\strong{Likelihood}}{\code{likelihood}, \code{me_variance},
#'     \code{lik_init}, \code{freq_band}, \code{system_priors},
#'     \code{tpf_options}, \code{filter_tunes},
#'     \code{heteroskedastic_shocks}, \code{stochastic_volatility},
#'     \code{plan}}
#'   \item{\strong{Mode-finding}}{\code{n_mode_iter}, \code{mode_method},
#'     \code{mode_n_starts}}
#'   \item{\strong{Sampler}}{\code{sampler}, \code{n_draws}, \code{n_warmup},
#'     \code{n_chains}, \code{n_particles}, \code{n_walkers},
#'     \code{analytic_grad}, \code{seed}, \code{checkpoint_dir},
#'     \code{resume}, \code{on_mismatch}, \code{...}}
#'   \item{\strong{Parallelism}}{\code{parallel}, \code{parallel_backend},
#'     \code{n_cores}}
#'   \item{\strong{Occasionally-binding constraints}}{\code{obc},
#'     \code{obc_max_inner}, \code{obc_ppf_reweight},
#'     \code{compute_smoother}}
#'   \item{\strong{Optimal policy}}{\code{run_ramsey}, \code{ramsey_order},
#'     \code{ramsey_n_periods}, \code{ramsey_burn_in}}
#'   \item{\strong{Diagnostics}}{\code{run_diag}}
#'   \item{\strong{Output}}{\code{output_dir}, \code{output_prefix},
#'     \code{verbose}}
#' }
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
#'   \code{"smc"}, \code{"nuts"}, \code{"dime"}, \code{"pmmh"} (random-walk
#'   MH over an unbiased particle likelihood), \code{"hmc"}, \code{"mala"},
#'   \code{"chees"}, \code{"dsmh"} or \code{"smc2"} -- every sampler of
#'   \code{\link{sampler_spec}}. Arguments a sampler does not take (e.g.
#'   \code{n_chains} for \code{"hmc"}) are left out, with a message.
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
#' @param analytic_grad For the gradient samplers (NUTS, HMC, MALA, ChEES), use
#'   the analytic gradient from \code{\link{make_posterior_grad}} instead of a
#'   finite-difference one. Its method is the \code{grad_method} option,
#'   default \code{"auto"}: the exact \code{"adjoint_solution"} for the
#'   Gaussian likelihood. Default \code{TRUE}: the exact gradient wherever the
#'   likelihood has one; a likelihood without one (the OBC filters,
#'   \code{"tpf"}, \code{"pskf"}, ...) keeps the numerical gradient without a
#'   warning. \code{FALSE} always finite-differences the log-posterior.
#'   Applies to the serial path and to the parallel multi-chain NUTS path
#'   (standard Gaussian models).
#' @param metric For \code{sampler = "nuts"}: the warmup mass-matrix
#'   adaptation. \code{"diagonal"} (default) starts from the inverse-Hessian
#'   diagonal at the mode and adapts a diagonal inverse mass equal to the
#'   warmup posterior variance (Stan's convention); \code{"warmup_dense"}
#'   ends warmup with a Ledoit-Wolf dense inverse mass; \code{"fisher_diag"}
#'   (opt-in) sets the diagonal inverse mass to
#'   \eqn{\sqrt{\mathrm{var}(x)/\mathrm{var}(\nabla \log p)}} from each
#'   window's draws and gradients (Seyboldt, Carlson & Carpenter 2026,
#'   arXiv:2603.18845); \code{"lowrank"} (opt-in) adapts a
#'   low-rank-plus-diagonal inverse mass from draws and gradients (same
#'   paper's Algorithm 1, Lao 2026 schedule). Applies to the serial NUTS path;
#'   the parallel multi-chain path always uses \code{"diagonal"} and warns
#'   otherwise. Ignored by the other samplers.
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
#'   noise diagonal. \code{NULL} (default) reads the \code{me_variance}
#'   option (0 unless set), as the estimation spec's
#'   \code{likelihood$me_variance} does.  Use a small positive value for
#'   stochastically singular models.
#' @param seed       Random seed for reproducibility (default \code{42L}): it
#'   seeds the mode stage (the caller's RNG stream is restored on exit) and
#'   re-seeds the sampler stage; the parallel samplers use it as their base
#'   seed.
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
#' @param ramsey_discount Planner discount factor for the Ramsey step.  When
#'   \code{NULL} (default) it is taken from the parameter \code{beta} or
#'   \code{betta}; if the model has neither, the Ramsey step is SKIPPED with a
#'   \code{dynhr_warn_ramsey_skipped} warning (\code{$ramsey} stays
#'   \code{NULL}) instead of aborting the whole estimation.
#' @param likelihood  Likelihood type: \code{"gaussian"} (Kalman filter, default),
#'   \code{"cumulant"}, \code{"whittle"}, \code{"tpf"} (Tempered Particle
#'   Filter), \code{"sv_rbpf"}, \code{"pskf"}, \code{"student_t"} (pass
#'   \code{student_df} in \code{...}), \code{"pruned"}, \code{"global_pf"}, or
#'   the OBC filters \code{"pkf"}, \code{"ppf"} and \code{"copf"} (which pick
#'   the filter of an OBC model) -- every likelihood of
#'   \code{\link{likelihood_spec}}.
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
#' @param stochastic_volatility  Call-level stochastic-volatility override.
#'   \code{NULL} (default) uses the \code{stochastic_volatility} block from
#'   the \code{.mod} file if present; supply one to attach or replace the
#'   SV specification for this call only, leaving the model object untouched.
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
#' @param checkpoint_dir Optional directory for streamed, restartable
#'   sampling (the samplers that support it are listed at
#'   \code{\link{run_posterior_estimation}}).
#' @param resume When \code{TRUE}, continue the chains saved in
#'   \code{checkpoint_dir}, adding \code{n_draws} draws.
#' @param on_mismatch What a \code{resume} does when the checkpoint's target
#'   differs from this call's, or a registered result change touches this
#'   run: \code{"refuse"} (default) or \code{"warn"} (resume; the result is
#'   marked). See \code{\link{run_estimation}}, section Checkpoints.
#' @param ...        Additional arguments forwarded to the chosen sampler
#'   (\code{student_df} and \code{obc_filter} go to the likelihood).
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
#'     \item{\code{resolved}}{What the sampler stage resolved at run time:
#'       \code{grad_method}, the analytic-gradient method it used (named by
#'       sampler, e.g. \code{"adjoint_solution"} for a requested
#'       \code{"auto"}), and \code{grad_method_requested}; \code{NULL} when
#'       no sampler ran. The run record copies it.}
#'     \item{\code{meta}}{Run metadata: mod_file, sampler, seed, timestamps, use_obc}
#'     \item{\code{run_record}}{A \code{dynhr_run_record}: resolved arguments,
#'       option snapshot, RNG state and provenance; replay it with
#'       \code{\link{dynhr_rerun}}.}
#'   }
#'
#' @seealso \code{\link{parse_mod}}, \code{\link{find_mode}}, \code{\link{dynhr_mcmc}},
#'   \code{\link{smc}}, \code{\link{nuts}}, \code{run_all_diagnostics},
#'   \code{write_llm_report}, \code{\link{run_estimation}} (the spec runner
#'   behind this wrapper), \code{\link{dynhr_rerun}}, \code{\link{dynhr_verify}}
#'
#' @examples
#' \donttest{
#' ## Model and data both ship with the package: nk_demo.mod is a textbook
#' ## three-equation New Keynesian model (nine estimated parameters) and
#' ## nk_demo_data.csv is 200 periods simulated from it at known parameters,
#' ## so the run has a right answer to be checked against.
#' mod_file <- system.file("extdata/models/nk_demo.mod", package = "dynhr")
#' est_data <- read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                  package = "dynhr"))
#'
#' fit <- run_full_estimation(
#'   mod_file    = mod_file,
#'   data        = est_data,
#'   obs_vars    = c("ygr", "infl", "intr"),
#'   output_dir  = tempdir(),
#'   sampler     = "smc",        # no mode required; returns log p(Y | M)
#'   n_particles = 500L,         # raise for a production run
#'   seed        = 1L,
#'   verbose     = FALSE
#' )
#'
#' summary(fit$chains)
#' fit$chains$log_marginal_lik
#' }
#' @export
run_full_estimation <- function(
    mod_file      = NULL,
    data          = NULL,
    obs_vars      = NULL,
    output_dir    = ".",
    output_prefix = "dynhr_est",
    sampler       = c("rwmh", "smc", "nuts", "dime", "pmmh", "hmc", "mala",
                      "chees", "dsmh", "smc2"),
    n_draws       = 10000L,
    n_warmup      = 5000L,
    n_chains      = 4L,
    n_particles   = 2000L,
    n_walkers     = NULL,
    parallel      = FALSE,
    parallel_backend = "mirai",
    n_cores       = NULL,
    analytic_grad = TRUE,
    metric        = c("diagonal", "warmup_dense", "fisher_diag", "lowrank"),
    n_mode_iter   = 10000L,
    mode_method   = "newrat",
    mode_n_starts = NULL,
    me_variance   = NULL,
    likelihood    = c("gaussian", "cumulant", "whittle", "tpf", "sv_rbpf",
                      "pskf", "student_t", "pruned", "global_pf", "pkf",
                      "ppf", "copf"),
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
    ramsey_discount = NULL,
    filter_tunes  = NULL,
    heteroskedastic_shocks = NULL,
    stochastic_volatility = NULL,
    tpf_options   = list(),
    plan          = NULL,
    checkpoint_dir = NULL,
    resume        = FALSE,
    on_mismatch   = "refuse",
    ...) {
  ## Retired spellings (0.9.3 renames): classed error naming the new one.
  .dynhr_reject_retired_args("run_full_estimation", ...names(),
                             names(sys.call()),
                             c(.dynhr_retired_count_args,
                               .dynhr_retired_data_args))
  ## Own the message epoch for this run: repeat-suppressed warnings
  ## (`.dynhr_warn(once = TRUE)`) are keyed within it and re-arm for the
  ## next run, and the close reports what it suppressed. A nested call
  ## inherits this epoch rather than opening a second one.
  .dynhr_run_epoch <- .dynhr_epoch("run_full_estimation")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)
  ## Run record (R/run-record.R): resolved args, option snapshot and RNG
  ## state at ENTRY -- before the body touches any argument or the RNG.
  .rr <- .dynhr_rr_begin("run_full_estimation", environment(), list(...))
  ## A thin wrapper. The arguments become an estimation spec
  ## (validate_spec() holds the cross-field checks) and the one spec runner
  ## runs its mode stage, sampler stage and outputs -- the same code as
  ## run_mode_finding() + run_posterior_estimation().
  spec <- as_estimation_spec(.rr$args, entry = "run_full_estimation")
  .run_estimation_impl(spec, rr = .rr)
}


# ============================================================================
# The spec runner
# ============================================================================

#' Run an estimation from its spec
#'
#' The single estimation runner: given a \code{\link{dynhr_estimation_spec}}
#' it runs the mode stage, then the sampler stage(s), then the outputs, and
#' returns the result classes of the familiar entry points.
#' \code{\link{run_mode_finding}}, \code{\link{run_posterior_estimation}} and
#' \code{\link{run_full_estimation}} are thin wrappers that build a spec from
#' their arguments (\code{\link{as_estimation_spec}}) and call the same code.
#'
#' @section Result form:
#' \code{spec$outputs$form} selects the result class: \code{"mode"} a
#' \code{dynhr_mode_result} (no sampler); \code{"posterior"} a
#' \code{dynhr_posterior_result} (the sampler(s) run from \code{mode$result}
#' when it is set, else from a fresh mode stage); \code{"full"} a
#' \code{dynhr_estimation_result} (mode, one sampler or none, then the
#' outputs: saved files, diagnostics, Ramsey, the OBC smoother). The default
#' \code{"auto"} picks \code{"mode"} without a sampler, \code{"posterior"}
#' with a precomputed mode result or a sampler sequence, else \code{"full"}.
#'
#' The prior-initialised samplers (\code{"smc"}, \code{"dsmh"},
#' \code{"dime"}, \code{"smc2"}) start from prior draws and use neither the
#' mode nor its proposal covariance, so when every sampler of the spec is one
#' of them the mode stage is not run -- unless an output is evaluated at the
#' mode (\code{outputs$diagnostics}, \code{outputs$ramsey}, the OBC
#' \code{outputs$smoother}). The result then has no mode: \code{$mode} of a
#' full result is \code{NULL}, and a posterior result's \code{mode_result}
#' carries the log-posterior but \code{theta_mode = NULL} and
#' \code{meta$mode_skipped = TRUE}. The draws are those of a run with a mode
#' stage (the sampler stage re-seeds). \code{mode$run} overrides this rule:
#' \code{"auto"} (the default) applies it, \code{"always"} runs the mode stage
#' anyway (for example to report the mode next to the draws), and
#' \code{"never"} does not run it -- an error (class
#' \code{dynhr_error_spec_mode_needed}) when a sampler or output uses the
#' mode.
#'
#' @section Options:
#' For the duration of the run the package option store holds the spec's
#' values: its option snapshot (\code{spec$options}) plus every result-changing
#' option taken from the typed field that owns it (for example
#' \code{power_posterior} from \code{likelihood$power_posterior},
#' \code{transform_params} from \code{mode$transform_params} during the mode
#' stage and from the sampler's field during sampling, \code{seed_base} from
#' \code{compute$seed}). Parallel daemons receive that store. The caller's
#' options are restored on exit, also on error.
#'
#' @section Seeds:
#' With \code{compute$seed} set, a run that includes the mode stage seeds it
#' and restores the caller's RNG stream on exit; the sampler stage re-seeds
#' with the same seed, and the parallel (mirai) samplers use it as their base
#' seed. Sampling a precomputed mode result (\code{mode$result}) seeds only
#' the sampler stage and leaves the RNG where it ends, as
#' \code{run_posterior_estimation()} does. \code{NULL} leaves the ambient RNG
#' stream untouched.
#'
#' @section Checkpoints:
#' With \code{compute$checkpoint_dir} a fresh run writes
#' \code{spec_integrity.rds} next to the chain files: the content hash of each
#' target/algorithm part (model, data, likelihood, mode, sampler without
#' \code{n_draws}) and the build and numerical environment. A resume
#' (\code{compute$resume = TRUE}) compares them: a different target part is
#' refused (error \code{dynhr_error_checkpoint_spec_mismatch}, naming the
#' differing parts and fields) unless \code{compute$on_mismatch = "warn"}. A
#' different dynhr version is looked up in the result-change registry (a
#' table of every change that alters results, with the components it
#' touches): when a registered change between the two versions touches a
#' component the continued chain uses (its likelihood, priors, sampler kernel,
#' analytic gradient, ...), the resume is refused as well (error
#' \code{dynhr_error_checkpoint_code_changed}, naming the changes and the
#' version to install) unless \code{compute$on_mismatch = "warn"}. Any other
#' code difference, or a different environment (R, platform, BLAS/LAPACK,
#' daemon count), warns: the continued chain is still a valid MCMC chain, but
#' not bit-identical to an uninterrupted run (equivalent in distribution). Any
#' such event is kept in \code{result$provenance$integrity} and the run
#' record, and \code{print()} shows it; an overridden refusal marks the result
#' as not clean. See \code{\link{dynhr_verify}} for checking a finished result
#' under the current build and environment.
#'
#' @param spec A \code{dynhr_estimation_spec} (from
#'   \code{\link{dynhr_estimation_spec}}, \code{\link{as_estimation_spec}},
#'   \code{update()} or \code{\link{read_spec}}).
#' @return A \code{dynhr_mode_result}, \code{dynhr_posterior_result} or
#'   \code{dynhr_estimation_result} (see Result form), with
#'   \code{$run_record} (schema 2: it carries the spec).
#' @seealso \code{\link{dynhr_estimation_spec}}, \code{\link{dynhr_rerun}},
#'   \code{\link{dynhr_verify}}, \code{\link{write_spec}};
#'   \code{vignette("estimation")}, section "Specs, records and reruns"
#' @examples
#' \donttest{
#' mod <- system.file("extdata/models/nk_demo.mod", package = "dynhr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr")))
#' spec <- dynhr_estimation_spec(mod, data = Y,
#'   mode = mode_spec(n_iter = 200L), sampler = NULL,
#'   compute = compute_spec(verbose = FALSE))
#' fit <- run_estimation(spec)
#' fit$theta_mode
#' }
#' @export
run_estimation <- function(spec) {
  .dynhr_run_epoch <- .dynhr_epoch("run_estimation")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)
  if (!inherits(spec, "dynhr_estimation_spec"))
    .dynhr_abort("run_estimation: `spec` must be a dynhr_estimation_spec ",
                 "(see dynhr_estimation_spec(), as_estimation_spec(), read_spec()).",
                 class = "dynhr_error_bad_argument")
  if (!identical(spec$spec_version, .spec_version))
    .dynhr_abort("run_estimation: spec_version ", format(spec$spec_version),
                 " is not supported (this dynhr reads version ", .spec_version,
                 ").", class = "dynhr_error_spec_version")
  ## before any stage runs: resuming needs a checkpoint on this machine
  resume_msg <- .spec_resume_problem(spec$compute)
  if (!is.null(resume_msg))
    .dynhr_abort("run_estimation: ", resume_msg,
                 class = c("dynhr_error_resume_no_checkpoint", "dynhr_error_spec_invalid"))
  .rr <- .dynhr_rr_begin("run_estimation", environment())
  .run_estimation_impl(spec, rr = .rr)
}

## The result form a spec asks for (outputs$form, "auto" resolved).
.est_form <- function(spec) {
  form <- spec$outputs$form
  if (!identical(form, "auto")) return(form)
  if (is.null(spec$sampler)) return("mode")
  if (!is.null(spec$mode$result) ||
      inherits(spec$sampler, "dynhr_sampler_sequence")) return("posterior")
  "full"
}

## The option store for one stage: the spec's snapshot plus each
## result-changing option from its typed home field (`samp`: the sampler
## spec whose fields apply; NULL during the mode stage). An option whose home
## this stage does not have stays at its registered default.
.est_stage_options <- function(spec, samp = NULL) {
  store <- spec$options
  attr(store, "set") <- NULL
  homes <- .spec_option_homes()
  for (op in names(homes)) {
    val <- NULL
    for (h in homes[[op]]) {
      cp  <- sub("\\$.*$", "", h)
      fld <- sub("^.*\\$", "", h)
      src <- if (identical(cp, "sampler")) samp else spec[[cp]]
      ## transform_params lives in mode and sampler: the stage's own wins
      if (identical(cp, "mode") && !is.null(samp) &&
          any(grepl("^sampler\\$", homes[[op]])) && fld %in% names(samp))
        next
      if (!is.null(src) && fld %in% names(src) && !is.null(src[[fld]]))
        val <- src[[fld]]
    }
    store[op] <- list(val)
  }
  store[!vapply(store, is.null, logical(1))]
}

## Replace the option store by `store` (restored by the caller).
.est_set_store <- function(store) {
  rm(list = ls(.dynhr_opts, all.names = TRUE), envir = .dynhr_opts)
  if (length(store)) list2env(store, envir = .dynhr_opts)
  invisible(NULL)
}

## The data matrix a spec's likelihood uses: the value, or the CSV (columns
## renamed by data_col_map, restricted to obs_vars); then the first_obs / nobs
## window.
.est_load_data <- function(spec) {
  obs_vars <- spec$obs_vars
  data <- spec$data$value
  if (is.null(data)) {
    data_raw <- utils::read.csv(spec$data$path)
    data_col_map <- spec$likelihood$data_col_map
    if (!is.null(data_col_map)) {
      for (mod_nm in names(data_col_map)) {
        dat_nm <- data_col_map[[mod_nm]]
        if (dat_nm %in% names(data_raw) && !(mod_nm %in% names(data_raw)))
          names(data_raw)[names(data_raw) == dat_nm] <- mod_nm
      }
    }
    if (!all(obs_vars %in% names(data_raw)))
      .dynhr_abort(sprintf("obs_vars not found in data: %s",
                           paste(setdiff(obs_vars, names(data_raw)), collapse = ", ")),
                   class = "dynhr_error_spec_invalid")
    data <- as.matrix(data_raw[, obs_vars])
  }
  if (is.null(colnames(data))) colnames(data) <- obs_vars
  lik <- spec$likelihood
  if (lik$first_obs > 1L || !is.null(lik$nobs)) {
    last <- if (is.null(lik$nobs)) nrow(data) else lik$first_obs + lik$nobs - 1L
    if (lik$first_obs > nrow(data) || last > nrow(data))
      .dynhr_abort("run_estimation: the likelihood sample (first_obs = ",
                   lik$first_obs, ", nobs = ", format(lik$nobs %||% "all"),
                   ") runs past the ", nrow(data), " data rows.",
                   class = "dynhr_error_spec_invalid")
    data <- data[seq.int(lik$first_obs, last), , drop = FALSE]
  }
  data
}

## Model, compiled model, data and the `solved` object the mode result keeps.
.est_inputs <- function(spec) {
  mp <- spec$model
  model <- mp$mod
  compiled <- mp$compiled %||% mp$solved$compiled %||%
    compile_model(model, verbose = FALSE, max_order = mp$max_order)
  solved <- mp$solved
  if (is.null(solved)) {
    solved <- structure(list(model = model, compiled = compiled),
                        class = c("dynhr_solved", "list"))
  } else if (is.null(solved$compiled)) {
    solved$compiled <- compiled
  }
  list(model = model, compiled = compiled, solved = solved,
       data = .est_load_data(spec))
}

## The log-posterior a spec defines, built by the runner's own builder: the
## mode stage run with optimise = FALSE (.est_mode_stage() then applies the
## plan / filter_tunes / heteroskedastic_shocks / stochastic_volatility
## overrides and .mod blocks, routes an OBC model to its filter, and builds
## the closure, but finds no mode and no proposal), under the spec's option
## store and OBC switch as .run_estimation_impl() sets them. Returns the
## stage's list (log_post_fn, result = the mode result without a mode,
## model, priors, data, obs_vars, ...). `spec` must be valid. dm_posterior()
## builds through this, so its closure is the runner's objective.
.est_build_posterior <- function(spec) {
  saved_opts <- as.list(.dynhr_opts)
  on.exit(.est_set_store(saved_opts), add = TRUE)
  if (isFALSE(spec$likelihood$obc)) {
    old_obc_off <- .mod_blocks_state$obc_off
    .mod_blocks_state$obc_off <- TRUE
    on.exit(.mod_blocks_state$obc_off <- old_obc_off, add = TRUE)
  }
  .est_set_store(.est_stage_options(spec, NULL))
  .est_mode_stage(spec, .est_inputs(spec), proposal = FALSE, optimise = FALSE)
}

## A spec's mode$result may be a dynhr_mode_ref (a run record's stand-in):
## only dynhr_rerun() can rebuild it.
.est_check_mode_result <- function(spec) {
  mr <- spec$mode$result
  if (inherits(mr, "dynhr_mode_ref"))
    .dynhr_abort("run_estimation: this spec's mode$result is a reference to a ",
                 "recorded mode run, not a mode result. Replay the record with ",
                 "dynhr_rerun(), or set mode$result to a dynhr_mode_result.",
                 class = "dynhr_error_rerun_needs_mode_result")
  if (!is.null(mr) && !inherits(mr, "dynhr_mode_result") &&
      !(is.list(mr) && !is.null(mr$theta_mode) && is.function(mr$log_post_fn)))
    .dynhr_abort("run_estimation: mode$result must be a dynhr_mode_result.",
                 class = "dynhr_error_spec_invalid")
  invisible(NULL)
}

## The runner. `rr`: the run record begun by the calling entry point.
.run_estimation_impl <- function(spec, rr = NULL) {
  .est_check_mode_result(spec)
  form <- .est_form(spec)
  cmp  <- spec$compute
  seed <- cmp$seed

  ## Scoped option store: the spec's values, restored on exit (also on error).
  saved_opts <- as.list(.dynhr_opts)
  on.exit(.est_set_store(saved_opts), add = TRUE)
  samplers <- .spec_sampler_list(spec$sampler)
  ## likelihood$obc = FALSE estimates an mcp model as linear on purpose: the
  ## posterior constructors need not warn that they leave the mcp tags out.
  if (isFALSE(spec$likelihood$obc)) {
    old_obc_off <- .mod_blocks_state$obc_off
    .mod_blocks_state$obc_off <- TRUE
    on.exit(.mod_blocks_state$obc_off <- old_obc_off, add = TRUE)
  }

  ## Checkpoint integrity (plan section 7): before any stage runs.
  integrity <- .est_checkpoint_integrity(spec)

  mode_result <- spec$mode$result
  mode_out <- NULL
  inp <- NULL
  if (is.null(mode_result)) {
    ## A run that finds the mode owns the seed: seeded here, and the caller's
    ## RNG stream is restored on exit.
    .local_seed(seed)
    .est_set_store(.est_stage_options(spec, NULL))
    inp <- .est_inputs(spec)
    ## a full run that samples nothing needs no proposal covariance; the
    ## prior-initialised samplers (SMC, DSMH, DIME, SMC2) need no mode at all
    ## unless an output asks for it, or mode$run says otherwise
    ## (.spec_mode_stage_needed())
    mode_out <- .est_mode_stage(spec, inp,
                                proposal = !(identical(form, "full") &&
                                               !length(samplers)),
                                optimise = .spec_mode_stage_needed(spec))
    mode_result <- mode_out$result
  }

  if (identical(form, "mode")) {
    result <- mode_result
  } else if (identical(form, "posterior")) {
    .est_set_store(.est_stage_options(spec, samplers[[1L]]))
    result <- .est_sampler_stage(.est_sampler_spec(spec, mode_result), mode_result)
  } else {
    result <- .est_full_result(spec, mode_out, samplers, rr)
  }

  if (!is.null(integrity)) result$provenance <- list(integrity = integrity)
  result$run_record <- if (!is.null(rr))
    .dynhr_rr_finish(rr, model = inp$model %||% mode_result$solved$model,
                     data = inp$data, spec = spec, integrity = integrity,
                     result = result)
  invisible(result)
}

## Output file prefix of a full run.
.est_prefix <- function(out) file.path(out$dir, out$prefix)

## The spec the sampler stage runs: when the mode stage was skipped (prior-
## initialised samplers) there is no mode to check, so its quality check is
## off.
.est_sampler_spec <- function(spec, mode_result) {
  if (isTRUE(mode_result$meta$mode_skipped))
    spec$mode <- .spec_build("mode", list(skip_check = TRUE), base = spec$mode)
  spec
}

## The dynhr_chains object and convergence of one sampler-stage method result.
.est_chains_object <- function(chain_res, sampler) {
  cl <- chain_res$chains
  if (length(cl) == 1L) {
    ## one chain: the sampler's own result (keeps post_logpost etc.)
    ch <- new_dynhr_chains(cl[[1L]], sampler)
    return(list(chains = ch,
                convergence = list(rhat = NULL, ess = NULL, combined = ch$chain)))
  }
  if (!is.null(chain_res$combined)) {
    ## .run_rwmh_batch(): combined chains + convergence already built
    return(list(chains = chain_res$combined, convergence = chain_res$convergence))
  }
  if (!length(cl)) return(list(chains = NULL, convergence = NULL))
  conv <- .compute_convergence(cl)
  st <- chain_res$chain_stats
  ch <- new_dynhr_chains(list(
    chain           = conv$combined,
    acceptance_rate = mean(st$accept_rate, na.rm = TRUE),
    sampler         = sampler, n_chains = length(cl),
    n_divergent     = if (!is.null(st$n_divergent)) sum(st$n_divergent, na.rm = TRUE),
    mean_treedepth  = if (!is.null(st$mean_treedepth))
                        mean(st$mean_treedepth, na.rm = TRUE),
    chain_list      = cl,
    chain_stats     = st), sampler)
  list(chains = ch, convergence = conv)
}

## The full-run form (formerly the tail of run_full_estimation()): the mode
## result is `mode_out` (.est_mode_stage()); runs the one sampler (if any),
## saves the outputs and runs the post-estimation extras.
.est_full_result <- function(spec, mode_out, samplers, rr) {
  if (is.null(mode_out))
    .dynhr_abort("run_estimation: outputs$form = \"full\" runs the mode stage; ",
                 "a precomputed mode$result needs form \"posterior\".",
                 class = "dynhr_error_spec_invalid")
  cmp <- spec$compute
  out <- spec$outputs
  lik <- spec$likelihood
  verbose <- cmp$verbose
  seed <- cmp$seed
  t_start <- Sys.time()
  .vcat <- function(...) if (verbose) .dynhr_cat(...)

  mode_result <- mode_out$result
  model      <- mode_out$model
  compiled   <- mode_out$compiled
  priors     <- mode_out$priors
  data       <- mode_out$data
  obs_vars   <- mode_out$obs_vars
  use_obc    <- mode_out$use_obc
  obc_specs  <- mode_out$obc_specs
  me_variance <- lik$me_variance
  theta_mode <- mode_result$theta_mode
  mode_res   <- mode_result$mode
  s1 <- if (length(samplers)) samplers[[1L]]

  ## Labels as run_full_estimation() reported them.
  rfe <- is.list(rr) && identical(rr$fn, "run_full_estimation")
  mod_file <- if (rfe && !is.null(rr$args$mod_file)) rr$args$mod_file else
    spec$model$path %||% model$mod_file %||% "(pre-parsed)"
  sampler_lbl <- if (!is.null(s1)) s1$method else if (rfe) rr$args$sampler[[1L]]
  n_draws  <- if (rfe) rr$args$n_draws  else s1$n_draws %||% 0L
  n_warmup <- if (rfe) rr$args$n_warmup else s1$n_warmup %||% 0L
  n_chains <- if (rfe) rr$args$n_chains else s1$n_chains %||% 1L

  prefix <- .est_prefix(out)
  if (isTRUE(out$save) && !is.null(mode_res)) {
    dir.create(out$dir, recursive = TRUE, showWarnings = FALSE)
    mode_path <- paste0(prefix, "_mode.rds")
    saveRDS(mode_res, mode_path)
    .vcat(sprintf("  Mode saved -> %s\n", mode_path))
  }

  # -------------------------------------------------------------------
  # Sampling: the sampler stage (run_posterior_estimation()'s code)
  # -------------------------------------------------------------------
  chains   <- NULL
  conv_res <- NULL
  post     <- NULL
  if (!is.null(s1)) {
    .est_set_store(.est_stage_options(spec, s1))
    post <- .est_sampler_stage(.est_sampler_spec(spec, mode_result), mode_result)
    co <- .est_chains_object(post$chains[[1L]], s1$method)
    chains   <- co$chains
    conv_res <- co$convergence
    if (isTRUE(out$save) && !is.null(chains)) {
      chains_path <- paste0(prefix, "_chains.rds")
      saveRDS(chains, chains_path)
      .vcat(sprintf("  Chains saved -> %s\n", chains_path))
    }
    if (!is.null(conv_res$rhat)) {
      .vcat(sprintf("  R-hat: max=%.3f  (>1.05: %d/%d)\n",
                    max(conv_res$rhat),
                    sum(conv_res$rhat > 1.05), length(conv_res$rhat)))
      .vcat(sprintf("  ESS:   median=%.0f  min=%.0f\n",
                    median(conv_res$ess), min(conv_res$ess)))
    }
  } else {
    .vcat("-- Sampling skipped (no sampler) --\n")
    ## A particle likelihood still gets its loglik-variance preflight (at
    ## the mode), as run_full_estimation() always reported it.
    post <- list(tpf_preflight = .est_pmcmc_preflight(
      mode_result$log_post_fn, theta_mode, mode_result$ctx, 1000L, verbose))
  }

  # -------------------------------------------------------------------
  # Ramsey policy workflow (optional)
  # -------------------------------------------------------------------
  ramsey_res <- NULL
  if (isTRUE(out$ramsey)) {
    .vcat("-- Ramsey policy workflow --\n")
    params_mode <- model$param_values
    for (nm in names(theta_mode)) {
      if (nm %in% names(params_mode)) params_mode[[nm]] <- theta_mode[[nm]]
    }
    ## 0.9.4: `ramsey_policy()` ABORTS when it can find no discount
    ## factor. Ramsey is an OPTIONAL extra, so resolve the discount here and
    ## skip the step with a classed warning naming the argument that fixes
    ## it, rather than losing the estimate.
    disc <- if (!is.null(out$ramsey_discount)) as.numeric(out$ramsey_discount)
            else .get_discount(params_mode)

    if (is.null(disc) || length(disc) != 1L || !is.finite(disc)) {
      .dynhr_warn(
        "run_full_estimation(run_ramsey = TRUE): skipping the Ramsey step. ",
        "The planner objective is discounted, but the model has no `beta` or ",
        "`betta` parameter and `ramsey_discount` was not supplied. ",
        "Pass ramsey_discount = <value> to run it; $ramsey stays NULL.",
        class = "dynhr_warn_ramsey_skipped")
    } else {
      ramsey_res <- ramsey_policy(
        model = model,
        compiled = compiled,
        params = params_mode,
        order = out$ramsey_order,
        n_periods = out$ramsey_n_periods,
        burn_in = out$ramsey_burn_in,
        discount = disc,
        verbose = FALSE
      )
    }
  }

  # -------------------------------------------------------------------
  # Diagnostics (optional)
  # -------------------------------------------------------------------
  diag_results <- NULL
  if (isTRUE(out$diagnostics) && !is.null(chains)) {
    .vcat("-- Diagnostics --\n")
    ## Solve at the posterior MODE through the same theta -> params -> steady
    ## state -> decision-rule pipeline the likelihood uses.
    params_diag <- apply_theta_to_params(model, theta_mode)
    ss_diag     <- solve_steady_state(model, compiled, params_diag, verbose = FALSE)
    dr_diag     <- NULL
    irf_diag    <- NULL
    if (!is.null(ss_diag) && isTRUE(ss_diag$converged)) {
      params_diag <- ss_diag$params %||% params_diag
      dr_diag <- solve_perturbation(model, compiled, ss_diag$ss, params_diag,
                                    verbose = FALSE)
      if (!is.null(dr_diag) && isTRUE(dr_diag$bk_satisfied))
        irf_diag <- compute_irfs(dr_diag, model, n_periods = 40L,
                                 params = params_diag)
      else
        .dynhr_warn("run_full_estimation(run_diag = TRUE): the posterior mode ",
                    "does not satisfy Blanchard-Kahn; diagnostics that need a ",
                    "decision rule will report ERROR.",
                    class = "dynhr_warn_diag_mode_solve")
    } else {
      .dynhr_warn("run_full_estimation(run_diag = TRUE): the steady state did ",
                  "not converge at the posterior mode; diagnostics that need ",
                  "a solution will report ERROR.",
                  class = "dynhr_warn_diag_mode_solve")
    }

    ## Provenance for the report: the orchestrator reads these off
    ## `draws` attributes.
    draws_diag <- chains$chain
    attr(draws_diag, "sampler")  <- sampler_lbl
    attr(draws_diag, "n_warmup") <- n_warmup
    attr(draws_diag, "seed")     <- seed

    diag_results <- run_all_diagnostics(
      model       = model,
      compiled    = compiled,
      dr          = dr_diag,
      ss          = if (!is.null(ss_diag)) ss_diag$ss else NULL,
      params      = params_diag,
      priors      = priors,
      data        = data,
      dates       = lik$dates,
      draws       = draws_diag,
      ## One post-warmup draw MATRIX per chain (D5 split R-hat / ESS need the
      ## chains separately; `chains$chain` is their row-bind).
      chains_list = if (!is.null(chains$chain_list))
                      lapply(chains$chain_list, function(ch) ch$chain %||% ch)
                    else list(chains$chain),
      irf         = irf_diag,
      obs_names   = obs_vars,
      theta_mode  = theta_mode,
      model_name  = if (!identical(mod_file, "(pre-parsed)")) basename(mod_file)
                    else out$prefix,
      ramsey_result = ramsey_res,
      verbose     = verbose
    )
    if (!is.null(diag_results))
      .vcat(sprintf("  %d diagnostics run\n", length(diag_results)))
  }

  # -------------------------------------------------------------------
  # OBC smoother + historical decomposition at mode (optional)
  # -------------------------------------------------------------------
  smoother_res <- NULL
  if (isTRUE(out$smoother) && use_obc) {
    .vcat("-- PKF smoother + historical decomposition at mode --\n")
    smoother_res <- .compute_smoother_at_mode(
      theta_mode, model, compiled, data, obs_vars,
      obc_specs, me_variance, lik$obc_max_inner
    )
    if (!is.null(smoother_res))
      .vcat(sprintf("  Smoothed states: %d x %d, Smoothed shocks: %d x %d, HD: %d shocks+constraint\n",
                    nrow(smoother_res$smoothed_states), ncol(smoother_res$smoothed_states),
                    nrow(smoother_res$smoothed_shocks), ncol(smoother_res$smoothed_shocks),
                    length(smoother_res$hd$contributions)))
  }

  # -------------------------------------------------------------------
  # PPF importance re-weighting (PKF adequacy check, opt-in)
  # -------------------------------------------------------------------
  ppf_reweight_res <- NULL
  if (use_obc && !is.null(chains) && isTRUE(out$obc_ppf_reweight)) {
    .vcat("-- PPF importance re-weighting (PKF adequacy check) --\n")
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
        .dynhr_warn("ppf_reweight_posterior failed: ", conditionMessage(e2),
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(ppf_reweight_res)) {
      if (ppf_reweight_res$ess_fraction < 0.5)
        .dynhr_warn(sprintf(
          "PKF may be inadequate for this dataset (ESS/n = %.2f < 0.50). ",
          ppf_reweight_res$ess_fraction),
          "Consider re-estimating with make_log_posterior_obc_ppf ",
          "(proposal='bootstrap' or 'copf').",
          call. = FALSE)
      .vcat(sprintf("  PPF reweight: %s\n", ppf_reweight_res$verdict))
    }
  }

  t_end   <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))
  .vcat(sprintf("  Outputs DONE in %.1f min\n", elapsed))

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
    tpf_preflight   = post$tpf_preflight,
    ## the gradient method(s) the sampler stage resolved (the run record
    ## copies them, .dynhr_rr_set_resolved())
    resolved        = post$resolved,
    ctx             = mode_result$ctx,
    meta            = list(
      mod_file    = mod_file,
      obs_vars    = obs_vars,
      sampler     = sampler_lbl,
      n_draws     = n_draws,
      n_warmup    = n_warmup,
      n_chains    = n_chains,
      seed        = seed,
      use_obc     = use_obc,
      run_ramsey  = isTRUE(out$ramsey),
      started_at  = format(t_start, "%Y-%m-%d %H:%M:%S"),
      elapsed_min = elapsed,
      dynhr_version = utils::packageVersion("dynhr")
    )
  )
  class(result) <- c("dynhr_estimation_result", "list")
  result
}


# ============================================================================
# Checkpoint integrity (plan section 7)
# ============================================================================

## (The result-change registry, .dynhr_result_changes_between() and the
## component tags a spec uses, .est_component_tags(), are in
## R/result-changes.R.)

.est_integrity_file <- function(dir) file.path(dir, "spec_integrity.rds")

## Class A: content hashes of the target / algorithm parts of a spec. The
## sampler hash leaves out n_draws (a resume adds draws on purpose).
.est_target_hashes <- function(spec) {
  sl <- lapply(.spec_sampler_list(spec$sampler), function(s) {
    s <- unclass(s)
    s$n_draws <- NULL
    s
  })
  list(model      = spec$hashes$model,
       data       = spec$hashes$data,
       likelihood = spec$hashes$likelihood,
       mode       = spec$hashes$mode,
       sampler    = .spec_hash(sl))
}

## Classes B (code) and C (numerical environment).
.est_integrity_provenance <- function(spec) {
  pv <- .dynhr_rr_provenance()
  cmp <- spec$compute
  list(version = pv$version, git_commit = pv$git_commit,
       r_version = pv$r_version, platform = pv$platform, os = pv$os,
       blas = pv$blas, lapack = pv$lapack, lapack_version = pv$lapack_version,
       n_cores = if (isTRUE(cmp$parallel)) .mirai_n_cores(cmp$n_cores) else 1L)
}

## Fresh checkpointed run: write the integrity file. Resume: compare it with
## the current spec and build (and act on) the verdict. Returns NULL when
## there is nothing to report, else a `dynhr_integrity` list.
.est_checkpoint_integrity <- function(spec) {
  cmp <- spec$compute
  dir <- cmp$checkpoint_dir
  if (is.null(dir)) return(NULL)
  now <- list(integrity_version = 1L,
              hashes     = .est_target_hashes(spec),
              provenance = .est_integrity_provenance(spec),
              spec       = .dynhr_rr_strip_spec(spec))
  f <- .est_integrity_file(dir)
  if (!isTRUE(cmp$resume)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(now, f)
    return(NULL)
  }
  if (!file.exists(f)) {
    .dynhr_warn("run_estimation: resuming from ", dir, ", which has no ",
                "spec_integrity.rds (a checkpoint written before dynhr's spec ",
                "runner): the target cannot be verified beyond the samplers' ",
                "own parameter check.",
                class = "dynhr_warning_checkpoint_unverified")
    return(structure(list(events = list(unverified = TRUE), checkpoint_dir = dir,
                          overridden = FALSE), class = "dynhr_integrity"))
  }
  was <- readRDS(f)
  events <- list()
  warn_mode <- identical(cmp$on_mismatch, "warn")

  ## ---- A: target / algorithm --------------------------------------------
  parts <- names(now$hashes)
  diff_parts <- parts[!vapply(parts, function(p)
    identical(now$hashes[[p]], was$hashes[[p]]), logical(1))]
  if (length(diff_parts)) {
    fl <- diff_specs(was$spec, now$spec)
    rx <- paste0("^(", paste(diff_parts, collapse = "|"), ")(\\[\\[\\d+\\]\\])?(\\$|$)")
    fl <- fl[grepl(rx, fl$path) & !grepl("\\$n_draws$", fl$path), , drop = FALSE]
    fields <- if (nrow(fl)) sprintf("%s: %s -> %s", fl$path, fl$a, fl$b)
              else character(0)
    msg <- paste0("the spec's ", paste(diff_parts, collapse = ", "),
                  " differ(s) from the checkpointed run",
                  if (length(fields)) paste0(" (", paste(fields, collapse = "; "), ")"),
                  ". The continued chain would target a different posterior or ",
                  "kernel, so the pooled chain would be invalid")
    if (!warn_mode)
      .dynhr_abort("run_estimation: refusing to resume from ", dir, ": ", msg,
                   ". Start a fresh checkpoint directory, or set ",
                   "compute$on_mismatch = \"warn\" to resume anyway (the ",
                   "result is then marked).",
                   class = "dynhr_error_checkpoint_spec_mismatch")
    .dynhr_warn("run_estimation: resuming from ", dir, " although ", msg,
                " (compute$on_mismatch = \"warn\"); the result is marked.",
                class = "dynhr_warning_checkpoint_spec_mismatch")
    events$target <- list(parts = diff_parts, fields = fields)
  }

  ## ---- B: code ----------------------------------------------------------
  pw <- was$provenance
  pn <- now$provenance
  if (!identical(pw$version, pn$version) ||
      !identical(pw$git_commit, pn$git_commit)) {
    ## The continued chain runs the sampler stage only: its target and its
    ## kernel are what a code change must not touch (the mode stage's own
    ## tags are the class-A mode hash's business).
    reg <- .dynhr_result_changes_between(pw$version, pn$version,
                                         .est_component_tags(spec, "sampler"))
    what <- sprintf("dynhr %s%s -> %s%s", format(pw$version),
                    if (!is.na(pw$git_commit %||% NA)) paste0(" (", substr(pw$git_commit, 1L, 12L), ")") else "",
                    format(pn$version),
                    if (!is.na(pn$git_commit %||% NA)) paste0(" (", substr(pn$git_commit, 1L, 12L), ")") else "")
    if (length(reg)) {
      msg <- paste0("the code changed (", what, ") through registered result ",
                    "change(s) touching this run: ", paste(reg, collapse = "; "),
                    ". The kernel or target changed mid-chain")
      if (!warn_mode)
        .dynhr_abort("run_estimation: refusing to resume from ", dir, ": ", msg,
                     ". Install dynhr ", format(pw$version), " (the version ",
                     "that wrote the checkpoint) to continue this chain, start ",
                     "a fresh checkpoint directory, or set ",
                     "compute$on_mismatch = \"warn\" to resume anyway (the ",
                     "result is then marked).",
                     class = "dynhr_error_checkpoint_code_changed")
      .dynhr_warn("run_estimation: resuming from ", dir, " although ", msg,
                  " (compute$on_mismatch = \"warn\"); the result is marked.",
                  class = "dynhr_warning_checkpoint_code_changed")
    } else {
      .dynhr_warn("run_estimation: resuming from ", dir, " with a different ",
                  "build (", what, "). No registered result change touches ",
                  "this run, so the chain stays a valid MCMC chain, but it is ",
                  "not bit-identical to an uninterrupted run (equivalent in ",
                  "distribution: means within MCSE).",
                  class = "dynhr_warning_checkpoint_code_changed")
    }
    events$code <- list(from = pw[c("version", "git_commit")],
                        to = pn[c("version", "git_commit")], registered = reg)
  }

  ## ---- C: numerical environment ------------------------------------------
  env_f <- c("r_version", "platform", "os", "blas", "lapack", "lapack_version",
             "n_cores")
  env_d <- env_f[!vapply(env_f, function(k) identical(pw[[k]], pn[[k]]), logical(1))]
  if (length(env_d)) {
    .dynhr_warn("run_estimation: resuming from ", dir, " in a different ",
                "numerical environment (",
                paste(sprintf("%s: %s -> %s", env_d,
                              vapply(env_d, function(k) paste(format(pw[[k]]), collapse = " "), ""),
                              vapply(env_d, function(k) paste(format(pn[[k]]), collapse = " "), "")),
                      collapse = "; "),
                "). Deterministic quantities agree to floating-point ",
                "precision, but the continued chain is not bit-identical to an ",
                "uninterrupted run (equivalent in distribution). The daemon ",
                "count alone does not change seeded draws.",
                class = "dynhr_warning_checkpoint_environment_changed")
    events$environment <- list(fields = env_d, from = pw[env_d], to = pn[env_d])
  }
  if (!length(events)) return(NULL)
  structure(list(events = events, checkpoint_dir = dir,
                 overridden = !is.null(events$target) ||
                   length(events$code$registered) > 0L),
            class = "dynhr_integrity")
}

## One-line summary of an integrity record (print banners).
.est_integrity_lines <- function(integ) {
  if (is.null(integ)) return(character(0))
  ev <- integ$events
  what <- c(
    if (!is.null(ev$target))
      paste0("resumed across a TARGET change (", paste(ev$target$parts, collapse = ", "),
             "; compute$on_mismatch = \"warn\")"),
    if (!is.null(ev$code))
      paste0("resumed across a code change (", format(ev$code$from$version), " -> ",
             format(ev$code$to$version),
             if (length(ev$code$registered))
               paste0("; ", length(ev$code$registered), " REGISTERED result ",
                      "change(s) touch this run; compute$on_mismatch = \"warn\""),
             ")"),
    if (!is.null(ev$environment))
      paste0("resumed in a different environment (",
             paste(ev$environment$fields, collapse = ", "), ")"),
    if (isTRUE(ev$unverified)) "resumed from an unverifiable checkpoint")
  c(sprintf("  !! INTEGRITY : %s", what[1L]),
    if (length(what) > 1L) sprintf("               %s", what[-1L]),
    if (isTRUE(integ$overridden))
      "               NOT a clean run: see $provenance$integrity")
}
