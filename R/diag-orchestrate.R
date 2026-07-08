## R/diag-orchestrate.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R; updated Phase-3+.
##
## run_all_diagnostics() orchestrator; run_model_diagnostics() convenience wrapper
## --------------------------------------------------------------------------

#' Run all applicable diagnostics
#'
#' Conditionally runs each diagnostic based on which inputs are provided.
#' Returns a named list of dynhr_diagnostic objects with class
#' \code{"dynhr_diagnostic_suite"}.
#'
#' When the first argument is a \code{dynhr_posterior_result} object (from
#' \code{\link{run_posterior_estimation}}), all required inputs are derived
#' automatically and the original \code{report} argument is available for
#' generating output documents.
#'
#' @param model           dynhr_mod or compatible model structure list.
#'   Alternatively, a \code{dynhr_posterior_result} object (all other args
#'   are then derived automatically).
#' @param compiled        Compiled/solved model object.
#' @param dr              Decision rules object.
#' @param ss              Named steady-state vector.
#' @param params          Named calibrated parameter vector.
#' @param priors          Prior specification (format depends on diagnostics).
#' @param data            Data matrix (T x n_obs).
#' @param draws           Posterior draws matrix (n_draws x n_params).
#' @param chains_list     List of draw matrices for multi-chain Gelman-Rubin.
#' @param irf             IRF data (named list of matrices or long data frame).
#' @param vd              Variance decomposition data.
#' @param hd              Historical decomposition (from historical_decomposition()).
#' @param shocks          Smoothed shocks matrix.
#' @param model_solve_fn  Function: theta -> moments vector (for D1, D3, D4).
#' @param prior_draw_fn   Function: () -> theta (for D4).
#' @param prior_density_fn Function: (x, param_name) -> density (for D6).
#' @param obs_names       Observable variable names.
#' @param param_names     Estimated parameter names.
#' @param moment_names    Moment names (for D1, D3, D4).
#' @param param_bounds    Matrix (n_par x 2) of lower/upper bounds (for D3).
#' @param data_moments    Named numeric vector of empirical moments (for D4, D9).
#' @param model_moments   List with \code{$sigma_y}, \code{$acf_y} (for D9).
#' @param models_list     Named list of model results for D14 (Bayes factor).
#' @param results_sub     Subsample draws for D16 stability check.
#' @param mode_multistart Output of run_mode_parallel() for D7 mode robustness.
#' @param dates           Date vector for time-series plots.
#' @param model_name      Human-readable model identifier used in plot titles.
#' @param param_names_c   Character vector of calibrated parameter names (for D26).
#' @param param_names_e   Character vector of estimated parameter names (for D26).
#' @param sigma_e         Shock covariance matrix (for D23 spectral ident).
#' @param abcd_solve_fn   Optional function: \code{theta -> list(A,B,C,D[,Sigma_e])}
#'   giving the ABCD state-space form for the Komunjer & Ng (2011) dynamic
#'   identification rank check (D37). When \code{NULL}, D37 is skipped (the
#'   check differentiates the state space w.r.t. theta, so it needs a
#'   re-solve closure, e.g. wrapping \code{solve_perturbation} +
#'   \code{build_dsge_state_space}).
#' @param ramsey_result   Optional \code{dynhr_ramsey_result2} enabling D24/D25.
#' @param mode_result     Optional \code{dynhr_mode_result} for prior-sensitivity diagnostic.
#' @param obc_specs       List of OBC specs (from \code{obc_parse_tags()}).
#' @param regime_defs     List of regime definitions for D28.
#' @param transition_matrix Square Markov transition matrix for D28.
#' @param d30_prec_bits   Bit precision for D30 (default 128).
#' @param d30_backend     Backend for D30: \code{"auto"}, \code{"Rmpfr"},
#'   \code{"JuliaCall"}, or \code{"base"}.
#' @param d20_weighting   Fisher-information weighting for D20:
#'   \code{"auto"} (default -- use the data sampling-covariance weighting
#'   \eqn{I = J'\,Var(\hat m)^{-1} J} over the var/autocovariance moments when
#'   \code{data}/\code{obs_names}/\code{dr}/\code{model} are available, else fall
#'   back to unweighted \eqn{J'J}), \code{"sampling"} (force the weighted form),
#'   or \code{"none"} (always unweighted). The sampling covariance is built from
#'   the data via a Newey-West HAC estimator divided by T, so \eqn{s_i=\theta/SE}
#'   is a genuine asymptotic t-ratio.
#' @param show_caption    Logical -- embed model name and data hash as a plot
#'   caption (default \code{TRUE}).
#' @param verbose         Logical -- print progress messages (default TRUE).
#' @param report          Report format when first arg is a
#'   \code{dynhr_posterior_result}: \code{"none"} (default, return results),
#'   \code{"html"}, \code{"pdf"}, \code{"md"}.  Ignored for the low-level API.
#' @param report_file     File path stem for the report (without extension).
#' @param output_dir      Directory for report outputs.
#' @param ...             Additional arguments (reserved).
#' @return A named list of \code{dynhr_diagnostic} objects with class
#'   \code{"dynhr_diagnostic_suite"}.
#' @noRd
run_all_diagnostics <- function(model            = NULL,
                                compiled         = NULL,
                                dr               = NULL,
                                dr2              = NULL,
                                ss               = NULL,
                                params           = NULL,
                                priors           = NULL,
                                data             = NULL,
                                draws            = NULL,
                                chains_list      = NULL,
                                irf              = NULL,
                                vd               = NULL,
                                hd               = NULL,
                                shocks           = NULL,
                                model_solve_fn   = NULL,
                                prior_draw_fn    = NULL,
                                prior_density_fn = NULL,
                                obs_names        = NULL,
                                param_names      = NULL,
                                moment_names     = NULL,
                                param_bounds     = NULL,
                                data_moments     = NULL,
                                model_moments    = NULL,
                                mode_multistart  = NULL,
                                draws_by_T       = NULL,
                                kps_sample_sizes = NULL,
                                kps_runner_fn    = NULL,
                                loglik_contrib_fn = NULL,
                                theta_mode       = NULL,
                                loglik_full_fn   = NULL,
                                estimated_names  = NULL,
                                dates            = NULL,
                                model_name       = NULL,
                                show_caption     = TRUE,
                                models_list     = NULL,
                                results_sub     = NULL,
                                param_names_c   = NULL,
                                param_names_e   = NULL,
                                sigma_e         = NULL,
                                abcd_solve_fn   = NULL,
                                solved          = NULL,
                                ramsey_result   = NULL,
                                mode_result     = NULL,
                                obc_specs       = NULL,
                                regime_defs     = NULL,
                                transition_matrix = NULL,
                                d30_prec_bits   = 128L,
                                d30_backend     = c("auto", "Rmpfr", "JuliaCall", "base"),
                                d20_weighting   = c("none", "auto", "sampling"),
                                halt_on_identification_fail = FALSE,
                                verbose         = TRUE,
                                report          = c("none", "html", "pdf", "md"),
                                report_file     = "diagnostic_report",
                                output_dir      = ".",
                                ...) {

  # ---- Auto-dispatch: if first arg is dynhr_posterior_result ----
  if (inherits(model, "dynhr_posterior_result")) {
    return(.run_diagnostics_from_posterior(
      posterior    = model,
      report       = match.arg(report),
      report_file  = report_file,
      output_dir   = output_dir,
      verbose      = verbose,
      halt_on_identification_fail = halt_on_identification_fail,
      ...
    ))
  }

  # Resolve match.arg() defaults BEFORE any lambda captures them.
  # match.arg() only works when called from a function where the argument is
  # a formal — it fails inside anonymous closures (no formals to inspect).
  d30_backend   <- match.arg(d30_backend)
  report        <- match.arg(report)
  d20_weighting <- match.arg(d20_weighting)

  results <- list()

  .msg <- function(txt) if (verbose) message("[dynhr] ", txt)

  # Per-diagnostic error isolation: one crashing diagnostic must not kill the
  # rest of the suite.  Pass the diagnostic call as a zero-arg closure.
  .safe_diag <- function(tag, expr_fn) {
    tryCatch(
      expr_fn(),
      error = function(e) {
        .msg(sprintf("%s FAILED: %s", tag, conditionMessage(e)))
        structure(
          list(status  = "error",
               pass    = NA,
               plots   = list(),
               summary = sprintf("%s failed: %s", tag, conditionMessage(e)),
               errored = TRUE),
          class = "dynhr_diagnostic"
        )
      }
    )
  }

  ## Build provenance descriptor used by all diagnostics
  meta <- diag_meta(model_name   = model_name,
                    data         = data,
                    show_caption = isTRUE(show_caption))

  # --- Expectations (from @dynhr:expectations block) ---
  if (!is.null(model)) {
    .msg("Expectations: checking @dynhr:expectations block...")
    results$expectations <- .safe_diag("Expectations", function()
      diag_expectations(model, data = data, params = params, irfs = irf, meta = meta))
  } else {
    .msg("Expectations: Skipped (model not provided)")
  }

  # --- Data diagnostic (if data provided) ---
  if (!is.null(data)) {
    .msg("Data: plotting observables...")
    results$data <- .safe_diag("Data", function()
      diag_data(data, obs_names = obs_names, dates = dates, meta = meta))
  }

  # --- Model structure summary ---
  .msg("Running model structure summary...")
  dr_for_struct <- if (!is.null(dr)) dr else compiled
  results$model_summary <- .safe_diag("ModelSummary", function()
    model_structure_summary(model, dr_for_struct, ss))

  # --- D0: static equation-system rank (redundant-equation pre-flight) ---
  if (!is.null(compiled) && !is.null(params) && !is.null(ss)) {
    .msg("D0: Static equation-system rank...")
    results$d0 <- .safe_diag("D0", function()
      d0_equation_rank(model, compiled, params, ss))
  }

  # --- D19: Second-order solution quality (Group A-3) ---
  if (!is.null(dr2) && inherits(dr2, "DecisionRules2")) {
    .msg("D19: Second-order solution quality...")
    results$d19 <- .safe_diag("D19", function()
      d19_second_order_accuracy(dr2 = dr2, model = model, ss = ss, params = params))
  } else {
    .msg("D19: Skipped (dr2 not provided; model solved at order < 2)")
  }

  # --- Pre-estimation ---
  # Define params_id at top level (used by D1-D30, not just when model_solve_fn available)
  params_id <- if (!is.null(params)) {
    if (!is.null(param_names) && all(param_names %in% names(params))) {
      params[param_names]
    } else {
      params
    }
  } else NULL

  if (!is.null(model_solve_fn) && !is.null(params)) {
    # D1 and D20 share the same numerical Jacobian of model_solve_fn at
    # params_id; compute it once here (~n_par+1 model solves) instead of twice.
    # If it fails, each diagnostic falls back to computing its own under
    # .safe_diag error isolation.
    J_shared <- tryCatch(.numerical_jacobian(model_solve_fn, params_id),
                         error = function(e) NULL)

    .msg("D1: Local identification...")
    results$d1 <- .safe_diag("D1", function()
      d1_local_identification(model_solve_fn, params_id, param_names, moment_names,
                              jacobian = J_shared, meta = meta))

    # D20 weighting: with data we can build the sampling covariance of the
    # moment estimator (Var(m_hat) = HAC long-run var / T) over the var/acv
    # second moments, and weight the Fisher information by it -- the principled
    # I = J' Var(m_hat)^-1 J on which s_i = theta/SE is a genuine t-ratio.
    # "auto" uses sampling when this is computable, else falls back to the
    # unweighted J'J (which D20 evaluates on model_solve_fn's own moments).
    d20_args <- list(model_solve_fn = model_solve_fn, theta = params_id,
                     param_names = param_names, moment_names = moment_names,
                     jacobian = J_shared, meta = meta)
    if (d20_weighting %in% c("auto", "sampling") &&
        !is.null(data) && !is.null(obs_names) && !is.null(dr) && !is.null(model)) {
      d20_sampling <- tryCatch({
        Yobs    <- as.matrix(data[, obs_names, drop = FALSE])
        T_obs   <- nrow(Yobs)
        dm_names <- names(.compute_data_moments(Yobs, max_lag = 4L))
        # Model-implied var/acv moments as a function of theta (aligns with
        # .compute_data_moments naming); shares its definition with D29.
        dm_solve_fn <- function(theta) {
          pp <- .update_params(params, theta)
          mm <- compute_moments(dr, model, params = pp)
          vc <- mm$var_cov; ac <- mm$autocorr; ml <- dim(ac)[3]
          out <- numeric(0)
          for (v in obs_names) out[paste0("var_", v)] <- vc[v, v]
          for (lag in seq_len(min(ml, 4L)))
            for (v in obs_names) out[paste0("acv", lag, "_", v)] <- ac[v, v, lag] * vc[v, v]
          out[dm_names]
        }
        J_dm  <- .numerical_jacobian(dm_solve_fn, params_id)
        # Var(m_hat): HAC long-run moment covariance divided by T.
        Omega <- .compute_moment_covariance(Yobs, max_lag = 4L, use_hac = TRUE) / T_obs
        dimnames(Omega) <- list(dm_names, dm_names)
        list(model_solve_fn = dm_solve_fn, theta = params_id,
             param_names = param_names, moment_names = dm_names,
             weighting = "sampling", moment_cov = Omega, jacobian = J_dm, meta = meta)
      }, error = function(e) { .msg(sprintf("D20 sampling weighting unavailable: %s", conditionMessage(e))); NULL })
      if (!is.null(d20_sampling)) d20_args <- d20_sampling
      else if (identical(d20_weighting, "sampling"))
        d20_args$weighting <- "sampling"   # let D20 note the fallback
    }
    .msg("D20: Identification strength (Fisher information)...")
    results$d20 <- .safe_diag("D20", function() do.call(d20_fisher_identification_strength, d20_args))

    .msg("D25: Higher-order identification (revisit)...")
    # params_id: estimated params only (matches param_names length).
    # params (full calibration) has more elements than param_names, causing
    # dimnames mismatch in D25's Jacobian column assignment.
    results$d25 <- .safe_diag("D25", function()
      d25_higher_order_identification(
        ramsey_result = ramsey_result, dr = dr, model = model, params = params_id,
        ss = ss, obs_mat = NULL, Sigma_e = sigma_e, param_names = param_names, verbose = verbose))

    if (isTRUE(halt_on_identification_fail) && identical(results$d1$pass, FALSE)) {
      .msg("D1 failed: halting remaining diagnostics (set halt_on_identification_fail=FALSE to continue).")
      class(results) <- "dynhr_diagnostic_suite"
      return(results)
    }
  } else {
    .msg("D1: Skipped (model_solve_fn or params not provided)")
    .msg("D20: Skipped (model_solve_fn or params not provided)")
  }

  # D2 (a re-labelled D23 stub) was removed; D23 below is the real spectral check.
  if (!is.null(dr) && !is.null(params)) {
    .msg("D23: Spectral identification check (Qu-Tkachenko)...")
    results$d23 <- .safe_diag("D23", function()
      d23_spectral_identification(
        dr = dr, model_solve_fn = model_solve_fn, theta = params_id,
        param_names = param_names, Sigma_e = sigma_e))

    .msg("D24: Global identification via KL divergence...")
    results$d24 <- .safe_diag("D24", function()
      d24_global_kl_identification(
        dr = dr, model = model, params = params_id, ramsey_result = ramsey_result,
        Sigma_e = sigma_e, param_names = param_names, verbose = verbose))
  } else {
    .msg("D23: Skipped (dr or params not provided)")
    .msg("D24: Skipped (dr or params not provided)")
  }

  if (!is.null(params) && !is.null(abcd_solve_fn)) {
    .msg("D37: Dynamic identification rank check (Komunjer-Ng)...")
    results$d37 <- .safe_diag("D37", function()
      d37_komunjer_ng(
        dr = dr, model = model, theta = params_id, param_names = param_names,
        model_solve_fn = abcd_solve_fn, Sigma_e = sigma_e, obs_vars = obs_names, meta = meta))
  } else {
    .msg("D37: Skipped (abcd_solve_fn not provided)")
  }

  if (!is.null(model_solve_fn) && !is.null(params) && !is.null(param_bounds)) {
    .msg("D3: Sensitivity analysis (Morris)...")
    results$d3 <- .safe_diag("D3", function()
      d3_sensitivity_morris(model_solve_fn, params_id, param_bounds, param_names, moment_names, meta = meta))
  } else {
    .msg("D3: Skipped (model_solve_fn, params, or param_bounds not provided)")
  }

  if (!is.null(model_solve_fn) && !is.null(params) && !is.null(obs_names)) {
    .msg("D22: Observable informativeness ranking...")
    results$d22 <- .safe_diag("D22", function()
      d22_observable_informativeness(
        model_solve_fn, params_id, obs_names = obs_names,
        param_names = param_names, moment_names = moment_names))
  } else {
    .msg("D22: Skipped (model_solve_fn, params, or obs_names not provided)")
  }

  if (!is.null(model_solve_fn) && !is.null(params) &&
      !is.null(param_names_c) && !is.null(param_names_e)) {
    .msg("D26: Calibration sensitivity diagnostic...")
    # Extract calibrated and estimated parameter subsets
    theta_full <- params
    theta_c <- theta_full[param_names_c]
    theta_e <- theta_full[param_names_e]
    results$d26 <- .safe_diag("D26", function()
      d26_calibration_sensitivity(
        model_solve_fn = model_solve_fn, theta_c = theta_c, theta_e = theta_e,
        param_names_c = param_names_c, param_names_e = param_names_e, meta = meta))
  } else {
    .msg("D26: Skipped (model_solve_fn, params, param_names_c, or param_names_e not provided)")
  }

  if (!is.null(model_solve_fn) && !is.null(prior_draw_fn) && !is.null(data_moments)) {
    .msg("D4: Prior predictive checks...")
    results$d4 <- .safe_diag("D4", function()
      d4_prior_predictive(model_solve_fn, prior_draw_fn, data_moments, moment_names, meta = meta))
  } else {
    .msg("D4: Skipped (model_solve_fn, prior_draw_fn, or data_moments not provided)")
  }

  # --- Estimation ---
  if (!is.null(draws)) {
    .msg("D5: MCMC convergence...")
    results$d5 <- .safe_diag("D5", function()
      d5_mcmc_convergence(draws, param_names, chains_list, meta = meta))
  } else {
    .msg("D5: Skipped (draws not provided)")
  }

  if (!is.null(draws) && !is.null(prior_density_fn)) {
    .msg("D6: Posterior vs prior (ridge density)...")
    results$d6 <- .safe_diag("D6", function()
      d6_posterior_vs_prior(draws, prior_density_fn, param_names, meta = meta))
  } else {
    .msg("D6: Skipped (draws or prior_density_fn not provided)")
  }

  if (!is.null(mode_multistart)) {
    .msg("D7: Mode-finding robustness...")
    results$d7 <- .safe_diag("D7", function()
      d7_mode_robustness(mode_multistart, meta = meta))
  } else {
    .msg("D7: Skipped (mode_multistart not provided)")
  }

  if (!is.null(draws_by_T) || !is.null(kps_runner_fn)) {
    .msg("D21: KPS posterior precision updating...")
    results$d21 <- .safe_diag("D21", function()
      d21_kps_precision_update(
        draws_by_T = draws_by_T, sample_sizes = kps_sample_sizes,
        kps_runner_fn = kps_runner_fn, param_names = param_names, meta = meta))
  } else {
    .msg("D21: Skipped (draws_by_T or kps_runner_fn not provided)")
  }

  # --- Post-estimation ---
  if (!is.null(irf)) {
    .msg("D8: IRF plausibility...")
    results$d8 <- .safe_diag("D8", function() d8_irf_plausibility(irf, metadata = meta))
  } else {
    .msg("D8: Skipped (irf not provided)")
  }

  if (!is.null(model_moments)) {
    .msg("D9: Moment matching...")
    results$d9 <- .safe_diag("D9", function() {
      if (is.list(data_moments) && !is.null(data_moments$sigma_y)) {
        d9_moment_matching(model_moments = model_moments,
                           data_moments = data_moments, obs_names = obs_names)
      } else {
        d9_moment_matching(model_moments = model_moments,
                           data = data, obs_names = obs_names)
      }
    })
  } else {
    .msg("D9: Skipped (model_moments not provided)")
  }

  if (!is.null(vd)) {
    .msg("D10: Variance decomposition...")
    results$d10 <- .safe_diag("D10", function() d10_variance_decomposition(vd))
  } else {
    .msg("D10: Skipped (vd not provided)")
  }

  if (!is.null(hd)) {
    .msg("D11: Historical decomposition...")
    results$d11 <- .safe_diag("D11", function() d11_historical_decomposition(hd, dates = dates, meta = meta))
  } else {
    .msg("D11: Skipped (hd not provided)")
  }

  if (!is.null(shocks)) {
    .msg("D12: Smoothed shocks...")
    results$d12 <- .safe_diag("D12", function() d12_smoothed_shocks(shocks, meta = meta))
  } else {
    .msg("D12: Skipped (shocks not provided)")
  }

  # --- Prior-sensitivity diagnostic (requires mode result) ---
  if (!is.null(solved) && !is.null(data) && !is.null(obs_names) && !is.null(mode_result)) {
    .msg("PriorSensitivity: Comparing informative vs flat priors...")
    me_var <- if (!is.null(data)) diag(cov(data)) * 0.07 else NULL
    results$prior_sensitivity <- .safe_diag("PriorSensitivity", function()
      diag_prior_sensitivity(solved = solved, data = data, obs_vars = obs_names,
                              mode_inf = mode_result, me_variance = me_var,
                              n_iter = 500L, verbose = FALSE))
  } else {
    .msg("PriorSensitivity: Skipped (solved, data, obs_names, or mode_result not provided)")
  }

  if (!is.null(dr) && !is.null(data)) {
    .msg("D13: Cross-equation restrictions...")
    results$d13 <- .safe_diag("D13", function()
      d13_cross_equation_restrictions(
        dr = dr, data = data, obs_names = obs_names,
        sigma_e = sigma_e, model = model, meta = meta))
  } else {
    .msg("D13: Cross-equation restrictions (placeholder -- dr or data not provided)...")
    results$d13 <- .safe_diag("D13", function() d13_cross_equation_restrictions())
  }

  if (!is.null(models_list)) {
    .msg("D14: Bayes factor model comparison...")
    results$d14 <- .safe_diag("D14", function() d14_bayes_factor(models_list, meta = meta))
  } else {
    .msg("D14: Skipped (models_list not provided)")
  }

  if (!is.null(dr) && !is.null(data)) {
    .msg("D15: DSGE-VAR tightness...")
    results$d15 <- .safe_diag("D15", function()
      d15_dsge_var(
        dr        = dr,
        data      = data,
        obs_names = obs_names,
        sigma_e   = sigma_e,
        model     = model,
        meta      = meta
      ))
  } else {
    .msg("D15: DSGE-VAR (placeholder -- dr or data not provided)...")
    results$d15 <- .safe_diag("D15", function() d15_dsge_var())
  }

  if (!is.null(draws) && !is.null(results_sub)) {
    .msg("D16: Subsample stability...")
    results$d16 <- d16_subsample_stability(draws, results_sub, param_names,
                                            meta = meta)
    # D34 reuses D16's regime draws (no extra estimation) and adds the
    # policy/private partition + super-exogeneity (operational Lucas critique).
    .msg("D34: Policy-partitioned invariance (Lucas)...")
    results$d34 <- .safe_diag("D34", function()
      d34_policy_invariance(model = model, draws = draws,
                            results_sub = results_sub, param_names = param_names,
                            meta = meta))
  } else {
    .msg("D16: Skipped (draws or results_sub not provided)")
  }

  .msg("D17: Narrative identification (requires hist_decomp + episodes)...")
  ## D17 requires hist_decomp, smoothed_shocks, and episode specs -- skip here.
  ## Call d17_narrative_identification() directly with full args.

  .msg("D18: Welfare plausibility...")
  results$d18 <- .safe_diag("D18", function()
    d18_welfare_plausibility(ramsey_result = ramsey_result))

  # --- Phase H: Frontier identification ---

  # D27: OBC / Piecewise-Linear Identification
  # Only applicable when the model has OBC/MCP constraints: an explicit
  # obc_specs, or MCP equation tags parseable from the model.
  obc_specs_d27 <- obc_specs
  if (is.null(obc_specs_d27) && !is.null(model)) {
    obc_specs_d27 <- tryCatch(obc_parse_tags(model), error = function(e) NULL)
  }
  if (!is.null(model) && !is.null(params) && length(obc_specs_d27) > 0L) {
    .msg("D27: OBC piecewise-linear identification...")
    results$d27 <- .safe_diag("D27", function()
      d27_obc_identification(model = model, dr = dr, params = params,
                              compiled = compiled, obc_specs = obc_specs_d27,
                              param_names = param_names, verbose = verbose))
  } else {
    .msg("D27: Skipped (model/params not provided, or no OBC/MCP constraints)")
  }

  # D28: Regime-Switching Identification
  if (!is.null(regime_defs) && !is.null(params)) {
    .msg("D28: Regime-switching identification...")
    results$d28 <- .safe_diag("D28", function()
      d28_regime_switching_identification(
        regime_defs = regime_defs, transition_matrix = transition_matrix,
        model = model, dr = dr, params = params, param_names = param_names, verbose = verbose))
  } else {
    .msg("D28: Skipped (regime_defs or params not provided)")
  }

  # D29: Data-Driven Constraints (requires data)
  if (!is.null(params) && !is.null(data) && !is.null(obs_names)) {
    .msg("D29: Data-driven constraints (Lanne-Luoto)...")
    d29_model_solve_fn <- NULL
    if (!is.null(dr) && !is.null(model)) {
      d29_model_solve_fn <- function(theta) {
        pp <- .update_params(params, theta)
        mm <- compute_moments(dr, model, params = pp)
        vc <- mm$var_cov; ac <- mm$autocorr; max_lag <- dim(ac)[3]
        out <- c()
        for (v in obs_names) out[paste0("var_", v)] <- vc[v, v]
        for (lag in seq_len(min(max_lag, 4L)))
          for (v in obs_names) out[paste0("acv", lag, "_", v)] <- ac[v, v, lag] * vc[v, v]
        out
      }
    }
    results$d29 <- .safe_diag("D29", function()
      d29_data_driven_constraints(
        data = data, model_solve_fn = d29_model_solve_fn %||% model_solve_fn,
        theta = params_id, param_names = param_names, verbose = verbose, meta = meta))
  } else {
    .msg("D29: Skipped (model_solve_fn, params, or data not provided)")
  }

  # D30: Arbitrary-Precision Rank Checks
  if (!is.null(model_solve_fn) && !is.null(params)) {
    .msg("D30: Arbitrary-precision rank checks...")
    results$d30 <- .safe_diag("D30", function()
      d30_arbitrary_precision_rank(
        model_solve_fn = model_solve_fn, theta = params_id,
        param_names = param_names, moment_names = moment_names,
        prec_bits = d30_prec_bits, backend = d30_backend, verbose = verbose))
  } else {
    .msg("D30: Skipped (model_solve_fn or params not provided)")
  }

  # --- D31 / D32: OBC scenario comparison + binding summary ---
  if (!is.null(obc_specs) && length(obc_specs) > 0L) {
    .msg("D31: OBC scenario comparison...")
    results$d31 <- .safe_diag("D31", function()
      d31_obc_scenario_comparison(
        sys = if (!is.null(dr)) dr$sys_mat else NULL, dr_slack = dr,
        obc_specs = obc_specs, obs_idx = if (!is.null(dr)) dr$obs_idx else NULL,
        model = model, constraint_name = obc_specs[[1]]$name %||% "OBC"))

    .msg("D32: OBC binding summary...")
    results$d32 <- .safe_diag("D32", function()
      d32_obc_binding_summary(
        sys = if (!is.null(dr)) dr$sys_mat else NULL, dr_slack = dr,
        obc_specs = obc_specs, obs_idx = if (!is.null(dr)) dr$obs_idx else NULL,
        model = model))
  } else {
    .msg("D31/D32: Skipped (obc_specs not provided; non-OBC model)")
  }

  # --- D35: misspecification softness (sandwich / IM-equality) ---
  # Needs a per-period log-likelihood; feeds the Passport's robust axis.
  if (!is.null(loglik_contrib_fn)) {
    th35 <- theta_mode
    if (is.null(th35) && !is.null(draws)) {
      dm <- as.matrix(draws)
      th35 <- stats::setNames(apply(dm, 2, stats::median), colnames(dm))
    }
    if (!is.null(th35)) {
      .msg("D35: Misspecification softness...")
      results$d35 <- .safe_diag("D35", function()
        d35_misspecification_softness(theta = th35,
                                      loglik_contrib_fn = loglik_contrib_fn,
                                      param_names = names(th35), meta = meta))
    }
  }

  # --- D41: KF innovation whiteness (per-observable var/acf1 z-stats) ---
  # Needs a decision rule + params + the filtered data. Skipped gracefully
  # (INFO, not ERROR) when the model's stationary P0 is not valid for this
  # dr/params combination (near-unit-root TT), mirroring kalman_filter's own
  # lik_init = "auto" eigenvalue guard (R/kalman-filter.R ~L845), since the
  # thin diagnostic filter in kf_innovation_diagnostics() only supports
  # lik_init = "stationary" (no exact-diffuse phase implemented there).
  if (!is.null(dr) && !is.null(model) && !is.null(params) && !is.null(data) &&
      !is.null(obs_names)) {
    state_idx_d41 <- dr$state_idx
    tt_d41 <- tryCatch(dr$ghx[state_idx_d41, , drop = FALSE], error = function(e) NULL)
    tt_evals_d41 <- if (!is.null(tt_d41) && nrow(tt_d41) > 0)
      tryCatch(eigen(tt_d41, symmetric = FALSE, only.values = TRUE)$values,
               error = function(e) NULL)
    else NULL
    near_unit_root_d41 <- is.null(tt_evals_d41) || any(Mod(tt_evals_d41) > 1 - 1e-6)

    if (!near_unit_root_d41) {
      .msg("D41: KF innovation whiteness...")
      results$d41 <- .safe_diag("D41", function()
        d41_innovation_whiteness(data, dr = dr, model = model, params = params,
                                 obs_vars = obs_names, meta = meta))
    } else {
      .msg("D41: Skipped (near-unit-root state transition; stationary P0 not valid)")
    }
  } else {
    .msg("D41: Skipped (dr, model, params, data, or obs_names not provided)")
  }

  # --- D36: calibration deepness (is each *fixed* deep param data-consistent?) ---
  # Needs a full-vector log-likelihood; feeds the Passport's calibrated axis.
  if (!is.null(loglik_full_fn) && !is.null(params)) {
    .msg("D36: Calibration deepness...")
    results$d36 <- .safe_diag("D36", function()
      d36_calibration_deepness(model = model, params = params,
                               loglik_fn = loglik_full_fn,
                               estimated = estimated_names, meta = meta))
  }

  # --- D33: structural-vs-reduced-form ("borrowed identification") ---
  # Needs only the model + @dynhr:deep map; feeds the Passport's structural axis.
  if (!is.null(model)) {
    .msg("D33: structural-vs-reduced-form...")
    results$d33 <- .safe_diag("D33", function()
      d33_structural_vs_reduced_form(model = model, meta = meta))
  }

  # --- Deep-Parameter Passport: synthesis over the identification block ---
  # Runs last so it can read the other diagnostics' results out of `results`.
  if (!is.null(model)) {
    .msg("Deep-Parameter Passport...")
    results$deep_passport <- .safe_diag("Deep-Parameter Passport", function()
      deep_parameter_passport(
        deep_spec = build_deep_spec(model = model),
        suite     = results,
        draws     = draws,
        meta      = meta))
  }

  # Sort results by group (A < B < C < D < ?) then importance rank (display
  # order only — does not rename or re-run any diagnostic function).
  .group_ord <- c(A = 1L, B = 2L, C = 3L, D = 4L, "?" = 5L)
  diag_nms   <- names(results)[vapply(results,
    function(r) inherits(r, "dynhr_diagnostic"), logical(1))]
  non_diag   <- setdiff(names(results), diag_nms)
  if (length(diag_nms) > 0L) {
    meta_list <- lapply(diag_nms, .diag_meta_for)
    sort_key  <- order(
      vapply(meta_list, function(m) .group_ord[[m$group]] %||% 5L, integer(1)),
      vapply(meta_list, function(m) m$importance,                    integer(1))
    )
    results <- c(results[non_diag], results[diag_nms[sort_key]])
  }

  # Set class
  class(results) <- "dynhr_diagnostic_suite"

  # Print summary
  if (verbose) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("dynhr diagnostic suite -- summary")
    message(paste(rep("=", 70), collapse = ""))

    for (nm in names(results)) {
      r <- results[[nm]]
      if (inherits(r, "dynhr_diagnostic")) {
        badge <- .badge_str(r)
        # One line per diagnostic: collapse embedded newlines, then truncate
        flat_summary  <- gsub("[[:space:]]+", " ", r$summary)
        short_summary <- substr(flat_summary, 1, 100)
        if (nchar(flat_summary) > 100) short_summary <- paste0(short_summary, "...")
        message(sprintf("  [%-5s] %s: %s", badge, nm, short_summary))
      }
    }

    n_pass <- sum(sapply(results, function(r) isTRUE(r$pass)), na.rm = TRUE)
    n_fail <- sum(sapply(results, function(r) identical(r$pass, FALSE)), na.rm = TRUE)
    n_err  <- sum(sapply(results, function(r) isTRUE(r$errored)), na.rm = TRUE)
    n_info <- sum(sapply(results, function(r) is.na(r$pass) && !isTRUE(r$errored)), na.rm = TRUE)
    message(sprintf("\nTotal: %d PASS, %d FAIL, %d ERROR, %d INFO",
                    n_pass, n_fail, n_err, n_info))
    message(paste(rep("=", 70), collapse = ""))
  }

  results
}


# ---------------------------------------------------------------------------
#' D41. Kalman-filter innovation whiteness
#'
#' Thin orchestrator wrapper around \code{\link{kf_innovation_diagnostics}}:
#' runs the per-observable standardized-innovation whiteness check on
#' \code{data} at \code{(dr, params)} and reports per-observable sample
#' variance / lag-1 acf z-stats, flagging \code{|z| > 4} (see
#' \code{kf_innovation_diagnostics}'s own conservative multiple-comparisons
#' threshold) as a warning-level (FAIL badge) finding.
#'
#' @param data      Observation matrix (\code{T x n_obs} or \code{n_obs x T};
#'   orientation is resolved against \code{obs_vars} the same way
#'   \code{\link{kf_innovation_diagnostics}} does).
#' @param dr        Decision rule (output of \code{\link{solve_perturbation}}).
#' @param model     Compiled/parsed model object.
#' @param params    Named numeric parameter vector.
#' @param obs_vars  Character vector of observed variable names.
#' @param me_variance Scalar measurement-error jitter (default 0), forwarded
#'   to \code{\link{kf_innovation_diagnostics}}.
#' @param meta      Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic}; \code{$result} is the raw
#'   \code{kf_innovation_diagnostics} object.
#' @noRd
d41_innovation_whiteness <- function(data, dr, model, params, obs_vars,
                                     me_variance = 0, meta = NULL) {
  tryCatch({
    ## data may be T x n_obs (the orchestrator's convention) or n_obs x T
    ## (kf_innovation_diagnostics's convention); match against obs_vars
    ## dimnames when available, else fall back on the n_obs x T assumption
    ## that kf_innovation_diagnostics itself uses internally.
    Y <- as.matrix(data)
    n_obs <- length(obs_vars)
    if (nrow(Y) != n_obs && ncol(Y) == n_obs) Y <- t(Y)

    diag_obj <- kf_innovation_diagnostics(Y, dr = dr, model = model,
                                          params = params, obs_vars = obs_vars,
                                          lik_init = "stationary",
                                          me_variance = me_variance)

    by_obs <- diag_obj$by_obs
    flagged <- by_obs$obs_var[which(abs(by_obs$z_var) > 4 | abs(by_obs$z_acf1) > 4)]
    pass <- length(flagged) == 0

    summary_txt <- sprintf(
      "D41 KF innovation whiteness: %d observable(s), max |z| = %.2f, %d flagged (|z| > 4)%s.",
      nrow(by_obs), diag_obj$joint$max_abs_z, diag_obj$joint$n_flagged,
      if (length(flagged)) paste0(" [", paste(flagged, collapse = ", "), "]") else "")

    badge <- if (pass) "PASS" else "FAIL"
    llm <- paste(c(
      sprintf("D41 | KF Innovation Whiteness | %s", badge),
      sprintf("  n_obs=%d max_abs_z=%.3f n_flagged=%d",
              nrow(by_obs), diag_obj$joint$max_abs_z, diag_obj$joint$n_flagged),
      sprintf("  per_obs: %s",
              paste(sprintf("%s var_z=%.3f(z=%.2f) acf1=%.3f(z=%.2f)",
                            by_obs$obs_var, by_obs$var_z, by_obs$z_var,
                            by_obs$acf1, by_obs$z_acf1), collapse = "; ")),
      if (length(flagged))
        sprintf("  flagged(|z|>4): %s", paste(flagged, collapse = ", ")),
      sprintf("  action: %s",
              if (pass)
                "Standardized innovations are consistent with white noise; no evidence of a likelihood-evaluation bug from this oracle."
              else "Innovations deviate from white noise (variance != 1 or nonzero lag-1 acf); check Sigma_e, ZZ/DD routing, steady state, and per-period tunes (me_extra/shock_scale) at this evaluation point.")
    ), collapse = "\n")

    .make_result(
      result = diag_obj, pass = pass, plots = list(),
      summary = summary_txt, llm_summary = llm)
  }, error = function(e) {
    .make_result(pass = NA, errored = TRUE,
                 summary = paste("D41 KF innovation whiteness: ERROR --",
                                 conditionMessage(e)))
  })
}


#' Convenience wrapper: run D8/D9/D10 on a StochSimulResult
#' @noRd
run_model_diagnostics <- function(model, dr = NULL, ss = NULL, verbose = TRUE) {
  results <- list()

  if (verbose) cat("\n-- Running model_structure_summary...\n")
  results$structure <- model_structure_summary(model, dr, ss)
  if (verbose) cat(results$structure$summary, "\n")

  irfs <- if (!is.null(dr$irfs)) dr$irfs else if (!is.null(dr$irf)) dr$irf else NULL
  if (!is.null(irfs)) {
    if (verbose) cat("\n-- Running d8_irf_plausibility...\n")
    results$d8 <- d8_irf_plausibility(irfs)
    if (verbose) cat(results$d8$summary, "\n")
  } else {
    if (verbose) cat("\n-- D8 skipped (no IRFs)\n")
  }

  if (!is.null(dr$moments)) {
    if (verbose) cat("\n-- Running d9_moment_matching...\n")
    results$d9 <- d9_moment_matching(dr)
    if (verbose) cat(results$d9$summary, "\n")
  } else {
    if (verbose) cat("\n-- D9 skipped (no moments)\n")
  }

  if (!is.null(dr$moments$var_decomp_pct)) {
    if (verbose) cat("\n-- Running d10_variance_decomposition...\n")
    results$d10 <- d10_variance_decomposition(dr)
    if (verbose) cat(results$d10$summary, "\n")
  } else {
    if (verbose) cat("\n-- D10 skipped (no variance decomposition)\n")
  }

  if (verbose) {
    cat("\n-- Diagnostic summary -----------------------------\n")
    for (nm in names(results)) {
      r <- results[[nm]]
      status <- if (is.na(r$pass)) "N/A" else if (r$pass) "PASS" else "FAIL"
      cat(sprintf("  %-20s  %s\n", nm, status))
    }
    cat("\n")
  }

  invisible(results)
}


# ============================================================================
# Helper: auto-derive diagnostics from dynhr_posterior_result
# ============================================================================

#' Run diagnostics from a dynhr_posterior_result object
#'
#' Derives all required inputs from the posterior result object and calls
#' \code{run_all_diagnostics()} with the expanded arguments.  Supports
#' optional report generation (HTML, PDF, or Markdown).
#'
#' @param posterior   A \code{dynhr_posterior_result} object.
#' @param report      Report format: \code{"none"} (default), \code{"html"},
#'   \code{"pdf"}, \code{"md"}.
#' @param report_file File path stem (without extension).
#' @param output_dir  Output directory.
#' @param verbose     Print progress.
#' @param ...         Additional arguments forwarded to \code{run_all_diagnostics()}.
#' @return A \code{dynhr_diagnostic_suite} object, invisibly.  When
#'   \code{report != "none"}, also writes the report file.
#' @noRd
.run_diagnostics_from_posterior <- function(posterior,
                                             report       = c("none", "html", "pdf", "md"),
                                             report_file  = "diagnostic_report",
                                             output_dir   = ".",
                                             verbose      = TRUE,
                                             halt_on_identification_fail = FALSE,
                                             bayesian_irf = FALSE,
                                             ...) {
  report <- match.arg(report)

  mode_res <- posterior$mode_result
  solved   <- mode_res$solved

  model      <- solved$model
  compiled   <- solved$compiled
  dr         <- solved$dr
  ss         <- if (!is.null(solved$ss)) solved$ss$values else NULL
  params     <- posterior$posterior_mean
  priors     <- mode_res$prior_spec
  data       <- mode_res$data
  obs_names  <- mode_res$obs_vars
  param_names <- priors$name
  draws      <- posterior$pooled_draws
  irf        <- posterior$posterior_irfs
  model_moments <- posterior$posterior_moments

  ## ---- D10 / D11 / D12 inputs (previously unreachable from a posterior) ----
  ## D10: variance decomposition, already computed at the posterior mean.
  vd <- if (!is.null(model_moments)) model_moments$var_decomp_pct else NULL

  ## D11/D12: run the Kalman smoother at the posterior mean to recover the
  ## smoothed shocks and the historical decomposition. Gaussian-linear only; for
  ## OBC / non-Gaussian likelihoods the linear smoother does not apply, so these
  ## stay NULL and D11/D12 are skipped exactly as before (graceful degradation).
  hd <- NULL; smoothed_shocks <- NULL
  use_obc_post <- !is.null(mode_res$obc_specs)
  if (!use_obc_post && !is.null(data) && !is.null(obs_names) &&
      !is.null(model) && !is.null(dr)) {
    sp_sm <- tryCatch(build_dsge_state_space(model, dr, obs_names, verbose = FALSE),
                      error = function(e) NULL)
    if (!is.null(sp_sm)) {
      ## data is T x n_obs (same orientation cov(data)/colMeans(data) use); the
      ## smoother takes Y as T x n_obs (nrow = T).
      Y_sm <- as.matrix(data[, obs_names, drop = FALSE])
      sm_out <- tryCatch(kalman_smoother(Y_sm, sp_sm), error = function(e) NULL)
      if (!is.null(sm_out)) {
        smoothed_shocks <- sm_out$smoothed_shocks
        hd <- tryCatch(historical_decomposition(smoothed_shocks, sp_sm),
                       error = function(e) NULL)
      }
    }
  }

  ## Posterior IRF credible bands (re-solves the model at many draws -- minutes;
  ## opt-in via bayesian_irf = TRUE). Stored as results$bayesian_irf below.
  bayesian_irf_result <- NULL
  if (isTRUE(bayesian_irf) && !is.null(draws) && !is.null(model) &&
      !is.null(compiled) && nrow(as.matrix(draws)) >= 10L) {
    bayesian_irf_result <- tryCatch(
      diag_bayesian_irf(model, compiled, draws), error = function(e) NULL)
  }

  # Build a model_solve_fn from the compiled model
  model_solve_fn <- NULL
  if (!is.null(compiled) && !is.null(dr)) {
    model_solve_fn <- function(theta) {
      pp <- .update_params(params, theta)
      moments <- compute_moments(dr, model, params = pp)
      sd_named <- moments$std_dev
      names(sd_named) <- paste0("sd_", names(sd_named))
      vd <- moments$var_decomp_pct
      vd_named <- as.numeric(vd)
      names(vd_named) <- paste0("vd_", 
                                 rep(rownames(vd), ncol(vd)), "_",
                                 rep(colnames(vd), each = nrow(vd)))
      c(sd_named, vd_named)
    }
  }

  # Prior density function for D6
  prior_density_fn <- NULL
  if (!is.null(priors)) {
    prior_density_fn <- function(x, param_name) {
      idx <- match(param_name, priors$name)
      if (is.na(idx)) return(dnorm(x, 0, 1))
      exp(log_prior_density(x, priors$distribution[idx],
                            priors$p1[idx], priors$p2[idx],
                            priors$lower[idx], priors$upper[idx]))
    }
  }

  # Prior draw function for D4
  prior_draw_fn <- NULL
  if (!is.null(priors)) {
    prior_sampler <- .smc_make_prior_sampler(priors)
    prior_draw_fn <- function() {
      setNames(vapply(names(prior_sampler),
                       function(nm) prior_sampler[[nm]](), numeric(1)),
               names(prior_sampler))
    }
  }

  # Data moments for D4, D9
  data_moments <- NULL
  if (!is.null(data)) {
    data_moments <- list(
      sigma_y = cov(data),
      mean_y  = colMeans(data, na.rm = TRUE)
    )
  }

  # Chains list for multi-chain convergence
  chains_list <- NULL
  if (!is.null(posterior$chains)) {
    chains_list <- list()
    for (m in names(posterior$chains)) {
      if (!is.null(posterior$chains[[m]]$chains)) {
        for (ch in posterior$chains[[m]]$chains) {
          if (!is.null(ch$chain))
            chains_list <- c(chains_list, list(ch$chain))
        }
      }
    }
    if (length(chains_list) == 0) chains_list <- NULL
  }

  # Determine model_name from file path
  model_name <- NULL
  if (!is.null(model$source_file) && !is.na(model$source_file)) {
    model_name <- tools::file_path_sans_ext(basename(model$source_file))
  }

  # Collect Ramsey result if present
  ramsey_result <- solved$ramsey

  # Run diagnostics with all derived args
  results <- run_all_diagnostics(
    model            = model,
    compiled         = compiled,
    dr               = dr,
    ss               = ss,
    params           = params,
    priors           = priors,
    data             = data,
    draws            = draws,
    chains_list      = chains_list,
    irf              = irf,
    vd               = vd,
    hd               = hd,
    shocks           = smoothed_shocks,
    model_moments    = model_moments,
    model_solve_fn   = model_solve_fn,
    prior_draw_fn    = prior_draw_fn,
    prior_density_fn = prior_density_fn,
    obs_names        = obs_names,
    param_names      = param_names,
    data_moments     = data_moments,
    solved           = solved,
    mode_result      = mode_res,
    ramsey_result    = ramsey_result,
    obc_specs        = mode_res$obc_specs,
    verbose          = verbose,
    halt_on_identification_fail = halt_on_identification_fail,
    ...
  )

  ## Attach posterior IRF credible bands (opt-in; NULL otherwise).
  if (!is.null(bayesian_irf_result))
    results$bayesian_irf <- bayesian_irf_result

  # ---- Report generation ----
  if (report != "none") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (report == "md") {
      report_path <- file.path(output_dir, paste0(report_file, ".md"))
      .write_diagnostic_md(results, report_path, model_name = model_name)
      if (verbose) message("[dynhr] Report written -> ", report_path)
    } else {
      # HTML / PDF: render via Quarto (report.qmd / report-pdf.qmd).
      # The templates read params$diags_rds (file path), so we write a temp RDS.
      qmd_name <- if (report == "pdf") "report-pdf.qmd" else "report.qmd"
      qmd_template <- system.file("templates", qmd_name, package = "dynhr")
      quarto_ok    <- requireNamespace("quarto",    quietly = TRUE)
      template_ok  <- nzchar(qmd_template) && file.exists(qmd_template)

      if (quarto_ok && template_ok) {
        # Write results to a temp RDS so the template can read it
        tmp_rds  <- tempfile(fileext = ".rds")
        on.exit(unlink(tmp_rds), add = TRUE)
        saveRDS(results, tmp_rds)

        qmd_out <- file.path(output_dir, paste0(report_file, ".qmd"))
        file.copy(qmd_template, qmd_out, overwrite = TRUE)
        quarto::quarto_render(
          input       = qmd_out,
          execute_params = list(
            diags_rds  = tmp_rds,
            model_name = model_name %||% ""
          ),
          quiet = !verbose
        )
        if (verbose)
          message("[dynhr] Quarto report -> ",
                  file.path(output_dir,
                            paste0(report_file,
                                   if (report == "pdf") ".pdf" else ".html")))
      } else {
        # Fallback: write markdown and warn clearly
        md_path <- file.path(output_dir, paste0(report_file, ".md"))
        .write_diagnostic_md(results, md_path, model_name = model_name)
        if (!quarto_ok) {
          warning(sprintf(
            "[dynhr] quarto package not installed; wrote Markdown to %s instead.",
            md_path))
        } else {
          warning(sprintf(
            "[dynhr] Quarto template '%s' not found; wrote Markdown to %s instead.",
            qmd_name, md_path))
        }
      }
    }
  }

  # ---- Executive summary ---------------------------------------------------
  if (report != "none") {
    summary_path <- file.path(output_dir, "executive_summary.md")
    write_executive_summary(results, file = summary_path, model_name = model_name)
    if (verbose) message("[dynhr] Executive summary -> ", summary_path)
  }

  invisible(results)
}


#' Write a diagnostic summary as Markdown
#' @noRd
.write_diagnostic_md <- function(results, file_path, model_name = NULL) {
  lines <- c()
  lines <- c(lines, "# dynhr Diagnostic Report")
  if (!is.null(model_name))
    lines <- c(lines, paste0("**Model:** ", model_name))
  lines <- c(lines, paste0("**Generated:** ", Sys.time()))
  lines <- c(lines, "")
  lines <- c(lines, "## Diagnostic Summary")
  lines <- c(lines, "")
  lines <- c(lines, "| Diagnostic | Status | Summary |")
  lines <- c(lines, "|---|---|---|")

  for (nm in names(results)) {
    r <- results[[nm]]
    if (!inherits(r, "dynhr_diagnostic")) next
    status <- if (isTRUE(r$errored))      "ERROR"
              else if (isTRUE(r$pass))    "PASS"
              else if (isFALSE(r$pass))   "FAIL"
              else                        "N/A"
    raw_sum <- if (!is.null(r$summary)) gsub("\n", " ", r$summary, fixed = TRUE) else ""
    safe_sum <- gsub("|", "\\|", raw_sum, fixed = TRUE)
    safe_nm  <- gsub("|", "\\|", nm, fixed = TRUE)
    lines <- c(lines, sprintf("| %s | %s | %s |", safe_nm, status, safe_sum))
  }
  lines <- c(lines, "")

  writeLines(lines, file_path)
}

