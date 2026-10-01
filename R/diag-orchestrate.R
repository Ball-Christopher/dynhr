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
#' @section Callback convention -- signal failure by RETURN VALUE, never by
#'   \code{stop()}:
#'
#' Every user-supplied callback here (\code{model_solve_fn},
#' \code{loglik_contrib_fn}, \code{loglik_full_fn}, \code{abcd_solve_fn},
#' \code{kps_runner_fn}, \code{d30_hp_solve_fn}, \code{prior_density_fn})
#' is called at MANY parameter vectors, some of which will be infeasible --
#' no steady state, a Blanchard-Kahn violation, a non-stationary draw, a
#' parameter outside its support. The contract is:
#'
#' \itemize{
#'   \item A log-likelihood / log-posterior callback returns \code{-Inf}.
#'   \item A moment / statistic callback returns \code{NA_real_} (a vector of
#'     \code{NA_real_} of the SAME LENGTH AND NAMES as a successful return --
#'     a short or unnamed vector is treated as a different moment set).
#'   \item Nothing calls \code{stop()} for an infeasible parameter vector.
#' }
#'
#' Reserve \code{stop()} for genuine programming errors (wrong argument type,
#' a missing model object) -- those SHOULD propagate. This mirrors Dynare's
#' \code{resol.m}, which returns integer \code{info} codes (20 = no steady
#' state, 21 = complex steady state, 30 = ergodic variance failure) rather
#' than throwing, precisely so the calling estimation routine can score the
#' draw and continue; Stan's \code{reject()} has the same
#' reject-this-state-as-\code{-Inf} semantics.
#'
#' A D3/D4 moment callback therefore looks like:
#' \preformatted{
#' moment_names <- names(my_moments(theta0))
#' na_out <- setNames(rep(NA_real_, length(moment_names)), moment_names)
#' model_solve_fn <- function(theta) {
#'   sol <- try_solve(theta)          # your own solve, returning NULL on failure
#'   if (is.null(sol)) return(na_out) # NOT stop("no steady state")
#'   my_moments(sol)
#' }
#' }
#' The diagnostics then report a non-finite Jacobian row / a rejected draw,
#' which is informative, instead of aborting the suite.
#'
#' One caveat for code OUTSIDE dynhr: general R optimisers are far less
#' tolerant of non-finite objectives than dynhr's own guarded call sites
#' (\code{nloptr} errors outright on \code{NA}). If you reuse one of these
#' callbacks under a raw \code{optim()}/\code{nloptr()} call, convert the
#' non-finite return into a large finite penalty at THAT boundary; dynhr's
#' internals already guard with \code{is.finite()}.
#'
#' @param model           dynhr_mod or compatible model structure list.
#'   Alternatively, a \code{dynhr_posterior_result} object (all other args
#'   are then derived automatically).
#' @param compiled        Compiled/solved model object.
#' @param dr              Decision rules object.
#' @param dr2             Optional second-order \code{DecisionRules2}; enables D19
#'   and D25 (D25 re-solves the model at order 2 via \code{model} + \code{compiled}).
#' @param ss              Named steady-state vector.
#' @param params          Named calibrated parameter vector.
#' @param priors          Prior specification (format depends on diagnostics).
#' @param data            Data matrix (T x n_obs).
#' @param draws           Posterior draws matrix (n_draws x n_params).
#' @param chains_list     List of draw matrices for multi-chain Gelman-Rubin.
#' @param irf             IRF data (named list of matrices or long data frame).
#' @param vd              Variance decomposition data.
#' @param hd              Historical decomposition (from historical_decomposition()).
#' @param shocks          Smoothed shocks: T x n matrix (columns = shock names), or a list with \code{smoothed_shocks} and \code{shock_sd} (per-period sd of the smoothed shock, used by D12 to standardise).
#' @param model_solve_fn  Function: theta -> moments vector (for D1, D3, D4).
#' @param prior_draw_fn   Function: () -> theta (for D4).
#' @param prior_density_fn Function: (x, param_name) -> PROPER prior density
#'   (for D6; not renormalised there). When NULL and \code{priors} is a
#'   prior_spec data.frame, it is built from \code{priors} with the
#'   estimation density \code{log_prior()}, normalised over [lower, upper].
#' @param loglik_contrib_fn Optional function: \code{theta -> numeric(T)} of
#'   PER-PERIOD log-likelihood contributions at \code{theta} (used by D21's
#'   precision-updating path and by the sample-split checks). Must follow the
#'   callback convention above: return \code{rep(-Inf, T)} -- or a length-T
#'   vector of \code{NA_real_} -- for an infeasible \code{theta}, never
#'   \code{stop()}.
#' @param theta_mode      Optional named numeric vector: the mode / point
#'   estimate at which the \code{loglik_*} callbacks are anchored. Defaults to
#'   \code{params} when NULL. Supply it when \code{params} is the posterior
#'   mean but the likelihood-based diagnostics should be evaluated at the mode.
#' @param loglik_full_fn  Optional function: a FULL named parameter vector ->
#'   scalar log-likelihood, re-solving the model at the vector it receives
#'   (a closure over a fixed decision rule makes every D36 profile flat).
#'   Enables D36. Same callback convention: non-finite return, no
#'   \code{stop()}.
#' @param estimated_names Optional character vector of estimated parameter
#'   names passed to D36 as its nuisance block (excluded from D36's default
#'   calibrated targets).
#' @param obs_names       Observable variable names.
#' @param param_names     Estimated parameter names.
#' @param moment_names    Moment names (for D1, D3, D4).
#' @param param_bounds    Matrix (n_par x 2) of lower/upper bounds (for D3).
#' @param data_moments    Named numeric vector of empirical moments (for D4, D9).
#' @param model_moments   List with \code{$sigma_y}, \code{$acf_y} (for D9).
#' @param models_list     Named list of model results for D14 (Bayes factor).
#' @param log_marglik_se  Optional named numeric vector of Monte Carlo standard
#'   errors of the log marginal likelihoods in \code{models_list} (SMC /
#'   bridge estimates); forwarded to D14, which reports INFO when the
#'   best-versus-runner-up gap is within \code{2} combined SEs.
#' @param results_sub     Subsample draws for D16 stability check.
#' @param mode_multistart Output of run_mode_parallel() for D7 mode robustness.
#' @param dates           Date vector for time-series plots.
#' @param model_name      Human-readable model identifier used in plot titles.
#' @param param_names_c   Character vector of calibrated parameter names (for D26).
#' @param param_names_e   Character vector of estimated parameter names (for D26).
#'   D26 uses the re-solving var/acv moments of \code{obs_names} (weighted by
#'   their sampling precision when \code{data} is given) when \code{model},
#'   \code{compiled} and \code{obs_names} are available, else
#'   \code{model_solve_fn}.
#' @param sigma_e         Shock covariance matrix (for D23 spectral ident).
#' @param abcd_solve_fn   Optional function: \code{theta -> list(A,B,C,D[,Sigma_e])}
#'   giving the ABCD state-space form for the Komunjer & Ng (2011) dynamic
#'   identification rank check (D37). When \code{NULL}, D37 is skipped (the
#'   check differentiates the state space w.r.t. theta, so it needs a
#'   re-solve closure, e.g. wrapping \code{solve_perturbation} +
#'   \code{build_dsge_state_space}).
#' @param ramsey_result   Optional \code{dynhr_ramsey_result2} (used by D18; D24/D25
#'   cannot re-solve a Ramsey problem and then report INFO).
#' @param mode_result     Optional \code{dynhr_mode_result} for prior-sensitivity diagnostic.
#' @param obc_specs       List of OBC specs (from \code{obc_parse_tags()}).
#' @param regime_defs     List of regime definitions for D28.
#' @param transition_matrix Square Markov transition matrix for D28.
#' @param d30_prec_bits   Bit precision for D30 (default 128).
#' @param d30_backend     Backend for D30: \code{"auto"}, \code{"Rmpfr"},
#'   or \code{"base"} (double baseline only).
#' @param d30_hp_solve_fn Optional mpfr-generic version of
#'   \code{model_solve_fn} for D30's end-to-end high-precision rank check.
#'   \code{NULL} (default) reports the double baseline only (pass = NA): a
#'   Jacobian computed in double gains nothing from mpfr.
#' @param d20_weighting   Fisher-information weighting for D20:
#'   \code{"none"} (default -- unweighted \eqn{J'J} on \code{model_solve_fn}'s
#'   moments; \eqn{s_i} then assumes unit moment variances and is not a
#'   t-ratio), \code{"auto"} (use the data sampling-covariance weighting
#'   \eqn{I = J'\,Var(\hat m)^{-1} J} over the var/autocovariance moments when
#'   \code{data}/\code{obs_names}/\code{dr}/\code{model} are available, else fall
#'   back to \code{"none"}), or \code{"sampling"} (as \code{"auto"}, noting the
#'   fallback when it is not computable). The sampling covariance is a Bartlett
#'   (Newey-West) HAC estimate with Andrews' plug-in bandwidth, divided by T, so
#'   \eqn{s_i=\theta/SE} is an asymptotic t-ratio.
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
                                log_marglik_se  = NULL,
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
                                d30_backend     = c("auto", "Rmpfr", "base"),
                                d30_hp_solve_fn = NULL,
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

  .msg <- function(txt) if (verbose) .dynhr_inform("[dynhr] ", txt)

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
  # `dr` lets D0 tell a unit-root singularity (WARN, expected) from a genuinely
  # redundant equation (FAIL) the way Dynare's model_diagnostics.m does.
  if (!is.null(compiled) && !is.null(params) && !is.null(ss)) {
    .msg("D0: Static equation-system rank...")
    results$d0 <- .safe_diag("D0", function()
      d0_equation_rank(model, compiled, params, ss, dr = dr))
  }

  # --- D40: near-unit root on the state transition -------------------------
  # `compiled` is deliberately NOT passed: the parameter sensitivities are two
  # extra model solves PER parameter, which does not belong in the default
  # suite. Call diag_near_unit_root() directly for those.
  if (!is.null(dr)) {
    .msg("D40: Near-unit root...")
    results$d40 <- .safe_diag("D40", function()
      diag_near_unit_root(model, compiled = NULL, theta = params, dr = dr,
                          n_obs = if (!is.null(data)) nrow(data) else NULL,
                          param_names = character(0), meta = meta))
  } else {
    .msg("D40: Skipped (dr not provided)")
  }

  # --- D19: Second-order solution quality (Group A-3) ---
  if (!is.null(dr2) && inherits(dr2, "DecisionRules2")) {
    .msg("D19: Second-order solution quality...")
    results$d19 <- .safe_diag("D19", function()
      d19_second_order_accuracy(dr2 = dr2, model = model, ss = ss,
                                params = params, meta = meta))
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
                         error = function(e) .dynhr_reraise_bug(e, NULL))

    .msg("D1: Local identification...")
    results$d1 <- .safe_diag("D1", function()
      d1_local_identification(model_solve_fn, params_id, param_names, moment_names,
                              jacobian = J_shared, meta = meta))

    # D20 weighting: with data we can build the sampling covariance of the
    # moment estimator (Var(m_hat) = HAC long-run var / T) over the var/acv
    # second moments, and weight the Fisher information by it -- the principled
    # I = J' Var(m_hat)^-1 J on which s_i = theta/SE is a genuine t-ratio.
    # The model side RE-SOLVES the model at each theta (needs `compiled`);
    # until 0.9.4 it reused the fixed `dr`, so only shock-covariance
    # parameters moved the moments and every structural parameter looked
    # unidentified. "auto" uses sampling when this is computable, else the
    # unweighted J'J on model_solve_fn's own moments.
    d20_args <- list(model_solve_fn = model_solve_fn, theta = params_id,
                     param_names = param_names, moment_names = moment_names,
                     jacobian = J_shared, meta = meta)
    d20_can_sample <- !is.null(data) && !is.null(obs_names) && !is.null(model) &&
      !is.null(compiled) && all(obs_names %in% colnames(data))
    if (d20_weighting %in% c("auto", "sampling") && d20_can_sample) {
      d20_args$model_solve_fn <- .d20_data_moment_fn(model, compiled, params,
                                                     obs_names, max_lag = 4L)
      d20_args$moment_names <- NULL
      d20_args$jacobian     <- NULL
      d20_args$weighting    <- "sampling"
      # Var(m_hat): Bartlett HAC (Andrews plug-in bandwidth) / T.
      d20_args$moment_cov   <- .d20_moment_sampling_cov(
        as.matrix(data[, obs_names, drop = FALSE]), max_lag = 4L)
    } else if (identical(d20_weighting, "sampling")) {
      d20_args$weighting <- "sampling"   # D20 notes the fallback to "none"
    }
    .msg("D20: Identification strength (Fisher information)...")
    results$d20 <- .safe_diag("D20", function() do.call(d20_fisher_identification_strength, d20_args))

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
    ## D23 differentiates the OBSERVABLES' spectral density, so it needs a
    ## theta -> state-space closure (model_solve_fn is theta -> moments and
    ## cannot serve). Build it from model + compiled + obs_names; Sigma_e is
    ## then re-derived at every theta (estimated stderrs get a derivative).
    ## Without those inputs D23 reports INFO (rank not checked).
    results$d23 <- .safe_diag("D23", function() {
      if (!is.null(model) && !is.null(compiled) && length(obs_names) > 0L) {
        d23_fn <- .d23_state_space_solve_fn(model, compiled, obs_names)
        d23_dr <- d23_fn(params_id)
        if (is.null(d23_dr))
          .dynhr_abort("the model does not solve at params.")
        d23_spectral_identification(
          dr = d23_dr, model_solve_fn = d23_fn, theta = params_id,
          param_names = param_names, Sigma_e = sigma_e, meta = meta)
      } else {
        d23_spectral_identification(
          dr = dr, model_solve_fn = model_solve_fn, theta = params_id,
          param_names = param_names, Sigma_e = sigma_e, meta = meta)
      }
    })

    .msg("D24: Global identification via KL divergence...")
    # D24 needs a re-solve closure: abcd_solve_fn (the D37 state space) is the
    # only one available here. It solves the competitive-equilibrium model, so
    # it is not used with a Ramsey DR (D24 then skips as unresolvable).
    results$d24 <- .safe_diag("D24", function()
      d24_global_kl_identification(
        dr = dr, model = model, params = params_id, ramsey_result = ramsey_result,
        Sigma_e = sigma_e, param_names = param_names,
        abcd_solve_fn = if (is.null(ramsey_result)) abcd_solve_fn else NULL,
        T_obs = if (!is.null(data)) NROW(data) else NULL,
        meta = meta, verbose = verbose))
  } else {
    .msg("D23: Skipped (dr or params not provided)")
    .msg("D24: Skipped (dr or params not provided)")
  }

  # D25 re-solves the model at order 2 at each perturbed theta (until 0.9.4 it
  # differentiated a FIXED dr, so every Jacobian was zero). It runs when the
  # model was solved at order >= 2 (dr2, or a DecisionRules2 dr); a Ramsey
  # solution cannot be re-solved here, so D25 then reports INFO.
  d25_has_o2 <- inherits(dr2, "DecisionRules2") || inherits(dr, "DecisionRules2")
  if (!is.null(params_id) && !is.null(ramsey_result)) {
    .msg("D25: Higher-order identification...")
    results$d25 <- .safe_diag("D25", function()
      d25_higher_order_identification(params = params_id,
                                      ramsey_result = ramsey_result, meta = meta))
  } else if (!is.null(params_id) && d25_has_o2 && !is.null(model) &&
             !is.null(compiled)) {
    .msg("D25: Higher-order identification (pruned order 2)...")
    results$d25 <- .safe_diag("D25", function()
      d25_higher_order_identification(
        params = params_id, model = model, compiled = compiled,
        obs_names = if (length(obs_names) > 0L) obs_names,
        Sigma_e = sigma_e, param_names = names(params_id),
        meta = meta, verbose = verbose))
  } else {
    .msg("D25: Skipped (needs params, model, compiled and an order >= 2 solution)")
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
        param_names = param_names, moment_names = moment_names, meta = meta))
  } else {
    .msg("D22: Skipped (model_solve_fn, params, or obs_names not provided)")
  }

  if (!is.null(model_solve_fn) && !is.null(params) &&
      !is.null(param_names_c) && !is.null(param_names_e)) {
    .msg("D26: Calibration sensitivity diagnostic...")
    # D26 needs moments that move with the CALIBRATED (structural) parameters.
    # run_model_diagnostics()'s model_solve_fn reuses the fixed `dr`, so every
    # calibrated column was exactly zero and D26 reported "robust". With the
    # model available, use the re-solving var/acv moment function (as D20),
    # weighted by the data moments' sampling precision when data are given.
    results$d26 <- .safe_diag("D26", function() {
      d26_fn <- model_solve_fn
      d26_W  <- NULL
      if (!is.null(model) && !is.null(compiled) && !is.null(obs_names)) {
        d26_fn <- .d20_data_moment_fn(model, compiled, params, obs_names, max_lag = 4L)
        if (!is.null(data) && all(obs_names %in% colnames(data))) {
          d26_Om <- .d20_moment_sampling_cov(
            as.matrix(data[, obs_names, drop = FALSE]), max_lag = 4L)
          d26_W <- crossprod(.robust_Omega_inv_sqrt(d26_Om))
          dimnames(d26_W) <- dimnames(d26_Om)
        }
      }
      d26_missing <- setdiff(c(param_names_c, param_names_e), names(params))
      if (length(d26_missing))
        .dynhr_abort("D26: not in `params`: ", paste(d26_missing, collapse = ", "), ".")
      ## Calibrate-vs-estimate ranking (Alegre Canton 2026): use each
      ## parameter's 90% prior interval as its plausible range where a prior
      ## exists (the paper's Section 5.4 choice); the rest keep D26's
      ## magnitude default. `n_obs` is deliberately NOT passed: d26_W is the
      ## inverse of the moment SAMPLING covariance (already O(1/T)), so the
      ## paper's sqrt(log n / n) weak-identification cutoff -- stated for a
      ## sample-size-free weight -- would pass trivially.
      d26_choice <- list()
      d26_rg <- .d26_prior_ranges(priors, c(param_names_c, param_names_e))
      if (length(d26_rg)) d26_choice$ranges <- d26_rg
      d26_calibration_sensitivity(
        model_solve_fn = d26_fn,
        theta_c = unlist(params[param_names_c]),
        theta_e = unlist(params[param_names_e]),
        param_names_c = param_names_c, param_names_e = param_names_e,
        weight = d26_W, meta = meta, calibration_choice = d26_choice)
    })
  } else {
    .msg("D26: Skipped (model_solve_fn, params, param_names_c, or param_names_e not provided)")
  }

  if (!is.null(model_solve_fn) && !is.null(prior_draw_fn) && !is.null(data_moments)) {
    .msg("D4: Prior predictive checks...")
    # The posterior path hands data_moments as list(sigma_y, mean_y) (the D9
    # form); its model_solve_fn returns sd_<var>/vd_<var>_<shock>. D4 needs a
    # named numeric vector, so compare the observables' sample SDs by name.
    # (Passing the list straight through skipped all 500 draws -> always NA.)
    d4_data_moments <- data_moments
    if (is.list(data_moments) && !is.null(data_moments$sigma_y)) {
      sig <- as.matrix(data_moments$sigma_y)
      d4_data_moments <- stats::setNames(sqrt(diag(sig)),
                                         paste0("sd_", colnames(sig)))
    }
    results$d4 <- .safe_diag("D4", function()
      d4_prior_predictive(model_solve_fn, prior_draw_fn, d4_data_moments,
                          moment_names, meta = meta))
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

  ## D6 needs the prior actually used in estimation; build it from `priors`
  ## (via log_prior(), truncation-normalised) when no density fn was given.
  if (is.null(prior_density_fn) && !is.null(priors) &&
      all(c("name", "distribution", "p1", "p2") %in% names(priors)))
    prior_density_fn <- .d6_prior_density_fn(priors)
  if (!is.null(draws) && !is.null(prior_density_fn)) {
    .msg("D6: Posterior vs prior (density overlay)...")
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
    # Raw data give the full sample-moment set (SD, ACF, cross-corr); the
    # posterior path's data_moments is only list(sigma_y = cov, mean_y), so it
    # is the fallback when no data matrix is available.
    results$d9 <- .safe_diag("D9", function() {
      ## Measurement error belongs in the MODEL SD: the data SD is the SD of
      ## y = y* + u, so the like-for-like model quantity carries Var(u) too.
      ## The posterior path knows the estimation's ME through `mode_result`.
      me_var_d9 <- .d9_me_var_from_mode(mode_result, obs_names)
      ## Posterior-predictive band (Faust-Gupta): INFO only, and only when we
      ## have draws AND can re-solve + simulate. One model solve + one
      ## T-length simulation per draw.
      ppc_fn_d9 <- .d9_ppc_sd_fn(model, compiled, obs_names, params,
                                 T_sim = if (!is.null(data)) nrow(as.matrix(data))
                                         else NULL)
      d9_args <- list(model_moments = model_moments, obs_names = obs_names,
                      me_var = me_var_d9, ppc_sd_fn = ppc_fn_d9,
                      ppc_draws = if (is.null(ppc_fn_d9)) NULL else draws,
                      meta = meta)
      if (!is.null(data)) d9_args$data <- data
      else d9_args$data_moments <- data_moments
      do.call(d9_moment_matching, d9_args)
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
    results$d11 <- .safe_diag("D11", function()
      d11_historical_decomposition(hd, dates = dates, var_names = obs_names,
                                   meta = meta))
  } else {
    .msg("D11: Skipped (hd not provided)")
  }

  if (!is.null(shocks)) {
    .msg("D12: Smoothed shocks...")
    ## A bare matrix carries no scale: standardise by the model's shock stds
    ## when the columns are the model's shocks (a list from the posterior path
    ## already carries the exact per-period sd).
    results$d12 <- .safe_diag("D12", function() {
      sd12 <- NULL
      cn12 <- colnames(shocks)
      if (is.matrix(shocks) && !is.null(model) && !is.null(cn12) &&
          all(cn12 %in% model$varexo_names)) {
        p12 <- apply_theta_to_params(model, params %||% numeric(0))
        sd12 <- sqrt(diag(as.matrix(.get_shock_cov(model, cn12, p12))))
        names(sd12) <- cn12
      }
      d12_smoothed_shocks(shocks, meta = meta, shock_sd = sd12, dates = dates)
    })
  } else {
    .msg("D12: Skipped (shocks not provided)")
  }

  # --- Prior-sensitivity diagnostic (requires mode result) ---
  if (!is.null(solved) && !is.null(data) && !is.null(obs_names) && !is.null(mode_result)) {
    .msg("PriorSensitivity: Comparing informative vs flat priors...")
    ## The flat run MUST use the SAME measurement-error variance the
    ## informative mode was fitted with -- otherwise the mode shift this
    ## diagnostic reports conflates the prior change with an ME change, and
    ## every parameter looks prior-driven. The 0.07*var heuristic survives only
    ## as the fallback for a mode result that carries no me_variance.
    ## The Kalman path's `me_variance` is a SCALAR (see kalman-filter.R, where
    ## it is used as `rep(me_variance, n_obs)` / `me_variance * diag(n_obs)`);
    ## the old per-observable `diag(cov(data)) * 0.07` vector was silently
    ## mis-recycled there and drove the flat loglik to -Inf. Average it.
    me_var <- if (!is.null(data)) mean(diag(cov(data))) * 0.07 else NULL
    me_var_ps <- .ps_me_variance(mode_result, me_var)
    results$prior_sensitivity <- .safe_diag("PriorSensitivity", function()
      diag_prior_sensitivity(solved = solved, data = data, obs_vars = obs_names,
                              mode_inf = mode_result, me_variance = me_var_ps,
                              n_iter = 500L, verbose = FALSE))
  } else {
    .msg("PriorSensitivity: Skipped (solved, data, obs_names, or mode_result not provided)")
  }

  if (!is.null(dr) && !is.null(data)) {
    .msg("D13: Cross-equation restrictions...")
    # Posterior path: hand D13 the measurement-error variances the estimation
    # actually used, so a model whose only misfit is unmodelled ME is not
    # penalised. mode_result$me_variance is the constant per-observable ME
    # (scalar = same for all); the time-varying me_extra is NOT included --
    # D13 compares stationary second moments.
    d13_me <- .d13_me_var_from_mode(mode_result, obs_names %||% colnames(data))
    results$d13 <- .safe_diag("D13", function()
      d13_cross_equation_restrictions(
        dr = dr, data = data, obs_names = obs_names,
        sigma_e = sigma_e, model = model, me_var = d13_me, meta = meta))
  } else {
    .msg("D13: Cross-equation restrictions (placeholder -- dr or data not provided)...")
    results$d13 <- .safe_diag("D13", function() d13_cross_equation_restrictions())
  }

  if (!is.null(models_list)) {
    .msg("D14: Bayes factor model comparison...")
    results$d14 <- .safe_diag("D14", function()
      d14_bayes_factor(models_list, log_marglik_se = log_marglik_se, meta = meta))
  } else {
    .msg("D14: Skipped (models_list not provided)")
  }

  if (!is.null(dr) && !is.null(data)) {
    .msg("D15: DSGE-VAR tightness...")
    results$d15 <- .safe_diag("D15", function()
      ## `data` reaches the orchestrator in LEVELS; d15_dsge_var() now demeans
      ## by the model steady state itself (demean = "steady_state" default),
      ## which is what DS04's zero-intercept dummy-observation prior assumes.
      d15_dsge_var(
        dr        = dr,
        data      = data,
        obs_names = obs_names,
        sigma_e   = sigma_e,
        model     = model,
        params    = params,
        meta      = meta,
        verbose   = verbose
      ))
  } else {
    .msg("D15: DSGE-VAR (placeholder -- dr or data not provided)...")
    results$d15 <- .safe_diag("D15", function() d15_dsge_var())
  }

  if (!is.null(draws) && !is.null(results_sub)) {
    .msg("D16: Subsample stability...")
    results$d16 <- .safe_diag("D16", function()
      d16_subsample_stability(draws, results_sub, param_names, meta = meta))
    # D34 reuses D16's regime draws (no extra estimation) and adds the
    # policy/private partition + super-exogeneity (operational Lucas critique).
    .msg("D34: Policy-partitioned invariance (Lucas)...")
    results$d34 <- .safe_diag("D34", function()
      d34_policy_invariance(model = model, draws = draws,
                            results_sub = results_sub, param_names = param_names,
                            meta = meta))
  } else {
    .msg("D16/D34: Skipped (draws or results_sub not provided)")
  }

  ## D17: episodes and categories come from the model's @dynhr:narratives /
  ## @dynhr:shock_categories blocks (parse_mod_file() stores them in
  ## model$metadata); rows of `hd` are matched to `dates`. Until 0.9.4 D17 was
  ## never run here, although the report templates describe it.
  meta_d17 <- if (is.list(model)) model$metadata else NULL
  narr_d17 <- if (is.list(meta_d17)) meta_d17$narratives else NULL
  if (is.list(hd) && !is.null(hd$contributions) && !is.null(dates) &&
      length(narr_d17) > 0L) {
    .msg("D17: Narrative identification...")
    results$d17 <- .safe_diag("D17", function()
      d17_narrative_identification(
        hist_decomp      = hd,
        dates            = dates,
        episodes         = narr_d17,
        shock_categories = meta_d17$shock_categories %||% list(),
        meta             = meta))
  } else {
    .msg("D17: Skipped (needs hd, dates and a @dynhr:narratives block in the model)")
  }

  .msg("D18: Welfare plausibility...")
  results$d18 <- .safe_diag("D18", function()
    d18_welfare_plausibility(ramsey_result = ramsey_result, meta = meta))

  # --- Phase H: Frontier identification ---

  # D27: OBC / Piecewise-Linear Identification
  # Only applicable when the model has OBC/MCP constraints: an explicit
  # obc_specs, or MCP equation tags parseable from the model.
  # D27 parses the tags itself (at every FD point, so a parameter-valued bound
  # moves with theta); here only check that an mcp tag is present. A malformed
  # tag then errors inside .safe_diag instead of silently skipping D27.
  has_mcp_d27 <- length(obc_specs) > 0L ||
    (!is.null(model) && any(vapply(model$equations, function(e)
      is.list(e) &&
        grepl("mcp", paste(c(e$tag_raw, e$tag), collapse = " "), fixed = TRUE),
      logical(1))))
  if (!is.null(model) && !is.null(params) && has_mcp_d27) {
    .msg("D27: OBC piecewise-linear identification...")
    results$d27 <- .safe_diag("D27", function()
      d27_obc_identification(model = model, params = params,
                              compiled = compiled, obc_specs = obc_specs,
                              param_names = param_names, obs_names = obs_names,
                              verbose = verbose, meta = meta))
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

  # D29: Data constraints (requires data). The model moments RE-SOLVE the
  # model at each theta (D20's helper, needs `compiled`); until 0.9.4 this
  # closure reused the fixed `dr` (and dr$Sigma_e), so the whole moment
  # Jacobian was zero and every parameter was flagged. Without `compiled`
  # the user's model_solve_fn is used (its moments are aligned by name).
  d29_obs_ok <- !is.null(data) && !is.null(obs_names) &&
    all(obs_names %in% colnames(data))
  d29_fn <- if (d29_obs_ok && !is.null(model) && !is.null(compiled)) {
    .d20_data_moment_fn(model, compiled, params, obs_names, max_lag = 4L)
  } else model_solve_fn
  # 0.9.4: Var(m-hat) is model-implied at theta (Christiano-Eichenbaum-
  # Trabandt) rather than a sample HAC, which over-rejected the S test.
  d29_acov <- if (d29_obs_ok && !is.null(model) && !is.null(compiled)) {
    acov_raw <- .d29_model_acov_fn(model, compiled, params, obs_names)
    function(K) acov_raw(params_id, K)
  } else NULL
  if (!is.null(params) && d29_obs_ok && !is.null(d29_fn)) {
    .msg("D29: Data constraints (Stock-Wright S test, INFO)...")
    results$d29 <- .safe_diag("D29", function()
      d29_data_driven_constraints(
        data = as.matrix(data[, obs_names, drop = FALSE]), model_solve_fn = d29_fn,
        theta = params_id, param_names = param_names, max_lag = 4L,
        model_acov_fn = d29_acov,
        verbose = verbose, meta = meta))
  } else {
    .msg("D29: Skipped (needs params, data with the obs_names columns, and compiled or model_solve_fn)")
  }

  # D30: Arbitrary-Precision Rank Checks
  if (!is.null(model_solve_fn) && !is.null(params)) {
    .msg("D30: Arbitrary-precision rank checks...")
    results$d30 <- .safe_diag("D30", function()
      d30_arbitrary_precision_rank(
        model_solve_fn = model_solve_fn, theta = params_id,
        param_names = param_names, moment_names = moment_names,
        hp_solve_fn = d30_hp_solve_fn,
        prec_bits = d30_prec_bits, backend = d30_backend, verbose = verbose,
        meta = meta))
  } else {
    .msg("D30: Skipped (model_solve_fn or params not provided)")
  }

  # --- D31 / D32: OBC scenario comparison + binding summary ---
  if (!is.null(obc_specs) && length(obc_specs) > 0L) {
    .msg("D31: OBC scenario comparison...")
    results$d31 <- .safe_diag("D31", function()
      # `dr` carries no system matrices, and may be solved at other parameter
      # values than `params` (posterior path): D31 re-solves both at `params`.
      d31_obc_scenario_comparison(
        obc_specs = obc_specs, model = model, compiled = compiled,
        params = params,
        constraint_name = {
          nm31 <- unlist(lapply(obc_specs, function(s) s$name))
          if (length(nm31) == length(obc_specs)) paste(nm31, collapse = " + ")
          else "OBC"
        },
        meta = meta))

    .msg("D32: OBC binding summary...")
    ## D32 auto-generates an IRF regime path, which needs the system matrices.
    ## No solver stores them on `dr` (dr$sys_mat was always NULL, so D32 was
    ## permanently a placeholder): rebuild them from compiled + ss + params.
    results$d32 <- .safe_diag("D32", function() {
      sys32 <- if (!is.null(dr)) dr$sys_mat else NULL
      ss32  <- if (is.list(ss)) ss$values else ss
      p32   <- params %||% (if (!is.null(model)) model$param_values)
      if (is.null(sys32) && !is.null(compiled) && !is.null(ss32) &&
          !is.null(p32))
        sys32 <- extract_system_matrices_fast(
          cache_system_structure(compiled), ss32, p32)
      Se32 <- sigma_e
      if (is.null(Se32) && !is.null(model) && !is.null(dr) && !is.null(p32))
        Se32 <- .get_shock_cov(model, dr$exo_names, p32)
      d32_obc_binding_summary(
        sys = sys32, dr_slack = dr, Sigma_e = Se32,
        obc_specs = obc_specs, obs_idx = if (!is.null(dr)) dr$obs_idx else NULL,
        model = model, meta = meta)
    })
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
  # Evaluated at `params`: with `compiled` the model is RE-SOLVED there (on the
  # posterior path `dr` is the calibration-time solution while `params` is the
  # posterior mean -- until 0.9.4 D41 mixed the two). The estimation's
  # likelihood settings come from `mode_result`. D41 itself reports INFO when
  # its thin filter cannot reproduce that likelihood (me_extra, shock_scale,
  # non-Gaussian, unit root without a stationary P0 -- the kalman_filter
  # lik_init = "auto" rule).
  if (!is.null(dr) && !is.null(model) && !is.null(params) && !is.null(data) &&
      !is.null(obs_names)) {
    .msg("D41: KF innovation whiteness...")
    mr41 <- mode_result
    lik41 <- if (!is.null(mr41$obc_specs)) "obc" else mr41$likelihood %||% "gaussian"
    results$d41 <- .safe_diag("D41", function()
      d41_innovation_whiteness(data, dr = dr, model = model, params = params,
                               obs_vars = obs_names, compiled = compiled,
                               me_variance = mr41$me_variance %||% 0,
                               me_extra = mr41$me_extra,
                               shock_scale = mr41$shock_scale,
                               lik_init = mr41$ctx$lik_init %||% "auto",
                               likelihood = lik41, meta = meta))
  } else {
    .msg("D41: Skipped (dr, model, params, data, or obs_names not provided)")
  }

  # --- D36: calibration deepness (is each *fixed* deep param data-consistent?) ---
  # Needs a full-vector log-likelihood; feeds the Passport's calibrated axis.
  # User-supplied only (never built here): it must RE-SOLVE the model at each
  # vector (e.g. sum(make_loglik_contrib(...)(p))) -- a fixed-dr closure makes
  # every profile flat, which D36 now reports as INFO, not PASS.
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

  # ---- Provenance ---------------------------------------------------------
  # A diagnostic report is a claim about a SPECIFIC model at a specific commit
  # against specific data. Before 0.9.4 the only provenance in the rendered
  # report was sessionInfo(), which runs in the Quarto subprocess and reports
  # whichever dynhr is INSTALLED -- not the one that produced these results.
  # Everything the suite can know about its own inputs is recorded here, in
  # the parent process, and `.report_meta()` renders it. Unknowns stay NA;
  # nothing here is inferred from the render environment.
  attr(results, "provenance") <- list(
    model_name  = model_name,
    model_file  = model$source_file %||% NA_character_,
    data_hash   = meta$data_hash %||% NA_character_,
    n_obs       = if (!is.null(data)) nrow(data) else NA_integer_,
    obs_vars    = obs_names,
    sample_span = .diag_sample_span(data, dates),
    frequency   = .diag_frequency(dates),
    seed        = attr(draws, "seed") %||% mode_result$seed %||% NA,
    sampler     = attr(draws, "sampler") %||% mode_result$sampler %||%
                  NA_character_,
    n_draws     = if (!is.null(draws)) nrow(draws) else NA_integer_,
    n_warmup    = attr(draws, "n_warmup") %||% mode_result$n_warmup %||% NA,
    n_chains    = if (!is.null(chains_list)) length(chains_list) else NA_integer_
  )

  # Print summary
  if (verbose) {
    .dynhr_inform("\n", paste(rep("=", 70), collapse = ""))
    .dynhr_inform("dynhr diagnostic suite -- summary")
    .dynhr_inform(paste(rep("=", 70), collapse = ""))

    for (nm in names(results)) {
      r <- results[[nm]]
      if (inherits(r, "dynhr_diagnostic")) {
        badge <- .badge_str(r)
        # One line per diagnostic: collapse embedded newlines, then truncate
        flat_summary  <- gsub("[[:space:]]+", " ", r$summary)
        short_summary <- substr(flat_summary, 1, 100)
        if (nchar(flat_summary) > 100) short_summary <- paste0(short_summary, "...")
        .dynhr_inform(sprintf("  [%-5s] %s: %s", badge, nm, short_summary))
      }
    }

    cnt <- .badge_counts(
      Filter(function(r) inherits(r, "dynhr_diagnostic"), results))
    .dynhr_inform(sprintf("\nTotal: %d PASS, %d WARN, %d FAIL, %d ERROR, %d INFO",
                    cnt[["PASS"]], cnt[["WARN"]], cnt[["FAIL"]],
                    cnt[["ERROR"]], cnt[["INFO"]]))
    .dynhr_inform(paste(rep("=", 70), collapse = ""))
  }

  results
}


#' Convenience wrapper: run D8/D9/D10 on a StochSimulResult
#' @noRd
run_model_diagnostics <- function(model, dr = NULL, ss = NULL, verbose = TRUE) {
  results <- list()

  if (verbose) .dynhr_cat("\n-- Running model_structure_summary...\n")
  results$structure <- model_structure_summary(model, dr, ss)
  if (verbose) .dynhr_cat(results$structure$summary, "\n")

  irfs <- if (!is.null(dr$irfs)) dr$irfs else if (!is.null(dr$irf)) dr$irf else NULL
  if (!is.null(irfs)) {
    if (verbose) .dynhr_cat("\n-- Running d8_irf_plausibility...\n")
    results$d8 <- d8_irf_plausibility(irfs)
    if (verbose) .dynhr_cat(results$d8$summary, "\n")
  } else {
    if (verbose) .dynhr_cat("\n-- D8 skipped (no IRFs)\n")
  }

  if (!is.null(dr$moments)) {
    if (verbose) .dynhr_cat("\n-- Running d9_moment_matching...\n")
    results$d9 <- d9_moment_matching(dr)
    if (verbose) .dynhr_cat(results$d9$summary, "\n")
  } else {
    if (verbose) .dynhr_cat("\n-- D9 skipped (no moments)\n")
  }

  if (!is.null(dr$moments$var_decomp_pct)) {
    if (verbose) .dynhr_cat("\n-- Running d10_variance_decomposition...\n")
    results$d10 <- d10_variance_decomposition(dr)
    if (verbose) .dynhr_cat(results$d10$summary, "\n")
  } else {
    if (verbose) .dynhr_cat("\n-- D10 skipped (no variance decomposition)\n")
  }

  if (verbose) {
    .dynhr_cat("\n-- Diagnostic summary -----------------------------\n")
    for (nm in names(results)) {
      r <- results[[nm]]
      status <- switch(.badge_str(r), INFO = "N/A", .badge_str(r))
      .dynhr_cat(sprintf("  %-20s  %s\n", nm, status))
    }
    .dynhr_cat("\n")
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

  ## ---- ONE point: the posterior-mean re-solve (0.9.4) ---------
  ## Until 0.9.4 this wrapper mixed two points: `dr` and `ss` were the
  ## CALIBRATION / mode-time solution off `solved`, while `params` was the
  ## posterior mean. Every diagnostic that took `dr`, `ss` and `params`
  ## together (D0's residual, D32's `sys`, D25/D31, ...) was evaluating an
  ## inconsistent model, and `model_solve_fn` below reused the fixed `dr`, so
  ## structural parameters moved NO moment: D1/D3/D4/D22/D26/D30 saw
  ## identically-zero Jacobians on this path.
  ##
  ## Re-solve once here, through the likelihood's own pipeline (steady state
  ## + first-order perturbation), and use that single point everywhere.
  ## `lik_init = "diffuse"` so a unit-root posterior mean is still described
  ## rather than rejected. When it does not solve we fall back to the
  ## calibration point and say so, rather than silently mislabelling it.
  post_sol <- NULL
  if (!is.null(model) && !is.null(compiled) && !is.null(params)) {
    post_sol <- .solve_dr_for_theta(model, compiled,
                                    cache_system_structure(compiled),
                                    params, new.env(parent = emptyenv()),
                                    lik_init = "diffuse")
    if (is.null(post_sol))
      .dynhr_warn("Diagnostics: the model does not solve at the posterior ",
                  "mean; falling back to the mode-time solution. Results are ",
                  "labelled posterior-mean but describe the mode.")
  }
  ## Full parameter vector at the posterior mean (structural params AND
  ## estimated `stderr <shock>` entries, which are shock-named and absent
  ## from model$param_values). Superset of `posterior$posterior_mean`, so
  ## run_all_diagnostics()'s `params[param_names]` indexing is unchanged,
  ## while D26's calibrated names and D32's `sys` rebuild now resolve.
  params_post <- params
  if (!is.null(post_sol)) {
    dr <- post_sol$dr
    ss <- if (inherits(dr$ys, "dynhr_steady")) dr$ys$values else dr$ys
    params_post <- post_sol$params
    for (nm in names(params))
      if (!nm %in% names(params_post)) params_post[[nm]] <- params[[nm]]
  }
  data       <- mode_res$data
  obs_names  <- mode_res$obs_vars
  param_names <- priors$name
  draws      <- posterior$pooled_draws
  ## Provenance (0.9.4 report review D3): run_all_diagnostics() reads the
  ## sampler settings off attributes of `draws`; the posterior result keeps
  ## them in $meta, so hand them over here rather than reporting "not
  ## recorded" for a run whose settings are known.
  if (!is.null(draws)) {
    pm <- posterior$meta %||% list()
    if (is.null(attr(draws, "sampler")) && !is.null(pm$methods))
      attr(draws, "sampler") <- paste(pm$methods, collapse = "+")
    if (is.null(attr(draws, "n_warmup")) && !is.null(pm$n_warmup))
      attr(draws, "n_warmup") <- pm$n_warmup
    if (is.null(attr(draws, "seed")) && !is.null(pm$seed))
      attr(draws, "seed") <- pm$seed
  }
  model_moments <- posterior$posterior_moments

  ## D8 input: IRFs AT THE POSTERIOR MEAN, re-solved here. Not
  ## posterior$posterior_irfs: run_posterior_estimation() builds those from the
  ## calibration-time `solved$dr` (ghx/ghu at the calibration, and its
  ## dr$Sigma_e snapshot overrides params), and drops estimated `stderr
  ## <shock>` values (shock-named, not in model$param_values) -- so they are
  ## the CALIBRATED IRFs labelled as posterior-mean ones. On nk_demo the e_g
  ## impact on g came out 0.300 (calibrated std) against a posterior mean of
  ## 0.287. An infeasible posterior mean leaves irf NULL (D8 skipped), never
  ## the stale collection.
  ##
  ## Since 0.9.4 the posterior-mean re-solve happens once above, so this
  ## reuses `post_sol` instead of solving a second time (same numbers as
  ## .d8_irfs_at_theta(), which is what it used to call).
  irf <- NULL
  if (!is.null(model) && !is.null(compiled)) {
    if (is.null(post_sol)) {
      .dynhr_inform(paste0("D8: the model does not solve at the posterior ",
                           "mean; IRF plausibility is skipped."),
                    verbose = verbose)
    } else {
      dr_irf <- post_sol$dr
      ## dr$Sigma_e is a solve-time snapshot of the CALIBRATED covariance and
      ## compute_irfs() would prefer it over params; params stay authoritative.
      dr_irf$Sigma_e <- NULL
      irf <- compute_irfs(dr_irf, model, n_periods = 40L,
                          params = post_sol$params)
    }
  }

  ## ---- D10 / D11 / D12 inputs (previously unreachable from a posterior) ----
  ## D10: variance decomposition, already computed at the posterior mean.
  ## Pass the whole compute_moments() list (not just $var_decomp_pct) so D10
  ## sees Sigma_e / var_cov and can flag the correlated-shock caveat.
  vd <- if (!is.null(model_moments$var_decomp_pct)) model_moments else NULL

  ## D11/D12: run the Kalman smoother at the posterior mean to recover the
  ## smoothed shocks and the historical decomposition. Gaussian-linear only; for
  ## OBC / non-Gaussian likelihoods the linear smoother does not apply, so these
  ## stay NULL and D11/D12 are skipped exactly as before (graceful degradation).
  hd <- NULL; smoothed_shocks <- NULL; shocks_d12 <- NULL
  use_obc_post <- !is.null(mode_res$obc_specs)
  if (!use_obc_post && !is.null(data) && !is.null(obs_names) &&
      !is.null(model) && !is.null(dr)) {
    ## Evaluate at the POSTERIOR MEAN. Until 0.9.4 this used `dr` and
    ## model$param_values -- the solution the user handed to mode-finding
    ## (the calibration), not the estimate -- so D11/D12 described the
    ## calibrated model whatever the data said. The likelihood settings of the
    ## estimation (ME, filter tunes, shock scaling, init) are carried too, or
    ## the smoothed shocks absorb the measurement error.
    lik_init_sm <- mode_res$ctx$lik_init %||% "auto"
    shock_scale_sm <- mode_res$shock_scale
    sp_sm <- tryCatch({
      if (is.null(params) || is.null(compiled)) {
        build_dsge_state_space(model, dr, obs_names, verbose = FALSE)
      } else {
        sol_pm <- .solve_dr_for_theta(model, compiled,
                                      cache_system_structure(compiled),
                                      params, new.env(parent = emptyenv()),
                                      lik_init = lik_init_sm,
                                      shock_scale = shock_scale_sm)
        if (is.null(sol_pm)) {
          .dynhr_inform("D11/D12: the model does not solve at the posterior ",
                        "mean; smoothed shocks skipped.")
          NULL
        } else {
          build_dsge_state_space(model, sol_pm$dr, obs_names, verbose = FALSE,
                                 params = sol_pm$params)
        }
      }
    }, error = function(e) .dynhr_reraise_bug(e, NULL))
    if (!is.null(sp_sm)) {
      ## data is T x n_obs (same orientation cov(data)/colMeans(data) use); the
      ## smoother takes Y as T x n_obs (nrow = T), in LEVELS -- the state space
      ## carries the observation intercept and subtracts it.
      ##
      ## Until 0.9.3 it did not: the smoother reached through a pre-built state
      ## space silently required DEVIATIONS, and this call handed it the raw
      ## level series, so D11/D12 ran on an un-demeaned likelihood. On nk_demo
      ## (observable steady states 0.5, 2, 4) that put the smoothed shocks out
      ## by 6x -- max |eps| 4.83 against 0.80 -- and the loglik at -33990.3
      ## against -757.6. Any historical decomposition produced by this path
      ## before 0.9.3 on a model with non-zero observable steady states is wrong.
      Y_sm <- as.matrix(data[, obs_names, drop = FALSE])
      sm_args <- list(me_variance = mode_res$me_variance %||% 0,
                      me_extra = mode_res$me_extra,
                      shock_scale = shock_scale_sm, lik_init = lik_init_sm)
      sm_out <- tryCatch({
        o <- do.call(.kalman_smoother_ss, c(list(Y_sm, sp_sm), sm_args))
        ## Per-period sd of the smoothed shock, same settings: D12 divides by
        ## it (dividing by sigma_j or by 1 mis-states the model's claim).
        o$shock_sd <- do.call(.d12_smoothed_shock_sd,
                              c(list(Y_sm, sp_sm), sm_args))
        o
      }, error = function(e) .dynhr_reraise_bug(e, NULL))
      if (!is.null(sm_out)) {
        smoothed_shocks <- sm_out$smoothed_shocks
        shocks_d12 <- list(smoothed_shocks = sm_out$smoothed_shocks,
                           shock_sd = sm_out$shock_sd)
        ## Pass the WHOLE smoother result: that supplies s_{0|T} (the
        ## "initial" component) and the smoothed states (the adding-up check
        ## D11 reports). The bare shock matrix defaulted s0 to zero, so the
        ## components did not add up to the data -- off by 1.07 on a
        ## two-shock AR(1) fixture whose data start away from steady state.
        hd <- tryCatch(historical_decomposition(sm_out, sp_sm),
                       error = function(e) .dynhr_reraise_bug(e, NULL))
      }
    }
  }

  ## Posterior IRF credible bands (re-solves the model at many draws -- minutes;
  ## opt-in via bayesian_irf = TRUE). Stored as results$bayesian_irf below.
  bayesian_irf_result <- NULL
  if (isTRUE(bayesian_irf) && !is.null(draws) && !is.null(model) &&
      !is.null(compiled) && nrow(as.matrix(draws)) >= 10L) {
    bayesian_irf_result <- tryCatch(
      diag_bayesian_irf(model, compiled, draws),
      error = function(e) .dynhr_reraise_bug(e, NULL))
  }

  # Build a model_solve_fn from the compiled model.
  #
  # 0.9.4: this used to close over the FIXED `dr` and only pass
  # `theta` down to compute_moments(), which rescales shocks but cannot move
  # ghx/ghu. Every structural column of the Jacobian was therefore EXACTLY
  # zero, and D1 reported "identified", D3 "negligible", D4 degenerate, D22
  # uninformative, D26 "robust" -- on every posterior-path run. It now
  # RE-SOLVES at each theta, the same way .d20_data_moment_fn() does (same
  # .solve_dr_for_theta() pipeline, same warm-start cache), at a cost of one
  # model solve per evaluation.
  model_solve_fn <- NULL
  if (!is.null(compiled) && !is.null(model) && !is.null(post_sol)) {
    ms_cache <- cache_system_structure(compiled)
    ms_state <- new.env(parent = emptyenv())
    ms_state$ss_warm <- NULL
    ms_moments <- function(sol) {
      moments  <- compute_moments(sol$dr, model, params = sol$params)
      sd_named <- moments$std_dev
      names(sd_named) <- paste0("sd_", names(sd_named))
      vd       <- moments$var_decomp_pct
      vd_named <- as.numeric(vd)
      names(vd_named) <- paste0("vd_",
                                 rep(rownames(vd), ncol(vd)), "_",
                                 rep(colnames(vd), each = nrow(vd)))
      c(sd_named, vd_named)
    }
    ms_names <- names(ms_moments(post_sol))
    ms_na    <- stats::setNames(rep(NA_real_, length(ms_names)), ms_names)
    model_solve_fn <- function(theta) {
      pp  <- .apply_theta_to_params(model, theta, params_post)
      sol <- .solve_dr_for_theta(model, compiled, ms_cache, pp, ms_state)
      ## Convention (0.9.4): a failing user solve RETURNS non-finite values,
      ## it does not stop -- the caller reports a non-finite Jacobian.
      if (is.null(sol)) return(ms_na)
      out <- ms_moments(sol)
      if (!identical(names(out), ms_names)) return(ms_na)
      out
    }
  }

  ## D3 bounds from the prior support (0.9.4): the posterior
  ## path never built `param_bounds`, so D3 was ALWAYS skipped there. Morris
  ## screening needs a finite box; take the prior's own [lower, upper] where
  ## finite and fall back to a +/- 3 sd interval around the prior mean. If any
  ## parameter yields no finite ordered pair, stay NULL (D3 skipped, as before).
  param_bounds <- NULL
  if (!is.null(priors) && !is.null(priors$name) && length(priors$name) > 0L) {
    lo <- as.numeric(priors$lower %||% rep(NA_real_, length(priors$name)))
    hi <- as.numeric(priors$upper %||% rep(NA_real_, length(priors$name)))
    mu <- as.numeric(priors$mean  %||% rep(NA_real_, length(priors$name)))
    sdv <- as.numeric(priors$std  %||% rep(NA_real_, length(priors$name)))
    wide <- is.finite(mu) & is.finite(sdv) & sdv > 0
    lo[!is.finite(lo) & wide] <- (mu - 3 * sdv)[!is.finite(lo) & wide]
    hi[!is.finite(hi) & wide] <- (mu + 3 * sdv)[!is.finite(hi) & wide]
    if (all(is.finite(lo)) && all(is.finite(hi)) && all(hi > lo)) {
      param_bounds <- cbind(lower = lo, upper = hi)
      rownames(param_bounds) <- priors$name
    }
  }

  # Prior density function for D6
  ## Same density the samplers use (log_prior), normalised over [lower, upper];
  ## NA for a name not in the spec (was: a silent N(0,1) stand-in).
  prior_density_fn <- NULL
  if (!is.null(priors)) prior_density_fn <- .d6_prior_density_fn(priors)

  # Prior draw function for D4
  prior_draw_fn <- NULL
  if (!is.null(priors)) {
    # .smc_make_prior_sampler() already returns function() -> named theta.
    # (It used to be indexed as a list of per-parameter samplers here, which
    # yielded an EMPTY draw every time, so D4 skipped all draws.)
    prior_draw_fn <- .smc_make_prior_sampler(priors)
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

  # Run diagnostics with all derived args.
  # `param_bounds` is added through the arg list (not as a named formal) so a
  # caller-supplied `param_bounds` in `...` still wins.
  diag_args <- list(
    model            = model,
    compiled         = compiled,
    dr               = dr,
    ss               = ss,
    params           = params_post,
    priors           = priors,
    data             = data,
    draws            = draws,
    chains_list      = chains_list,
    irf              = irf,
    vd               = vd,
    hd               = hd,
    shocks           = shocks_d12,
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
  if (is.null(diag_args$param_bounds) && !is.null(param_bounds))
    diag_args$param_bounds <- param_bounds
  results <- do.call(run_all_diagnostics, diag_args)

  ## Attach posterior IRF credible bands (opt-in; NULL otherwise).
  if (!is.null(bayesian_irf_result))
    results$bayesian_irf <- bayesian_irf_result

  # ---- Report generation ----
  if (report != "none") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (report == "md") {
      report_path <- file.path(output_dir, paste0(report_file, ".md"))
      .write_diagnostic_md(results, report_path, model_name = model_name)
      if (verbose) .dynhr_inform("[dynhr] Report written -> ", report_path)
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
          .dynhr_inform("[dynhr] Quarto report -> ",
                  file.path(output_dir,
                            paste0(report_file,
                                   if (report == "pdf") ".pdf" else ".html")))
      } else {
        # Fallback: write markdown and warn clearly
        md_path <- file.path(output_dir, paste0(report_file, ".md"))
        .write_diagnostic_md(results, md_path, model_name = model_name)
        if (!quarto_ok) {
          .dynhr_warn(sprintf(
            "[dynhr] quarto package not installed; wrote Markdown to %s instead.",
            md_path))
        } else {
          .dynhr_warn(sprintf(
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
    if (verbose) .dynhr_inform("[dynhr] Executive summary -> ", summary_path)
  }

  invisible(results)
}


## ---------------------------------------------------------------------------
## D9 wiring helpers (0.9.4)
## ---------------------------------------------------------------------------

## Measurement-error VARIANCES for D9, off the estimation's mode result.
## `me_variance` is the scalar/vector the likelihood used; `me_extra` is the
## n_obs x T per-period EXTRA variance, whose per-observable MEAN is added
## (a single SD comparison cannot represent a time-varying variance; the mean
## is the variance of the ME process averaged over the sample). NULL when the
## estimation carried no measurement error, so D9 is unchanged.
.d9_me_var_from_mode <- function(mode_result, obs_names) {
  if (is.null(mode_result) || is.null(obs_names)) return(NULL)
  mv <- mode_result$me_variance
  out <- NULL
  if (!is.null(mv) && any(is.finite(mv) & mv > 0)) out <- mv
  ex <- mode_result$me_extra
  if (!is.null(ex) && is.matrix(ex) && any(is.finite(ex) & ex > 0)) {
    rm_ <- rowMeans(ex, na.rm = TRUE)
    if (is.null(names(rm_)) && length(rm_) == length(obs_names))
      names(rm_) <- obs_names
    base <- stats::setNames(rep(0, length(obs_names)), obs_names)
    if (!is.null(out)) {
      if (length(out) == 1L && is.null(names(out))) base[] <- as.numeric(out)
      else if (!is.null(names(out))) {
        hit <- intersect(obs_names, names(out)); base[hit] <- out[hit]
      } else if (length(out) == length(obs_names)) base[] <- as.numeric(out)
    }
    hit <- intersect(obs_names, names(rm_))
    base[hit] <- base[hit] + rm_[hit]
    out <- base
  }
  out
}

## theta -> sample SDs of a T_sim-length path simulated at theta, for D9's
## posterior-predictive band. NULL (band skipped) when the model cannot be
## re-solved or T_sim is unknown.
##
## Follows the package callback convention: an infeasible draw returns a
## zero-length value that .d9_ppc_band() drops -- it never stops.
.d9_ppc_sd_fn <- function(model, compiled, obs_names, params, T_sim) {
  if (is.null(model) || is.null(compiled) || is.null(obs_names) ||
      is.null(T_sim) || !is.finite(T_sim) || T_sim < 20L) return(NULL)
  cache <- cache_system_structure(compiled)
  state <- new.env(parent = emptyenv())
  state$ss_warm <- NULL
  base_params <- params
  function(theta) {
    pp  <- .apply_theta_to_params(model, theta, base_params)
    sol <- .solve_dr_for_theta(model, compiled, cache, pp, state)
    if (is.null(sol)) return(numeric(0))
    sim <- simulate_model(sol$dr, n_periods = as.integer(T_sim), model = model)
    df  <- if (is.list(sim) && !is.null(sim$data)) sim$data else sim
    df  <- as.data.frame(df)
    hit <- intersect(obs_names, names(df))
    if (length(hit) == 0L) return(numeric(0))
    vapply(df[hit], function(cl) stats::sd(as.numeric(cl)), numeric(1))
  }
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
    ## .badge_str() is the single source of truth, so the markdown report gets
    ## the 0.9.4 WARN level too; "N/A" is this renderer's label for INFO.
    status <- switch(.badge_str(r), INFO = "N/A", .badge_str(r))
    raw_sum <- if (!is.null(r$summary)) gsub("\n", " ", r$summary, fixed = TRUE) else ""
    safe_sum <- gsub("|", "\\|", raw_sum, fixed = TRUE)
    safe_nm  <- gsub("|", "\\|", nm, fixed = TRUE)
    lines <- c(lines, sprintf("| %s | %s | %s |", safe_nm, status, safe_sum))
  }
  lines <- c(lines, "")

  writeLines(lines, file_path)
}



#' Plausible-range widths for D26's calibrate-vs-estimate ranking: the 90%
#' interval of each parameter's prior (5% to 95% quantile of the law
#' `.lp_dist1()` scores, incl. Dynare p3/p4 and truncation to [lower, upper]).
#' Parameters without a (proper) prior are left out, so D26 falls back to its
#' magnitude default for them.
#' @noRd
.d26_prior_ranges <- function(priors, names_wanted) {
  if (!is.data.frame(priors) || !all(c("name", "distribution", "p1", "p2") %in%
                                      names(priors)))
    return(numeric(0))
  out <- numeric(0)
  for (i in which(priors$name %in% names_wanted)) {
    d  <- .normalize_dist(priors$distribution[i])
    if (!d %in% c("normal", "beta", "gamma", "inv_gamma", "inv_gamma1",
                  "inv_gamma2", "uniform")) next
    p3 <- if ("p3" %in% names(priors)) priors$p3[i] else NA_real_
    p4 <- if ("p4" %in% names(priors)) priors$p4[i] else NA_real_
    if (!.diag_prior_law_proper(d, priors$p1[i], priors$p2[i], p3, p4)) next
    law <- .prior_law(d, priors$p1[i], priors$p2[i], p3, p4,
                      label = priors$name[i])
    lo <- if ("lower" %in% names(priors)) priors$lower[i] else -Inf
    hi <- if ("upper" %in% names(priors)) priors$upper[i] else Inf
    lo <- if (is.finite(lo)) max(lo, law$lo) else law$lo
    hi <- if (is.finite(hi)) min(hi, law$hi) else law$hi
    ## Quantiles of the law truncated to [lo, hi] by inverse CDF.
    plo <- if (is.finite(lo) && lo > law$lo) law$p(lo) else 0
    phi <- if (is.finite(hi) && hi < law$hi) law$p(hi) else 1
    w <- law$q(plo + 0.95 * (phi - plo)) - law$q(plo + 0.05 * (phi - plo))
    if (is.finite(w) && w > 0) out[[priors$name[i]]] <- w
  }
  out
}
