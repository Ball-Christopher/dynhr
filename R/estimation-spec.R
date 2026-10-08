## R/estimation-spec.R
## ---------------------------------------------------------------------------
## The estimation spec (E5 part C1): one versioned object that fully determines
## an estimation run -- model, data, likelihood, mode finding, sampler(s),
## compute settings, post-estimation outputs and package options.
##
## Layout
##   * `.spec_schema`: the ONE table of sub-spec fields (type, default, the
##     global option a default comes from, one-line doc). The constructors,
##     update(), the flat-argument converter, the YAML/JSON document layer
##     (R/spec-serialize.R) and the Rd field lists (@eval) all read it.
##   * Sub-specs: likelihood_spec(), mode_spec(), sampler_spec(method),
##     compute_spec(), outputs_spec() -- validated field lists with print and
##     update() methods.
##   * dynhr_estimation_spec(): assembles them with the model, data and the
##     option snapshot, hashes every component and runs validate_spec(), the
##     single home of cross-field checks.
##   * as_estimation_spec(): from the flat arguments of run_full_estimation(),
##     run_posterior_estimation() and run_mode_finding(), from a run record
##     (schema 1), from an estimation result, or from a .mod file's
##     estimation(...) command.
##
## Result-changing registered options (`changes_results = TRUE` in
## `.dynhr_option_registry`) are TYPED fields of the sub-spec that uses them;
## their default is read from the option store when the sub-spec is built, so
## dynhr_set_options() seeds new specs but never changes an existing one. All
## other options are kept in `spec$options`, a snapshot.
## ---------------------------------------------------------------------------

.spec_version <- 1L

.spec_likelihood_types <- c("gaussian", "cumulant", "whittle", "tpf", "pskf",
                            "student_t", "pkf", "ppf", "copf", "pruned",
                            "sv_rbpf", "global_pf")
## Unbiased-but-noisy particle likelihoods: PMMH targets (pmmh()'s own list).
.spec_pmmh_likelihoods <- c("tpf", "ppf", "copf", "sv_rbpf")
## Noisy likelihoods: the flat run_mode_finding() does not accept them (its
## `likelihood` choices exclude them); an estimation spec does -- its mode
## stage runs the deterministic optimiser on one noisy evaluation per point --
## and validate_spec() warns about it.
.spec_noisy_mode_likelihoods <- c("tpf", "ppf", "copf", "global_pf")
## Samplers the spec runner does not run, with the message it refuses with.
.spec_unsupported_samplers <- "smc2"
.spec_unsupported_sampler_msg <- function(method)
  paste0("sampler \"", method, "\" is not run by the spec runner yet; call ",
         method, "() directly.")
## Likelihood types that exist only as OBC filters.
.spec_obc_likelihoods <- c("pkf", "ppf", "copf")

.spec_samplers <- c("rwmh", "pmmh", "nuts", "hmc", "mala", "chees", "smc",
                    "dsmh", "dime", "smc2")

## One schema entry. `option`: the registered option the default is read
## from; `option_mode = "always"` reads its effective value (registry default
## when unset -- then `default` equals the registry default), `"set"` reads it
## only when explicitly set and otherwise uses `default` (for options whose
## registry default belongs to a different code path, e.g. mode_method).
## `text = FALSE`: the value cannot go into YAML/JSON (RDS only).
.spec_f <- function(type, default = NULL, doc, choices = NULL, option = NULL,
                    option_mode = "always", nullable = FALSE, check = NULL,
                    text = TRUE) {
  list(type = type, default = default, doc = doc, choices = choices,
       option = option, option_mode = option_mode, nullable = nullable,
       check = check, text = text)
}

.spec_chk_nonneg <- function(x)
  if (any(!is.finite(x)) || any(x < 0)) "must be finite and >= 0"
.spec_chk_pos <- function(x)
  if (any(!is.finite(x)) || any(x <= 0)) "must be finite and > 0"
## The NUTS runner raises any serial time limit below this to it
## (max(30L, timeout)), so a smaller value would be silently replaced.
.spec_nuts_min_timeout <- 30L
.spec_chk_timeout <- function(x)
  if (any(!is.finite(x)) || any(x < .spec_nuts_min_timeout))
    paste0("must be >= ", .spec_nuts_min_timeout, " (seconds): the NUTS ",
           "runner never uses a shorter limit and would silently raise it")
.spec_chk_finite <- function(x)
  if (any(!is.finite(x))) "must be finite"
.spec_chk_band <- function(x)
  if (length(x) != 2L || any(!is.finite(x)) || x[1L] >= x[2L])
    "must be two finite numbers c(lo, hi) with lo < hi"
.spec_chk_named_chr <- function(x)
  if (!is.character(x) || is.null(names(x)) || any(!nzchar(names(x))))
    "must be a named character vector"
## (The match against the prior's parameter names happens in the mode stage,
## .rmf_check_theta_init(), where the prior spec is known.)
.spec_chk_theta <- function(x) {
  if (!is.numeric(x) || !is.null(dim(x)))
    return("must be a named numeric vector")
  nm <- names(x)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)))
    return(paste0("must be NAMED: its names are matched against the prior ",
                  "spec's parameter names"))
  if (anyDuplicated(nm))
    return(paste0("has duplicated name(s): ",
                  paste(unique(nm[duplicated(nm)]), collapse = ", ")))
  if (any(!is.finite(x)))
    return(paste0("has non-finite value(s) for: ",
                  paste(nm[!is.finite(x)], collapse = ", "),
                  " (NA/NaN/Inf starts are not allowed)"))
  NULL
}
.spec_chk_matrix <- function(x)
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x))
    "must be a square numeric matrix"
.spec_chk_ramsey_order <- function(x)
  if (!x %in% 1:2) "must be 1 or 2"
.spec_chk_pruned_order <- function(x)
  if (!x %in% 2:3) "must be 2 or 3"

.spec_schema <- list(
  likelihood = list(
    type = .spec_f("chr1", "gaussian", choices = .spec_likelihood_types,
      doc = "Likelihood: gaussian (Kalman filter), cumulant, whittle, tpf, pskf, student_t, pkf, ppf, copf, pruned, sv_rbpf or global_pf"),
    me_variance = .spec_f("num", 0, option = "me_variance",
      check = .spec_chk_nonneg,
      doc = "Measurement-error variance: one value, or one per observable"),
    lik_init = .spec_f("chr1", "auto",
      choices = c("auto", "stationary", "diffuse", "kappa", "fixed_unknown"),
      doc = "Kalman filter P0 initialisation (fixed_unknown: IRIS's GLS-estimated unit-root level; gaussian only)"),
    filter_method = .spec_f("chr1", "auto",
      choices = c("auto", "dare", "chandrasekhar", "standard", "reference",
                  "univariate", "univariate_ss"),
      doc = "Kalman filter recursion of the posterior's value path, as the method argument of kalman_filter() (gaussian only; see make_log_posterior())"),
    singular_obs = .spec_f("chr1", "reject", choices = c("reject", "skip"),
      doc = "Informative singular-observation policy of the Gaussian value path: reject (log-likelihood -Inf when the filter would skip a component whose forecast variance is ~0 but whose innovation is not; default) or skip (Dynare convention, finite likelihood + warning)"),
    freq_band = .spec_f("num", c(0, pi), check = .spec_chk_band,
      doc = "Whittle frequency band c(lo, hi) in radians (whittle only)"),
    pskf_cdf = .spec_f("chr1", "accurate", check = function(x) .pskf_cdf_check(x),
      option = "pskf_cdf",
      doc = "PSKF CDF evaluation: accurate (default and only value; pskf only)"),
    power_posterior = .spec_f("num1", 1, option = "power_posterior",
      check = .spec_chk_pos,
      doc = "Power-posterior tempering exponent in (0, 1] applied to the log-likelihood"),
    gradient_policy = .spec_f("chr1", "auto",
      choices = c("auto", "analytic", "numerical"),
      doc = "Gradient path of make_posterior_grad(): auto, analytic or numerical"),
    student_df = .spec_f("num1", NULL, nullable = TRUE,
      doc = "Student-t degrees of freedom (> 2); required for student_t"),
    pruned_order = .spec_f("int1", 2L, check = .spec_chk_pruned_order,
      doc = "Pruned state-space order (2 or 3) for pruned and tpf"),
    first_obs = .spec_f("int1", 1L, check = .spec_chk_pos,
      doc = "First data row used by the likelihood (Dynare first_obs)"),
    nobs = .spec_f("int1", NULL, nullable = TRUE, check = .spec_chk_pos,
      doc = "Number of data rows used from first_obs on; NULL = all (Dynare nobs)"),
    system_priors = .spec_f("any", NULL, nullable = TRUE,
      doc = "System-prior specification (system_prior_spec) or NULL"),
    tpf_options = .spec_f("list", list(),
      doc = "Options for the tempered particle filter (make_log_posterior_tpf)"),
    obc = .spec_f("lgl1", NULL, nullable = TRUE,
      doc = "OBC filter: NULL auto-detects mcp tags, TRUE forces it, FALSE disables it"),
    obc_max_inner = .spec_f("int1", 10L, check = .spec_chk_pos,
      doc = "Maximum inner PKF iterations per period (OBC only)"),
    obc_filter = .spec_f("chr1", "pkf", choices = c("pkf", "ppf", "copf"),
      doc = "OBC filter type: pkf, or the particle filters ppf / copf (OBC only)"),
    obc_options = .spec_f("list", list(),
      doc = "Options for the OBC particle filters"),
    filter_tunes = .spec_f("any", NULL, nullable = TRUE,
      doc = "filter_tunes override: NULL uses the .mod block, FALSE disables it, or a filter_tunes_spec"),
    heteroskedastic_shocks = .spec_f("any", NULL, nullable = TRUE,
      doc = "heteroskedastic_shocks override (NULL uses the .mod block)"),
    stochastic_volatility = .spec_f("any", NULL, nullable = TRUE,
      doc = "stochastic_volatility override (NULL uses the .mod block)"),
    plan = .spec_f("any", NULL, nullable = TRUE,
      doc = "A dynhr_plan bundling filter tunes and shock scales"),
    sample_start = .spec_f("chr1", NULL, nullable = TRUE,
      doc = "Sample start date for date-literal plan periods"),
    ms_spec = .spec_f("any", NULL, nullable = TRUE,
      doc = "Markov-switching shock-variance spec (ms_dsge_spec) or NULL"),
    ms_struct_spec = .spec_f("any", NULL, nullable = TRUE,
      doc = "Structural Markov-switching spec (ms_struct_spec) or NULL"),
    ms_collapse = .spec_f("chr1", "gpb2", choices = c("gpb2", "gpb3", "imm"),
      doc = "Markov-switching filter collapse depth"),
    infeasible_penalty = .spec_f("num1", NULL, nullable = TRUE,
      doc = "Log-posterior value returned at infeasible draws; NULL = -Inf"),
    data_col_map = .spec_f("any", NULL, nullable = TRUE,
      check = .spec_chk_named_chr,
      doc = "Named character vector mapping model observables to CSV column names"),
    dates = .spec_f("any", NULL, nullable = TRUE,
      doc = "Date vector for the data rows (diagnostic labels)"),
    extra = .spec_f("list", list(),
      doc = "Further arguments for the log-posterior constructor (make_log_posterior's ...)")
  ),
  mode = list(
    method = .spec_f("chr1", "newrat", option = "mode_method",
      option_mode = "set",
      doc = "Optimiser sequence (see find_mode), e.g. newrat, cmaes_newrat, combined"),
    n_iter = .spec_f("int1", 10000L, option = "nm_maxit", option_mode = "set",
      check = .spec_chk_pos, doc = "Optimiser iteration budget"),
    n_starts = .spec_f("int1", NULL, nullable = TRUE, check = .spec_chk_pos,
      doc = "Parallel multi-start count; NULL = one per daemon"),
    transform_params = .spec_f("lgl1", TRUE, option = "transform_params",
      doc = "Search in unconstrained eta-space"),
    proposal_cov = .spec_f("chr1", "diagonal", choices = c("diagonal", "full"),
      option = "proposal_cov_method",
      doc = "Proposal-covariance strategy: diagonal or full"),
    exact_hessian = .spec_f("lgl1", FALSE, option = "use_exact_hessian",
      doc = "Build the proposal from the analytic posterior Hessian"),
    analytic_h0 = .spec_f("lgl1", TRUE, option = "use_analytic_hess",
      doc = "newrat / cmaes_newrat on the standard Gaussian Kalman likelihood: seed the initial curvature H0 from the analytic posterior Hessian (FALSE: csminwel's default H0); options$use_analytic_hess overrides it"),
    perturb_scale = .spec_f("num1", 0.5, option = "perturb_scale",
      check = .spec_chk_pos,
      doc = "Multi-start perturbation scale (fraction of the prior sd)"),
    theta_init = .spec_f("any", NULL, nullable = TRUE, check = .spec_chk_theta,
      doc = "Named starting vector; NULL = INITVAL / prior means"),
    run = .spec_f("chr1", "auto", choices = c("auto", "always", "never"),
      doc = "Run the mode stage? auto: unless every sampler is prior-initialised (smc, dsmh, dime, smc2) and no output needs the mode; always: also then (e.g. to report the mode); never: not at all (an error when a sampler or output needs the mode)"),
    skip_check = .spec_f("lgl1", FALSE,
      doc = "Skip the mode-quality check before sampling"),
    options = .spec_f("list", list(),
      doc = "Extra optimiser options (run_mode_finding's mode_options)"),
    extra = .spec_f("list", list(),
      doc = "run_mode_finding's ...: forwarded to the log-posterior constructor and the optimiser"),
    result = .spec_f("any", NULL, nullable = TRUE, text = FALSE,
      doc = "A precomputed dynhr_mode_result; when set the mode stage is not run (RDS only)")
  ),
  sampler = list(
    method = .spec_f("chr1", "rwmh", choices = .spec_samplers,
      doc = "Sampler: rwmh, pmmh, nuts, hmc, mala, chees, smc, dsmh, dime or smc2"),
    n_draws = .spec_f("int1", 10000L, check = .spec_chk_pos,
      doc = "Post-warmup draws to keep (>= 1; for dime: post-warmup iterations per walker, so the draws kept are n_draws x n_walkers)"),
    n_warmup = .spec_f("int1", 5000L, check = .spec_chk_nonneg,
      doc = "Warmup / burn-in draws, discarded"),
    n_chains = .spec_f("int1", 4L, check = .spec_chk_pos,
      doc = "Number of chains"),
    n_particles = .spec_f("int1", 2000L, check = .spec_chk_pos,
      doc = "Particle count (SMC / DSMH cloud, SMC2 outer particles, PMMH preflight)"),
    n_walkers = .spec_f("int1", NULL, nullable = TRUE, check = .spec_chk_pos,
      doc = "DIME ensemble walkers; NULL = max(5 * n_par, 20)"),
    transform_params = .spec_f("lgl1", TRUE, option = "transform_params",
      doc = "Sample in unconstrained eta-space"),
    metric = .spec_f("chr1", "diagonal",
      doc = "Mass matrix / metric (choices depend on the method)"),
    adapt = .spec_f("chr1", "independent", choices = c("independent", "pooled"),
      option = "nuts_adapt",
      doc = "Parallel multi-chain NUTS warmup: independent or pooled"),
    analytic_grad = .spec_f("lgl1", TRUE,
      doc = "Use the exact (analytic) posterior gradient wherever the likelihood has one; where it has none (OBC, tpf, pskf, ...) the numerical gradient is used without a warning. FALSE: always the numerical (finite-difference) gradient"),
    grad_method = .spec_f("chr1", "auto", choices = c("auto", "hybrid", "implicit", "adjoint", "adjoint_solution"),
      option = "grad_method",
      doc = "Analytic gradient method: auto (default; resolved at build time, see make_posterior_grad), hybrid, implicit, adjoint or adjoint_solution"),
    adapt_cov = .spec_f("lgl1", FALSE, option = "rwmh_adapt_cov",
      doc = "RWMH Haario adaptive proposal covariance"),
    n_blocks = .spec_f("int1", 1L, option = "rwmh_n_blocks",
      check = .spec_chk_pos,
      doc = "RWMH randomized parameter blocks"),
    monge_alpha = .spec_f("num1", 1, option = "monge_alpha",
      check = .spec_chk_nonneg,
      doc = "Monge-metric softness alpha (metric = monge)"),
    timeout = .spec_f("int1", 300L, check = .spec_chk_timeout,
      doc = "Serial NUTS wall-clock limit in seconds (>= 30: the runner never uses a shorter limit)"),
    Sigma_prop = .spec_f("any", NULL, nullable = TRUE, check = .spec_chk_matrix,
      doc = "Proposal covariance; NULL = the mode stage's"),
    extra = .spec_f("list", list(),
      doc = "Further arguments for the sampler function")
  ),
  compute = list(
    parallel = .spec_f("lgl1", FALSE,
      doc = "Run chains / starts / Hessians on a mirai daemon pool"),
    backend = .spec_f("chr1", "mirai", choices = "mirai",
      doc = "Parallel backend"),
    n_cores = .spec_f("int1", NULL, nullable = TRUE, check = .spec_chk_pos,
      doc = "Daemon count; NULL = auto"),
    seed = .spec_f("int1", 42L, nullable = TRUE, option = "seed_base",
      doc = "RNG seed; NULL leaves the ambient RNG stream untouched"),
    verbose = .spec_f("lgl1", TRUE, doc = "Print progress"),
    checkpoint_dir = .spec_f("chr1", NULL, nullable = TRUE,
      doc = "Directory for streamed, restartable sampling; NULL = none"),
    resume = .spec_f("lgl1", FALSE,
      doc = "Continue the chains saved in checkpoint_dir"),
    on_mismatch = .spec_f("chr1", "refuse", choices = c("refuse", "warn"),
      doc = "Resume when the checkpoint's target (model, data, likelihood, mode or sampler) differs: refuse, or warn and mark the result")
  ),
  outputs = list(
    form = .spec_f("chr1", "auto", choices = c("auto", "mode", "posterior", "full"),
      doc = "Result class of run_estimation(): mode (dynhr_mode_result), posterior (dynhr_posterior_result), full (dynhr_estimation_result) or auto"),
    dir = .spec_f("chr1", ".", doc = "Directory for saved outputs"),
    prefix = .spec_f("chr1", "dynhr_est", doc = "Prefix for saved file names"),
    save = .spec_f("lgl1", FALSE,
      doc = "Save the mode and chains as RDS files in dir"),
    diagnostics = .spec_f("lgl1", FALSE, doc = "Run the diagnostic battery"),
    stoch_simul = .spec_f("lgl1", FALSE,
      doc = "IRFs and moments at the posterior mean"),
    smoother = .spec_f("lgl1", FALSE,
      doc = "OBC smoother and historical decomposition at the mode (OBC only)"),
    obc_ppf_reweight = .spec_f("lgl1", FALSE,
      doc = "PPF importance re-weighting of the PKF draws (OBC only)"),
    ramsey = .spec_f("lgl1", FALSE,
      doc = "Run the Ramsey policy workflow at the mode"),
    ramsey_order = .spec_f("int1", 1L, check = .spec_chk_ramsey_order,
      doc = "Ramsey perturbation order (1 or 2)"),
    ramsey_n_periods = .spec_f("int1", 400L, check = .spec_chk_pos,
      doc = "Ramsey welfare simulation length"),
    ramsey_burn_in = .spec_f("int1", 100L, check = .spec_chk_nonneg,
      doc = "Ramsey welfare simulation burn-in"),
    ramsey_discount = .spec_f("num1", NULL, nullable = TRUE,
      check = .spec_chk_finite,
      doc = "Planner discount factor; NULL = the beta / betta parameter")
  )
)

## Which sampler fields each method takes (in schema order, after `method`).
.spec_sampler_table <- list(
  rwmh  = c("n_draws", "n_warmup", "n_chains", "transform_params",
            "adapt_cov", "n_blocks", "Sigma_prop", "extra"),
  pmmh  = c("n_draws", "n_warmup", "n_chains", "n_particles",
            "transform_params", "adapt_cov", "n_blocks", "Sigma_prop", "extra"),
  nuts  = c("n_draws", "n_warmup", "n_chains", "transform_params", "metric",
            "adapt", "analytic_grad", "grad_method", "timeout", "Sigma_prop",
            "extra"),
  hmc   = c("n_draws", "n_warmup", "transform_params", "metric", "analytic_grad",
            "grad_method", "extra"),
  mala  = c("n_draws", "n_warmup", "transform_params", "metric",
            "monge_alpha", "analytic_grad", "grad_method", "Sigma_prop",
            "extra"),
  chees = c("n_draws", "n_warmup", "transform_params", "analytic_grad",
            "grad_method", "Sigma_prop", "extra"),
  smc   = c("n_particles", "extra"),
  dsmh  = c("n_particles", "extra"),
  dime  = c("n_draws", "n_warmup", "n_walkers", "extra"),
  smc2  = c("n_particles", "extra")
)
.spec_sampler_metrics <- list(
  nuts = c("diagonal", "hessian", "warmup_dense", "whittle_fim", "lowrank",
           "fisher_diag"),
  mala = c("diagonal", "hessian", "monge"),
  hmc  = c("diagonal", "warmup_dense"))
## Method-specific defaults that differ from the schema's.
.spec_sampler_defaults <- list(smc2 = list(n_particles = 200L))
## Samplers that stream to a checkpoint (serial path; parallel: rwmh/nuts).
.spec_checkpoint_serial   <- c("rwmh", "pmmh", "nuts", "mala", "chees", "dime")
.spec_checkpoint_parallel <- c("rwmh", "pmmh", "nuts")

.spec_components <- c("likelihood", "mode", "sampler", "compute", "outputs")
.spec_component_doc <- c(
  likelihood = "Likelihood and sample (likelihood_spec())",
  mode       = "Posterior mode finding (mode_spec())",
  sampler    = "Posterior sampler (sampler_spec()); null = no sampling",
  compute    = "Execution: parallelism, seed, checkpointing (compute_spec())",
  outputs    = "Post-estimation extras; they change no draws (outputs_spec())")

.spec_class <- function(component)
  c(paste0("dynhr_", component, "_spec"), "dynhr_subspec")

.spec_component_of <- function(x) {
  for (cp in .spec_components)
    if (inherits(x, paste0("dynhr_", cp, "_spec"))) return(cp)
  NA_character_
}

## Schema entries of one component (sampler: restricted to `method`).
.spec_component_schema <- function(component, method = NULL) {
  sch <- .spec_schema[[component]]
  if (identical(component, "sampler")) {
    keep <- c("method", .spec_sampler_table[[method]])
    sch <- sch[keep]
    if (!is.null(.spec_sampler_metrics[[method]]))
      sch$metric$choices <- .spec_sampler_metrics[[method]]
    for (nm in names(.spec_sampler_defaults[[method]]))
      sch[[nm]]$default <- .spec_sampler_defaults[[method]][[nm]]
  }
  sch
}

## Registered options that are typed spec fields, with their homes.
.spec_option_homes <- function() {
  out <- list()
  for (cp in .spec_components) {
    sch <- .spec_schema[[cp]]
    for (nm in names(sch)) {
      op <- sch[[nm]]$option
      if (!is.null(op)) out[[op]] <- c(out[[op]], paste0(cp, "$", nm))
    }
  }
  out
}

.spec_result_options <- function() {
  reg <- .dynhr_option_registry
  names(reg)[vapply(reg, function(e) isTRUE(e$changes_results), logical(1))]
}

# ---------------------------------------------------------------------------
# Option-sourced defaults
# ---------------------------------------------------------------------------

## A scoped option source: while `$snap` is non-NULL (an option snapshot such
## as a run record's `options`, with attribute "set"), `.spec_opt()` reads it
## instead of the live option store. Used to rebuild a spec under the options
## a recorded run was made with.
.spec_opt_source <- new.env(parent = emptyenv())
.spec_opt_source$snap <- NULL

.spec_opt <- function(name, fallback, mode = "always") {
  snap <- .spec_opt_source$snap
  reg_default <- .dynhr_option_registry[[name]]$default
  if (!is.null(snap)) {
    set <- attr(snap, "set")
    if (name %in% set || (identical(mode, "always") && name %in% names(snap)))
      return(snap[[name]])
    return(if (identical(mode, "always")) reg_default else fallback)
  }
  if (exists(name, envir = .dynhr_opts, inherits = FALSE))
    return(get(name, envir = .dynhr_opts))
  if (identical(mode, "always")) reg_default else fallback
}

.spec_with_snapshot <- function(snap, fn) {
  old <- .spec_opt_source$snap
  .spec_opt_source$snap <- snap
  on.exit(.spec_opt_source$snap <- old, add = TRUE)
  fn()
}

.spec_default <- function(def) {
  if (is.null(def$option)) return(def$default)
  .spec_opt(def$option, def$default, def$option_mode)
}

# ---------------------------------------------------------------------------
# Field checking and sub-spec construction
# ---------------------------------------------------------------------------

.spec_bad <- function(path, what, class = "dynhr_error_spec_invalid")
  .dynhr_abort("estimation spec: `", path, "` ", what, ".", class = class)

## Coerce one field value to its schema type and check it.
.spec_check_value <- function(v, def, path) {
  if (is.null(v)) {
    if (isTRUE(def$nullable)) return(NULL)
    .spec_bad(path, "must not be NULL")
  }
  scalar <- function() length(v) == 1L && !is.na(v)
  switch(def$type,
    lgl1 = {
      if (!is.logical(v) || !scalar()) .spec_bad(path, "must be TRUE or FALSE")
      v <- as.vector(v, "logical")
    },
    int1 = {
      if (!is.numeric(v) || !scalar() || !is.finite(v) || v != round(v) ||
          abs(v) > .Machine$integer.max)
        .spec_bad(path, "must be a whole number")
      v <- as.vector(v, "integer")
    },
    num1 = {
      if (!is.numeric(v) || !scalar()) .spec_bad(path, "must be one number")
      v <- as.vector(v, "double")
    },
    chr1 = {
      if (!is.character(v) || !scalar()) .spec_bad(path, "must be one string")
      v <- as.vector(v, "character")
    },
    num = {
      if (!is.numeric(v) || !length(v) || anyNA(v) || is.object(v))
        .spec_bad(path, "must be a numeric vector without NA")
      storage.mode(v) <- "double"
    },
    list = {
      if (!is.list(v) || is.object(v))
        .spec_bad(path, "must be a plain list")
      nm <- names(v)
      if (length(v) && (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)))
        .spec_bad(path, "must be a list with unique, non-empty names")
      if (!length(v)) v <- list()   # canonical empty list (no names attribute)
    },
    any = NULL)
  if (!is.null(def$choices) && !v %in% def$choices)
    .spec_bad(path, paste0("must be one of ",
                           paste0("\"", def$choices, "\"", collapse = ", "),
                           " (got \"", v, "\")"))
  if (is.function(def$check)) {
    msg <- def$check(v)
    if (!is.null(msg)) .spec_bad(path, msg)
  }
  v
}

## Build (or, with `base`, update) one sub-spec from a named list of field
## values. A field absent from `fields` keeps its `base` value, else its
## default; a NULL value for a non-nullable field means "the default".
.spec_build <- function(component, fields, method = NULL, base = NULL) {
  fields <- as.list(fields)
  nms <- names(fields)
  if (length(fields) && (is.null(nms) || any(!nzchar(nms))))
    .dynhr_abort("estimation spec: every field of ", component,
                 "_spec() must be named.", class = "dynhr_error_spec_invalid")
  if (anyDuplicated(nms))
    .dynhr_abort("estimation spec: duplicated field(s) ",
                 paste(unique(nms[duplicated(nms)]), collapse = ", "),
                 " in ", component, "_spec().",
                 class = "dynhr_error_spec_invalid")
  if (identical(component, "sampler")) {
    method <- tolower(method %||% fields$method %||% base$method %||% "rwmh")
    if (!is.character(method) || length(method) != 1L ||
        !method %in% .spec_samplers)
      .dynhr_abort("estimation spec: unknown sampler method \"",
                   paste(method, collapse = " "), "\". Valid: ",
                   paste(.spec_samplers, collapse = ", "), ".",
                   class = "dynhr_error_spec_invalid")
    fields$method <- method
    nms <- names(fields)
  }
  sch <- .spec_component_schema(component, method)
  unknown <- setdiff(nms, names(sch))
  if (length(unknown)) {
    where <- if (identical(component, "sampler"))
      sprintf("sampler \"%s\"", method) else paste0(component, "_spec()")
    other <- if (identical(component, "sampler")) {
      takes <- names(Filter(function(f) any(unknown %in% f), .spec_sampler_table))
      if (length(takes)) paste0(" (",
                                paste(unknown[unknown %in% unlist(.spec_sampler_table)],
                                      collapse = ", "),
                                " belong(s) to: ", paste(takes, collapse = ", "), ")")
    }
    .dynhr_abort("estimation spec: unknown or irrelevant field(s) for ", where,
                 ": ", paste(unknown, collapse = ", "), other,
                 ". Valid fields: ", paste(names(sch), collapse = ", "), ".",
                 class = "dynhr_error_spec_unknown_field")
  }
  out <- list()
  for (nm in names(sch)) {
    def <- sch[[nm]]
    if (nm %in% nms) {
      v <- fields[[nm]]
      if (is.null(v) && !isTRUE(def$nullable)) v <- .spec_default(def)
    } else if (!is.null(base) && nm %in% names(base)) {
      v <- base[[nm]]
    } else {
      v <- .spec_default(def)
    }
    out[nm] <- list(.spec_check_value(v, def, paste0(component, "$", nm)))
  }
  class(out) <- .spec_class(component)
  if (identical(component, "likelihood")) .spec_check_likelihood(out)
  out
}

## Within-likelihood checks: estimation_context()'s own validation (the
## likelihood spec is a thin layer over it) plus the run_full_estimation /
## run_mode_finding argument conflicts.
.spec_check_likelihood <- function(x) {
  if (!is.null(x$plan) && !is.null(x$filter_tunes))
    .dynhr_abort("estimation spec: supply either likelihood$plan or ",
                 "likelihood$filter_tunes, not both.",
                 class = "dynhr_error_spec_invalid")
  if (!is.null(x$plan) && !is.null(x$heteroskedastic_shocks))
    .dynhr_abort("estimation spec: supply either likelihood$plan or ",
                 "likelihood$heteroskedastic_shocks, not both.",
                 class = "dynhr_error_spec_invalid")
  if (x$power_posterior > 1)
    .dynhr_warn("estimation spec: likelihood$power_posterior > 1 produces a ",
                "'cold' (over-confident) posterior. This is valid but ",
                "unusual; set it <= 1 for standard generalised-Bayes tempering.",
                class = "dynhr_warning_spec_cold_posterior")
  .spec_likelihood_context(x)
  invisible(x)
}

## The estimation_context() a likelihood spec describes (me_extra and
## shock_scale are resolved later, against the model and data).
.spec_likelihood_context <- function(x) {
  estimation_context(
    me_variance     = x$me_variance,
    likelihood      = x$type,
    lik_init        = x$lik_init,
    filter_method   = x$filter_method,
    singular_obs    = x$singular_obs,
    freq_band       = x$freq_band,
    system_priors   = x$system_priors,
    tpf_options     = x$tpf_options,
    obc_options     = x$obc_options,
    gradient_policy = x$gradient_policy,
    plan            = x$plan,
    sample_start    = x$sample_start,
    ms_spec         = x$ms_spec,
    ms_struct_spec  = x$ms_struct_spec,
    ms_collapse     = x$ms_collapse,
    student_df      = x$student_df,
    pruned_order    = x$pruned_order)
}

.spec_rd_escape <- function(s) gsub("([%{}\\\\])", "\\\\\\1", s)

.spec_rd_default <- function(def) {
  shown <- function(v) paste0("\\code{", .spec_rd_escape(
    paste(deparse(v, width.cutoff = 500L), collapse = " ")), "}")
  if (is.null(def$option)) return(shown(def$default))
  if (identical(def$option_mode, "set"))
    return(paste0("the \\code{", def$option, "} option when set, else ",
                  shown(def$default)))
  paste0("the \\code{", def$option, "} option (registered default ",
         shown(def$default), ")")
}

## Rd for a component's fields, from the schema (used through @eval).
.spec_rd_fields <- function(component) {
  sch <- .spec_schema[[component]]
  if (identical(component, "sampler")) sch <- sch[setdiff(names(sch), "method")]
  items <- vapply(names(sch), function(nm) {
    def <- sch[[nm]]
    ch <- if (!is.null(def$choices) && length(def$choices) > 1L)
      paste0(" One of ", paste0("\\code{\"", def$choices, "\"}", collapse = ", "),
             ".") else ""
    sprintf("  \\item{\\code{%s}}{%s.%s Default: %s.}", nm,
            .spec_rd_escape(def$doc), ch, .spec_rd_default(def))
  }, character(1))
  extra <- if (identical(component, "sampler")) {
    rows <- vapply(names(.spec_sampler_table), function(m)
      sprintf("  \\item{\\code{\"%s\"}}{%s%s}", m,
              paste0("\\code{", .spec_sampler_table[[m]], "}", collapse = ", "),
              if (!is.null(.spec_sampler_metrics[[m]]))
                paste0("; \\code{metric} one of ",
                       paste0("\\code{\"", .spec_sampler_metrics[[m]], "\"}",
                              collapse = ", ")) else ""),
      character(1))
    c("@section Fields by method:",
      "A method takes only the fields listed for it; any other field is an error.",
      "\\describe{", rows, "}")
  }
  c("@section Fields:", "\\describe{", items, "}", extra)
}

# ---------------------------------------------------------------------------
# Sub-spec constructors
# ---------------------------------------------------------------------------

#' Likelihood specification
#'
#' The likelihood part of a \code{\link{dynhr_estimation_spec}}: which
#' likelihood is evaluated, on which rows of the data, with which measurement
#' error, filter overrides and occasionally-binding-constraint treatment. It is
#' a thin layer over \code{\link{estimation_context}} (whose checks it runs)
#' that also holds the likelihood inputs that are otherwise only arguments of
#' \code{\link{run_full_estimation}} / \code{\link{run_mode_finding}}.
#'
#' Fields whose default names an option read that option when the spec is
#' built (see \code{\link{dynhr_set_options}}): a later change of the option
#' does not change an existing spec.
#'
#' @param type Likelihood type (a field; see Fields).
#' @param ... Further fields, by name. An unknown field is an error.
#' @eval .spec_rd_fields("likelihood")
#' @return A \code{dynhr_likelihood_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}, \code{\link{mode_spec}},
#'   \code{\link{sampler_spec}}, \code{\link{compute_spec}},
#'   \code{\link{outputs_spec}}
#' @examples
#' likelihood_spec("gaussian", me_variance = 1e-4)
#' @export
likelihood_spec <- function(type = "gaussian", ...) {
  .spec_build("likelihood", c(list(type = type), list(...)))
}

#' Mode-finding specification
#'
#' The posterior-mode stage of a \code{\link{dynhr_estimation_spec}}.
#' Fields whose default names an option read that option when the spec is
#' built.
#'
#' @param ... Fields, by name. An unknown field is an error.
#' @eval .spec_rd_fields("mode")
#' @return A \code{dynhr_mode_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @examples
#' mode_spec(method = "cmaes_newrat", n_iter = 2000L)
#' @export
mode_spec <- function(...) .spec_build("mode", list(...))

#' Sampler specification
#'
#' The posterior-sampling stage of a \code{\link{dynhr_estimation_spec}}.
#' Each method takes its own set of fields (see \emph{Fields by method}); a
#' field the method does not use is an error rather than silently ignored.
#' \code{"pmmh"} is random-walk Metropolis-Hastings over an unbiased particle
#' likelihood (see \code{\link{pmmh}}).
#'
#' @param method Sampler: \code{"rwmh"}, \code{"pmmh"}, \code{"nuts"},
#'   \code{"hmc"}, \code{"mala"}, \code{"chees"}, \code{"smc"},
#'   \code{"dsmh"}, \code{"dime"} or \code{"smc2"} (case-insensitive).
#' @param ... Fields of that method, by name.
#' @eval .spec_rd_fields("sampler")
#' @return A \code{dynhr_sampler_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @examples
#' sampler_spec("nuts", n_draws = 2000L, n_warmup = 1000L, metric = "lowrank")
#' @export
sampler_spec <- function(method = "rwmh", ...) {
  .spec_build("sampler", list(...), method = method)
}

#' Compute specification
#'
#' Execution settings of a \code{\link{dynhr_estimation_spec}}: parallelism,
#' the seed, progress output and checkpointing.
#'
#' @param ... Fields, by name. An unknown field is an error.
#' @eval .spec_rd_fields("compute")
#' @return A \code{dynhr_compute_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @examples
#' compute_spec(parallel = TRUE, n_cores = 4L, seed = 1L)
#' @export
compute_spec <- function(...) .spec_build("compute", list(...))

#' Outputs specification
#'
#' Post-estimation extras of a \code{\link{dynhr_estimation_spec}}: saved
#' files, diagnostics, stochastic simulation, the OBC smoother and the Ramsey
#' step. None of them changes a draw.
#'
#' @param ... Fields, by name. An unknown field is an error.
#' @eval .spec_rd_fields("outputs")
#' @return A \code{dynhr_outputs_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @examples
#' outputs_spec(diagnostics = TRUE, save = TRUE, dir = tempdir())
#' @export
outputs_spec <- function(...) .spec_build("outputs", list(...))

#' @rdname likelihood_spec
#' @param object A sub-spec to update.
#' @details \code{update(object, ...)} changes the named fields and keeps
#'   every other field as it is (it does not re-read the options); a
#'   \code{NULL} value resets a non-nullable field to its default. For a
#'   sampler spec, a different \code{method} starts a fresh spec of that
#'   method from the given fields.
#' @importFrom stats update
#' @export
update.dynhr_subspec <- function(object, ...) {
  cp <- .spec_component_of(object)
  new <- list(...)
  if (identical(cp, "sampler") && !is.null(new$method) &&
      !identical(tolower(new$method), object$method))
    return(.spec_build("sampler", new, method = new$method))
  .spec_build(cp, new, method = object$method, base = object)
}

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

.spec_fmt_value <- function(v, width = 60L) {
  if (is.null(v)) return("NULL")
  if (is.function(v)) return("<function>")
  if (is.matrix(v)) return(sprintf("<%s matrix %d x %d>", typeof(v), nrow(v), ncol(v)))
  if (is.data.frame(v)) return(sprintf("<data.frame %d x %d>", nrow(v), ncol(v)))
  if (is.object(v) && !is.atomic(v)) return(paste0("<", class(v)[1L], ">"))
  if (is.list(v)) {
    if (!length(v)) return("list()")
    s <- paste0("list(", paste(if (is.null(names(v))) seq_along(v) else names(v),
                                collapse = ", "), ")")
  } else if (is.atomic(v)) {
    if (!length(v)) return(paste0(typeof(v), "(0)"))
    el <- if (is.character(v)) paste0("\"", v, "\"") else format(v, digits = 7L)
    if (!is.null(names(v))) el <- paste0(names(v), "=", el)
    s <- paste(el, collapse = ", ")
    if (length(v) > 1L) s <- paste0("c(", s, ")")
  } else {
    s <- paste0("<", typeof(v), ">")
  }
  if (nchar(s) > width) paste0(substr(s, 1L, width - 3L), "...") else s
}

#' @rdname likelihood_spec
#' @param x A sub-spec.
#' @export
format.dynhr_subspec <- function(x, ...) {
  cp <- .spec_component_of(x)
  w  <- max(nchar(names(x)))
  c(sprintf("<%s_spec>%s", cp,
            if (identical(cp, "sampler")) paste0(" ", x$method) else ""),
    sprintf("  %-*s : %s", w, names(x), vapply(x, .spec_fmt_value, character(1))))
}

#' @rdname likelihood_spec
#' @export
print.dynhr_subspec <- function(x, ...) {
  writeLines(format(x, ...))
  invisible(x)
}

#' @rdname likelihood_spec
#' @export
as.list.dynhr_subspec <- function(x, ...) unclass(x)[names(x)]

# ---------------------------------------------------------------------------
# The estimation spec
# ---------------------------------------------------------------------------

## Model part: list(mod, path, prior_spec, max_order, compiled, solved).
## `compiled` / `solved` are caches (RDS only; not hashed); the source path is
## the parsed model's `source_file` when that file exists.
.spec_model_part <- function(model, compiled = NULL, max_order = 1L,
                             prior_spec = NULL) {
  solved <- NULL
  if (inherits(model, "dynhr_model")) {
    compiled   <- compiled %||% model$compiled
    prior_spec <- prior_spec %||% model$prior_spec
    model      <- model$model
  } else if (inherits(model, "dynhr_solved") ||
             (is.list(model) && !is.object(model) &&
              inherits(model$model, "dynhr_mod"))) {
    ## a dynhr_solved, or a bare list(model =, compiled =) standing in for one
    solved   <- model
    compiled <- compiled %||% model$compiled
    model    <- model$model
  } else if (is.character(model) && length(model) == 1L) {
    if (!file.exists(model))
      .dynhr_abort("estimation spec: model file not found: ", model,
                   class = "dynhr_error_spec_missing_file")
    model <- parse_mod(normalizePath(model), verbose = FALSE)
  }
  if (!inherits(model, "dynhr_mod"))
    .dynhr_abort("estimation spec: `model` must be a .mod path, a parsed ",
                 "dynhr_mod, a dynhr_solved or a dynhr_model.",
                 class = "dynhr_error_spec_invalid")
  path <- model$source_file
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path))
    path <- NULL
  ## A prior spec identical to the model's own is not stored. (A model with
  ## no estimated_params block has no own prior spec to compare with.)
  ep <- model$estimated_params
  if (!is.null(prior_spec) && !is.null(ep) && NROW(ep) > 0L) {
    own <- extract_prior_spec(model, verbose = FALSE)
    if (identical(own, prior_spec)) prior_spec <- NULL
  }
  if (!is.null(compiled) && !inherits(compiled, "dynhr_compiled"))
    .dynhr_abort("estimation spec: `compiled` must be a dynhr_compiled ",
                 "object from compile_model().", class = "dynhr_error_spec_invalid")
  if (!is.null(compiled)) max_order <- compiled$max_order %||% max_order
  list(mod = model, path = path, prior_spec = prior_spec,
       max_order = .spec_check_value(max_order, .spec_f("int1", doc = ""),
                                     "model$max_order"),
       compiled = compiled, solved = solved)
}

## Data part: list(value = matrix | NULL, path = CSV path | NULL).
.spec_data_part <- function(data) {
  if (is.null(data)) return(list(value = NULL, path = NULL))
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data))
      .dynhr_abort("estimation spec: data file not found: ", data,
                   class = "dynhr_error_spec_missing_file")
    return(list(value = NULL, path = normalizePath(data)))
  }
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    .dynhr_abort("estimation spec: `data` must be a numeric T x n_obs matrix, ",
                 "a data.frame or a CSV path.", class = "dynhr_error_spec_invalid")
  list(value = data, path = NULL)
}

## Sub-spec argument: an object of the right class, or a named list of fields.
.spec_as_sub <- function(x, component) {
  if (inherits(x, paste0("dynhr_", component, "_spec"))) return(x)
  if (is.list(x) && !is.object(x))
    return(.spec_build(component, x, method = if (identical(component, "sampler"))
                         x$method %||% "rwmh"))
  .dynhr_abort("estimation spec: `", component, "` must be a ", component,
               "_spec() or a named list of its fields.",
               class = "dynhr_error_spec_invalid")
}

## Sampler argument: NULL / FALSE (no sampling), a sampler spec, a named list
## of fields, or an unnamed list of those (a sequence, run in order).
.spec_as_sampler <- function(x) {
  if (is.null(x) || isFALSE(x)) return(NULL)
  if (inherits(x, "dynhr_sampler_sequence")) return(x)
  if (is.list(x) && !is.object(x) && length(x) && is.null(names(x))) {
    seq_ <- lapply(x, .spec_as_sub, component = "sampler")
    if (length(seq_) == 1L) return(seq_[[1L]])
    return(structure(seq_, class = "dynhr_sampler_sequence"))
  }
  .spec_as_sub(x, "sampler")
}

## The samplers of a spec as a plain list (0, 1 or more specs).
.spec_sampler_list <- function(s) {
  if (is.null(s)) list()
  else if (inherits(s, "dynhr_sampler_sequence")) unclass(s)
  else list(s)
}

## Option snapshot: the non-result-changing registered options plus any set
## unregistered ones, from the live store (or the scoped snapshot).
.spec_options_snapshot <- function() {
  snap <- .spec_opt_source$snap
  eff  <- if (is.null(snap)) dynhr_get_options(effective = TRUE) else snap
  attr(eff, "set") <- NULL
  eff[setdiff(names(eff), .spec_result_options())]
}

## Apply named option values: result-changing options go to their typed
## fields, the rest into the snapshot. `parts` is the list of spec parts.
.spec_apply_options <- function(parts, opts) {
  if (is.null(opts)) return(parts)
  if (!is.list(opts) || (length(opts) &&
        (is.null(names(opts)) || any(!nzchar(names(opts))))))
    .dynhr_abort("estimation spec: `options` must be a named list.",
                 class = "dynhr_error_spec_invalid")
  homes <- .spec_option_homes()
  for (nm in names(opts)) {
    val <- opts[[nm]]
    if (nm %in% names(homes)) {
      applied <- FALSE
      for (h in homes[[nm]]) {
        cp <- sub("\\$.*$", "", h); fld <- sub("^.*\\$", "", h)
        if (identical(cp, "sampler")) {
          sl <- .spec_sampler_list(parts$sampler)
          for (i in seq_along(sl)) if (fld %in% names(sl[[i]])) {
            sl[[i]] <- .spec_build("sampler", setNames(list(val), fld),
                                   method = sl[[i]]$method, base = sl[[i]])
            applied <- TRUE
          }
          if (length(sl) == 1L) parts$sampler <- sl[[1L]]
          else if (length(sl) > 1L)
            parts$sampler <- structure(sl, class = "dynhr_sampler_sequence")
        } else {
          parts[[cp]] <- .spec_build(cp, setNames(list(val), fld),
                                     base = parts[[cp]])
          applied <- TRUE
        }
      }
      if (!applied)
        .dynhr_abort("estimation spec: option `", nm, "` is a typed field (",
                     paste(homes[[nm]], collapse = ", "), ") that this spec ",
                     "does not have (sampler: ",
                     paste(vapply(.spec_sampler_list(parts$sampler),
                                  `[[`, "", "method"), collapse = ", "),
                     ").", class = "dynhr_error_spec_unknown_field")
    } else if (nm %in% names(.dynhr_option_registry)) {
      parts$options[nm] <- list(if (is.null(val))
        .dynhr_option_registry[[nm]]$default else val)
    } else if (is.null(val)) {
      parts$options[[nm]] <- NULL
    } else {
      parts$options[nm] <- list(val)
    }
  }
  parts
}

# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------

## Canonical form for hashing: two identical() values must serialise to the
## same bytes, which they need not do as built -- attributes can be stored in
## a different order (a data.frame rebuilt from YAML), NA can carry another
## payload, -0 and 0 differ in their bits, and strings can carry different
## encoding marks. Closures, environments and language objects are kept as
## they are.
.spec_canon <- function(x) {
  if (is.null(x) || is.function(x) || is.environment(x) || isS4(x) ||
      is.language(x) || !typeof(x) %in% c("logical", "integer", "double",
                                           "character", "list"))
    return(x)
  at <- attributes(x)
  y <- x
  attributes(y) <- NULL
  if (is.list(y)) {
    for (i in seq_along(y)) if (!is.null(y[[i]])) y[i] <- list(.spec_canon(y[[i]]))
  } else if (is.double(y)) {
    y[which(is.na(y) & !is.nan(y))] <- NA_real_
    y[which(y == 0)] <- 0
  } else if (is.character(y)) {
    y <- enc2utf8(y)
  }
  if (length(at)) {
    at <- at[order(names(at), method = "radix")]
    attributes(y) <- lapply(at, .spec_canon)
  }
  y
}

## Content hash of an R object: sha256 (digest, Suggests) or md5, of the
## version-2 XDR serialisation WITHOUT its header (which carries the writing R
## version), so the hash depends on content only. Version 2 also writes ALTREP
## vectors expanded, so 1:3 and c(1L, 2L, 3L) hash alike.
.spec_hash <- function(x) {
  raw <- serialize(.spec_canon(x), NULL, xdr = TRUE, version = 2L)
  raw <- raw[-seq_len(14L)]
  if (requireNamespace("digest", quietly = TRUE))
    return(paste0("sha256:", digest::digest(raw, algo = "sha256",
                                            serialize = FALSE)))
  f <- tempfile("dynhr-spec-")
  on.exit(unlink(f), add = TRUE)
  writeBin(raw, f)
  paste0("md5:", unname(tools::md5sum(f)))
}

## Hash of a file's bytes (sha256 via digest, else md5), NA when absent.
.spec_file_hash <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE))
    return(paste0("sha256:", digest::digest(file = path, algo = "sha256")))
  paste0("md5:", unname(tools::md5sum(path)))
}

.spec_model_hash <- function(mp) {
  mod <- mp$mod
  mod$source_file <- NA_character_
  .spec_hash(list(mod = mod, prior_spec = mp$prior_spec, max_order = mp$max_order))
}

## A sub-spec as hashed: plain lists; a precomputed mode result is reduced to
## the pieces the sampler consumes (its closures are not content).
.spec_hashable <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "dynhr_sampler_sequence")) return(lapply(unclass(x), unclass))
  x <- unclass(x)
  ## mode$run at its default "auto" is left out, so a spec keeps the hash
  ## (and a checkpoint its resumability) it had before the field existed
  if (identical(x[["run"]], "auto")) x["run"] <- NULL
  ## likewise mode$analytic_h0 at its default TRUE
  if (isTRUE(x[["analytic_h0"]])) x["analytic_h0"] <- NULL
  ## likelihood$pskf_cdf is inert unless the likelihood is pskf: left out
  ## there, so non-PSKF specs keep the hash they had before the field existed
  if (!is.null(x[["type"]]) && !identical(x[["type"]], "pskf")) x["pskf_cdf"] <- NULL
  if (!is.null(x$result))
    x$result <- list(theta_mode = x$result$theta_mode,
                     Sigma_prop = x$result$Sigma_prop,
                     V_mode     = x$result$V_mode)
  x
}

.spec_hashes <- function(parts) {
  h <- list(
    model      = .spec_model_hash(parts$model),
    data       = .spec_hash(list(
      value = if (!is.null(parts$data$value)) .spec_hash(parts$data$value),
      file  = if (!is.null(parts$data$path)) .spec_file_hash(parts$data$path))),
    likelihood = .spec_hash(.spec_hashable(parts$likelihood)),
    mode       = .spec_hash(.spec_hashable(parts$mode)),
    sampler    = .spec_hash(.spec_hashable(parts$sampler)),
    compute    = .spec_hash(.spec_hashable(parts$compute)),
    outputs    = .spec_hash(.spec_hashable(parts$outputs)),
    options    = .spec_hash(parts$options))
  h$spec <- .spec_hash(c(list(spec_version = .spec_version,
                             obs_vars = parts$obs_vars), h))
  h
}

## Assemble parts into a spec: hashes, then validate_spec().
.spec_assemble <- function(parts, validate = TRUE) {
  spec <- structure(list(
    spec_version = .spec_version,
    model      = parts$model,
    data       = parts$data,
    obs_vars   = parts$obs_vars,
    likelihood = parts$likelihood,
    mode       = parts$mode,
    sampler    = parts$sampler,
    compute    = parts$compute,
    outputs    = parts$outputs,
    options    = parts$options,
    hashes     = .spec_hashes(parts)
  ), class = "dynhr_estimation_spec")
  if (validate) validate_spec(spec)
  spec
}

.spec_parts <- function(spec) unclass(spec)[c("model", "data", "obs_vars",
  "likelihood", "mode", "sampler", "compute", "outputs", "options")]

.spec_resolve_obs_vars <- function(obs_vars, mod) {
  if (is.null(obs_vars) || !length(obs_vars))
    obs_vars <- mod$obs_vars %||% mod$varobs_names
  if (is.null(obs_vars) || !length(obs_vars))
    .dynhr_abort("estimation spec: `obs_vars` not supplied and the model ",
                 "declares no `varobs`; pass obs_vars= or add a varobs line ",
                 "to the .mod.", class = "dynhr_error_spec_invalid")
  if (!is.character(obs_vars) || anyNA(obs_vars) || anyDuplicated(obs_vars))
    .dynhr_abort("estimation spec: `obs_vars` must be unique names.",
                 class = "dynhr_error_spec_invalid")
  as.vector(obs_vars, "character")
}

#' Estimation specification: one object that determines an estimation run
#'
#' Collects everything an estimation needs -- the model, the data, the
#' likelihood, the mode-finding and sampling stages, compute settings,
#' post-estimation outputs and the package options -- in one versioned,
#' validated object (\code{spec_version = 1}). It can be printed, compared
#' (\code{\link{diff_specs}}), edited (\code{update()}), and written to and
#' read from YAML, JSON or RDS (\code{\link{write_spec}},
#' \code{\link{read_spec}}).
#'
#' @section Options:
#' Every registered option that can change a result (see
#' \code{\link{dynhr_set_options}}) is a typed field of the sub-spec that uses
#' it, and its default is the option's value when that sub-spec is built:
#' \code{me_variance} and \code{power_posterior} (likelihood);
#' \code{mode_method}, \code{nm_maxit} (both only when set),
#' \code{transform_params}, \code{proposal_cov_method},
#' \code{use_exact_hessian}, \code{use_analytic_hess} (field
#' \code{analytic_h0}) and \code{perturb_scale} (mode);
#' \code{transform_params}, \code{nuts_adapt}, \code{grad_method},
#' \code{rwmh_adapt_cov}, \code{rwmh_n_blocks} and \code{monge_alpha}
#' (sampler); \code{seed_base} (compute \code{seed}). The other options
#' (diagnostics, guards and I/O) are stored in \code{spec$options}, a snapshot.
#' So a spec, once built, does not depend on the session's options.
#'
#' @section Hashes:
#' \code{spec$hashes} holds a content hash of each part (\code{model},
#' \code{data}, \code{likelihood}, \code{mode}, \code{sampler},
#' \code{compute}, \code{outputs}, \code{options}) and of the whole
#' (\code{spec}): sha256 when the suggested package \pkg{digest} is
#' installed, else md5, prefixed with the algorithm. The model hash ignores
#' where the model was read from; cached compiled objects are not hashed.
#'
#' @param model A \code{.mod} path, a parsed \code{dynhr_mod}
#'   (\code{\link{parse_mod}}), a \code{dynhr_solved}
#'   (\code{\link{solve_model}}) or a \code{\link{dynhr_model}}: the spec
#'   that object holds is returned, with the components given here applied
#'   as \code{update()} applies them (a missing argument keeps the held
#'   component; \code{sampler = NULL} removes the sampling stage).
#' @param data A \eqn{T \times n_{obs}} numeric matrix or data.frame, or a CSV
#'   path (read when the estimation runs).
#' @param obs_vars Observable names; \code{NULL} = the model's \code{varobs}.
#' @param likelihood A \code{\link{likelihood_spec}} or a named list of its
#'   fields.
#' @param mode A \code{\link{mode_spec}} or a named list of its fields.
#' @param sampler A \code{\link{sampler_spec}}, a named list of its fields
#'   (with \code{method}), an unnamed list of those (a sequence run in
#'   order), or \code{NULL} (mode finding only).
#' @param compute A \code{\link{compute_spec}} or a named list of its fields.
#' @param outputs A \code{\link{outputs_spec}} or a named list of its fields.
#' @param options \code{NULL} (default) snapshots the current package options.
#'   A named list sets options for this spec: result-changing options are
#'   applied to their typed fields (after the sub-specs are built), the others
#'   are stored in the snapshot.
#' @param max_order Derivative order to compile the model to (when no compiled
#'   model is supplied).
#' @return A \code{dynhr_estimation_spec}: a list with \code{spec_version},
#'   \code{model} (\code{mod}, \code{path}, \code{prior_spec},
#'   \code{max_order}, and the caches \code{compiled} and \code{solved}),
#'   \code{data} (\code{value} or \code{path}), \code{obs_vars}, the five
#'   sub-specs, \code{options} and \code{hashes}.
#' @seealso \code{\link{as_estimation_spec}}, \code{\link{validate_spec}},
#'   \code{\link{write_spec}}, \code{\link{diff_specs}}
#' @examples
#' mod <- system.file("extdata/models/nk_demo.mod", package = "dynhr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr")))
#' spec <- dynhr_estimation_spec(mod, data = Y,
#'   sampler = sampler_spec("nuts", n_draws = 1000L, n_warmup = 500L))
#' spec
#' spec2 <- update(spec, sampler = list(n_chains = 8L))
#' diff_specs(spec, spec2)
#' @export
dynhr_estimation_spec <- function(model, data = NULL, obs_vars = NULL,
                                  likelihood = likelihood_spec(),
                                  mode = mode_spec(),
                                  sampler = sampler_spec("rwmh"),
                                  compute = compute_spec(),
                                  outputs = outputs_spec(),
                                  options = NULL,
                                  max_order = 1L) {
  if (inherits(model, "dynhr_model")) {
    ## The dynhr_model HOLDS a spec (its single source of truth): return it,
    ## with the components given here replacing / editing its own as
    ## update() does (a missing argument keeps the held component).
    args <- list(model$spec, data = data, obs_vars = obs_vars,
                 options = options)
    if (!missing(likelihood)) args$likelihood <- likelihood
    if (!missing(mode))       args$mode       <- mode
    if (!missing(sampler))    args$sampler    <- if (is.null(sampler)) FALSE else sampler
    if (!missing(compute))    args$compute    <- compute
    if (!missing(outputs))    args$outputs    <- outputs
    return(do.call(update.dynhr_estimation_spec, args))
  }
  mp <- .spec_model_part(model, max_order = max_order)
  parts <- list(
    model      = mp,
    data       = .spec_data_part(data),
    obs_vars   = .spec_resolve_obs_vars(obs_vars, mp$mod),
    likelihood = .spec_as_sub(likelihood, "likelihood"),
    mode       = .spec_as_sub(mode, "mode"),
    sampler    = .spec_as_sampler(sampler),
    compute    = .spec_as_sub(compute, "compute"),
    outputs    = .spec_as_sub(outputs, "outputs"),
    options    = .spec_options_snapshot())
  parts <- .spec_apply_options(parts, options)
  .spec_assemble(parts)
}

#' @rdname dynhr_estimation_spec
#' @param object A \code{dynhr_estimation_spec}.
#' @param ... Not used (an error when non-empty).
#' @details \code{update(spec, ...)} edits a spec component by component:
#'   \code{likelihood}, \code{mode}, \code{compute}, \code{outputs} and
#'   \code{sampler} take a named list of fields to change (the other fields
#'   keep their values -- options are not re-read) or a complete sub-spec that
#'   replaces the component; a named \code{sampler} list is applied to every
#'   sampler of a sequence that has those fields, \code{sampler = FALSE}
#'   removes the sampling stage; \code{options} is merged into the snapshot,
#'   with result-changing options routed to their typed fields. The result is
#'   validated again.
#' @export
update.dynhr_estimation_spec <- function(object, model = NULL, data = NULL,
                                         obs_vars = NULL, likelihood = NULL,
                                         mode = NULL, sampler = NULL,
                                         compute = NULL, outputs = NULL,
                                         options = NULL, ...) {
  if (length(list(...)))
    .dynhr_abort("update.dynhr_estimation_spec: unknown argument(s) ",
                 paste(names(list(...)), collapse = ", "), ".",
                 class = "dynhr_error_spec_unknown_field")
  .spec_assemble(.spec_update_parts(object, model = model, data = data,
                                    obs_vars = obs_vars,
                                    likelihood = likelihood, mode = mode,
                                    sampler = sampler, compute = compute,
                                    outputs = outputs, options = options))
}

## The parts of `object` with update()'s edits applied (not yet assembled or
## validated; the dynhr_model holds a spec without data, which it assembles
## unvalidated). `lenient_obs`: an empty obs_vars on a model without varobs
## stays empty instead of erroring (a dynhr_model that has no data yet).
## `model_part`: a model part built by the caller (the dynhr_model compiles
## a new model eagerly), used in place of `model`.
.spec_update_parts <- function(object, model = NULL, data = NULL,
                               obs_vars = NULL, likelihood = NULL,
                               mode = NULL, sampler = NULL, compute = NULL,
                               outputs = NULL, options = NULL,
                               lenient_obs = FALSE, model_part = NULL) {
  parts <- .spec_parts(object)
  if (!is.null(model) || !is.null(model_part)) {
    parts$model <- model_part %||%
      .spec_model_part(model, max_order = parts$model$max_order)
    if (is.null(obs_vars)) obs_vars <- parts$obs_vars
  }
  if (!is.null(data)) parts$data <- .spec_data_part(data)
  if (isTRUE(lenient_obs) && !length(obs_vars) && !is.null(obs_vars) &&
      !length(parts$model$mod$obs_vars %||% parts$model$mod$varobs_names))
    obs_vars <- NULL
  if (!is.null(obs_vars))
    parts$obs_vars <- .spec_resolve_obs_vars(obs_vars, parts$model$mod)
  for (cp in c("likelihood", "mode", "compute", "outputs")) {
    new <- switch(cp, likelihood = likelihood, mode = mode, compute = compute,
                  outputs = outputs)
    if (is.null(new)) next
    if (inherits(new, paste0("dynhr_", cp, "_spec"))) {
      parts[[cp]] <- new
    } else if (is.list(new) && !is.object(new)) {
      parts[[cp]] <- .spec_build(cp, new, base = parts[[cp]])
    } else {
      .dynhr_abort("update: `", cp, "` must be a named list of fields or a ",
                   cp, "_spec().", class = "dynhr_error_spec_invalid")
    }
  }
  if (isFALSE(sampler)) {
    parts["sampler"] <- list(NULL)
  } else if (!is.null(sampler)) {
    named_fields <- is.list(sampler) && !is.object(sampler) &&
      length(sampler) && !is.null(names(sampler))
    cur <- .spec_sampler_list(parts$sampler)
    same_method <- is.null(sampler$method) ||
      (length(cur) == 1L && identical(tolower(sampler$method), cur[[1L]]$method))
    if (named_fields && length(cur) && same_method) {
      flds <- sampler[setdiff(names(sampler), "method")]
      if (length(cur) == 1L) {
        cur[[1L]] <- .spec_build("sampler", flds, method = cur[[1L]]$method,
                                 base = cur[[1L]])
      } else {
        unk <- setdiff(names(flds), unlist(lapply(cur, names), use.names = FALSE))
        if (length(unk))
          .dynhr_abort("update: sampler field(s) ", paste(unk, collapse = ", "),
                       " are not fields of any sampler in this spec (",
                       paste(vapply(cur, `[[`, "", "method"), collapse = ", "),
                       ").", class = "dynhr_error_spec_unknown_field")
        for (i in seq_along(cur)) {
          f_i <- flds[intersect(names(flds), names(cur[[i]]))]
          if (length(f_i))
            cur[[i]] <- .spec_build("sampler", f_i, method = cur[[i]]$method,
                                    base = cur[[i]])
        }
      }
      parts$sampler <- if (length(cur) == 1L) cur[[1L]] else
        structure(cur, class = "dynhr_sampler_sequence")
    } else {
      parts$sampler <- .spec_as_sampler(sampler)
    }
  }
  .spec_apply_options(parts, options)
}

# ---------------------------------------------------------------------------
# .spec_resume_problem(): compute$resume = TRUE with a checkpoint_dir that
# holds no checkpoint. A property of the machine, not of the spec, so it is
# checked when run_estimation() starts (before the mode stage), not in
# validate_spec(): a spec written with resume = TRUE must still load
# elsewhere.
.spec_resume_problem <- function(cmp) {
  if (!isTRUE(cmp$resume) || is.null(cmp$checkpoint_dir)) return(NULL)
  meta <- .ckpt_paths(cmp$checkpoint_dir)$meta
  if (file.exists(meta)) return(NULL)
  paste0("compute$resume = TRUE, but compute$checkpoint_dir \"",
         cmp$checkpoint_dir, "\" holds no checkpoint (no ", basename(meta),
         "): there is no run to resume. Run once with resume = FALSE first, ",
         "or point checkpoint_dir at an existing checkpoint.")
}

# validate_spec(): the single home of cross-field checks
# ---------------------------------------------------------------------------

.spec_has_mcp <- function(mod) {
  any(vapply(mod$equations, function(e) {
    tg <- e$tag
    !is.null(tg) && !is.na(tg) && grepl("^\\s*mcp\\s*=", tg, perl = TRUE)
  }, logical(1L)))
}

## Is the OBC filter active? likelihood$obc, else auto-detected mcp tags.
.spec_obc_active <- function(spec) {
  obc <- spec$likelihood$obc
  if (isTRUE(obc)) return(TRUE)
  if (isFALSE(obc)) return(FALSE)
  .spec_has_mcp(spec$model$mod)
}

## Data column names after the data_col_map renaming (NULL when unknown).
.spec_data_columns <- function(spec) {
  if (!is.null(spec$data$value)) return(colnames(spec$data$value))
  if (is.null(spec$data$path)) return(NULL)
  cols <- names(utils::read.csv(spec$data$path, nrows = 1L))
  map <- spec$likelihood$data_col_map
  for (mod_nm in names(map)) {
    dat_nm <- map[[mod_nm]]
    if (dat_nm %in% cols && !(mod_nm %in% cols)) cols[cols == dat_nm] <- mod_nm
  }
  cols
}

# ---------------------------------------------------------------------------
# `extra` keys that shadow a typed field
# ---------------------------------------------------------------------------
## The `extra` lists are forwarded as `...`: likelihood$extra and mode$extra
## to the log-posterior constructor, sampler$extra to the sampler. A key that
## names a typed field (or a result-changing option) under its own or another
## name would bypass that field: likelihood$extra = list(power = 0.5) tempers
## the run while the spec records power_posterior = 1 -- and the mode stage's
## gradient, built from the typed fields, targets another posterior.

## Argument names of the constructors that set a typed field under another
## name. "" = no single field holds it (the spec builds that value itself).
.spec_extra_aliases <- c(
  power       = "likelihood$power_posterior",
  likelihood  = "likelihood$type",
  me_extra    = "likelihood$filter_tunes",
  shock_scale = "likelihood$heteroskedastic_shocks",
  ctx         = "")

## The typed field an `extra` key of `component` shadows, as "component$field"
## ("" when only the spec itself may set it), or NULL. `lik_type`: the spec's
## likelihood type (`order` is the TPF's pruned order); `method`: the sampler.
.spec_extra_shadow <- function(key, component, lik_type = NULL, method = NULL) {
  if (key %in% names(.spec_extra_aliases)) return(.spec_extra_aliases[[key]])
  if (identical(key, "order") && identical(lik_type, "tpf") &&
      component %in% c("likelihood", "mode"))
    return("likelihood$pruned_order")
  own <- switch(component,
    likelihood = c(setdiff(names(.spec_schema$likelihood), "extra")),
    ## mode$extra also reaches the log-posterior constructor
    mode = c(setdiff(names(.spec_schema$mode), "extra"),
             setdiff(names(.spec_schema$likelihood), "extra")),
    sampler = c("method", setdiff(.spec_sampler_table[[method]], "extra")))
  if (key %in% own) {
    cp <- if (identical(component, "mode") &&
              !key %in% names(.spec_schema$mode)) "likelihood" else component
    return(paste0(cp, "$", key))
  }
  homes <- .spec_option_homes()
  if (key %in% .spec_result_options() && !is.null(homes[[key]])) {
    h <- homes[[key]]
    same <- h[startsWith(h, paste0(component, "$"))]
    return(if (length(same)) same[[1L]] else h[[1L]])
  }
  NULL
}

## Every shadowing key of a spec (or of its parts): a data.frame with the
## `extra` list it sits in, the key and the typed field it shadows.
.spec_extra_shadows <- function(x) {
  lik_type <- x$likelihood$type
  sl <- .spec_sampler_list(x$sampler)
  src <- c(list(list(where = "likelihood$extra", cp = "likelihood",
                     ex = x$likelihood$extra, method = NULL),
                list(where = "mode$extra", cp = "mode", ex = x$mode$extra,
                     method = NULL)),
           lapply(seq_along(sl), function(i) list(
             where = if (length(sl) > 1L) sprintf("sampler[[%d]]$extra", i)
                     else "sampler$extra",
             cp = "sampler", ex = sl[[i]]$extra, method = sl[[i]]$method)))
  rows <- lapply(src, function(s) {
    tg <- vapply(names(s$ex), function(k)
      .spec_extra_shadow(k, s$cp, lik_type, s$method) %||% NA_character_,
      character(1))
    keep <- !is.na(tg)
    data.frame(where = rep(s$where, sum(keep)), key = names(s$ex)[keep],
               target = unname(tg[keep]), stringsAsFactors = FALSE)
  })
  do.call(rbind, c(list(data.frame(where = character(0), key = character(0),
                                   target = character(0),
                                   stringsAsFactors = FALSE)), rows))
}

## Can field `fld` of sub-spec `s` (component `cp`) take `v` from an `extra`
## key? Yes when it holds its default or already `v`.
.spec_route_settable <- function(s, cp, fld, v) {
  def <- .spec_component_schema(cp, s$method)[[fld]]
  !is.null(def) && (identical(s[[fld]], v) ||
                      identical(s[[fld]], .spec_default(def)))
}

## Route one `extra` value to its typed field `tgt` ("component$field").
## `idx`: the samplers a sampler field may go to. Returns list(done, parts, sl).
.spec_route_one <- function(parts, sl, tgt, v, idx) {
  no <- list(done = FALSE, parts = parts, sl = sl)
  if (is.null(tgt) || !grepl("$", tgt, fixed = TRUE)) return(no)
  cp  <- sub("\\$.*$", "", tgt)
  fld <- sub("^.*\\$", "", tgt)
  if (identical(cp, "sampler")) {
    idx <- idx[vapply(idx, function(i) fld %in% names(sl[[i]]), logical(1))]
    if (!length(idx) || !all(vapply(idx, function(i)
      .spec_route_settable(sl[[i]], "sampler", fld, v), logical(1))))
      return(no)
    for (i in idx)
      sl[[i]] <- .spec_build("sampler", stats::setNames(list(v), fld),
                             method = sl[[i]]$method, base = sl[[i]])
  } else {
    if (!cp %in% c("likelihood", "mode", "compute") ||
        !.spec_route_settable(parts[[cp]], cp, fld, v))
      return(no)
    parts[[cp]] <- .spec_build(cp, stats::setNames(list(v), fld),
                               base = parts[[cp]])
  }
  list(done = TRUE, parts = parts, sl = sl)
}

## as_estimation_spec() / the dynhr_model bridge: move each `extra` key that
## shadows a typed field to that field. The field takes the key's value when
## it still holds its default (or the same value); a key whose field was set to
## something else, or that no field holds (ctx, me_extra, shock_scale), stays
## in `extra` for validate_spec() to reject with a pointer to the field.
.spec_route_extras <- function(parts) {
  lik_type <- parts$likelihood$type
  sl <- .spec_sampler_list(parts$sampler)
  for (cp_src in c("likelihood", "mode")) {
    ex <- parts[[cp_src]]$extra
    moved <- character(0)
    for (key in names(ex)) {
      r <- .spec_route_one(parts, sl, .spec_extra_shadow(key, cp_src, lik_type),
                           ex[[key]], idx = seq_along(sl))
      if (r$done) {
        parts <- r$parts; sl <- r$sl; moved <- c(moved, key)
      }
    }
    if (length(moved))
      parts[[cp_src]] <- .spec_build(cp_src, list(extra = ex[setdiff(names(ex), moved)]),
                                     base = parts[[cp_src]])
  }
  for (i in seq_along(sl)) {
    ex <- sl[[i]]$extra
    moved <- character(0)
    for (key in names(ex)) {
      r <- .spec_route_one(parts, sl,
                           .spec_extra_shadow(key, "sampler", lik_type, sl[[i]]$method),
                           ex[[key]], idx = i)
      if (r$done) {
        parts <- r$parts; sl <- r$sl; moved <- c(moved, key)
      }
    }
    if (length(moved))
      sl[[i]] <- .spec_build("sampler", list(extra = ex[setdiff(names(ex), moved)]),
                             method = sl[[i]]$method, base = sl[[i]])
  }
  if (length(sl) == 1L) parts$sampler <- sl[[1L]]
  else if (length(sl) > 1L)
    parts$sampler <- structure(sl, class = "dynhr_sampler_sequence")
  parts
}

## The prior-initialised samplers: they draw their starting cloud / ensemble
## from the prior and use neither the mode nor its proposal covariance.
.spec_prior_init_samplers <- c("smc", "dsmh", "dime", "smc2")

## Does the runner run the mode stage? mode$run = "always" / "never" says so
## outright; "auto" (the default) runs it unless .spec_mode_stage_uses()
## finds nothing that uses the mode. A precomputed mode$result: no stage.
.spec_mode_stage_needed <- function(spec) {
  if (!is.null(spec$mode$result)) return(FALSE)
  switch(spec$mode$run %||% "auto",
         always = TRUE,
         never  = FALSE,
         length(.spec_mode_stage_uses(spec)) > 0L)
}

## What in `spec` uses the posterior mode (character; empty when nothing
## does): no sampler (the run is the mode stage), outputs$form = "mode", a
## sampler that starts from the mode (all but the prior-initialised SMC /
## DSMH / DIME / SMC2), or an output evaluated at the mode (diagnostics,
## Ramsey, the OBC smoother).
.spec_mode_stage_uses <- function(spec) {
  sl <- .spec_sampler_list(spec$sampler)
  if (!length(sl)) return("no sampler (the run is the mode stage)")
  out <- spec$outputs
  meth <- vapply(sl, `[[`, character(1), "method")
  meth <- unique(meth[!meth %in% .spec_prior_init_samplers])
  c(if (identical(out$form, "mode")) "outputs$form = \"mode\"",
    if (length(meth)) paste0("sampler \"", meth, "\""),
    if (isTRUE(out$diagnostics)) "outputs$diagnostics",
    if (isTRUE(out$ramsey)) "outputs$ramsey",
    if (isTRUE(out$smoother) && .spec_obc_active(spec)) "outputs$smoother")
}

## Run a runtime rule that signals its own error and re-signal that error as a
## spec error (same message), so validate_spec()'s errors keep one class.
.spec_reclass <- function(expr, class = NULL)
  withCallingHandlers(expr, error = function(e)
    .dynhr_abort("validate_spec: ", conditionMessage(e),
                 class = c(class, "dynhr_error_spec_invalid")))

## Can filter_tunes add observables to this spec's run? (They enlarge the
## observable set the likelihood sees, so a per-observable length cannot be
## checked against spec$obs_vars then.)
.spec_tunes_possible <- function(spec) {
  lik <- spec$likelihood
  if (!is.null(lik$plan)) return(TRUE)
  ft <- lik$filter_tunes
  if (isFALSE(ft)) return(FALSE)
  if (!is.null(ft)) return(TRUE)
  tn <- spec$model$mod$filter_tunes
  df <- if (is.data.frame(tn)) tn else if (is.list(tn)) tn$tunes
  is.data.frame(df) && nrow(df) > 0L
}

## The planner objective text the Ramsey workflow works from: the parsed
## `planner_objective(...)` of the model ("" when it has none).
.spec_planner_objective <- function(mod)
  mod$planner_objective$text %||% ""

## The prior spec the mode stage uses (the spec's own, else the model's);
## NULL when the model estimates nothing.
.spec_prior_spec <- function(spec) {
  if (!is.null(spec$model$prior_spec)) return(spec$model$prior_spec)
  ep <- spec$model$mod$estimated_params
  if (is.null(ep) || NROW(ep) == 0L) return(NULL)
  extract_prior_spec(spec$model$mod, verbose = FALSE)
}

## ---- sampler arguments, as the runner passes them -----------------------
## The sampler function a spec method calls, so its arguments can be judged
## by the sampler's own checker (and the names it gets are checked against the
## function's formals). pmmh is RWMH over a particle likelihood; smc runs on a
## mirai pool when compute$parallel is set.
.spec_sampler_fn <- function(method, parallel = FALSE)
  switch(method,
         rwmh = , pmmh = rwmh,
         nuts = dynhr_nuts, hmc = dynhr_hmc, mala = dynhr_mala,
         chees = dynhr_chees, dime = run_dime,
         smc = if (isTRUE(parallel)) run_smc_mirai else dynhr_smc,
         dsmh = dynhr_dsmh, smc2 = dynhr_smc2)

## The function a sampler's PARALLEL (mirai) path calls when that path is
## not the serial sampler with more cores: multi-chain NUTS and DIME. NULL
## for the others (their parallel paths take the same arguments, or are
## covered by .spec_sampler_fn).
.spec_parallel_sampler_fn <- function(method)
  switch(method, nuts = run_nuts_mirai, dime = run_dime_mirai, NULL)

## sampler$extra names the parallel function cannot take. The runner forwards
## the others and refuses these (so they are never silently dropped);
## validate_spec() refuses them when the spec is built.
.spec_parallel_extra_unsupported <- function(method, extra) {
  fn <- .spec_parallel_sampler_fn(method)
  if (is.null(fn) || !length(extra) || is.null(names(extra))) return(character(0))
  setdiff(names(extra), names(formals(fn)))
}

## Does the runner take the parallel path for this sampler spec?
.spec_takes_parallel_path <- function(spec, samp)
  isTRUE(spec$compute$parallel) &&
    (identical(samp$method, "dime") ||
       (identical(samp$method, "nuts") && isTRUE(samp$n_chains > 1L)))

## NUTS takes these metrics by name (a fixed dense one arrives as a matrix).
.spec_nuts_adapted_metrics <- c("diagonal", "warmup_dense", "lowrank",
                                "fisher_diag")

## The named list of arguments the spec runner hands the sampler function of
## `samp` (a sampler_spec): the typed fields under the function's own argument
## names, then sampler$extra (which wins, as in the runner's do.call). Only
## what a spec states is listed; the runner's own additions (starting point,
## gradient, metric matrices, checkpoint, ...) are not judged here. rwmh's
## n_draws counts burn-in (the runner passes n_draws + n_warmup), DIME's
## walkers / iterations are run_dime's n_chain / n_iter, and NUTS takes only
## the adapted metrics by name.
.spec_sampler_call_args <- function(spec, samp) {
  m  <- samp$method
  nd <- samp$n_draws
  nb <- samp$n_warmup
  base <- switch(m,
    rwmh = , pmmh = list(n_draws = nd + nb, n_burn = nb,
                         adapt_cov = samp$adapt_cov, n_blocks = samp$n_blocks,
                         Sigma_prop = samp$Sigma_prop),
    nuts  = list(n_draws = nd, n_warmup = nb,
                 metric = if (samp$metric %in% .spec_nuts_adapted_metrics)
                            samp$metric),
    hmc   = list(n_draws = nd, n_warmup = nb, metric = samp$metric),
    mala  = , chees = list(n_draws = nd, n_warmup = nb),
    dime  = list(n_chain = samp$n_walkers, n_iter = nd, n_burn = nb),
    smc   = , dsmh = , smc2 = list(n_particles = samp$n_particles))
  base <- Filter(Negate(is.null), base)
  ex <- samp$extra
  if (!length(ex)) return(base)
  if (is.null(names(ex)) || any(!nzchar(names(ex)))) return(c(base, ex))
  base[names(ex)] <- ex
  base
}

## Every problem the sampler's own argument checker finds in the spec's
## sampler (character(0) if none); the condition classes to abort with are
## kept in the result's "class_names" attribute. `n_par` is the number of
## estimated parameters, NULL when unknown.
.spec_sampler_problems <- function(spec, samp, n_par = NULL) {
  m    <- samp$method
  args <- .spec_sampler_call_args(spec, samp)
  par  <- isTRUE(spec$compute$parallel)
  fn   <- .spec_sampler_fn(m, par)
  fn_name <- switch(m, rwmh = , pmmh = "rwmh", nuts = "dynhr_nuts",
                    hmc = "dynhr_hmc", mala = "dynhr_mala",
                    chees = "dynhr_chees", dime = "run_dime",
                    smc = if (par) "run_smc_mirai" else "dynhr_smc",
                    dsmh = "dynhr_dsmh", smc2 = "dynhr_smc2")
  ## names that are not arguments of the sampler function (the checkers of
  ## the SMC family do not look at names)
  unknown <- function() .mcmc_check_args(args, fn, fn_name, list(), n_par)
  p <- switch(m,
    rwmh  = .rwmh_args_problem(args, n_par),
    pmmh  = c(.rwmh_args_problem(args, n_par),
              .pmmh_args_problem(c(list(n_particles = samp$n_particles),
                                   args[intersect("methods", names(args))]))),
    nuts  = .nuts_args_problem(args, n_par),
    hmc   = .hmc_args_problem(args, n_par),
    mala  = .mala_args_problem(args, n_par),
    chees = .chees_args_problem(args, n_par),
    dime  = .dime_args_problem(args, n_par),
    smc   = c(unknown(), .smc_args_problem(args, n_par)),
    dsmh  = c(unknown(), .dsmh_args_problem(args, n_par)),
    smc2  = c(unknown(), .smc2_args_problem(args, n_par)))
  if (.spec_takes_parallel_path(spec, samp)) {
    uns <- .spec_parallel_extra_unsupported(m, samp$extra)
    if (length(uns))
      p <- c(p, stats::setNames(sprintf(paste0(
        "`%s` is not supported on the parallel (mirai) %s path, which would ",
        "ignore it; set compute$parallel = FALSE or drop it."),
        uns, toupper(m)), rep("dynhr_error_inapplicable_argument", length(uns))))
  }
  nm <- names(p)
  cls <- c(if (!is.null(nm)) nm[nzchar(nm)],
           switch(m, dsmh = "dynhr_error_dsmh_args",
                  smc = , smc2 = "dynhr_error_smc_args",
                  pmmh = NULL, "dynhr_error_invalid_argument"),
           if (any(grepl("^unknown argument", p))) "dynhr_error_unknown_argument")
  p <- unique(unname(p))
  attr(p, "class_names") <- unique(cls)
  p
}

## The nearest existing ancestor of `path` (path itself when it exists).
.spec_existing_ancestor <- function(path) {
  p <- path
  while (!file.exists(p)) {
    up <- dirname(p)
    if (identical(up, p)) return(p)
    p <- up
  }
  p
}

#' Validate an estimation spec
#'
#' The single place where an estimation spec's fields are checked against each
#' other. Each sub-spec is checked field by field when it is built; this adds
#' the checks that span components. Errors have class
#' \code{dynhr_error_spec_invalid}; warnings are classed
#' \code{dynhr_warning_spec_*} (and \code{dynhr_warning_metric_ignored}).
#'
#' Checks (E = error, W = warning):
#' \itemize{
#'   \item E: data present; observables resolved and present among the data
#'     columns (after \code{data_col_map}); the \code{first_obs} /
#'     \code{nobs} window inside the data.
#'   \item E: an OBC model (mcp tags, or \code{obc = TRUE}) with a likelihood
#'     other than \code{"gaussian"} (auto-switched to the OBC filter),
#'     \code{"pkf"}, \code{"ppf"} or \code{"copf"}; OBC-only likelihoods or a
#'     non-default \code{obc_filter} without OBC.
#'   \item E (class \code{dynhr_error_spec_shadowed_field}): a
#'     \code{likelihood$extra}, \code{mode$extra} or \code{sampler$extra} key
#'     that shadows a typed field or a result-changing option -- e.g.
#'     \code{power} or \code{power_posterior} (\code{likelihood$power_posterior}),
#'     \code{me_variance}, \code{grad_method}, \code{transform_params},
#'     \code{ctx} -- which would bypass the field. The error names the field
#'     to set. (\code{\link{as_estimation_spec}} and the
#'     \code{dynhr_model} bridge move such keys to their field.)
#'   \item W: mode finding on a noisy particle likelihood (\code{"tpf"},
#'     \code{"ppf"}, \code{"copf"}, \code{"global_pf"}) -- when the mode stage
#'     runs (it is skipped for prior-initialised samplers, see
#'     \code{\link{run_estimation}}). The spec runner does run it (the
#'     optimiser then works on a noisy objective); the flat
#'     \code{\link{run_mode_finding}} does not accept these likelihoods.
#'   \item E: arguments that would only fail after an expensive stage:
#'     \code{mode$method} not a mode-finding optimiser; unknown
#'     \code{mode$options} keys; \code{likelihood$extra} / \code{mode$extra}
#'     names the log-posterior constructor (and, for \code{mode$extra}, the
#'     optimiser) does not take; a whittle \code{freq_band} outside
#'     \eqn{[0, \pi]}; a \code{me_variance} vector of the wrong length (or a
#'     per-observable vector the likelihood cannot use);
#'     \code{pruned_order = 3} on a model compiled below \code{max_order = 2};
#'     \code{mode$theta_init} names or bounds the mode stage would refuse;
#'     \code{likelihood$dates} not one per sample row; \code{outputs$ramsey}
#'     on a model without a planner objective (or a free-instrument system, or
#'     \code{ramsey_order = 2} below \code{max_order = 2}); \code{outputs$save}
#'     with a \code{dir} that is (or sits under) a file, or a \code{prefix}
#'     with a path separator. W: \code{outputs$ramsey} without a discount
#'     factor (\code{beta} / \code{betta} or \code{ramsey_discount}), and
#'     \code{likelihood$dates} when diagnostics are off (it only labels them).
#'   \item E (class \code{dynhr_error_spec_unsupported}): a sampler the spec
#'     runner does not run (\code{"smc2"}; call \code{\link{dynhr_smc2}}
#'     directly).
#'   \item E (class \code{dynhr_error_spec_mode_needed}):
#'     \code{mode$run = "never"} while something uses the mode: no sampler,
#'     \code{outputs$form = "mode"}, a sampler other than the
#'     prior-initialised \code{"smc"}, \code{"dsmh"}, \code{"dime"},
#'     \code{"smc2"}, or an output evaluated at the mode
#'     (\code{outputs$diagnostics}, \code{outputs$ramsey}, the OBC
#'     \code{outputs$smoother}).
#'   \item E: \code{metric = "monge"} / \code{"whittle_fim"} unless the
#'     \code{allow_monge_metric} / \code{allow_whittle_fim_metric} option is
#'     \code{TRUE} in \code{spec$options}; W for \code{"monge"} when allowed.
#'   \item W: parallel multi-chain NUTS with a fixed metric (\code{"hessian"},
#'     \code{"whittle_fim"}), which that path does not use.
#'   \item E: parallel multi-chain NUTS with \code{adapt = "pooled"} and
#'     \code{metric = "warmup_dense"} or a checkpoint.
#'   \item E: \code{checkpoint_dir} with a sampler (or parallel path) that
#'     cannot stream to it; \code{resume = TRUE} without \code{checkpoint_dir}.
#'     (\code{resume = TRUE} on a directory that holds no checkpoint is a
#'     property of the machine, not the spec: \code{\link{run_estimation}}
#'     refuses it before any stage runs.)
#'   \item E: a sampler argument the sampler itself would refuse (or that
#'     would run a frozen chain, a collapsed particle cloud or a
#'     non-monotone tempering ladder): each sampler has one argument checker,
#'     applied here to the arguments the runner would pass it (the typed
#'     fields under the sampler function's own names, then
#'     \code{sampler$extra}) and listing every problem -- e.g. NUTS
#'     \code{max_treedepth = 0}, SMC \code{n_mh_steps} / \code{ess_target} /
#'     a \code{lambda_schedule} that is not increasing in (0, 1], RWMH
#'     \code{scale}, \code{target_rate} or a \code{Sigma_prop} that is not
#'     positive definite or not of the estimated parameters' order, fewer DIME
#'     walkers than parameters plus one, an unknown \code{sampler$extra} name.
#'     Condition classes: \code{dynhr_error_invalid_argument} (the MCMC
#'     samplers), \code{dynhr_error_smc_args},
#'     \code{dynhr_error_unknown_argument} (an unknown name). One rule stays
#'     at run time: a \code{"dsmh"} sampler needs \code{n_obs},
#'     \code{lambda1} or \code{lambda_schedule}, which the runner takes from
#'     the data.
#'   \item E (class \code{dynhr_error_dsmh_args}): a \code{"dsmh"} sampler
#'     whose \code{n_particles} is not a multiple of \code{n_groups} or
#'     \code{n_strata} (set through \code{sampler$extra}; defaults 10 and 20),
#'     or another \code{dynhr_dsmh()} argument rule -- caught when the spec is
#'     built, not after the mode stage.
#'   \item W: \code{"pmmh"} on a likelihood that is not an unbiased particle
#'     likelihood; E: \code{"smc2"} on a likelihood other than \code{"tpf"} /
#'     \code{"sv_rbpf"}.
#'   \item W: the OBC smoother or PPF re-weighting requested without OBC.
#'   \item E: an \code{outputs$form} the stages cannot produce:
#'     \code{"mode"} with a sampler or a precomputed \code{mode$result};
#'     \code{"posterior"} without a sampler; \code{"full"} with a sampler
#'     sequence or a precomputed \code{mode$result}.
#' }
#'
#' @param spec A \code{dynhr_estimation_spec}.
#' @return \code{spec}, invisibly.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @export
validate_spec <- function(spec) {
  bad <- function(...) .dynhr_abort("validate_spec: ", ...,
                                    class = "dynhr_error_spec_invalid")
  if (!inherits(spec, "dynhr_estimation_spec"))
    bad("`spec` is not a dynhr_estimation_spec.")
  if (!identical(spec$spec_version, .spec_version))
    .dynhr_abort("validate_spec: spec_version ", format(spec$spec_version),
                 " is not supported (this dynhr reads version ", .spec_version,
                 ").", class = "dynhr_error_spec_version")
  for (cp in c("likelihood", "mode", "compute", "outputs"))
    if (!inherits(spec[[cp]], paste0("dynhr_", cp, "_spec")))
      bad("`", cp, "` is not a ", cp, "_spec().")
  samplers <- .spec_sampler_list(spec$sampler)
  for (s in samplers)
    if (!inherits(s, "dynhr_sampler_spec")) bad("`sampler` is not a sampler_spec().")
  lik <- spec$likelihood
  cmp <- spec$compute
  out <- spec$outputs
  opts <- spec$options

  ## ---- data and observables (run_full_estimation / run_mode_finding) ----
  if (is.null(spec$data$value) && is.null(spec$data$path))
    bad("no data: supply `data` (a matrix, data.frame or CSV path).")
  if (!length(spec$obs_vars)) bad("no observables (`obs_vars`).")
  cols <- .spec_data_columns(spec)
  if (!is.null(cols)) {
    miss <- setdiff(spec$obs_vars, cols)
    if (length(miss))
      bad("obs_vars not found in the data columns: ",
          paste(miss, collapse = ", "), ".")
  }
  if (!is.null(spec$data$value)) {
    nr <- nrow(spec$data$value)
    last <- if (is.null(lik$nobs)) nr else lik$first_obs + lik$nobs - 1L
    if (lik$first_obs > nr || last > nr)
      bad("the likelihood sample (first_obs = ", lik$first_obs, ", nobs = ",
          format(lik$nobs %||% "all"), ") runs past the ", nr, " data rows.")
  }

  ## ---- OBC x likelihood (run_mode_finding's rule) ------------------------
  obc_on <- .spec_obc_active(spec)
  if (obc_on && !lik$type %in% c("gaussian", .spec_obc_likelihoods))
    bad("the model has occasionally-binding constraints (mcp tags or ",
        "likelihood$obc = TRUE), which only the OBC/PKF filter (or PPF / ",
        "COPF via likelihood$obc_filter) can evaluate; likelihood$type = \"",
        lik$type, "\" is incompatible with it. Use type = \"gaussian\" ",
        "(switched to the OBC filter) or set likelihood$obc = FALSE to ",
        "estimate the model as linear.")
  if (!obc_on && lik$type %in% .spec_obc_likelihoods)
    bad("likelihood$type = \"", lik$type, "\" is an OBC filter, but the ",
        "model has no occasionally-binding constraints (no mcp tags and ",
        "likelihood$obc is not TRUE).")
  if (!obc_on && !identical(lik$obc_filter, "pkf"))
    bad("likelihood$obc_filter = \"", lik$obc_filter, "\" needs an OBC model.")
  eff_type <- if (obc_on) lik$obc_filter else lik$type

  ## ---- `extra` keys that shadow a typed field --------------
  sh <- .spec_extra_shadows(spec)
  if (nrow(sh)) {
    what <- ifelse(nzchar(sh$target),
                   paste0("set ", sh$target, " instead"),
                   "the spec builds that value from its likelihood fields")
    .dynhr_abort("validate_spec: ",
                 paste0(sh$where, " key `", sh$key, "` shadows a typed field (",
                        what, ")", collapse = "; "),
                 ". An `extra` entry would bypass the field: the spec, its ",
                 "hash and the mode stage's gradient would describe another ",
                 "target than the run evaluates.",
                 class = c("dynhr_error_spec_shadowed_field",
                           "dynhr_error_spec_invalid"))
  }

  ## ---- mode$run = "never" with something that uses the mode --------------
  if (is.null(spec$mode$result) && identical(spec$mode$run, "never")) {
    uses <- .spec_mode_stage_uses(spec)
    if (length(uses))
      .dynhr_abort("validate_spec: mode$run = \"never\", but ",
                   paste(uses, collapse = ", "), " use(s) the posterior mode. ",
                   "Only the prior-initialised samplers (",
                   paste(.spec_prior_init_samplers, collapse = ", "),
                   ") run without it; set mode$run = \"auto\" or supply ",
                   "mode$result.",
                   class = c("dynhr_error_spec_mode_needed",
                             "dynhr_error_spec_invalid"))
  }

  ## ---- mode stage on a noisy likelihood ----------------------------------
  if (is.null(spec$mode$result) && eff_type %in% .spec_noisy_mode_likelihoods &&
      .spec_mode_stage_needed(spec))
    .dynhr_warn("validate_spec: mode finding on the noisy particle likelihood \"",
                eff_type, "\": the deterministic optimisers assume a fixed ",
                "objective, but each evaluation here is a fresh noisy estimate, ",
                "so the mode and its curvature are noisy too (the flat ",
                "run_mode_finding() does not accept this likelihood). Prefer a ",
                "deterministic likelihood for the mode, or a prior-initialised ",
                "sampler (SMC, DSMH, DIME) with mode$run = \"never\".",
                class = "dynhr_warning_spec_noisy_mode")

  ## ---- per sampler -------------------------------------------------------
  ## analytic_grad (default TRUE) means "the exact gradient wherever the
  ## likelihood has one". Where it has none (OBC, or a likelihood / gradient
  ## policy / per-observable me_variance without one) the runner keeps the
  ## numerical gradient without a warning -- the default must stay quiet, and
  ## TRUE cannot be told from the default -- and its verbose output names the
  ## gradient each stage used.
  par <- isTRUE(cmp$parallel)
  ## the number of estimated parameters, for the rules that need it
  n_par_est <- if (length(samplers)) {
    pr <- .spec_prior_spec(spec)
    if (!is.null(pr)) length(pr$name)
  }
  for (s in samplers) {
    m <- s$method
    if (identical(s$metric, "monge")) {
      if (!isTRUE(opts$allow_monge_metric))
        bad("metric = \"monge\" is experimental and disabled by default: the ",
            "position-dependent Monge metric collapses the proposal step on ",
            "sharp / near-unit-root posteriors. Prefer \"hessian\" or ",
            "\"diagonal\"; to use it anyway set options = ",
            "list(allow_monge_metric = TRUE).")
      .dynhr_warn("validate_spec: metric = \"monge\" can collapse the proposal ",
                  "step on sharp / near-unit-root posteriors; check mixing / ESS.",
                  class = "dynhr_warning_spec_monge_metric")
    }
    if (identical(s$metric, "whittle_fim") && !isTRUE(opts$allow_whittle_fim_metric))
      bad("metric = \"whittle_fim\" is experimental and disabled by default; ",
          "to use it anyway set options = list(allow_whittle_fim_metric = TRUE).")
    par_nuts <- identical(m, "nuts") && par && s$n_chains > 1L
    if (par_nuts && s$metric %in% c("hessian", "whittle_fim"))
      .dynhr_warn(sprintf(paste0(
        "validate_spec: metric = \"%s\" is not supported on the parallel ",
        "(mirai) NUTS path; the chains use the adapted \"diagonal\" metric. ",
        "Set compute$parallel = FALSE to use \"%s\"."), s$metric, s$metric),
        class = "dynhr_warning_metric_ignored")
    if (par_nuts && identical(s$adapt, "pooled")) {
      if (identical(s$metric, "warmup_dense"))
        bad("adapt = \"pooled\" supports metric = \"diagonal\", ",
            "\"fisher_diag\" or \"lowrank\", not \"warmup_dense\".")
      if (!is.null(cmp$checkpoint_dir))
        bad("adapt = \"pooled\" does not support a checkpoint (compute$checkpoint_dir).")
    }
    if (!is.null(cmp$checkpoint_dir)) {
      ok <- if (par && !m %in% c("mala", "chees", "hmc")) m %in% .spec_checkpoint_parallel
            else m %in% .spec_checkpoint_serial
      if (!ok)
        bad("compute$checkpoint_dir is set, but sampler \"", m, "\"",
            if (par) " on the parallel path", " cannot stream to a checkpoint ",
            "(serial: ", paste(.spec_checkpoint_serial, collapse = ", "),
            "; parallel: ", paste(.spec_checkpoint_parallel, collapse = ", "), ").")
    }
    if (identical(m, "pmmh") && !eff_type %in% .spec_pmmh_likelihoods)
      .dynhr_warn("validate_spec: sampler \"pmmh\" on likelihood \"", eff_type,
                  "\", which is not an unbiased particle likelihood (",
                  paste(.spec_pmmh_likelihoods, collapse = ", "), "); this is ",
                  "plain random-walk Metropolis-Hastings.",
                  class = "dynhr_warning_spec_pmmh_likelihood")
    ## the sampler's own argument rules (one checker per sampler, shared with
    ## the sampler's entry), applied to the arguments the runner would pass
    sp <- .spec_sampler_problems(spec, s, n_par_est)
    if (length(sp))
      .dynhr_abort("validate_spec: sampler \"", m, "\": ",
                   paste(sp, collapse = " "),
                   if (identical(m, "dsmh"))
                     " (n_groups / n_strata are set through sampler$extra.)",
                   class = c(attr(sp, "class_names"),
                             "dynhr_error_spec_invalid"))
    if (identical(m, "smc2") && !lik$type %in% c("tpf", "sv_rbpf"))
      bad("sampler \"smc2\" needs likelihood$type \"tpf\" or \"sv_rbpf\" ",
          "(got \"", lik$type, "\").")
    if (m %in% .spec_unsupported_samplers)
      .dynhr_abort("validate_spec: ", .spec_unsupported_sampler_msg(m),
                   class = c("dynhr_error_spec_unsupported",
                             "dynhr_error_spec_invalid"))
  }
  if (isTRUE(cmp$resume) && is.null(cmp$checkpoint_dir))
    bad("compute$resume = TRUE needs compute$checkpoint_dir.")

  ## ---- fields whose rule needs the model, the data or the file system ----
  mode_runs <- is.null(spec$mode$result) && .spec_mode_stage_needed(spec)
  n_rows <- if (!is.null(spec$data$value)) {
    nr <- nrow(spec$data$value)
    last <- if (is.null(lik$nobs)) nr else lik$first_obs + lik$nobs - 1L
    last - lik$first_obs + 1L
  }

  ## mode$method: one of the optimisers the dispatcher implements
  if (!spec$mode$method %in% .mode_finding_methods)
    bad("mode$method = \"", spec$mode$method, "\" is not a mode-finding ",
        "optimiser. Valid: ", paste(.mode_finding_methods, collapse = ", "), ".")

  ## mode$options: the mode stage reads a fixed set of keys, so another is a typo
  bad_opt <- setdiff(names(spec$mode$options), .mode_option_keys)
  if (length(bad_opt))
    bad("mode$options has key(s) the mode stage does not read: ",
        paste0("`", bad_opt, "`", collapse = ", "), ". Valid: ",
        paste(.mode_option_keys, collapse = ", "), ".")

  ## likelihood$extra / mode$extra: the log-posterior constructor's `...` takes
  ## the names of some likelihood; mode$extra also takes the optimiser's
  ## arguments. (OBC models build another constructor and ignore both.)
  if (!obc_on) {
    known <- unique(unlist(.mlp_dots_by_likelihood(), use.names = FALSE))
    ## (the constructor's own error class is kept, so a caller that handles
    ## dynhr_error_unknown_argument sees the same condition earlier)
    bad_arg <- function(...)
      .dynhr_abort("validate_spec: ", ...,
                   class = c("dynhr_error_unknown_argument",
                             "dynhr_error_spec_invalid"))
    lik_bad <- setdiff(names(lik$extra), c(known, "verbose"))
    if (length(lik_bad))
      bad_arg("likelihood$extra has argument(s) the log-posterior constructor ",
          "does not take: ", paste0("`", lik_bad, "`", collapse = ", "),
          ". It takes (by likelihood): ",
          paste(sort(known), collapse = ", "), ".")
    mode_bad <- setdiff(names(spec$mode$extra),
                        c(known, .rmf_optimiser_arg_names(), "verbose"))
    if (length(mode_bad))
      bad_arg("mode$extra has argument(s) neither the log-posterior ",
          "constructor nor the optimiser takes: ",
          paste0("`", mode_bad, "`", collapse = ", "), ".")
  }

  ## whittle freq_band: 0 <= lo < hi <= pi (make_log_posterior_whittle's rule)
  if (identical(eff_type, "whittle")) {
    fb <- lik$freq_band
    if (!.whittle_freq_band_ok(fb))
      bad("likelihood$freq_band must be c(lo, hi) with 0 <= lo < hi <= pi ",
          "(radians); got c(", paste(format(fb), collapse = ", "), ").")
  }

  ## me_variance: a scalar, or one per observable (the constructor's rule)
  if (!obc_on && length(lik$me_variance) > 1L && !.spec_tunes_possible(spec))
    .spec_reclass(.kf_me_variance(
      lik$me_variance, spec$obs_vars,
      sprintf("likelihood$me_variance (likelihood = \"%s\")", lik$type),
      allow_vector = lik$type %in% c("gaussian", "pskf") &&
        is.null(lik$ms_spec) && is.null(lik$ms_struct_spec)),
      class = "dynhr_error_me_variance")

  ## pruned_order = 3 needs derivatives of order 3 (solve_perturbation's rule:
  ## orders 2 and 3 need a model compiled to max_order >= 2)
  if (is.null(spec$mode$result) && eff_type %in% c("pruned", "tpf") &&
      identical(lik$pruned_order, 3L) && spec$model$max_order < 2L)
    bad("likelihood$pruned_order = 3 needs a model compiled to max_order >= ",
        "2 (this spec's model$max_order is ", spec$model$max_order, "): the ",
        "order-3 solve would fail at every draw and the mode stage would ",
        "'converge' at -Inf. Build the spec with max_order = 2 (or a ",
        "compiled model of that order).")

  ## mode$theta_init: exactly the estimated parameters, inside their bounds
  ## (run_mode_finding's rule)
  if (!is.null(spec$mode$theta_init) && mode_runs) {
    pr <- .spec_prior_spec(spec)
    if (!is.null(pr))
      .spec_reclass(.rmf_check_theta_init(spec$mode$theta_init, pr))
  }

  ## likelihood$dates: used by the diagnostics only, one per sample row
  if (!is.null(lik$dates)) {
    if (!is.null(n_rows) && length(lik$dates) != n_rows)
      bad("likelihood$dates has ", length(lik$dates), " entries but the ",
          "likelihood sample has ", n_rows, " data rows: give one date per ",
          "row of the sample (first_obs / nobs window).")
    if (!isTRUE(out$diagnostics))
      .dynhr_warn("validate_spec: likelihood$dates labels the diagnostics' ",
                  "time axis only and outputs$diagnostics is FALSE, so it has ",
                  "no effect.", class = "dynhr_warning_spec_output_ignored")
  }

  ## outputs$ramsey: the Ramsey workflow needs a planner objective in the
  ## model, a square system, derivatives of its order and a discount factor
  if (isTRUE(out$ramsey)) {
    mod <- spec$model$mod
    if (!nzchar(trimws(.spec_planner_objective(mod))))
      bad("outputs$ramsey = TRUE, but the model has no planner objective ",
          "(add `planner_objective(...);` to the .mod file): the sampler ",
          "would run to the end and the Ramsey step then fail with \"No ",
          "planner objective provided\".")
    if (length(mod$equations) < length(mod$var_names))
      bad("outputs$ramsey = TRUE, but the model is a free-instrument system (",
          length(mod$equations), " equations, ", length(mod$var_names),
          " endogenous variables): ramsey_policy() needs a square system.")
    if (out$ramsey_order >= 2L && spec$model$max_order < 2L)
      bad("outputs$ramsey_order = ", out$ramsey_order, " needs a model ",
          "compiled to max_order >= 2 (this spec's model$max_order is ",
          spec$model$max_order, ").")
    if (is.null(out$ramsey_discount) &&
        is.null(.get_discount(mod$param_values)))
      .dynhr_warn("validate_spec: outputs$ramsey = TRUE, but the model has no ",
                  "`beta` / `betta` parameter and outputs$ramsey_discount is ",
                  "NULL: the Ramsey step will be skipped after estimation. ",
                  "Set outputs$ramsey_discount.",
                  class = "dynhr_warning_spec_ramsey_discount")
  }

  ## outputs$save: dir must be (or be creatable as) a directory; prefix is a
  ## file-name prefix, not a path
  if (isTRUE(out$save)) {
    anc <- .spec_existing_ancestor(out$dir)
    if (!dir.exists(anc))
      bad("outputs$save = TRUE, but outputs$dir \"", out$dir, "\" is, or sits ",
          "under, an existing file (\"", anc, "\"), not a directory.")
    if (grepl("[/\\\\]", out$prefix))
      bad("outputs$prefix \"", out$prefix, "\" contains a path separator; ",
          "it is a file-name prefix (put the directory in outputs$dir).")
  }

  ## ---- result form (run_estimation()) ------------------------------------
  form <- out$form
  has_mr <- !is.null(spec$mode$result)
  if (identical(form, "mode") && (length(samplers) || has_mr))
    bad("outputs$form = \"mode\" runs the mode stage only: the spec must have ",
        "no sampler (sampler = NULL) and no precomputed mode$result.")
  if (identical(form, "posterior") && !length(samplers))
    bad("outputs$form = \"posterior\" needs a sampler.")
  if (identical(form, "full") && length(samplers) > 1L)
    bad("outputs$form = \"full\" runs one sampler; a sampler sequence needs ",
        "form \"posterior\".")
  if (identical(form, "full") && has_mr)
    bad("outputs$form = \"full\" runs the mode stage itself; a precomputed ",
        "mode$result needs form \"posterior\".")

  ## ---- outputs that need OBC ---------------------------------------------
  if (!obc_on && (isTRUE(out$smoother) || isTRUE(out$obc_ppf_reweight)))
    .dynhr_warn("validate_spec: outputs$",
                if (isTRUE(out$smoother)) "smoother" else "obc_ppf_reweight",
                " has an effect only on an OBC model; this model has none.",
                class = "dynhr_warning_spec_output_ignored")
  invisible(spec)
}

# ---------------------------------------------------------------------------
# Printing, summaries, diffs
# ---------------------------------------------------------------------------

## Fields of a sub-spec that differ from the schema default.
.spec_nondefault <- function(x) {
  cp  <- .spec_component_of(x)
  sch <- .spec_component_schema(cp, x$method)
  x <- unclass(x)
  keep <- vapply(names(x), function(nm)
    !identical(x[[nm]], sch[[nm]]$default), logical(1))
  x[keep & names(x) != "method"]
}

.spec_fmt_fields <- function(x) {
  if (!length(x)) return("(defaults)")
  paste(sprintf("%s=%s", names(x), vapply(x, .spec_fmt_value, character(1),
                                          width = 40L)), collapse = ", ")
}

.spec_short_hash <- function(h)
  if (length(h) != 1L || is.na(h)) "none" else sub("^(\\w+):(.{12}).*$", "\\1:\\2", h)

#' @rdname dynhr_estimation_spec
#' @param x A \code{dynhr_estimation_spec}.
#' @export
format.dynhr_estimation_spec <- function(x, ...) {
  mp <- x$model
  model_txt <- paste0(if (!is.null(mp$path)) basename(mp$path) else "<in-memory dynhr_mod>",
                      sprintf(" (%d endo, %d exo, %d params)",
                              length(mp$mod$var_names), length(mp$mod$varexo_names),
                              length(mp$mod$param_names)),
                      if (!is.null(mp$prior_spec)) " + custom prior_spec")
  data_txt <- if (!is.null(x$data$value))
    sprintf("%d x %d matrix", nrow(x$data$value), ncol(x$data$value))
  else if (!is.null(x$data$path)) basename(x$data$path) else "none"
  sl <- .spec_sampler_list(x$sampler)
  samp_txt <- if (!length(sl)) "none (mode finding only)" else
    paste(vapply(sl, function(s) paste0(s$method, ": ",
                                        .spec_fmt_fields(.spec_nondefault(s))),
                 character(1)), collapse = "\n               then ")
  opt_set <- x$options[vapply(names(x$options), function(nm) {
    reg <- .dynhr_option_registry[[nm]]
    is.null(reg) || !identical(reg$default, x$options[[nm]])
  }, logical(1))]
  c(sprintf("<dynhr_estimation_spec> version %d  [%s]", x$spec_version,
            .spec_short_hash(x$hashes$spec)),
    paste0("  model      : ", model_txt),
    paste0("  data       : ", data_txt, "  (obs: ", paste(x$obs_vars, collapse = ", "), ")"),
    paste0("  likelihood : ", x$likelihood$type, "  ",
           .spec_fmt_fields(.spec_nondefault(x$likelihood)[
             setdiff(names(.spec_nondefault(x$likelihood)), "type")])),
    paste0("  mode       : ", if (!is.null(x$mode$result)) "precomputed result; " else "",
           .spec_fmt_fields(.spec_nondefault(x$mode)[
             setdiff(names(.spec_nondefault(x$mode)), "result")])),
    paste0("  sampler    : ", samp_txt),
    paste0("  compute    : ", .spec_fmt_fields(.spec_nondefault(x$compute))),
    paste0("  outputs    : ", .spec_fmt_fields(.spec_nondefault(x$outputs))),
    paste0("  options    : ", .spec_fmt_fields(opt_set)))
}

#' @rdname dynhr_estimation_spec
#' @export
print.dynhr_estimation_spec <- print.dynhr_subspec

## Flatten a spec to path -> value (field level; model and data by hash).
.spec_flatten <- function(spec) {
  out <- list(
    "model"      = spec$hashes$model,
    "model$path" = spec$model$path,
    "data"       = spec$hashes$data,
    "data$path"  = spec$data$path,
    "obs_vars"   = spec$obs_vars)
  for (cp in c("likelihood", "mode", "compute", "outputs")) {
    x <- unclass(spec[[cp]])
    for (nm in names(x)) out[paste0(cp, "$", nm)] <- list(x[[nm]])
  }
  sl <- .spec_sampler_list(spec$sampler)
  if (!length(sl)) out["sampler"] <- list(NULL)
  for (i in seq_along(sl)) {
    pre <- if (length(sl) == 1L) "sampler$" else sprintf("sampler[[%d]]$", i)
    x <- unclass(sl[[i]])
    for (nm in names(x)) out[paste0(pre, nm)] <- list(x[[nm]])
  }
  for (nm in names(spec$options)) out[paste0("options$", nm)] <- list(spec$options[[nm]])
  out
}

#' @rdname dynhr_estimation_spec
#' @details \code{summary(spec)} returns every field as a table (class
#'   \code{dynhr_spec_summary}): \code{path}, \code{value} and whether the
#'   value is the schema default.
#' @export
summary.dynhr_estimation_spec <- function(object, ...) {
  fl <- .spec_flatten(object)
  is_def <- vapply(names(fl), function(p) {
    if (!grepl("^(likelihood|mode|sampler|compute|outputs)(\\[\\[\\d+\\]\\])?\\$", p))
      return(NA)
    cp  <- sub("(\\[\\[\\d+\\]\\])?\\$.*$", "", p)
    fld <- sub("^.*\\$", "", p)
    meth <- if (identical(cp, "sampler")) fl[[sub("\\$[^$]*$", "$method", p)]]
    def <- .spec_component_schema(cp, meth)[[fld]]
    identical(fl[[p]], def$default)
  }, logical(1))
  structure(data.frame(path = names(fl),
                       value = vapply(fl, .spec_fmt_value, character(1)),
                       default = unname(is_def), stringsAsFactors = FALSE,
                       row.names = NULL),
            class = c("dynhr_spec_summary", "data.frame"))
}

#' @export
print.dynhr_spec_summary <- function(x, ...) {
  w <- max(nchar(x$path))
  cat(sprintf("%-*s  %s%s", w, x$path, x$value,
              ifelse(x$default %in% FALSE, "   *", "")), sep = "\n")
  cat("(* = differs from the schema default)\n")
  invisible(x)
}

#' @rdname dynhr_estimation_spec
#' @export
as.list.dynhr_estimation_spec <- function(x, ...) {
  x <- unclass(x)
  for (cp in c("likelihood", "mode", "compute", "outputs")) x[[cp]] <- unclass(x[[cp]])
  if (!is.null(x$sampler))
    x$sampler <- if (inherits(x$sampler, "dynhr_sampler_sequence"))
      lapply(unclass(x$sampler), unclass) else unclass(x$sampler)
  x
}

#' @export
print.dynhr_sampler_sequence <- function(x, ...) {
  cat(sprintf("<sampler sequence: %s>\n",
              paste(vapply(x, `[[`, "", "method"), collapse = " -> ")))
  for (s in x) print(s)
  invisible(x)
}

#' Compare two estimation specs field by field
#'
#' @param a,b Two \code{dynhr_estimation_spec} objects.
#' @return A data.frame (class \code{dynhr_spec_diff}) with one row per field
#'   that differs: \code{path}, and the formatted values \code{a} and
#'   \code{b}. The model and the data are compared by content hash. Zero rows
#'   when the specs agree.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @export
diff_specs <- function(a, b) {
  if (!inherits(a, "dynhr_estimation_spec") || !inherits(b, "dynhr_estimation_spec"))
    .dynhr_abort("diff_specs: `a` and `b` must be dynhr_estimation_spec objects.",
                 class = "dynhr_error_bad_argument")
  fa <- .spec_flatten(a)
  fb <- .spec_flatten(b)
  paths <- unique(c(names(fa), names(fb)))
  miss <- function(f, p) if (p %in% names(f)) .spec_fmt_value(f[[p]]) else "<absent>"
  diffp <- paths[!vapply(paths, function(p)
    (p %in% names(fa)) == (p %in% names(fb)) && identical(fa[[p]], fb[[p]]),
    logical(1))]
  structure(data.frame(path = diffp,
                       a = vapply(diffp, miss, character(1), f = fa),
                       b = vapply(diffp, miss, character(1), f = fb),
                       stringsAsFactors = FALSE, row.names = NULL),
            class = c("dynhr_spec_diff", "data.frame"))
}

#' @export
print.dynhr_spec_diff <- function(x, ...) {
  if (!nrow(x)) {
    cat("<spec diff> no differences\n")
  } else {
    cat(sprintf("<spec diff> %d field(s) differ\n", nrow(x)))
    w <- max(nchar(x$path))
    cat(sprintf("  %-*s : %s -> %s", w, x$path, x$a, x$b), sep = "\n")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# as_estimation_spec(): from flat arguments, records, results and .mod files
# ---------------------------------------------------------------------------

## Where each argument of the three entry points lands in the spec. `...`
## goes to one `extra` list. (Used by as_estimation_spec() and its tests.)
.spec_flat_map <- list(
  run_full_estimation = c(
    mod_file = "model$path", data = "data", obs_vars = "obs_vars",
    output_dir = "outputs$dir", output_prefix = "outputs$prefix",
    sampler = "sampler$method", n_draws = "sampler$n_draws",
    n_warmup = "sampler$n_warmup", n_chains = "sampler$n_chains",
    n_particles = "sampler$n_particles", n_walkers = "sampler$n_walkers",
    parallel = "compute$parallel", parallel_backend = "compute$backend",
    n_cores = "compute$n_cores", analytic_grad = "sampler$analytic_grad",
    metric = "sampler$metric", n_mode_iter = "mode$n_iter",
    mode_method = "mode$method", mode_n_starts = "mode$n_starts",
    me_variance = "likelihood$me_variance", likelihood = "likelihood$type",
    lik_init = "likelihood$lik_init", filter_method = "likelihood$filter_method",
    singular_obs = "likelihood$singular_obs",
    freq_band = "likelihood$freq_band",
    system_priors = "likelihood$system_priors", seed = "compute$seed",
    run_diag = "outputs$diagnostics", verbose = "compute$verbose",
    model = "model$mod", compiled = "model$compiled",
    data_col_map = "likelihood$data_col_map", dates = "likelihood$dates",
    obc = "likelihood$obc", obc_max_inner = "likelihood$obc_max_inner",
    compute_smoother = "outputs$smoother",
    obc_ppf_reweight = "outputs$obc_ppf_reweight", run_ramsey = "outputs$ramsey",
    ramsey_order = "outputs$ramsey_order",
    ramsey_n_periods = "outputs$ramsey_n_periods",
    ramsey_burn_in = "outputs$ramsey_burn_in",
    ramsey_discount = "outputs$ramsey_discount",
    filter_tunes = "likelihood$filter_tunes",
    heteroskedastic_shocks = "likelihood$heteroskedastic_shocks",
    stochastic_volatility = "likelihood$stochastic_volatility",
    tpf_options = "likelihood$tpf_options", plan = "likelihood$plan",
    checkpoint_dir = "compute$checkpoint_dir", resume = "compute$resume",
    on_mismatch = "compute$on_mismatch",
    "..." = "sampler$extra"),
  run_mode_finding = c(
    solved = "model$mod", data = "data", obs_vars = "obs_vars",
    n_iter = "mode$n_iter", method = "mode$method",
    me_variance = "likelihood$me_variance", likelihood = "likelihood$type",
    pruned_order = "likelihood$pruned_order",
    data_col_map = "likelihood$data_col_map", mode_options = "mode$options",
    posterior_options = "likelihood$extra", parallel = "compute$parallel",
    n_cores = "compute$n_cores", n_starts = "mode$n_starts",
    proposal_cov_method = "mode$proposal_cov",
    transform_params = "mode$transform_params",
    filter_tunes = "likelihood$filter_tunes",
    heteroskedastic_shocks = "likelihood$heteroskedastic_shocks",
    stochastic_volatility = "likelihood$stochastic_volatility",
    plan = "likelihood$plan", use_exact_hessian = "mode$exact_hessian",
    verbose = "compute$verbose", theta_init = "mode$theta_init",
    "..." = "mode$extra"),
  run_posterior_estimation = c(
    mode_result = "mode$result", methods = "sampler$method",
    n_warmup = "sampler$n_warmup", n_draws = "sampler$n_draws",
    n_chains = "sampler$n_chains", n_particles = "sampler$n_particles",
    n_walkers = "sampler$n_walkers", parallel = "compute$parallel",
    parallel_backend = "compute$backend", n_cores = "compute$n_cores",
    analytic_grad = "sampler$analytic_grad", Sigma_prop = "sampler$Sigma_prop",
    run_stoch_simul = "outputs$stoch_simul",
    skip_mode_finding_check = "mode$skip_check",
    nuts_timeout_seconds = "sampler$timeout",
    transform_params = "sampler$transform_params",
    rwmh_adapt_cov = "sampler$adapt_cov", rwmh_n_blocks = "sampler$n_blocks",
    metric = "sampler$metric", monge_alpha = "sampler$monge_alpha",
    checkpoint_dir = "compute$checkpoint_dir", resume = "compute$resume",
    on_mismatch = "compute$on_mismatch",
    verbose = "compute$verbose", seed = "compute$seed",
    "..." = "sampler$extra")
)

## Keys of run_mode_finding()'s posterior_options that are likelihood fields.
.spec_rmf_posterior_keys <- c(lik_init = "lik_init", filter_method = "filter_method",
  singular_obs = "singular_obs",
  freq_band = "freq_band",
  system_priors = "system_priors", infeasible_penalty = "infeasible_penalty",
  obc_filter = "obc_filter", max_inner = "obc_max_inner",
  tpf_options = "tpf_options", obc_options = "obc_options")

## Formals of `entry` that the entry point resolves with match.arg().
.spec_matcharg <- list(run_full_estimation = c("sampler", "metric", "likelihood"),
                       run_mode_finding = "likelihood",
                       run_posterior_estimation = "metric")

## Resolve a flat argument list against `entry`'s formals: every formal gets
## its supplied value or its current default; match.arg formals collapse to
## one choice; names that are not formals are the `...` arguments.
.spec_flat_resolve <- function(x, entry) {
  f    <- get(entry, envir = asNamespace("dynhr"), inherits = FALSE)
  fmls <- setdiff(names(formals(f)), "...")
  nms  <- names(x)
  if (length(x) && (is.null(nms) || any(!nzchar(nms))))
    .dynhr_abort("as_estimation_spec: the argument list for ", entry,
                 "() must be fully named.", class = "dynhr_error_bad_argument")
  args <- list()
  defaults <- list()
  for (nm in fmls) {
    d <- .dynhr_rr_current_default(f, nm)
    if (isTRUE(d$ok)) defaults[nm] <- list(d$value)
    args[nm] <- list(if (nm %in% nms) x[[nm]] else if (isTRUE(d$ok)) d$value)
    if (nm %in% .spec_matcharg[[entry]]) {
      choices <- eval(formals(f)[[nm]], envir = asNamespace("dynhr"))
      v <- args[[nm]]
      args[nm] <- list(if (identical(v, choices)) choices[[1L]] else
        match.arg(v, choices))
      defaults[nm] <- list(choices[[1L]])
    }
  }
  list(args = args, dots = x[setdiff(nms, fmls)], defaults = defaults)
}

## Say which supplied (non-default) arguments a converter could not place.
.spec_note_dropped <- function(entry, a, defaults, dropped) {
  dropped <- unique(dropped)
  dropped <- dropped[!vapply(dropped, function(nm)
    identical(a[[nm]], defaults[[nm]]), logical(1))]
  if (length(dropped))
    .dynhr_inform("as_estimation_spec: ", entry, "() argument(s) ",
                  paste(dropped, collapse = ", "), " do not apply to the ",
                  "chosen sampler(s) and are not part of the spec (",
                  entry, "() ignored them too).",
                  class = "dynhr_message_spec_ignored_args")
  invisible(NULL)
}

## A sampler spec from candidate fields: keeps the fields `method` takes.
## Returns list(spec, dropped = names of non-NULL candidates it does not take).
.spec_sampler_from <- function(method, cand) {
  takes <- .spec_sampler_table[[method]]
  keep  <- intersect(names(cand), takes)
  list(spec = .spec_build("sampler", cand[keep], method = method),
       dropped = setdiff(names(cand), takes))
}

## Metric of a run_posterior_estimation() / run_full_estimation() call as the
## given method used it (hmc / mala silently used "diagonal" for the metrics
## they do not implement).
.spec_method_metric <- function(method, metric) {
  ch <- .spec_sampler_metrics[[method]]
  if (is.null(ch)) return(NULL)
  if (metric %in% ch) metric else "diagonal"
}

.spec_parts_rfe <- function(a, dots) {
  mp <- if (!is.null(a$model)) {
    .spec_model_part(a$model, compiled = a$compiled)
  } else {
    if (is.null(a$mod_file))
      .dynhr_abort("as_estimation_spec: provide `mod_file` or `model`.",
                   class = "dynhr_error_spec_invalid")
    .spec_model_part(a$mod_file, compiled = a$compiled)
  }
  ## likelihood = "pkf" / "ppf" / "copf" names the OBC filter (E5 C3: the
  ## widened choices); an explicit obc_filter in `...` wins.
  obc_filter <- dots$obc_filter %||%
    (if (a$likelihood %in% .spec_obc_likelihoods) a$likelihood else "pkf")
  dots$obc_filter <- NULL
  student_df <- dots$student_df
  dots$student_df <- NULL
  lik_extra <- if (identical(a$likelihood, "sv_rbpf"))
    list(n_particles = a$n_particles) else list()
  lik <- likelihood_spec(a$likelihood,
    me_variance = a$me_variance, lik_init = a$lik_init,
    filter_method = a$filter_method,
    singular_obs = a$singular_obs,
    freq_band = a$freq_band, system_priors = a$system_priors,
    tpf_options = a$tpf_options, obc = a$obc, obc_max_inner = a$obc_max_inner,
    obc_filter = obc_filter, filter_tunes = a$filter_tunes,
    heteroskedastic_shocks = a$heteroskedastic_shocks,
    stochastic_volatility = a$stochastic_volatility, plan = a$plan,
    data_col_map = a$data_col_map, dates = a$dates, student_df = student_df,
    extra = lik_extra)
  mode <- mode_spec(method = a$mode_method, n_iter = a$n_mode_iter,
                    n_starts = a$mode_n_starts)
  samp <- NULL
  dropped <- character(0)
  if (a$n_draws > 0) {
    cand <- list(n_draws = a$n_draws, n_warmup = a$n_warmup,
                 n_chains = a$n_chains, n_particles = a$n_particles,
                 n_walkers = a$n_walkers, analytic_grad = a$analytic_grad,
                 metric = .spec_method_metric(a$sampler, a$metric) %||% a$metric,
                 extra = dots)
    r <- .spec_sampler_from(a$sampler, cand)
    samp <- r$spec
    dropped <- r$dropped
    ## sv_rbpf: n_particles is the RB-PF particle count (likelihood$extra)
    if (identical(a$likelihood, "sv_rbpf")) dropped <- setdiff(dropped, "n_particles")
  }
  list(parts = list(
    model = mp, data = .spec_data_part(a$data),
    obs_vars = .spec_resolve_obs_vars(a$obs_vars, mp$mod),
    likelihood = lik, mode = mode, sampler = samp,
    compute = compute_spec(parallel = a$parallel, backend = a$parallel_backend,
                           n_cores = a$n_cores, seed = a$seed,
                           verbose = a$verbose,
                           checkpoint_dir = a$checkpoint_dir,
                           resume = a$resume %||% FALSE,
                           on_mismatch = a$on_mismatch %||% "refuse"),
    outputs = outputs_spec(form = "full",
      dir = a$output_dir, prefix = a$output_prefix,
      save = TRUE, diagnostics = a$run_diag, smoother = a$compute_smoother,
      obc_ppf_reweight = a$obc_ppf_reweight, ramsey = a$run_ramsey,
      ramsey_order = a$ramsey_order, ramsey_n_periods = a$ramsey_n_periods,
      ramsey_burn_in = a$ramsey_burn_in, ramsey_discount = a$ramsey_discount),
    options = .spec_options_snapshot()),
    dropped = dropped)
}

.spec_parts_rmf <- function(a, dots, max_order = 1L) {
  if (is.null(a$solved))
    .dynhr_abort("as_estimation_spec: run_mode_finding() needs `solved`.",
                 class = "dynhr_error_spec_invalid")
  mp <- .spec_model_part(a$solved, max_order = max_order)
  po <- a$posterior_options %||% list()
  lik_f <- list()
  for (k in intersect(names(po), names(.spec_rmf_posterior_keys)))
    lik_f[.spec_rmf_posterior_keys[[k]]] <- list(po[[k]])
  lik_extra <- po[setdiff(names(po), names(.spec_rmf_posterior_keys))]
  if (!is.null(dots$student_df)) {
    lik_f$student_df <- dots$student_df
    dots$student_df <- NULL
  }
  lik <- do.call(likelihood_spec, c(list(a$likelihood,
    me_variance = a$me_variance, pruned_order = a$pruned_order,
    data_col_map = a$data_col_map, filter_tunes = a$filter_tunes,
    heteroskedastic_shocks = a$heteroskedastic_shocks,
    stochastic_volatility = a$stochastic_volatility, plan = a$plan,
    extra = lik_extra), lik_f))
  mode <- mode_spec(method = a$method, n_iter = a$n_iter, n_starts = a$n_starts,
                    transform_params = a$transform_params,
                    proposal_cov = a$proposal_cov_method,
                    exact_hessian = a$use_exact_hessian,
                    theta_init = a$theta_init, options = a$mode_options,
                    extra = dots)
  list(parts = list(
    model = mp, data = .spec_data_part(a$data),
    obs_vars = .spec_resolve_obs_vars(a$obs_vars, mp$mod),
    likelihood = lik, mode = mode, sampler = NULL,
    compute = compute_spec(parallel = a$parallel, n_cores = a$n_cores,
                           seed = NULL, verbose = a$verbose),
    outputs = outputs_spec(form = "mode"),
    options = .spec_options_snapshot()),
    dropped = character(0))
}

## Parts for the mode side of a mode result that carries no run record.
.spec_parts_from_mode_result <- function(mr) {
  ctx <- if (inherits(mr$ctx, "dynhr_estimation_context")) mr$ctx else
    ctx_from_mode_result(mr)
  obc <- !is.null(mr$obc_specs)
  lik <- likelihood_spec(
    if (obc) "gaussian" else ctx$likelihood,
    me_variance = ctx$me_variance, lik_init = ctx$lik_init,
    filter_method = ctx$filter_method %||% "auto",
    singular_obs = ctx$singular_obs %||% "reject",
    freq_band = ctx$freq_band, system_priors = ctx$system_priors,
    tpf_options = ctx$tpf_options, obc_options = ctx$obc_options,
    obc_filter = if (obc) ctx$likelihood else "pkf",
    gradient_policy = ctx$gradient_policy, student_df = ctx$student_df,
    pruned_order = ctx$pruned_order, ms_spec = ctx$ms_spec,
    ms_struct_spec = ctx$ms_struct_spec, ms_collapse = ctx$ms_collapse)
  mp <- .spec_model_part(mr$solved)
  list(model = mp, data = .spec_data_part(mr$data),
       obs_vars = .spec_resolve_obs_vars(mr$obs_vars, mp$mod),
       likelihood = lik, mode = mode_spec())
}

## Sampler / compute / outputs parts of a run_posterior_estimation() call on
## top of `base` (the mode side's parts).
.spec_parts_rpe <- function(a, dots, base, mode_result = NULL) {
  methods <- tolower(as.character(a$methods))
  n <- length(methods)
  rec <- function(x) if (is.null(x)) NULL else rep_len(x, n)
  nw <- rec(a$n_warmup); nd <- rec(a$n_draws); nc <- rec(a$n_chains)
  np <- rec(a$n_particles)
  eff_type <- if (.spec_has_mcp(base$model$mod) && !isFALSE(base$likelihood$obc))
    base$likelihood$obc_filter else base$likelihood$type
  specs <- list()
  dropped <- character(0)
  used <- character(0)
  for (i in seq_len(n)) {
    m <- methods[[i]]
    if (!m %in% .spec_samplers)
      .dynhr_abort("as_estimation_spec: unknown method \"", toupper(m),
                   "\" for run_posterior_estimation().",
                   class = "dynhr_error_spec_invalid")
    if (identical(m, "rwmh") && eff_type %in% .spec_pmmh_likelihoods) m <- "pmmh"
    cand <- list(n_warmup = nw[[i]], n_draws = nd[[i]], n_chains = nc[[i]],
                 n_particles = np[[i]], n_walkers = a$n_walkers,
                 analytic_grad = a$analytic_grad,
                 metric = .spec_method_metric(m, a$metric) %||% a$metric,
                 transform_params = a$transform_params,
                 adapt_cov = a$rwmh_adapt_cov, n_blocks = a$rwmh_n_blocks,
                 monge_alpha = a$monge_alpha, timeout = a$nuts_timeout_seconds,
                 Sigma_prop = a$Sigma_prop, extra = dots)
    r <- .spec_sampler_from(m, cand)
    specs[[i]] <- r$spec
    dropped <- c(dropped, r$dropped)
    used <- c(used, intersect(names(cand), .spec_sampler_table[[m]]))
  }
  argname <- c(n_walkers = "n_walkers", analytic_grad = "analytic_grad",
               metric = "metric", adapt_cov = "rwmh_adapt_cov",
               n_blocks = "rwmh_n_blocks", monge_alpha = "monge_alpha",
               timeout = "nuts_timeout_seconds", Sigma_prop = "Sigma_prop",
               n_particles = "n_particles", n_chains = "n_chains",
               n_draws = "n_draws", n_warmup = "n_warmup",
               transform_params = "transform_params")
  never <- setdiff(dropped, c(used, "extra"))
  mode <- base$mode
  mode <- .spec_build("mode", list(result = mode_result,
                                   skip_check = a$skip_mode_finding_check),
                      base = mode)
  list(parts = c(base[c("model", "data", "obs_vars", "likelihood")], list(
    mode = mode,
    sampler = if (n == 1L) specs[[1L]] else
      structure(specs, class = "dynhr_sampler_sequence"),
    compute = compute_spec(parallel = a$parallel, backend = a$parallel_backend,
                           n_cores = a$n_cores, seed = a$seed,
                           verbose = a$verbose,
                           checkpoint_dir = a$checkpoint_dir,
                           resume = a$resume,
                           on_mismatch = a$on_mismatch %||% "refuse"),
    outputs = outputs_spec(form = "posterior", stoch_simul = a$run_stoch_simul),
    options = .spec_options_snapshot())),
    dropped = unname(argname[never]))
}

## Parts from flat arguments of `entry` (under the current option source).
.spec_parts_flat <- function(x, entry) {
  r <- .spec_flat_resolve(x, entry)
  a <- r$args
  out <- switch(entry,
    run_full_estimation = .spec_parts_rfe(a, r$dots),
    run_mode_finding    = .spec_parts_rmf(a, r$dots),
    run_posterior_estimation = {
      mr <- a$mode_result
      if (is.null(mr))
        .dynhr_abort("as_estimation_spec: run_posterior_estimation() needs ",
                     "`mode_result`.", class = "dynhr_error_spec_invalid")
      mrr <- mr$run_record
      base <- if (inherits(mrr, "dynhr_run_record") &&
                  identical(mrr$fn, "run_mode_finding"))
        .spec_parts_record(mrr)$parts else .spec_parts_from_mode_result(mr)
      .spec_parts_rpe(a, r$dots, base, mode_result = mr)
    })
  .spec_note_dropped(entry, a, r$defaults, out$dropped)
  ## `...` / posterior_options keys that name a typed field go to that field
  out$parts <- .spec_route_extras(out$parts)
  out
}

## Parts from a run record, each stage under the option snapshot it ran with.
.spec_parts_record <- function(rec) {
  if (!inherits(rec, "dynhr_run_record"))
    .dynhr_abort("as_estimation_spec: not a dynhr_run_record.",
                 class = "dynhr_error_no_run_record")
  sv <- rec$schema_version
  if (!is.numeric(sv) || length(sv) != 1L || sv > .dynhr_run_record_schema)
    .dynhr_abort("as_estimation_spec: run record schema version ", format(sv),
                 " is not supported.", class = "dynhr_error_run_record_schema")
  if (!rec$fn %in% .dynhr_rr_entry_points)
    .dynhr_abort("as_estimation_spec: the record names an unknown entry point.",
                 class = "dynhr_error_run_record_schema")
  if (identical(rec$fn, "run_estimation"))
    return(list(parts = .spec_parts(rec$spec %||% rec$args$spec), dropped = character(0)))
  args <- rec$args
  .spec_with_snapshot(rec$options, function() {
    out <- switch(rec$fn,
      run_mode_finding = {
        r <- .spec_flat_resolve(args, "run_mode_finding")
        .spec_parts_rmf(r$args, r$dots,
                        max_order = rec$rebuild$solved_compiled$max_order %||% 1L)
      },
      run_full_estimation = {
        if (!is.null(rec$inputs$model)) args["model"] <- list(rec$inputs$model)
        if (!is.null(rec$inputs$data))  args["data"]  <- list(rec$inputs$data)
        r <- .spec_flat_resolve(args, "run_full_estimation")
        .spec_parts_rfe(r$args, r$dots)
      },
      run_posterior_estimation = {
        ref <- args$mode_result
        if (!inherits(ref, "dynhr_mode_ref") ||
            !inherits(ref$run_record, "dynhr_run_record"))
          .dynhr_abort("as_estimation_spec: this run_posterior_estimation() ",
                       "record has no replayable mode run record.",
                       class = "dynhr_error_rerun_needs_mode_result")
        base <- .spec_parts_record(ref$run_record)$parts
        args["mode_result"] <- list(NULL)
        r <- .spec_flat_resolve(args, "run_posterior_estimation")
        .spec_parts_rpe(r$args, r$dots, base, mode_result = NULL)
      })
    out$parts <- .spec_route_extras(out$parts)
    out
  })
}

#' Build an estimation spec from arguments, a run record, a result or a .mod
#'
#' \code{as_estimation_spec()} converts the existing ways of describing an
#' estimation into a \code{\link{dynhr_estimation_spec}}:
#' \describe{
#'   \item{A flat argument list}{(a named list, with \code{entry} naming the
#'     function) of \code{\link{run_full_estimation}},
#'     \code{\link{run_posterior_estimation}} or
#'     \code{\link{run_mode_finding}}: exactly the arguments that function
#'     takes, the unnamed-in-the-signature ones (\code{...}) included.
#'     Arguments not supplied take the function's defaults; \code{NULL}
#'     arguments that the function resolves from an option read the option.
#'     A \code{run_posterior_estimation()} list becomes a spec whose mode
#'     stage is the supplied \code{mode_result} (and, when that result has a
#'     run record, its settings); several \code{methods} become a sampler
#'     sequence; \code{"RWMH"} on a particle likelihood becomes
#'     \code{"pmmh"}. Arguments that do not apply to the chosen sampler (and
#'     that the function ignores) are left out, with a message.}
#'   \item{A \code{dynhr_run_record}}{(schema 1), or an estimation result
#'     carrying one: the recorded arguments, with each stage's option-sourced
#'     fields taken from that stage's recorded option snapshot.}
#'   \item{A \code{.mod} path}{Its \code{estimation(...)} options become
#'     defaults, mapped conservatively: \code{mh_replic} and \code{mh_drop}
#'     (draws and warmup), \code{mh_nblocks} (chains),
#'     \code{posterior_sampling_method} (random-walk MH, \code{hssmc},
#'     \code{dsmh}), \code{mode_compute} (4 = \code{"newrat"}, 9 =
#'     \code{"cmaes"}), \code{first_obs} and \code{nobs} (the likelihood
#'     sample), \code{lik_init} (1 = stationary, 3 = diffuse) and
#'     \code{diffuse_filter}, \code{prefilter = 0} and \code{order = 1} (no-ops),
#'     and a CSV \code{datafile}. Display-only options are ignored; every other
#'     option is listed in a message. The arguments in \code{...} (as for
#'     \code{\link{dynhr_estimation_spec}}) override the mapped values: a
#'     named list is merged over them, a complete sub-spec replaces them.}
#' }
#'
#' @param x A named argument list, a \code{dynhr_run_record}, an estimation
#'   result, a \code{.mod} path, a \code{dynhr_estimation_spec} (returned
#'   validated) or a \code{\link{dynhr_model}} (the spec it holds, validated).
#' @param ... Method arguments: \code{entry} for a list; for a \code{.mod}
#'   path, \code{data}, \code{obs_vars}, \code{likelihood}, \code{mode},
#'   \code{sampler}, \code{compute}, \code{outputs} and \code{options}.
#' @return A \code{dynhr_estimation_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}, \code{\link{dynhr_rerun}}
#' @examples
#' mod <- system.file("extdata/models/nk_demo.mod", package = "dynhr")
#' dat <- system.file("extdata/models/nk_demo_data.csv", package = "dynhr")
#' as_estimation_spec(list(mod_file = mod, data = dat,
#'                         obs_vars = c("ygr", "infl", "intr"),
#'                         sampler = "nuts", n_draws = 1000L),
#'                    entry = "run_full_estimation")
#' @export
as_estimation_spec <- function(x, ...) UseMethod("as_estimation_spec")

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.default <- function(x, ...) {
  .dynhr_abort("as_estimation_spec: cannot build a spec from an object of ",
               "class ", paste(class(x), collapse = "/"), ".",
               class = "dynhr_error_bad_argument")
}

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.dynhr_estimation_spec <- function(x, ...) validate_spec(x)

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.dynhr_model <- function(x, ...) {
  if (length(list(...)))
    .dynhr_abort("as_estimation_spec: a dynhr_model takes no further ",
                 "arguments; edit its spec with update(dm, ...) or ",
                 "dynhr_estimation_spec(dm, ...).",
                 class = "dynhr_error_bad_argument")
  validate_spec(x$spec)
}

#' @rdname as_estimation_spec
#' @param entry For a list: \code{"run_full_estimation"},
#'   \code{"run_posterior_estimation"} or \code{"run_mode_finding"}.
#' @export
as_estimation_spec.list <- function(x, entry = c("run_full_estimation",
                                                 "run_posterior_estimation",
                                                 "run_mode_finding"), ...) {
  entry <- match.arg(entry)
  .spec_assemble(.spec_parts_flat(x, entry)$parts)
}

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.dynhr_run_record <- function(x, ...) {
  .spec_assemble(.spec_parts_record(x)$parts)
}

.spec_from_result <- function(x) {
  if (!inherits(x$run_record, "dynhr_run_record"))
    .dynhr_abort("as_estimation_spec: this result carries no run record.",
                 class = "dynhr_error_no_run_record")
  as_estimation_spec.dynhr_run_record(x$run_record)
}

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.dynhr_mode_result <- function(x, ...) .spec_from_result(x)

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.dynhr_posterior_result <- as_estimation_spec.dynhr_mode_result

#' @rdname as_estimation_spec
#' @export
as_estimation_spec.dynhr_estimation_result <- as_estimation_spec.dynhr_mode_result

## Dynare estimation() options with no bearing on results.
.spec_dynare_display_opts <- c("nograph", "nodisplay", "graph_format", "tex",
  "noprint", "console", "nocorr", "plot_priors", "mode_check",
  "mode_check_neighbourhood_size", "mode_check_symmetric_plots",
  "mode_check_number_of_points", "conf_sig", "mh_conf_sig", "silent_optimizer",
  "nodiagnostic")

.spec_unquote <- function(v)
  if (is.character(v)) gsub("^['\"]|['\"]$", "", v) else v

## Map a .mod's estimation(...) options to spec field lists.
## Returns list(likelihood, mode, sampler (list, or FALSE = no sampling),
## data (path or NULL), unmapped (character)).
.spec_map_dynare <- function(opts, mod_dir) {
  lik <- list(); mode <- list(); samp <- list(); data <- NULL
  unmapped <- character(0)
  num1 <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  nm <- names(opts)
  if (is.null(nm)) nm <- character(0)
  method <- "rwmh"
  if ("posterior_sampling_method" %in% nm) {
    psm <- .spec_unquote(opts$posterior_sampling_method)
    method <- switch(psm, random_walk_metropolis_hastings = "rwmh",
                     hssmc = "smc", dsmh = "dsmh", NA_character_)
    if (is.na(method)) {
      unmapped <- c(unmapped, "posterior_sampling_method")
      method <- "rwmh"
    }
  }
  samp$method <- method
  if ("mh_replic" %in% nm && num1(opts$mh_replic)) {
    rep_ <- as.integer(opts$mh_replic)
    if (rep_ == 0L) {
      samp <- FALSE
    } else if ("n_draws" %in% .spec_sampler_table[[method]]) {
      drop <- if ("mh_drop" %in% nm && num1(opts$mh_drop)) opts$mh_drop else 0.5
      samp$n_warmup <- as.integer(floor(drop * rep_))
      samp$n_draws  <- rep_ - samp$n_warmup
    } else {
      unmapped <- c(unmapped, "mh_replic")
    }
  } else {
    if ("mh_replic" %in% nm) unmapped <- c(unmapped, "mh_replic")
    if ("mh_drop" %in% nm) unmapped <- c(unmapped, "mh_drop")
  }
  for (k in intersect(c("mh_nblocks", "mh_nblck"), nm)) {
    if (is.list(samp) && num1(opts[[k]]) &&
        "n_chains" %in% .spec_sampler_table[[method]])
      samp$n_chains <- as.integer(opts[[k]])
    else if (!isFALSE(samp)) unmapped <- c(unmapped, k)
  }
  if ("mode_compute" %in% nm) {
    mc <- opts$mode_compute
    meth <- if (num1(mc)) switch(as.character(mc), "4" = "newrat",
                                 "9" = "cmaes", NA_character_) else NA_character_
    if (is.na(meth)) unmapped <- c(unmapped, "mode_compute") else mode$method <- meth
  }
  if ("first_obs" %in% nm) {
    if (num1(opts$first_obs)) lik$first_obs <- as.integer(opts$first_obs)
    else unmapped <- c(unmapped, "first_obs")
  }
  if ("nobs" %in% nm) {
    if (num1(opts$nobs)) lik$nobs <- as.integer(opts$nobs)
    else unmapped <- c(unmapped, "nobs")
  }
  if ("lik_init" %in% nm) {
    li <- if (num1(opts$lik_init)) switch(as.character(opts$lik_init),
            "1" = "stationary", "3" = "diffuse", NA_character_) else NA_character_
    if (is.na(li)) unmapped <- c(unmapped, "lik_init") else lik$lik_init <- li
  }
  if ("diffuse_filter" %in% nm) {
    if (isTRUE(opts$diffuse_filter)) lik$lik_init <- "diffuse"
    else unmapped <- c(unmapped, "diffuse_filter")
  }
  if ("prefilter" %in% nm && !identical(opts$prefilter, 0))
    unmapped <- c(unmapped, "prefilter")
  if ("order" %in% nm && !identical(opts$order, 1))
    unmapped <- c(unmapped, "order")
  if ("datafile" %in% nm) {
    df <- .spec_unquote(as.character(opts$datafile))
    cand <- file.path(mod_dir, c(df, paste0(df, ".csv")))
    cand <- cand[file.exists(cand) & grepl("\\.csv$", cand, ignore.case = TRUE)]
    if (length(cand)) data <- cand[[1L]] else unmapped <- c(unmapped, "datafile")
  }
  handled <- c("posterior_sampling_method", "mh_replic", "mh_drop",
               "mh_nblocks", "mh_nblck", "mode_compute", "first_obs", "nobs",
               "lik_init", "diffuse_filter", "prefilter", "order", "datafile",
               .spec_dynare_display_opts)
  unmapped <- c(unmapped, setdiff(nm, handled))
  list(likelihood = lik, mode = mode, sampler = samp, data = data,
       unmapped = unique(unmapped))
}

#' @rdname as_estimation_spec
#' @param data,obs_vars,likelihood,mode,sampler,compute,outputs,options For a
#'   \code{.mod} path: as in \code{\link{dynhr_estimation_spec}}, overriding
#'   the values mapped from the file's \code{estimation(...)} command.
#' @export
as_estimation_spec.character <- function(x, data = NULL, obs_vars = NULL,
                                         likelihood = NULL, mode = NULL,
                                         sampler = NULL, compute = NULL,
                                         outputs = NULL, options = NULL, ...) {
  if (length(x) != 1L || !file.exists(x))
    .dynhr_abort("as_estimation_spec: `x` must be the path of an existing ",
                 ".mod file.", class = "dynhr_error_spec_missing_file")
  path <- normalizePath(x)
  mod <- parse_mod(path, verbose = FALSE)
  est <- Filter(function(cm) identical(cm$name, "estimation"), mod$commands)
  opts <- if (length(est)) est[[1L]]$options else list()
  mp <- .spec_map_dynare(opts, dirname(path))
  if (length(mp$unmapped))
    .dynhr_inform("as_estimation_spec: estimation() option(s) with no dynhr ",
                  "spec equivalent, not applied: ",
                  paste(mp$unmapped, collapse = ", "), ".",
                  class = "dynhr_message_spec_unmapped_mod_options")
  merge <- function(mapped, user) {
    if (is.null(user)) return(mapped)
    if (inherits(user, "dynhr_subspec") || inherits(user, "dynhr_sampler_sequence") ||
        isFALSE(user) || !is.list(user) || is.null(names(user)))
      return(user)
    if (isFALSE(mapped)) mapped <- list()
    for (nm in names(user)) mapped[nm] <- list(user[[nm]])
    mapped
  }
  samp <- merge(mp$sampler, sampler)
  if (is.list(samp) && !inherits(samp, "dynhr_subspec") && !is.null(names(samp)) &&
      !is.null(sampler$method) && !identical(tolower(sampler$method), mp$sampler$method)) {
    ## a user-chosen method: keep only the mapped fields it takes
    keep <- c("method", .spec_sampler_table[[tolower(sampler$method)]])
    samp <- samp[intersect(names(samp), keep)]
  }
  dynhr_estimation_spec(mod,
    data = data %||% mp$data, obs_vars = obs_vars,
    likelihood = merge(mp$likelihood, likelihood),
    mode = merge(mp$mode, mode),
    sampler = samp,
    compute = merge(list(), compute),
    outputs = merge(list(), outputs),
    options = options)
}
