## R/result-changes.R
## ---------------------------------------------------------------------------
## The result-change registry (E5 C3, plan section 7, class B "code").
##
## A table of every dynhr change that alters results -- a likelihood, a
## sampler kernel, a solution, a prior law -- with the version that shipped it
## and the COMPONENT TAGS it touches. An estimation spec declares the tags it
## uses (.est_component_tags()), so for two builds the registry says whether
## the code between them changed anything THIS run depends on:
##
##   * resuming a checkpoint across such a change is refused (the kernel or
##     target changed mid-chain) unless compute$on_mismatch = "warn";
##   * dynhr_rerun() warns, naming the changes and the version to install for
##     a bit-identical replay; any other build difference is only a message;
##   * dynhr_verify() lists them and expects the results to differ.
##
## Adding an entry: every fix whose NEWS entry says results change gets a row
## here, in the same commit, tagged with the narrowest tags that cover it (add
## a tag to .dynhr_component_tag_doc -- and to .est_component_tags() when a
## spec can use it -- if none fits). A tag no spec ever declares (girf, osr,
## ...) still documents the change; it just never blocks an estimation.
## Changes that only turned a crash or a loud error into a result, or moved
## values at the floating-point level (< 1e-10 relative), are not registered.
## ---------------------------------------------------------------------------

## What each component tag means: when a spec uses it. Tags of the form
## "likelihood:<type>", "obc:<filter>", "mode:<method>" and "sampler:<method>"
## are generated from the spec as well and may appear in registry rows.
.dynhr_component_tag_doc <- c(
  ## target: model, priors, data, likelihood (every stage uses them)
  prior_rules       = "every estimation (the estimated_params prior law and bounds)",
  bk_rank           = "every estimation (the order-1 Blanchard-Kahn rank condition, a singular QZ block Z11)",
  ghu_solve         = "every estimation (the order-1 shock loading ghu solve)",
  prior_uniform     = "a uniform prior",
  prior_stderr      = "an estimated shock stderr / var prior",
  prior_corr_skew   = "an estimated corr or skew row in estimated_params",
  shock_default_var = "a shock that the shocks block does not list",
  mod_parser        = "the .mod parse (documenting only: a spec keeps the parsed model by value and its content hash is class A)",
  deriv_special     = "model equations using normcdf / normpdf / erf / steady_state()",
  obs_trends        = "observation_trends / deterministic_trends",
  hetero_shocks     = "heteroskedastic_shocks (model block, override or plan)",
  perturbation_ho   = "a higher-order (pruned) solution",
  ss_linear         = "model(linear)",
  kalman            = "the Gaussian Kalman-filter likelihood",
  dare_missing      = "the Kalman likelihood with measurement error and missing data",
  whittle_lik       = "the Whittle likelihood",
  student_t         = "the Student-t filter likelihood",
  pskf_order1       = "the order-1 PSKF likelihood",
  pskf_mean         = "the PSKF shock mean (correlated skewed shocks)",
  pskf_cdf          = "the PSKF multivariate-normal CDF evaluation (pruning compensation, Phi_2 .. Phi_7)",
  obc_pwl           = "the piecewise-linear OBC solution (any OBC model)",
  obc_filter        = "the OBC Kalman (PKF) likelihood",
  obc_ppf           = "the OBC particle filters (PPF / COPF)",
  obc_mcp           = "mcp-tagged complementarity constraints",
  ms_struct         = "the structural Markov-switching filter",
  ms_regime_ss      = "structural Markov switching with regime steady states",
  ramsey_mult       = "a Ramsey model (lagged-multiplier state columns)",
  ramsey_foc        = "a Ramsey model (FOCs through model-locals)",
  planner_discount  = "a Ramsey / welfare step (planner discount)",
  ## mode stage
  grad_intercept    = "the analytic posterior gradient (measurement-constant scores)",
  grad_power        = "the analytic gradient / Hessian under power_posterior or system priors",
  grad_eqmap        = "the analytic solution-derivative gradient (implicit / adjoint)",
  grad_init         = "the posterior gradient's filter initialisation (lik_init) and hybrid FD base value",
  kalman_singular_F = "the Kalman filter's singular-F test (univariate fallback) and steady-state lock",
  kalman_diffuse    = "the exact-diffuse Kalman initialisation (unit-root models, lik_init = diffuse / auto)",
  whittle_grad      = "the analytic Whittle gradient",
  exact_hessian_opt = "mode$exact_hessian = TRUE (the use_exact_hessian option)",
  proposal_bounds   = "the Step-6 proposal covariance at a mode on a prior bound",
  proposal_hessian  = "the Step-6 proposal Hessian (Sigma_prop / V_mode)",
  rfe_options       = "a full run (run_full_estimation) with non-default mode / RWMH options",
  ## sampler stage
  nuts_mass         = "NUTS / HMC / ChEES diagonal mass adaptation",
  hmc_jitter        = "the HMC step-size jitter",
  chees_mass        = "ChEES checkpoint resume (inverse mass)",
  dsmh_streams      = "DSMH chain random streams",
  smc               = "the SMC sampler (target, stage-0 draws, THAMES)",
  parallel_seed     = "parallel chains with a seed other than 42",
  sampler_repeat    = "a sampler sequence that repeats a method",
  ## outputs, and tools outside estimation (documenting only)
  chain_summaries   = "a sampler stage (the chains' ESS / R-hat summaries and the THAMES marginal-likelihood SE)",
  posterior_outputs = "IRFs and moments at the posterior mean (outputs$stoch_simul)",
  diagnostics       = "the diagnostic battery (outputs$diagnostics)",
  moments_irfs      = "compute_moments() / compute_irfs() and variance decompositions outside estimation",
  theta_names       = "a posterior / gradient closure called directly with an unnamed or reordered theta (numDeriv, optim)",
  simulation        = "simulate_model() and the order-2 / order-3 simulators",
  cond_forecast     = "conditional forecasts",
  girf              = "compute_girf() order-2 GIRFs",
  obc_pf            = "the perfect-foresight OBC path solver",
  ramsey_obc        = "Ramsey with occasionally-binding constraints",
  osr               = "optimal simple rules",
  discretion        = "discretionary_policy()",
  benchmark         = "dynhr_benchmark()")

## The registry. One row per result-changing change: `version` (the release
## that shipped it), `components` (tags) and `description`.
##
## 0.9.3.8-0.9.3.33 had no per-version NEWS; their rows are backfilled from
## the commit messages (git log fdf271b8..2c6974b0) and the consolidated
## 0.9.4 NEWS. 0.9.3.8-0.9.3.12 (the messaging layer: conditions, cat() onto
## the message stream) and the report-rendering releases 0.9.3.30 / 0.9.3.32
## change no results and have no rows; 0.9.3.29 (run_diag = TRUE crashed
## after sampling) turned a crash into a result and is not registered.
.dynhr_result_changes <- local({
  rows <- list(
    ## ---- the D0-D41 adversarial diagnostics refresh (0.9.3.13-0.9.3.21) ----
    list("0.9.3.13", "diagnostics",
         "D0 scales rows / columns before its SVD (false FAIL on badly scaled equations) and checks the steady-state residual; D1 uses a finite-difference-aware rank tolerance (exactly unidentified models passed); D3 divides downward Morris steps by the signed step"),
    list("0.9.3.15", "diagnostics",
         "D4's posterior path is no longer always NA; D5 computes rank-normalised split R-hat and bulk / tail ESS matching posterior 1.7.0 (its ESS was ~3x too low, a scale-only mismatch passed); D6 no longer renormalises the prior over the posterior range (overlap inflated 15x) and builds it from log_prior()"),
    list("0.9.3.16", "diagnostics",
         "D7 no longer counts failed starts as basins; D8 re-solves at the posterior mean (its IRFs were calibration values); D9's posterior path compared 0 observables; D10 fractions read 100x small; D11 runs on historical_decomposition() output; D12 standardises by the per-period smoothed-shock sd (its gate failed 100% of correct models)"),
    list("0.9.3.17", "diagnostics",
         "D13 / D15 Gamma_k gain the missing Z T^(k-1) R Sigma_e D' term and no longer transpose VAR(p >= 2) lag blocks; D14 compares the best model with the runner-up (was the weakest competitor); D16 uses a nested-subsample Hausman test (false-FAILed 84% of stable models); D17 is wired into run_all_diagnostics(); D18 adds the per-period gap"),
    list("0.9.3.18", "diagnostics",
         "D19 compares order-1 and order-2 IRFs for the same shock vector (gap 4.5 against 0.02); D20, D23 use the finite-difference-aware rank (unidentified models passed); D23 (FAIL on every run) and D24 (always degenerate) re-solve per theta; D21 / D22 match parameters and observables by name"),
    list("0.9.3.19", "diagnostics",
         "D25 and D29 re-solve per theta (their Jacobians were zero); D27 re-solves the binding policy at each FD point (FAIL for every OBC model); D28 solves the ergodic regime distribution directly (power iteration was wrong for persistent chains); D26 computes the miscalibration derivative; D30 runs a genuine mpfr chain (JuliaCall backend removed)"),
    list("0.9.3.20", "diagnostics",
         "D31 runs (it never did) and compares constrained with unconstrained paths; D32 decodes the regime bitfield per constraint; D33 evaluates calibration maps in the allowlist sandbox; D34 Cochran-Q / nested Hausman tests; D35 uses a HAC score meat (plain OPG was blind to serial correlation); D36 LR intervals"),
    list("0.9.3.21", "diagnostics",
         "D37 uses Komunjer-Ng Prop. 3-S in the singular case (Delta_U was dropped: false PASS) and the FD-aware rank; D38 a unit-invariant log-parameter spectrum; D39 separates seven point classes and uses the shared prior sampler (its private copy had the IG1 / IG2 bug); D40 starts from the model initval; D41 re-solves at params with the estimation's measurement error"),
    ## ---- 0.9.4 bug-fix waves WS1-WS5 (0.9.3.23-0.9.3.28) ------------------
    list("0.9.3.23", c("posterior_outputs", "diagnostics"),
         "the posterior path re-solves the decision rule at the posterior mean: posterior_irfs, posterior_moments and the diagnostics' model_solve_fn used the calibration / mode-time rule (no structural parameter moved them; estimated shock stderrs were dropped)"),
    list("0.9.3.23", "bk_rank",
         "solve_perturbation() reports a singular QZ block Z11 (the Blanchard-Kahn rank condition) as bk_satisfied = FALSE with a classed warning; its pseudo-inverse rule was accepted by the likelihood before"),
    list("0.9.3.23", "obc_pwl",
         "OBC bounds are levels everywhere: the OBC-binding and Boehl paths disagreed when the steady state was non-zero"),
    list("0.9.3.23", c("ss_linear", "benchmark"),
         "model(linear) solves the affine steady state A y = -c instead of zero-filling non-zero constants (every level -- measurement constants, OBC bounds, forecasts -- was off by the constant): the SW2007 dynhr_benchmark() log posterior is -874.6253133934 (was -1968.70 with ctrend, constepinf, constelab and conster dropped; fingerprint re-recorded in 0.9.3.27)"),
    list("0.9.3.23", "ms_struct",
         "the ergodic regime distribution (Markov-switching specs, regime Ramsey) is one SVD null-space solver, replacing power iterations that never left the uniform start for persistent chains"),
    list("0.9.3.23", c("moments_irfs", "posterior_outputs", "diagnostics"),
         "order-1 variance decompositions (compute_moments(), D10) Cholesky-orthogonalise correlated shocks in declared order, as Dynare does; both IRF APIs use one shock-scale helper with params winning over dr$Sigma_e"),
    list("0.9.3.23", c("perturbation_ho", "moments_irfs"),
         "pruned order-2 lag autocovariances (compute_moments_order2(), pruned_ss_moments()) use the exact augmented-system formula: the lag-1 autocorrelation of an AR(1) state came out rho^3"),
    list("0.9.3.23", "simulation",
         "simulate_model() honours shock correlations; simulate_model_order2(pruning = FALSE) is a genuinely unpruned recursion"),
    list("0.9.3.23", "planner_discount",
         "ramsey_policy() honours a betta discount; order-1 unconditional welfare is NA with a reason; welfare_se is stored (n_periods default 2000)"),
    list("0.9.3.24", "chain_summaries",
         "one ESS / R-hat implementation (split, rank-normalised, folded R-hat; bulk / tail ESS matching posterior 1.7.0) for the posterior runner's convergence summary, sampler_diagnostics() and SBC: the old shared helper counted lag 0 twice (ESS ~3x too low) and pooled chains by concatenation"),
    list("0.9.3.24", c("chain_summaries", "smc"),
         "thames_mdd(se_method = \"iid\") includes the truncation indicator's variance (mean SE / across-seed sd 0.41, now 0.93)"),
    list("0.9.3.24", c("prior_rules", "smc"),
         "log_prior() and log_prior_density() (and the prior gradient) share one density per distribution, IG1 / IG2 with sd = Inf using the mean-preserving limits; the prior sampler (SMC, SMC2, SBC, D39) draws truncated normal / gamma / inverse gamma by inverse CDF instead of clamping, which put point masses at the bounds"),
    list("0.9.3.25", "diagnostics",
         "a WARN level (pass = TRUE, warn = TRUE); D5 WARNs at R-hat > 1.01, FAILs at >= 1.05 and WARNs when bulk / tail ESS < 100 x n_chains (the flat 1000 target is gone); D0 separates unit-root singularity (WARN) from redundant equations (FAIL); D40 flags near-unit roots and is wired in; D20's strength table gates only under weighting = \"sampling\"; D29 uses the model-implied moment covariance"),
    list("0.9.3.26", "diagnostics",
         "D8's default benchmarks check response signs only, over periods 1-8; D9's model sds include measurement-error variances; D41's variance test is kurtosis-robust (t5 innovations reject 0.071 against 0.251 at nominal 5%); D15 demeans by the decision-rule steady state and adds the lambda = Inf endpoint"),
    list("0.9.3.27", c("obc_pwl", "obc_filter"),
         "the one-period OBC regime check and the LCP solver (Lemke / Newton) compare deviations with deviation-form bounds (they used the level bounds: Lemke returned an empty spell, Newton bound every period on a non-zero-steady-state model)"),
    list("0.9.3.27", c("perturbation_ho", "moments_irfs"),
         "compute_moments_order2() variance decompositions Cholesky-orthogonalise correlated shocks like order 1 (22% of the variance was unattributed)"),
    list("0.9.3.27", "simulation",
         "simulate_model_order3(pruning = FALSE) is a genuinely unpruned cubic recursion (it used the order-1 state rule)"),
    list("0.9.3.28", "prior_corr_skew",
         "estimated corr a, b and skew <shock> rows were applied as the first shock's stderr, so the estimated correlation / skewness never reached the likelihood; they now set the correlation / skewness (overriding shocks-block expressions)"),
    list("0.9.3.28", "simulation",
         "simulate_model() for skewed shocks, simulate_model_order2 / order3 and the SBC data generator draw the PSKF likelihood's joint closed-skew-normal law (the order-2/3 simulators drew Gaussian shocks; simulate_model ignored corr for skewed shocks); the order-3 Gaussian path honours corr"),
    list("0.9.3.28", "obc_filter",
         "the PKF binding check uses the shared level-to-deviation bound conversion"),
    ## ---- prior sensitivity (0.9.3.31, 0.9.3.33) ----------------------------
    list("0.9.3.31", "diagnostics",
         "prior sensitivity builds its flat priors from the parsed estimated_params instead of re-parsing the .mod text (trailing % comments with commas were read as bounds) and reuses the informative fit's me_variance (a vector fallback gave loglik -Inf)"),
    list("0.9.3.33", "diagnostics",
         "prior sensitivity starts its flat-prior search at the informative mode (from the support midpoints NZSIM stopped short and flagged 30 of 68 parameters as prior-driven); a search ending below its start is INFO"),
    ## ---- 0.9.3.34 onwards (per-version NEWS archived in dev/) ----------------
    list("0.9.3.34", "shock_default_var",
         "a shock missing from the shocks block no longer borrows its variance from a same-named sig_/stderr_/sigma_ parameter (it has zero variance, as in Dynare)"),
    list("0.9.3.35", "whittle_lik",
         "the Whittle likelihood counts each ordinate once and compares the periodogram with the spectral density (not 2*pi times it): Whittle posteriors were too wide and the shock scale biased"),
    list("0.9.3.35", "student_t",
         "the Student-t filter adds the measurement-error term K me K' to the state-covariance update"),
    list("0.9.3.36", c("mod_parser", "deriv_special"),
         "the .mod parser no longer truncates equations or swallows malformed # lines; derivatives of x^steady_state(y), erf and 2/3-argument normcdf / normpdf fixed"),
    list("0.9.3.37", "mod_parser",
         "predetermined_variables re-timing through model-locals and lagged auxiliaries; EXPECTATION(k)(expr)"),
    list("0.9.3.38", "cond_forecast",
         "conditional forecasts use the full shock covariance"),
    list("0.9.3.38", "ramsey_foc",
         "Ramsey first-order conditions keep derivatives through model-locals; the generated Ramsey .mod is written at full precision"),
    list("0.9.3.38", "perturbation_ho",
         "models with no state variables get their order-2 to order-5 shock terms"),
    list("0.9.3.38", "ss_linear",
         "model(linear) with a singular static Jacobian is solved by minimum-norm least squares; an inconsistent system reports converged = FALSE"),
    list("0.9.3.38", "discretion",
         "discretionary_policy() requires a discount instead of using 0.99 silently"),
    list("0.9.3.39", c("prior_rules", "smc"),
         "estimated_params follows Dynare: bounds are the prior support within [LB, UB], p3/p4 are a generalised-beta support / shift / truncation, INITVAL starts the mode search; one prior sampler for SMC and the other samplers"),
    list("0.9.3.39", "prior_stderr",
         "an estimated stderr prior is renamed to sig_X only when the shock's stderr expression uses sig_X"),
    list("0.9.3.40", "grad_power",
         "the analytic posterior gradient and exact Hessian include power_posterior and system priors (gradient mode finding targeted another posterior)"),
    list("0.9.3.40", "diagnostics",
         "chain_diagnostics() uses rank-normalised split R-hat and Geyer ESS; prior CDFs and quantiles honour p3/p4"),
    list("0.9.3.42", "smc",
         "SMC tempers power_posterior and system priors and draws stage 0 from the parameter prior; THAMES re-indexes resampled particles; the adaptive lambda targets the combined-weight ESS"),
    list("0.9.3.43", "perturbation_ho",
         "order-4/5 sigma terms when a current shock enters the model nonlinearly"),
    list("0.9.3.45", "sampler_repeat",
         "a repeated sampler method is pooled once (its draws were counted twice)"),
    list("0.9.3.47", c("obs_trends", "mod_parser"),
         "observation_trends are applied to the data (they were ignored); var(log) is implemented"),
    list("0.9.3.50", "nuts_mass",
         "the diagonal mass adaptation of NUTS / HMC / ChEES sets the INVERSE mass to the warmup variance (it set the mass)"),
    list("0.9.3.51", "chees_mass",
         "ChEES checkpoints store the inverse mass (a resumed chain drifted by ~1e-15)"),
    list("0.9.3.56", c("hetero_shocks", "mod_parser"),
         "integer heteroskedastic_shocks periods are offset by first_obs > 1; model-locals with a lead or lag are re-timed"),
    list("0.9.3.57", "hmc_jitter",
         "dynhr_hmc() jitters its step size by +/-20% by default"),
    list("0.9.3.59", "obc_mcp",
         "MCP upper-bound and two-sided constraints (sign of the Fischer-Burmeister Jacobian)"),
    list("0.9.3.66", "diagnostics",
         "kalman_smoother(method = \"auto\") restarts on the univariate recursion when the innovation covariance is singular"),
    list("0.9.3.71", "osr",
         "osr() no longer imposes a hidden [0.1, 10] box and penalises unit-root loss variables"),
    list("0.9.3.76", "prior_uniform",
         "uniform_pdf P1, P2 are the mean and standard deviation, as in Dynare (they were read as bounds)"),
    list("0.9.3.78", c("mod_parser", "ramsey_obc"),
         "Dynare command options with lists are no longer corrupted; ramsey_constraints are imposed"),
    list("0.9.3.81", "girf",
         "compute_girf()'s analytic order-2 GIRF for horizons >= 2"),
    list("0.9.3.82", c("obc_pf", "ramsey_obc", "planner_discount"),
         "pf_newton_solve() releases binding OBC periods; planner_discount from the .mod is honoured"),
    list("0.9.3.83", "ms_struct",
         "the structural Markov-switching filter scores y_t with the regime in force at t"),
    list("0.9.3.84", "whittle_grad",
         "the analytic Whittle gradient no longer drops a finite-differenced parameter's slot"),
    list("0.9.3.86", c("dsmh_streams", "sampler:dsmh"),
         "DSMH chains run on per-chain seeded streams (seeded draws change)"),
    list("0.9.3.88", "ramsey_mult",
         "solve_perturbation() zeroes a state's ghx column only for a literal v = 0 equation (the Ramsey lagged multiplier kept its commitment)"),
    list("0.9.3.89", c("mod_parser", "obs_trends"),
         "trend_var / deflators are detrended at parse; deterministic_trends merge into observation_trends; var(...) with nested parentheses parses"),
    list("0.9.3.90", c("ms_struct", "ms_regime_ss"),
         "structural Markov switching with regime-dependent steady states: the exact constant k_s, and the filter's state and observation intercepts"),
    list("0.9.3.92", "obc_pwl",
         "piecewise-linear OBC paths for binding spells longer than one period and interacting constraints; OBC bounds replace the right equation"),
    list("0.9.3.93", "obc_filter",
         "the OBC Kalman filter is Dynare 7's OccBin piecewise-linear filter (OBC likelihoods change)"),
    list("0.9.3.94", "obc_ppf",
         "the OBC particle filters use the OccBin multi-period solution per particle and propagate all-missing periods"),
    list("0.9.3.96", "dare_missing",
         "kalman_filter(method = \"dare\") handles missing observations when me_variance > 0 (it returned -Inf)"),
    list("0.9.3.96", "pskf_mean",
         "the PSKF shock mean under correlated skewed shocks is the exact closed-skew-normal mean"),
    list("0.9.3.100", "benchmark",
         "dynhr_benchmark() reports the acceptance rate (it was always NA)"),
    list("0.9.3.101", "exact_hessian_opt",
         "run_mode_finding() honours dynhr_set_options(use_exact_hessian = TRUE)"),
    list("0.9.3.104", "pskf_order1",
         "the order-1 PSKF likelihood uses the exact contemporaneous state space (it used the wrong law whenever an observable loads a current shock)"),
    list("0.9.3.105", "grad_intercept",
         "the analytic gradient keeps the likelihood score of parameters that only shift an observable's steady state (every grad_method)"),
    list("0.9.3.106", "grad_eqmap",
         "analytic solution-derivative gradients (implicit / adjoint / adjoint_solution) use the system's compound-LHS equation ordering"),
    list("0.9.3.107", "parallel_seed",
         "parallel RWMH / NUTS / SMC / DIME chains honour the seed (they always used seed_base = 42)"),
    list("0.9.3.107", "rfe_options",
         "run_full_estimation() honours transform_params, rwmh_*, proposal_cov_method and use_exact_hessian"),
    list("0.9.3.109", "grad_eqmap",
         "grad_method option default is now \"auto\" (adjoint_solution for the Gaussian likelihood): sampler gradients with analytic_grad = TRUE change"),
    list("0.9.3.110", c("grad_eqmap", "hetero_shocks"),
         "cumulant adjoint_solution with estimated shock stds; shock_scale on a unit-root model no longer gives loglik = -Inf"),
    list("0.9.3.111", "grad_power",
         "the mode-stage analytic gradient receives the same target as the objective (power_posterior / likelihood extras); shadowing extras are mapped to typed fields"),
    list("0.9.3.115", "grad_init",
         "make_posterior_grad() honours lik_init (diffuse near-unit roots had the \"auto\" slope); hybrid takes its FD base value from the objective (rounding-level change; ~1e8 error on singular-F models such as art_zlb_mcp, where analytic methods are now refused)"),
    list("0.9.3.117", "grad_init",
         "the sampler-stage (NUTS/MALA/ChEES, serial, mirai, pooled) and parallel mode gradients honour lik_init; posterior_hessian() refuses a non-stationary P0 in force"),
    list("0.9.3.118", "pskf_cdf",
         "the PSKF likelihood (also pskf_smoother and the HANK truncation filter) evaluates the pruning compensation with the deterministic CDF up to dimension 5 (was Mendell-Elston from dimension 3) and checks 4-/5-dimensional Miwa CDFs, falling back to a separation-of-variables rule: up to 12-18 nats with several skewed shocks, 0.14 nat with one (option pskf_cdf = \"fast\" restores the old evaluation); round-off variances no longer enter the CSN pseudoinverse as directions"),
    list("0.9.3.119", "kalman_singular_F",
         "the singular-F univariate fallback follows Dynare's rule (rcond(F) < kalman_tol and a diag or correlation-rcond below it): small-scale well-conditioned F (art_zlb_mcp: 1227.49 -> 3280.37) no longer drops components; the steady-state lock is relative to max|P| (logliks move at the 1e-7 level on small-P models); lik_init reaches the parallel mode portfolio, newrat H0 seed, at-mode Hessian, Sigma_prop and SBC gradients"),
    list("0.9.3.121", c("kalman_singular_F", "grad_init"),
         "Lyapunov/P0 solvers (R and C++) converge on a relative tolerance (P0 last bits move; small-scale models were up to 6e-6 nats off); the Kalman score and steady-state adjoint lock where the forward filter locks (hybrid kept an early-locked shock-std score, 1e-5 relative off, on small-P models); missing-data / shock_scale / me_extra branches apply the singular-F rule (were up to ~500 nats off on relatively singular F); parallel mode finding passes system_priors to the daemons' objective and gradient"),
    list("0.9.3.122", "pskf_cdf",
         "every PSKF multivariate-normal CDF of dimension 3-7 uses a deterministic C++ tilted separation-of-variables lattice rule with log-scale error below max(1e-5, 1e-7|log p|) (checked/plain Miwa and the R lattice were 1e-3 to 1e-1 nat off in orthant tails, non-finite in deep tails); dimension-2 CDFs on the log scale (lost all mass below -8, floored at log p = -36); correlations below 1e-3 no longer zeroed; pskf_cdf = \"fast\" unchanged"),
    list("0.9.3.123", c("kalman", "kalman_diffuse", "grad_init"),
         "the exact-diffuse phase tests the unit-free F_inf / P_inf on their own (Dynare's form): with innovation variances above ~1e10 in the data's units the diffuse phase ended early (10-30 nats); the diffuse smoother and adjoint gradient follow"),
    list("0.9.3.123", c("kalman", "grad_init", "student_t", "ms_struct"),
         "lik_init = \"auto\" accepts a stationary P0 that is PSD up to round-off relative to its own scale (was an absolute -1e-8): singular-but-valid large-scale P0s no longer switch to the diffuse / kappa init (18 nats); a user P0 is validated relative to its scale"),
    list("0.9.3.123", c("perturbation_ho", "likelihood:pruned"),
         "the pruned augmented-state covariance is solved block-balanced (order 2 and 3): the first-order block of the order-3 covariance was up to 1.4e-6 relative off when shock scales exceed 1"),
    list("0.9.3.123", "likelihood:cumulant",
         "the GMM weight matrix uses a per-moment ridge (correlation form) and a relative lag truncation: the weighted cumulant likelihood no longer depends on the data's units (was 2.6-3.0 nats off at scale 1e-4)"),
    list("0.9.3.123", "likelihood:tpf",
         "the TPF initial-cloud Cholesky jitter is relative to each state's variance (fixed-seed logliks move ~1.5e-10 relative)"),
    list("0.9.3.124", c("obc_pwl", "obc_filter", "obc_ppf"),
         "OccBin regime decisions (PKF, PPF/COPF, boehl_solve_regime_path, occbin simulation) use Dynare's plain inequalities with a round-off band relative to |x|+|b| and to the tagged equation's term magnitudes, not an absolute 1e-8: bound violations below 1e-8 in level units now bind; changes only near-kink cases or small-unit models (COPF was 3.2e-2 nats off at scale 1e-4); Dynare-parity goldens unchanged"),
    list("0.9.3.126", c("obc_pwl", "obc_mcp"),
         "solve_obc_lcp (Lemke and Newton) decides relative to max|q|, max|M| and the relative OccBin round-off band (was absolute 1e-8 / 1e-6): a violation below 1e-8 in level units now binds; occbin_solve_path's Newton stops per equation relative to its term magnitudes (small-unit paths stayed at the steady state); results at natural scale unchanged"),
    list("0.9.3.126", "obc_ppf",
         "the PPF/COPF initial particle cloud falls back in the model's units (one-period shock covariance without a stationary P_0; eigen square root of a PSD-singular P_0), was diag(1e-6) / diag(1e-3): logliks change only for singular-P_0 or unit-root slack models (up to 71 nats at scale 1e-4)"),
    list("0.9.3.127", c("grad_eqmap", "grad_init"),
         "make_posterior_grad()'s own default grad_method is now \"auto\" (was \"hybrid\"), so mode finding, the Hessian, SBC and the parallel mode / NUTS entry points take the exact gradient: hybrid missed the 1e-6 exactness gate on every benchmark model (8e-4 .. 2.5e-2 relative; 0 instead of 236.5 at a bound on NZSIM) and the exact method was 1.5-6.6x faster on 6 of 7 models (replication/grad_bench_2026-09-28)"),
    list("0.9.3.128", c("obc_pf", "ramsey_obc"),
         "pf_newton_solve / ramsey_obc_pf regime decisions use Dynare's plain inequalities with a round-off band relative to |y|+|b| and the relaxed equation's term magnitudes (was absolute tol), and Newton also requires every equation's residual within tol times its term magnitude (max|R| < tol kept): small-unit models no longer stay all-slack (paths 20-99% off at 1e-6 x natural units); natural-scale results bit-identical (0 of 9,104 decisions moved)"),
    list("0.9.3.128", c("sampler:nuts", "diagnostics"),
         "dynhr_nuts()$acceptance_rate (and NUTS chain_stats$accept_rate, sampler_diagnostics) is the mean post-warmup NUTS acceptance statistic (new $accept_stat / $accept_stats / $move_rate); it was 1 - mean(treedepth == 0), identically 1; draws unchanged"),
    list("0.9.3.130", "proposal_bounds",
         "Step-6 proposal covariance at a mode on a prior bound (the central FD stencil left the support): bound parameters get their one-sided inward curvature and are decoupled, the rest the curvature with them held at the bound (was: non-finite Hessian -> .make_pd prior variances for every coupled parameter); the exact-Hessian path decouples likewise; run_mode_finding() returns $Sigma_prop_eta there. Interior modes unchanged"),
    list("0.9.3.131", "proposal_bounds",
         "transform_params samplers use run_mode_finding()'s $Sigma_prop_eta at a mode on a prior bound (RWMH/CPM proposal, parallel NUTS proposal, serial NUTS/ChEES initial mass; was the delta method ~1e6+ too wide -> RWMH froze); the estimate runner's proposal is build_sigma_prop's bound-aware Step 6"),
    list("0.9.3.131", c("sampler:nuts", "sampler:chees"),
         "serial NUTS / ChEES initial eta-space mass is the eta precision d^2 / Var(theta) (was 1 / (Var(theta) d^2): off by d^4 per parameter, ~1e8 for a log-transformed parameter near 0.01); warmup starts from the right metric"),
    list("0.9.3.132", c("sampler:nuts", "sampler:mala"),
         "metric = \"hessian\" under transform_params handed the THETA-space Sigma_prop over as the eta metric (off by d_i d_j per entry); now D^-1 Sigma D^-1 (Sigma_prop_eta at a bound mode)"),
    list("0.9.3.132", "sampler:rwmh",
         "the CPM path (tpf with cpm_rho_u) drove a THETA-space random walk with an ETA-space covariance; rwmh_cpm() now takes the transform and walks in eta"),
    list("0.9.3.132", "sampler:nuts",
         "parallel (mirai) and pooled NUTS with analytic_grad and transform_params applied the eta chain rule to the gradient twice (an eta closure fed theta); the daemons now pass the theta gradient -- draws were valid (the energy used the true target) but the trajectories were inefficient"),
    list("0.9.3.133", "sampler:nuts",
         "parallel NUTS takes its initial eta mass / chain dispersal from the same covariance as the serial path (sampler's Sigma_prop > Sigma_prop_eta at a bound mode > the mode's scaled Sigma_prop; was the unscaled V_mode first); metric = \"whittle_fim\" now genuinely falls back to the eta-converted Hessian metric (was always diagonal) and a degenerate FIM is no longer used as an identity metric"),
    list("0.9.3.134", "proposal_hessian",
         "Step-6 proposal covariance (Sigma_prop, V_mode; proposal_cov(method = \"full\")) is the symmetrised central difference of the mode stage's exact gradient when it has one (gaussian/cumulant/whittle, no OBC, non-hybrid; 2n gradient calls; was the 2n(n+1)-evaluation log-posterior stencil): more accurate (fs2000 4e-6 vs 8e-4 relative), Sigma_prop moves ~1e-3 relative; the parallel Step-6 posterior now includes the system prior (was dropped on the pool)"),
    list("0.9.3.135", c("exact_hessian_opt", "mode:newrat"),
         "posterior_hessian() (t2_method loop / contract_once, the default) looked up the second-order solution blocks by param_names position instead of structural-parameter position: whenever an estimated shock stderr preceded a structural parameter the exact Hessian was wrong (SW2007: 13% Frobenius, stderr block 33%, csigl diagonal 22.7x); affects use_exact_hessian / mode$exact_hessian proposals and the newrat analytic H0 seed"),
    list("0.9.3.136", "grad_eqmap",
         "exact gradient faster, result change at FD-noise level: dSigma_e/dtheta is the exact derivative of .get_shock_cov() (stderr / variance / corr priority chain, stats::D() on shocks-block expressions; FD kept only where D() is not exact) -- the gradient moves <= 1e-10 relative on small models, 2e-9 on NZSIM (old FD rounding amplified by 1/s^2); the parameter Jacobian is emitted sparse and contracted over its non-zeros (<= 1e-14)"),
    list("0.9.3.137", c("grad_eqmap", "likelihood:cumulant", "likelihood:pruned", "sampler:nuts"),
         "dSigma_e/dtheta in the order-2 solution derivatives / order-2 adjoint (cumulant gradient, pruned-SS chain, method-of-moments analytic Jacobian) and in the Whittle-FIM NUTS metric is the exact derivative of .get_shock_cov() (was central FD; FD kept only where stats::D() is not exact): gradients / Jacobian / FIM move <= ~2e-10 relative"),
    list("0.9.3.138", c("sampler:nuts", "sampler:hmc", "sampler:chees", "sampler:mala"),
         "fused log-posterior + gradient (attr(make_posterior_grad(...), \"logpost_grad\")) pair-checked against log_post_fn at the start (1e-10 relative): one evaluation per new leapfrog point (was two gradients + one log-posterior), gradient carried; chains move at rounding level (SW2007 NUTS 3.6e-10 after 12 iterations); the unfused path is unchanged"),
    list("0.9.3.139", c("grad_power", "sampler:nuts", "sampler:hmc", "sampler:chees", "sampler:mala"),
         "the runner's sampler gradient (serial NUTS / ChEES / MALA / HMC and the parallel NUTS daemons, independent and pooled) includes the context's system prior and power_posterior, so it is the gradient of the sampled log-posterior and the fused path applies (was: the gradient without the system prior); fusion is dropped with an infeasibility penalty; HMC takes the runner's analytic gradient (sampler_spec(\"hmc\") gains analytic_grad / grad_method); the unfused path carries the gradient (identical draws, half the gradient calls; NUTS lowrank / fisher_diag n_grad_evals drop)"),
    list("0.9.4", c("sampler:nuts", "sampler:hmc", "sampler:mala", "sampler:chees", "grad_eqmap"),
         "the runners' gradient samplers (NUTS / HMC / MALA / ChEES in run_posterior_estimation / run_full_estimation / run_estimation) default to analytic_grad = TRUE: the exact gradient wherever the likelihood has one (was a finite-difference gradient unless requested); default-argument draws change; analytic_grad = FALSE reproduces the old behaviour"),
    list("0.9.4", "mode:newrat",
         "dynhr_set_options(use_analytic_hess = FALSE) is honoured (mode spec field analytic_h0; the newrat analytic H0 seed); it was accepted and ignored"),
    list("0.9.4", "theta_names",
         "every theta-taking closure (make_posterior / make_log_posterior and every likelihood factory, make_posterior_grad and its logpost_grad companion, make_transformed_logpost / make_transformed_grad, make_loglik_contrib) reads theta by name: an UNNAMED theta of length nrow(prior_spec) is taken in prior_spec$name order (was: mapped onto nothing -- logprior 0 and the calibrated model's likelihood; sw2007 at the mode -2383 instead of -874.6), a permuted named theta maps by name (the gradient relabelled it positionally), and a wrong-length or mis-named theta is an error of class dynhr_error_theta_names. Named prior-order calls, and so every spec-driven run, are unchanged"),
    list("0.9.4", "likelihood:cumulant",
         "order-3/4 cumulant log-likelihood and adjoint gradient no longer abort with solve.default(V): computationally singular when eigen(hx) returns a numerically singular eigenvector matrix (an exactly repeated eigenvalue, e.g. fs2000's double structural zero at alp = 0.34 / 0.35601 / 0.40): those points use a doubling solve (matches a dense Kronecker solve to machine precision); order 4 no longer inverts blockdiag(Sigma_x, Sigma_e), so rank-deficient Sigma_x models evaluate (fs2000 errored at every theta); previously finite values unchanged"),
    list("0.9.4", "smc",
         "smc(phi_schedule = ) without approx_loglik_fn (likelihood tempering) is used as lambda_schedule (was: ignored, adaptive tempering ran silently); dynhr_smc_model_tempered() withholds phi_schedule from Stage 1; dynhr_set_options(x = NULL) un-sets an option so its registered default applies (the do.call(dynhr_set_options, old) restore idiom used to store NULL)"),
    list("0.9.4", "ghu_solve",
         "order-1 shock loading ghu = -(f_0 + f_+ P)^{-1} f_u is solved by plain LU (tol = 0, Dynare parity) however ill-conditioned; an rcond < 2.2e-16 used to route it to a silent truncated-SVD pseudo-inverse (a wrong ghu under bk_satisfied = TRUE; nk_small near-unit-root corner loglik -668108 -> -11337); an exactly singular system now sets bk_satisfied = FALSE with dynhr_warning_ghu_singular"),
    list("0.9.4", "ms_struct",
         "ms_solve per-regime ghu and fixed-point step use plain LU (tol = 0): ill-conditioned but nonsingular systems are solved instead of pseudo-inverted or declared diverged"),
    list("0.9.4", "pskf_cdf",
         "accurate-mode bivariate normal log-CDF uses Genz's BVND (C++) when p >= 1e-3 (log-scale quadrature kept below) and logcdf_ME_r's block dispatch runs in C++ (bit-identical): each bivariate log-CDF moves <= 1e-13, PSKF loglik <= 3e-12 relative (the general path's own one-ulp floor); pskf_cdf = \"fast\" unchanged; ~2.4x faster on the Reiter-HANK probe"),
    list("0.9.4", "pskf_cdf",
         "accurate-mode q = 3 normal log-CDF blocks use Genz's exact trivariate TVN (C++, mvn_logcdf3_cpp) when p >= 1e-3, det(correlation) >= 1e-8 and the adaptive rule converges (the lattice otherwise; the same rule in the C++ and R dispatch): each such call moves by the lattice's error (~3.6e-7) to within ~1e-15 of mvtnorm TVPACK; PSKF loglik moves ~1e-5 to 1e-4 toward the exact value (Reiter-HANK probe +3.8e-5, now 3e-13 from a TVPACK reference); pskf_cdf = \"fast\" unchanged; ~1.85x faster on the Reiter-HANK probe"))
  data.frame(
    version     = vapply(rows, `[[`, "", 1L),
    components  = I(lapply(rows, `[[`, 2L)),
    description = vapply(rows, `[[`, "", 3L),
    stringsAsFactors = FALSE)
})

## A version string as a package_version, or NULL when it is not one.
.dynhr_as_version <- function(v) {
  if (inherits(v, "numeric_version")) return(v)
  if (!is.character(v) || length(v) != 1L || is.na(v)) return(NULL)
  pv <- package_version(v, strict = FALSE)
  if (is.na(pv)) NULL else pv
}

## Registered result changes between two dynhr versions that touch the
## `components` a run uses (tags from .est_component_tags(); NULL = unknown,
## every change counts). A change belongs to the interval when it shipped
## after the older version and no later than the newer one -- in either
## direction: going back across a fix changes results as much as going
## forward. An unreadable older version counts every change up to the newer
## one. Returns one line per change, "<version> [<tags hit>] <description>",
## oldest first; character(0) when none. `registry`: a table shaped like
## .dynhr_result_changes (tests pass synthetic ones).
.dynhr_result_changes_between <- function(v_old, v_new, components,
                                          registry = .dynhr_result_changes) {
  a <- .dynhr_as_version(as.character(v_old %||% NA_character_))
  b <- .dynhr_as_version(as.character(v_new %||% NA_character_))
  if (is.null(b)) return(character(0))
  if (!is.null(a) && a == b) return(character(0))
  reg <- registry
  if (!nrow(reg)) return(character(0))
  rv  <- package_version(reg$version)
  lo  <- if (is.null(a)) NULL else if (a < b) a else b
  hi  <- if (is.null(a) || a < b) b else a
  inside <- rv <= hi & (if (is.null(lo)) TRUE else rv > lo)
  hit <- lapply(reg$components, function(tg)
    if (is.null(components)) tg else intersect(tg, components))
  keep <- which(inside & lengths(hit) > 0L)
  keep <- keep[order(rv[keep])]
  vapply(keep, function(i)
    sprintf("%s [%s] %s", reg$version[[i]], paste(hit[[i]], collapse = ", "),
            reg$description[[i]]), character(1))
}

# ---------------------------------------------------------------------------
# The component tags a spec uses
# ---------------------------------------------------------------------------

## Mode methods whose serial stage uses the analytic posterior gradient, and
## the likelihoods make_posterior_grad() serves there (run-mode-finding.R).
.est_grad_mode_methods <- c("combined", "nelder", "cmaes_jade", "newrat",
                            "cmaes_newrat")
.est_grad_likelihoods  <- c("gaussian", "cumulant", "whittle")

## The component tags an estimation spec uses. Target tags (model, priors,
## data, likelihood) are always included; `stages` adds those of the mode
## stage, the sampler stage and the outputs. A checkpoint resume passes
## "sampler" only (the continued chain runs no mode stage; its target is
## covered by the target tags); dynhr_rerun() and dynhr_verify() pass all.
.est_component_tags <- function(spec, stages = c("mode", "sampler", "outputs")) {
  lik <- spec$likelihood
  mod <- spec$model$mod
  sl  <- .spec_sampler_list(spec$sampler)
  obc_on <- .spec_obc_active(spec)
  eff <- if (obc_on) lik$obc_filter else lik$type
  pw_or_sp <- !isTRUE(all.equal(lik$power_posterior, 1)) ||
    !is.null(lik$system_priors)
  ep <- mod$estimated_params
  prior_col <- if (is.data.frame(ep)) as.character(ep$prior) else character(0)
  shocks_listed <- mod$shocks$variances$name
  eq_txt <- unlist(lapply(mod$equations, `[[`, "text"), use.names = FALSE)
  po <- mod$planner_objective
  ramsey_model <- is.list(po) && length(po$text) == 1L && !is.na(po$text) &&
    nzchar(po$text)
  hs <- mod$heteroskedastic_shocks$scales
  me_pos <- any(lik$me_variance > 0)
  data_na <- if (is.null(spec$data$value)) NA else anyNA(spec$data$value)

  tg <- c(
    paste0("likelihood:", eff),
    if (obc_on) c(paste0("obc:", lik$obc_filter), "obc_pwl",
                  if (identical(lik$obc_filter, "pkf")) "obc_filter" else "obc_ppf"),
    if (.spec_has_mcp(mod)) "obc_mcp",
    "prior_rules", "bk_rank", "ghu_solve",
    if (any(grepl("uniform", prior_col, fixed = TRUE))) "prior_uniform",
    if (is.data.frame(ep) && any(ep$type %in% c("stderr", "var"))) "prior_stderr",
    if (is.data.frame(ep) && any(ep$type %in% c("corr", "skew"))) "prior_corr_skew",
    if (length(setdiff(mod$varexo_names, shocks_listed))) "shock_default_var",
    if (any(grepl("normcdf|normpdf|erf|steady_state\\(", eq_txt))) "deriv_special",
    if (!is.null(mod$observation_trends)) "obs_trends",
    if ((is.data.frame(hs) && nrow(hs) > 0L) || !is.null(lik$heteroskedastic_shocks) ||
        !is.null(lik$plan)) "hetero_shocks",
    if (eff %in% c("pruned", "tpf", "global_pf") ||
        isTRUE(spec$model$max_order > 1L)) "perturbation_ho",
    if (isTRUE(mod$model_options$linear)) "ss_linear",
    if (identical(eff, "gaussian")) c("kalman", "kalman_singular_F",
      if (isTRUE(lik$lik_init %in% c("auto", "diffuse"))) "kalman_diffuse"),
    if (identical(eff, "gaussian") && me_pos && !isFALSE(data_na)) "dare_missing",
    if (identical(eff, "whittle")) "whittle_lik",
    if (identical(eff, "student_t")) "student_t",
    if (identical(eff, "pskf")) c("pskf_order1", "pskf_mean", "pskf_cdf"),
    if (!is.null(lik$ms_struct_spec)) c("ms_struct", "ms_regime_ss"),
    if (ramsey_model) c("ramsey_mult", "ramsey_foc", "planner_discount"))

  ## the analytic-gradient tags of one gradient use with `grad_method`
  grad_tags <- function(grad_method)
    c("grad_intercept", "grad_init", if (pw_or_sp) "grad_power",
      if (!identical(grad_method, "hybrid")) "grad_eqmap",
      if (identical(eff, "whittle")) "whittle_grad")
  full_form <- identical(.est_form(spec), "full")

  if ("mode" %in% stages) {
    md <- spec$mode
    tg <- c(tg, paste0("mode:", md$method))
    if (!obc_on && md$method %in% .est_grad_mode_methods &&
        lik$type %in% .est_grad_likelihoods &&
        !isFALSE(md$options$use_analytic_grad))
      tg <- c(tg, grad_tags(.dynhr_option_registry[["grad_method"]]$default))
    if (isTRUE(md$exact_hessian)) tg <- c(tg, "exact_hessian_opt")
    tg <- c(tg, "proposal_bounds", "proposal_hessian")
    if (full_form && (!isTRUE(md$transform_params) || isTRUE(md$exact_hessian) ||
                      !identical(md$proposal_cov, "diagonal")))
      tg <- c(tg, "rfe_options")
  }
  if ("sampler" %in% stages && length(sl)) {
    ms <- vapply(sl, `[[`, "", "method")
    tg <- c(tg, paste0("sampler:", ms))
    if (any(ms %in% c("nuts", "hmc", "chees"))) tg <- c(tg, "nuts_mass")
    if ("hmc" %in% ms) tg <- c(tg, "hmc_jitter")
    if ("chees" %in% ms) tg <- c(tg, "chees_mass")
    if ("dsmh" %in% ms) tg <- c(tg, "dsmh_streams")
    if ("smc" %in% ms) tg <- c(tg, "smc")
    if (anyDuplicated(ms)) tg <- c(tg, "sampler_repeat")
    grad_ok <- !obc_on && .ctx_allows_analytic_gradient(
      list(gradient_policy = lik$gradient_policy, likelihood = lik$type,
           me_variance = lik$me_variance))
    for (s in sl) {
      if (grad_ok && isTRUE(s$analytic_grad))
        tg <- c(tg, grad_tags(s$grad_method %||% "auto"))
      if (full_form && (isFALSE(s$transform_params) || isTRUE(s$adapt_cov) ||
                        isTRUE((s$n_blocks %||% 1L) != 1L)))
        tg <- c(tg, "rfe_options")
    }
    seed <- spec$compute$seed
    if (isTRUE(spec$compute$parallel) && !is.null(seed) && !identical(as.integer(seed), 42L) &&
        any(ms %in% c("rwmh", "pmmh", "nuts", "smc", "dime")))
      tg <- c(tg, "parallel_seed")
  }
  if ("outputs" %in% stages) {
    out <- spec$outputs
    ## the chains' ESS / R-hat / THAMES summaries are computed after sampling
    ## from all draws, so they belong to the outputs: a checkpoint resume
    ## ("sampler" only) is not refused over them
    if (length(sl)) tg <- c(tg, "chain_summaries")
    if (isTRUE(out$stoch_simul)) tg <- c(tg, "posterior_outputs")
    if (isTRUE(out$diagnostics)) tg <- c(tg, "diagnostics")
    if (isTRUE(out$ramsey)) tg <- c(tg, "ramsey_mult", "ramsey_foc", "planner_discount")
  }
  unique(tg)
}
