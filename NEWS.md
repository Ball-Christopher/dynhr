# dynhr 0.9.4

This is a large correctness, performance and API release, following a full
package review. It consolidates development versions 0.9.3.8 to 0.9.3.139.
The dev version of each change is given in parentheses. The per-version
notes for 0.9.3.34 onwards are archived in the development repository as
`dev/NEWS-0.9.3-dev.md` (not shipped with the package); 0.9.3.8 to 0.9.3.33
(the messaging layer, the diagnostics refresh and the report fixes) had no
per-version notes and are summarised here from the commit history.

Many fixes change results: likelihoods, gradients, samplers, proposals, OBC
paths, priors, moments and diagnostics. **Estimates made with 0.9.3.7 or
earlier should be re-run**, and in particular anything using PSKF, OBC
filters, structural Markov switching, the Whittle likelihood, the analytic
gradient, bounded `estimated_params` rows, or estimated `corr` / `skew`
rows. The machine-readable list is the result-change
registry, `dynhr:::.dynhr_result_changes`: each row gives the version and
the components it touches. Estimation results now carry a run record, and
`dynhr_verify()` / `dynhr_rerun()` name the registered changes that separate
a result from the installed build (see *New features*).

The release carries breaking API changes; read the first section before
upgrading.

## Breaking changes

### Messages, warnings and progress output (0.9.3.8-0.9.3.12)

- **Progress output is now on the message stream, not stdout.** 587 `cat()`
  calls that printed progress inside computations (the `verbose = TRUE`
  text of `solve_model()`, the samplers, ...) now emit messages, so
  **`capture.output()` / `sink()` no longer capture them**; use
  `capture.output(type = "message")`, `suppressMessages()` or the verbosity
  level. Objects printing themselves (`print` / `summary` / `format`
  methods, and report printers such as `verify_steady_state()` and
  `hmc_summary()`) still write to stdout (0.9.3.9, 0.9.3.10). The moments
  tables of `stoch_simul(verbose = TRUE)` are now one message (0.9.3.12),
  and `profile_ci()` prints one line per grid point instead of a dot
  (0.9.3.11).
- **Warnings and messages are classed conditions** (`dynhr_warning`,
  `dynhr_message`, `dynhr_error`, plus a subclass where given). They are
  base R conditions, so `suppressWarnings()`, `tryCatch()` and
  `expect_warning()` work as before, and message text is unchanged
  (0.9.3.8).
- **One verbosity level:** `dynhr_set_verbosity()` / `dynhr_verbosity()`,
  one of `"silent"`, `"error"`, `"warn"`, `"info"` (default) or `"debug"`.
  The default keeps everything that printed before; `"warn"` mutes progress
  and `"silent"` also drops warnings. Errors are never suppressed. The level
  can override a per-call `verbose = TRUE`, `verbose = FALSE` still does not
  silence warnings, and the level reaches mirai daemons (0.9.3.8, 0.9.3.9).
- **"Warn once" is per run, not per session.** The old latches were never
  reset, so a second run in one session emitted none of those warnings. Now
  12 entry points (`run_full_estimation()`, `run_posterior_estimation()`,
  `run_mode_finding()`, `run_estimation_passport()`, `dynhr_sbc()`,
  `dynhr_smc2()`, `dynhr_benchmark()`, `dynhr_model()`, `solve_model()`,
  `method_of_moments()`, `forecast_backtest()`, `conditional_forecast()`)
  each open a run: a repeated warning is shown once per run, and the run
  reports how many repeats it suppressed when it ends (0.9.3.8).

### Renamed or removed

- **`mcmc()` is now `dynhr_mcmc()`** (it masked `coda::mcmc`), and **`tune()`
  is now `filter_tune()`** (it masked `e1071::tune` and the tidymodels `tune`
  package). There are no aliases; the `"mcmc"` class is unchanged (0.9.3.75).
- **`dm_test()` is now `diebold_mariano_test()`**, class
  `dynhr_diebold_mariano_test`. No alias (0.9.3.129).
- **`hank_filter_shocks()` takes `me_variance`** instead of `me_sd`
  (0.9.3.58).
- **`boehl_spell_trajectory()` is removed** and `obc_simulate()` has a new
  signature (0.9.3.92).
- **`obc_guess_verify()` loses its observation pre-pass and `obs_tol`**;
  `regime_path_init` now only seeds period 1, and the OBC filters gain
  `horizon = 200` (0.9.3.93).
- **D30 drops its `"JuliaCall"` backend**: `run_all_diagnostics(d30_backend
  =)` takes `"auto"`, `"Rmpfr"` or `"base"` (0.9.3.19).

### The `dynhr_model` / `dm_*` redesign (0.9.3.129)

- A `dynhr_model` now holds a `dynhr_estimation_spec` (`dm$spec`). Its fields
  (`dm$model`, `dm$data`, `dm$obs_vars`, `dm$priors`, ...) are read-only
  through `$` / `[[`; edit them with `update(dm, ...)`, which also clears
  every dependent cached result.
- `dm_mode()` and `dm_sample()` run their stage through `run_estimation()`,
  so the spec defaults apply (4 chains, 10000 / 5000 draws). `$mode` is a
  `dynhr_mode_result` (log posterior at `$mode$mode$logpost`) and a new `$fit`
  slot holds the posterior result. `dm_sample()` loses `theta0`, and unknown
  `...` names are an error. Neither verb needs `dm_posterior()` first.
- `dm_posterior()` builds the spec runner's objective, so the `.mod`'s
  heteroskedastic_shocks / filter_tunes / stochastic-volatility / OBC blocks
  are now applied (they were ignored).

### Priors and `.mod` semantics

- **Uniform priors follow Dynare:** `uniform_pdf, P1, P2` is now **mean and
  standard deviation** (support P1 +/- sqrt(3) P2), not the bounds. Write
  bounds as `NAME, uniform_pdf, , , LB, UB;`. Giving P1/P2 together with
  P3/P4, giving neither, or a non-positive sd is an error; the silent [0, 1]
  default is gone. `dynhr_warning_uniform_parameterisation` flags rows that
  were probably written as bounds (0.9.3.76).
- **`estimated_params` follows Dynare's rules.** `LB`/`UB` bounds are now
  honoured (they were ignored). `p3`/`p4` are no longer truncation bounds for
  every shape: BETA is a generalised beta on `[p3, p4]`, GAMMA and the
  inverse gammas take `p3` as a shift, NORMAL is truncated, UNIFORM uses
  `[p3, p4]`. **A short-form row `NAME, SHAPE, p1, p2, lo, hi;` that meant
  truncation now means a generalised or shifted prior; convert it to the
  long form `NAME, p1, lo, hi, SHAPE, p1, p2;`** as the bundled models were
  (0.9.3.39).
- `INITVAL` / `estimated_params_init` now give `run_mode_finding()`'s default
  start (0.9.3.39).
- A shock missing from the shocks block gets zero variance, as in Dynare; it
  no longer borrows a same-named `sig_`/`stderr_`/`sigma_` parameter
  (0.9.3.34).
- Integer `heteroskedastic_shocks` periods are now offset by an integer
  `first_obs > 1` (0.9.3.56).

### Changed defaults

- **The exact gradient is the default everywhere.** The `grad_method` option
  defaults to `"auto"` (0.9.3.109) and so does `make_posterior_grad()`'s own
  formal, `run_nuts_mirai()` and the parallel mode tasks (0.9.3.127): mode
  finding, the Hessian and SBC change. `dynhr_set_options(grad_method =
  "hybrid")` restores the old sampler behaviour.
- **The estimation runners' gradient samplers use the exact gradient by
  default** (0.9.4). NUTS / HMC / MALA / ChEES in `run_posterior_estimation()`,
  `run_full_estimation()` and `run_estimation()` default to
  `analytic_grad = TRUE` (the spec field `sampler_spec(...)$analytic_grad`
  too); before, they used a finite-difference gradient unless asked. The
  exact gradient is used wherever the likelihood has one; other likelihoods
  keep the numerical gradient without a warning (the former
  `dynhr_warning_spec_analytic_grad_unavailable` is gone). Default-argument
  draws change; `analytic_grad = FALSE` reproduces the old behaviour. The
  low-level `nuts()` / `dynhr_nuts()` without `grad_fn` still use finite
  differences.
- `dynhr_hmc()` jitters its step size by +/-20% (`step_jitter = 0.2`);
  `step_jitter = 0` reproduces the old chain (0.9.3.57).
- PSKF pruning compensation: `offset_miwa_qmax` defaults to 5 (0.9.3.118).
- `kalman_smoother(method = "auto")` restarts on the univariate recursion
  when F is singular (0.9.3.66).
- `me_variance` defaults to `NULL` (the option) in `run_full_estimation()`,
  `run_mode_finding()` and `dynhr_model()` (0.9.3.111).
- `run_mode_finding(use_exact_hessian = NULL)` (was `FALSE`, which masked the
  option) (0.9.3.101).
- `run_mcmc_mirai()`, `run_nuts_mirai()` and `run_mode_mirai()` take
  `seed_base = NULL`, resolved from `options(seed_base)` (42 when unset)
  (0.9.3.114).

### Changed return values

- `dynhr_nuts()$acceptance_rate` is now the mean post-warmup NUTS acceptance
  statistic; it used to be identically 1 (0.9.3.128).
- `run_estimation()` and its wrappers skip the mode stage when every sampler is
  prior-initialised (SMC, DSMH, DIME, SMC²) and nothing needs the mode:
  `$mode` is then `NULL` (0.9.3.111).
- `osr()`: `result$loss_weights` is the weight matrix (0.9.3.71).
- Structural MS: `c_const` is now the steady-state-gap constant k_s, a named
  length-n_endo vector (0.9.3.90).
- Run-record hashes carry their algorithm as a prefix (`"sha256:"` /
  `"md5:"`) (0.9.3.102).
- Diagnostic results have a WARN level (`pass = TRUE, warn = TRUE`),
  rendered by every printer and report; several diagnostics' thresholds and
  badges changed (see *Diagnostics*) (0.9.3.25).
- `ramsey_policy()`: order-1 unconditional welfare is `NA` with a reason, and
  `welfare_se` is stored (`n_periods` default 2000) (0.9.3.23).

### Now an error or refused (was silently ignored or wrong)

- **Theta is read by name everywhere (0.9.4).** Every theta-taking closure
  (`make_log_posterior()` / `make_posterior()` for every likelihood,
  `make_posterior_grad()` and its fused companion, `make_transformed_logpost()` /
  `make_transformed_grad()`, `make_loglik_contrib()`) maps theta through one
  helper: an UNNAMED theta of length `nrow(prior_spec)` is taken in
  `prior_spec$name` order (as `numDeriv` / `optim` pass it), a named theta in
  any order maps by name, and a wrong length or names that do not match
  `prior_spec$name` (including a partial named subset) raise
  `dynhr_error_theta_names`; a length-0 theta is the explicit "evaluate at the
  calibrated parameters" call and passes through. Before, an unnamed theta
  mapped onto nothing -- the closure silently scored the calibrated model with
  logprior 0 (SW2007: -2383 against -874.6) -- and the gradient closures
  relabelled a permuted named theta by POSITION, evaluating it at the wrong
  point. An unnamed `theta_init` to the mode finder now keeps its names and
  the prior-bound box constraints (they were dropped).

- A `.mod` expression calling anything outside the sandbox allowlist aborts
  the parse (`dynhr_error_unsafe_mod_expression`) instead of yielding `NA`
  (0.9.3.24); see *Security*.
- `.mod` parsing: leftover tokens and unknown characters abort
  (`dynhr_error_mod_syntax`) instead of truncating the equation (0.9.3.36);
  undeclared `predetermined_variables` and nested `EXPECTATION` abort
  (0.9.3.37); `adl()` is an error (0.9.3.72); an unbalanced growth model
  aborts (`dynhr_error_unbalanced_growth`, stricter than Dynare) (0.9.3.89).
- Student-t likelihood: `student_df <= 2` aborts; `me_extra`, `shock_scale`,
  `known_shocks`, `filter_tunes` and heteroskedastic shocks are refused
  (`dynhr_error_student_t_unsupported`) (0.9.3.35).
- A trended model (`observation_trends`) is refused by every path except the
  Gaussian Kalman filter (0.9.3.47).
- A vector `me_variance` is refused where it is not supported
  (`dynhr_error_me_variance_vector`) (0.9.3.96).
- `discretionary_policy()` and the Ramsey / welfare functions error when no
  discount is found (`dynhr_error_discount_missing`,
  `dynhr_error_no_discount`); they used a silent 0.99 (0.9.3.38, 0.9.3.82).
  `osr()` driven by `planner_objective` requires `discount=` (0.9.3.78).
- `osr()` no longer imposes a hidden [0.1, 10] box, and a unit-root loss
  variable is penalised instead of counting as zero loss (0.9.3.71).
- Estimation validation: OBC with the cumulant likelihood, an mcp model with
  a Kalman-only likelihood, an unknown method, `checkpoint_dir` with SMC, and
  `resume` without a directory are errors (0.9.3.107); spec `extra` keys that
  shadow a typed field are refused (`dynhr_error_spec_shadowed_field`)
  (0.9.3.111).
- Analytic gradients refuse `lik_init = "kappa"` and singular-F models
  (`dynhr_error_grad_lik_init`, `dynhr_error_grad_singular_F`);
  `posterior_hessian()` refuses a non-stationary P0 (0.9.3.115, 0.9.3.117).
- Resuming a checkpoint with a different target, or across a registered
  result change, is refused unless `on_mismatch = "warn"` (0.9.3.107,
  0.9.3.108).
- Parallel workers running a different build from the session abort
  (`dynhr_error_worker_version_skew`): **run `R CMD INSTALL` after changing
  the source, before any parallel run** (0.9.3.41).
- Programming errors (subscript out of bounds, object not found, ...) are
  re-raised instead of becoming a silent numerical fallback, and malformed
  user text (e.g. a `planner_objective` naming an undefined symbol) errors
  (0.9.3.87, 0.9.3.97-0.9.3.99).
- **Arguments that used to pass silently through `...` (0.9.4).**
  `make_log_posterior()` / `make_posterior()` reject a `...` name that no
  likelihood takes (`dynhr_error_unknown_argument`; `verbose` is accepted
  and ignored). A name another likelihood takes (e.g. `order = 99L` with the
  Gaussian likelihood), and `freq_band`, `pruned_order`, `student_df` or
  `lik_init` set away from their defaults for a likelihood that does not use
  them, are ignored with a warning (`dynhr_warning_inapplicable_argument`).
  `run_mode_finding()` forwards its `...` to this constructor, so it rejects
  the same names. `run_posterior_estimation()` and `run_full_estimation()`
  reject the retired spellings `nburn` / `ndraws` / `nchains` /
  `nparticles` / `nwalkers` (and `run_full_estimation()` also `Y`,
  `obs_names`, `observables`, `me_var`, `me_sd`) with
  `dynhr_error_retired_argument`, naming the current argument; they used to
  fail deep inside a sampler with an unclassed "unused argument" (and
  `me_var` was partial-matched to `me_variance`).
- `smc(phi_schedule = )` without `approx_loglik_fn` (likelihood tempering)
  used to be ignored, so the run was adaptive; it is now an alias for
  `lambda_schedule` (a one-time `dynhr_message_schedule_alias` says so), and
  giving both with different values is an error
  (`dynhr_error_schedule_conflict`). `smc_model_tempered()` keeps
  `phi_schedule` for its bridge stage only (0.9.4).
- `dynhr_set_options(name = NULL)` un-sets the option. The documented
  restore idiom `old <- dynhr_set_options(x = v); do.call(dynhr_set_options,
  old)` returns `NULL` for an option that was unset, and used to STORE that
  `NULL`, which then shadowed the registered default (e.g. `grad_method`
  read as `NULL`). The options live in dynhr's own store, not in R's
  `options()` (`getOption("dynhr.allow_monge_metric")` is always `NULL`);
  read them with `dynhr_get_options(effective = TRUE)`.
  `allow_monge_metric` still gates `metric = "monge"` (0.9.4).

### Seeded output changes

Seeded chains differ from 0.9.3.7 even where the target is unchanged: NUTS /
HMC / ChEES (mass adaptation, 0.9.3.50; initial mass, 0.9.3.131), gradient
samplers with `analytic_grad = TRUE` (gradient default, 0.9.3.109), HMC (step
jitter, 0.9.3.57), DSMH (per-chain streams, 0.9.3.86), and parallel chains,
which now honour `seed` (0.9.3.107).

## Results that change

Sizes are those measured in the dev notes; see the registry for the
component tags.

### Security

- **Reading a `.mod` file could execute arbitrary R code.** Calibration,
  initval / histval, shocks, verbatim and `@#define` expressions were
  evaluated in a `baseenv()` child, so `system()` / `unlink()` in a `.mod`
  ran at parse time. One allowlist sandbox now covers the parser and D33
  (0.9.3.24), and also `filter_tunes` / `heteroskedastic_shocks` values,
  `steady_state_model`, shocks-block expressions and OBC / OccBin / MCP
  bounds, checked once at parse so MCMC draws pay nothing extra (0.9.3.34).

### Posterior outputs

- **The posterior path reused the calibration / mode decision rule**, so
  `posterior_irfs()`, `posterior_moments()` and the diagnostics' re-solve
  saw calibrated dynamics. They now re-solve at the posterior mean
  (0.9.3.23).
- **`run_full_estimation(run_diag = TRUE)` crashed after sampling**, every
  time. It now solves at the posterior mode and passes the solution, priors,
  observables, IRFs and per-chain draws to `run_all_diagnostics()`
  (0.9.3.29).

### Kalman filter and initialisation

- **Singular-F test** used an absolute 1e-10 pivot cut, so a
  well-conditioned F at small scale went univariate: art_zlb_mcp scored
  1227.49 against the dense 3280.37, and rescaling the data changed the
  answer. It now follows Dynare 7.1's rcond rule (0.9.3.119). The
  missing-data, `shock_scale` and `me_extra` branches apply the same rule
  (they were 249 and 503 nats off) (0.9.3.121).
- **Steady-state lock and Lyapunov / P0 tolerances are relative** (the
  absolute lock left 0.66 nat at scale 1e-3) (0.9.3.119, 0.9.3.121).
- **Exact-diffuse phase ended early at large data scales** (10.1 nats at
  scale 1e4, 28.8 at 1e6); `lik_init = "auto"` switched a valid singular P0
  to diffuse / kappa (18 nats) (0.9.3.123).
- **`method = "dare"` returned -Inf on any missing observation when
  `me_variance > 0`**, which also hit every `return_ll_contrib = TRUE` call
  (0.9.3.96).
- **Student-t filter and `kf_innovation_diagnostics()` dropped `K me K'`**
  (0.6-1.7 nats at `me_variance = 0.1`) (0.9.3.35).
- `shock_scale` on a unit-root model gave -Inf; it is now supported through
  the univariate filter (0.9.3.110).
- `kalman_smoother()` silently dropped zero-variance components (7.9e-3
  adding-up residual); `historical_decomposition()`'s coherence checks had a
  `max(1, .)` floor that hid small-scale incoherence (0.9.3.66, 0.9.3.123).
- Pruned state space: the augmented covariance is solved block-balanced
  (1.4e-6 relative at scale 100); the TPF initial-cloud jitter is relative to
  each state's variance (0.9.3.123).
- **Conditional forecasts ignored shock correlations** (0.9.3.38).

### Gradients and Hessians

- **Compound-LHS equations broke the analytic solution-derivative
  gradients** (implicit / adjoint / adjoint_solution, orders 1 and 2):
  SW2007 `crhoa` -36.9 against 88.97, max relative error 21.8; fs2000 order-2
  and pruned gradients off by up to 21x (0.9.3.106).
- **Parameters that only shift an observable's steady state**
  (`y_obs = y + mu`) returned only their prior score under every
  `grad_method`: SW2007 `constepinf` -15.97 against -55.28. Gradient mode
  finding could stall (0.9.3.105).
- **The old default `hybrid` gradient** missed the 1e-6 exactness gate on
  every benchmark model (8e-4 to 2.5e-2 relative) and returned 0 instead of
  236.5 for an NZSIM parameter at its bound. The exact gradient is now the
  default; fixture modes moved by up to 2.5e-6 in theta and every optimiser
  reaches the same, slightly higher mode (0.9.3.127).
- **`posterior_hessian()`'s exact Hessian misplaced its second-order
  blocks** whenever an estimated stderr preceded a structural parameter:
  SW2007 13% Frobenius error, 22.7x on the `csigl` diagonal. Affects
  `use_exact_hessian` proposals and the newrat H0 seed (0.9.3.135).
- **Gradient and exact Hessian ignored `power_posterior` and system
  priors** (0.9.3.40); the mode stage's gradient also ignored Whittle
  `debias` and cumulant orders / weights (0.9.3.111); parallel mode finding
  and the parallel Step-6 posterior dropped `system_priors` (0.9.3.121,
  0.9.3.134); the runner's sampler gradients omitted the system prior
  (0.9.3.139).
- **`lik_init` did not reach the gradient**: every gradient was the
  `"auto"` one (d/drho 7.747 against the objective's 9.565 under
  `"diffuse"`) (0.9.3.115, 0.9.3.117, 0.9.3.119). Hybrid now takes its FD
  base value from the objective, and the gradient recursions lock where the
  filter locks (0.9.3.121).
- **Step-6 proposal Hessian** is the central difference of the exact
  gradient (fs2000 relative error 4e-6 against 8e-4) (0.9.3.134); at a mode
  on a prior bound it no longer steps across the bound and hands prior
  variances to coupled parameters (0.9.3.130).
- Cumulant orders within 1:2 get an exact gradient (was FD, ~1e-5 relative)
  (0.9.3.115); the shock-covariance derivative is exact (FD-noise-level
  change) (0.9.3.136, 0.9.3.137).
- `dynhr_set_options(use_exact_hessian = TRUE)` was never honoured
  (0.9.3.101).

### Samplers and proposals

- **The default diagonal mass adaptation in NUTS, HMC and ChEES was
  inverted**, squaring the target's conditioning: on a 10-d correlated
  Gaussian the posterior means were 19 MCSE off and the covariance 48% off
  (0.9.3.50).
- **Theta/eta errors under `transform_params`:** the serial NUTS / ChEES
  initial mass was off by d^4 (~1e8) (0.9.3.131); the dense `"hessian"`
  metric was off by d_i d_j (~1e3), the CPM path walked in the wrong space,
  and parallel NUTS applied the eta chain rule twice (0.9.3.132); parallel
  NUTS started from a different metric than serial, and the `whittle_fim`
  fallback was always diagonal (0.9.3.133).
- **At a mode on a prior bound, RWMH froze** (acceptance 0) because the
  delta-method proposal was ~1e6 too wide; samplers now use the bound-aware
  `$Sigma_prop_eta` (acceptance 0.33) (0.9.3.131).
- **SMC sampled a different posterior**: it ignored `power_posterior`,
  dropped system priors in the pruned / PSKF / Student-t factories, and drew
  bounded betas from a different law (sd 0.173 against 0.158) (0.9.3.39,
  0.9.3.42). THAMES paired draws with the wrong particle and ignored
  ellipsoid mass outside the support (~+0.05 nats per boundary parameter),
  and model-tempered SMC gave Inf / NaN on infeasible M0 particles
  (0.9.3.42).
- **Parallel chains ignored `seed`** (seeds 31 and 99 gave identical draws),
  and `run_full_estimation()` ignored `transform_params`, `rwmh_*`,
  `proposal_cov_method` and `use_exact_hessian` (0.9.3.107).
- **`dm_*`:** the default RWMH chain was frozen (acceptance 0, 1 unique row
  in 60; now 0.66), and SMC diagnostics treated the weighted cloud as
  equal-weight (means off by up to 0.65 sd) (0.9.3.129).
- A repeated sampler method was pooled twice (0.9.3.45); sequences with
  different draw counts crashed the convergence step (0.9.3.108).
- **One ESS / R-hat implementation** (split, rank-normalised, folded R-hat;
  bulk / tail ESS, matching `posterior` 1.7.0) replaces three; the old
  shared helper double-counted lag 0 (ESS ~3x too low) and pooled chains by
  concatenation (0.9.3.24). `chain_diagnostics()` used the wrong
  autocorrelation pairing, capped ESS at N and had no rank-normalisation
  (0.9.3.40).
- `thames_mdd(se_method = "iid")` dropped the truncation indicator's
  variance (mean SE / across-seed sd 0.41, now 0.93) (0.9.3.24).

### PSKF and skewed shocks

- **The order-1 PSKF used the wrong law whenever an observable loads a
  contemporaneous shock** (nearly every model): with zero skew it did not
  reduce to the Kalman filter (AR(1) -25.11 against -11.38). PSKF estimates
  from earlier versions are biased (AR(1) mode (0.84, 0.21) against
  (0.74, 0.27)); re-run them. The conditional-forecast PSKF path had the same
  bug (0.9.3.104).
- **Correlated skewed shocks were mis-centred** (true mean 0.36 against 0.76
  in the example) (0.9.3.96).
- **Pruning compensation and CDF accuracy:** 9.3 / 18.4 / 26.6 nats off an
  unbiased particle filter at T = 50 / 100 / 200 with two skewed shocks, now
  within 0.02 nat (0.9.3.118); the 3-7 dimensional CDFs are a deterministic
  C++ evaluator (the previous one was 1e-3 to 1e-1 nat off in orthant
  tails), dimension-2 CDFs no longer floor at log p = -36, and small
  correlations are no longer zeroed (0.9.3.122).
- **`pskf_smoother()`** used an absolute 1e-10 jitter (1e-6 relative gain
  error) and silently dropped the skewness correction on error (0.9.3.113);
  its smoothed means at t < T were off by up to 0.6 posterior sd and are now
  exact given the retained latents (0.9.3.116). Linearly dependent
  noise-free observables no longer give -Inf or a Lapack error (0.9.3.116).

### OBC, OccBin and perfect foresight

- **Piecewise-linear OBC paths were wrong for spells longer than one
  period** (next period assumed slack; interacting constraints missed): a
  3-period ZLB spell was off by 0.019, a 6-period spell by 0.18. Affects
  `boehl_*`, `solve_obc_lcp`, `compute_irfs_obc`, `obc_simulate`,
  `occbin_solve_path(method = "pwlinear")` and `ramsey_obc_pwlinear`. OBC
  specs could also replace the wrong equation when .mod order and
  declaration order differ (0.9.3.92).
- **The OBC Kalman filter is now Dynare 7's OccBin piecewise-linear filter**:
  the old one gave 665.63 where Dynare gives 860.786. OBC posteriors,
  smoother, decomposition and inversion filter change, and the likelihood is
  no longer path-dependent on the previous draw (0.9.3.93).
- **The OBC particle filters (PPF / COPF)** now solve each particle's
  multi-period regime sequence: COPF gave -18010 where the PKF gives 394.76,
  and bootstrap and COPF differed by 32.8 nats. All-missing periods no
  longer slip the time index (0.9.3.94).
- **`pf_newton_solve()` never released a binding period**, so spells could
  come out too long; `ramsey_obc_pwlinear()` kept a one-period ZLB episode
  binding for the whole horizon (0.9.3.82).
- **MCP upper-bound and two-sided constraints never solved** (the failure
  was loud) (0.9.3.59).
- **OBC bounds are levels everywhere.** The Boehl and OBC-binding paths
  disagreed when the steady state was non-zero (0.9.3.23); the regime and
  LCP solvers compared deviations with level bounds (Lemke returned an empty
  spell, Newton bound every period on a non-zero-steady-state toy)
  (0.9.3.27); the PKF binding check shares the bound conversion (0.9.3.28).
- **Regime decisions are scale-free** (plain inequalities with a relative
  round-off band, as in Dynare) in the PWL check, PKF, PPF / COPF, Lemke LCP,
  `occbin_solve_path()`, `pf_newton_solve()` and `ramsey_obc_pf()`. At small
  units a 0.5% bound violation did not bind, Lemke returned z = 0, and
  perfect-foresight paths were 20-99% off. Natural-scale Dynare-parity cases
  are unchanged (0.9.3.124, 0.9.3.126, 0.9.3.128).

### Markov switching

- **The structural MS filter scored y_t with the previous regime's
  observation rows** (0.9.3.83), and **ignored regime-dependent steady
  states** (off by 24-8374 nats against exact enumeration; `c_const` was
  wrong) (0.9.3.90). Structural log-likelihoods, regime probabilities,
  smoothed states and `ms_smoothed_fit_struct()` change; re-run structural
  MS estimates. Reduced-form MS filters are unchanged.
- The ergodic regime distribution is one SVD solver, replacing two
  power-iteration copies, which were wrong for persistent chains (0.9.3.19,
  0.9.3.23).

### Moments, simulation and IRFs

- **Pruned order-2 lag autocovariances were wrong** (the lag-1
  autocorrelation of an AR(1) state came out rho^3) (0.9.3.23);
  `pruned_ss_moments3()` gains lag autocovariances (0.9.3.27).
- **Variance decompositions orthogonalise correlated shocks** by Cholesky
  in declared order, as Dynare does, at order 1 (0.9.3.23) and in
  `compute_moments_order2()`, where 22% of the variance was unattributed
  (0.9.3.27).
- **Simulators:** `simulate_model()` honours shock correlations (0.9.3.23),
  and for skewed shocks too (0.9.3.28); the order-2 / order-3 simulators
  and the SBC data generator drew Gaussian shocks for skewed models and now
  use the PSKF likelihood's joint closed-skew-normal law; the order-3
  Gaussian path
  honours `corr` (0.9.3.28). `simulate_model_order2(pruning = FALSE)`
  (0.9.3.23) and `simulate_model_order3(pruning = FALSE)` (which used the
  order-1 state rule) are genuinely unpruned (0.9.3.27).
- Both IRF APIs use one shock-scale helper, with `params` winning
  (0.9.3.23).

### Whittle, cumulant and GMM

- **Order-4 cumulant values before 0.9.2.0004 were corrupt; the identity
  weight is not unit-free.** Before 0.9.2.0004 the model's order-4 cumulant
  block was silently recycled against the sample's (250 entries against 16)
  and the sample fourth cumulant was divided by T+1: on a linear 3-observable
  model the 1:4 log-likelihood read -1591 where the exact identity criterion
  is -61941. Cumulant order-4 posteriors or evidence from before 0.9.2.0004
  are invalid. The default `cumulant_weight = "identity"` criterion is in the
  data's units (the order-k block scales as sigma^(2k)), so it is not
  comparable across units or models; for evidence pass
  `weight_matrix = estimate_gmm_weight_matrix(..., method = "analytic")`,
  which now builds the order-3/4 blocks (Gaussian long-run variances via
  Isserlis' theorem) instead of erroring (0.9.4). On that model the order-4
  term becomes -0.37 nats.
- **Cumulant likelihood at repeated eigenvalues (0.9.4).** The order-3/4
  cumulant log-likelihood and its adjoint gradient aborted with
  `solve.default(V): computationally singular` wherever `eigen(hx)` returned a
  numerically singular eigenvector matrix -- an exactly repeated eigenvalue,
  e.g. fs2000's double structural zero at scattered `alp` values (0.34,
  0.35601, 0.40). Those points now use a doubling solve (machine precision
  against a dense Kronecker solve); previously finite values are
  bit-identical. Order 4 no longer inverts blockdiag(Sigma_x, Sigma_e), so a
  rank-deficient Sigma_x (fs2000: rank 2 of 4) evaluates -- it errored at
  every theta.

- **The Whittle likelihood was wrong twice:** half-weighted ordinates, and
  (with `debias = FALSE`) a 2 pi scale error that biased the shock std down
  by sqrt(2 pi). Whittle posteriors were far too wide (0.9.3.35).
- The analytic Whittle gradient errored whenever a finite-differenced
  parameter preceded an analytic one (0.9.3.84); the cumulant
  `adjoint_solution` gradient crashed with estimated shock stds (0.9.3.110).
- The cumulant GMM weight matrix depended on the data's units (2.6-3.0 nats
  at scale 1e-4) (0.9.3.123); `method_of_moments(weighting = "diagonal")`
  weights were off by 2 at scale 1e-4 (0.9.3.124).

### HANK

- **`hank_determinacy()` uses the sequence-space winding-number
  criterion**; the old `rcond(H_U)` test called an indeterminate operator
  (e.g. a Taylor rule with phi < 1) determinate (0.9.3.54).

### Priors

- **Estimated `corr a, b` and `skew <shock>` rows were applied as the first
  shock's stderr**, so the estimated correlation or skewness never reached
  the likelihood (0.9.3.28).
- `log_prior()` and `log_prior_density()` share one density per
  distribution; IG1 / IG2 with sd = Inf use the mean-preserving limits; the
  prior sampler draws truncated normal / gamma / inverse gamma by inverse
  CDF instead of clamping (which put point masses at the bounds) (0.9.3.24).
- An estimated `stderr eps_X` prior was renamed to `sig_X` whenever such a
  parameter existed, even when the shocks block did not use it, so the
  estimated std had no effect on the likelihood (0.9.3.39).
- One prior sampler draws exactly what the density scores (SMC, SMC², SBC,
  DIME, prior predictive, prior sensitivity) (0.9.3.39); prior CDFs and
  quantiles honour `p3`/`p4` (0.9.3.40).
- **Modes and draws produced before 0.9.3.39 may lie outside the
  now-enforced `estimated_params` support.** The `LB` / `UB` bounds were
  ignored before, so an old mode with a parameter below its `LB` now has
  log posterior `-Inf` (the SW2007 stress variant's old mode had
  `ctou = 0.00094` against an `LB` of 0.005). Prior draws hand-rolled from
  the untruncated laws also differ from the bounded prior: SW2007's
  indeterminate share of the prior is 3.1% under the untruncated laws and
  0.45% under the bounded prior (`crpi >= 1`). Re-draw from the prior with
  the package's sampler, and re-find old modes, before comparing.

### Diagnostics

**All 41 D-series diagnostics were adversarially reviewed and fixed**
(0.9.3.13-0.9.3.22); about half had been silently degenerate (zero
Jacobians from a fixed decision rule, never wired, or placeholder paths).
Expect different badges. The largest changes:

- **Never ran or always degenerate:** D11 on `historical_decomposition()`
  output, D17 (not wired into `run_all_diagnostics()`), D23 (FAIL on every
  run), D24, D25, D29 (zero Jacobian) and D31 (0.9.3.16-0.9.3.20); D27 was
  FAIL for every OBC model (0.9.3.19).
- **Identification:** D1, D20, D23 and D37 use a finite-difference-aware
  rank, so exactly unidentified models no longer pass; D37 gave a false PASS
  in the singular case (0.9.3.13, 0.9.3.18, 0.9.3.21).
- **Posterior path:** D4 was always NA, D8's IRFs were calibration values,
  and D9 compared 0 observables (0.9.3.15, 0.9.3.16).
- **False verdicts:** D12's gate failed 100% of correct models, D16's
  false-FAILed 84% of stable ones, D14 compared the best model with the
  weakest competitor instead of the runner-up, D6's overlap was inflated
  15x, and D5's ESS was ~3x too low (0.9.3.15-0.9.3.17).
- **Wrong formulas:** D13 / D15's `Gamma_k` missed a term and transposed
  VAR(p >= 2) lag blocks; D19 compared order-1 and order-2 IRFs for
  different shock vectors (4.5 against 0.02); D35's plain OPG was blind to
  serial correlation (0.9.3.17, 0.9.3.18, 0.9.3.20).
- **New thresholds, checked against the literature** (0.9.3.25,
  0.9.3.26): D5 WARNs at R-hat > 1.01 and FAILs at >= 1.05, and WARNs when
  bulk or tail ESS is below 100 x n_chains (the flat 1000 target is gone);
  D0 separates unit-root singularity (WARN) from redundant equations
  (FAIL); D40 flags near-unit roots and is wired in; D8 checks response
  signs only, over periods 1-8, by default; D9 includes measurement-error
  variances; D20's strength table gates only with `weighting = "sampling"`;
  D29 uses the model-implied moment covariance (size 0.06 against 0.21);
  D41's variance test is kurtosis-robust (t5 innovations reject 0.071
  against 0.251 at nominal 5%); D15 demeans by the decision-rule steady
  state and adds the lambda = Inf endpoint.
- **Prior sensitivity** no longer re-parses the `.mod` text (trailing
  `%` comments with commas were read as bounds) (0.9.3.31), and starts its
  flat-prior search at the informative mode: on NZSIM the search from the
  support midpoints stopped short and flagged 30 of 68 parameters as
  prior-driven. A search that ends below its start is reported as INFO
  (0.9.3.33).
- D26's calibrate-vs-estimate ranking uses 90% prior intervals as plausible
  ranges (0.9.3.64); D41 and D37 use relative tests (0.9.3.124).
- `dynhr_sbc()` reports failed replications (the parallel path lost every
  message) (0.9.3.87); `dynhr_benchmark()` always reported
  `accept_rate = NA` (0.9.3.100); `run_mode_finding()$V_mode` was always NULL
  (0.9.3.45).

### Solvers and parser

- **The first-order shock loading `ghu` could be silently wrong (0.9.4).**
  `solve_perturbation()` solved `ghu = -(f_0 + f_+ P)^{-1} f_u` with R's
  `solve()`, which refuses any system with rcond below 2.2e-16, and fell back
  to a silent truncated-SVD pseudo-inverse. The result was a rank-deficient
  `ghu` with O(1) model-equation residuals, reported with
  `bk_satisfied = TRUE`. At a near-unit-root corner of a small NK model
  (rcond 7e-19) the `g` and `z` rows were (-0.10, 0.28, -0.29) and
  (0.11, -0.29, 0.86) instead of unit vectors, the Kalman fallback skipped 653
  of 900 observations, and optimisers were drawn to a spurious likelihood far
  above the genuine mode. `ghu` is now solved by plain LU (as Dynare's
  backslash) however ill-conditioned; only an exactly singular system falls
  back, with a `dynhr_warning_ghu_singular` warning, `ghu_singular = TRUE` and
  `bk_satisfied = FALSE`. `ms_solve()`'s per-regime `ghu` and its fixed-point
  step (which declared ill-conditioned steps "diverged") use the same rule.
- **The parser could turn a typo into a different model**: `y = 2*z) +
  5*z(-1);` parsed as `y = 2*z` (0.9.3.36).
- **`predetermined_variables`** gave wrong solutions with a `#` local or a
  `k(-1)` term (0.9.3.37).
- **Ramsey FOCs lost every derivative through a `#` model-local**; the
  generated Ramsey `.mod` rounded to 7 digits (0.9.3.38).
- **`solve_perturbation()` zeroed the ghx column of states** whose own rows
  were zero, e.g. the Ramsey lagged multiplier (a discretion-like rule)
  (0.9.3.88).
- **Order-4/5 sigma terms were wrong** when a current shock entered
  nonlinearly (`ghss2` 0.0507 against 0.0048) (0.9.3.43); models with no
  states got zero higher-order terms at orders 2-5 (0.9.3.38).
- **`model(linear)` constants were zero-filled** in the steady state; it now
  solves the affine system. The SW2007 benchmark log posterior is
  -874.6253133934 (it was -1968.70 with `ctrend`, `constepinf`, `constelab`
  and `conster` dropped) (0.9.3.23, 0.9.3.27). With a singular Jacobian it
  reported convergence with a non-zero residual; ordered QZ accepted a
  failing `dgges` INFO (0.9.3.38).
- `solve_perturbation()` flags a singular Z11 as `bk_satisfied = FALSE` with
  a classed warning (0.9.3.23).
- `ramsey_policy()` honours `betta` (0.9.3.23);
  `run_full_estimation(ramsey_discount =)` warns and skips on models with
  neither `beta` nor `betta` instead of aborting the estimation (0.9.3.27).
- **`compute_girf()`'s analytic order-2 GIRF was wrong for h >= 2** (up to
  75% of the peak on rbc2shock) (0.9.3.81).
- Dynare command options with lists (`instruments=(i,tau)`) were corrupted,
  block options with lists hid the block, and `var(...)` with nested
  parentheses was dropped (0.9.3.78, 0.9.3.82, 0.9.3.89).
- `planner_discount` is honoured (0.9.3.82).

## New features

### Estimation specs, runner and reproducibility

- **`dynhr_estimation_spec()`** built from `likelihood_spec()`,
  `mode_spec()`, `sampler_spec()`, `compute_spec()` and `outputs_spec()`,
  with validation, `print`, `update()` and content hashes;
  `validate_spec()`, `as_estimation_spec()` (from flat arguments, a run
  record or a .mod), `diff_specs()`, and `write_spec()` / `read_spec()` in
  YAML, JSON and RDS (0.9.3.103). New fields `mode$run` (0.9.3.114) and
  `likelihood$pskf_cdf` (0.9.3.120).
- **`run_estimation(spec)`** is the single code path;
  `run_full_estimation()`, `run_posterior_estimation()` and
  `run_mode_finding()` are thin wrappers (0.9.3.107).
  `run_full_estimation()` accepts every sampler and likelihood and gains
  `checkpoint_dir`, `resume` and `on_mismatch` (0.9.3.108).
- **Run records and replay:** results carry `$run_record` (arguments, data
  and model, every option, RNG state, provenance, hashes); `dynhr_rerun()`
  replays bit-identically (0.9.3.101). Records carry the full spec
  (schema 2) and the resolved gradient method (0.9.3.107, 0.9.3.110,
  0.9.3.114).
- `run_mode_finding(theta_init =)` sets the optimiser start (0.9.3.33).
- **`dynhr_verify(result)`** reports PASS/FAIL for target, code and
  environment, with deterministic and distributional checks (0.9.3.108).
- **Result-change registry** and checkpoint integrity checks (0.9.3.107,
  0.9.3.108).
- **Option registry:** every global option has a documented default, and
  `dynhr_get_options(effective = TRUE)` lists the values in force
  (0.9.3.101).
- **`dm_*` verbs over the spec**, `update(dm, ...)`, and
  `at = c("params", "mode", "posterior_mean")` on `dm_solve()`, `dm_irf()`
  and `dm_forecast()` (0.9.3.129).

### Samplers

- `sbc_uniformity_test()` and a new `sbc_ranks()` are exported for bespoke SBC
  harnesses that run their own sampler (0.9.4). Ranks lie in `0:L`, where
  `L = length(seq.int(1L, n_keep, by = thin))` is the number of thinned draws;
  that `L` is what `sbc_uniformity_test(L = )` expects (`floor(n_keep / thin)`
  can be one smaller).
- `dynhr_dsmh()`: Dynamic Striated Metropolis-Hastings, a port of Dynare
  7.1's `dsmh.m`, also `run_posterior_estimation(methods = "DSMH")`
  (0.9.3.68); parallel via `n_cores` (0.9.3.86, 0.9.3.87).
- Pooled multi-chain NUTS warmup, `dynhr_set_options(nuts_adapt =
  "pooled")` (0.9.3.77).
- NUTS metrics `"lowrank"` (0.9.3.49) and `"fisher_diag"` (0.9.3.57), passed
  through `run_*_estimation(metric =)` including the parallel path
  (0.9.3.57).
- NUTS `$accept_stat`, `$accept_stats` and `$move_rate` (0.9.3.128);
  `rwmh_cpm(transform =)` (0.9.3.132); `sampler_spec("hmc")` gains
  `analytic_grad` / `grad_method` (0.9.3.139).
- `grad_method = "auto"` picks the exact method per likelihood and reports
  it (`[NUTS] gradient method: auto -> adjoint_solution`) (0.9.3.109,
  0.9.3.110). `make_posterior_grad()` and `posterior_hessian()` gain
  `lik_init`, `power` and `system_priors` (0.9.3.40, 0.9.3.115, 0.9.3.117).

### Likelihoods and filters

- `kalman_filter()` reports `$diagnostics$n_dropped_informative` and
  `dropped_informative_by_period`, and warns (class
  `dynhr_warning_dropped_observations`, once per estimation run) when the
  univariate filter skips a component with a zero forecast variance whose
  innovation is not negligible -- the likelihood then omits real data. The
  loglik is unchanged. The standard-to-univariate fallback message no longer
  claims exact agreement when components were skipped, and `n_dropped` now
  also counts skips under the steady-state lock (`univariate_ss`) (0.9.4).
- Per-observable `me_variance` (a length-n_obs vector) in `kalman_filter()`,
  `kalman_smoother()`, PSKF and `make_log_posterior()` (0.9.3.96).
- `kalman_smoother(method = c("auto", "durbin-koopman", "univariate"))`
  (0.9.3.66).
- IMM filtering for Markov switching, `collapse = "imm"`, reduced-form
  (0.9.3.60) and structural (0.9.3.83).
- `pskf_cdf = "fast"` restores the old PSKF CDF evaluation (0.9.3.118).

### Dynare parity

- `observation_trends` and `var(log)` (0.9.3.47); `trend_var`,
  `log_trend_var`, `var(deflator=)`, `var(log_deflator=)` and
  `deterministic_trends` (0.9.3.89). No Dynare construct remains on the
  "unsupported" list (0.9.3.89).
- Complementarity tags, `shock_paths`, dated `heteroskedastic_shocks`, and
  leads / lags of `#` locals (0.9.3.56); `model_replace` / `model_remove` /
  `var_remove`, `diff()`, unquoted equation tags, and a macro processor with
  its own interpreter (never evaluated as R) (0.9.3.72);
  `EXPECTATION(k)(expr)` (0.9.3.37).
- `perfect_foresight_controlled_paths` (0.9.3.73).
- Optimal policy: `optim_weights` (with cross terms), `osr_params_bounds`,
  `osr` from `planner_objective` alone, multi-instrument
  `discretionary_policy()`, `ramsey_constraints` (0.9.3.71, 0.9.3.78), and
  `ramsey_obc_pwlinear()` with `ramsey_constraints` (0.9.3.92).
- D0 names each redundant-equation relation, as Dynare's
  `model_diagnostics` does (0.9.3.70).

### Diagnostic reports (0.9.3.30-0.9.3.32)

After an adversarial review of the HTML / PDF reports:
- The templates read one spooled contract (badges, counts, colours,
  explanations, actions, provenance) written by `write_report()`, with
  output escaped per target format.
- A methodology table covers every emitted diagnostic, with the new rules.
- Reports gain a provenance section (package version, `GIT_COMMIT`, model
  file, data dimensions and hash, observables, sampler settings), also on
  the posterior-result and `run_diag = TRUE` paths.
- WARN level in the PDF; severity-ordered sections opening on the first
  non-empty level; per-diagnostic anchors; no Google Fonts or MathJax; the
  "All diagnostics" double render is gone (22.8 MB -> ~10 MB on the scale
  check).
- Large models render: plot heights are capped, many-panel figures are
  split into pages (PDF 5 rows / 26 levels, HTML 12 / 60) up to 8 pages,
  and fixed-width callout text wraps (0.9.3.31, 0.9.3.32).
- The temporary directory is cleaned unless `dynhr.report.keep_tmp` is set.

### Analysis tools

- `variance_decomposition_nonlinear()` for pruned order-2 / order-3
  solutions (0.9.3.48), with exact polynomial Shapley above 12 shocks, and
  `historical_decomposition_nonlinear()` (0.9.3.65).
- `lp_estimand()`: what a local projection recovers from a solved model
  (0.9.3.80).
- `svar_ica()` and `match_irfs_svar_ica()`: indirect inference with a
  non-Gaussian SVAR (0.9.3.74).
- `mom_se_bounds()` and `method_of_moments(se_bounds = TRUE)` (0.9.3.62).
- `posterior_log_scores()`, `psis_lfo()`, `mdd_modified_harmonic_mean()` and
  `marginal_likelihoods()` (0.9.3.52).
- D26 ranks which parameters to calibrate and which to estimate (0.9.3.44);
  `run_all_diagnostics(log_marglik_se =)` gives D14 a 2-SE INFO band
  (0.9.3.25, 0.9.3.26); D9 adds an INFO posterior-predictive band (0.9.3.26).
- `posterior::as_draws*()` and `coda::as.mcmc*()` methods, plus
  `as.data.frame()`, `coef()`, `vcov()`, `logLik()`, `nobs()` (0.9.3.45).
- `spectral_density(normalise =)` (0.9.3.43); `dynhr_sitrep()` (0.9.3.41).

### HANK

- `hank_filter_shocks()`: sequence-space least-squares shock filter
  (0.9.3.55).
- `hank_limited_info_md()` / `hank_limited_info_ssj()`: limited-information
  estimation of a heterogeneous-agent block (0.9.3.53).
- `hank_nonlinearity_scan()` and `hank_residual_audit()` (0.9.3.61);
  `hank_steady_state_scan()` (0.9.3.70).
- Behavioural (non-FIRE) expectations in the sequence space,
  `expectations = list(theta =, gamma =, type =)` (0.9.3.95).

### Robustness

- Catch-all error handlers (about 40 in 0.9.3.87, then 187 in three sweeps)
  re-raise programming errors while keeping numerical fallbacks (0.9.3.87,
  0.9.3.97-0.9.3.99).
- `.mod` blocks a call does not apply (`filter_tunes`,
  `heteroskedastic_shocks`, `stochastic_volatility`, `mcp`) warn
  (`dynhr_warning_mod_blocks_ignored`) (0.9.3.111).
- Seeded functions restore the caller's `.Random.seed` (0.9.3.41).
- `dynhr_set_options(use_analytic_hess = FALSE)` is honoured (new mode spec
  field `analytic_h0`, the newrat analytic H0 seed); it was accepted and
  ignored (0.9.4).
- `geigen` (Suggests) is optional in the BK-distance / BK-wall tools: without
  it they fall back to LAPACK zggev via QZ (0.9.4).
- The shipped SBC coverage matrix (`inst/extdata/sbc_matrix.rds`, rendered by
  `vignette("sbc-matrix")`) is regenerated with the 2026-09 re-certification
  (0.9.4).

## Performance

Installed-build timings from the dev notes (Apple M5):

| | SW2007 | NZSIM |
|---|---|---|
| exact gradient vs old hybrid default (0.9.3.127) | 14.5 vs 80.4 ms | 45 vs 268 ms |
| exact gradient, 0.9.3.133 -> 0.9.3.134 | 14.2 -> 12.1 ms | 50.0 -> 26.4 ms |
| Step-6 proposal Hessian (0.9.3.134) | 7.36 s -> 0.68 s | 34.3 s -> 4.9 s |
| NUTS per leaf, fused lp + gradient (0.9.3.138) | 18.0 -> 8.3 ms | 43.3 -> 19.7 ms |

- PSKF accurate-mode normal CDFs (0.9.4): bivariate terms use Genz's BVND and
  trivariate terms Genz's (2004) exact TVN (Plackett reduction, adaptive
  Gauss-Kronrod, |error in p| ~1e-15) in C++ when p >= 1e-3 (TVN also needs
  det(correlation) >= 1e-8); the log-scale quadrature and the lattice are kept
  for the tails and near-singular blocks, and the block-factorisation dispatch
  runs in C++. `hank_reiter_pskf_loglik` at n_a = 16, T = 200: about 4.4x faster
  per evaluation (BVND + dispatch 2.4x, ~340 -> 144 ms; then TVN 1.85x,
  measured A/B on a loaded machine). The slowdown since
  0.9.0 (~88 ms) came from 0.9.3.104 (the corrected law keeps 3 skew
  dimensions, not 1) and 0.9.3.122 (bivariate quadrature, 1e-12 snap). The
  BVND / dispatch step moves likelihoods only at round-off (<= 3e-12
  relative); the TVN removes the lattice's per-call error (~3.6e-7), so PSKF
  logliks move ~1e-5 to 1e-4 TOWARD the exact value (the probe +3.8e-5, now
  3e-13 from a TVPACK reference).
- The exact gradient is 1.5-6.6x faster than hybrid on 6 of 7 benchmark
  models (on nk_hs2016 hybrid is 1.3x faster but 1e-3 off) (0.9.3.127).
- Unfused samplers carry the gradient: identical draws, half the gradient
  calls (e.g. 936 -> 469) (0.9.3.139).
- Exact shock-covariance derivative: one `.o2sd_dSigma_e` pass 6.8 ms ->
  0.03-0.8 ms (0.9.3.137); sparse parameter Jacobian (0.9.3.136).
- The adjoint Kalman filter solves its Lyapunov equations by doubling (was
  an O(n^6) Kronecker LU) (0.9.3.105).
- PSKF's C++ MVN CDF: the accurate filter costs ~3.7x the fast setting,
  against 5-11x for 0.9.3.118's evaluator (0.9.3.122).
- Cumulant `adjoint_solution` gradient up to ~100x faster than `implicit`
  (rbc2shock orders 1:3: 22 vs 1921 ms) (0.9.3.110).
- OBC particle filters group particles by regime guess, ~60x faster than
  solving one by one (0.9.3.94).
- Exact Shapley for 20 shocks at order 2: 0.7 s instead of 49 s of
  permutation MC (0.9.3.65); IMM is 1.5-2.3x faster than GPB(2) (0.9.3.60).
- Pooled NUTS warmup: 3.3-3.9x (`lowrank`) and 1.7-3.4x (`diagonal`)
  effective samples per gradient on a badly scaled d = 20 Gaussian
  (0.9.3.77); the `lowrank` metric is >100x the default at condition number
  1e4 (0.9.3.49).

## Validation

- **SBC re-certification** (frozen 0.9.3.87 build, AR(1)): chees x Gaussian is
  now **certified** (200 replications, min p = 0.441); hmc x Gaussian
  (min p = 0.510) and nuts x Gaussian (min p = 0.707) re-certified;
  nuts x Whittle is calibrated at T = 200 (min p = 0.346; it was < 1e-4
  before the Whittle fix) but stays *characterized* (0.9.3.125).
- **Dynare 7.1 parity:** OccBin PKF likelihood 860.78638394718621 and every
  per-period density (0.9.3.93); `occbin_solver` spells to 1e-9-1e-12
  (0.9.3.92); `perfect_foresight_solver(lmmcp)` to 2e-8 (0.9.3.82);
  `ramsey_constraints` paths to 5e-8 (0.9.3.78); Ramsey ghx/ghu to 2e-13
  (0.9.3.88); OSR loss to 1e-10 and objective 3.57149415801871 (0.9.3.71,
  0.9.3.78); multi-instrument discretion to 1e-7 (0.9.3.71); growth models
  to 1e-8 and trended log-likelihood to 1e-9 (0.9.3.47, 0.9.3.89);
  `var(log)` to 1e-11 (0.9.3.47); order-4/5 `k_order_solver` terms
  (0.9.3.43); controlled paths to 1e-7 (0.9.3.73); uniform priors
  bit-for-bit (0.9.3.76); `planner_discount` to 1e-8 (0.9.3.82).
- **Scale-equivariance audit:** every likelihood and smoother was put
  through loglik(c) + N log c = constant for data / shock scales c in
  1e-4..1e4; the absolute tolerances that broke it were made relative
  (0.9.3.119-0.9.3.128). Gaussian, Student-t, sv_rbpf, global_pf, Kim / IMM /
  structural MS, HANK, OBC PKF / PPF and PSKF were checked (0.9.3.123).
- **Independent oracles:** analytic gradients against Richardson FD on 7
  models (0.9.3.127); the exact Hessian against FD of the exact gradient to
  7.5e-8 (0.9.3.135); PSKF against a grid filter and an unbiased particle
  filter (0.9.3.118), the MVN CDF against a quadrature oracle to log p = -300
  (0.9.3.122); structural MS against exact regime-path enumeration
  (0.9.3.83, 0.9.3.90); Whittle against the exact Gaussian likelihood to
  ~0.2 nats (0.9.3.35); SV RB-PF against exact grid and rejection-sampling
  likelihoods (0.9.3.67).
- **Diagnostics:** the D5 estimators match `posterior` 1.7.0 (0.9.3.15);
  D15's lambda = Inf endpoint matches a from-scratch VAR likelihood to 1e-9
  (0.9.3.26); the corr / skew prior fix is checked against a closed-form
  AR(1) likelihood, a 601 x 601 quadrature and an independent rejection
  sampler (0.9.3.28).
- **Equivalence:** `dm_*` equals `run_estimation()` bit for bit on the same
  spec (0.9.3.129), and the thin wrappers are bit-identical to the old
  serial runners (0.9.3.107).

## Known limitations

- The OccBin PKF is about 6x slower per evaluation than the old filter
  (0.11 s against 0.018 s at T = 80) (0.9.3.93).
- PSKF: a genuine pruning error of 0.08-0.13 nat remains in heavier
  configurations (3 skewed shocks into 2 states) (0.9.3.118); the smoother
  is ~0.05 sd from the exact posterior at the default `max_q = 5` (0.9.3.116).
- Markov switching: GPB(2) can fail badly where GPB(3) does not, GPB(3) is
  not exact in general, and the one-step `collapse_max` can miss multi-period
  losses; IMM can be much worse than GPB(2) when the regimes'
  previous-period states differ a lot (0.9.3.83, 0.9.3.90).
- `observation_trends` is honoured only by the Gaussian Kalman path; mode
  finding falls back to FD gradients for trended models (0.9.3.47).
- A vector `me_variance` is refused by the other likelihoods, MS, OBC, the
  particle filters and `make_posterior_grad()` (0.9.3.96).
- The Student-t likelihood refuses `me_extra`, `shock_scale`,
  `known_shocks`, `filter_tunes` and heteroskedastic shocks (0.9.3.35).
- `lik_init = "kappa"` has no analytic gradient kernel (0.9.3.115).
- With a system prior, SMC's `log_marginal_lik` includes the system prior's
  unknown normalising constant (0.9.3.42).
- Behavioural expectations apply to linear paths only; the horizon-varying
  variant is not implemented (0.9.3.95).
- `hank_limited_info_md()`: the analytic chi-squared J test is oversized;
  prefer `se = "bootstrap"` (0.9.3.53). `svar_ica(method = "dcov")` is
  O(T^2) per evaluation; use `"fastica"` for long samples (0.9.3.74).
- `fisher_diag` is not a uniform improvement and stays opt-in (0.9.3.57).
- Not every error is classed yet: `stop()` calls get a `dynhr_error` class
  only where a file was already being changed. The default verbosity stays
  `"info"` for this release (0.9.3.8-0.9.3.12).
- Not done (measured): the Kalman adjoint runs at ~1.4x the forward filter
  (0.9.3.137); the per-session JIT warm-up is ~0.5 s once per R process
  (0.9.3.134).

# dynhr 0.9.3.7

**`dr$Sigma_e` was silently ignored by the Kalman path, and there was no way
to inject a covariance deliberately.**

Follow-up to the 0.9.3.6 report, from the real model behind it. That model had
TWO zero shock rows, and the adapter working around them patched
`dr$Sigma_e` -- which changed nothing. `solve_perturbation()` populates that
field (from its own `Sigma_e` argument when given, otherwise from the shocks
block) and `compute_irfs()` / `compute_moments()` honour it, so patching it is
a natural thing to try. `kalman_filter()` and `build_dsge_state_space()` both
recomputed from `.get_shock_cov(model, exo, params)` and ignored it, with no
override argument anywhere. The zero rows survived into `ss$Sigma_e` and
surfaced far downstream as a decomposition that would not add up.

- **`build_dsge_state_space(..., Sigma_e = )`** is the injection point: an
  explicit covariance overrides everything, unambiguously.

- **A disagreeing `dr$Sigma_e` is now reported.** The warning says it is not
  used, why, how to inject one deliberately, and that `compute_irfs()` /
  `compute_moments()` *do* honour it -- that asymmetry being the whole trap. An
  agreeing `dr$Sigma_e` says nothing, so the ordinary path is silent.

- **`ss$zero_variance_shocks`** records which shocks carry no variance, and
  `historical_decomposition()`'s incoherence warning names them. This is
  recorded rather than warned about at construction, because a per-shock
  `stderr 0` is legitimate -- it is what makes a deterministic `known_shocks`
  or `shock_means` injection meaningful. It surfaces only when something
  actually fails to add up, which is the partial-zero case 0.9.3.6's all-zero
  guard correctly stays quiet about.

**`params` remains authoritative, and that is deliberate.** Preferring
`dr$Sigma_e` was implemented first and the Markov-switching P = I oracle
rejected it within minutes: that test evaluates a fixed-regime likelihood by
reusing ONE solved decision rule while passing RESCALED `params`, and the stale
`dr$Sigma_e` made the filter ignore the rescaling -- 737 nats of error across
7 failures. Solve-once, evaluate-at-many-params is exactly what estimation
does, so `params` has to win. `kalman_filter()` is unchanged byte for byte;
there is now a regression test pinning the reuse pattern and the reason.

# dynhr 0.9.3.6

## A model with no shock variance is now called out at the source

Reported as a historical-decomposition add-up failure: a two-state linear model
with cross-loaded exact observations (`y1 = x`, `y2 = x + z`) and
`lik_init = "kappa"` gave `adding_up_residual` 7.9e-3, while the independent,
duplicate and single-observable variants balanced to 1e-17. The natural reading
-- that the cross-loaded geometry breaks the singular-innovation path -- is not
what was happening.

**The reproducer's model has no `shocks;` block, so `Sigma_e` is entirely
zero.** The model has no stochastic structure at all. With `Q = 0` the whole
path is determined by `s_0`, and under a kappa prior `s_{0|T}` is computed as
`P_{0|0} r_0` with `P_{0|0} = 1e6 * I` -- a large number times a small one --
so kappa's round-off lands directly in the initial state and `s_1` stops
equalling `T s_0`. The cross-loaded geometry only decides whether that
round-off is *visible*; the other three variants are degenerate in ways that
happen to hide it.

Give the same model a `shocks;` block and it balances to 1.1e-10 under kappa
and 1.7e-16 on the default initialisation, dropping nothing at all -- the
singular-innovation path is not even entered.

`kalman_filter()` and `kalman_smoother()` now warn when EVERY shock has zero
variance, naming the missing `shocks;` block as the usual cause. A per-shock
`stderr 0` stays legitimate and silent -- it is what makes a deterministic
`known_shocks` or `shock_means` injection meaningful; all of them being zero is
the different thing.

## `adding_up_residual` is necessary but not sufficient, and now says so

The decomposition propagates its components with `T_mat`/`R_mat` while the
add-up check rebuilds the path from the smoother's states separately -- so what
can that comparison actually see? Measured, by corrupting the smoother's output
so the answer is known by construction:

| corruption (nk_demo)               | `adding_up` | `transition` |
|------------------------------------|-------------|--------------|
| shock at t = 1                     | 2.56        | 2.57         |
| shock mid-sample                   | 2.56        | 2.57         |
| state mid-sample                   | 1.00        | 1.00         |
| shock at the LAST period           | **5.9e-15** | 2.57         |
| all shocks x1.5 in the LAST period | **5.9e-15** | 0.82         |
| state at the LAST period           | **5.9e-15** | 1.00         |

The contemporaneous `ghu %*% eps_t` term appears identically on both sides of
the adding-up comparison and CANCELS, so a shock error is visible only through
its propagated (t+1 onward) effect, damped by T -- and the final period is not
checked at all.

`historical_decomposition()` therefore also returns `$transition_residual`,
`$transition_residual_by_period`, `$transition_worst_period` and
`$transition_ok`: the direct question, does the smoother's own output satisfy
its own transition, `s_t = T s_{t-1} + R eps_t`, period by period. It has
neither blind spot and catches all six corruptions. **Look at it first** when a
decomposition will not add up -- it separates "the decomposition is wrong" from
"its INPUTS are incoherent", and the per-period vector localises the latter. On
the reported case it puts the entire error in period 1, which is what
identified the initial state as the culprit.

This also answers the report's second acceptance criterion directly: the
decomposition no longer returns a non-additive result as if it were valid --
it warns, says which of the two checks failed, and where.

Nothing about the computed contributions changes.

## Cross-loaded exact observations, pinned as coherent

Three observables each loading both states, none with measurement error, the
third an exact combination of the other two -- a predictable component dropped
in every period (160 of 160), and with ragged edges the dropped SET varying
period to period (130 of 160). Transition residual 6.7e-16 and 1.3e-15, with
the smoothed states exact against the simulated truth. Both are regression
tests now, as is the reporter's own four-variant reproducer.

# dynhr 0.9.3.5

**A successful Cholesky is not a test that the innovation covariance is
invertible, and every multivariate path was using it as one.**

Reported from a shock- and historical-decomposition matching exercise:
`me_variance = 0` gave adding-up residuals of ~2e6 and `me_variance = 1e-12`
about 4. The decomposition was reporting the problem faithfully -- its
residual tracks the smoother's own transition residual, and the smoother's
states and shocks had stopped being consistent with each other.

**The defect.** `chol()` can factorise a matrix that is singular to round-off
and return a garbage pivot. On a stochastically singular system -- more
observables than shocks, or an observable that is an exact combination of
others -- that produces a finite, badly wrong answer instead of a detected
failure. Measured on a two-observable / one-shock fixture with
`me_variance = 0`: the innovation covariance had `rcond` 1.8e-17 and a
NEGATIVE determinant, `chol()` succeeded, and

```
kalman_filter()   +292.7      <- garbage, and too HIGH
kalman_smoother()  -10.8      <- dropped the component in 19 of 20 periods
univariate filter  -26.8      <- correct
```

The smoother's error sat entirely in the ONE period where `chol()` happened to
succeed: its smoothed states were exact while its smoothed shocks were
inconsistent with them by 0.18. A likelihood that is too high is the dangerous
direction -- an optimiser walks straight into it.

**The fix**, in the three places that inverted F: after a successful `chol()`,
test the pivots. The i-th squared diagonal of the Cholesky factor IS that
observable's conditional variance -- the same quantity the univariate filter
has always skipped on, and the one `.smoother_informative_obs()` already used
to pick the informative subset. It was simply never consulted when `chol()`
succeeded. Sites: `.kf_step()` (the `standard` and `dare` paths), the C++
`kalman_standard_loop_cpp` fast path, and the smoother's forward pass.

**After the fix, every path agrees.** All five filter methods and the smoother
return -26.828497 on that fixture, each multivariate path detecting the
singularity and rerouting to the univariate filter, which drops the
uninformative component instead of inverting through it. The decomposition
residual goes from 2.9e-1 to 8.9e-16 at `me_variance = 0`, and from 3.4e-4 to
2.7e-15 at 1e-12.

**Where the boundary sits.** A component whose conditional variance is below
`kalman_tol` (default 1e-10, relative to F's scale) is dropped and the answer
is exact. Above it the component genuinely carries information and the
accuracy is the conditioning limit of inverting F, about `eps / me_variance` --
a numerical fact rather than a defect, and one that
`historical_decomposition()`'s adding-up check now surfaces rather than
leaving to be discovered downstream.

**One test changed meaning.** `test-kalman-smoother-na.R`'s "loglik decreases
when data are removed" used a fixture that is an AR(1) plus an ALIAS -- two
observables, one shock -- and blanked one of them. It passed only because the
smoother was inverting through the singular F, so the "information" being
removed was round-off. With an exact alias, removing either series alone costs
nothing, because the other still pins the state; only removing both does. The
test now asserts that, which is sharper than what it replaced: the two
one-column runs agree with the full run to the bit and with each other.

# dynhr 0.9.3.4

A fourth report against the filtering surface, on historical decomposition.
The decomposition function was correct; its integration contract was neither
documented nor checked, and one input combination was genuinely wrong.

**Bug fix: `pre_sample` put the decomposition out of phase.** With a backfill,
`kalman_smoother()` returns its series trimmed to the caller's sample, but
`smoothed_initial` is still `s_{0|T}` for the PADDED sample -- the state `k`
periods earlier. `smoothed_initial_state()` handed that to
`historical_decomposition()` as the anchor, so the initial-condition
trajectory ran `k` periods out of phase and the components stopped adding up
(measured 3.8e-2 on the two-shock fixture at `k = 2`), with every dimension
still correct. It now returns the last pre-sample row, which is the period
immediately before the returned rows.

**The contract is now enforced, not just described.** A correctly-SIZED but
wrongly-ORDERED `s0`, or shock columns in a different order from
`ss$shock_names`, used to be accepted in silence. Since the adding-up residual
is exactly the size of the initial-condition error, on a model with large
states -- or a kappa-initialised smoother, whose `s_{0|T}` carries the
arbitrary prior -- that silence surfaces as a residual of 1e8-1e9 with nothing
to say why. Now: columns and named vectors are matched to the state space by
NAME and reordered, a mismatch is an error naming both sets, and a transposed
matrix is refused with the orientation it wanted.

**`adding_up_residual` is always present.** It used to appear only when
`smoothed_states` was supplied -- so the one call shape that can be silently
wrong was also the one with no diagnostic at all. It is now a number when
there is something to check against and `NA` with `$adding_up_note` when there
is not, alongside `$adding_up_relative`, `$adding_up_ok` and a `tol` argument
(default 1e-8, relative to the path's scale). Exceeding it warns and names the
likely causes.

**Documented contract**, in `?historical_decomposition`: orientation (rows are
periods; `kalman_filter()`'s state matrices are the other way round and need
`t()`), the state at the first contribution period (row 1 loads the PRE-SAMPLE
state, which is why the `initial` column exists), pre- versus post-transition
(pre-: the state entering the period plus that period's shock), and how
`shock_means` / `known_shocks` enter (through the smoother's output, in their
own shock's column, needing nothing here).

**Verified against IRIS `simulate(..., 'contributions', true)`** to 2e-16 on a
two-shock linear model: each isolated-shock column, the initial-condition
column, and the total. Two IRIS columns have no dynhr counterpart by
construction, and the docs now say so: a measurement-shock column (dynhr's
`me_variance` is observation noise, not a structural shock, so compare against
IRIS's structural columns plus its init column rather than its grand total),
and a nonlinear column that is identically zero for a linear model.
`$has_nonlinear_column` and `$has_residual_column` report this.

# dynhr 0.9.3.3

A third report against the filtering surface, this one a feature request
rather than a defect: an IRIS-compatible deterministic shock-mean path, and an
explicit state-timing contract on the result.

**Verified against IRIS itself.** With `shock_means` and the default
`shock_timing = "dated"`, dynhr reproduces IRIS's
`filter(m, d, range, 'vary=', j)` to machine precision -- on all three of its
outputs and on the smoothed shocks -- with **no adapter shift of any kind**
(checked against IRIS Toolbox Release 20180308 under Octave; the reference
values are transcribed into `test-shock-means.R` at full precision, since IRIS
is not a dependency). The names line up one-for-one: IRIS `'predict'` /
`'filter'` / `'smooth'` are dynhr `predicted_states` / `updated_states` /
`smoothed_states`.

If you previously needed a one-period shift to line the two up, the cause was
the MECHANISM, not the timing. `known_shocks` is an exact observation of the
shock; an IRIS `vary` tune sets the mean and leaves the shock random, so its
smoothed shock is *revised away* from the injected value -- 1.3516 against an
injected 1.5 on the fixture -- and no shift of an exactly-pinned path can
reproduce that. `shock_means` is the matching statement. IRIS says as much in
its own source: "The std dev of the tuned shocks remain unchanged and hence
the filtered shocks can differ from its tunes".


A third report against the filtering surface, this one a feature request
rather than a defect: an IRIS-compatible deterministic shock-mean path, and an
explicit state-timing contract on the result.

## Filtering

- **`shock_means` on `kalman_filter()` and `kalman_smoother()`: deterministic
  shock MEANS.** An
  `n_exo x T` matrix of mean shifts (`NA` and `0` both mean "no shift here", so
  the matrix shape `known_shocks` uses works unchanged), with rows matched by
  name.

  **This is a different statement from `known_shocks`, and the difference is
  the point of having both.** `known_shocks` says the REALISATION is known,
  `eps = v`: the shock stops being random, its variance is used up, and where
  it has a prior density the value enters the likelihood. `shock_means` says
  the MEAN is known and the shock **keeps its variance**: nothing is observed,
  nothing about `m` is estimated, and the likelihood gains no term. It is a
  deterministic input, entering the transition and measurement constants as
  `R m_t` and `D m_t` before the update.

  The two coincide exactly for a shock with no prior variance -- knowing the
  mean and knowing the realisation are then the same statement -- and the
  tests pin both halves: identical states and log-likelihood there, materially
  different where the shock still has variance. Agreement alone would prove
  nothing, since an argument that was quietly ignored would also produce it.

  No method routing is involved. The mean path splits off as a deterministic
  trajectory subtracted from the data and added back to the reported states,
  so every `method` evaluates it identically (`known_shocks`, by contrast, can
  only be expressed on the augmented state and forces `method = "univariate"`),
  and zero-variance shocks, missing observations, `a0`/`P0` and both diffuse
  recursions are non-events.

- **`shock_timing`** selects how the columns are read: `"dated"` (default)
  makes column `t` the shock dated `t`, entering `s_t` and `y_t` -- the dating
  `known_shocks` and `shock_scale` already use -- while
  `"transition_next"` reads column `t` as driving the transition OUT of period
  `t`, so it lands on `s_{t+1}`. The second is the one-period adapter shift a
  caller comparing against a package with the other convention would otherwise
  apply by hand; it is now stated in the call.

- **`updated_states` and `predicted_states`**, with the timing in the names:
  `updated_states[, t]` is `s_{t|t}` and `predicted_states[, t]` is `s_{t|t-1}`
  (with `s_{1|0}` from `a0`). `filtered_states` is the same matrix as
  `updated_states`, kept under the name the rest of the package uses. The
  prediction needs nothing extra from the recursions: `E[eps_t] = 0` in the
  deviation system makes `s_{t|t-1} = T s_{t-1|t-1}` exact.

  A unit pulse in `shock_means` with no data to update on therefore traces the
  model's own impulse response exactly -- `updated_states[, t] = T^(t-1) R e_j`
  under `"dated"`, with the two paths coinciding because there is nothing to
  update with. That identity is what a cross-package shock-response comparison
  reduces to once the timing is fixed, and it is asserted to 1e-12 without
  needing the other package present.

- **One deliberate asymmetry between the two entry points.** A known shock's
  REALISATION is fixed, so `kalman_smoother()` reports the injected value back
  as itself; a MEAN leaves the shock random, so `smoothed_shocks` reports
  `m_t + u_{t|T}` -- the mean plus the smoothed deviation around it. A smoother
  that returned the mean unrevised would be ignoring the data; one that
  ignored the mean would be ignoring the input. The two coincide on a
  zero-variance shock, where there is no deviation left to revise.

- `kalman_smoother()` also gains `updated_states` and `predicted_states` under
  the same names, the latter being the mean counterpart of `predicted_cov` --
  computed all along and simply not returned. Note the orientation differs
  from the filter's by long-standing convention: rows are periods here.

- `$diagnostics$shock_means` reports the timing, the number of shifted cells
  and which shocks they belong to, on both entry points. `loglik_type` stays
  `"marginal"`: an input conditions nothing.

# dynhr 0.9.3.2

The follow-up to 0.9.3.1, from a second report against the same surface. Four
behavioural gaps, one of which was a silent wrong answer.

**Two existing calls return different numbers, both deliberately.**

1. `kalman_filter(known_shocks = ...)` on a shock with **zero prior variance**
   used to ignore the injected value completely. The injected row's forecast
   variance is exactly zero, so the sequential filter's `F > kalman_tol` guard
   skipped it and the shock never reached the transition -- no warning, and a
   filtered path that stayed on its unforced trajectory. On the one-state
   fixture in `test-statespace-semantics.R` the filtered state was
   `0 0 0 0 0 0` where the answer is `0 1 0.8 0.64 0.512 0.4096`, and the
   log-likelihood was -10.257 instead of -13.741. If you have conditioned on a
   deterministic shock (`stderr 0`, or a `shock_scale` of zero at that period),
   re-run it.

2. `kalman_smoother()` on a **unit-root model with `lik_init = "auto"`** now
   runs the exact diffuse recursion instead of substituting `P0 = 1e6 * I`.
   The smoothed states move by ~1e-8 (the fallback was accurate); the
   log-likelihood moves by a lot, because the kappa one carried an arbitrary
   additive constant and could not be compared with anything. It now equals
   `kalman_filter(lik_init = "diffuse", method = "univariate")` to machine
   precision. `lik_init = "kappa"` still asks for the old prior explicitly.

Everything else in this release is additive.

## Filtering and smoothing

- **Deterministic known shocks.** A shock with no prior variance is
  deterministic BEFORE conditioning, so an injected value is an INPUT, not an
  observation carrying no information. It is now applied as a mean shift with
  no covariance update -- which is exactly what conditioning on a degenerate
  component is, and exactly the limit of the ordinary update as the variance
  vanishes. A point mass carries no density, so the joint and the conditional
  likelihood coincide there and the filter agrees with the smoother to machine
  precision rather than up to a `log p(eps = v)` term.

- **`kalman_filter()` reports the state hand-off: `final_state`,
  `final_cov`.** `s_{T|T}` and `Var(s_T | y_{1:T})`, named, in the same
  convention `a0` / `P0` take -- so a split sample can be filtered in two
  calls and the prediction-error decomposition holds:

  ```r
  f1 <- kalman_filter(y[1:k, ],     dr, m, p, obs_vars = obs)
  f2 <- kalman_filter(y[(k+1):T, ], dr, m, p, obs_vars = obs,
                      a0 = f1$final_state, P0 = f1$final_cov)
  f1$loglik + f2$loglik           # == the unsplit loglik
  ```

  Pinned to 1e-9 for `standard`, `univariate` and `dare` in
  `test-statespace-semantics.R`. `final_cov` is `NULL` for
  `method = "chandrasekhar"`, which propagates low-rank increments precisely
  so that P is never formed -- reporting NULL is the honest answer, and "auto"
  does not choose that method below `n_state = 100`.

  Before this, `a0` and `P0` existed but there was no way to get a covariance
  out of the filter to put into them.

- **Exact diffuse smoothing (`lik_init = "diffuse"` on
  `kalman_smoother()`).** A new sequential (Koopman-Durbin) smoother on the
  augmented state `[s_{t-1}; eps_t]` -- `R/smoother-diffuse.R`. The augmented
  form is what makes it tractable: dynhr's shocks enter the observation
  equation directly, and the textbook diffuse recursions assume uncorrelated
  noise, which the augmentation restores. It also makes the disturbance
  smoother free -- `E[eps_t | y]` is a block of the smoothed state, so there
  is no second recursion to keep consistent with the first.

  The recursions are derived in the file as the `kappa -> Inf` expansion of
  the ordinary ones, not quoted, and every `L` is a rank-one update, so the
  `N` recursion costs O(nb^2) per observable rather than O(nb^3).

  Certified against `.gls_smoother()` with zero prior precision on the diffuse
  coordinates -- a projection, no recursion, no kappa: **3e-15** on states,
  shocks and the pre-sample initial state, where the old fallback managed
  5e-8. The likelihood matches the filter's exact diffuse value to 12 digits.
  Missing observations do not downgrade it: "missing", "diffuse" and "exactly
  predictable" are three independent decisions per (period, observable) in a
  sequential filter.

- **Structured run diagnostics: `$diagnostics` on both entry points.** Same
  field names on the filter and the smoother, so a parity harness reads either
  without a special case: `method_requested` / `method_used`,
  `lik_init_requested` / `lik_init_used`, `routing` (a data frame of
  `from` / `to` / `reason`, one row per automatic reroute or fallback),
  `diffuse_periods`, `missing_by_period` / `n_missing`, `dropped_by_period` /
  `n_dropped`, `known_shocks` and `loglik_type`.

  `missing` and `dropped` are counted separately on purpose: "not observed"
  and "observed but exactly predictable" shorten the conditioning set for
  different reasons, and conflating them hides a stochastic singularity.
  `loglik_type` is `"marginal"`, `"joint"` (an injected shock with a prior
  density) or `"conditional"` (all injections deterministic; also the
  smoother's convention).

- **`pre_sample` now pads `known_shocks` too.** The backfill rows are
  prepended before the known-shock matrix is validated, so an `n_exo x T`
  matrix on the caller's sample used to be rejected as the wrong width. It is
  padded with `NA` (= unknown) instead. The pre-sample series are also read
  after the deterministic trajectory is added back, so a backfill run together
  with `known_shocks` reports the padded periods with the injected shocks in
  them.

## Testing

- `test-statespace-semantics.R` (new): the four requests plus the requested
  A-F regression matrix over (missing data, zero covariance, unit root, known
  shocks). The matrix asserts **superposition** rather than a pinned number --
  injecting a deterministic shock must equal subtracting its trajectory from
  the data and injecting nothing, which is linearity of the state space and
  needs no oracle.

- `test-kalman-smoother-exact.R` block (d) is inverted: it recorded that the
  smoother had NO exact-diffuse initialisation. It now pins that the exact
  answer IS the `kappa -> Inf` limit -- along the RIGHT approximating sequence
  (`P_star + kappa * P_inf`, diffuse only in the unit-root direction). The old
  `kappa * I` fallback is diffuse in every direction, a different prior whose
  limit is 9.3e-4 away and flat in kappa.

# dynhr 0.9.3.1

A filtering and smoothing release. No export is added or removed and every
addition is an optional argument on a function that already existed, so no
existing call signature changes.

**One existing call DOES return a different number, deliberately.**
`lik_init = "diffuse"` with missing observations used to warn and silently
downgrade to `"kappa"`: -44.097 where the exact diffuse log-likelihood is
-36.271 on the local-level fixture with three gaps. That is the fix in this
release rather than a side effect of it -- but if you have run a
diffuse-initialised filter on data with gaps, you were getting the kappa
answer, and you will now get a different (correct) one.

Every other new argument is inert when omitted, and that is asserted rather
than assumed: `a0 = 0`, `pre_sample = 0` and an all-`NA` `known_shocks` each
reproduce the 0.9.3 result at `tolerance = 0`.

The rest closes gaps that had been worked around by hand: no way to set the
initial state, no way to backfill the latent history before the sample, and no
way to tell the filter about a shock you already know.

## Filtering and smoothing

- **Known historical shocks can be injected: `known_shocks` on
  `kalman_filter()` and `kalman_smoother()`.** An `n_exo x T` matrix carrying
  the value where a shock is known and `NA` where it is not -- the
  `NA`-as-unknown convention `data` already uses, and the `n_exo x T` shape
  `shock_scale` already uses, with rows matched by name. For an announced
  policy change, a measured intervention, or a judgemental adjustment;
  `known_shocks_sd` (filter) makes the injection soft rather than exact.

  The two entry points reach it by different mechanisms, because their
  recursions differ. The filter's univariate path already runs on the
  augmented state `[s_{t-1}; eps_t]`, so the shocks ARE state components there
  and a known shock is an exact observation of one -- passing `known_shocks`
  routes to `method = "univariate"`. The smoother instead splits the known
  shock off as the deterministic part of the system, smooths the remainder,
  and adds it back. Both were checked against closed forms that use no
  recursion at all (`helper-gls-smoother.R`): the smoother matches a
  projection oracle to 2e-16, and the filter matches an analytic joint density
  to 1e-10.

  **Semantics, because it is a modelling choice rather than a detail:** the
  filter treats the known value as an observation and reports the JOINT
  `log p(y, eps = v)`, so a known shock informs that shock's own standard
  error. The smoother conditions. The two differ by exactly `log p(eps = v)`
  -- verified to 0.0e+00 -- so subtract that term if you want the conditional
  likelihood.


- **`kalman_filter()` and `kalman_smoother()` take `a0` and `P0`.** The
  initial state mean and covariance were hard-coded -- zero, and whatever
  `lik_init` implied -- so there was no way to carry a state across a sample
  split, condition on a known history, or hand the smoother a prior of your
  own. The internals had always taken them; only the public shape was missing.
  `a0` is in **deviations from the steady state** (the convention the states
  are reported in; `data` is in levels) and is matched BY NAME when named.
  `P0` must be symmetric and positive semi-definite, may be a scalar for a
  multiple of the identity, and is reported back as `lik_init = "user"`.

  The identity that pins this: splitting a sample at `k`, filtering the first
  block, and starting the second from `(filtered_states[k, ],
  filtered_cov[, , k])` reproduces the whole-sample log-likelihood -- verified
  to 5e-13 on `nk_demo`.

  Two traps were found while wiring it, both of the "accepted, validated,
  reported, then silently dropped" kind: `kalman_standard_loop_cpp()` takes no
  initial state (it starts at zero unconditionally), so the default fast path
  for a stationary model now falls through to the R loop whenever `a0` is
  non-zero; and the smoother's backward pass read a hard zero for
  \eqn{s_{0|0}}, which would have ignored `a0` in `smoothed_initial` while the
  forward pass honoured it.

- **`lik_init = "diffuse"` now works with missing observations.** It used to
  warn and downgrade to `"kappa"`, which answers a different question: on a
  local-level model with three gaps, -44.097 against the exact -36.271. The
  exact-diffuse recursion and missing data meet in the sequential
  (Koopman-Durbin univariate) filter, which skips a missing observable one at
  a time -- as it already did for every other initialisation -- so the
  capability was implemented and unreachable rather than absent. `method =
  "auto"` now routes the combination there; an explicit multivariate method
  reroutes with one notice rather than downgrading.

  The difference is not stylistic: the exact-diffuse filtered states are
  invariant to the diffuse scale to the last bit, while the `"kappa"` states
  keep moving as kappa grows (5.4e-5 from 1e4 to 1e6, 5.4e-7 from 1e6 to 1e8)
  -- an approximating sequence, not the answer. Complete-data results are
  unchanged, and on a single unit root the two diffuse conventions still agree
  to 5e-12.

- **`kalman_smoother(pre_sample = k)`** backfills `k` periods of latent
  history before the first observation and returns them in
  `presample_states` / `presample_shocks` / `presample_cov`, chronological,
  with every other series still aligned to `data`. No new recursion: an
  all-missing period is predict-only, so this is the ordinary backward pass
  over `k` padded rows -- the mechanism that has always produced the single
  `smoothed_initial` period, and the last backfilled row reproduces it
  exactly. The log-likelihood is unchanged and the in-sample states move by
  2e-14. `me_extra` and `shock_scale` are padded with it, so their per-period
  columns stay aligned.

- **How far the smoother's kappa fallback sits from the exact diffuse answer
  is now measured, not assumed.** `helper-gls-smoother.R` computes the
  smoother by projection instead of recursion -- the diffuse case is a zero on
  the prior precision, not a separate algorithm -- and is certified against
  the existing smoother on a stationary model to 3e-16 before being used
  where the answer is unknown. On a unit-root fixture the kappa fallback's
  smoothed states are within **5e-8 absolute on a state of scale 3.3**, about
  eight significant figures. The error falls as `1/kappa` and then RISES
  again as round-off takes over (5.0e-6, 5.0e-8, 4.8e-9, 2.8e-7 at kappa
  1e4/1e6/1e8/1e10), so raising `.DIFFUSE_SCALE` is not the fix -- but nor is
  the fallback the crude approximation it looked like. An exact diffuse
  SMOOTHER remains unimplemented; the tests now state its acceptance
  criterion. Where the diffuse treatment genuinely matters is the
  log-likelihood's unbounded kappa-dependent constant, and
  `kalman_filter(lik_init = "diffuse")` already handles that exactly.

- Multi-period pre-sample backfill needs no manual padding: the argument
  above wraps it. On `nk_demo` with four padded periods
  the in-sample states are unchanged to 2e-14, the log-likelihood is identical
  (missing rows contribute nothing), the transition identity holds across the
  pad boundary to 2e-16, and the last padded row reproduces
  `smoothed_initial`. Exact for a stationary model, since `P_0` is then the
  unconditional covariance; on a unit-root model it inherits the smoother's
  kappa fallback, which awaits an exact diffuse smoother.

# dynhr 0.9.3

This release adds four new estimation and solution capabilities — global
(projection) solutions with their own likelihood, Markov-switching DSGE
filtering and smoothing, mixed-frequency observation blocks, and
moment-based estimation.

It also carries a substantial correctness pass on the parts of the package
you reach for when fitting a model to data. Three of those are worth
singling out, because each was silent:

- `run_full_estimation()` built its sampler proposal from the **prior**
  rather than the posterior curvature — a structurally unreachable branch —
  which mixed poorly on most models and froze the chain outright on a
  well-identified one.
- `kalman_smoother()` never subtracted the **observation intercept**, so it
  silently required data in deviations and disagreed with `kalman_filter()`
  on the same series by tens of thousands of log points — and the historical
  decomposition behind diagnostics D11/D12 was running on exactly that
  mismatch.
- A parameter named **`sigma_e`** was deleted from the calibration before it
  was read, leaving it `NA` while compilation and solving continued.

Measurement error is now one noise model across every filter,
`kalman_smoother()` takes the same arguments as `kalman_filter()`, and every
filtering and smoothing entry point takes observables in **levels**.

The release carries breaking API changes; read the first section before
upgrading.

## Breaking changes

- **Three empty pass-through aliases are removed**, with no deprecation
  shims: use `solve_model()` for `solve_dsge()`, `run_mode_finding()` for
  `estimate_mode()`, and `run_posterior_estimation()` for
  `estimate_posterior()`.
- **Argument spellings are unified across every exported likelihood,
  filter and forecast function.** The data argument is `data` (was `Y`),
  the observable-name argument is `obs_vars` (was `obs_names` or
  `observables`), and the measurement-error argument is `me_variance` (was
  `me_var` or `me_sd`) — matching `make_posterior()` and
  `run_full_estimation()`. Matrix orientation is unchanged and is now
  stated per function.
- **Two of those renames change what you pass, not just the name.**
  `hank_loglik_ar()`, `hank_loglik_ar_grad()`,
  `hank_loglik_ar_structural_grad()` and `hank_simulate_aggregate()` now
  take a VARIANCE where they took a standard deviation — pass `me_sd^2`.
  `make_log_posterior_hank()`, `make_posterior_grad_hank_ar()` and
  `hank_ar_target()` took both `me_var` and an independently overridable
  `me_sd`; they now take the single `me_variance`, so their two likelihood
  branches always describe the same measurement error.
  `hank_loglik_ar_grad()$me` is still the score with respect to the
  standard deviation.
- **`run_posterior_estimation()`'s count arguments follow the
  sampler-level convention**: `nburn`/`ndraws`/`nchains`/`nparticles`/
  `nwalkers` are now `n_warmup`/`n_draws`/`n_chains`/`n_particles`/
  `n_walkers`, matching `run_full_estimation()` and `mcmc()`.
- **The inert `mh_scale` argument is removed** from
  `make_log_posterior_tpf()` and `dynhr_smc2()`'s `likelihood_args`, which
  now rejects unrecognised keys instead of silently dropping them. It had
  done nothing since the tempered particle filter's mutation step was
  corrected to hold the ancestor state fixed.
- **`dynhr_smc2()` returns a `dynhr_chains` object** like every other
  sampler entry point, so `print()`, `summary()` and `plot()` work on it.
  All previous fields are retained.
- **`kalman_smoother()`'s second argument is a decision-rule object, and
  every filtering and smoothing entry point takes observables in LEVELS.**
  Passing a pre-built `dsge_ss` — the original signature — is an error that
  names its replacement; the deprecated `ss = ` alias is gone with it. The
  data convention changed with the shape: a state space now carries its own
  observation intercept in the `d` field `new_dsge_ss()` has documented all
  along, so `kalman_smoother()`, `realtime_decomposition()`,
  `forecast_backtest()` and `conditional_forecast()`'s anchoring data are all
  in levels and the model's steady state is subtracted for you. If your
  series are already in deviations, pass `d = 0` to `kalman_smoother()`.
  Two entry points to the same recursion silently requiring different data
  was the defect; one release of polymorphic shim would have kept it alive
  in a second form, so it is settled here instead.

## New features

### Global (projection) solutions

- **`solve_global()` is usable as an estimation target.**
  `make_log_posterior(likelihood = "global_pf")` runs a bootstrap particle
  filter over a projection solution, so the model is never linearised. The
  estimate is unbiased for the marginal likelihood, making `pmmh()` over it
  a valid pseudo-marginal sampler. `global_pf_sbc()` is the matching
  rank-uniformity certification.
- **The default collocation domain is measured, not guessed.** It is
  derived from the model's own shock process and unconditional state
  dispersion, so an AR(1) with a large innovation gets a wider grid
  automatically; a four-fixture cover study set the default. `solve_global()`
  also validates its model class structurally and fails loud rather than
  silently returning a bad approximation.
- **`euler_errors()` is model-agnostic** — it reads the Euler equations from
  the parsed model instead of assuming RBC parameter names — and
  **`den_haan_marcet()`** is new. Both work on perturbation solutions too,
  so they can be used to decide whether a global solve is needed at all.
- `simulate()` on a `GlobalSolution` had an off-by-one in shock timing: the
  shock drawn in period `t` was applied to the wrong period. Fixed.

### Markov-switching DSGE

- **Filtering, smoothing and IRFs across regimes**: `ms_kim_filter()`,
  `ms_kim_smoother()`, `ms_kim_smoother_struct()` (structural switching) and
  `ms_irf()`.
- **GPB(3) collapse** (`collapse = "gpb3"`) keeps the pair
  `(s_{t-1}, s_t)` and collapses over two lags instead of one. It is
  strictly weaker as an approximation than GPB(2) and markedly more accurate
  when regimes are persistent: against an all-regime-path enumeration
  oracle, smoothed regime probabilities improved from 1.00 away to 2e-15 and
  states from 25 sd to 3e-15 sd on the GPB(2) breakdown draw. It is not a
  pointwise improvement — on 7 of 48 grid draws its error is up to 2.1x
  GPB(2)'s, at absolute levels below 3e-5 — and costs about 1.75x at
  `T = 300`, `h = 2`. `"gpb2"` remains the default and is bit-identical to
  before.
- The filter's covariance update is now Joseph-form, and a collapse
  diagnostic is available via `return_collapse_diag`.

### Mixed-frequency observations

- **`obs_aggregation`** declares an observable as the *k*-period temporal
  aggregate of a higher-frequency model variable, so a monthly model can be
  estimated on quarterly data without leaving the monthly frequency. Four
  aggregators: `flow_sum`, `flow_mean`, `stock_end` and `triangle`
  (Mariano–Murasawa). Implemented as fixed-weight state augmentation
  (Harvey 1989 §6.3), so `ZZ` and `TT` stay constant, the filter's hot path
  and both C++ kernels are untouched, and a model without `obs_aggregation`
  returns a byte-identical log-likelihood. `mf_augment_state_space()`,
  `mf_expand_observations()` and `mf_aggregation_weights()` expose the
  machinery directly.

### Estimation

- **`method_of_moments()`** — GMM and SMM by moment matching over the
  model-implied autocovariance structure, with identity / optimal /
  Newey–West / diagonal weighting, optional two-step, and analytic moment
  Jacobians where dynhr has them.
- **`forecast_backtest()`** — recursive expanding-window out-of-sample
  scoring with CRPS, log score, PIT and interval coverage; re-estimate at
  each origin or roll a single fit forward.
- **Delayed acceptance**: `mcmc(..., screen_fn = )` evaluates a cheap
  approximate likelihood first and only runs the expensive one on proposals
  that survive. The two-stage acceptance ratio keeps the exact posterior
  invariant.
- **Resumable chains**: `mcmc(checkpoint_dir =, resume =, flush_every =)`
  writes a checksummed, atomically-written state pack carrying the position,
  log-posterior, adaptation state and `.Random.seed`, so a resumed chain is
  statistically identical to the uninterrupted run. `mcmc_chain_state()`,
  `mcmc_chain_save()`, `mcmc_chain_restore()` and `mcmc_chain_extend()` are
  the sampler-agnostic primitives; restore refuses a tampered pack or one
  saved mid-adaptation.
- **`dynhr_model()` and the `dm_*()` verbs** (`dm_solve`, `dm_posterior`,
  `dm_mode`, `dm_sample`, `dm_diagnostics`, `dm_forecast`, `dm_irf`,
  `dm_test`) — one pipeline object carrying model, compiled, steady state,
  decision rules, data and priors, instead of threading six arguments
  through every call. Every vignette pipeline routed through the object
  returns bit-identical numbers to the functional path.

### Model input and output

- **`write_mod()`** serialises a parsed model back to Dynare `.mod` source.
  A `parse_mod()` -> `write_mod()` -> `parse_mod()` round trip is a cheap
  check that dynhr read a file the way you meant it.
- **`histval` blocks are parsed** into `model$histval` (lag-indexed), and
  **`smoother2histval()`** builds that history from a completed smoother run
  — the standard way to start a forecast or counterfactual from where the
  data left off.
- **`shock_groups` blocks are parsed** into `model$shock_groups`, consumed
  by the decomposition functions.

### Heterogeneous agents

- **`hank_ks_aggregate_risk()`** — Krusell–Smith with genuine aggregate
  risk, plus `hank_ks_risk_irf()` for generalised impulse responses and
  `hank_ks_ergodic_mean()` for the aggregate precautionary term without a
  full simulation. `hank_tfp_chain()` builds the aggregate productivity
  chain.

### Decompositions

- **Shock decompositions add up exactly.** `historical_decomposition()` and
  its OBC variant gain an `initial` column for the contribution of the
  initial state, which is not zero unless the sample starts at the steady
  state; the columns now reproduce the data to machine precision. Both take
  `shock_groups`. **`realtime_decomposition()`** re-runs the decomposition
  across data vintages, so a given quarter's story can be tracked as it was
  revised.

## Correctness

### Measurement error is one noise model everywhere

Measurement error was implemented inconsistently across filters: on several
paths it entered the forecast covariance only, acting as a regulariser
rather than as observation noise, which made those likelihoods disagree with
the exact Kalman filter by O(`me_variance`).

- The **multivariate Kalman filter**, the **Markov-switching filters** and
  the **SV Rao-Blackwellised particle filter** (`kf_step()` and its compiled
  kernel) now all implement the true i.i.d. law: `me_variance` enters the
  forecast covariance AND the Joseph state-covariance update.
- **`dynhr_sbc()`'s data-generating process now adds measurement error**, so
  the DGP and the likelihood describe the same model. An SBC on a filter
  with `me_variance > 0` against a noiseless DGP was certifying a
  mis-specification.

### Filtering and smoothing

- **`kalman_filter(method = "chandrasekhar")` is exact again**, and is
  re-admitted to `method = "auto"` above `n_state = 100`.
- **`kalman_smoother()`'s state pass is exact** under dynhr's timing
  convention.
- **`pkf_smoother_obc()`** uses the correlated-noise disturbance smoother.
- The Kim smoothers gained a joint-probability regime pass
  (`regime_pass = "joint"`).
- **`kalman_smoother()` rejected systems `kalman_filter()` handled.** Reported
  for a unit-root, `shock_scale`d model with a non-positive innovation
  covariance. The two functions treated a singular `F` by different
  mechanisms: `kalman_filter()` falls back to the univariate (Koopman-Durbin)
  filter, which skips any component whose conditional variance is below
  `kalman_tol` -- the correct treatment, since such a component is predictable
  exactly and carries no information -- while `kalman_smoother()` added JITTER
  on an absolute ladder (`1e-8` ... `1e-2`, then an unguarded
  `chol(F + 0.1 I)`). `F` is not O(1): a unit-root smoother starts from
  `P = 1e6 I` and `shock_scale` multiplies `Q` on top, so the ladder was
  either far too small -- and the unguarded rung threw, which is the reported
  rejection -- or it "worked" and silently corrupted the result. Switching one
  shock off via `shock_scale` (the `u_k = 0` idiom for forcing a series to its
  observed value, and what a hard `filter_tunes` tune does underneath) moved
  the smoother's log-likelihood to **-8.5e+09** where the filter returned
  **-63.7**. The smoother now makes the filter's decision: a zero-variance
  component is dropped for that period, exactly as it already treats a
  *missing* observable, and the update proceeds on the informative subset,
  whose `F` is positive definite by construction. The value becomes -70.0,
  and the smoother's offset from the filter is now the same constant whether
  or not a shock is switched off. Dropped components are reported, naming how
  many periods and components and stating that the log-likelihood is not
  comparable with an undropped run.
- **`kalman_smoother()` silently required data in DEVIATIONS.**
  `kalman_filter()` takes raw level data and subtracts the observation
  intercept `d = dr$ys[obs_vars]`; the smoother never did, and
  `build_dsge_state_space()` carries no steady state at all, so the
  requirement was unstated and its violation silent. On any model whose
  observables have non-zero steady states -- which is most of them --
  filtering and smoothing the same series disagreed wildly: on the bundled
  `nk_demo` (observable steady states 0.5, 2 and 4) the filter returned
  -757.6 and the smoother -33990.3. Hand-demeaning the data closes the gap to
  2e-13, confirming the intercept was the whole of it. The state space now
  carries the intercept (`build_dsge_state_space()$d`, rescaled by `sum(w)`
  under a mixed-frequency aggregator exactly as `kalman_filter()` rescales
  its own), every entry point subtracts it, and `d = 0` is the explicit
  escape hatch for data already in deviations.
- **The historical decomposition behind D11/D12 was running on levels
  through the deviations-only path.** `run_all_diagnostics(posterior)` built
  a state space, handed it the raw observable columns and smoothed them, so
  on any model with non-zero observable steady states the smoothed shocks
  absorbed the level offset: on `nk_demo`, `max |eps|` 4.83 against 0.80 —
  six times too large — and a log-likelihood of -33990.3 against -757.6. D12
  reports the mean of each smoothed shock and passes it at `|mean| < 0.1`, so
  the diagnostic was reporting the bug as a model failure. **Any historical
  decomposition or smoothed-shock series produced through
  `run_all_diagnostics()` before this release, on a model whose observables
  have non-zero steady states, should be recomputed.** The same latent defect
  is closed in `realtime_decomposition()`, `forecast_backtest()` (whose
  predictive mean now carries the intercept back so it is scored against the
  realised level) and `conditional_forecast()`, whose Gaussian branch alone
  took deviations while its `tpf` and `pskf` branches demeaned for
  themselves — the same call needed different data depending on
  `ctx$likelihood`.

- **`kalman_smoother()` takes the same arguments as `kalman_filter()`.** Every
  other filter/smoother entry point -- `kalman_filter()`, `ms_kim_filter()`,
  `ms_kim_smoother()`, `kf_innovation_diagnostics()` -- takes
  `(data, dr, model, params, obs_vars, me_variance, ...)`. The Gaussian
  smoother took `(data, ss, Q, me_extra, shock_scale)`, so the obvious call by
  analogy after filtering failed and the only signpost to the required
  `build_dsge_state_space()` step was a single `@param` line; `?kalman_smoother`
  had no example and no `\\seealso`. It now accepts the filter's arguments,
  including **`me_variance`**, which it previously lacked entirely -- so a
  model filtered with measurement error could not be smoothed under the same
  noise model. Filter and smoother now agree to 1e-8 across `me_variance`
  values. The pre-built state space is no longer accepted: see the breaking
  changes above for the one-line migration.
- **`kalman_smoother(lik_init = )`.** `kalman_filter()` refuses
  `lik_init = "auto"` together with `shock_scale` on a nonstationary model and
  instructs the caller to pass `"kappa"` or `"stationary"` explicitly -- which
  the smoother had no way to accept, so the two could not be made comparable
  even in principle. `"stationary"` now errors on a unit root instead of
  silently returning a kappa-initialised answer. There is still no exact
  *diffuse* initialisation in the smoother.
- **The same defect is fixed in `conditional_forecast()`**, in four places.
  Its internal Kalman pass added a `1e-10` ridge on every period (so it was
  never an unregularised filter) and fell back to an unguarded
  `chol(F + 1e-6 I)`; and its three condition-system solves used unguarded
  `chol(M + 1e-12 I)`. Those Gram matrices lose rank exactly when the
  conditions are collinear, over-specified, or routed through a switched-off
  shock -- normal things to ask for. All four now drop uninformative
  components or take the minimum-norm solution via a relative-cutoff
  pseudo-inverse, and warn rather than throwing.

### Parsing

- **A parameter named `sigma_e` was silently dropped** (reported against
  0.9.1). `remove_blocks()` strips Dynare *command* statements before the
  calibration is read, and its keyword list ends with `Sigma_e` — Dynare's
  shock-covariance assignment — but the loop matched case-INSENSITIVELY. A
  user's `sigma_e = 1;` was therefore deleted as if it were that command, and
  the parameter survived declared but `NA`, after which `compile_model()` and
  the solvers happily proceeded on an invalid calibration behind a warning.
  Dynare identifiers are case-sensitive and its command is spelled with a
  capital S, so lower-case `sigma_e` is an ordinary parameter name; it is now
  matched case-sensitively. `Sigma_e` was the only entry in that list whose
  lower-case form is a legal user identifier — the rest (`stoch_simul`,
  `steady`, `check`, …) are genuine Dynare reserved words and stay
  case-insensitive.

### Sampling

- **The one-call estimation API proposed from the PRIOR, not the posterior.**
  `run_full_estimation()` and the estimation runner both built the RWMH
  proposal as `if (!is.null(mode_res$V_mode)) ... else diag(prior_spec$std^2)`,
  but the mode result on that path comes from the optimiser core
  `.run_mode_finding()`, which returns only the mode and its convergence
  record and never sets `V_mode`. The condition was structurally unreachable,
  so every proposal was a prior-variance diagonal and the posterior curvature
  the mode-finder had just located was silently discarded. On a model whose
  prior and posterior sit at a similar scale this merely mixed poorly
  (`fs2000`: 10.2% acceptance); on a well-identified one it froze the chain
  outright — 0% acceptance and exactly zero posterior variance, every draw
  equal to the mode. The Hessian-based proposal logic
  (`.make_pd` regularisation plus eigen-basis capping at the prior scale) is
  now shared with `run_mode_finding()` rather than duplicated, and the
  samplers compute the Hessian at the mode and use it. Measured after the
  fix: the nine-parameter `nk_demo` fixture goes 0% -> 34.3% acceptance with
  every posterior mean within one standard deviation of the values the data
  were simulated at, and `fs2000` goes 10.2% -> 29.4%. If the Hessian genuinely
  cannot be evaluated the proposal still degrades to the prior diagonal, but
  now **warns**: a silent version of that fallback is what hid this for three
  releases. `test-sampler-proposal.R` pins the property no test had asserted —
  that the chain moves at all.

### Numerics

- **The SV RB-PF agrees between R and C++ past the volatility overflow
  point.** The two diverged without bound once a volatility particle left
  the useful double range. The overflow was not the cause: the forecast
  covariance and its Cholesky factor are bit-identical in both, but
  `chol2inv()` and `arma::inv_sympd()` differ by one ulp in `F^-1`, and when
  the observation block is perfectly informative the exact Kalman gain is
  the identity, so the Joseph factors are exactly zero. That one ulp made
  them entirely rounding noise, the state covariance rounding noise squared,
  and the next period's inverse amplified it without limit. Both kernels now
  snap a Joseph entry lying within a few ulps of the magnitudes that
  cancelled to the exact zero it approximates; `kf_step()` also rejects a
  non-finite forecast covariance before `chol()`, as the compiled kernel
  already did. Verified over 864 parameter/seed/length combinations:
  0 divergent, worst relative gap 4.1e-16.
- **Order-3 cumulant moments were projected onto only the `(i,i,k)` slice**
  of the third-moment tensor. Fixed, along with three further cumulant/GMM
  defects (the analytic GMM weight matrix's lag handling among them).
- **Seeded particle-filter closures no longer reset the caller's RNG
  stream** — a seeded closure used to reseed the global stream on every
  evaluation, silently correlating an outer sampler's own draws.
- **`power` (power-posterior tempering) reaches every likelihood, exactly
  once**, and is shipped to parallel workers.
- **Stationary initial covariances (PSKF, TPF) come from the real Lyapunov
  solution** rather than a truncated series.
- **Order-2 shortcuts no longer discard shock correlations**: `.linear_dr2`
  dropped the off-diagonal of `Sigma_e`.
- **`kalman_filter()`'s singularity fallback is conditional and audible**
  instead of silent, and `.safe_inv()` truncates on a relative singular-value
  cutoff rather than an absolute one.
- **`hank_het_block()` fails loud on reducible income chains** and gains
  `dist_init` — at `p_un = p_nu = 0` the non-participation state is a closed
  class, and a uniform-seeded power iteration stranded about a third of the
  mass there.
- **`dynhr_set_options()` values now reach mirai daemons**, and parallel
  workers receive the host's option state.
- The cumulant gradient no longer returns a silent `NaN` for an explosive or
  non-stationary draw.

### API and structure

- Two exported functions failed on every call and had no test:
  `ramsey_obc_pwlinear()` built its shock sequence transposed, and
  `diag_prior_sensitivity()` referenced an undeclared argument. Both fixed,
  with smoke tests added for every previously untested export.
- Five `print()`/`summary()` methods were written but never registered, so
  from an installed package they fell through to `print.default()` and
  dumped the whole object. All registered, with a structural guard against
  recurrence.
- `run_posterior_estimation()` with `n_chains >= 2` crashed on any
  one-parameter model. Fixed.
- The HANK result classes share one compact `print.hank_block()`, so
  printing a block no longer dumps the stationary distribution. Nineteen
  internal-but-exported oracles are marked `@keywords internal`.
- One discrete-Lyapunov solver (`solve_lyapunov()`); the removed
  direct-Kronecker variant returned a false `NaN` on highly non-normal
  stable matrices and a negative variance for a scalar explosive root.
- `ast_to_string()` under-parenthesised `a - (b - c)`.

# dynhr 0.9.2

The headline addition is **SMC² (`dynhr_smc2()`)**: sequential Monte Carlo
over parameters driven by the package's particle-filter likelihoods — the
order-3 tempered particle filter and the stochastic-volatility
Rao-Blackwellized filter, both SBC-certified in this cycle. Supporting it,
the tempered particle filter received a substantial correctness and
performance pass, and the HANK household line gains three-state (E/U/N)
labour status, asset-indexed transfer incidence, and a complete
distribution-Jacobian family.

dynhr now requires **R >= 4.3.0** (declared via `Depends` and
`SystemRequirements: C++20`): the compiled two/three-asset solvers use
`std::barrier`, which needs GCC >= 11 — Rtools43 on Windows. This
formalizes the 0.9.1.1 hotfix (`CXX_STD = CXX20` in both Makevars, so
Windows toolchains no longer fall back to C++17 and fail at
`#include <barrier>`).

## New features

- **`dynhr_smc2()` — SMC² over parameters.** Likelihood-tempered SMC on
  theta where each evaluation is an unbiased particle-filter estimate
  (`likelihood = "tpf"` or `"sv_rbpf"`), with pseudo-marginal mutation
  moves: every theta-particle carries its stored likelihood estimate, the
  incumbent is never re-evaluated, and fresh filter randomness is drawn
  only at proposals — so the final-stage marginal is the exact posterior
  and the telescoped evidence estimate is valid. Filter-level tuning
  routes through `likelihood_args =`; `parallel = TRUE` evaluates
  theta-particles on a mirai daemon pool. Oracle-tested against exact-KF
  SMC on degenerate-volatility and linear-observation models.
- **Order-3 tempered particle filter** (`make_log_posterior_tpf(order = 3)`):
  the particle transition implements the pruned third-order recursion
  exactly, verified against the reference simulator and the order-3
  Gaussian pruned-KF, and **certified calibrated by a full R = 100
  rank-uniformity SBC** (`tpf_order3_sbc()`, the new full-certification
  tier under `DYNHR_SBC_FULL`, joining `sv_rbpf_sbc()`).
- **Three-state (E/U/N) household labour status**
  (`hank_employment_income3()`): all six transition rates are
  Jacobian-ready inputs (direct N<->E flows default to zero — a
  documented, overridable restriction), non-employed income is
  calibratable, and the three stocks plus six gross flows are selectable,
  ND-verified outputs. With the participation margins zeroed the
  two-state model is recovered (chain, income and shares bit-identically;
  the jointly-solved policies to machine precision).
- **Asset-indexed transfer incidence**: `hank_het_block(Tr_incidence =)`
  accepts an `n_e x n_a` matrix `omega(e, a)` over beginning-of-period
  states (the transfer stays lump-sum), enabling wealth-correlated
  incidence schedules. The block reports `Omega_ss` and a new `"Omega"`
  Jacobian output (the fiscal outlay aggregate); cash-on-hand positivity
  is asserted at the steady state. Vector and default incidence are
  bit-identical to 0.9.1.
- **Distribution-Jacobian family completed**: new
  `hank_het2_dist_jacobian()` / `hank_het3_dist_jacobian()` (+ `_nd`
  numerical oracles) give the full cross-sectional distribution response
  via the fake-news algorithm, including the `theta_coll` (collateral
  LTV) column; `hank_model_dist_irf()` gains the `"het2"`/`"het3"` arms.
  The three-asset backward sweep is now a single compiled call with a
  persistent worker pool, and the expectation stream is shared across
  inputs — both bit-identical to the R paths they replace.
- `hank_model()` fails loud on structural GE singularity, naming the
  offending unknown/target instead of an opaque LAPACK error.

## Tempered particle filter: correctness and performance

All four items below change `likelihood = "tpf"` values relative to 0.9.1;
fixed-seed pins were re-computed.

- **Observation-equation fidelity fix.** The TPF previously used a
  linearized observation equation, dropping the quadratic/cubic
  observation-row tensors the pruned model's data actually carry — on
  models with nonlinear observables it evaluated the wrong model's
  likelihood (found by SBC certification as a decisively miscalibrated
  posterior). The full observation reconstruction is now used at both
  orders and in both backends. Found in passing: the compiled order-2
  kernel had a Kronecker index-layout bug in its `ghxu` term, active only
  for multi-shock models with `n_e != n_s`; R and C++ backends now agree
  on a full-filter multi-shock run and a permanent parity test guards it.
- **Mutation-kernel unbiasedness fix.** The RWMH mutation step previously
  random-walked the ancestor state with a likelihood-only acceptance
  ratio — a kernel that biased the likelihood estimator upward,
  increasingly with `n_mh`. Mutation now refreshes only the period-t
  shock with the ancestor state fixed (Herbst-Schorfheide 2019); the
  estimator is verified unbiased against the exact true-ME likelihood
  and the order-3 SBC certification was re-run under the corrected
  kernel. `mh_scale` is now inert (kept for API compatibility).
- **Missing-data support.** Fully-missing periods now propagate the
  particle cloud (previously the period was skipped entirely, gluing the
  gap's endpoints together — the likelihood of a different model);
  partially-missing periods evaluate the observed elements only, matching
  the univariate Kalman filter's per-element NA convention. Fully
  observed datasets are unchanged.
- **Stationary particle initialization by default** (`burn_in_init = 50`):
  the cloud starts from (approximately) the full pruned stationary joint
  instead of zeroed higher-order layers, removing a short-sample
  persistence bias identified by SBC. Pass `burn_in_init = 0` to
  reproduce 0.9.1 behavior (required with correlated-pseudo-marginal
  `U_list` evaluation, which has no burn-in slots).
- **~8x faster evaluation** at estimation-scale settings: the per-particle
  `kronecker()` loops were replaced with bit-identical row-indexed
  products and per-period invariants hoisted out of the mutation loop.
  A full order-3 SBC certification now takes ~40 minutes, down from ~10
  hours.

## Fixes

- Direct `kalman_filter()` calls with `me_variance > 0` no longer re-run a
  150-step R-level Riccati detector on every call (a ~13x per-call
  overhead; estimation closures were unaffected) — the result is memoized
  and the iteration exits on convergence.
- `sv_rbpf` no longer zero-weights a particle when `inv_sympd()` fails at
  the PD boundary despite a successful Cholesky (the inverse is recovered
  from the factor), and a period where every particle fails now warns
  instead of silently returning `-Inf`.
- `run_full_estimation(likelihood = "tpf")` no longer errors when
  `tpf_options` carries orchestration-level keys (e.g. `cpm_rho_u`).
- `sbc_uniformity_test()` takes the true rank support `L` explicitly (the
  observed maximum can under-count bins) and its tail-asymmetry statistic
  is exactly centered under the null for unequal edge bins.
- Rank-histogram plotting, preflight warning messages, and several Rd/
  documentation defects (including two README function names that did not
  exist) were corrected.

# dynhr 0.9.1

The headline addition is a first-class **stochastic-volatility (SV) estimation
layer** — latent, per-shock time-varying volatility on the shocks of a linear
DSGE, estimated end to end. No named peer (Dynare, gEcon, MacroModelling.jl)
ships a first-class SV declaration/estimation path. The other changes in this
cycle are package hygiene, bug fixes, and incremental improvements.

## New features

- **Stochastic-volatility-on-shocks estimation.** Declare independent AR(1)
  log-variance processes on any subset of a model's shocks
  (`stochastic_volatility;` mod-file block; `stochastic_volatility()` /
  `sv_entry()` constructors) and estimate their hyperparameters jointly with
  the structural parameters. Conditional on a volatility path the model is
  exactly linear-Gaussian with `shock_scale = exp(h_t/2)`, so estimation uses
  a **Rao-Blackwellized particle filter** — particle-filtering only the
  low-dimensional log-variance states and integrating the DSGE states
  analytically via the Kalman recursion — giving an unbiased marginal
  likelihood that plugs directly into `pmmh()`. New API surface:
  `likelihood = "sv_rbpf"` in `make_log_posterior()`, `run_mode_finding()`,
  and `run_full_estimation()`; a `stochastic_volatility=` argument on the
  runners; and the exported single-step Kalman primitives `kf_step()` /
  `kf_stationary_init()`. Scope (v1): stationary models, order-1 (linear)
  solutions only — structural SV (volatility perceived in the decision
  rules) remains available via the order-2 pruned/TPF likelihoods.
- **Cross-machine benchmark.** `dynhr_benchmark()` runs a fixed, realistic
  estimation workload — Smets & Wouters (2007), 36 estimated parameters, 7
  observables, 160 quarters — through RWMH at a sweep of core counts and
  reports per-draw/per-second throughputs (`draws_per_sec`, `us_per_draw`,
  `speedup`, `efficiency`) alongside the system detail
  (`dynhr_system_info()`: hardware, OS, R build, BLAS/LAPACK, dynhr build)
  needed to compare machines. The workload is fingerprinted by its log
  posterior at the published mode; two results are comparable only if the
  fingerprint and the BLAS agree. The model, data, published mode and mode
  Hessian ship in `inst/extdata/models` (AER replication deposit
  openicpsr-116269-V1, BSD-3-Clause/CC BY 4.0; see `sw2007_SOURCE.md`).

## Fixes

- **Gradient-based estimation of pruned-state-space models.**
  `run_posterior_estimation()` / `find_mode()` now use the analytic
  (order-2) / semi-analytic (order-3) pruned-state-space gradient for
  `likelihood = "pruned"` instead of silently falling back to numerical
  finite differences — a stale gate had excluded it although the gradient
  was implemented and tested. This makes structural-SV order-2 models
  NUTS/HMC-able. Validated against `numDeriv` finite differences.
- **Three-asset fake-news accounting.** The forward derivative now streams
  the actual plus/minus policy and transition legs through the Young
  operator, preserving mass and all three asset first moments at active
  bounds (previously reconstructed and clipped a second symmetric policy
  pair). The analytic consumption Jacobian is assembled from the exact
  household budget, including the internal income-distribution response
  and the multiplicative foreign-price terms, so the six reported output
  Jacobians satisfy the aggregate linear budget identity to roundoff.
  Fixed-`Pi` price inputs (including transfer blocks, where the precondition
  is that income depends on the income state alone — satisfied by
  `y = w*e + Tr*omega` for any `Tr`) now carry an exact zero
  income-distribution response, avoiding an `O(N)/step` summation residue
  on large grids; a transfer block previously fell through to the
  accumulated path and carried a residue five orders larger.

# dynhr 0.9.0

This is a feature release. The headline addition since 0.8.1 is a complete
**heterogeneous-agent (HANK) estimation and welfare stack**, alongside an
**exact-Hessian curvature toolkit** and new **chain / determinacy diagnostics**.
Everything from 0.8.1 is retained and unchanged.

## Highlights

### Heterogeneous-agent (HANK) estimation

dynhr can now solve, estimate, and do welfare analysis on heterogeneous-agent
New-Keynesian models end to end, not just linearise them.

- **Solving.** Krusell–Smith / one-asset steady states and sequence-space
  Jacobians (`hank_ks_steady()`, `hank_ks_model()`, `hank_egm_solve()`,
  `hank_het_jacobian()`), a global nonlinear transition solver
  (`hank_td_nonlinear()`), a finite-Reiter linearisation path
  (`hank_finite_solve()` / `hank_reiter_statespace()`), and a
  **discount-heterogeneity mixture economy** (`hank_mixture_ks_steady()`,
  `hank_mixture_ks_assemble()`, `hank_mixture_ks_model()`).
- **Estimating.** `hank_mixture_joint_logpost()` evaluates a joint posterior
  over three information channels — aggregate macro dynamics (Kalman),
  the stationary cross-sectional wealth **level**, and the cross-sectional
  **response** to a price shock. Recommended posterior routes are the exact
  grid (`hank_mixture_sbc()`) or a Laplace approximation
  (`hank_mixture_laplace()`), **not** a diagonal random-walk sampler — the
  channels are strongly non-diagonal. `hank_mixture_emulator()` provides a
  distribution-agnostic surrogate, and `hank_mixture_sbc()` ships
  simulation-based-calibration certification (marginal and joint
  test-quantity ranks).
- **Welfare.** `hank_welfare_posterior()`, `hank_value_transition()`,
  `hank_mixture_welfare_pool()`, `hank_cev()`, and
  `hank_welfare_channels()` decompose consumption-equivalent welfare
  (interest- vs labour-income incidence) across the wealth distribution.
- **Identification.** The headline scientific finding baked into these tools:
  the discount-rate *spread* in a mixture economy is identified by the
  stationary wealth **level** (à la cstwMPC), not by the price-shock response
  (`hank_partial_id_level_response()`, `hank_reweight_level_metric()`).

### Exact-Hessian curvature and gradients

- **`posterior_hessian()`** is now exported, with four second-order-term
  methods — `t2_method = "loop"`, `"contract_once"`, `"hvp_solution"`, and
  `"adjoint_solution"` (the last two avoid forming the state-space second
  derivative; `"adjoint_solution"` is exact and finite-difference-free) — plus
  a `check_mode` guard that diagnoses a Σ_ε mismatch at the mode.
- `posterior_hessian_fd_grad()`, `laplace_log_marglik()`, and
  `make_posterior_grad(grad_method = "adjoint_solution")` round out the
  curvature stack.

### Diagnostics and determinacy

- **`chain_diagnostics()`** — MCMC chain summaries (R-hat, ESS, etc.).
- **`bk_distance()`** / `solution_pencil_spectrum()` — Blanchard–Kahn
  determinacy distance and the generalised-eigenvalue spectrum of the
  solution pencil.
- `kf_innovation_diagnostics()` — Kalman innovation whiteness checks.

### Pathological-DSGE estimation robustness

- `check_hessian_conditioning()`, `fd_safe_hessian()`, `profile_ci()`,
  `run_estimation_passport()`, and `make_loglik_contrib()` support inference
  on weakly-identified / ill-conditioned posteriors.

### Other additions

- Order-3 pruned state space: `pruned_state_space3()`, `pruned_ss_moments3()`,
  `pruned_ss_loglik3()`.
- `compute_fourth_cumulant(method = "closed_form")` — opt-in chain-exact /
  contemp-only fourth-cumulant trace (default remains `"window"`).
- A canonical occasionally-binding-constraint regime-path object
  (`new_obc_regime_path()` / `as_obc_regime_path()`).

## Two- and three-asset HANK households

- **Two-asset (liquid/illiquid) household**: liquid `b` at return `rb`,
  illiquid `a` at return `ra` with a convex adjustment cost
  (`hank_egm2_solve()` / `hank_het2_block()`, EGM over the `(Vb, Va)` pair,
  R + ~7x-faster C++ backends), a joint `(e, b, a)` distribution
  (`hank_forward_operator2()`, `hank_aggregate2()`), a fake-news Jacobian
  (`hank_het2_jacobian()` / `hank_het2_jacobian_nd()`), and GE composition
  (`hank_het2_block_spec()`, `hank_twoasset_steady()`, `hank_twoasset_model()`).
  With `ra > rb` this produces **wealthy hand-to-mouth** households
  (liquid-constrained while holding substantial illiquid wealth) — the
  Kaplan-Moll-Violante mechanism a one-asset HANK cannot represent. Ported
  from and validated against the sequence-jacobian reference implementation.
  Estimation works through the existing `hank_state_space()` /
  `hank_kalman_loglik()` / `hank_loglik_ar()` likelihood stack with no new
  code once a block's Jacobian feeds the block-DAG. Calibration gates:
  `hank_twoasset_grid_check()` (grid adequacy), `hank_theta_boundary_check()`
  (sequence-space horizon adequacy), `hank_twoasset_htm_stats()` (reports
  floor-mass and policy-constrained mass separately), and
  `hank_euler2_residual()` (FOC-residual oracle). Two-asset blocks are
  rejected loudly by every one-asset entry point and vice versa.
- **Three-asset household** (domestic liquid `d`, gross foreign `f`,
  domestic capital `a`): `hank_het3_block()` / `hank_egm3_solve()`,
  `hank_het3_jacobian()` / `hank_het3_jacobian_nd()`, `hank_euler3_residual()`
  (the Stage-4 budget/Euler diagnostic), and `hank_het3_manifest()` (a
  reproducibility record: version/commit, backend,
  grid, iterations, convergence gaps, wall time, peak RSS). Warm continuation
  across a calibration sweep via `Vd_init`/`Vf_init`/`Va_init` and `relax`
  (defaults reproduce the cold-start solve exactly); `hank_egm3_regrid_values()`
  interpolates marginal values across grid refinements. Non-convergence
  returns a `hank_het3_block_failed` object (with `strict = FALSE`) reporting
  iterations, `tol`, `relax` and both convergence gaps, rather than erroring
  or silently propagating into downstream Jacobians. `threads` controls
  cross-platform parallelism over the household loop; output is
  bit-identical at every thread count. The foreign valuation channel `px`
  (the price of foreign claims) is a first-class sequence-space input on
  `hank_het3_block()`, `hank_td3_nonlinear()`, and both Jacobians. Adjustment
  resources (`chi`, `phi`, aggregated as `CHI`/`PHI`) are reported so the
  household budget and the goods-market resource identity close exactly. A
  fix that forms consumption from the budget after fixing assets (matching
  the two-asset kernel's convention) changed the fixed point and is also
  roughly a **10x speedup**; every three-asset performance figure recorded
  before it is superseded.
- **Discrete-adjustment two-asset household** (fixed-cost/taste-shock,
  internal): `hank_egm2d_solve()`, `hank_het2d_block()` (per-branch
  policies, adjust probability `P`, `(V, Vb, Va)` envelope,
  `hank_forward_operator2d()`), `hank_td2d_nonlinear()` (the nonlinear
  transition), `hank_het2d_jacobian()` / `hank_het2d_jacobian_nd()`,
  `hank_het2d_block_spec()` (`kind = "het2d"`).
  A KiwiSaver-style locked-contribution channel (`phi_contrib`) resolves a
  fixed-cost participation trap in the stationary distribution. Remains
  internal (not exported) pending the full stack; C++ kernels are a
  documented follow-up.
- **Bond pricing and revaluation**: `hank_bond_block()` / `hank_bond_ss()` —
  a geometric (delta-coupon) bond, the cleanest revaluation instrument
  (fixed coupons make realized-return movements pure revaluation).
  `hank_td2_reval_decompose()` runs a household block on full /
  cash-flow-only / revaluation-only return paths to separate wealth from
  cash-flow incidence of a rate change.
- **Uniform lump-sum transfer `Tr`** on both one- and two-asset household
  blocks (`y = w*e + Tr`), a first-class Jacobian input with `Tr_path` on
  `hank_td_nonlinear()` / `hank_td2_nonlinear()`. Backward compatible:
  `Tr = 0` is the default and reproduces prior behavior exactly.
  `hank_impc()` now differences around the transfer-inclusive income
  definition.

## Debt-side channels: collateral, repricing, and the Fisher effect

- **Collateral-linked borrowing**: households may borrow beyond the
  unsecured liquid floor by pledging illiquid wealth,
  `b' >= b_grid[1] - theta_coll * a'`, at loan-to-value `theta_coll`, on
  the one-asset block. `theta_coll = 0` reproduces the pre-collateral
  solver exactly; a time-varying `theta_t` is supported, with a `theta_1`
  argument re-basing the date-1 distribution when the initial LTV differs
  from the block's steady-state value.
- **Staggered debt repricing**: `hank_repricing_block()` — the effective
  rate on the debt stock reprices by a fraction `phi_r` per period
  (`phi_r = 1` is instant repricing / the identity).
- **Borrowing wedge**: `r_minus` is a second aggregate rate input on
  `hank_het_block(r_minus = )`, letting deposit and borrowing rates diverge
  and reprice on different schedules. Wedge-unaware routines
  (`hank_impc()`, `hank_mpc()`, `hank_euler_residual()`,
  `hank_het_dist_jacobian()`, `hank_sam_reiter_linearize()`) now refuse a
  wedge-carrying block loudly rather than silently pricing the debt side at
  the saving rate.
- **Fisher channel**: `hank_fisher_block()` (anticipated inflation as an
  exact real-rate block) and `hank_liquid_reval_d0()` (the unanticipated
  date-1 nominal-stock revaluation, applied only to the liquid axis). The
  aggregate consumption response to surprise inflation follows the net
  nominal position of the household block, so per-group incidence is the
  robust object to report.

## General sequence-space engine and estimation

- **`hank_model()`** generalizes the Krusell-Smith-specific GE solve into
  an arbitrary directed-acyclic-graph of blocks: `hank_simple_block()`
  (analytic or finite-difference block Jacobian) and
  `hank_het_block_spec()` (fake-news), composed by the chain rule into GE
  Jacobians `H_U`/`H_Z` and solved via `hank_model_irf()`. `hank_ks_model()`
  builds Krusell-Smith through it, and `hank_nk_hank()` assembles a
  one-asset New Keynesian HANK (heterogeneous households, Taylor rule,
  Fisher equation, output-gap NKPC, bond clearing).
  `hank_model_nonlinear_irf()` is a general nonlinear perfect-foresight
  transition solver for any `hank_model` object. GE factorizations are
  cached and reused (keyed on `H_U`, bit-identical results); disable via
  `options(dynhr.hank_ge_cache = FALSE)`.
- **Exact-AR(1) stacked-covariance likelihood**: `hank_loglik_ar()` /
  `hank_loglik_aggregate_ar()` / `hank_autocov_ar()` replace the
  truncated-MA approximation's dropped shock-persistence tail with an
  exact closed-form stacked autocovariance — material (hundreds of
  log-points) near a unit root. `hank_theta_boundary_check()` diagnoses
  when the solve horizon is too short for a given persistence and warns
  (`check_boundary = FALSE` to silence). `hank_loglik_ar_grad()` supplies
  the analytic/semi-analytic gradient (sigma exact, measurement-error
  exact, persistence semi-analytic) and `make_posterior_grad_hank_ar()`
  wires it into a `make_log_posterior_hank()`-shaped closure. A `cache =`
  argument on `hank_loglik_ar()` / `hank_loglik_aggregate_ar()` memoizes
  per-shock autocovariance slabs (~15x faster warm evaluations);
  `make_log_posterior_hank(likelihood = "exact_ar")` enables it
  automatically. Note: `run_mode_finding()` / `run_full_estimation()` do
  not yet accept `likelihood = "exact_ar"` directly.
- **Kalman state-space bridge**: `hank_ma_state_space()` packs the
  aggregate MA representation into a `dsge_ss` state space, and
  `hank_loglik_ss()` evaluates it through dynhr's ordinary Kalman engine —
  giving a linearized HANK model access to the full estimation stack
  (smoother, frequency-domain likelihoods, priors, samplers, SBC).
- **Wealth heterogeneity axis in mixtures**: `hank_mixture_ks_steady_hetinc()`
  accepts optional per-type `eis` and `amin` (borrowing constraint on the
  shared asset grid), completing the wealth-heterogeneity path through
  general equilibrium (`hank_mixture_block_spec()`, `hank_mixture_dist()`);
  entries without `eis`/`amin` are unaffected.

## Diagnostics and robustness

- **`hank_determinacy()`** reports the conditioning of the GE Jacobian
  `H_U` and its determinacy verdict; `hank_model_irf()` warns on a
  near-singular `H_U` instead of returning a silent garbage IRF.
- **`hank_impc()`** — the intertemporal MPC matrix (Auclert-Rognlie-Straub
  intertemporal Keynesian cross) via fake-news. **`hank_mpc()`** now also
  reports MPC `by_income`. **`hank_distribution_stats()`** /
  **`hank_gini()`** — wealth/consumption Gini, top wealth shares,
  hand-to-mouth fraction, wealth percentiles.
- **`validate_hank_block()`** — a preflight diagnostic for a
  `hank_het_block` (or a mixture's block list): grid monotonicity, the
  shared Markov transition contract, forward-operator row-stochasticity,
  distribution validity, policy feasibility, consumption positivity,
  stored-aggregate identity checks, stationary-solver convergence, and
  optional R/C++ backend parity.
- **macOS Accelerate/vecLib compatibility fix**: a Hermitian eigendecomposition
  segfault under vecLib-backed BLAS (affecting the Whittle likelihood and
  Fisher-information paths) is fixed by an internal real-embedding
  eigensolver that is byte-identical on real input.
- Input-contract hardening across the HANK primitives: `hank_stationary_dist()`
  no longer reports false convergence on a degenerate (NaN) iterate,
  `hank_forward_operator()` / `hank_egm_solve()` validate the Markov/grid
  contract, and `hank_aggregate()` rejects dimension mismatches it would
  otherwise silently recycle.

## Other additions

- **`mode_trust_region()`** — a deterministic dogleg trust-region Newton
  posterior-mode optimizer (with an eigenvalue-clamp Hessian modification),
  complementing the stochastic/global finders and the `newrat` default.
- **`new_obc_regime_path()` / `as_obc_regime_path()`** — a canonical
  occasionally-binding-constraint regime-path object unifying the varied
  return shapes of dynhr's OBC solvers.
- Bugfix: `make_posterior_grad(likelihood = "cumulant")` used the wrong
  solve order (and consequently the wrong sign on some parameters) when
  `cumulant_orders` excluded 3 and 4; the default `cumulant_orders = 1:4`
  was unaffected. `make_posterior_grad()` now errors, instead of silently
  returning a Gaussian gradient, on likelihoods with no gradient path
  (`pskf`, `student_t`, `tpf`, `ppf`, `copf`).
- Foundational HANK / sequence-space Jacobian module (Auclert, Bardóczy,
  Rognlie & Straub 2021): household primitives `hank_income_rouwenhorst()`,
  `hank_asset_grid()`, `hank_egm_solve()`, `hank_euler_residual()`,
  `hank_forward_operator()`, `hank_stationary_dist()`, `hank_aggregate()`;
  the fake-news Jacobian `hank_het_jacobian()` (verified against the
  brute-force `hank_het_jacobian_nd()` to machine precision); GE via
  `hank_ks_steady()`, `hank_ks_linear_irf()`, `hank_ks_nonlinear_irf()`;
  the estimation bridge `hank_ma_coefficients()`, `hank_autocov()`,
  `hank_loglik_aggregate()`, `hank_simulate_aggregate()`; and diagnostics
  `hank_mpc()`, `hank_plot_jacobian()`, `hank_plot_irf()`.

# dynhr 0.8.1

This is the **first public release**. dynhr parses Dynare-style `.mod` files,
solves DSGE models by perturbation (orders 1-3) or global projection, and
estimates them with a wide range of (Bayesian) likelihood-based methods,
backed by an extensive model-diagnostics suite. Converted from a source-script
toolkit into a proper R package: `R CMD INSTALL`/`R CMD check` clean, parser
and solver behavior parity-tested against Dynare/Dynare.jl, vignettes
(`vignette("solving")`, `vignette("diagnostics")`), and a benchmarking harness
(`inst/benchmarks/`).

## Model solving

- Perturbation solving to third order, including `steady_state_model` blocks
  with derived parameters (analytic first- and second-order sensitivity when
  well-conditioned, with a finite-difference fallback) and pruned state
  spaces at order 2/3. `compile_model(param_deriv = "auto"/"on"/"off"/
  "second")` controls how much analytic parameter-derivative machinery is
  compiled (`"second"` enables the exact-Hessian codegen; `"off"` uses the FD
  fallback).
- A **global/projection solver**: Chebyshev collocation with Coleman time
  iteration and Gauss-Hermite expectations, for models with strong
  nonlinearity or occasionally-binding behaviour that perturbation cannot
  capture.
- **Markov-switching DSGE**: shock-variance switching (`ms_dsge_spec`,
  `ms_kim_filter` — Kim-Nelson GPB(2) filter with a Hamilton collapse) and a
  structural switching solver (`solve_ms_perturbation`, regime-coupled
  first-order perturbation).
- `simulate_model()` stochastically simulates a solved model at its compiled
  order (`simulate_model_order2()` / `simulate_model_order3()` for the
  pruned second-/third-order recursions); `compute_irfs()` computes impulse
  responses. `compute_moments_order2()` — deterministic order-2 unconditional
  moments (mean, full covariance, autocorrelations, per-shock variance
  decomposition) from the augmented pruned state space, without simulation
  noise.

## State-space filtering and likelihoods

- **`kalman_filter()`**: `method = "standard"/"dare"/"chandrasekhar"/
  "univariate"` (the last a Koopman-Durbin state-augmented sequential filter
  that handles singular innovation covariances by construction and is the
  automatic fallback whenever another method's `F` is singular or explodes),
  exact-diffuse initialization for unit-root/local-level models
  (`lik_init = "auto"/"diffuse"/"kappa"`), and `return_ll_contrib = TRUE` for
  per-period log-likelihood contributions. **`me_variance` now defaults to
  `0`** (previously `1e-8`) now that the univariate filter removes the need
  for a jitter crutch. **`kalman_smoother()`** returns per-period
  `filtered_cov`/`predicted_cov`/`smoothed_cov` state covariance arrays.
  `build_dsge_state_space()` constructs the `dsge_ss` object from a solved
  model; `new_dsge_ss()` / `ss_convert_timing()` give it an explicit
  lagged/current timing tag and convert between the two conventions.
- **Likelihood families**, all reachable via `make_log_posterior(likelihood
  = )` / `estimation_context()`: Gaussian; `"whittle"` (exact complex-spectral
  multivariate likelihood with analytic gradient, `debias = TRUE` by default
  — the Sykulski et al. 2019 expected-periodogram correction); `"cumulant"`
  (orders up to 4, with an analytic gradient); `"pskf"` — the Closed
  Skew-Normal Kalman filter (Guljanov, Mutschler & Trede 2026) for models
  with skew-normal shocks (`skew SHOCKNAME = expr` in the `shocks;` block;
  `pskf_smoother(method = "csn")` for the backward pass); `"tpf"` — the
  tempered particle filter (C++ `tpf_run_period_cpp` backend,
  `tpf_loglik_sd_preflight()` for a pre-run SD/particle-count check,
  `tpf_options` on `run_full_estimation()`); and the OBC (occasionally-
  binding-constraint) family — `kalman_filter_obc()` / `kalman_filter_obc_pkf()` /
  `kalman_filter_obc_inversion()` (exactly-identified deterministic
  inversion filter) and `ppf_likelihood()` / `make_log_posterior_obc_ppf()`
  (the bootstrap Piecewise Particle Filter), with a conditionally-optimal
  proposal (`proposal = c("bootstrap", "copf")`) and post-hoc importance
  reweighting (`ppf_reweight_posterior()`). `kalman_filter_student_t()`
  gives multivariate-t measurement/forecast errors for heavy-tailed data.
  A `power_posterior` tempering exponent is wired through all five main
  likelihood factories (Gaussian, cumulant, Whittle, PSKF, TPF).
- Shocks-block expressions (`stderr_expr`, `variance_expr`, `corr_expr`,
  `skew_expr`) are all re-evaluated against the current parameter vector on
  every draw, so estimated stderrs/correlations/skewness track theta rather
  than freezing at their parse-time snapshot.
- Correlated-pseudo-marginal proposals for particle likelihoods:
  `rwmh_cpm()` and `tpf_options$cpm_rho_u`.

## Analytic gradients and Hessians

- Tangent- and adjoint-mode Kalman-filter gradients (driving
  `make_posterior_grad()`) support time-varying `me_extra` (per-period
  measurement-error additions) and `shock_scale` (per-period shock-std
  scaling), and — for `steady_state_model`-derived parameters — the full
  first- and second-order sensitivity chain through the derived parameter.
  Analytic gradients are also available for the Whittle and cumulant
  likelihoods.
- **`posterior_hessian()`** — the exact second-order posterior Hessian,
  gated behind `compile_model(param_deriv = "second")`.
- Performance-sensitive gradient/Hessian kernels have C++ ports used
  automatically when available (`dynhr.use_rcpp` option controls the
  Kalman steady-state C++ path; set `FALSE` for the bit-exact R fallback).

## Bayesian estimation and sampling

- Samplers: RWMH, NUTS, SMC, and DIME, plus geometry-aware additions —
  `dynhr_mala()` (Laplace-MALA / simplified-manifold MALA), dense
  (non-diagonal) mass-matrix HMC/NUTS with `softabs_metric()` for indefinite
  Hessians, a gradient-only `monge_metric_fn()`, `whittle_fim()` (a
  frequency-domain Fisher-information metric), and `dynhr_chees()`
  (ChEES-HMC adaptive trajectory length, a tree-free alternative to NUTS).
- **`estimation_context()`** bundles the per-estimation options (likelihood,
  `lik_init`, `me_extra`, `shock_scale`, `freq_band`, `system_priors`,
  `tpf_options`, gradient policy); `run_mode_finding()`,
  `run_full_estimation()`, and `run_posterior_estimation()` drive mode-
  finding and sampling on top of it, with parallel chains/particles via
  `mirai` and live progress bars (`progress = ` on `run_mcmc_mirai()` /
  `run_mode_mirai()`).
- **`dynhr_plan()`** (an IRIS-"plan"-style judgment object) bundles
  `plan_tune()` (in-sample hard/soft tunes), `plan_condition()`
  (out-of-sample Waggoner-Zha conditioning), and `plan_scale_shock()`
  (heteroskedastic shock scaling); `conditional_forecast()` and
  `bayesian_conditional_forecast()` (posterior-draw conditional forecasting)
  both accept a `plan =`.
- GMM estimation with an analytic block-diagonal optimal weight matrix
  (alongside a Newey-West HAC estimator); `score_forecast()` for
  proper-scoring-rule forecast evaluation (CRPS, energy score, variogram
  score).
- Marginal likelihood / model comparison: `thames_mdd()` (the truncated
  harmonic-mean THAMES estimator, Metodiev et al. 2024) and
  `smc_model_tempered()` (two-stage SMC model tempering, Mlikota-Schorfheide
  2024).
- **`robust_confidence_set()`** — weak-identification-robust inference via
  the Andrews-Mikusheva (2015) LM/score test with test inversion.
- **`dynhr_sbc()`** — simulation-based calibration across
  `sampler = c("rwmh", "nuts", "smc", "dime")` and the likelihood families
  above (`nuts` is refused with the noisy `"tpf"` likelihood).

## Diagnostics

- **`run_all_diagnostics()`** runs a broad battery (D0 through D36),
  including: `$d0`, a static-Jacobian rank preflight for local
  well-posedness; variance/historical decomposition and smoothed shocks at
  the posterior mean (with opt-in `bayesian_irf = TRUE` posterior IRF
  credible bands); spectral identification with the full complex Hermitian
  Gram matrix (`spectral = c("exact", "companion")`); and system priors
  including `sp_spectral_peak()` (a spectral-density-peak-frequency prior
  feature).
- **Deep-parameter diagnostics**: an `@dynhr:deep` `.mod` metadata block
  classifies each parameter (`extract_mod_metadata()`, `build_deep_spec()`)
  as a deep primitive or reduced-form/auxiliary quantity; four new
  diagnostics assess structural-vs-reduced-form confounding
  ("borrowed identification"), policy-partitioned invariance (an
  operational Lucas-critique test, with cheap Laplace draws via
  `deep_laplace_draws()`), misspecification softness (a sandwich/
  information-matrix-equality check), and calibration deepness/tension for
  parameters that are calibrated rather than estimated. **
  `deep_parameter_passport()`** synthesizes all of this into a per-parameter
  A-F scorecard. `kalman_filter(return_ll_contrib = TRUE)` supplies the
  per-period likelihood contributions these diagnostics consume.
- `@dynhr:expectations` `.mod` metadata declares model-specific checks
  (`data_mean`, `data_sd`, `data_ratio`, `param_range`, `irf_sign`),
  evaluated by `diag_expectations()`. `run_diagnostics()` / `write_report()`
  are the public aliases for `run_all_diagnostics()` / `write_llm_report()`.

## Welfare and optimal policy

- `ramsey_model(order = 2)` — full second-order (augmented-system) Ramsey
  optimal policy. `osr()` (optimal simple rules) gains `order = ` and
  `planner_objective = ` to minimize expected welfare directly at order 2
  instead of an order-invariant quadratic variance loss.
  `conditional_welfare()` supports `method = "analytic"` (fast deterministic Taylor
  approximation), `"stochastic"` (Monte Carlo), or `"deterministic"` (the
  zero-shock path) for the order-2 state-dependent risk correction, and
  correctly conditions on a supplied `initial_state`.
