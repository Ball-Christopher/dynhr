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
