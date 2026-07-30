# dynhr 0.9.1

A feature release. The headline additions since 0.9.0 are **multi-asset
heterogeneous-agent households**, **stochastic volatility on shocks**, and a
**cross-machine benchmark**. Everything in 0.9.0 is retained.

## Multi-asset HANK households

0.9.0 shipped the one-asset (Krusell-Smith) household. 0.9.1 adds two- and
three-asset households with per-asset adjustment costs, each with the full
sequence-space toolchain.

- **Two-asset** (liquid / illiquid): `hank_het2_block()`, `hank_het2_jacobian()`,
  `hank_td2_nonlinear()`, `hank_egm2_solve()`, `hank_forward_operator2()`, plus
  a deposit variant (`hank_het2d_block()`, `hank_het2d_jacobian()`) and a
  ready-made GE model (`hank_twoasset_model()`, `hank_twoasset_steady()`).
- **Three-asset** (domestic / foreign / illiquid, with a foreign price `px`
  entering both the purchase and payoff legs of the budget):
  `hank_het3_block()`, `hank_het3_jacobian()`, `hank_td3_nonlinear()`,
  `hank_egm3_solve()`. `hank_het3_jacobian_checkpoint()` and
  `hank_het3_jacobian_spot()` support long-running builds on large grids.
- **Numerical-differentiation oracles** ship alongside the analytic Jacobians
  (`hank_het2_jacobian_nd()`, `hank_het3_jacobian_nd()`) rather than being
  test-only, so a user can verify an analytic Jacobian on their own block.
- **Reproducibility**: `hank_het_manifest()` / `hank_het2_manifest()` /
  `hank_het3_manifest()` and the matching `*_fingerprint()` functions record and
  hash the exact block a result came from.
- Building blocks for richer environments: `hank_bond_block()`,
  `hank_fisher_block()`, `hank_repricing_block()`, `hank_employment_income()`,
  `hank_incidence_earnings()`, and a search-and-matching Reiter linearisation
  (`hank_sam_reiter_statespace()`).

## HANK estimation with an exact sequence-space score

- `hank_loglik_ar()` and `hank_loglik_ar_grad()` give an exact-AR likelihood and
  its analytic gradient through the structural block
  (`hank_loglik_ar_structural_grad()`), with `make_posterior_grad_hank_ar()`
  wiring it into the gradient samplers.
- `hank_run_estimation()` / `hank_run_mode_finding()` drive the pipeline;
  `hank_model_dtheta()` propagates parameter derivatives through the model DAG.

## Stochastic volatility on shocks

Declare independent AR(1) log-variance processes on any subset of a model's
shocks with `stochastic_volatility()` / `sv_entry()`, and estimate them with a
Rao-Blackwellised particle filter (`make_log_posterior_sv_rbpf()`). The
degenerate-volatility limit reproduces the exact Kalman likelihood.
`sv_rbpf_sbc()` provides simulation-based-calibration certification, and
`kf_step()` / `kf_stationary_init()` expose the single-step filter primitives.

## Benchmarking

`dynhr_benchmark()` runs a fixed Smets & Wouters (2007) estimation workload --
36 estimated parameters, 7 observables, 160 quarters, started at the published
posterior mode -- through random-walk Metropolis at a sweep of core counts. It
reports per-chain and aggregate throughputs, all normalised per draw or per
second so runs with different draw counts stay comparable, together with a
workload fingerprint. `dynhr_system_info()` records CPU, RAM, OS, R build and
the BLAS/LAPACK actually linked; two results are only comparable if those agree.

The Smets-Wouters model, data, published mode and mode Hessian ship in
`inst/extdata/models` from the AEA replication deposit openicpsr-116269-V1
(BSD-3-Clause for code, CC BY 4.0 for data); see `sw2007_SOURCE.md`.

## Newly exported helpers

Previously internal, exported on request from downstream users:
`make_log_posterior()`, `extract_prior_spec()`, `log_prior()`,
`solve_steady_state()`, `solve_lyapunov()`, `apply_theta_to_params()`,
`build_param_transform()`, `make_transformed_logpost()`,
`make_transformed_grad()`.

## Other additions

- `bk_wall_transform()` -- a reparameterisation that removes the divergences
  gradient samplers suffer near the Blanchard-Kahn determinacy boundary.
- `mdd_calibration()` -- marginal-data-density estimators calibrated against
  analytic and quadrature ground truth.

## Fixes

- README: `theoretical_moments()` and `ramsey_regime_markov()` were named in the
  0.9.0 capability map but do not exist; the functions are `compute_moments()`
  and `ramsey_regime_independent()`.
- Non-ASCII characters in roxygen documentation replaced with ASCII or `\eqn{}`
  markup, so the PDF manual builds.

# dynhr 0.9.0

A feature release. The headline addition since 0.8.1 is a complete
**heterogeneous-agent (HANK) estimation and welfare stack**, alongside an
**exact-Hessian curvature toolkit** and new **chain / determinacy diagnostics**.
Everything in 0.8.1 is retained.

## Heterogeneous-agent (HANK) estimation

dynhr can now solve, estimate, and do welfare analysis on heterogeneous-agent
New-Keynesian models end to end, using the sequence-space Jacobian.

- **Solving** — Krusell–Smith / one-asset steady states and sequence-space
  Jacobians (`hank_ks_steady()`, `hank_ks_model()`, `hank_egm_solve()`,
  `hank_het_jacobian()`), a global nonlinear transition solver
  (`hank_td_nonlinear()`), a finite-Reiter linearisation path
  (`hank_finite_solve()`, `hank_reiter_statespace()`), and a
  discount-heterogeneity **mixture economy** (`hank_mixture_ks_steady()`,
  `hank_mixture_ks_assemble()`, `hank_mixture_ks_model()`).
- **Estimating** — `hank_mixture_joint_logpost()` evaluates a joint posterior
  over three information channels: aggregate **macro** dynamics (Kalman), the
  stationary cross-sectional wealth **level**, and the cross-sectional
  **response** to a shock. The recommended posterior routes are the exact grid
  (`hank_mixture_sbc()`) or a Laplace approximation (`hank_mixture_laplace()`) —
  not a diagonal random-walk sampler, as the channels are strongly
  non-diagonal. `hank_mixture_emulator()` gives a distribution-agnostic
  surrogate, and `hank_mixture_sbc()` ships simulation-based-calibration
  certification.
- **Welfare** — `hank_welfare_posterior()`, `hank_value_transition()`,
  `hank_mixture_welfare_pool()`, `hank_cev()`, and `hank_welfare_channels()`
  decompose consumption-equivalent welfare (interest- vs labour-income
  incidence) across the wealth distribution.
- **Identification** — the headline result baked into these tools: the
  discount-rate *spread* in a mixture economy is identified by the stationary
  wealth **level** (à la cstwMPC), not the price-shock response
  (`hank_partial_id_level_response()`, `hank_reweight_level_metric()`).

See `vignette("hank")` for a worked example.

## Exact-Hessian curvature and gradients

- **`posterior_hessian()`** is now exported, with four second-order-term methods
  (`t2_method = "loop"`, `"contract_once"`, `"hvp_solution"`,
  `"adjoint_solution"`; the last two avoid forming the state-space second
  derivative, and `"adjoint_solution"` is exact and finite-difference-free), plus
  a `check_mode` guard.
- `posterior_hessian_fd_grad()`, `laplace_log_marglik()`, and
  `make_posterior_grad(grad_method = "adjoint_solution")`.

## Diagnostics and determinacy

- **`chain_diagnostics()`** — MCMC chain summaries (split-R-hat, ESS, MCSE).
- **`bk_distance()`** / `solution_pencil_spectrum()` — Blanchard–Kahn
  determinacy distance and the solution-pencil generalised-eigenvalue spectrum.
- `kf_innovation_diagnostics()` — Kalman innovation whiteness checks.

## Pathological-DSGE estimation robustness

- `check_hessian_conditioning()`, `fd_safe_hessian()`, `profile_ci()`,
  `run_estimation_passport()`, and `make_loglik_contrib()` for inference on
  weakly-identified / ill-conditioned posteriors.

## Other additions

- Order-3 pruned state space: `pruned_state_space3()`, `pruned_ss_moments3()`,
  `pruned_ss_loglik3()`.
- `compute_fourth_cumulant(method = "closed_form")` (opt-in; default `"window"`).
- A canonical occasionally-binding-constraint regime-path object
  (`new_obc_regime_path()` / `as_obc_regime_path()`).

As in 0.8.1, this is a curated build for distribution: the test suite,
model/data fixtures, replication material, and third-party tooling are
maintained separately and are not part of the released package.


# dynhr 0.8.1

First public release.

`dynhr` is a self-contained R toolkit for medium-scale DSGE models: it parses
Dynare `.mod` files; solves by perturbation (orders 1–5) with occasionally-
binding-constraint methods (OccBin / MCP / LCP / Boehl); runs Bayesian
estimation (random-walk Metropolis, SMC, HMC/NUTS, with mode-finding and prior
tooling); and reports an identification → convergence → fit → narrative
diagnostic battery. The solver and filter are validated against Dynare and
Dynare.jl.

This is a curated build for distribution: the test suite, model/data fixtures,
replication material, and third-party tooling are maintained separately and are
not part of the released package.
