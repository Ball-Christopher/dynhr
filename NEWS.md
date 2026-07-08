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
