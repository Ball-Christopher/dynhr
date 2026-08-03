# dynhr

A self-contained R toolkit for parsing, solving, estimating, and diagnosing
medium-scale DSGE models. It reads Dynare-format `.mod` files and provides, in
pure R (with optional Rcpp/Armadillo acceleration):

- **Solving** — steady state, perturbation to **orders 1–5** (deterministic and
  stochastic/`sigma` corrections), IRFs, moments, and stochastic simulation.
- **Heterogeneous agents (HANK)** — solve, estimate, and do welfare analysis on
  heterogeneous-agent New-Keynesian models: Krusell–Smith steady states and
  sequence-space Jacobians (`hank_ks_steady()`, `hank_het_jacobian()`), a
  discount-heterogeneity mixture economy (`hank_mixture_ks_model()`), a joint
  posterior over macro / cross-sectional-level / cross-sectional-response
  channels (`hank_mixture_joint_logpost()`, with exact-grid or Laplace posterior
  routes and SBC certification), and consumption-equivalent welfare
  (`hank_welfare_posterior()`). See `vignette("hank")`.
- **Filtering & smoothing** — Kalman filter (standard, Chandrasekhar, DARE
  oracle) and smoother, plus a piecewise Kalman filter for occasionally-binding
  constraints. (The default filter adds **no** measurement-error variance
  (`me_variance = 0`) for exact Dynare parity, falling back automatically to the
  univariate filter on a singular innovation covariance; see `?kalman_filter`.)
- **Estimation** — random-walk Metropolis-Hastings, sequential Monte Carlo, and
  NUTS, with mode-finding (Nelder-Mead, CMA-ES, JADE), prior tooling, and an
  exact-Hessian curvature stack (`posterior_hessian()` with finite-difference-
  free adjoint second-order terms, `laplace_log_marglik()`, `profile_ci()`).
- **Stochastic volatility on shocks** — declare AR(1) log-variance processes on
  any subset of a model's shocks (`stochastic_volatility()`) and estimate them
  with a Rao-Blackwellised particle filter (`make_log_posterior_sv_rbpf()`).
- **Benchmarking** — `dynhr_benchmark()` runs a fixed Smets-Wouters (2007)
  estimation workload across a sweep of core counts and reports normalised
  throughputs plus full system information, so two machines can be compared.
- **Occasionally-binding constraints** — OccBin/MCP/LCP and Boehl-style solvers,
  plus Ramsey/OSR/discretionary optimal-policy machinery.
- **Diagnostics** — an identification → convergence → fit → narrative battery
  that renders a Markdown report, plus MCMC chain summaries
  (`chain_diagnostics()`) and Blanchard–Kahn determinacy distance
  (`bk_distance()`).

The solver and filter are validated against **Dynare 7.0** and **Dynare.jl**,
and the high-order sigma terms against a closed-form ground-truth model. The
parity test suite and its fixtures are maintained separately from this released
package.

## Install

```r
# install.packages("pak")
pak::pak("Ball-Christopher/dynhr")

library(dynhr)
```

(`devtools::install_github()` is deprecated as of devtools 2.5.0; if you
prefer not to use pak, `remotes::install_github("Ball-Christopher/dynhr")`
still works.)

The package compiles a small amount of C++ (`src/`, via Rcpp + RcppArmadillo,
using the C++20 standard); a C++ toolchain is required to install from
source — on Windows that means **Rtools43 or newer** (i.e. R >= 4.3), on
macOS the Xcode command-line tools. Julia is **optional** and only needed for
the Dynare.jl interop.

## Quick start

```r
library(dynhr)

mod      <- parse_mod(system.file("extdata/models/rbc.mod", package = "dynhr"))
compiled <- compile_model(mod)
steady   <- solve_steady(compiled, mod$param_values)
dr       <- solve_perturbation(mod, compiled, steady$values, mod$param_values)

print(dr)                                              # decision rules
sim <- simulate_model(dr, n_periods = 200, model = mod)
```

Higher orders (set `max_order` when compiling, then pass `order` to the solver):

```r
compiled <- compile_model(mod, max_order = 3L)
dr3      <- solve_perturbation(mod, compiled, steady$values,
                               mod$param_values, order = 3L)
```

End-to-end Bayesian estimation (mode-finding + sampling + diagnostics) is driven
by `run_full_estimation()`.

## Where things are

- `R/` — package source
- `src/` — Rcpp/Armadillo backends (folded Faà-di-Bruno compose, Kalman steady
  state, sparse MCP solve), each with a pure-R fallback toggled by
  `options(dynhr.use_rcpp = )`
- `inst/extdata/models/` — a few reference DSGE models used by the examples,
  plus the Smets-Wouters (2007) model, data and published mode that
  `dynhr_benchmark()` runs (provenance and licensing in `sw2007_SOURCE.md`)
- `inst/templates/` — report templates for the diagnostic battery

## A note on AI and reliability

AI tools were used extensively in the development of this package. The code has
been tested thoroughly throughout development, but it remains **experimental**
and may contain errors — **use at your own risk**, and validate results against
a trusted reference for any consequential use. `dynhr` is part of the author's
ongoing experimentation with AI-assisted development tools, and feedback and bug
reports are welcome.

## Reference: capability map and function index

`dynhr` is a self-contained R implementation of the full DSGE modelling
workflow — a Dynare alternative that reads the same `.mod` files but runs
entirely in R (with optional Rcpp/Armadillo kernels). This section is a complete
map of the public API, grouped by task, so the package can be understood and
used from this page alone.

### Model input
- `parse_mod(file_or_text)` — parse a Dynare-format `.mod` file (or a string)
  into a model object. Supports `var`, `varexo`, `parameters`, `model`,
  `initval` / `steady_state_model`, `shocks`, `estimated_params`, and OBC tags.
- `compile_model(model, max_order = 1L)` — compile symbolic derivatives (to
  order 5) and fast evaluators; raise `max_order` for higher-order perturbation.

### Solving
- `solve_steady(compiled, params)` — steady state (analytic `steady_state_model`
  seed or numerical Newton).
- `solve_perturbation(model, compiled, ss, params, order = 1L)` — perturbation
  decision rules to **orders 1–5** (deterministic + stochastic `sigma` terms).
- `compute_irfs(dr, model, n_periods = 40L)`, `simulate_model()`,
  `compute_moments()`, `stoch_simul()` — IRFs, stochastic simulation,
  model-implied moments.
- Determinacy: `bk_distance()` (signed Blanchard–Kahn distance to the
  indeterminacy / no-solution boundary), `solution_pencil_spectrum()`
  (generalised-eigenvalue spectrum of the solution pencil).
- Higher-order pruned state space: `pruned_state_space()` /
  `pruned_ss_moments()` / `pruned_ss_loglik()` (order 2) and
  `pruned_state_space3()` / `pruned_ss_moments3()` / `pruned_ss_loglik3()`
  (order 3).

### Occasionally-binding constraints (OBC / ZLB)
- Solvers: `occbin_solve_path()`, `mcp_solve_path()` (mixed complementarity),
  `boehl_solve_regime_path()`.
- IRFs / decompositions: `compute_irfs_obc()`, `historical_decomposition_obc()`.
- Filtering under OBC: `kalman_filter_obc_pkf()` (piecewise Kalman filter) and
  particle filters (`make_log_posterior_obc_ppf()`, `make_log_posterior_tpf()`).
- Canonical regime-path object: `new_obc_regime_path()` / `as_obc_regime_path()`.

### Filtering & smoothing
- `kalman_filter(Y, dr, model, params, obs_vars, method = "auto")` — Gaussian
  likelihood; methods `"auto"`, `"standard"`, `"chandrasekhar"`, `"dare"`,
  `"univariate"`; stationary / exact-diffuse / auto initialisation (`lik_init`).
  Default `me_variance = 0` for exact Dynare parity, with an automatic
  univariate fallback on a singular innovation covariance.
- `kalman_smoother()` — RTS smoother.
- Skewed / particle likelihoods: `make_log_posterior_pskf_order2()` (pruned
  skewed Kalman filter), `make_log_posterior_tpf()` (tempered particle filter).

### Bayesian estimation
- `prior_spec(model)` — extract priors from the `estimated_params` block.
- `make_posterior(model, data, prior_spec, obs_vars, compiled, me_variance = 0)`
  — build a log-posterior closure (Gaussian, cumulant, or Whittle likelihood).
- `find_mode(log_post_fn, theta_init, prior_spec, method = "newrat")` — posterior
  mode; `method` ∈ `"newrat"` (csminwel, = Dynare `mode_compute 4`), `"cmaes"`,
  `"nelder"`, `"jade"`, `"combined"`.
- Samplers: `mcmc()` (random-walk Metropolis), `smc()` (sequential Monte Carlo,
  returns a log-marginal-likelihood estimate), `nuts()` (plus MALA / HMC / CHEES
  through `run_full_estimation`).
- `run_full_estimation(...)` — one-call pipeline (mode → sample → diagnostics),
  multi-chain with Gelman–Rubin convergence.

### Curvature, model evidence, and weak identification
- `posterior_hessian(log_post, theta, t2_method = "adjoint_solution")` — exact
  analytic Hessian. `t2_method` ∈ `"loop"`, `"contract_once"`, `"hvp_solution"`,
  `"adjoint_solution"`; the last two never form the state-space second
  derivative, and `"adjoint_solution"` is exact and finite-difference-free.
- `laplace_log_marglik()` — Laplace model evidence from the mode plus Hessian.
- `make_posterior_grad(grad_method = "adjoint_solution")` — analytic score for
  gradient samplers.
- `check_hessian_conditioning()`, `fd_safe_hessian()`, `profile_ci()`,
  `run_estimation_passport()` — weak-identification / ill-conditioning tooling.
- `chain_diagnostics(draws)` — split-R-hat, bulk/tail ESS, and Monte-Carlo
  standard errors from a draws matrix or chain list.

### Heterogeneous agents (HANK)
- Household / income: `hank_income_rouwenhorst()`, `hank_asset_grid()`,
  `hank_egm_solve()`, `hank_stationary_dist()`, `hank_mpc()`.
- GE & sequence-space Jacobians: `hank_ks_steady()`, `hank_ks_model()`,
  `hank_het_jacobian()`, `hank_td_nonlinear()` (global nonlinear transition),
  `hank_reiter_statespace()` (finite Reiter linearisation).
- Discount-heterogeneity mixture economy: `hank_mixture_ks_steady()`,
  `hank_mixture_ks_assemble()`, `hank_mixture_ks_model()`,
  `hank_mixture_agg_irf()`.
- Estimation: `hank_mixture_joint_logpost()` — a joint posterior over three
  channels (aggregate **macro** dynamics, the stationary cross-sectional wealth
  **level**, and the cross-sectional **response** to a shock). Recommended
  posterior routes: the exact grid (`hank_mixture_sbc()`) or a Laplace
  approximation (`hank_mixture_laplace()`), not a diagonal random-walk sampler.
  `hank_mixture_emulator()` is a distribution-agnostic surrogate;
  `hank_mixture_sbc()` ships simulation-based-calibration certification.
- Multi-asset households: liquid/illiquid two-asset blocks (`hank_het2_block()`,
  `hank_het2_jacobian()`, `hank_td2_nonlinear()`) and a three-asset block with
  domestic, foreign and illiquid claims plus per-asset adjustment costs
  (`hank_het3_block()`, `hank_het3_jacobian()`, `hank_td3_nonlinear()`). Both
  carry an exact numerical-differentiation oracle (`*_jacobian_nd()`) and a
  reproducibility fingerprint / manifest (`hank_het3_fingerprint()`,
  `hank_het3_manifest()`).
- Welfare: `hank_welfare_posterior()`, `hank_cev()`, `hank_value_transition()`,
  `hank_welfare_channels()`, `hank_mixture_welfare_pool()`.
- Identification result baked into the tools: in a mixture economy the discount
  *spread* is identified by the stationary wealth **level** (à la cstwMPC), not
  by the price-shock response (`hank_partial_id_level_response()`).

### Optimal policy
- Ramsey: `ramsey_model()` (augmented FOC system, Bodenstein–Guerrieri),
  `ramsey_nn1()` ((n, n+1) approximation, Gross–Hansen), `ramsey_obc_pf()` /
  `ramsey_obc_pwlinear()` (with OBC), `ramsey_regime_deterministic()` /
  `ramsey_regime_independent()` (regime-dependent).
- `osr()` — optimal simple rules; `discretionary_policy()` — Markov-perfect
  discretion; `nash_ramsey_cooperative()` / `nash_ramsey_openloop()` — policy
  games.
- `opp_sufficient_stats()` — sufficient-statistics optimal policy (Barnichon &
  Mesters 2023).
- Welfare: `welfare_compute()`, `welfare_ce_diff()`, `welfare_cost_of_rule()`,
  `welfare_decompose()`.

### Diagnostics
- `run_diagnostics()` — the D1–D30 battery (identification → convergence → fit →
  narrative), rendered by `write_report()`.
- Standalone helpers: `chain_diagnostics()`, `bk_distance()`,
  `kf_innovation_diagnostics()`, `solution_pencil_spectrum()`.

### Benchmarking
- `dynhr_benchmark()` — run a fixed Smets-Wouters (2007) estimation workload
  (36 estimated parameters, 7 observables, 160 quarters) through random-walk
  Metropolis at a sweep of core counts. Reports per-chain and aggregate
  throughputs, all normalised per draw or per second so runs with different
  draw counts stay comparable, plus a workload fingerprint.
- `dynhr_system_info()` — CPU, RAM, OS, R build and the BLAS/LAPACK actually
  linked. Two benchmark results are only comparable if these agree.

### Cookbook

Solve and simulate:

```r
library(dynhr)
m   <- parse_mod(system.file("extdata/models/rbc.mod", package = "dynhr"))
cm  <- compile_model(m, max_order = 2L)
ss  <- solve_steady(cm, m$param_values)
dr  <- solve_perturbation(m, cm, ss$values, m$param_values, order = 2L)
sim <- simulate_model(dr, n_periods = 200, model = m)
irf <- compute_irfs(dr, model = m, n_periods = 40L)
```

Estimate (mode → exact-Hessian curvature → NUTS → diagnostics):

```r
priors <- prior_spec(m)
theta0 <- setNames(priors$mean, priors$name)
lp     <- make_posterior(m, data = Y, prior_spec = priors,
                         obs_vars = obs, compiled = cm, me_variance = 0)
mode   <- find_mode(log_post_fn = lp, theta_init = theta0, prior_spec = priors)
H      <- posterior_hessian(lp, mode$theta_mode, t2_method = "adjoint_solution")
chains <- nuts(log_post_fn = lp, theta0 = mode$theta_mode,
               n_draws = 4000L, n_warmup = 2000L)
chain_diagnostics(chains$chain)
```

HANK mixture (discount heterogeneity):

```r
inc   <- hank_income_rouwenhorst(rho = 0.966, sigma = 0.5, n = 7)
ag    <- hank_asset_grid(amax = 200, n = 500, amin = 0)
mks   <- hank_mixture_ks_steady(ag, inc$Pi, inc$e, betas = c(0.985, 0.96),
                                omega = c(0.6, 0.4), eis = 1,
                                alpha = 0.36, delta = 0.025)
model <- hank_mixture_ks_model(mks, T_h = 300L)
irf   <- hank_model_irf(model, 0.01 * 0.9^(0:299))
# joint posterior over macro / level / response channels:
#   ?hank_mixture_joint_logpost  and  vignette("hank")
```

## Licence

MIT — see [LICENSE](LICENSE).
