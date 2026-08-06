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
  channels (`hank_mixture_joint_logpost()`, with exact-grid or Laplace
  posterior routes and SBC certification), and consumption-equivalent welfare
  (`hank_welfare_posterior()`).
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
  with a Rao-Blackwellised particle filter (`make_log_posterior_sv_rbpf()`),
  SBC-certified via `sv_rbpf_sbc()`.
- **Benchmarking** — `dynhr_benchmark()` runs a fixed Smets-Wouters (2007)
  estimation workload across a sweep of core counts and reports normalised
  throughputs plus full system information, so two machines can be compared.
- **Occasionally-binding constraints** — OccBin/MCP/LCP and Boehl-style solvers,
  plus Ramsey/OSR/discretionary optimal-policy machinery.
- **Diagnostics** — an identification → convergence → fit → narrative battery
  that renders a Markdown report, plus MCMC chain summaries
  (`chain_diagnostics()`) and Blanchard–Kahn determinacy distance
  (`bk_distance()`).

A parity test suite validates the solver and filter against **Dynare 7.0** and
**Dynare.jl** golden files, and the high-order sigma terms against a closed-form
ground-truth model.

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
macOS the Xcode command-line tools. Julia, Dynare, and Octave are **optional**
and only needed to regenerate parity goldens.

### Build configuration that affects PERFORMANCE (and reproducibility)

dynhr supplies no compiler flags of its own. `src/Makevars` is a single line —
`PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` — so optimisation level, BLAS
and LAPACK all come from **your R installation's** configuration, not from the
package. That makes a few otherwise-invisible choices matter, especially if you
are timing the heterogeneous-agent (HANK) routines or comparing runs across
machines.

**1. `devtools::load_all()` compiles at `-O0`. Never benchmark on it.**
This is the single most common way to get badly wrong numbers from this
package: `load_all()` builds unoptimised objects into `src/`, and they *persist*
and can be picked up by a subsequent install. A HANK Jacobian measured this way
has been observed ~6x slower than the same code installed normally. Before any
timing, and after any `load_all()`:

```sh
R CMD INSTALL --preclean .   # --preclean is what discards the -O0 objects
```

`R CMD config CXX20FLAGS` should show `-O2` (the R default). Correctness is
unaffected either way — only speed.

**2. Check which BLAS you are actually linked against.** On macOS, R can be
configured to use Apple's Accelerate (vecLib) instead of the reference BLAS,
which is substantially faster for the dense linear algebra in the solvers and
Kalman filters. It is a symlink, and it is easy not to know which one you have:

```sh
ls -l "$(R RHOME)/lib/libRblas.dylib"      # -> libRblas.vecLib.dylib if Accelerate
otool -L "$(R RHOME)/library/dynhr/libs/dynhr.so" | grep -i accelerate
```

On Linux the analogue is whether R is linked against OpenBLAS/MKL or the
reference BLAS (`sessionInfo()` reports it). Two builds of *identical* dynhr
source can differ severalfold in wall time on this alone, so state it when
reporting timings.

**3. macOS arm64 needs a gfortran whose runtime matches R's `FLIBS`.**
R 4.6 arm64 expects gfortran 14.2 (`/opt/gfortran`); a mismatched one produces
link errors or, worse, a package that loads but misbehaves. `otool -L` on the
installed `.so` should show `libgfortran.5.dylib` from that prefix.

**4. Build and check in a UTF-8 locale** (`LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8`),
or `R CMD build`/`check` can fail on encoded characters in the sources.

**5. Worker threads for the compiled HANK kernels** are resolved as: an explicit
`threads` argument, then `getOption("dynhr.hank_threads")`, then a
machine-derived default. Output is **bit-identical at every thread count** — it
is purely a throughput knob. Two things worth knowing: `R CMD check` sets
`_R_CHECK_LIMIT_CORES_`, which clamps the resolved count to 2 (so a "slow" gate
run may simply be running two-wide); and the machine-derived default is tuned
for the three-asset kernel, so on **two-asset** problems an explicit smaller
count is often faster than the default. To see what a long run actually
resolved to:

```r
options(dynhr.hank_report_threads = TRUE)   # reports on CHANGE, not per call
```

**6. Reproducing a specific build.** An installed dynhr records the commit it
was built from in `inst/GIT_COMMIT`: line 1 is the 40-hex SHA, line 2 is
`version: <DESCRIPTION Version>`, so commit-and-version agreement can be checked
from the installed files alone, with no git and no network.
`hank_het3_manifest()` surfaces it as `git_commit` for run manifests.

```r
readLines(system.file("GIT_COMMIT", package = "dynhr"))
```

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
by `run_full_estimation()`. See the vignettes for worked examples.

## Vignettes

- `vignette("dynhr")` — getting started
- `vignette("solving")` — steady state, perturbation, IRFs, moments, Kalman
- `vignette("estimation")` — priors, mode-finding, MCMC/SMC/NUTS, exact-Hessian curvature
- `vignette("hank")` — heterogeneous-agent (HANK) solving, mixture estimation, welfare
- `vignette("diagnostics")` — the diagnostic battery and custom expectations
- `vignette("mod-conversion")` — Dynare `.mod` compatibility notes
- `vignette("mod-syntax")` — dynhr `.mod` syntax reference
- `vignette("sbc-matrix")` — the SBC (Simulation-Based Calibration) coverage matrix

## Where things are

- `R/` — package source
- `src/` — Rcpp/Armadillo backends (folded Faà-di-Bruno compose, Kalman steady
  state, sparse MCP solve), each with a pure-R fallback toggled by
  `options(dynhr.use_rcpp = )`
- `inst/extdata/models/` — reference DSGE models for examples and tests
- `inst/extdata/golden/` — Dynare/Dynare.jl reference outputs for parity tests
- `inst/pipelines/` — full estimation pipeline scripts
- `inst/julia/`, `inst/octave/` — scripts that regenerate golden files
- `tests/testthat/` — unit and parity tests

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
  through `run_full_estimation`), `dynhr_smc2()` (SMC^2: outer theta-tempering
  around an inner noisy-but-unbiased particle-filter likelihood — `tpf` or
  `sv_rbpf`).
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
- Employment margin: `hank_employment_income()` (two-state E/U) and
  `hank_employment_income3()` (three-state E/U/N, with a separate
  not-in-labour-force state and its own job-finding/separation rates).
- Income incidence for a heterogeneous block: `hank_incidence_earnings()`
  (normalises an income profile against a block's grid/transition matrix).
- GE & sequence-space Jacobians: `hank_ks_steady()`, `hank_ks_model()`,
  `hank_het_jacobian()`, `hank_td_nonlinear()` (global nonlinear transition),
  `hank_reiter_statespace()` (finite Reiter linearisation).
- Distribution Jacobians (perturbation of the cross-sectional distribution
  itself): `hank_het_dist_jacobian()` / `hank_het2_dist_jacobian()` /
  `hank_het3_dist_jacobian()` (one/two/three-asset), each with a
  `_nd()` finite-difference check counterpart, plus
  `hank_mixture_dist_jacobian()` for the mixture economy.
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
- Welfare: `hank_welfare_posterior()`, `hank_cev()`, `hank_value_transition()`,
  `hank_welfare_channels()`, `hank_mixture_welfare_pool()`.
- Identification result baked into the tools: in a mixture economy the discount
  *spread* is identified by the stationary wealth **level** (à la cstwMPC), not
  by the price-shock response (`hank_partial_id_level_response()`).

### Optimal policy
- Ramsey: `ramsey_model()` (augmented FOC system, Bodenstein–Guerrieri),
  `ramsey_nn1()` ((n, n+1) approximation, Gross–Hansen), `ramsey_obc_pf()` /
  `ramsey_obc_pwlinear()` (with OBC), `ramsey_regime_deterministic()` /
  `ramsey_regime_independent()` (regime-dependent: deterministic switch /
  Markov-switching).
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
