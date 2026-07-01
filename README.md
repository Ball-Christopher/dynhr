# dynhr

A self-contained R toolkit for parsing, solving, estimating, and diagnosing
medium-scale DSGE models. It reads Dynare-format `.mod` files and provides, in
pure R (with optional Rcpp/Armadillo acceleration):

- **Solving** — steady state, perturbation to **orders 1–5** (deterministic and
  stochastic/`sigma` corrections), IRFs, moments, and stochastic simulation.
- **Filtering & smoothing** — Kalman filter (standard, Chandrasekhar, DARE
  oracle) and smoother, plus a piecewise Kalman filter for occasionally-binding
  constraints. (Note: the default Kalman filter adds 1e-8*I measurement-error
  variance as a positive-definiteness safeguard; see `?kalman_filter` for how to
  disable it for exact Dynare parity.)
- **Estimation** — random-walk Metropolis-Hastings, sequential Monte Carlo, and
  NUTS, with mode-finding (Nelder-Mead, CMA-ES, JADE) and prior tooling.
- **Occasionally-binding constraints** — OccBin/MCP/LCP and Boehl-style solvers,
  plus Ramsey/OSR/discretionary optimal-policy machinery.
- **Diagnostics** — an identification → convergence → fit → narrative battery
  that renders a Markdown report.

The solver and filter are validated against **Dynare 7.0** and **Dynare.jl**,
and the high-order sigma terms against a closed-form ground-truth model. The
parity test suite and its fixtures are maintained separately from this released
package.

## Install

```r
# install.packages("devtools")
devtools::install_github("Ball-Christopher/dynhr")

library(dynhr)
```

The package compiles a small amount of C++ (`src/`, via Rcpp + RcppArmadillo);
a C++ toolchain (Rtools on Windows, or Xcode command-line tools on macOS) is
required to install from source. Julia is **optional** and only needed for the
Dynare.jl interop.

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
- `inst/extdata/models/` — a few reference DSGE models used by the examples
- `inst/templates/` — report templates for the diagnostic battery

## A note on AI and reliability

AI tools were used extensively in the development of this package. The code has
been tested thoroughly throughout development, but it remains **experimental**
and may contain errors — **use at your own risk**, and validate results against
a trusted reference for any consequential use. `dynhr` is part of the author's
ongoing experimentation with AI-assisted development tools, and feedback and bug
reports are welcome.

## Licence

MIT — see [LICENSE](LICENSE).
