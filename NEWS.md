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
