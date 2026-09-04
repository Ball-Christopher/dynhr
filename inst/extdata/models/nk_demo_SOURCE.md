# `nk_demo` — provenance and licence

`nk_demo.mod` and `nk_demo_data.csv` are the bundled worked example used by the
README quick start and by the estimation and diagnostics vignettes.

## Licence

**Both files were written for dynhr and are covered by the package's MIT
licence.** Neither is copied, adapted or derived from any third-party model
file, replication package or software distribution.

## The model

A three-equation New Keynesian model — IS curve, Phillips curve, Taylor rule —
augmented with two AR(1) driving processes and three measurement equations.

The economics is textbook. The canonical derivation is in

- Gali, J. (2015). *Monetary Policy, Inflation, and the Business Cycle*, 2nd
  ed., ch. 3. Princeton University Press.
- Woodford, M. (2003). *Interest and Prices*, ch. 4. Princeton University Press.

Those are cited as references **for the economics only**. The equations here
were written directly from that standard form in dynhr's own notation; no text
was taken from either book, from any software distribution's example suite, or
from any author's replication code. A linearised three-equation New Keynesian
model is common property of the literature — what a licence could attach to is
a particular *file*, and this file is ours.

Parameter values are round-number conventional calibrations (`beta = 0.99`,
`phi_pi = 1.5`, and so on) chosen for this demonstration, not taken from any
published estimate.

## The data

`nk_demo_data.csv` is **synthetic data simulated from `nk_demo.mod` itself** by
dynhr. It is not, and does not approximate, any real economic series.

It is reproducible from the shipped package:

```r
library(dynhr)
m  <- parse_mod(system.file("extdata/models/nk_demo.mod", package = "dynhr"))
cp <- compile_model(m, verbose = FALSE)
ss <- solve_steady(cp, m$param_values)
dr <- solve_perturbation(m, cp, ss$values, m$param_values)

set.seed(20260904L)
sim <- simulate_model(dr, n_periods = 600L, model = m)   # deviations from SS
obs <- c("ygr", "infl", "intr")
lev <- sweep(sim[, obs, drop = FALSE], 2, dr$ys[obs], "+")   # -> levels
Y   <- as.data.frame(lev[401:600, , drop = FALSE])           # 400 burn-in dropped
identical(round(Y, 6), read.csv(system.file("extdata/models/nk_demo_data.csv",
                                            package = "dynhr")))
```

200 quarterly observations of three observables, rounded to six decimals:

| column | meaning | steady state | sample mean | sample sd |
|---|---|---|---|---|
| `ygr` | quarterly output growth, % | 0.5 | 0.494 | 0.999 |
| `infl` | annualised inflation, % | 2.0 | 1.642 | 1.414 |
| `intr` | annualised nominal rate, % | 4.0 | 3.369 | 2.088 |

Sample means differ from the steady state by the usual amount for 200 draws of
a persistent process; that is a property of the draw, not an error.

## Why this file exists

The 0.9.3 cut needed a bundled estimable model whose provenance is beyond
question. The previous candidate for that role, `fs2000.mod`, could not be
used: this repository's own `inst/extdata/models/README.md` recorded the intent
to source it *"from the Dynare distribution"*, its variable declaration matches
Dynare's `fs2000.mod` exactly, and its 192-observation dataset matches Dynare's
`fsdat_simul`. Dynare is **GPL-3**, which cannot be redistributed inside an
MIT-licensed package. `nk_demo` removes the question rather than answering it.

## Scope of the claim

This is a **demonstration fixture**, not a replication target and not a
parity oracle. It exists so that a reader who has just installed the package
can run a complete parse -> compile -> steady state -> mode-finding -> sampling
pass against real-looking data in a few seconds. Because the data are simulated
from the model at known parameter values, posterior mass should sit near those
values — which makes it a convenient smoke test, but it says nothing about the
package's accuracy on real data. Dynare cross-checking lives in the
replication suite; the benchmark workload is `sw2007` (see
`sw2007_SOURCE.md`).
