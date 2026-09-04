# Bundled models — provenance and licence

Every `.mod` file and dataset shipped in the released package is listed here.
Nothing ships without a row. `sw2007` and `nk_demo` have their own detailed
files (`sw2007_SOURCE.md`, `nk_demo_SOURCE.md`); the rest are covered below.

| file | origin | licence |
|---|---|---|
| `nk_demo.mod`, `nk_demo_data.csv` | written for dynhr; data simulated from the model | MIT (this package) — see `nk_demo_SOURCE.md` |
| `rbc.mod` | written for dynhr | MIT (this package) |
| `rbc2shock.mod` | written for dynhr | MIT (this package) |
| `nk_2obc.mod` | written for dynhr | MIT (this package) |
| `nk_zlb_dynare.mod` | written for dynhr; cleared at the 0.8.1 public cut — see the note below | MIT (this package) |
| `sw2007*` | AEA/ICPSR replication package openicpsr-116269-V1 | BSD-3-Clause (code) / CC BY 4.0 (data) — see `sw2007_SOURCE.md` |

## The self-authored models

`rbc.mod` and `rbc2shock.mod` are the standard King–Plosser–Rebelo real
business cycle model — Cobb–Douglas production, CRRA utility, an AR(1) TFP
process, and in the two-shock variant an AR(1) preference shock. `nk_2obc.mod`
is a minimal linearised three-equation New Keynesian model carrying two
simultaneously-bindable occasionally-binding constraints, written as a test
fixture for the multi-constraint regime logic.

These were written in dynhr's own notation from the standard published form of
each model. The economics is common property of the macro literature
(King, Plosser & Rebelo 1988; Galí 2015 for the New Keynesian block); what a
licence could attach to is a particular *file*, and these files are ours.

## Note on `nk_zlb_dynare.mod`

The name records what the file is *for* — cross-checking dynhr's
occasionally-binding-constraint solvers against Dynare's OccBin implementation,
so it is written in the OccBin `mcp` tag syntax Dynare accepts. It is a
three-equation linearised New Keynesian model with a zero lower bound on the
nominal rate: the canonical worked example of that literature, and the same
three equations as `nk_2obc.mod` with one constraint instead of two.

Its header cites Guerrieri & Iacoviello (2015) and Giovannini, Pfeiffer &
Ratto (2021) as the **method** references for OccBin, not as sources of text.

**Cleared as part of the 0.8.1 public release review (maintainer, 2026-09-04).**
This file was one of exactly three models in the first curated public
artifact — the 0.8.1 cut `d8f3b00` shipped `nk_2obc.mod`,
`nk_zlb_dynare.mod` and `rbc.mod` — and has shipped in every public release
since (0.8.1, 0.9.0, 0.9.1, 0.9.2). Its provenance was settled at that cut;
this record simply writes down what was previously only implicit.

Consistent with that: it carries dynhr's own `@dynhr-model` metadata header
and is the same three equations as the self-authored `nk_2obc.mod` with one
constraint instead of two. Note for anyone re-checking that syntactic
similarity to a Dynare test model is expected and is not evidence of copying
— the `mcp` tag syntax is fixed by Dynare's parser, and a three-equation NK
model with a ZLB has very little room for expressive variation.

## For anyone adding a model

Add a row here before shipping it, and do not ship a file whose provenance you
cannot state. In particular, **do not take `.mod` files from the Dynare
distribution, its test suite, or its example directory**: Dynare is GPL-3 and
incompatible with this package's MIT licence. This was not hypothetical — the
0.9.3 cut nearly shipped `fs2000.mod` and its dataset, which this repository's
own development notes record as intended to be sourced "from the Dynare
distribution", whose variable declaration matches Dynare's `fs2000.mod` exactly
and whose 192-observation dataset matches Dynare's `fsdat_simul`. It was
replaced by `nk_demo`, written from scratch with data simulated from the model
itself.

Published *equations* are not licensable — transcribing a model from a journal
article or textbook into your own file is fine, and is what the self-authored
models above are. Copying somebody's *file* is not.
