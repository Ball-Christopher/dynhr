# Smets & Wouters (2007) — provenance and licence

These four files are the benchmark problem used by `dynhr_benchmark()`.

## Source

Smets, Frank and Rafael Wouters (2007), "Shocks and Frictions in US Business
Cycles: A Bayesian DSGE Approach", *American Economic Review* 97(3), 586–606.

Obtained from the AEA/ICPSR replication package **openicpsr-116269-V1**,
"Data and Code for: Shocks and Frictions in US Business Cycles".
Copyright 2007 American Economic Association.

## Licence

The replication package's own `LICENSE.txt` states:

- **Modified BSD Licence (BSD-3-Clause)** — applies to all code, scripts,
  programs and software.
- **Creative Commons Attribution 4.0 International (CC BY 4.0)** — applies to
  databases, images, tables, text and any other objects.

Both are permissive and redistributable with attribution, which is what this
file provides. `sw2007.mod` is code (BSD-3-Clause); `sw2007_data.csv`,
`sw2007_mode.csv` and `sw2007_hessian.csv` are data (CC BY 4.0).

NOTE for anyone re-checking this: the **CC BY-NC 3.0** notice one sees on the
openICPSR landing page applies to *ICPSR's own catalogue metadata records* (the
`*_ddi_2.5.xml` / `*_oai_pmh.xml` descriptors), NOT to the contents of the
deposit. The deposit's `LICENSE.txt` above governs these files. The
non-commercial clause is therefore not inherited, and nothing here conflicts
with dynhr's MIT licence.

## Files

| file | what | derivation |
|---|---|---|
| `sw2007.mod` | the model | `usmodel.mod`, byte-for-byte apart from line endings (the deposit's CRLF normalised to LF on checkout). No model content, parameter value or comment altered; the BSD notice requirement is met by this file plus the copyright line above |
| `sw2007_data.csv` | 230 quarterly observations of the 7 observables `dy dc dinve labobs pinfobs dw robs` | `usmodel_data.mat` (the file Dynare's `datafile=usmodel_data` actually reads), converted MAT -> CSV with no numerical transformation |
| `sw2007_mode.csv` | the published posterior mode, 36 values | `usmodel_mode.mat$xparam1`, in the order of the `.mod`'s `estimated_params` block |
| `sw2007_hessian.csv` | the mode Hessian, 36x36 | `usmodel_mode.mat$hh` |

The `.xls` copy of the data in the deposit is a working spreadsheet with
multiple header rows and intermediate columns; the `.mat` is the machine-read
series and is what was converted.

## Scope of the claim

`dynhr_benchmark()` uses this as a **fixed, realistic workload**. It is NOT a
Dynare-parity replication: the benchmark applies `first_obs = 71` but does not
reproduce Dynare's `presample = 4` or `lik_init = 2`, so the log posterior it
reports will not equal Dynare's. That number is recorded only as a
reproducibility fingerprint — two machines running the benchmark must agree on
it to the last digit, which is what makes their timings comparable. Actual
Dynare cross-checking of Smets-Wouters lives in the replication suite.
