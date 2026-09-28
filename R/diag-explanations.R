## R/diag-explanations.R
## --------------------------------------------------------------------------
## Methodology text for every diagnostic the suite can emit.
##
## Why this lives in the package and not in the Quarto templates
## ------------------------------------------------------------
## Until 0.9.4 this text was duplicated BY HAND in inst/templates/report.qmd
## and inst/templates/report-pdf.qmd, and nobody updated either one when a
## diagnostic changed. The 2026-09-17 report review found the two copies had
## drifted apart AND away from the code: the D5 entry still quoted the
## pre-0.9.4 flat "Bulk-ESS > 1000, Tail-ESS > 400" rule that the diagnostic's
## own roxygen calls a silent bug, D8 still said "sign and magnitude" when the
## shipped defaults are sign-only, D14 claimed SMC gives an EXACT log marginal
## likelihood (it gives a noisy unbiased estimate -- which is why the 2-SE gate
## exists), and 10 of the 17 diagnostics in a real posterior run had no
## explanation at all.
##
## One table, keyed exactly like `.DIAG_META`, serialised into report-meta.rds
## by `.report_meta()`. `test-report-094-explanations.R` asserts the two tables
## have identical key sets, so a new diagnostic cannot ship without text.
##
## Each entry: $what (what the diagnostic computes), $look_for (how to read
## the badge and the numbers), $citation.
## --------------------------------------------------------------------------

.DIAG_EXPLANATIONS <- list(

  # ---- Group A: pre-solve / structure -------------------------------------
  model_summary = list(
    what = paste0(
      "Model structure inventory: counts of endogenous variables, exogenous ",
      "shocks, equations and states, the Blanchard-Kahn condition at the ",
      "calibration point, and the dimensions of the policy functions ghx ",
      "(states) and ghu (shocks)."),
    look_for = paste0(
      "BK satisfied, and the state/shock counts matching what the .mod file ",
      "declares. A BK failure means there is no unique stable solution and ",
      "every downstream diagnostic is meaningless."),
    citation = "Blanchard & Kahn (1980) Econometrica."),

  data = list(
    what = paste0(
      "Observable data diagnostic: plots each observed series, reports its ",
      "mean, standard deviation, first-order autocorrelation and the count ",
      "of missing observations, and records the sample span."),
    look_for = paste0(
      "Series that are demeaned/detrended as the measurement equations ",
      "assume, no unintended level shifts, and a missing-value pattern the ",
      "filter can handle. A near-unit-root observable against a stationary ",
      "model is the classic source of a D9 failure."),
    citation = "Pfeifer (2013), 'A Guide to Specifying Observation Equations'."),

  expectations = list(
    what = paste0(
      "Checks the user-declared quantitative expectations in the ",
      "`@dynhr:expectations` block of the .mod file: data means and ratios, ",
      "parameter ranges, and IRF signs."),
    look_for = paste0(
      "All checks PASS. A failed `data_mean` or `data_ratio` check says the ",
      "data are inconsistent with what the model was written to reproduce; a ",
      "failed `param_range` check usually means an implausible mode."),
    citation = "User-declared; see ?diag_expectations."),

  d0 = list(
    what = paste0(
      "D0 static equation-system rank: the pre-flight check that the static ",
      "Jacobian has full rank at the steady state, i.e. that no equation is ",
      "redundant and no variable is left unpinned. Mirrors Dynare's ",
      "model_diagnostics.m."),
    look_for = paste0(
      "PASS = full rank. WARN = singular BUT the solved system has a unit ",
      "root, which makes the singularity expected (Dynare reports the same). ",
      "FAIL = a genuinely redundant or missing equation: read the listed ",
      "collinear equations and unpinned variables."),
    citation = "Dynare model_diagnostics.m; Adjemian et al. (2011) Dynare WP 1."),

  d19 = list(
    what = paste0(
      "D19 second-order solution quality: checks the order-2 perturbation ",
      "coefficients (ghxx/ghuu/ghxu and the shift correction) for finiteness ",
      "and magnitude before the solution is used for simulation."),
    look_for = paste0(
      "PASS = the second-order terms are finite and small relative to the ",
      "first-order ones. Large or non-finite terms mean the expansion point ",
      "is a poor approximation; consider pruning or a different ",
      "parameterisation."),
    citation = "Schmitt-Grohe & Uribe (2004) JEDC; Andreasen et al. (2018) REStud."),

  d40 = list(
    what = paste0(
      "D40 near-unit-root check: eigenvalues of the state transition block ",
      "ghx[state_idx, ], their half-lives log(0.5)/log(|lambda|), the ",
      "participation of each state in the dominant root, and central finite ",
      "differences d(spectral radius)/d(theta)."),
    look_for = paste0(
      "PASS = every half-life is short relative to the sample T. WARN = a ",
      "root is merely long-lived (half-life > halflife_frac x T): the ",
      "stationary P0 is very diffuse and moment-based diagnostics (D9, D10) ",
      "are unreliable. FAIL = a unit or explosive root: use ",
      "lik_init = 'auto'/'diffuse' and do not read unconditional moments at ",
      "all."),
    citation = "Durbin & Koopman (2012) ch. 5 (diffuse initialisation)."),

  # ---- Group B: identification & sensitivity ------------------------------
  d1 = list(
    what = paste0(
      "D1 local identification (Iskrev 2010): rank of the Jacobian of the ",
      "model-implied moments with respect to the estimated parameters, at ",
      "the calibration point."),
    look_for = paste0(
      "PASS = full column rank, so the parameters are locally identified by ",
      "that moment set. FAIL = rank deficient; the listed weak parameters lie ",
      "in the null space. Add moments, calibrate one of them, or ",
      "reparameterise."),
    citation = "Iskrev (2010) JME; Komunjer & Ng (2011) Econometrica."),

  d23 = list(
    what = paste0(
      "D23 spectral identification (Qu & Tkachenko 2012): rank of the ",
      "frequency-domain information matrix built from the spectral density ",
      "of the observables over the whole frequency band."),
    look_for = paste0(
      "PASS = full rank, so the parameters are identified from the second ",
      "moments of the observables at every frequency. FAIL = rank deficient; ",
      "add observables, restrict priors, or supply a dr-solve function as ",
      "`model_solve_fn` so the derivative is taken through the solution."),
    citation = "Qu & Tkachenko (2012) Quantitative Economics."),

  d24 = list(
    what = paste0(
      "D24 global identification via KL divergence (Qu & Tkachenko 2017): ",
      "profiles the Kullback-Leibler distance between the spectral densities ",
      "at theta and at perturbed parameter values, looking for distant ",
      "observationally-equivalent points that a local rank test cannot see."),
    look_for = paste0(
      "PASS = the KL profile rises away from theta in every direction. FAIL = ",
      "a flat or multi-modal profile: a distinct parameter vector implies the ",
      "same spectrum, so the posterior mode is not unique."),
    citation = "Qu & Tkachenko (2017) Journal of Econometrics."),

  d37 = list(
    what = paste0(
      "D37 dynamic identification rank check (Komunjer & Ng 2011): the rank ",
      "condition on the (A, B, C, D) state-space representation, which ",
      "handles minimality and observational equivalence up to a similarity ",
      "transform."),
    look_for = paste0(
      "PASS = the Delta matrix has full column rank, so theta is identified ",
      "from the transfer function. FAIL = rank deficient: the null space ",
      "names the unidentified parameter directions."),
    citation = "Komunjer & Ng (2011) Econometrica."),

  d20 = list(
    what = paste0(
      "D20 identification strength from the Fisher information matrix: ",
      "decomposes strength into a sensitivity part and a collinearity part ",
      "and reports a Cramer-Rao-based strength index s_i per parameter."),
    look_for = paste0(
      "Low |s_i| parameters are weakly identified even when the rank ",
      "condition passes -- prime candidates for calibration or ",
      "reparameterisation. A rank-deficient Fisher matrix makes every CRLB in ",
      "the table unreliable, not just the flagged rows."),
    citation = "Andrle (2010) IMF WP; Ratto & Iskrev (2011) Dynare WP 8."),

  d25 = list(
    what = paste0(
      "D25 higher-order identification: repeats the rank test using the ",
      "pruned second-order solution, so that information carried by third ",
      "moments and by the risk-adjustment terms enters the Jacobian."),
    look_for = paste0(
      "PASS = the order-2 moment set identifies parameters the order-1 set ",
      "does not. If it reports no gain, recheck the D1 baseline rank first -- ",
      "a deficient baseline makes the comparison meaningless."),
    citation = "Mutschler (2015) JEDC (identification at higher order)."),

  d3 = list(
    what = paste0(
      "D3 Morris elementary-effects screening: ranks parameters by mu*, the ",
      "mean absolute elementary effect, on each model moment. Because moments ",
      "carry different units the ranking uses the RELATIVE mu* -- each moment ",
      "column divided by its own largest entry, so 1 is the most influential ",
      "parameter for that moment."),
    look_for = paste0(
      "A parameter whose relative mu* is near zero for EVERY moment is ",
      "insensitive to the whole moment set and likely unidentified from it. A ",
      "near-zero ABSOLUTE mu* on one moment means nothing on its own -- it ",
      "may simply be a small-scale moment."),
    citation = "Morris (1991) Technometrics; Campolongo et al. (2007) EMS."),

  d22 = list(
    what = paste0(
      "D22 observable informativeness: maps local sensitivity contributions ",
      "from moments back to the observables that carry them, giving a ",
      "dominant-observable map per parameter."),
    look_for = paste0(
      "Use the map to decide which series identify which parameter, and where ",
      "two observables are redundant. A parameter with no dominant observable ",
      "is identified only through cross-equation restrictions."),
    citation = "Iskrev (2019) European Economic Review."),

  d26 = list(
    what = paste0(
      "D26 calibration sensitivity (Iskrev 2019b): perturbs each CALIBRATED ",
      "parameter over a grid and records the elasticity of the estimated ",
      "parameters, and of their precision, with respect to it. It also ranks ",
      "every split of the parameters into estimated and calibrated sets by the ",
      "Alegre Canton (2026) sensitivity statistic K: the worst-case ",
      "first-order bias of the estimates (or of a chosen object of interest) ",
      "per normalised calibration error, among the locally identified splits."),
    look_for = paste0(
      "WARN = an estimate moves more than one-for-one with a calibrated value ",
      "(|elasticity| > elasticity_tol). FAIL = precision collapses or ",
      "identification is lost somewhere on the grid. Either way: estimate the ",
      "listed calibrated parameters, or report results across calibrations. ",
      "The calibration-choice line (INFO, never changes the badge) names the ",
      "least-sensitive split and the current split's K for comparison: a ",
      "current K far above the minimum means a different set of parameters ",
      "should be fixed."),
    citation = paste0(
      "Iskrev (2019b), 'On the sensitivity of estimates to calibration'; ",
      "Alegre Canton (2026), 'Choosing What to Calibrate and What to Estimate ",
      "in Structural Models', arXiv:2606.25688.")),

  d4 = list(
    what = paste0(
      "D4 prior predictive check: draws from the joint prior, solves the ",
      "model at each draw and compares the resulting distribution of model ",
      "moments to the observed data moments, as a prior-predictive p-value."),
    look_for = paste0(
      "Every moment's PPP inside [0.025, 0.975]. A p-value at either extreme ",
      "says the prior itself is inconsistent with the data, before a single ",
      "posterior draw is taken -- tighten or widen the relevant priors."),
    citation = "Gelman, Meng & Stern (1996) Statistica Sinica; Fernandez-Villaverde & Guerron-Quintana (2020) JEDC."),

  d30 = list(
    what = paste0(
      "D30 arbitrary-precision rank checks: recomputes the identification ",
      "singular values in extended precision (Rmpfr) so that a rank verdict ",
      "sitting near the double-precision noise floor can be resolved."),
    look_for = paste0(
      "Inspect `sv_comparison`: a singular value that moves materially ",
      "between double and extended precision was never a reliable rank ",
      "verdict. Values stable across precisions can be trusted."),
    citation = "Higham (2002) Accuracy and Stability of Numerical Algorithms."),

  d27 = list(
    what = paste0(
      "D27 OBC / piecewise-linear identification: re-runs the identification ",
      "check in the presence of an occasionally-binding constraint, where the ",
      "regime indicator itself carries information about the parameters."),
    look_for = paste0(
      "PASS = the piecewise parameters are identified given the observed ",
      "binding frequency. FAIL = identification depends on episodes the ",
      "sample barely contains."),
    citation = "Guerrieri & Iacoviello (2015) JME."),

  d28 = list(
    what = paste0(
      "D28 regime-switching identification: checks identification of the ",
      "regime-dependent parameters and of the transition matrix, given how ",
      "often each regime is visited in the sample."),
    look_for = paste0(
      "PASS = every regime is visited often enough to identify its own ",
      "parameters. FAIL = a rarely-visited regime whose parameters are driven ",
      "by the prior."),
    citation = "Hamilton (1989) Econometrica; Farmer, Waggoner & Zha (2011) JoE."),

  d33 = list(
    what = paste0(
      "D33 structural-vs-reduced-form gap ('borrowed identification'): asks ",
      "whether a deep parameter is identified in its own right or only ",
      "through a reduced-form combination declared via `reduced_form=`."),
    look_for = paste0(
      "PASS = the deep parameters have their own identification. FAIL = a ",
      "structural interpretation is being borrowed from a reduced-form ",
      "coefficient; declare the map, or stop reading the parameter ",
      "structurally."),
    citation = "Lucas (1976) Carnegie-Rochester; Nakamura & Steinsson (2018) JEP."),

  # ---- Group C: estimation / MCMC -----------------------------------------
  d5 = list(
    what = paste0(
      "D5 MCMC convergence: rank-normalised split R-hat, Bulk-ESS and ",
      "Tail-ESS (Vehtari et al. 2021), trace and rank-overlay plots, the ",
      "lag-1 ACF, and -- for NUTS -- BFMI, divergences and treedepth ",
      "saturation."),
    look_for = paste0(
      "FAIL when any R-hat is >= 1.05 or is not computable (a stuck or ",
      "constant chain). WARN when an R-hat is in the 1.01-1.05 band, or when ",
      "Bulk- or Tail-ESS falls below the target. The ESS target is Vehtari et ",
      "al.'s own scaling rule, 100 x n_chains (a single chain is split, so ",
      "the floor is 200) -- NOT the flat, chain-count-independent ESS floor ",
      "dynhr used before 0.9.4, which silently assumed ten chains and ",
      "appears nowhere in the literature. There is no ESS FAIL tier: a low ",
      "ESS means the Monte Carlo error is larger than you want, not that the ",
      "draws come from the wrong distribution."),
    citation = "Vehtari, Gelman, Simpson, Carpenter & Burkner (2021) Bayesian Analysis."),

  d7 = list(
    what = paste0(
      "D7 mode-finding robustness: clusters the end points of every ",
      "optimisation start into basins and reports whether they agree on one ",
      "posterior mode."),
    look_for = paste0(
      "One basin containing all starts = robust. Several basins = the ",
      "posterior may be multimodal; re-run from more starting points, or use ",
      "an SMC sampler rather than trusting a single mode."),
    citation = "Gelman & Rubin (1992) Statistical Science."),

  d6 = list(
    what = paste0(
      "D6 prior-vs-posterior updating: overlays the prior and posterior ",
      "density of every estimated parameter and reports the overlap ",
      "coefficient, the fraction of shared density area."),
    look_for = paste0(
      "Overlap well below 0.80 = the data are informative about that ",
      "parameter. Overlap near 1 = the posterior is the prior; the parameter ",
      "is not identified by this sample, whatever its credible interval ",
      "looks like."),
    citation = "Geweke (1992) Bayesian Statistics 4; Koop (2003) Bayesian Econometrics."),

  prior_sensitivity = list(
    what = paste0(
      "Prior sensitivity: re-estimates (or re-weights) under flat/uniform ",
      "priors and compares the resulting posterior locations with the ",
      "informative-prior run. Each flat prior is uniform over the ",
      "informative prior's own bounds when both are finite, and otherwise ",
      "over its 0.5%-99.5% quantile range intersected with those bounds; ",
      "the support used per parameter is reported alongside the shifts. The ",
      "flat-prior search STARTS AT THE INFORMATIVE MODE, not at the prior ",
      "means (for a uniform prior, the midpoint of its support): the question ",
      "is local -- does the mode move when the prior is flattened? -- and a ",
      "far-away start would turn it into a global search whose failure to ",
      "converge would be misreported as prior-drivenness."),
    look_for = paste0(
      "Estimates that barely move are data-driven. Estimates that move a long ",
      "way are prior-driven and should be reported as such, with the flat-",
      "prior result alongside. An INFO badge means the flat search never beat ",
      "its own starting log-posterior, so it did not converge and the ",
      "comparison says nothing about the priors."),
    citation = "Muller (2012) Journal of Monetary Economics."),

  d21 = list(
    what = paste0(
      "D21 Bayesian KPS precision updating: regresses posterior precision on ",
      "sample size across nested subsamples and reports the slope per ",
      "parameter."),
    look_for = paste0(
      "A slope near 1 (precision growing with T) implies strong ",
      "identification. A slope near 0 means adding data does not sharpen the ",
      "posterior -- weak or non-identification."),
    citation = "Koop, Pesaran & Smith (2013) JBES."),

  # ---- Group D: post-estimation -------------------------------------------
  d8 = list(
    what = paste0(
      "D8 IRF plausibility: checks impulse responses against benchmarks, ",
      "either the user's `@dynhr:benchmarks` block or the shipped defaults ",
      "(monetary, technology and UIP shocks on output, inflation and the real ",
      "exchange rate)."),
    look_for = paste0(
      "The shipped defaults check the response SIGN ONLY, over horizons 1-8 ",
      "(`sign_window = c(1, 8)`); peak magnitude and peak timing windows are ",
      "OPT-IN (`peak_range` / `peak_horizon`, both NULL by default). The peak ",
      "value and horizon printed for each benchmark are reported OUTPUT, not ",
      "a restriction that was tested. SKIP = no shock in the model matched ",
      "the benchmark's shock-name list. A failed sign check usually points at ",
      "the policy rule or the shock calibration."),
    citation = "Christiano, Eichenbaum & Evans (1999) Handbook of Macroeconomics."),

  d9 = list(
    what = paste0(
      "D9 moment matching: compares each observable's model-implied standard ",
      "deviation (at the posterior mode, INCLUDING measurement error -- ",
      "`model_sd`; `model_sd_state` is the state-only figure) with its ",
      "empirical standard deviation, and optionally reports a posterior-",
      "predictive band."),
    look_for = paste0(
      "Each row carries a status: ok (ratio inside [0.5, 2.0]), marginal ",
      "(0.5-0.6 or 1.7-2.0), outside, or non-finite. The badge is a ",
      "70 %-rule: PASS when at least 70 % of the FINITE ratios are inside the ",
      "band, so a single outside row is not a failure. Non-finite ratios are ",
      "excluded from the denominator, not counted as failures. The posterior-",
      "predictive band, when computed, is reported as INFO and does not gate."),
    citation = "Schorfheide (2000) JAE; Smets & Wouters (2007) AER."),

  d10 = list(
    what = paste0(
      "D10 variance decomposition: the share of each observable's forecast-",
      "error variance attributable to each structural shock, unconditionally ",
      "and by horizon."),
    look_for = paste0(
      "Informational. One shock accounting for nearly all of a series' ",
      "variance usually signals misspecification or a missing propagation ",
      "mechanism. For open economies, foreign shocks are typically 20-40 % of ",
      "output variance."),
    citation = "Sims (1980) Econometrica; Kamber, Morley & Wong (2018) REStat."),

  d11 = list(
    what = paste0(
      "D11 historical decomposition: splits the observed path of each series ",
      "into contributions from each structural shock plus the initial-state ",
      "effect, using the Kalman smoother."),
    look_for = paste0(
      "Informational. Check that recognisable historical episodes are ",
      "attributed to plausible shocks, and that the contributions add up to ",
      "the observed series (a residual that does not close means the ",
      "decomposition is being read outside its sample)."),
    citation = "Banbura, Giannone & Reichlin (2010) JAE."),

  d12 = list(
    what = paste0(
      "D12 smoothed-shock diagnostics: tests whether the smoothed structural ",
      "shocks are serially uncorrelated (Ljung-Box), zero-mean and ",
      "approximately normal, and reports their per-period standard deviation ",
      "against the model's."),
    look_for = paste0(
      "Ljung-Box p > 0.05 for every shock. Serial correlation or a non-zero ",
      "mean says the model is missing a persistent component and the filter ",
      "is absorbing it into the shocks."),
    citation = "Ljung & Box (1978) Biometrika."),

  d41 = list(
    what = paste0(
      "D41 Kalman-filter innovation whiteness: tests the one-step-ahead ",
      "prediction errors (standardised by their filtered variance) for serial ",
      "correlation, for a zero mean, and for unit variance -- the properties ",
      "the Gaussian likelihood assumes."),
    look_for = paste0(
      "PASS = white, zero-mean, unit-variance innovations, so the likelihood ",
      "is correctly specified for this data. Serial correlation in the ",
      "innovations means the state space is missing dynamics; a variance far ",
      "from 1 means the measurement-error or shock scaling is wrong. Unlike ",
      "D12 this reads the FILTER's own residuals, so it is valid even where ",
      "unconditional moments are not."),
    citation = "Harvey (1989) Forecasting, Structural Time Series and the Kalman Filter, ch. 5."),

  d13 = list(
    what = paste0(
      "D13 cross-equation restrictions (Del Negro & Schorfheide 2004): ",
      "compares the DSGE-implied reduced-form VAR, obtained from the state ",
      "space by Yule-Walker, with an unrestricted OLS VAR on the same data."),
    look_for = paste0(
      "PASS = small Frobenius relative difference (< 0.5) and a Wald p-value ",
      "above 0.05, i.e. the restrictions are not rejected. A rejection says ",
      "the model's cross-equation implications, not just its parameters, are ",
      "at odds with the data."),
    citation = "Del Negro & Schorfheide (2004, 2009) JBES."),

  d15 = list(
    what = paste0(
      "D15 DSGE-VAR: estimates a VAR whose prior is centred on the DSGE ",
      "restrictions with tightness lambda, and picks the lambda* that ",
      "maximises the marginal likelihood over the grid."),
    look_for = paste0(
      "lambda = Inf is a first-class grid point and the top of the default ",
      "grid: lambda* = Inf means the cross-equation restrictions bind exactly ",
      "-- the DSGE-implied VAR beats every loosened alternative, the best ",
      "available outcome. lambda* > 2 = restrictions well supported. ",
      "lambda* < 0.5 = the data want a much looser VAR, i.e. the DSGE ",
      "restrictions are rejected."),
    citation = "Del Negro & Schorfheide (2006) JBES."),

  d29 = list(
    what = paste0(
      "D29 data constraints: sampling-weighted minimum-distance t-ratios ",
      "(|theta| / SE) showing how tightly the data's variance and ",
      "autocovariance moments pin each parameter, plus the Stock-Wright S ",
      "statistic."),
    look_for = paste0(
      "Small t-ratios name the parameters the moments barely constrain: add ",
      "informative observables or calibrate them. D20 with ",
      "weighting = \"sampling\" is the same information in rank form. The ",
      "Stock-Wright S statistic is REPORTED, not gated."),
    citation = "Stock & Wright (2000) Econometrica; Christiano, Trabandt & Walentin (2010) JEDC."),

  d16 = list(
    what = paste0(
      "D16 subsample stability: re-estimates on sample splits and flags ",
      "parameters whose credible intervals do not overlap across the splits."),
    look_for = paste0(
      "Overlapping intervals = stable structure. Non-overlap flags a ",
      "structural break or genuinely time-varying behaviour, and makes the ",
      "full-sample posterior an average of two regimes."),
    citation = "Lubik & Schorfheide (2004) AER; Fernandez-Villaverde & Rubio-Ramirez (2008) REStud."),

  d14 = list(
    what = paste0(
      "D14 Bayes factor comparison: log marginal likelihoods across the ",
      "supplied models, the natural-log and log10 Bayes factors against the ",
      "best model, and the prior/posterior model probabilities."),
    look_for = paste0(
      "Kass-Raftery / Jeffreys bands on log10 BF: > 2 decisive, 1-2 strong, ",
      "0.5-1 substantial. PASS requires at least 'strong' evidence against ",
      "EVERY competitor, so a tie for best is a failure. Marginal likelihoods ",
      "are themselves Monte Carlo estimates -- SMC and bridge sampling give a ",
      "noisy unbiased estimate, not an exact value -- so when ",
      "`log_marglik_se` is supplied and the best-vs-runner-up gap falls ",
      "inside `se_multiple` (default 2) combined standard errors, the badge ",
      "is downgraded to INFO: there is no winner to report. That noise check ",
      "is layered on top of the Kass-Raftery scale, not a replacement for it."),
    citation = "Kass & Raftery (1995) JASA; Jeffreys (1961) Theory of Probability."),

  model_comparison = list(
    what = paste0(
      "Multi-model comparison table: side-by-side fit, log marginal ",
      "likelihood and parameter estimates for two or more estimated dynhr ",
      "models."),
    look_for = paste0(
      "Informational. Read the relative log Bayes factors with the same ",
      "simulation-noise caution as D14, and check the models were estimated ",
      "on identical data before comparing evidences at all."),
    citation = "Geweke (2005) Contemporary Bayesian Econometrics and Statistics."),

  posterior_predictive = list(
    what = paste0(
      "Posterior predictive checks: draws from the posterior, solves and ",
      "simulates at each draw, and compares the distribution of model-implied ",
      "moments with the data moments via posterior-predictive p-values."),
    look_for = paste0(
      "PPP near 0.5 = the data moment is typical of what the fitted model ",
      "generates. PPP < 0.05 = the model systematically under-produces that ",
      "moment; > 0.95 = it over-produces. Unlike a prior-predictive failure, ",
      "this is a statement about the fitted model, not the prior."),
    citation = "Gelman, Meng & Stern (1996); Herbst & Schorfheide (2016) Bayesian Estimation of DSGE Models."),

  bayesian_irf = list(
    what = paste0(
      "Bayesian IRFs: impulse responses recomputed at every posterior draw, ",
      "summarised as pointwise posterior medians with credible bands."),
    look_for = paste0(
      "Informational. A band that straddles zero over the whole horizon means ",
      "the sign of that response is not identified by the data, however clean ",
      "the mode IRF looks."),
    citation = "Sims & Zha (1999) Econometrica."),

  d17 = list(
    what = paste0(
      "D17 narrative identification: checks that the smoothed shocks satisfy ",
      "user-declared sign and magnitude restrictions for named historical ",
      "episodes from the `@dynhr:narratives` block."),
    look_for = paste0(
      "PASS = the shock story is consistent with the narrative. WEAK = the ",
      "shock has the right sign but is below the declared threshold. FAIL = ",
      "the smoothed shock contradicts the narrative, which usually means the ",
      "shock is soaking up something else."),
    citation = "Antolin-Diaz & Rubio-Ramirez (2018) AER."),

  shock_dominant = list(
    what = paste0(
      "Dominant-shock summary: for each observable, the shock with the ",
      "largest variance share and the size of that share."),
    look_for = paste0(
      "Informational. A single shock dominating every observable is the ",
      "classic sign of a model with one effective propagation channel."),
    citation = "Sims (1980) Econometrica."),

  d18 = list(
    what = paste0(
      "D18 welfare plausibility: compares the welfare implied by the ",
      "estimated model with the Ramsey / optimal-policy benchmark and reports ",
      "the gap as a percentage of steady-state welfare."),
    look_for = paste0(
      "FAIL when |gap / steady_welfare| x 100 exceeds `gap_threshold`. A very ",
      "large gap usually means the Ramsey problem is mis-specified (check the ",
      "discounting of the FOCs) rather than that policy is wildly ",
      "suboptimal."),
    citation = "Lucas (1987) Models of Business Cycles; Schmitt-Grohe & Uribe (2007) JME."),

  d31 = list(
    what = paste0(
      "D31 OBC scenario comparison: solves the constrained and unconstrained ",
      "paths and reports where the bound binds and by how much the paths ",
      "differ."),
    look_for = paste0(
      "FAIL = the constrained path violates a bound, or the complementarity ",
      "solve did not converge. Bounds are compared in LEVELS (steady state + ",
      "deviation): check the units of `bound`, and that a bound given as a ",
      "parameter name carries the value you intended."),
    citation = "Guerrieri & Iacoviello (2015) JME (OccBin)."),

  d32 = list(
    what = paste0(
      "D32 OBC binding periods: counts and dates the periods in which the ",
      "occasionally-binding constraint is active along the smoothed or ",
      "simulated path."),
    look_for = paste0(
      "Informational. A constraint that binds in almost every period, or in ",
      "none, means the OBC machinery is not doing what the model intends -- ",
      "and in the 'almost never' case the OBC parameters are unidentified."),
    citation = "Guerrieri & Iacoviello (2015) JME."),

  d34 = list(
    what = paste0(
      "D34 policy-partitioned invariance (operational Lucas critique): ",
      "re-estimates the PRIVATE-block parameters under perturbed policy ",
      "parameters and asks whether they stay put, as a structural ",
      "interpretation requires them to."),
    look_for = paste0(
      "PASS = the private-block estimates are invariant to the policy block. ",
      "FAIL = they move with policy, so they are not policy-invariant deep ",
      "parameters and counterfactual policy experiments are not valid."),
    citation = "Lucas (1976) Carnegie-Rochester Conference Series."),

  d35 = list(
    what = paste0(
      "D35 misspecification softness: compares sandwich standard errors with ",
      "Hessian-only standard errors (the information-matrix equality) to see ",
      "how much of the reported precision survives misspecification."),
    look_for = paste0(
      "Sandwich and Hessian SEs close together = the IM equality roughly ",
      "holds. A sandwich SE much larger than the Hessian SE means the ",
      "reported credible intervals are too narrow -- the parameter is 'soft'."),
    citation = "White (1982) Econometrica; Muller (2013) Econometrica."),

  d36 = list(
    what = paste0(
      "D36 calibration deepness: tests whether the CALIBRATED parameters are ",
      "in tension with the data, by checking how much the fit would improve ",
      "if each were freed."),
    look_for = paste0(
      "FAIL = the data want a calibrated value materially different from the ",
      "one imposed; estimate that parameter or justify the calibration ",
      "against the same data."),
    citation = "Iskrev (2019b); Canova & Sala (2009) JME."),

  deep_passport = list(
    what = paste0(
      "Deep-parameter passport: synthesises D33-D36 into a per-parameter ",
      "grade covering own identification, policy invariance, softness under ",
      "misspecification, and calibration tension."),
    look_for = paste0(
      "Read it as the answer to 'which of these parameters may I interpret ",
      "structurally, and which are reduced-form coefficients wearing a ",
      "structural name?'. A low grade does not invalidate the fit; it ",
      "invalidates the counterfactual."),
    citation = "Lucas (1976); Nakamura & Steinsson (2018) JEP.")
)


#' Explanation text for a diagnostic, with a generic fallback
#'
#' @param nm Diagnostic name as used in the results list.
#' @return List with `what`, `look_for`, `citation`.
#' @noRd
.diag_explanation_for <- function(nm) {
  e <- .DIAG_EXPLANATIONS[[nm]]
  if (!is.null(e)) return(e)
  list(
    what = paste0(
      "No methodology note is registered for '", nm, "' in ",
      "R/diag-explanations.R."),
    look_for = paste0(
      "Read the diagnostic's own summary text and its help page (?", nm,
      "). The badge still follows the package-wide convention: PASS / WARN ",
      "(soft failure) / FAIL / ERROR (did not run) / INFO (reported, not ",
      "gated)."),
    citation = "")
}
