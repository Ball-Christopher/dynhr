# dynhr 0.9.3.6

## A model with no shock variance is now called out at the source

Reported as a historical-decomposition add-up failure: a two-state linear model
with cross-loaded exact observations (`y1 = x`, `y2 = x + z`) and
`lik_init = "kappa"` gave `adding_up_residual` 7.9e-3, while the independent,
duplicate and single-observable variants balanced to 1e-17. The natural reading
-- that the cross-loaded geometry breaks the singular-innovation path -- is not
what was happening.

**The reproducer's model has no `shocks;` block, so `Sigma_e` is entirely
zero.** The model has no stochastic structure at all. With `Q = 0` the whole
path is determined by `s_0`, and under a kappa prior `s_{0|T}` is computed as
`P_{0|0} r_0` with `P_{0|0} = 1e6 * I` -- a large number times a small one --
so kappa's round-off lands directly in the initial state and `s_1` stops
equalling `T s_0`. The cross-loaded geometry only decides whether that
round-off is *visible*; the other three variants are degenerate in ways that
happen to hide it.

Give the same model a `shocks;` block and it balances to 1.1e-10 under kappa
and 1.7e-16 on the default initialisation, dropping nothing at all -- the
singular-innovation path is not even entered.

`kalman_filter()` and `kalman_smoother()` now warn when EVERY shock has zero
variance, naming the missing `shocks;` block as the usual cause. A per-shock
`stderr 0` stays legitimate and silent -- it is what makes a deterministic
`known_shocks` or `shock_means` injection meaningful; all of them being zero is
the different thing.

## `adding_up_residual` is necessary but not sufficient, and now says so

The decomposition propagates its components with `T_mat`/`R_mat` while the
add-up check rebuilds the path from the smoother's states separately -- so what
can that comparison actually see? Measured, by corrupting the smoother's output
so the answer is known by construction:

| corruption (nk_demo)               | `adding_up` | `transition` |
|------------------------------------|-------------|--------------|
| shock at t = 1                     | 2.56        | 2.57         |
| shock mid-sample                   | 2.56        | 2.57         |
| state mid-sample                   | 1.00        | 1.00         |
| shock at the LAST period           | **5.9e-15** | 2.57         |
| all shocks x1.5 in the LAST period | **5.9e-15** | 0.82         |
| state at the LAST period           | **5.9e-15** | 1.00         |

The contemporaneous `ghu %*% eps_t` term appears identically on both sides of
the adding-up comparison and CANCELS, so a shock error is visible only through
its propagated (t+1 onward) effect, damped by T -- and the final period is not
checked at all.

`historical_decomposition()` therefore also returns `$transition_residual`,
`$transition_residual_by_period`, `$transition_worst_period` and
`$transition_ok`: the direct question, does the smoother's own output satisfy
its own transition, `s_t = T s_{t-1} + R eps_t`, period by period. It has
neither blind spot and catches all six corruptions. **Look at it first** when a
decomposition will not add up -- it separates "the decomposition is wrong" from
"its INPUTS are incoherent", and the per-period vector localises the latter. On
the reported case it puts the entire error in period 1, which is what
identified the initial state as the culprit.

This also answers the report's second acceptance criterion directly: the
decomposition no longer returns a non-additive result as if it were valid --
it warns, says which of the two checks failed, and where.

Nothing about the computed contributions changes.

## Cross-loaded exact observations, pinned as coherent

Three observables each loading both states, none with measurement error, the
third an exact combination of the other two -- a predictable component dropped
in every period (160 of 160), and with ragged edges the dropped SET varying
period to period (130 of 160). Transition residual 6.7e-16 and 1.3e-15, with
the smoothed states exact against the simulated truth. Both are regression
tests now, as is the reporter's own four-variant reproducer.

# dynhr 0.9.3.5

**A successful Cholesky is not a test that the innovation covariance is
invertible, and every multivariate path was using it as one.**

Reported from a shock- and historical-decomposition matching exercise:
`me_variance = 0` gave adding-up residuals of ~2e6 and `me_variance = 1e-12`
about 4. The decomposition was reporting the problem faithfully -- its
residual tracks the smoother's own transition residual, and the smoother's
states and shocks had stopped being consistent with each other.

**The defect.** `chol()` can factorise a matrix that is singular to round-off
and return a garbage pivot. On a stochastically singular system -- more
observables than shocks, or an observable that is an exact combination of
others -- that produces a finite, badly wrong answer instead of a detected
failure. Measured on a two-observable / one-shock fixture with
`me_variance = 0`: the innovation covariance had `rcond` 1.8e-17 and a
NEGATIVE determinant, `chol()` succeeded, and

```
kalman_filter()   +292.7      <- garbage, and too HIGH
kalman_smoother()  -10.8      <- dropped the component in 19 of 20 periods
univariate filter  -26.8      <- correct
```

The smoother's error sat entirely in the ONE period where `chol()` happened to
succeed: its smoothed states were exact while its smoothed shocks were
inconsistent with them by 0.18. A likelihood that is too high is the dangerous
direction -- an optimiser walks straight into it.

**The fix**, in the three places that inverted F: after a successful `chol()`,
test the pivots. The i-th squared diagonal of the Cholesky factor IS that
observable's conditional variance -- the same quantity the univariate filter
has always skipped on, and the one `.smoother_informative_obs()` already used
to pick the informative subset. It was simply never consulted when `chol()`
succeeded. Sites: `.kf_step()` (the `standard` and `dare` paths), the C++
`kalman_standard_loop_cpp` fast path, and the smoother's forward pass.

**After the fix, every path agrees.** All five filter methods and the smoother
return -26.828497 on that fixture, each multivariate path detecting the
singularity and rerouting to the univariate filter, which drops the
uninformative component instead of inverting through it. The decomposition
residual goes from 2.9e-1 to 8.9e-16 at `me_variance = 0`, and from 3.4e-4 to
2.7e-15 at 1e-12.

**Where the boundary sits.** A component whose conditional variance is below
`kalman_tol` (default 1e-10, relative to F's scale) is dropped and the answer
is exact. Above it the component genuinely carries information and the
accuracy is the conditioning limit of inverting F, about `eps / me_variance` --
a numerical fact rather than a defect, and one that
`historical_decomposition()`'s adding-up check now surfaces rather than
leaving to be discovered downstream.

**One test changed meaning.** `test-kalman-smoother-na.R`'s "loglik decreases
when data are removed" used a fixture that is an AR(1) plus an ALIAS -- two
observables, one shock -- and blanked one of them. It passed only because the
smoother was inverting through the singular F, so the "information" being
removed was round-off. With an exact alias, removing either series alone costs
nothing, because the other still pins the state; only removing both does. The
test now asserts that, which is sharper than what it replaced: the two
one-column runs agree with the full run to the bit and with each other.

# dynhr 0.9.3.4

A fourth report against the filtering surface, on historical decomposition.
The decomposition function was correct; its integration contract was neither
documented nor checked, and one input combination was genuinely wrong.

**Bug fix: `pre_sample` put the decomposition out of phase.** With a backfill,
`kalman_smoother()` returns its series trimmed to the caller's sample, but
`smoothed_initial` is still `s_{0|T}` for the PADDED sample -- the state `k`
periods earlier. `smoothed_initial_state()` handed that to
`historical_decomposition()` as the anchor, so the initial-condition
trajectory ran `k` periods out of phase and the components stopped adding up
(measured 3.8e-2 on the two-shock fixture at `k = 2`), with every dimension
still correct. It now returns the last pre-sample row, which is the period
immediately before the returned rows.

**The contract is now enforced, not just described.** A correctly-SIZED but
wrongly-ORDERED `s0`, or shock columns in a different order from
`ss$shock_names`, used to be accepted in silence. Since the adding-up residual
is exactly the size of the initial-condition error, on a model with large
states -- or a kappa-initialised smoother, whose `s_{0|T}` carries the
arbitrary prior -- that silence surfaces as a residual of 1e8-1e9 with nothing
to say why. Now: columns and named vectors are matched to the state space by
NAME and reordered, a mismatch is an error naming both sets, and a transposed
matrix is refused with the orientation it wanted.

**`adding_up_residual` is always present.** It used to appear only when
`smoothed_states` was supplied -- so the one call shape that can be silently
wrong was also the one with no diagnostic at all. It is now a number when
there is something to check against and `NA` with `$adding_up_note` when there
is not, alongside `$adding_up_relative`, `$adding_up_ok` and a `tol` argument
(default 1e-8, relative to the path's scale). Exceeding it warns and names the
likely causes.

**Documented contract**, in `?historical_decomposition`: orientation (rows are
periods; `kalman_filter()`'s state matrices are the other way round and need
`t()`), the state at the first contribution period (row 1 loads the PRE-SAMPLE
state, which is why the `initial` column exists), pre- versus post-transition
(pre-: the state entering the period plus that period's shock), and how
`shock_means` / `known_shocks` enter (through the smoother's output, in their
own shock's column, needing nothing here).

**Verified against IRIS `simulate(..., 'contributions', true)`** to 2e-16 on a
two-shock linear model: each isolated-shock column, the initial-condition
column, and the total. Two IRIS columns have no dynhr counterpart by
construction, and the docs now say so: a measurement-shock column (dynhr's
`me_variance` is observation noise, not a structural shock, so compare against
IRIS's structural columns plus its init column rather than its grand total),
and a nonlinear column that is identically zero for a linear model.
`$has_nonlinear_column` and `$has_residual_column` report this.

# dynhr 0.9.3.3

A third report against the filtering surface, this one a feature request
rather than a defect: an IRIS-compatible deterministic shock-mean path, and an
explicit state-timing contract on the result.

**Verified against IRIS itself.** With `shock_means` and the default
`shock_timing = "dated"`, dynhr reproduces IRIS's
`filter(m, d, range, 'vary=', j)` to machine precision -- on all three of its
outputs and on the smoothed shocks -- with **no adapter shift of any kind**
(checked against IRIS Toolbox Release 20180308 under Octave; the reference
values are transcribed into `test-shock-means.R` at full precision, since IRIS
is not a dependency). The names line up one-for-one: IRIS `'predict'` /
`'filter'` / `'smooth'` are dynhr `predicted_states` / `updated_states` /
`smoothed_states`.

If you previously needed a one-period shift to line the two up, the cause was
the MECHANISM, not the timing. `known_shocks` is an exact observation of the
shock; an IRIS `vary` tune sets the mean and leaves the shock random, so its
smoothed shock is *revised away* from the injected value -- 1.3516 against an
injected 1.5 on the fixture -- and no shift of an exactly-pinned path can
reproduce that. `shock_means` is the matching statement. IRIS says as much in
its own source: "The std dev of the tuned shocks remain unchanged and hence
the filtered shocks can differ from its tunes".


A third report against the filtering surface, this one a feature request
rather than a defect: an IRIS-compatible deterministic shock-mean path, and an
explicit state-timing contract on the result.

## Filtering

- **`shock_means` on `kalman_filter()` and `kalman_smoother()`: deterministic
  shock MEANS.** An
  `n_exo x T` matrix of mean shifts (`NA` and `0` both mean "no shift here", so
  the matrix shape `known_shocks` uses works unchanged), with rows matched by
  name.

  **This is a different statement from `known_shocks`, and the difference is
  the point of having both.** `known_shocks` says the REALISATION is known,
  `eps = v`: the shock stops being random, its variance is used up, and where
  it has a prior density the value enters the likelihood. `shock_means` says
  the MEAN is known and the shock **keeps its variance**: nothing is observed,
  nothing about `m` is estimated, and the likelihood gains no term. It is a
  deterministic input, entering the transition and measurement constants as
  `R m_t` and `D m_t` before the update.

  The two coincide exactly for a shock with no prior variance -- knowing the
  mean and knowing the realisation are then the same statement -- and the
  tests pin both halves: identical states and log-likelihood there, materially
  different where the shock still has variance. Agreement alone would prove
  nothing, since an argument that was quietly ignored would also produce it.

  No method routing is involved. The mean path splits off as a deterministic
  trajectory subtracted from the data and added back to the reported states,
  so every `method` evaluates it identically (`known_shocks`, by contrast, can
  only be expressed on the augmented state and forces `method = "univariate"`),
  and zero-variance shocks, missing observations, `a0`/`P0` and both diffuse
  recursions are non-events.

- **`shock_timing`** selects how the columns are read: `"dated"` (default)
  makes column `t` the shock dated `t`, entering `s_t` and `y_t` -- the dating
  `known_shocks` and `shock_scale` already use -- while
  `"transition_next"` reads column `t` as driving the transition OUT of period
  `t`, so it lands on `s_{t+1}`. The second is the one-period adapter shift a
  caller comparing against a package with the other convention would otherwise
  apply by hand; it is now stated in the call.

- **`updated_states` and `predicted_states`**, with the timing in the names:
  `updated_states[, t]` is `s_{t|t}` and `predicted_states[, t]` is `s_{t|t-1}`
  (with `s_{1|0}` from `a0`). `filtered_states` is the same matrix as
  `updated_states`, kept under the name the rest of the package uses. The
  prediction needs nothing extra from the recursions: `E[eps_t] = 0` in the
  deviation system makes `s_{t|t-1} = T s_{t-1|t-1}` exact.

  A unit pulse in `shock_means` with no data to update on therefore traces the
  model's own impulse response exactly -- `updated_states[, t] = T^(t-1) R e_j`
  under `"dated"`, with the two paths coinciding because there is nothing to
  update with. That identity is what a cross-package shock-response comparison
  reduces to once the timing is fixed, and it is asserted to 1e-12 without
  needing the other package present.

- **One deliberate asymmetry between the two entry points.** A known shock's
  REALISATION is fixed, so `kalman_smoother()` reports the injected value back
  as itself; a MEAN leaves the shock random, so `smoothed_shocks` reports
  `m_t + u_{t|T}` -- the mean plus the smoothed deviation around it. A smoother
  that returned the mean unrevised would be ignoring the data; one that
  ignored the mean would be ignoring the input. The two coincide on a
  zero-variance shock, where there is no deviation left to revise.

- `kalman_smoother()` also gains `updated_states` and `predicted_states` under
  the same names, the latter being the mean counterpart of `predicted_cov` --
  computed all along and simply not returned. Note the orientation differs
  from the filter's by long-standing convention: rows are periods here.

- `$diagnostics$shock_means` reports the timing, the number of shifted cells
  and which shocks they belong to, on both entry points. `loglik_type` stays
  `"marginal"`: an input conditions nothing.

# dynhr 0.9.3.2

The follow-up to 0.9.3.1, from a second report against the same surface. Four
behavioural gaps, one of which was a silent wrong answer.

**Two existing calls return different numbers, both deliberately.**

1. `kalman_filter(known_shocks = ...)` on a shock with **zero prior variance**
   used to ignore the injected value completely. The injected row's forecast
   variance is exactly zero, so the sequential filter's `F > kalman_tol` guard
   skipped it and the shock never reached the transition -- no warning, and a
   filtered path that stayed on its unforced trajectory. On the one-state
   fixture in `test-statespace-semantics.R` the filtered state was
   `0 0 0 0 0 0` where the answer is `0 1 0.8 0.64 0.512 0.4096`, and the
   log-likelihood was -10.257 instead of -13.741. If you have conditioned on a
   deterministic shock (`stderr 0`, or a `shock_scale` of zero at that period),
   re-run it.

2. `kalman_smoother()` on a **unit-root model with `lik_init = "auto"`** now
   runs the exact diffuse recursion instead of substituting `P0 = 1e6 * I`.
   The smoothed states move by ~1e-8 (the fallback was accurate); the
   log-likelihood moves by a lot, because the kappa one carried an arbitrary
   additive constant and could not be compared with anything. It now equals
   `kalman_filter(lik_init = "diffuse", method = "univariate")` to machine
   precision. `lik_init = "kappa"` still asks for the old prior explicitly.

Everything else in this release is additive.

## Filtering and smoothing

- **Deterministic known shocks.** A shock with no prior variance is
  deterministic BEFORE conditioning, so an injected value is an INPUT, not an
  observation carrying no information. It is now applied as a mean shift with
  no covariance update -- which is exactly what conditioning on a degenerate
  component is, and exactly the limit of the ordinary update as the variance
  vanishes. A point mass carries no density, so the joint and the conditional
  likelihood coincide there and the filter agrees with the smoother to machine
  precision rather than up to a `log p(eps = v)` term.

- **`kalman_filter()` reports the state hand-off: `final_state`,
  `final_cov`.** `s_{T|T}` and `Var(s_T | y_{1:T})`, named, in the same
  convention `a0` / `P0` take -- so a split sample can be filtered in two
  calls and the prediction-error decomposition holds:

  ```r
  f1 <- kalman_filter(y[1:k, ],     dr, m, p, obs_vars = obs)
  f2 <- kalman_filter(y[(k+1):T, ], dr, m, p, obs_vars = obs,
                      a0 = f1$final_state, P0 = f1$final_cov)
  f1$loglik + f2$loglik           # == the unsplit loglik
  ```

  Pinned to 1e-9 for `standard`, `univariate` and `dare` in
  `test-statespace-semantics.R`. `final_cov` is `NULL` for
  `method = "chandrasekhar"`, which propagates low-rank increments precisely
  so that P is never formed -- reporting NULL is the honest answer, and "auto"
  does not choose that method below `n_state = 100`.

  Before this, `a0` and `P0` existed but there was no way to get a covariance
  out of the filter to put into them.

- **Exact diffuse smoothing (`lik_init = "diffuse"` on
  `kalman_smoother()`).** A new sequential (Koopman-Durbin) smoother on the
  augmented state `[s_{t-1}; eps_t]` -- `R/smoother-diffuse.R`. The augmented
  form is what makes it tractable: dynhr's shocks enter the observation
  equation directly, and the textbook diffuse recursions assume uncorrelated
  noise, which the augmentation restores. It also makes the disturbance
  smoother free -- `E[eps_t | y]` is a block of the smoothed state, so there
  is no second recursion to keep consistent with the first.

  The recursions are derived in the file as the `kappa -> Inf` expansion of
  the ordinary ones, not quoted, and every `L` is a rank-one update, so the
  `N` recursion costs O(nb^2) per observable rather than O(nb^3).

  Certified against `.gls_smoother()` with zero prior precision on the diffuse
  coordinates -- a projection, no recursion, no kappa: **3e-15** on states,
  shocks and the pre-sample initial state, where the old fallback managed
  5e-8. The likelihood matches the filter's exact diffuse value to 12 digits.
  Missing observations do not downgrade it: "missing", "diffuse" and "exactly
  predictable" are three independent decisions per (period, observable) in a
  sequential filter.

- **Structured run diagnostics: `$diagnostics` on both entry points.** Same
  field names on the filter and the smoother, so a parity harness reads either
  without a special case: `method_requested` / `method_used`,
  `lik_init_requested` / `lik_init_used`, `routing` (a data frame of
  `from` / `to` / `reason`, one row per automatic reroute or fallback),
  `diffuse_periods`, `missing_by_period` / `n_missing`, `dropped_by_period` /
  `n_dropped`, `known_shocks` and `loglik_type`.

  `missing` and `dropped` are counted separately on purpose: "not observed"
  and "observed but exactly predictable" shorten the conditioning set for
  different reasons, and conflating them hides a stochastic singularity.
  `loglik_type` is `"marginal"`, `"joint"` (an injected shock with a prior
  density) or `"conditional"` (all injections deterministic; also the
  smoother's convention).

- **`pre_sample` now pads `known_shocks` too.** The backfill rows are
  prepended before the known-shock matrix is validated, so an `n_exo x T`
  matrix on the caller's sample used to be rejected as the wrong width. It is
  padded with `NA` (= unknown) instead. The pre-sample series are also read
  after the deterministic trajectory is added back, so a backfill run together
  with `known_shocks` reports the padded periods with the injected shocks in
  them.

## Testing

- `test-statespace-semantics.R` (new): the four requests plus the requested
  A-F regression matrix over (missing data, zero covariance, unit root, known
  shocks). The matrix asserts **superposition** rather than a pinned number --
  injecting a deterministic shock must equal subtracting its trajectory from
  the data and injecting nothing, which is linearity of the state space and
  needs no oracle.

- `test-kalman-smoother-exact.R` block (d) is inverted: it recorded that the
  smoother had NO exact-diffuse initialisation. It now pins that the exact
  answer IS the `kappa -> Inf` limit -- along the RIGHT approximating sequence
  (`P_star + kappa * P_inf`, diffuse only in the unit-root direction). The old
  `kappa * I` fallback is diffuse in every direction, a different prior whose
  limit is 9.3e-4 away and flat in kappa.

# dynhr 0.9.3.1

A filtering and smoothing release. No export is added or removed and every
addition is an optional argument on a function that already existed, so no
existing call signature changes.

**One existing call DOES return a different number, deliberately.**
`lik_init = "diffuse"` with missing observations used to warn and silently
downgrade to `"kappa"`: -44.097 where the exact diffuse log-likelihood is
-36.271 on the local-level fixture with three gaps. That is the fix in this
release rather than a side effect of it -- but if you have run a
diffuse-initialised filter on data with gaps, you were getting the kappa
answer, and you will now get a different (correct) one.

Every other new argument is inert when omitted, and that is asserted rather
than assumed: `a0 = 0`, `pre_sample = 0` and an all-`NA` `known_shocks` each
reproduce the 0.9.3 result at `tolerance = 0`.

The rest closes gaps that had been worked around by hand: no way to set the
initial state, no way to backfill the latent history before the sample, and no
way to tell the filter about a shock you already know.

## Filtering and smoothing

- **Known historical shocks can be injected: `known_shocks` on
  `kalman_filter()` and `kalman_smoother()`.** An `n_exo x T` matrix carrying
  the value where a shock is known and `NA` where it is not -- the
  `NA`-as-unknown convention `data` already uses, and the `n_exo x T` shape
  `shock_scale` already uses, with rows matched by name. For an announced
  policy change, a measured intervention, or a judgemental adjustment;
  `known_shocks_sd` (filter) makes the injection soft rather than exact.

  The two entry points reach it by different mechanisms, because their
  recursions differ. The filter's univariate path already runs on the
  augmented state `[s_{t-1}; eps_t]`, so the shocks ARE state components there
  and a known shock is an exact observation of one -- passing `known_shocks`
  routes to `method = "univariate"`. The smoother instead splits the known
  shock off as the deterministic part of the system, smooths the remainder,
  and adds it back. Both were checked against closed forms that use no
  recursion at all (`helper-gls-smoother.R`): the smoother matches a
  projection oracle to 2e-16, and the filter matches an analytic joint density
  to 1e-10.

  **Semantics, because it is a modelling choice rather than a detail:** the
  filter treats the known value as an observation and reports the JOINT
  `log p(y, eps = v)`, so a known shock informs that shock's own standard
  error. The smoother conditions. The two differ by exactly `log p(eps = v)`
  -- verified to 0.0e+00 -- so subtract that term if you want the conditional
  likelihood.


- **`kalman_filter()` and `kalman_smoother()` take `a0` and `P0`.** The
  initial state mean and covariance were hard-coded -- zero, and whatever
  `lik_init` implied -- so there was no way to carry a state across a sample
  split, condition on a known history, or hand the smoother a prior of your
  own. The internals had always taken them; only the public shape was missing.
  `a0` is in **deviations from the steady state** (the convention the states
  are reported in; `data` is in levels) and is matched BY NAME when named.
  `P0` must be symmetric and positive semi-definite, may be a scalar for a
  multiple of the identity, and is reported back as `lik_init = "user"`.

  The identity that pins this: splitting a sample at `k`, filtering the first
  block, and starting the second from `(filtered_states[k, ],
  filtered_cov[, , k])` reproduces the whole-sample log-likelihood -- verified
  to 5e-13 on `nk_demo`.

  Two traps were found while wiring it, both of the "accepted, validated,
  reported, then silently dropped" kind: `kalman_standard_loop_cpp()` takes no
  initial state (it starts at zero unconditionally), so the default fast path
  for a stationary model now falls through to the R loop whenever `a0` is
  non-zero; and the smoother's backward pass read a hard zero for
  \eqn{s_{0|0}}, which would have ignored `a0` in `smoothed_initial` while the
  forward pass honoured it.

- **`lik_init = "diffuse"` now works with missing observations.** It used to
  warn and downgrade to `"kappa"`, which answers a different question: on a
  local-level model with three gaps, -44.097 against the exact -36.271. The
  exact-diffuse recursion and missing data meet in the sequential
  (Koopman-Durbin univariate) filter, which skips a missing observable one at
  a time -- as it already did for every other initialisation -- so the
  capability was implemented and unreachable rather than absent. `method =
  "auto"` now routes the combination there; an explicit multivariate method
  reroutes with one notice rather than downgrading.

  The difference is not stylistic: the exact-diffuse filtered states are
  invariant to the diffuse scale to the last bit, while the `"kappa"` states
  keep moving as kappa grows (5.4e-5 from 1e4 to 1e6, 5.4e-7 from 1e6 to 1e8)
  -- an approximating sequence, not the answer. Complete-data results are
  unchanged, and on a single unit root the two diffuse conventions still agree
  to 5e-12.

- **`kalman_smoother(pre_sample = k)`** backfills `k` periods of latent
  history before the first observation and returns them in
  `presample_states` / `presample_shocks` / `presample_cov`, chronological,
  with every other series still aligned to `data`. No new recursion: an
  all-missing period is predict-only, so this is the ordinary backward pass
  over `k` padded rows -- the mechanism that has always produced the single
  `smoothed_initial` period, and the last backfilled row reproduces it
  exactly. The log-likelihood is unchanged and the in-sample states move by
  2e-14. `me_extra` and `shock_scale` are padded with it, so their per-period
  columns stay aligned.

- **How far the smoother's kappa fallback sits from the exact diffuse answer
  is now measured, not assumed.** `helper-gls-smoother.R` computes the
  smoother by projection instead of recursion -- the diffuse case is a zero on
  the prior precision, not a separate algorithm -- and is certified against
  the existing smoother on a stationary model to 3e-16 before being used
  where the answer is unknown. On a unit-root fixture the kappa fallback's
  smoothed states are within **5e-8 absolute on a state of scale 3.3**, about
  eight significant figures. The error falls as `1/kappa` and then RISES
  again as round-off takes over (5.0e-6, 5.0e-8, 4.8e-9, 2.8e-7 at kappa
  1e4/1e6/1e8/1e10), so raising `.DIFFUSE_SCALE` is not the fix -- but nor is
  the fallback the crude approximation it looked like. An exact diffuse
  SMOOTHER remains unimplemented; the tests now state its acceptance
  criterion. Where the diffuse treatment genuinely matters is the
  log-likelihood's unbounded kappa-dependent constant, and
  `kalman_filter(lik_init = "diffuse")` already handles that exactly.

- Multi-period pre-sample backfill needs no manual padding: the argument
  above wraps it. On `nk_demo` with four padded periods
  the in-sample states are unchanged to 2e-14, the log-likelihood is identical
  (missing rows contribute nothing), the transition identity holds across the
  pad boundary to 2e-16, and the last padded row reproduces
  `smoothed_initial`. Exact for a stationary model, since `P_0` is then the
  unconditional covariance; on a unit-root model it inherits the smoother's
  kappa fallback, which awaits an exact diffuse smoother.

# dynhr 0.9.3

This release adds four new estimation and solution capabilities — global
(projection) solutions with their own likelihood, Markov-switching DSGE
filtering and smoothing, mixed-frequency observation blocks, and
moment-based estimation.

It also carries a substantial correctness pass on the parts of the package
you reach for when fitting a model to data. Three of those are worth
singling out, because each was silent:

- `run_full_estimation()` built its sampler proposal from the **prior**
  rather than the posterior curvature — a structurally unreachable branch —
  which mixed poorly on most models and froze the chain outright on a
  well-identified one.
- `kalman_smoother()` never subtracted the **observation intercept**, so it
  silently required data in deviations and disagreed with `kalman_filter()`
  on the same series by tens of thousands of log points — and the historical
  decomposition behind diagnostics D11/D12 was running on exactly that
  mismatch.
- A parameter named **`sigma_e`** was deleted from the calibration before it
  was read, leaving it `NA` while compilation and solving continued.

Measurement error is now one noise model across every filter,
`kalman_smoother()` takes the same arguments as `kalman_filter()`, and every
filtering and smoothing entry point takes observables in **levels**.

The release carries breaking API changes; read the first section before
upgrading.

## Breaking changes

- **Three empty pass-through aliases are removed**, with no deprecation
  shims: use `solve_model()` for `solve_dsge()`, `run_mode_finding()` for
  `estimate_mode()`, and `run_posterior_estimation()` for
  `estimate_posterior()`.
- **Argument spellings are unified across every exported likelihood,
  filter and forecast function.** The data argument is `data` (was `Y`),
  the observable-name argument is `obs_vars` (was `obs_names` or
  `observables`), and the measurement-error argument is `me_variance` (was
  `me_var` or `me_sd`) — matching `make_posterior()` and
  `run_full_estimation()`. Matrix orientation is unchanged and is now
  stated per function.
- **Two of those renames change what you pass, not just the name.**
  `hank_loglik_ar()`, `hank_loglik_ar_grad()`,
  `hank_loglik_ar_structural_grad()` and `hank_simulate_aggregate()` now
  take a VARIANCE where they took a standard deviation — pass `me_sd^2`.
  `make_log_posterior_hank()`, `make_posterior_grad_hank_ar()` and
  `hank_ar_target()` took both `me_var` and an independently overridable
  `me_sd`; they now take the single `me_variance`, so their two likelihood
  branches always describe the same measurement error.
  `hank_loglik_ar_grad()$me` is still the score with respect to the
  standard deviation.
- **`run_posterior_estimation()`'s count arguments follow the
  sampler-level convention**: `nburn`/`ndraws`/`nchains`/`nparticles`/
  `nwalkers` are now `n_warmup`/`n_draws`/`n_chains`/`n_particles`/
  `n_walkers`, matching `run_full_estimation()` and `mcmc()`.
- **The inert `mh_scale` argument is removed** from
  `make_log_posterior_tpf()` and `dynhr_smc2()`'s `likelihood_args`, which
  now rejects unrecognised keys instead of silently dropping them. It had
  done nothing since the tempered particle filter's mutation step was
  corrected to hold the ancestor state fixed.
- **`dynhr_smc2()` returns a `dynhr_chains` object** like every other
  sampler entry point, so `print()`, `summary()` and `plot()` work on it.
  All previous fields are retained.
- **`kalman_smoother()`'s second argument is a decision-rule object, and
  every filtering and smoothing entry point takes observables in LEVELS.**
  Passing a pre-built `dsge_ss` — the original signature — is an error that
  names its replacement; the deprecated `ss = ` alias is gone with it. The
  data convention changed with the shape: a state space now carries its own
  observation intercept in the `d` field `new_dsge_ss()` has documented all
  along, so `kalman_smoother()`, `realtime_decomposition()`,
  `forecast_backtest()` and `conditional_forecast()`'s anchoring data are all
  in levels and the model's steady state is subtracted for you. If your
  series are already in deviations, pass `d = 0` to `kalman_smoother()`.
  Two entry points to the same recursion silently requiring different data
  was the defect; one release of polymorphic shim would have kept it alive
  in a second form, so it is settled here instead.

## New features

### Global (projection) solutions

- **`solve_global()` is usable as an estimation target.**
  `make_log_posterior(likelihood = "global_pf")` runs a bootstrap particle
  filter over a projection solution, so the model is never linearised. The
  estimate is unbiased for the marginal likelihood, making `pmmh()` over it
  a valid pseudo-marginal sampler. `global_pf_sbc()` is the matching
  rank-uniformity certification.
- **The default collocation domain is measured, not guessed.** It is
  derived from the model's own shock process and unconditional state
  dispersion, so an AR(1) with a large innovation gets a wider grid
  automatically; a four-fixture cover study set the default. `solve_global()`
  also validates its model class structurally and fails loud rather than
  silently returning a bad approximation.
- **`euler_errors()` is model-agnostic** — it reads the Euler equations from
  the parsed model instead of assuming RBC parameter names — and
  **`den_haan_marcet()`** is new. Both work on perturbation solutions too,
  so they can be used to decide whether a global solve is needed at all.
- `simulate()` on a `GlobalSolution` had an off-by-one in shock timing: the
  shock drawn in period `t` was applied to the wrong period. Fixed.

### Markov-switching DSGE

- **Filtering, smoothing and IRFs across regimes**: `ms_kim_filter()`,
  `ms_kim_smoother()`, `ms_kim_smoother_struct()` (structural switching) and
  `ms_irf()`.
- **GPB(3) collapse** (`collapse = "gpb3"`) keeps the pair
  `(s_{t-1}, s_t)` and collapses over two lags instead of one. It is
  strictly weaker as an approximation than GPB(2) and markedly more accurate
  when regimes are persistent: against an all-regime-path enumeration
  oracle, smoothed regime probabilities improved from 1.00 away to 2e-15 and
  states from 25 sd to 3e-15 sd on the GPB(2) breakdown draw. It is not a
  pointwise improvement — on 7 of 48 grid draws its error is up to 2.1x
  GPB(2)'s, at absolute levels below 3e-5 — and costs about 1.75x at
  `T = 300`, `h = 2`. `"gpb2"` remains the default and is bit-identical to
  before.
- The filter's covariance update is now Joseph-form, and a collapse
  diagnostic is available via `return_collapse_diag`.

### Mixed-frequency observations

- **`obs_aggregation`** declares an observable as the *k*-period temporal
  aggregate of a higher-frequency model variable, so a monthly model can be
  estimated on quarterly data without leaving the monthly frequency. Four
  aggregators: `flow_sum`, `flow_mean`, `stock_end` and `triangle`
  (Mariano–Murasawa). Implemented as fixed-weight state augmentation
  (Harvey 1989 §6.3), so `ZZ` and `TT` stay constant, the filter's hot path
  and both C++ kernels are untouched, and a model without `obs_aggregation`
  returns a byte-identical log-likelihood. `mf_augment_state_space()`,
  `mf_expand_observations()` and `mf_aggregation_weights()` expose the
  machinery directly.

### Estimation

- **`method_of_moments()`** — GMM and SMM by moment matching over the
  model-implied autocovariance structure, with identity / optimal /
  Newey–West / diagonal weighting, optional two-step, and analytic moment
  Jacobians where dynhr has them.
- **`forecast_backtest()`** — recursive expanding-window out-of-sample
  scoring with CRPS, log score, PIT and interval coverage; re-estimate at
  each origin or roll a single fit forward.
- **Delayed acceptance**: `mcmc(..., screen_fn = )` evaluates a cheap
  approximate likelihood first and only runs the expensive one on proposals
  that survive. The two-stage acceptance ratio keeps the exact posterior
  invariant.
- **Resumable chains**: `mcmc(checkpoint_dir =, resume =, flush_every =)`
  writes a checksummed, atomically-written state pack carrying the position,
  log-posterior, adaptation state and `.Random.seed`, so a resumed chain is
  statistically identical to the uninterrupted run. `mcmc_chain_state()`,
  `mcmc_chain_save()`, `mcmc_chain_restore()` and `mcmc_chain_extend()` are
  the sampler-agnostic primitives; restore refuses a tampered pack or one
  saved mid-adaptation.
- **`dynhr_model()` and the `dm_*()` verbs** (`dm_solve`, `dm_posterior`,
  `dm_mode`, `dm_sample`, `dm_diagnostics`, `dm_forecast`, `dm_irf`,
  `dm_test`) — one pipeline object carrying model, compiled, steady state,
  decision rules, data and priors, instead of threading six arguments
  through every call. Every vignette pipeline routed through the object
  returns bit-identical numbers to the functional path.

### Model input and output

- **`write_mod()`** serialises a parsed model back to Dynare `.mod` source.
  A `parse_mod()` -> `write_mod()` -> `parse_mod()` round trip is a cheap
  check that dynhr read a file the way you meant it.
- **`histval` blocks are parsed** into `model$histval` (lag-indexed), and
  **`smoother2histval()`** builds that history from a completed smoother run
  — the standard way to start a forecast or counterfactual from where the
  data left off.
- **`shock_groups` blocks are parsed** into `model$shock_groups`, consumed
  by the decomposition functions.

### Heterogeneous agents

- **`hank_ks_aggregate_risk()`** — Krusell–Smith with genuine aggregate
  risk, plus `hank_ks_risk_irf()` for generalised impulse responses and
  `hank_ks_ergodic_mean()` for the aggregate precautionary term without a
  full simulation. `hank_tfp_chain()` builds the aggregate productivity
  chain.

### Decompositions

- **Shock decompositions add up exactly.** `historical_decomposition()` and
  its OBC variant gain an `initial` column for the contribution of the
  initial state, which is not zero unless the sample starts at the steady
  state; the columns now reproduce the data to machine precision. Both take
  `shock_groups`. **`realtime_decomposition()`** re-runs the decomposition
  across data vintages, so a given quarter's story can be tracked as it was
  revised.

## Correctness

### Measurement error is one noise model everywhere

Measurement error was implemented inconsistently across filters: on several
paths it entered the forecast covariance only, acting as a regulariser
rather than as observation noise, which made those likelihoods disagree with
the exact Kalman filter by O(`me_variance`).

- The **multivariate Kalman filter**, the **Markov-switching filters** and
  the **SV Rao-Blackwellised particle filter** (`kf_step()` and its compiled
  kernel) now all implement the true i.i.d. law: `me_variance` enters the
  forecast covariance AND the Joseph state-covariance update.
- **`dynhr_sbc()`'s data-generating process now adds measurement error**, so
  the DGP and the likelihood describe the same model. An SBC on a filter
  with `me_variance > 0` against a noiseless DGP was certifying a
  mis-specification.

### Filtering and smoothing

- **`kalman_filter(method = "chandrasekhar")` is exact again**, and is
  re-admitted to `method = "auto"` above `n_state = 100`.
- **`kalman_smoother()`'s state pass is exact** under dynhr's timing
  convention.
- **`pkf_smoother_obc()`** uses the correlated-noise disturbance smoother.
- The Kim smoothers gained a joint-probability regime pass
  (`regime_pass = "joint"`).
- **`kalman_smoother()` rejected systems `kalman_filter()` handled.** Reported
  for a unit-root, `shock_scale`d model with a non-positive innovation
  covariance. The two functions treated a singular `F` by different
  mechanisms: `kalman_filter()` falls back to the univariate (Koopman-Durbin)
  filter, which skips any component whose conditional variance is below
  `kalman_tol` -- the correct treatment, since such a component is predictable
  exactly and carries no information -- while `kalman_smoother()` added JITTER
  on an absolute ladder (`1e-8` ... `1e-2`, then an unguarded
  `chol(F + 0.1 I)`). `F` is not O(1): a unit-root smoother starts from
  `P = 1e6 I` and `shock_scale` multiplies `Q` on top, so the ladder was
  either far too small -- and the unguarded rung threw, which is the reported
  rejection -- or it "worked" and silently corrupted the result. Switching one
  shock off via `shock_scale` (the `u_k = 0` idiom for forcing a series to its
  observed value, and what a hard `filter_tunes` tune does underneath) moved
  the smoother's log-likelihood to **-8.5e+09** where the filter returned
  **-63.7**. The smoother now makes the filter's decision: a zero-variance
  component is dropped for that period, exactly as it already treats a
  *missing* observable, and the update proceeds on the informative subset,
  whose `F` is positive definite by construction. The value becomes -70.0,
  and the smoother's offset from the filter is now the same constant whether
  or not a shock is switched off. Dropped components are reported, naming how
  many periods and components and stating that the log-likelihood is not
  comparable with an undropped run.
- **`kalman_smoother()` silently required data in DEVIATIONS.**
  `kalman_filter()` takes raw level data and subtracts the observation
  intercept `d = dr$ys[obs_vars]`; the smoother never did, and
  `build_dsge_state_space()` carries no steady state at all, so the
  requirement was unstated and its violation silent. On any model whose
  observables have non-zero steady states -- which is most of them --
  filtering and smoothing the same series disagreed wildly: on the bundled
  `nk_demo` (observable steady states 0.5, 2 and 4) the filter returned
  -757.6 and the smoother -33990.3. Hand-demeaning the data closes the gap to
  2e-13, confirming the intercept was the whole of it. The state space now
  carries the intercept (`build_dsge_state_space()$d`, rescaled by `sum(w)`
  under a mixed-frequency aggregator exactly as `kalman_filter()` rescales
  its own), every entry point subtracts it, and `d = 0` is the explicit
  escape hatch for data already in deviations.
- **The historical decomposition behind D11/D12 was running on levels
  through the deviations-only path.** `run_all_diagnostics(posterior)` built
  a state space, handed it the raw observable columns and smoothed them, so
  on any model with non-zero observable steady states the smoothed shocks
  absorbed the level offset: on `nk_demo`, `max |eps|` 4.83 against 0.80 —
  six times too large — and a log-likelihood of -33990.3 against -757.6. D12
  reports the mean of each smoothed shock and passes it at `|mean| < 0.1`, so
  the diagnostic was reporting the bug as a model failure. **Any historical
  decomposition or smoothed-shock series produced through
  `run_all_diagnostics()` before this release, on a model whose observables
  have non-zero steady states, should be recomputed.** The same latent defect
  is closed in `realtime_decomposition()`, `forecast_backtest()` (whose
  predictive mean now carries the intercept back so it is scored against the
  realised level) and `conditional_forecast()`, whose Gaussian branch alone
  took deviations while its `tpf` and `pskf` branches demeaned for
  themselves — the same call needed different data depending on
  `ctx$likelihood`.

- **`kalman_smoother()` takes the same arguments as `kalman_filter()`.** Every
  other filter/smoother entry point -- `kalman_filter()`, `ms_kim_filter()`,
  `ms_kim_smoother()`, `kf_innovation_diagnostics()` -- takes
  `(data, dr, model, params, obs_vars, me_variance, ...)`. The Gaussian
  smoother took `(data, ss, Q, me_extra, shock_scale)`, so the obvious call by
  analogy after filtering failed and the only signpost to the required
  `build_dsge_state_space()` step was a single `@param` line; `?kalman_smoother`
  had no example and no `\\seealso`. It now accepts the filter's arguments,
  including **`me_variance`**, which it previously lacked entirely -- so a
  model filtered with measurement error could not be smoothed under the same
  noise model. Filter and smoother now agree to 1e-8 across `me_variance`
  values. The pre-built state space is no longer accepted: see the breaking
  changes above for the one-line migration.
- **`kalman_smoother(lik_init = )`.** `kalman_filter()` refuses
  `lik_init = "auto"` together with `shock_scale` on a nonstationary model and
  instructs the caller to pass `"kappa"` or `"stationary"` explicitly -- which
  the smoother had no way to accept, so the two could not be made comparable
  even in principle. `"stationary"` now errors on a unit root instead of
  silently returning a kappa-initialised answer. There is still no exact
  *diffuse* initialisation in the smoother.
- **The same defect is fixed in `conditional_forecast()`**, in four places.
  Its internal Kalman pass added a `1e-10` ridge on every period (so it was
  never an unregularised filter) and fell back to an unguarded
  `chol(F + 1e-6 I)`; and its three condition-system solves used unguarded
  `chol(M + 1e-12 I)`. Those Gram matrices lose rank exactly when the
  conditions are collinear, over-specified, or routed through a switched-off
  shock -- normal things to ask for. All four now drop uninformative
  components or take the minimum-norm solution via a relative-cutoff
  pseudo-inverse, and warn rather than throwing.

### Parsing

- **A parameter named `sigma_e` was silently dropped** (reported against
  0.9.1). `remove_blocks()` strips Dynare *command* statements before the
  calibration is read, and its keyword list ends with `Sigma_e` — Dynare's
  shock-covariance assignment — but the loop matched case-INSENSITIVELY. A
  user's `sigma_e = 1;` was therefore deleted as if it were that command, and
  the parameter survived declared but `NA`, after which `compile_model()` and
  the solvers happily proceeded on an invalid calibration behind a warning.
  Dynare identifiers are case-sensitive and its command is spelled with a
  capital S, so lower-case `sigma_e` is an ordinary parameter name; it is now
  matched case-sensitively. `Sigma_e` was the only entry in that list whose
  lower-case form is a legal user identifier — the rest (`stoch_simul`,
  `steady`, `check`, …) are genuine Dynare reserved words and stay
  case-insensitive.

### Sampling

- **The one-call estimation API proposed from the PRIOR, not the posterior.**
  `run_full_estimation()` and the estimation runner both built the RWMH
  proposal as `if (!is.null(mode_res$V_mode)) ... else diag(prior_spec$std^2)`,
  but the mode result on that path comes from the optimiser core
  `.run_mode_finding()`, which returns only the mode and its convergence
  record and never sets `V_mode`. The condition was structurally unreachable,
  so every proposal was a prior-variance diagonal and the posterior curvature
  the mode-finder had just located was silently discarded. On a model whose
  prior and posterior sit at a similar scale this merely mixed poorly
  (`fs2000`: 10.2% acceptance); on a well-identified one it froze the chain
  outright — 0% acceptance and exactly zero posterior variance, every draw
  equal to the mode. The Hessian-based proposal logic
  (`.make_pd` regularisation plus eigen-basis capping at the prior scale) is
  now shared with `run_mode_finding()` rather than duplicated, and the
  samplers compute the Hessian at the mode and use it. Measured after the
  fix: the nine-parameter `nk_demo` fixture goes 0% -> 34.3% acceptance with
  every posterior mean within one standard deviation of the values the data
  were simulated at, and `fs2000` goes 10.2% -> 29.4%. If the Hessian genuinely
  cannot be evaluated the proposal still degrades to the prior diagonal, but
  now **warns**: a silent version of that fallback is what hid this for three
  releases. `test-sampler-proposal.R` pins the property no test had asserted —
  that the chain moves at all.

### Numerics

- **The SV RB-PF agrees between R and C++ past the volatility overflow
  point.** The two diverged without bound once a volatility particle left
  the useful double range. The overflow was not the cause: the forecast
  covariance and its Cholesky factor are bit-identical in both, but
  `chol2inv()` and `arma::inv_sympd()` differ by one ulp in `F^-1`, and when
  the observation block is perfectly informative the exact Kalman gain is
  the identity, so the Joseph factors are exactly zero. That one ulp made
  them entirely rounding noise, the state covariance rounding noise squared,
  and the next period's inverse amplified it without limit. Both kernels now
  snap a Joseph entry lying within a few ulps of the magnitudes that
  cancelled to the exact zero it approximates; `kf_step()` also rejects a
  non-finite forecast covariance before `chol()`, as the compiled kernel
  already did. Verified over 864 parameter/seed/length combinations:
  0 divergent, worst relative gap 4.1e-16.
- **Order-3 cumulant moments were projected onto only the `(i,i,k)` slice**
  of the third-moment tensor. Fixed, along with three further cumulant/GMM
  defects (the analytic GMM weight matrix's lag handling among them).
- **Seeded particle-filter closures no longer reset the caller's RNG
  stream** — a seeded closure used to reseed the global stream on every
  evaluation, silently correlating an outer sampler's own draws.
- **`power` (power-posterior tempering) reaches every likelihood, exactly
  once**, and is shipped to parallel workers.
- **Stationary initial covariances (PSKF, TPF) come from the real Lyapunov
  solution** rather than a truncated series.
- **Order-2 shortcuts no longer discard shock correlations**: `.linear_dr2`
  dropped the off-diagonal of `Sigma_e`.
- **`kalman_filter()`'s singularity fallback is conditional and audible**
  instead of silent, and `.safe_inv()` truncates on a relative singular-value
  cutoff rather than an absolute one.
- **`hank_het_block()` fails loud on reducible income chains** and gains
  `dist_init` — at `p_un = p_nu = 0` the non-participation state is a closed
  class, and a uniform-seeded power iteration stranded about a third of the
  mass there.
- **`dynhr_set_options()` values now reach mirai daemons**, and parallel
  workers receive the host's option state.
- The cumulant gradient no longer returns a silent `NaN` for an explosive or
  non-stationary draw.

### API and structure

- Two exported functions failed on every call and had no test:
  `ramsey_obc_pwlinear()` built its shock sequence transposed, and
  `diag_prior_sensitivity()` referenced an undeclared argument. Both fixed,
  with smoke tests added for every previously untested export.
- Five `print()`/`summary()` methods were written but never registered, so
  from an installed package they fell through to `print.default()` and
  dumped the whole object. All registered, with a structural guard against
  recurrence.
- `run_posterior_estimation()` with `n_chains >= 2` crashed on any
  one-parameter model. Fixed.
- The HANK result classes share one compact `print.hank_block()`, so
  printing a block no longer dumps the stationary distribution. Nineteen
  internal-but-exported oracles are marked `@keywords internal`.
- One discrete-Lyapunov solver (`solve_lyapunov()`); the removed
  direct-Kronecker variant returned a false `NaN` on highly non-normal
  stable matrices and a negative variance for a scalar explosive root.
- `ast_to_string()` under-parenthesised `a - (b - c)`.

# dynhr 0.9.2

The headline addition is **SMC² (`dynhr_smc2()`)**: sequential Monte Carlo
over parameters driven by the package's particle-filter likelihoods — the
order-3 tempered particle filter and the stochastic-volatility
Rao-Blackwellized filter, both SBC-certified in this cycle. Supporting it,
the tempered particle filter received a substantial correctness and
performance pass, and the HANK household line gains three-state (E/U/N)
labour status, asset-indexed transfer incidence, and a complete
distribution-Jacobian family.

dynhr now requires **R >= 4.3.0** (declared via `Depends` and
`SystemRequirements: C++20`): the compiled two/three-asset solvers use
`std::barrier`, which needs GCC >= 11 — Rtools43 on Windows. This
formalizes the 0.9.1.1 hotfix (`CXX_STD = CXX20` in both Makevars, so
Windows toolchains no longer fall back to C++17 and fail at
`#include <barrier>`).

## New features

- **`dynhr_smc2()` — SMC² over parameters.** Likelihood-tempered SMC on
  theta where each evaluation is an unbiased particle-filter estimate
  (`likelihood = "tpf"` or `"sv_rbpf"`), with pseudo-marginal mutation
  moves: every theta-particle carries its stored likelihood estimate, the
  incumbent is never re-evaluated, and fresh filter randomness is drawn
  only at proposals — so the final-stage marginal is the exact posterior
  and the telescoped evidence estimate is valid. Filter-level tuning
  routes through `likelihood_args =`; `parallel = TRUE` evaluates
  theta-particles on a mirai daemon pool. Oracle-tested against exact-KF
  SMC on degenerate-volatility and linear-observation models.
- **Order-3 tempered particle filter** (`make_log_posterior_tpf(order = 3)`):
  the particle transition implements the pruned third-order recursion
  exactly, verified against the reference simulator and the order-3
  Gaussian pruned-KF, and **certified calibrated by a full R = 100
  rank-uniformity SBC** (`tpf_order3_sbc()`, the new full-certification
  tier under `DYNHR_SBC_FULL`, joining `sv_rbpf_sbc()`).
- **Three-state (E/U/N) household labour status**
  (`hank_employment_income3()`): all six transition rates are
  Jacobian-ready inputs (direct N<->E flows default to zero — a
  documented, overridable restriction), non-employed income is
  calibratable, and the three stocks plus six gross flows are selectable,
  ND-verified outputs. With the participation margins zeroed the
  two-state model is recovered (chain, income and shares bit-identically;
  the jointly-solved policies to machine precision).
- **Asset-indexed transfer incidence**: `hank_het_block(Tr_incidence =)`
  accepts an `n_e x n_a` matrix `omega(e, a)` over beginning-of-period
  states (the transfer stays lump-sum), enabling wealth-correlated
  incidence schedules. The block reports `Omega_ss` and a new `"Omega"`
  Jacobian output (the fiscal outlay aggregate); cash-on-hand positivity
  is asserted at the steady state. Vector and default incidence are
  bit-identical to 0.9.1.
- **Distribution-Jacobian family completed**: new
  `hank_het2_dist_jacobian()` / `hank_het3_dist_jacobian()` (+ `_nd`
  numerical oracles) give the full cross-sectional distribution response
  via the fake-news algorithm, including the `theta_coll` (collateral
  LTV) column; `hank_model_dist_irf()` gains the `"het2"`/`"het3"` arms.
  The three-asset backward sweep is now a single compiled call with a
  persistent worker pool, and the expectation stream is shared across
  inputs — both bit-identical to the R paths they replace.
- `hank_model()` fails loud on structural GE singularity, naming the
  offending unknown/target instead of an opaque LAPACK error.

## Tempered particle filter: correctness and performance

All four items below change `likelihood = "tpf"` values relative to 0.9.1;
fixed-seed pins were re-computed.

- **Observation-equation fidelity fix.** The TPF previously used a
  linearized observation equation, dropping the quadratic/cubic
  observation-row tensors the pruned model's data actually carry — on
  models with nonlinear observables it evaluated the wrong model's
  likelihood (found by SBC certification as a decisively miscalibrated
  posterior). The full observation reconstruction is now used at both
  orders and in both backends. Found in passing: the compiled order-2
  kernel had a Kronecker index-layout bug in its `ghxu` term, active only
  for multi-shock models with `n_e != n_s`; R and C++ backends now agree
  on a full-filter multi-shock run and a permanent parity test guards it.
- **Mutation-kernel unbiasedness fix.** The RWMH mutation step previously
  random-walked the ancestor state with a likelihood-only acceptance
  ratio — a kernel that biased the likelihood estimator upward,
  increasingly with `n_mh`. Mutation now refreshes only the period-t
  shock with the ancestor state fixed (Herbst-Schorfheide 2019); the
  estimator is verified unbiased against the exact true-ME likelihood
  and the order-3 SBC certification was re-run under the corrected
  kernel. `mh_scale` is now inert (kept for API compatibility).
- **Missing-data support.** Fully-missing periods now propagate the
  particle cloud (previously the period was skipped entirely, gluing the
  gap's endpoints together — the likelihood of a different model);
  partially-missing periods evaluate the observed elements only, matching
  the univariate Kalman filter's per-element NA convention. Fully
  observed datasets are unchanged.
- **Stationary particle initialization by default** (`burn_in_init = 50`):
  the cloud starts from (approximately) the full pruned stationary joint
  instead of zeroed higher-order layers, removing a short-sample
  persistence bias identified by SBC. Pass `burn_in_init = 0` to
  reproduce 0.9.1 behavior (required with correlated-pseudo-marginal
  `U_list` evaluation, which has no burn-in slots).
- **~8x faster evaluation** at estimation-scale settings: the per-particle
  `kronecker()` loops were replaced with bit-identical row-indexed
  products and per-period invariants hoisted out of the mutation loop.
  A full order-3 SBC certification now takes ~40 minutes, down from ~10
  hours.

## Fixes

- Direct `kalman_filter()` calls with `me_variance > 0` no longer re-run a
  150-step R-level Riccati detector on every call (a ~13x per-call
  overhead; estimation closures were unaffected) — the result is memoized
  and the iteration exits on convergence.
- `sv_rbpf` no longer zero-weights a particle when `inv_sympd()` fails at
  the PD boundary despite a successful Cholesky (the inverse is recovered
  from the factor), and a period where every particle fails now warns
  instead of silently returning `-Inf`.
- `run_full_estimation(likelihood = "tpf")` no longer errors when
  `tpf_options` carries orchestration-level keys (e.g. `cpm_rho_u`).
- `sbc_uniformity_test()` takes the true rank support `L` explicitly (the
  observed maximum can under-count bins) and its tail-asymmetry statistic
  is exactly centered under the null for unequal edge bins.
- Rank-histogram plotting, preflight warning messages, and several Rd/
  documentation defects (including two README function names that did not
  exist) were corrected.

# dynhr 0.9.1

The headline addition is a first-class **stochastic-volatility (SV) estimation
layer** — latent, per-shock time-varying volatility on the shocks of a linear
DSGE, estimated end to end. No named peer (Dynare, gEcon, MacroModelling.jl)
ships a first-class SV declaration/estimation path. The other changes in this
cycle are package hygiene, bug fixes, and incremental improvements.

## New features

- **Stochastic-volatility-on-shocks estimation.** Declare independent AR(1)
  log-variance processes on any subset of a model's shocks
  (`stochastic_volatility;` mod-file block; `stochastic_volatility()` /
  `sv_entry()` constructors) and estimate their hyperparameters jointly with
  the structural parameters. Conditional on a volatility path the model is
  exactly linear-Gaussian with `shock_scale = exp(h_t/2)`, so estimation uses
  a **Rao-Blackwellized particle filter** — particle-filtering only the
  low-dimensional log-variance states and integrating the DSGE states
  analytically via the Kalman recursion — giving an unbiased marginal
  likelihood that plugs directly into `pmmh()`. New API surface:
  `likelihood = "sv_rbpf"` in `make_log_posterior()`, `run_mode_finding()`,
  and `run_full_estimation()`; a `stochastic_volatility=` argument on the
  runners; and the exported single-step Kalman primitives `kf_step()` /
  `kf_stationary_init()`. Scope (v1): stationary models, order-1 (linear)
  solutions only — structural SV (volatility perceived in the decision
  rules) remains available via the order-2 pruned/TPF likelihoods.
- **Cross-machine benchmark.** `dynhr_benchmark()` runs a fixed, realistic
  estimation workload — Smets & Wouters (2007), 36 estimated parameters, 7
  observables, 160 quarters — through RWMH at a sweep of core counts and
  reports per-draw/per-second throughputs (`draws_per_sec`, `us_per_draw`,
  `speedup`, `efficiency`) alongside the system detail
  (`dynhr_system_info()`: hardware, OS, R build, BLAS/LAPACK, dynhr build)
  needed to compare machines. The workload is fingerprinted by its log
  posterior at the published mode; two results are comparable only if the
  fingerprint and the BLAS agree. The model, data, published mode and mode
  Hessian ship in `inst/extdata/models` (AER replication deposit
  openicpsr-116269-V1, BSD-3-Clause/CC BY 4.0; see `sw2007_SOURCE.md`).

## Fixes

- **Gradient-based estimation of pruned-state-space models.**
  `run_posterior_estimation()` / `find_mode()` now use the analytic
  (order-2) / semi-analytic (order-3) pruned-state-space gradient for
  `likelihood = "pruned"` instead of silently falling back to numerical
  finite differences — a stale gate had excluded it although the gradient
  was implemented and tested. This makes structural-SV order-2 models
  NUTS/HMC-able. Validated against `numDeriv` finite differences.
- **Three-asset fake-news accounting.** The forward derivative now streams
  the actual plus/minus policy and transition legs through the Young
  operator, preserving mass and all three asset first moments at active
  bounds (previously reconstructed and clipped a second symmetric policy
  pair). The analytic consumption Jacobian is assembled from the exact
  household budget, including the internal income-distribution response
  and the multiplicative foreign-price terms, so the six reported output
  Jacobians satisfy the aggregate linear budget identity to roundoff.
  Fixed-`Pi` price inputs (including transfer blocks, where the precondition
  is that income depends on the income state alone — satisfied by
  `y = w*e + Tr*omega` for any `Tr`) now carry an exact zero
  income-distribution response, avoiding an `O(N)/step` summation residue
  on large grids; a transfer block previously fell through to the
  accumulated path and carried a residue five orders larger.

# dynhr 0.9.0

This is a feature release. The headline addition since 0.8.1 is a complete
**heterogeneous-agent (HANK) estimation and welfare stack**, alongside an
**exact-Hessian curvature toolkit** and new **chain / determinacy diagnostics**.
Everything from 0.8.1 is retained and unchanged.

## Highlights

### Heterogeneous-agent (HANK) estimation

dynhr can now solve, estimate, and do welfare analysis on heterogeneous-agent
New-Keynesian models end to end, not just linearise them.

- **Solving.** Krusell–Smith / one-asset steady states and sequence-space
  Jacobians (`hank_ks_steady()`, `hank_ks_model()`, `hank_egm_solve()`,
  `hank_het_jacobian()`), a global nonlinear transition solver
  (`hank_td_nonlinear()`), a finite-Reiter linearisation path
  (`hank_finite_solve()` / `hank_reiter_statespace()`), and a
  **discount-heterogeneity mixture economy** (`hank_mixture_ks_steady()`,
  `hank_mixture_ks_assemble()`, `hank_mixture_ks_model()`).
- **Estimating.** `hank_mixture_joint_logpost()` evaluates a joint posterior
  over three information channels — aggregate macro dynamics (Kalman),
  the stationary cross-sectional wealth **level**, and the cross-sectional
  **response** to a price shock. Recommended posterior routes are the exact
  grid (`hank_mixture_sbc()`) or a Laplace approximation
  (`hank_mixture_laplace()`), **not** a diagonal random-walk sampler — the
  channels are strongly non-diagonal. `hank_mixture_emulator()` provides a
  distribution-agnostic surrogate, and `hank_mixture_sbc()` ships
  simulation-based-calibration certification (marginal and joint
  test-quantity ranks).
- **Welfare.** `hank_welfare_posterior()`, `hank_value_transition()`,
  `hank_mixture_welfare_pool()`, `hank_cev()`, and
  `hank_welfare_channels()` decompose consumption-equivalent welfare
  (interest- vs labour-income incidence) across the wealth distribution.
- **Identification.** The headline scientific finding baked into these tools:
  the discount-rate *spread* in a mixture economy is identified by the
  stationary wealth **level** (à la cstwMPC), not by the price-shock response
  (`hank_partial_id_level_response()`, `hank_reweight_level_metric()`).

### Exact-Hessian curvature and gradients

- **`posterior_hessian()`** is now exported, with four second-order-term
  methods — `t2_method = "loop"`, `"contract_once"`, `"hvp_solution"`, and
  `"adjoint_solution"` (the last two avoid forming the state-space second
  derivative; `"adjoint_solution"` is exact and finite-difference-free) — plus
  a `check_mode` guard that diagnoses a Σ_ε mismatch at the mode.
- `posterior_hessian_fd_grad()`, `laplace_log_marglik()`, and
  `make_posterior_grad(grad_method = "adjoint_solution")` round out the
  curvature stack.

### Diagnostics and determinacy

- **`chain_diagnostics()`** — MCMC chain summaries (R-hat, ESS, etc.).
- **`bk_distance()`** / `solution_pencil_spectrum()` — Blanchard–Kahn
  determinacy distance and the generalised-eigenvalue spectrum of the
  solution pencil.
- `kf_innovation_diagnostics()` — Kalman innovation whiteness checks.

### Pathological-DSGE estimation robustness

- `check_hessian_conditioning()`, `fd_safe_hessian()`, `profile_ci()`,
  `run_estimation_passport()`, and `make_loglik_contrib()` support inference
  on weakly-identified / ill-conditioned posteriors.

### Other additions

- Order-3 pruned state space: `pruned_state_space3()`, `pruned_ss_moments3()`,
  `pruned_ss_loglik3()`.
- `compute_fourth_cumulant(method = "closed_form")` — opt-in chain-exact /
  contemp-only fourth-cumulant trace (default remains `"window"`).
- A canonical occasionally-binding-constraint regime-path object
  (`new_obc_regime_path()` / `as_obc_regime_path()`).

## Two- and three-asset HANK households

- **Two-asset (liquid/illiquid) household**: liquid `b` at return `rb`,
  illiquid `a` at return `ra` with a convex adjustment cost
  (`hank_egm2_solve()` / `hank_het2_block()`, EGM over the `(Vb, Va)` pair,
  R + ~7x-faster C++ backends), a joint `(e, b, a)` distribution
  (`hank_forward_operator2()`, `hank_aggregate2()`), a fake-news Jacobian
  (`hank_het2_jacobian()` / `hank_het2_jacobian_nd()`), and GE composition
  (`hank_het2_block_spec()`, `hank_twoasset_steady()`, `hank_twoasset_model()`).
  With `ra > rb` this produces **wealthy hand-to-mouth** households
  (liquid-constrained while holding substantial illiquid wealth) — the
  Kaplan-Moll-Violante mechanism a one-asset HANK cannot represent. Ported
  from and validated against the sequence-jacobian reference implementation.
  Estimation works through the existing `hank_state_space()` /
  `hank_kalman_loglik()` / `hank_loglik_ar()` likelihood stack with no new
  code once a block's Jacobian feeds the block-DAG. Calibration gates:
  `hank_twoasset_grid_check()` (grid adequacy), `hank_theta_boundary_check()`
  (sequence-space horizon adequacy), `hank_twoasset_htm_stats()` (reports
  floor-mass and policy-constrained mass separately), and
  `hank_euler2_residual()` (FOC-residual oracle). Two-asset blocks are
  rejected loudly by every one-asset entry point and vice versa.
- **Three-asset household** (domestic liquid `d`, gross foreign `f`,
  domestic capital `a`): `hank_het3_block()` / `hank_egm3_solve()`,
  `hank_het3_jacobian()` / `hank_het3_jacobian_nd()`, `hank_euler3_residual()`
  (the Stage-4 budget/Euler diagnostic), and `hank_het3_manifest()` (a
  reproducibility record: version/commit, backend,
  grid, iterations, convergence gaps, wall time, peak RSS). Warm continuation
  across a calibration sweep via `Vd_init`/`Vf_init`/`Va_init` and `relax`
  (defaults reproduce the cold-start solve exactly); `hank_egm3_regrid_values()`
  interpolates marginal values across grid refinements. Non-convergence
  returns a `hank_het3_block_failed` object (with `strict = FALSE`) reporting
  iterations, `tol`, `relax` and both convergence gaps, rather than erroring
  or silently propagating into downstream Jacobians. `threads` controls
  cross-platform parallelism over the household loop; output is
  bit-identical at every thread count. The foreign valuation channel `px`
  (the price of foreign claims) is a first-class sequence-space input on
  `hank_het3_block()`, `hank_td3_nonlinear()`, and both Jacobians. Adjustment
  resources (`chi`, `phi`, aggregated as `CHI`/`PHI`) are reported so the
  household budget and the goods-market resource identity close exactly. A
  fix that forms consumption from the budget after fixing assets (matching
  the two-asset kernel's convention) changed the fixed point and is also
  roughly a **10x speedup**; every three-asset performance figure recorded
  before it is superseded.
- **Discrete-adjustment two-asset household** (fixed-cost/taste-shock,
  internal): `hank_egm2d_solve()`, `hank_het2d_block()` (per-branch
  policies, adjust probability `P`, `(V, Vb, Va)` envelope,
  `hank_forward_operator2d()`), `hank_td2d_nonlinear()` (the nonlinear
  transition), `hank_het2d_jacobian()` / `hank_het2d_jacobian_nd()`,
  `hank_het2d_block_spec()` (`kind = "het2d"`).
  A KiwiSaver-style locked-contribution channel (`phi_contrib`) resolves a
  fixed-cost participation trap in the stationary distribution. Remains
  internal (not exported) pending the full stack; C++ kernels are a
  documented follow-up.
- **Bond pricing and revaluation**: `hank_bond_block()` / `hank_bond_ss()` —
  a geometric (delta-coupon) bond, the cleanest revaluation instrument
  (fixed coupons make realized-return movements pure revaluation).
  `hank_td2_reval_decompose()` runs a household block on full /
  cash-flow-only / revaluation-only return paths to separate wealth from
  cash-flow incidence of a rate change.
- **Uniform lump-sum transfer `Tr`** on both one- and two-asset household
  blocks (`y = w*e + Tr`), a first-class Jacobian input with `Tr_path` on
  `hank_td_nonlinear()` / `hank_td2_nonlinear()`. Backward compatible:
  `Tr = 0` is the default and reproduces prior behavior exactly.
  `hank_impc()` now differences around the transfer-inclusive income
  definition.

## Debt-side channels: collateral, repricing, and the Fisher effect

- **Collateral-linked borrowing**: households may borrow beyond the
  unsecured liquid floor by pledging illiquid wealth,
  `b' >= b_grid[1] - theta_coll * a'`, at loan-to-value `theta_coll`, on
  the one-asset block. `theta_coll = 0` reproduces the pre-collateral
  solver exactly; a time-varying `theta_t` is supported, with a `theta_1`
  argument re-basing the date-1 distribution when the initial LTV differs
  from the block's steady-state value.
- **Staggered debt repricing**: `hank_repricing_block()` — the effective
  rate on the debt stock reprices by a fraction `phi_r` per period
  (`phi_r = 1` is instant repricing / the identity).
- **Borrowing wedge**: `r_minus` is a second aggregate rate input on
  `hank_het_block(r_minus = )`, letting deposit and borrowing rates diverge
  and reprice on different schedules. Wedge-unaware routines
  (`hank_impc()`, `hank_mpc()`, `hank_euler_residual()`,
  `hank_het_dist_jacobian()`, `hank_sam_reiter_linearize()`) now refuse a
  wedge-carrying block loudly rather than silently pricing the debt side at
  the saving rate.
- **Fisher channel**: `hank_fisher_block()` (anticipated inflation as an
  exact real-rate block) and `hank_liquid_reval_d0()` (the unanticipated
  date-1 nominal-stock revaluation, applied only to the liquid axis). The
  aggregate consumption response to surprise inflation follows the net
  nominal position of the household block, so per-group incidence is the
  robust object to report.

## General sequence-space engine and estimation

- **`hank_model()`** generalizes the Krusell-Smith-specific GE solve into
  an arbitrary directed-acyclic-graph of blocks: `hank_simple_block()`
  (analytic or finite-difference block Jacobian) and
  `hank_het_block_spec()` (fake-news), composed by the chain rule into GE
  Jacobians `H_U`/`H_Z` and solved via `hank_model_irf()`. `hank_ks_model()`
  builds Krusell-Smith through it, and `hank_nk_hank()` assembles a
  one-asset New Keynesian HANK (heterogeneous households, Taylor rule,
  Fisher equation, output-gap NKPC, bond clearing).
  `hank_model_nonlinear_irf()` is a general nonlinear perfect-foresight
  transition solver for any `hank_model` object. GE factorizations are
  cached and reused (keyed on `H_U`, bit-identical results); disable via
  `options(dynhr.hank_ge_cache = FALSE)`.
- **Exact-AR(1) stacked-covariance likelihood**: `hank_loglik_ar()` /
  `hank_loglik_aggregate_ar()` / `hank_autocov_ar()` replace the
  truncated-MA approximation's dropped shock-persistence tail with an
  exact closed-form stacked autocovariance — material (hundreds of
  log-points) near a unit root. `hank_theta_boundary_check()` diagnoses
  when the solve horizon is too short for a given persistence and warns
  (`check_boundary = FALSE` to silence). `hank_loglik_ar_grad()` supplies
  the analytic/semi-analytic gradient (sigma exact, measurement-error
  exact, persistence semi-analytic) and `make_posterior_grad_hank_ar()`
  wires it into a `make_log_posterior_hank()`-shaped closure. A `cache =`
  argument on `hank_loglik_ar()` / `hank_loglik_aggregate_ar()` memoizes
  per-shock autocovariance slabs (~15x faster warm evaluations);
  `make_log_posterior_hank(likelihood = "exact_ar")` enables it
  automatically. Note: `run_mode_finding()` / `run_full_estimation()` do
  not yet accept `likelihood = "exact_ar"` directly.
- **Kalman state-space bridge**: `hank_ma_state_space()` packs the
  aggregate MA representation into a `dsge_ss` state space, and
  `hank_loglik_ss()` evaluates it through dynhr's ordinary Kalman engine —
  giving a linearized HANK model access to the full estimation stack
  (smoother, frequency-domain likelihoods, priors, samplers, SBC).
- **Wealth heterogeneity axis in mixtures**: `hank_mixture_ks_steady_hetinc()`
  accepts optional per-type `eis` and `amin` (borrowing constraint on the
  shared asset grid), completing the wealth-heterogeneity path through
  general equilibrium (`hank_mixture_block_spec()`, `hank_mixture_dist()`);
  entries without `eis`/`amin` are unaffected.

## Diagnostics and robustness

- **`hank_determinacy()`** reports the conditioning of the GE Jacobian
  `H_U` and its determinacy verdict; `hank_model_irf()` warns on a
  near-singular `H_U` instead of returning a silent garbage IRF.
- **`hank_impc()`** — the intertemporal MPC matrix (Auclert-Rognlie-Straub
  intertemporal Keynesian cross) via fake-news. **`hank_mpc()`** now also
  reports MPC `by_income`. **`hank_distribution_stats()`** /
  **`hank_gini()`** — wealth/consumption Gini, top wealth shares,
  hand-to-mouth fraction, wealth percentiles.
- **`validate_hank_block()`** — a preflight diagnostic for a
  `hank_het_block` (or a mixture's block list): grid monotonicity, the
  shared Markov transition contract, forward-operator row-stochasticity,
  distribution validity, policy feasibility, consumption positivity,
  stored-aggregate identity checks, stationary-solver convergence, and
  optional R/C++ backend parity.
- **macOS Accelerate/vecLib compatibility fix**: a Hermitian eigendecomposition
  segfault under vecLib-backed BLAS (affecting the Whittle likelihood and
  Fisher-information paths) is fixed by an internal real-embedding
  eigensolver that is byte-identical on real input.
- Input-contract hardening across the HANK primitives: `hank_stationary_dist()`
  no longer reports false convergence on a degenerate (NaN) iterate,
  `hank_forward_operator()` / `hank_egm_solve()` validate the Markov/grid
  contract, and `hank_aggregate()` rejects dimension mismatches it would
  otherwise silently recycle.

## Other additions

- **`mode_trust_region()`** — a deterministic dogleg trust-region Newton
  posterior-mode optimizer (with an eigenvalue-clamp Hessian modification),
  complementing the stochastic/global finders and the `newrat` default.
- **`new_obc_regime_path()` / `as_obc_regime_path()`** — a canonical
  occasionally-binding-constraint regime-path object unifying the varied
  return shapes of dynhr's OBC solvers.
- Bugfix: `make_posterior_grad(likelihood = "cumulant")` used the wrong
  solve order (and consequently the wrong sign on some parameters) when
  `cumulant_orders` excluded 3 and 4; the default `cumulant_orders = 1:4`
  was unaffected. `make_posterior_grad()` now errors, instead of silently
  returning a Gaussian gradient, on likelihoods with no gradient path
  (`pskf`, `student_t`, `tpf`, `ppf`, `copf`).
- Foundational HANK / sequence-space Jacobian module (Auclert, Bardóczy,
  Rognlie & Straub 2021): household primitives `hank_income_rouwenhorst()`,
  `hank_asset_grid()`, `hank_egm_solve()`, `hank_euler_residual()`,
  `hank_forward_operator()`, `hank_stationary_dist()`, `hank_aggregate()`;
  the fake-news Jacobian `hank_het_jacobian()` (verified against the
  brute-force `hank_het_jacobian_nd()` to machine precision); GE via
  `hank_ks_steady()`, `hank_ks_linear_irf()`, `hank_ks_nonlinear_irf()`;
  the estimation bridge `hank_ma_coefficients()`, `hank_autocov()`,
  `hank_loglik_aggregate()`, `hank_simulate_aggregate()`; and diagnostics
  `hank_mpc()`, `hank_plot_jacobian()`, `hank_plot_irf()`.

# dynhr 0.8.1

This is the **first public release**. dynhr parses Dynare-style `.mod` files,
solves DSGE models by perturbation (orders 1-3) or global projection, and
estimates them with a wide range of (Bayesian) likelihood-based methods,
backed by an extensive model-diagnostics suite. Converted from a source-script
toolkit into a proper R package: `R CMD INSTALL`/`R CMD check` clean, parser
and solver behavior parity-tested against Dynare/Dynare.jl, vignettes
(`vignette("solving")`, `vignette("diagnostics")`), and a benchmarking harness
(`inst/benchmarks/`).

## Model solving

- Perturbation solving to third order, including `steady_state_model` blocks
  with derived parameters (analytic first- and second-order sensitivity when
  well-conditioned, with a finite-difference fallback) and pruned state
  spaces at order 2/3. `compile_model(param_deriv = "auto"/"on"/"off"/
  "second")` controls how much analytic parameter-derivative machinery is
  compiled (`"second"` enables the exact-Hessian codegen; `"off"` uses the FD
  fallback).
- A **global/projection solver**: Chebyshev collocation with Coleman time
  iteration and Gauss-Hermite expectations, for models with strong
  nonlinearity or occasionally-binding behaviour that perturbation cannot
  capture.
- **Markov-switching DSGE**: shock-variance switching (`ms_dsge_spec`,
  `ms_kim_filter` — Kim-Nelson GPB(2) filter with a Hamilton collapse) and a
  structural switching solver (`solve_ms_perturbation`, regime-coupled
  first-order perturbation).
- `simulate_model()` stochastically simulates a solved model at its compiled
  order (`simulate_model_order2()` / `simulate_model_order3()` for the
  pruned second-/third-order recursions); `compute_irfs()` computes impulse
  responses. `compute_moments_order2()` — deterministic order-2 unconditional
  moments (mean, full covariance, autocorrelations, per-shock variance
  decomposition) from the augmented pruned state space, without simulation
  noise.

## State-space filtering and likelihoods

- **`kalman_filter()`**: `method = "standard"/"dare"/"chandrasekhar"/
  "univariate"` (the last a Koopman-Durbin state-augmented sequential filter
  that handles singular innovation covariances by construction and is the
  automatic fallback whenever another method's `F` is singular or explodes),
  exact-diffuse initialization for unit-root/local-level models
  (`lik_init = "auto"/"diffuse"/"kappa"`), and `return_ll_contrib = TRUE` for
  per-period log-likelihood contributions. **`me_variance` now defaults to
  `0`** (previously `1e-8`) now that the univariate filter removes the need
  for a jitter crutch. **`kalman_smoother()`** returns per-period
  `filtered_cov`/`predicted_cov`/`smoothed_cov` state covariance arrays.
  `build_dsge_state_space()` constructs the `dsge_ss` object from a solved
  model; `new_dsge_ss()` / `ss_convert_timing()` give it an explicit
  lagged/current timing tag and convert between the two conventions.
- **Likelihood families**, all reachable via `make_log_posterior(likelihood
  = )` / `estimation_context()`: Gaussian; `"whittle"` (exact complex-spectral
  multivariate likelihood with analytic gradient, `debias = TRUE` by default
  — the Sykulski et al. 2019 expected-periodogram correction); `"cumulant"`
  (orders up to 4, with an analytic gradient); `"pskf"` — the Closed
  Skew-Normal Kalman filter (Guljanov, Mutschler & Trede 2026) for models
  with skew-normal shocks (`skew SHOCKNAME = expr` in the `shocks;` block;
  `pskf_smoother(method = "csn")` for the backward pass); `"tpf"` — the
  tempered particle filter (C++ `tpf_run_period_cpp` backend,
  `tpf_loglik_sd_preflight()` for a pre-run SD/particle-count check,
  `tpf_options` on `run_full_estimation()`); and the OBC (occasionally-
  binding-constraint) family — `kalman_filter_obc()` / `kalman_filter_obc_pkf()` /
  `kalman_filter_obc_inversion()` (exactly-identified deterministic
  inversion filter) and `ppf_likelihood()` / `make_log_posterior_obc_ppf()`
  (the bootstrap Piecewise Particle Filter), with a conditionally-optimal
  proposal (`proposal = c("bootstrap", "copf")`) and post-hoc importance
  reweighting (`ppf_reweight_posterior()`). `kalman_filter_student_t()`
  gives multivariate-t measurement/forecast errors for heavy-tailed data.
  A `power_posterior` tempering exponent is wired through all five main
  likelihood factories (Gaussian, cumulant, Whittle, PSKF, TPF).
- Shocks-block expressions (`stderr_expr`, `variance_expr`, `corr_expr`,
  `skew_expr`) are all re-evaluated against the current parameter vector on
  every draw, so estimated stderrs/correlations/skewness track theta rather
  than freezing at their parse-time snapshot.
- Correlated-pseudo-marginal proposals for particle likelihoods:
  `rwmh_cpm()` and `tpf_options$cpm_rho_u`.

## Analytic gradients and Hessians

- Tangent- and adjoint-mode Kalman-filter gradients (driving
  `make_posterior_grad()`) support time-varying `me_extra` (per-period
  measurement-error additions) and `shock_scale` (per-period shock-std
  scaling), and — for `steady_state_model`-derived parameters — the full
  first- and second-order sensitivity chain through the derived parameter.
  Analytic gradients are also available for the Whittle and cumulant
  likelihoods.
- **`posterior_hessian()`** — the exact second-order posterior Hessian,
  gated behind `compile_model(param_deriv = "second")`.
- Performance-sensitive gradient/Hessian kernels have C++ ports used
  automatically when available (`dynhr.use_rcpp` option controls the
  Kalman steady-state C++ path; set `FALSE` for the bit-exact R fallback).

## Bayesian estimation and sampling

- Samplers: RWMH, NUTS, SMC, and DIME, plus geometry-aware additions —
  `dynhr_mala()` (Laplace-MALA / simplified-manifold MALA), dense
  (non-diagonal) mass-matrix HMC/NUTS with `softabs_metric()` for indefinite
  Hessians, a gradient-only `monge_metric_fn()`, `whittle_fim()` (a
  frequency-domain Fisher-information metric), and `dynhr_chees()`
  (ChEES-HMC adaptive trajectory length, a tree-free alternative to NUTS).
- **`estimation_context()`** bundles the per-estimation options (likelihood,
  `lik_init`, `me_extra`, `shock_scale`, `freq_band`, `system_priors`,
  `tpf_options`, gradient policy); `run_mode_finding()`,
  `run_full_estimation()`, and `run_posterior_estimation()` drive mode-
  finding and sampling on top of it, with parallel chains/particles via
  `mirai` and live progress bars (`progress = ` on `run_mcmc_mirai()` /
  `run_mode_mirai()`).
- **`dynhr_plan()`** (an IRIS-"plan"-style judgment object) bundles
  `plan_tune()` (in-sample hard/soft tunes), `plan_condition()`
  (out-of-sample Waggoner-Zha conditioning), and `plan_scale_shock()`
  (heteroskedastic shock scaling); `conditional_forecast()` and
  `bayesian_conditional_forecast()` (posterior-draw conditional forecasting)
  both accept a `plan =`.
- GMM estimation with an analytic block-diagonal optimal weight matrix
  (alongside a Newey-West HAC estimator); `score_forecast()` for
  proper-scoring-rule forecast evaluation (CRPS, energy score, variogram
  score).
- Marginal likelihood / model comparison: `thames_mdd()` (the truncated
  harmonic-mean THAMES estimator, Metodiev et al. 2024) and
  `smc_model_tempered()` (two-stage SMC model tempering, Mlikota-Schorfheide
  2024).
- **`robust_confidence_set()`** — weak-identification-robust inference via
  the Andrews-Mikusheva (2015) LM/score test with test inversion.
- **`dynhr_sbc()`** — simulation-based calibration across
  `sampler = c("rwmh", "nuts", "smc", "dime")` and the likelihood families
  above (`nuts` is refused with the noisy `"tpf"` likelihood).

## Diagnostics

- **`run_all_diagnostics()`** runs a broad battery (D0 through D36),
  including: `$d0`, a static-Jacobian rank preflight for local
  well-posedness; variance/historical decomposition and smoothed shocks at
  the posterior mean (with opt-in `bayesian_irf = TRUE` posterior IRF
  credible bands); spectral identification with the full complex Hermitian
  Gram matrix (`spectral = c("exact", "companion")`); and system priors
  including `sp_spectral_peak()` (a spectral-density-peak-frequency prior
  feature).
- **Deep-parameter diagnostics**: an `@dynhr:deep` `.mod` metadata block
  classifies each parameter (`extract_mod_metadata()`, `build_deep_spec()`)
  as a deep primitive or reduced-form/auxiliary quantity; four new
  diagnostics assess structural-vs-reduced-form confounding
  ("borrowed identification"), policy-partitioned invariance (an
  operational Lucas-critique test, with cheap Laplace draws via
  `deep_laplace_draws()`), misspecification softness (a sandwich/
  information-matrix-equality check), and calibration deepness/tension for
  parameters that are calibrated rather than estimated. **
  `deep_parameter_passport()`** synthesizes all of this into a per-parameter
  A-F scorecard. `kalman_filter(return_ll_contrib = TRUE)` supplies the
  per-period likelihood contributions these diagnostics consume.
- `@dynhr:expectations` `.mod` metadata declares model-specific checks
  (`data_mean`, `data_sd`, `data_ratio`, `param_range`, `irf_sign`),
  evaluated by `diag_expectations()`. `run_diagnostics()` / `write_report()`
  are the public aliases for `run_all_diagnostics()` / `write_llm_report()`.

## Welfare and optimal policy

- `ramsey_model(order = 2)` — full second-order (augmented-system) Ramsey
  optimal policy. `osr()` (optimal simple rules) gains `order = ` and
  `planner_objective = ` to minimize expected welfare directly at order 2
  instead of an order-invariant quadratic variance loss.
  `conditional_welfare()` supports `method = "analytic"` (fast deterministic Taylor
  approximation), `"stochastic"` (Monte Carlo), or `"deterministic"` (the
  zero-shock path) for the order-2 state-dependent risk correction, and
  correctly conditions on a supplied `initial_state`.
