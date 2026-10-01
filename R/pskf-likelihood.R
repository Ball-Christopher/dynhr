## R/pskf-likelihood.R
## --------------------------------------------------------------------------
## Pruned Skewed Kalman Filter (PSKF) likelihood for first-order DSGE models
## with skew-normal shocks.
##
## Reference: Guljanov, Mutschler & Trede (2026), "Pruned Skewed Kalman Filter
##   and Smoother with Application to DSGE Models," JEDC Vol. 187 (Dynare WP
##   #78). Reference implementation: github.com/gguljanov/pruned-skewed-kalman
##   (R_codes/skalman_filter.R).
##
## The shocks follow a Closed Skew-Normal (CSN) distribution:
##   x ~ CSN(mu, Sigma, Gamma, nu, Delta)
##   pdf ∝ N(x; mu, Sigma) * Phi_q(Gamma(x-mu) - nu; 0, Delta)
##            / Phi_q(-nu; 0, Delta + Gamma Sigma Gamma')
## Gamma = 0 collapses to Gaussian.  Closed under affine maps and conditioning,
## which enables an exact Kalman-filter analog.
##
## CRITICAL CONVENTIONS:
##  1. ghu EXCLUDES Sigma_e (it is the loading matrix R in x_t = G x_{t-1} + R e_t).
##     Sigma_eta = RR Sigma_e RR'; Sigma_eps = DD Sigma_e DD' + me_variance*I.
##  4. Sigma_eps is assembled as ONE matrix (HH + me_variance*I).
##  6. mu_eta is mean-corrected so E[eta] = 0: mu_eta = -RR E[e], with E[e] the
##     mean of the JOINT CSN shock law (.csn_shock_mean). For uncorrelated
##     skewed shocks E[e_i] = sigma_i delta_i sqrt(2/pi); under correlation the
##     per-shock formula is wrong (fixed 0.9.4 -- it used to be applied always).
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Deterministic multivariate-normal CDF
## ---------------------------------------------------------------------------
## logcdf_ME_r(x, S) computes log Phi_q(x; 0, S).
##
## Strategy:
##   q = 1: pnorm(log.p = TRUE) -- exact.
##   q = 2: .mvn_logcdf2() -- Genz's BVND in C++ when p >= 1e-3 (0.9.4; <= ~1e-13
##          in log p), otherwise the conditional 1-D integral on the LOG scale
##          around its (log-concave) mode; relative accuracy ~1e-10 down to
##          log p ~ -1e3.
##   q = 3 (0.9.4): mvn_logcdf3_cpp() -- Genz's (2004) exact
##          trivariate normal (Plackett reduction, adaptive Gauss-Kronrod;
##          |error in p| <= ~1e-15) when p >= 1e-3, det(correlation) >= 1e-8
##          and the adaptive rule converged; otherwise the lattice below.
##          Same function, same inputs, in the C++ and the R dispatch.
##   3 <= q <= miwa_qmax: .mvn_logcdf_sov() -- the C++ evaluator
##          mvn_logcdf_cpp() (src/mvn_cdf.cpp): Genz separation of variables
##          with Genz-Bretz variable reordering and Botev minimax tilting,
##          on fixed rank-1 lattices with fixed shifts (no RNG), log scale,
##          lattice size raised until the shift-spread error estimate is
##          below 0.1-0.2 x max(1e-5, 1e-7 |log p|) (continuous hand-over
##          between sizes, see mvn_cdf.cpp). Since 0.9.3.122; measured
##          against independent 1- and 2-factor quadrature oracles (160
##          problems, q = 3..7, log p in [-1, -300]) and mvtnorm Genz-Bretz
##          (40, q = 3..5): every |error| < max(1e-5, 1e-7 |log p|), worst
##          0.39 of it. PSKF calls cost ~0.1 / 0.3 / 1.4 / 3 / 10 ms at
##          q = 3 / 4 / 5 / 6 / 7 (-O2). It replaced the
##          checked-Miwa / R-lattice evaluator (0.9.3.118) (1e-3 to 2.5e-2 nat off in
##          moderate-to-deep tails) and plain Miwa(128) for q = 3 and 6-7.
##          Falls back to plain Miwa(128), then Mendell-Elston, only when C
##          is numerically singular (a Cholesky pivot <= 1e-10).
##   The snap-to-zero below uses |rho| < 1e-12 (exact block factorisation):
##   the C++ evaluator has no tiny-correlation pocket, and the 1e-3 snap of
##   the Miwa path costs up to |b_i b_j| * 1e-3 nat in log p in the tails.
##
##   q > miwa_qmax: Mendell-Elston (1974) sequential
##          conditioning approximation. First-order moment-matching;
##          absolute error ~1e-3 (near-diagonal S) to ~0.05 (high
##          correlations); measured against mvtnorm (40 random 3-7 dim
##          problems) median |log-CDF error| 0.016, max 0.86 (before 0.9.4.6:
##          median 1.1, max 6.5, from a sign / double-shrink bug). The checked
##          lattice / Genz evaluators above remain far more accurate (~1e-6
##          or better), which is why they are the default. See inline
##          accuracy note below.
##
## Why not GenzBretz?  GenzBretz (mvtnorm default) uses randomised QMC --
## non-reproducible across calls, which breaks gradient-based mode-finding and
## SBC.  Miwa() and the fixed-lattice rule are deterministic.  A tryCatch
## wraps the Miwa call so that non-PD edge-cases fall back to ME without
## crashing.
##
## miwa_qmax: largest dimension evaluated with the deterministic evaluators
## before falling back to Mendell-Elston.  An exact 2-D grid-filter oracle
## showed that the multi-shock "pruning" bias formerly seen (-7.3 nats at
## T=12 for alpha=(+2,-2)) was Mendell-Elston evaluation error at q>5, NOT
## discarded skew mass: swapping ME for an exact Phi_q at unchanged
## cut_tol=0.01 collapsed the gap to |0.005|.  (That error was mostly the
## sign / double-shrink bug in the ME evaluator fixed in 0.9.4.6; with the
## fix, an uncapped max_q = Inf run is 0.01-0.03 nat from the grid.)
## Keeping q inside the exact range via rank-capped pruning (see dim_red4_r
## max_q) remains the default, since the exact evaluators are still more
## accurate than ME.
#' @noRd
logcdf_ME_r <- function(x, S, miwa_qmax = 5L, use_cpp = TRUE) {
  q <- length(x)
  if (q == 0L) return(0)

  ## q = 1: exact
  if (q == 1L) {
    b <- x[1] / sqrt(S[1, 1])
    return(pnorm(b, log.p = TRUE))
  }

  ## Accurate path in C++ (0.9.4): mvn_logcdf_dispatch_cpp() runs
  ## the snap / block factorisation / dispatch below operation for operation
  ## (bit-identical results) and returns NA whenever this R code would reach
  ## a fallback evaluator (bivariate quadrature below p = 1e-3, Miwa,
  ## Mendell-Elston) or meets a non-finite input -- the R path then runs
  ## unchanged. It removes the ~30 us of interpreter overhead per call.
  ## use_cpp = FALSE skips the C++ fast paths (this one, the bivariate Genz
  ## rule in .mvn_logcdf2 and the trivariate one in .mvn_logcdf_sov) -- the
  ## general path, for the equivalence tests
  ## (test-fix-0929-pskf-cdf-fast-path.R, test-fix-0929-pskf-tvn.R).
  if (use_cpp) {
    v_cpp <- mvn_logcdf_dispatch_cpp(as.numeric(x), as.matrix(S),
                                     as.numeric(miwa_qmax))
    if (!is.na(v_cpp)) return(v_cpp)
  }

  ## ---- SNAP + BLOCK FACTORIZATION (Miwa-pocket fix) ---------------
  ## mvtnorm's Miwa algorithm has an instability pocket for TINY-but-nonzero
  ## correlations: measured on a well-conditioned 3x3 with mixed-sign
  ## rho ~ 1e-5..1e-3 it returned Phi_3 = 1.0856 (> 1!) at steps = 128 and was
  ## still 0.04 absolute off at steps = 512, while rho = 0, 1e-7 and 0.2 are
  ## all ~1e-6 accurate (Reiter-HANK PSKF investigation -- this made the
  ## filter's likelihood IMPROPER).  Fix:
  ## (a) standardize to correlation form and SNAP |rho| < 1e-3 to exactly 0
  ##     (error bound: |dPhi/drho| = phi_2 <= 1/(2*pi) per pair, so <= ~1.6e-4
  ##     per snapped pair in PROBABILITY -- but |d log Phi / d rho| grows like
  ##     |b_i b_j| in the tails; hence the accurate path, which does not use
  ##     Miwa, snaps only |rho| < 1e-12);
  ## (b) factor the CDF over the connected components of the snapped
  ##     correlation graph (EXACT given the snap) -- singletons/pairs then use
  ##     the exact pnorm/bivariate paths and the q >= 3 evaluator only sees
  ##     well-coupled blocks.
  sdv <- sqrt(pmax(diag(as.matrix(S)), .Machine$double.eps))
  Cm  <- as.matrix(S) / outer(sdv, sdv)
  off <- row(Cm) != col(Cm)
  Cm[off & abs(Cm) < 1e-12] <- 0
  x <- as.numeric(x) / sdv
  S <- Cm
  adj  <- Cm != 0
  comp <- integer(q); n_comp <- 0L
  if (all(adj)) {
    n_comp <- 1L                   # fully coupled (the usual case): one block
  } else {
    for (s0 in seq_len(q)) {
      if (comp[s0] > 0L) next
      n_comp <- n_comp + 1L
      stack <- s0
      while (length(stack) > 0L) {
        v <- stack[[1L]]; stack <- stack[-1L]
        if (comp[v] > 0L) next
        comp[v] <- n_comp
        stack <- c(stack, which(adj[v, ] & comp == 0L))
      }
    }
  }
  if (n_comp > 1L) {
    out <- 0
    for (cc in seq_len(n_comp)) {
      ii  <- which(comp == cc)
      out <- out + logcdf_ME_r(x[ii], Cm[ii, ii, drop = FALSE], miwa_qmax,
                               use_cpp = use_cpp)
    }
    return(out)
  }

  ## WEAK COUPLINGS. A block of q >= 4 needs the lattice evaluator or
  ## Mendell-Elston, whose errors (>= 1e-7 nat; median 0.016 for ME) dwarf the effect
  ## of a correlation below 1e-8 (|d log Phi / d rho| ~ |b_i b_j|). Couplings
  ## that weak are round-off, not structure: the covariances the PSKF passes
  ## here are built as Delta + Gamma Sigma Gamma' with Gamma entries in the
  ## hundreds, so cross-block entries of 1e-12 .. 1e-11 (correlation scale)
  ## appear where the exact value is zero. One such entry above the 1e-12 snap
  ## joins blocks that are each inside the exact range (measured: a 6-row
  ## stack of blocks 4 + 1 + 1 became one block of 6, went to Mendell-Elston
  ## and gave a mean-offset factor of 0.657 where the closed form for the
  ## independent row is 0.3505; a 4-block of two independent pairs went to
  ## the lattice, whose reordering flips with the round-off). Zero them and
  ## split again; blocks of 2 and 3 keep the 1e-12 snap because their
  ## evaluators are exact to 1e-13. mvn_logcdf_dispatch_cpp applies the same
  ## rule.
  if (q >= 4L) {
    weak <- off & S != 0 & abs(S) < 1e-8
    if (any(weak)) {
      S[weak] <- 0
      return(logcdf_ME_r(x, S, miwa_qmax, use_cpp = use_cpp))
    }
  }

  ## q = 2
  if (q == 2L) return(.mvn_logcdf2(x[1], x[2], S[1, 2], use_cpp = use_cpp))

  ## 3 <= q <= miwa_qmax: the C++ lattice evaluator. Numerically singular S
  ## (C++ NA) and Miwa's impossible values fall through to Miwa / ME.
  if (q <= miwa_qmax) {
    val <- .mvn_logcdf_sov(x, S, use_cpp = use_cpp)
    if (!is.na(val)) return(val)
    if (requireNamespace("mvtnorm", quietly = TRUE)) {
      val <- .mvn_miwa_prob(x, S, 128L)
      if (!is.na(val)) return(log(max(val, .Machine$double.eps)))
    }
    ## Fall through to ME on error / impossible value (non-PD S, etc.)
  }

  ## q > miwa_qmax (or both evaluators failed): Mendell-Elston sequential
  ## conditioning.
  ##
  ## ACCURACY NOTE (Mendell & Elston 1974 first-order moment-matching):
  ## ME approximates Phi_q(x; 0, S) by sequential univariate conditioning:
  ##   Phi_q ≈ prod_{j=1}^{q} Phi_1(b_j(x, S)), where b_j is the conditional
  ##   mean adjustment at each step.
  ## Error bounds (Mendell & Elston 1974; empirical benchmarks vs Miwa):
  ##   - Near-diagonal S (max |rho_ij| < 0.1): absolute error ~1e-3 per period.
  ##   - Moderate correlations (max |rho_ij| ~ 0.3-0.5): absolute error ~0.01.
  ##   - High correlations (max |rho_ij| > 0.8): absolute error ~0.03-0.05.
  ##   These translate to per-period loglik errors of O(abs_err/Phi_q) in the
  ##   CSN loglik correction terms.
  ## MEASURED against mvtnorm (40 random 3-7 dimensional cases): median
  ## absolute log-CDF error 0.016, max 0.86 in the tails. Before 0.9.4.6 the
  ## conditional-mean shift had the wrong sign, which produced errors of
  ## several nats (median 1.1, max 6.5 on the same cases). The branch is
  ## reached with the default max_q = 5 too: the pre-pruning stack of
  ## max_q + n_exo rows can form a coupled block wider than miwa_qmax, and
  ## the lattice / Miwa evaluators can decline a near-singular block.
  ##
  ## Algorithm (Mendell & Elston 1974, corrected for covariance -- not correlation -- form):
  ## For j = 1 .. q-1:
  ##   standardised bound: bj = b[j] / sqrt(S[j,j])
  ##   P_j = Phi(bj)  (marginal probability for dim j)
  ##   Mills ratio lambda = phi(bj)/Phi(bj)
  ##   mean shift for k > j:  b[k] += S[j,k]/sqrt(S[j,j]) * lambda
  ##     (E[X_k | X_j <= b_j] = -S[j,k]/sqrt(S[j,j]) * lambda, so the bound
  ##     for X_k rises; the sign was flipped before 0.9.4.6)
  ##   variance shrinkage for k,l > j:
  ##     S[k,l] -= S[j,k]*S[j,l]/S[j,j] * lambda*(bj + lambda)
  ## Final: P_q = Phi(b[q]/sqrt(S[q,q]))
  ## Product = prod P_j * P_q

  b  <- as.numeric(x)
  CS <- as.matrix(S)
  log_p <- 0

  for (j in seq_len(q - 1L)) {
    sjj <- CS[j, j]
    sj  <- sqrt(max(sjj, .Machine$double.eps))
    bj  <- b[j] / sj
    lPj <- pnorm(bj, log.p = TRUE)
    log_p <- log_p + lPj

    lphi_j       <- dnorm(bj, log = TRUE)
    lambda        <- exp(lphi_j - lPj)         # phi(bj) / Phi(bj)
    delta_factor  <- lambda * (bj + lambda)     # variance shrinkage factor

    idx <- (j + 1L):q
    for (k in idx) {
      ## Mean shift
      b[k] <- b[k] + CS[j, k] / sj * lambda
      ## Covariance shrinkage: the upper triangle once, then mirrored (looping
      ## l over all of idx decremented every off-diagonal twice before 0.9.4.6)
      for (l in idx[idx >= k]) {
        CS[k, l] <- CS[k, l] - CS[j, k] * CS[j, l] / sjj * delta_factor
        CS[l, k] <- CS[k, l]
      }
    }
  }

  bq    <- b[q] / sqrt(max(CS[q, q], .Machine$double.eps))
  log_p <- log_p + pnorm(bq, log.p = TRUE)

  log_p
}


## Phi_q(x; 0, S) by mvtnorm's Miwa rule with `steps` grid points, or NA on an
## error or an impossible value. IMPOSSIBLE-VALUE guard (Miwa-pocket fix,
## part c): Miwa can return probabilities > 1 or <= 0 WITHOUT erroring on
## pathological inputs -- those are treated like an error so the caller falls
## through instead of corrupting the loglik.
#' @noRd
.mvn_miwa_prob <- function(x, S, steps) {
  val <- tryCatch(
    mvtnorm::pmvnorm(upper = as.numeric(x), sigma = as.matrix(S),
                     algorithm = mvtnorm::Miwa(steps = steps))[1L],
    error = function(e) NA_real_
  )
  if (is.na(val) || !is.finite(val) || val <= 0 || val > 1 + 1e-8)
    return(NA_real_)
  min(val, 1)
}


## log Phi_2(h1, h2; rho) for standardized bounds, on the LOG scale:
##   Phi_2 = int_{-Inf}^{h1} phi(u) Phi((h2 - rho u) / sr) du,  sr^2 = 1 - rho^2,
## with h1 the smaller bound. The log-integrand g(u) is strictly concave
## (g'' <= -1), so the mass sits within ~40 of its maximiser um on
## (-Inf, h1] (g <= g(um) - (u - um)^2 / 2 away from it): exp(g - g(um)) is
## integrated over [um - 40, min(h1, um + 40)], split at um +- 8, and g(um)
## added back. No lower limit of -8 (the pre-0.9.3.122 path lost the mass below it
## -- all of it for h1 = 0, h2 = -10, rho = -0.9) and no floor at
## log(.Machine$double.eps). |rho| is capped at 0.9999 as before.
## When p >= 1e-3 the quadrature is skipped for Genz's BVND
## (mvn_logcdf2_cpp, src/mvn_cdf.cpp): the same value to <= ~1e-13 in log p
## (6.6e-14 measured, the quadrature's own tolerance) at ~1/900 of the cost.
#' @noRd
.mvn_logcdf2 <- function(h1, h2, rho, use_cpp = TRUE) {
  if (is.na(h1) || is.na(h2) || is.na(rho)) return(NA_real_)
  if (h1 == -Inf || h2 == -Inf) return(-Inf)
  if (h1 == Inf) return(pnorm(h2, log.p = TRUE))
  if (h2 == Inf) return(pnorm(h1, log.p = TRUE))
  if (h2 < h1) { tmp <- h1; h1 <- h2; h2 <- tmp }
  rho <- max(-0.9999, min(0.9999, rho))
  ## Fast exact path (0.9.4): Genz's BVND in C++ (mvn_logcdf2_cpp),
  ## absolute error ~1e-16 in p, returned only for p >= 1e-3 where that is
  ## <= ~1e-13 in log p -- below the quadrature's own rel.tol. NA (deeper
  ## tail) falls through to the log-scale quadrature.
  if (use_cpp) {
    v_cpp <- mvn_logcdf2_cpp(h1, h2, rho)
    if (!is.na(v_cpp)) return(v_cpp)
  }
  sr  <- sqrt(1 - rho^2)
  g   <- function(u) dnorm(u, log = TRUE) +
                     pnorm((h2 - rho * u) / sr, log.p = TRUE)
  gp  <- function(u) {
    z <- (h2 - rho * u) / sr
    -u - rho / sr * exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))
  }
  if (gp(h1) >= 0) {
    um <- h1
  } else {
    lo <- h1 - 1
    while (gp(lo) < 0) lo <- h1 - 2 * (h1 - lo)
    um <- stats::uniroot(gp, c(lo, h1), tol = 1e-12)$root
  }
  gm  <- g(um)
  f   <- function(u) exp(g(u) - gm)
  br  <- c(um - 40, um - 8, um, um + 8, um + 40)
  br  <- pmin(br, h1)
  tot <- 0
  for (k in 1:4) {
    if (br[k + 1L] > br[k])
      tot <- tot + stats::integrate(f, br[k], br[k + 1L], rel.tol = 1e-11,
                                    abs.tol = 0, subdivisions = 200L)$value
  }
  gm + log(tot)
}


## log Phi_q(b; 0, C) by the C++ lattice evaluator mvn_logcdf_cpp()
## (src/mvn_cdf.cpp: Genz separation of variables, Genz-Bretz reordering,
## Botev minimax tilting, fixed rank-1 lattices and shifts, log scale).
## Deterministic; NA when C is numerically singular (a Cholesky pivot
## <= 1e-10 in correlation form). Bounds of +Inf are marginalised out.
#' @noRd
.mvn_logcdf_sov <- function(b, C, use_cpp = TRUE) {
  b <- as.numeric(b)
  C <- as.matrix(C)
  if (anyNA(b) || anyNA(C)) return(NA_real_)
  if (any(b == -Inf)) return(-Inf)
  fin <- b < Inf
  if (!all(fin)) {
    if (!any(fin)) return(0)
    b <- b[fin]
    C <- C[fin, fin, drop = FALSE]
  }
  if (length(b) == 1L) return(pnorm(b / sqrt(C[1L, 1L]), log.p = TRUE))
  if (length(b) == 2L) {
    sd <- sqrt(diag(C))
    return(.mvn_logcdf2(b[1L] / sd[1L], b[2L] / sd[2L], C[1L, 2L] / (sd[1L] * sd[2L]),
                        use_cpp = use_cpp))
  }
  ## q = 3 (0.9.4): Genz's exact TVN (mvn_logcdf3_cpp) when its
  ## threshold rule accepts -- p >= 1e-3, det(correlation) >= 1e-8, adaptive
  ## rule converged; NA otherwise, and the lattice runs. The C++ dispatch
  ## calls the same function on the same inputs (bit-identical).
  if (length(b) == 3L && use_cpp) {
    v3 <- mvn_logcdf3_cpp(b, C)
    if (!is.na(v3)) return(v3)
  }
  mvn_logcdf_cpp(b, C)[1L]
}


## ---------------------------------------------------------------------------
## Pruning: dim_red4_r
## ---------------------------------------------------------------------------
## Keeps only the rows of (Gamma, nu) whose maximum absolute correlation with
## any other row of Gamma Sigma Gamma' + Delta is >= cut_tol.
## Without pruning the skewness dimension q grows by n_exo every period;
## pruning keeps q bounded (typically 1-3 for small DSGE).
##
## max_q: HARD rank-based cap on the retained skew dimension.  After the
## cut_tol threshold filter, if more than max_q rows survive, only the max_q
## rows with the LARGEST skew-vs-state correlation are kept (in original
## order).  Rationale: the Phi_q evaluator is deterministic-accurate to
## ~1e-6 only for q <= miwa_qmax (5 in the filter's likelihood terms; the C++
## lattice evaluator since 0.9.3.122, Miwa before); beyond that the
## Mendell-Elston approximation is used, whose per-period error (median
## 0.016 nat after the 0.9.4.6 fix) accumulates over periods.  Rank-capping q
## at 5 keeps every CDF call inside the exact-evaluator range; the discarded
## low-correlation skew mass costs at most ~0.3 nat at T=12 on the worst
## measured fixture.  Historical note: the -7.3 nat multi-shock bias
## (T=12, 2 shocks, alpha=(+2,-2)) that motivated the cap was the ME
## sign / double-shrink bug fixed in 0.9.4.6; uncapped (max_q = Inf) the
## error is now 0.01-0.03 nat from an exact grid on that fixture, so the
## cap is a precision / cost safeguard rather than a bias correction.
##
## Mean offset of a CSN(0, Sigma, Gamma, nu, Delta) relative to its Gaussian
## location:  E[X] - mu = Sigma Gamma' g,
##   g_j = phi(-nu_j; V_jj) * Phi_{q-1}(cond_j) / Phi_q(-nu; V),
##   V = Delta + Gamma Sigma Gamma',
## (gradient of the log-normaliser wrt nu; verified against rejection-sampling
## MC). miwa_qmax is the range of
## the Phi evaluations (logcdf_ME_r); Mendell-Elston beyond it.
## DEFAULT 5 since 0.9.3.118 (it was 2, i.e. Mendell-Elston for every Phi of
## dimension >= 3, on the argument that the errors cancel in the ratio -- they
## do not). Measured against INDEPENDENT oracles (exact grid filter;
## Rao-Blackwellised particle filter on the half-normal selection latents).
## The former ME errors below were dominated by the sign / double-shrink bug
## fixed in 0.9.4.6; "ME" is the fixed Mendell-Elston evaluator
## (miwa_qmax = 2), "before" the buggy one:
##   skewed AR(1), alpha 1.5, T = 20 (test-pskf-smoother.R fixture):
##     ME 3e-4 nat off the grid likelihood (before 0.138), miwa_qmax 5: 2e-4;
##   2 states, 2 skew shocks alpha (3, -2), measurement-error var 1e-3,
##     T = 100: ME 0.002 nat off the particle filter (before -18.4),
##     miwa_qmax 5: 0.07 MCSE;
##   1 state, 2 skew shocks (4, -3), T = 60: before -2.3 nat, miwa_qmax 5 -0.07.
## The accurate evaluators are used throughout: they are ~1e-6 accurate per
## call versus ME's 0.016 median (0.86 max) log-CDF error.
## The smoother passes 7 (dim_red4_r / .pskf_filter offset_miwa_qmax): there
## the compensation is carried back to earlier periods (the buggy
## ME-evaluated g put t < T smoothed means up to 2.4 posterior sd off).
#' @noRd
.csn_mean_offset <- function(Gamma, nu, Delta, Sigma) {
  as.numeric(Sigma %*% t(Gamma) %*% .csn_offset_g(Gamma, nu, Delta, Sigma))
}

## The latent-space factor g of .csn_mean_offset (offset = Sigma Gamma' g);
## all zeros when any entry is non-finite (the offset is then dropped, as
## before).
#' @noRd
.csn_offset_g <- function(Gamma, nu, Delta, Sigma, miwa_qmax = 5L) {
  q <- nrow(Gamma)
  if (q == 0L) return(numeric(0))
  V <- Delta + Gamma %*% Sigma %*% t(Gamma)
  V <- (V + t(V)) / 2
  logZ <- logcdf_ME_r(-nu, V, miwa_qmax = miwa_qmax)
  g <- numeric(q)
  for (j in seq_len(q)) {
    Vjj <- V[j, j]
    if (Vjj <= 0) next
    lphi <- dnorm(-nu[j], 0, sqrt(Vjj), log = TRUE)
    if (q == 1L) {
      lcond <- 0
    } else {
      mcond <- (-nu[-j]) - V[-j, j] / Vjj * (-nu[j])
      Scond <- V[-j, -j, drop = FALSE] -
               V[-j, j, drop = FALSE] %*% V[j, -j, drop = FALSE] / Vjj
      Scond <- (Scond + t(Scond)) / 2
      lcond <- logcdf_ME_r(mcond, Scond, miwa_qmax = miwa_qmax)
    }
    g[j] <- exp(lphi + lcond - logZ)
  }
  if (!all(is.finite(g))) return(rep(0, q))
  g
}

## Returns list(Gamma = q'xp, nu = q', Delta = q'xq', mu_shift = p-vector,
## keep = integer(q'), lambda = numeric(q)) with q' <= q. keep holds the
## retained rows' indices in the input stack (increasing); lambda is the cut
## expressed in the latent space of the input stack,
## lambda = g_before - (g_after placed at the keep rows), so that
## mu_shift = Sigma Gamma' lambda = Cov(x, W) lambda. The smoother
## (pskf_smoother) needs both: keep to track each retained skew latent back
## to its birth period, lambda to carry the cut's compensation back to the
## periods where the dropped latents were still alive. offset_miwa_qmax is
## the deterministic-evaluator range of the Phi evaluations inside the
## compensation (.csn_offset_g, logcdf_ME_r; default 5).
## mu_shift is the FIRST-MOMENT COMPENSATION for the cut:
## deleting a skew row removes that dimension's contribution to the CSN mean
## (first-order in its skew-state correlation), so the caller must add
## mu_shift = offset(before) - offset(after) to the Gaussian location.
## Without it, saturated pruning (rank cap binding every period, e.g. a
## persistent single-shock model with |alpha| large) accumulates a systematic
## state-mean drift that makes the likelihood IMPROPER (one-step predictive
## densities integrating to 0.03-0.97) -- caught by the alpha_z SBC on the
## Reiter HANK state space (coverage collapsed to 3-6%).
#' @noRd
dim_red4_r <- function(Gamma, nu, Delta, Sigma, cut_tol = 0.01, max_q = 5L,
                       offset_miwa_qmax = 5L, keep = NULL) {
  q <- nrow(Gamma)
  no_shift <- rep(0, ncol(Gamma))
  if (q == 0L)
    return(list(Gamma = Gamma, nu = nu, Delta = Delta, mu_shift = no_shift,
                keep = integer(0), lambda = numeric(0)))
  ## keep (optional): the retained rows, increasing indices into the stack.
  ## The selection below is then skipped and these rows are used as given;
  ## the mean compensation for the cut is computed exactly as for a selected
  ## set. This is the frozen-selection branch that finite-difference
  ## gradients evaluate (see .pskf_filter, keep_override).
  if (is.null(keep)) {
    keep_idx <- .dim_red4_select(Gamma, Delta, Sigma, cut_tol, max_q)
  } else {
    keep_idx <- as.integer(keep)
  }
  if (length(keep_idx) == q) {
    return(list(Gamma = Gamma, nu = nu, Delta = Delta, mu_shift = no_shift,
                keep = seq_len(q), lambda = rep(0, q)))
  }

  SGt      <- Sigma %*% t(Gamma)
  g_before <- .csn_offset_g(Gamma, nu, Delta, Sigma,
                            miwa_qmax = offset_miwa_qmax)
  if (length(keep_idx) == 0L) {
    return(list(
      Gamma    = matrix(0, nrow = 0L, ncol = ncol(Gamma)),
      nu       = numeric(0),
      Delta    = matrix(0, nrow = 0L, ncol = 0L),
      mu_shift = as.numeric(SGt %*% g_before),
      keep     = integer(0),
      lambda   = g_before
    ))
  }
  G_k <- Gamma[keep_idx, , drop = FALSE]
  n_k <- nu[keep_idx]
  D_k <- Delta[keep_idx, keep_idx, drop = FALSE]
  lambda <- g_before
  lambda[keep_idx] <- lambda[keep_idx] -
    .csn_offset_g(G_k, n_k, D_k, Sigma, miwa_qmax = offset_miwa_qmax)
  list(
    Gamma    = G_k,
    nu       = n_k,
    Delta    = D_k,
    mu_shift = as.numeric(SGt %*% lambda),
    keep     = as.integer(keep_idx),
    lambda   = lambda
  )
}

## The kept rows of dim_red4_r's natural selection (increasing indices into
## the stack): the cut_tol threshold, the collinearity guard, the max_q cap.
#' @noRd
.dim_red4_select <- function(Gamma, Delta, Sigma, cut_tol, max_q) {
  q <- nrow(Gamma)
  ## Pruning criterion (reference dim_red4): the correlation between each
  ## skewness dimension and the STATE,
  ##   corr(skew_i, x_j) = (Gamma Sigma)[i, j] /
  ##                       sqrt((Delta + Gamma Sigma Gamma')[i, i] * Sigma[j, j])
  ## A skew dimension whose correlation with every state is ~0 contributes an
  ## (almost) constant Phi factor to the density -- it cancels between the
  ## numerator and normaliser and can be dropped.  Keep rows with
  ## max_j |corr| >= cut_tol.
  ## (A previous version measured correlations AMONG skew rows, which dropped
  ## the dominant incoming-shock dimension every period -- caught by an exact
  ## grid-filter oracle.)
  GS       <- Gamma %*% Sigma                                   # q x n_state
  cov_full <- Delta + GS %*% t(Gamma)                            # q x q
  d_skew   <- sqrt(pmax(diag(cov_full), .Machine$double.eps))
  d_state  <- sqrt(pmax(diag(Sigma),    .Machine$double.eps))
  corr_sx  <- abs(GS) / outer(d_skew, d_state)                   # q x n_state
  ## row maxima (max.col picks the argmax column; NaN rows are not ordered
  ## by it, so a non-finite matrix takes the plain apply() route)
  max_corr <- if (all(is.finite(corr_sx)) && ncol(corr_sx) > 0L)
    corr_sx[cbind(seq_len(q), max.col(corr_sx, ties.method = "first"))]
  else apply(corr_sx, 1, max)

  keep_idx <- which(max_corr >= cut_tol)
  ## COLLINEARITY guard (Reiter-HANK investigation): when the skew
  ## rows are propagated images of the SAME shock direction (persistent
  ## single-shock models), corr among skew dims -> 1 within a few periods.
  ## Near-singular V = Delta + G S G' breaks the Phi_q evaluators (Miwa
  ## errors -> ME fallback whose ~0.03 ABSOLUTE error on ~1e-8 tail CDFs is
  ## large relative to the CDF itself, hence in the top-minus-bottom log
  ## difference -> IMPROPER likelihood, one-step predictive integrals
  ## 0.04-1e14; measured with the ME evaluator before its 0.9.4.6 sign fix). Nearly-duplicate
  ## constraints are nearly-free to drop UNDER MEAN COMPENSATION (below), so
  ## iteratively drop the weaker row of any pair with |corr| > 0.9 -- this
  ## keeps V numerically nonsingular and the deterministic evaluators (the
  ## C++ lattice rule since 0.9.3.122) away from their singular-C fallback.
  if (length(keep_idx) > 1L) {
    repeat {
      Vk <- cov_full[keep_idx, keep_idx, drop = FALSE]
      dk <- d_skew[keep_idx]
      Ck <- abs(Vk) / outer(dk, dk)
      diag(Ck) <- 0
      mx <- which(Ck == max(Ck), arr.ind = TRUE)[1L, ]
      if (Ck[mx[1L], mx[2L]] <= 0.9 || length(keep_idx) <= 1L) break
      pair <- keep_idx[c(mx[1L], mx[2L])]
      drop_row <- pair[which.min(max_corr[pair])]
      keep_idx <- setdiff(keep_idx, drop_row)
    }
  }
  ## Rank-based cap: keep the max_q rows with the largest skew-vs-state
  ## correlation (original row order preserved for reproducibility).
  if (is.finite(max_q) && length(keep_idx) > max_q) {
    keep_idx <- keep_idx[order(max_corr[keep_idx],
                               decreasing = TRUE)[seq_len(max_q)]]
    keep_idx <- sort(keep_idx)
  }
  keep_idx
}


## ---------------------------------------------------------------------------
## Assemble CSN shock parameters from model and current params
## ---------------------------------------------------------------------------
## Returns list with (n_chi = dim of the contemporaneous state, see
## .pskf_order1_statespace):
##   Sigma_eta : n_chi x n_chi     = RR Sigma_e RR'  (RR = ghu[chi_idx, ])
##   Sigma_eps : n_obs x n_obs     = me_variance * I (or diag(me)); DD = 0
##   Gamma_eta : q_eta x n_chi     (skewness loading through RR)
##   nu_eta    : q_eta vector (= 0)
##   Delta_eta : q_eta x q_eta (= I)
##   mu_eta    : n_chi vector (mean correction, = 0 when all alpha_i = 0)
##   mu_eps    : n_obs vector   (= 0)
##   alpha     : n_exo named vector of shape parameters
##   TT, ZZ, chi_idx, state_pos, obs_pos : the state space to filter with
##
## DERIVATION (CSN linear-transform convention; file-header convention 1):
##   The shock vector in model space is e ~ skewNormal(alpha_i, sigma_i) for
##   each component i.  Each scalar skew-normal can be written as:
##     e_i = mu_i + sigma_i * Z_i,
##       where Z_i ~ CSN(0, 1, alpha_i, 0, 1) (unit skew-normal).
##   Collecting all shocks: e ~ CSN(mu_e, Sigma_e, Gamma_e, 0, I_{n_exo}) where
##     mu_e    = -E[e]  (joint CSN mean, .csn_shock_mean; = -sigma_i delta_i sqrt(2/pi)
##               per shock only when no skewed shock is correlated)
##     Sigma_e = diag(sigma_i^2)  (ignoring cross-correlations for skewness)
##     Gamma_e = diag(alpha_i / sigma_i)   (the shape scales by 1/sigma_i)
##   NOTE: cross-shock correlations enter Sigma_e normally; the CSN skewness
##   parameterisation for the state-space applies the affine transform:
##     eta = RR e   =>  eta ~ CSN(RR mu_e, RR Sigma_e RR', Gamma_e Sigma_e RR' (RR Sigma_e RR')^{-1}, ...)
##   By the CSN linear-map closure (Dominguez-Molina et al. 2003):
##     If x ~ CSN(mu, S, G, nu, D)  then  A x + b ~ CSN(A mu + b, A S A',
##       G S A' (A S A')^{-1}, nu, D + G S G' - G S A'(A S A')^{-1} A S G')
##   Applied here with A = RR, x = e (the shock), x ~ CSN(mu_e, Sigma_e, Gamma_e, 0_q, I_q):
##     mu_eta   = RR mu_e            (n_state)
##     S_eta    = RR Sigma_e RR'     (n_state x n_state)
##     Gamma_1  = Gamma_e Sigma_e RR' (S_eta)^{-1}   [the "new" skewness loading]
##     nu_1     = 0   [unchanged]
##     Delta_1  = I + Gamma_e Sigma_e Gamma_e' - Gamma_e Sigma_e RR' S_eta^{-1} RR Sigma_e Gamma_e'
##              = I + Gamma_e [Sigma_e - Sigma_e RR'(RR Sigma_e RR')^{-1} RR Sigma_e] Gamma_e'
##              = I + Gamma_e Sigma_e_{perp} Gamma_e'
##   where Sigma_e_{perp} = Sigma_e - Sigma_e RR' (RR Sigma_e RR')^{-1} RR Sigma_e
##   is the Schur complement (residual variance of e not explained by eta = RR e).
##   For a square invertible RR (n_state = n_exo case) Sigma_e_{perp} = 0 =>
##   Delta_1 = I.  More generally Delta_1 >= I.
##
##   SIMPLIFICATION used here (first-order PSKF v0):
##   We store the skewness loading in state space as:
##     Gamma_eta = Gamma_e Sigma_e RR' S_eta^{-1}   (q_eta x n_state)
##   and correspondingly Delta_eta = I + Gamma_e Sigma_e_{perp} Gamma_e'.
##   This is exact by CSN closure.
##   When n_state = n_exo and RR is square invertible:
##     Gamma_eta = diag(alpha_i/sigma_i) * Sigma_e * RR' * (RR Sigma_e RR')^{-1}
##   which simplifies to diag(alpha_i) * (RR)^{-1} when Sigma_e = diag(sigma_i^2).
##
## STATE-SPACE FORM (0.9.4): the lift is applied to the CONTEMPORANEOUS
## order-1 system built by .pskf_order1_statespace() -- state chi_t = the
## model variables [states; non-state observables] at t, RR = ghu[K, ], and a
## noise-free selection observation (DD = 0, so Sigma_eps = diag(me) only).
## The returned list additionally carries TT, ZZ (and the index bookkeeping)
## so every consumer filters the SAME system. See .pskf_order1_statespace for
## why the former (TT = ghx[state, ], ZZ = ghx[obs, ], eps = DD e) assembly
## was the wrong law whenever DD != 0.
#' @noRd
.get_csn_shock_params <- function(model, exo_names, obs_vars, dr, params,
                                   me_variance = 0) {
  ssm <- .pskf_order1_statespace(dr, obs_vars)
  RR  <- ssm$RR                              # n_chi x n_exo
  DD  <- ssm$DD                              # n_obs x n_exo (zero)

  Sigma_e <- .get_shock_cov(model, exo_names, params)  # n_exo x n_exo
  alpha   <- .get_shock_skewness(model, exo_names, params)  # n_exo named vector

  ## Correlated + skewed shocks: the FULL joint CSN.
  ##
  ## The joint shock law is  e ~ CSN(mu_e, Sigma_e, Gamma_e, 0, I_{n_exo})  with
  ##   Gamma_e = diag(alpha_i / sigma_i),   Delta_e = I,   seed cov = full Sigma_e
  ## (off-diagonals included).  This is the construction pinned numerically in
  ## Three properties were verified by brute force
  ## (1e6-draw CSN rejection sampler, 2-shock alpha = (+2, -2)):
  ##   (a) rho_12 = 0  collapses Sigma_e to diagonal, so the build is BIT-IDENTICAL
  ##       to the old per-shock-independent path (verified delta == 0).
  ##   (b) Delta_e = I is always PD, so the joint law is a valid CSN for any
  ##       admissible Sigma_e (no PD trap; the naive coupling
  ##       Delta_e[i,j] = alpha_i alpha_j rho_ij is NON-PD and is NOT used).
  ##   (c) the off-diagonal coupling is carried entirely by the seed Sigma_e:
  ##       the cross-shock CO-SKEWNESS  E[(e_i-Ee_i)^2 (e_j-Ee_j)]  flips sign
  ##       with the sign of rho_12  (measured +0.018 at rho=+0.5, -0.060 at
  ##       rho=-0.5, ~0 at rho=0).  The truncation latents inherit the seed
  ##       correlation through the CSN stochastic representation
  ##       [e; U] ~ N(., [[Sigma_e, Sigma_e Gamma_e']; [Gamma_e Sigma_e, .]]),
  ##       e | (U >= 0):  off-diagonal Sigma_e ==> correlated U_i, U_j ==>
  ##       nonzero cross-coskewness, exactly the coupling the guard refused.
  ##       Nothing in Delta_e changes; Sigma_e does ALL the work.
  ##
  ## The linear-map lift below already consumes the full Sigma_e (Sigma_eta,
  ## Gamma_eta, and the Schur term Sigma_e_perp all use Sigma_e, not its
  ## diagonal), so removing the guard is sufficient.  The SBC DGP draw
  ## (validate-sbc.R) is replaced by the MATCHING joint CSN rejection sampler so
  ## the DGP and likelihood share the same law (consistency).

  c(.csn_state_noise_lift(RR, DD, Sigma_e, alpha, me_variance),
    ssm[c("TT", "ZZ", "chi_idx", "state_pos", "obs_pos")])
}


## ---------------------------------------------------------------------------
## Order-1 PSKF state space in CONTEMPORANEOUS form
## ---------------------------------------------------------------------------
## dynhr's decision rule is in the LAGGED convention
##   s_t = TT s_{t-1} + RR e_t,      y_t = ZZ s_{t-1} + DD e_t
## (TT = ghx[state, ], ZZ = ghx[obs, ], RR/DD = ghu rows). .pskf_filter wants
##   x_t = TT x_{t-1} + eta_t,       y_t = ZZ x_t + eps_t,  eta INDEPENDENT of eps.
## Feeding it (ghx[state, ], ghx[obs, ], eta = RR e, eps = DD e) -- the
## pre-0.9.4 assembly -- re-reads the lagged system with x_t := s_{t-1}, which
## makes eta_t = RR e_{t-1} and eps_{t-1} = DD e_{t-1} the SAME shock while the
## filter treats them as independent: the wrong law whenever DD != 0 (every
## model in which an observable loads a contemporaneous shock). On the AR(1)
## x = 0.7 x(-1) + e, sd 0.3, observed x, it gave -25.114 against the exact
## (and Gaussian-KF) -11.379 at ZERO skewness.
##
## Exact contemporaneous form used instead: carry the model variables that are
## either states or observables, chi_t = y_t[K] with K = c(state_idx,
## setdiff(obs_idx, state_idx)) (states first, so the first n_s columns of
## ghx[K, ] are exactly the state columns):
##   chi_t = G chi_{t-1} + H e_t,   G = [ghx[K, ], 0],  H = ghu[K, ]
##   y_t   = S chi_t + u_t,         S = selection,      u_t = measurement error
## The observation is noise-free apart from the measurement error, so the CSN
## lift sees eta = H e and Sigma_eps = diag(me) only (DD = 0). The Gaussian
## stationary init P0 = Lyap(G, H Sigma_e H') is the stationary law of chi_0,
## whose state block is the stationary s_0 the Gaussian KF starts from, so at
## alpha = 0 the likelihood equals kalman_filter()'s exactly. (Same idea as
## the order-2 path below, which carries the raw innovation in the state.)
##
## Returns list(TT, ZZ, RR, DD, chi_idx, state_pos, obs_pos): chi_idx are the
## endo indices of chi, state_pos the positions of the states in chi
## (= 1..n_s), obs_pos the positions of the observables in chi.
#' @noRd
.pskf_order1_statespace <- function(dr, obs_vars) {
  state_idx <- dr$state_idx
  obs_idx   <- match(obs_vars, dr$endo_names)
  if (any(is.na(obs_idx)))
    .dynhr_abort("PSKF: observed variable(s) not in the model: ",
                 paste(obs_vars[is.na(obs_idx)], collapse = ", "),
                 class = "dynhr_error_pskf")
  n_s   <- length(state_idx)
  K     <- c(state_idx, setdiff(obs_idx, state_idx))
  n_chi <- length(K)
  n_obs <- length(obs_idx)

  TT <- matrix(0, n_chi, n_chi)
  TT[, seq_len(n_s)] <- dr$ghx[K, , drop = FALSE]
  RR <- dr$ghu[K, , drop = FALSE]
  obs_pos <- match(obs_idx, K)
  ZZ <- matrix(0, n_obs, n_chi)
  ZZ[cbind(seq_len(n_obs), obs_pos)] <- 1
  DD <- matrix(0, n_obs, ncol(RR))

  list(TT = TT, ZZ = ZZ, RR = RR, DD = DD, chi_idx = K,
       state_pos = seq_len(n_s), obs_pos = obs_pos)
}


#' Scale-free Moore-Penrose pseudoinverse of a symmetric PSD matrix
#'
#' Used for \eqn{(RR \Sigma_e RR')^{+}} in \code{.csn_state_noise_lift}.
#' The rank decision is made on the CORRELATION form
#' \eqn{D^{-1/2} S D^{-1/2}} (\eqn{D = diag(S)}) rather than on \code{S}
#' itself, so a state whose variance is genuinely small but nonzero (a shock
#' with \code{stderr = 1e-5} alongside one with \code{stderr = 1}) is not
#' truncated as a null direction the way \code{MASS::ginv} would. Rows/columns
#' whose variance is zero or at the rounding level of the matrix's own scale
#' (\code{diag(S) <= n * eps * max(diag(S))}) ARE null directions and are
#' dropped outright.
#'
#' @param S    Symmetric positive-semidefinite matrix.
#' @param rtol Relative eigenvalue cutoff on the correlation form
#'   (default \code{sqrt(.Machine$double.eps)}, matching \code{MASS::ginv}).
#' @return Matrix of the same dimension as \code{S}.
#' @noRd
.csn_sym_pinv <- function(S, rtol = sqrt(.Machine$double.eps)) {
  n <- nrow(S)
  P <- matrix(0, n, n)
  if (n == 0L) return(P)
  ## pmax(): a round-off-negative diagonal (-1e-20 on an exactly-known state
  ## direction) is a zero-variance coordinate, not a NaN-with-warning.
  ## ABSOLUTE FLOOR: a variance at the rounding level of the matrix's
  ## own scale, n * eps * max diag(S) (>= n * eps * lambda_max / n), is a
  ## known direction too. Without it a round-off residual on an exactly-known
  ## state (P - K F K' left 1e-312, or -- positive -- 1e-17 against O(1)
  ## variances) passed `d > 0`, its correlation-form row was ~e_i (the
  ## round-off covariances divided by a tiny d), so it was kept as a genuine
  ## direction with pseudoinverse entry 1 / 1e-312 = Inf. The rank decision
  ## among the kept coordinates stays scale-free (correlation form).
  dS <- diag(S)
  dS[!is.finite(dS)] <- 0
  floor_abs <- n * .Machine$double.eps * max(dS, 0)
  d <- sqrt(pmax(dS, 0))
  keep <- dS > floor_abs
  if (!any(keep)) return(P)
  ds <- d[keep]
  C  <- S[keep, keep, drop = FALSE] / outer(ds, ds)   # correlation form
  C  <- 0.5 * (C + t(C))
  e  <- eigen(C, symmetric = TRUE)
  ev_max <- max(e$values)
  if (!is.finite(ev_max) || ev_max <= 0) return(P)
  pos <- e$values > rtol * ev_max
  if (!any(pos)) return(P)
  V <- e$vectors[, pos, drop = FALSE]
  P[keep, keep] <- (V %*% ((1 / e$values[pos]) * t(V))) / outer(ds, ds)
  0.5 * (P + t(P))
}

## Inverse of a symmetric positive-semidefinite S: a pivoted-Cholesky inverse
## when S is unambiguously nonsingular, .csn_sym_pinv(S) otherwise. Where the
## Cholesky route is taken the result equals the pseudoinverse (the inverse IS
## the pseudoinverse of a nonsingular matrix) to round-off.
##
## The route is decided on PIVOTS of the correlation form, never on chol()
## merely succeeding (a garbage pivot from a matrix singular to round-off can
## survive chol() and give a wildly wrong inverse). The Cholesky route needs
##   1. every diagonal above .csn_sym_pinv's rounding-level floor (so no row
##      would be dropped as a known direction),
##   2. a full-rank PIVOTED factorisation of the correlation form C, and
##   3. n * ||C^{-1}||_F < max_cond. Since lambda_max(C) <= tr(C) = n and
##      1 / lambda_min(C) <= ||C^{-1}||_F, this bounds the condition number of
##      C from above; the bound sits far below 1 / rtol, so every matrix that
##      takes this route is also one whose eigenvalues all clear
##      .csn_sym_pinv's rank cutoff.
## A matrix failing any test goes to .csn_sym_pinv unchanged, so singular and
## near-singular inputs (a singular initial covariance, exactly known state
## directions) keep the eigen rule.
#' @noRd
.csn_sym_inv <- function(S, rtol = sqrt(.Machine$double.eps), fast = TRUE,
                         max_cond = 1e5) {
  n <- nrow(S)
  if (!fast || n == 0L) return(.csn_sym_pinv(S, rtol))
  dS <- diag(S)
  if (!all(is.finite(dS)) || min(dS) <= n * .Machine$double.eps * max(dS))
    return(.csn_sym_pinv(S, rtol))
  d <- sqrt(dS)
  C <- S / outer(d, d)
  C <- 0.5 * (C + t(C))
  if (!all(is.finite(C))) return(.csn_sym_pinv(S, rtol))
  ## The LAPACK pivoted factorisation warns (and reports rank < n) on a
  ## rank-deficient input; that case is exactly the fallback.
  R <- suppressWarnings(chol(C, pivot = TRUE))
  if (!identical(attr(R, "rank"), n)) return(.csn_sym_pinv(S, rtol))
  piv <- attr(R, "pivot")
  Ci  <- matrix(0, n, n)
  Ci[piv, piv] <- chol2inv(R)
  if (!all(is.finite(Ci)) || n * sqrt(sum(Ci^2)) >= min(max_cond, 0.25 / rtol))
    return(.csn_sym_pinv(S, rtol))
  P <- Ci / outer(d, d)
  0.5 * (P + t(P))
}


## log P(U >= 0) for U ~ N(0, V): the zero-orthant probability. By symmetry it
## equals log Phi_q(0; 0, V). Closed forms (Sheppard 1899; Plackett 1954) for
## q <= 3 -- exact and deterministic -- and the accurate logcdf_ME_r() path
## (C++ lattice evaluator) for 4 <= q <= 5 (Mendell-Elston beyond, see
## logcdf_ME_r).
#' @noRd
.csn_log_orthant0 <- function(V) {
  q <- nrow(V)
  if (q == 0L) return(0)
  if (q == 1L) return(log(0.5))
  if (q <= 3L) {
    d <- sqrt(diag(V))
    R <- V / outer(d, d)
    if (q == 2L) return(log(0.25 + asin(R[1L, 2L]) / (2 * pi)))
    return(log(0.125 + (asin(R[1L, 2L]) + asin(R[1L, 3L]) +
                          asin(R[2L, 3L])) / (4 * pi)))
  }
  logcdf_ME_r(rep(0, q), V)
}


#' Mean of the package's joint skew-normal shock law
#'
#' The skewed shocks are \eqn{e \sim CSN_{n,q}(0, \Sigma_e, \Gamma_e, 0, I_q)}
#' (closed skew-normal, Gonzalez-Farias, Dominguez-Molina and Gupta 2004), with
#' one truncation latent per skewed shock: \eqn{\Gamma_e} has the rows
#' \eqn{(\alpha_j/\sigma_j) e_j'} and the seed covariance is the FULL
#' \eqn{\Sigma_e}. This is NOT Azzalini's multivariate \eqn{SN(\Omega, \alpha)}
#' (which has a single latent, \eqn{q = 1}); the two coincide only for one
#' skewed shock. The stochastic representation is \eqn{e = x | U \ge 0} with
#' \eqn{[x; U] \sim N(0, [[\Sigma_e, \Sigma_e\Gamma_e'], [\Gamma_e\Sigma_e,
#' V]])}, \eqn{V = I + \Gamma_e \Sigma_e \Gamma_e'}, and its mean is
#' \deqn{E[e] = \Sigma_e \Gamma_e' g, \quad g_j = \phi(0; V_{jj})
#'   \Phi_{q-1}(0; V_{-j|j}) / \Phi_q(0; V),}
#' the gradient of the log normaliser with respect to the truncation point.
#'
#' When no skewed shock is correlated with any other shock, \eqn{V} is
#' diagonal, \eqn{g_j = \sqrt{2/\pi}/\sqrt{1+\alpha_j^2}} and the formula
#' reduces to the per-shock \eqn{\sigma_j \delta_j \sqrt{2/\pi}}; that case
#' returns the per-shock expression itself (bit-identical to the pre-0.9.4
#' code). Under correlation the per-shock formula is wrong twice: \eqn{V}'s
#' off-diagonals \eqn{\alpha_i \alpha_j \rho_{ij}} change \eqn{g}, and the
#' off-diagonal \eqn{\Sigma_e} gives a skewed shock's correlates -- even
#' Gaussian ones -- a nonzero mean.
#'
#' @param Sigma_e n_exo x n_exo seed covariance.
#' @param alpha   length-n_exo shape vector (0 = Gaussian shock).
#' @return length-n_exo numeric vector \eqn{E[e]}.
#' @noRd
.csn_shock_mean <- function(Sigma_e, alpha) {
  alpha <- as.numeric(alpha)
  n     <- length(alpha)
  Sigma_e <- as.matrix(Sigma_e)
  sigma <- sqrt(diag(Sigma_e))
  delta <- alpha / sqrt(1 + alpha^2)
  diag_mean <- sigma * delta * sqrt(2 / pi)
  ## A zero-stderr shock cannot be skewed (alpha/sigma undefined): no latent.
  sk <- which(alpha != 0 & is.finite(alpha) & sigma > 0)
  if (length(sk) == 0L) return(diag_mean)
  cross <- Sigma_e[, sk, drop = FALSE]
  cross[cbind(sk, seq_along(sk))] <- 0
  if (all(cross == 0)) return(diag_mean)

  q   <- length(sk)
  Gam <- matrix(0, q, n)
  Gam[cbind(seq_len(q), sk)] <- alpha[sk] / sigma[sk]
  V <- diag(q) + Gam %*% Sigma_e %*% t(Gam)
  V <- 0.5 * (V + t(V))
  log_Z <- .csn_log_orthant0(V)
  g <- numeric(q)
  for (j in seq_len(q)) {
    l_phi  <- dnorm(0, 0, sqrt(V[j, j]), log = TRUE)
    l_cond <- 0
    if (q > 1L) {
      Sc <- V[-j, -j, drop = FALSE] -
        V[-j, j, drop = FALSE] %*% V[j, -j, drop = FALSE] / V[j, j]
      l_cond <- .csn_log_orthant0(0.5 * (Sc + t(Sc)))
    }
    g[j] <- exp(l_phi + l_cond - log_Z)
  }
  if (!all(is.finite(g)))
    .dynhr_abort(".csn_shock_mean(): non-finite skew-normal mean; the shock ",
                 "covariance is not usable.", class = "dynhr_error_pskf")
  as.numeric(Sigma_e %*% t(Gam) %*% g)
}


#' Lift a skewed shock law through a linear state-space loading (generic core)
#'
#' The generic CSN state-noise construction shared by the DSGE path
#' (\code{.get_csn_shock_params}, which extracts \code{RR}/\code{DD} from a
#' perturbation solution) and non-DSGE linear state spaces (e.g. the Reiter
#' HANK adapter \code{hank_reiter_pskf_loglik}): given state loading
#' \code{eta = RR e}, observation feedthrough \code{DD e}, seed covariance
#' \code{Sigma_e} and per-shock skewness \code{alpha}, build the CSN
#' parameters of \code{eta} consumed by \code{.pskf_filter}. All derivation
#' notes live at the (single) call site above.
#' @noRd
.csn_state_noise_lift <- function(RR, DD, Sigma_e, alpha, me_variance = 0) {
  n_exo <- ncol(RR)
  n_obs <- nrow(DD)

  ## State noise covariance (header convention 1: ghu excludes Sigma_e)
  Sigma_eta <- RR %*% Sigma_e %*% t(RR)   # n_state x n_state

  ## Obs noise covariance (header convention 4: assembled as ONE matrix). me_variance is
  ## a scalar (H = me I) or one variance per observable (H = diag(me)).
  H_me <- if (length(me_variance) == 1L) me_variance * diag(n_obs)
          else diag(as.numeric(me_variance), nrow = n_obs)
  Sigma_eps <- DD %*% Sigma_e %*% t(DD) + H_me

  ## Mean correction (header convention 6): mu_eta is chosen s.t. E[eta] = 0
  ## (preserves the model steady state). E[e] is the mean of the JOINT CSN
  ## law, which only equals the per-shock sigma_i delta_i sqrt(2/pi) when no
  ## skewed shock is correlated with any other shock -- see .csn_shock_mean().
  sigma_e   <- sqrt(diag(Sigma_e))           # n_exo
  mean_e    <- .csn_shock_mean(Sigma_e, alpha)   # E[e], n_exo
  mu_eta    <- -as.numeric(RR %*% mean_e)    # n_state (sign: E[eta] = RR*E[e] -> correct to 0)
  mu_eps    <- rep(0, n_obs)                 # no mean shift on obs

  ## CSN skewness loading in state space (see DERIVATION above)
  ## q_eta = n_exo (one skewness dimension per shock)
  Gamma_e   <- diag(alpha / sigma_e, nrow = n_exo)  # n_exo x n_exo

  ## ---- Sigma_eta^{-1}: ONE rank-revealing pseudoinverse for BOTH terms -----
  ##
  ## Gamma_eta    = Gamma_e Sigma_e RR' S_eta^+                (n_exo x n_state)
  ## Sigma_e_perp = Sigma_e - Sigma_e RR' S_eta^+ RR Sigma_e   (Schur complement)
  ## Delta_eta    = I + Gamma_e Sigma_e_perp Gamma_e'
  ##
  ## Both used to be written as tryCatch(solve(Sigma_eta), <fallback>), with
  ## DIFFERENT fallbacks: ginv for Gamma_eta, but a bare identity for
  ## Delta_eta. That construction was wrong twice over:
  ##
  ##  (1) base::solve() only errors when rcond(Sigma_eta) < .Machine$double.eps
  ##      (~2.2e-16). A NEARLY singular Sigma_eta -- rcond ~1e-14, which is the
  ##      routine case when n_state > n_exo and the unshocked states are only
  ##      weakly excited -- does NOT throw; it returns a wildly amplified
  ##      "inverse". So the fallbacks fired almost never, and precisely in the
  ##      regime they were written for the code returned garbage instead.
  ##  (2) the identity fallback for Delta_eta is only correct when the Schur
  ##      complement vanishes, i.e. when rank(RR) = n_exo. That precondition
  ##      was documented but never checked.
  ##
  ## Using the SAME Moore-Penrose pseudoinverse in both places removes the
  ## inconsistency and needs no fallback at all: when rank(RR) = n_exo the
  ## Schur complement comes out numerically zero and Delta_eta = I falls out
  ## automatically; when it does not, the (correct, nonzero) Schur complement
  ## is carried. Rank is decided on the CORRELATION form of Sigma_eta so that
  ## states with a genuinely small-but-nonzero variance are not mistaken for
  ## null directions (a plain MASS::ginv, which thresholds on the largest
  ## absolute singular value, does exactly that).
  S_eta_pinv <- .csn_sym_pinv(Sigma_eta)

  Sigma_e_RRt <- Sigma_e %*% t(RR)                    # n_exo x n_state
  Gamma_eta   <- Gamma_e %*% Sigma_e_RRt %*% S_eta_pinv

  Sigma_e_perp <- Sigma_e - Sigma_e_RRt %*% S_eta_pinv %*% t(Sigma_e_RRt)
  Sigma_e_perp <- 0.5 * (Sigma_e_perp + t(Sigma_e_perp))
  Delta_eta    <- diag(n_exo) + Gamma_e %*% Sigma_e_perp %*% t(Gamma_e)
  Delta_eta    <- 0.5 * (Delta_eta + t(Delta_eta))

  ## Explicit precondition check (previously only a comment). Sigma_e_perp is
  ## a Schur complement of the PSD matrix [[Sigma_e, Sigma_e RR'],
  ## [RR Sigma_e, Sigma_eta]], hence PSD, so Delta_eta >= I always. A
  ## violation means the pseudoinverse rank decision was wrong (Sigma_eta is
  ## not the covariance of RR e, or Sigma_e is not PSD) and the CSN lift is
  ## not trustworthy -- reject the draw rather than filter with it.
  if (n_exo > 0L) {
    if (!all(is.finite(Delta_eta)) || !all(is.finite(Gamma_eta)))
      stop(".csn_state_noise_lift(): non-finite CSN skewness parameters; ",
           "Sigma_eta is not usable.")
    ev_min <- min(eigen(Delta_eta, symmetric = TRUE, only.values = TRUE)$values)
    tol_D  <- 1e-8 * max(1, max(abs(Delta_eta)))
    if (!is.finite(ev_min) || ev_min < 1 - tol_D)
      stop(sprintf(paste0(".csn_state_noise_lift(): Delta_eta is not >= I ",
                          "(min eigenvalue %.6g). The Schur complement ",
                          "Sigma_e - Sigma_e RR' (RR Sigma_e RR')^+ RR Sigma_e ",
                          "must be PSD; it is not, so the rank of Sigma_eta ",
                          "was mis-determined."), ev_min))
  }

  ## nu_eta = 0 (standard CSN for skew-normal shocks)
  nu_eta <- rep(0, n_exo)

  list(
    Sigma_eta = Sigma_eta,
    Sigma_eps = Sigma_eps,
    Gamma_eta = Gamma_eta,   # n_exo x n_state
    nu_eta    = nu_eta,
    Delta_eta = Delta_eta,   # n_exo x n_exo
    mu_eta    = mu_eta,
    mu_eps    = mu_eps,
    alpha     = alpha
  )
}


## ---------------------------------------------------------------------------
## Degenerate innovation covariances
## ---------------------------------------------------------------------------
## .pskf_independent_obs(Omega): order-preserving, scale-free selection of the
## linearly independent observables of a PSD innovation covariance Omega. It
## is an in-order Cholesky factorisation of the CORRELATION form
## D^{-1/2} Omega D^{-1/2} that skips row j when its pivot -- the conditional
## variance of observable j given the rows already kept, relative to its
## marginal variance -- is <= rtol (default sqrt(eps), the same scale-free
## cutoff as .csn_sym_pinv). Zero-variance rows are skipped outright. The
## pivots are TESTED, never trusted to chol(): chol() succeeds on matrices
## singular to round-off with a garbage pivot (the memory-noted "garbage pivot
## gave loglik too HIGH" failure).
## Returns list(keep = kept row indices (increasing), L = lower-triangular
## factor with Omega[keep, keep] = L L').
#' @noRd
.pskf_independent_obs <- function(Omega, rtol = sqrt(.Machine$double.eps)) {
  m    <- nrow(Omega)
  d    <- sqrt(pmax(diag(Omega), 0))
  keep <- integer(0)
  Lc   <- matrix(0, 0L, 0L)   # lower Cholesky factor of the kept correlation block
  for (j in seq_len(m)) {
    if (!(d[j] > 0)) next
    if (length(keep) == 0L) {
      w   <- numeric(0)
      piv <- 1
    } else {
      c_j <- Omega[keep, j] / (d[keep] * d[j])
      w   <- forwardsolve(Lc, c_j)
      piv <- 1 - sum(w^2)
    }
    if (piv > rtol) {
      k    <- length(keep)
      Lc   <- rbind(cbind(Lc, matrix(0, k, 1L)), c(w, sqrt(piv)))
      keep <- c(keep, j)
    }
  }
  list(keep = keep, L = d[keep] * Lc)
}

## .pskf_dependent_obs_consistent(): do the dependent observables' innovations
## equal the combination of the independent ones that Omega implies? For a
## dependent row r the residual e_r = v_r - Omega[r, K] Omega[K, K]^{-1} v_K
## has (to the rank rule's resolution) zero variance, so the data are
## consistent iff e_r = 0. Tolerance: round-off in v (rtol times the size of
## the terms y, ZZ mu, mu_eps it was formed from, v_scale) plus the rank
## rule's resolution sqrt(rtol * Omega_rr) -- the largest conditional sd
## .pskf_independent_obs treats as zero.
#' @noRd
.pskf_dependent_obs_consistent <- function(Omega, v, v_scale, ind,
                                           rtol = sqrt(.Machine$double.eps)) {
  keep <- ind$keep
  dep  <- setdiff(seq_along(v), keep)
  if (length(dep) == 0L) return(TRUE)
  if (length(keep) == 0L) {
    r       <- v[dep]
    r_scale <- v_scale[dep]
  } else {
    A <- t(backsolve(t(ind$L),
                     forwardsolve(ind$L, Omega[keep, dep, drop = FALSE])))
    r       <- v[dep] - as.numeric(A %*% v[keep])
    r_scale <- v_scale[dep] + as.numeric(abs(A) %*% v_scale[keep])
  }
  tol <- rtol * r_scale + sqrt(rtol * pmax(diag(Omega)[dep], 0))
  all(abs(r) <= tol)
}


## ---------------------------------------------------------------------------
## PSKF filter recursion
## ---------------------------------------------------------------------------
## Implements the CSN Kalman filter (PSKF) following skalman_filter.R /
## the published recursion.
##
## State space (first-order):
##   x_t = G x_{t-1} + eta_t,   eta_t ~ CSN(mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta)
##   y_t = F x_t + eps_t,        eps_t ~ N(mu_eps, Sigma_eps)
##
## At each step the filter maintains the CSN state distribution:
##   x_{t|t-1} ~ CSN(mu_pred, Sigma_pred, Gamma_pred, nu_pred, Delta_pred)
## with q-dimensional skewness parameter growing by q_eta per step, then
## pruned back by dim_red4_r.
##
## Returns scalar total log-likelihood when store_path=FALSE (default).
## When store_path=TRUE returns a list: list(ll, mu_pred_path, Sigma_pred_path,
## mu_filt_path, Sigma_filt_path, K_gauss_path, ZZ_path,
## Gamma_filt_path, nu_filt_path, Delta_filt_path,
## Gamma_pred_path, nu_pred_path, Delta_pred_path, keep_path, lambda_path) —
## arrays needed by the CSN backward pass in pskf_smoother(). keep_path[[t]]
## indexes the rows of the period-t pre-prune skew stack [U_{t-1}; u_t] that
## survived dim_red4_r (so row i of Gamma_pred_path[[t]] is latent
## keep_path[[t]][i] of that stack); lambda_path[[t]] is dim_red4_r's
## latent-space cut compensation on that stack (mu_shift = Sigma Gamma' lambda).
## Linearly dependent observables (singular Omega, e.g. duplicated noise-free
## rows) are handled exactly: see .pskf_independent_obs.
## offset_miwa_qmax (default 5): deterministic-evaluator range of the Phi
## evaluations in the pruning compensation (dim_red4_r); pskf_smoother raises
## it to 7. The likelihood's CDF terms and the compensation use the accurate
## evaluators of logcdf_ME_r (log-scale bivariate quadrature, C++ lattice rule
## for q >= 3); see .csn_mean_offset and logcdf_ME_r for their errors.
## max_q (default 5): hard rank-based cap on the retained skew dimension,
## chosen so every Phi_q call stays inside the deterministic-evaluator range
## (q <= 5).  The historical multi-shock bias (-7.3 nats at T=12,
## alpha=(+2,-2)) was ENTIRELY Mendell-Elston Phi_q evaluation error at
## q > 5, not discarded skew mass -- specifically the sign / double-shrink
## bug fixed in 0.9.4.6; capping q at 5 cut it to |gap| <= 0.14 nat on the
## worst measured fixture.  With the ME fix, max_q = Inf (cut_tol-threshold
## pruning only, ME for q > 5) is 0.01-0.03 nat from an exact grid on the
## uncapped multishock fixture, so the cap is now a cost / precision
## safeguard; it stays the default because the exact evaluators are more
## accurate than ME (~1e-6 vs 0.016 median).
##
## keep_override (default NULL): a FROZEN pruning selection -- a list of T
## integer vectors, element t the rows of the period-t pre-prune stack to
## keep (the keep_path of an earlier run). dim_red4_r then skips its
## selection and uses those rows, still computing the mean compensation for
## the cut. The stack at t has length(keep_override[[t - 1]]) + nrow(Gamma_eta)
## rows whatever theta is, so a path recorded at one theta is a valid
## selection at any other theta of the same model: the skew dimension is
## fixed by the path, not by the data or the parameters. Supplying the path
## the natural selection would have made changes nothing (same operations).
## store_keep = TRUE returns the loglik with attr(, "keep_path") (only when it
## is finite); the value itself is unchanged.
##
## Inside a frozen-selection scope (.pskf_fd_gradient, .pskf_freeze_wrap)
## a call without keep_override and without store_path takes the selection
## recorded by the scope's first evaluation, so every point of a
## finite-difference stencil differentiates the same smooth branch.
#' @noRd
.pskf_filter <- function(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                          Delta_eta, mu_eps, Sigma_eps,
                          cut_tol = 0.01,
                          max_q = 5L,
                          offset_miwa_qmax = 5L,
                          store_path = FALSE,
                          fast_solve = TRUE,
                          use_cpp = TRUE,
                          keep_override = NULL,
                          store_keep = FALSE) {
  q_eta <- nrow(as.matrix(Gamma_eta))
  if (!is.null(keep_override) &&
      !.pskf_keep_override_ok(keep_override, ncol(Y), q_eta))
    .dynhr_abort(".pskf_filter(): `keep_override` must be a list of ",
                 ncol(Y), " increasing integer vectors, element t indexing ",
                 "the rows of the period-t pre-prune skew stack ",
                 "(length(keep_override[[t - 1]]) + ", q_eta, " rows).",
                 class = "dynhr_error_bad_argument")
  st <- .pskf_freeze_state
  slot <- 0L
  if (is.null(keep_override) && !store_path && !store_keep && st$active) {
    st$i <- st$i + 1L
    slot <- st$i
    rec  <- if (slot <= length(st$paths)) st$paths[[slot]] else NULL
    if (is.null(rec)) {
      store_keep <- TRUE
    } else if (.pskf_keep_override_ok(rec, ncol(Y), q_eta)) {
      keep_override <- rec
    } else {
      slot <- 0L   # recorded for a different filter call: run unfrozen
    }
  }
  ll <- .pskf_filter_run(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                         Delta_eta, mu_eps, Sigma_eps, cut_tol, max_q,
                         offset_miwa_qmax, store_path, fast_solve, use_cpp,
                         keep_override, store_keep)
  if (slot > 0L && is.null(keep_override)) {
    kp <- attr(ll, "keep_path")
    if (!is.null(kp)) st$paths[[slot]] <- kp
    attr(ll, "keep_path") <- NULL
  }
  ll
}

## A keep path that fits the recursion: T increasing integer vectors, each
## inside the stack its predecessor leaves (length(prev) + q_eta rows).
#' @noRd
.pskf_keep_override_ok <- function(keep, n_T, q_eta) {
  if (!is.list(keep) || length(keep) != n_T) return(FALSE)
  q_prev <- 0L
  for (t in seq_len(n_T)) {
    k <- keep[[t]]
    q_new <- q_prev + q_eta
    if (!is.numeric(k) || anyNA(k) || any(k != round(k)) ||
        any(k < 1) || any(k > q_new) || is.unsorted(k, strictly = TRUE))
      return(FALSE)
    q_prev <- length(k)
  }
  TRUE
}

## ---------------------------------------------------------------------------
## Frozen-selection finite differences
## ---------------------------------------------------------------------------
## The pruning selection (dim_red4_r) is piecewise constant in theta, so the
## PSKF loglik is piecewise smooth with jumps of up to ~2e-3 nat where the
## kept set changes (measured on a 2-shock fixture). A finite-difference
## stencil straddling such a switch differentiates the jump: a central
## difference returns jump / (2h), 1-10 at h = 1e-4 and ~100 at the samplers'
## relative step 1e-5 -- as large as, or far larger than, the gradient itself.
## Inside a frozen-selection scope every PSKF filter call replays the
## selection recorded by the scope's FIRST evaluation, so all stencil points
## lie on one smooth branch and the difference is the derivative of that
## branch (the loglik itself, where no switch lies between the points).
##
## State: `active`; `paths`, the recorded keep paths in filter-call order
## (one per .pskf_filter call of an evaluation, so a posterior that filters
## twice is frozen per call); `i`, the call counter within an evaluation.
.pskf_freeze_state <- new.env(parent = emptyenv())
.pskf_freeze_state$active <- FALSE
.pskf_freeze_state$paths  <- list()
.pskf_freeze_state$i      <- 0L

## Open a scope; returns the previous state for .pskf_freeze_close (scopes
## nest: an inner gradient records and replays its own selection). `paths`
## seeds the scope with selections recorded elsewhere (.pskf_record_centre,
## e.g. in the main process for a daemon's stencil): every evaluation then
## replays them from the start instead of recording its own.
#' @noRd
.pskf_freeze_open <- function(paths = list()) {
  st <- .pskf_freeze_state
  prev <- list(active = st$active, paths = st$paths, i = st$i)
  st$active <- TRUE
  st$paths  <- paths
  st$i      <- 0L
  prev
}

## The selections a PSKF evaluation of fn at theta makes, one keep path per
## .pskf_filter call in call order (an empty list when fn runs no PSKF filter
## or its loglik is not finite), for seeding scopes elsewhere. The caller's
## RNG stream is left as it was, so a stochastic fn (a particle filter)
## evaluated here does not move anything drawn afterwards.
#' @noRd
.pskf_record_centre <- function(fn, theta) {
  ge  <- globalenv()
  had <- exists(".Random.seed", envir = ge, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = ge, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = ge)
    } else if (exists(".Random.seed", envir = ge, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = ge)
    }
  }, add = TRUE)
  prev <- .pskf_freeze_open()
  on.exit(.pskf_freeze_close(prev), add = TRUE, after = FALSE)
  .pskf_freeze_wrap(fn)(theta)
  .pskf_freeze_state$paths
}

#' @noRd
.pskf_freeze_close <- function(prev) {
  st <- .pskf_freeze_state
  st$active <- prev$active
  st$paths  <- prev$paths
  st$i      <- prev$i
  invisible(NULL)
}

## f wrapped so that each call is one evaluation of the scope: the filter-call
## counter restarts, so the k-th filter call of every evaluation replays the
## k-th recorded path. The first evaluation in the scope records.
#' @noRd
.pskf_freeze_wrap <- function(f) {
  force(f)
  function(x, ...) {
    .pskf_freeze_state$i <- 0L
    f(x, ...)
  }
}

## TRUE once the scope has recorded a selection (a PSKF filter ran to a
## finite loglik in it).
#' @noRd
.pskf_freeze_recorded <- function() {
  length(.pskf_freeze_state$paths) > 0L
}

## For a stencil whose first evaluation is not the centre (central
## differences): when that evaluation recorded a PSKF selection, discard it
## and record at the centre instead (one evaluation of fz, the wrapped
## function, at theta). Returns TRUE when the centre selected differently --
## the caller then re-evaluates its first stencil point on the centre's
## branch (an unchanged selection means that point already is on it). A
## function that runs no PSKF filter costs nothing extra.
#' @noRd
.pskf_freeze_anchor <- function(fz, theta) {
  if (!.pskf_freeze_recorded()) return(FALSE)
  first <- .pskf_freeze_state$paths
  .pskf_freeze_state$paths <- list()
  fz(theta)
  !identical(.pskf_freeze_state$paths, first)
}

## Frozen-selection finite-difference gradient of a scalar function f (e.g. a
## PSKF log-posterior) at theta, with the selection recorded at theta itself.
## method "central" (2d + 1 evaluations) or "forward" (d + 1); h is the
## absolute step (scalar or length d). lower / upper (scalar or length d):
## a central step that would leave [lower, upper] becomes the one-sided
## difference on the inside. A non-finite f at any stencil point gives NA
## for that coordinate. For a function that runs no PSKF filter this is the
## plain finite difference.
#' @noRd
.pskf_fd_gradient <- function(f, theta, h = 1e-4,
                              method = c("central", "forward"),
                              lower = -Inf, upper = Inf) {
  method <- match.arg(method)
  d <- length(theta)
  h <- rep_len(h, d)
  lower <- rep_len(lower, d)
  upper <- rep_len(upper, d)
  prev <- .pskf_freeze_open()
  on.exit(.pskf_freeze_close(prev), add = TRUE)
  fz <- .pskf_freeze_wrap(f)
  f0 <- fz(theta)
  g <- numeric(d)
  for (k in seq_len(d)) {
    e <- replace(numeric(d), k, h[k])
    up_ok <- theta[k] + h[k] <= upper[k]
    dn_ok <- theta[k] - h[k] >= lower[k]
    if (method == "central" && up_ok && dn_ok) {
      g[k] <- (fz(theta + e) - fz(theta - e)) / (2 * h[k])
    } else if (up_ok) {
      g[k] <- (fz(theta + e) - f0) / h[k]
    } else {
      g[k] <- (f0 - fz(theta - e)) / h[k]
    }
  }
  g[!is.finite(g)] <- NA_real_
  names(g) <- names(theta)
  g
}

## Theta-space gradient closure for mode finding with a PSKF posterior: the
## frozen-selection central difference of log_post_fn (value or list with
## $logpost) with optim()'s default step (ndeps = 1e-3), one-sided inside
## the prior bounds (prior_spec$lower / $upper, less the 1e-8 margin the
## box-constrained optimisers keep). Given to the L-BFGS-B polish in place
## of optim()'s internal finite differences, which cannot be frozen.
#' @noRd
.pskf_mode_grad_fn <- function(log_post_fn, prior_spec, h = 1e-3) {
  force(log_post_fn)
  lo <- stats::setNames(prior_spec$lower, prior_spec$name)
  hi <- stats::setNames(prior_spec$upper, prior_spec$name)
  lo[is.na(lo)] <- -Inf
  hi[is.na(hi)] <- Inf
  lo <- lo + 1e-8
  hi <- hi - 1e-8
  f <- function(th) {
    r <- log_post_fn(th)
    if (is.list(r)) r$logpost else r
  }
  function(theta) {
    l <- unname(lo[names(theta)]); u <- unname(hi[names(theta)])
    if (length(l) != length(theta)) l <- u <- NA_real_
    l[is.na(l)] <- -Inf
    u[is.na(u)] <- Inf
    .pskf_fd_gradient(f, theta, h = h, lower = l, upper = u)
  }
}

## The optimiser gradient a mode-finding run gets when it has no analytic
## one: for likelihood = "pskf" and a method with an L-BFGS-B polish stage
## (combined, nelder, cmaes_jade) the frozen-selection difference
## (.pskf_mode_grad_fn), otherwise NULL (the optimiser differences
## internally). Shared by run_mode_finding() and run_mode_mirai()'s chains.
#' @noRd
.pskf_mode_polish_grad <- function(likelihood, method, log_post_fn, prior_spec) {
  if (identical(likelihood, "pskf") &&
      length(method) == 1L && method %in% c("combined", "nelder", "cmaes_jade"))
    .pskf_mode_grad_fn(log_post_fn, prior_spec)
  else NULL
}

## The recursion behind .pskf_filter (arguments as there, all required).
#' @noRd
.pskf_filter_run <- function(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                              Delta_eta, mu_eps, Sigma_eps,
                              cut_tol, max_q, offset_miwa_qmax, store_path,
                              fast_solve, use_cpp, keep_override, store_keep) {
  ## Y: n_obs x T matrix
  n_obs   <- nrow(Y)
  n_T     <- ncol(Y)
  n_state <- nrow(TT)   # = ncol(TT)
  ## a scalar mu_eps means "the same mean for every observable"; expand it so
  ## the per-period row subsetting (missing / dependent rows) is well defined
  if (length(mu_eps) == 1L && n_obs > 1L) mu_eps <- rep(mu_eps, n_obs)

  ## ---- Initialisation -------------------------------------------------------
  ## Start from the Gaussian stationary distribution (Lyapunov initialization).
  ## The skewness is zero at t=0: Gamma = [], nu = [], Delta = [].
  ## Solve P0 from P0 = TT P0 TT' + Sigma_eta with the package's shared
  ## Lyapunov solver (R/stochsimul-monolith.R): doubling algorithm, RELATIVE
  ## convergence test, NaN on a non-stationary TT.
  ##
  ## The former local .lyapunov_solve_r() ran at most 500 PLAIN fixed-point
  ## steps against an ABSOLUTE 1e-10 tolerance and, on non-convergence,
  ## returned the last iterate SILENTLY. Fixed-point convergence is geometric
  ## at rate rho(TT)^2, so for rho = 0.999 the 500th iterate still carries
  ## ~0.999^1000 = 37% of the stationary variance missing -- a badly wrong P0
  ## that fed a finite, plausible-looking loglik. And for a non-stationary TT
  ## (no stationary P0 exists) it returned a divergent iterate, which the old
  ## `!is.finite` guard then replaced by Sigma_eta -- fabricating a finite
  ## likelihood for an infeasible draw. A NaN P0 now rejects the draw.
  P0 <- solve_lyapunov(TT, Sigma_eta)
  if (is.null(P0) || any(!is.finite(P0))) {
    if (store_path) return(list(ll = -Inf))
    return(-Inf)
  }

  ## The recursion in C++ (pskf_filter_cpp, src/pskf_filter.cpp): the same
  ## operations as the loop below, no R callbacks, with or without the path
  ## storage the smoother needs. It answers only
  ## for the accurate CDF evaluators; wherever the loop
  ## below would leave the code the C++ mirrors (a non-finite matrix, a CDF
  ## the C++ evaluators decline, a failed factorisation) it reports status 1
  ## and the R loop below runs from the start. use_cpp = FALSE is the
  ## reference R recursion.
  if (use_cpp &&
      is.numeric(Y) && length(mu_eps) == n_obs &&
      length(mu_eta) == n_state && is.matrix(Sigma_eps) &&
      identical(dim(Sigma_eps), c(n_obs, n_obs)) &&
      identical(dim(as.matrix(ZZ)), c(n_obs, n_state)) &&
      identical(dim(as.matrix(Sigma_eta)), c(n_state, n_state)) &&
      ncol(as.matrix(Gamma_eta)) == n_state &&
      length(nu_eta) == nrow(as.matrix(Gamma_eta)) &&
      identical(dim(as.matrix(Delta_eta)), rep(length(nu_eta), 2L))) {
    res <- pskf_filter_cpp(
      matrix(as.double(Y), n_obs, n_T), as.matrix(TT), as.matrix(ZZ),
      as.numeric(mu_eta), as.matrix(Sigma_eta), as.matrix(Gamma_eta),
      as.numeric(nu_eta), as.matrix(Delta_eta), as.numeric(mu_eps),
      as.matrix(Sigma_eps), P0, as.numeric(cut_tol), as.numeric(max_q),
      as.numeric(offset_miwa_qmax), isTRUE(fast_solve), isTRUE(store_path),
      keep_override, isTRUE(store_keep))
    if (res$status == 0L) {
      if (!store_path) {
        if (isTRUE(store_keep) && is.finite(res$ll))
          return(structure(res$ll, keep_path = res$keep_path))
        return(res$ll)
      }
      if (!is.finite(res$ll)) return(list(ll = res$ll))
      res$status <- NULL
      return(res)
    }
  }

  mu_filt    <- rep(0, n_state)
  Sigma_filt <- P0
  Gamma_filt <- matrix(0, nrow = 0L, ncol = n_state)
  nu_filt    <- numeric(0)
  Delta_filt <- matrix(0, nrow = 0L, ncol = 0L)

  ## Log-likelihood accumulator
  ll <- 0

  ## Loop invariants
  tTT   <- t(TT)
  I_n   <- diag(n_state)
  tGe   <- t(Gamma_eta)

  ## Pre-compute log(2*pi) constant
  ll_const <- -0.5 * n_obs * log(2 * pi)

  ## Per-period path storage (when store_path = TRUE)
  if (store_path) {
    mu_pred_path    <- vector("list", n_T)
    Sigma_pred_path <- vector("list", n_T)
    mu_filt_path    <- vector("list", n_T)
    Sigma_filt_path <- vector("list", n_T)
    K_gauss_path    <- vector("list", n_T)
    ZZ_path         <- vector("list", n_T)
    ## CSN skewness path (needed for CSN backward smoother)
    Gamma_filt_path <- vector("list", n_T)
    nu_filt_path    <- vector("list", n_T)
    Delta_filt_path <- vector("list", n_T)
    Gamma_pred_path <- vector("list", n_T)
    nu_pred_path    <- vector("list", n_T)
    Delta_pred_path <- vector("list", n_T)
    keep_path       <- vector("list", n_T)
    lambda_path     <- vector("list", n_T)
  } else if (store_keep) {
    keep_path       <- vector("list", n_T)
  }

  for (t in seq_len(n_T)) {
    ## ---- PREDICTION STEP ----------------------------------------------------
    ## Gaussian part: standard Kalman prediction
    ##   mu_{t|t-1}    = TT mu_{t-1|t-1} + mu_eta
    ##   Sigma_{t|t-1} = TT Sigma_{t-1|t-1} TT' + Sigma_eta
    mu_pred    <- as.numeric(TT %*% mu_filt) + mu_eta
    TT_Sfilt   <- TT %*% Sigma_filt
    Sigma_pred <- TT_Sfilt %*% tTT + Sigma_eta

    ## CSN part: by the affine-map + sum closure of the CSN distribution.
    ## Stacking the current filtered skewness (Gamma_filt, nu_filt, Delta_filt)
    ## with the incoming shock skewness (Gamma_eta, nu_eta, Delta_eta).
    ##
    ## DERIVATION of the prediction Gamma/nu/Delta (from skalman_filter.R):
    ## Let x_{t-1|t-1} ~ CSN(m, S, G, nu, D) and eta_t ~ CSN(m_e, S_e, G_e, nu_e, D_e)
    ## independent. Then x_t = TT x_{t-1} + eta_t gives (by independence + closure):
    ##   x_t ~ CSN(TT m + m_e, TT S TT' + S_e, Gamma_stacked, nu_stacked, Delta_stacked)
    ## where:
    ##   Gamma_stacked = [G_filt * S * TT' * S_pred^{-1}; G_eta * S_eta * S_pred^{-1}]
    ##   (each block has its own covariance projected onto S_pred = TT S TT' + S_e)
    ## and nu_stacked = [nu_filt; nu_eta] (vertical concatenation)
    ## and Delta_stacked is a 4-block matrix:
    ##   [D_filt + G_filt(S - S*TT'*S_pred^{-1}*TT*S)*G_filt',  -G_filt*S*TT'*S_pred^{-1}*S_eta*G_eta']
    ##   [-G_eta*S_eta*S_pred^{-1}*TT*S*G_filt',     D_eta + G_eta*(S_eta - S_eta*S_pred^{-1}*S_eta)*G_eta']
    ## (This is the covariance formula from the Schur complement of S_pred in the
    ## joint covariance of (x_{t-1}, eta_t) under the linear map x_t = TT x_{t-1} + eta_t.)

    S_pred     <- Sigma_pred   # n_state x n_state
    q_filt     <- nrow(Gamma_filt)
    q_eta      <- nrow(Gamma_eta)
    q_new      <- q_filt + q_eta

    if (q_new == 0L) {
      ## Pure Gaussian: no skewness
      Gamma_pred <- matrix(0, nrow = 0L, ncol = n_state)
      nu_pred    <- numeric(0)
      Delta_pred <- matrix(0, nrow = 0L, ncol = 0L)
    } else {
      ## Sigma_pred^+ for the CSN projections Cov(U, x) Var(x)^+. Sigma_pred
      ## is routinely SINGULAR: the contemporaneous order-1 state
      ## (.pskf_order1_statespace) carries observables that are exact linear
      ## functions of the states and the shocks, noise-free observation makes
      ## some state directions known exactly, and the order-2 state carries
      ## Kronecker duplicates. For a genuine joint Gaussian (x, U), Cov(U, x)
      ## lies in the range of Var(x), so the Moore-Penrose pseudoinverse gives
      ## the EXACT conditional (the null directions have Cov(U, x) = 0 too).
      ## The former solve() with a `+ 1e-8 I` fallback on error was not exact
      ## there: on an AR(2) observed without error (a lagged state known
      ## exactly) it was 4.5e-6 nats off the closed-form skew-normal likelihood,
      ## and a near-singular (but not error-raising) Sigma_pred went through
      ## solve() unregularised. Same scale-free rank rule as the lift
      ## (.csn_sym_pinv).
      S_pred_inv <- .csn_sym_inv(S_pred, fast = fast_solve)

      ## Shared projections (each computed once): with the symmetric
      ## Sigma_filt, Sigma_filt TT' = (TT Sigma_filt)'.
      if (q_filt > 0L) {
        Sfilt_TT_t <- t(TT_Sfilt)                       # Sigma_filt TT'
        M_filt     <- Sfilt_TT_t %*% S_pred_inv         # n_state x n_state
        G_filt_block <- Gamma_filt %*% M_filt           # q_filt x n_state
      }
      if (q_eta > 0L) {
        M_eta       <- Sigma_eta %*% S_pred_inv         # n_state x n_state
        G_eta_block <- Gamma_eta %*% M_eta              # q_eta  x n_state
      }

      ## Stacked Gamma (q_new x n_state)
      if (q_filt == 0L) {
        Gamma_pred <- G_eta_block
      } else if (q_eta == 0L) {
        Gamma_pred <- G_filt_block
      } else {
        Gamma_pred <- rbind(G_filt_block, G_eta_block)
      }

      ## Stacked nu
      nu_pred <- c(nu_filt, nu_eta)

      ## Delta 4-block (Schur complement of S_pred in the joint CSN skewness
      ## covariance):
      ##   D11 = D_filt + G_filt (S_filt - S_filt TT' S_pred^{-1} TT S_filt) G_filt'
      ##   D22 = D_eta  + G_eta  (S_eta  - S_eta  S_pred^{-1}       S_eta)   G_eta'
      ##   D12 = -G_filt S_filt TT' S_pred^{-1} S_eta G_eta'
      ##   D21 = D12'
      if (q_filt > 0L) {
        Schur_filt <- M_filt %*% TT_Sfilt
        D11 <- Delta_filt +
               Gamma_filt %*% (Sigma_filt - Schur_filt) %*% t(Gamma_filt)
      }
      if (q_eta > 0L) {
        Schur_eta <- Sigma_eta - M_eta %*% Sigma_eta
        D22 <- Delta_eta + Gamma_eta %*% Schur_eta %*% tGe
      }
      if (q_filt == 0L) {
        Delta_pred <- D22
      } else if (q_eta == 0L) {
        Delta_pred <- D11
      } else {
        D12 <- -(G_filt_block %*% Sigma_eta) %*% tGe
        Delta_pred <- rbind(
          cbind(D11, D12),
          cbind(t(D12), D22)
        )
      }
    }

    ## ---- PRUNE ---------------------------------------------------------------
    ## Keep q bounded Pruning after prediction before update.
    ## keep_t: which rows of the pre-prune stack [U_{t-1}; u_t] survive (the
    ## exact smoother tracks each retained skew latent back to its birth).
    keep_t   <- integer(0)
    lambda_t <- numeric(0)
    if (q_new > 0L) {
      pruned <- dim_red4_r(Gamma_pred, nu_pred, Delta_pred, Sigma_pred, cut_tol,
                           max_q = max_q,
                           offset_miwa_qmax = offset_miwa_qmax,
                           keep = keep_override[[t]])
      Gamma_pred <- pruned$Gamma
      nu_pred    <- pruned$nu
      Delta_pred <- pruned$Delta
      keep_t     <- pruned$keep
      lambda_t   <- pruned$lambda
      ## first-moment compensation for the cut skew mass (see dim_red4_r):
      ## without this, saturated pruning drifts the state mean systematically
      mu_pred <- mu_pred + pruned$mu_shift
    }
    q_pred <- nrow(Gamma_pred)
    if (store_keep) keep_path[[t]] <- keep_t

    ## ---- OBSERVED ROWS -------------------------------------------------------
    ## Partial observation: update only on available obs (mirrors Gaussian KF)
    y_t <- Y[, t]
    obs_mask <- !is.na(y_t)
    if (all(obs_mask)) {
      ZZ_t        <- ZZ
      Sigma_eps_t <- Sigma_eps
      mu_eps_t    <- mu_eps
    } else {
      ZZ_t        <- ZZ[obs_mask, , drop = FALSE]
      Sigma_eps_t <- Sigma_eps[obs_mask, obs_mask, drop = FALSE]
      mu_eps_t    <- mu_eps[obs_mask]
      y_t         <- y_t[obs_mask]
    }

    ## ---- DEGENERATE (LINEARLY DEPENDENT) OBSERVATIONS -------------------------
    ## Gaussian innovation:
    ##   v_t   = y_t - ZZ mu_{t|t-1} - mu_eps
    ##   Omega = ZZ Sigma_{t|t-1} ZZ' + Sigma_eps  (innovation covariance)
    ## Omega is SINGULAR when an observable is (to the scale-free rank rule of
    ## .pskf_independent_obs) an exact linear function of the others: e.g. two
    ## noise-free observations of the same state, or a noise-free observable
    ## of a zero-variance state. chol() is NOT a singularity test: on an Omega
    ## singular to round-off it can succeed with a garbage pivot (and the old
    ## solve(Omega) then failed with an unclassed Lapack error), and when it
    ## fails the draw used to be rejected (-Inf) even for data that satisfy the
    ## degeneracy. The exact treatment: a dependent observable is a
    ## deterministic function of the independent ones given the state, so
    ##   * if its innovation equals the implied combination of the independent
    ##     innovations (consistent data) it carries no further information and
    ##     is dropped -- the likelihood is the density of the first-in-order
    ##     linearly independent observables, i.e. EXACTLY the likelihood with
    ##     the dependent row removed;
    ##   * otherwise the data have zero density under the model: ll = -Inf.
    if (length(y_t) > 0L) {
      Zmu_t <- as.numeric(ZZ_t %*% mu_pred)
      v_t   <- as.numeric(y_t - Zmu_t - mu_eps_t)
      Omega <- ZZ_t %*% Sigma_pred %*% t(ZZ_t) + Sigma_eps_t
      if (!all(is.finite(Omega)) || !all(is.finite(v_t))) { ll <- -Inf; break }
      ind <- .pskf_independent_obs(Omega)
      if (length(ind$keep) < length(v_t)) {
        v_scale <- abs(y_t) + abs(Zmu_t) + abs(rep_len(mu_eps_t, length(v_t)))
        if (!.pskf_dependent_obs_consistent(Omega, v_t, v_scale, ind)) {
          ll <- -Inf; break
        }
        k           <- ind$keep
        ZZ_t        <- ZZ_t[k, , drop = FALSE]
        Sigma_eps_t <- as.matrix(Sigma_eps_t)[k, k, drop = FALSE]
        y_t         <- y_t[k]
        v_t         <- v_t[k]
        Omega       <- Omega[k, k, drop = FALSE]
      }
    }

    if (length(y_t) == 0L) {
      ## Nothing (independent) observed: prediction becomes filtered
      mu_filt    <- mu_pred
      Sigma_filt <- Sigma_pred
      Gamma_filt <- Gamma_pred
      nu_filt    <- nu_pred
      Delta_filt <- Delta_pred
      ## No likelihood contribution; store degenerate path entry (K=0, no obs)
      if (store_path) {
        mu_pred_path[[t]]    <- mu_pred
        Sigma_pred_path[[t]] <- Sigma_pred
        mu_filt_path[[t]]    <- mu_filt
        Sigma_filt_path[[t]] <- Sigma_filt
        K_gauss_path[[t]]    <- matrix(0, n_state, 0L)
        ZZ_path[[t]]         <- matrix(0, 0L, n_state)
        Gamma_filt_path[[t]] <- Gamma_filt
        nu_filt_path[[t]]    <- nu_filt
        Delta_filt_path[[t]] <- Delta_filt
        Gamma_pred_path[[t]] <- Gamma_pred
        nu_pred_path[[t]]    <- nu_pred
        Delta_pred_path[[t]] <- Delta_pred
        keep_path[[t]]       <- keep_t
        lambda_path[[t]]     <- lambda_t
      }
      next
    }
    n_obs_t <- length(y_t)

    ## ---- UPDATE STEP ---------------------------------------------------------
    ## Omega is now the kept block, whose correlation-form pivots all exceed
    ## sqrt(eps) (tested above), so chol() cannot meet a near-zero pivot and
    ## no jitter is ever added (the former chol(Omega + 1e-8 I) fabricated a
    ## finite loglik on infeasible draws). The rcond check of solve() is
    ## disabled (tol = 0) for the same reason: a badly SCALED Omega that is
    ## well conditioned in correlation form is not singular. (The gain keeps
    ## the explicit-inverse form: a noise-free observation of a state then
    ## gets the exact gain 1; a sqrt-factor route leaves a 1-ulp residual
    ## variance that the scale-free pseudoinverse of the next Sigma_pred
    ## treats as genuine -- it broke the AR(2) closed-form test.)
    Omega_chol    <- chol(Omega)
    log_det_Omega <- 2 * sum(log(diag(Omega_chol)))
    Omega_inv_v   <- backsolve(Omega_chol, forwardsolve(t(Omega_chol), v_t))

    ## Gaussian log-likelihood contribution: log N(v_t; 0, Omega)
    ll_gauss <- -0.5 * (n_obs_t * log(2 * pi) + log_det_Omega +
                        sum(v_t * Omega_inv_v))

    ## Kalman gain (Gaussian)
    K_gauss <- Sigma_pred %*% t(ZZ_t) %*% solve(Omega, tol = 0)

    ## CSN update:
    ## Gamma and Delta UNCHANGED; only nu shifts:
    ##   nu_{t|t} = nu_{t|t-1} - K_skewed * v_t
    ## where K_skewed = Gamma_pred %*% K_gauss  (q_pred x n_obs_t)
    if (q_pred > 0L) {
      K_skewed <- Gamma_pred %*% K_gauss
      nu_upd   <- nu_pred - as.numeric(K_skewed %*% v_t)
    } else {
      nu_upd   <- nu_pred
    }

    ## Standard Gaussian state update
    mu_upd    <- mu_pred    + as.numeric(K_gauss %*% v_t)
    I_KZ      <- I_n - K_gauss %*% ZZ_t
    ## Joseph form. Algebraically identical to the short form
    ## `I_KZ %*% Sigma_pred` at the exact Kalman gain, but it is a sum of two
    ## explicitly symmetric PSD terms, so it stays symmetric and PSD under
    ## round-off. The short form is neither: it is not symmetric even in exact
    ## arithmetic as written (only S_upd = S_pred - K Omega K' is), and the
    ## asymmetry was fed straight into D_top = Delta_pred + Gamma_pred
    ## Sigma_upd Gamma_pred' -- symmetrised there, but only after the CDF
    ## argument had already been built from an asymmetric covariance -- and
    ## into the next period's prediction. Explicit symmetrisation on top costs
    ## one n_state^2 add and removes the drift entirely.
    Sigma_upd <- I_KZ %*% Sigma_pred %*% t(I_KZ) +
                 K_gauss %*% Sigma_eps_t %*% t(K_gauss)
    Sigma_upd <- 0.5 * (Sigma_upd + t(Sigma_upd))

    ## ---- LOGLIK CSN CORRECTION -----------------------------------------------
    ## The full CSN loglik is:
    ##   log p(y_t | Y_{1:t-1}) = log N(y_t; ZZ mu_pred + mu_eps, Omega)
    ##                           + log Phi_q(nu_adj_top;  0, D_top)
    ##                           - log Phi_q(-nu_pred;    0, D_bot)
    ## where:
    ##   D_bot = Delta_pred + Gamma_pred Sigma_pred Gamma_pred'
    ##           (the normalisation constant of the CSN predictive distribution)
    ##   The "top" arguments come from the CSN after conditioning on y_t:
    ##   nu_adj_top = nu_{t|t} = nu_pred - K_skewed v_t   (already computed above)
    ##   D_top      = Delta_pred + Gamma_pred Sigma_upd Gamma_pred'
    ##                (skewness covariance with UPDATED state uncertainty)
    ##
    ## NOTE: when alpha = 0 for all shocks, Gamma_pred = 0 (all rows zero),
    ## q_pred = 0 after pruning, so both Phi_q terms are log(1) = 0 and the
    ## loglik reduces exactly to the Gaussian log N term.  The two-CDF
    ## correction vanishes ALGEBRAICALLY, ensuring the zero-alpha reduction
    ## is exact.

    if (q_pred > 0L) {
      ## Bottom: CSN normalisation constant of predictive distribution
      D_bot <- Delta_pred + Gamma_pred %*% Sigma_pred %*% t(Gamma_pred)
      ## Ensure symmetry (numerical noise)
      D_bot <- 0.5 * (D_bot + t(D_bot))

      ## Top: CSN normalisation constant of updated distribution
      D_top <- Delta_pred + Gamma_pred %*% Sigma_upd %*% t(Gamma_pred)
      D_top <- 0.5 * (D_top + t(D_top))

      ## log Phi_q(-nu_pred; 0, D_bot) -- denominator
      ll_cdf_bot <- logcdf_ME_r(-nu_pred, D_bot)

      ## log Phi_q(-nu_{t|t}; 0, D_top) -- numerator: the normalisation
      ## constant of the UPDATED (posterior) CSN is Phi_q(-nu_upd; 0, D_top),
      ## and Phi_q(K_skew v - nu_pred) == Phi_q(-nu_upd).  (Sign verified
      ## against an exact single-period quadrature oracle: with +nu_upd the
      ## skew direction inverts and the loglik is systematically wrong.)
      ll_cdf_top <- logcdf_ME_r(-nu_upd, D_top)

      ll_skew <- ll_cdf_top - ll_cdf_bot
    } else {
      ## Pure Gaussian: both CDF terms are log(1) = 0, correction = 0.
      ll_skew <- 0
    }

    ll <- ll + ll_gauss + ll_skew

    ## ---- Store filtered state for next period --------------------------------
    mu_filt    <- mu_upd
    Sigma_filt <- Sigma_upd
    Gamma_filt <- Gamma_pred   # unchanged (update rule)
    nu_filt    <- nu_upd
    Delta_filt <- Delta_pred   # unchanged

    ## ---- Save path for smoother ----------------------------------------------
    if (store_path) {
      mu_pred_path[[t]]    <- mu_pred
      Sigma_pred_path[[t]] <- Sigma_pred
      mu_filt_path[[t]]    <- mu_filt
      Sigma_filt_path[[t]] <- Sigma_filt
      K_gauss_path[[t]]    <- K_gauss
      ZZ_path[[t]]         <- ZZ_t
      ## CSN skewness path
      Gamma_filt_path[[t]] <- Gamma_filt     # = Gamma_pred (update rule: unchanged)
      nu_filt_path[[t]]    <- nu_filt         # = nu_upd
      Delta_filt_path[[t]] <- Delta_filt     # = Delta_pred (update rule: unchanged)
      Gamma_pred_path[[t]] <- Gamma_pred
      nu_pred_path[[t]]    <- nu_pred
      Delta_pred_path[[t]] <- Delta_pred
      keep_path[[t]]       <- keep_t
      lambda_path[[t]]     <- lambda_t
    }
  }

  if (store_path) {
    list(
      ll              = ll,
      mu_pred_path    = mu_pred_path,
      Sigma_pred_path = Sigma_pred_path,
      mu_filt_path    = mu_filt_path,
      Sigma_filt_path = Sigma_filt_path,
      K_gauss_path    = K_gauss_path,
      ZZ_path         = ZZ_path,
      Gamma_filt_path = Gamma_filt_path,
      nu_filt_path    = nu_filt_path,
      Delta_filt_path = Delta_filt_path,
      Gamma_pred_path = Gamma_pred_path,
      nu_pred_path    = nu_pred_path,
      Delta_pred_path = Delta_pred_path,
      keep_path       = keep_path,
      lambda_path     = lambda_path
    )
  } else if (store_keep && is.finite(ll)) {
    structure(ll, keep_path = keep_path)
  } else {
    ll
  }
}


## ---------------------------------------------------------------------------
## .pskf_orient_data: resolve the T x n_obs / n_obs x T orientation
## ---------------------------------------------------------------------------
#' Orient an observation matrix as n_obs x T
#'
#' Accepts either \code{T x n_obs} (columns = observables) or
#' \code{n_obs x T}. A non-square matrix is oriented by its shape (a matrix
#' whose column count equals \code{length(obs_vars)} is transposed). A square
#' matrix is ambiguous by shape, so it is oriented by its dimnames: column
#' names equal to \code{obs_vars} mean \code{T x n_obs}, row names equal to
#' \code{obs_vars} mean \code{n_obs x T}. A square matrix that cannot be
#' oriented from its names (unnamed, or named both ways) is an error rather
#' than a guess.
#'
#' @param data     Numeric matrix (or data frame).
#' @param obs_vars Character vector of observed variable names.
#' @param fn       Calling function name, for the error message.
#' @return Numeric matrix, n_obs x T.
#' @noRd
.pskf_orient_data <- function(data, obs_vars, fn) {
  n_obs <- length(obs_vars)
  if (nrow(data) == ncol(data) && ncol(data) == n_obs) {
    col_ok <- !is.null(colnames(data)) && identical(colnames(data), obs_vars)
    row_ok <- !is.null(rownames(data)) && identical(rownames(data), obs_vars)
    if (col_ok && !row_ok) return(t(data))
    if (row_ok && !col_ok) return(as.matrix(data))
    .dynhr_abort(sprintf(
      paste0("%s: `data` is square (%d x %d), so its orientation cannot be ",
             "inferred from its shape. Give it column names equal to ",
             "`obs_vars` (T x n_obs) or row names equal to `obs_vars` ",
             "(n_obs x T), but not both."), fn, n_obs, n_obs),
      call. = FALSE)
  }
  if (ncol(data) == n_obs) t(data) else data
}

## ---------------------------------------------------------------------------
## make_log_posterior_pskf: factory function (mirrors make_log_posterior_tpf)
## ---------------------------------------------------------------------------
#' Create a PSKF log-posterior evaluator
#'
#' Factory called once before MCMC. Returns a closure function(theta) ->
#' list(logpost, loglik, logprior) using the Pruned Skewed Kalman Filter.
#'
#' @param model       dynhr_mod (from parse_mod())
#' @param data        Observation matrix (T x n_obs), columns = obs_vars
#' @param prior_spec  Prior specification (from extract_prior_spec())
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled (from compile_model())
#' @param me_variance Measurement error variance (default 0)
#' @param system_priors Named list of system-prior closures, or NULL
#' @param cut_tol     Pruning tolerance (default 0.01)
#' @param max_q       Hard cap on the retained skew dimension (default 5,
#'   the deterministic Phi_q evaluator range; see .pskf_filter). Inf = uncapped
#'   (Mendell-Elston for q > 5; within ~0.03 nat of an exact grid on the
#'   measured multi-shock fixture, but less accurate than the capped default).
#' @param power       Power-posterior exponent applied to the LIKELIHOOD only
#'   (prior and system priors stay untempered). \code{NULL} (default) resolves
#'   the \code{power_posterior} package option ONCE, at factory time, so every
#'   draw evaluated by the returned closure uses the same tempering -- an
#'   option flipped mid-chain can no longer silently change the target
#'   distribution between draws.
#' @param pskf_cdf    "accurate" (the only value; see .pskf_cdf_settings); NULL
#'   (default) resolves the \code{pskf_cdf} package option once, at factory
#'   time.
#' @param ...         Ignored (for interface compatibility)
#' @return function(theta) -> list(logpost, loglik, logprior)
#' @noRd
make_log_posterior_pskf <- function(model, data, prior_spec, obs_vars,
                                     compiled, me_variance = 0,
                                     system_priors = NULL,
                                     cut_tol = 0.01,
                                     max_q = 5L,
                                     power = NULL,
                                     pskf_cdf = NULL,
                                     ...) {
  ## Resolve the power-posterior exponent ONCE, here, rather than on every
  ## evaluation: the closure's target must not change under the caller's feet.
  power <- .dynhr_opt("power_posterior", power, default = 1)
  cdf   <- .pskf_cdf_settings(pskf_cdf)
  ## Scalar or per-observable H = diag(me), validated at build time (the
  ## per-draw tryCatch would turn a bad value into -Inf at every theta).
  me_variance <- .kf_me_variance(me_variance, obs_vars,
                                 "make_log_posterior_pskf")
  ## Data orientation: either T x n_obs or n_obs x T -> Y is n_obs x T
  Y <- .pskf_orient_data(data, obs_vars, "make_log_posterior_pskf")

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  ## Adapter over the shared closure builder (R/posterior-closure.R). Specific
  ## to this branch: the CSN state-space assembly + .pskf_filter (loglik hook)
  ## and the system-prior convention -- the PSKF pair takes a BARE LIST of
  ## `function(dr, params)` closures, summed with a short-circuit on the first
  ## non-finite value, rather than the `system_prior_spec` objects the
  ## Kalman/pruned branches take, and (mode "extra") keeps the result out of
  ## $logprior. Both are pinned by test-posterior-closure-parity.R.
  .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      dr <- sol$dr
      obs_idx <- match(obs_vars, dr$endo_names)
      if (any(is.na(obs_idx))) return(NULL)

      ## Observation mean: d_obs (includes SS value and any DR offset)
      d_obs <- dr$ys[obs_vars]
      ## Demean Y
      Y_dm <- Y - d_obs   # n_obs x T (broadcasting over columns)

      ## --- CSN shock parameters + the contemporaneous state space ---
      ## (Estimated shock skewness arrives through `params`.) TT/ZZ come from
      ## .pskf_order1_statespace() via csn: the lagged-convention pair
      ## (ghx[state, ], ghx[obs, ]) with eps = DD e is NOT the model's law
      ## when DD != 0 -- see the note at .pskf_order1_statespace.
      csn <- tryCatch(
        .get_csn_shock_params(model, dr$exo_names, obs_vars, dr, params,
                              me_variance),
        error = function(e) .dynhr_reraise_bug(e, NULL)
      )
      if (is.null(csn)) return(NULL)

      ## --- Run PSKF ---
      ll <- tryCatch(
        .pskf_filter(
          Y         = Y_dm,
          TT        = csn$TT,
          ZZ        = csn$ZZ,
          mu_eta    = csn$mu_eta,
          Sigma_eta = csn$Sigma_eta,
          Gamma_eta = csn$Gamma_eta,
          nu_eta    = csn$nu_eta,
          Delta_eta = csn$Delta_eta,
          mu_eps    = csn$mu_eps,
          Sigma_eps = csn$Sigma_eps,
          cut_tol   = cut_tol,
          max_q     = max_q,
          offset_miwa_qmax = cdf$offset_miwa_qmax
        ),
        error = function(e) .dynhr_reraise_bug(e, -Inf)
      )
      if (!is.finite(ll)) return(NULL)
      list(loglik = ll)
    },
    power             = power,
    needs_me_floor    = FALSE,
    system_prior      = system_priors,
    system_prior_fn   = .pskf_system_prior_sum,
    system_prior_mode = "extra")
}


## The PSKF CDF-evaluation setting (package option `pskf_cdf`, resolved once
## per factory). The only value is "accurate": pruning compensation with Phi
## evaluated up to dimension 5, every Phi of dimension 2 by log-scale
## quadrature and of dimension 3-5 by the C++ lattice evaluator -- see
## logcdf_ME_r. A cheaper Mendell-Elston / plain-Miwa setting ("fast") existed
## up to 0.9.4.23; it ran the interpreted R filter and was slower than the C++
## default as well as less accurate, and was removed in 0.9.4.24.
## Returns list(offset_miwa_qmax).
#' @noRd
.pskf_cdf_settings <- function(pskf_cdf = NULL) {
  .pskf_cdf_check(.dynhr_opt("pskf_cdf", pskf_cdf))
  list(offset_miwa_qmax = 5L)
}

## Validator for a `pskf_cdf` value (argument, package option or spec field):
## aborts unless the value is "accurate"; NULL otherwise, so it doubles as a
## spec-field check. "fast" gets its own message: it was accepted up to
## 0.9.4.23 and must not be silently mapped to "accurate".
#' @noRd
.pskf_cdf_check <- function(pskf_cdf) {
  if (identical(pskf_cdf, "fast"))
    .dynhr_abort("pskf_cdf = \"fast\" was removed in dynhr 0.9.4.24. ",
                 "Use \"accurate\" (the default): it runs the C++ PSKF filter, ",
                 "so it is both faster and more accurate than the removed ",
                 "Mendell-Elston / plain-Miwa setting.",
                 class = "dynhr_error_pskf")
  if (!is.character(pskf_cdf) || length(pskf_cdf) != 1L ||
      !identical(pskf_cdf, "accurate"))
    .dynhr_abort("pskf_cdf must be \"accurate\", not ",
                 paste(format(pskf_cdf), collapse = ", "), ".",
                 class = "dynhr_error_pskf")
  invisible(NULL)
}


## System-prior evaluator for the PSKF pair.
##
## Unlike the Kalman / pruned / whittle / cumulant branches (which take a
## `system_prior_spec` and go through `.eval_system_priors()`), both PSKF
## factories accept a BARE LIST of `function(dr, params)` closures and sum
## them, stopping at the first non-finite partial sum. Kept as its own hook
## rather than unified: the two shapes are different public contracts, and
## test-posterior-closure-parity.R pins both.
#' @noRd
.pskf_system_prior_sum <- function(spec, theta, sol, params, res) {
  lsp <- 0
  if (!is.null(spec) && length(spec) > 0L) {
    for (fn in spec) {
      v   <- tryCatch(fn(sol$dr, params), error = function(e) .dynhr_reraise_bug(e, -Inf))
      lsp <- lsp + v
      if (!is.finite(lsp)) break
    }
  }
  lsp
}


## ===========================================================================
## PSKF ON THE PRUNED ORDER-2 STATE SPACE
## ===========================================================================
## The order-2 AFVRR augmented system (R/pruned-state-space.R) is LINEAR in the
## augmented state xi_t = [x1; x2; x1(x)x1]:
##   xi_{t+1} = Tlin xi_t + c_drift + G  r_t
##   y_t      = Dxi  xi_t + d_y     + Gv r_t
## with the single "raw innovation" r_t = [eps; eps(x)x1; x1(x)eps; eps(x)eps]
## driving BOTH the state (via G) and the observation (via Gv) -- i.e. the state
## and observation noise are CORRELATED (Cov = G Cr0 Gv' =: SS != 0, verified
## ~23 on rbc2shock).  .pskf_filter has no cross-covariance slot; it assumes
## eta_t (state) INDEPENDENT of eps_t (obs).
##
## RESOLUTION (exact, keeps .pskf_filter verbatim): carry the raw innovation
## r_t IN the augmented state.  Define chi_t = [xi_t; r_t] (dim d + Dr).  Then
##   chi_{t+1} = [ Tlin  G ] chi_t + [ 0 ] r_{t+1}
##               [  0    0 ]         [ I ]
##   y_t       = [ Dxi  Gv ] chi_t                     (NOISE-FREE observation!)
## The observation feedthrough Gv r_t is now a deterministic map of the state
## block r_t, so the state noise (r_{t+1}, loading [0; I]) is uncorrelated with
## the -- now zero -- observation noise.  This is EXACTLY the .pskf_filter form
##   chi_t = TT chi_{t-1} + eta_t,   y_t = ZZ chi_t + eps_t (Sigma_eps = 0),
## with TT = [[Tlin, G],[0,0]], ZZ = [Dxi, Gv], RR = [0; I_Dr], DD = 0.
## The innovation Omega = ZZ Sigma_pred ZZ' stays PD because the r-block of the
## state carries the full Cr0 (Omega picks up Gv Cr0 Gv' = HH > 0).
##
## SKEW LIFT: the raw-innovation covariance is Cr0 (dim Dr); the CSN skew lives
## ONLY on its first n_u components (the eps block).  Feed the generic lift
## .csn_state_noise_lift(RR, DD, Sigma_e = Cr0, alpha = [alpha_shocks; 0..0]):
## the eps rows carry Gamma_e = diag(alpha/sigma_e) and the mean correction; the
## Gaussian/quadratic augmentation rows (eps(x)x1, x1(x)eps, eps(x)eps) all have
## alpha = 0, so their Gamma_e rows and mean corrections vanish -- exactly the
## "original shocks carry the skew, augmentation blocks are Gaussian."
##
## TIMING / DEMEANING (for exact Gaussian-limit parity with pruned_ss_loglik):
## .pskf_filter inits mu_filt = 0, Sigma_filt = P0 (Lyapunov fixed point of
## (TT, Sigma_eta)) then does ONE predict before the first update, so its t=1
## prediction is the augmented STATIONARY (mean 0, cov P0).  pruned_ss_loglik
## instead inits directly at the stationary (mu0, Sxi0).  To align, we work in
## DEVIATIONS from the stationary mean (state mean 0, matching .pskf_filter's
## init exactly) and subtract the Gaussian stationary observation mean
## d_y + Dxi mu0 (mu0 = (I - Tlin)^{-1} c_drift) from Y.  Under skew the filter's
## own mu_eta keeps E[eta] = 0, so subtracting the Gaussian mean stays correct.
## Verified: alpha=0 parity vs pruned_ss_loglik = O(1e-13) (test-pskf-order2.R).

## Assemble the augmented CSN order-2 state space from a pruned_ss object.
## Returns list(TT, ZZ, RR, DD, Cr0, obs_mean, n_u, Dr) for the demeaned filter.
#' @noRd
.pskf_order2_augment <- function(pss, obs_vars) {
  sys     <- pss$sys
  d       <- sys$d
  n_u     <- sys$n_u
  obs_idx <- match(obs_vars, pss$endo_names)
  if (any(is.na(obs_idx)))
    stop(".pskf_order2_augment: obs_vars not in pss: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  n_obs <- length(obs_vars)

  Tlin <- sys$Tlin                        # d x d
  G    <- sys$G                           # d x Dr
  Dr   <- ncol(G)
  Dxi  <- sys$Dxi[obs_idx, , drop = FALSE]  # n_obs x d
  Gv   <- sys$Gv[obs_idx,  , drop = FALSE]  # n_obs x Dr

  ## Stationary raw-innovation covariance Cr0 (a = 0, P = Sigma_x)
  Sigma_x <- solve_lyapunov(sys$hx, sys$hu %*% sys$Sigma_e %*% t(sys$hu))
  Cr0     <- .order2_cov_r(numeric(sys$n_s), Sigma_x, sys$Sigma_e)

  da <- d + Dr
  TT <- matrix(0, da, da)
  TT[seq_len(d), seq_len(d)]        <- Tlin
  TT[seq_len(d), (d + 1L):da]       <- G
  ZZ <- matrix(0, n_obs, da)
  ZZ[, seq_len(d)]                  <- Dxi
  ZZ[, (d + 1L):da]                 <- Gv
  RR <- rbind(matrix(0, d, Dr), diag(Dr))   # da x Dr  (noise into r block only)
  DD <- matrix(0, n_obs, Dr)                 # no direct obs feedthrough

  ## Gaussian stationary observation mean (demean target):
  ##   d_y = ys + 0.5*ghss + c_v ;  mu0 = (I - Tlin)^{-1} c_drift
  ##   obs_mean = Dxi mu0 + d_y
  c_drift  <- sys$cc + sys$c_u
  mu0      <- as.numeric(solve(diag(d) - Tlin, c_drift))
  d_y      <- pss$ys[obs_vars] + 0.5 * sys$ghss[obs_idx] + sys$c_v[obs_idx]
  obs_mean <- as.numeric(Dxi %*% mu0) + d_y

  list(TT = TT, ZZ = ZZ, RR = RR, DD = DD, Cr0 = Cr0,
       obs_mean = obs_mean, n_u = n_u, Dr = Dr)
}


#' Create a PSKF log-posterior evaluator on the pruned ORDER-2 state space
#'
#' Factory (mirror of \code{make_log_posterior_pskf}) that runs the
#' closed-skew-normal (CSN) Kalman filter on the AFVRR (2018) pruned
#' second-order augmented state space, so skewed shocks propagate through the
#' second-order dynamics.  The augmented system is linear in the augmented
#' state, so the SAME CSN machinery (\code{.pskf_filter} + the shock-CSN lift
#' \code{.csn_state_noise_lift}) applies once the raw innovation \eqn{r_t} is
#' carried in the state (decorrelating the state/observation noise; see the
#' derivation note in this file).
#'
#' At \code{alpha = 0} for every shock this reduces EXACTLY (to \eqn{O(10^{-13})})
#' to the Gaussian pruned-order-2 filter \code{\link{pruned_ss_loglik}}.
#'
#' @param model       dynhr_mod (from \code{parse_mod}).
#' @param data        Observation matrix (T x n_obs or n_obs x T). A square
#'   matrix (T equal to the number of observables) is oriented by its
#'   dimnames: column names equal to \code{obs_vars} mean T x n_obs, row
#'   names equal to \code{obs_vars} mean n_obs x T; an unnamed square matrix
#'   is an error.
#' @param prior_spec  Prior specification (from \code{extract_prior_spec}).
#' @param obs_vars    Character vector of observed variable names.
#' @param compiled    dynhr_compiled (from \code{compile_model}, order >= 2).
#' @param me_variance Measurement-error variance (default 0): a scalar, or
#'   one variance per observable (\eqn{H = diag(me)}).
#' @param system_priors Named list of system-prior closures, or NULL.
#' @param cut_tol     Pruning tolerance (default 0.01; see \code{.pskf_filter}).
#' @param max_q       Hard cap on the retained skew dimension (default 5, the
#'   deterministic Phi_q evaluator range).
#' @param power       Power-posterior exponent applied to the LIKELIHOOD only
#'   (prior and system priors stay untempered). \code{NULL} (default) resolves
#'   the \code{power_posterior} package option ONCE, at factory time, so every
#'   draw evaluated by the returned closure uses the same tempering -- an
#'   option flipped mid-chain can no longer silently change the target
#'   distribution between draws.
#' @param pskf_cdf    Evaluation of the multivariate-normal CDFs of the filter:
#'   only \code{"accurate"} is accepted (\code{"fast"} was removed in 0.9.4.24;
#'   see \code{\link{dynhr_set_options}}).
#'   \code{NULL} (default) resolves the \code{pskf_cdf} package option once,
#'   at factory time.
#' @param ...         Ignored (interface compatibility).
#' @return function(theta) -> list(logpost, loglik, logprior)
#' @export
make_log_posterior_pskf_order2 <- function(model, data, prior_spec, obs_vars,
                                            compiled, me_variance = 0,
                                            system_priors = NULL,
                                            cut_tol = 0.01,
                                            max_q = 5L,
                                            power = NULL,
                                            pskf_cdf = NULL,
                                            ...) {
  ## Resolved ONCE at factory time -- see make_log_posterior_pskf().
  power <- .dynhr_opt("power_posterior", power, default = 1)
  cdf   <- .pskf_cdf_settings(pskf_cdf)
  me_variance <- .kf_me_variance(me_variance, obs_vars,
                                 "make_log_posterior_pskf_order2")
  ## Data orientation: either T x n_obs or n_obs x T -> Y is n_obs x T
  Y <- .pskf_orient_data(data, obs_vars, "make_log_posterior_pskf_order2")

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  ## Adapter over the shared closure builder (R/posterior-closure.R); the
  ## order-2 twin of make_log_posterior_pskf's. The order-2 lift, the pruned-SS
  ## object and the augmented CSN state space go in the SOLVE hook (they are
  ## what this branch solves for); the skew lift and the filter go in the
  ## loglik hook. Same bare-list system-prior contract, evaluated against dr2.
  .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    solve_fn = function(model, compiled, sys_cache, ss, params, theta) {
      ## Order-1 solve first: BK + stationarity guard, and the order-2 input.
      s1 <- .posterior_solve1(model, compiled, sys_cache, ss, params,
                              "spectral")
      if (is.null(s1)) return(NULL)
      dr2 <- tryCatch(
        solve_perturbation_order2(model, compiled, ss, params, s1$dr,
                                  verbose = FALSE),
        error = function(e) .dynhr_reraise_bug(e, NULL)
      )
      if (is.null(dr2)) return(NULL)
      pss <- tryCatch(pruned_state_space(dr2, model, params),
                      error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(pss)) return(NULL)
      aug <- tryCatch(.pskf_order2_augment(pss, obs_vars),
                      error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(aug)) return(NULL)
      list(dr = dr2, pss = pss, aug = aug)
    },
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      aug <- sol$aug
      ## Per-shock skewness lifted onto the raw innovation r_t: only the first
      ## n_u components (the eps block) carry alpha; augmentation blocks = 0.
      alpha_shocks <- .get_shock_skewness(model, sol$dr$exo_names, params)
      alpha_aug    <- c(as.numeric(alpha_shocks), rep(0, aug$Dr - aug$n_u))

      csn <- tryCatch(
        .csn_state_noise_lift(aug$RR, aug$DD, aug$Cr0, alpha_aug, me_variance),
        error = function(e) .dynhr_reraise_bug(e, NULL)
      )
      if (is.null(csn)) return(NULL)

      ## Demean Y by the Gaussian stationary observation mean (see timing note)
      Y_dm <- Y - aug$obs_mean

      ll <- tryCatch(
        .pskf_filter(
          Y         = Y_dm,
          TT        = aug$TT,
          ZZ        = aug$ZZ,
          mu_eta    = csn$mu_eta,
          Sigma_eta = csn$Sigma_eta,
          Gamma_eta = csn$Gamma_eta,
          nu_eta    = csn$nu_eta,
          Delta_eta = csn$Delta_eta,
          mu_eps    = csn$mu_eps,
          Sigma_eps = csn$Sigma_eps,
          cut_tol   = cut_tol,
          max_q     = max_q,
          offset_miwa_qmax = cdf$offset_miwa_qmax
        ),
        error = function(e) .dynhr_reraise_bug(e, -Inf)
      )
      if (!is.finite(ll)) return(NULL)
      list(loglik = ll)
    },
    power             = power,
    needs_me_floor    = FALSE,
    system_prior      = system_priors,
    system_prior_fn   = .pskf_system_prior_sum,
    system_prior_mode = "extra")
}
