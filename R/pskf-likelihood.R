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
## CRITICAL CONVENTIONS (brief landmines):
##  1. ghu EXCLUDES Sigma_e (it is the loading matrix R in x_t = G x_{t-1} + R e_t).
##     Sigma_eta = RR Sigma_e RR'; Sigma_eps = DD Sigma_e DD' + me_variance*I.
##  4. Sigma_eps is assembled as ONE matrix (HH + me_variance*I).
##  6. mu_eta is mean-corrected so E[eta] = 0: for skew-normal with shape alpha
##     and std sigma, E[eta] = sigma * delta * sqrt(2/pi), so mu_eta = -E[eta].
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Deterministic multivariate-normal CDF
## ---------------------------------------------------------------------------
## logcdf_ME_r(x, S) computes log Phi_q(x; 0, S).
##
## Strategy:
##   q = 1: pnorm() -- exact.
##   q = 2: stats::integrate() over the conditional 1-D integral -- deterministic,
##          accurate to ~1e-10.
##   3 <= q <= 5: mvtnorm::pmvnorm(algorithm = Miwa()) when mvtnorm is available.
##          Miwa (2004) is fully deterministic (no random seed); accuracy ~1e-10.
##          Falls back to Mendell-Elston if mvtnorm is unavailable.
##   q > 5: Mendell-Elston (1974) sequential conditioning approximation.
##          First-order moment-matching; absolute error ~1e-3 (near-diagonal S)
##          to ~0.05 (high correlations).  Miwa's O(n_pts^q) cost is impractical
##          for q > 5 so ME is the fallback.  See inline accuracy note below.
##
## Why Miwa over GenzBretz?  GenzBretz (mvtnorm default) uses randomised QMC --
## non-reproducible across calls, which breaks gradient-based mode-finding and
## SBC.  Miwa() IS deterministic and exact for q <= 5.  A tryCatch wraps the
## Miwa call so that non-PD edge-cases fall back to ME without crashing.
##
## miwa_qmax: largest dimension evaluated with the deterministic-exact Miwa
## algorithm before falling back to Mendell-Elston.  The 2026-07-03 pruning-
## bias investigation (scratchpad/pskf-*.R; memory note
## pskf-multishock-pruning-bias) proved with an exact 2-D grid-filter oracle
## that the entire multi-shock "pruning" bias (-7.3 nats at T=12 for
## alpha=(+2,-2)) was Mendell-Elston evaluation error at q>5, NOT discarded
## skew mass: swapping ME for an exact Phi_q at unchanged cut_tol=0.01
## collapsed the gap to |0.005|.  Keeping q inside the Miwa-exact range via
## rank-capped pruning (see dim_red4_r max_q) is therefore the fix.
#' @noRd
logcdf_ME_r <- function(x, S, miwa_qmax = 5L) {
  q <- length(x)
  if (q == 0L) return(0)

  ## q = 1: exact
  if (q == 1L) {
    b <- x[1] / sqrt(S[1, 1])
    return(pnorm(b, log.p = TRUE))
  }

  ## ---- SNAP + BLOCK FACTORIZATION (2026-07 Miwa-pocket fix) ---------------
  ## mvtnorm's Miwa algorithm has an instability pocket for TINY-but-nonzero
  ## correlations: measured on a well-conditioned 3x3 with mixed-sign
  ## rho ~ 1e-5..1e-3 it returned Phi_3 = 1.0856 (> 1!) at steps = 128 and was
  ## still 0.04 absolute off at steps = 512, while rho = 0, 1e-7 and 0.2 are
  ## all ~1e-6 accurate (scratchpad m5_trace_cdf.R, Reiter-HANK PSKF
  ## investigation -- this made the filter's likelihood IMPROPER).  Fix:
  ## (a) standardize to correlation form and SNAP |rho| < 1e-3 to exactly 0
  ##     (error bound: |dPhi/drho| = phi_2 <= 1/(2*pi) per pair, so <= ~1.6e-4
  ##     per snapped pair -- strictly smaller than Miwa's own pocket error);
  ## (b) factor the CDF over the connected components of the snapped
  ##     correlation graph (EXACT given the snap) -- singletons/pairs then use
  ##     the exact pnorm/bivariate paths and Miwa only sees well-coupled
  ##     blocks.
  sdv <- sqrt(pmax(diag(as.matrix(S)), .Machine$double.eps))
  Cm  <- as.matrix(S) / outer(sdv, sdv)
  off <- row(Cm) != col(Cm)
  Cm[off & abs(Cm) < 1e-3] <- 0
  x <- as.numeric(x) / sdv
  S <- Cm
  adj  <- Cm != 0
  comp <- integer(q); n_comp <- 0L
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
  if (n_comp > 1L) {
    out <- 0
    for (cc in seq_len(n_comp)) {
      ii  <- which(comp == cc)
      out <- out + logcdf_ME_r(x[ii], Cm[ii, ii, drop = FALSE], miwa_qmax)
    }
    return(out)
  }

  ## q = 2: exact via stats::integrate (deterministic, accurate ~1e-10)
  if (q == 2L) {
    ## Standardise bounds to correlation form
    sd1  <- sqrt(S[1, 1])
    sd2  <- sqrt(S[2, 2])
    h1   <- x[1] / sd1
    h2   <- x[2] / sd2
    rho  <- S[1, 2] / (sd1 * sd2)
    rho  <- max(-0.9999, min(0.9999, rho))   # guard |rho| = 1
    sr   <- sqrt(1 - rho^2)
    ## P(X1<=h1, X2<=h2) = integral_{-inf}^{h1} phi(u) Phi((h2 - rho*u)/sr) du
    f <- function(u) dnorm(u) * pnorm((h2 - rho * u) / sr)
    val <- tryCatch(
      integrate(f, lower = -8, upper = h1, rel.tol = 1e-8,
                subdivisions = 100L)$value,
      error = function(e) {
        ## Fallback: product approximation Phi(h1)*Phi(h2) (lower bound)
        pnorm(h1) * pnorm(h2)
      }
    )
    val <- max(val, .Machine$double.eps)
    return(log(val))
  }

  ## 3 <= q <= 5: use deterministic Miwa algorithm via mvtnorm if available.
  ## Miwa (2004) is exact for any q; O(n_pts^q) cost is practical for q <= 5.
  ## A tryCatch guards against non-PD edge cases; those fall through to ME.
  if (q <= miwa_qmax && requireNamespace("mvtnorm", quietly = TRUE)) {
    val <- tryCatch(
      mvtnorm::pmvnorm(upper = as.numeric(x), sigma = as.matrix(S),
                       algorithm = mvtnorm::Miwa())[1L],
      error = function(e) NA_real_
    )
    ## IMPOSSIBLE-VALUE guard (Miwa-pocket fix, part c): Miwa can return
    ## probabilities > 1 or <= 0 WITHOUT erroring on pathological inputs --
    ## treat those like an error and fall through to the deterministic ME
    ## approximation instead of corrupting the loglik.
    if (!is.na(val) && is.finite(val) && val > 0 && val <= 1 + 1e-8) {
      return(log(max(min(val, 1), .Machine$double.eps)))
    }
    ## Fall through to ME on error / impossible value (non-PD S, etc.)
  }

  ## q > 5 (or Miwa unavailable / errored): Mendell-Elston sequential conditioning.
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
  ## MEASURED (2026-07-03, scratchpad/step3-evaluators.R): under strong
  ## multi-shock skew (2 shocks, |alpha| = 2-3) the CSN correction arguments
  ## sit deep in the orthant tails and the ME error per call reaches SEVERAL
  ## NATS at q = 6-16, accumulating to -7.3 nats by T=12 for alpha=(+2,-2)
  ## (sign follows the sign of the off-diagonal correlations of S).
  ## The filter therefore rank-caps q at max_q = 5 (dim_red4_r), so this
  ## branch is only reached when the caller explicitly raises max_q (or
  ## mvtnorm is unavailable); the loglik is then approximate and biased
  ## under strong skew -- see memory note pskf-multishock-pruning-bias.
  ##
  ## Algorithm (Mendell & Elston 1974, corrected for covariance -- not correlation -- form):
  ## For j = 1 .. q-1:
  ##   standardised bound: bj = b[j] / sqrt(S[j,j])
  ##   P_j = Phi(bj)  (marginal probability for dim j)
  ##   Mills ratio lambda = phi(bj)/Phi(bj)
  ##   mean shift for k > j:  b[k] -= S[j,k]/sqrt(S[j,j]) * lambda
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
      b[k] <- b[k] - CS[j, k] / sj * lambda
      ## Covariance shrinkage
      for (l in idx) {
        CS[k, l] <- CS[k, l] - CS[j, k] * CS[j, l] / sjj * delta_factor
        CS[l, k] <- CS[k, l]
      }
    }
  }

  bq    <- b[q] / sqrt(max(CS[q, q], .Machine$double.eps))
  log_p <- log_p + pnorm(bq, log.p = TRUE)

  log_p
}


## ---------------------------------------------------------------------------
## Pruning: dim_red4_r
## ---------------------------------------------------------------------------
## Keeps only the rows of (Gamma, nu) whose maximum absolute correlation with
## any other row of Gamma Sigma Gamma' + Delta is >= cut_tol.
## Without pruning the skewness dimension q grows by n_exo every period;
## pruning keeps q bounded (typically 1-3 for small DSGE, brief Landmine 2).
##
## max_q: HARD rank-based cap on the retained skew dimension.  After the
## cut_tol threshold filter, if more than max_q rows survive, only the max_q
## rows with the LARGEST skew-vs-state correlation are kept (in original
## order).  Rationale (2026-07-03 investigation, scratchpad/pskf-*.R): the
## Phi_q evaluator is deterministic-exact (Miwa) only for q <= 5; beyond that
## the Mendell-Elston approximation's error GROWS with the per-period loglik
## contribution and accumulates (measured -7.3 nats at T=12, 2 shocks,
## alpha=(+2,-2)).  Rank-capping q at 5 keeps every CDF call inside the
## exact-evaluator range; the discarded low-correlation skew mass costs far
## less (|gap| <= ~0.3 nat at T=12 on the worst measured fixture) than the
## ME error it avoids.
##
## Mean offset of a CSN(0, Sigma, Gamma, nu, Delta) relative to its Gaussian
## location:  E[X] - mu = Sigma Gamma' g,
##   g_j = phi(-nu_j; V_jj) * Phi_{q-1}(cond_j) / Phi_q(-nu; V),
##   V = Delta + Gamma Sigma Gamma',
## (gradient of the log-normaliser wrt nu; verified against rejection-sampling
## MC, 2026-07 pruning-mean-drift investigation). The conditional CDFs use the
## cheap ME path (miwa_qmax = 2): errors largely cancel in the ratio and the
## offset is itself a correction term -- Miwa-5 here would dominate filter
## cost under saturated pruning (one cut per period).
#' @noRd
.csn_mean_offset <- function(Gamma, nu, Delta, Sigma) {
  q <- nrow(Gamma)
  if (q == 0L) return(rep(0, ncol(Gamma)))
  V <- Delta + Gamma %*% Sigma %*% t(Gamma)
  V <- (V + t(V)) / 2
  logZ <- logcdf_ME_r(-nu, V, miwa_qmax = 2L)
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
      lcond <- logcdf_ME_r(mcond, Scond, miwa_qmax = 2L)
    }
    g[j] <- exp(lphi + lcond - logZ)
  }
  if (!all(is.finite(g))) return(rep(0, ncol(Gamma)))
  as.numeric(Sigma %*% t(Gamma) %*% g)
}

## Returns list(Gamma = q'xp, nu = q', Delta = q'xq', mu_shift = p-vector)
## with q' <= q.  mu_shift is the FIRST-MOMENT COMPENSATION for the cut:
## deleting a skew row removes that dimension's contribution to the CSN mean
## (first-order in its skew-state correlation), so the caller must add
## mu_shift = offset(before) - offset(after) to the Gaussian location.
## Without it, saturated pruning (rank cap binding every period, e.g. a
## persistent single-shock model with |alpha| large) accumulates a systematic
## state-mean drift that makes the likelihood IMPROPER (one-step predictive
## densities integrating to 0.03-0.97) -- caught by the alpha_z SBC on the
## Reiter HANK state space (2026-07; coverage collapsed to 3-6%).
#' @noRd
dim_red4_r <- function(Gamma, nu, Delta, Sigma, cut_tol = 0.01, max_q = 5L) {
  q <- nrow(Gamma)
  no_shift <- rep(0, ncol(Gamma))
  if (q == 0L)
    return(list(Gamma = Gamma, nu = nu, Delta = Delta, mu_shift = no_shift))

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
  max_corr <- apply(corr_sx, 1, max)

  keep_idx <- which(max_corr >= cut_tol)
  ## COLLINEARITY guard (2026-07, Reiter-HANK investigation): when the skew
  ## rows are propagated images of the SAME shock direction (persistent
  ## single-shock models), corr among skew dims -> 1 within a few periods.
  ## Near-singular V = Delta + G S G' breaks the Phi_q evaluators (Miwa
  ## errors -> ME fallback whose ~0.03 ABSOLUTE error on ~1e-8 tail CDFs is
  ## tens of nats in the top-minus-bottom log difference -> IMPROPER
  ## likelihood, one-step predictive integrals 0.04-1e14). Nearly-duplicate
  ## constraints are nearly-free to drop UNDER MEAN COMPENSATION (below), so
  ## iteratively drop the weaker row of any pair with |corr| > 0.995 -- this
  ## keeps V numerically nonsingular and the Miwa path exact.
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
  if (length(keep_idx) == q) {
    return(list(Gamma = Gamma, nu = nu, Delta = Delta, mu_shift = no_shift))
  }
  off_before <- .csn_mean_offset(Gamma, nu, Delta, Sigma)
  if (length(keep_idx) == 0L) {
    return(list(
      Gamma    = matrix(0, nrow = 0L, ncol = ncol(Gamma)),
      nu       = numeric(0),
      Delta    = matrix(0, nrow = 0L, ncol = 0L),
      mu_shift = off_before
    ))
  }
  G_k <- Gamma[keep_idx, , drop = FALSE]
  n_k <- nu[keep_idx]
  D_k <- Delta[keep_idx, keep_idx, drop = FALSE]
  list(
    Gamma    = G_k,
    nu       = n_k,
    Delta    = D_k,
    mu_shift = off_before - .csn_mean_offset(G_k, n_k, D_k, Sigma)
  )
}


## ---------------------------------------------------------------------------
## Assemble CSN shock parameters from model and current params
## ---------------------------------------------------------------------------
## Returns list with:
##   Sigma_eta : n_state x n_state = RR Sigma_e RR'
##   Sigma_eps : n_obs   x n_obs   = DD Sigma_e DD' + me_variance * I
##   Gamma_eta : q_eta x n_state   (skewness loading through RR)
##   nu_eta    : q_eta vector (= 0)
##   Delta_eta : q_eta x q_eta (= I)
##   mu_eta    : n_state vector (mean correction, = 0 when all alpha_i = 0)
##   mu_eps    : n_obs vector   (= 0)
##   alpha     : n_exo named vector of shape parameters
##
## DERIVATION (CSN linear-transform convention, brief Landmine 1):
##   The shock vector in model space is e ~ skewNormal(alpha_i, sigma_i) for
##   each component i.  Each scalar skew-normal can be written as:
##     e_i = mu_i + sigma_i * Z_i,
##       where Z_i ~ CSN(0, 1, alpha_i, 0, 1) (unit skew-normal).
##   Collecting all shocks: e ~ CSN(mu_e, Sigma_e, Gamma_e, 0, I_{n_exo}) where
##     mu_e    = -E[e_i] diagonal corrections (= -sigma_i delta_i sqrt(2/pi))
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
#' @noRd
.get_csn_shock_params <- function(model, exo_names, obs_vars, dr, params,
                                   me_variance = 0) {
  ## Extract state-space matrices (mirrors kalman_filter convention)
  state_idx <- dr$state_idx
  obs_idx   <- match(obs_vars, dr$endo_names)
  RR <- dr$ghu[state_idx, , drop = FALSE]   # n_state x n_exo
  DD <- dr$ghu[obs_idx,   , drop = FALSE]   # n_obs   x n_exo

  n_state <- nrow(RR)
  n_exo   <- ncol(RR)
  n_obs   <- nrow(DD)

  Sigma_e <- .get_shock_cov(model, exo_names, params)  # n_exo x n_exo
  alpha   <- .get_shock_skewness(model, exo_names, params)  # n_exo named vector

  ## Correlated + skewed shocks (Tier 10 item 3): the FULL joint CSN.
  ##
  ## The joint shock law is  e ~ CSN(mu_e, Sigma_e, Gamma_e, 0, I_{n_exo})  with
  ##   Gamma_e = diag(alpha_i / sigma_i),   Delta_e = I,   seed cov = full Sigma_e
  ## (off-diagonals included).  This is the construction pinned numerically in
  ## Phase 0 of the v2-csn brief.  Three properties were verified by brute force
  ## (1e6-draw CSN rejection sampler, 2-shock alpha = (+2, -2)):
  ##   (a) rho_12 = 0  collapses Sigma_e to diagonal, so the build is BIT-IDENTICAL
  ##       to the old per-shock-independent path (A3.5; verified delta == 0).
  ##   (b) Delta_e = I is always PD, so the joint law is a valid CSN for any
  ##       admissible Sigma_e (no PD trap; the brief's naive coupling
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
  ## the DGP and likelihood share the same law (A3.3 consistency).

  .csn_state_noise_lift(RR, DD, Sigma_e, alpha, me_variance)
}


#' Scale-free Moore-Penrose pseudoinverse of a symmetric PSD matrix
#'
#' Used for \eqn{(RR \Sigma_e RR')^{+}} in \code{.csn_state_noise_lift}.
#' The rank decision is made on the CORRELATION form
#' \eqn{D^{-1/2} S D^{-1/2}} (\eqn{D = diag(S)}) rather than on \code{S}
#' itself, so a state whose variance is genuinely small but nonzero (a shock
#' with \code{stderr = 1e-5} alongside one with \code{stderr = 1}) is not
#' truncated as a null direction the way \code{MASS::ginv} would. Exactly
#' zero-variance rows/columns ARE null directions and are dropped outright.
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
  d <- sqrt(diag(S))
  d[!is.finite(d)] <- 0
  keep <- d > 0
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

  ## State noise covariance (Landmine 1: ghu excludes Sigma_e)
  Sigma_eta <- RR %*% Sigma_e %*% t(RR)   # n_state x n_state

  ## Obs noise covariance (Landmine 4: assembled as ONE matrix)
  Sigma_eps <- DD %*% Sigma_e %*% t(DD) + me_variance * diag(n_obs)

  ## Mean correction (Landmine 6): E[e_i] = sigma_i * delta_i * sqrt(2/pi)
  ## where delta_i = alpha_i / sqrt(1 + alpha_i^2).
  ## mu_eta is chosen s.t. E[eta] = 0 (preserves model steady state).
  sigma_e   <- sqrt(diag(Sigma_e))           # n_exo  (ignores off-diagonal for mean)
  delta_e   <- alpha / sqrt(1 + alpha^2)     # n_exo
  mean_e    <- sigma_e * delta_e * sqrt(2 / pi)  # E[e_i], n_exo
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
## PSKF filter recursion
## ---------------------------------------------------------------------------
## Implements the CSN Kalman filter (PSKF) following skalman_filter.R /
## the brief's recursion spec.
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
## Gamma_pred_path, nu_pred_path, Delta_pred_path) — per-period arrays
## needed by the CSN backward pass in pskf_smoother().
## max_q (default 5): hard rank-based cap on the retained skew dimension,
## chosen so every Phi_q call stays inside the deterministic-EXACT Miwa range
## (q <= 5).  Verified 2026-07-03 (scratchpad/pskf-*.R + memory note
## pskf-multishock-pruning-bias): the historical multi-shock bias
## (-7.3 nats at T=12, alpha=(+2,-2)) was ENTIRELY Mendell-Elston Phi_q
## evaluation error at q > 5, not discarded skew mass; capping q at 5 cuts
## it to |gap| <= 0.14 nat on the worst measured fixture.  max_q = Inf
## restores the old uncapped behavior (cut_tol-threshold pruning only, ME
## for q > 5) -- NOT recommended for multi-shock skew.
#' @noRd
.pskf_filter <- function(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                          Delta_eta, mu_eps, Sigma_eps,
                          cut_tol = 0.01,
                          max_q = 5L,
                          store_path = FALSE) {
  ## Y: n_obs x T matrix
  n_obs   <- nrow(Y)
  n_T     <- ncol(Y)
  n_state <- nrow(TT)   # = ncol(TT)

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

  mu_filt    <- rep(0, n_state)
  Sigma_filt <- P0
  Gamma_filt <- matrix(0, nrow = 0L, ncol = n_state)
  nu_filt    <- numeric(0)
  Delta_filt <- matrix(0, nrow = 0L, ncol = 0L)

  ## Log-likelihood accumulator
  ll <- 0

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
  }

  for (t in seq_len(n_T)) {
    ## ---- PREDICTION STEP ----------------------------------------------------
    ## Gaussian part: standard Kalman prediction
    ##   mu_{t|t-1}    = TT mu_{t-1|t-1} + mu_eta
    ##   Sigma_{t|t-1} = TT Sigma_{t-1|t-1} TT' + Sigma_eta
    mu_pred    <- as.numeric(TT %*% mu_filt) + mu_eta
    Sigma_pred <- TT %*% Sigma_filt %*% t(TT) + Sigma_eta

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
      ## Cholesky of Sigma_pred (needed for solve)
      ## Use tryCatch to handle near-singularity gracefully
      S_pred_inv <- tryCatch(solve(S_pred), error = function(e) {
        S_pred_reg <- S_pred + 1e-8 * diag(n_state)
        solve(S_pred_reg)
      })

      ## Build the stacked Gamma (q_new x n_state)
      if (q_filt == 0L) {
        ## Only shock skewness
        G_eta_block <- Gamma_eta %*% Sigma_eta %*% S_pred_inv   # q_eta x n_state
        Gamma_pred  <- G_eta_block
      } else if (q_eta == 0L) {
        ## Only filtered skewness
        G_filt_block <- Gamma_filt %*% Sigma_filt %*% t(TT) %*% S_pred_inv
        Gamma_pred   <- G_filt_block
      } else {
        G_filt_block <- Gamma_filt %*% Sigma_filt %*% t(TT) %*% S_pred_inv  # q_filt x n_state
        G_eta_block  <- Gamma_eta  %*% Sigma_eta  %*% S_pred_inv             # q_eta  x n_state
        Gamma_pred   <- rbind(G_filt_block, G_eta_block)                     # q_new  x n_state
      }

      ## Stacked nu
      nu_pred <- c(nu_filt, nu_eta)

      ## Delta 4-block (conditional covariance structure):
      ## Using the general formula:
      ##   D11 = D_filt + G_filt (S_filt - S_filt TT' S_pred^{-1} TT S_filt) G_filt'
      ##   D22 = D_eta  + G_eta  (S_eta  - S_eta  S_pred^{-1}       S_eta)   G_eta'
      ##   D12 = -G_filt S_filt TT' S_pred^{-1} S_eta G_eta'
      ##   D21 = D12'
      ## The block structure is the Schur complement of S_pred in the
      ## joint CSN skewness covariance.

      if (q_filt == 0L) {
        ## Only eta block
        Schur_eta <- Sigma_eta - Sigma_eta %*% S_pred_inv %*% Sigma_eta
        Delta_pred <- Delta_eta + Gamma_eta %*% Schur_eta %*% t(Gamma_eta)
      } else if (q_eta == 0L) {
        ## Only filt block
        TT_Sfilt <- TT %*% Sigma_filt
        Schur_filt <- Sigma_filt %*% t(TT) %*% S_pred_inv %*% TT_Sfilt
        D11 <- Delta_filt + Gamma_filt %*%
               (Sigma_filt - Schur_filt) %*% t(Gamma_filt)
        Delta_pred <- D11
      } else {
        ## Full 4-block
        TT_Sfilt <- TT %*% Sigma_filt
        ## Schur contributions:
        ## S_filt - S_filt TT' S_pred^{-1} TT S_filt
        Sfilt_TT_t    <- Sigma_filt %*% t(TT)   # n_state x n_state
        Schur_filt    <- Sfilt_TT_t %*% S_pred_inv %*% TT %*% Sigma_filt
        D11 <- Delta_filt +
               Gamma_filt %*% (Sigma_filt - Schur_filt) %*% t(Gamma_filt)

        ## S_eta - S_eta S_pred^{-1} S_eta
        Schur_eta <- Sigma_eta - Sigma_eta %*% S_pred_inv %*% Sigma_eta
        D22 <- Delta_eta + Gamma_eta %*% Schur_eta %*% t(Gamma_eta)

        ## Off-diagonal: -G_filt S_filt TT' S_pred^{-1} S_eta G_eta'
        D12 <- -Gamma_filt %*% Sfilt_TT_t %*%
               S_pred_inv %*% Sigma_eta %*% t(Gamma_eta)

        ## Assemble block matrix
        Delta_pred <- rbind(
          cbind(D11, D12),
          cbind(t(D12), D22)
        )
      }
    }

    ## ---- PRUNE ---------------------------------------------------------------
    ## Keep q bounded (brief Landmine 2). Pruning after prediction before update.
    if (q_new > 0L) {
      pruned <- dim_red4_r(Gamma_pred, nu_pred, Delta_pred, Sigma_pred, cut_tol,
                           max_q = max_q)
      Gamma_pred <- pruned$Gamma
      nu_pred    <- pruned$nu
      Delta_pred <- pruned$Delta
      ## first-moment compensation for the cut skew mass (see dim_red4_r):
      ## without this, saturated pruning drifts the state mean systematically
      mu_pred <- mu_pred + pruned$mu_shift
    }
    q_pred <- nrow(Gamma_pred)

    ## ---- HANDLE MISSING OBSERVATIONS ----------------------------------------
    y_t <- Y[, t]
    obs_mask <- !is.na(y_t)
    if (!all(obs_mask)) {
      ## Partial observation: update only on available obs (mirrors Gaussian KF)
      if (!any(obs_mask)) {
        ## All missing: prediction becomes filtered
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
        }
        next
      }
      ZZ_t      <- ZZ[obs_mask,   , drop = FALSE]
      Sigma_eps_t <- Sigma_eps[obs_mask, obs_mask, drop = FALSE]
      mu_eps_t  <- mu_eps[obs_mask]
      y_t       <- y_t[obs_mask]
      n_obs_t   <- sum(obs_mask)
    } else {
      ZZ_t        <- ZZ
      Sigma_eps_t <- Sigma_eps
      mu_eps_t    <- mu_eps
      n_obs_t     <- n_obs
    }

    ## ---- UPDATE STEP ---------------------------------------------------------
    ## Gaussian innovation:
    ##   v_t   = y_t - ZZ mu_{t|t-1} - mu_eps
    ##   Omega = ZZ Sigma_{t|t-1} ZZ' + Sigma_eps  (innovation covariance)
    v_t   <- as.numeric(y_t - ZZ_t %*% mu_pred - mu_eps_t)
    Omega <- ZZ_t %*% Sigma_pred %*% t(ZZ_t) + Sigma_eps_t

    ## Cholesky of Omega for logdet and solve. chol() fails only when Omega is
    ## NOT positive-definite -- a singular/indefinite innovation covariance,
    ## which is an infeasible draw (e.g. no measurement error on a rank-deficient
    ## observable block). The former code jittered it (chol(Omega + 1e-8 I)) and
    ## proceeded, silently fabricating a finite loglik that the sampler then
    ## ACCEPTED (a wrongly-accepted, corrupted draw). Reject it instead.
    Omega_chol <- tryCatch(chol(Omega), error = function(e) NULL)
    if (is.null(Omega_chol)) { ll <- -Inf; break }
    log_det_Omega <- 2 * sum(log(diag(Omega_chol)))
    Omega_inv_v   <- backsolve(Omega_chol, forwardsolve(t(Omega_chol), v_t))

    ## Gaussian log-likelihood contribution: log N(v_t; 0, Omega)
    ll_gauss <- -0.5 * (n_obs_t * log(2 * pi) + log_det_Omega +
                        sum(v_t * Omega_inv_v))

    ## Kalman gain (Gaussian)
    K_gauss <- Sigma_pred %*% t(ZZ_t) %*% solve(Omega)

    ## CSN update (brief's update step):
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
    I_KZ      <- diag(n_state) - K_gauss %*% ZZ_t
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
    ## (brief mandatory test, brief p.2) is exact.

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
    Gamma_filt <- Gamma_pred   # unchanged (brief update rule)
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
      Delta_pred_path = Delta_pred_path
    )
  } else {
    ll
  }
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
#' @param cut_tol     Pruning tolerance (default 0.01; brief Landmine 2)
#' @param max_q       Hard cap on the retained skew dimension (default 5,
#'   the Miwa-exact Phi_q range; see .pskf_filter). Inf = old uncapped
#'   behavior (Mendell-Elston for q > 5; biased under strong multi-shock skew).
#' @param power       Power-posterior exponent applied to the LIKELIHOOD only
#'   (prior and system priors stay untempered). \code{NULL} (default) resolves
#'   the \code{power_posterior} package option ONCE, at factory time, so every
#'   draw evaluated by the returned closure uses the same tempering -- an
#'   option flipped mid-chain can no longer silently change the target
#'   distribution between draws.
#' @param ...         Ignored (for interface compatibility)
#' @return function(theta) -> list(logpost, loglik, logprior)
#' @noRd
make_log_posterior_pskf <- function(model, data, prior_spec, obs_vars,
                                     compiled, me_variance = 0,
                                     system_priors = NULL,
                                     cut_tol = 0.01,
                                     max_q = 5L,
                                     power = NULL,
                                     ...) {
  ## Resolve the power-posterior exponent ONCE, here, rather than on every
  ## evaluation: the closure's target must not change under the caller's feet.
  power <- .dynhr_opt("power_posterior", power, default = 1)
  ## Validate data orientation: need T x n_obs
  if (ncol(data) == length(obs_vars)) {
    Y <- t(data)   # -> n_obs x T
  } else {
    Y <- data       # assume already n_obs x T
  }

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
      ## --- Assemble state space ---
      state_idx <- dr$state_idx
      obs_idx   <- match(obs_vars, dr$endo_names)
      if (any(is.na(obs_idx))) return(NULL)

      TT <- dr$ghx[state_idx, , drop = FALSE]
      ZZ <- dr$ghx[obs_idx,   , drop = FALSE]
      ## Observation mean: d_obs (includes SS value and any DR offset)
      d_obs <- dr$ys[obs_vars]
      ## Demean Y
      Y_dm <- Y - d_obs   # n_obs x T (broadcasting over columns)

      ## --- CSN shock parameters ---
      ## (Estimated shock skewness arrives through `params`; for the v0 wiring
      ## the shocks block provides a fixed alpha.)
      csn <- tryCatch(
        .get_csn_shock_params(model, dr$exo_names, obs_vars, dr, params,
                              me_variance),
        error = function(e) NULL
      )
      if (is.null(csn)) return(NULL)

      ## --- Run PSKF ---
      ll <- tryCatch(
        .pskf_filter(
          Y         = Y_dm,
          TT        = TT,
          ZZ        = ZZ,
          mu_eta    = csn$mu_eta,
          Sigma_eta = csn$Sigma_eta,
          Gamma_eta = csn$Gamma_eta,
          nu_eta    = csn$nu_eta,
          Delta_eta = csn$Delta_eta,
          mu_eps    = csn$mu_eps,
          Sigma_eps = csn$Sigma_eps,
          cut_tol   = cut_tol,
          max_q     = max_q
        ),
        error = function(e) -Inf
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
      v   <- tryCatch(fn(sol$dr, params), error = function(e) -Inf)
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
## brief's "original shocks carry the skew, augmentation blocks are Gaussian."
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
#' @param data        Observation matrix (T x n_obs or n_obs x T).
#' @param prior_spec  Prior specification (from \code{extract_prior_spec}).
#' @param obs_vars    Character vector of observed variable names.
#' @param compiled    dynhr_compiled (from \code{compile_model}, order >= 2).
#' @param me_variance Measurement-error variance (default 0).
#' @param system_priors Named list of system-prior closures, or NULL.
#' @param cut_tol     Pruning tolerance (default 0.01; see \code{.pskf_filter}).
#' @param max_q       Hard cap on the retained skew dimension (default 5, the
#'   Miwa-exact Phi_q range).
#' @param power       Power-posterior exponent applied to the LIKELIHOOD only
#'   (prior and system priors stay untempered). \code{NULL} (default) resolves
#'   the \code{power_posterior} package option ONCE, at factory time, so every
#'   draw evaluated by the returned closure uses the same tempering -- an
#'   option flipped mid-chain can no longer silently change the target
#'   distribution between draws.
#' @param ...         Ignored (interface compatibility).
#' @return function(theta) -> list(logpost, loglik, logprior)
#' @export
make_log_posterior_pskf_order2 <- function(model, data, prior_spec, obs_vars,
                                            compiled, me_variance = 0,
                                            system_priors = NULL,
                                            cut_tol = 0.01,
                                            max_q = 5L,
                                            power = NULL,
                                            ...) {
  ## Resolved ONCE at factory time -- see make_log_posterior_pskf().
  power <- .dynhr_opt("power_posterior", power, default = 1)
  ## Validate data orientation: need n_obs x T
  if (ncol(data) == length(obs_vars)) {
    Y <- t(data)   # -> n_obs x T
  } else {
    Y <- data       # assume already n_obs x T
  }

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
        error = function(e) NULL
      )
      if (is.null(dr2)) return(NULL)
      pss <- tryCatch(pruned_state_space(dr2, model, params),
                      error = function(e) NULL)
      if (is.null(pss)) return(NULL)
      aug <- tryCatch(.pskf_order2_augment(pss, obs_vars),
                      error = function(e) NULL)
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
        error = function(e) NULL
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
          max_q     = max_q
        ),
        error = function(e) -Inf
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
