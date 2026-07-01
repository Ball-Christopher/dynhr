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
#' @noRd
logcdf_ME_r <- function(x, S) {
  q <- length(x)
  if (q == 0L) return(0)

  ## q = 1: exact
  if (q == 1L) {
    b <- x[1] / sqrt(S[1, 1])
    return(pnorm(b, log.p = TRUE))
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
  if (q <= 5L && requireNamespace("mvtnorm", quietly = TRUE)) {
    val <- tryCatch(
      mvtnorm::pmvnorm(upper = as.numeric(x), sigma = as.matrix(S),
                       algorithm = mvtnorm::Miwa())[1L],
      error = function(e) NA_real_
    )
    if (!is.na(val)) {
      return(log(max(val, .Machine$double.eps)))
    }
    ## Fall through to ME on error (non-PD S, etc.)
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
  ##   CSN loglik correction terms. For the PSKF at q<=3 (typical after pruning
  ##   with cut_tol=0.01), q>5 is rarely encountered in practice; q>5 arises
  ##   only in long filter runs on very high-dimensional models before pruning.
  ## Recommendation: for q > 5, the loglik is approximate with error ~0.01-0.05
  ## absolute per CDF call. Use cut_tol >= 0.01 to keep q small. If q > 5 is
  ## frequent, install mvtnorm and increase cut_tol to reduce q to <= 5.
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
## Returns list(Gamma = q'xp, nu = q', Delta = q'xq') with q' <= q.
#' @noRd
dim_red4_r <- function(Gamma, nu, Delta, Sigma, cut_tol = 0.01) {
  q <- nrow(Gamma)
  if (q == 0L) return(list(Gamma = Gamma, nu = nu, Delta = Delta))

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
  if (length(keep_idx) == q) {
    return(list(Gamma = Gamma, nu = nu, Delta = Delta))
  }
  if (length(keep_idx) == 0L) {
    return(list(
      Gamma = matrix(0, nrow = 0L, ncol = ncol(Gamma)),
      nu    = numeric(0),
      Delta = matrix(0, nrow = 0L, ncol = 0L)
    ))
  }
  list(
    Gamma = Gamma[keep_idx, , drop = FALSE],
    nu    = nu[keep_idx],
    Delta = Delta[keep_idx, keep_idx, drop = FALSE]
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

  ## Gamma_eta = Gamma_e Sigma_e RR' (RR Sigma_e RR')^{-1}   (n_exo x n_state)
  ## If Sigma_eta is nearly singular, use pseudoinverse (models with
  ## n_state > n_exo have redundant states that don't get shocked).
  tryCatch({
    S_eta_inv  <- solve(Sigma_eta)
    Gamma_eta  <- Gamma_e %*% Sigma_e %*% t(RR) %*% S_eta_inv
  }, error = function(e) {
    ## Pseudoinverse fallback for singular Sigma_eta (degenerate models)
    S_eta_pinv <- MASS::ginv(Sigma_eta)
    Gamma_eta  <<- Gamma_e %*% Sigma_e %*% t(RR) %*% S_eta_pinv
  })

  ## Delta_eta = I + Gamma_e Sigma_e_perp Gamma_e'
  ## Sigma_e_perp = Sigma_e - Sigma_e RR' S_eta^{-1} RR Sigma_e (Schur complement)
  tryCatch({
    Sigma_e_perp <- Sigma_e - Sigma_e %*% t(RR) %*% solve(Sigma_eta) %*% RR %*% Sigma_e
    Delta_eta    <- diag(n_exo) + Gamma_e %*% Sigma_e_perp %*% t(Gamma_e)
  }, error = function(e) {
    Delta_eta <<- diag(n_exo)  # fallback: square invertible case
  })

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
#' @noRd
.pskf_filter <- function(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                          Delta_eta, mu_eps, Sigma_eps,
                          cut_tol = 0.01,
                          store_path = FALSE) {
  ## Y: n_obs x T matrix
  n_obs   <- nrow(Y)
  n_T     <- ncol(Y)
  n_state <- nrow(TT)   # = ncol(TT)

  ## ---- Initialisation -------------------------------------------------------
  ## Start from the Gaussian stationary distribution (Lyapunov initialization).
  ## The skewness is zero at t=0: Gamma = [], nu = [], Delta = [].
  ## Solve P0 from P0 = TT P0 TT' + Sigma_eta via dlyap (iterative fallback).
  P0 <- .lyapunov_solve_r(TT, Sigma_eta)
  ## Fallback if Lyapunov solution failed
  if (is.null(P0) || any(!is.finite(P0)))
    P0 <- Sigma_eta

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
      pruned <- dim_red4_r(Gamma_pred, nu_pred, Delta_pred, Sigma_pred, cut_tol)
      Gamma_pred <- pruned$Gamma
      nu_pred    <- pruned$nu
      Delta_pred <- pruned$Delta
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
    Sigma_upd <- I_KZ %*% Sigma_pred   # Joseph form below is more stable but slower
    ## Joseph form for numerical stability (optional, use when needed):
    ## Sigma_upd <- I_KZ %*% Sigma_pred %*% t(I_KZ) +
    ##              K_gauss %*% Sigma_eps_t %*% t(K_gauss)

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
## Discrete Lyapunov equation solver (pure R, for Sigma_eta initialisation)
## ---------------------------------------------------------------------------
## Solve  P = A P A' + Q  by Schur decomposition (Bartels-Stewart).
## Falls back to fixed-point iteration for small matrices.
#' @noRd
.lyapunov_solve_r <- function(A, Q, max_iter = 500L, tol = 1e-10) {
  n <- nrow(A)
  ## Try a simple fixed-point iteration: P_{k+1} = A P_k A' + Q
  P <- Q
  for (i in seq_len(max_iter)) {
    P_new <- A %*% P %*% t(A) + Q
    if (max(abs(P_new - P)) < tol) return(P_new)
    P <- P_new
  }
  ## Did not converge (unit root / non-stationary): return last iterate
  P
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
#' @param ...         Ignored (for interface compatibility)
#' @return function(theta) -> list(logpost, loglik, logprior)
#' @noRd
make_log_posterior_pskf <- function(model, data, prior_spec, obs_vars,
                                     compiled, me_variance = 0,
                                     system_priors = NULL,
                                     cut_tol = 0.01,
                                     ...) {
  ## Validate data orientation: need T x n_obs
  if (ncol(data) == length(obs_vars)) {
    Y <- t(data)   # -> n_obs x T
  } else {
    Y <- data       # assume already n_obs x T
  }

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  sys_cache <- cache_system_structure(compiled)

  ## Warm-start cache for steady-state solve
  ss_warm <- NULL

  function(theta) {
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    params <- .apply_theta_to_params(model, theta)

    ## Also update estimated shock skewness from theta (v1 wiring;
    ## for v0 the shocks block provides fixed alpha)
    ## --- Steady state ---
    ss_result <- solve_steady_state(model, compiled, params,
                                    y0 = ss_warm, verbose = FALSE)
    if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
      if (!is.null(ss_warm))
        ss_result <- solve_steady_state(model, compiled, params,
                                        verbose = FALSE)
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
        ss_warm <<- NULL
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
      }
    }
    ss_warm <<- ss_result$ss

    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
    dr  <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Stationarity guard (mirrors gaussian path)
    ns <- length(dr$state_idx)
    ev <- dr$eigenvalues
    spectral_radius <- if (!is.null(ev) && length(ev) >= ns)
      max(Mod(ev[seq_len(ns)]))
    else
      max(Mod(eigen(dr$ghx[dr$state_idx, , drop = FALSE],
                    only.values = TRUE)$values))
    if (spectral_radius >= 1)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## --- Assemble state space ---
    state_idx <- dr$state_idx
    obs_idx   <- match(obs_vars, dr$endo_names)
    if (any(is.na(obs_idx)))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    TT <- dr$ghx[state_idx, , drop = FALSE]
    ZZ <- dr$ghx[obs_idx,   , drop = FALSE]
    ## Observation mean: d_obs (includes SS value and any DR offset)
    d_obs <- dr$ys[obs_vars]
    ## Demean Y
    Y_dm <- Y - d_obs   # n_obs x T (broadcasting over columns)

    ## --- CSN shock parameters ---
    exo_names <- dr$exo_names
    csn <- tryCatch(
      .get_csn_shock_params(model, exo_names, obs_vars, dr, params, me_variance),
      error = function(e) NULL
    )
    if (is.null(csn))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

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
        cut_tol   = cut_tol
      ),
      error = function(e) -Inf
    )
    if (!is.finite(ll))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## --- System priors ---
    lsp <- 0
    if (!is.null(system_priors) && length(system_priors) > 0L) {
      for (fn in system_priors) {
        v <- tryCatch(fn(dr, params), error = function(e) -Inf)
        lsp <- lsp + v
        if (!is.finite(lsp)) break
      }
    }

    ## power-posterior: temper the LIKELIHOOD only; lp (prior) and lsp (system
    ## prior) are prior-side and untempered.
    logpost <- lp + .dynhr_opt("power_posterior", default = 1) * ll + lsp
    list(logpost = logpost, loglik = ll, logprior = lp)
  }
}
