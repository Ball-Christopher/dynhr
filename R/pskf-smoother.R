## R/pskf-smoother.R
## --------------------------------------------------------------------------
## PSKF smoother: exact CSN posterior moments from the CSN forward pass.
##
## Reference: Guljanov, Mutschler & Trede (2026), "Pruned Skewed Kalman Filter
##   and Smoother with Application to DSGE Models," JEDC Vol. 187 (Dynare WP
##   #78). Reference implementation: github.com/gguljanov/pruned-skewed-kalman.
##
## IMPLEMENTATION (method = "csn"):
##   Step 1 (forward): .pskf_filter(store_path = TRUE) stores the per-period
##     filtered / predicted CSN parameters and, per period, which rows of the
##     pre-prune skew stack survived dim_red4_r (keep_path).
##   Step 2 (Gaussian part): standard RTS on (mu, Sigma) gives the moments
##     (s_G, P_G) of x_t | Y_{1:T} with the skew truncations dropped.
##   Step 3 (exact skew correction): every skew latent is a Gaussian variable
##     truncated at nu; the filter's time-T CSN says the latents U_T it
##     retains satisfy W = U_T - E[U_T | Y] >= nu_T with W ~ N(0, V_T).
##     Therefore, for every t,
##       E[x_t | Y]   = s_G,t + C_t g_T,   Var(x_t | Y) = P_G,t + C_t H_T C_t',
##     C_t = Cov(x_t, U_T | Y) (Gaussian part), (g_T, H_T) the gradient and
##     Hessian of log Phi_q(z; 0, V_T) at z = -nu_T. C_t is propagated by an
##     RTS pass on the augmented state (x_t, U_t) -- see pskf_smoother().
##   Exact given the filter's retained latents; when dim_red4_r prunes, the
##   smoother inherits the pruned filter's approximation (measured in
##   tests/testthat/test-fix-0926-pskf-exact-smoother.R).
##
##   History: the earlier backward pass shifted nu only (Gamma / Delta held
##   at their filtered values) -- first-order, ignoring latents born after t;
##   off an importance-sampled exact posterior by ~0.3 posterior sd at t < T
##   on a 2-state / 2-shock fixture. Its smoothed / filtered
##   covariances were the Gaussian-part covariances; they are now the exact
##   CSN posterior covariances.
##
##   When Gamma is all-zero (pure Gaussian), the correction vanishes -- exact
##   Gaussian RTS (method = "gaussian" also does this).
##
## METHOD = "gaussian" (legacy):
##   Gaussian RTS backward pass on the CSN forward-pass moments WITHOUT the
##   CSN skewness correction (means and covariances are the Gaussian part).
##   Exact at alpha=0; approximate at alpha != 0 (mean gap ~4.7e-3 for
##   alpha=2). Kept as a cheap baseline.
##
## Entry point: pskf_smoother()
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## .csn_logcdf_derivs / .csn_hazard_rate
##
## Analytic gradient (and optionally Hessian) of L(z) = log Phi_q(z; 0, V).
## With phi_i the N(0, V_ii) density at z_i and Phi_{q-1}(i) the conditional
## CDF of the other coordinates given coordinate i at z_i,
##   dL/dz_i = phi_i(z_i) Phi_{q-1}(z_{-i} - V_{-i,i} z_i / V_ii;
##                                   V_{-i,-i} - V_{-i,i} V_{i,-i} / V_ii) / Phi_q
## (exact; for q = 2 the conditional CDF is a univariate pnorm). Second
## derivatives of Phi_q: for i != j
##   d2 Phi / dz_i dz_j = phi_2(z_i, z_j; V_{ij,ij}) Phi_{q-2}(conditional on i, j)
## and, differentiating dPhi/dz_i through both its arguments,
##   d2 Phi / dz_i^2 = -(z_i / V_ii) dPhi/dz_i
##                     - sum_{k != i} (V_ki / V_ii) d2 Phi / dz_i dz_k,
## so Hess L = (Hess Phi) / Phi - g g'. Every Phi_m is logcdf_ME_r(), the
## filter's evaluator (exact pnorm / log-scale quadrature for m <= 2, the C++
## lattice evaluator mvn_logcdf_cpp for 3 <= m <= miwa_qmax; Miwa(128)
## / checked Miwa before -- Mendell-Elston beyond).
## miwa_qmax defaults to 7 here
## (the filter's likelihood path uses 5): these calls are made once per
## period by the smoother, the C++ evaluator costs ~3 ms at dim 6 and
## ~10 ms at dim 7 on PSKF calls (it raises the lattice size until the error
## estimate meets max(1e-5, 1e-7 |log p|); Miwa(128) was ~12 ms at dim 7),
## and
## Mendell-Elston is less accurate in the orthant tails that skewed
## posteriors live in (median |log-CDF error| 0.016, max 0.86 against
## mvtnorm after its 0.9.4.6 sign / double-shrink fix; before the fix a
## q = 6 unpruned fixture put ME-based smoothed means hundreds of MCSE off
## an exact posterior). The exact evaluators (~1e-6) are used where they
## are affordable.
##
## The former q >= 2 branch of .csn_hazard_rate was a central finite
## difference with an ABSOLUTE step 1e-5 (its comment claimed the step
## scaled with sqrt(D_ii)): at a latent scale sd(D_ii) ~ 1e-5 the step is one
## sd and the "gradient" is garbage.
##
## Uses: the truncated-normal moments of W ~ N(0, V) conditioned on W >= nu
## (the CSN latent representation) are E[W | .] = V g and
## Var(W | .) = V + V H V, with g, H the gradient / Hessian at z = -nu.
## Returns list(g = numeric(q), H = q x q matrix or NULL).
#' @noRd
.csn_logcdf_derivs <- function(z, V, hessian = FALSE, miwa_qmax = 7L) {
  q <- length(z)
  if (q == 0L)
    return(list(g = numeric(0), H = if (hessian) matrix(0, 0L, 0L)))
  V <- as.matrix(V)
  z <- as.numeric(z)

  if (q == 1L) {
    s  <- sqrt(max(V[1L, 1L], .Machine$double.eps))
    a  <- z / s
    r  <- exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))  # Mills: phi/Phi
    g  <- r / s
    H  <- if (hessian) matrix(-r * (a + r) / s^2, 1L, 1L)
    return(list(g = g, H = H))
  }

  L0 <- logcdf_ME_r(z, V, miwa_qmax = miwa_qmax)
  sd <- sqrt(diag(V))
  g  <- numeric(q)
  for (i in seq_len(q)) {
    mi <- z[-i] - V[-i, i] / V[i, i] * z[i]
    Si <- V[-i, -i, drop = FALSE] -
          V[-i, i, drop = FALSE] %*% V[i, -i, drop = FALSE] / V[i, i]
    Si <- 0.5 * (Si + t(Si))
    g[i] <- exp(dnorm(z[i], 0, sd[i], log = TRUE) + logcdf_ME_r(mi, Si, miwa_qmax = miwa_qmax) - L0)
  }
  if (!hessian) return(list(g = g, H = NULL))

  ## A = (Hess Phi) / Phi
  A <- matrix(0, q, q)
  for (i in seq_len(q - 1L)) {
    for (j in (i + 1L):q) {
      ij  <- c(i, j)
      Vij <- V[ij, ij]
      zij <- z[ij]
      det_ij <- Vij[1L, 1L] * Vij[2L, 2L] - Vij[1L, 2L]^2
      Vij_inv <- matrix(c(Vij[2L, 2L], -Vij[1L, 2L], -Vij[1L, 2L], Vij[1L, 1L]),
                        2L, 2L) / det_ij
      lphi2 <- -log(2 * pi) - 0.5 * log(det_ij) -
               0.5 * sum(zij * (Vij_inv %*% zij))
      if (q > 2L) {
        B  <- V[-ij, ij, drop = FALSE] %*% Vij_inv
        mc <- z[-ij] - as.numeric(B %*% zij)
        Sc <- V[-ij, -ij, drop = FALSE] - B %*% V[ij, -ij, drop = FALSE]
        Sc <- 0.5 * (Sc + t(Sc))
        lc <- logcdf_ME_r(mc, Sc, miwa_qmax = miwa_qmax)
      } else {
        lc <- 0
      }
      A[i, j] <- A[j, i] <- exp(lphi2 + lc - L0)
    }
  }
  for (i in seq_len(q)) {
    A[i, i] <- -(z[i] / V[i, i]) * g[i] - sum(V[-i, i] / V[i, i] * A[i, -i])
  }
  H <- A - tcrossprod(g)
  list(g = g, H = 0.5 * (H + t(H)))
}

## Hazard-rate vector h = d/dz log Phi_q(z; 0, D) at z = -nu (the CSN mean
## offset of the latent is D h); kept as the entry point used by the
## filtered/smoothed mean corrections.
#' @noRd
.csn_hazard_rate <- function(nu, D) {
  .csn_logcdf_derivs(-as.numeric(nu), D)$g
}


## ---------------------------------------------------------------------------
## .pskf_csn_derivs_checked
##
## (g, H) of log Phi_q(z; 0, D) at z = -nu (.csn_logcdf_derivs), refusing
## non-finite values. The earlier call sites wrapped the hazard in
## tryCatch(error = function(e) rep(0, q)), which would have silently dropped
## the whole skewness correction (returning the Gaussian-part mean) on any
## failure. A non-finite hazard is an explicit classed error instead.
#' @noRd
.pskf_csn_derivs_checked <- function(nu, D, t, what) {
  dv <- .csn_logcdf_derivs(-as.numeric(nu), D, hessian = TRUE)
  if (!all(is.finite(dv$g)) || !all(is.finite(dv$H))) {
    .dynhr_abort(
      "pskf_smoother(): non-finite CSN hazard rate for the ", what,
      " state at t = ", t, " (skew dimension q = ", length(nu), ").",
      class = "dynhr_pskf_smoother_error")
  }
  dv
}

## Exact mean / covariance offsets of x ~ CSN(mu, Sigma, Gamma, nu, Delta)
## relative to its Gaussian part N(mu, Sigma): with D = Delta + Gamma Sigma
## Gamma' and (g, H) the gradient / Hessian of log Phi_q(z; 0, D) at -nu,
##   E[x] - mu = Sigma Gamma' g,   Var(x) - Sigma = Sigma Gamma' H Gamma Sigma.
## NULL when q = 0 (no correction).
#' @noRd
.pskf_csn_moment_correction <- function(Sigma, Gamma, nu, Delta, t, what) {
  if (nrow(Gamma) == 0L) return(NULL)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  D  <- Delta + Gamma %*% Sigma %*% t(Gamma)
  D  <- 0.5 * (D + t(D))
  dv <- .pskf_csn_derivs_checked(nu, D, t, what)
  SG <- Sigma %*% t(Gamma)
  V  <- SG %*% dv$H %*% t(SG)
  list(mean = as.numeric(SG %*% dv$g), cov = 0.5 * (V + t(V)))
}


## ---------------------------------------------------------------------------
## pskf_smoother
##
## Run the PSKF forward pass (store_path=TRUE) and then apply the backward
## pass to produce smoothed state estimates.
##
## Two methods:
##   method = "csn"     -- exact CSN posterior means AND covariances (given
##                         the filter's retained skew latents; see the file
##                         header). Reduces to Gaussian RTS at alpha=0.
##   method = "gaussian" -- Gaussian RTS backward pass only; no CSN skewness
##                         correction. Legacy v1 method; preserved for backward
##                         compatibility and as a cheap baseline.
##
## Arguments:
##   Y         n_obs x T matrix of observations (NA = missing).
##   TT        n_state x n_state transition matrix.
##   ZZ        n_obs x n_state observation matrix.
##   mu_eta    length-n_state mean-correction for state noise (E[eta]=0 after).
##   Sigma_eta n_state x n_state state noise covariance.
##   Gamma_eta q_eta x n_state skewness loading matrix for state noise.
##   nu_eta    length-q_eta skewness location for state noise.
##   Delta_eta q_eta x q_eta skewness covariance for state noise.
##   mu_eps    length-n_obs measurement noise mean.
##   Sigma_eps n_obs x n_obs measurement noise covariance.
##   cut_tol   Pruning tolerance passed to .pskf_filter (default 0.01).
##   max_q     Skew-dimension rank cap passed to .pskf_filter (default 5).
##             Since the pruning mean-compensation fix, a cut
##             shifts the Gaussian location -- pass cut_tol = 0, max_q = Inf
##             when an exactly-Gaussian forward mean recursion is required
##             (e.g. the method = "gaussian" backward-compat oracle).
##   method    "csn" (default) or "gaussian". See above.
##
## Returns: list(
##   smoothed_means  T x n_state matrix of E[x_t | Y_{1:T}].
##   smoothed_covs   n_state x n_state x T array of Var(x_t | Y_{1:T})
##                   (method = "gaussian": the Gaussian-part covariances).
##   filtered_means  T x n_state matrix of E[x_t | Y_{1:t}].
##   filtered_covs   n_state x n_state x T array of Var(x_t | Y_{1:t})
##                   (method = "gaussian": the Gaussian-part covariances).
##   loglik          scalar log p(y_{1:T}) from the PSKF forward pass.
## )
##
## Note: at T=1 there is no backward pass; smoothed == filtered.
##
## Backward compatibility: all existing positional/named calls with 11 arguments
##   (Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, mu_eps,
##    Sigma_eps, cut_tol) work unchanged and return the same fields plus
##   the CSN correction in smoothed_means.  Callers that assume the Gaussian
##   baseline can pass method = "gaussian".
## @noRd
pskf_smoother <- function(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                           Delta_eta, mu_eps, Sigma_eps, cut_tol = 0.01,
                           method = c("csn", "gaussian"), max_q = 5L) {
  method <- match.arg(method)

  ## Ensure Y is a matrix (n_obs x T)
  if (!is.matrix(Y)) Y <- matrix(Y, nrow = 1L)
  n_obs   <- nrow(Y)
  n_T     <- ncol(Y)
  n_state <- nrow(TT)

  ## ---- Forward pass (with path storage) ------------------------------------
  fwd <- .pskf_filter(
    Y         = Y,
    TT        = TT,
    ZZ        = ZZ,
    mu_eta    = mu_eta,
    Sigma_eta = Sigma_eta,
    Gamma_eta = Gamma_eta,
    nu_eta    = nu_eta,
    Delta_eta = Delta_eta,
    mu_eps    = mu_eps,
    Sigma_eps = Sigma_eps,
    cut_tol   = cut_tol,
    max_q     = max_q,
    ## pruning compensation with the deterministic Phi evaluator up to dim 7
    ## (the pre-prune stack at the default max_q = 5 with two skew shocks;
    ## dims 3-7: the C++ lattice evaluator, see logcdf_ME_r). The backward
    ## pass carries the compensation (lambda_path) to earlier periods, where
    ## the former likelihood-path ME evaluation (dim > 2; before its 0.9.4.6
    ## sign fix) was measured up to 2.4 posterior sd off (see
    ## .csn_mean_offset). Cost: ~3 ms at dim 6,
    ## ~10 ms at dim 7 on PSKF calls -- hence the cap. When a connected Phi
    ## block of dimension 6-7 occurs, the returned loglik can differ slightly
    ## from .pskf_filter()'s default (offset_miwa_qmax = 5: ME there, whose
    ## log-CDF error is ~0.016 median).
    offset_miwa_qmax = 7L,
    store_path = TRUE
  )

  ## An infeasible forward pass (non-stationary TT: no Lyapunov initial
  ## covariance; or data with zero density under a degenerate innovation
  ## covariance, e.g. a noise-free observable of a zero-variance state that
  ## does not equal its prediction) returns ll = -Inf with a missing or
  ## truncated path. There is no filtered distribution to smooth, so say so
  ## rather than failing on a NULL path.
  if (!is.finite(fwd$ll)) {
    .dynhr_abort(
      "pskf_smoother(): the PSKF forward pass is infeasible (loglik = -Inf): ",
      "either TT has no stationary initial covariance, or the data have zero ",
      "density under a singular innovation covariance ZZ Sigma_pred ZZ' + ",
      "Sigma_eps (a linearly dependent / exactly predictable noise-free ",
      "observable that the data contradict).",
      class = "dynhr_pskf_smoother_error")
  }

  loglik          <- fwd$ll
  mu_pred_path    <- fwd$mu_pred_path
  Sigma_pred_path <- fwd$Sigma_pred_path
  mu_filt_path    <- fwd$mu_filt_path
  Sigma_filt_path <- fwd$Sigma_filt_path
  Gamma_filt_path <- fwd$Gamma_filt_path
  nu_filt_path    <- fwd$nu_filt_path
  Delta_filt_path <- fwd$Delta_filt_path
  Gamma_pred_path <- fwd$Gamma_pred_path
  Delta_pred_path <- fwd$Delta_pred_path
  keep_path       <- fwd$keep_path
  lambda_path     <- fwd$lambda_path

  ## ---- Gaussian RTS backward pass ------------------------------------------
  ## Gaussian part of the posterior: moments of x_t | Y_{1:T} with the skew
  ## truncations dropped.
  s_smooth <- matrix(0, n_T, n_state)   # smoothed means (T x n_state)
  P_smooth <- array(0, dim = c(n_state, n_state, n_T))
  s_smooth[n_T, ]   <- mu_filt_path[[n_T]]
  P_smooth[, , n_T] <- Sigma_filt_path[[n_T]]

  ## Backward sweep: t = T-1 down to 1 (empty when n_T == 1)
  for (step in seq_len(n_T - 1L)) {
    t <- n_T - step   # reverse order: T-1, T-2, ..., 1

    mu_f    <- mu_filt_path[[t]]
    Sigma_f <- Sigma_filt_path[[t]]
    mu_p    <- mu_pred_path[[t + 1L]]
    Sigma_p <- Sigma_pred_path[[t + 1L]]

    ## RTS gain: G_t = Cov(x_t, x_{t+1} | Y_{1:t}) Var(x_{t+1} | Y_{1:t})^+
    ##                 = Sigma_f TT' Sigma_p^+.
    ## Sigma_p = Sigma_{t+1|t} is routinely SINGULAR here, for the same reasons
    ## as in the filter's prediction step (see .pskf_filter): the
    ## contemporaneous order-1 state carries observables that are exact linear
    ## functions of states and shocks, noise-free observation makes state
    ## directions known exactly (an AR(2) observed without error knows its
    ## lagged state), and the order-2 state carries Kronecker duplicates. For
    ## the genuine joint Gaussian (x_t, x_{t+1}) | Y_{1:t} the cross-covariance
    ## Sigma_f TT' annihilates null(Sigma_p) (a null direction v has
    ## Var(v' x_{t+1}) = 0, hence Sigma_f TT' v = Cov(x_t, v' x_{t+1}) = 0),
    ## i.e. its rows lie in range(Sigma_p), so the Moore-Penrose
    ## pseudoinverse gives the EXACT conditional mean and covariance:
    ## G_t Sigma_p G_t' = Sigma_f TT' Sigma_p^+ TT Sigma_f, the exact Schur
    ## complement.
    ## The former `solve(Sigma_p + 1e-10 I)` was not exact: it is an
    ## ABSOLUTE jitter, so every eigen-direction lambda of Sigma_p got the
    ## relative gain error 1e-10 / lambda -- 1e-6 at DSGE scale (shock
    ## stderr 0.01, lambda ~ 1e-4), O(1) at stderr 1e-5 -- and its SVD
    ## fallback (reached through tryCatch only when solve() errored) used yet
    ## another rank rule. Same scale-free rank rule as the filter
    ## (.csn_sym_pinv).
    Sigma_p_sym <- 0.5 * (Sigma_p + t(Sigma_p))   # enforce symmetry
    G_t <- Sigma_f %*% t(TT) %*% .csn_sym_pinv(Sigma_p_sym)

    ## Gaussian smoother update (mean and covariance)
    s_smooth[t, ] <- mu_f + as.numeric(G_t %*% (s_smooth[t + 1L, ] - mu_p))
    P_s_next      <- P_smooth[, , t + 1L]
    diff_P        <- P_s_next - Sigma_p_sym
    P_smooth[, , t] <- Sigma_f + G_t %*% diff_P %*% t(G_t)

    ## Enforce symmetry of smoothed covariance
    P_smooth[, , t] <- 0.5 * (P_smooth[, , t] + t(P_smooth[, , t]))
  }

  ## Filtered moments (Gaussian part; the CSN correction is added below).
  s_filt <- matrix(0, n_T, n_state)
  P_filt <- array(0, dim = c(n_state, n_state, n_T))
  for (t in seq_len(n_T)) {
    s_filt[t, ]   <- mu_filt_path[[t]]
    P_filt[, , t] <- Sigma_filt_path[[t]]
  }

  if (method == "csn") {
    ## ---- Exact CSN correction ----------------------------------------------
    ## Latent representation: every skew-normal shock row carries a Gaussian
    ## latent u = Gamma_eta (eta - mu_eta) + e, e ~ N(0, Delta_eta), and the
    ## CSN law is the Gaussian law conditioned on u >= nu_eta. The filter's
    ## period-T law x_T | Y ~ CSN(mu_T, Sigma_T, Gamma_T, nu_T, Delta_T) is
    ## exactly this statement for the latents U_T it retains: in the Gaussian
    ## part, W = U_T - E[U_T | Y] ~ N(0, V_T), V_T = Delta_T + Gamma_T
    ## Sigma_T Gamma_T', and the truncation is W >= nu_T. Writing
    ## x_t = s_G,t + B W + (independent of W) with B V_T = C_t =
    ## Cov(x_t, U_T | Y), and using E[W | W >= nu_T] = V_T g_T and
    ## Var(W | W >= nu_T) = V_T + V_T H_T V_T ((g_T, H_T) = gradient / Hessian
    ## of log Phi_q(z; 0, V_T) at z = -nu_T, .csn_logcdf_derivs), for EVERY t
    ##   E[x_t | Y]   = s_G,t + C_t g_T,
    ##   Var(x_t | Y) = P_G,t + C_t H_T C_t'.
    ## At t = T, C_T = Sigma_T Gamma_T' (the filtered CSN moments).
    ##
    ## C_t: RTS on the AUGMENTED Gaussian state z_t = (x_t, U_t), U_t the
    ## latents retained at t. z is Markov (z_{t+1} = F z_t + independent
    ## noise; F keeps the surviving old latents, keep_path), so
    ##   Cov(z_t, U_T | Y) = J_t Cov(z_{t+1}, U_T | Y),
    ##   J_t = Cov(z_t, z_{t+1} | Y_{1:t}) Var(z_{t+1} | Y_{1:t})^+,
    ## with every block available from the filter's paths:
    ##   Var(z_t | Y_{1:t}) = [Sigma_t, Sigma_t Gamma_t'; Gamma_t Sigma_t, V_t],
    ##   Cov(z_t, x_{t+1})  = Var(z_t)[, x] TT',
    ##   Cov(z_t, u)        = Var(z_t)[, k] for a surviving old latent k,
    ##                        0 for a latent born at t + 1.
    ## The pseudoinverse is exact for the same reason as the RTS gain above
    ## (cross-covariance rows lie in the range of the joint covariance).
    for (t in seq_len(n_T)) {
      corr <- .pskf_csn_moment_correction(
        P_filt[, , t], Gamma_filt_path[[t]], nu_filt_path[[t]],
        Delta_filt_path[[t]], t, "filtered")
      if (is.null(corr)) next
      s_filt[t, ]   <- s_filt[t, ] + corr$mean
      P_filt[, , t] <- P_filt[, , t] + corr$cov
    }

    ## PRUNING. When dim_red4_r cuts latents at t + 1, the filter replaces
    ## their truncation by a first-moment shift of x_{t+1}:
    ## mu_shift = Cov(x_{t+1}, W) lambda_{t+1} (W the pre-prune stack
    ## [U_t; u_{t+1}], lambda its latent-space cut compensation). The model
    ## this approximation describes is "W's mean is shifted by
    ## Var(W) lambda", which also shifts every z_t correlated with W. The
    ## Gaussian RTS pass above treats mu_shift as an exogenous input to
    ## x_{t+1} and so reproduces that model for every period >= t + 1 (the
    ## data after t see W only through x_{t+1}), but not before: there the
    ## missing term is Cov(z_t, W | Y_{1:t}) lambda (only the OLD latents of
    ## W, the columns of Var(z_t | Y_{1:t}) at U_t, are correlated with z_t),
    ## carried further back by the same augmented gains J_t as the U_T term
    ## (Cov(z_t, z_s | Y_{1:s}) = J_t ... J_{s-1} Var(z_s | Y_{1:s})). So the
    ## total mean correction of z_t is the vector
    ##   v_t = J_t v_{t+1} + Var(z_t | Y_{1:t})[, U_t] lambda_{t+1}[old],
    ##   v_T = Cov(z_T, U_T | Y) g_T.
    ## With no cut this is exactly M_t g_T. The cut term is exact for the
    ## pruned filter's own approximating model (not for the true posterior);
    ## its variance effect is not added (smoothed covariances then carry the
    ## U_T term only).
    q_T   <- nrow(Gamma_filt_path[[n_T]])
    Sig_T <- Sigma_filt_path[[n_T]]
    Gam_T <- Gamma_filt_path[[n_T]]
    V_T   <- Delta_filt_path[[n_T]] + Gam_T %*% Sig_T %*% t(Gam_T)
    V_T   <- 0.5 * (V_T + t(V_T))
    ## M = Cov(z_t, U_T | Y), rows z_t = (x_t, U_t); v = mean correction of z_t
    M <- rbind(Sig_T %*% t(Gam_T), V_T)
    if (q_T > 0L) {
      dv <- .pskf_csn_derivs_checked(nu_filt_path[[n_T]], V_T, n_T,
                                     "smoothed")
      v  <- as.numeric(M %*% dv$g)
    } else {
      v  <- rep(0, n_state)
    }
    for (t in n_T:1L) {
      if (t < n_T) {
        Sig_t <- Sigma_filt_path[[t]]
        Gam_t <- Gamma_filt_path[[t]]
        q_t   <- nrow(Gam_t)
        V_t   <- Delta_filt_path[[t]] + Gam_t %*% Sig_t %*% t(Gam_t)
        Pz_t  <- rbind(cbind(Sig_t, Sig_t %*% t(Gam_t)),
                       cbind(Gam_t %*% Sig_t, V_t))
        Sig_p <- Sigma_pred_path[[t + 1L]]
        Gam_p <- Gamma_pred_path[[t + 1L]]
        V_p   <- Delta_pred_path[[t + 1L]] + Gam_p %*% Sig_p %*% t(Gam_p)
        Pz_p  <- rbind(cbind(Sig_p, Sig_p %*% t(Gam_p)),
                       cbind(Gam_p %*% Sig_p, V_p))
        Pz_p  <- 0.5 * (Pz_p + t(Pz_p))
        keep  <- keep_path[[t + 1L]]
        old   <- keep <= q_t
        C_u   <- matrix(0, n_state + q_t, length(keep))
        C_u[, old] <- Pz_t[, n_state + keep[old], drop = FALSE]
        C_z   <- cbind(Pz_t[, seq_len(n_state), drop = FALSE] %*% t(TT), C_u)
        J_t   <- C_z %*% .csn_sym_pinv(Pz_p)
        M     <- J_t %*% M
        v     <- as.numeric(J_t %*% v)
        lam   <- lambda_path[[t + 1L]]
        if (q_t > 0L && length(lam) > 0L && any(lam[seq_len(q_t)] != 0)) {
          v <- v + as.numeric(Pz_t[, n_state + seq_len(q_t), drop = FALSE] %*%
                              lam[seq_len(q_t)])
        }
      }
      s_smooth[t, ] <- s_smooth[t, ] + v[seq_len(n_state)]
      if (q_T > 0L) {
        C_t  <- M[seq_len(n_state), , drop = FALSE]
        P_st <- P_smooth[, , t] + C_t %*% dv$H %*% t(C_t)
        P_smooth[, , t] <- 0.5 * (P_st + t(P_st))
      }
    }
  }

  list(
    smoothed_means  = s_smooth,
    smoothed_covs   = P_smooth,
    filtered_means  = s_filt,
    filtered_covs   = P_filt,
    loglik          = loglik
  )
}
