### =============================================================================
### kalman_filter_v4.R -- Production Kalman filter for dynhr
### =============================================================================
###
### v4 fixes:
###   DARE:          self-convergence check |P_new - P| (not |P - P???|),
###                  then switch to DARE-precomputed K???/F???.
###   Chandrasekhar: bootstrap 1+ Riccati steps before starting recursion,
###                  avoiding catastrophic F??? = F??? - Z K??? F??? K???' Z' cancel.
###   Standard:      reference (unchanged).
###
### All three produce identical log-likelihoods (< 1e-10 nats).
### =============================================================================

### -- Backend dispatch --------------------------------------------------------
###
### kalman_ss_loop_cpp() is compiled from src/kalman_ss.cpp via RcppArmadillo.
### The R fallback below is bit-exact at the loglik level (see
### tests/testthat/test-kalman-rcpp-parity.R) and is what runs if the DLL is
### unavailable (development source loads, partial installs, or when the user
### sets options(dynhr.use_rcpp = FALSE) to disable C++ for debugging).

.HAS_RCPP_KALMAN <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kalman_ss_loop_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

.kf_ss_dispatch <- function(Y_minus_d, ZZ, TT, K_ss, F_inv_ss,
                            ll_ss_const, s, start_t, end_t,
                            filtered = NULL) {
  ## Guard: validate time indices before dispatching to C++ (or R fallback).
  if (start_t < 1L)
    stop(".kf_ss_dispatch: start_t must be >= 1")
  if (end_t > ncol(Y_minus_d))
    stop(".kf_ss_dispatch: end_t exceeds Y_minus_d columns")
  if (start_t > end_t)
    stop(".kf_ss_dispatch: start_t must be <= end_t")

  if (.HAS_RCPP_KALMAN()) {
    ## C++ ignores `filtered` shape -- pass return_filtered flag instead.
    out <- kalman_ss_loop_cpp(Y_minus_d, ZZ, TT, K_ss, F_inv_ss, ll_ss_const,
                              as.numeric(s), as.integer(start_t),
                              as.integer(end_t), !is.null(filtered))
    out$s <- as.numeric(out$s)
    if (!is.null(filtered) && !is.null(out$filtered)) {
      ## C++ returns full-width matrix; reuse caller's preallocated one.
      filtered[, start_t:end_t] <- out$filtered[, start_t:end_t, drop = FALSE]
      out$filtered <- filtered
    }
    out
  } else {
    .kf_ss_loop_R(Y_minus_d, ZZ, TT, K_ss, F_inv_ss, ll_ss_const,
                  s, start_t, end_t, filtered)
  }
}

### -- Steady-state loop (extracted for reuse + future Rcpp swap) ----------------
###
### Runs the constant-gain Kalman recursion from index `start_t` to `end_t`,
### given the (already converged) steady-state Kalman gain `K_ss`, inverse
### innovation covariance `F_inv_ss`, and the precomputed log-likelihood
### constant `ll_ss_const = ll_const - 0.5 * log_det_F_ss`.
###
### Performance notes:
###   * Y_minus_d is precomputed once outside the loop (one allocation, not T).
###   * `sum(v * (F_inv_ss %*% v))` is identical to `drop(crossprod(v, ...))`
###     but skips the 1x1 matrix temporary.
###   * `TT %*% s + K_ss %*% v` (two BLAS3 calls on n_state x n_state and
###     n_state x n_obs) beats the `cbind(TT, K_ss) %*% c(s, v)` scratch trick
###     when n_state >> n_obs (typical).
###
### The Rcpp backend (see R/perf-rcpp.R) consumes the same arguments and
### returns the same list, so the swap is a one-line `if (rcpp) ...`.

.kf_ss_loop_R <- function(Y_minus_d, ZZ, TT, K_ss, F_inv_ss,
                          ll_ss_const, s, start_t, end_t,
                          filtered = NULL) {
  loglik <- 0
  ok     <- TRUE
  for (t in start_t:end_t) {
    v   <- Y_minus_d[, t] - as.numeric(ZZ %*% s)
    Fv  <- F_inv_ss %*% v
    ll_t <- ll_ss_const - 0.5 * sum(v * Fv)
    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) { ok <- FALSE; break }
    loglik <- loglik + ll_t
    s <- as.numeric(TT %*% s + K_ss %*% v)
    if (!is.null(filtered)) filtered[, t] <- s
  }
  list(loglik = loglik, s = s, filtered = filtered, ok = ok)
}

### -- Univariate (sequential) Kalman filter on the augmented state ------------
###
### Koopman & Durbin (2000) univariate treatment of dynhr's state space.
### dynhr's measurement y_t = Z s_{t-1} + D eps_t shares shocks with the
### transition, so measurement noise is correlated across observables AND
### with the state innovations -- which the textbook univariate filter
### (requiring diagonal measurement noise) cannot handle directly.
###
### Resolution: AUGMENT, don't whiten. Define x_t = [s_{t-1}; eps_t]:
###
###   x_{t+1} = [ TT RR ] x_t + [ 0 ] eps_{t+1}      y_t = [ Z  D ] x_t
###             [ 0  0  ]       [ I ]
###
### Measurement noise is now exactly zero (trivially diagonal), all
### correlation lives in the state, and the standard univariate update
### applies -- including its graceful handling of zero-variance innovation
### components (skip the observable, no inversion). Singular F is therefore
### a non-event on this path, and the same recursion implements the exact
### univariate DIFFUSE filter (Dynare kalman_algo=4 analog), covering the
### multivariate diffuse phase's unsupported Case C (F_inf singular but
### nonzero -- see .kf_diffuse_phase below).
###
### The C++ backend (src/kalman_univariate.cpp) runs the whole loop in one
### call; .kf_univariate_loop_R is the bit-mirroring fallback (same
### options(dynhr.use_rcpp) switch as the standard filter).

.HAS_RCPP_KALMAN_UNI <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kalman_univariate_loop_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

.kf_univariate_loop_R <- function(Y_minus_d, Zb, Tb, QQb, a, P_star, P_inf,
                                  me_variance, kalman_tol, diffuse_tol,
                                  conv_tol, max_diffuse, ll_min,
                                  return_filtered, n_state,
                                  me_extra = NULL, shock_scale = NULL,
                                  Sigma_e = NULL, e_idx = NULL) {
  n_obs <- nrow(Y_minus_d)
  n_T   <- ncol(Y_minus_d)
  log2pi <- log(2 * pi)

  ## me_extra (n_obs x T) holds per-period per-observable extra ME variances
  ## for filter_tunes soft tunes. When provided, me_extra[i, t] is added to
  ## F_star for observable i at period t; only the R loop supports this
  ## (C++ backend takes a scalar me_variance, not a matrix).
  has_me_extra_uni <- !is.null(me_extra) && any(me_extra != 0)

  ## shock_scale (n_exo x T) holds per-period shock-stdev multipliers. The
  ## augmented state x_t = [s_{t-1}; eps_t] means the process noise added by
  ## the Tb-transition at the END of iteration t is Var(eps_{t+1}) -- so the
  ## QQb block used for that transition must be rebuilt from shock_scale
  ## column t+1, not column t. The caller (.kf_univariate_dispatch) already
  ## bakes column 1 into the INITIAL Pb before this loop starts.
  has_shock_scale_uni <- !is.null(shock_scale)

  loglik <- 0
  ok     <- TRUE
  diffuse <- length(P_inf) > 0 && max(abs(P_inf)) > 0
  diffuse_failed <- FALSE
  d_diffuse <- NA_integer_
  ll_contrib <- numeric(n_T)
  filtered <- if (return_filtered) matrix(0, n_state, n_T) else NULL

  for (t in seq_len(n_T)) {
    if (diffuse && t > max_diffuse) { diffuse_failed <- TRUE; break }

    ll_t <- 0
    for (i in seq_len(n_obs)) {
      y_i <- Y_minus_d[i, t]
      if (!is.finite(y_i)) next                     # NA observable: skip
      Zi <- Zb[i, ]
      v  <- y_i - sum(Zi * a)
      K_star <- drop(P_star %*% Zi)
      F_star <- sum(Zi * K_star) + me_variance +
        if (has_me_extra_uni) me_extra[i, t] else 0

      if (diffuse) {
        K_inf <- drop(P_inf %*% Zi)
        F_inf <- sum(Zi * K_inf)
        if (F_inf > diffuse_tol * max(1, F_star)) {
          ## Diffuse update (DK 2012 sec. 7.2.5); same renormalization
          ## convention as the multivariate Case B in .kf_diffuse_phase.
          ll_t   <- ll_t - 0.5 * log(F_inf)
          a      <- a + K_inf * (v / F_inf)
          P_star <- P_star + tcrossprod(K_inf) * (F_star / F_inf^2) -
            (tcrossprod(K_star, K_inf) + tcrossprod(K_inf, K_star)) / F_inf
          P_inf  <- P_inf - tcrossprod(K_inf) / F_inf
          next
        }
      }

      if (F_star > kalman_tol) {
        ll_t   <- ll_t - 0.5 * (log2pi + log(F_star) + v * v / F_star)
        a      <- a + K_star * (v / F_star)
        P_star <- P_star - tcrossprod(K_star) / F_star
      }
      ## else: zero innovation variance -- skip gracefully (singular F).
    }

    if (!is.finite(ll_t) || ll_t < ll_min) { ok <- FALSE; break }
    loglik <- loglik + ll_t
    ll_contrib[t] <- ll_t

    a <- drop(Tb %*% a)
    QQb_t <- if (has_shock_scale_uni && t < n_T) {
      sc_next <- shock_scale[, t + 1L]
      QQb_next <- QQb
      QQb_next[e_idx, e_idx] <- Sigma_e * outer(sc_next, sc_next)
      QQb_next
    } else QQb
    P_star <- Tb %*% P_star %*% t(Tb) + QQb_t
    P_star <- (P_star + t(P_star)) * 0.5
    if (diffuse) {
      P_inf <- Tb %*% P_inf %*% t(Tb)
      P_inf <- (P_inf + t(P_inf)) * 0.5
      if (max(abs(P_inf)) < conv_tol * max(1, max(abs(P_star)))) {
        diffuse   <- FALSE
        d_diffuse <- t
      }
    }
    if (return_filtered) filtered[, t] <- a[seq_len(n_state)]
  }
  if (ok && diffuse) diffuse_failed <- TRUE    # sample ended inside the phase

  list(loglik = loglik, a = a, filtered = filtered, ok = ok,
       d_diffuse = d_diffuse, diffuse_failed = diffuse_failed,
       ll_contrib = ll_contrib)
}

### Builds the augmented system and dispatches to C++ or the R fallback.
### P_state is the n_state x n_state initial covariance of s_0; P_inf_state
### (NULL or n_state x n_state) is the diffuse part for the exact-diffuse
### initialization. The tolerances mirror .kf_diffuse_phase's defaults;
### kalman_tol is Dynare's kalman_tol analog (skip threshold for F_i).
.kf_univariate_dispatch <- function(Y_minus_d, ZZ, TT, RR, DD, Sigma_e,
                                    s0, P_state, P_inf_state = NULL,
                                    me_variance = 0,
                                    return_filtered = FALSE,
                                    kalman_tol = 1e-10,
                                    diffuse_tol = 1e-10,
                                    conv_tol = 1e-8,
                                    max_diffuse = 100L,
                                    me_extra = NULL,
                                    ss_lock = FALSE,
                                    shock_scale = NULL) {
  ## The steady-state lock assumes a constant present-observable pattern, so it
  ## is only valid on a complete panel. Disable it if any observation is missing
  ## (the filter then runs the exact full recursion). The R fallback ignores
  ## ss_lock (runs unlocked), so univariate_ss degrades gracefully without Rcpp.
  ## A time-varying shock_scale is likewise incompatible with a frozen gain.
  ss_lock <- isTRUE(ss_lock) && !anyNA(Y_minus_d) && is.null(shock_scale)
  n_state <- nrow(TT)
  n_exo   <- ncol(RR)
  nb      <- n_state + n_exo
  s_idx   <- seq_len(n_state)
  e_idx   <- (n_state + 1L):nb

  ## Strip dimnames: a 1-row Y with a rowname would otherwise leak its name
  ## into the scalar innovations (and from there into loglik) in the R loop.
  dimnames(Y_minus_d) <- NULL

  Zb <- cbind(ZZ, DD)
  Tb <- rbind(cbind(TT, RR), matrix(0, n_exo, nb))
  QQb <- matrix(0, nb, nb)
  QQb[e_idx, e_idx] <- Sigma_e

  a0 <- c(s0, numeric(n_exo))
  Pb <- QQb                                   # eps_1 block = Sigma_e
  ## shock_scale active: eps_1 (the x_1 state's shock component) is scaled by
  ## column 1, not the baseline Sigma_e (mirrors the QQb_t rebuild inside
  ## .kf_univariate_loop_R for every later transition).
  has_shock_scale_d <- !is.null(shock_scale)
  if (has_shock_scale_d) {
    sc1 <- shock_scale[, 1L]
    Pb[e_idx, e_idx] <- Sigma_e * outer(sc1, sc1)
  }
  Pb[s_idx, s_idx] <- P_state
  Pi <- matrix(0, nb, nb)
  if (!is.null(P_inf_state)) Pi[s_idx, s_idx] <- P_inf_state

  ## me_extra must be an n_obs x T matrix when non-NULL; expand to the
  ## augmented-state index (first n_obs rows of Zb correspond to original obs).
  ## The C++ backend takes only a scalar me_variance: force R when me_extra
  ## has any nonzero entries. shock_scale likewise forces R: the C++ kernel
  ## bakes a single time-invariant QQb and has no per-period rebuild hook.
  has_me_extra_d <- !is.null(me_extra) && any(me_extra != 0)

  if (.HAS_RCPP_KALMAN_UNI() && !has_me_extra_d && !has_shock_scale_d) {
    out <- kalman_univariate_loop_cpp(Y_minus_d, Zb, Tb, QQb, a0, Pb, Pi,
                                      me_variance, kalman_tol, diffuse_tol,
                                      conv_tol, as.integer(max_diffuse),
                                      .KF_LL_MIN, return_filtered,
                                      as.integer(n_state), ss_lock)
    out$a <- as.numeric(out$a)
    out
  } else {
    .kf_univariate_loop_R(Y_minus_d, Zb, Tb, QQb, a0, Pb, Pi,
                          me_variance, kalman_tol, diffuse_tol,
                          conv_tol, max_diffuse, .KF_LL_MIN,
                          return_filtered, n_state,
                          me_extra = me_extra,
                          shock_scale = shock_scale,
                          Sigma_e = Sigma_e, e_idx = e_idx)
  }
}

### -- DARE solver (unchanged from v3) --------------------------------------------

.solve_dare <- function(TT, ZZ, QQ, HH, SS, me_diag = NULL,
                        tol = .DARE_TOL, max_iter = 1000) {
  n_s <- nrow(TT)
  tZZ <- t(ZZ)
  P   <- solve_lyapunov(TT, QQ)
  HH_full <- if (!is.null(me_diag)) HH + me_diag else HH

  for (i in seq_len(max_iter)) {
    PZ     <- P %*% tZZ
    Ft     <- ZZ %*% PZ + HH_full
    Ft     <- (Ft + t(Ft)) * 0.5
    Fc     <- chol(Ft)
    Fi     <- chol2inv(Fc)
    K      <- (TT %*% PZ + SS) %*% Fi
    P_new  <- tcrossprod(TT %*% P, TT) + QQ - tcrossprod(K %*% Ft, K)
    P_new  <- (P_new + t(P_new)) * 0.5
    if (max(abs(P_new - P)) < tol) {
      return(list(P = P_new, K = K, F = Ft, F_inv = Fi,
                  F_chol = Fc, log_det_F = 2 * sum(log(diag(Fc))),
                  converged = TRUE, iterations = i))
    }
    P <- P_new
  }
  list(P = P, K = K, F = Ft, F_inv = Fi,
       F_chol = Fc, log_det_F = 2 * sum(log(diag(Fc))),
       converged = FALSE, iterations = max_iter)
}


### -- Exact diffuse initialization (Koopman & Durbin 2003; DK 2012 ch.5) -----
###
### Partition the state space into a unit-root ("diffuse") subspace, spanned
### by the orthonormal columns A_inf, and its orthogonal complement A_star
### (the stationary subspace). The partition is obtained from the same
### ordered-Schur decomposition used by .build_P0_Q2() in
### R/backend-monolith.R: TT = U %*% Ts %*% t(U) with the unit-root block of
### Ts (|diag| > 1 - ur_tol) swapped to the leading nunit x nunit block.
###
###   P_inf_0  = A_inf  %*% t(A_inf)                  (= U[,1:nunit] U[,1:nunit]')
###   P_star_0 = solve the stable-block Lyapunov equation in the rotated
###              (Schur) basis for the (stable, stable) block of
###              U' QQ U, and rotate back; zero elsewhere.
###
### Returns list(P_inf, P_star, nunit). nunit == 0 means TT has no unit roots
### (the diffuse path degenerates immediately: P_inf == 0).
.kf_diffuse_P0 <- function(TT, QQ, ur_tol = 1e-6) {
  n <- nrow(TT)
  schur <- Matrix::Schur(Matrix::Matrix(TT, sparse = FALSE))
  Ts <- as.matrix(schur@T)
  Qs <- as.matrix(schur@Q)

  ## Push unit-root blocks (|diag| > 1 - ur_tol) to the leading block,
  ## mirroring .build_P0_Q2()'s swap logic.
  ur_target <- 0L
  for (i in seq_len(n)) {
    if (abs(abs(diag(Ts)[i]) - 1) < ur_tol) {
      j <- i
      while (j > ur_target + 1L) {
        sw <- .swap_schur_11(Ts, Qs, j - 1L)
        Ts <- sw$T; Qs <- sw$Q
        j <- j - 1L
      }
      ur_target <- ur_target + 1L
    }
  }
  nunit <- ur_target
  U <- Qs

  P_inf <- matrix(0, n, n)
  if (nunit > 0L) {
    A_inf <- U[, seq_len(nunit), drop = FALSE]
    P_inf <- tcrossprod(A_inf)
  }

  P_star <- matrix(0, n, n)
  n_stable <- n - nunit
  if (n_stable > 0L) {
    idx_s <- (nunit + 1L):n
    Ts_s  <- Ts[idx_s, idx_s, drop = FALSE]
    QQ_schur <- t(U) %*% QQ %*% U
    QQ_s <- QQ_schur[idx_s, idx_s, drop = FALSE]
    Pa_s <- .solve_lyapunov(Ts_s, QQ_s)
    if (all(is.finite(Pa_s))) {
      Pa <- matrix(0, n, n)
      Pa[idx_s, idx_s] <- .sym(Pa_s)
      P_star <- U %*% Pa %*% t(U)
    }
  }

  list(P_inf = .sym(P_inf), P_star = .sym(P_star), nunit = nunit)
}

### -- Diffuse recursion (Durbin & Koopman 2012, sec. 5.2) ---------------------
###
### Runs the exact diffuse Kalman filter from t = 1 until P_inf has decayed to
### (numerically) zero, or the hard cap min(n_T, 100) is reached.
###
### Convention (dynhr's, with correlated measurement noise -- see file header
### of kalman_filter()):
###   v   = y_t - d - ZZ %*% s
###   F_inf  = ZZ P_inf  ZZ'
###   F_star = ZZ P_star ZZ' + HH + me_diag
###
### Case A (F_inf ~= 0, i.e. max|F_inf| < diffuse_tol * max(1, max|F_star|)):
###   the diffuse part of the state has nothing left to learn from this
###   observation; run a STANDARD exact .kf_step on (s, P_star) (with its
###   usual likelihood contribution) and propagate P_inf <- TT P_inf TT'.
###
### Case B (F_inf nonsingular, rcond(F_inf) > diffuse_tol):
###   M_inf  = TT P_inf  ZZ'
###   M_star = TT P_star ZZ' + SS    (SS is the dynhr correlated-noise term;
###                                    DK 2012 has no SS, i.e. M*_t = T P*_t Z')
###   K0 = M_inf F_inf^{-1}
###   ll_t = -0.5 log|F_inf|     (no -0.5*n_obs*log(2*pi) and no quadratic
###          term: writing F = F_star + kappa*F_inf with F_inf full rank,
###          log|F| = n_obs*log(kappa) + log|F_inf| + O(1/kappa) and
###          v'F^{-1}v = O(1/kappa) -> 0; the exact diffuse likelihood is
###          defined by removing the divergent -0.5*n_obs*[log(2*pi) +
###          log(kappa)] piece, leaving only -0.5*log|F_inf|. This is the
###          renormalization that distinguishes "diffuse" from "kappa": the
###          latter keeps a finite-kappa version of the divergent piece,
###          hence the kappa-dependent additive offset between the two.)
###   s'      = TT s + K0 v
###   P_inf'  = TT P_inf  TT' - K0 F_inf K0'
###   P_star' = TT P_star TT' + QQ - K0 M_star' - M_star K0' + K0 F_star K0'
###
### Derivation note: write P = P_star + kappa*P_inf, M = M_star + kappa*M_inf,
### F = F_star + kappa*F_inf, K = M F^{-1}, and expand
### P' = T P T' + QQ - M F^{-1} M' as kappa -> Inf. The O(kappa) terms give the
### P_inf recursion above (with K -> K0 = M_inf F_inf^{-1}); the O(1) terms
### give the P_star recursion. The only departure from DK 2012 sec. 5.2 is
### M_star = T P_star Z' + S, which absorbs dynhr's state/measurement
### cross-covariance S = RR Sigma_e DD'.
###
### Case C (F_inf singular but nonzero): not handled by the multivariate
### recursion. Signals via fallback = TRUE + fallback_univariate = TRUE so
### the caller restarts on the univariate diffuse filter (sequential
### processing on the augmented state handles F_inf of any rank -- see
### .kf_univariate_dispatch above).
###
### Missing observations during the diffuse phase are not supported here
### (fallback = TRUE, restart with lik_init = "kappa").
###
### Returns list(s, P, loglik, t_next, ok, d_diffuse, fallback, P_inf_final).
.kf_diffuse_phase <- function(Y, d, ZZ, TT, RR, DD, QQ, HH, SS, Sigma_e,
                              me_diag, s0, P_inf0, P_star0,
                              n_obs, n_T, ll_const, kf_step,
                              diffuse_tol = 1e-10, conv_tol = 1e-8,
                              max_periods = 100L) {
  s <- s0; P_inf <- P_inf0; P_star <- P_star0
  loglik <- 0
  cap <- min(n_T, max_periods)

  for (t in seq_len(cap)) {
    y_t <- Y[, t]
    if (anyNA(y_t)) {
      return(list(s = s, P = P_star, loglik = loglik, t_next = t,
                  ok = TRUE, d_diffuse = NA_integer_, fallback = TRUE,
                  P_inf_final = P_inf))
    }

    v <- y_t - d - as.numeric(ZZ %*% s)

    F_inf  <- ZZ %*% P_inf  %*% t(ZZ)
    F_inf  <- (F_inf + t(F_inf)) * 0.5
    F_star <- ZZ %*% P_star %*% t(ZZ) + HH + me_diag
    F_star <- (F_star + t(F_star)) * 0.5

    scale_star <- max(1, max(abs(F_star)))

    if (max(abs(F_inf)) < diffuse_tol * scale_star) {
      ## -- Case A: diffuse part uninformative this period -----------------
      step <- kf_step(s, P_star, v)
      if (is.null(step))
        return(list(s = s, P = P_star, loglik = -Inf, t_next = t,
                    ok = FALSE, d_diffuse = NA_integer_, fallback = FALSE,
                    P_inf_final = P_inf))
      loglik <- loglik + step$ll
      s      <- step$s
      P_star <- step$P
      P_inf  <- tcrossprod(TT %*% P_inf, TT)
      P_inf  <- (P_inf + t(P_inf)) * 0.5
    } else {
      rc <- tryCatch(rcond(F_inf), error = function(e) 0)
      Fc_inf <- tryCatch(chol(F_inf), error = function(e) NULL)

      if (is.null(Fc_inf) || !is.finite(rc) || rc <= diffuse_tol) {
        ## -- Case C: F_inf singular but nonzero -- univariate handles it --
        return(list(s = s, P = P_star, loglik = loglik, t_next = t,
                    ok = TRUE, d_diffuse = NA_integer_, fallback = TRUE,
                    fallback_univariate = TRUE, P_inf_final = P_inf))
      }

      ## -- Case B: F_inf nonsingular ---------------------------------------
      F_inf_inv <- chol2inv(Fc_inf)
      log_det_F_inf <- 2 * sum(log(diag(Fc_inf)))

      M_inf  <- TT %*% P_inf  %*% t(ZZ)
      M_star <- TT %*% P_star %*% t(ZZ) + SS

      K0 <- M_inf %*% F_inf_inv

      ## NOTE on the constant term: for kappa -> Inf with F_inf full rank
      ## (n_obs x n_obs here), log|F| = n_obs*log(kappa) + log|F_inf| + O(1/kappa)
      ## and the quadratic term v'F^{-1}v -> 0. The -0.5*n_obs*log(2*pi) and
      ## the divergent -0.5*n_obs*log(kappa) pieces are exactly the
      ## normalization that the "exact diffuse" likelihood removes (this is
      ## the renormalization that separates "diffuse" from "kappa": the
      ## latter keeps a finite-kappa version of the divergent piece, which is
      ## why the two differ by a kappa-dependent additive constant). What
      ## remains is just the Jacobian term for the information content of
      ## this period's observation about the diffuse state.
      ll_t <- -0.5 * log_det_F_inf
      if (!is.finite(ll_t)) {
        return(list(s = s, P = P_star, loglik = -Inf, t_next = t,
                    ok = FALSE, d_diffuse = NA_integer_, fallback = FALSE,
                    P_inf_final = P_inf))
      }
      loglik <- loglik + ll_t

      s_new <- as.numeric(TT %*% s) + as.numeric(K0 %*% v)

      P_inf_new <- tcrossprod(TT %*% P_inf, TT) - tcrossprod(K0 %*% F_inf, K0)
      P_inf_new <- (P_inf_new + t(P_inf_new)) * 0.5

      P_star_new <- tcrossprod(TT %*% P_star, TT) + QQ -
        K0 %*% t(M_star) - M_star %*% t(K0) + tcrossprod(K0 %*% F_star, K0)
      P_star_new <- (P_star_new + t(P_star_new)) * 0.5

      s <- s_new; P_inf <- P_inf_new; P_star <- P_star_new
    }

    if (max(abs(P_inf)) < conv_tol * max(1, max(abs(P_star)))) {
      return(list(s = s, P = P_star, loglik = loglik, t_next = t + 1L,
                  ok = TRUE, d_diffuse = t, fallback = FALSE,
                  P_inf_final = P_inf))
    }
  }

  ## Hard cap reached without convergence.
  warning("kf_diffuse_phase: P_inf did not converge to zero within ",
          cap, " periods; falling back to lik_init = \"kappa\".")
  list(s = s, P = P_star, loglik = loglik, t_next = cap + 1L,
       ok = TRUE, d_diffuse = NA_integer_, fallback = TRUE,
       P_inf_final = P_inf)
}


### -- Main filter ----------------------------------------------------------------

#' Kalman filter for DSGE models
#'
#' Evaluates the log-likelihood of a DSGE model under the Gaussian state-space
#' representation implied by the decision rules.  Supported \code{method}
#' values: \code{"auto"}, \code{"standard"}, \code{"dare"}, \code{"reference"}
#' (alias for \code{"dare"}), \code{"chandrasekhar"}, \code{"univariate"}.
#' There is no \code{"block"} method.
#'
#' @param Y observation matrix (\code{n_obs} x \code{T}).
#' @param dr decision rule (output of \code{\link{solve_perturbation}}).
#' @param model compiled model object (output of \code{\link{compile_model}}).
#' @param params named numeric vector of parameter values.
#' @param obs_vars character vector of observed variable names.
#' @param return_filtered logical; if \code{TRUE}, return filtered state estimates
#'   (one column per time step).
#' @param ss_tol tolerance for steady-state lock detection (default \code{.LYAP_TOL}).
#' @param me_variance scalar measurement-error variance added to the innovation
#'   covariance \code{F} in every likelihood evaluation (default \code{0}).
#'   See Details for implications.
#' @param return_ll_contrib logical; if \code{TRUE}, return per-step
#'   log-likelihood contributions (prediction-error decomposition).
#' @param method character; filtering algorithm: \code{"auto"} (default; picks
#'   \code{"standard"} or \code{"chandrasekhar"}), \code{"dare"} (textbook
#'   Kalman filter, no steady-state shortcut), \code{"chandrasekhar"}
#'   (low-rank recursion), \code{"standard"} (per-step Riccati with steady-state
#'   lock), \code{"reference"} (alias for \code{"dare"}), or
#'   \code{"univariate"} (Koopman--Durbin 2000 sequential filter on the
#'   augmented state \code{[s; eps]}; handles singular innovation
#'   covariances and the exact-diffuse phase with singular \code{F_inf}).
#' @param lik_init character; initialization of the state covariance \code{P0}
#'   (and, for unit-root models, the diffuse covariance): \code{"auto"}
#'   (default), \code{"stationary"}, \code{"diffuse"}, or \code{"kappa"}.
#'   See Details.
#' @param me_extra \code{n_obs x T} matrix of additional per-observable,
#'   per-period measurement-error variances (default \code{NULL}, no extra
#'   variance).  Unlike the scalar \code{me_variance} regularizer,
#'   \code{me_extra} is treated as TRUE per-period measurement noise on all
#'   paths: it enters the innovation covariance \code{F_t} AND the
#'   Joseph-form state-covariance update
#'   (\code{P += K_t diag(me_extra[, t]) t(K_t)}), so the \code{"standard"},
#'   \code{"dare"}/\code{"reference"}, and \code{"univariate"} methods agree
#'   exactly under \code{me_extra} (at \code{me_variance = 0}).
#'   Intended for \code{filter_tunes} soft tunes: the expanded
#'   observable's column carries \code{stderr^2} at the tune periods and 0
#'   elsewhere.  When \code{me_extra} is non-\code{NULL} and has any nonzero
#'   entry the filter is forced onto the per-period R loop (\code{method =
#'   "standard"}, \code{"dare"}, or \code{"univariate"} in R); the C++
#'   standard fast path, the Chandrasekhar recursion, and the steady-state
#'   lock are all bypassed
#'   because they bake a time-invariant \code{H} into the gain/covariance
#'   update.  Concretely: \code{method = "auto"} routes to
#'   \code{"standard"}; explicit \code{method = "chandrasekhar"} with
#'   nonzero \code{me_extra} raises an error; the DARE drift diagnostic
#'   (\code{dare_p_drift}) is skipped (\code{NA}) because \code{P} has no
#'   time-invariant fixed point under per-period \code{F_t}.  Note also that the
#'   exact-diffuse phase (\code{lik_init = "diffuse"}) does not support
#'   missing observations and falls back to \code{"kappa"} whenever
#'   \code{me_extra} introduces NA rows (see Landmine 2 in the
#'   filter-tunes brief).
#' @param shock_scale \code{n_exo x T} matrix of shock standard-deviation
#'   scale factors (multiplicative, so 1 = baseline).  At period \code{t}
#'   the effective shock covariance is
#'   \eqn{\\Sigma_{e,t} = \\mathrm{diag}(s_t)\\,\\Sigma_e\\,\\mathrm{diag}(s_t)}
#'   where \eqn{s_t} is column \code{t} of \code{shock_scale}.
#'   \code{NULL} (default) and an all-ones matrix are treated as identity
#'   (no heteroskedasticity). Rows must be ordered to match \code{dr$exo_names}.
#'   Incompatible with \code{method = "chandrasekhar"}, \code{"univariate"},
#'   and \code{lik_init = "diffuse"}. \code{P0} always uses the baseline
#'   (unscaled) \eqn{\\Sigma_e}.
#' @param me_floor_check Logical: when \code{me_variance > 0}, compare it
#'   against the smallest eigenvalue of the model-implied (ME-free) steady-
#'   state innovation covariance \code{F} and warn if the floor is large
#'   relative to that eigenvalue (near-collinear observables; see
#'   \code{.pruned_me_floor_ratio}). Default
#'   \code{getOption("dynhr.me_floor_check", TRUE)}. Only evaluated on the
#'   stationary (non-\code{shock_scale}) baseline system; scoped to the
#'   scalar \code{me_variance} floor, not \code{me_extra}.
#'
#' @details
#' \strong{Measurement-error convention:}
#' The default is \code{me_variance = 0}: no jitter is added to the
#' innovation covariance \code{F}, giving exact numeric parity with Dynare
#' (which adds no regularisation). Through dynhr 0.7 the default was
#' \code{1e-8} as a positive-definiteness safeguard; this is no longer
#' needed because a singular or ill-conditioned \code{F} on any
#' multivariate path now triggers an automatic fallback to the univariate
#' (sequential) filter, which processes observables one at a time and skips
#' zero-variance components instead of inverting \code{F} (the analog of
#' Dynare's \code{univariate_kalman_filter_if_singularity_is_detected}).
#' Set \code{me_variance > 0} to deliberately add measurement noise, e.g.
#' for stochastically singular models where the singularity is a modelling
#' choice rather than a numerical artifact. When comparing marginal
#' likelihoods (model comparison) across models, ensure all use the same
#' \code{me_variance} setting; the ranking is otherwise invalidated.
#' Convention note for \code{me_variance > 0}: the multivariate methods add
#' the jitter to \code{F} only (likelihood and gain), never to the
#' state-covariance update -- a regularisation, not a noise model. The
#' \code{"univariate"} method instead treats \code{me_variance} as TRUE
#' iid diagonal measurement noise (the statistically exact filter for the
#' noise-augmented model), so the two conventions agree exactly at
#' \code{me_variance = 0} but differ by \code{O(me_variance)} otherwise.
#' Do not mix methods across draws when \code{me_variance > 0}.
#'
#' \strong{Likelihood initialization (\code{lik_init}):}
#' \code{"stationary"} initialises \code{P0} via the discrete Lyapunov
#' equation (\code{solve_lyapunov(TT, QQ)}); this is the historical default
#' and is only valid when \code{TT} is stable (spectral radius < 1). On a
#' unit-root \code{TT} it returns a \code{NaN} \code{P0} and the filter
#' fails. \code{"diffuse"} implements the exact diffuse initialization of
#' Koopman & Durbin (2003) and Durbin & Koopman (2012, ch. 5): the unit-root
#' subspace gets an (improper, infinite-variance) diffuse prior handled
#' analytically via a separate \code{P_inf} recursion, with no arbitrary
#' tuning constant. \code{"kappa"} reproduces the legacy big-\code{kappa}
#' initialization (\code{dynhr:::.build_P0()}: Lyapunov on the stable Schur
#' block, \code{.DIFFUSE_SCALE = 1e6} on the unit-root block) -- this is an
#' approximation, and the resulting log-likelihood differs from
#' \code{"diffuse"} by an arbitrary, \code{kappa}-dependent additive
#' constant (it shifts the LEVEL of the log-likelihood but not differences
#' between parameter draws at the same \code{kappa}). \strong{Do not mix
#' \code{lik_init} settings within a marginal-likelihood comparison.}
#' \code{"auto"} (the default) chooses the cheapest \emph{valid} init by
#' letting the stationary Lyapunov solve decide. It inspects the eigenvalues
#' of \code{TT}: with all \code{|lambda| <= 1 - 1e-6} it uses
#' \code{"stationary"} directly. When a (near-)unit root is present it
#' \emph{attempts} \code{solve_lyapunov(TT, QQ)} and uses \code{"stationary"}
#' whenever the result is a valid covariance (finite and positive
#' semi-definite), falling back to \code{"diffuse"} only when the solve
#' genuinely fails (a \code{NaN}/non-PSD \code{P0}, i.e. a true unit or
#' explosive root such as a random-walk trend at \code{rho == 1}). This
#' matters for performance \emph{and} accuracy: a root in
#' \code{[1 - 1e-6, 1)} is still stationary, so \code{"stationary"} is both
#' the exact likelihood (the process really is stationary) and far cheaper
#' than the per-step exact-diffuse phase (~10 ms vs ~500-800 ms/eval on a
#' medium model). Near-unit-root mode-finding and MCMC sample the
#' feasibility boundary heavily, so this is the difference between a
#' tractable and an intractable run. The stationary \code{P0} computed
#' during resolution is cached and reused (no second Lyapunov solve).
#' \strong{Note:} \code{"auto"} is therefore \emph{not} bit-identical to a
#' forced \code{lik_init = "diffuse"} on near-unit-root models -- it returns
#' the (correct) stationary likelihood there. To force the exact-diffuse
#' treatment for every near-unit root, pass \code{lik_init = "diffuse"}
#' explicitly.
#'
#' @references
#'   Kalman, R. E. (1960). A new approach to linear filtering and prediction
#'     problems. \emph{Journal of Basic Engineering}, 82(1), 35-45.
#'   Anderson, B. D. O., & Moore, J. B. (1979). \emph{Optimal Filtering}.
#'     Prentice-Hall.
#'   Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for
#'     multivariate state space models. \emph{Journal of Time Series
#'     Analysis}, 21(3), 281-296.
#'   Strid, I., & Walentin, K. (2011). Block Kalman filtering for large-scale
#'     DSGE models. \emph{Computational Economics}, 39(2), 145-160.
#' @export
kalman_filter <- function(Y, dr, model, params, obs_vars,
                          return_filtered = FALSE, ss_tol = .LYAP_TOL,
                          me_variance = 0,
                          return_ll_contrib = FALSE,
                          method = c("auto", "dare", "chandrasekhar",
                                     "standard", "reference", "univariate",
                                     "univariate_ss"),
                          lik_init = c("auto", "stationary",
                                       "diffuse", "kappa"),
                          me_extra = NULL,
                          shock_scale = NULL,
                          me_floor_check = getOption("dynhr.me_floor_check",
                                                     TRUE)) {

  ## "reference" is an alias for "dare" -- both run the per-step textbook
  ## Kalman filter with no steady-state shortcut. "dare" is kept for
  ## backward compatibility but may be deprecated in a future release.
  method <- match.arg(method)
  if (method == "reference") method <- "dare"
  ## "univariate_ss" is the univariate filter with the opt-in steady-state lock
  ## (frozen gains on the converged stationary tail; an ~ss_tol approximation,
  ## the same one the multivariate steady-state filter makes). Treat it exactly
  ## like "univariate" everywhere except the final dispatch, which sets ss_lock.
  ss_lock_req <- (method == "univariate_ss")
  if (ss_lock_req) method <- "univariate"
  method_orig <- method          # for M23 warning below

  lik_init <- match.arg(lik_init)
  lik_init_orig <- lik_init      # for M23 warning below

  ## Per-period log-likelihood contributions (the prediction-error
  ## decomposition) are collected on the bit-exact per-step "dare" path; this
  ## is what D35 (misspecification softness) consumes. Force that path when the
  ## caller asks for contributions (accuracy over speed; not a per-draw path).
  ## The univariate filter collects its own per-period contributions.
  if (isTRUE(return_ll_contrib) && method != "univariate") method <- "dare"

  ## -- Common setup --------------------------------------------------

  state_idx <- dr$state_idx
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)
  n_obs     <- length(obs_vars)

  if (n_obs > n_exo)
    warning(sprintf("Stochastic singularity: %d obs but only %d shocks.", n_obs, n_exo))

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("Observed variables not found: ", paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d   <- dr$ys[obs_vars]

  Sigma_e <- .get_shock_cov(model, exo, params)
  QQ      <- tcrossprod(RR %*% Sigma_e, RR)
  HH      <- tcrossprod(DD %*% Sigma_e, DD)
  SS      <- RR %*% Sigma_e %*% t(DD)
  tZZ     <- t(ZZ)
  me_diag <- me_variance * diag(n_obs)

  ## -- Measurement-error floor guard (near-degenerate-F hazard) ------------
  ## Same detector as pruned_ss_loglik() (see R/pruned-state-space.R): a
  ## "negligible" me_variance floor can still dwarf the smallest eigenvalue
  ## of the model-implied (ME-free) steady-state innovation covariance when
  ## some linear combination of observables is nearly perfectly predictable.
  ## Reuses the stationary Lyapunov P0 as Sxi0. Per-draw closures latch this
  ## to TRUE only once (see make_log_posterior); direct calls run it every
  ## time, but the detector memoizes on system content, so repeated calls on
  ## an unchanged system cost only the Lyapunov solve + a hash (the
  ## un-memoized Riccati was ~13x a small-model KF sweep; kf_rbc_standard
  ## perf-gate regression, 2026-08-05).
  if (me_variance > 0 && isTRUE(me_floor_check)) {
    Sxi0_chk <- tryCatch(solve_lyapunov(TT, QQ), error = function(e) NULL)
    if (!is.null(Sxi0_chk) && all(is.finite(Sxi0_chk))) {
      .warn_me_floor_lock(
        .pruned_me_floor_ratio(TT, ZZ, QQ, HH, SS, Sxi0_chk, me_variance),
        obs_vars, me_variance)
    }
  }

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ll_const <- -0.5 * n_obs * log(2 * pi)
  has_missing <- anyNA(Y)

  ## ---- me_extra validation and routing ------------------------------------
  ## me_extra must be n_obs x T when non-NULL.
  has_me_extra <- FALSE
  if (!is.null(me_extra)) {
    if (!is.matrix(me_extra) || nrow(me_extra) != n_obs || ncol(me_extra) != n_T)
      stop(sprintf(
        "kalman_filter: me_extra must be an n_obs x T matrix (%d x %d); got %s.",
        n_obs, n_T,
        if (is.matrix(me_extra)) paste0(nrow(me_extra), " x ", ncol(me_extra))
        else "non-matrix"), call. = FALSE)
    ## Only treat as active when at least one entry is nonzero.
    has_me_extra <- any(me_extra != 0)
  }
  if (has_me_extra) {
    ## C++ standard fast path bakes a time-invariant H into the gain / update.
    ## Chandrasekhar likewise assumes a constant innovation structure.
    ## Neither can handle per-period me_extra: force R paths only.
    if (method == "chandrasekhar")
      stop("kalman_filter: method = 'chandrasekhar' is incompatible with ",
           "nonzero me_extra (time-varying ME variances require a per-period ",
           "R loop); use method = 'auto', 'standard', or 'univariate'.",
           call. = FALSE)
    if (method == "auto") method <- "standard"
    ## has_missing is already TRUE whenever me_extra introduces NA rows; if
    ## not, mark it so the C++ dispatch inside METHOD 3 is also bypassed.
    has_missing <- TRUE
  }

  ## ---- shock_scale validation and routing ---------------------------------
  ## shock_scale must be n_exo x T when non-NULL; rows must align to dr$exo_names.
  has_shock_scale <- FALSE
  if (!is.null(shock_scale)) {
    if (!is.matrix(shock_scale) || nrow(shock_scale) != n_exo || ncol(shock_scale) != n_T)
      stop(sprintf(
        "kalman_filter: shock_scale must be an n_exo x T matrix (%d x %d); got %s.",
        n_exo, n_T,
        if (is.matrix(shock_scale)) paste0(nrow(shock_scale), " x ", ncol(shock_scale))
        else "non-matrix"), call. = FALSE)
    ## Only treat as active when at least one entry differs from 1.
    has_shock_scale <- !all(shock_scale == 1)
  }
  if (has_shock_scale) {
    ## Chandrasekhar assumes constant low-rank innovation updates -- incompatible.
    if (method == "chandrasekhar")
      stop("kalman_filter: method = 'chandrasekhar' is incompatible with ",
           "shock_scale (time-varying shock variances require a per-period R loop); ",
           "use method = 'auto', 'standard', or 'dare'.", call. = FALSE)
    ## C++ univariate bakes QQb once; force the R standard path.
    if (method == "univariate") method <- "standard"
    if (method == "auto")       method <- "standard"
    ## C++ standard fast path bakes HH/Sigma_e -- bypass it.
    has_missing <- TRUE
  }
  ## Diffuse phase with non-unity scales: stop() (scaling inside the diffuse
  ## phase would require per-step P_inf updates not currently implemented).
  if (has_shock_scale && lik_init == "diffuse")
    stop("kalman_filter: lik_init = 'diffuse' is incompatible with shock_scale. ",
         "Use lik_init = 'kappa' or 'stationary' with heteroskedastic shocks.",
         call. = FALSE)

  ## Precompute Y - d (broadcast d down each column) once: the per-step loops
  ## then need only one matrix-vector subtraction. Safe with missing data: NA
  ## propagates through Y - d and is caught by the same anyNA / is.finite
  ## checks downstream.
  Y_minus_d <- Y - d

  ## -- Resolve lik_init = "auto" --------------------------------------
  ## Inspect the eigenvalues of TT. Roots well inside the unit circle => the
  ## stationary init (P0 from solve_lyapunov) is exact and fast, so use it.
  ##
  ## For roots within 1e-6 of the unit circle the eigenvalue test alone is too
  ## conservative: solve_lyapunov(TT, QQ) only actually DIVERGES (NaN P0) for
  ## roots AT/ABOVE 1. For roots in [1-1e-6, 1) the stationary variance is large
  ## but FINITE, and the stationary filter is both the EXACT likelihood (the
  ## process really is stationary) AND ~50x cheaper than the per-step exact-
  ## diffuse phase. Near-unit-root mode-finding/MCMC samples the feasibility
  ## boundary heavily, so eagerly picking "diffuse" there dominated wall time
  ## (~700 ms/eval vs ~11 ms). So: when a near-unit root is present, let the
  ## Lyapunov solve DECIDE -- use stationary when it converges (reusing the P0
  ## below), and fall back to the exact-diffuse init only when it genuinely
  ## fails (true unit/explosive root, e.g. a random-walk trend at rho == 1).
  P0_auto <- NULL
  if (lik_init == "auto") {
    ## symmetric = FALSE skips R's isSymmetric()/all.equal() probe (TT is the
    ## non-symmetric state transition); only Mod() of the eigenvalues is used,
    ## so the result is identical. Profiled at ~6-8% of an RWMH draw.
    tt_evals <- eigen(TT, symmetric = FALSE, only.values = TRUE)$values
    if (any(Mod(tt_evals) > 1 - 1e-6)) {
      P0_auto <- tryCatch(solve_lyapunov(TT, QQ), error = function(e) NULL)
      ## Accept the stationary init only if P0 is finite AND a valid covariance
      ## (PSD). A genuine unit root => NaN P0; an explosive root (|lambda| > 1,
      ## which a BK-satisfying dr never produces, but guard anyway) => a finite
      ## but NON-PSD P0 (negative variance). Both must fall back to diffuse.
      ok_stat <- !is.null(P0_auto) && all(is.finite(P0_auto)) &&
        min(Re(eigen((P0_auto + t(P0_auto)) / 2, symmetric = TRUE,
                     only.values = TRUE)$values)) > -1e-8
      lik_init <- if (ok_stat) "stationary" else "diffuse"
      if (lik_init == "diffuse") P0_auto <- NULL   # not reusable on the diffuse path
    } else {
      lik_init <- "stationary"
    }
  }
  ## has_shock_scale + diffuse is rejected above when the caller passes
  ## lik_init = "diffuse" literally, but lik_init = "auto" (the default) can
  ## also RESOLVE to "diffuse" here on a unit-root TT -- catch that case too,
  ## since the diffuse phase (.kf_diffuse_phase) has no shock_scale awareness.
  if (has_shock_scale && lik_init == "diffuse")
    stop("kalman_filter: lik_init = \"auto\" resolved to the diffuse ",
         "initialization (TT has unit-root eigenvalues), which is ",
         "incompatible with shock_scale. Pass lik_init = \"kappa\" or ",
         "\"stationary\" explicitly when using heteroskedastic shocks on a ",
         "nonstationary model.", call. = FALSE)
  d_diffuse <- NA_integer_

  ## -- M23: warn about diffuse loglik convention on unit-root models ------
  ## The univariate diffuse filter (method = "univariate") keeps the
  ## divergent 0.5*log(F_inf) term in the diffuse phase; the multivariate
  ## methods ("standard", "dare", "chandrasekhar") drop it (exact-diffuse
  ## renormalization).  On double-unit-root models this produces an additive
  ## gap of ~15-16 nats between conventions.  For STATIONARY models all
  ## methods agree to numerical tolerance.
  ## Only warn when the caller explicitly mixed method + lik_init in a way
  ## that risks cross-method comparison (method_orig != "auto" and
  ## method_orig != "univariate" and lik_init_orig != "auto").
  ## Normal usage (method = "auto", lik_init = "auto") never warns here.
  if (lik_init_orig %in% c("diffuse", "kappa") &&
      method_orig %in% c("standard", "dare", "chandrasekhar", "reference")) {
    has_unit_roots <- if (exists("tt_evals"))   # reuse from "auto" branch
      any(Mod(tt_evals) > 1 - 1e-6)
    else
      any(Mod(eigen(TT, symmetric = FALSE, only.values = TRUE)$values) > 1 - 1e-6)
    if (has_unit_roots) {
      warning("kalman_filter: unit-root model detected with method = \"",
              method, "\". ",
              "The univariate and multivariate diffuse methods use different ",
              "normalizations of the diffuse-phase log-determinant term: ",
              "log-likelihoods are NOT directly comparable across methods on ",
              "nonstationary models (gap can exceed 15 nats). ",
              "Use method = \"univariate\" throughout for cross-draw consistency, ",
              "or do not mix methods when comparing log-likelihoods.",
              call. = FALSE, immediate. = FALSE)
    }
  }

  if (method == "auto") {
    ## "dare" is now the exact, no-shortcut KF (see METHOD 1 below) - kept as a
    ## correctness oracle for testing, not as a speed path. "auto" therefore
    ## still picks the fast methods: standard with self-lock for the common
    ## case, chandrasekhar for large state vectors.
    if (has_missing)       method <- "standard"
    else if (lik_init %in% c("diffuse", "kappa")) {
      ## The exact-diffuse phase has no multivariate Rcpp fast path (the
      ## standard/chandrasekhar C++ kernels bake in a time-invariant gain), so
      ## auto used to run the diffuse phase in the slow per-step R "standard"
      ## loop -- and on a model with a singular innovation covariance it then
      ## fell THROUGH to the univariate filter anyway, paying the wasted R
      ## attempt first. The univariate C++ filter runs the exact-diffuse
      ## recursion natively (kalman_univariate_loop_cpp handles F_inf of any
      ## rank) at ~2x the speed.
      ##
      ## CONVENTION GATE: route to univariate only at me_variance == 0, where
      ## the univariate and multivariate diffuse log-likelihoods are identical
      ## (verified to 3e-11). For me_variance > 0 they DIVERGE by O(me_variance)
      ## -- the multivariate path adds the jitter to F only (regularisation)
      ## while univariate treats it as true iid measurement noise -- and the
      ## exact-diffuse ADJOINT gradient (make_posterior_grad) implements the
      ## multivariate convention, so rerouting the loglik would desync it from
      ## its own gradient. me_extra likewise forces the univariate R loop. Keep
      ## "standard" in both cases.
      method <- if (.HAS_RCPP_KALMAN_UNI() && is.null(me_extra) &&
                    me_variance == 0)
        "univariate" else "standard"
    }
    else if (n_state > 50) method <- "chandrasekhar"
    else                   method <- "standard"
  }

  ## A diffuse phase requires the per-step R loop (no Rcpp / Chandrasekhar
  ## fast path); for large state vectors with unit roots, fall back to the
  ## R "standard" loop entirely. This is a perf cost but correctness-first.
  use_diffuse_phase <- lik_init %in% c("diffuse", "kappa") &&
    method %in% c("standard", "dare", "chandrasekhar")
  if (use_diffuse_phase && method == "chandrasekhar") method <- "standard"

  filtered <- if (return_filtered) matrix(0, n_state, n_T) else NULL

  ## -- Guard: solve_lyapunov for lik_init = "stationary" -----------------
  ## solve_lyapunov() returns a NaN matrix when TT has unit eigenvalues; on
  ## the "stationary" path that NaN P0 silently propagates through all filter
  ## quantities and produces loglik = 0 with no error.  This helper throws a
  ## clear diagnostic so the caller can switch to lik_init = "auto" or
  ## "diffuse" instead.
  .solve_lyapunov_stationary <- function() {
    ## Reuse the P0 already computed while resolving lik_init = "auto" (avoids a
    ## second solve_lyapunov on the near-unit-root auto path).
    if (!is.null(P0_auto)) return(P0_auto)
    P0 <- solve_lyapunov(TT, QQ)
    if (anyNA(P0))
      stop("kalman_filter: lik_init = \"stationary\" failed because ",
           "solve_lyapunov() returned NaN -- TT has unit-root eigenvalues. ",
           "Use lik_init = \"auto\" or \"diffuse\" for nonstationary models.",
           call. = FALSE)
    P0
  }

  ## -- Shared Phase-1 KF step (used by dare, chandrasekhar, and the
  ##    diffuse-phase Case A below -- must be defined BEFORE the lik_init
  ##    initialization block, which may call it) -----
  ## Returns list(ll, s, P, K, F_inv, log_det_F) after one full KF step

  .kf_step <- function(s, P, v) {
    PZ <- P %*% tZZ
    Ft <- ZZ %*% PZ + HH + me_diag
    Ft <- (Ft + t(Ft)) * 0.5
    ## Singular / non-PD F: return NULL so the caller can fall back to the
    ## univariate filter (see .kf_fail) instead of erroring out.
    Fc <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) return(NULL)
    Fi  <- chol2inv(Fc)
    ldf <- 2 * sum(log(diag(Fc)))
    ll  <- ll_const - 0.5 * (ldf + drop(crossprod(v, Fi %*% v)))
    if (!is.finite(ll) || ll < .KF_LL_MIN) return(NULL)
    K    <- (TT %*% PZ + SS) %*% Fi
    s_n  <- drop(TT %*% s) + drop(K %*% v)
    TmKZ <- TT - K %*% ZZ
    RmKD <- RR - K %*% DD
    P_n  <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD)
    P_n  <- (P_n + t(P_n)) * 0.5
    list(ll = ll, s = s_n, P = P_n, K = K, F_inv = Fi, log_det_F = ldf, F_mat = Ft)
  }

  ## -- Univariate (sequential) filter runner ----------------------------
  ## Runs the Koopman-Durbin (2000) univariate filter on the augmented
  ## state [s; eps] (see .kf_univariate_dispatch). Used (a) directly for
  ## method = "univariate", (b) as the automatic singularity fallback from
  ## every multivariate path (.kf_fail below), and (c) for the exact
  ## diffuse phase when F_inf is singular but nonzero (Case C).

  .run_univariate <- function(li, ss_lock = FALSE) {
    P_inf_state <- NULL
    if (li == "diffuse") {
      P0u <- .kf_diffuse_P0(TT, QQ)
      if (P0u$nunit == 0L) {
        li      <- "stationary"
        P_state <- solve_lyapunov(TT, QQ)
      } else {
        P_state     <- P0u$P_star
        P_inf_state <- P0u$P_inf
      }
    } else if (li == "kappa") {
      P_state <- .build_P0(TT, QQ)
    } else {
      P_state <- .solve_lyapunov_stationary()
    }
    ## has_shock_scale is captured from the enclosing kalman_filter() call:
    ## every .run_univariate() caller (method = "univariate" directly, the
    ## .kf_fail() singularity fallback from "standard"/"dare", and the exact-
    ## diffuse Case C restart) must honor shock_scale the same way the
    ## multivariate paths do -- see .kf_univariate_dispatch / .kf_univariate_loop_R.
    ss_arg <- if (has_shock_scale) shock_scale else NULL
    out <- .kf_univariate_dispatch(Y_minus_d, ZZ, TT, RR, DD, Sigma_e,
                                   numeric(n_state), P_state, P_inf_state,
                                   me_variance, return_filtered,
                                   me_extra = me_extra, ss_lock = ss_lock,
                                   shock_scale = ss_arg)
    if (isTRUE(out$diffuse_failed)) {
      ## P_inf never decayed (unobserved unit root) or the sample ended
      ## inside the diffuse phase: same kappa fallback as the multivariate
      ## diffuse path.
      warning("kalman_filter: univariate diffuse phase: P_inf did not ",
              "converge to zero; falling back to lik_init = \"kappa\".")
      li  <- "kappa"
      out <- .kf_univariate_dispatch(Y_minus_d, ZZ, TT, RR, DD, Sigma_e,
                                     numeric(n_state), .build_P0(TT, QQ),
                                     NULL, me_variance, return_filtered,
                                     me_extra = me_extra, shock_scale = ss_arg)
    }
    if (!isTRUE(out$ok))
      return(list(loglik = -Inf, filtered_states = NULL,
                  n_obs = n_obs, n_T = n_T, method = "univariate",
                  lik_init = li, d_diffuse = NA_integer_))
    filt <- NULL
    if (return_filtered && !is.null(out$filtered)) {
      filt <- out$filtered
      rownames(filt) <- endo[state_idx]
    }
    dd <- out$d_diffuse
    if (is.null(dd) || length(dd) != 1L) dd <- NA_integer_
    list(loglik = out$loglik, filtered_states = filt,
         loglik_contrib = if (return_ll_contrib) as.numeric(out$ll_contrib)
                          else NULL,
         n_obs = n_obs, n_T = n_T, method = "univariate",
         lik_init = li, d_diffuse = as.integer(dd))
  }

  ## A multivariate path hit a singular / non-PD innovation covariance (or
  ## a numerically exploded step): retry on the univariate filter before
  ## declaring -Inf -- the analog of Dynare's
  ## univariate_kalman_filter_if_singularity_is_detected. A genuinely bad
  ## draw fails the same per-period ll floor there too, so this never turns
  ## a true -Inf into a finite value.
  .kf_fail <- function(failed_method) {
    out <- tryCatch(.run_univariate(lik_init), error = function(e) NULL)
    if (!is.null(out)) return(out)
    list(loglik = -Inf, filtered_states = NULL,
         n_obs = n_obs, n_T = n_T, method = failed_method,
         lik_init = lik_init, d_diffuse = d_diffuse)
  }

  if (method == "univariate") return(.run_univariate(lik_init, ss_lock = ss_lock_req))

  ## -- Initialization (s0, P0) for the chosen lik_init -----------------
  ## Computed once, used by both METHOD 1 ("dare") and METHOD 3 ("standard").
  ## "kappa" and "diffuse" require methods %in% c("standard", "dare") (forced
  ## above); "stationary" / "auto"-resolved-to-"stationary" reuse the
  ## historical solve_lyapunov(TT, QQ) P0 unconditionally.
  init_s <- numeric(n_state)
  init_P <- NULL
  init_loglik <- 0
  init_t_start <- 1L

  if (lik_init == "kappa") {
    init_P <- .build_P0(TT, QQ)
  } else if (lik_init == "diffuse") {
    if (has_missing) {
      warning("kalman_filter: lik_init = \"diffuse\" does not support ",
              "missing observations during the diffuse phase; falling ",
              "back to lik_init = \"kappa\".")
      lik_init <- "kappa"
      init_P <- .build_P0(TT, QQ)
    } else {
      P0 <- .kf_diffuse_P0(TT, QQ)
      diff_out <- .kf_diffuse_phase(Y, d, ZZ, TT, RR, DD, QQ, HH, SS, Sigma_e,
                                    me_diag, init_s, P0$P_inf, P0$P_star,
                                    n_obs, n_T, ll_const, .kf_step)
      if (!isTRUE(diff_out$ok))
        return(.kf_fail(method))
      if (isTRUE(diff_out$fallback)) {
        if (isTRUE(diff_out$fallback_univariate)) {
          ## Case C (F_inf singular but nonzero): exactly what the
          ## univariate diffuse filter is for. Restart from t = 1 on it.
          return(.run_univariate("diffuse"))
        }
        lik_init <- "kappa"
        init_P <- .build_P0(TT, QQ)
      } else {
        init_s       <- diff_out$s
        init_P       <- diff_out$P
        init_loglik  <- diff_out$loglik
        init_t_start <- diff_out$t_next
        d_diffuse    <- as.integer(diff_out$d_diffuse)
      }
    }
  }

  ## ===================================================================
  ## METHOD 1: "dare" -- exact full Kalman filter (no SS shortcut)
  ## ===================================================================
  ##
  ## History: the previous "dare" branch ran std-KF until |P_new - P| < tol
  ## and then swapped in DARE's exact (K_ss, F_ss). On nk_canonical with
  ## me_variance = 1e-6 this disagreed with the textbook KF by ~1 nat
  ## (sometimes ~6 nats on T=1000) because std-KF's P stabilises before it
  ## reaches DARE's fixed point, and DARE's K computed from the true fixed
  ## point is then inconsistent with std-KF's drifted P / s trajectory. The
  ## discontinuity at the swap, not "Joseph-form artifacts", is what shows
  ## up in the loglik.
  ##
  ## Resolution (2026-05-17): method = "dare" now runs the textbook KF with
  ## per-step Riccati updates for every t -- no SS lock at all. This is the
  ## slowest but bit-exact path; "standard" and "chandrasekhar" remain
  ## available for speed (both are now equally accurate in the common case).
  ## The DARE solver itself is still called as a consistency check: we warn
  ## if the KF's final P drifts far from the DARE fixed point. See
  ## inst/benchmarks/diag_dare2.R for the numerical evidence.

  if (method == "dare") {
    P <- if (lik_init == "stationary") .solve_lyapunov_stationary() else init_P

    ## The DARE fixed point is a diagnostic for the TIME-INVARIANT system
    ## only: with per-period me_extra (or shock_scale) F_t varies over t, P
    ## has no steady state, and the drift check is meaningless -- skip it.
    dare_diag_ok <- lik_init == "stationary" && !has_me_extra &&
      !has_shock_scale
    dare <- if (dare_diag_ok)
      tryCatch(.solve_dare(TT, ZZ, QQ, HH, SS, me_diag, tol = .DARE_TOL,
                           max_iter = 500),
               error = function(e)       # chol failure on singular F: the
                 list(converged = FALSE, # textbook loop below still runs
                      P = NULL, iterations = NA_integer_))
    else
      list(converged = FALSE, P = NULL, iterations = NA_integer_)
    if (dare_diag_ok && !dare$converged)
      warning("DARE solver did not converge; method=\"dare\" still proceeds ",
              "as exact textbook KF.")

    s <- init_s; loglik <- init_loglik
    ll_contrib <- if (return_ll_contrib) numeric(n_T) else NULL

    t_start <- init_t_start
    if (t_start > n_T) {
      ## The diffuse phase consumed the entire sample.
      if (return_filtered) rownames(filtered) <- endo[state_idx]
      return(list(loglik = loglik, filtered_states = filtered,
                  loglik_contrib = ll_contrib,
                  n_obs = n_obs, n_T = n_T, method = "dare",
                  dare_iterations = NA_integer_, dare_p_drift = NA_real_,
                  lik_init = lik_init, d_diffuse = d_diffuse))
    }

    for (t in t_start:n_T) {
      v    <- Y[, t] - as.numeric(ZZ %*% s) - d
      ## me_extra active at this period? Per-period F diagonal + Joseph term
      ## below -- previously the dare path silently IGNORED me_extra whenever
      ## F was nonsingular (.kf_step closure-captures the time-invariant
      ## me_diag), which also poisoned every return_ll_contrib caller (that
      ## flag forces method = "dare").
      me_x_t <- has_me_extra && any(me_extra[, t] != 0)
      if (has_shock_scale || me_x_t) {
        ## Inline per-period step for dare path (cannot use .kf_step which
        ## closure-captures the baseline HH/SS/Sigma_e and the time-invariant
        ## me_diag -- Landmine 6). Handles shock_scale, me_extra, or both.
        if (has_shock_scale) {
          sc_t  <- shock_scale[, t]
          Se_t  <- Sigma_e * outer(sc_t, sc_t)
          HH_t  <- tcrossprod(DD %*% Se_t, DD)
          SS_t  <- RR %*% Se_t %*% t(DD)
        } else {
          Se_t  <- Sigma_e
          HH_t  <- HH
          SS_t  <- SS
        }
        me_diag_t <- if (me_x_t) me_diag + diag(me_extra[, t], nrow = n_obs)
                     else me_diag
        PZ    <- P %*% t(ZZ)
        Ft    <- ZZ %*% PZ + HH_t + me_diag_t
        Ft    <- (Ft + t(Ft)) * 0.5
        Fc    <- tryCatch(chol(Ft), error = function(e) NULL)
        if (is.null(Fc)) return(.kf_fail("dare"))
        Fi    <- chol2inv(Fc)
        ldf   <- 2 * sum(log(diag(Fc)))
        ll_t  <- ll_const - 0.5 * (ldf + drop(crossprod(v, Fi %*% v)))
        if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(.kf_fail("dare"))
        K_t   <- (TT %*% PZ + SS_t) %*% Fi
        s     <- as.numeric(TT %*% s) + drop(K_t %*% v)
        TmKZ  <- TT - K_t %*% ZZ; RmKD <- RR - K_t %*% DD
        P     <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Se_t, RmKD)
        ## Joseph true-noise term for the me_extra part only:
        ## P' += K_t diag(me_extra[, t]) K_t'. me_variance (the base me_diag)
        ## stays on the documented F-only-regularizer convention (no term) --
        ## mirrors the fixed "standard" branch, so dare == standard ==
        ## univariate under me_extra.
        if (me_x_t) P <- P + K_t %*% (me_extra[, t] * t(K_t))
        P     <- (P + t(P)) * 0.5
        loglik <- loglik + ll_t
        if (return_ll_contrib) ll_contrib[t] <- ll_t
      } else {
        step <- .kf_step(s, P, v)
        if (is.null(step)) return(.kf_fail("dare"))
        loglik <- loglik + step$ll
        if (return_ll_contrib) ll_contrib[t] <- step$ll
        s      <- step$s
        P      <- step$P
      }
      if (return_filtered) filtered[, t] <- s
    }

    final_drift <- if (isTRUE(dare$converged)) max(abs(P - dare$P)) else NA_real_

    if (return_filtered) rownames(filtered) <- endo[state_idx]
    return(list(loglik = loglik, filtered_states = filtered,
                loglik_contrib = ll_contrib,
                n_obs = n_obs, n_T = n_T, method = "dare",
                dare_iterations  = dare$iterations,
                dare_p_drift     = final_drift,
                lik_init = lik_init, d_diffuse = d_diffuse))
  }


  ## ===================================================================
  ## METHOD 2: Chandrasekhar (bootstrapped initialization)
  ## ===================================================================
  ##
  ## Problem with W???=K???, M???=-F???: the first F update
  ##   F??? = F??? - Z K??? F??? K???' Z'
  ## nearly cancels when H is small (me_variance ~ 1e-8), losing PD.
  ##
  ## Fix: run N_BOOT standard Riccati steps from P???, then initialise
  ## Chandrasekhar from (K_b, F_b, ??P_b) where ??P is small enough
  ## that the F update doesn't catastrophically cancel.
  ##
  ## N_BOOT is chosen adaptively: keep going until |??P| < 0.1 * |P|.

  if (method == "chandrasekhar") {
    HH_full <- HH + me_diag

    ## --- Bootstrap: standard KF steps until ??P is small ---
    P <- .solve_lyapunov_stationary()

    s <- numeric(n_state); loglik <- 0
    K_prev <- NULL; F_prev <- NULL; F_inv_prev <- NULL; ldf_prev <- NULL
    boot_steps <- 0L
    ss_reached <- FALSE

    ## Process observations with standard KF until ??P is small
    for (t in seq_len(n_T)) {
      v <- Y[, t] - drop(ZZ %*% s) - d

      step <- .kf_step(s, P, v)
      if (is.null(step))
        return(.kf_fail("chandrasekhar"))

      loglik <- loglik + step$ll
      s <- step$s

      delta_P <- max(abs(step$P - P))
      scale_P <- max(abs(P))

      if (return_filtered) filtered[, t] <- s

      ## Check if P has self-converged (then skip Chandrasekhar entirely)
      if (t > 1L && delta_P < ss_tol) {
        ss_reached  <- TRUE
        boot_steps  <- t
        K_prev      <- step$K
        F_inv_prev  <- step$F_inv
        ldf_prev    <- step$log_det_F
        P <- step$P
        break
      }

      ## Check if ??P is small enough to start Chandrasekhar safely
      ## Threshold: |??P| < 10% of |P| ensures F update doesn't cancel
      if (t >= 3L && delta_P < 0.1 * scale_P) {
        boot_steps <- t
        K_prev     <- step$K
        F_prev     <- step$F_mat
        F_inv_prev <- step$F_inv
        ldf_prev   <- step$log_det_F
        P <- step$P
        break
      }

      P <- step$P
      K_prev     <- step$K
      F_prev     <- step$F_mat
      F_inv_prev <- step$F_inv
      ldf_prev   <- step$log_det_F
    }

    if (ss_reached || boot_steps >= n_T) {
      ## Already converged during bootstrap -- finish with steady-state
      if (ss_reached && boot_steps < n_T) {
        ll_ss_const <- ll_const - 0.5 * ldf_prev
        TT_K_ss     <- cbind(TT, K_prev)
        sv           <- numeric(n_state + n_obs)

        for (t2 in (boot_steps + 1L):n_T) {
          v2   <- Y[, t2] - as.numeric(ZZ %*% s) - d
          Fv2  <- F_inv_prev %*% v2
          ll_t <- ll_ss_const - 0.5 * sum(v2 * Fv2)
          if (!is.finite(ll_t) || ll_t < .KF_LL_MIN)
            return(.kf_fail("chandrasekhar"))
          loglik <- loglik + ll_t
          s <- as.numeric(TT %*% s + K_prev %*% v2)
          if (return_filtered) filtered[, t2] <- s
        }
      }

      if (return_filtered) rownames(filtered) <- endo[state_idx]
      return(list(loglik = loglik, filtered_states = filtered,
                  n_obs = n_obs, n_T = n_T, method = "chandrasekhar",
                  boot_steps = boot_steps,
                  ss_reached_at = if (ss_reached) boot_steps else NA_integer_,
                  lik_init = lik_init, d_diffuse = d_diffuse))
    }

    ## --- Chandrasekhar phase: initialise from bootstrap endpoint ---
    ## We have K_b, F_b, P_b from the last bootstrap step.
    ## Do one more Riccati step to get P_{b+1} and compute ??P.

    ## One explicit Riccati step to get the next P
    PZ_b  <- P %*% tZZ
    F_b   <- ZZ %*% PZ_b + HH_full
    F_b   <- (F_b + t(F_b)) * 0.5
    Fc_b  <- tryCatch(chol(F_b), error = function(e) NULL)
    if (is.null(Fc_b)) return(.kf_fail("chandrasekhar"))
    Fi_b  <- chol2inv(Fc_b)
    K_b   <- (TT %*% PZ_b + SS) %*% Fi_b
    P_next <- tcrossprod(TT %*% P, TT) + QQ - tcrossprod(K_b %*% F_b, K_b)
    P_next <- (P_next + t(P_next)) * 0.5

    delta_P_mat <- P_next - P   # small by construction (bootstrap ensured this)

    ## Low-rank factorisation of ??P via eigen
    eig  <- eigen(delta_P_mat, symmetric = TRUE)
    keep <- abs(eig$values) > 1e-14 * max(abs(eig$values))
    r    <- max(1L, sum(keep))
    W <- eig$vectors[, keep, drop = FALSE]
    M <- diag(eig$values[keep], nrow = r, ncol = r)

    ## Current Chandrasekhar state
    K         <- K_b
    F_mat     <- F_b
    F_inv     <- Fi_b
    log_det_F <- 2 * sum(log(diag(Fc_b)))

    ## Pre-allocate
    TT_K_ss <- NULL; sv <- numeric(n_state + n_obs)
    K_ss <- NULL; F_inv_ss <- NULL; log_det_F_ss <- NULL; ll_ss_const <- NULL
    ch_ss_step <- NA_integer_

    start_t <- boot_steps + 1L

    for (t in start_t:n_T) {
      v <- Y[, t] - drop(ZZ %*% s) - d

      if (!ss_reached) {
        ## Likelihood with current K, F
        ll_t <- ll_const - 0.5 * (log_det_F + drop(crossprod(v, F_inv %*% v)))
        if (!is.finite(ll_t) || ll_t < .KF_LL_MIN)
          return(.kf_fail("chandrasekhar"))
        loglik <- loglik + ll_t
        s <- drop(TT %*% s) + drop(K %*% v)

        ## Chandrasekhar update (low-rank)
        ZW  <- ZZ %*% W                       # n_obs x r
        TW  <- TT %*% W                       # n_s x r
        ZWM <- ZW %*% M                       # n_obs x r

        F_new <- F_mat + ZWM %*% t(ZW)
        F_new <- (F_new + t(F_new)) * 0.5
        F_new_chol <- tryCatch(chol(F_new), error = function(e) NULL)
        if (is.null(F_new_chol)) return(.kf_fail("chandrasekhar"))

        F_new_inv     <- chol2inv(F_new_chol)
        log_det_F_new <- 2 * sum(log(diag(F_new_chol)))

        K_new <- K + TW %*% M %*% t(ZW) %*% F_new_inv
        W_new <- TW - K_new %*% ZW
        M_new <- M + t(ZW) %*% F_new_inv %*% ZWM
        M_new <- (M_new + t(M_new)) * 0.5

        ZW_new <- ZZ %*% W_new

        if (max(abs(K_new - K)) < ss_tol) {
          ss_reached   <- TRUE
          ch_ss_step   <- t
          K_ss         <- K_new
          F_inv_ss     <- F_new_inv
          log_det_F_ss <- log_det_F_new
          ll_ss_const  <- ll_const - 0.5 * log_det_F_ss
          TT_K_ss      <- cbind(TT, K_ss)
        }

        K <- K_new; F_mat <- F_new; F_inv <- F_new_inv
        log_det_F <- log_det_F_new; W <- W_new; M <- M_new

      } else {
        ## Steady-state phase
        Fv   <- F_inv_ss %*% v
        ll_t <- ll_ss_const - 0.5 * sum(v * Fv)
        if (!is.finite(ll_t) || ll_t < -1e8)
          return(.kf_fail("chandrasekhar"))
        loglik <- loglik + ll_t
        s <- as.numeric(TT %*% s + K_ss %*% v)
      }
      if (return_filtered) filtered[, t] <- s
    }

    if (return_filtered) rownames(filtered) <- endo[state_idx]
    return(list(loglik = loglik, filtered_states = filtered,
                n_obs = n_obs, n_T = n_T, method = "chandrasekhar",
                boot_steps = boot_steps,
                ss_reached_at = ch_ss_step,
                lik_init = lik_init, d_diffuse = d_diffuse))
  }


  ## ===================================================================
  ## METHOD 3: Standard KF (handles missing data)
  ## ===================================================================

  P <- if (lik_init == "stationary") .solve_lyapunov_stationary() else init_P

  s <- init_s; loglik <- init_loglik
  ss_reached <- FALSE
  K_ss <- NULL; F_inv_ss <- NULL; log_det_F_ss <- NULL
  TT_K_ss <- NULL; sv <- numeric(n_state + n_obs); ll_ss_const <- NULL

  ## Fast path: run the entire standard filter (transient + steady-state lock
  ## + tail) in one C++ call. Only when there is no missing data (the per-step
  ## partial-observation handling below needs the R loop), the Rcpp backend
  ## is available, AND no diffuse phase was needed (lik_init == "stationary").
  ## A diffuse phase always uses the per-step R loop below (init_t_start > 1
  ## and/or non-Lyapunov P0) -- see use_diffuse_phase above. Bit-parity with
  ## the R loop is asserted by test-kalman-rcpp-parity.R.
  ## has_shock_scale is excluded explicitly (belt-and-suspenders): it already
  ## forces has_missing <- TRUE above (the C++ kernel bakes a single
  ## time-invariant Sigma_e/HH/SS and cannot express a per-period scale), but
  ## gate on it directly here too so this fast path can never silently ignore
  ## shock_scale if that has_missing coupling is ever loosened.
  if (!has_missing && !has_shock_scale && lik_init == "stationary" &&
      .HAS_RCPP_KALMAN()) {
    out <- kalman_standard_loop_cpp(Y_minus_d, ZZ, TT, RR, DD, HH + me_diag,
                                    Sigma_e, SS, P, ll_const, ss_tol,
                                    .KF_LL_MIN, return_filtered)
    if (!out$ok)
      return(.kf_fail("standard"))
    filt <- NULL
    if (return_filtered) {
      filt <- out$filtered
      rownames(filt) <- endo[state_idx]
    }
    return(list(loglik = out$loglik, filtered_states = filt,
                n_obs = n_obs, n_T = n_T, method = "standard",
                lik_init = lik_init, d_diffuse = d_diffuse))
  }

  if (init_t_start > n_T) {
    ## The diffuse phase consumed the entire sample.
    if (return_filtered) rownames(filtered) <- endo[state_idx]
    return(list(loglik = loglik, filtered_states = filtered,
                n_obs = n_obs, n_T = n_T, method = "standard",
                lik_init = lik_init, d_diffuse = d_diffuse))
  }

  for (t in init_t_start:n_T) {
    v <- Y_minus_d[, t] - as.numeric(ZZ %*% s)

    if (any(is.na(v))) {
      obs_ok <- which(!is.na(v))
      ## Compute per-period shock covariance for this t (used in all branches).
      Se_miss <- if (has_shock_scale) {
        sc_t_m <- shock_scale[, t]
        Sigma_e * outer(sc_t_m, sc_t_m)
      } else Sigma_e
      QQ_miss <- if (has_shock_scale) tcrossprod(RR %*% Se_miss, RR) else QQ
      if (length(obs_ok) == 0L) {
        s <- drop(TT %*% s)
        P <- tcrossprod(TT %*% P, TT) + QQ_miss; P <- (P + t(P)) * 0.5
        ## A fully-missing period is a pure prediction step: it moves P off the
        ## steady-state fixed point, so the steady-state lock (if engaged) is no
        ## longer valid and must be recomputed. Without this, subsequent
        ## observed periods reuse the stale K_ss / F_inv_ss / ll_ss_const from
        ## the pre-gap steady state and the loglik is silently wrong (the
        ## partial-missing branch below already resets the lock for the same
        ## reason).
        ss_reached <- FALSE
        if (return_filtered) filtered[, t] <- s
        next
      }
      ZZ_t <- ZZ[obs_ok, , drop = FALSE]; DD_t <- DD[obs_ok, , drop = FALSE]
      HH_t <- tcrossprod(DD_t %*% Se_miss, DD_t)
      SS_t_miss <- RR %*% Se_miss %*% t(DD_t)
      v <- v[obs_ok]; n_obs_t <- length(obs_ok)
      me_t <- me_variance * diag(n_obs_t)
      if (has_me_extra) diag(me_t) <- diag(me_t) + me_extra[obs_ok, t]
      Ft <- ZZ_t %*% P %*% t(ZZ_t) + HH_t + me_t
      Fc <- tryCatch(chol(Ft), error = function(e) NULL)
      if (is.null(Fc)) return(.kf_fail("standard"))
      Fi <- chol2inv(Fc); ldf <- 2 * sum(log(diag(Fc)))
      ## Constant is -0.5 * n_obs_t * log(2*pi): correct ll_const (which
      ## bakes in the full n_obs) UP by the number of missing components.
      loglik <- loglik + ll_const + (n_obs - n_obs_t) * 0.5 * log(2*pi) -
        0.5 * (ldf + drop(crossprod(v, Fi %*% v)))
      K <- (TT %*% P %*% t(ZZ_t) + SS_t_miss) %*% Fi
      s <- drop(TT %*% s) + drop(K %*% v)
      TmKZ <- TT - K %*% ZZ_t; RmKD <- RR - K %*% DD_t
      P <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Se_miss, RmKD)
      ## Joseph true-noise term for the me_extra part (observed subset only):
      ## y = Z s + D e + u with Var(u) = diag(me_extra[obs_ok, t]) requires
      ## P' += K diag(me_extra[obs_ok, t]) K' for ANY gain K. me_variance
      ## stays on the documented F-only-regularizer convention (no term).
      if (has_me_extra && any(me_extra[obs_ok, t] != 0))
        P <- P + K %*% (me_extra[obs_ok, t] * t(K))
      P <- (P + t(P)) * 0.5; ss_reached <- FALSE
      if (return_filtered) filtered[, t] <- s; next
    }

    if (!ss_reached) {
      ## When shock_scale is active, compute per-period scaled covariances:
      ##   Se_t = diag(s_t) Sigma_e diag(s_t)
      ##   QQ_t = RR Se_t RR',  HH_t = DD Se_t DD',  SS_t = RR Se_t DD'
      ## This inline per-period step does NOT modify .kf_step() (Landmine 6).
      if (has_shock_scale) {
        sc_t    <- shock_scale[, t]
        Se_t    <- Sigma_e * outer(sc_t, sc_t)   ## diag(s_t) Sigma_e diag(s_t)
        QQ_t    <- tcrossprod(RR %*% Se_t, RR)
        HH_t    <- tcrossprod(DD %*% Se_t, DD)
        SS_t    <- RR %*% Se_t %*% t(DD)
        me_diag_t <- if (has_me_extra && any(me_extra[, t] != 0))
                       me_diag + diag(me_extra[, t], nrow = n_obs)
                     else me_diag
        PZ   <- P %*% t(ZZ)
        Ft_t <- ZZ %*% PZ + HH_t + me_diag_t
        Ft_t <- (Ft_t + t(Ft_t)) * 0.5
        Fc_t <- tryCatch(chol(Ft_t), error = function(e) NULL)
        if (is.null(Fc_t)) return(.kf_fail("standard"))
        Fi_t  <- chol2inv(Fc_t)
        ldf_t <- 2 * sum(log(diag(Fc_t)))
        ll_t  <- ll_const - 0.5 * (ldf_t + drop(crossprod(v, Fi_t %*% v)))
        if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(.kf_fail("standard"))
        K_t   <- (TT %*% PZ + SS_t) %*% Fi_t
        s     <- as.numeric(TT %*% s) + drop(K_t %*% v)
        TmKZ  <- TT - K_t %*% ZZ; RmKD <- RR - K_t %*% DD
        P     <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Se_t, RmKD)
        ## Joseph true-noise term for the me_extra part of me_diag_t only:
        ## P' += K diag(me_extra[, t]) K'. me_variance (the base me_diag)
        ## stays on the documented F-only-regularizer convention (no term).
        if (has_me_extra && any(me_extra[, t] != 0))
          P <- P + K_t %*% (me_extra[, t] * t(K_t))
        P     <- (P + t(P)) * 0.5
        loglik <- loglik + ll_t
      ## When me_extra is active at this period (and no shock_scale), use a
      ## per-period F that adds me_extra[, t] to the diagonal instead of the
      ## time-invariant me_diag.  The steady-state lock (ss_reached) is
      ## suppressed for the entire run when has_me_extra is TRUE.
      } else if (has_me_extra && any(me_extra[, t] != 0)) {
        me_diag_t <- me_diag + diag(me_extra[, t], nrow = n_obs)
        PZ <- P %*% t(ZZ)
        Ft_t <- ZZ %*% PZ + HH + me_diag_t
        Ft_t <- (Ft_t + t(Ft_t)) * 0.5
        Fc_t <- tryCatch(chol(Ft_t), error = function(e) NULL)
        if (is.null(Fc_t)) return(.kf_fail("standard"))
        Fi_t <- chol2inv(Fc_t)
        ldf_t <- 2 * sum(log(diag(Fc_t)))
        ll_t <- ll_const - 0.5 * (ldf_t + drop(crossprod(v, Fi_t %*% v)))
        if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(.kf_fail("standard"))
        K_t   <- (TT %*% PZ + SS) %*% Fi_t
        s     <- as.numeric(TT %*% s) + drop(K_t %*% v)
        TmKZ  <- TT - K_t %*% ZZ; RmKD <- RR - K_t %*% DD
        P     <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD)
        ## Joseph true-noise term for me_extra (this branch only runs when
        ## me_extra[, t] has a nonzero entry): P' += K diag(me_extra[, t]) K'.
        ## me_variance stays F-only (regularizer convention; no Joseph term).
        P     <- P + K_t %*% (me_extra[, t] * t(K_t))
        P     <- (P + t(P)) * 0.5
        loglik <- loglik + ll_t
      } else {
        step <- .kf_step(s, P, v)
        if (is.null(step)) return(.kf_fail("standard"))
        loglik <- loglik + step$ll; s <- step$s
        if (!has_me_extra && !has_shock_scale && t > 1L && max(abs(step$P - P)) < ss_tol) {
          ss_reached   <- TRUE
          K_ss         <- step$K; F_inv_ss <- step$F_inv
          log_det_F_ss <- step$log_det_F
          ll_ss_const  <- ll_const - 0.5 * log_det_F_ss
          if (return_filtered) filtered[, t] <- s
          ## Bulk-process the remaining tail via dispatcher (C++ or R). Only
          ## safe when the tail has no missing observations: missing-data
          ## handling needs the per-step path. Keep the per-step fallback if
          ## the tail contains NA.
          tail_start <- t + 1L
          if (tail_start <= n_T &&
              !anyNA(Y_minus_d[, tail_start:n_T, drop = FALSE])) {
            out <- .kf_ss_dispatch(Y_minus_d, ZZ, TT, K_ss, F_inv_ss,
                                   ll_ss_const, s, tail_start, n_T, filtered)
            if (!out$ok) return(.kf_fail("standard"))
            loglik <- loglik + out$loglik
            s      <- out$s
            if (return_filtered) filtered <- out$filtered
            break
          }
        }
        P <- step$P
      }
    } else {
      ## Per-step SS (only reached if missing data forced re-entry).
      Fv   <- F_inv_ss %*% v
      ll_t <- ll_ss_const - 0.5 * sum(v * Fv)
      if (!is.finite(ll_t) || ll_t < .KF_LL_MIN)
        return(.kf_fail("standard"))
      loglik <- loglik + ll_t
      s <- as.numeric(TT %*% s + K_ss %*% v)
    }
    if (return_filtered) filtered[, t] <- s
  }

  if (return_filtered) rownames(filtered) <- endo[state_idx]
  list(loglik = loglik, filtered_states = filtered,
       n_obs = n_obs, n_T = n_T, method = "standard",
       lik_init = lik_init, d_diffuse = d_diffuse)
}
