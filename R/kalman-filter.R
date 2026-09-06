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
                                  Sigma_e = NULL, e_idx = NULL,
                                  det_rows = NULL, det_cols = NULL) {
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

  ## Deterministic known shocks: rows whose target state component has (at
  ## this period) no prior variance at all, so the ordinary sequential update
  ## is degenerate. See .kf_univariate_dispatch for why this is a mean shift.
  n_det <- length(det_rows)
  det_applied <- if (n_det) matrix(FALSE, n_det, n_T) else NULL
  is_det_row <- logical(n_obs)
  if (n_det) is_det_row[det_rows] <- TRUE
  ## Per-period count of observation components skipped because their forecast
  ## variance was (numerically) zero -- i.e. exactly predictable from what has
  ## already been processed. Reported in kalman_filter()'s diagnostics so a
  ## parity harness can see it without parsing a warning.
  n_skipped <- integer(n_T)

  loglik <- 0
  ok     <- TRUE
  diffuse <- length(P_inf) > 0 && max(abs(P_inf)) > 0
  diffuse_failed <- FALSE
  d_diffuse <- NA_integer_
  ll_contrib <- numeric(n_T)
  filtered <- if (return_filtered) matrix(0, n_state, n_T) else NULL

  for (t in seq_len(n_T)) {
    if (diffuse && t > max_diffuse) { diffuse_failed <- TRUE; break }

    ## Deterministic injections first: the period's observables load eps_t
    ## through D, so their innovations have to be taken against a state that
    ## already carries the shocks the caller supplied.
    for (k in seq_len(n_det)) {
      y_k <- Y_minus_d[det_rows[k], t]
      if (!is.finite(y_k)) next                   # shock unknown this period
      col <- det_cols[k]
      F_k <- P_star[col, col] + me_variance +
        if (has_me_extra_uni) me_extra[det_rows[k], t] else 0
      if (F_k <= kalman_tol) {
        ## P[col, ] is identically zero (a PSD matrix with a zero diagonal
        ## entry has a zero row), so nothing but this mean moves, and the row
        ## is skipped by the F_star guard below with a zero innovation.
        a[col] <- y_k
        det_applied[k, t] <- TRUE
      }
    }

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
      } else if (!is_det_row[i]) {
        ## Zero innovation variance: skip gracefully (singular F), and record
        ## it. A deterministic known-shock row is not counted -- it was
        ## applied above, not dropped.
        n_skipped[t] <- n_skipped[t] + 1L
      }
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
       ll_contrib = ll_contrib, det_applied = det_applied, P = P_star,
       n_skipped = n_skipped)
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
                                    shock_scale = NULL,
                                    known = NULL) {
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

  ## ---- Known shocks: one extra observation row per known shock ------------
  ## Row j is [0_{n_state}, e_j'], so it observes eps_j directly. The data row
  ## carries the known value and NA elsewhere, and the loop's existing
  ## `if (!is.finite(y_i)) next` handles the unknown periods -- the injection
  ## is only "on" in the periods the caller named.
  ##
  ## The measurement variance for these rows goes through me_extra, not the
  ## scalar me_variance, so that a HARD injection (variance 0) can coexist with
  ## me_variance > 0 on the real observables. That forces the R loop, which is
  ## correct: the C++ kernel takes a single scalar and could not express the
  ## two different noise levels.
  det_rows <- NULL; det_cols <- NULL; det_possible <- FALSE
  if (!is.null(known)) {
    n_k  <- length(known$idx)
    n_obs_real <- nrow(Y_minus_d)
    Zk   <- matrix(0, n_k, nb)
    for (i in seq_len(n_k)) Zk[i, n_state + known$idx[i]] <- 1
    Zb   <- rbind(Zb, Zk)
    me_full <- matrix(me_variance, nrow(Y_minus_d), ncol(Y_minus_d))
    if (!is.null(me_extra)) me_full <- me_full + me_extra
    me_extra <- rbind(me_full, known$sd^2)
    me_variance <- 0
    Y_minus_d <- rbind(Y_minus_d, unname(known$values))
    ## The injected rows are NA wherever the shock is unknown, so the frozen
    ## gain the lock assumes does not exist.
    ss_lock <- FALSE

    ## ---- DETERMINISTIC known shocks (zero prior variance) ---------------
    ## A shock with `stderr 0` is deterministic BEFORE conditioning, and the
    ## observation row above cannot express that: F = Z P Z' + H is exactly
    ## zero for it, so the sequential filter's `if (F_star > kalman_tol)`
    ## guard skipped the row and the supplied value never reached the state
    ## -- a silently ignored input, not a refused one.
    ##
    ## Zero prior variance does not mean "no information", it means the
    ## component IS its mean; conditioning on eps_j = v is then a MEAN SHIFT
    ## with no covariance update at all (the row and column of P are zero, so
    ## nothing else moves), and it carries no density -- a point mass has no
    ## dimension to integrate over. That is exactly the limit of the ordinary
    ## update as the variance goes to zero: the filtered path converges to the
    ## deterministic one, and the joint density's prior factor log p(eps = v)
    ## drops out. So a deterministic injection reports the SAME number under
    ## either semantics, and agrees with kalman_smoother() to machine
    ## precision rather than up to a prior term.
    ##
    ## `det_rows`/`det_cols` name the observation row and the augmented-state
    ## column of every candidate; the override is applied in the loop, where
    ## the prior variance for the period is known. It must run BEFORE the
    ## real observables of that period: y_t loads eps_t through D, so a
    ## trailing override would leave the period's innovations computed
    ## against a shock the caller had already told us.
    det_rows <- n_obs_real + seq_len(n_k)
    det_cols <- n_state + known$idx
    ## Is any candidate actually deterministic? The eps block is rebuilt from
    ## Sigma_e at every transition, so its prior variance at period t is
    ## Sigma_e[j,j] scaled by shock_scale[j,t]^2 -- known statically, which is
    ## what lets the C++ kernel (which has no override) be ruled out up front
    ## rather than discovered mid-recursion.
    sv <- diag(as.matrix(Sigma_e))[known$idx]
    sv_mat <- matrix(sv, n_k, ncol(known$values))
    if (!is.null(shock_scale))
      sv_mat <- sv_mat * shock_scale[known$idx, , drop = FALSE]^2
    det_possible <- any(!is.na(known$values) & known$sd == 0 &
                        sv_mat <= kalman_tol)
  }
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

  if (.HAS_RCPP_KALMAN_UNI() && !has_me_extra_d && !has_shock_scale_d &&
      !det_possible) {
    out <- kalman_univariate_loop_cpp(Y_minus_d, Zb, Tb, QQb, a0, Pb, Pi,
                                      me_variance, kalman_tol, diffuse_tol,
                                      conv_tol, as.integer(max_diffuse),
                                      .KF_LL_MIN, return_filtered,
                                      as.integer(n_state), ss_lock)
    out$a <- as.numeric(out$a)
    out$P_state <- out$P[s_idx, s_idx, drop = FALSE]
    out$P <- NULL
    out
  } else {
    out <- .kf_univariate_loop_R(Y_minus_d, Zb, Tb, QQb, a0, Pb, Pi,
                          me_variance, kalman_tol, diffuse_tol,
                          conv_tol, max_diffuse, .KF_LL_MIN,
                          return_filtered, n_state,
                          me_extra = me_extra,
                          shock_scale = shock_scale,
                          Sigma_e = Sigma_e, e_idx = e_idx,
                          det_rows = det_rows, det_cols = det_cols)
    out$P_state <- out$P[s_idx, s_idx, drop = FALSE]
    out$P <- NULL
    out
  }
}

### -- Deterministic shock MEANS ---------------------------------------------
###
### A different statement from `known_shocks`, and the difference is the whole
### point of having both:
###
###   known_shocks  eps_{j,t} = v          the REALISATION is known. The shock
###                                        stops being random: its variance is
###                                        used up, and (with a prior density)
###                                        the value informs the likelihood.
###   shock_means   E[eps_{j,t}] = m       the MEAN is known and the shock
###                                        keeps its variance. Nothing is
###                                        observed and nothing is estimated
###                                        about m; it is an INPUT.
###
### The second is IRIS's `filter(..., 'vary', J)`: a deterministic path added
### to the transition and measurement constants before the Kalman update. It
### is not expressible as an observation of eps, and forcing it through the
### augmented-state route would both change the covariance and add a density
### term that does not belong.
###
### It IS expressible as a data adjustment, because the system is linear.
### With eps_t = m_t + u_t, u_t ~ N(0, Sigma_e):
###   s_t = T s_{t-1} + R eps_t = s^det_t + s~_t,  s^det_t = T s^det_{t-1} + R m_t
###   y_t = Z s_{t-1} + D eps_t = y^det_t          + y~_t,  y^det_t = Z s^det_{t-1} + D m_t
### so subtracting y^det from the data leaves an ORDINARY filtering problem in
### s~ with no input at all, and s^det is added back to the reported states
### afterwards. No covariance changes (a mean shift moves no second moment),
### no likelihood term is added, and every method works unmodified -- no
### routing to the univariate filter, no new recursion, and zero-variance
### shocks and missing observations are non-events.
###
### TIMING. `shock_means[, t]` is the shock DATED t by default, exactly as
### `known_shocks` and `shock_scale` date their columns: it enters s_t and y_t.
### `shock_timing = "transition_next"` instead reads column t as driving the
### transition OUT of period t, so it enters s_{t+1} and y_{t+1} -- one column
### to the right, which is the adapter shift a caller comparing against a
### package with the other convention would otherwise apply by hand.
.kf_shock_means <- function(shock_means, shock_timing, shock_names, n_T,
                            what = "kalman_filter") {
  if (is.null(shock_means)) return(NULL)
  n_e <- length(shock_names)
  M <- as.matrix(shock_means)
  if (nrow(M) != n_e || ncol(M) != n_T)
    stop(sprintf("%s: `shock_means` must be n_exo x T (%d x %d); got %d x %d.",
                 what, n_e, n_T, nrow(M), ncol(M)), call. = FALSE)
  if (!is.null(rownames(M))) {
    if (!setequal(rownames(M), shock_names))
      stop(sprintf("%s: `shock_means` rownames do not match the model's shocks (%s).",
                   what, paste(shock_names, collapse = ", ")), call. = FALSE)
    M <- M[shock_names, , drop = FALSE]
  }
  if (any(!is.na(M) & !is.finite(M)))
    stop(what, ": `shock_means` entries must be finite or NA.", call. = FALSE)
  ## NA means "no shift here", i.e. zero -- the same shape a caller builds for
  ## `known_shocks`, where NA means "unknown". A mean of zero IS no shift, so
  ## the two readings coincide and neither surprises.
  M[is.na(M)] <- 0
  if (identical(shock_timing, "transition_next")) {
    ## Column t drives the transition out of t, so it lands on t+1. The last
    ## column would land outside the sample; it is dropped, which is what
    ## "the sample ends" means, not an error.
    M <- cbind(matrix(0, n_e, 1L), M[, -n_T, drop = FALSE])
  }
  if (all(M == 0)) return(NULL)
  M
}

### The deterministic trajectory implied by a mean path: y^det (subtracted from
### the data) and s^det (added back to the reported states). Column t of
### `s_det` is s^det_t; `s_det_prev` carries s^det_{t-1}, which is what the
### observation row and the one-step prediction load.
.kf_det_path <- function(M, TT, RR, ZZ, DD) {
  n_T <- ncol(M); n_s <- nrow(TT); n_o <- nrow(ZZ)
  s_det <- matrix(0, n_s, n_T)
  y_det <- matrix(0, n_o, n_T)
  s_prev <- numeric(n_s)
  for (t in seq_len(n_T)) {
    y_det[, t] <- drop(ZZ %*% s_prev) + drop(DD %*% M[, t])
    s_prev     <- drop(TT %*% s_prev) + drop(RR %*% M[, t])
    s_det[, t] <- s_prev
  }
  list(s_det = s_det, y_det = y_det)
}


### -- Routing table for $diagnostics ----------------------------------------
###
### Built on EVERY filter call, including the per-draw ones inside an MCMC
### chain, so it goes through structure() rather than data.frame(): the latter
### measured 0.065 ms against a 1.92 ms filter call on the rbc fixture (3.4%
### of a small-model likelihood evaluation, which is exactly the kind of R
### overhead the per-draw path is bound by), the former ~10x less.
.kf_routing_df <- function(rows) {
  n <- length(rows)
  if (!n) rows <- list()
  structure(list(from   = vapply(rows, `[[`, "", "from"),
                 to     = vapply(rows, `[[`, "", "to"),
                 reason = vapply(rows, `[[`, "", "reason")),
            class = "data.frame",
            row.names = if (n) seq_len(n) else integer(0))
}


### -- Known historical shocks ------------------------------------------------
###
### `known_shocks` is n_exo x T: the VALUE where a shock is known, NA where it
### is not. NA-as-unknown is the same convention `data` already uses, and the
### n_exo x T shape is the one `shock_scale` already uses, so a caller who has
### met either has met this.
###
### The mechanism costs almost nothing because of an accident of the univariate
### filter: it runs on the AUGMENTED state x_t = [s_{t-1}; eps_t] (see
### .kf_univariate_dispatch), so the shocks ARE state components there. A known
### shock is then an exact observation of a component that already exists -- an
### extra observation row Z = [0, e_j'] with zero measurement error -- and the
### periods where it is unknown are NA rows, which the loop already skips one
### observable at a time. No new recursion, and no change to the ones that are
### there.
###
### `known_shocks_sd` makes the injection SOFT: the value is then an observation
### with that standard deviation rather than an exact constraint, which is what
### a judgemental adjustment ("about 0.7, but I would not die for it") actually
### is. NA or 0 means exact.
.kf_known_shocks <- function(known_shocks, known_shocks_sd, shock_names, n_T,
                             what = "kalman_filter") {
  if (is.null(known_shocks)) return(NULL)
  n_e <- length(shock_names)
  K <- as.matrix(known_shocks)
  if (nrow(K) != n_e || ncol(K) != n_T)
    stop(sprintf("%s: `known_shocks` must be n_exo x T (%d x %d); got %d x %d.",
                 what, n_e, n_T, nrow(K), ncol(K)), call. = FALSE)
  if (!is.null(rownames(K))) {
    if (!setequal(rownames(K), shock_names))
      stop(sprintf("%s: `known_shocks` rownames do not match the model's shocks (%s).",
                   what, paste(shock_names, collapse = ", ")), call. = FALSE)
    K <- K[shock_names, , drop = FALSE]
  }
  if (any(is.finite(K) & !is.finite(K)) || any(!is.na(K) & !is.finite(K)))
    stop(what, ": `known_shocks` entries must be finite or NA.", call. = FALSE)

  S <- NULL
  if (!is.null(known_shocks_sd)) {
    S <- as.matrix(known_shocks_sd)
    if (length(S) == 1L) S <- matrix(as.numeric(S), n_e, n_T)
    if (nrow(S) != n_e || ncol(S) != n_T)
      stop(sprintf("%s: `known_shocks_sd` must be n_exo x T (%d x %d) or a scalar; got %d x %d.",
                   what, n_e, n_T, nrow(S), ncol(S)), call. = FALSE)
    if (!is.null(rownames(S)) && setequal(rownames(S), shock_names))
      S <- S[shock_names, , drop = FALSE]
    if (any(!is.na(S) & (!is.finite(S) | S < 0)))
      stop(what, ": `known_shocks_sd` must be non-negative and finite (or NA).",
           call. = FALSE)
    S[is.na(S)] <- 0
  }
  rows <- which(rowSums(!is.na(K)) > 0L)     # only shocks that are ever known
  if (!length(rows)) return(NULL)
  list(idx = rows, values = K[rows, , drop = FALSE],
       sd = if (is.null(S)) matrix(0, length(rows), n_T) else S[rows, , drop = FALSE],
       names = shock_names[rows])
}


### -- User-supplied initial condition (a0 / P0) ------------------------------
###
### Both entry points hard-coded a zero-mean state and a lik_init-derived P0,
### so there was no way to start the recursion anywhere else -- no way to carry
### a state across a sample split, to condition on a known history, or to hand
### the smoother a prior of your own. The internal machinery always took them
### (.kf_univariate_dispatch(s0, P_state, P_inf_state), .kf_diffuse_phase(s0,
### ...)); only the public shape was missing.
###
### `a0` is in DEVIATIONS from the steady state, like every state quantity in
### this package -- the observables are in levels and the filter subtracts
### `dr$ys[obs_vars]`, but the STATE vector it reports is deviations. Passing a
### level here is the same class of error the smoother's missing intercept was.
###
### Matched BY NAME whenever `a0` is named: a same-length vector in a different
### order is a reorder request, and silently accepting it by position would
### load each state into the wrong equation while every dimension still checked
### out (the rule .fbt_data() already applies to observables).
.kf_init_mean <- function(a0, state_names, n_state, what = "kalman_filter") {
  if (is.null(a0)) return(numeric(n_state))
  a0 <- unlist(a0, use.names = TRUE)
  if (!is.numeric(a0) || anyNA(a0) || any(!is.finite(a0)))
    stop(what, ": `a0` must be finite and numeric.", call. = FALSE)
  if (!is.null(names(a0)) && !is.null(state_names)) {
    miss <- setdiff(state_names, names(a0))
    if (length(miss))
      stop(sprintf("%s: `a0` is named but does not cover every state: %s. The states are: %s.",
                   what, paste(miss, collapse = ", "),
                   paste(state_names, collapse = ", ")), call. = FALSE)
    extra <- setdiff(names(a0), state_names)
    if (length(extra))
      stop(sprintf("%s: `a0` names entries that are not states: %s.",
                   what, paste(extra, collapse = ", ")), call. = FALSE)
    return(as.numeric(a0[state_names]))
  }
  if (length(a0) != n_state)
    stop(sprintf(paste0("%s: `a0` has %d entr%s but the model has %d state(s)",
                        "%s. Name it to be matched by name."),
                 what, length(a0), if (length(a0) == 1L) "y" else "ies", n_state,
                 if (is.null(state_names)) ""
                 else paste0(" (", paste(state_names, collapse = ", "), ")")),
         call. = FALSE)
  as.numeric(a0)
}

### P0 validation. Symmetrised on the way through (round-off in a
### user-constructed covariance is expected); a genuinely asymmetric or
### negative-definite matrix is refused rather than silently repaired.
.kf_init_cov <- function(P0, state_names, n_state, what = "kalman_filter") {
  if (is.null(P0)) return(NULL)
  if (length(P0) == 1L && is.numeric(P0)) P0 <- diag(as.numeric(P0), n_state)
  P0 <- as.matrix(P0)
  if (!is.numeric(P0) || nrow(P0) != n_state || ncol(P0) != n_state)
    stop(sprintf("%s: `P0` must be %d x %d (or a scalar for a multiple of the identity); got %d x %d.",
                 what, n_state, n_state, nrow(P0), ncol(P0)), call. = FALSE)
  if (anyNA(P0) || any(!is.finite(P0)))
    stop(what, ": `P0` must be finite.", call. = FALSE)
  if (!is.null(rownames(P0)) && !is.null(state_names)) {
    if (!setequal(rownames(P0), state_names))
      stop(what, ": `P0` has dimnames that do not match the state names.",
           call. = FALSE)
    P0 <- P0[state_names, state_names, drop = FALSE]
  }
  asym <- max(abs(P0 - t(P0)))
  if (asym > 1e-8 * max(1, max(abs(P0))))
    stop(sprintf("%s: `P0` is not symmetric (max |P0 - t(P0)| = %.3g).",
                 what, asym), call. = FALSE)
  P0 <- (P0 + t(P0)) * 0.5
  ev <- tryCatch(min(eigen(P0, symmetric = TRUE, only.values = TRUE)$values),
                 error = function(e) NA_real_)
  if (is.finite(ev) && ev < -1e-8 * max(1, max(abs(P0))))
    stop(sprintf("%s: `P0` is not positive semi-definite (smallest eigenvalue %.3g).",
                 what, ev), call. = FALSE)
  dimnames(P0) <- if (is.null(state_names)) NULL else list(state_names, state_names)
  P0
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
#' @param data observation matrix (\code{n_obs} x \code{T}).
#' @param dr decision rule (output of \code{\link{solve_perturbation}}).
#' @param model compiled model object (output of \code{\link{compile_model}}).
#' @param params named numeric vector of parameter values.
#' @param obs_vars character vector of observed variable names.
#' @param return_filtered logical; if \code{TRUE}, return filtered state estimates
#'   (one column per time step).
#' @param ss_tol tolerance for steady-state lock detection (default \code{.LYAP_TOL}).
#' @param me_variance scalar variance of iid Gaussian measurement error on
#'   every observable (default \code{0}). It enters the innovation covariance
#'   \code{F} AND the state-covariance (Joseph) update on every method --
#'   the exact likelihood of the noise-augmented model. See Details.
#' @param return_ll_contrib logical; if \code{TRUE}, return per-step
#'   log-likelihood contributions (prediction-error decomposition).
#' @param method character; filtering algorithm: \code{"auto"} (default; picks
#'   \code{"standard"} or \code{"chandrasekhar"}), \code{"dare"} (textbook
#'   Kalman filter, no steady-state shortcut), \code{"chandrasekhar"}
#'   (Morf--Sidhu--Kailath low-rank increment recursion; requires
#'   \code{lik_init = "stationary"}, a fully observed panel and no
#'   \code{me_extra} / \code{shock_scale}, and errors otherwise --
#'   \code{"auto"} selects it only when \code{n_state > 100}, the measured
#'   crossover against the C++ \code{"standard"} loop),
#'   \code{"standard"} (per-step Riccati with steady-state
#'   lock), \code{"reference"} (alias for \code{"dare"}), or
#'   \code{"univariate"} (Koopman--Durbin 2000 sequential filter on the
#'   augmented state \code{[s; eps]}; handles singular innovation
#'   covariances and the exact-diffuse phase with singular \code{F_inf}).
#' @param lik_init character; initialization of the state covariance \code{P0}
#'   (and, for unit-root models, the diffuse covariance): \code{"auto"}
#'   (default), \code{"stationary"}, \code{"diffuse"}, or \code{"kappa"}.
#'   See Details. Supplying \code{P0} overrides this and is reported back as
#'   \code{lik_init = "user"}.
#' @param a0 Initial state mean, length \code{n_state}, in \strong{deviations
#'   from the steady state} -- the convention the filter's own
#'   \code{filtered_states} are in. \code{NULL} (default) starts at the steady
#'   state, i.e. a vector of zeros. Matched BY NAME when named (a same-length
#'   vector in a different order is a reorder request, not a relabelling); the
#'   state names are the \code{rownames} of \code{filtered_states}. Composes
#'   with every \code{lik_init}, including \code{"diffuse"}.
#' @param P0 Initial state covariance, \code{n_state x n_state} (or a scalar
#'   for a multiple of the identity). \code{NULL} (default) takes the
#'   covariance implied by \code{lik_init}. Must be symmetric and positive
#'   semi-definite; matched by \code{dimnames} when it has them. Mutually
#'   exclusive with \code{lik_init = "diffuse"}, which builds its own
#'   \code{(P_inf, P_star)} split. A supplied \code{P0} is not the Lyapunov
#'   prior the Chandrasekhar increment recursion is initialised from, so
#'   \code{method = "auto"} routes to \code{"standard"} and explicit
#'   \code{method = "chandrasekhar"} is refused.
#' @param me_extra \code{n_obs x T} matrix of additional per-observable,
#'   per-period measurement-error variances (default \code{NULL}, no extra
#'   variance).  Like the scalar \code{me_variance}, \code{me_extra} is TRUE
#'   per-period measurement noise on all paths: it enters the innovation
#'   covariance \code{F_t} AND the Joseph-form state-covariance update
#'   (\code{P += K_t \%*\% (diag(me_variance + me_extra[, t]) \%*\% t(K_t))}), so the
#'   \code{"standard"}, \code{"dare"}/\code{"reference"}, and
#'   \code{"univariate"} methods agree exactly under any combination of
#'   \code{me_variance} and \code{me_extra}.
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
#'   exact-diffuse phase with missing observations (including NA rows
#'   introduced by \code{me_extra}) runs on the sequential univariate filter,
#'   which handles it exactly -- see \code{lik_init} in Details.
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
#' @param obs_aggregation Optional named list declaring one or more
#'   observables as TEMPORAL AGGREGATES of a higher-frequency model variable,
#'   e.g. \code{list(gdp_q = list(of = "gdp_m", type = "flow_sum", k = 3L))}
#'   ("the observable \code{gdp_q} is the 3-period sum of the model variable
#'   \code{gdp_m}").  \code{NULL} (default) falls back to
#'   \code{model$obs_aggregation}, and when that is also \code{NULL} the
#'   filter is byte-identical to the non-aggregated one.  Names must appear in
#'   \code{obs_vars}; \code{of} must be a model variable; \code{type} is one of
#'   \code{"flow_sum"}, \code{"flow_mean"}, \code{"stock_end"} or
#'   \code{"triangle"} (see \code{\link{mf_aggregation_weights}}).
#'   Implemented by fixed-weight state augmentation (Harvey 1989;
#'   Mariano--Murasawa 2003): the aggregator's \code{m-1} lags become extra
#'   states, so \code{ZZ} and \code{TT} stay constant.  Supply the data at the
#'   HIGH frequency with the aggregate observed at the last period of each
#'   window and \code{NA} in between (\code{\link{mf_expand_observations}}
#'   builds such a column); the existing missing-data path handles the gaps.
#'   The steady-state offset of an aggregated row is scaled by the sum of its
#'   weights, so the data must be on the aggregate's own scale.
#'   Consumed by the Gaussian Kalman likelihood only: the other likelihoods in
#'   \code{\link{make_log_posterior}} (particle filters, pruned/PSKF,
#'   \code{whittle}, \code{cumulant}, \code{student_t}) and the
#'   Markov-switching Kim filter reject it, and the ANALYTIC score
#'   (\code{make_posterior_grad}) has no aggregation awareness -- sample a
#'   mixed-frequency posterior with a gradient-free sampler.
#' @param known_shocks Known historical shock values: an \code{n_exo x T}
#'   matrix carrying the value where a shock is known and \code{NA} where it is
#'   not -- the \code{NA}-as-unknown convention \code{data} uses, and the
#'   \code{n_exo x T} shape \code{shock_scale} uses. Rows are matched BY NAME
#'   when the matrix has rownames. \code{NULL} (default), or a matrix that is
#'   all \code{NA}, is a no-op.
#'
#'   Use it for a shock you actually know: an announced policy change, a
#'   measured intervention, a judgemental adjustment carried over from another
#'   exercise.
#'
#'   \strong{Semantics.} The known value is treated as an OBSERVATION, so the
#'   reported log-likelihood is the \strong{joint} \eqn{\log p(y, \varepsilon =
#'   v)} -- the known shock is data, and it informs that shock's own standard
#'   error. It therefore differs from the conditional \eqn{\log p(y \mid
#'   \varepsilon = v)} by exactly \eqn{\log p(\varepsilon = v)}; subtract that
#'   term if the conditional is what you want. \code{\link{kalman_smoother}}
#'   conditions, so its \code{loglik} is the conditional one.
#'
#'   \strong{Deterministic shocks.} When the injected shock has \emph{no}
#'   prior variance (a zero diagonal of \code{Sigma_e}, or a
#'   \code{shock_scale} of zero at that period) the value is applied as a
#'   deterministic INPUT: the state mean shifts and the covariance does not
#'   move, because a component with zero variance IS its mean. A point mass
#'   carries no density, so the two semantics coincide there -- the reported
#'   log-likelihood is both the joint and the conditional one, and it agrees
#'   with \code{\link{kalman_smoother}} exactly rather than up to a prior
#'   term. It is also the limit of the positive-variance case: as the variance
#'   goes to zero the filtered path converges to the deterministic one.
#'   (Before 0.9.3.2 such a value was silently ignored -- the injected row's
#'   forecast variance was exactly zero, so the sequential filter skipped it.)
#'
#'   Only the univariate (Koopman--Durbin) filter can express this: it runs on
#'   the augmented state \eqn{[s_{t-1}; \varepsilon_t]}, where the shocks are
#'   state components and a known shock is an exact observation of one. Passing
#'   \code{known_shocks} therefore routes to \code{method = "univariate"}.
#' @param known_shocks_sd Optional standard errors for \code{known_shocks},
#'   same shape (or a scalar). The injection is then SOFT -- the value is an
#'   observation with that standard deviation rather than an exact constraint,
#'   which is what a judgemental adjustment usually is. \code{NULL} or
#'   \code{NA} means exact.
#' @param shock_means Deterministic shock MEANS: an \code{n_exo x T} matrix of
#'   mean shifts, \code{NULL} (default) for none. \code{NA} and \code{0} both
#'   mean "no shift here", so the same matrix shape used for
#'   \code{known_shocks} works. Rows are matched BY NAME when the matrix has
#'   rownames.
#'
#'   \strong{This is a different statement from \code{known_shocks}, and the
#'   difference is the point of having both.} \code{known_shocks} says the
#'   REALISATION is known, \eqn{\varepsilon_{j,t} = v}: the shock stops being
#'   random, its variance is used up, and (where it has a prior density) the
#'   value enters the likelihood. \code{shock_means} says the MEAN is known,
#'   \eqn{E[\varepsilon_{j,t}] = m}, and the shock \strong{keeps its
#'   variance}: nothing is observed, nothing about \eqn{m} is estimated, and
#'   the likelihood gains no term. It is a deterministic INPUT -- the analogue
#'   of a \code{vary} / plan path in other state-space packages, entering the
#'   transition and measurement constants as \eqn{R m_t} and \eqn{D m_t}
#'   before the update.
#'
#'   Because the system is linear the two coincide exactly for a shock with no
#'   prior variance (a \code{stderr 0} shock, or a \code{shock_scale} of zero
#'   at that period): knowing the mean and knowing the realisation are then the
#'   same statement, and both entry points return the same states and the same
#'   log-likelihood. They are pinned against each other by test.
#'
#'   No method routing is involved: the mean path splits off as a deterministic
#'   trajectory that is subtracted from the data and added back to the reported
#'   states, so every \code{method} evaluates it identically, and zero-variance
#'   shocks, missing observations and \code{a0}/\code{P0} are non-events.
#' @param shock_timing How to read the columns of \code{shock_means}:
#'   \describe{
#'     \item{\code{"dated"}}{(default) column \eqn{t} is the shock DATED
#'       \eqn{t}. It enters \eqn{s_t} and \eqn{y_t}, the same dating
#'       \code{known_shocks} and \code{shock_scale} use for their columns.}
#'     \item{\code{"transition_next"}}{column \eqn{t} drives the transition
#'       OUT of period \eqn{t}, so it enters \eqn{s_{t+1}} and
#'       \eqn{y_{t+1}} -- the whole matrix shifted one column right, with the
#'       last column falling outside the sample. This is the convention to pass
#'       when adapting a path written for a package that dates its input by the
#'       transition rather than by the shock.}
#'   }
#'   Ignored when \code{shock_means} is \code{NULL}.
#' @param me_floor_check Logical: when \code{me_variance > 0}, compare it
#'   against the smallest eigenvalue of the model-implied (ME-free) steady-
#'   state innovation covariance \code{F} and warn when the assumed
#'   measurement-noise variance is large relative to that eigenvalue -- i.e.
#'   when the observation noise, not the model, dominates some near-collinear
#'   combination of observables (see
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
#' multivariate path triggers a \emph{conditional} automatic fallback to the
#' univariate (sequential) filter, which processes observables one at a time
#' and skips zero-variance components instead of inverting \code{F} (the
#' analog of Dynare's
#' \code{univariate_kalman_filter_if_singularity_is_detected}); see
#' \strong{Singularity-fallback contract} below.
#' Set \code{me_variance > 0} to deliberately add measurement noise, e.g.
#' for stochastically singular models where the singularity is a modelling
#' choice rather than a numerical artifact. When comparing marginal
#' likelihoods (model comparison) across models, ensure all use the same
#' \code{me_variance} setting; the ranking is otherwise invalidated.
#' \code{me_variance > 0} is a genuine noise model on EVERY method: the
#' observation equation becomes \code{y_t = ZZ x_t + DD e_t + u_t} with
#' \code{u_t ~ N(0, me_variance I)}, so the variance enters both \code{F}
#' and the state-covariance update
#' (\code{P += K \%*\% (me_variance * t(K))}). All methods --
#' \code{"standard"}, \code{"dare"}, \code{"chandrasekhar"},
#' \code{"univariate"} -- and \code{\link{kalman_smoother}} therefore
#' evaluate the SAME likelihood at any \code{me_variance}, and it equals the
#' brute-force joint-Gaussian projection of the whole sample. (Before dynhr 0.9.2.x the multivariate
#' methods added \code{me_variance} to \code{F} only -- a regulariser --
#' which understated the state uncertainty and biased the log-likelihood
#' UPWARD by \code{O(me_variance)} per period against the package's own
#' data-generating process; \code{me_variance = 0} results are unchanged.)
#'
#' \strong{Singularity-fallback contract:}
#' When a multivariate method (\code{"standard"}, \code{"dare"},
#' \code{"chandrasekhar"}) hits a singular / non-positive-definite innovation
#' covariance \code{F}, a failed DARE/Chandrasekhar step, or a non-finite
#' per-period contribution, the filter retries on the \code{"univariate"}
#' filter \strong{only when the two paths evaluate the SAME likelihood},
#' namely when all of
#' \itemize{
#'   \item \code{me_variance == 0} (a conservative gate retained from the
#'     pre-F3-D regulariser convention: the two paths now evaluate the same
#'     likelihood at any \code{me_variance}, but the gate still refuses to
#'     switch estimator mid-chain),
#'   \item \code{me_extra} is \code{NULL} or all zero,
#'   \item the resolved \code{lik_init} is not \code{"diffuse"} (the
#'     univariate diffuse filter keeps the divergent \code{0.5*log(F_inf)}
#'     term that the multivariate exact-diffuse recursion renormalises away;
#'     the two conventions can differ by more than 15 nats on unit-root
#'     models)
#' }
#' hold.  In that case a \code{warning()} is emitted \strong{once per
#' \code{kalman_filter()} call} naming the method that failed, and the
#' returned \code{$method} is \code{"univariate"}.
#' \strong{Otherwise the draw is rejected: the function returns
#' \code{loglik = -Inf} and \code{$method} names the failed multivariate
#' method}, which is what Dynare does without
#' \code{univariate_kalman_filter_if_singularity_is_detected}.  The point is
#' that an MCMC chain in which only \emph{some} draws trip the fallback must
#' not silently sample a mixture of two different likelihood definitions.
#' (One routing decision is deliberately \emph{not} a fallback and is not
#' covered by this contract: with \code{lik_init = "diffuse"} and an
#' \code{F_inf} that is singular but nonzero -- "Case C" -- the multivariate
#' exact-diffuse recursion is undefined and the univariate diffuse filter is
#' the algorithm for that model, not a substitute for it.)
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
#' @return A list with
#'   \describe{
#'     \item{\code{loglik}}{the log-likelihood.}
#'     \item{\code{updated_states}, \code{predicted_states}}{the two state
#'       paths, \code{n_state x T}, with the timing in the names:
#'       \code{updated_states[, t]} is \eqn{s_{t|t}} (conditioned on
#'       \eqn{y_{1:t}}) and \code{predicted_states[, t]} is \eqn{s_{t|t-1}}
#'       (conditioned on \eqn{y_{1:t-1}}, with \eqn{s_{1|0}} taken from
#'       \code{a0}). Both only when \code{return_filtered = TRUE}.}
#'     \item{\code{filtered_states}}{the SAME matrix as
#'       \code{updated_states}, under the name the rest of the package uses.}
#'     \item{\code{loglik_contrib}}{per-period contributions, when asked for.}
#'     \item{\code{final_state}, \code{final_cov}}{the state hand-off:
#'       \eqn{s_{T|T}} and \eqn{Var(s_T \mid y_{1:T})}, named, in the same
#'       convention \code{a0} / \code{P0} take. \code{final_cov} is
#'       \code{NULL} for \code{method = "chandrasekhar"}, which never forms
#'       \eqn{P}, and after a failed run.}
#'     \item{\code{method}, \code{lik_init}, \code{d_diffuse},
#'       \code{n_obs}, \code{n_T}}{what ran.}
#'     \item{\code{diagnostics}}{machine-readable record of the run:
#'       \code{method_requested} / \code{method_used},
#'       \code{lik_init_requested} / \code{lik_init_used}, \code{routing}
#'       (a data frame of \code{from} / \code{to} / \code{reason} rows, one
#'       per automatic reroute or fallback), \code{diffuse_periods},
#'       \code{missing_by_period} / \code{n_missing} (observations that were
#'       absent), \code{dropped_by_period} / \code{n_dropped} (observations
#'       that were present but exactly predictable, hence carried no
#'       information), \code{known_shocks} (how many injected cells, and how
#'       many of them were deterministic) and \code{loglik_type}
#'       (\code{"marginal"}, \code{"joint"} or \code{"conditional"}).
#'       \code{\link{kalman_smoother}} returns the same fields.}
#'   }
#'
#' @section State timing:
#' The state space is written in LAGGED form -- \eqn{s_t = T s_{t-1} +
#' R \varepsilon_t} and \eqn{y_t = Z s_{t-1} + D \varepsilon_t} -- so
#' \eqn{s_t} is the model's state variables dated \eqn{t}, and
#' \code{updated_states[, t]} is those variables conditioned on data through
#' \eqn{t}. \code{predicted_states} is the one-step-ahead pair,
#' \eqn{s_{t|t-1}}; it needs nothing extra from the recursion, because
#' \eqn{E[\varepsilon_t] = 0} in the deviation system makes
#' \eqn{s_{t|t-1} = T s_{t-1|t-1}} exactly (plus \eqn{R m_t} when a
#' \code{shock_means} path is supplied).
#'
#' A unit pulse in \code{shock_means} with no data to update on therefore
#' traces the model's own impulse response exactly: with
#' \code{shock_timing = "dated"} and a pulse in column 1,
#' \code{updated_states[, t]} is \eqn{T^{t-1} R e_j}, and the two paths
#' coincide because there is nothing to update with. That identity is the
#' bridge to a package that reports shock responses on a different grid --
#' compare against it once, choose the \code{shock_timing} that lines up, and
#' the convention is then stated in the call rather than applied by hand.
#'
#' @section Matching another package's shock timing:
#' Two packages can inject the same deterministic shock path and produce
#' responses one period apart, purely because they index the input
#' differently. \code{shock_timing} exists so that difference is declared in
#' the call instead of applied by hand to the output.
#'
#' \strong{For IRIS, the answer is the default and no shift is needed.}
#' Verified against IRIS Toolbox Release 20180308 under Octave: with
#' \code{shock_means} and \code{shock_timing = "dated"}, dynhr reproduces
#' \code{filter(m, d, range, 'vary=', j)} to machine precision -- on all three
#' of IRIS's outputs and on the smoothed shocks -- while
#' \code{"transition_next"} is measurably wrong for it. The names line up
#' one-for-one:
#'
#' \tabular{ll}{
#'   IRIS \code{'output=', 'predict'} \tab \code{predicted_states} \cr
#'   IRIS \code{'output=', 'filter'}  \tab \code{updated_states}   \cr
#'   IRIS \code{'output=', 'smooth'}  \tab \code{smoothed_states} (from
#'     \code{\link{kalman_smoother}}) \cr
#'   IRIS \code{'vary='} shock tunes  \tab \code{shock_means}
#' }
#'
#' If you previously needed a one-period shift to line the two up, the cause
#' was the MECHANISM rather than the timing: \code{known_shocks} is an exact
#' observation of the shock, whereas an IRIS \code{vary} tune sets the mean
#' and leaves the shock random -- so its smoothed shock is revised away from
#' the injected value, and no shift of an exactly-pinned path can reproduce
#' that. \code{shock_means} is the matching statement; use it and drop the
#' adapter.
#'
#' For any other reference implementation, three lines settle it, with no data
#' needed:
#'
#' \preformatted{
#' Yna <- Y; Yna[] <- NA                  # same shape, no observations
#' M   <- matrix(0, n_exo, T, dimnames = list(exo_names, NULL))
#' M["e_j", 1] <- 1                       # unit pulse in the FIRST column
#' kalman_filter(Yna, dr, model, params, obs_vars = obs,
#'               return_filtered = TRUE, shock_means = M)$updated_states[, 1:4]
#' }
#'
#' With no observations there is nothing to update with, so the filtered path
#' IS the deterministic path, and the deterministic path of a unit pulse is
#' the model's impulse response: under \code{"dated"} that is
#' \eqn{T^{t-1} R e_j}, with the IMPACT in column 1 -- the same column as the
#' pulse -- and \code{predicted_states} identical to it. Produce the other
#' package's response to the same pulse and compare: impact in the same period
#' means keep \code{"dated"}; impact one period later means
#' \code{"transition_next"}, which is exactly the same as shifting the matrix
#' one column right yourself.
#'
#' Do not infer the answer from a HISTORICAL filtered path, where an offset in
#' the shock input and an offset in the state output look alike. The pulse
#' experiment fixes the input convention with no data in play, and
#' \code{updated_states} / \code{predicted_states} then fix the output
#' convention by name (\eqn{s_{t|t}} against \eqn{s_{t|t-1}}).
#'
#' @section Filtering a split sample:
#' \code{final_state} / \code{final_cov} are exactly what \code{a0} /
#' \code{P0} take, so a sample can be filtered in two calls and the
#' prediction-error decomposition holds to numerical tolerance:
#'
#' \preformatted{
#' f1 <- kalman_filter(y[1:k, ], dr, model, params, obs_vars = obs)
#' f2 <- kalman_filter(y[(k+1):T, ], dr, model, params, obs_vars = obs,
#'                     a0 = f1$final_state, P0 = f1$final_cov)
#' f1$loglik + f2$loglik   # == the unsplit loglik
#' }
#'
#' Both are in DEVIATIONS from the steady state (the data is in levels and the
#' filter subtracts \code{dr$ys[obs_vars]} itself), and \code{a0} is
#' \eqn{s_0} -- the state one period BEFORE the first row of \code{data}, so
#' the hand-off lines up without an off-by-one. Supplying \code{P0} sets
#' \code{lik_init = "user"}, which is reported back in \code{$lik_init}.
#' For latent history BEFORE the first observation, use
#' \code{kalman_smoother(pre_sample = k)} instead: that estimates the padded
#' periods from the data that follows them, which \code{a0} cannot do.
#'
#' @seealso \code{\link{kalman_smoother}}, \code{\link{make_posterior}},
#'   \code{\link{kf_innovation_diagnostics}}
#' @examples
#' model    <- parse_mod(system.file("extdata/models/rbc.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, verbose = FALSE)
#' steady   <- solve_steady(compiled, model$param_values,
#'                          endo_names = model$var_names,
#'                          exo_names  = model$varexo_names, verbose = FALSE)
#' dr <- solve_perturbation(model, compiled, steady$values,
#'                          model$param_values, verbose = FALSE)
#'
#' ## Simulated data standing in for observations. rbc.mod has ONE shock, so
#' ## one observable keeps the innovation covariance non-singular.
#' set.seed(1)
#' paths <- simulate_model(dr, n_periods = 100L, model = model, burn_in = 20L)
#' Y <- as.matrix(paths[, "y", drop = FALSE])
#'
#' kf <- kalman_filter(Y, dr, model, model$param_values, obs_vars = "y",
#'                     me_variance = 1e-6)
#' kf$loglik
#' @export
kalman_filter <- function(data, dr, model, params, obs_vars,
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
                          obs_aggregation = NULL,
                          a0 = NULL,
                          P0 = NULL,
                          known_shocks = NULL,
                          known_shocks_sd = NULL,
                          shock_means = NULL,
                          shock_timing = c("dated", "transition_next"),
                          me_floor_check = getOption("dynhr.me_floor_check",
                                                     TRUE)) {
  shock_timing <- match.arg(shock_timing)

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

  ## ---- Mixed-frequency / temporal aggregation ----------------------------
  ## Resolve the aggregation spec BEFORE the observation rows are cut: an
  ## aggregated observable's NAME (e.g. "gdp_q") is not a model variable, the
  ## higher-frequency variable it aggregates (`of`, e.g. "gdp_m") is. With no
  ## spec (the corpus case) .mf_resolve() returns NULL and every line below is
  ## the pre-aggregation one, byte for byte -- including the `obs_vars` used
  ## as the row selector. See R/mixed-frequency.R for the algebra.
  mf <- .mf_resolve(obs_aggregation %||% model$obs_aggregation, obs_vars,
                    known = endo)
  obs_base <- if (is.null(mf)) obs_vars else mf$base

  obs_idx <- match(obs_base, endo)
  if (any(is.na(obs_idx)))
    stop("Observed variables not found: ", paste(obs_base[is.na(obs_idx)], collapse = ", "))

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d   <- dr$ys[obs_base]
  state_names_out <- endo[state_idx]

  if (!is.null(mf)) {
    ## Fixed-weight state augmentation: carry the m-1 lags each aggregator
    ## needs as extra states, so ZZ and TT stay CONSTANT and the recursions
    ## (and both C++ kernels) are untouched. The aggregate is observed only
    ## every k-th period; the NA periods take the existing missing-data path.
    aug <- .mf_augment_matrices(TT, RR, ZZ, DD, mf$w_list, obs_vars)
    TT <- aug$TT; RR <- aug$RR; ZZ <- aug$ZZ; DD <- aug$DD
    state_names_out <- c(state_names_out, aug$aug_names)
    n_state <- n_state + aug$n_aug
    ## An aggregate's steady state is sum(w) times the underlying variable's
    ## (3 * ys for a 3-period flow sum, 1 * ys for a mean or an end-of-period
    ## stock), so the demeaning offset has to be rescaled with it.
    d <- d * mf$scale
  }

  ## -- User-supplied initial condition -----------------------------------
  ## Validated here, once TT/n_state/state_names_out are final (the
  ## mixed-frequency augmentation above adds states, so a0/P0 must cover them
  ## too -- name them and the match is checked for you).
  a0_vec <- .kf_init_mean(a0, state_names_out, n_state)
  P0_mat <- .kf_init_cov(P0, state_names_out, n_state)
  if (!is.null(P0_mat)) {
    if (lik_init_orig == "diffuse")
      stop("kalman_filter: `P0` and lik_init = \"diffuse\" are two different ",
           "initialisations -- the exact-diffuse recursion builds its own ",
           "(P_inf, P_star) split and has nothing to do with a supplied P0. ",
           "Pass one or the other. `a0` composes with either.", call. = FALSE)
    ## "user" is an INTERNAL lik_init value: it exists so the two
    ## `if (lik_init == "stationary")` sites below cannot overwrite P0 with the
    ## Lyapunov solution, and so that $lik_init reports what was actually run.
    lik_init <- "user"
  }

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

  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  if (nrow(data) != n_obs) data <- t(data)
  n_T <- ncol(data)
  known <- .kf_known_shocks(known_shocks, known_shocks_sd, exo, n_T)
  mean_path <- .kf_shock_means(shock_means, shock_timing, exo, n_T)

  ## Two deterministic statements about the same shock in the same period are
  ## a contradiction, not a composition: `known_shocks` fixes the realisation,
  ## so a mean for it is already spoken for. Different shocks, or different
  ## periods, compose fine and are left alone.
  if (!is.null(mean_path) && !is.null(known)) {
    clash <- !is.na(known$values) & mean_path[known$idx, , drop = FALSE] != 0
    if (any(clash)) {
      hit <- which(apply(clash, 1L, any))
      stop(sprintf(paste0("kalman_filter: `shock_means` and `known_shocks` both ",
                          "specify %s. A known shock's REALISATION is fixed, so ",
                          "its mean is already determined -- pass one or the ",
                          "other for a given shock and period."),
                   paste(known$names[hit], collapse = ", ")), call. = FALSE)
    }
  }

  ## A low-frequency series stored on the high-frequency grid must be spaced
  ## in multiples of its aggregation length; a misaligned column would
  ## otherwise be filtered as if it were high-frequency, silently.
  if (!is.null(mf)) .mf_check_pattern(data, mf, obs_vars)

  ll_const <- -0.5 * n_obs * log(2 * pi)
  has_missing <- anyNA(data)

  ## ---- Structured run diagnostics (R4) -----------------------------------
  ## Every routing decision below appends a (from, to, reason) row here, and
  ## .kf_result() turns them into $diagnostics$routing. The point is that a
  ## parity harness should not have to parse warning text to learn that
  ## method = "auto" ran the univariate filter, or that a singular F sent a
  ## multivariate path somewhere else: the run describes itself.
  ##
  ## Appended positionally (no `<<-`) so the record stays a plain local.
  route_log <- list()

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
    if (method == "auto") {
      method <- "standard"
      route_log[[length(route_log) + 1L]] <-
        c(from = "auto", to = "standard",
          reason = "me_extra (per-period measurement error) needs the R loop")
    }
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
    if (method %in% c("univariate", "auto")) {
      route_log[[length(route_log) + 1L]] <-
        c(from = method, to = "standard",
          reason = "shock_scale (per-period shock variances) needs the R loop")
      method <- "standard"
    }
    ## C++ standard fast path bakes HH/Sigma_e -- bypass it.
    has_missing <- TRUE
  }
  ## Diffuse phase with non-unity scales: stop() (scaling inside the diffuse
  ## phase would require per-step P_inf updates not currently implemented).
  if (has_shock_scale && lik_init == "diffuse")
    stop("kalman_filter: lik_init = 'diffuse' is incompatible with shock_scale. ",
         "Use lik_init = 'kappa' or 'stationary' with heteroskedastic shocks.",
         call. = FALSE)

  ## ---- Known-shock metadata (R4) -----------------------------------------
  ## Which injected cells are DETERMINISTIC -- zero prior variance, hence a
  ## mean shift carrying no density (see .kf_univariate_dispatch). The same
  ## static test the dispatch uses to rule out the C++ kernel, computed here
  ## so the result can report what was applied and under which likelihood
  ## convention, without re-deriving it from the recursion.
  known_meta <- NULL
  if (!is.null(known)) {
    sv <- diag(as.matrix(Sigma_e))[known$idx]
    sv_mat <- matrix(sv, length(known$idx), n_T)
    if (has_shock_scale)
      sv_mat <- sv_mat * shock_scale[known$idx, , drop = FALSE]^2
    cell  <- !is.na(known$values)
    hard  <- cell & known$sd == 0
    known_meta <- list(names = known$names,
                       n_cells = sum(cell),
                       n_deterministic = sum(hard & sv_mat <= .KF_ZERO_VAR_TOL),
                       n_soft = sum(cell & known$sd > 0))
  }

  ## Precompute Y - d (broadcast d down each column) once: the per-step loops
  ## then need only one matrix-vector subtraction. Safe with missing data: NA
  ## propagates through Y - d and is caught by the same anyNA / is.finite
  ## checks downstream.
  ## ---- Deterministic shock means: split the trajectory off ---------------
  ## See .kf_shock_means. Subtracting y^det leaves an ordinary filtering
  ## problem; s^det is added back in .kf_result(). Runs AFTER the
  ## mixed-frequency augmentation, so it uses the same TT/RR/ZZ/DD the
  ## recursions do.
  ##
  ## The subtraction is applied to `data` and NOT to `Y_minus_d` alone,
  ## because Y_minus_d is not the only thing the recursions read: the dare
  ## loop forms its own innovation as `data[, t] - ZZ %*% s - d`, and
  ## .kf_diffuse_phase() takes `data` and `d` separately. Adjusting only the
  ## precomputed matrix left both of those evaluating the UNSHIFTED model --
  ## caught by the cross-method test, which is why it asserts every method
  ## rather than trusting one.
  det_path <- NULL
  if (!is.null(mean_path)) {
    det_path <- .kf_det_path(mean_path, TT, RR, ZZ, DD)
    data     <- data - det_path$y_det
  }

  Y_minus_d <- data - d
  ## Per-period count of observations that are simply absent, kept apart from
  ## the components a singular F drops: "not there" and "carries no
  ## information" are different reasons for a shorter conditioning set, and a
  ## parity harness has to be able to tell them apart.
  missing_by_period <- as.integer(colSums(is.na(Y_minus_d)))

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
  ## Skipped when the diffuse + missing-data reroute below is going to fire:
  ## that path emits a more specific notice which already carries this caveat,
  ## and two warnings saying "use univariate" about a call that is ABOUT to use
  ## univariate is noise.
  if (lik_init_orig %in% c("diffuse", "kappa") &&
      !(has_missing && lik_init_orig == "diffuse") &&
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
    if (has_missing && lik_init == "diffuse") {
      ## The exact-diffuse recursion and missing observations meet only in the
      ## SEQUENTIAL (univariate) filter. .kf_diffuse_phase() processes the
      ## observation vector as a block and bails out on any NA -- correctly, it
      ## has no way to drop a component -- whereas the univariate loop skips a
      ## missing observable one at a time, which is exactly right and is what
      ## it already does for every other initialisation.
      ##
      ## This used to fall through to "standard", which then downgraded
      ## lik_init to "kappa" and returned a DIFFERENT likelihood: on the
      ## local-level fixture with three gaps, -44.097 against the exact
      ## -36.271, and filtered states that keep moving as kappa grows while
      ## the diffuse ones are scale-free to the last bit.
      method <- "univariate"
      route_log[[length(route_log) + 1L]] <-
        c(from = "auto", to = "univariate",
          reason = paste("exact diffuse initialisation with missing",
                         "observations: only the sequential filter can skip",
                         "a component inside the diffuse phase"))
    }
    else if (has_missing) {
      method <- "standard"
      route_log[[length(route_log) + 1L]] <-
        c(from = "auto", to = "standard",
          reason = "missing observations need the per-step loop")
    }
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
      ## ROUTING GATE (conservative): route to univariate only at
      ## me_variance == 0. Since F3-D both paths implement the SAME true-ME
      ## law, so this is no longer a correctness requirement -- it is kept so
      ## that turning on a measurement-error floor cannot silently change
      ## which ALGORITHM (and hence which round-off/steady-state behaviour)
      ## a diffuse-init likelihood runs through. me_extra likewise forces the
      ## univariate R loop. Keep "standard" in both cases.
      method <- if (.HAS_RCPP_KALMAN_UNI() && is.null(me_extra) &&
                    me_variance == 0)
        "univariate" else "standard"
      route_log[[length(route_log) + 1L]] <-
        c(from = "auto", to = method,
          reason = if (method == "univariate")
            "diffuse initialisation: the sequential filter runs it natively"
          else paste("diffuse initialisation with measurement error or",
                     "me_extra: kept on the multivariate per-step loop"))
    }
    ## F4-B: the Chandrasekhar branch is exact again (see METHOD 2), so it is
    ## back in "auto" -- but only above the MEASURED crossover. The increment
    ## recursion runs in R at O(n_state^2 n_obs) per step, the "standard"
    ## filter in C++ at O(n_state^3); on this machine (T = 200, random dense
    ## systems, R CMD INSTALL build) chandrasekhar/standard wall clock was
    ##   n_state   50    100    200    300      (n_obs = 3)
    ##   ratio   1.85   1.01   0.64   0.33
    ##   n_state  100    200                    (n_obs = 7)
    ##   ratio   1.11   0.61
    ## i.e. the old n_state > 50 rule was a ~1.9x PESSIMIZATION and the
    ## crossover sits near n_state = 100 (weakly dependent on n_obs).
    else if (n_state > 100) {
      method <- "chandrasekhar"
      route_log[[length(route_log) + 1L]] <-
        c(from = "auto", to = "chandrasekhar",
          reason = sprintf("n_state = %d > 100 (measured crossover)", n_state))
    } else {
      method <- "standard"
      route_log[[length(route_log) + 1L]] <-
        c(from = "auto", to = "standard", reason = "default for n_state <= 100")
    }
  }

  ## A user-supplied P0 is not the Lyapunov prior the Chandrasekhar increment
  ## recursion is initialised from, and that branch refuses it loudly below --
  ## so route around it. Deliberately OUTSIDE the `if (method == "auto")`
  ## chain above: dropping it in between two `else if` arms silently re-parses
  ## the tail of the chain as `if (user) ... else if (n_state > 100) ... else
  ## "standard"`, which sends EVERY call to "standard" no matter what the
  ## earlier arms decided.
  ## Known shocks are observations of eps, and eps is only a state component on
  ## the univariate path (x_t = [s_{t-1}; eps_t]). The multivariate recursions
  ## have nowhere to put the constraint, so route rather than refuse -- the same
  ## call made for the exact-diffuse phase with missing data.
  if (!is.null(known) && method != "univariate") {
    if (method_orig != "auto")
      warning("kalman_filter: `known_shocks` observes eps directly, which only ",
              "the univariate (Koopman-Durbin) filter can express -- it runs on ",
              "the augmented state [s_{t-1}; eps_t], where the shocks ARE state ",
              "components. Ignoring method = \"", method, "\" and reporting ",
              "method = \"univariate\".", call. = FALSE)
    route_log[[length(route_log) + 1L]] <-
      c(from = method, to = "univariate",
        reason = paste("known_shocks observes eps directly, which only the",
                       "augmented-state (univariate) filter can express"))
    method <- "univariate"
  }

  ## ...but only when "auto" chose it. An EXPLICIT method = "chandrasekhar"
  ## falls through to that branch's own scope guard, which refuses a non-
  ## Lyapunov initialisation loudly -- this file's rule is to refuse rather
  ## than quietly answer a different question.
  if (lik_init == "user" && method == "chandrasekhar" && method_orig == "auto") {
    route_log[[length(route_log) + 1L]] <-
      c(from = "chandrasekhar", to = "standard",
        reason = "a user-supplied P0 is not the Lyapunov prior the Chandrasekhar recursion starts from")
    method <- "standard"
  }

  ## A diffuse phase requires the per-step R loop (no Rcpp / Chandrasekhar
  ## fast path); for large state vectors with unit roots, fall back to the
  ## R "standard" loop entirely. This is a perf cost but correctness-first.
  use_diffuse_phase <- lik_init %in% c("diffuse", "kappa") &&
    method %in% c("standard", "dare", "chandrasekhar")
  if (use_diffuse_phase && method == "chandrasekhar") {
    route_log[[length(route_log) + 1L]] <-
      c(from = "chandrasekhar", to = "standard",
        reason = "a diffuse phase needs the per-step loop")
    method <- "standard"
  }

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

  ## -- Result assembler: one shape for every exit --------------------------
  ## Nine `return`s used to build the result list by hand, which is how the
  ## per-path fields drifted apart in the first place. Everything that is the
  ## same everywhere is assembled here instead, so `$diagnostics` has one
  ## schema no matter which path ran, and a new field cannot be added to some
  ## exits and forgotten at others.
  ##
  ## `final_state`/`final_cov` are the state hand-off: (s_{T|T}, Var(s_T|y))
  ## in the SAME convention `a0`/`P0` take, so filtering a split sample is
  ## kf2 <- kalman_filter(y2, ..., a0 = kf1$final_state, P0 = kf1$final_cov).
  ## `final_cov` is NULL only where the path never forms P (Chandrasekhar
  ## propagates low-rank increments instead) or where the run failed.
  ##
  ## `route_log` and `known_meta` are read from the enclosing frame at CALL
  ## time, so every exit sees the routing decisions that actually happened.
  .kf_result <- function(loglik, filtered_states, method_used, lik_init_used,
                         d_diffuse, final_state = NULL, final_cov = NULL,
                         loglik_contrib = NULL, dropped = NULL,
                         known_applied = NULL, fallback = NULL,
                         extra = list()) {
    ## ---- Timing contract, and the deterministic add-back ----------------
    ## Two state paths are reported, and the names say which is which rather
    ## than leaving "filtered" to be interpreted:
    ##   updated_states[, t]   = s_{t|t}    conditioned on y_1..t
    ##   predicted_states[, t] = s_{t|t-1}  conditioned on y_1..t-1
    ## `filtered_states` is the SAME matrix as `updated_states`, kept because
    ## it is the name the rest of the package and every existing caller uses.
    ##
    ## The prediction needs no extra state from the recursions: E[eps_t] = 0 in
    ## the deviation system, so s_{t|t-1} = T s_{t-1|t-1} exactly, with
    ## s_{1|0} = T a0. A deterministic mean path adds R m_t on top, which is
    ## precisely s^det_t - T s^det_{t-1}; so predicting from the STOCHASTIC
    ## filtered path and then adding s^det gives both series in one step.
    predicted_states <- NULL
    if (!is.null(filtered_states)) {
      prev <- cbind(a0_vec, filtered_states[, -n_T, drop = FALSE])
      predicted_states <- TT %*% prev
      if (!is.null(det_path)) {
        filtered_states  <- filtered_states  + det_path$s_det
        predicted_states <- predicted_states + det_path$s_det
      }
      rownames(filtered_states)  <- state_names_out
      rownames(predicted_states) <- state_names_out
    }
    if (!is.null(final_state)) {
      final_state <- as.numeric(final_state)
      ## The hand-off is the TOTAL state, so a second call started from it
      ## carries the deterministic level already accumulated and restarts its
      ## own mean path from zero.
      if (!is.null(det_path)) final_state <- final_state + det_path$s_det[, n_T]
      names(final_state) <- state_names_out
    }
    if (!is.null(final_cov)) {
      final_cov <- as.matrix(final_cov)
      dimnames(final_cov) <- list(state_names_out, state_names_out)
    }
    routes <- route_log
    ## `fallback` is one row, or several: a singularity fallback into the
    ## univariate filter can be followed by that filter's own diffuse -> kappa
    ## fallback, and both belong in the record.
    if (!is.null(fallback))
      routes <- c(routes, if (is.list(fallback)) fallback else list(fallback))
    routing <- .kf_routing_df(routes)
    dd <- if (is.null(d_diffuse) || length(d_diffuse) != 1L) NA_integer_
          else as.integer(d_diffuse)
    if (is.null(dropped)) dropped <- integer(n_T)
    ## Likelihood convention. Only the injected shocks can move it off
    ## "marginal": a shock observed with a prior density contributes it
    ## (JOINT), a deterministic one is a point mass and contributes nothing
    ## (CONDITIONAL). A mixture of the two is joint in the ones that have a
    ## density -- reported as "joint" so it is never read as conditional.
    ll_type <- if (is.null(known_meta)) "marginal"
      else if (known_meta$n_deterministic == known_meta$n_cells) "conditional"
      else "joint"
    c(list(loglik = loglik, filtered_states = filtered_states,
           updated_states = filtered_states,
           predicted_states = predicted_states,
           loglik_contrib = loglik_contrib,
           n_obs = n_obs, n_T = n_T, method = method_used,
           lik_init = lik_init_used, d_diffuse = dd,
           final_state = final_state, final_cov = final_cov,
           diagnostics = list(
             method_requested   = method_orig,
             method_used        = method_used,
             lik_init_requested = lik_init_orig,
             lik_init_used      = lik_init_used,
             routing            = routing,
             diffuse_periods    = if (is.na(dd)) integer(0) else seq_len(dd),
             missing_by_period  = missing_by_period,
             n_missing          = sum(missing_by_period),
             dropped_by_period  = as.integer(dropped),
             n_dropped          = sum(as.integer(dropped)),
             known_shocks       = if (is.null(known_meta)) NULL
                                  else c(known_meta,
                                         list(n_applied = known_applied)),
             ## A deterministic mean path is an INPUT, so it is reported as
             ## one: it conditions nothing and adds no density term, and
             ## `loglik_type` stays "marginal" for it.
             shock_means        = if (is.null(mean_path)) NULL
                                  else list(timing = shock_timing,
                                            n_cells = sum(mean_path != 0),
                                            names = exo[rowSums(mean_path != 0) > 0]),
             loglik_type        = ll_type)),
      extra)
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
    ## TRUE measurement-noise law (F3-D): y_t = ZZ x_t + DD e_t + u_t with
    ## Var(u_t) = me_variance * I requires P' += K me_variance I K' for ANY
    ## gain K. Before F3-D `me_variance` entered F only (a regulariser), which
    ## made the multivariate paths disagree with the univariate filter, the
    ## smoother and the DARE fixed point by O(me_variance).
    if (me_variance != 0) P_n <- P_n + me_variance * tcrossprod(K)
    P_n  <- (P_n + t(P_n)) * 0.5
    list(ll = ll, s = s_n, P = P_n, K = K, F_inv = Fi, log_det_F = ldf, F_mat = Ft)
  }

  ## -- Univariate (sequential) filter runner ----------------------------
  ## Runs the Koopman-Durbin (2000) univariate filter on the augmented
  ## state [s; eps] (see .kf_univariate_dispatch). Used (a) directly for
  ## method = "univariate", (b) as the automatic singularity fallback from
  ## every multivariate path (.kf_fail below), and (c) for the exact
  ## diffuse phase when F_inf is singular but nonzero (Case C).

  .run_univariate <- function(li, ss_lock = FALSE, fallback = NULL) {
    P_inf_state <- NULL
    if (li == "user") {
      P_state <- P0_mat
    } else if (li == "diffuse") {
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
    ## The steady-state lock freezes the gain on the converged tail, which is
    ## reasoning about the STATIONARY prior; a supplied P0 is not that, so the
    ## lock is dropped rather than assumed to still apply.
    out <- .kf_univariate_dispatch(Y_minus_d, ZZ, TT, RR, DD, Sigma_e,
                                   a0_vec, P_state, P_inf_state,
                                   me_variance, return_filtered,
                                   me_extra = me_extra,
                                   ss_lock = ss_lock && is.null(P0_mat),
                                   shock_scale = ss_arg, known = known)
    if (isTRUE(out$diffuse_failed)) {
      ## P_inf never decayed (unobserved unit root) or the sample ended
      ## inside the diffuse phase: same kappa fallback as the multivariate
      ## diffuse path.
      warning("kalman_filter: univariate diffuse phase: P_inf did not ",
              "converge to zero; falling back to lik_init = \"kappa\".")
      fb_diffuse <- c(from = "diffuse", to = "kappa",
                      reason = "P_inf did not converge to zero within the sample")
      fallback <- c(if (is.null(fallback)) list()
                    else if (is.list(fallback)) fallback else list(fallback),
                    list(fb_diffuse))
      li  <- "kappa"
      out <- .kf_univariate_dispatch(Y_minus_d, ZZ, TT, RR, DD, Sigma_e,
                                     a0_vec, .build_P0(TT, QQ),
                                     NULL, me_variance, return_filtered,
                                     me_extra = me_extra, shock_scale = ss_arg,
                                     known = known)
    }
    if (!isTRUE(out$ok))
      return(.kf_result(-Inf, NULL, "univariate", li, NA_integer_,
                        fallback = fallback))
    filt <- NULL
    if (return_filtered && !is.null(out$filtered)) {
      filt <- out$filtered
      rownames(filt) <- state_names_out
    }
    .kf_result(out$loglik, filt, "univariate", li, out$d_diffuse,
               final_state = out$a[seq_len(n_state)],
               final_cov   = out$P_state,
               loglik_contrib = if (return_ll_contrib)
                 as.numeric(out$ll_contrib) else NULL,
               dropped = out$n_skipped,
               known_applied = if (is.null(out$det_applied)) NULL
                               else sum(out$det_applied),
               fallback = fallback)
  }

  ## -- A5: the singularity-fallback CONTRACT ----------------------------
  ##
  ## A multivariate path hit a singular / non-PD innovation covariance (or a
  ## numerically exploded step): retry on the univariate filter before
  ## declaring -Inf -- the analog of Dynare's
  ## univariate_kalman_filter_if_singularity_is_detected.
  ##
  ## But the univariate filter does NOT always evaluate the same likelihood
  ## as the multivariate paths, and switching conventions mid-chain is worse
  ## than rejecting the draw: an MCMC chain in which SOME draws trip the
  ## fallback would then be sampling a mixture of two different densities,
  ## silently. The two conventions coincide only when ALL of the following
  ## hold (each condition documented elsewhere in this file):
  ##   * me_variance == 0 -- a CONSERVATIVE gate retained from the pre-F3-D
  ##     regulariser convention. Since F3-D the multivariate paths treat
  ##     me_variance as TRUE iid measurement noise, exactly like the
  ##     univariate filter (test-kf-true-me.R pins agreement to 1e-10 at
  ##     me_variance in {1e-3, 1e-2}), so the fallback would now be sound at
  ##     me_variance > 0 too; it is left closed so this wave changes the
  ##     LIKELIHOOD only, not the failure policy;
  ##   * me_extra inactive -- gated for the same reason (per-observable,
  ##     per-period ME, and the me_extra routing block forces the R loops);
  ##   * lik_init has not resolved to the EXACT-DIFFUSE init -- the
  ##     univariate diffuse filter keeps the divergent 0.5*log(F_inf) term
  ##     that the multivariate exact-diffuse recursion renormalises away
  ##     (see the M23 warning above: the gap exceeds 15 nats on double-unit-
  ##     root models). "kappa" and "stationary" are fine: both paths build
  ##     the same P0.
  ## Read lik_init LAZILY (not captured at closure-creation time): the
  ## diffuse block below may downgrade lik_init to "kappa" before the
  ## METHOD 1/3 loops run, and the fallback is legitimate again after that.
  ##
  ## Otherwise the draw is rejected with loglik = -Inf -- which is exactly
  ## what Dynare does when
  ## univariate_kalman_filter_if_singularity_is_detected is not set.
  ##
  ## A genuinely bad draw fails the same per-period ll floor inside the
  ## univariate filter too, so a permitted fallback never turns a true -Inf
  ## into a finite value.
  ##
  ## `.kf_fallback_warned` is a LOCAL latch of this kalman_filter() call (not
  ## a package-level env): one warning per call, no cross-call state.
  .kf_fallback_warned <- FALSE

  .kf_fail <- function(failed_method) {
    fb <- c(from = failed_method, to = "univariate",
            reason = paste("singular / non-positive-definite innovation",
                           "covariance (or a non-finite step)"))
    hard_fail <- .kf_result(-Inf, NULL, failed_method, lik_init, d_diffuse)
    if (!(me_variance == 0 && !has_me_extra && lik_init != "diffuse"))
      return(hard_fail)
    out <- tryCatch(.run_univariate(lik_init, fallback = fb),
                    error = function(e) NULL)
    if (is.null(out)) return(hard_fail)
    if (!.kf_fallback_warned) {
      .kf_fallback_warned <<- TRUE
      warning("kalman_filter: method = \"", failed_method, "\" hit a ",
              "singular / non-positive-definite innovation covariance ",
              "(or a non-finite step); the log-likelihood was evaluated ",
              "with the UNIVARIATE (Koopman-Durbin) filter instead. The two ",
              "agree exactly under the current settings (me_variance = 0, ",
              "no me_extra, lik_init = \"", lik_init, "\"), but $method is ",
              "reported as \"univariate\" -- expect this on some draws only.",
              call. = FALSE)
    }
    out
  }

  if (method == "univariate") return(.run_univariate(lik_init, ss_lock = ss_lock_req))

  ## -- Initialization (s0, P0) for the chosen lik_init -----------------
  ## Computed once, used by both METHOD 1 ("dare") and METHOD 3 ("standard").
  ## "kappa" and "diffuse" require methods %in% c("standard", "dare") (forced
  ## above); "stationary" / "auto"-resolved-to-"stationary" reuse the
  ## historical solve_lyapunov(TT, QQ) P0 unconditionally.
  init_s <- a0_vec
  init_P <- NULL
  init_loglik <- 0
  init_t_start <- 1L

  if (lik_init == "user") {
    init_P <- P0_mat
  } else if (lik_init == "kappa") {
    init_P <- .build_P0(TT, QQ)
  } else if (lik_init == "diffuse") {
    if (has_missing) {
      ## The caller asked for a multivariate method explicitly (method = "auto"
      ## already routes this combination to the univariate filter). Route it
      ## anyway rather than downgrading: the sequential filter evaluates the
      ## EXACT diffuse likelihood with gaps, which is what was asked for, while
      ## the old "kappa" downgrade answered a different question -- on the
      ## local-level fixture with three gaps, -44.097 instead of -36.271.
      warning("kalman_filter: method = \"", method, "\" cannot run the ",
              "exact-diffuse recursion with missing observations (the ",
              "multivariate diffuse phase has no way to drop one component of ",
              "the observation vector). Using the UNIVARIATE (Koopman-Durbin) ",
              "filter, which does it exactly, and reporting method = ",
              "\"univariate\". On a model with TWO OR MORE unit roots the two ",
              "diffuse conventions differ by an additive constant (see the ",
              "unit-root note in ?kalman_filter), so do not compare this ",
              "log-likelihood with one from a complete-data multivariate run; ",
              "pass method = \"univariate\" throughout instead.",
              call. = FALSE)
      return(.run_univariate("diffuse"))
    } else {
      P0 <- .kf_diffuse_P0(TT, QQ)
      diff_out <- .kf_diffuse_phase(data, d, ZZ, TT, RR, DD, QQ, HH, SS, Sigma_e,
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
      if (return_filtered) rownames(filtered) <- state_names_out
      return(.kf_result(loglik, filtered, "dare", lik_init, d_diffuse,
                        final_state = s, final_cov = init_P,
                        loglik_contrib = ll_contrib,
                        extra = list(dare_iterations = NA_integer_,
                                     dare_p_drift = NA_real_)))
    }

    for (t in t_start:n_T) {
      v    <- data[, t] - as.numeric(ZZ %*% s) - d
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
        ## Joseph true-noise term for the FULL measurement-error diagonal
        ## (base me_variance + this period's me_extra):
        ## P' += K_t diag(me_variance + me_extra[, t]) K_t'.
        me_vec_t <- rep(me_variance, n_obs)
        if (me_x_t) me_vec_t <- me_vec_t + me_extra[, t]
        if (any(me_vec_t != 0)) P <- P + K_t %*% (me_vec_t * t(K_t))
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

    if (return_filtered) rownames(filtered) <- state_names_out
    return(.kf_result(loglik, filtered, "dare", lik_init, d_diffuse,
                      final_state = s, final_cov = P,
                      loglik_contrib = ll_contrib,
                      extra = list(dare_iterations = dare$iterations,
                                   dare_p_drift    = final_drift)))
  }


  ## ===================================================================
  ## METHOD 2: Chandrasekhar recursions (Morf-Sidhu-Kailath; Herbst 2015)
  ## ===================================================================
  ##
  ## For a TIME-INVARIANT system started at the stationary covariance
  ## P_1 = P_0 = solve_lyapunov(TT, QQ) the Riccati INCREMENT
  ##   dP_t = P_{t+1} - P_t
  ## has rank <= n_obs for every t, so the filter can propagate the low-rank
  ## factors (W_t, M_t) of dP_t = W_t M_t W_t' (n_state x n_obs and
  ## n_obs x n_obs) instead of the n_state x n_state matrix P_t. With
  ##   F_t = ZZ P_t ZZ' + HH + me_diag,   K_t = (TT P_t ZZ' + SS) F_t^{-1}
  ## the EXACT recursions are (ZW = ZZ W_t, TW = TT W_t):
  ##   F_{t+1} = F_t + ZW M_t ZW'
  ##   K_{t+1} = K_t + (TW - K_t ZW) M_t ZW' F_{t+1}^{-1}
  ##   W_{t+1} = TW - K_{t+1} ZW
  ##   M_{t+1} = M_t + M_t ZW' F_t^{-1} ZW M_t        <-- OLD F_t, M on BOTH sides
  ## and the likelihood contributions are the standard filter's,
  ##   ll_t = ll_const - 0.5 (log|F_t| + v_t' F_t^{-1} v_t).
  ## Initialisation needs no bootstrap: with P_1 = P_0 the Riccati gives
  ## dP_1 = TT P_0 TT' + QQ - K_1 F_1 K_1' - P_0 = -K_1 F_1 K_1', i.e.
  ##   W_1 = K_1,  M_1 = -F_1  (exactly rank n_obs).
  ##
  ## F4-B (2026-09-04) rewrote this branch. The previous implementation
  ## bootstrapped with standard Riccati steps and an eigen-factorisation of
  ## dP, and was wrong in four independent ways:
  ##   (D1) when the bootstrap loop finished the whole sample without its
  ##        |dP| < 0.1|P| break firing (short samples), boot_steps stayed 0
  ##        and the Chandrasekhar phase restarted at t = 1 -- every
  ##        observation's contribution was counted TWICE (rbc, T = 8: 23.1
  ##        nats off the exact likelihood);
  ##   (D2) the gain update used K_{t+1} = K_t + TW M ZW' F_{t+1}^{-1},
  ##        dropping the -K_t ZW M ZW' F_{t+1}^{-1} term implied by
  ##        K_{t+1} F_{t+1} = K_t F_t + TW M ZW';
  ##   (D3) the increment update used M_{t+1} = M_t + ZW' F_{t+1}^{-1} ZW M_t
  ##        -- the NEW F instead of F_t, and only ONE factor of M_t (the
  ##        correct form is quadratic in M_t, and is what keeps M_t
  ##        symmetric);
  ##   (D4) the eigen-factorisation of dP at the bootstrap endpoint kept up
  ##        to n_state directions with a relative 1e-14 cut, so the carried
  ##        rank (and hence the recursion) depended on round-off; combined
  ##        with (D2)/(D3) this is what drove F to lose positive-definiteness
  ##        on fs2000 at me_variance = 1e-4 (chol failure -> -Inf).
  ## The recursions above are exact, so no bootstrap, no eigen and no
  ## near-cancellation heuristic is needed.

  if (method == "chandrasekhar") {
    ## -- Scope guards: refuse LOUDLY rather than return a wrong number. -----
    ## me_extra / shock_scale are already rejected upstream (both make the
    ## innovation structure time-varying, which breaks the increment
    ## recursion); missing observations and a non-stationary initialisation
    ## do the same, so they must not fall through to a silent -Inf.
    if (anyNA(data))
      stop("kalman_filter: method = 'chandrasekhar' is incompatible with ",
           "missing observations (the increment recursion assumes a ",
           "time-invariant observation equation); use method = 'auto', ",
           "'standard' or 'univariate'.", call. = FALSE)
    if (lik_init != "stationary")
      stop("kalman_filter: method = 'chandrasekhar' requires lik_init = ",
           "\"stationary\" (the increment recursion is initialised from the ",
           "Lyapunov P0); got lik_init = \"", lik_init, "\". Use method = ",
           "'standard' or 'dare' for a diffuse/kappa initialisation.",
           call. = FALSE)

    HH_full <- HH + me_diag
    P0      <- .solve_lyapunov_stationary()

    ## --- Exact initialisation: F_1, K_1 at P_1 = P_0, W_1 = K_1, M_1 = -F_1
    PZ0   <- P0 %*% tZZ
    F_mat <- ZZ %*% PZ0 + HH_full
    F_mat <- (F_mat + t(F_mat)) * 0.5
    Fc    <- tryCatch(chol(F_mat), error = function(e) NULL)
    if (is.null(Fc)) return(.kf_fail("chandrasekhar"))
    F_inv     <- chol2inv(Fc)
    log_det_F <- 2 * sum(log(diag(Fc)))
    K         <- (TT %*% PZ0 + SS) %*% F_inv
    W         <- K
    M         <- -F_mat

    ## a0 composes with the increment recursion (only P0 is constrained to the
    ## Lyapunov prior, which the scope guard above enforces).
    s <- a0_vec; loglik <- 0
    ch_ss_step <- NA_integer_

    for (t in seq_len(n_T)) {
      v    <- Y_minus_d[, t] - as.numeric(ZZ %*% s)
      ll_t <- ll_const - 0.5 * (log_det_F + drop(crossprod(v, F_inv %*% v)))
      if (!is.finite(ll_t) || ll_t < .KF_LL_MIN)
        return(.kf_fail("chandrasekhar"))
      loglik <- loglik + ll_t
      s <- as.numeric(TT %*% s) + drop(K %*% v)
      if (return_filtered) filtered[, t] <- s
      if (t == n_T) break

      ## Shared low-rank blocks. WM is n_state x n_obs, ZWM is n_obs x n_obs.
      ZW  <- ZZ %*% W
      TW  <- TT %*% W
      WM  <- W %*% M
      ZWM <- ZZ %*% WM

      ## Steady-state lock: dP_t = W M W' is exactly the quantity the
      ## standard filter compares against ss_tol (max|P_{t+1} - P_t|), and
      ## the frozen (K_t, F_t) are the same ones -- so the two methods lock
      ## at the same period with the same gain.
      if (t > 1L && max(abs(tcrossprod(WM, W))) < ss_tol) {
        ch_ss_step  <- t
        ll_ss_const <- ll_const - 0.5 * log_det_F
        tail_start  <- t + 1L
        out <- .kf_ss_dispatch(Y_minus_d, ZZ, TT, K, F_inv, ll_ss_const,
                               s, tail_start, n_T, filtered)
        if (!out$ok) return(.kf_fail("chandrasekhar"))
        loglik <- loglik + out$loglik
        s      <- out$s
        if (return_filtered) filtered <- out$filtered
        break
      }

      ## --- Chandrasekhar increment update -------------------------------
      F_new <- F_mat + ZWM %*% t(ZW)
      F_new <- (F_new + t(F_new)) * 0.5
      Fc_n  <- tryCatch(chol(F_new), error = function(e) NULL)
      if (is.null(Fc_n)) return(.kf_fail("chandrasekhar"))
      F_new_inv     <- chol2inv(Fc_n)
      log_det_F_new <- 2 * sum(log(diag(Fc_n)))

      tZWM  <- t(ZWM)                       # = M ZW'
      K_new <- K + (TW - K %*% ZW) %*% (tZWM %*% F_new_inv)
      W_new <- TW - K_new %*% ZW
      ## M_{t+1} = M + M ZW' F_t^{-1} ZW M -- quadratic in M, hence symmetric.
      M_new <- M + tZWM %*% F_inv %*% ZWM
      M_new <- (M_new + t(M_new)) * 0.5

      K <- K_new; F_mat <- F_new; F_inv <- F_new_inv
      log_det_F <- log_det_F_new; W <- W_new; M <- M_new
    }

    if (return_filtered) rownames(filtered) <- state_names_out
    ## No `final_cov`: this recursion propagates the low-rank increments
    ## (W, M) precisely so that P is never formed. Ask for the state hand-off
    ## with method = "standard" (or let "auto" pick it) rather than paying
    ## O(n_state^2) per step to rebuild what this method exists to avoid.
    return(.kf_result(loglik, filtered, "chandrasekhar", lik_init, d_diffuse,
                      final_state = s, final_cov = NULL,
                      extra = list(boot_steps = 0L,
                                   ss_reached_at = ch_ss_step)))
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
  ## `a0` gate: kalman_standard_loop_cpp() takes no initial state -- it starts
  ## the recursion at zero, unconditionally. Without this test a non-zero `a0`
  ## would be accepted, validated, reported, and then SILENTLY DROPPED on the
  ## default fast path for every stationary model. Fall through to the R loop,
  ## which reads init_s.
  if (!has_missing && !has_shock_scale && lik_init == "stationary" &&
      all(a0_vec == 0) && .HAS_RCPP_KALMAN()) {
    out <- kalman_standard_loop_cpp(Y_minus_d, ZZ, TT, RR, DD, HH + me_diag,
                                    Sigma_e, SS, P, ll_const, ss_tol,
                                    .KF_LL_MIN, return_filtered,
                                    rep(me_variance, n_obs))
    if (!out$ok)
      return(.kf_fail("standard"))
    filt <- NULL
    if (return_filtered) {
      filt <- out$filtered
      rownames(filt) <- state_names_out
    }
    return(.kf_result(out$loglik, filt, "standard", lik_init, d_diffuse,
                      final_state = out$s, final_cov = out$P))
  }

  if (init_t_start > n_T) {
    ## The diffuse phase consumed the entire sample.
    if (return_filtered) rownames(filtered) <- state_names_out
    return(.kf_result(loglik, filtered, "standard", lik_init, d_diffuse,
                      final_state = s, final_cov = P))
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
      ## Joseph true-noise term for the FULL ME diagonal (observed subset
      ## only): y = Z s + D e + u with Var(u) = diag(me_variance +
      ## me_extra[obs_ok, t]) requires P' += K diag(.) K' for ANY gain K.
      me_vec_t <- rep(me_variance, n_obs_t)
      if (has_me_extra) me_vec_t <- me_vec_t + me_extra[obs_ok, t]
      if (any(me_vec_t != 0)) P <- P + K %*% (me_vec_t * t(K))
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
        ## Joseph true-noise term for the FULL ME diagonal me_diag_t
        ## (base me_variance + this period's me_extra).
        me_vec_t <- rep(me_variance, n_obs)
        if (has_me_extra) me_vec_t <- me_vec_t + me_extra[, t]
        if (any(me_vec_t != 0)) P <- P + K_t %*% (me_vec_t * t(K_t))
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
        ## Joseph true-noise term for the FULL ME diagonal (base me_variance
        ## + me_extra[, t]): P' += K diag(.) K'.
        P     <- P + K_t %*% ((me_variance + me_extra[, t]) * t(K_t))
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

  if (return_filtered) rownames(filtered) <- state_names_out
  .kf_result(loglik, filtered, "standard", lik_init, d_diffuse,
             final_state = s, final_cov = P)
}
