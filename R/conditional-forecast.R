## R/conditional-forecast.R
## --------------------------------------------------------------------------
## Waggoner-Zha (1999) conditional forecasting on a linear state-space.
##
## State-space convention (same as smoother-monolith.R / build_dsge_state_space):
##   s_h = T s_{h-1} + R eps_h          [n_state x 1]
##   y_h = Z s_{h-1} + D eps_h          [n_obs   x 1]
##
## where s_0 = s_{T|T} (filtered terminal state) and deviations are from
## steady state (all variables are in deviation form throughout).
##
## Core references:
##   Waggoner, D. F. & Zha, T. (1999). Conditional Forecasts in Dynamic
##     Multivariate Models. Review of Economics and Statistics, 81(4), 639-651.
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Internal: build the stacked A and M matrices
##
## For the H-period forecast horizon the stacked system is:
##   Y_stack = M s_0 + A eps_stack
##
## where Y_stack = [y_1; y_2; ...; y_H] (n_obs*H x 1),
##       eps_stack = [eps_1; eps_2; ...; eps_H] (n_shock*H x 1),
##       s_0 = terminal state (n_state x 1).
##
## M[h, :] = Z T^{h-1}   (n_obs x n_state block for horizon h)
##
## A is block lower-triangular:
##   A[(h-1)*n_obs + 1 : h*n_obs, (j-1)*n_shock + 1 : j*n_shock]
##     = Z T^{h-1-j} R    for j < h
##     = D                for j == h
##     = 0                for j > h
## ---------------------------------------------------------------------------
.build_forecast_stack <- function(ss, H) {
  T_mat  <- ss$T_mat
  R_mat  <- ss$R_mat
  Z_mat  <- ss$Z_mat
  D_mat  <- ss$D_mat
  n_s    <- ss$n_state
  n_obs  <- ss$n_obs
  n_shk  <- ss$n_shock

  ## Pre-compute T^0, T^1, ..., T^{H-1}
  T_pow <- vector("list", H)
  T_pow[[1]] <- diag(n_s)       ## T^0 = I
  for (h in 2:H) T_pow[[h]] <- T_pow[[h - 1L]] %*% T_mat

  ## M matrix (n_obs*H x n_state)
  M <- matrix(0, n_obs * H, n_s)
  for (h in seq_len(H)) {
    rows <- (h - 1L) * n_obs + seq_len(n_obs)
    M[rows, ] <- Z_mat %*% T_pow[[h]]     ## Z T^h (T^0 = I, so row 1 = Z)
  }

  ## A matrix (n_obs*H x n_shock*H), block lower-triangular
  A <- matrix(0, n_obs * H, n_shk * H)
  for (h in seq_len(H)) {
    rows <- (h - 1L) * n_obs + seq_len(n_obs)
    ## Diagonal block: D (contemporaneous shock -> obs)
    cols_h <- (h - 1L) * n_shk + seq_len(n_shk)
    A[rows, cols_h] <- D_mat
    ## Sub-diagonal blocks: Z T^{h-1-j} R for j = 1, ..., h-1
    for (j in seq_len(h - 1L)) {
      cols_j <- (j - 1L) * n_shk + seq_len(n_shk)
      ## power index: h-1-j; T_pow index is h-j (1-based, T_pow[[1]] = T^0)
      A[rows, cols_j] <- Z_mat %*% T_pow[[h - j]] %*% R_mat
    }
  }

  list(A = A, M = M, T_pow = T_pow)
}


## ---------------------------------------------------------------------------
## Internal: hard anticipated conditioning (Waggoner-Zha 1999)
##
## Minimum-variance (most likely) shock path under Q_stack = diag(Q_diag),
## computed in whitened coordinates A_w = A_c W, W = diag(sqrt(Q_diag)):
##   eps_star = W A_w' (A_w A_w')^{-1} b_c
##
## Conditional covariance for draws (eta space, mapped back via W):
##   Sigma_c = I - A_w' (A_w A_w')^{-1} A_w
##
## Returns point forecast paths and optionally a draw matrix.
##
## @param A_full     n_obs*H x n_shock*H stacked response matrix
## @param M          n_obs*H x n_state   state-to-obs propagator
## @param s0         n_state x 1 terminal state
## @param cond_rows  Integer vector indexing rows of A_full that are conditioned
## @param b_cond     Numeric vector of condition values minus deterministic mean
## @param free_cols  Integer vector indexing COLUMNS of A_full to use (free shocks)
## @param n_draws    Number of draws (0 = point only)
## @param Q_diag     Shock variance diagonal (length n_shock*H or scalar 1)
## @return List: eps_star (n_shock*H), paths (H x n_obs), draws (list), shock_paths
## ---------------------------------------------------------------------------
.hard_anticipated <- function(A_full, M, s0, cond_rows, b_cond,
                              free_cols, n_draws, Q_diag, ss, H) {
  n_obs  <- ss$n_obs
  n_shk  <- ss$n_shock

  ## Restrict to free shocks
  A_c <- A_full[cond_rows, free_cols, drop = FALSE]

  ## Check row rank
  n_cond <- length(cond_rows)
  if (n_cond > length(free_cols)) {
    stop(sprintf(
      "conditional_forecast: %d conditions but only %d free shock-periods -- system is over-determined.",
      n_cond, length(free_cols)), call. = FALSE)
  }
  rk <- qr(A_c)$rank
  if (rk < n_cond) {
    stop(sprintf(
      "conditional_forecast: A_cond has rank %d < %d conditions -- the restriction is infeasible with the chosen free_shocks.",
      rk, n_cond), call. = FALSE)
  }

  ## Q-weighted minimum-variance solution (Waggoner-Zha): with stacked shock
  ## covariance Q_stack = diag(Q_diag), whiten eta = W^{-1} eps with
  ## W = diag(sqrt(Q_diag[free])) and solve in eta space, so the conditional
  ## mean is the most likely shock path under the model's shock covariance
  ## and conditioned directions get exactly zero draw variance for any Q.
  ## Q_diag = NULL keeps the legacy unit-metric minimum-norm solution.
  n_free <- length(free_cols)
  w   <- if (!is.null(Q_diag)) sqrt(Q_diag[free_cols]) else rep(1, n_free)
  A_w <- sweep(A_c, 2L, w, "*")

  ## eta_star = A_w' (A_w A_w')^{-1} b_cond; eps_star = W eta_star
  AtA  <- tcrossprod(A_w)          ## n_cond x n_cond
  AtA  <- (AtA + t(AtA)) * 0.5
  ## Minimum-norm solve, not a ridge. `AtA` is rank-deficient whenever the
  ## conditions are collinear or over-specified -- forcing two series that the
  ## model cannot independently hit, or conditioning through a shock that has
  ## been switched off (weight 0) -- which is a normal thing for a user to ask
  ## for, not an error. The old `chol(AtA + 1e-12 * diag(n_cond))` was both
  ## ABSOLUTE (meaningless when AtA is not O(1)) and UNGUARDED, so it threw.
  ## .safe_inv() truncates on a RELATIVE singular-value cutoff and returns the
  ## minimum-norm solution, warning when it does.
  AtA_inv <- .safe_inv(AtA, warn_label = "conditional_forecast: condition system")

  eps_free_star <- w * (t(A_w) %*% (AtA_inv %*% b_cond))   ## n_free x 1

  eps_star <- numeric(ncol(A_full))
  eps_star[free_cols] <- eps_free_star

  ## Recover paths
  ## Y_stack is [y_1; y_2; ...; y_H] with each y_h of length n_obs,
  ## so reshape with byrow = TRUE: row h gets elements [(h-1)*n_obs+1 .. h*n_obs].
  Y_stack_point  <- as.numeric(A_full %*% eps_star) + as.numeric(M %*% s0)
  paths_point    <- matrix(Y_stack_point, nrow = H, ncol = n_obs, byrow = TRUE)
  colnames(paths_point) <- ss$obs_names_fcst

  ## Shock paths (H x n_shock): eps_star is [eps_1; ...; eps_H] stacked
  shock_paths <- matrix(eps_star, nrow = n_shk, ncol = H, byrow = FALSE)
  shock_paths <- t(shock_paths)
  colnames(shock_paths) <- ss$shock_names

  draws_out <- NULL
  if (n_draws > 0L) {
    ## Conditional covariance in whitened eta space:
    ##   Sigma_c = I_{n_free} - A_w' (A_w A_w')^{-1} A_w
    ## (exactly singular along conditioned directions); draws map back to
    ## eps space via W. Use SVD for a numerically stable square-root.
    Proj    <- t(A_w) %*% AtA_inv %*% A_w  ## n_free x n_free
    Sigma_c <- diag(n_free) - Proj
    Sigma_c <- (Sigma_c + t(Sigma_c)) * 0.5

    sv <- svd(Sigma_c)
    ## Clip tiny/negative eigenvalues
    sv$d <- pmax(sv$d, 0)
    Sigma_sqrt <- sv$u %*% diag(sqrt(sv$d), nrow = length(sv$d)) %*% t(sv$u)

    draws_list <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      z <- rnorm(n_free)
      eps_d <- numeric(ncol(A_full))
      eps_d[free_cols] <- as.numeric(eps_free_star) +
        w * as.numeric(Sigma_sqrt %*% z)
      Y_d <- as.numeric(A_full %*% eps_d) + as.numeric(M %*% s0)
      d_mat <- matrix(Y_d, nrow = H, ncol = n_obs, byrow = TRUE)
      colnames(d_mat) <- ss$obs_names_fcst
      draws_list[[i]] <- d_mat
    }
    draws_out <- draws_list
  }

  list(paths_point = paths_point, shock_paths = shock_paths, draws = draws_out)
}


## ---------------------------------------------------------------------------
## Internal: hard unanticipated conditioning
##
## Period-by-period: at each horizon h, agent is surprised; the shock at
## period h (restricted to free_shocks) is chosen to satisfy conditions
## exactly.  The state propagates forward using the realised shocks.
##
## @param ss          State-space list
## @param s0          Terminal state (n_state x 1)
## @param H           Horizon
## @param cond_df     Data.frame with cols: horizon (1..H), var_idx, value
## @param free_shk_idx Integer indices of free shocks (into 1..n_shock)
## @param n_draws      Number of draws
## @param Q            n_shock x n_shock shock covariance (NULL = unit metric)
## @return Same structure as .hard_anticipated
## ---------------------------------------------------------------------------
.hard_unanticipated <- function(ss, s0, H, cond_df, free_shk_idx, n_draws,
                                Q = NULL) {
  T_mat  <- ss$T_mat
  R_mat  <- ss$R_mat
  Z_mat  <- ss$Z_mat
  D_mat  <- ss$D_mat
  n_s    <- ss$n_state
  n_obs  <- ss$n_obs
  n_shk  <- ss$n_shock

  paths_point <- matrix(0, H, n_obs)
  colnames(paths_point) <- ss$obs_names_fcst
  shock_paths <- matrix(0, H, n_shk)
  colnames(shock_paths) <- ss$shock_names

  ## Q-weighted solution, same whitening as .hard_anticipated: with
  ## eta = W^{-1} eps, W = diag(sqrt(diag(Q))), the minimum-norm eta is the
  ## most likely shock under N(0, Q). Q = NULL keeps the unit metric.
  w_all  <- if (!is.null(Q)) sqrt(diag(Q)) else rep(1, n_shk)
  w_free <- w_all[free_shk_idx]

  s_prev <- s0
  for (h in seq_len(H)) {
    ## Unconditional predicted obs and state
    y_pred_mean <- as.numeric(Z_mat %*% s_prev)  ## no shock yet

    ## Conditions at this horizon
    h_cond <- cond_df[cond_df$horizon == h, , drop = FALSE]

    eps_h <- numeric(n_shk)
    if (nrow(h_cond) > 0L) {
      ## Restricted obs equation: y_h = y_pred_mean + D eps_h
      ## For conditioned variables i: D[i, free_shk_idx] eps_free = b_i
      cond_idx <- h_cond$var_idx     ## indices into 1..n_obs
      b_h      <- h_cond$value - y_pred_mean[cond_idx]

      D_c <- D_mat[cond_idx, free_shk_idx, drop = FALSE]
      n_c <- nrow(D_c)
      n_f <- length(free_shk_idx)

      rk <- qr(D_c)$rank
      if (rk < n_c)
        stop(sprintf(
          "conditional_forecast (unanticipated): at horizon %d, rank(%d) < %d conditions -- infeasible with chosen free_shocks.",
          h, rk, n_c), call. = FALSE)

      D_w <- sweep(D_c, 2L, w_free, "*")
      DDt <- tcrossprod(D_w)
      DDt <- (DDt + t(DDt)) * 0.5
      ## See the note at the condition-system solve above: DDt loses rank when
      ## a conditioning shock is switched off (w_free = 0) or the conditions
      ## outnumber the free shocks. Minimum-norm, not a ridge.
      DDt_inv <- .safe_inv(DDt, warn_label = "conditional_forecast: shock system")
      eps_free <- w_free * as.numeric(t(D_w) %*% (DDt_inv %*% b_h))
      eps_h[free_shk_idx] <- eps_free
    }

    y_h <- y_pred_mean + as.numeric(D_mat %*% eps_h)
    paths_point[h, ] <- y_h
    shock_paths[h, ] <- eps_h

    ## Propagate state
    s_prev <- as.numeric(T_mat %*% s_prev) + as.numeric(R_mat %*% eps_h)
  }

  ## Draws: same structure but with noise on unconstrained shock components
  draws_out <- NULL
  if (n_draws > 0L) {
    draws_list <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      s_d   <- s0
      paths_d <- matrix(0, H, n_obs)
      for (h in seq_len(H)) {
        y_pred_mean_d <- as.numeric(Z_mat %*% s_d)
        h_cond <- cond_df[cond_df$horizon == h, , drop = FALSE]

        ## Draw all shocks from N(0, Q)
        eps_d <- w_all * rnorm(n_shk)

        if (nrow(h_cond) > 0L) {
          cond_idx <- h_cond$var_idx
          b_h      <- h_cond$value - y_pred_mean_d[cond_idx]

          ## Project out in whitened eta space: set free shocks to satisfy
          ## conditions, keep null-space component for the draw
          D_c <- D_mat[cond_idx, free_shk_idx, drop = FALSE]
          D_w <- sweep(D_c, 2L, w_free, "*")
          DDt <- tcrossprod(D_w)
          DDt_inv <- .safe_inv((DDt + t(DDt)) * 0.5,
                               warn_label = "conditional_forecast: shock system")
          ## Particular solution (whitened)
          eta_particular <- t(D_w) %*% (DDt_inv %*% b_h)
          ## Null-space component (keeps random variation on free shocks
          ## orthogonal to the conditions)
          Proj <- t(D_w) %*% DDt_inv %*% D_w   ## n_free x n_free
          eta_free <- ifelse(w_free > 0, eps_d[free_shk_idx] / w_free, 0)
          eta_null <- (diag(length(free_shk_idx)) - Proj) %*% eta_free
          eps_d[free_shk_idx] <- w_free * as.numeric(eta_particular + eta_null)
        }

        y_d <- y_pred_mean_d + as.numeric(D_mat %*% eps_d)
        paths_d[h, ] <- y_d
        s_d <- as.numeric(T_mat %*% s_d) + as.numeric(R_mat %*% eps_d)
      }
      colnames(paths_d) <- ss$obs_names_fcst
      draws_list[[i]] <- paths_d
    }
    draws_out <- draws_list
  }

  list(paths_point = paths_point, shock_paths = shock_paths, draws = draws_out)
}


## ---------------------------------------------------------------------------
## Internal: soft conditioning (forward KF with pseudo-observations)
##
## Treats conditions as noisy observations of the forecast-period observables.
## Runs kalman_smoother over the forecast horizon with:
##   - Y_fcst: H x n_obs matrix, NA everywhere except conditioned (var, horizon) pairs
##   - me_extra: n_obs x H extra measurement-error variances (stderr^2 at
##     conditioned pairs, 0 elsewhere)
##
## This is exactly the filter_tunes pattern applied forward.
##
## @param Q  n_shock x n_shock shock covariance (NULL = identity; correct only
##   when every shock has stderr 1, since ghu holds unit-shock responses)
## ---------------------------------------------------------------------------
.soft_forecast <- function(ss, s0, P0, H, cond_df, n_draws, Q = NULL) {
  n_obs <- ss$n_obs
  n_shk <- ss$n_shock
  n_s   <- ss$n_state

  ## Build Y_fcst and me_extra
  Y_fcst   <- matrix(NA_real_, H, n_obs)
  me_extra <- matrix(0, n_obs, H)

  for (i in seq_len(nrow(cond_df))) {
    h   <- cond_df$horizon[i]
    idx <- cond_df$var_idx[i]
    Y_fcst[h, idx]   <- cond_df$value[i]
    se <- cond_df$stderr[i]
    if (!is.na(se) && se > 0) {
      me_extra[idx, h] <- se^2
    }
    ## Hard soft conditions: stderr = 0 -> very small variance for numerical stability
    if (is.na(se) || se == 0) me_extra[idx, h] <- 0
  }

  ## Run a forward KF over the forecast horizon with supplied initial state.
  ## We need to replicate the KF forward pass from kalman_smoother but
  ## starting from s0 / P0 (not from the Lyapunov steady state).
  T_mat  <- ss$T_mat
  R_mat  <- ss$R_mat
  Z_mat  <- ss$Z_mat
  D_mat  <- ss$D_mat
  if (is.null(Q)) Q <- diag(n_shk)

  RQR <- R_mat %*% Q %*% t(R_mat)
  DQD <- D_mat %*% Q %*% t(D_mat)
  RQD <- R_mat %*% Q %*% t(D_mat)

  has_me_extra <- any(me_extra != 0)

  s_filt <- matrix(0, H, n_s)
  P_filt <- array(0, c(n_s, n_s, H))
  s_pred <- matrix(0, H, n_s)
  P_pred <- array(0, c(n_s, n_s, H))

  s_tt <- s0
  P_tt <- P0

  for (h in seq_len(H)) {
    ## Predict
    s_tp <- as.numeric(T_mat %*% s_tt)
    P_tp <- T_mat %*% P_tt %*% t(T_mat) + RQR
    s_pred[h, ]   <- s_tp
    P_pred[, , h] <- P_tp

    ## Innovation (against s_{h-1} = s_tt, lagged-state convention)
    y_pred_h <- as.numeric(Z_mat %*% s_tt)
    v_h      <- Y_fcst[h, ] - y_pred_h

    obs_ok <- which(!is.na(v_h))
    n_ok   <- length(obs_ok)

    if (n_ok == 0L) {
      s_tt <- s_tp
      P_tt <- P_tp
    } else {
      v_ok  <- v_h[obs_ok]
      Zt    <- Z_mat[obs_ok, , drop = FALSE]
      DQDt  <- tcrossprod(D_mat[obs_ok, , drop = FALSE] %*% Q,
                          D_mat[obs_ok, , drop = FALSE])
      RQDt  <- RQD[, obs_ok, drop = FALSE]

      F_t <- Zt %*% P_tt %*% t(Zt) + DQDt
      if (has_me_extra) diag(F_t) <- diag(F_t) + me_extra[obs_ok, h]
      F_t <- (F_t + t(F_t)) * 0.5

      ## Singular F: drop the zero-variance components, do NOT ridge them.
      ## This used to read
      ##     F_ch <- tryCatch(chol(F_t + 1e-10 * diag(n_ok)), ...)
      ##     if (is.null(F_ch)) F_ch <- chol(F_t + 1e-6 * diag(n_ok))
      ## which carried three problems. The 1e-10 ridge was added on EVERY
      ## period, so this pass was never an unregularised Kalman filter; both
      ## ridges were ABSOLUTE while F_t is not O(1) (a unit-root state space
      ## with `shock_scale` puts it many orders higher); and the fallback was
      ## UNGUARDED, so when 1e-6 was not enough it threw. Same defect as the
      ## one fixed in kalman_smoother() -- see test-smoother-singular-F.R.
      F_ch <- tryCatch(chol(F_t), error = function(e) NULL)
      if (is.null(F_ch)) {
        keep <- .smoother_informative_obs(F_t)
        if (!any(keep)) {
          ## No component carries information this period: predict-only.
          s_tt <- s_tp
          P_tt <- P_tp
          s_filt[h, ]   <- s_tt
          P_filt[, , h] <- P_tt
          next
        }
        obs_ok <- obs_ok[keep]
        n_ok   <- length(obs_ok)
        v_ok   <- v_ok[keep]
        Zt     <- Z_mat[obs_ok, , drop = FALSE]
        DQDt   <- tcrossprod(D_mat[obs_ok, , drop = FALSE] %*% Q,
                             D_mat[obs_ok, , drop = FALSE])
        RQDt   <- RQD[, obs_ok, drop = FALSE]
        F_t    <- Zt %*% P_tt %*% t(Zt) + DQDt
        if (has_me_extra) diag(F_t) <- diag(F_t) + me_extra[obs_ok, h]
        F_t    <- (F_t + t(F_t)) * 0.5
        F_ch   <- tryCatch(chol(F_t), error = function(e) NULL)
        if (is.null(F_ch))
          stop("conditional_forecast: the innovation covariance at horizon ",
               h, " is not positive definite even after dropping every ",
               "zero-variance observation component.", call. = FALSE)
      }
      F_inv <- chol2inv(F_ch)

      K_t  <- (T_mat %*% P_tt %*% t(Zt) + RQDt) %*% F_inv
      s_tt <- s_tp + as.numeric(K_t %*% v_ok)
      P_tt <- P_tp - K_t %*% F_t %*% t(K_t)
      P_tt <- 0.5 * (P_tt + t(P_tt))
    }

    s_filt[h, ]   <- s_tt
    P_filt[, , h] <- P_tt
  }

  ## Recover paths from filtered states: the filtered state at h gives
  ## s_{h|conds}, and y_h = Z s_{h-1|conds}.
  paths_point <- matrix(0, H, n_obs)
  s_cur <- s0
  for (h in seq_len(H)) {
    paths_point[h, ] <- as.numeric(Z_mat %*% s_cur)
    s_cur <- s_filt[h, ]
  }
  colnames(paths_point) <- ss$obs_names_fcst

  ## Shock paths: recover from filtered states
  R_pinv <- MASS::ginv(R_mat)
  shock_paths <- matrix(0, H, n_shk)
  colnames(shock_paths) <- ss$shock_names
  s_cur2 <- s0
  for (h in seq_len(H)) {
    residual <- s_filt[h, ] - as.numeric(T_mat %*% s_cur2)
    shock_paths[h, ] <- as.numeric(R_pinv %*% residual)
    s_cur2 <- s_filt[h, ]
  }

  ## Draws via simulation smoother (simplified: draw from P_filt conditional)
  draws_out <- NULL
  if (n_draws > 0L) {
    draws_list <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      ## Simple: perturb the filtered mean with conditional covariance
      paths_d <- matrix(0, H, n_obs)
      s_d <- s0
      s_d_prev <- s0
      for (h in seq_len(H)) {
        ## Draw state perturbation from P_filt[,,h]
        P_h  <- P_filt[, , h]
        sv_h <- svd(P_h)
        sv_h$d <- pmax(sv_h$d, 0)
        Prt  <- sv_h$u %*% diag(sqrt(sv_h$d), nrow = length(sv_h$d)) %*% t(sv_h$u)
        s_draw_h <- s_filt[h, ] + as.numeric(Prt %*% rnorm(n_s))
        paths_d[h, ] <- as.numeric(Z_mat %*% s_d_prev)
        s_d_prev <- s_draw_h
      }
      colnames(paths_d) <- ss$obs_names_fcst
      draws_list[[i]] <- paths_d
    }
    draws_out <- draws_list
  }

  list(paths_point = paths_point, shock_paths = shock_paths, draws = draws_out)
}


## ---------------------------------------------------------------------------
## Internal: unconditional forecast (plain propagation from s0)
## Used for the no-conditions case and internal consistency checks.
## ---------------------------------------------------------------------------
.unconditional_forecast <- function(ss, s0, H) {
  T_mat <- ss$T_mat
  Z_mat <- ss$Z_mat
  n_obs <- ss$n_obs

  paths <- matrix(0, H, n_obs)
  colnames(paths) <- ss$obs_names_fcst
  s_prev <- s0
  for (h in seq_len(H)) {
    paths[h, ] <- as.numeric(Z_mat %*% s_prev)
    s_prev <- as.numeric(T_mat %*% s_prev)
  }
  paths
}


## ---------------------------------------------------------------------------
## Internal: extract terminal state s_{T|T} and P_{T|T} from Y
##
## Dispatch table on ctx$likelihood:
##   NULL / "gaussian" / "whittle" -> standard Kalman smoother (bit-identical
##     to the original code path).
##   "tpf"  -> Tempered Particle Filter; particle mean (x1 + x2) at T.
##             Requires dr to be a DecisionRules2 object and ctx$me_variance > 0.
##   "pskf" -> PSKF smoother; smoothed mean at T via Gaussian RTS backward pass.
##   "pkf"  -> OBC PKF forward pass; filtered state s_{T|T} and P_{T|T}.
##
## @param data     T x n_obs observation matrix (rows = time).
## @param ss_raw   State-space list from build_dsge_state_space().
## @param Q        n_shock x n_shock shock covariance.
## @param ctx      estimation_context or NULL.
## @param model    Parsed model (for tpf/pskf/pkf paths).
## @param dr       Decision rules (DecisionRules or DecisionRules2).
## @param obs_vars  Character vector of observable names.
## @param compiled  Compiled model (dynhr_compiled) — required for pkf path.
## @return list(s0 = n_state numeric, P0 = n_state x n_state matrix).
## ---------------------------------------------------------------------------
.extract_terminal_state <- function(data, ss_raw, Q, ctx, model, dr, obs_vars,
                                    compiled = NULL) {

  lik <- if (is.null(ctx)) "gaussian" else ctx$likelihood

  ## Dispatch table — "pkf", "ppf", "copf" added for OBC models (Tier 10 item 5;
  ## Tier 15 §C: ppf/copf keys added so particle-filter-estimated OBC models use
  ## the correct terminal-state summary rather than silently falling back to pkf).
  dispatch <- list(
    gaussian = "gaussian",
    whittle  = "gaussian",   ## Whittle uses the same smoother for s0
    tpf      = "tpf",
    pskf     = "pskf",
    pkf      = "pkf",
    ppf      = "ppf",
    copf     = "copf"
  )
  path <- dispatch[[lik]]
  if (is.null(path))
    stop(
      "conditional_forecast: unrecognised ctx$likelihood '", lik, "'. ",
      "Valid values: gaussian, whittle, tpf, pskf, pkf, ppf, copf.",
      call. = FALSE
    )

  ## ---- Gaussian / Whittle path (original code) ----------------------------
  if (identical(path, "gaussian")) {
    ## Levels in: `ss_raw$d` is the observation intercept and the smoother
    ## subtracts it, the same demeaning the tpf and pskf paths below do
    ## explicitly. Before 0.9.3 this branch alone took deviations, so the same
    ## conditional_forecast() call needed different data depending on
    ## ctx$likelihood.
    sm <- .kalman_smoother_ss(data, ss_raw, Q = Q)
    s0 <- as.numeric(sm$filtered_states[nrow(sm$filtered_states), ])
    P0 <- sm$P_filt_last
    return(list(s0 = s0, P0 = P0))
  }

  ## ---- TPF path -----------------------------------------------------------
  if (identical(path, "tpf")) {
    if (!inherits(dr, "DecisionRules2"))
      stop(
        "conditional_forecast (ctx$likelihood = \"tpf\"): the TPF terminal-state ",
        "path requires a second-order decision-rules object.\n",
        "Call solve_perturbation_order2() and pass the result as 'dr'.",
        call. = FALSE
      )
    me_var <- ctx$me_variance
    if (!is.numeric(me_var) || length(me_var) != 1L ||
        !is.finite(me_var) || me_var <= 0)
      stop(
        "conditional_forecast (ctx$likelihood = \"tpf\"): ctx$me_variance must ",
        "be a finite positive scalar (the TPF tempering instrument).",
        call. = FALSE
      )

    ## Merge tpf_options (ctx wins over defaults)
    tpf_opts <- if (!is.null(ctx$tpf_options) && length(ctx$tpf_options) > 0L)
      ctx$tpf_options else list()
    n_particles <- tpf_opts$n_particles %||% 500L
    ess_target  <- tpf_opts$ess_target  %||% 0.5
    n_mh        <- tpf_opts$n_mh        %||% 1L
    ## No mh_scale: removed from the TPF stack 2026-09-02 (inert since the
    ## mutation step was corrected to hold the ancestor state fixed).
    tpf_seed    <- tpf_opts$seed        %||% NULL

    ## State-space from the order-2 DR
    state_idx <- dr$state_idx
    obs_idx   <- match(obs_vars, dr$endo_names)
    if (any(is.na(obs_idx)))
      stop("conditional_forecast (tpf): some obs_names not found in dr$endo_names.",
           call. = FALSE)

    ZZ       <- dr$ghx[obs_idx,   , drop = FALSE]
    DD       <- dr$ghu[obs_idx,   , drop = FALSE]
    d_obs    <- dr$ys[obs_vars]
    ghss_obs <- 0.5 * dr$ghss[obs_idx]

    n_s  <- length(state_idx)
    n_2s <- 2L * n_s
    N    <- as.integer(n_particles)

    ## Cholesky of Sigma_e
    Sigma_e <- ss_raw$Sigma_e
    L_e <- tryCatch(t(chol(Sigma_e)), error = function(e) {
      eg   <- eigen(Sigma_e, symmetric = TRUE)
      vals <- pmax(eg$values, 0)
      eg$vectors %*% diag(sqrt(vals), nrow = length(vals))
    })

    ## Lyapunov P0 for initial particle cloud
    TT_s <- dr$ghx[state_idx, , drop = FALSE]
    RR_s <- dr$ghu[state_idx, , drop = FALSE]
    QQ_s <- tcrossprod(RR_s %*% Sigma_e, RR_s)
    P0_tpf <- tryCatch({
      Pk <- QQ_s
      for (iter in seq_len(500L)) {
        Pk_new <- TT_s %*% Pk %*% t(TT_s) + QQ_s
        if (max(abs(Pk_new - Pk)) < 1e-12 * (1 + max(abs(Pk_new)))) break
        Pk <- Pk_new
      }
      Pk
    }, error = function(e) QQ_s)

    if (!is.null(tpf_seed)) set.seed(tpf_seed)

    if (is.null(P0_tpf) || !is.finite(max(abs(P0_tpf)))) {
      particles <- matrix(0, nrow = n_2s, ncol = N)
    } else {
      L_P0    <- tryCatch(t(chol(P0_tpf + diag(1e-12, n_s))),
                          error = function(e) diag(sqrt(diag(P0_tpf) + 1e-12), n_s))
      x1_init <- L_P0 %*% matrix(rnorm(n_s * N), nrow = n_s)
      particles <- rbind(x1_init, matrix(0, nrow = n_s, ncol = N))
    }

    ## Run TPF forward pass over all T periods; keep final particles
    Y_mat <- t(data)    ## n_obs x T (columns = periods)
    T_obs <- ncol(Y_mat)
    for (t in seq_len(T_obs)) {
      y_t <- Y_mat[, t]
      if (any(!is.finite(y_t))) next
      res <- tpf_run_period(
        particles   = particles,
        y_t         = y_t,
        dr2         = dr,
        Sigma_e     = Sigma_e,
        L_e         = L_e,
        ZZ          = ZZ,
        DD          = DD,
        d_obs       = d_obs,
        ghss_obs    = ghss_obs,
        me_variance = me_var,
        ess_target  = ess_target,
        n_mh        = n_mh,
        use_rcpp    = .HAS_RCPP_TPF(),
        backend     = if (.HAS_RCPP_TPF_PERIOD()) "cpp" else "R"
      )
      particles <- res$particles
    }

    ## Particle mean: x1 + x2 components (both n_s x N)
    x1_final <- particles[seq_len(n_s), , drop = FALSE]
    x2_final <- particles[seq_len(n_s) + n_s, , drop = FALSE]
    s0 <- as.numeric(rowMeans(x1_final + x2_final))

    ## Covariance from particle cloud
    state_cloud <- x1_final + x2_final   ## n_s x N
    cloud_centred <- state_cloud - s0
    P0 <- tcrossprod(cloud_centred) / (N - 1L)

    return(list(s0 = s0, P0 = P0))
  }

  ## ---- PSKF path ----------------------------------------------------------
  if (identical(path, "pskf")) {
    me_var <- if (!is.null(ctx$me_variance)) ctx$me_variance else 0

    state_idx <- dr$state_idx
    obs_idx   <- match(obs_vars, dr$endo_names)
    if (any(is.na(obs_idx)))
      stop("conditional_forecast (pskf): some obs_names not found in dr$endo_names.",
           call. = FALSE)

    TT <- dr$ghx[state_idx, , drop = FALSE]
    ZZ <- dr$ghx[obs_idx,   , drop = FALSE]

    ## CSN shock parameters
    exo_names <- dr$exo_names
    csn <- tryCatch(
      .get_csn_shock_params(model, exo_names, obs_vars, dr,
                            model$param_values, me_var),
      error = function(e) {
        stop(
          "conditional_forecast (pskf): failed to assemble CSN shock parameters.\n",
          "Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )

    ## Demean Y (PSKF works on deviations from SS)
    d_obs <- dr$ys[obs_vars]
    Y_dm  <- t(data) - d_obs   ## n_obs x T

    sm <- tryCatch(
      pskf_smoother(
        Y         = Y_dm,
        TT        = TT,
        ZZ        = ZZ,
        mu_eta    = csn$mu_eta,
        Sigma_eta = csn$Sigma_eta,
        Gamma_eta = csn$Gamma_eta,
        nu_eta    = csn$nu_eta,
        Delta_eta = csn$Delta_eta,
        mu_eps    = csn$mu_eps,
        Sigma_eps = csn$Sigma_eps
      ),
      error = function(e) {
        stop(
          "conditional_forecast (pskf): pskf_smoother() failed.\n",
          "Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )

    ## PSKF timing: its observation equation is y_t = ZZ x_t (current-state),
    ## while the forecaster consumes s_T in the lagged convention
    ## (y_{T+1} = Z s_T).  Matching E[y_{T+1} | y_{1:T}] across the two
    ## (observationally equivalent) representations requires one transition:
    ## s0 = TT x_T (E[eta] = 0 by mu_eta construction) and
    ## P0 = TT P_T TT' + Sigma_eta.
    T_obs <- nrow(sm$smoothed_means)
    xT    <- as.numeric(sm$smoothed_means[T_obs, ])
    PT    <- sm$smoothed_covs[, , T_obs]
    s0    <- as.numeric(TT %*% xT)
    P0    <- TT %*% PT %*% t(TT) + csn$Sigma_eta
    return(list(s0 = s0, P0 = P0))
  }

  ## ---- OBC PKF path -------------------------------------------------------
  if (identical(path, "pkf")) {
    ## Need compiled to call cache_system_structure -> extract_system_matrices_fast
    ## -> .solve_from_system -> obc_ensure_policy.  Either passed via compiled= or
    ## available from model (re-compile as fallback; ~50 ms overhead).
    if (is.null(compiled)) {
      message(
        "conditional_forecast (pkf): 'compiled' not supplied; re-compiling model.\n",
        "Pass compiled = solved$compiled to avoid this overhead."
      )
      compiled <- compile_model(model, verbose = FALSE)
    }

    me_var <- ctx$me_variance %||% 1e-8

    ## OBC specs: use pre-parsed ones from ctx if available, else parse now.
    specs <- ctx$obc_specs %||% obc_parse_tags(model)

    obs_idx <- match(obs_vars, model$var_names)
    if (any(is.na(obs_idx)))
      stop("conditional_forecast (pkf): some obs_names not found in model$var_names.",
           call. = FALSE)

    ## Ensure the compiled object has lead_lag_incidence attached.
    if (is.null(compiled$lead_lag_incidence) &&
        !is.null(compiled$model$lead_lag_incidence))
      compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

    ## Solve steady state and extract system matrices at current param values.
    params <- model$param_values
    ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
    if (is.null(ss_result) || !ss_result$converged)
      stop("conditional_forecast (pkf): steady-state solver did not converge.",
           call. = FALSE)

    sys_cache <- cache_system_structure(compiled)
    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys       <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)

    ## Solve for the slack-regime decision rules.
    dr_slack <- tryCatch(
      .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE),
      error = function(e) stop(
        "conditional_forecast (pkf): perturbation solver failed.\n",
        "Original error: ", conditionMessage(e), call. = FALSE
      )
    )
    if (is.null(dr_slack) || !isTRUE(dr_slack$bk_satisfied))
      stop("conditional_forecast (pkf): BK condition not satisfied.", call. = FALSE)

    ## me-floor hazard guard (see R/obc-regime.R .obc_warn_me_floor_lock() and
    ## R/pruned-state-space.R); single-shot call site, no closure latch needed.
    .obc_warn_me_floor_lock(
      dr_slack, model, params, obs_vars, obs_idx, me_var,
      check = isTRUE(getOption("dynhr.me_floor_check", TRUE)))

    ## Seed regime cache with the slack policy.
    regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
    obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)

    ## Ensure Y is n_obs x T.
    Y_mat <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else
              if (nrow(data) == length(obs_vars)) data else t(data)

    ## PKF forward pass: collect filtered state and P_{T|T}.
    kf <- kalman_filter_obc_pkf(
      Y_mat, dr_slack, regime_cache, sys,
      model, params, obs_vars, specs,
      obs_idx         = obs_idx,
      me_variance     = me_var,
      return_filtered = TRUE,
      return_P_last   = TRUE,
      return_store    = FALSE
    )

    if (is.null(kf) || !is.finite(kf$loglik))
      stop("conditional_forecast (pkf): PKF forward pass returned non-finite loglik.",
           call. = FALSE)

    n_T <- ncol(Y_mat)
    s0  <- as.numeric(kf$filtered_states[, n_T])
    P0  <- kf$P_last %||% diag(1e-6, length(s0))

    return(list(s0 = s0, P0 = P0))
  }

  ## ---- OBC PPF / COPF path (Tier 15 §C) ----------------------------------
  ## Both "ppf" and "copf" summarise the terminal particle cloud to (mean, cov)
  ## exactly as the TPF path does, giving a Gaussian summary for the downstream
  ## Waggoner-Zha / soft-KF forecaster (which is purely Gaussian-state).
  ##
  ## The difference between "ppf" and "copf" is only the proposal argument
  ## forwarded to ppf_likelihood(); the terminal-state summary is identical.
  if (path %in% c("ppf", "copf")) {
    if (is.null(compiled)) {
      message(
        "conditional_forecast (", path, "): 'compiled' not supplied; re-compiling model.\n",
        "Pass compiled = solved$compiled to avoid this overhead."
      )
      compiled <- compile_model(model, verbose = FALSE)
    }

    me_var <- ctx$me_variance %||% 1e-4

    ## OBC specs
    specs <- ctx$obc_specs %||% obc_parse_tags(model)

    obs_idx <- match(obs_vars, model$var_names)
    if (any(is.na(obs_idx)))
      stop("conditional_forecast (", path, "): some obs_names not found in model$var_names.",
           call. = FALSE)

    ## Ensure compiled has lead_lag_incidence
    if (is.null(compiled$lead_lag_incidence) &&
        !is.null(compiled$model$lead_lag_incidence))
      compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

    ## Steady state + system matrices at current param values
    params <- model$param_values
    ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
    if (is.null(ss_result) || !ss_result$converged)
      stop("conditional_forecast (", path, "): steady-state solver did not converge.",
           call. = FALSE)

    sys_cache <- cache_system_structure(compiled)
    params    <- ss_result$params %||% params
    sys       <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)

    dr_slack <- tryCatch(
      .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE),
      error = function(e) stop(
        "conditional_forecast (", path, "): perturbation solver failed.\n",
        "Original error: ", conditionMessage(e), call. = FALSE
      )
    )
    if (is.null(dr_slack) || !isTRUE(dr_slack$bk_satisfied))
      stop("conditional_forecast (", path, "): BK condition not satisfied.", call. = FALSE)

    ## Read PPF/COPF options from ctx$obc_options (with safe defaults)
    obc_opts     <- ctx$obc_options %||% list()
    N_particles  <- as.integer(obc_opts$N_particles %||% 2000L)
    proposal_arg <- if (path == "copf") "copf" else "bootstrap"
    regime_guess <- obc_opts$regime_guess %||% "ancestor"
    ppf_seed     <- obc_opts$seed %||% NULL

    ## Build a FRESH regime_cache (theta-dependent matrices at current params)
    regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
    ## For COPF: obc_ensure_policy needs copf_args; we add them after Sigma_e inversion.
    ## For bootstrap: no copf_args needed (NULL).
    copf_args <- NULL
    if (path == "copf") {
      Sigma_e_ppf <- .get_shock_cov(model, dr_slack$exo_names, params)
      ch_Se       <- tryCatch(chol(Sigma_e_ppf), error = function(e2) NULL)
      if (!is.null(ch_Se))
        copf_args <- list(Sigma_e     = Sigma_e_ppf,
                          Sigma_e_inv = chol2inv(ch_Se),
                          me_variance = me_var)
    }
    obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx, copf_args)

    ## Ensure Y is n_obs x T
    Y_mat <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else
              if (nrow(data) == length(obs_vars)) data else t(data)

    ## Run PPF/COPF forward pass over the full history, returning the terminal
    ## particle cloud.
    pf <- ppf_likelihood(
      Y_mat, dr_slack, regime_cache, sys,
      model, params, obs_vars, specs,
      obs_idx          = obs_idx,
      N                = N_particles,
      me_variance      = me_var,
      proposal         = proposal_arg,
      regime_guess     = regime_guess,
      return_particles = TRUE,
      seed             = ppf_seed
    )

    if (is.null(pf) || !is.finite(pf$loglik))
      stop("conditional_forecast (", path, "): PPF forward pass returned non-finite loglik.",
           call. = FALSE)

    particles <- pf$particles   ## n_state x N_particles
    N_p       <- ncol(particles)

    ## Summarise terminal cloud to (mean, empirical cov) — same as TPF path.
    s0        <- rowMeans(particles)
    cloud_c   <- particles - s0
    P0        <- tcrossprod(cloud_c) / (N_p - 1L)

    return(list(s0 = s0, P0 = P0))
  }

  stop("conditional_forecast: unhandled likelihood path '", path, "'.", call. = FALSE)
}


## ---------------------------------------------------------------------------
#' Waggoner-Zha (1999) conditional forecasts
#'
#' Computes conditional forecasts for a DSGE model by imposing linear
#' restrictions on future observables, finding the (minimum-norm) shock path
#' that delivers the conditions.
#'
#' @param model       Parsed model object from \code{parse_mod()}.
#' @param dr          Decision rules from \code{solve_perturbation()}.
#' @param data        \code{T x n_obs} matrix of observables in
#'   \strong{levels}, used to anchor the forecast (the filtered terminal
#'   state \code{s_{T|T}} is extracted internally, subtracting the model's
#'   steady state as \code{\link{kalman_filter}()} does). May also be a list
#'   output from \code{\link{kalman_smoother}()} (then its
#'   \code{filtered_states} and \code{P_filt_last} fields are used directly).
#'   Note that the FORECAST PATHS and the \code{conditions} remain in
#'   deviations from steady state: only the anchoring data is in levels, and
#'   only because every filtering entry point now is.
#' @param conditions  Data frame with columns:
#'   \describe{
#'     \item{\code{var}}{Character: observable name.}
#'     \item{\code{horizon}}{Integer: forecast horizon (1 = one period ahead).}
#'     \item{\code{value}}{Numeric: target value (in model deviation units).}
#'     \item{\code{stderr}}{Numeric (optional, \code{NA} or absent for hard
#'       conditions): standard deviation of the conditioning noise for soft
#'       conditions.}
#'   }
#'   Pass \code{NULL} or a zero-row data frame for an unconditional forecast.
#' @param horizon     Integer: total forecast horizon H.
#' @param obs_vars    Character vector: names of observables in \code{data}
#'   (column order). Defaults to \code{colnames(data)} if available.
#' @param type        \code{"anticipated"} (default): the whole shock path is
#'   chosen at time T (Waggoner-Zha stacked solution). \code{"unanticipated"}:
#'   agents are surprised each period; the conditioning shock is chosen
#'   period-by-period.
#' @param method      \code{"hard"} (default): conditions are imposed exactly
#'   (zero measurement error). \code{"soft"}: conditions are treated as noisy
#'   pseudo-observations with measurement error \code{stderr^2}.
#' @param n_draws     Integer: number of draws from the conditional shock
#'   distribution (default 0 = point forecast only).
#' @param free_shocks Character vector (or \code{NULL}): names of shocks allowed
#'   to move to satisfy the conditions. Defaults to all shocks. Restricting
#'   \code{free_shocks} implements the Waggoner-Zha "attribution" practice
#'   where only named shocks are assigned to deliver the conditions.
#' @param Q           \code{n_shock x n_shock} shock covariance matrix (default
#'   \code{NULL}: use \code{Sigma_e} from the model's \code{shocks;} block, as
#'   returned in \code{build_dsge_state_space()$Sigma_e}; \code{ghu} holds
#'   unit-shock responses, so the identity would be correct only when every
#'   shock has \code{stderr 1}).
#' @param ctx         Optional \code{\link{estimation_context}} object. When
#'   supplied, dispatches the terminal-state extraction on
#'   \code{ctx$likelihood}: \code{"gaussian"} and \code{"whittle"} use the
#'   standard Kalman smoother (default, bit-identical to \code{ctx = NULL});
#'   \code{"tpf"} runs the Tempered Particle Filter (requires \code{dr} to be
#'   a \code{DecisionRules2} object from \code{solve_perturbation_order2()} and
#'   \code{ctx$me_variance > 0}); \code{"pskf"} uses the PSKF smoother
#'   (smoothed mean at T).
#' @param compiled    Optional \code{dynhr_compiled} from
#'   \code{\link{compile_model}}.  Required when \code{ctx$likelihood = "pkf"}
#'   (OBC piecewise-linear filter); if \code{NULL} the model is recompiled
#'   on-the-fly (slower).
#' @param ...         Currently unused.
#'
#' @return An object of class \code{"dynhr_cfcst"} with:
#'   \describe{
#'     \item{\code{paths_point}}{H x n_obs matrix of point-forecast observable paths.}
#'     \item{\code{paths_uncond}}{H x n_obs unconditional forecast paths.}
#'     \item{\code{shock_paths}}{H x n_shock matrix of implied shock paths.}
#'     \item{\code{draws}}{List of H x n_obs matrices (length \code{n_draws}),
#'       or \code{NULL} if \code{n_draws = 0}.}
#'     \item{\code{conditions}}{The \code{conditions} data frame (copy).}
#'     \item{\code{horizon}}{H.}
#'     \item{\code{type}, \code{method}}{As supplied.}
#'     \item{\code{obs_names}, \code{shock_names}}{Character vectors.}
#'   }
#'
#' @param plan Optional \code{\link{dynhr_plan}} carrying out-of-sample
#'   conditions (see \code{\link{plan_condition}}); mutually exclusive with
#'   \code{conditions}/\code{type}/\code{method}.
#'
#' @references
#'   Waggoner, D. F. & Zha, T. (1999). Conditional Forecasts in Dynamic
#'     Multivariate Models. \emph{Review of Economics and Statistics}, 81(4),
#'     639-651.
#'
#' @seealso \code{\link{plan_condition}}, \code{\link{bayesian_conditional_forecast}},
#'   \code{\link{kalman_smoother}}
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
#' set.seed(1)
#' Y <- as.matrix(simulate_model(dr, n_periods = 100L, model = model,
#'                               burn_in = 20L)[, "y", drop = FALSE])
#'
#' ## Hard-condition output growth for the first two forecast quarters
#' conds <- data.frame(var = "y", horizon = 1:2, value = c(0.01, 0.008))
#' cf <- conditional_forecast(model, dr, Y, conditions = conds,
#'                            horizon = 8L, obs_vars = "y")
#' cf
#' head(cf$shock_paths)        # the shock path that delivers the conditions
#'
#' @export
## ---------------------------------------------------------------------------
conditional_forecast <- function(model, dr, data, conditions = NULL,
                                 horizon = 8L,
                                 obs_vars = NULL,
                                 type   = c("anticipated", "unanticipated"),
                                 method = c("hard", "soft"),
                                 n_draws = 0L,
                                 free_shocks = NULL,
                                 Q = NULL,
                                 plan = NULL,
                                 ctx = NULL,
                                 compiled = NULL,
                                 ...) {
  ## ---- Resolve plan= if supplied ----
  if (!is.null(plan)) {
    if (!inherits(plan, "dynhr_plan"))
      stop("conditional_forecast: 'plan' must be a dynhr_plan object.", call. = FALSE)
    if (!is.null(conditions))
      stop("conditional_forecast: supply either 'plan' or 'conditions', not both.",
           call. = FALSE)
    ## Warn if plan has shock_scale entries (out-of-sample scaling unsupported)
    if (length(plan$shock_scales) > 0L)
      warning(
        "conditional_forecast: plan contains shock_scale entries, but out-of-sample ",
        "per-horizon shock scaling is not supported in this version. The shock_scale ",
        "entries are ignored on the forecast side. Pass a custom Q matrix to ",
        "conditional_forecast() to control forecast-side shock covariance.",
        call. = FALSE)
    cond_df <- plan_to_conditions(plan)
    conditions <- cond_df
    ## Override type and method from plan attributes if conditions are present
    if (nrow(cond_df) > 0L) {
      type   <- attr(cond_df, "type")   %||% type[1L]
      method <- attr(cond_df, "method") %||% method[1L]
    }
  }

  type   <- match.arg(type)
  method <- match.arg(method)
  H      <- as.integer(horizon)
  stopifnot(H >= 1L)

  ## ---- Resolve observable names ----
  ## Y may be a raw matrix/data.frame or a kalman_smoother() result list.
  smoother_result <- NULL
  if (is.list(data) && !is.data.frame(data) && !is.null(data$filtered_states)) {
    smoother_result <- data
    data <- NULL
    if (is.null(obs_vars))
      stop("conditional_forecast: when Y is a kalman_smoother() result, supply obs_names explicitly.",
           call. = FALSE)
  } else {
    data <- as.matrix(data)
    if (is.null(obs_vars)) {
      if (!is.null(colnames(data))) obs_vars <- colnames(data)
      else stop("conditional_forecast: obs_names must be supplied when Y has no colnames.",
                call. = FALSE)
    }
  }

  ## ---- Build state space ----
  ss_raw <- build_dsge_state_space(model, dr, obs_vars, verbose = FALSE)
  ## Attach obs_names for internal use
  ss_raw$obs_names_fcst <- obs_vars

  ## Default shock covariance: Sigma_e from the shocks; block (ghu holds
  ## unit-shock responses). Threaded through the terminal-state smoother run
  ## and every conditioning path (anticipated, unanticipated, soft).
  if (is.null(Q)) Q <- ss_raw$Sigma_e

  ## ---- Extract terminal state s_{T|T} and P_{T|T} ----
  if (!is.null(smoother_result)) {
    s0 <- as.numeric(smoother_result$filtered_states[nrow(smoother_result$filtered_states), ])
    P0 <- if (!is.null(smoother_result$P_filt_last)) smoother_result$P_filt_last
          else diag(ss_raw$n_state) * 1e-6   ## fallback: near-zero uncertainty
  } else {
    ## Dispatch terminal-state extraction on ctx$likelihood
    ts <- .extract_terminal_state(data, ss_raw, Q, ctx, model, dr, obs_vars,
                                  compiled = compiled)
    s0 <- ts$s0
    P0 <- ts$P0
  }

  ## ---- Unconditional forecast ----
  paths_uncond <- .unconditional_forecast(ss_raw, s0, H)

  ## ---- Handle null / empty conditions -> unconditional ----
  no_conds <- is.null(conditions) ||
    (is.data.frame(conditions) && nrow(conditions) == 0L)

  if (no_conds) {
    result <- structure(
      list(
        paths_point  = paths_uncond,
        paths_uncond = paths_uncond,
        shock_paths  = matrix(0, H, ss_raw$n_shock,
                              dimnames = list(NULL, ss_raw$shock_names)),
        draws        = NULL,
        conditions   = data.frame(var = character(0), horizon = integer(0),
                                  value = numeric(0), stringsAsFactors = FALSE),
        horizon      = H,
        type         = type,
        method       = method,
        obs_names    = obs_vars,
        shock_names  = ss_raw$shock_names
      ),
      class = "dynhr_cfcst"
    )
    return(result)
  }

  ## ---- Validate conditions ----
  if (!is.data.frame(conditions))
    stop("conditional_forecast: 'conditions' must be a data.frame.", call. = FALSE)
  req_cols <- c("var", "horizon", "value")
  miss <- setdiff(req_cols, names(conditions))
  if (length(miss))
    stop(sprintf("conditional_forecast: 'conditions' is missing columns: %s.",
                 paste(miss, collapse = ", ")), call. = FALSE)

  bad_var <- setdiff(conditions$var, obs_vars)
  if (length(bad_var))
    stop(sprintf("conditional_forecast: condition variable(s) not in obs_names: %s.",
                 paste(bad_var, collapse = ", ")), call. = FALSE)

  bad_h <- which(conditions$horizon < 1L | conditions$horizon > H)
  if (length(bad_h))
    stop(sprintf("conditional_forecast: condition horizon(s) out of range [1, %d]: %s.",
                 H, paste(conditions$horizon[bad_h], collapse = ", ")), call. = FALSE)

  ## Add var_idx and normalise stderr
  cond_df <- conditions
  cond_df$var_idx <- match(cond_df$var, obs_vars)
  if (!"stderr" %in% names(cond_df)) cond_df$stderr <- NA_real_
  cond_df$stderr <- as.numeric(cond_df$stderr)

  ## ---- free_shocks ----
  all_shocks <- ss_raw$shock_names
  if (is.null(free_shocks)) {
    free_shk_idx <- seq_along(all_shocks)
  } else {
    bad_shk <- setdiff(free_shocks, all_shocks)
    if (length(bad_shk))
      stop(sprintf("conditional_forecast: free_shocks not in model: %s.",
                   paste(bad_shk, collapse = ", ")), call. = FALSE)
    free_shk_idx <- match(free_shocks, all_shocks)
  }

  ## ---- Dispatch ----
  if (method == "soft") {
    ## Soft: use forward KF with pseudo-observations
    ## (type distinction is implicit: unanticipated path embedded in KF)
    res_inner <- .soft_forecast(ss_raw, s0, P0, H, cond_df, n_draws, Q = Q)

  } else {
    ## Hard conditioning
    if (type == "anticipated") {
      ## Build stacked matrices
      stk <- .build_forecast_stack(ss_raw, H)
      A_full <- stk$A
      M      <- stk$M

      ## Deterministic mean: M s0 for each obs-horizon pair
      mean_stack <- as.numeric(M %*% s0)   ## n_obs*H

      ## Build row index into A_full for each condition
      n_obs_ss <- ss_raw$n_obs
      cond_rows <- (cond_df$horizon - 1L) * n_obs_ss + cond_df$var_idx
      b_cond    <- cond_df$value - mean_stack[cond_rows]

      ## free_cols: columns of A_full corresponding to free shocks across all H
      free_cols <- as.integer(outer((free_shk_idx - 1L), (0L:(H - 1L)) * ss_raw$n_shock,
                                    "+") + 1L)
      ## Flatten to sorted vector
      free_cols <- sort(as.integer(free_cols))

      Q_diag <- if (!is.null(Q)) rep(diag(Q), H) else NULL

      res_inner <- .hard_anticipated(A_full, M, s0, cond_rows, b_cond,
                                     free_cols, n_draws, Q_diag, ss_raw, H)
    } else {
      ## Unanticipated: period-by-period
      res_inner <- .hard_unanticipated(ss_raw, s0, H, cond_df, free_shk_idx,
                                       n_draws, Q = Q)
    }
  }

  result <- structure(
    list(
      paths_point  = res_inner$paths_point,
      paths_uncond = paths_uncond,
      shock_paths  = res_inner$shock_paths,
      draws        = res_inner$draws,
      conditions   = conditions,
      horizon      = H,
      type         = type,
      method       = method,
      obs_names    = obs_vars,
      shock_names  = all_shocks
    ),
    class = "dynhr_cfcst"
  )
  result
}


## ---------------------------------------------------------------------------
#' Print method for conditional forecast results
#'
#' @param x   A \code{dynhr_cfcst} object.
#' @param ...  Unused.
#' @export
## ---------------------------------------------------------------------------
print.dynhr_cfcst <- function(x, ...) {
  cat(sprintf(
    "dynhr conditional forecast  [type=%s, method=%s, horizon=%d]\n",
    x$type, x$method, x$horizon))
  cat(sprintf("  Observables : %s\n", paste(x$obs_names, collapse = ", ")))
  cat(sprintf("  Shocks      : %s\n", paste(x$shock_names, collapse = ", ")))
  n_cond <- if (is.data.frame(x$conditions)) nrow(x$conditions) else 0L
  cat(sprintf("  Conditions  : %d\n", n_cond))
  if (!is.null(x$draws))
    cat(sprintf("  Draws       : %d\n", length(x$draws)))
  cat("\nPoint forecast paths (obs in deviation from steady state):\n")
  print(round(x$paths_point, 6))
  invisible(x)
}


## ---------------------------------------------------------------------------
#' Plot method for conditional forecast results
#'
#' Plots conditioned observable paths against the unconditional forecast,
#' optionally with quantile fan from draws.
#'
#' @param x      A \code{dynhr_cfcst} object.
#' @param vars   Character vector of observable names to plot. Defaults to all.
#' @param probs  Numeric vector of quantile probabilities for the fan (default
#'   \code{c(0.1, 0.9)}). Ignored if no draws.
#' @param ...    Additional arguments passed to \code{matplot()}.
#' @return \code{x}, invisibly.
#' @export
## ---------------------------------------------------------------------------
plot.dynhr_cfcst <- function(x, vars = NULL, probs = c(0.1, 0.9), ...) {
  if (is.null(vars)) vars <- x$obs_names
  bad <- setdiff(vars, x$obs_names)
  if (length(bad))
    stop(sprintf("plot.dynhr_cfcst: variable(s) not in forecast: %s.",
                 paste(bad, collapse = ", ")), call. = FALSE)
  idx <- match(vars, x$obs_names)
  n_v <- length(idx)
  H   <- x$horizon

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  n_cols <- min(n_v, 3L)
  n_rows <- ceiling(n_v / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1), mgp = c(1.5, 0.5, 0))

  for (v in vars) {
    vi   <- match(v, x$obs_names)
    pt   <- x$paths_point[, vi]
    unc  <- x$paths_uncond[, vi]
    ylim <- range(c(pt, unc), na.rm = TRUE)

    ## Quantile fan from draws
    if (!is.null(x$draws) && length(x$draws) > 0L) {
      draw_mat <- sapply(x$draws, function(d) d[, vi])
      q_lo <- apply(draw_mat, 1, quantile, probs = min(probs), na.rm = TRUE)
      q_hi <- apply(draw_mat, 1, quantile, probs = max(probs), na.rm = TRUE)
      ylim <- range(c(ylim, q_lo, q_hi), na.rm = TRUE)
    }

    plot(seq_len(H), pt, type = "n", ylim = ylim,
         xlab = "Horizon", ylab = "Dev. from SS", main = v,
         cex.main = 0.9, ...)

    if (!is.null(x$draws) && length(x$draws) > 0L) {
      polygon(c(seq_len(H), rev(seq_len(H))), c(q_lo, rev(q_hi)),
              col = "#2166AC30", border = NA)
    }

    lines(seq_len(H), unc,  col = "#888888", lty = 2, lwd = 1)
    lines(seq_len(H), pt,   col = "#2166AC", lwd = 1.5)

    ## Mark conditioned points
    if (is.data.frame(x$conditions) && nrow(x$conditions) > 0L) {
      cv <- x$conditions[x$conditions$var == v, , drop = FALSE]
      if (nrow(cv) > 0L)
        points(cv$horizon, cv$value, pch = 19, col = "#B2182B", cex = 0.9)
    }
  }
  invisible(x)
}
