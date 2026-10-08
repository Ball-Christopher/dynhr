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
  for (h in seq_len(H - 1L) + 1L) T_pow[[h]] <- T_pow[[h - 1L]] %*% T_mat

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
## Internal: law of the FREE shocks under the full shock covariance Q
##
## The conditioning solvers pick the most likely shock path under
## eps ~ N(0, Q). Only the free shocks move; the others are held at a given
## value eps_n (zero on the point path, a draw in the unanticipated draws).
## Given eps_n, the free shocks are Gaussian with
##   mean  B eps_n,          B = Q_fn Q_nn^+
##   cov   S = Q_ff - Q_fn Q_nn^+ Q_nf      (Schur complement)
## so the whitening is eps_f = B eps_n + L eta, eta ~ N(0, I), L L' = S.
## Whitening with sqrt(diag(Q)) instead (the old code) dropped every
## shock correlation -- the Q = I bug class one level down. With all shocks
## free this is exactly the GLS solution eps* = Q A'(A Q A')^+ b.
##
## @param Q        n_shock x n_shock covariance, or NULL (unit metric)
## @param free_idx Integer indices of the free shocks, in the order the
##   caller stacks them
## @param n_shk    Number of shocks
## @return list(L = |free| x |free| square root of S, B = |free| x |nonfree|,
##   L_n = square root of Q_nn, nonfree_idx)
## ---------------------------------------------------------------------------
.cf_sqrt_psd <- function(S) {
  S <- (S + t(S)) * 0.5
  n <- nrow(S)
  if (n == 0L) return(matrix(0, 0L, 0L))
  off <- S
  diag(off) <- 0
  if (all(off == 0)) return(diag(sqrt(pmax(diag(S), 0)), nrow = n))
  eg <- eigen(S, symmetric = TRUE)
  d  <- pmax(eg$values, 0)
  eg$vectors %*% (sqrt(d) * t(eg$vectors))
}

.cf_free_shock_law <- function(Q, free_idx, n_shk) {
  nonfree_idx <- setdiff(seq_len(n_shk), free_idx)
  if (is.null(Q)) Q <- diag(n_shk)
  Q <- as.matrix(Q)
  if (!identical(dim(Q), c(n_shk, n_shk)))
    stop(sprintf("conditional_forecast: Q must be %d x %d (one row/column per shock).",
                 n_shk, n_shk), call. = FALSE)
  Q <- (Q + t(Q)) * 0.5
  Q_ff <- Q[free_idx, free_idx, drop = FALSE]
  if (length(nonfree_idx) == 0L) {
    return(list(L = .cf_sqrt_psd(Q_ff),
                B = matrix(0, length(free_idx), 0L),
                L_n = matrix(0, 0L, 0L), nonfree_idx = nonfree_idx))
  }
  Q_fn <- Q[free_idx, nonfree_idx, drop = FALSE]
  Q_nn <- Q[nonfree_idx, nonfree_idx, drop = FALSE]
  if (all(Q_fn == 0)) {
    B <- matrix(0, length(free_idx), length(nonfree_idx))
    S <- Q_ff
  } else {
    B <- Q_fn %*% MASS::ginv(Q_nn)
    S <- Q_ff - B %*% t(Q_fn)
  }
  list(L = .cf_sqrt_psd(S), B = B, L_n = .cf_sqrt_psd(Q_nn),
       nonfree_idx = nonfree_idx)
}


## ---------------------------------------------------------------------------
## Internal: refuse a hard condition the solved shocks do not actually meet.
## The minimum-norm solve truncates tiny singular values, so a condition that
## is reachable only through a numerically negligible shock would otherwise be
## silently missed. `achieved` and `target` are the delivered and requested
## (mean-adjusted) condition values.
## ---------------------------------------------------------------------------
.cf_check_exact <- function(achieved, target, label, where = "") {
  miss <- max(abs(achieved - target))
  if (!is.finite(miss) || miss > 1e-8 * (1 + max(abs(target))))
    .dynhr_abort(sprintf(
      "%s%s: the hard conditions cannot be met exactly by the free shocks under the shock covariance (largest miss %.3g) -- the restriction is infeasible.",
      label, where, miss), class = "dynhr_error_cf_infeasible")
  invisible(TRUE)
}


## ---------------------------------------------------------------------------
## Internal: the c CONDITIONED rows of the stacked response matrix
##
## Row (obs i, horizon h) of A is [Z_i T^{h-1-j} R for j < h ; D_i at j = h],
## and the matching row of M is Z_i T^{h-1}. This builds ONLY those rows (and
## only the free-shock columns, period-major in ascending shock order), so
## memory is O(c * n_free * H) instead of O(n_obs * n_shock * H^2). Equal to
## .build_forecast_stack()$A[cond_rows, free_cols] and ($M %*% s0)[cond_rows].
##
## @return list(A_c = c x (n_free*H), m0 = c-vector Z_i T^{h-1} s0)
## ---------------------------------------------------------------------------
.cf_conditioned_rows <- function(ss, s0, H, cond_df, free_idx) {
  T_mat <- ss$T_mat; R_mat <- ss$R_mat; Z_mat <- ss$Z_mat; D_mat <- ss$D_mat
  n_f <- length(free_idx)
  n_c <- nrow(cond_df)
  A_c <- matrix(0, n_c, n_f * H)
  m0  <- numeric(n_c)
  for (i in unique(cond_df$var_idx)) {
    ks   <- which(cond_df$var_idx == i)
    hmax <- max(cond_df$horizon[ks])
    ## v_p = Z_i T^p; Gs[p + 1, ] = v_p R (response of y_{p+1+j} to eps_j)
    Gs   <- matrix(0, hmax, ss$n_shock)
    vs0  <- numeric(hmax)
    v    <- Z_mat[i, , drop = FALSE]
    for (p in seq_len(hmax)) {
      Gs[p, ] <- v %*% R_mat
      vs0[p]  <- sum(v * s0)
      if (p < hmax) v <- v %*% T_mat
    }
    for (k in ks) {
      h <- cond_df$horizon[k]
      m0[k] <- vs0[h]
      if (h > 1L)
        A_c[k, seq_len((h - 1L) * n_f)] <-
          as.numeric(t(Gs[(h - 1L):1L, free_idx, drop = FALSE]))
      A_c[k, (h - 1L) * n_f + seq_len(n_f)] <- D_mat[i, free_idx]
    }
  }
  list(A_c = A_c, m0 = m0)
}

## Observable paths from a shock matrix (n_shock x H) by the state recursion.
.cf_paths_from_shocks <- function(ss, s0, H, eps_mat) {
  paths  <- matrix(0, H, ss$n_obs)
  s_prev <- s0
  for (h in seq_len(H)) {
    paths[h, ] <- as.numeric(ss$Z_mat %*% s_prev) +
      as.numeric(ss$D_mat %*% eps_mat[, h])
    s_prev <- as.numeric(ss$T_mat %*% s_prev) +
      as.numeric(ss$R_mat %*% eps_mat[, h])
  }
  colnames(paths) <- ss$obs_names_fcst
  paths
}


## ---------------------------------------------------------------------------
## Internal: hard anticipated conditioning (Waggoner-Zha 1999)
##
## Minimum-variance (most likely) shock path under the free-shock covariance
## S_stack = I_H (x) L L' (L from .cf_free_shock_law(); shocks are iid over
## time so the square root is block-diagonal and is applied per period), in
## whitened coordinates eta (eps_free = W eta, W = I_H (x) L):
##   B = A_c W                      (c x F; only the c conditioned rows)
##   t(B) P = U R                   (thin rank-revealing QR, U F x c)
##   eta*   = U R^{-T} P' b
## Draws: eta_j = eta* + z - U U' z, z ~ N(0, I_F) -- the exact Gaussian
## conditional; no F x F matrix and no n_obs*H x n_shock*H stack is formed.
##
## @param Q  n_shock x n_shock shock covariance, or NULL (unit metric)
## @param free_idx sorted indices of the free shocks
## @param m_path optional n_shock x H matrix of deterministic shock MEANS
##   (IRIS `vary`): eps = m + u, u ~ N(0, Q). The conditions are solved for
##   u around the mean-inclusive baseline; shock_paths reports m + u.
## @return List: paths_point (H x n_obs), shock_paths (H x n_shock), draws
## ---------------------------------------------------------------------------
.hard_anticipated <- function(ss, s0, H, cond_df, free_idx, Q, n_draws,
                              m_path = NULL) {
  n_shk  <- ss$n_shock
  n_f    <- length(free_idx)
  n_free <- n_f * H
  n_cond <- nrow(cond_df)

  if (n_cond > n_free) {
    .dynhr_abort(sprintf(
      "conditional_forecast: %d conditions but only %d free shock-periods -- system is over-determined.",
      n_cond, n_free), class = "dynhr_error_cf_infeasible")
  }

  L <- if (!is.null(Q)) .cf_free_shock_law(Q, free_idx, n_shk)$L else diag(n_f)
  cr  <- .cf_conditioned_rows(ss, s0, H, cond_df, free_idx)
  A_c <- cr$A_c
  b_cond <- cond_df$value - cr$m0
  if (!is.null(m_path)) {
    ## Baseline path of the deterministic mean: y^m = M s0 + A m.
    base <- .cf_paths_from_shocks(ss, s0, H, m_path)
    b_cond <- cond_df$value -
      base[cbind(as.integer(cond_df$horizon), cond_df$var_idx)]
  }
  m_use <- if (is.null(m_path)) matrix(0, n_shk, H) else m_path

  ## B = A_c (I_H (x) L), block by period
  B <- A_c
  if (!is.null(Q)) {
    for (j in seq_len(H)) {
      blk <- (j - 1L) * n_f + seq_len(n_f)
      B[, blk] <- A_c[, blk, drop = FALSE] %*% L
    }
  }
  free_eps <- function(eta) L %*% matrix(eta, nrow = n_f, ncol = H)
  apply_W <- function(eta) {
    out <- matrix(0, n_shk, H)
    out[free_idx, ] <- free_eps(eta)
    out
  }

  eta_star <- numeric(n_free)
  U <- NULL
  if (n_cond > 0L) {
    qrB <- qr(t(B))
    ## Feasibility is a property of the WHITENED system: a shock with zero
    ## variance cannot move.
    if (qrB$rank < n_cond) {
      .dynhr_abort(sprintf(
        "conditional_forecast: the whitened condition system has rank %d < %d conditions -- the restriction is infeasible with the chosen free_shocks and shock covariance (a condition reachable only through zero-variance shocks cannot be met).",
        qrB$rank, n_cond), class = "dynhr_error_cf_infeasible")
    }
    U  <- qr.Q(qrB)[, seq_len(n_cond), drop = FALSE]
    Rm <- qr.R(qrB)[seq_len(n_cond), seq_len(n_cond), drop = FALSE]
    y  <- forwardsolve(t(Rm), b_cond[qrB$pivot[seq_len(n_cond)]])
    eta_star <- as.numeric(U %*% y)
  }
  eps_star <- apply_W(eta_star)
  if (n_cond > 0L)
    .cf_check_exact(as.numeric(A_c %*% as.numeric(free_eps(eta_star))),
                    b_cond, "conditional_forecast")

  eps_star <- eps_star + m_use
  paths_point <- .cf_paths_from_shocks(ss, s0, H, eps_star)
  shock_paths <- t(eps_star)
  colnames(shock_paths) <- ss$shock_names

  draws_out <- NULL
  if (n_draws > 0L) {
    draws_out <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      z <- rnorm(n_free)
      if (!is.null(U)) z <- z - as.numeric(U %*% crossprod(U, z))
      draws_out[[i]] <- .cf_paths_from_shocks(ss, s0, H,
                                              apply_W(eta_star + z) + m_use)
    }
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
## @param m_path       optional n_shock x H deterministic shock means (see
##   .hard_anticipated); the solved shock is the random part around m_h.
## @return Same structure as .hard_anticipated
## ---------------------------------------------------------------------------
.hard_unanticipated <- function(ss, s0, H, cond_df, free_shk_idx, n_draws,
                                Q = NULL, m_path = NULL) {
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

  ## Q-weighted solution, same whitening as .hard_anticipated: given the
  ## non-free shocks eps_n, eps_free = B eps_n + L eta with L L' the FULL
  ## conditional covariance of the free shocks (.cf_free_shock_law()), and the
  ## minimum-norm eta is the most likely shock under N(0, Q). Q = NULL keeps
  ## the unit metric.
  law         <- .cf_free_shock_law(Q, free_shk_idx, n_shk)
  L_free      <- law$L
  B_free      <- law$B
  L_non       <- law$L_n
  nonfree_idx <- law$nonfree_idx

  if (is.null(m_path)) m_path <- matrix(0, n_shk, H)

  s_prev <- s0
  for (h in seq_len(H)) {
    m_h <- m_path[, h]
    ## Unconditional predicted obs and state (the mean's own contribution
    ## D m_h is part of the baseline)
    y_pred_mean <- as.numeric(Z_mat %*% s_prev) +
      as.numeric(D_mat %*% m_h)

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

      ## Whitened feasibility (see .hard_anticipated): zero-variance free
      ## shocks cannot deliver a condition.
      D_w <- D_c %*% L_free
      rk <- qr(D_w)$rank
      if (rk < n_c)
        .dynhr_abort(sprintf(
          "conditional_forecast (unanticipated): at horizon %d, the whitened system has rank %d < %d conditions -- infeasible with the chosen free_shocks and shock covariance.",
          h, rk, n_c), class = "dynhr_error_cf_infeasible")
      DDt <- tcrossprod(D_w)
      DDt <- (DDt + t(DDt)) * 0.5
      ## DDt loses rank only when the conditions outnumber the free shocks
      ## (refused above); the solve is minimum-norm, not a ridge.
      DDt_inv <- .safe_inv(DDt, warn_label = "conditional_forecast: shock system")
      eps_free <- as.numeric(L_free %*% (t(D_w) %*% (DDt_inv %*% b_h)))
      .cf_check_exact(as.numeric(D_c %*% eps_free), b_h,
                      "conditional_forecast (unanticipated)",
                      sprintf(" at horizon %d", h))
      eps_h[free_shk_idx] <- eps_free
    }

    y_h <- y_pred_mean + as.numeric(D_mat %*% eps_h)
    paths_point[h, ] <- y_h
    shock_paths[h, ] <- eps_h + m_h

    ## Propagate state with the TOTAL shock m + u
    s_prev <- as.numeric(T_mat %*% s_prev) +
      as.numeric(R_mat %*% (eps_h + m_h))
  }

  ## Draws: same structure but with noise on unconstrained shock components
  draws_out <- NULL
  if (n_draws > 0L) {
    draws_list <- vector("list", n_draws)
    for (i in seq_len(n_draws)) {
      s_d   <- s0
      paths_d <- matrix(0, H, n_obs)
      for (h in seq_len(H)) {
        m_h <- m_path[, h]
        y_pred_mean_d <- as.numeric(Z_mat %*% s_d) +
          as.numeric(D_mat %*% m_h)
        h_cond <- cond_df[cond_df$horizon == h, , drop = FALSE]

        ## Draw all shocks from N(0, Q): the non-free block from its
        ## marginal, the free block from its law given the non-free draw,
        ## eps_f = B eps_n + L eta (correlations kept in both pieces).
        z     <- rnorm(n_shk)
        eps_d <- numeric(n_shk)
        if (length(nonfree_idx))
          eps_d[nonfree_idx] <- as.numeric(L_non %*% z[nonfree_idx])
        mu_free  <- as.numeric(B_free %*% eps_d[nonfree_idx])
        eta_free <- z[free_shk_idx]

        if (nrow(h_cond) > 0L) {
          cond_idx <- h_cond$var_idx
          ## The condition must hold for the TOTAL shock, so the non-free
          ## draw's and the conditional mean's loadings come off the target.
          b_h      <- h_cond$value - y_pred_mean_d[cond_idx] -
            as.numeric(D_mat[cond_idx, , drop = FALSE] %*% eps_d) -
            as.numeric(D_mat[cond_idx, free_shk_idx, drop = FALSE] %*% mu_free)

          ## Project out in whitened eta space: set free shocks to satisfy
          ## conditions, keep null-space component for the draw
          D_c <- D_mat[cond_idx, free_shk_idx, drop = FALSE]
          D_w <- D_c %*% L_free
          DDt <- tcrossprod(D_w)
          DDt_inv <- .safe_inv((DDt + t(DDt)) * 0.5,
                               warn_label = "conditional_forecast: shock system")
          ## Particular solution (whitened)
          eta_particular <- t(D_w) %*% (DDt_inv %*% b_h)
          ## Null-space component (keeps random variation on free shocks
          ## orthogonal to the conditions)
          Proj <- t(D_w) %*% DDt_inv %*% D_w   ## n_free x n_free
          eta_null <- (diag(length(free_shk_idx)) - Proj) %*% eta_free
          eta_free <- as.numeric(eta_particular + eta_null)
        }
        eps_d[free_shk_idx] <- mu_free + as.numeric(L_free %*% eta_free)

        y_d <- y_pred_mean_d + as.numeric(D_mat %*% eps_d)
        paths_d[h, ] <- y_d
        s_d <- as.numeric(T_mat %*% s_d) +
          as.numeric(R_mat %*% (eps_d + m_h))
      }
      colnames(paths_d) <- ss$obs_names_fcst
      draws_list[[i]] <- paths_d
    }
    draws_out <- draws_list
  }

  list(paths_point = paths_point, shock_paths = shock_paths, draws = draws_out)
}


## ---------------------------------------------------------------------------
## Internal: soft conditioning (exact Gaussian conditional laws)
##
## Each condition is a noisy pseudo-observation of a forecast-period
## observable,
##   value_k = y_{h_k, i_k} + u_k,   u_k ~ N(0, stderr_k^2)   (0 for hard rows).
## s0 is taken as known, exactly as in the hard solvers. Only the FREE shocks
## move; the others are fixed at zero, and the free shocks have the law
## N(0, S) with S the Schur complement given the non-free shocks at zero
## (.cf_free_shock_law()), eps_free = L eta, eta ~ N(0, I).
## Paths are always generated THROUGH the model transition from s0 with the
## shocks in hand, so the contemporaneous D eps_h term is in every path and
## every draw is a model-consistent trajectory
## (y_h = Z s_{h-1} + D eps_h along its own state path).
##
## type = "anticipated": all conditions are known at the origin. The stacked
##   shock vector has Gaussian posterior
##     eta | conds ~ N(A_w' G^+ b,  I - A_w' G^+ A_w),
##     A_w = A_c W,  G = A_w A_w' + diag(stderr^2),  b = value - (M s0)_c,
##   the point path is the posterior mean and each draw is a joint draw of the
##   whole sequence. stderr -> 0 gives the hard anticipated (Waggoner-Zha)
##   solution.
##
## type = "unanticipated": period by period, agents do not foresee later
##   conditions. eps_h is drawn from its exact conditional law given ONLY the
##   conditions at horizon h and the realised path up to h-1:
##     eta_h | conds_h ~ N(D_w' G_h^+ b_h,  I - D_w' G_h^+ D_w),
##     D_w = D_c L,  G_h = D_w D_w' + diag(stderr^2),
##     b_h = value - Z_c s_{h-1}.
##   The state then advances with the realised shock. The point path is the
##   same recursion with the conditional means; draws are sequential exact
##   conditional draws. stderr -> 0 gives .hard_unanticipated().
##
## @param free_shk_idx Indices of the shocks allowed to move
## @param type "anticipated" or "unanticipated"
## @param Q  n_shock x n_shock shock covariance (NULL = identity; correct only
##   when every shock has stderr 1, since ghu holds unit-shock responses)
## @param m_path optional n_shock x H deterministic shock means; the
##   conditioning acts on the random part u around the mean-inclusive baseline
##   and the reported/propagated shock is m + u.
## ---------------------------------------------------------------------------
.soft_forecast <- function(ss, s0, H, cond_df, n_draws, Q = NULL,
                           free_shk_idx = seq_len(ss$n_shock),
                           type = "anticipated", m_path = NULL) {
  n_obs <- ss$n_obs
  n_shk <- ss$n_shock
  T_mat <- ss$T_mat
  R_mat <- ss$R_mat
  Z_mat <- ss$Z_mat
  D_mat <- ss$D_mat

  free_shk_idx <- sort(free_shk_idx)
  n_f   <- length(free_shk_idx)
  L     <- .cf_free_shock_law(Q, free_shk_idx, n_shk)$L
  se    <- cond_df$stderr
  se[is.na(se)] <- 0
  m_use <- if (is.null(m_path)) matrix(0, n_shk, H) else m_path

  ## Advance the model through the transition with a shock matrix
  ## (n_shock x H), returning the observable paths.
  path_from_shocks <- function(eps_mat) {
    paths  <- matrix(0, H, n_obs)
    s_prev <- s0
    for (h in seq_len(H)) {
      paths[h, ] <- as.numeric(Z_mat %*% s_prev) +
        as.numeric(D_mat %*% eps_mat[, h])
      s_prev <- as.numeric(T_mat %*% s_prev) +
        as.numeric(R_mat %*% eps_mat[, h])
    }
    colnames(paths) <- ss$obs_names_fcst
    paths
  }

  eps_from_eta <- function(eta_vec, n_periods_eta) {
    ## eta_vec is period-major over the free shocks
    eta_mat <- matrix(eta_vec, nrow = n_f, ncol = n_periods_eta)
    out <- matrix(0, n_shk, n_periods_eta)
    out[free_shk_idx, ] <- L %*% eta_mat
    out
  }

  if (type == "anticipated") {
    stk <- .build_forecast_stack(ss, H)
    cond_rows <- (cond_df$horizon - 1L) * n_obs + cond_df$var_idx
    b_cond <- cond_df$value - as.numeric(stk$M %*% s0)[cond_rows]
    if (!is.null(m_path))
      b_cond <- b_cond - as.numeric(stk$A %*% as.numeric(m_path))[cond_rows]
    free_cols <- sort(as.integer(outer(free_shk_idx - 1L,
                                       (0L:(H - 1L)) * n_shk, "+") + 1L))
    A_w <- stk$A[cond_rows, free_cols, drop = FALSE] %*% kronecker(diag(H), L)
    G   <- tcrossprod(A_w) + diag(se^2, nrow = length(se))
    G   <- (G + t(G)) * 0.5
    G_inv <- if (length(se)) {
      .safe_inv(G, warn_label = "conditional_forecast: soft condition system")
    } else G
    eta_mean <- as.numeric(crossprod(A_w, G_inv %*% b_cond))
    n_eta    <- length(eta_mean)
    C_sqrt <- NULL
    if (n_draws > 0L) {
      Sigma_c <- diag(n_eta) - crossprod(A_w, G_inv %*% A_w)
      C_sqrt  <- .cf_sqrt_psd(Sigma_c)
    }
    eps_pt <- eps_from_eta(eta_mean, H) + m_use
    draw_eps <- function() eps_from_eta(eta_mean + as.numeric(C_sqrt %*% rnorm(n_eta)), H) + m_use
  } else {
    ## Sequential recursion; `draw = FALSE` gives the conditional-mean path.
    run_seq <- function(draw) {
      eps_mat <- matrix(0, n_shk, H)
      s_prev  <- s0
      for (h in seq_len(H)) {
        k   <- which(cond_df$horizon == h)
        eta <- if (draw) rnorm(n_f) else numeric(n_f)
        if (length(k)) {
          idx <- cond_df$var_idx[k]
          b_h <- cond_df$value[k] -
            as.numeric(Z_mat[idx, , drop = FALSE] %*% s_prev) -
            as.numeric(D_mat[idx, , drop = FALSE] %*% m_use[, h])
          D_w <- D_mat[idx, free_shk_idx, drop = FALSE] %*% L
          G   <- tcrossprod(D_w) + diag(se[k]^2, nrow = length(k))
          G   <- (G + t(G)) * 0.5
          G_inv <- .safe_inv(G, warn_label = "conditional_forecast: soft condition system")
          eta_mean <- as.numeric(crossprod(D_w, G_inv %*% b_h))
          eta <- if (draw) {
            Sigma_c <- diag(n_f) - crossprod(D_w, G_inv %*% D_w)
            eta_mean + as.numeric(.cf_sqrt_psd(Sigma_c) %*% eta)
          } else eta_mean
        }
        eps_mat[free_shk_idx, h] <- as.numeric(L %*% eta)
        s_prev <- as.numeric(T_mat %*% s_prev) +
          as.numeric(R_mat %*% (eps_mat[, h] + m_use[, h]))
      }
      eps_mat + m_use
    }
    eps_pt   <- run_seq(FALSE)
    draw_eps <- function() run_seq(TRUE)
  }

  paths_point <- path_from_shocks(eps_pt)
  shock_paths <- t(eps_pt)
  colnames(shock_paths) <- ss$shock_names

  draws_out <- NULL
  if (n_draws > 0L) {
    draws_out <- vector("list", n_draws)
    for (i in seq_len(n_draws)) draws_out[[i]] <- path_from_shocks(draw_eps())
  }

  list(paths_point = paths_point, shock_paths = shock_paths, draws = draws_out)
}


## ---------------------------------------------------------------------------
## Internal: unconditional forecast (plain propagation from s0)
## Used for the no-conditions case and internal consistency checks.
## ---------------------------------------------------------------------------
.unconditional_forecast <- function(ss, s0, H, m_path = NULL) {
  T_mat <- ss$T_mat
  Z_mat <- ss$Z_mat
  n_obs <- ss$n_obs

  paths <- matrix(0, H, n_obs)
  colnames(paths) <- ss$obs_names_fcst
  s_prev <- s0
  for (h in seq_len(H)) {
    paths[h, ] <- as.numeric(Z_mat %*% s_prev)
    s_prev <- as.numeric(T_mat %*% s_prev)
    if (!is.null(m_path)) {
      ## Deterministic contribution of the shock means (IRIS `vary`).
      paths[h, ] <- paths[h, ] + as.numeric(ss$D_mat %*% m_path[, h])
      s_prev <- s_prev + as.numeric(ss$R_mat %*% m_path[, h])
    }
  }
  paths
}

## ---------------------------------------------------------------------------
## Internal: resolve conditional_forecast()'s shock_means into the sample part
## and the horizon part, both DATED (timing already applied).
##
## `shock_means` is an n_exo x (T + H) matrix (rows by name when named; NA/0 =
## no shift), or a named numeric vector = a constant mean for every period
## (unnamed shocks get 0). `shock_timing` is read over the WHOLE T + H axis,
## exactly as kalman_filter() reads it, so that under "transition_next" the
## sample's last column lands on forecast period 1. n_T = 0 is the
## smoother-result case: the sample part is already inside the supplied
## filtered state and the matrix covers the horizon only.
## @return list(sample = n_exo x T matrix or NULL, horizon = n_exo x H matrix
##   or NULL). NULL horizon = no deterministic forecast contribution.
## ---------------------------------------------------------------------------
.cf_resolve_shock_means <- function(shock_means, shock_timing, shock_names,
                                    n_T, H) {
  if (is.null(shock_means)) return(list(sample = NULL, horizon = NULL))
  n_e <- length(shock_names)
  if (is.null(dim(shock_means))) {
    v <- shock_means
    if (!is.numeric(v) || is.null(names(v)) || any(!nzchar(names(v))))
      .dynhr_abort("conditional_forecast: a vector `shock_means` must be a NAMED numeric vector (shock name -> constant mean); use an n_exo x (T + horizon) matrix for a time-varying path.",
                   class = "dynhr_error_bad_argument")
    bad <- setdiff(names(v), shock_names)
    if (length(bad))
      .dynhr_abort(sprintf("conditional_forecast: `shock_means` names not in the model's shocks: %s.",
                           paste(bad, collapse = ", ")),
                   class = "dynhr_error_bad_argument")
    full <- numeric(n_e); names(full) <- shock_names
    full[names(v)] <- v
    shock_means <- matrix(full, n_e, n_T + H, dimnames = list(shock_names, NULL))
  }
  M <- .kf_shock_means(shock_means, shock_timing, shock_names, n_T + H,
                       what = "conditional_forecast")
  if (is.null(M)) return(list(sample = NULL, horizon = NULL))
  dimnames(M) <- list(shock_names, NULL)
  list(sample  = if (n_T > 0L && any(M[, seq_len(n_T)] != 0))
                   M[, seq_len(n_T), drop = FALSE] else NULL,
       horizon = if (any(M[, n_T + seq_len(H)] != 0))
                   M[, n_T + seq_len(H), drop = FALSE] else NULL)
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
## @param shock_means_sample,a0,P0  Sample part of conditional_forecast()'s
##   shock_means (already dated, n_exo x T) and the filter's initial state /
##   covariance. Only the Gaussian/Whittle smoother can carry them; every other
##   dispatch branch refuses (classed) rather than silently ignore them.
## @return list(s0 = n_state numeric, P0 = n_state x n_state matrix).
## ---------------------------------------------------------------------------
.extract_terminal_state <- function(data, ss_raw, Q, ctx, model, dr, obs_vars,
                                    compiled = NULL, shock_means_sample = NULL,
                                    a0 = NULL, P0 = NULL,
                                    shock_means_given = FALSE,
                                    lik_init = "auto") {

  lik <- if (is.null(ctx)) "gaussian" else ctx$likelihood
  ## Only the Gaussian Kalman smoother subtracts a deterministic trend; the
  ## particle / OBC terminal-state branches would filter untrended data.
  if (!lik %in% c("gaussian", "whittle"))
    .refuse_obs_trends(model, sprintf("conditional_forecast(ctx$likelihood = \"%s\")", lik),
                       dr = dr)
  if ((shock_means_given || !is.null(shock_means_sample) || !is.null(a0) ||
       !is.null(P0)) && !lik %in% c("gaussian", "whittle"))
    .dynhr_abort(sprintf(
      "conditional_forecast: `shock_means` / `a0` / `P0` are carried only by the Gaussian Kalman terminal-state filter; ctx$likelihood = \"%s\" cannot take them (silently dropping them would forecast from the wrong state).",
      lik), class = "dynhr_error_bad_argument")
  if (!identical(lik_init, "auto") && !lik %in% c("gaussian", "whittle"))
    .dynhr_abort(sprintf(
      "conditional_forecast: `lik_init` initialises the Gaussian Kalman terminal-state filter only; ctx$likelihood = \"%s\" has its own terminal-state run. Leave lik_init = \"auto\".",
      lik), class = "dynhr_error_bad_argument")

  ## Dispatch table — "pkf", "ppf", "copf" added for OBC models (the
  ## ppf/copf keys exist so particle-filter-estimated OBC models use
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
    sm <- .kalman_smoother_ss(data, ss_raw, Q = Q,
                              shock_means = shock_means_sample,
                              shock_timing = "dated", a0 = a0, P0 = P0,
                              lik_init = lik_init)
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
    }, error = function(e) .dynhr_reraise_bug(e, QQ_s))

    ## Seeded for the particle cloud only; the caller's RNG stream is restored
    ## when this function returns.
    .local_seed(tpf_seed)

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

    obs_idx   <- match(obs_vars, dr$endo_names)
    if (any(is.na(obs_idx)))
      stop("conditional_forecast (pskf): some obs_names not found in dr$endo_names.",
           call. = FALSE)

    ## CSN shock parameters AND the contemporaneous state space they belong
    ## to (csn$TT / csn$ZZ; see .pskf_order1_statespace in
    ## R/pskf-likelihood.R). The lagged-convention pair (ghx[state, ],
    ## ghx[obs, ]) with eps = DD e is not the model's law when DD != 0.
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
        TT        = csn$TT,
        ZZ        = csn$ZZ,
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

    ## PSKF timing: the contemporaneous state chi_t CONTAINS s_t (its first
    ## n_s coordinates, csn$state_pos), and the forecaster consumes s_T in the
    ## lagged convention (y_{T+1} = Z s_T + D e_{T+1}), so the initial
    ## condition is the time-T (smoothed == filtered) state block itself --
    ## no extra transition. (The former s0 = TT x_T, P0 = TT P_T TT' +
    ## Sigma_eta was the matching step for the old x_t := s_{t-1} reading.)
    T_obs <- nrow(sm$smoothed_means)
    sp    <- csn$state_pos
    s0    <- as.numeric(sm$smoothed_means[T_obs, sp])
    P0    <- sm$smoothed_covs[sp, sp, T_obs]
    P0    <- matrix(P0, length(sp), length(sp))
    return(list(s0 = s0, P0 = P0))
  }

  ## ---- OBC PKF path -------------------------------------------------------
  if (identical(path, "pkf")) {
    ## Need compiled to call cache_system_structure -> extract_system_matrices_fast
    ## -> .solve_from_system -> obc_ensure_policy.  Either passed via compiled= or
    ## available from model (re-compile as fallback; ~50 ms overhead).
    if (is.null(compiled)) {
      .dynhr_inform(
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
    ## (no-op for non-SSM-parameter models).
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

  ## ---- OBC PPF / COPF path ----------------------------------
  ## Both "ppf" and "copf" summarise the terminal particle cloud to (mean, cov)
  ## exactly as the TPF path does, giving a Gaussian summary for the downstream
  ## Waggoner-Zha / soft-KF forecaster (which is purely Gaussian-state).
  ##
  ## The difference between "ppf" and "copf" is only the proposal argument
  ## forwarded to ppf_likelihood(); the terminal-state summary is identical.
  if (path %in% c("ppf", "copf")) {
    if (is.null(compiled)) {
      .dynhr_inform(
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
#'     \item{\code{value}}{Finite numeric: target value (in model deviation
#'       units). \code{NA} is refused, not dropped.}
#'     \item{\code{stderr}}{Numeric (optional): standard deviation of the
#'       conditioning noise for soft conditions. \code{NA}, absent or \code{0}
#'       means a hard condition; a negative or non-finite value is an error.}
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
#'   pseudo-observations with measurement error \code{stderr^2}; the point
#'   path uses the conditional-mean shocks and each draw's shocks are pushed
#'   through the model transition from the terminal state (so every draw is a
#'   model-consistent trajectory). With \code{type =
#'   "anticipated"} the whole shock sequence is drawn jointly given all soft
#'   conditions; with \code{type = "unanticipated"} each period's shock is
#'   drawn from its exact conditional law given only that period's conditions
#'   and the realised path so far. As \code{stderr} tends to 0 each variant
#'   reproduces the corresponding hard solution. \code{free_shocks} is honoured:
#'   the other shocks are held at zero.
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
#' @param shock_means Deterministic shock MEANS (IRIS \code{vary}): every
#'   shock is \eqn{\varepsilon = m + u}, \eqn{u \sim N(0, Q)}, \eqn{m}
#'   deterministic. Either an \code{n_exo x (nrow(data) + horizon)} matrix
#'   (rows by name when named; \code{NA} or 0 = no shift, as in
#'   \code{\link{kalman_filter}}), or a NAMED numeric vector giving a constant
#'   mean for every period (unnamed shocks get 0). This is how a
#'   trend-with-drift model (no static steady state) is forecast: solve the
#'   decision rule with the drift switched off, carry the drift as the mean of
#'   its shock, and start the growth state on its balanced-growth path with
#'   \code{a0}. The sample columns enter the terminal-state filter; the horizon
#'   columns enter the forecast baseline through the transition. Conditions
#'   (hard or soft, anticipated or not, with draws) are solved for the random
#'   part \eqn{u} around that mean-inclusive baseline; \code{shock_paths}
#'   reports the TOTAL \eqn{m + u}, and draws add \eqn{u}-noise only. When
#'   \code{data} is a \code{kalman_smoother()} result the sample part is
#'   already inside it, so \code{shock_means} then covers the horizon only
#'   (\code{horizon} columns). Refused (classed error) under a non-Gaussian
#'   \code{ctx$likelihood} (tpf, pskf, pkf, ppf, copf).
#' @param shock_timing How to read the columns of \code{shock_means}: as in
#'   \code{\link{kalman_filter}}, \code{"dated"} (column t enters period t) or
#'   \code{"transition_next"} (column t drives the transition out of t, so it
#'   lands on t + 1; the last sample column is the first forecast period).
#' @param a0,P0 Initial state mean / covariance of the terminal-state filter,
#'   as in \code{\link{kalman_smoother}} (e.g. the balanced-growth start of a
#'   drift model). Not allowed with a \code{kalman_smoother()} result as
#'   \code{data} or a non-Gaussian \code{ctx$likelihood}.
#' @param units \code{"auto"} (default), \code{"deviations"} or
#'   \code{"levels"}. \code{"levels"} (any model): conditions are read, and
#'   \code{paths_point}, \code{paths_uncond} and \code{draws} returned, as
#'   steady state + deterministic trend + deviation. \code{"deviations"}:
#'   deviations from steady state; refused on a trend model. \code{"auto"}:
#'   \code{"levels"} on a trend model, \code{"deviations"} otherwise. The
#'   resolved choice is \code{result$units}.
#' @param lik_init Initialisation of the Gaussian terminal-state filter, as in
#'   \code{\link{kalman_smoother}}: \code{"auto"} (default), \code{"stationary"},
#'   \code{"kappa"}, \code{"diffuse"} or \code{"fixed_unknown"} (IRIS's
#'   default: the unit-root level estimated by GLS). The terminal state
#'   \eqn{s_{T|T}} under \code{"fixed_unknown"} equals the exact diffuse one
#'   once the unit-root level is identified. Refused (classed error) with a
#'   \code{kalman_smoother()} result as \code{data} or a non-Gaussian
#'   \code{ctx$likelihood}.
#' @param ...         Currently unused.
#'
#' @section Trend models:
#' A model with \code{observation_trends} and/or a balanced-growth decision
#' rule (\code{solve_model(steady_options = list(growth = TRUE))}, pass the
#' returned \code{model} and \code{dr}) is supported with the Gaussian
#' terminal-state filter. The terminal state comes from the trend-honouring
#' smoother, the forecast is computed in deviations, and the deterministic
#' trend is added back at the sample's period index continued past the last
#' observation (\code{first_obs + T + h - 1}). On such models
#' \code{paths_point}, \code{paths_uncond} and \code{draws} are LEVELS
#' (steady-state intercept + slope * index + deviation), and \code{conditions}
#' are given in levels too. A \code{kalman_smoother()} result as \code{data}
#' works (T is its number of rows; a growth model needs its
#' \code{growth_path}). Models without trends keep the deviation convention.
#' Non-Gaussian terminal-state likelihoods refuse (classed
#' \code{dynhr_error_observation_trends_unsupported}). Composes with
#' \code{shock_means}.
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
#'     \item{\code{type}, \code{method}}{As used.}
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
                                 shock_means = NULL,
                                 shock_timing = c("dated", "transition_next"),
                                 a0 = NULL, P0 = NULL,
                                 units = c("auto", "deviations", "levels"),
                                 lik_init = c("auto", "stationary", "kappa",
                                              "diffuse", "fixed_unknown"),
                                 ...) {
  ## Own the message epoch for this run: repeat-suppressed warnings
  ## (`.dynhr_warn(once = TRUE)`) are keyed within it and re-arm for the
  ## next run, and the close reports what it suppressed. A nested call
  ## inherits this epoch rather than opening a second one.
  .dynhr_run_epoch <- .dynhr_epoch("conditional_forecast")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)
  ## Deterministic trends (observation_trends and/or a balanced-growth dr):
  ## the Gaussian terminal-state filter subtracts them; the forecast is
  ## computed in deviations and the trend is added back below. The model and
  ## the dr must agree about balanced growth.
  .check_growth_pair(model, dr)
  lik_init <- match.arg(lik_init)
  ## ---- Resolve plan= if supplied ----
  if (!is.null(plan)) {
    if (!inherits(plan, "dynhr_plan"))
      stop("conditional_forecast: 'plan' must be a dynhr_plan object.", call. = FALSE)
    if (!is.null(conditions))
      stop("conditional_forecast: supply either 'plan' or 'conditions', not both.",
           call. = FALSE)
    ## Warn if plan has shock_scale entries (out-of-sample scaling unsupported)
    if (length(plan$shock_scales) > 0L)
      .dynhr_warn(
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
  shock_timing <- match.arg(shock_timing)
  units <- match.arg(units)
  is_trend_model <- isTRUE(model$balanced_growth) || .has_obs_trends(model, dr)
  if (identical(units, "deviations") && is_trend_model)
    .dynhr_abort("conditional_forecast: units = \"deviations\" is meaningless on a trend model (observation_trends or balanced growth): deviations from an arbitrary unit-root level carry no information. Use units = \"levels\" (or \"auto\").",
                 class = "dynhr_error_bad_argument")
  if (identical(units, "auto"))
    units <- if (is_trend_model) "levels" else "deviations"
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

  ## ---- Deterministic shock means (IRIS `vary`) ----
  ## eps = m + u. The sample part of m (and a0 / P0) goes into the terminal
  ## filter; the horizon part shifts the forecast baseline.
  n_T_sample <- if (is.null(smoother_result)) nrow(data) else 0L
  if (!is.null(smoother_result) && (!is.null(a0) || !is.null(P0)))
    .dynhr_abort("conditional_forecast: `a0` / `P0` initialise the terminal-state filter, which was already run when `data` is a kalman_smoother() result; pass them to kalman_smoother() instead.",
                 class = "dynhr_error_bad_argument")
  if (!is.null(smoother_result) && !identical(lik_init, "auto"))
    .dynhr_abort("conditional_forecast: `lik_init` initialises the terminal-state filter, which was already run when `data` is a kalman_smoother() result; pass it to kalman_smoother() instead.",
                 class = "dynhr_error_bad_argument")
  sm_split <- .cf_resolve_shock_means(shock_means, shock_timing,
                                      ss_raw$shock_names, n_T_sample, H)
  m_path <- sm_split$horizon

  ## ---- Extract terminal state s_{T|T} and P_{T|T} ----
  if (!is.null(smoother_result)) {
    if (isTRUE(model$balanced_growth) && is.null(smoother_result$growth_path))
      .dynhr_abort("conditional_forecast: `model` is a balanced-growth model but the kalman_smoother() result in `data` carries no `growth_path`: it was run without the balanced-growth dr, so its terminal state is not a deviation from the growth path. Re-run kalman_smoother() with the dr and model from solve_model().",
                   class = "dynhr_error_balanced_growth_mismatch")
    s0 <- as.numeric(smoother_result$filtered_states[nrow(smoother_result$filtered_states), ])
    P0 <- if (!is.null(smoother_result$P_filt_last)) smoother_result$P_filt_last
          else diag(ss_raw$n_state) * 1e-6   ## fallback: near-zero uncertainty
  } else {
    ## Dispatch terminal-state extraction on ctx$likelihood
    ts <- .extract_terminal_state(data, ss_raw, Q, ctx, model, dr, obs_vars,
                                  compiled = compiled,
                                  shock_means_sample = sm_split$sample,
                                  a0 = a0, P0 = P0,
                                  shock_means_given = !is.null(shock_means),
                                  lik_init = lik_init)
    s0 <- ts$s0
    P0 <- ts$P0
  }

  ## ---- Deterministic trend over the horizon ----
  ## Period T + h of the sample's trend index (first_obs + t - 1, continued):
  ## the observable's level is d + slope * index + deviation. NULL when the
  ## model carries no trend (deviations in, deviations out, as before).
  trend_off <- NULL
  if (identical(units, "levels")) {
    n_T_idx <- if (is.null(smoother_result)) nrow(data)
               else nrow(smoother_result$filtered_states)
    idx <- ss_raw$obs_trend_first_obs - 1 + n_T_idx + seq_len(H)
    slopes <- if (is.null(ss_raw$obs_trend)) numeric(length(obs_vars))
              else as.numeric(ss_raw$obs_trend)
    trend_off <- outer(idx, slopes) +
      matrix(as.numeric(ss_raw$d), H, length(obs_vars), byrow = TRUE)
    colnames(trend_off) <- obs_vars
  }
  add_trend <- function(P) if (is.null(trend_off)) P else P + trend_off

  ## ---- Unconditional forecast ----
  paths_uncond <- add_trend(.unconditional_forecast(ss_raw, s0, H, m_path = m_path))

  ## ---- Handle null / empty conditions -> unconditional ----
  no_conds <- is.null(conditions) ||
    (is.data.frame(conditions) && nrow(conditions) == 0L)

  if (no_conds) {
    result <- structure(
      list(
        paths_point  = paths_uncond,
        paths_uncond = paths_uncond,
        shock_paths  = if (is.null(m_path))
                         matrix(0, H, ss_raw$n_shock,
                                dimnames = list(NULL, ss_raw$shock_names))
                       else `colnames<-`(t(m_path), ss_raw$shock_names),
        draws        = NULL,
        conditions   = data.frame(var = character(0), horizon = integer(0),
                                  value = numeric(0), stringsAsFactors = FALSE),
        horizon      = H,
        type         = type,
        method       = method,
        obs_names    = obs_vars,
        shock_names  = ss_raw$shock_names,
        units        = units
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

  cond_horizon <- conditions$horizon
  if (!is.numeric(cond_horizon) || any(!is.finite(cond_horizon)) ||
      any(cond_horizon != round(cond_horizon)))
    .dynhr_abort("conditional_forecast: condition horizons must be finite integer-valued numbers (no NA, no fractions).",
                 class = "dynhr_error_bad_argument")
  bad_h <- which(cond_horizon < 1L | cond_horizon > H)
  if (length(bad_h))
    .dynhr_abort(sprintf("conditional_forecast: condition horizon(s) out of range [1, %d]: %s.",
                 H, paste(cond_horizon[bad_h], collapse = ", ")),
                 class = "dynhr_error_bad_argument")

  if (!is.numeric(conditions$value) || any(!is.finite(conditions$value)))
    .dynhr_abort("conditional_forecast: condition values must be finite numbers (an NA value would silently drop the condition; remove the row instead).",
                 class = "dynhr_error_bad_argument")
  if ("stderr" %in% names(conditions)) {
    se_in <- conditions$stderr
    if (!is.numeric(se_in) && !all(is.na(se_in)))
      .dynhr_abort("conditional_forecast: condition stderr must be numeric (NA or 0 = hard, finite >= 0 = soft noise).",
                   class = "dynhr_error_bad_argument")
    if (any(!is.na(se_in) & (!is.finite(se_in) | se_in < 0)))
      .dynhr_abort("conditional_forecast: condition stderr must be NA (hard), 0 (hard) or a finite non-negative number.",
                   class = "dynhr_error_bad_argument")
  }

  ## Add var_idx and normalise stderr
  cond_df <- conditions
  cond_df$var_idx <- match(cond_df$var, obs_vars)
  if (!"stderr" %in% names(cond_df)) cond_df$stderr <- NA_real_
  cond_df$stderr <- as.numeric(cond_df$stderr)
  ## Conditions are in LEVELS; the solvers work in deviations.
  if (!is.null(trend_off))
    cond_df$value <- cond_df$value -
      trend_off[cbind(as.integer(cond_df$horizon), cond_df$var_idx)]

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
    ## Soft: exact Gaussian conditional laws (joint posterior of the stacked
    ## shock path for "anticipated", period-by-period for "unanticipated").
    res_inner <- .soft_forecast(ss_raw, s0, H, cond_df, n_draws, Q = Q,
                                free_shk_idx = free_shk_idx, type = type,
                                m_path = m_path)

  } else {
    ## Hard conditioning
    if (type == "anticipated") {
      res_inner <- .hard_anticipated(ss_raw, s0, H, cond_df,
                                     sort(free_shk_idx), Q, n_draws,
                                     m_path = m_path)
    } else {
      ## Unanticipated: period-by-period
      res_inner <- .hard_unanticipated(ss_raw, s0, H, cond_df, free_shk_idx,
                                       n_draws, Q = Q, m_path = m_path)
    }
  }

  result <- structure(
    list(
      paths_point  = add_trend(res_inner$paths_point),
      paths_uncond = paths_uncond,
      shock_paths  = res_inner$shock_paths,
      draws        = if (is.null(res_inner$draws)) NULL
                     else lapply(res_inner$draws, add_trend),
      conditions   = conditions,
      horizon      = H,
      type         = type,
      method       = method,
      obs_names    = obs_vars,
      shock_names  = all_shocks,
      units        = units
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
  u <- x$units %||% "deviations"
  cat(sprintf("  Units       : %s\n", u))
  cat(if (identical(u, "levels"))
        "\nPoint forecast paths (levels: steady state + trend + deviation):\n"
      else "\nPoint forecast paths (obs in deviation from steady state):\n")
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
         xlab = "Horizon", ylab = if (identical(x$units, "levels")) "Level" else "Dev. from SS", main = v,
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
