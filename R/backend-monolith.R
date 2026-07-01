###############################################################################
#  dynhr_backend.R
#  General-purpose state-space toolkit for dynhr.
#  Dependencies: base R + Matrix (for Schur decomposition in build_state_space).
###############################################################################

# =============================================================================
# 0.  HELPERS
# =============================================================================

## `%||%` is now defined once in R/utils-pipe.R (loaded first via Collate).
## Local definition removed during phase-0 packaging.

.sym <- function(X) (X + t(X)) / 2

.solve_lyapunov <- function(A, Q, tol = .LYAP_TOL) {
  n <- nrow(A)
  if (n == 1L) return(matrix(Q / (1 - A[1, 1]^2), 1, 1))
  M <- diag(n * n) - kronecker(A, A)
  # Check conditioning before solve; return NaN on singular system
  if (rcond(M) < .Machine$double.eps) return(matrix(NaN, n, n))
  P_vec <- solve(M, as.vector(Q))
  matrix(P_vec, n, n)
}

.pinv <- function(M, tol = .PINV_TOL) {
  s <- svd(M)
  pos <- s$d > max(tol * s$d[1], .Machine$double.eps)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% diag(1 / s$d[pos], nrow = sum(pos)) %*%
    t(s$u[, pos, drop = FALSE])
}

#' Robust Cholesky decomposition with diagonal fallback
#'
#' Computes t(chol(Sigma)) for a covariance matrix Sigma. If Sigma is not
#' positive definite, falls back to diag(sqrt(diag(Sigma))).
#' @param Sigma Covariance matrix
#' @param n     Dimension
#' @return Lower triangular factor (n x n)
#' @noRd
.robust_chol <- function(Sigma, n) {
  eig_vals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig_vals) > .Machine$double.eps * max(eig_vals) * n) {
    t(chol(Sigma))
  } else {
    diag(sqrt(pmax(diag(Sigma), 0)), nrow = n)
  }
}

.smoother_gain <- function(P_filt, F_mat, P_pred) {
  P_filt %*% t(F_mat) %*% .pinv(P_pred)
}

.swap_schur_11 <- function(T, Q, k) {
  n <- nrow(T)
  a <- T[k, k];  b <- T[k + 1, k + 1]
  if (abs(a - b) < 1e-15) {
    T[k, k] <- b;  T[k + 1, k + 1] <- a
    return(list(T = T, Q = Q))
  }
  t12 <- T[k, k + 1]
  r   <- sqrt(t12^2 + (a - b)^2)
  cc  <- t12 / r;  ss <- (a - b) / r
  G <- diag(n)
  G[k, k] <- cc;      G[k, k + 1] <- ss
  G[k + 1, k] <- -ss;  G[k + 1, k + 1] <- cc
  T <- crossprod(G, T) %*% G
  Q <- Q %*% G
  T[k + 1, k] <- 0
  list(T = T, Q = Q)
}

.build_P0 <- function(F_mat, GQG, diffuse_scale = .DIFFUSE_SCALE, tol = .LYAP_TOL) {
  n <- nrow(F_mat)
  P <- diag(1, n)
  for (iter in 1:500) {
    P_new <- F_mat %*% P %*% t(F_mat) + GQG
    if (max(abs(P_new - P)) < tol) return(.sym(P_new))
    if (any(!is.finite(P_new)) || max(abs(P_new)) > 1e20) break
    P <- P_new
  }
  evals <- abs(eigen(F_mat, only.values = TRUE)$values)
  is_unit <- evals > (1 - 1e-6)
  stable_idx <- which(!is_unit)
  P0 <- diag(diffuse_scale, n)
  if (length(stable_idx) > 0) {
    F_s   <- F_mat[stable_idx, stable_idx, drop = FALSE]
    GQG_s <- GQG[stable_idx, stable_idx, drop = FALSE]
    P_s   <- .solve_lyapunov(F_s, GQG_s)
    if (all(is.finite(P_s))) P0[stable_idx, stable_idx] <- .sym(P_s)
  }
  P0
}

.build_P0_Q2 <- function(F_mat, GQG, const = NULL, ur_tol = 1e-6) {
  ## ur_tol kept as 1e-6 — threshold for unit-root detection is model-specific.
  n <- nrow(F_mat)
  if (is.null(const)) const <- rep(0, n)
  
  schur <- Matrix::Schur(Matrix::Matrix(F_mat, sparse = FALSE))
  Ts <- as.matrix(schur@T)
  Qs <- as.matrix(schur@Q)
  
  ur_target <- 0L
  for (i in seq_len(n)) {
    if (abs(diag(Ts)[i] - 1) < ur_tol) {
      j <- i
      while (j > ur_target + 1L) {
        sw <- .swap_schur_11(Ts, Qs, j - 1L)
        Ts <- sw$T;  Qs <- sw$Q
        j <- j - 1L
      }
      ur_target <- ur_target + 1L
    }
  }
  nunit <- ur_target
  U <- Qs
  
  Pa0 <- matrix(0, n, n)
  n_stable <- n - nunit
  if (n_stable > 0L) {
    idx_s <- (nunit + 1L):n
    Ts_s  <- Ts[idx_s, idx_s, drop = FALSE]
    U_inv <- t(U)
    GQG_schur <- U_inv %*% GQG %*% U
    GQG_s <- GQG_schur[idx_s, idx_s, drop = FALSE]
    Pa0[idx_s, idx_s] <- .sym(.solve_lyapunov(Ts_s, GQG_s))
  }
  P0 <- U %*% Pa0 %*% t(U)
  
  Q2 <- if (nunit > 0L) U[, 1:nunit, drop = FALSE] else NULL
  
  x0_unc <- rep(0, n)
  if (n_stable > 0L && any(const != 0)) {
    idx_s <- (nunit + 1L):n
    Ka_s <- as.vector(t(U) %*% const)[idx_s]
    alpha0 <- rep(0, n)
    alpha0[idx_s] <- solve(diag(n_stable) - Ts[idx_s, idx_s, drop = FALSE], Ka_s)
    x0_unc <- as.vector(U %*% alpha0)
  }
  
  list(P0 = .sym(P0), Q2 = Q2, x0_unc = x0_unc, nunit = nunit,
       eigenvalues = diag(Ts))
}


# =============================================================================
# 1.  dynhr_kalman_smoother
# =============================================================================

dynhr_kalman_smoother <- function(Y, F_mat, G_mat, H_mat, Q_mat,
                                  R_mat   = NULL,
                                  const   = NULL,
                                  d_const = NULL,
                                  x0      = NULL,
                                  P0      = NULL,
                                  init_Q2 = NULL) {
  
  TT    <- nrow(Y)
  n_s   <- nrow(F_mat)
  n_obs <- ncol(Y)
  n_shk <- ncol(G_mat)
  
  if (is.null(R_mat))   R_mat   <- matrix(0, n_obs, n_obs)
  if (is.null(const))   const   <- rep(0, n_s)
  if (is.null(d_const)) d_const <- rep(0, n_obs)
  if (is.null(x0))      x0      <- rep(0, n_s)
  
  GQG <- G_mat %*% Q_mat %*% t(G_mat)
  if (is.null(P0)) P0 <- .build_P0(F_mat, GQG)
  
  use_gls   <- !is.null(init_Q2)
  nunit_gls <- if (use_gls) ncol(init_Q2) else 0L
  init_est  <- NULL
  
  .forward_pass <- function(x0_use, accumulate_gls = FALSE) {
    x_pred <- matrix(0, TT, n_s)
    x_filt <- matrix(0, TT, n_s)
    P_pred <- array(0, dim = c(n_s, n_s, TT))
    P_filt <- array(0, dim = c(n_s, n_s, TT))
    loglik <- 0
    x_tt <- x0_use
    P_tt <- P0
    
    if (accumulate_gls) {
      Q2_filt   <- init_Q2
      sumMtSiM  <- matrix(0, nunit_gls, nunit_gls)
      sumMtSiv  <- numeric(nunit_gls)
    }
    
    for (t in seq_len(TT)) {
      x_tp <- as.numeric(F_mat %*% x_tt + const)
      P_tp <- .sym(F_mat %*% P_tt %*% t(F_mat) + GQG)
      x_pred[t, ]   <- x_tp
      P_pred[, , t] <- P_tp
      
      if (accumulate_gls) Q2_pred <- F_mat %*% Q2_filt
      
      obs_avail <- !is.na(Y[t, ])
      n_avail   <- sum(obs_avail)
      
      if (n_avail > 0L) {
        H_t <- H_mat[obs_avail, , drop = FALSE]
        R_t <- R_mat[obs_avail, obs_avail, drop = FALSE]
        d_t <- d_const[obs_avail]
        y_t <- Y[t, obs_avail]
        
        v_t <- y_t - H_t %*% x_tp - d_t
        S_t <- .sym(H_t %*% P_tp %*% t(H_t) + R_t)
        # Use pseudoinverse if S_t is near-singular
        S_inv <- if (rcond(S_t) > .Machine$double.eps) solve(S_t) else .pinv(S_t)
        K_t   <- P_tp %*% t(H_t) %*% S_inv
        
        x_tt <- as.numeric(x_tp + K_t %*% v_t)
        P_tt <- .sym(P_tp - K_t %*% H_t %*% P_tp)
        
        log_det_S <- determinant(S_t, logarithm = TRUE)$modulus[1]
        loglik <- loglik - 0.5 * (n_avail * log(2 * pi) + log_det_S +
                                    as.numeric(t(v_t) %*% S_inv %*% v_t))
        
        if (accumulate_gls) {
          M_t <- H_t %*% Q2_pred
          SiM <- S_inv %*% M_t
          sumMtSiM <- sumMtSiM + crossprod(M_t, SiM)
          sumMtSiv <- sumMtSiv + as.numeric(crossprod(M_t, S_inv %*% v_t))
          Q2_filt  <- Q2_pred - K_t %*% M_t
        }
      } else {
        x_tt <- x_tp
        P_tt <- P_tp
        if (accumulate_gls) Q2_filt <- Q2_pred
      }
      
      x_filt[t, ]   <- x_tt
      P_filt[, , t] <- P_tt
    }
    
    result <- list(x_pred = x_pred, x_filt = x_filt,
                   P_pred = P_pred, P_filt = P_filt, loglik = loglik)
    if (accumulate_gls) {
      result$sumMtSiM <- sumMtSiM
      result$sumMtSiv <- sumMtSiv
    }
    result
  }
  
  pass1 <- .forward_pass(x0, accumulate_gls = use_gls)
  
  if (use_gls) {
    init_est <- as.numeric(solve(pass1$sumMtSiM, pass1$sumMtSiv))
    x0_star  <- x0 + as.numeric(init_Q2 %*% init_est)
    pass2    <- .forward_pass(x0_star, accumulate_gls = FALSE)
    x_pred <- pass2$x_pred;  x_filt <- pass2$x_filt
    P_pred <- pass2$P_pred;  P_filt <- pass2$P_filt
    loglik <- pass2$loglik
  } else {
    x_pred <- pass1$x_pred;  x_filt <- pass1$x_filt
    P_pred <- pass1$P_pred;  P_filt <- pass1$P_filt
    loglik <- pass1$loglik
  }
  
  x_smooth <- matrix(0, TT, n_s)
  P_smooth <- array(0, dim = c(n_s, n_s, TT))
  x_smooth[TT, ] <- x_filt[TT, ]
  P_smooth[, , TT] <- P_filt[, , TT]
  
  for (t in (TT - 1):1) {
    J_t <- .smoother_gain(P_filt[, , t], F_mat, P_pred[, , t + 1])
    x_smooth[t, ] <- x_filt[t, ] +
      J_t %*% (x_smooth[t + 1, ] - x_pred[t + 1, ])
    P_smooth[, , t] <- .sym(
      P_filt[, , t] + J_t %*% (P_smooth[, , t + 1] - P_pred[, , t + 1]) %*% t(J_t)
    )
  }
  
  GtG_inv_Gt <- .pinv(G_mat)
  eps_smooth  <- matrix(0, TT, n_shk)
  for (t in seq_len(TT - 1)) {
    residual <- x_smooth[t + 1, ] - F_mat %*% x_smooth[t, ] - const
    eps_smooth[t, ] <- GtG_inv_Gt %*% residual
  }
  eps_smooth[TT, ] <- NA_real_
  
  list(
    smoothed_states  = x_smooth,
    smoothed_shocks  = eps_smooth,
    filtered_states  = x_filt,
    predicted_states = x_pred,
    smoothed_P       = P_smooth,
    loglik           = loglik,
    init_est         = init_est
  )
}

# =============================================================================
# 2.  build_state_space
# =============================================================================

build_state_space <- function(spec, sigma_u_overrides = NULL) {
  
  core     <- spec$core_states
  st_ident <- spec$state_identities %||% list()
  derived  <- spec$derived_variables %||% list()
  meas     <- spec$measurements
  m_ident  <- spec$meas_identities %||% list()
  
  core_names <- names(core)
  ident_names <- names(st_ident)
  
  I1_states <- core_names[vapply(core, function(s) s$type == "I1", logical(1))]
  lag_names <- paste0(I1_states, "_lag")
  state_names <- c(ident_names, core_names, lag_names)
  n_s <- length(state_names)
  
  expand_weights <- function(wt) {
    out <- numeric(0)
    for (nm in names(wt)) {
      if (nm %in% names(derived)) {
        dv <- derived[[nm]]
        for (dnm in names(dv)) {
          val <- wt[nm] * dv[dnm]
          if (dnm %in% names(out)) { out[dnm] <- out[dnm] + val
          } else { out[dnm] <- val }
        }
      } else {
        if (nm %in% names(out)) { out[nm] <- out[nm] + wt[nm]
        } else { out[nm] <- wt[nm] }
      }
    }
    out
  }
  
  # -- F --------------------------------------------------------------
  F_mat <- matrix(0, n_s, n_s, dimnames = list(state_names, state_names))
  for (nm in core_names) {
    s <- core[[nm]]
    if (s$type == "I0") {
      F_mat[nm, nm] <- s$rho
    } else {
      F_mat[nm, nm] <- 1 + s$rho
      F_mat[nm, paste0(nm, "_lag")] <- -s$rho
    }
  }
  for (nm in I1_states) F_mat[paste0(nm, "_lag"), nm] <- 1
  
  for (nm in ident_names) {
    wt <- expand_weights(st_ident[[nm]])
    wy <- rep(0, n_s); names(wy) <- state_names
    for (snm in names(wt)) if (snm %in% state_names) wy[snm] <- wt[snm]
    F_mat[nm, ] <- as.numeric(crossprod(wy, F_mat))
  }
  
  # -- G --------------------------------------------------------------
  G_mat <- matrix(0, n_s, length(core_names),
                  dimnames = list(state_names, paste0("e_", core_names)))
  for (j in seq_along(core_names)) G_mat[core_names[j], j] <- 1
  
  for (nm in ident_names) {
    wt <- expand_weights(st_ident[[nm]])
    wy <- rep(0, n_s); names(wy) <- state_names
    for (snm in names(wt)) if (snm %in% state_names) wy[snm] <- wt[snm]
    G_mat[nm, ] <- as.numeric(crossprod(wy, G_mat))
  }
  
  # -- Q --------------------------------------------------------------
  sigmas_core <- vapply(core, function(s) s$sigma %||% 1, numeric(1))
  Q_mat <- diag(sigmas_core^2, nrow = length(core_names))
  dimnames(Q_mat) <- list(paste0("e_", core_names), paste0("e_", core_names))
  
  # -- K (constant) ---------------------------------------------------
  K <- rep(0, n_s); names(K) <- state_names
  for (nm in core_names) {
    s <- core[[nm]]
    if (s$type == "I0") K[nm] <- (1 - s$rho) * (s$ss %||% 0)
    else                K[nm] <- (1 - s$rho) * (s$drift %||% 0)
  }
  for (nm in ident_names) {
    wt <- expand_weights(st_ident[[nm]])
    wy <- rep(0, n_s); names(wy) <- state_names
    for (snm in names(wt)) if (snm %in% state_names) wy[snm] <- wt[snm]
    K[nm] <- sum(wy * K)
  }
  
  # -- Measurement system: Am z = Bm DD x + Cm u ---------------------
  obs_names <- names(meas)
  n_obs <- length(obs_names)
  
  tvar_names <- c(state_names, names(derived))
  n_tvar <- length(tvar_names)
  
  DD <- matrix(0, n_tvar, n_s, dimnames = list(tvar_names, state_names))
  for (nm in state_names) DD[nm, nm] <- 1
  for (dnm in names(derived)) {
    dv <- derived[[dnm]]
    for (snm in names(dv)) DD[dnm, snm] <- dv[snm]
  }
  
  shock_names <- paste0("u_", sub("_$", "", obs_names))
  sigma_u_vec <- numeric(n_obs); names(sigma_u_vec) <- shock_names
  for (i in seq_along(obs_names)) {
    sigma_u_vec[i] <- meas[[obs_names[i]]]$sigma_u %||% 0
  }
  
  message("  sigma_u_vec: ", paste(names(sigma_u_vec), round(sigma_u_vec, 2), 
                                   sep="=", collapse=", "))
  
  if (!is.null(sigma_u_overrides)) {
    valid <- intersect(names(sigma_u_overrides), shock_names)
    sigma_u_vec[valid] <- sigma_u_overrides[valid]
  }
  
  
  message("  sigma_u_vec after overrides: ", paste(names(sigma_u_vec), round(sigma_u_vec, 2), 
                                   sep="=", collapse=", "))
  
  
  Am <- matrix(0, n_obs, n_obs, dimnames = list(obs_names, obs_names))
  Bm <- matrix(0, n_obs, n_tvar, dimnames = list(obs_names, tvar_names))
  Cm <- matrix(0, n_obs, n_obs, dimnames = list(obs_names, shock_names))
  
  identity_meas <- character(0)
  
  for (i in seq_along(obs_names)) {
    nm <- obs_names[i]
    m  <- meas[[nm]]
    m_type <- m$type %||% (if (!is.null(m$gap_terms)) "gap" else "direct")
    
    Am[nm, nm] <- 1
    
    if (m_type == "identity") {
      identity_meas <- c(identity_meas, nm)
      wt <- m$weights
      for (wnm in names(wt)) {
        if (wnm %in% names(m_ident)) {
          exp_to <- m_ident[[wnm]]$expand_to
          for (enm in names(exp_to))
            Am[nm, enm] <- Am[nm, enm] - wt[wnm] * exp_to[enm]
        } else {
          Am[nm, wnm] <- Am[nm, wnm] - wt[wnm]
        }
      }
    } else {
      state_nm <- m$state
      if (!is.null(state_nm) && state_nm %in% tvar_names) Bm[nm, state_nm] <- 1
      Cm[nm, shock_names[i]] <- sigma_u_vec[i]
      if (!is.null(m$gap_terms)) {
        for (gt_nm in names(m$gap_terms)) {
          gt_coef <- m$gap_terms[[gt_nm]]
          Am[nm, gt_nm] <- Am[nm, gt_nm] - gt_coef
          gt_state <- sub("_$", "", gt_nm)
          if (gt_state %in% tvar_names)
            Bm[nm, gt_state] <- Bm[nm, gt_state] - gt_coef
        }
      }
    }
  }
  
  Am_inv <- solve(Am)
  H_mat  <- Am_inv %*% Bm %*% DD
  HH     <- Am_inv %*% Cm
  R_mat  <- HH %*% t(HH)
  d_const <- rep(0, n_obs); names(d_const) <- obs_names
  dimnames(H_mat) <- list(obs_names, state_names)
  
  # -- P0, Q2, x0_unc ------------------------------------------------
  GQG <- G_mat %*% Q_mat %*% t(G_mat)
  pq  <- .build_P0_Q2(F_mat, GQG, const = K)
  
  names(pq$x0_unc)  <- state_names
  dimnames(pq$P0)   <- list(state_names, state_names)
  if (!is.null(pq$Q2)) rownames(pq$Q2) <- state_names
  
  list(
    F_mat = F_mat, G_mat = G_mat, H_mat = H_mat,
    Q_mat = Q_mat, R_mat = R_mat,
    const = K, d_const = d_const,
    P0 = pq$P0, Q2 = pq$Q2, x0_unc = pq$x0_unc,
    state_names = state_names, obs_names = obs_names,
    n_state = n_s, n_obs = n_obs, n_shock = length(core_names),
    identity_meas = identity_meas,
    nunit = pq$nunit, eigenvalues = pq$eigenvalues
  )
}

# --- End of dynhr_backend.R -------------------------------------------------
