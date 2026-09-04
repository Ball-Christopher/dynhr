###############################################################################
#  dynhr_backend.R
#  General-purpose state-space toolkit for dynhr.
#  Dependencies: base R + Matrix (for the Schur decomposition in .build_P0_Q2).
###############################################################################

# =============================================================================
# 0.  HELPERS
# =============================================================================

## `%||%` is now defined once in R/utils-pipe.R (loaded first via Collate).
## Local definition removed during phase-0 packaging.

.sym <- function(X) (X + t(X)) / 2

## NOTE (D2, 2026-09-02): .solve_lyapunov() used to be a SECOND, direct-
## kronecker Lyapunov solver defined here.  It is now a one-line alias of the
## package's single solve_lyapunov() (R/solve-helpers.R), which keeps the same
## NaN-on-nonstationary contract but also converges on highly non-normal but
## stable transition matrices where the kron rcond gate falsely fired.


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


# --- End of dynhr_backend.R -------------------------------------------------
