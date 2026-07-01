## R/solve-perturbation-sylvester.R
## --------------------------------------------------------------------------
## Compact generalized-Sylvester solver for the deterministic Kronecker
## systems that arise at every perturbation order:
##
##     A_L * X + fp * X * hx^{⊗k} = RHS,        X is n × ns^k.
##
## The dense reference solver (.solve_kron_direct, order4 file) forms the full
## ns^k × ns^k matrix hx^{⊗k} and runs a real Schur factorisation on it — an
## O(ns^{3k}) operation plus an O(ns^{2k}) dense object.  This routine instead
## triangularises only the small n × n pencil (A_L, fp) and diagonalises only
## the small ns × ns state matrix hx, applying hx^{⊗k} implicitly by mode-wise
## tensor contraction.  Cost drops to roughly O(n^3 + n^2 · k · ns^{k+1}); the
## giant Kronecker object is never formed.
##
## Correctness notes (these are the traps the earlier per-order attempts hit):
##   * COMPLEX QZ.  qz() on real (A_L, fp) returns a *real* quasi-triangular
##     pencil with 2×2 blocks for complex-conjugate eigenpairs; the old per-row
##     back-substitution treated every diagonal as 1×1 and silently dropped the
##     2×2 coupling.  We feed complex inputs so S, T are strictly upper
##     triangular and the scalar back-substitution is exact.
##   * APPLY DIRECTION.  Solving y·(αI + β·hx^{⊗k}) = b via hx = V Λ V^{-1}
##     needs  b̃ = b·V^{⊗k},  then divide by (α + β·λ-grid),  then ·(V^{-1})^{⊗k}.
##     The order-5 helper had V and V^{-1} swapped.
##   * COMPLEX-SAFE CONTRACTION.  Mode application uses t(M) (plain transpose),
##     never crossprod()/Conj — eigenvectors of a non-symmetric hx are complex.
##   * DIAGONALISABILITY GUARD.  If V is ill-conditioned (defective hx) the
##     eigen route is unstable, so we fall back to the dense Schur solver, which
##     needs no diagonalisability.  The dense path therefore remains the safety
##     net for every solve.
## --------------------------------------------------------------------------


#' Apply the k-fold Kronecker product M^{⊗k} to a row vector y (length ns^k).
#'
#' Computes  y · (M ⊗ M ⊗ ... ⊗ M)  (k factors) by reshaping y as a k-mode
#' tensor and contracting each mode with M.  Complex-safe: uses t(M), not
#' crossprod(), so it is correct for the complex eigenvectors of a
#' non-symmetric state matrix.
#'
#' Layout: y is column-major over the grid (mode 1 fastest), matching
#' array(y, dim = rep(ns, k)) and the .quint_names / odometer conventions
#' used by the forcing assemblers.
#'
#' @param y Numeric or complex vector of length ns^k.
#' @param M ns × ns matrix (real or complex).
#' @param k Kronecker power (k >= 1).
#' @return Vector of length ns^k: y · M^{⊗k}.
#' @noRd
.apply_kronk <- function(y, M, k) {
  ns <- nrow(M)
  if (k <= 0L) return(y)
  tM <- t(M)                                   # hoisted: one transpose, not k
  if (k == 1L) return(as.vector(tM %*% matrix(y, ns, 1L)))

  dimk <- rep(ns, k)
  cols <- ns^(k - 1L)
  Tk   <- array(y, dim = dimk)
  for (mode in seq_len(k)) {
    perm <- c(mode, seq_len(k)[-mode])
    Tk <- aperm(array(tM %*% matrix(aperm(Tk, perm), ns, cols), dim = dimk),
                order(perm))
  }
  as.vector(Tk)
}


#' Build the eigenvalue grid λ_{j1}·λ_{j2}·…·λ_{jk} over the ns^k Kronecker
#' index set, in the same column-major (mode-1-fastest) order as .apply_kronk.
#' @noRd
.kron_lambda_grid <- function(lambda, k) {
  g <- lambda
  if (k >= 2L) for (a in 2:k) g <- as.vector(outer(g, lambda))
  g
}


#' Compact generalized-Sylvester solve:  A_L·X + fp·X·hx^{⊗k} = RHS.
#'
#' Drop-in replacement for \code{.solve_kron_direct} with identical sign and
#' layout conventions (X is n × ns^k; RHS already carries its sign).  Avoids
#' forming hx^{⊗k}.
#'
#' The eigenbasis route diagonalises hx, so its accuracy degrades like
#' kappa(V)^k when the eigenvector matrix V is non-orthogonal (typical for
#' Blanchard-Kahn state transitions).  We therefore wrap it in iterative
#' refinement: the eigen factors are reused to solve the correction equation
#' A_L·dX + fp·dX·hx^{⊗k} = residual, recovering full accuracy in a few cheap
#' steps.  A final residual check GUARANTEES correctness — if refinement cannot
#' drive the residual below \code{accept_tol}, we fall back to the dense Schur
#' solver (which needs no diagonalisability).  The dense path is thus always the
#' safety net, while the common (well-conditioned) case stays on the fast path.
#'
#' @param A_L     n × n effective feedback matrix.
#' @param fp      n × n f_plus matrix.
#' @param hx      ns × ns state-transition matrix.
#' @param k       Kronecker power (k >= 1).
#' @param RHS     n × ns^k right-hand side.
#' @param cond_tol Ceiling on kappa(V)^k; above it we go straight to the dense
#'   solver rather than waste refinement steps (default 1e13).
#' @param tol     Target relative residual for early exit (default 1e-11).
#' @param accept_tol Relative residual that must be met to accept the compact
#'   result; otherwise fall back to dense (default 1e-8).
#' @param maxit   Maximum iterative-refinement steps (default 8).
#' @param verbose Logical.
#' @return n × ns^k real matrix X.
#' @noRd
.solve_kron_compact <- function(A_L, fp, hx, k, RHS, cond_tol = 1e13,
                                tol = 1e-11, accept_tol = 1e-8, maxit = 8L,
                                verbose = FALSE) {
  n  <- nrow(A_L)
  ns <- nrow(hx)
  m  <- ns^k
  stopifnot(ncol(RHS) == m)

  if (k < 1L) stop(".solve_kron_compact: k must be >= 1")
  if (ns == 0L) return(matrix(0, n, 0L))

  dense <- function() .solve_kron_direct(A_L, fp, hx, k, RHS)

  # QZ must be available; otherwise use the dense reference solver.
  if (!requireNamespace("QZ", quietly = TRUE)) {
    if (verbose) cat("  QZ unavailable; using dense Schur solver.\n")
    return(dense())
  }

  # --- Eigendecomposition of the small state matrix, with stability guard ---
  # The eigenbasis solve amplifies roundoff by ~kappa(V)^k (each of the k
  # Kronecker modes contributes one factor of kappa(V)).  Iterative refinement
  # recovers accuracy only while that amplification stays below 1/eps, so we
  # guard on kappa(V)^k, not kappa(V).  Clustered/near-defective eigenvalues
  # (common in tightly-calibrated models) blow this up — those fall back to the
  # dense Schur solver, which is unconditionally stable.
  eg    <- eigen(hx)
  V     <- eg$vectors
  lam   <- eg$values
  V_inv <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(V_inv) ||
      !all(is.finite(as.numeric(rbind(Re(V_inv), Im(V_inv))))) ||
      kappa(V)^k > cond_tol) {
    if (verbose) cat(sprintf(
      "  kappa(V)^k = %.2e exceeds %.2e; using dense Schur solver.\n",
      kappa(V)^k, cond_tol))
    return(dense())
  }

  # --- Complex generalized Schur of (A_L, fp): A_L = Q S Z^H, fp = Q T Z^H ---
  # Force complex so S, T are strictly upper triangular (no 2×2 blocks).
  qzr <- tryCatch(QZ::qz(A_L + 0i, fp + 0i), error = function(e) NULL)
  if (is.null(qzr)) {
    if (verbose) cat("  Complex QZ failed; using dense Schur solver.\n")
    return(dense())
  }
  S <- qzr$S; T <- qzr$T; Q <- qzr$Q; Z <- qzr$Z
  Qh <- Conj(t(Q))
  lam_grid <- .kron_lambda_grid(lam, k)        # length m, complex

  # Solve A_L·X + fp·X·hx^{⊗k} = B for given B, via the triangular (S,T) pencil
  # and the eigenbasis of hx.  Returns a real n × m matrix.
  solve_once <- function(B) {
    G <- Qh %*% B                               # S W + T W C = G,  W = Z^H X
    W <- matrix(0 + 0i, n, m)
    for (i in n:1) {
      rhs_i <- G[i, ]
      if (i < n) {
        for (j in (i + 1L):n) {
          wj <- W[j, ]
          rhs_i <- rhs_i - S[i, j] * wj - T[i, j] * .apply_kronk(wj, hx, k)
        }
      }
      a <- S[i, i]; b <- T[i, i]
      bt <- .apply_kronk(rhs_i, V, k)           # rhs_i · V^{⊗k}
      W[i, ] <- .apply_kronk(bt / (a + b * lam_grid), V_inv, k)  # · (V^{-1})^{⊗k}
    }
    Re(Z %*% W)
  }

  # Apply the LHS operator  A_L·X + fp·X·hx^{⊗k}  (for residual computation).
  apply_lhs <- function(X) {
    XK <- matrix(0, n, m)
    for (rr in seq_len(n)) XK[rr, ] <- Re(.apply_kronk(X[rr, ], hx, k))
    A_L %*% X + fp %*% XK
  }

  scale  <- max(1, max(abs(RHS)))
  X      <- solve_once(RHS)
  r_best <- Inf
  X_best <- X
  for (it in 0:maxit) {
    Rr <- RHS - apply_lhs(X)
    rn <- max(abs(Rr))
    if (rn < r_best) { r_best <- rn; X_best <- X }
    if (rn <= tol * scale || it == maxit) break
    X <- X + solve_once(Rr)                      # iterative refinement
  }

  if (r_best > accept_tol * scale) {
    if (verbose) cat(sprintf(
      "  compact residual %.2e > accept %.2e; using dense Schur solver.\n",
      r_best, accept_tol * scale))
    return(dense())
  }
  rownames(X_best) <- rownames(A_L)
  X_best
}
