## R/solve-helpers.R
## --------------------------------------------------------------------------
## Small numerical utilities used by the perturbation solver: safe SVD-based
## pseudoinverse, etc.
##
## Phase-1 split from perturbation-monolith.R (no logic changes).
## --------------------------------------------------------------------------

# ---- Utility: safe matrix inverse via SVD ----

#' SVD-based pseudoinverse with fallback
#'
#' @param M     Square or rectangular matrix
#' @param tol   Singular value threshold below which values are treated as zero
#' @return Pseudoinverse of M
#' @noRd
.safe_inv <- function(M, tol = 1e-10) {
  s <- svd(M)
  ## Numerically identical to ifelse(s$d > tol, 1/s$d, 0) but without the
  ## ifelse overhead; and scale rows of t(u) by d_inv directly instead of
  ## forming a diagonal matrix (diag(d) %*% X == d * X), saving an allocation
  ## and a matrix multiply.
  d_inv <- 1 / s$d
  d_inv[s$d <= tol] <- 0
  s$v %*% (d_inv * t(s$u))
}
