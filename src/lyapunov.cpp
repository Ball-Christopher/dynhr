// lyapunov.cpp -- the doubling iteration behind solve_lyapunov(), in C++.
//
// Mirrors the R loop in solve_lyapunov() (R/solve-helpers.R) step for step:
//   X_{k+1} = X_k + A_k X_k A_k',   A_{k+1} = A_k A_k,
// stopping when max|X_{k+1} - X_k| <= tol * max|X_{k+1}| (RELATIVE, so the
// solution is scale-equivariant), and giving up (converged = FALSE) on any
// non-finite entry so the R caller can run its stability gate and the
// vec/kronecker fallback. Only the iteration lives here: the zero-B shortcut,
// the spectral-radius gate and the fallback stay in R.
//
// The matrix products keep the R association ((A X) A') and go straight to
// dgemm, exactly as R's %*% does. Armadillo's own operator* would multiply
// matrices of order <= 4 with an emulated loop whose summation order differs
// from the BLAS routine, which moves the result in the last bit and, through
// a likelihood, the last bits of a posterior -- so the BLAS call is made
// explicitly and the solution is bit-identical to the R loop.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// C = A B through the BLAS (dgemm), no small-matrix shortcut.
static void blas_mm(const arma::mat& A, const arma::mat& B, arma::mat& C) {
  const arma::blas_int m = arma::blas_int(A.n_rows);
  const arma::blas_int n = arma::blas_int(B.n_cols);
  const arma::blas_int k = arma::blas_int(A.n_cols);
  C.set_size(A.n_rows, B.n_cols);
  if (C.n_elem == 0) return;
  const char t = 'N';
  const double one = 1.0, zero = 0.0;
  arma::blas::gemm<double>(&t, &t, &m, &n, &k, &one, A.memptr(), &m,
                           B.memptr(), &k, &zero, C.memptr(), &m);
}

// [[Rcpp::export]]
List lyapunov_doubling_cpp(const arma::mat& A, const arma::mat& B,
                           int max_iter, double tol) {
  arma::mat X = B;
  arma::mat A_pow = A;
  arma::mat AX, X_new, At, A_next;
  bool converged = false;
  int iters = 0;
  for (int iter = 0; iter < max_iter; ++iter) {
    iters = iter + 1;
    blas_mm(A_pow, X, AX);
    At = A_pow.t();
    blas_mm(AX, At, X_new);
    X_new += X;
    if (!X_new.is_finite()) break;
    const double diff = arma::abs(X_new - X).max();
    if (!std::isfinite(diff)) break;
    if (diff <= tol * arma::abs(X_new).max()) {
      // The R loop returns the PREVIOUS iterate X here, not X_new.
      converged = true;
      break;
    }
    blas_mm(A_pow, A_pow, A_next);
    A_pow = A_next;
    if (!A_pow.is_finite()) break;
    X = X_new;
  }
  return List::create(_["X"] = X, _["converged"] = converged,
                      _["iter"] = iters);
}
