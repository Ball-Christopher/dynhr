// mcp_sparse_solve.cpp -- Sparse linear solver for the MCP stacked system
//
// Provides one exported function, mcp_sparse_solve_cpp(), used by
// mcp_solve_path() when "backend = 'rcpp'" is selected.
// The pure-R fallback (Matrix::solve on dgCMatrix) is the default.
//
// Uses Armadillo's spsolve() which auto-detects the best available solver:
//   - SuperLU (if installed and ARMA_USE_SUPERLU is defined)
//   - Built-in LAPACK-based sparse solver (always available)
//
// To enable SuperLU on this system, install libsuperlu and set:
//   Sys.setenv("PKG_LIBS" = "-lsuperlu")
// before compilation.  Without SuperLU, the built-in solver is used
// (still faster than R's Matrix::sparseQR for most systems).
//
// Ref: https://gallery.rcpp.org/articles/armadillo-with-superlu/
//
// Args: triplet-form sparse matrix (i, j, x vectors) + dense RHS vector.
// Returns: solution vector.

// SuperLU note: to enable SuperLU (2-5x faster for sparse systems),
// install libsuperlu and add to src/Makevars.win:
//   PKG_LIBS += -lsuperlu
// then uncomment the following line.  Without SuperLU, Armadillo uses
// its built-in LAPACK-based sparse solver (still faster than Matrix::sparseQR).
// #define ARMA_USE_SUPERLU 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::_;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

// [[Rcpp::export]]
NumericVector mcp_sparse_solve_cpp(IntegerVector i_idx,
                                   IntegerVector j_idx,
                                   NumericVector x_vals,
                                   NumericVector rhs,
                                   int n_dim) {
  // ----- Input validation -----
  if (n_dim <= 0)
    Rcpp::stop("mcp_sparse_solve_cpp: n_dim must be >= 1");
  if (i_idx.size() != j_idx.size())
    Rcpp::stop("mcp_sparse_solve_cpp: i_idx and j_idx must have same length");
  if (i_idx.size() != x_vals.size())
    Rcpp::stop("mcp_sparse_solve_cpp: i_idx and x_vals must have same length");
  if (rhs.size() != n_dim)
    Rcpp::stop("mcp_sparse_solve_cpp: rhs must have length n_dim");

  arma::uword n = static_cast<arma::uword>(n_dim);
  arma::uword nnz = static_cast<arma::uword>(i_idx.size());

  // Build sparse matrix from triplets (1-based R indices → 0-based C++)
  arma::umat locations(2, nnz);
  arma::vec values(nnz);

  for (arma::uword k = 0; k < nnz; ++k) {
    int row = i_idx[k];
    int col = j_idx[k];
    if (row < 1 || row > n_dim)
      Rcpp::stop("mcp_sparse_solve_cpp: row index %d out of range [1, %d]", row, n_dim);
    if (col < 1 || col > n_dim)
      Rcpp::stop("mcp_sparse_solve_cpp: col index %d out of range [1, %d]", col, n_dim);
    locations(0, k) = static_cast<arma::uword>(row - 1);  // row (0-based)
    locations(1, k) = static_cast<arma::uword>(col - 1);  // col (0-based)
    values[k] = x_vals[k];
  }

  arma::sp_mat A(locations, values, n, n);
  arma::vec b = arma::vec(rhs.begin(), n);

  // Solve: A * x = b using the default solver.
  // Armadillo auto-selects: SuperLU (if available) → built-in LAPACK.
  // We use the two-argument form (returns bool, throws no exceptions).
  arma::vec x;
  bool success = arma::spsolve(x, A, b);

  if (!success) {
    // Fall back to iterative least-squares for rank-deficient systems
    success = arma::spsolve(x, A, b, "lsqr");
  }

  if (!success) {
    Rcpp::stop("mcp_sparse_solve_cpp: solve failed (system may be singular)");
  }

  return Rcpp::wrap(x);
}
