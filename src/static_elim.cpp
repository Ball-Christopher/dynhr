// static_elim.cpp -- C++ static-variable QR elimination for the first-order
// perturbation solver (dynhr).
//
// Replaces the hot R block in .solve_from_system() (Step 3): a full QR of the
// static-variable columns of the reordered Jacobian, followed by transforming
// the four system matrices by Q'. The pure-R path (qr/qr.Q + four t(Q) %*% .)
// stays in R/solve-perturbation.R as a bit-parity fallback.
//
// IMPORTANT: this uses *non-pivoted* QR (arma::qr), matching R's default qr()
// (LINPACK dqrdc2) which does no column pivoting in the full-rank case. The
// pivot is therefore the identity; the R caller sets inv_piv <- seq_len(n_s).
// Householder sign-convention differences between LAPACK and dqrdc2 cancel
// through the downstream static-recovery math, so ghx/ghu are unchanged
// (asserted by test-static-elim-parity.R). Column *pivoting*, by contrast,
// does change the result, which is exactly why we avoid it here.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// [[Rcpp::export]]
List qr_static_transform_cpp(const arma::mat& f_static,
                             const arma::mat& f_minus_r,
                             const arma::mat& f_zero_r,
                             const arma::mat& f_plus_r,
                             const arma::mat& f_exo_r) {
  arma::mat Q, R;
  bool ok = arma::qr(Q, R, f_static);   // full (complete) QR; Q is n x n
  if (!ok)
    Rcpp::stop("qr_static_transform_cpp: arma::qr failed to factorize f_static");

  const arma::mat Qt = Q.t();
  return List::create(
    _["Qf_minus"] = Qt * f_minus_r,
    _["Qf_zero"]  = Qt * f_zero_r,
    _["Qf_plus"]  = Qt * f_plus_r,
    _["Qf_exo"]   = Qt * f_exo_r
  );
}
