// qr_implicit.cpp -- static-variable QR elimination that never forms Q.
//
// qr_static_transform_cpp() (static_elim.cpp) builds the full n x n Q with
// dorgqr and multiplies Q' into four n x k blocks: O(n^3) for Q plus
// O(n^2 k) for the products. This variant factors the n x n_s static block
// with dgeqrf and applies the Householder reflectors directly to each block
// with a blocked reflector (dlarft + dlarfb) -- the qr.qty() route -- so the
// cost is O(n n_s^2 + n n_s k) and Q is never materialised. The reflectors are the same ones the full-Q
// routine uses, so Q'F agrees to round-off.
//
// Plain LAPACK on column-major arrays; no Armadillo in this translation unit
// (its LAPACK prototypes clash with the explicit ones, see ordered_qz.cpp).

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <vector>
#include <algorithm>
#ifndef FCONE
# define FCONE
#endif

using Rcpp::NumericMatrix;
using Rcpp::List;
using Rcpp::_;

// Q' C for the reflectors held in `qr` (dgeqrf output), applied as ONE block
// reflector H = I - V T V': Q' = H' = I - V T' V', two matrix products on the
// n x n_s panel V rather than n_s rank-one updates. dormqr switches to its
// unblocked, BLAS-2 path whenever n_s is below its block size (32), which is
// exactly the common case of a model with a few dozen static variables; the
// explicit dlarft/dlarfb pair keeps the work in BLAS-3 for every n_s.
static NumericMatrix apply_qt(const std::vector<double>& qr, int n, int ns,
                              const std::vector<double>& T,
                              const NumericMatrix& B,
                              std::vector<double>& work) {
  NumericMatrix C = Rcpp::clone(B);
  const int k = C.ncol();
  if (k == 0 || n == 0 || ns == 0) return C;
  const char side = 'L', trans = 'T', direct = 'F', storev = 'C';
  int ldwork = k;
  F77_CALL(dlarfb)(&side, &trans, &direct, &storev, &n, &k, &ns, qr.data(),
                   &n, T.data(), &ns, C.begin(), &n, work.data(), &ldwork
                   FCONE FCONE FCONE FCONE);
  return C;
}

// [[Rcpp::export]]
List qr_static_transform_implicit_cpp(const NumericMatrix& f_static,
                                      const NumericMatrix& f_minus_r,
                                      const NumericMatrix& f_zero_r,
                                      const NumericMatrix& f_plus_r,
                                      const NumericMatrix& f_exo_r) {
  const int n  = f_static.nrow();
  const int ns = f_static.ncol();
  std::vector<double> qr(f_static.begin(), f_static.end());
  std::vector<double> tau(ns > 0 ? ns : 1);
  int info = 0, lwork = -1;
  double wq = 0.0;
  F77_CALL(dgeqrf)(&n, &ns, qr.data(), &n, tau.data(), &wq, &lwork, &info);
  lwork = static_cast<int>(wq);
  std::vector<double> work(lwork > 1 ? lwork : 1);
  F77_CALL(dgeqrf)(&n, &ns, qr.data(), &n, tau.data(), work.data(), &lwork,
                   &info);
  if (info != 0)
    Rcpp::stop("qr_static_transform_implicit_cpp: dgeqrf failed to factorize f_static");

  // Triangular factor of the block reflector.
  std::vector<double> T((size_t) ns * ns, 0.0);
  if (ns > 0) {
    const char direct = 'F', storev = 'C';
    F77_CALL(dlarft)(&direct, &storev, &n, &ns, qr.data(), &n, tau.data(),
                     T.data(), &ns FCONE FCONE);
  }

  const int kmax = std::max(std::max(f_minus_r.ncol(), f_zero_r.ncol()),
                            std::max(f_plus_r.ncol(), f_exo_r.ncol()));
  std::vector<double> w((size_t) std::max(kmax, 1) * std::max(ns, 1));

  return List::create(
    _["Qf_minus"] = apply_qt(qr, n, ns, T, f_minus_r, w),
    _["Qf_zero"]  = apply_qt(qr, n, ns, T, f_zero_r,  w),
    _["Qf_plus"]  = apply_qt(qr, n, ns, T, f_plus_r,  w),
    _["Qf_exo"]   = apply_qt(qr, n, ns, T, f_exo_r,   w)
  );
}
