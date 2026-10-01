// pskf_lapack.cpp -- the exact LAPACK routines behind the R-level pieces the
// C++ PSKF recursion mirrors: eigen(symmetric = TRUE) (dsyevr) and
// chol(pivot = TRUE) + chol2inv() (dpstrf, dpotri). Calling the same
// routines with the same arguments as R gives the same rounding, so the
// pseudoinverse / Cholesky inverse of the recursion agree with the R path to
// the last bits rather than to the conditioning of the matrix.
//
// No Armadillo here: its own LAPACK prototypes clash with explicit ones (see
// ordered_qz.cpp), so this file only sees plain arrays.

#define USE_FC_LEN_T
#include <R.h>
#include <R_ext/RS.h>
#include <vector>
#include "pskf_lapack.h"
#include <R_ext/Lapack.h>   // dsyevr / dpstrf / dpotri with FCLEN hidden lengths
#ifndef FCONE
# define FCONE
#endif

// R's La_rs (eigen(x, symmetric = TRUE)): dsyevr, jobz "V", range "A",
// uplo "L", abstol 0; values returned in DECREASING order, vectors' columns
// reversed to match.
bool pskf_sym_eigen(int n, const double* A, double* values, double* vectors) {
  std::vector<double> a(A, A + (size_t) n * n), w(n), z((size_t) n * n);
  std::vector<int> isuppz(2 * (size_t) n);
  const char jobz = 'V', range = 'A', uplo = 'L';
  const double vl = 0.0, vu = 0.0, abstol = 0.0;
  const int il = 0, iu = 0;
  int m = 0, info = 0, lwork = -1, liwork = -1, iwork_q = 0;
  double work_q = 0.0;
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, a.data(), &n, &vl, &vu, &il, &iu,
                   &abstol, &m, w.data(), z.data(), &n, isuppz.data(), &work_q,
                   &lwork, &iwork_q, &liwork, &info FCONE FCONE FCONE);
  if (info != 0) return false;
  lwork = (int) work_q;
  liwork = iwork_q;
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, a.data(), &n, &vl, &vu, &il, &iu,
                   &abstol, &m, w.data(), z.data(), &n, isuppz.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE
                   FCONE);
  if (info != 0) return false;
  for (int j = 0; j < n; ++j) {
    values[j] = w[n - 1 - j];
    for (int i = 0; i < n; ++i)
      vectors[i + (size_t) j * n] = z[i + (size_t) (n - 1 - j) * n];
  }
  return true;
}

// chol(C, pivot = TRUE) then chol2inv(): out(piv[a], piv[b]) = the inverse of
// C[piv, piv] scattered back. Returns false when the pivoted factorisation is
// rank deficient (rank < n) or LAPACK reports an error.
bool pskf_pivchol_inverse(int n, const double* C, double* out) {
  std::vector<double> a((size_t) n * n, 0.0), work(2 * (size_t) n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i <= j; ++i) a[i + (size_t) j * n] = C[i + (size_t) j * n];
  std::vector<int> piv(n);
  const char uplo = 'U';
  double tol = -1.0;   // non-const: R's prototype takes double*
  int rank = 0, info = 0;
  F77_CALL(dpstrf)(&uplo, &n, a.data(), &n, piv.data(), &rank, &tol,
                   work.data(), &info FCONE);
  if (info < 0 || rank != n) return false;
  F77_CALL(dpotri)(&uplo, &n, a.data(), &n, &info FCONE);
  if (info != 0) return false;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i <= j; ++i) {
      const double v = a[i + (size_t) j * n];
      const size_t pi = piv[i] - 1, pj = piv[j] - 1;
      out[pi + pj * n] = v;
      out[pj + pi * n] = v;
    }
  return true;
}
