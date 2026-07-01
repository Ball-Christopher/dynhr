// ordered_qz.cpp -- robust ordered generalized Schur (QZ) decomposition.
//
// Calls LAPACK dgges directly with a custom SORT callback that places stable
// eigenvalues (|alpha/beta| < 1 + 1e-6) in the top-left, exactly mirroring
// Dynare's mjdgges and the threshold used by dynhr's previous QZ::qz() +
// QZ::qz.dtgsen() path (select = |eig| <= 1 + 1e-6).
//
// Motivation: QZ::qz.dtgsen (the post-hoc eigenvalue reorder in the QZ package,
// v0.2-4) has an input-INDEPENDENT memory-corruption bug -- hammering it on a
// fixed, well-conditioned pencil segfaults after enough calls (missing PROTECT /
// workspace under-allocation in the package C glue). Ordering *during* the
// decomposition via dgges's SELCTG callback avoids dtgsen entirely. The R path
// (QZ::qz + qz.dtgsen) is retained as a fallback in R/solve-qz.R for builds
// without the compiled DLL.
//
// Returns ALPHAR/ALPHAI/BETA so the caller's eigenvalue classification and
// Blanchard-Kahn counting logic are unchanged (dgges computes the same Schur
// form and generalized eigenvalues as QZ::qz; only the ordering routine differs).

// This file calls LAPACK's dgges_ directly (declared below) and uses Armadillo
// only as a matrix container -- no Armadillo LAPACK-backed ops. Suppress
// Armadillo's own LAPACK prototypes so its dgges_ declaration (which carries
// hidden Fortran string-length args) does not conflict with ours. Newer
// RcppArmadillo (>= 15) + clang turn that clash into a hard error.
#define ARMA_DONT_USE_LAPACK

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// LAPACK dgges (double real generalized Schur). Declared explicitly so we do not
// depend on Armadillo's internal LAPACK wrappers (signatures vary across
// versions). R ships LAPACK (linked via -lRlapack).
typedef int (*dgges_selctg)(const double*, const double*, const double*);

extern "C" {
  void dgges_(const char* JOBVSL, const char* JOBVSR, const char* SORT,
              dgges_selctg SELCTG, const int* N,
              double* A, const int* LDA, double* B, const int* LDB,
              int* SDIM, double* ALPHAR, double* ALPHAI, double* BETA,
              double* VSL, const int* LDVSL, double* VSR, const int* LDVSR,
              double* WORK, const int* LWORK, int* BWORK, int* INFO);
}

// Stability threshold for the SORT callback. dgges's SELCTG takes no extra
// arguments, so the criterium is passed via a translation-unit-local variable.
// dgges runs synchronously and single-threaded from our call, so this is safe.
static double qz_critmod = 1.0 + 1e-6;

// Select "stable" eigenvalues: |alpha/beta| < critmod, written without division
// so beta == 0 (infinite eigenvalue) correctly classifies as unstable (returns 0).
static int qz_selctg(const double* ar, const double* ai, const double* beta) {
  const double lhs = (*ar) * (*ar) + (*ai) * (*ai);
  const double rhs = qz_critmod * qz_critmod * (*beta) * (*beta);
  return (lhs < rhs) ? 1 : 0;
}

// [[Rcpp::export]]
List ordered_qz_cpp(const arma::mat& E, const arma::mat& D,
                    double critmod = 1.0000010) {
  const int N = static_cast<int>(E.n_rows);
  List out;
  if (N == 0) {
    return List::create(_["ok"] = true,
                        _["S"] = arma::mat(0, 0), _["T"] = arma::mat(0, 0),
                        _["Q"] = arma::mat(0, 0), _["Z"] = arma::mat(0, 0),
                        _["ALPHAR"] = arma::vec(), _["ALPHAI"] = arma::vec(),
                        _["BETA"] = arma::vec(), _["sdim"] = 0, _["info"] = 0);
  }
  if (static_cast<int>(E.n_cols) != N ||
      static_cast<int>(D.n_rows) != N || static_cast<int>(D.n_cols) != N)
    Rcpp::stop("ordered_qz_cpp: E and D must be square and conformable");

  qz_critmod = critmod;

  // dgges overwrites A and B with the (quasi-triangular) Schur factors S and T.
  // The pencil orientation matches QZ::qz(E, D): A = E, B = D.
  arma::mat A = E;   // -> S on exit
  arma::mat B = D;   // -> T on exit
  arma::mat VSL(N, N), VSR(N, N);
  arma::vec ALPHAR(N), ALPHAI(N), BETA(N);
  std::vector<int> BWORK(N);

  const char JOBVSL = 'V', JOBVSR = 'V', SORT = 'S';
  int SDIM = 0, INFO = 0;

  // Workspace query.
  double work_query = 0.0;
  int LWORK = -1;
  dgges_(&JOBVSL, &JOBVSR, &SORT, &qz_selctg, &N,
         A.memptr(), &N, B.memptr(), &N, &SDIM,
         ALPHAR.memptr(), ALPHAI.memptr(), BETA.memptr(),
         VSL.memptr(), &N, VSR.memptr(), &N,
         &work_query, &LWORK, BWORK.data(), &INFO);
  if (INFO != 0)
    return List::create(_["ok"] = false, _["info"] = INFO, _["stage"] = "query");

  LWORK = static_cast<int>(work_query);
  if (LWORK < 8 * N + 16) LWORK = 8 * N + 16;   // LAPACK minimum guard
  std::vector<double> WORK(LWORK);

  dgges_(&JOBVSL, &JOBVSR, &SORT, &qz_selctg, &N,
         A.memptr(), &N, B.memptr(), &N, &SDIM,
         ALPHAR.memptr(), ALPHAI.memptr(), BETA.memptr(),
         VSL.memptr(), &N, VSR.memptr(), &N,
         WORK.data(), &LWORK, BWORK.data(), &INFO);

  // INFO in [1, N]: QZ iteration failed. INFO == N+1: reordering failed but
  // eigenvalues are still correct. INFO == N+2: after reordering, roundoff made
  // a selected pair no longer satisfy the criterion (cosmetic). Treat <=0 and
  // N+1/N+2 as usable; hard failures (1..N) signal the caller to fall back.
  bool ok = (INFO == 0) || (INFO == N + 1) || (INFO == N + 2);

  return List::create(_["ok"] = ok, _["info"] = INFO, _["sdim"] = SDIM,
                      _["S"] = A, _["T"] = B, _["Q"] = VSL, _["Z"] = VSR,
                      _["ALPHAR"] = ALPHAR, _["ALPHAI"] = ALPHAI,
                      _["BETA"] = BETA);
}
