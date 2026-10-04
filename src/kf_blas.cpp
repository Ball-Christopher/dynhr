// kf_blas.cpp -- see kf_blas.h. No Armadillo in this translation unit.

// Fortran hidden string-length arguments: USE_FC_LEN_T must be defined before
// any R header so that Rconfig.h defines FC_LEN_T and R_ext/RS.h defines FCONE
// (Writing R Extensions, "Fortran character strings").
#define USE_FC_LEN_T
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include "kf_blas.h"

#ifndef FCONE
# define FCONE
#endif

void kf_blas_trmm_right_lower(int n, int m, const double* L, int ldl,
                              double* B, int ldb) {
  if (n <= 0 || m <= 0) return;
  const char side = 'R', uplo = 'L', trans = 'N', diag = 'N';
  const double one = 1.0;
  F77_CALL(dtrmm)(&side, &uplo, &trans, &diag, &n, &m, &one, L, &ldl, B, &ldb
                  FCONE FCONE FCONE FCONE);
}

void kf_blas_syr2k_lower(int n, int k, const double* A, int lda,
                         const double* B, int ldb, double* C, int ldc) {
  if (n <= 0 || k <= 0) return;
  const char uplo = 'L', trans = 'N';
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsyr2k)(&uplo, &trans, &n, &k, &one, A, &lda, B, &ldb, &zero, C,
                   &ldc FCONE FCONE);
}

void kf_blas_syrk_lower_acc(int n, int k, const double* A, int lda,
                            double* C, int ldc) {
  if (n <= 0 || k <= 0) return;
  const char uplo = 'L', trans = 'N';
  const double one = 1.0;
  F77_CALL(dsyrk)(&uplo, &trans, &n, &k, &one, A, &lda, &one, C, &ldc
                  FCONE FCONE);
}
