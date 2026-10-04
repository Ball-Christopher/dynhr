// kf_blas.h -- three BLAS level-3 calls the structured Kalman kernel needs
// (triangular multiply, symmetric rank-k and rank-2k updates). They live in
// their own translation unit because R_ext/BLAS.h and Armadillo's BLAS
// prototypes must not meet in one file (see ordered_qz.cpp for the same
// constraint on LAPACK). Plain pointers, column-major, no Armadillo types.

#ifndef DYNHR_KF_BLAS_H
#define DYNHR_KF_BLAS_H

// B (n x m) := B * L, L (m x m) lower triangular, non-unit diagonal.
void kf_blas_trmm_right_lower(int n, int m, const double* L, int ldl,
                              double* B, int ldb);

// C (n x n, LOWER triangle only) := A * B' + B * A', A and B n x k.
// The strictly upper triangle of C is not referenced.
void kf_blas_syr2k_lower(int n, int k, const double* A, int lda,
                         const double* B, int ldb, double* C, int ldc);

// C (n x n, LOWER triangle only) += A * A', A n x k.
void kf_blas_syrk_lower_acc(int n, int k, const double* A, int lda,
                            double* C, int ldc);

#endif
