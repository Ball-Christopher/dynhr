// pskf_lapack.h -- LAPACK wrappers shared with pskf_filter.cpp.
#ifndef DYNHR_PSKF_LAPACK_H
#define DYNHR_PSKF_LAPACK_H

// eigen(symmetric = TRUE) of the n x n column-major A (lower triangle read):
// values in decreasing order, vectors n x n with matching columns.
bool pskf_sym_eigen(int n, const double* A, double* values, double* vectors);

// chol(C, pivot = TRUE) followed by chol2inv(): the inverse of C (upper
// triangle read), n x n column-major, in the ORIGINAL ordering. False when
// the pivoted factorisation has rank < n.
bool pskf_pivchol_inverse(int n, const double* C, double* out);

#endif
