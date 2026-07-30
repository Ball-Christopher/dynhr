// hank_ar_score.cpp -- compact contraction weights for the exact-AR
// sequence-space score (R/hank-kalman-ar-grad.R).
//
// The exact-AR likelihood assembles the stacked covariance by a FIXED index
// gather from the compact autocovariance array G (n_obs x n_obs x T_data,
// ~19k entries): S[r, c] = G[Sidx[r, c]]. The gather is linear, so for any
// parameter phi, dS/dphi = gather(dG/dphi) and
//
//     d loglik/d phi = -0.5 [ tr(Sinv dS) - v' dS v ]
//                    = -0.5 sum_k w[k] (dG/dphi)[k],
//     w[k] = sum_{(r,c): Sidx[r,c] = k} ( Sinv[r,c] - v[r] v[c] ).
//
// So ONE pass over the N^2 stacked entries collapses both contractions onto
// G-space, after which each parameter costs a ~19k dot product instead of a
// fresh N x N gather plus an N^2 trace and an N^2 matvec. This is the
// "compact accumulation" the paper brief (notes/dynhr_brief_seqspace_
// gradients.md 4b, briefs/20-seqspace-gradient-scope.md 5) identified as the
// only place the real speedup lives -- and which is 3x SLOWER when written in
// R with rowsum(), because grouping millions of entries into thousands of
// groups is interpreter-bound there and trivial here.
//
// The R fallback (.hank_ar_score_weights(use_cpp = FALSE)) computes the same
// quantity via rowsum(); parity is asserted in test-hank-kalman-ar-grad.R.

#include <Rcpp.h>

// Block-Toeplitz pass of the autocovariance-kernel adjoint
// (R/hank-kalman-ar-grad.R::.hank_ar_slab_adjoint):
//
//   dL/dPsi[t,] = sum_{d <= t-1} B_d Psi[t-d,] + sum_{d <= q-t} B_d' Psi[t+d,]
//
// Both folds of every B_d, in one pass. `Bvec` is n_obs^2 x q with column d+1
// holding vec(B_d); `Psi` is q x n_obs. The q x q x n_obs^2 shape makes this
// the dominant cost of the adjoint (~q^2 n_obs^2 FMAs), and in R it is q BLAS
// calls on skinny matrices -- measured 1.4x SLOWER than the finite-difference
// rho path it is meant to replace. Here it is one loop nest.
//
// B_d decays geometrically in d (the rho^|k-d| / rho^(k+d) weight kernels), so
// the all-zero skip below is a real saving once the tail underflows.
// [[Rcpp::export]]
Rcpp::NumericMatrix hank_ar_slab_adjoint_psi_cpp(const Rcpp::NumericMatrix& Bvec,
                                                 const Rcpp::NumericMatrix& Psi) {
  const int q = Psi.nrow();
  const int n_obs = Psi.ncol();
  const int nn = n_obs * n_obs;
  if (Bvec.nrow() != nn)
    Rcpp::stop("hank_ar_slab_adjoint_psi_cpp: nrow(Bvec) must be ncol(Psi)^2.");
  if (Bvec.ncol() != q)
    Rcpp::stop("hank_ar_slab_adjoint_psi_cpp: ncol(Bvec) must be nrow(Psi).");

  Rcpp::NumericMatrix GPsi(q, n_obs);
  double* g = &GPsi[0];
  const double* p = &Psi[0];
  const double* B = &Bvec[0];

  for (int d = 0; d < q; ++d) {
    const double* Bd = B + (std::size_t)d * nn;
    bool nonzero = false;
    for (int e = 0; e < nn && !nonzero; ++e) nonzero = (Bd[e] != 0.0);
    if (!nonzero) continue;                     // rho^d underflowed
    for (int j = 0; j < q - d; ++j) {
      const int t = j + d;
      for (int b = 0; b < n_obs; ++b) {
        const double p_jb = p[j + (std::size_t)b * q];
        double acc_b = 0.0;                     // (B_d' Psi[t,])[b]
        for (int a = 0; a < n_obs; ++a) {
          const double Bab = Bd[a + b * n_obs];
          // fold 1: GPsi[t,] += B_d Psi[j,]
          g[t + (std::size_t)a * q] += Bab * p_jb;
          // fold 2: GPsi[j,] += B_d' Psi[t,]
          acc_b += Bab * p[t + (std::size_t)a * q];
        }
        g[j + (std::size_t)b * q] += acc_b;
      }
    }
  }
  return GPsi;
}


// [[Rcpp::export]]
Rcpp::NumericVector hank_ar_score_weights_cpp(const Rcpp::IntegerMatrix& Sidx,
                                              const Rcpp::NumericMatrix& Sinv,
                                              const Rcpp::NumericVector& v,
                                              int n_g) {
  const int n = Sidx.nrow();
  if (Sidx.ncol() != n)
    Rcpp::stop("hank_ar_score_weights_cpp: Sidx must be square.");
  if (Sinv.nrow() != n || Sinv.ncol() != n)
    Rcpp::stop("hank_ar_score_weights_cpp: Sinv must match dim(Sidx).");
  if (v.size() != n)
    Rcpp::stop("hank_ar_score_weights_cpp: length(v) must match nrow(Sidx).");
  if (n_g <= 0)
    Rcpp::stop("hank_ar_score_weights_cpp: n_g must be positive.");

  Rcpp::NumericVector w(n_g);
  double* wp = &w[0];
  const int* ip = &Sidx[0];
  const double* sp = &Sinv[0];
  const double* vp = &v[0];

  for (int c = 0; c < n; ++c) {
    const double vc = vp[c];
    const int off = c * n;
    for (int r = 0; r < n; ++r) {
      const int k = ip[off + r] - 1;               // Sidx is 1-based
      if (k < 0 || k >= n_g)
        Rcpp::stop("hank_ar_score_weights_cpp: Sidx entry out of range "
                   "[1, n_g] -- gather index and G dimensions disagree.");
      wp[k] += sp[off + r] - vp[r] * vc;
    }
  }
  return w;
}
