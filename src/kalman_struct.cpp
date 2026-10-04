// kalman_struct.cpp -- structured standard Kalman filter for dynhr
//
// Structured standard Kalman filter (complete data, zero a0, stationary
// initial covariance -- the case kalman_standard_loop_cpp() in kalman_ss.cpp
// handles).
//
// The state-space matrices a solved DSGE model produces are dense as stored
// but numerically sparse: the static-elimination / QZ pipeline leaves entries
// at round-off (1e-16 .. 1e-14 of the matrix scale) wherever the structural
// coefficient is zero, and the shock blocks RR / DD / Sigma_e carry exact
// zeros. This kernel works on a THRESHOLDED view of TT, ZZ, RR and DD (an
// entry is dropped when |x| <= zero_tol * max|matrix|; a matrix holding a
// non-finite entry keeps every non-zero one) and uses it in three ways:
//   * ZZ, TT, DD are compressed sparse column, so Z s, P Z', Z (P Z'),
//     T (P Z'), K Z and K DD cost their non-zeros instead of a dense product;
//   * states whose TT AND ZZ columns are both empty never feed the next
//     state or the next innovation, so the Joseph term A P A' (A = TT - K ZZ)
//     is restricted to the remaining "active" states;
//   * Sigma_e diagonal (the usual case): B Sigma_e B' (B = RR - K DD) is one
//     rank-k update over the shocks with a positive variance.
// The covariance update A P A' uses P = L + L' (L = lower triangle of P with
// halved diagonal, exact) so that A P A' = M A' + A M' with M = A L: one
// triangular multiply and one rank-2k update, 1.5 n^3 flops instead of 2 n^3,
// and the result is exactly symmetric, so no symmetrising pass is needed.
//
// The recursion is the Joseph form of the dense kernel, operation for
// operation, over the thresholded matrices: it differs from the dense kernel
// by round-off (re-associated sums) plus the effect of the dropped entries,
// which sit at or below the accuracy of the solved matrices themselves.
// zero_tol = 0 drops only exact zeros. The dense kernel stays as the
// reference (R option dynhr.kf_zero_tol < 0 selects it).
//
// Both kernels here also run the simple covariance update (upd = 1,
// P' = T P T' + R Sigma R' - K F K', the algebraic equivalent of the Joseph
// form for the optimal gain) with the risk indicators of kf_guard.h; the
// general-loop kernel at the end of this file handles missing data, a given
// initial state / covariance, shock_scale and me_extra.

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <chrono>
#include "kf_blas.h"
#include "kf_guard.h"
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// Defined in kalman_ss.cpp (Dynare's badly-conditioned-F rule).
bool kf_general_F_singular(const arma::mat& Rc, const arma::mat& Ft,
                           const arma::mat& Fi, double kalman_tol);

namespace {

struct KfCsc {
  std::vector<arma::uword> ptr;   // size nc + 1
  std::vector<arma::uword> idx;   // row index of each stored entry
  std::vector<double>      val;
  arma::uword nr = 0, nc = 0;
};

// Drop threshold for one matrix: zero_tol * max|x|, or 0 (exact zeros only)
// when any entry is non-finite -- a NaN / Inf must reach the likelihood
// checks, never be thresholded away.
double kf_struct_thr(const arma::mat& M, double zero_tol) {
  double mx = 0.0;
  for (arma::uword k = 0; k < M.n_elem; ++k) {
    const double a = std::fabs(M[k]);
    if (!std::isfinite(a)) return 0.0;
    if (a > mx) mx = a;
  }
  return zero_tol * mx;
}

// `!(|a| <= thr)` keeps NaN; with thr == 0 it drops exactly the zeros.
KfCsc kf_struct_csc(const arma::mat& M, double thr) {
  KfCsc S;
  S.nr = M.n_rows; S.nc = M.n_cols;
  S.ptr.assign(S.nc + 1, 0);
  for (arma::uword j = 0; j < S.nc; ++j) {
    for (arma::uword i = 0; i < S.nr; ++i) {
      const double a = M(i, j);
      if (!(std::fabs(a) <= thr)) {
        S.idx.push_back(i);
        S.val.push_back(a);
      }
    }
    S.ptr[j + 1] = S.idx.size();
  }
  return S;
}

// CSC of the transpose (= compressed sparse ROW of the original).
KfCsc kf_struct_transpose(const KfCsc& S) {
  KfCsc R;
  R.nr = S.nc; R.nc = S.nr;
  R.ptr.assign(R.nc + 1, 0);
  for (arma::uword e = 0; e < S.idx.size(); ++e) ++R.ptr[S.idx[e] + 1];
  for (arma::uword j = 0; j < R.nc; ++j) R.ptr[j + 1] += R.ptr[j];
  R.idx.resize(S.idx.size());
  R.val.resize(S.val.size());
  std::vector<arma::uword> next(R.ptr.begin(), R.ptr.end() - 1);
  for (arma::uword j = 0; j < S.nc; ++j)
    for (arma::uword e = S.ptr[j]; e < S.ptr[j + 1]; ++e) {
      const arma::uword pos = next[S.idx[e]]++;
      R.idx[pos] = j;
      R.val[pos] = S.val[e];
    }
  return R;
}

arma::mat kf_struct_dense(const KfCsc& S) {
  arma::mat M(S.nr, S.nc, arma::fill::zeros);
  for (arma::uword j = 0; j < S.nc; ++j)
    for (arma::uword e = S.ptr[j]; e < S.ptr[j + 1]; ++e)
      M(S.idx[e], j) = S.val[e];
  return M;
}

// out += S * x  (S sparse, x dense vector)
inline void kf_struct_gaxpy(const KfCsc& S, const double* x, double* out) {
  for (arma::uword j = 0; j < S.nc; ++j) {
    const double xj = x[j];
    for (arma::uword e = S.ptr[j]; e < S.ptr[j + 1]; ++e)
      out[S.idx[e]] += S.val[e] * xj;
  }
}

inline bool kf_struct_col_nonempty(const KfCsc& S, arma::uword j) {
  return S.ptr[j + 1] > S.ptr[j];
}

// out_i = (Z B Z')_ii for every row i of Z, from the LOWER triangle of the
// symmetric B; `Zrow` is the CSC of Z' (column i = row i of Z).
void kf_struct_zbz_diag(const KfCsc& Zrow, const arma::mat& Blow,
                        arma::vec& out) {
  out.set_size(Zrow.nc);
  for (arma::uword i = 0; i < Zrow.nc; ++i) {
    double acc = 0.0;
    for (arma::uword e1 = Zrow.ptr[i]; e1 < Zrow.ptr[i + 1]; ++e1) {
      const arma::uword r1 = Zrow.idx[e1];
      for (arma::uword e2 = Zrow.ptr[i]; e2 < Zrow.ptr[i + 1]; ++e2) {
        const arma::uword r2 = Zrow.idx[e2];
        acc += Zrow.val[e1] * Zrow.val[e2] *
               (r1 >= r2 ? Blow(r1, r2) : Blow(r2, r1));
      }
    }
    out[i] = acc;
  }
}

}  // namespace

// [[Rcpp::export]]
List kalman_standard_struct_loop_cpp(const arma::mat& Y_minus_d,
                                     const arma::mat& ZZ,
                                     const arma::mat& TT,
                                     const arma::mat& RR,
                                     const arma::mat& DD,
                                     const arma::mat& HH_full,
                                     const arma::mat& Sigma_e,
                                     const arma::mat& SS,
                                     arma::mat P,
                                     double ll_const,
                                     double ss_tol,
                                     double ll_min,
                                     bool return_filtered,
                                     const arma::vec& me_diag_vec,
                                     double kalman_tol,
                                     double zero_tol,
                                     int upd = 0,
                                     double guard_piv = -1.0,
                                     double guard_r2 = -1.0,
                                     double guard_ret = -1.0,
                              double guard_amp = -1.0) {
  using arma::uword;
  const uword n = TT.n_rows;
  const uword p = ZZ.n_rows;
  const uword q = RR.n_cols;
  const uword n_T = Y_minus_d.n_cols;
  if (!(zero_tol >= 0.0) || !std::isfinite(zero_tol))
    Rcpp::stop("kalman_standard_struct_loop_cpp: zero_tol must be a finite number >= 0");
  if (TT.n_cols != n || ZZ.n_cols != n || RR.n_rows != n || DD.n_rows != p ||
      DD.n_cols != q || Sigma_e.n_rows != q || Sigma_e.n_cols != q ||
      HH_full.n_rows != p || HH_full.n_cols != p || SS.n_rows != n ||
      SS.n_cols != p || P.n_rows != n || P.n_cols != n ||
      Y_minus_d.n_rows != p)
    Rcpp::stop("kalman_standard_struct_loop_cpp: inconsistent matrix dimensions");

  // ---- once per call: thresholded sparse views and the active sets --------
  const KfCsc Tsp = kf_struct_csc(TT, kf_struct_thr(TT, zero_tol));
  const KfCsc Zsp = kf_struct_csc(ZZ, kf_struct_thr(ZZ, zero_tol));
  const KfCsc Dsp = kf_struct_csc(DD, kf_struct_thr(DD, zero_tol));
  const KfCsc Rsp = kf_struct_csc(RR, kf_struct_thr(RR, zero_tol));
  const KfCsc Zrow = kf_struct_transpose(Zsp);       // column i = row i of ZZ
  const arma::mat TTc = kf_struct_dense(Tsp);
  const arma::mat RRc = kf_struct_dense(Rsp);

  // Simple covariance update P' = T P T' + R Sigma R' - K F K' (upd = 1): the
  // algebraic equivalent of the Joseph recursion below for the optimal gain,
  // with F carrying the measurement-error diagonal. T P T' runs on the sparse
  // T (nnz * n flops twice), K F K' = (K/2) N' + N (K/2)' with N = K F is one
  // rank-2p update.
  const bool simple = (upd == 1);
  KfGuard guard(guard_piv, guard_r2, guard_ret, guard_amp);
  const KfCsc Trow = kf_struct_transpose(Tsp);       // column i = row i of TT
  arma::mat QQc, M1, Bs, Xm;
  if (simple) {
    QQc = (RRc * Sigma_e) * RRc.t();
    M1.set_size(n, n); Bs.set_size(n, n); Xm.set_size(n, n);
  }

  // States that neither feed the transition nor are observed: the matching
  // column of A = TT - K ZZ is identically zero.
  std::vector<uword> J;
  for (uword j = 0; j < n; ++j)
    if (kf_struct_col_nonempty(Tsp, j) || kf_struct_col_nonempty(Zsp, j))
      J.push_back(j);
  const uword m = J.size();

  // Sigma_e diagonal with non-negative entries -> rank-k form over the
  // shocks with a positive variance; anything else takes the general product.
  bool sig_diag = true;
  for (uword c = 0; c < q && sig_diag; ++c)
    for (uword r = 0; r < q; ++r)
      if (r != c && Sigma_e(r, c) != 0.0) { sig_diag = false; break; }
  if (sig_diag)
    for (uword k = 0; k < q; ++k)
      if (!(Sigma_e(k, k) >= 0.0)) { sig_diag = false; break; }
  std::vector<uword> Q;
  for (uword k = 0; k < q; ++k) {
    const bool act = kf_struct_col_nonempty(Rsp, k) ||
                     kf_struct_col_nonempty(Dsp, k);
    if (!act) continue;
    if (sig_diag && Sigma_e(k, k) == 0.0) continue;
    Q.push_back(k);
  }
  const uword qa = Q.size();
  arma::vec sqrtS(qa, arma::fill::zeros);
  arma::mat Sq;
  if (sig_diag) {
    for (uword kk = 0; kk < qa; ++kk)
      sqrtS[kk] = std::sqrt(Sigma_e(Q[kk], Q[kk]));
  } else {
    Sq.set_size(qa, qa);
    for (uword b = 0; b < qa; ++b)
      for (uword a = 0; a < qa; ++a) Sq(a, b) = Sigma_e(Q[a], Q[b]);
  }
  const bool has_me = me_diag_vec.n_elem > 0 &&
                      arma::any(arma::abs(me_diag_vec) > 0.0);
  const bool me_nonneg = !has_me || arma::all(me_diag_vec >= 0.0);

  arma::vec s(n, arma::fill::zeros);
  double loglik = 0.0;
  bool ok = true;
  bool ss_reached = false;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n, n_T);
  arma::mat K_ss, F_inv_ss;
  double ll_ss_const = 0.0;

  arma::vec v, zs, Fiv, Kv, s_n, tmp(n);
  arma::mat PZ(n, p), Ft, Rc, Fi, TPZ(n, p), K, Cm(n, n);
  arma::mat Aa(n, m), Ma(n, m), Lm(m, m), Bm(n, qa), Ks;

  for (uword t = 0; t < n_T; ++t) {
    zs.zeros(p);
    kf_struct_gaxpy(Zsp, s.memptr(), zs.memptr());
    v = Y_minus_d.col(t) - zs;

    if (!ss_reached) {
      // PZ = P ZZ'  (column i = sum over the non-zeros of row i of ZZ)
      PZ.zeros();
      for (uword i = 0; i < p; ++i) {
        double* pz = PZ.colptr(i);
        for (uword e = Zrow.ptr[i]; e < Zrow.ptr[i + 1]; ++e) {
          const double w = Zrow.val[e];
          const double* pc = P.colptr(Zrow.idx[e]);
          for (uword r = 0; r < n; ++r) pz[r] += w * pc[r];
        }
      }
      // Ft = ZZ PZ + HH, symmetrised
      Ft.set_size(p, p);
      for (uword c = 0; c < p; ++c) {
        const double* pzc = PZ.colptr(c);
        for (uword i = 0; i < p; ++i) {
          double acc = 0.0;
          for (uword e = Zrow.ptr[i]; e < Zrow.ptr[i + 1]; ++e)
            acc += Zrow.val[e] * pzc[Zrow.idx[e]];
          Ft(i, c) = acc;
        }
      }

      Ft += HH_full;
      Ft = 0.5 * (Ft + Ft.t());
      if (!arma::chol(Rc, Ft)) { ok = false; break; }
      if (!arma::inv_sympd(Fi, Ft)) { ok = false; break; }
      if (kf_general_F_singular(Rc, Ft, Fi, kalman_tol)) { ok = false; break; }
      const double ldf = 2.0 * arma::sum(arma::log(Rc.diag()));
      Fiv = Fi * v;
      const double ll = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;

      // K = (TT PZ + SS) Fi
      TPZ.zeros();
      for (uword c = 0; c < p; ++c) {
        const double* pzc = PZ.colptr(c);
        double* tc = TPZ.colptr(c);
        for (uword k = 0; k < n; ++k) {
          const double b = pzc[k];
          for (uword e = Tsp.ptr[k]; e < Tsp.ptr[k + 1]; ++e)
            tc[Tsp.idx[e]] += Tsp.val[e] * b;
        }
      }
      s_n.zeros(n);
      kf_struct_gaxpy(Tsp, s.memptr(), s_n.memptr());

      TPZ += SS;
      K = TPZ * Fi;
      Kv = K * v;
      s_n += Kv;

      if (simple) {
        guard.note_F(Rc, Ft, Fi);
        if (guard.tripped) break;
        // M1 = T P (sparse T times dense P)
        M1.zeros();
        for (uword j = 0; j < n; ++j) {
          const double* pj = P.colptr(j);
          double* mj = M1.colptr(j);
          for (uword k = 0; k < n; ++k) {
            const double w = pj[k];
            for (uword e = Tsp.ptr[k]; e < Tsp.ptr[k + 1]; ++e)
              mj[Tsp.idx[e]] += Tsp.val[e] * w;
          }
        }
        // Bs (lower) = M1 T' + R Sigma R'
        for (uword i = 0; i < n; ++i) {
          double* bi = Bs.colptr(i);
          for (uword r = i; r < n; ++r) bi[r] = QQc(r, i);
          for (uword e = Trow.ptr[i]; e < Trow.ptr[i + 1]; ++e) {
            const double w = Trow.val[e];
            const double* mk = M1.colptr(Trow.idx[e]);
            for (uword r = i; r < n; ++r) bi[r] += w * mk[r];
          }
        }
        kf_struct_zbz_diag(Zrow, Bs, guard.pending);
        // Xm (lower) = 0.5 (K N' + N K'),  N = K F = TPZ
        Ks = K;
        Ks *= 0.5;
        kf_blas_syr2k_lower(static_cast<int>(n), static_cast<int>(p),
                            Ks.memptr(), static_cast<int>(n),
                            TPZ.memptr(), static_cast<int>(n),
                            Xm.memptr(), static_cast<int>(n));
        for (uword c = 0; c < n; ++c) {
          for (uword r = c; r < n; ++r) Cm(r, c) = Bs(r, c) - Xm(r, c);
        }
        for (uword c = 0; c < n; ++c)
          for (uword r = c + 1; r < n; ++r) Cm(c, r) = Cm(r, c);
        // diagonal-only guard read (lower triangle of Bs holds the base)
        guard.note_P(Bs, Cm);
        if (guard.tripped) break;
      } else {
      // A = TT - K ZZ over the active columns; B = RR - K DD over the
      // active shocks (scaled by the shock std when Sigma_e is diagonal).
      for (uword jj = 0; jj < m; ++jj) {
        const uword j = J[jj];
        tmp.zeros();
        double* tp = tmp.memptr();
        for (uword e = Zsp.ptr[j]; e < Zsp.ptr[j + 1]; ++e) {
          const double w = Zsp.val[e];
          const double* kc = K.colptr(Zsp.idx[e]);
          for (uword r = 0; r < n; ++r) tp[r] += w * kc[r];
        }
        double* ac = Aa.colptr(jj);
        const double* tc = TTc.colptr(j);
        for (uword r = 0; r < n; ++r) ac[r] = tc[r] - tp[r];
      }
      for (uword kk = 0; kk < qa; ++kk) {
        const uword k = Q[kk];
        tmp.zeros();
        double* tp = tmp.memptr();
        for (uword e = Dsp.ptr[k]; e < Dsp.ptr[k + 1]; ++e) {
          const double w = Dsp.val[e];
          const double* kc = K.colptr(Dsp.idx[e]);
          for (uword r = 0; r < n; ++r) tp[r] += w * kc[r];
        }
        double* bc = Bm.colptr(kk);
        const double* rc = RRc.colptr(k);
        const double sc = sig_diag ? sqrtS[kk] : 1.0;
        for (uword r = 0; r < n; ++r) bc[r] = sc * (rc[r] - tp[r]);
      }

      // P' = A P A' + B Sigma B' (+ K diag(me) K'), lower triangle, mirrored.
      Cm.zeros();
      if (m > 0) {
        for (uword b = 0; b < m; ++b) {
          double* lc = Lm.colptr(b);
          for (uword a = 0; a < m; ++a)
            lc[a] = (a > b) ? P(J[a], J[b])
                            : (a == b ? 0.5 * P(J[a], J[a]) : 0.0);
        }
        Ma = Aa;
        kf_blas_trmm_right_lower(static_cast<int>(n), static_cast<int>(m),
                                 Lm.memptr(), static_cast<int>(m),
                                 Ma.memptr(), static_cast<int>(n));
        kf_blas_syr2k_lower(static_cast<int>(n), static_cast<int>(m),
                            Ma.memptr(), static_cast<int>(n),
                            Aa.memptr(), static_cast<int>(n),
                            Cm.memptr(), static_cast<int>(n));
      }
      if (qa > 0) {
        if (sig_diag) {
          kf_blas_syrk_lower_acc(static_cast<int>(n), static_cast<int>(qa),
                                 Bm.memptr(), static_cast<int>(n),
                                 Cm.memptr(), static_cast<int>(n));
        } else {
          const arma::mat G = (Bm * Sq) * Bm.t();
          for (uword c = 0; c < n; ++c)
            for (uword r = c; r < n; ++r) Cm(r, c) += G(r, c);
        }
      }
      if (has_me) {
        if (me_nonneg) {
          Ks = K;
          for (uword i = 0; i < p; ++i)
            Ks.col(i) *= std::sqrt(me_diag_vec[i]);
          kf_blas_syrk_lower_acc(static_cast<int>(n), static_cast<int>(p),
                                 Ks.memptr(), static_cast<int>(n),
                                 Cm.memptr(), static_cast<int>(n));
        } else {
          const arma::mat G = (K * arma::diagmat(me_diag_vec)) * K.t();
          for (uword c = 0; c < n; ++c)
            for (uword r = c; r < n; ++r) Cm(r, c) += G(r, c);
        }
      }
      for (uword c = 0; c < n; ++c)
        for (uword r = c + 1; r < n; ++r) Cm(c, r) = Cm(r, c);

      }
      s = s_n;
      if (t > 0 && arma::abs(Cm - P).max() < ss_tol * arma::abs(Cm).max()) {
        ss_reached  = true;
        K_ss        = K;
        F_inv_ss    = Fi;
        ll_ss_const = ll_const - 0.5 * ldf;
      }
      P.swap(Cm);
    } else {
      Fiv = F_inv_ss * v;
      const double ll = ll_ss_const - 0.5 * arma::dot(v, Fiv);
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      s_n.zeros(n);
      kf_struct_gaxpy(Tsp, s.memptr(), s_n.memptr());
      Kv = K_ss * v;
      s_n += Kv;
      s = s_n;
    }

    if (return_filtered) filtered.col(t) = s;
  }

  return List::create(_["loglik"]     = loglik,
                      _["s"]          = s,
                      _["P"]          = P,
                      _["filtered"]   = return_filtered ? Rcpp::wrap(filtered)
                                                        : R_NilValue,
                      _["ok"]         = ok,
                      _["ss_reached"] = ss_reached,
                      _["n_active_states"] = static_cast<int>(m),
                      _["n_active_shocks"] = static_cast<int>(qa),
                      _["nnz_TT"]     = static_cast<int>(Tsp.idx.size()),
                      _["nnz_ZZ"]     = static_cast<int>(Zsp.idx.size()),
                      _["sigma_diag"] = sig_diag,
                      _["guard_tripped"] = guard.tripped,
                      _["ind_piv"]    = guard.piv,
                      _["ind_r2"]     = guard.r2,
                      _["ind_ret"]    = guard.ret,
                      _["ind_amp"]    = guard.amp,
                      _["ind_bad"]    = guard.bad);
}


// How fast is this process's BLAS relative to a plain compiled loop? Returns
// t(dgemm, 40 x 40 x 40) / t(the same product as a column-oriented axpy
// loop), each best-of-N, measured once per process and cached. An optimised
// BLAS (Accelerate, OpenBLAS, MKL) returns ~0.05-0.15; a reference BLAS
// (R's bundled libRblas / Rblas.dll, "Matrix products: default") returns
// ~0.4-1. kalman_filter() uses it to pick between the dense-BLAS kernel and
// the structured kernel, whose sparse loops beat only a slow dgemm.
// [[Rcpp::export]]
double kalman_blas_probe_cpp() {
  static const double ratio = []() {
    const arma::uword n = 40;
    arma::mat A(n, n), B(n, n), C(n, n), D(n, n);
    for (arma::uword j = 0; j < n; ++j)
      for (arma::uword i = 0; i < n; ++i) {
        A(i, j) = 1.0 / static_cast<double>(1 + i + 2 * j);
        B(i, j) = 1.0 / static_cast<double>(3 + 2 * i + j);
      }
    using clock = std::chrono::steady_clock;
    const int reps = 40, trials = 9;
    double t_blas = 1e300, t_loop = 1e300, sink = 0.0;
    C = A * B;                                   // warm the BLAS / caches
    for (int tr = 0; tr < trials; ++tr) {
      auto t0 = clock::now();
      for (int r = 0; r < reps; ++r) { C = A * B; sink += C(r % n, 0); }
      auto t1 = clock::now();
      for (int r = 0; r < reps; ++r) {
        D.zeros();
        for (arma::uword j = 0; j < n; ++j) {
          double* dj = D.colptr(j);
          for (arma::uword l = 0; l < n; ++l) {
            const double b = B(l, j);
            const double* al = A.colptr(l);
            for (arma::uword i = 0; i < n; ++i) dj[i] += al[i] * b;
          }
        }
        sink += D(r % n, 0);
      }
      auto t2 = clock::now();
      t_blas = std::min(t_blas, std::chrono::duration<double>(t1 - t0).count());
      t_loop = std::min(t_loop, std::chrono::duration<double>(t2 - t1).count());
    }
    if (!std::isfinite(sink)) return 1.0;
    return t_blas / t_loop;
  }();
  return ratio;
}


// ---------------------------------------------------------------------------
// Structured general standard filter: the per-step loop of
// kalman_standard_general_loop_cpp() (kalman_ss.cpp) -- missing observations,
// a given initial state / covariance / start period, per-period shock scales
// and extra measurement-error variances -- on the thresholded sparse views of
// the structured kernel above. A period with observed rows O updates with the
// observed rows of ZZ / DD only; the sparse view of those rows (observed-row
// CSC of ZZ and DD, the active states and shocks) is built once per DISTINCT
// missingness pattern in the call and reused by every period that shows it.
// Branch for branch the steady-state lock, the lock release on a missing
// period, the ll floor and the ll constant correction are the dense general
// loop's; the update is the Joseph form (upd = 0) or the simple form (upd = 1)
// exactly as in the complete-data kernel.
// ---------------------------------------------------------------------------

namespace {

struct KfPat {
  std::vector<arma::uword> O;      // observed rows
  KfCsc Zo, Zr, Do;                // Zo: p_o x n (CSC), Zr = Zo', Do: p_o x q
  arma::mat Dd;                    // dense Do
  std::vector<arma::uword> J;      // states feeding T or observed in O
  std::vector<arma::uword> Q;      // shocks entering R or D_o (variance > 0 if diagonal)
  arma::mat HHo, SSo;              // baseline blocks restricted to O
};

KfCsc kf_struct_restrict_rows(const KfCsc& S, const std::vector<arma::uword>& rowmap,
                              arma::uword n_new) {
  const arma::uword none = static_cast<arma::uword>(-1);
  KfCsc R;
  R.nr = n_new; R.nc = S.nc;
  R.ptr.assign(S.nc + 1, 0);
  for (arma::uword j = 0; j < S.nc; ++j) {
    for (arma::uword e = S.ptr[j]; e < S.ptr[j + 1]; ++e) {
      const arma::uword r = rowmap[S.idx[e]];
      if (r == none) continue;
      R.idx.push_back(r);
      R.val.push_back(S.val[e]);
    }
    R.ptr[j + 1] = R.idx.size();
  }
  return R;
}

// Cm (LOWER) = T P T' + QQ, with T sparse (Tsp, Trow = its transpose).
void kf_struct_tpt_plus(const KfCsc& Tsp, const KfCsc& Trow, const arma::mat& P,
                        const arma::mat& QQ, arma::mat& M1, arma::mat& Cm) {
  using arma::uword;
  const uword n = P.n_rows;
  M1.zeros();
  for (uword j = 0; j < n; ++j) {
    const double* pj = P.colptr(j);
    double* mj = M1.colptr(j);
    for (uword k = 0; k < n; ++k) {
      const double w = pj[k];
      for (uword e = Tsp.ptr[k]; e < Tsp.ptr[k + 1]; ++e)
        mj[Tsp.idx[e]] += Tsp.val[e] * w;
    }
  }
  for (uword i = 0; i < n; ++i) {
    double* bi = Cm.colptr(i);
    for (uword r = i; r < n; ++r) bi[r] = QQ(r, i);
    for (uword e = Trow.ptr[i]; e < Trow.ptr[i + 1]; ++e) {
      const double w = Trow.val[e];
      const double* mk = M1.colptr(Trow.idx[e]);
      for (uword r = i; r < n; ++r) bi[r] += w * mk[r];
    }
  }
}

}  // namespace

// [[Rcpp::export]]
List kalman_standard_general_struct_loop_cpp(const arma::mat& Y_minus_d,
                                             const arma::mat& ZZ,
                                             const arma::mat& TT,
                                             const arma::mat& RR,
                                             const arma::mat& DD,
                                             const arma::mat& HH,
                                             const arma::mat& QQ,
                                             const arma::mat& Sigma_e,
                                             const arma::mat& SS,
                                             arma::vec s,
                                             arma::mat P,
                                             int t_start,
                                             double init_loglik,
                                             double ll_const,
                                             double ss_tol,
                                             double ll_min,
                                             bool return_filtered,
                                             const arma::vec& me_vec,
                                             const arma::mat& me_extra,
                                             const arma::mat& shock_scale,
                                             double kalman_tol,
                                             double zero_tol,
                                             int upd = 0,
                                             double guard_piv = -1.0,
                                             double guard_r2 = -1.0,
                                             double guard_ret = -1.0,
                              double guard_amp = -1.0) {
  using arma::uword;
  const uword n = TT.n_rows;
  const uword p = ZZ.n_rows;
  const uword q = RR.n_cols;
  const uword n_T = Y_minus_d.n_cols;
  const bool has_mx = me_extra.n_elem > 0;
  const bool has_sc = shock_scale.n_elem > 0;
  if (!(zero_tol >= 0.0) || !std::isfinite(zero_tol))
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: zero_tol must be a finite number >= 0");
  if (TT.n_cols != n || ZZ.n_cols != n || RR.n_rows != n || DD.n_rows != p ||
      DD.n_cols != q || Sigma_e.n_rows != q || Sigma_e.n_cols != q ||
      HH.n_rows != p || HH.n_cols != p || SS.n_rows != n || SS.n_cols != p ||
      QQ.n_rows != n || QQ.n_cols != n || Y_minus_d.n_rows != p)
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: inconsistent matrix dimensions");
  if (has_mx && (me_extra.n_rows != p || me_extra.n_cols != n_T))
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: me_extra must be n_obs x T");
  if (has_sc && (shock_scale.n_rows != q || shock_scale.n_cols != n_T))
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: shock_scale must be n_exo x T");
  if (me_vec.n_elem != p)
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: me_vec must have length n_obs");
  if (s.n_elem != n || P.n_rows != n || P.n_cols != n)
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: s / P do not match TT");
  if (t_start < 1)
    Rcpp::stop("kalman_standard_general_struct_loop_cpp: t_start must be >= 1");

  const KfCsc Tsp = kf_struct_csc(TT, kf_struct_thr(TT, zero_tol));
  const KfCsc Zsp = kf_struct_csc(ZZ, kf_struct_thr(ZZ, zero_tol));
  const KfCsc Dsp = kf_struct_csc(DD, kf_struct_thr(DD, zero_tol));
  const KfCsc Rsp = kf_struct_csc(RR, kf_struct_thr(RR, zero_tol));
  const KfCsc Trow = kf_struct_transpose(Tsp);
  const KfCsc Zrow_all = kf_struct_transpose(Zsp);
  const arma::mat RRc = kf_struct_dense(Rsp);
  const arma::mat TTc = kf_struct_dense(Tsp);

  const bool simple = (upd == 1);
  KfGuard guard(guard_piv, guard_r2, guard_ret, guard_amp);

  bool sig_diag = true;
  for (uword c = 0; c < q && sig_diag; ++c)
    for (uword r = 0; r < q; ++r)
      if (r != c && Sigma_e(r, c) != 0.0) { sig_diag = false; break; }
  if (sig_diag)
    for (uword k = 0; k < q; ++k)
      if (!(Sigma_e(k, k) >= 0.0)) { sig_diag = false; break; }
  // A per-period scale keeps a diagonal Sigma diagonal; a non-finite scale
  // does not -- the general product then carries it to the likelihood checks.
  if (has_sc && sig_diag && !shock_scale.is_finite()) sig_diag = false;

  const arma::mat QQc = (RRc * Sigma_e) * RRc.t();
  const double log2pi = std::log(2.0 * M_PI);

  std::map<std::string, KfPat> cache;
  const uword none = static_cast<uword>(-1);
  auto make_pat = [&](const std::string& key, const std::vector<uword>& O) {
    KfPat pt;
    pt.O = O;
    const uword po = O.size();
    std::vector<uword> rowmap(p, none);
    for (uword a = 0; a < po; ++a) rowmap[O[a]] = a;
    pt.Zo = kf_struct_restrict_rows(Zsp, rowmap, po);
    pt.Do = kf_struct_restrict_rows(Dsp, rowmap, po);
    pt.Zr = kf_struct_transpose(pt.Zo);
    pt.Dd = kf_struct_dense(pt.Do);
    for (uword j = 0; j < n; ++j)
      if (kf_struct_col_nonempty(Tsp, j) || kf_struct_col_nonempty(pt.Zo, j))
        pt.J.push_back(j);
    for (uword k = 0; k < q; ++k) {
      if (!(kf_struct_col_nonempty(Rsp, k) || kf_struct_col_nonempty(pt.Do, k)))
        continue;
      if (sig_diag && Sigma_e(k, k) == 0.0) continue;
      pt.Q.push_back(k);
    }
    pt.HHo.set_size(po, po);
    pt.SSo.set_size(n, po);
    for (uword b = 0; b < po; ++b) {
      for (uword a = 0; a < po; ++a) pt.HHo(a, b) = HH(O[a], O[b]);
      for (uword r = 0; r < n; ++r) pt.SSo(r, b) = SS(r, O[b]);
    }
    return cache.emplace(key, std::move(pt)).first;
  };

  double loglik = init_loglik;
  bool ok = true, ss_reached = false;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n, n_T);
  arma::mat K_ss, F_inv_ss;
  double ll_ss_const = 0.0;

  arma::vec y, v, zs, v_o, s_n, Kv, me_t, tmp(n), sqrtS;
  arma::mat Se_t, QQ_t, PZ, Ft, Rc, Fi, N, K, Cm(n, n), Bs(n, n), Xm(n, n);
  arma::mat M1(n, n), Aa, Ma, Lm, Bm, Ks, Sq, Dsc;
  std::vector<uword> O;
  std::string key(p, '0');

  for (uword t = static_cast<uword>(t_start - 1); t < n_T; ++t) {
    y = Y_minus_d.col(t);
    O.clear();
    for (uword i = 0; i < p; ++i) {
      const bool obs = !std::isnan(y[i]);
      key[i] = obs ? '1' : '0';
      if (obs) O.push_back(i);
    }
    const uword n_ok = O.size();
    zs.zeros(p);
    kf_struct_gaxpy(Zsp, s.memptr(), zs.memptr());
    v = y - zs;

    if (has_sc) {
      const arma::vec sc = shock_scale.col(t);
      Se_t = Sigma_e % (sc * sc.t());
      QQ_t = (RRc * Se_t) * RRc.t();
    }
    const arma::mat& Se = has_sc ? Se_t : Sigma_e;
    const arma::mat& QQe = has_sc ? QQ_t : QQc;

    if (n_ok == 0) {
      ss_reached = false;
      s_n.zeros(n);
      kf_struct_gaxpy(Tsp, s.memptr(), s_n.memptr());
      s = s_n;
      kf_struct_tpt_plus(Tsp, Trow, P, QQe, M1, Cm);
      guard.pending.reset();
      for (uword c = 0; c < n; ++c)
        for (uword r = c + 1; r < n; ++r) Cm(c, r) = Cm(r, c);
      P = Cm;
      if (return_filtered) filtered.col(t) = s;
      continue;
    }

    const bool complete = (n_ok == p);
    if (!complete) ss_reached = false;

    if (ss_reached) {                          // complete, locked
      const arma::vec Fiv = F_inv_ss * v;
      const double ll = ll_ss_const - 0.5 * arma::dot(v, Fiv);
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      s_n.zeros(n);
      kf_struct_gaxpy(Tsp, s.memptr(), s_n.memptr());
      Kv = K_ss * v;
      s_n += Kv;
      s = s_n;
      if (return_filtered) filtered.col(t) = s;
      continue;
    }

    auto it = cache.find(key);
    if (it == cache.end()) it = make_pat(key, O);
    const KfPat& pt = it->second;
    const uword po = n_ok;

    me_t.set_size(po);
    for (uword a = 0; a < po; ++a)
      me_t[a] = me_vec[pt.O[a]] + (has_mx ? me_extra(pt.O[a], t) : 0.0);

    // PZ = P Z_o'
    PZ.set_size(n, po);
    PZ.zeros();
    for (uword a = 0; a < po; ++a) {
      double* pz = PZ.colptr(a);
      for (uword e = pt.Zr.ptr[a]; e < pt.Zr.ptr[a + 1]; ++e) {
        const double w = pt.Zr.val[e];
        const double* pc = P.colptr(pt.Zr.idx[e]);
        for (uword r = 0; r < n; ++r) pz[r] += w * pc[r];
      }
    }
    // Ft = Z_o PZ + D_o Se D_o' + diag(me_t)
    Ft.set_size(po, po);
    for (uword c = 0; c < po; ++c) {
      const double* pzc = PZ.colptr(c);
      for (uword a = 0; a < po; ++a) {
        double acc = 0.0;
        for (uword e = pt.Zr.ptr[a]; e < pt.Zr.ptr[a + 1]; ++e)
          acc += pt.Zr.val[e] * pzc[pt.Zr.idx[e]];
        Ft(a, c) = acc;
      }
    }
    if (has_sc) {
      Dsc = (pt.Dd * Se) * pt.Dd.t();
      Ft += Dsc;
    } else {
      Ft += pt.HHo;
    }
    for (uword a = 0; a < po; ++a) Ft(a, a) += me_t[a];
    Ft = 0.5 * (Ft + Ft.t());
    if (!arma::chol(Rc, Ft)) { ok = false; break; }
    if (!arma::inv_sympd(Fi, Ft)) { ok = false; break; }
    if (kf_general_F_singular(Rc, Ft, Fi, kalman_tol)) { ok = false; break; }
    const double ldf = 2.0 * arma::sum(arma::log(Rc.diag()));
    v_o.set_size(po);
    for (uword a = 0; a < po; ++a) v_o[a] = v[pt.O[a]];
    const double quad = arma::dot(v_o, Fi * v_o);
    if (complete) {
      const double ll = ll_const - 0.5 * (ldf + quad);
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
    } else {
      // ll_const bakes in the full n_obs: correct it up by the number of
      // missing components (no ll_min floor on a partial period, as in the
      // dense loop).
      loglik += ll_const + 0.5 * static_cast<double>(p - n_ok) * log2pi -
        0.5 * (ldf + quad);
    }
    if (simple) {
      guard.note_F(Rc, Ft, Fi, pt.O.data());
      if (guard.tripped) break;
    }

    // N = T PZ + S_o
    N.set_size(n, po);
    N.zeros();
    for (uword c = 0; c < po; ++c) {
      const double* pzc = PZ.colptr(c);
      double* nc = N.colptr(c);
      for (uword k = 0; k < n; ++k) {
        const double b = pzc[k];
        for (uword e = Tsp.ptr[k]; e < Tsp.ptr[k + 1]; ++e)
          nc[Tsp.idx[e]] += Tsp.val[e] * b;
      }
    }
    if (has_sc) N += (RRc * Se) * pt.Dd.t();
    else        N += pt.SSo;
    K = N * Fi;
    s_n.zeros(n);
    kf_struct_gaxpy(Tsp, s.memptr(), s_n.memptr());
    Kv = K * v_o;
    s_n += Kv;

    if (simple) {
      kf_struct_tpt_plus(Tsp, Trow, P, QQe, M1, Bs);
      kf_struct_zbz_diag(Zrow_all, Bs, guard.pending);
      Ks = K;
      Ks *= 0.5;
      kf_blas_syr2k_lower(static_cast<int>(n), static_cast<int>(po),
                          Ks.memptr(), static_cast<int>(n),
                          N.memptr(), static_cast<int>(n),
                          Xm.memptr(), static_cast<int>(n));
      for (uword c = 0; c < n; ++c)
        for (uword r = c; r < n; ++r) Cm(r, c) = Bs(r, c) - Xm(r, c);
      for (uword c = 0; c < n; ++c)
        for (uword r = c + 1; r < n; ++r) Cm(c, r) = Cm(r, c);
      guard.note_P(Bs, Cm);
      if (guard.tripped) break;
    } else {
      const uword m = pt.J.size();
      const uword qa = pt.Q.size();
      Cm.zeros();
      Aa.set_size(n, m);
      for (uword jj = 0; jj < m; ++jj) {
        const uword j = pt.J[jj];
        tmp.zeros();
        double* tp = tmp.memptr();
        for (uword e = pt.Zo.ptr[j]; e < pt.Zo.ptr[j + 1]; ++e) {
          const double w = pt.Zo.val[e];
          const double* kc = K.colptr(pt.Zo.idx[e]);
          for (uword r = 0; r < n; ++r) tp[r] += w * kc[r];
        }
        double* ac = Aa.colptr(jj);
        const double* tc = TTc.colptr(j);
        for (uword r = 0; r < n; ++r) ac[r] = tc[r] - tp[r];
      }
      Bm.set_size(n, qa);
      if (sig_diag) sqrtS.set_size(qa);
      for (uword kk = 0; kk < qa; ++kk) {
        const uword k = pt.Q[kk];
        tmp.zeros();
        double* tp = tmp.memptr();
        for (uword e = pt.Do.ptr[k]; e < pt.Do.ptr[k + 1]; ++e) {
          const double w = pt.Do.val[e];
          const double* kc = K.colptr(pt.Do.idx[e]);
          for (uword r = 0; r < n; ++r) tp[r] += w * kc[r];
        }
        double sc1 = 1.0;
        if (sig_diag) sc1 = sqrtS[kk] = std::sqrt(Se(k, k));
        double* bc = Bm.colptr(kk);
        const double* rc = RRc.colptr(k);
        for (uword r = 0; r < n; ++r) bc[r] = sc1 * (rc[r] - tp[r]);
      }
      if (m > 0) {
        Lm.set_size(m, m);
        for (uword b = 0; b < m; ++b) {
          double* lc = Lm.colptr(b);
          for (uword a = 0; a < m; ++a)
            lc[a] = (a > b) ? P(pt.J[a], pt.J[b])
                            : (a == b ? 0.5 * P(pt.J[a], pt.J[a]) : 0.0);
        }
        Ma = Aa;
        kf_blas_trmm_right_lower(static_cast<int>(n), static_cast<int>(m),
                                 Lm.memptr(), static_cast<int>(m),
                                 Ma.memptr(), static_cast<int>(n));
        kf_blas_syr2k_lower(static_cast<int>(n), static_cast<int>(m),
                            Ma.memptr(), static_cast<int>(n),
                            Aa.memptr(), static_cast<int>(n),
                            Cm.memptr(), static_cast<int>(n));
      }
      if (qa > 0) {
        if (sig_diag) {
          kf_blas_syrk_lower_acc(static_cast<int>(n), static_cast<int>(qa),
                                 Bm.memptr(), static_cast<int>(n),
                                 Cm.memptr(), static_cast<int>(n));
        } else {
          Sq.set_size(qa, qa);
          for (uword b = 0; b < qa; ++b)
            for (uword a = 0; a < qa; ++a) Sq(a, b) = Se(pt.Q[a], pt.Q[b]);
          const arma::mat G = (Bm * Sq) * Bm.t();
          for (uword c = 0; c < n; ++c)
            for (uword r = c; r < n; ++r) Cm(r, c) += G(r, c);
        }
      }
      if (arma::any(me_t != 0.0)) {
        if (arma::all(me_t >= 0.0)) {
          Ks = K;
          for (uword a = 0; a < po; ++a) Ks.col(a) *= std::sqrt(me_t[a]);
          kf_blas_syrk_lower_acc(static_cast<int>(n), static_cast<int>(po),
                                 Ks.memptr(), static_cast<int>(n),
                                 Cm.memptr(), static_cast<int>(n));
        } else {
          const arma::mat G = (K * arma::diagmat(me_t)) * K.t();
          for (uword c = 0; c < n; ++c)
            for (uword r = c; r < n; ++r) Cm(r, c) += G(r, c);
        }
      }
      for (uword c = 0; c < n; ++c)
        for (uword r = c + 1; r < n; ++r) Cm(c, r) = Cm(r, c);
    }

    s = s_n;
    // The lock engages on a complete, time-invariant period only (never in a
    // run carrying shock_scale or me_extra), as in the dense loop.
    if (complete && !has_mx && !has_sc && t > 0 &&
        arma::abs(Cm - P).max() < ss_tol * arma::abs(Cm).max()) {
      ss_reached  = true;
      K_ss        = K;
      F_inv_ss    = Fi;
      ll_ss_const = ll_const - 0.5 * ldf;
    }
    P = Cm;
    if (return_filtered) filtered.col(t) = s;
  }

  return List::create(_["loglik"]     = loglik,
                      _["s"]          = s,
                      _["P"]          = P,
                      _["filtered"]   = return_filtered ? Rcpp::wrap(filtered)
                                                        : R_NilValue,
                      _["ok"]         = ok,
                      _["ss_reached"] = ss_reached,
                      _["n_patterns"] = static_cast<int>(cache.size()),
                      _["guard_tripped"] = guard.tripped,
                      _["ind_piv"]    = guard.piv,
                      _["ind_r2"]     = guard.r2,
                      _["ind_ret"]    = guard.ret,
                      _["ind_amp"]    = guard.amp,
                      _["ind_bad"]    = guard.bad);
}
