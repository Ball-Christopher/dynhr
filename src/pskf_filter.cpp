// pskf_filter.cpp -- the PSKF per-period recursion in C++.
//
// pskf_filter_cpp() is the body of .pskf_filter() (R/pskf-likelihood.R),
// operation for operation; with store_path = TRUE it also returns the
// per-period arrays the smoother's backward pass consumes (same names and
// shapes as the R loop's lists). The recursion: Kalman prediction, CSN
// stacking, the S_pred inverse (pivot-tested Cholesky, eigen pseudoinverse
// otherwise: .csn_sym_inv / .csn_sym_pinv), pruning with the mean
// compensation (dim_red4_r, .csn_offset_g), missing and linearly dependent
// observables (.pskf_independent_obs, .pskf_dependent_obs_consistent), the
// update and the two CDF terms through the C++ dispatch of mvn_cdf.cpp.
//
// The R code stays the reference. Wherever it would leave code this file
// mirrors -- a non-finite matrix, a CDF the C++ evaluators decline
// (bivariate quadrature below p = 1e-3, Miwa, a numerically singular lattice
// call), a failed factorisation -- the recursion throws Fallback and the
// function returns status 1: the caller then runs the R recursion from the
// start, so results never depend on which path answered. The only CDF
// evaluator reproduced beyond mvn_cdf.cpp is Mendell-Elston (blocks larger
// than miwa_qmax), copied from logcdf_ME_r including its update order.
//
// keep_override (a list of per-period kept-row vectors, or NULL) freezes the
// pruning selection: dim_red4 uses the given rows instead of selecting, and
// still computes the mean compensation (.pskf_filter, keep_override).
// store_keep = TRUE (without store_path) adds keep_path to the result.
//
// Returns list(status, ll[, <path lists>]): status 0 = ll is the
// log-likelihood (possibly -Inf; the path lists are then absent), status 1 =
// fall back to R.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include "mvn_cdf.h"
#include "pskf_lapack.h"
// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;
using arma::vec;
using arma::uword;

namespace {

struct Fallback {};

// (A B) C in R's left-to-right order (Armadillo may pick the other order, at
// a different rounding)
inline mat mul3(const mat& A, const mat& B, const mat& C) {
  const mat AB = A * B;
  return AB * C;
}

const double kRtol = 1.4901161193847656e-08;  // sqrt(.Machine$double.eps)

// ---- log Phi_q: logcdf_ME_r(x, S, miwa_qmax, check = TRUE) ------------------
double logcdf_full(const std::vector<double>& x, const std::vector<double>& S,
                   int q, double miwa_qmax) {
  if (q == 0) return 0.0;
  for (int i = 0; i < q; ++i) if (!std::isfinite(x[i])) throw Fallback();
  for (size_t i = 0; i < S.size(); ++i) if (!std::isfinite(S[i])) throw Fallback();
  if (q == 1) {
    const double b = x[0] / std::sqrt(S[0]);
    if (!std::isfinite(b)) throw Fallback();
    return R::pnorm(b, 0.0, 1.0, 1, 1);
  }
  double v;
  if (logcdf_dispatch(x, S, q, miwa_qmax, v)) return v;

  // the R path after the dispatch declined: snap, connected blocks, ...
  std::vector<double> sdv(q), Cm(q * q), xs(q);
  for (int i = 0; i < q; ++i)
    sdv[i] = std::sqrt(std::max(S[i + i * q], DBL_EPSILON));
  bool all_adj = true;
  for (int j = 0; j < q; ++j)
    for (int i = 0; i < q; ++i) {
      double c = S[i + j * q] / (sdv[i] * sdv[j]);
      if (i != j && std::fabs(c) < 1e-12) c = 0.0;
      Cm[i + j * q] = c;
      if (c == 0.0) all_adj = false;
    }
  for (int i = 0; i < q; ++i) xs[i] = x[i] / sdv[i];
  if (!all_adj) {
    std::vector<int> comp(q, 0);
    int n_comp = 0;
    for (int s0 = 0; s0 < q; ++s0) {
      if (comp[s0] > 0) continue;
      ++n_comp;
      std::vector<int> stack(1, s0);
      while (!stack.empty()) {
        const int v0 = stack.front();
        stack.erase(stack.begin());
        if (comp[v0] > 0) continue;
        comp[v0] = n_comp;
        for (int w = 0; w < q; ++w)
          if (Cm[v0 + w * q] != 0.0 && comp[w] == 0) stack.push_back(w);
      }
    }
    if (n_comp > 1) {
      double tot = 0.0;
      for (int cc = 1; cc <= n_comp; ++cc) {
        std::vector<int> ii;
        for (int i = 0; i < q; ++i) if (comp[i] == cc) ii.push_back(i);
        const int m = ii.size();
        std::vector<double> xb(m), Sb(m * m);
        for (int a = 0; a < m; ++a) {
          xb[a] = xs[ii[a]];
          for (int b = 0; b < m; ++b) Sb[a + b * m] = Cm[ii[a] + ii[b] * q];
        }
        tot = tot + logcdf_full(xb, Sb, m, miwa_qmax);
      }
      return tot;
    }
  }
  // one connected block the C++ evaluators declined: only q > miwa_qmax
  // (Mendell-Elston) is reproduced here
  if (q >= 4) {
    // weak couplings in a block that needs the lattice or Mendell-Elston:
    // zero them and split again (logcdf_ME_r, same rule and threshold)
    bool weak = false;
    std::vector<double> C2(Cm);
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < q; ++i)
        if (i != j && Cm[i + j * q] != 0.0 && std::fabs(Cm[i + j * q]) < 1e-8) {
          C2[i + j * q] = 0.0;
          weak = true;
        }
    if (weak) return logcdf_full(xs, C2, q, miwa_qmax);
  }
  if (q == 2 || (double) q <= miwa_qmax) throw Fallback();
  std::vector<double> b(xs), CS(Cm);
  double log_p = 0.0;
  for (int j = 0; j < q - 1; ++j) {
    const double sjj = CS[j + j * q];
    const double sj = std::sqrt(std::max(sjj, DBL_EPSILON));
    const double bj = b[j] / sj;
    const double lPj = R::pnorm(bj, 0.0, 1.0, 1, 1);
    log_p = log_p + lPj;
    const double lphi_j = R::dnorm(bj, 0.0, 1.0, 1);
    const double lambda = std::exp(lphi_j - lPj);
    const double delta_factor = lambda * (bj + lambda);
    for (int k = j + 1; k < q; ++k) {
      // E[X_k | X_j <= b_j] = -S_jk / sqrt(S_jj) * lambda, so the bound for
      // X_k rises by S_jk / sqrt(S_jj) * lambda.
      b[k] = b[k] + CS[j + k * q] / sj * lambda;
      // upper triangle once, then mirrored (l from j + 1 decremented every
      // off-diagonal twice)
      for (int l = k; l < q; ++l) {
        CS[k + l * q] = CS[k + l * q] - CS[j + k * q] * CS[j + l * q] / sjj *
                        delta_factor;
        CS[l + k * q] = CS[k + l * q];
      }
    }
  }
  const double bq = b[q - 1] / std::sqrt(std::max(CS[q * q - 1], DBL_EPSILON));
  return log_p + R::pnorm(bq, 0.0, 1.0, 1, 1);
}

double logcdf_vec(const vec& x, const mat& S, double miwa_qmax) {
  const int q = x.n_elem;
  std::vector<double> xv(x.begin(), x.end()), Sv(S.begin(), S.end());
  return logcdf_full(xv, Sv, q, miwa_qmax);
}

// ---- .csn_offset_g -----------------------------------------------------------
vec csn_offset_g(const mat& Gamma, const vec& nu, const mat& Delta,
                 const mat& Sigma, double miwa_qmax) {
  const uword q = Gamma.n_rows;
  if (q == 0) return vec();
  mat V = Delta + mul3(Gamma, Sigma, Gamma.t());
  V = (V + V.t()) / 2.0;
  const double logZ = logcdf_vec(-nu, V, miwa_qmax);
  vec g(q, arma::fill::zeros);
  for (uword j = 0; j < q; ++j) {
    const double Vjj = V(j, j);
    if (Vjj <= 0) continue;
    const double lphi = R::dnorm(-nu[j], 0.0, std::sqrt(Vjj), 1);
    double lcond;
    if (q == 1) {
      lcond = 0.0;
    } else {
      vec mcond(q - 1);
      mat Scond(q - 1, q - 1);
      uword a = 0;
      for (uword i = 0; i < q; ++i) {
        if (i == j) continue;
        mcond[a] = (-nu[i]) - V(i, j) / Vjj * (-nu[j]);
        uword b = 0;
        for (uword k = 0; k < q; ++k) {
          if (k == j) continue;
          Scond(a, b) = V(i, k) - (V(i, j) * V(j, k)) / Vjj;
          ++b;
        }
        ++a;
      }
      Scond = (Scond + Scond.t()) / 2.0;
      lcond = logcdf_vec(mcond, Scond, miwa_qmax);
    }
    g[j] = std::exp(lphi + lcond - logZ);
  }
  if (!g.is_finite()) return vec(q, arma::fill::zeros);
  return g;
}

// ---- .csn_sym_pinv / .csn_sym_inv --------------------------------------------
mat csn_sym_pinv(const mat& S) {
  const uword n = S.n_rows;
  mat P(n, n, arma::fill::zeros);
  if (n == 0) return P;
  vec dS = S.diag();
  double mx = 0.0;
  for (uword i = 0; i < n; ++i) {
    if (!std::isfinite(dS[i])) dS[i] = 0.0;
    mx = std::max(mx, dS[i]);
  }
  const double floor_abs = n * DBL_EPSILON * mx;
  std::vector<uword> keep;
  for (uword i = 0; i < n; ++i) if (dS[i] > floor_abs) keep.push_back(i);
  if (keep.empty()) return P;
  const uword m = keep.size();
  vec ds(m);
  for (uword a = 0; a < m; ++a) ds[a] = std::sqrt(std::max(dS[keep[a]], 0.0));
  mat C(m, m);
  for (uword a = 0; a < m; ++a)
    for (uword b = 0; b < m; ++b)
      C(a, b) = S(keep[a], keep[b]) / (ds[a] * ds[b]);
  C = 0.5 * (C + C.t());
  if (!C.is_finite()) throw Fallback();
  vec ev(m);
  mat V(m, m);
  if (!pskf_sym_eigen((int) m, C.memptr(), ev.memptr(), V.memptr()))
    throw Fallback();
  const double ev_max = ev[0];   // decreasing order, as R's eigen()
  if (!std::isfinite(ev_max) || ev_max <= 0) return P;
  std::vector<uword> pos;
  for (uword i = 0; i < m; ++i) if (ev[i] > kRtol * ev_max) pos.push_back(i);
  if (pos.empty()) return P;
  mat Vp(m, pos.size());
  vec inv_ev(pos.size());
  for (uword c = 0; c < pos.size(); ++c) {
    Vp.col(c) = V.col(pos[c]);
    inv_ev[c] = 1.0 / ev[pos[c]];
  }
  mat Vt = Vp.t();
  Vt.each_col() %= inv_ev;
  const mat Q = Vp * Vt;
  for (uword a = 0; a < m; ++a)
    for (uword b = 0; b < m; ++b)
      P(keep[a], keep[b]) = Q(a, b) / (ds[a] * ds[b]);
  return 0.5 * (P + P.t());
}

mat csn_sym_inv(const mat& S, bool fast) {
  const uword n = S.n_rows;
  if (!fast || n == 0) return csn_sym_pinv(S);
  const vec dS = S.diag();
  if (!dS.is_finite() || dS.min() <= n * DBL_EPSILON * dS.max())
    return csn_sym_pinv(S);
  const vec d = arma::sqrt(dS);
  mat C = S / (d * d.t());
  C = 0.5 * (C + C.t());
  if (!C.is_finite()) return csn_sym_pinv(S);
  mat Ci(n, n);
  if (!pskf_pivchol_inverse((int) n, C.memptr(), Ci.memptr()))
    return csn_sym_pinv(S);
  const double max_cond = std::min(1e5, 0.25 / kRtol);
  if (!Ci.is_finite() || n * std::sqrt(arma::accu(Ci % Ci)) >= max_cond)
    return csn_sym_pinv(S);
  mat P = Ci / (d * d.t());
  return 0.5 * (P + P.t());
}

// ---- dim_red4_r ----------------------------------------------------------------
struct Pruned {
  mat Gamma;
  vec nu;
  mat Delta;
  vec mu_shift;
  std::vector<int> keep;  // 1-based rows of the pre-prune stack that survived
  vec lambda;             // latent-space cut compensation on that stack
};

// The natural selection of dim_red4_r (.dim_red4_select): cut_tol threshold,
// collinearity guard, max_q rank cap. 0-based, increasing.
std::vector<uword> dim_red4_select(const mat& Gamma, const mat& Delta,
                                   const mat& Sigma, double cut_tol,
                                   double max_q) {
  const uword q = Gamma.n_rows, p = Gamma.n_cols;
  const mat GS = Gamma * Sigma;
  const mat cov_full = Delta + GS * Gamma.t();
  vec d_skew(q), d_state(p);
  for (uword i = 0; i < q; ++i)
    d_skew[i] = std::sqrt(std::max(cov_full(i, i), DBL_EPSILON));
  for (uword j = 0; j < p; ++j)
    d_state[j] = std::sqrt(std::max(Sigma(j, j), DBL_EPSILON));
  vec max_corr(q, arma::fill::zeros);
  for (uword i = 0; i < q; ++i) {
    double m = -arma::datum::inf;
    for (uword j = 0; j < p; ++j)
      m = std::max(m, std::fabs(GS(i, j)) / (d_skew[i] * d_state[j]));
    max_corr[i] = m;
  }
  std::vector<uword> keep;
  for (uword i = 0; i < q; ++i) if (max_corr[i] >= cut_tol) keep.push_back(i);
  if (keep.size() > 1) {
    for (;;) {
      const uword m = keep.size();
      if (m <= 1) break;
      // Ck = |V_k| / outer(d_k, d_k) with a zero diagonal; the first (column-
      // major) position of its maximum, as which(Ck == max(Ck), arr.ind)[1, ]
      double best = -arma::datum::inf;
      uword br = 0, bc = 0;
      for (uword c = 0; c < m; ++c)
        for (uword r = 0; r < m; ++r) {
          const double v = (r == c) ? 0.0
              : std::fabs(cov_full(keep[r], keep[c])) /
                (d_skew[keep[r]] * d_skew[keep[c]]);
          if (v > best) { best = v; br = r; bc = c; }
        }
      if (best <= 0.9) break;
      const uword pr = keep[br], pc = keep[bc];
      const uword drop_row = (max_corr[pc] < max_corr[pr]) ? pc : pr;
      keep.erase(std::find(keep.begin(), keep.end(), drop_row));
    }
  }
  if (std::isfinite(max_q) && (double) keep.size() > max_q) {
    std::vector<uword> ord(keep);
    std::stable_sort(ord.begin(), ord.end(), [&](uword a, uword b) {
      return max_corr[a] > max_corr[b];
    });
    ord.resize((size_t) max_q);
    std::sort(ord.begin(), ord.end());
    keep = ord;
  }
  return keep;
}

// dim_red4_r. keep_ov (0-based, increasing; may be null): a frozen selection
// used in place of dim_red4_select; the mean compensation is computed as for
// a selected set.
Pruned dim_red4(const mat& Gamma, const vec& nu, const mat& Delta,
                const mat& Sigma, double cut_tol, double max_q,
                double offset_miwa_qmax, const std::vector<uword>* keep_ov) {
  const uword q = Gamma.n_rows, p = Gamma.n_cols;
  Pruned out;
  {
    const mat GS = Gamma * Sigma;
    const mat cov_full = Delta + GS * Gamma.t();
    if (!GS.is_finite() || !cov_full.is_finite()) throw Fallback();
  }
  const std::vector<uword> keep = keep_ov
      ? *keep_ov
      : dim_red4_select(Gamma, Delta, Sigma, cut_tol, max_q);
  if (keep.size() == q) {
    out.Gamma = Gamma;
    out.nu = nu;
    out.Delta = Delta;
    out.mu_shift = vec(p, arma::fill::zeros);
    out.keep.resize(q);
    for (uword i = 0; i < q; ++i) out.keep[i] = (int) i + 1;
    out.lambda = vec(q, arma::fill::zeros);
    return out;
  }
  const mat SGt = Sigma * Gamma.t();
  const vec g_before = csn_offset_g(Gamma, nu, Delta, Sigma, offset_miwa_qmax);
  if (keep.empty()) {
    out.Gamma = mat(0, p);
    out.nu = vec();
    out.Delta = mat(0, 0);
    out.mu_shift = SGt * g_before;
    out.lambda = g_before;
    return out;
  }
  const uword k = keep.size();
  mat G_k(k, p), D_k(k, k);
  vec n_k(k);
  for (uword a = 0; a < k; ++a) {
    G_k.row(a) = Gamma.row(keep[a]);
    n_k[a] = nu[keep[a]];
    for (uword b = 0; b < k; ++b) D_k(a, b) = Delta(keep[a], keep[b]);
  }
  vec lambda = g_before;
  const vec g_after = csn_offset_g(G_k, n_k, D_k, Sigma, offset_miwa_qmax);
  for (uword a = 0; a < k; ++a) lambda[keep[a]] = lambda[keep[a]] - g_after[a];
  out.Gamma = G_k;
  out.nu = n_k;
  out.Delta = D_k;
  out.mu_shift = SGt * lambda;
  out.lambda = lambda;
  for (uword a = 0; a < k; ++a) out.keep.push_back((int) keep[a] + 1);
  return out;
}

// ---- .pskf_independent_obs / .pskf_dependent_obs_consistent -------------------
struct IndepObs {
  std::vector<uword> keep;
  mat L;
};

IndepObs independent_obs(const mat& Omega) {
  const uword m = Omega.n_rows;
  vec d(m);
  for (uword i = 0; i < m; ++i) d[i] = std::sqrt(std::max(Omega(i, i), 0.0));
  IndepObs r;
  mat Lc(0, 0);
  for (uword j = 0; j < m; ++j) {
    if (!(d[j] > 0)) continue;
    const uword k = r.keep.size();
    vec w;
    double piv = 1.0;
    if (k > 0) {
      vec c_j(k);
      for (uword a = 0; a < k; ++a)
        c_j[a] = Omega(r.keep[a], j) / (d[r.keep[a]] * d[j]);
      w = arma::solve(arma::trimatl(Lc), c_j, arma::solve_opts::fast);
      piv = 1.0 - arma::accu(w % w);
    }
    if (piv > kRtol) {
      mat Ln(k + 1, k + 1, arma::fill::zeros);
      if (k > 0) {
        Ln.submat(0, 0, k - 1, k - 1) = Lc;
        for (uword a = 0; a < k; ++a) Ln(k, a) = w[a];
      }
      Ln(k, k) = std::sqrt(piv);
      Lc = Ln;
      r.keep.push_back(j);
    }
  }
  const uword k = r.keep.size();
  r.L.set_size(k, k);
  for (uword a = 0; a < k; ++a) r.L.row(a) = d[r.keep[a]] * Lc.row(a);
  return r;
}

bool dependent_consistent(const mat& Omega, const vec& v, const vec& v_scale,
                          const IndepObs& ind) {
  const uword m = v.n_elem;
  std::vector<bool> is_keep(m, false);
  for (uword a : ind.keep) is_keep[a] = true;
  std::vector<uword> dep;
  for (uword i = 0; i < m; ++i) if (!is_keep[i]) dep.push_back(i);
  if (dep.empty()) return true;
  const uword nd = dep.size(), nk = ind.keep.size();
  vec r(nd), r_scale(nd);
  if (nk == 0) {
    for (uword a = 0; a < nd; ++a) { r[a] = v[dep[a]]; r_scale[a] = v_scale[dep[a]]; }
  } else {
    mat Okd(nk, nd);
    for (uword a = 0; a < nk; ++a)
      for (uword b = 0; b < nd; ++b) Okd(a, b) = Omega(ind.keep[a], dep[b]);
    const mat X = arma::solve(arma::trimatl(ind.L), Okd, arma::solve_opts::fast);
    const mat Y = arma::solve(arma::trimatu(ind.L.t()), X, arma::solve_opts::fast);
    const mat A = Y.t();  // nd x nk
    vec vk(nk), sk(nk);
    for (uword a = 0; a < nk; ++a) { vk[a] = v[ind.keep[a]]; sk[a] = v_scale[ind.keep[a]]; }
    const vec Av = A * vk, Asv = arma::abs(A) * sk;
    for (uword a = 0; a < nd; ++a) {
      r[a] = v[dep[a]] - Av[a];
      r_scale[a] = v_scale[dep[a]] + Asv[a];
    }
  }
  for (uword a = 0; a < nd; ++a) {
    const double tol = kRtol * r_scale[a] +
                       std::sqrt(kRtol * std::max(Omega(dep[a], dep[a]), 0.0));
    if (!(std::fabs(r[a]) <= tol)) return false;
  }
  return true;
}

// R's solve(A, tol = 0): dgesv against the identity
bool lu_inverse(const mat& A, mat& out) {
  arma::blas_int n = A.n_rows, nrhs = n, info = 0;
  mat W = A;
  out.eye(n, n);
  std::vector<arma::blas_int> ipiv(n);
  arma::lapack::gesv(&n, &nrhs, W.memptr(), &n, ipiv.data(), out.memptr(), &n, &info);
  return info == 0;
}

}  // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::List pskf_filter_cpp(const arma::mat& Y, const arma::mat& TT,
                           const arma::mat& ZZ, const arma::vec& mu_eta,
                           const arma::mat& Sigma_eta,
                           const arma::mat& Gamma_eta, const arma::vec& nu_eta,
                           const arma::mat& Delta_eta, const arma::vec& mu_eps,
                           const arma::mat& Sigma_eps, const arma::mat& P0,
                           double cut_tol, double max_q,
                           double offset_miwa_qmax, bool fast_solve,
                           bool store_path, SEXP keep_override,
                           bool store_keep) {
  const uword n_obs = Y.n_rows, n_T = Y.n_cols, n = TT.n_rows;
  double ll = 0.0;
  // keep_override: NULL, or a list of n_T 1-based kept-row vectors (validated
  // by the R caller, .pskf_keep_override_ok)
  const bool frozen = !Rf_isNull(keep_override);
  std::vector<std::vector<uword>> keep_ov;
  if (frozen) {
    const Rcpp::List kl(keep_override);
    if ((uword) kl.size() != n_T) Rcpp::stop("keep_override: length != T");
    keep_ov.resize(n_T);
    for (uword t = 0; t < n_T; ++t) {
      const Rcpp::NumericVector kt(kl[t]);
      for (R_xlen_t a = 0; a < kt.size(); ++a)
        keep_ov[t].push_back((uword) kt[a] - 1);
    }
  }
  const bool keep_out = store_path || store_keep;
  // per-period path (store_path = TRUE), same names and shapes as the R loop
  Rcpp::List mu_pred_path(n_T), Sigma_pred_path(n_T), mu_filt_path(n_T),
      Sigma_filt_path(n_T), K_gauss_path(n_T), ZZ_path(n_T),
      Gamma_filt_path(n_T), nu_filt_path(n_T), Delta_filt_path(n_T),
      Gamma_pred_path(n_T), nu_pred_path(n_T), Delta_pred_path(n_T),
      keep_path(n_T), lambda_path(n_T);
  auto as_vec = [](const vec& v) {
    return Rcpp::NumericVector(v.begin(), v.end());
  };
  auto as_ivec = [](const std::vector<int>& v) {
    return Rcpp::IntegerVector(v.begin(), v.end());
  };
  try {
    vec mu_filt(n, arma::fill::zeros);
    mat Sigma_filt = P0;
    mat Gamma_filt(0, n);
    vec nu_filt;
    mat Delta_filt(0, 0);
    const mat tTT = TT.t();
    const mat I_n = arma::eye(n, n);
    const uword q_eta = Gamma_eta.n_rows;
    const double log2pi = std::log(2.0 * arma::datum::pi);

    for (uword t = 0; t < n_T; ++t) {
      // prediction
      vec mu_pred = TT * mu_filt + mu_eta;
      const mat TT_Sfilt = TT * Sigma_filt;
      const mat Sigma_pred = TT_Sfilt * tTT + Sigma_eta;
      if (!Sigma_pred.is_finite()) throw Fallback();
      const uword q_filt = Gamma_filt.n_rows, q_new = q_filt + q_eta;
      std::vector<int> keep_t;
      vec lambda_t;
      mat Gamma_pred(0, n);
      vec nu_pred;
      mat Delta_pred(0, 0);
      if (q_new > 0) {
        const mat S_pred_inv = csn_sym_inv(Sigma_pred, fast_solve);
        mat M_filt, M_eta, G_filt_block, G_eta_block, D11, D22;
        if (q_filt > 0) {
          M_filt = TT_Sfilt.t() * S_pred_inv;
          G_filt_block = Gamma_filt * M_filt;
        }
        if (q_eta > 0) {
          M_eta = Sigma_eta * S_pred_inv;
          G_eta_block = Gamma_eta * M_eta;
        }
        if (q_filt == 0) Gamma_pred = G_eta_block;
        else if (q_eta == 0) Gamma_pred = G_filt_block;
        else Gamma_pred = arma::join_cols(G_filt_block, G_eta_block);
        nu_pred = arma::join_cols(nu_filt, nu_eta);
        if (q_filt > 0) {
          const mat Schur_filt = M_filt * TT_Sfilt;
          D11 = Delta_filt + mul3(Gamma_filt, Sigma_filt - Schur_filt, Gamma_filt.t());
        }
        if (q_eta > 0) {
          const mat Schur_eta = Sigma_eta - M_eta * Sigma_eta;
          D22 = Delta_eta + mul3(Gamma_eta, Schur_eta, Gamma_eta.t());
        }
        if (q_filt == 0) {
          Delta_pred = D22;
        } else if (q_eta == 0) {
          Delta_pred = D11;
        } else {
          const mat D12 = -(G_filt_block * Sigma_eta) * Gamma_eta.t();
          Delta_pred = arma::join_cols(arma::join_rows(D11, D12),
                                       arma::join_rows(D12.t(), D22));
        }
        if (!Gamma_pred.is_finite() || !Delta_pred.is_finite() ||
            !nu_pred.is_finite())
          throw Fallback();
        // pruning with the first-moment compensation
        Pruned pr = dim_red4(Gamma_pred, nu_pred, Delta_pred, Sigma_pred,
                             cut_tol, max_q, offset_miwa_qmax,
                             frozen ? &keep_ov[t] : nullptr);
        Gamma_pred = pr.Gamma;
        nu_pred = pr.nu;
        Delta_pred = pr.Delta;
        mu_pred = mu_pred + pr.mu_shift;
        keep_t = pr.keep;
        lambda_t = pr.lambda;
      }
      const uword q_pred = Gamma_pred.n_rows;
      if (keep_out) keep_path[t] = as_ivec(keep_t);

      // observed rows
      std::vector<uword> obs;
      for (uword i = 0; i < n_obs; ++i)
        if (!std::isnan(Y(i, t))) obs.push_back(i);
      mat ZZ_t(obs.size(), n), Sigma_eps_t(obs.size(), obs.size());
      vec y_t(obs.size()), mu_eps_t(obs.size());
      for (uword a = 0; a < obs.size(); ++a) {
        ZZ_t.row(a) = ZZ.row(obs[a]);
        y_t[a] = Y(obs[a], t);
        mu_eps_t[a] = mu_eps[obs[a]];
        for (uword b = 0; b < obs.size(); ++b)
          Sigma_eps_t(a, b) = Sigma_eps(obs[a], obs[b]);
      }

      vec v_t;
      mat Omega;
      if (y_t.n_elem > 0) {
        const vec Zmu_t = ZZ_t * mu_pred;
        v_t = y_t - Zmu_t - mu_eps_t;
        Omega = mul3(ZZ_t, Sigma_pred, ZZ_t.t()) + Sigma_eps_t;
        if (!Omega.is_finite() || !v_t.is_finite())
          return Rcpp::List::create(Rcpp::_["status"] = 0,
                                    Rcpp::_["ll"] = R_NegInf);
        const IndepObs ind = independent_obs(Omega);
        if (ind.keep.size() < v_t.n_elem) {
          const vec v_scale = arma::abs(y_t) + arma::abs(Zmu_t) +
                              arma::abs(mu_eps_t);
          if (!dependent_consistent(Omega, v_t, v_scale, ind))
            return Rcpp::List::create(Rcpp::_["status"] = 0,
                                      Rcpp::_["ll"] = R_NegInf);
          const uword k = ind.keep.size();
          mat ZZk(k, n), Sek(k, k), Omk(k, k);
          vec yk(k), vk(k);
          for (uword a = 0; a < k; ++a) {
            ZZk.row(a) = ZZ_t.row(ind.keep[a]);
            yk[a] = y_t[ind.keep[a]];
            vk[a] = v_t[ind.keep[a]];
            for (uword b = 0; b < k; ++b) {
              Sek(a, b) = Sigma_eps_t(ind.keep[a], ind.keep[b]);
              Omk(a, b) = Omega(ind.keep[a], ind.keep[b]);
            }
          }
          ZZ_t = ZZk; Sigma_eps_t = Sek; y_t = yk; v_t = vk; Omega = Omk;
        }
      }

      if (y_t.n_elem == 0) {
        mu_filt = mu_pred;
        Sigma_filt = Sigma_pred;
        Gamma_filt = Gamma_pred;
        nu_filt = nu_pred;
        Delta_filt = Delta_pred;
        if (store_path) {
          mu_pred_path[t] = as_vec(mu_pred);
          Sigma_pred_path[t] = Sigma_pred;
          mu_filt_path[t] = as_vec(mu_filt);
          Sigma_filt_path[t] = Sigma_filt;
          K_gauss_path[t] = mat(n, 0);
          ZZ_path[t] = mat(0, n);
          Gamma_filt_path[t] = Gamma_filt;
          nu_filt_path[t] = as_vec(nu_filt);
          Delta_filt_path[t] = Delta_filt;
          Gamma_pred_path[t] = Gamma_pred;
          nu_pred_path[t] = as_vec(nu_pred);
          Delta_pred_path[t] = Delta_pred;
          keep_path[t] = as_ivec(keep_t);
          lambda_path[t] = as_vec(lambda_t);
        }
        continue;
      }
      const uword n_obs_t = y_t.n_elem;

      // update
      mat Oc;
      if (!arma::chol(Oc, Omega)) throw Fallback();
      double log_det = 0.0;
      for (uword i = 0; i < n_obs_t; ++i) log_det += std::log(Oc(i, i));
      log_det *= 2.0;
      const vec fw = arma::solve(arma::trimatl(Oc.t()), v_t,
                                 arma::solve_opts::fast);
      const vec Oiv = arma::solve(arma::trimatu(Oc), fw, arma::solve_opts::fast);
      const double ll_gauss = -0.5 * (n_obs_t * log2pi + log_det +
                                      arma::accu(v_t % Oiv));
      mat Om_inv;
      if (!lu_inverse(Omega, Om_inv)) throw Fallback();
      const mat K = mul3(Sigma_pred, ZZ_t.t(), Om_inv);
      vec nu_upd = nu_pred;
      if (q_pred > 0) nu_upd = nu_pred - (Gamma_pred * K) * v_t;
      const vec mu_upd = mu_pred + K * v_t;
      const mat I_KZ = I_n - K * ZZ_t;
      mat Sigma_upd = mul3(I_KZ, Sigma_pred, I_KZ.t()) + mul3(K, Sigma_eps_t, K.t());
      Sigma_upd = 0.5 * (Sigma_upd + Sigma_upd.t());

      double ll_skew = 0.0;
      if (q_pred > 0) {
        mat D_bot = Delta_pred + mul3(Gamma_pred, Sigma_pred, Gamma_pred.t());
        D_bot = 0.5 * (D_bot + D_bot.t());
        mat D_top = Delta_pred + mul3(Gamma_pred, Sigma_upd, Gamma_pred.t());
        D_top = 0.5 * (D_top + D_top.t());
        const double ll_cdf_bot = logcdf_vec(-nu_pred, D_bot, 5.0);
        const double ll_cdf_top = logcdf_vec(-nu_upd, D_top, 5.0);
        ll_skew = ll_cdf_top - ll_cdf_bot;
      }
      ll = ll + ll_gauss + ll_skew;

      mu_filt = mu_upd;
      Sigma_filt = Sigma_upd;
      Gamma_filt = Gamma_pred;
      nu_filt = nu_upd;
      Delta_filt = Delta_pred;
      if (store_path) {
        mu_pred_path[t] = as_vec(mu_pred);
        Sigma_pred_path[t] = Sigma_pred;
        mu_filt_path[t] = as_vec(mu_filt);
        Sigma_filt_path[t] = Sigma_filt;
        K_gauss_path[t] = K;
        ZZ_path[t] = ZZ_t;
        Gamma_filt_path[t] = Gamma_filt;
        nu_filt_path[t] = as_vec(nu_filt);
        Delta_filt_path[t] = Delta_filt;
        Gamma_pred_path[t] = Gamma_pred;
        nu_pred_path[t] = as_vec(nu_pred);
        Delta_pred_path[t] = Delta_pred;
        keep_path[t] = as_ivec(keep_t);
        lambda_path[t] = as_vec(lambda_t);
      }
    }
  } catch (const Fallback&) {
    return Rcpp::List::create(Rcpp::_["status"] = 1, Rcpp::_["ll"] = NA_REAL);
  }
  if (!store_path && store_keep)
    return Rcpp::List::create(Rcpp::_["status"] = 0, Rcpp::_["ll"] = ll,
                              Rcpp::_["keep_path"] = keep_path);
  if (!store_path)
    return Rcpp::List::create(Rcpp::_["status"] = 0, Rcpp::_["ll"] = ll);
  return Rcpp::List::create(
      Rcpp::_["status"] = 0, Rcpp::_["ll"] = ll,
      Rcpp::_["mu_pred_path"] = mu_pred_path,
      Rcpp::_["Sigma_pred_path"] = Sigma_pred_path,
      Rcpp::_["mu_filt_path"] = mu_filt_path,
      Rcpp::_["Sigma_filt_path"] = Sigma_filt_path,
      Rcpp::_["K_gauss_path"] = K_gauss_path,
      Rcpp::_["ZZ_path"] = ZZ_path,
      Rcpp::_["Gamma_filt_path"] = Gamma_filt_path,
      Rcpp::_["nu_filt_path"] = nu_filt_path,
      Rcpp::_["Delta_filt_path"] = Delta_filt_path,
      Rcpp::_["Gamma_pred_path"] = Gamma_pred_path,
      Rcpp::_["nu_pred_path"] = nu_pred_path,
      Rcpp::_["Delta_pred_path"] = Delta_pred_path,
      Rcpp::_["keep_path"] = keep_path,
      Rcpp::_["lambda_path"] = lambda_path);
}
