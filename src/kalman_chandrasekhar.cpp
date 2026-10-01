// kalman_chandrasekhar.cpp -- the Chandrasekhar increment recursion of
// kalman_filter() (METHOD 2 in R/kalman-filter.R) in C++, for the complete-
// data, stationary-initialisation case.
//
// Mirrors the R loop one-for-one:
//   F_{t+1} = F_t + ZW M ZW'
//   K_{t+1} = K_t + (TW - K_t ZW) M ZW' F_{t+1}^{-1}
//   W_{t+1} = TW - K_{t+1} ZW
//   M_{t+1} = M_t + M ZW' F_t^{-1} ZW M
// initialised at P_1 = P0 (the Lyapunov prior), W_1 = K_1, M_1 = -F_1, with
// the same relative steady-state lock (max|dP| < ss_tol * max|P|); after the
// lock the constant-gain tail runs here too. P is never formed beyond the
// running tracker the lock test needs, so the cost per step is
// O(n_state^2 n_obs), against O(n_state^3) for the standard filter.
//
// Any non-positive-definite F, non-finite or implausibly low period
// log-likelihood returns ok = FALSE and the R caller reports the failure.

#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// [[Rcpp::export]]
List kalman_chandrasekhar_loop_cpp(const arma::mat& Y_minus_d,
                                   const arma::mat& ZZ,
                                   const arma::mat& TT,
                                   const arma::mat& SS,
                                   const arma::mat& HH_full,
                                   const arma::mat& P0,
                                   arma::vec s,
                                   double ll_const,
                                   double ss_tol,
                                   double ll_min,
                                   bool return_filtered) {
  const arma::uword n_state = TT.n_rows;
  const arma::uword n_T     = Y_minus_d.n_cols;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n_state, n_T);
  double loglik = 0.0;
  int ss_step = 0;
  bool ok = true;

  auto result = [&]() {
    return List::create(_["loglik"]   = loglik,
                        _["s"]        = s,
                        _["filtered"] = return_filtered ? Rcpp::wrap(filtered)
                                                        : R_NilValue,
                        _["ok"]       = ok,
                        _["ss_step"]  = ss_step);
  };

  const arma::mat ZZt = ZZ.t();
  arma::mat Rc, F_inv, F_new_inv;

  // Exact initialisation at P_1 = P0.
  const arma::mat PZ0 = P0 * ZZt;
  arma::mat F_mat = ZZ * PZ0 + HH_full;
  F_mat = 0.5 * (F_mat + F_mat.t());
  if (!arma::chol(Rc, F_mat) || !arma::inv_sympd(F_inv, F_mat)) {
    ok = false; return result();
  }
  double log_det_F = 2.0 * arma::sum(arma::log(Rc.diag()));
  arma::mat K = (TT * PZ0 + SS) * F_inv;
  arma::mat W = K;
  arma::mat M = -F_mat;
  arma::mat P_ch = P0;

  arma::vec v;
  arma::mat ZW, TW, WM, ZWM, dP, F_new, tZWM, K_new, W_new, M_new, tmp;

  for (arma::uword t = 0; t < n_T; ++t) {
    v = Y_minus_d.col(t) - ZZ * s;
    double ll = ll_const - 0.5 * (log_det_F + arma::dot(v, F_inv * v));
    if (!std::isfinite(ll) || ll < ll_min) { ok = false; return result(); }
    loglik += ll;
    s = TT * s + K * v;
    if (return_filtered) filtered.col(t) = s;
    if (t + 1 == n_T) break;

    ZW  = ZZ * W;
    TW  = TT * W;
    WM  = W * M;
    ZWM = ZZ * WM;

    // Steady-state lock (1-based period t+1 > 1 in the R loop).
    dP = WM * W.t();
    P_ch += dP;
    if (t > 0 && arma::abs(dP).max() < ss_tol * arma::abs(P_ch).max()) {
      ss_step = static_cast<int>(t) + 1;
      const double ll_ss_const = ll_const - 0.5 * log_det_F;
      for (arma::uword u = t + 1; u < n_T; ++u) {
        v = Y_minus_d.col(u) - ZZ * s;
        const double llu = ll_ss_const - 0.5 * arma::dot(v, F_inv * v);
        if (!std::isfinite(llu) || llu < ll_min) { ok = false; return result(); }
        loglik += llu;
        s = TT * s + K * v;
        if (return_filtered) filtered.col(u) = s;
      }
      break;
    }

    // Increment update.
    F_new = F_mat + ZWM * ZW.t();
    F_new = 0.5 * (F_new + F_new.t());
    if (!arma::chol(Rc, F_new) || !arma::inv_sympd(F_new_inv, F_new)) {
      ok = false; return result();
    }
    const double log_det_F_new = 2.0 * arma::sum(arma::log(Rc.diag()));

    tZWM  = ZWM.t();                       // M ZW'
    tmp   = TW - K * ZW;
    K_new = K + tmp * (tZWM * F_new_inv);
    W_new = TW - K_new * ZW;
    M_new = M + tZWM * F_inv * ZWM;
    M_new = 0.5 * (M_new + M_new.t());

    K = K_new; F_mat = F_new; F_inv = F_new_inv;
    log_det_F = log_det_F_new; W = W_new; M = M_new;
  }
  return result();
}
