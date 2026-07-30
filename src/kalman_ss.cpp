// kalman_ss.cpp -- C++ steady-state Kalman recursion for dynhr
//
// Provides one exported function, kalman_ss_loop_cpp(), used by
// kalman_filter()'s steady-state branch when the DLL is available.
// The pure-R fallback (.kf_ss_loop_R) is kept in R/kalman-filter.R and
// gives identical results to the bit; see test-kalman-rcpp-parity.R.
//
// Args mirror the R helper one-for-one to keep the switch trivial.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// [[Rcpp::export]]
List kalman_ss_loop_cpp(const arma::mat& Y_minus_d,
                        const arma::mat& ZZ,
                        const arma::mat& TT,
                        const arma::mat& K_ss,
                        const arma::mat& F_inv_ss,
                        double ll_ss_const,
                        arma::vec s,
                        int start_t,
                        int end_t,
                        bool return_filtered) {
  double loglik = 0.0;
  const arma::uword n_state = s.n_elem;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n_state, Y_minus_d.n_cols);
  bool ok = true;

  // Convert to 0-based; loop runs t in [t0, t1).
  const int t0 = start_t - 1;
  const int t1 = end_t;
  if (t0 < 0)
    Rcpp::stop("kalman_ss_loop_cpp: start_t must be >= 1");
  if (t1 > static_cast<int>(Y_minus_d.n_cols))
    Rcpp::stop("kalman_ss_loop_cpp: end_t exceeds Y_minus_d columns");
  if (t0 >= t1)
    Rcpp::stop("kalman_ss_loop_cpp: start_t must be <= end_t");

  for (int t = t0; t < t1; ++t) {
    arma::vec v   = Y_minus_d.col(t) - ZZ * s;
    arma::vec Fv  = F_inv_ss * v;
    double    llt = ll_ss_const - 0.5 * arma::dot(v, Fv);
    if (!std::isfinite(llt) || llt < -1e8) { ok = false; break; }
    // NOTE: the -1e8 threshold here mirrors .KF_LL_MIN in R/dynhr-package.R;
    // C++ cannot directly reference R constants. Keep in sync manually.
    loglik += llt;
    s = TT * s + K_ss * v;
    if (return_filtered) filtered.col(t) = s;
  }

  return List::create(_["loglik"]   = loglik,
                      _["s"]        = s,
                      _["filtered"] = return_filtered ? Rcpp::wrap(filtered)
                                                      : R_NilValue,
                      _["ok"]       = ok);
}


// Full standard Kalman filter (transient Riccati + steady-state lock + tail)
// in one C++ call. Mirrors METHOD 3 in R/kalman-filter.R exactly for the
// no-missing-data fast path; the R loop there still handles missing data and
// the options(dynhr.use_rcpp = FALSE) debug fallback. Running the whole filter
// here avoids ~n_transient R-interpreter round-trips per likelihood eval,
// which dominate the per-draw MCMC cost on small state vectors (the
// n_state x n_state matrices are tiny, so the cost is call overhead, not BLAS).
//
// HH_full must already include the measurement-error ridge (HH + me_diag).
// Bit-parity with the R path is covered by test-kalman-rcpp-parity.R.
//
// [[Rcpp::export]]
List kalman_standard_loop_cpp(const arma::mat& Y_minus_d,
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
                              bool return_filtered) {
  const arma::uword n_state = TT.n_rows;
  const arma::uword n_T     = Y_minus_d.n_cols;
  arma::vec s(n_state, arma::fill::zeros);
  double loglik = 0.0;
  bool ok = true;
  bool ss_reached = false;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n_state, n_T);
  const arma::mat ZZt = ZZ.t();
  arma::mat  K_ss, F_inv_ss;
  double     ll_ss_const = 0.0;

  // Work buffers hoisted out of the loop. Armadillo's operator= reuses the
  // existing allocation when the result dimensions are unchanged, so after the
  // first iteration these incur no heap traffic -- the transient branch runs
  // on tiny matrices where per-iteration malloc/free and expression-template
  // temporaries dominate (it is allocation-bound, not BLAS-bound). The
  // arithmetic and its association are kept identical to the original compound
  // expressions, so the result is bit-for-bit unchanged (test-kalman-rcpp-
  // parity.R asserts this).
  arma::vec v, Fiv, Kv, s_n;
  arma::mat PZ, Ft, Rc, Fi, TPZ, K, KZ, KD, TmKZ, RmKD, P_n, tmpA, tmpB;

  for (arma::uword t = 0; t < n_T; ++t) {
    v = Y_minus_d.col(t) - ZZ * s;

    if (!ss_reached) {
      PZ = P * ZZt;                       // n_state x n_obs
      Ft = ZZ * PZ;                       // n_obs x n_obs (= ZZ*P*ZZt)
      Ft += HH_full;
      Ft = 0.5 * (Ft + Ft.t());
      if (!arma::chol(Rc, Ft)) { ok = false; break; }   // upper: Ft = Rc'Rc
      // Non-throwing form: chol() success does not imply inv_sympd() success
      // (see kalman_adjoint.cpp) -- degrade gracefully instead of throwing.
      if (!arma::inv_sympd(Fi, Ft)) { ok = false; break; }
      double ldf = 2.0 * arma::sum(arma::log(Rc.diag()));
      Fiv = Fi * v;
      double ll = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      // K = (TT*PZ + SS) * Fi
      TPZ = TT * PZ;
      TPZ += SS;
      K = TPZ * Fi;                       // n_state x n_obs
      // s_n = TT*s + K*v
      s_n = TT * s;
      Kv  = K * v;
      s_n += Kv;
      // TmKZ = TT - K*ZZ ;  RmKD = RR - K*DD
      KZ = K * ZZ;   TmKZ = TT - KZ;
      KD = K * DD;   RmKD = RR - KD;
      // P_n = TmKZ*P*TmKZ.t() + RmKD*Sigma_e*RmKD.t()
      tmpA = TmKZ * P;
      P_n  = tmpA * TmKZ.t();
      tmpB = RmKD * Sigma_e;
      P_n += tmpB * RmKD.t();
      P_n = 0.5 * (P_n + P_n.t());
      s = s_n;
      if (t > 0 && arma::abs(P_n - P).max() < ss_tol) {
        ss_reached  = true;
        K_ss        = K;
        F_inv_ss    = Fi;
        ll_ss_const = ll_const - 0.5 * ldf;
      }
      P = P_n;
    } else {
      Fiv = F_inv_ss * v;
      double ll = ll_ss_const - 0.5 * arma::dot(v, Fiv);
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      s_n = TT * s;
      Kv  = K_ss * v;
      s_n += Kv;
      s = s_n;
    }

    if (return_filtered) filtered.col(t) = s;
  }

  return List::create(_["loglik"]     = loglik,
                      _["s"]          = s,
                      _["filtered"]   = return_filtered ? Rcpp::wrap(filtered)
                                                        : R_NilValue,
                      _["ok"]         = ok,
                      _["ss_reached"] = ss_reached);
}
