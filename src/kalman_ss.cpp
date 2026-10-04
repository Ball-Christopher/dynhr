// kalman_ss.cpp -- C++ steady-state Kalman recursion for dynhr
//
// Provides one exported function, kalman_ss_loop_cpp(), used by
// kalman_filter()'s steady-state branch when the DLL is available.
// The pure-R fallback (.kf_ss_loop_R) is kept in R/kalman-filter.R and
// gives identical results to the bit; see test-kalman-rcpp-parity.R.
//
// Args mirror the R helper one-for-one to keep the switch trivial.

#include <RcppArmadillo.h>
#include <limits>
#include <algorithm>
#include "kf_guard.h"
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
// complete-data, zero-a0, stationary-init case; everything else in that loop
// (missing data, a0 / P0, shock_scale, me_extra) runs on
// kalman_standard_general_loop_cpp() below, and the R loop is the
// options(dynhr.use_rcpp = FALSE) reference for both. Running the whole filter
// here avoids ~n_transient R-interpreter round-trips per likelihood eval,
// which dominate the per-draw MCMC cost on small state vectors (the
// n_state x n_state matrices are tiny, so the cost is call overhead, not BLAS).
//
// HH_full must already include the measurement-error diagonal (HH + me_diag);
// me_diag_vec carries that same diagonal SEPARATELY because the Joseph
// covariance update needs it as TRUE observation noise (P' += K diag(me) K').
// Pass an empty vector for me_diag_vec when there is no measurement error.
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
                              bool return_filtered,
                              const arma::vec& me_diag_vec,
                              double kalman_tol,
                              int upd = 0,
                              double guard_piv = -1.0,
                              double guard_r2 = -1.0,
                              double guard_ret = -1.0,
                              double guard_amp = -1.0) {
  const arma::uword n_state = TT.n_rows;
  const arma::uword n_T     = Y_minus_d.n_cols;
  // upd = 1: simple covariance update P' = T P T' + R Sigma R' - K F K'
  // (algebraically the Joseph recursion below for the optimal gain, with F
  // carrying the measurement-error diagonal); see kf_guard.h for the guard.
  const bool simple = (upd == 1);
  KfGuard guard(guard_piv, guard_r2, guard_ret, guard_amp);
  arma::mat QQs, Xs, base_s;
  if (simple) QQs = RR * Sigma_e * RR.t();
  arma::vec s(n_state, arma::fill::zeros);
  double loglik = 0.0;
  bool ok = true;
  bool ss_reached = false;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n_state, n_T);
  const arma::mat ZZt = ZZ.t();
  const bool has_me = me_diag_vec.n_elem > 0 &&
                      arma::any(arma::abs(me_diag_vec) > 0.0);
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
      // ...and a successful chol() is NOT a sufficient test: on a
      // stochastically singular system it can succeed with a pivot at
      // round-off and return a finite, badly wrong likelihood (measured
      // +292.7 against a correct -26.8). Dynare 7.1's rule
      // (kalman_filter.m, badly_conditioned_F): univariate only if
      //   rcond(F) < tol && (any(diag F < tol) || rcond(corr F) < tol),
      // rcond = exact reciprocal 1-norm condition from F and its inverse.
      // Mirrors .kf_F_singular() in R/kalman-filter.R; keep the two in step.
      {
        const arma::vec piv = Rc.diag();
        const arma::vec dF  = Ft.diag();
        if (!piv.is_finite() || piv.min() <= 0.0 || !dF.is_finite() ||
            dF.min() <= 0.0 || !Fi.is_finite()) { ok = false; break; }
        const double rc = 1.0 / (arma::norm(Ft, 1) * arma::norm(Fi, 1));
        if (rc < kalman_tol) {
          if (dF.min() < kalman_tol) { ok = false; break; }
          const arma::vec sg = arma::sqrt(dF);
          const arma::mat SG = sg * sg.t();
          const double rcc = 1.0 / (arma::norm(Ft / SG, 1) *
                                    arma::norm(Fi % SG, 1));
          if (rcc < kalman_tol) { ok = false; break; }
        }
      }
      double ldf = 2.0 * arma::sum(arma::log(Rc.diag()));
      Fiv = Fi * v;
      double ll = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      if (simple) {
        guard.note_F(Rc, Ft, Fi);
        if (guard.tripped) break;
      }
      // K = (TT*PZ + SS) * Fi
      TPZ = TT * PZ;
      TPZ += SS;
      K = TPZ * Fi;                       // n_state x n_obs
      // s_n = TT*s + K*v
      s_n = TT * s;
      Kv  = K * v;
      s_n += Kv;
      if (simple) {
        // P' = T P T' + R Sigma R' - K F K', with K F = TPZ.
        base_s = TT * P * TT.t();
        base_s += QQs;
        guard.pending = arma::sum((ZZ * base_s) % ZZ, 1);
        Xs = K * TPZ.t();
        P_n = base_s - 0.5 * (Xs + Xs.t());
        P_n = 0.5 * (P_n + P_n.t());
        guard.note_P(base_s, P_n);
        if (guard.tripped) break;
      } else {
      // TmKZ = TT - K*ZZ ;  RmKD = RR - K*DD
      KZ = K * ZZ;   TmKZ = TT - KZ;
      KD = K * DD;   RmKD = RR - KD;
      // P_n = TmKZ*P*TmKZ.t() + RmKD*Sigma_e*RmKD.t()
      tmpA = TmKZ * P;
      P_n  = tmpA * TmKZ.t();
      tmpB = RmKD * Sigma_e;
      P_n += tmpB * RmKD.t();
      // TRUE measurement-noise law (F3-D): P' += K diag(me) K'.
      if (has_me) P_n += (K * arma::diagmat(me_diag_vec)) * K.t();
      P_n = 0.5 * (P_n + P_n.t());
      }
      s = s_n;
      // RELATIVE lock, as in the R loop (kalman_filter): an absolute ss_tol
      // froze the gain at t = 2 on a model whose P is ~1e-15.
      if (t > 0 && arma::abs(P_n - P).max() < ss_tol * arma::abs(P_n).max()) {
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
                      // Var(s_T | y_1..T): the state hand-off a split sample
                      // needs as its P0 (kalman_filter's `final_cov`). Frozen
                      // at the steady-state value once the lock engages, which
                      // is the same approximation the loglik already makes.
                      _["P"]          = P,
                      _["filtered"]   = return_filtered ? Rcpp::wrap(filtered)
                                                        : R_NilValue,
                      _["ok"]         = ok,
                      _["ss_reached"] = ss_reached,
                      _["guard_tripped"] = guard.tripped,
                      _["ind_piv"]    = guard.piv,
                      _["ind_r2"]     = guard.r2,
                      _["ind_ret"]    = guard.ret,
                      _["ind_amp"]    = guard.amp,
                      _["ind_bad"]    = guard.bad);
}


// Dynare's badly-conditioned-F rule, as in kalman_standard_loop_cpp() above
// and .kf_F_singular() in R/kalman-filter.R (keep the three in step). `Rc` is
// the upper Cholesky factor of `Ft`, `Fi` its inverse. TRUE = treat as
// singular (the caller hands the evaluation to the univariate filter).
// External linkage: kalman_struct.cpp's structured kernel applies the same rule.
bool kf_general_F_singular(const arma::mat& Rc, const arma::mat& Ft,
                           const arma::mat& Fi, double kalman_tol) {
  const arma::vec piv = Rc.diag();
  const arma::vec dF  = Ft.diag();
  if (!piv.is_finite() || piv.min() <= 0.0 || !dF.is_finite() ||
      dF.min() <= 0.0 || !Fi.is_finite()) return true;
  const double rc = 1.0 / (arma::norm(Ft, 1) * arma::norm(Fi, 1));
  if (rc < kalman_tol) {
    if (dF.min() < kalman_tol) return true;
    const arma::vec sg = arma::sqrt(dF);
    const arma::mat SG = sg * sg.t();
    const double rcc = 1.0 / (arma::norm(Ft / SG, 1) *
                              arma::norm(Fi % SG, 1));
    if (rcc < kalman_tol) return true;
  }
  return false;
}


// General standard Kalman filter: the whole of kalman_filter()'s per-step
// "standard" loop in one C++ call, for every case the complete-data kernel
// above does not take --
//   * missing observations (NaN entries of Y_minus_d): a partially observed
//     period updates with the observed rows of ZZ / DD / the ME diagonal
//     only, a fully missing period is a pure prediction step;
//   * a given initial state s and covariance P (a0, a user P0, the kappa
//     prior, or the hand-off of the exact-diffuse phase, which also supplies
//     t_start > 1 and that phase's loglik as init_loglik);
//   * per-period shock scales (shock_scale, n_exo x T; empty = none) and
//     per-period extra measurement-error variances (me_extra, n_obs x T;
//     empty = none).
// The steady-state lock is the R loop's: it engages only on a complete,
// time-invariant period (never in a run with shock_scale or me_extra), at
// t > 1 under the relative tolerance, and a missing period releases it --
// the covariance then leaves its fixed point, so the Riccati recursion runs
// again until it re-converges and may re-lock. Each branch mirrors the
// corresponding branch of the R loop (the reference behind
// options(dynhr.use_rcpp = FALSE)), including which steps apply the ll_min
// floor. The missing-observation pattern is read off the NaNs here, once per
// period: it costs O(n_obs) per period against the O(n_state^3) step.
//
// HH, QQ and SS are the BASELINE (unscaled, ME-free) blocks; me_vec is the
// length-n_obs measurement-error diagonal (zeros when there is none).
//
// [[Rcpp::export]]
List kalman_standard_general_loop_cpp(const arma::mat& Y_minus_d,
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
                                      int upd = 0,
                                      double guard_piv = -1.0,
                                      double guard_r2 = -1.0,
                                      double guard_ret = -1.0,
                              double guard_amp = -1.0) {
  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_T     = Y_minus_d.n_cols;
  const bool has_mx = me_extra.n_elem > 0;
  const bool has_sc = shock_scale.n_elem > 0;
  if (Y_minus_d.n_rows != n_obs)
    Rcpp::stop("kalman_standard_general_loop_cpp: Y_minus_d must have n_obs rows");
  if (has_mx && (me_extra.n_rows != n_obs || me_extra.n_cols != n_T))
    Rcpp::stop("kalman_standard_general_loop_cpp: me_extra must be n_obs x T");
  if (has_sc && (shock_scale.n_rows != Sigma_e.n_rows ||
                 shock_scale.n_cols != n_T))
    Rcpp::stop("kalman_standard_general_loop_cpp: shock_scale must be n_exo x T");
  if (me_vec.n_elem != n_obs)
    Rcpp::stop("kalman_standard_general_loop_cpp: me_vec must have length n_obs");
  if (s.n_elem != n_state || P.n_rows != n_state || P.n_cols != n_state)
    Rcpp::stop("kalman_standard_general_loop_cpp: s / P do not match TT");
  if (t_start < 1)
    Rcpp::stop("kalman_standard_general_loop_cpp: t_start must be >= 1");

  double loglik = init_loglik;
  bool ok = true;
  bool ss_reached = false;
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n_state, n_T);
  const arma::mat ZZt = ZZ.t();
  const arma::mat HH_full = HH + arma::diagmat(me_vec);
  const bool has_me = arma::any(arma::abs(me_vec) > 0.0);
  const double log2pi = std::log(2.0 * M_PI);
  arma::mat  K_ss, F_inv_ss;
  double     ll_ss_const = 0.0;

  // upd = 1: simple covariance update at every observed period (see
  // kalman_standard_loop_cpp); the prediction-only step is the same in both.
  const bool simple = (upd == 1);
  KfGuard guard(guard_piv, guard_r2, guard_ret, guard_amp);
  arma::mat base_s, Xs, QQ_t;
  arma::vec y, v, Fiv, Kv, s_n, me_t;
  arma::mat PZ, Ft, Rc, Fi, TPZ, K, KZ, KD, TmKZ, RmKD, P_n, tmpA, tmpB;
  arma::mat Se_t, HH_t, SS_t;
  arma::uvec obs_ok;

  for (arma::uword t = static_cast<arma::uword>(t_start - 1); t < n_T; ++t) {
    y = Y_minus_d.col(t);
    arma::uword n_ok = 0;
    for (arma::uword i = 0; i < n_obs; ++i) if (!std::isnan(y[i])) ++n_ok;
    v = y - ZZ * s;

    // This period's shock covariance (scaled when shock_scale is supplied).
    if (has_sc) {
      const arma::vec sc = shock_scale.col(t);
      Se_t = Sigma_e % (sc * sc.t());
    }
    const arma::mat& Se = has_sc ? Se_t : Sigma_e;

    if (n_ok < n_obs) {
      // ---- missing observation(s): drop them, release the lock ----------
      ss_reached = false;
      if (n_ok == 0) {
        guard.pending.reset();
        s = TT * s;
        if (has_sc) P = TT * P * TT.t() + RR * Se * RR.t();
        else        P = TT * P * TT.t() + QQ;
        P = 0.5 * (P + P.t());
        if (return_filtered) filtered.col(t) = s;
        continue;
      }
      obs_ok = arma::find(y == y);               // NaN != NaN; +-Inf kept, as is.na()
      const arma::mat ZZ_o = ZZ.rows(obs_ok);
      const arma::mat DD_o = DD.rows(obs_ok);
      const arma::vec v_o  = v.elem(obs_ok);
      me_t = me_vec.elem(obs_ok);
      if (has_mx) {
        const arma::vec mx = me_extra.col(t);
        me_t += mx.elem(obs_ok);
      }
      Ft = ZZ_o * P * ZZ_o.t() + DD_o * Se * DD_o.t() + arma::diagmat(me_t);
      Ft = 0.5 * (Ft + Ft.t());
      if (!arma::chol(Rc, Ft)) { ok = false; break; }
      if (!arma::inv_sympd(Fi, Ft)) { ok = false; break; }
      if (kf_general_F_singular(Rc, Ft, Fi, kalman_tol)) { ok = false; break; }
      const double ldf = 2.0 * arma::sum(arma::log(Rc.diag()));
      // ll_const bakes in the full n_obs; correct it up by the number of
      // missing components (the R loop applies no ll_min floor here).
      loglik += ll_const + 0.5 * static_cast<double>(n_obs - n_ok) * log2pi -
        0.5 * (ldf + arma::dot(v_o, Fi * v_o));
      const arma::mat N_o = TT * P * ZZ_o.t() + RR * Se * DD_o.t();
      K = N_o * Fi;
      s = TT * s + K * v_o;
      if (simple) {
        guard.note_F(Rc, Ft, Fi, obs_ok.memptr());
        if (guard.tripped) break;
        QQ_t = has_sc ? arma::mat(RR * Se * RR.t()) : QQ;
        base_s = TT * P * TT.t() + QQ_t;
        guard.pending = arma::sum((ZZ * base_s) % ZZ, 1);
        Xs = K * N_o.t();
        P_n = base_s - 0.5 * (Xs + Xs.t());
        P = 0.5 * (P_n + P_n.t());
        guard.note_P(base_s, P);
        if (guard.tripped) break;
      } else {
        TmKZ = TT - K * ZZ_o;
        RmKD = RR - K * DD_o;
        P_n  = TmKZ * P * TmKZ.t() + RmKD * Se * RmKD.t();
        if (arma::any(me_t != 0.0)) P_n += (K * arma::diagmat(me_t)) * K.t();
        P = 0.5 * (P_n + P_n.t());
      }
      if (return_filtered) filtered.col(t) = s;
      continue;
    }

    if (ss_reached) {
      Fiv = F_inv_ss * v;
      const double ll = ll_ss_const - 0.5 * arma::dot(v, Fiv);
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      s_n = TT * s;
      Kv  = K_ss * v;
      s_n += Kv;
      s = s_n;
      if (return_filtered) filtered.col(t) = s;
      continue;
    }

    const bool mx_t = has_mx && arma::any(me_extra.col(t) != 0.0);
    if (has_sc || mx_t) {
      // ---- time-varying period: per-period Q / H, never locks -----------
      me_t = me_vec;
      if (has_mx) me_t += me_extra.col(t);
      if (has_sc) { HH_t = DD * Se * DD.t(); SS_t = RR * Se * DD.t(); }
      PZ = P * ZZt;
      Ft = ZZ * PZ + (has_sc ? HH_t : HH) + arma::diagmat(me_t);
      Ft = 0.5 * (Ft + Ft.t());
      if (!arma::chol(Rc, Ft)) { ok = false; break; }
      if (!arma::inv_sympd(Fi, Ft)) { ok = false; break; }
      if (kf_general_F_singular(Rc, Ft, Fi, kalman_tol)) { ok = false; break; }
      const double ldf = 2.0 * arma::sum(arma::log(Rc.diag()));
      const double ll = ll_const - 0.5 * (ldf + arma::dot(v, Fi * v));
      if (!std::isfinite(ll) || ll < ll_min) { ok = false; break; }
      loglik += ll;
      const arma::mat N_t = TT * PZ + (has_sc ? SS_t : SS);
      K = N_t * Fi;
      s = TT * s + K * v;
      if (simple) {
        guard.note_F(Rc, Ft, Fi);
        if (guard.tripped) break;
        QQ_t = has_sc ? arma::mat(RR * Se * RR.t()) : QQ;
        base_s = TT * P * TT.t() + QQ_t;
        guard.pending = arma::sum((ZZ * base_s) % ZZ, 1);
        Xs = K * N_t.t();
        P_n = base_s - 0.5 * (Xs + Xs.t());
        P = 0.5 * (P_n + P_n.t());
        guard.note_P(base_s, P);
        if (guard.tripped) break;
      } else {
        TmKZ = TT - K * ZZ;
        RmKD = RR - K * DD;
        P_n  = TmKZ * P * TmKZ.t() + RmKD * Se * RmKD.t();
        if (arma::any(me_t != 0.0)) P_n += (K * arma::diagmat(me_t)) * K.t();
        P = 0.5 * (P_n + P_n.t());
      }
      if (return_filtered) filtered.col(t) = s;
      continue;
    }

    // ---- complete, time-invariant period: the complete-data kernel's
    //      step, operation for operation -------------------------------
    PZ = P * ZZt;
    Ft = ZZ * PZ;
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
    TPZ = TT * PZ;
    TPZ += SS;
    K = TPZ * Fi;
    s_n = TT * s;
    Kv  = K * v;
    s_n += Kv;
    if (simple) {
      guard.note_F(Rc, Ft, Fi);
      if (guard.tripped) break;
      base_s = TT * P * TT.t() + QQ;
      guard.pending = arma::sum((ZZ * base_s) % ZZ, 1);
      Xs = K * TPZ.t();
      P_n = base_s - 0.5 * (Xs + Xs.t());
      P_n = 0.5 * (P_n + P_n.t());
      guard.note_P(base_s, P_n);
      if (guard.tripped) break;
    } else {
    KZ = K * ZZ;   TmKZ = TT - KZ;
    KD = K * DD;   RmKD = RR - KD;
    tmpA = TmKZ * P;
    P_n  = tmpA * TmKZ.t();
    tmpB = RmKD * Sigma_e;
    P_n += tmpB * RmKD.t();
    if (has_me) P_n += (K * arma::diagmat(me_vec)) * K.t();
    P_n = 0.5 * (P_n + P_n.t());
    }
    s = s_n;
    if (!has_mx && !has_sc && t > 0 &&
        arma::abs(P_n - P).max() < ss_tol * arma::abs(P_n).max()) {
      ss_reached  = true;
      K_ss        = K;
      F_inv_ss    = Fi;
      ll_ss_const = ll_const - 0.5 * ldf;
    }
    P = P_n;
    if (return_filtered) filtered.col(t) = s;
  }

  return List::create(_["loglik"]     = loglik,
                      _["s"]          = s,
                      _["P"]          = P,
                      _["filtered"]   = return_filtered ? Rcpp::wrap(filtered)
                                                        : R_NilValue,
                      _["ok"]         = ok,
                      _["ss_reached"] = ss_reached,
                      _["guard_tripped"] = guard.tripped,
                      _["ind_piv"]    = guard.piv,
                      _["ind_r2"]     = guard.r2,
                      _["ind_ret"]    = guard.ret,
                      _["ind_amp"]    = guard.amp,
                      _["ind_bad"]    = guard.bad);
}
