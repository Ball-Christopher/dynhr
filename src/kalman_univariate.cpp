// kalman_univariate.cpp -- C++ univariate (sequential) Kalman filter for dynhr
//
// Koopman & Durbin (2000) univariate treatment, run on the AUGMENTED state
// x_t = [s_{t-1}; eps_t] so that dynhr's correlated measurement noise
// (y_t = Z s_{t-1} + D eps_t shares shocks with the transition) becomes a
// measurement with EXACTLY zero noise:
//
//   x_{t+1} = Tb x_t + [0; I] eps_{t+1},   Tb = [ TT RR ]
//                                               [ 0  0  ]
//   y_t     = Zb x_t,                      Zb = [ Z  D ]
//
// Observables are processed one at a time, so a singular innovation
// covariance F never needs inverting: zero-variance components are skipped.
// The same recursion covers the exact-diffuse initialization phase (the
// Dynare kalman_algo=4 analog) for F_inf of ANY rank -- including the
// multivariate diffuse filter's unsupported "Case C" (F_inf singular but
// nonzero, see .kf_diffuse_phase in R/kalman-filter.R).
//
// One exported function, consumed by .kf_univariate_dispatch() in
// R/kalman-filter.R. The pure-R fallback .kf_univariate_loop_R() mirrors it
// one-for-one; parity is asserted in test-kalman-univariate.R.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// [[Rcpp::export]]
List kalman_univariate_loop_cpp(const arma::mat& Y_minus_d,
                                const arma::mat& Zb,
                                const arma::mat& Tb,
                                const arma::mat& QQb,
                                arma::vec a,
                                arma::mat P_star,
                                arma::mat P_inf,
                                double me_variance,
                                double kalman_tol,
                                double diffuse_tol,
                                double conv_tol,
                                int max_diffuse,
                                double ll_min,
                                bool return_filtered,
                                int n_state,
                                bool ss_lock) {
  const arma::uword n_obs = Y_minus_d.n_rows;
  const arma::uword n_T   = Y_minus_d.n_cols;
  const double log2pi = std::log(2.0 * arma::datum::pi);

  double loglik       = 0.0;
  bool   ok           = true;
  bool   diffuse      = P_inf.n_elem > 0 && arma::abs(P_inf).max() > 0.0;
  bool   diffuse_failed = false;
  int    d_diffuse    = NA_INTEGER;
  arma::vec ll_contrib(n_T, arma::fill::zeros);
  // Per-period count of observation components skipped because their
  // forecast variance was (numerically) zero (see the skip branch below).
  arma::ivec n_skipped(n_T, arma::fill::zeros);
  arma::mat filtered;
  if (return_filtered) filtered.zeros(n_state, n_T);

  // Work buffers hoisted out of the loop (same rationale as in
  // kalman_ss.cpp: tiny matrices, allocation-bound not BLAS-bound). KKt is the
  // rank-1 outer product reused by the per-observation covariance downdate so
  // the hot stationary-tail path allocates no nb x nb temporary per observable.
  arma::vec K_star, K_inf;
  arma::mat KKt;

  // Time-update block structure. Tb = [[TT RR],[0 0]] (the bottom n_exo rows
  // are exactly zero), so  Tb*P*Tb' = [[G*P*G', 0],[0, 0]]  with G = the top
  // n_state rows of Tb, and QQb = blkdiag(0, Sigma_e) then fills the eps block
  // with the constant Sigma_e. Computing only the n_state-block G*P*G' instead
  // of the full nb*nb*nb product halves the (dominant) time-update flops and
  // yields the same matrix bit-for-bit (the dropped products are exact zeros).
  // G and the Sigma_e block are loop-invariant -- hoist them.
  const arma::uword ns    = static_cast<arma::uword>(n_state);
  const arma::uword nb    = Tb.n_rows;
  const arma::uword n_exo = nb - ns;
  const arma::mat   G     = Tb.head_rows(ns);              // [TT RR], ns x nb
  arma::mat Sig_e;
  if (n_exo > 0) Sig_e = QQb.submat(ns, ns, nb - 1, nb - 1);   // = Sigma_e

  // -- Steady-state lock (opt-in via ss_lock) -----------------------------
  // Once the start-of-period P_star converges (stationary tail), the within-
  // period sequence of per-observable Kalman gains is fixed. Freeze that
  // sequence and skip ALL P updates thereafter -- the same approximation the
  // multivariate steady-state filter (kalman_standard_loop_cpp) makes. Engaged
  // only after the diffuse phase; the caller gates ss_lock on complete data so
  // the present-observable pattern is constant period to period.
  bool ss_engaged = false, ss_capture = false, have_prev = false;
  arma::mat P_prev;
  std::vector<arma::vec>   g_seq;     // frozen Kalman gains K_i
  std::vector<double>      invF_seq;  // 1 / F_i
  std::vector<double>      cst_seq;   // 0.5 * (log2pi + log F_i)
  std::vector<arma::uword> idx_seq;   // observable index i

  for (arma::uword t = 0; t < n_T; ++t) {
    if (ss_engaged) {
      // Locked fast path: frozen gains, live innovations, no P update.
      double ll_t = 0.0;
      for (std::size_t s = 0; s < idx_seq.size(); ++s) {
        const double y_i = Y_minus_d(idx_seq[s], t);
        const double v   = y_i - arma::dot(Zb.row(idx_seq[s]), a);
        ll_t -= cst_seq[s] + 0.5 * v * v * invF_seq[s];
        a    += g_seq[s] * (v * invF_seq[s]);
      }
      if (!std::isfinite(ll_t) || ll_t < ll_min) { ok = false; break; }
      loglik += ll_t; ll_contrib(t) = ll_t;
      arma::vec a_top = G * a; a.head(ns) = a_top;
      if (n_exo > 0) a.tail(n_exo).zeros();
      if (return_filtered) filtered.col(t) = a.head(ns);
      continue;
    }
    if (diffuse && static_cast<int>(t) >= max_diffuse) {
      // Mirrors .kf_diffuse_phase's hard cap: P_inf never decayed, the
      // diffuse subspace is (numerically) unobservable. The caller falls
      // back to lik_init = "kappa".
      diffuse_failed = true;
      break;
    }

    double ll_t = 0.0;
    for (arma::uword i = 0; i < n_obs; ++i) {
      const double y_i = Y_minus_d(i, t);
      if (!std::isfinite(y_i)) continue;            // NA observable: skip
      const arma::rowvec Zi = Zb.row(i);
      const double v = y_i - arma::dot(Zi, a);
      K_star = P_star * Zi.t();
      const double F_star = arma::dot(Zi, K_star) + me_variance;

      if (diffuse) {
        K_inf = P_inf * Zi.t();
        const double F_inf = arma::dot(Zi, K_inf);
        if (F_inf > diffuse_tol * std::max(1.0, F_star)) {
          // Diffuse update (DK 2012 sec. 7.2.5). No 2*pi and no quadratic
          // term: same renormalization convention as the multivariate
          // Case B in .kf_diffuse_phase (see the NOTE there).
          ll_t   -= 0.5 * std::log(F_inf);
          a      += K_inf * (v / F_inf);
          P_star += (K_inf * K_inf.t()) * (F_star / (F_inf * F_inf))
                  - (K_star * K_inf.t() + K_inf * K_star.t()) / F_inf;
          P_inf  -= (K_inf * K_inf.t()) / F_inf;
          continue;
        }
        // F_inf ~ 0 for this observable: fall through to the standard
        // scalar update on P_star (the per-observable analog of Case A).
      }

      if (F_star > kalman_tol) {
        ll_t   -= 0.5 * (log2pi + std::log(F_star) + v * v / F_star);
        a      += K_star * (v / F_star);
        // Sherman-Morrison rank-1 covariance downdate, reusing KKt to avoid an
        // nb x nb allocation per observable. Bit-identical to
        //   P_star -= (K_star * K_star.t()) / F_star;
        // (same products, same scalar division, same subtraction order).
        KKt = K_star * K_star.t();
        KKt /= F_star;
        P_star -= KKt;
        if (ss_capture) {            // record the converged gain for this obs
          idx_seq.push_back(i);
          g_seq.push_back(K_star);
          invF_seq.push_back(1.0 / F_star);
          cst_seq.push_back(0.5 * (log2pi + std::log(F_star)));
        }
      }
      else {
        // (numerically) zero innovation variance -- the observable is
        // an exact linear combination of already-processed information.
        // Skip it gracefully (no inversion), per Koopman & Durbin (2000).
        // This is what makes singular F a non-event on this path. Counted so
        // kalman_filter() can report it as structured diagnostics rather than
        // leaving the caller to infer it.
        n_skipped(t) += 1;
      }
    }

    if (!std::isfinite(ll_t) || ll_t < ll_min) { ok = false; break; }
    // NOTE: ll_min mirrors .KF_LL_MIN in R/dynhr-package.R (passed in).
    loglik += ll_t;
    ll_contrib(t) = ll_t;

    // Time update (block-structured; see the hoist above). Bit-identical to
    //   a      = Tb * a;
    //   P_star = Tb * P_star * Tb.t() + QQb;
    // because Tb's bottom n_exo rows are exactly zero (their products are exact
    // zeros) and QQb is blkdiag(0, Sigma_e). The 0.5*(.+.t()) symmetrisation is
    // kept identical so the result matches the R fallback to the bit.
    {
      arma::vec a_top = G * a;
      a.head(ns) = a_top;
      if (n_exo > 0) a.tail(n_exo).zeros();
      arma::mat M = G * P_star * G.t();          // = (Tb P Tb')[0:ns, 0:ns]
      P_star.zeros();
      P_star.submat(0, 0, ns - 1, ns - 1) = M;
      if (n_exo > 0) P_star.submat(ns, ns, nb - 1, nb - 1) = Sig_e;
      P_star = 0.5 * (P_star + P_star.t());
    }
    if (diffuse) {
      arma::mat Mi = G * P_inf * G.t();          // = (Tb P_inf Tb')[0:ns, 0:ns]
      P_inf.zeros();
      P_inf.submat(0, 0, ns - 1, ns - 1) = Mi;
      P_inf = 0.5 * (P_inf + P_inf.t());
      if (arma::abs(P_inf).max() <
          conv_tol * std::max(1.0, arma::abs(P_star).max())) {
        diffuse   = false;
        d_diffuse = static_cast<int>(t) + 1;  // 1-based, as in R
      }
    }
    // Steady-state lock bookkeeping (stationary, post-diffuse tail only). When
    // the start-of-period covariance stops changing, the gain sequence captured
    // this period is frozen and the lock engages from the next period.
    if (ss_lock && !diffuse) {
      if (ss_capture) {
        ss_engaged = true;            // gains captured this period; lock onward
        ss_capture = false;
      } else if (have_prev &&
                 arma::abs(P_star - P_prev).max() <
                   conv_tol * std::max(1.0, arma::abs(P_star).max())) {
        // RELATIVE tolerance: a near-unit-root model has a large stationary
        // covariance, so an absolute test never converges (cf. the diffuse
        // collapse check above and the solve_lyapunov relative-tol fix).
        ss_capture = true;            // capture the gain sequence next period
      }
      P_prev    = P_star;             // start-of-next-period covariance
      have_prev = true;
    }
    if (return_filtered) filtered.col(t) = a.head(n_state);
  }
  if (ok && diffuse) diffuse_failed = true;  // sample ended inside the phase

  return List::create(_["loglik"]         = loglik,
                      _["a"]              = a,
                      // Final start-of-next-period covariance of the AUGMENTED
                      // state [s_T; eps_{T+1}]; the caller slices the state
                      // block out of it (kalman_filter's `final_cov`).
                      _["P"]              = P_star,
                      _["filtered"]       = return_filtered
                                              ? Rcpp::wrap(filtered)
                                              : R_NilValue,
                      _["ok"]             = ok,
                      _["d_diffuse"]      = d_diffuse,
                      _["diffuse_failed"] = diffuse_failed,
                      _["ll_contrib"]     = ll_contrib,
                      _["n_skipped"]      = n_skipped);
}
