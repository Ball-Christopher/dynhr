// src/sv_rbpf.cpp
// ---------------------------------------------------------------------------
// Compiled kernel for the SV-on-shocks Rao-Blackwellized particle filter
// (R reference: .sv_rbpf_loglik in R/sv-rbpf.R). One call runs the full
// T-period sweep: AR(1) log-volatility propagation, the exact per-particle
// Kalman step conditional on the particle's shock scale, log-mean-exp
// marginal-likelihood increments, and systematic resampling.
//
// RNG: uses R's RNG stream (R::rnorm / R::unif_rand) in EXACTLY the R
// reference's draw order (init h by particle-major column fill; per period
// the same fill; ONE uniform per period for the resampler), so a given
// set.seed() produces the same volatility particles as the R path and the
// two implementations agree to numerical precision (parity-tested in
// test-sv-rbpf.R).
//
// Per-particle failure (non-PD forecast covariance) => weight 0 for that
// particle; a fully infeasible cloud => -Inf, matching the R reference.
//
// MEASUREMENT ERROR (F4-A, 2026-09-04): me_diag is TRUE i.i.d. observation
// noise, not an F-only regulariser -- it enters Ft AND the Joseph term
// (Pn += K diag(me) K'), exactly as R/kf-step.R and kalman_filter() do.
//
// RETURN CONTRACT (B5, 2026-09-02): a LIST, not a bare double --
//   $loglik      scalar log-likelihood (-Inf on a fully infeasible cloud)
//   $fail_period 1-based period at which EVERY particle failed, or 0.
// The all-particles-failed diagnostic used to be an Rcpp::warning() raised
// from inside this kernel. Under options(warn = 2) R promotes a warning to an
// error, and R's error mechanism is a longjmp: it unwinds PAST every C++
// destructor on the stack, so any RNGScope/Armadillo cleanup between here and
// the .Call boundary is skipped (a stale .Random.seed, leaked buffers) --
// precisely on the runs a user asked to be strict. The kernel therefore only
// REPORTS the period; R/sv-rbpf.R raises the warning after the call returns.
// ---------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List sv_rbpf_loglik_cpp(const arma::mat& Y,        // n_obs x T
                          const arma::mat& TT,       // n_s x n_s
                          const arma::mat& ZZ,       // n_obs x n_s
                          const arma::mat& RR,       // n_s x n_e
                          const arma::mat& DD,       // n_obs x n_e
                          const arma::mat& Sigma_e,  // n_e x n_e
                          const arma::vec& d,        // n_obs (may be zeros)
                          const arma::mat& P0,       // n_s x n_s
                          const arma::uvec& sv_idx1, // 1-based positions
                          const arma::vec& mu,
                          const arma::vec& rho,
                          const arma::vec& seta,
                          const int n_particles,
                          const arma::vec& me_diag)  // n_obs or empty
{
  const int n_T   = Y.n_cols;
  const int n_obs = Y.n_rows;
  const int n_s   = TT.n_rows;
  const int n_e   = RR.n_cols;
  const int n_sv  = sv_idx1.n_elem;
  const int N     = n_particles;
  const double neg_inf = -std::numeric_limits<double>::infinity();
  const double ll_const = -0.5 * n_obs * std::log(2.0 * M_PI);
  // A few ulps of headroom for the Joseph cancellation snap below; the same
  // 8 * .Machine$double.eps R/kf-step.R uses.
  const double eps_c = 8.0 * std::numeric_limits<double>::epsilon();

  arma::uvec sv_idx = sv_idx1 - 1;             // 0-based

  const bool has_me = me_diag.n_elem > 0;
  arma::mat ME(n_obs, n_obs, arma::fill::zeros);
  if (has_me) ME.diag() = me_diag;

  // ---- init: h ~ stationary AR(1) marginal, (s, P) shared baseline --------
  arma::vec sd0 = seta / arma::sqrt(1.0 - arma::square(rho));
  arma::mat h(n_sv, N);
  for (int i = 0; i < N; ++i)                  // R column-major fill order
    for (int k = 0; k < n_sv; ++k)
      h(k, i) = R::rnorm(mu(k), sd0(k));

  arma::mat s(n_s, N, arma::fill::zeros);
  arma::cube P(n_s, n_s, N);
  for (int i = 0; i < N; ++i) P.slice(i) = P0;

  arma::vec scale_full(n_e, arma::fill::ones);
  arma::vec ll_t(N);
  arma::mat s_new(n_s, N);
  arma::cube P_new(n_s, n_s, N);
  arma::mat h_res(n_sv, N);
  arma::uvec idx(N);

  double loglik = 0.0;

  for (int t = 0; t < n_T; ++t) {
    // -- propagate volatility (bootstrap proposal), R draw order ------------
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < n_sv; ++k)
        h(k, i) = mu(k) + rho(k) * (h(k, i) - mu(k)) + R::rnorm(0.0, seta(k));

    const arma::vec y_t = Y.col(t);

    // -- exact conditional Kalman step per particle --------------------------
    for (int i = 0; i < N; ++i) {
      for (int k = 0; k < n_sv; ++k)
        scale_full(sv_idx(k)) = std::exp(h(k, i) / 2.0);

      arma::mat Se = Sigma_e % (scale_full * scale_full.t());
      arma::mat HH = DD * Se * DD.t();
      arma::mat SS = RR * Se * DD.t();

      const arma::mat& Pi = P.slice(i);
      arma::mat PZ = Pi * ZZ.t();
      arma::mat Ft = ZZ * PZ + HH;
      if (has_me) Ft += ME;
      Ft = 0.5 * (Ft + Ft.t());

      // NON-FINITE GUARD (F4-A follow-up). An extreme volatility particle can
      // overflow exp(h/2) to +Inf (the SBC's inv_gamma sigma_eta prior has a
      // heavy tail), and Se = Inf then produces NaN in HH = DD Se DD' wherever
      // DD has a structural zero (0 * Inf). Symmetrising does NOT remove a
      // NaN, and arma::is_symmetric() compares entries with exact !=, which is
      // ALWAYS true for NaN -- so arma::chol() printed
      //   "warning: chol(): given matrix is not symmetric"
      // to stderr before correctly returning false. The likelihood was never
      // wrong (the particle takes the weight-0 path either way); the warning
      // was pure stderr noise, unsuppressable from R, that surfaced ~5 times
      // in the R = 100 SBC certification. Fail the particle BEFORE chol sees
      // the matrix -- the same point at which R/kf-step.R's
      // tryCatch(chol(Ft)) fails ("the leading minor of order 1 is not
      // positive"), so R/C++ parity is unchanged.
      arma::mat Fc;
      if (!Ft.is_finite()) {                   // NaN/Inf: weight 0, carry state
        ll_t(i) = neg_inf;
        s_new.col(i) = s.col(i);
        P_new.slice(i) = Pi;
        continue;
      }
      if (!arma::chol(Fc, Ft)) {               // non-PD: weight 0, carry state
        ll_t(i) = neg_inf;
        s_new.col(i) = s.col(i);
        P_new.slice(i) = Pi;
        continue;
      }
      // Non-throwing 2-arg form: chol() succeeding does not guarantee
      // inv_sympd() succeeds (different LAPACK paths disagree right at the
      // PD boundary — the 0.9.0.0001 KF-kernel crash class). Rather than
      // degrade a particle we already have a valid Cholesky factor for,
      // recover Fi from Fc directly: arma::chol(Fc, Ft) returns Fc UPPER
      // triangular with Ft = Fc.t() * Fc, so Ft^{-1} = Rinv * Rinv.t()
      // where Rinv = inv(trimatu(Fc)) = Fc^{-1}. Only fall back to the
      // weight-0 path if THAT also fails (near-singular Fc).
      arma::mat Fi;
      if (!arma::inv_sympd(Fi, Ft)) {
        arma::mat Rinv;
        if (!arma::inv(Rinv, arma::trimatu(Fc))) {
          ll_t(i) = neg_inf;
          s_new.col(i) = s.col(i);
          P_new.slice(i) = Pi;
          continue;
        }
        Fi = Rinv * Rinv.t();
      }
      double ldf = 2.0 * arma::accu(arma::log(Fc.diag()));

      arma::vec v = y_t - ZZ * s.col(i) - d;
      arma::vec Fiv = Fi * v;
      double ll = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
      if (!std::isfinite(ll)) {
        ll_t(i) = neg_inf;
        s_new.col(i) = s.col(i);
        P_new.slice(i) = Pi;
        continue;
      }
      ll_t(i) = ll;

      arma::mat K = (TT * PZ + SS) * Fi;
      s_new.col(i) = TT * s.col(i) + K * v;
      arma::mat KZ = K * ZZ;
      arma::mat KD = K * DD;
      arma::mat TmKZ = TT - KZ;
      arma::mat RmKD = RR - KD;
      // CANCELLATION SNAP (see R/kf-step.R for the full argument). On a
      // perfectly-observed step K = I exactly, so both Joseph factors are
      // exactly zero and P' = 0. In floating point they come out at ~1 ulp of
      // the terms that cancelled and P' is rounding noise squared, which the
      // next period's F^-1 amplifies without limit. Ft and its Cholesky
      // factor are BIT-IDENTICAL to R's here (same LAPACK dpotrf); the one
      // ulp that separated the two implementations entered through
      // arma::inv_sympd() vs R's chol2inv(), became a 100%-relative
      // difference in TT - K ZZ, and ended as a factor 2-4 in a |ll| ~ 1e30
      // particle weight past the SV volatility overflow point. Snapping an
      // entry that is within a few ulps of the magnitudes that formed it
      // replaces noise with the value it approximates, and both
      // implementations snap to the same exact zero.
      for (arma::uword c = 0; c < TmKZ.n_cols; ++c)
        for (arma::uword r = 0; r < TmKZ.n_rows; ++r)
          if (std::abs(TmKZ(r, c)) <=
              eps_c * std::max(std::abs(TT(r, c)), std::abs(KZ(r, c))))
            TmKZ(r, c) = 0.0;
      for (arma::uword c = 0; c < RmKD.n_cols; ++c)
        for (arma::uword r = 0; r < RmKD.n_rows; ++r)
          if (std::abs(RmKD(r, c)) <=
              eps_c * std::max(std::abs(RR(r, c)), std::abs(KD(r, c))))
            RmKD(r, c) = 0.0;
      arma::mat Pn = TmKZ * Pi * TmKZ.t() + RmKD * Se * RmKD.t();
      // TRUE i.i.d. measurement-noise law (F4-A): y_t = ZZ x_t + DD e_t + u_t
      // with Var(u_t) = diag(me_diag) requires Pn += K diag(me_diag) K' for
      // ANY gain K. Before F4-A me_diag entered Ft only (an F-regulariser),
      // which made the RB-PF disagree with kalman_filter(method =
      // "univariate") -- and with R/kf-step.R -- by O(me_variance).
      if (has_me) Pn += K * ME * K.t();
      // Symmetrised on the way out (as R/kf-step.R:142 does), so the Pi that
      // feeds next period's Ft is exactly symmetric -- 0.5*(A + A.t()) is
      // bit-exact in IEEE-754 because addition is commutative. The Joseph and
      // K ME K' terms therefore cannot be the source of an asymmetric chol
      // argument; only a NaN can be, which the guard above catches.
      P_new.slice(i) = 0.5 * (Pn + Pn.t());
    }

    // -- log-mean-exp increment ----------------------------------------------
    double m = ll_t.max();
    if (!std::isfinite(m)) {
      // Every particle infeasible: report the period and stop. The R wrapper
      // turns fail_period > 0 into the warning this used to raise here.
      return Rcpp::List::create(Rcpp::_["loglik"]      = neg_inf,
                                Rcpp::_["fail_period"] = t + 1);
    }
    arma::vec w_un = arma::exp(ll_t - m);
    loglik += m + std::log(arma::mean(w_un));

    // -- systematic resample (ONE uniform, matching the R reference) --------
    double u0 = R::unif_rand();
    arma::vec cumw = arma::cumsum(w_un / arma::accu(w_un));
    cumw(N - 1) = 1.0;                         // guard fp drift
    int ii = 0;
    for (int j = 0; j < N; ++j) {
      double pos = (u0 + j) / N;
      while (cumw(ii) < pos) ++ii;
      idx(j) = ii;
    }

    for (int j = 0; j < N; ++j) {
      h_res.col(j) = h.col(idx(j));
      s.col(j) = s_new.col(idx(j));
      P.slice(j) = P_new.slice(idx(j));
    }
    h = h_res;
  }

  return Rcpp::List::create(Rcpp::_["loglik"]      = loglik,
                            Rcpp::_["fail_period"] = 0);
}
