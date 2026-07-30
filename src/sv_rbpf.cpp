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
// ---------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double sv_rbpf_loglik_cpp(const arma::mat& Y,        // n_obs x T
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

      arma::mat Fc;
      if (!arma::chol(Fc, Ft)) {               // non-PD: weight 0, carry state
        ll_t(i) = neg_inf;
        s_new.col(i) = s.col(i);
        P_new.slice(i) = Pi;
        continue;
      }
      arma::mat Fi = arma::inv_sympd(Ft);
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
      arma::mat TmKZ = TT - K * ZZ;
      arma::mat RmKD = RR - K * DD;
      arma::mat Pn = TmKZ * Pi * TmKZ.t() + RmKD * Se * RmKD.t();
      P_new.slice(i) = 0.5 * (Pn + Pn.t());
    }

    // -- log-mean-exp increment ----------------------------------------------
    double m = ll_t.max();
    if (!std::isfinite(m)) return neg_inf;
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

  return loglik;
}
