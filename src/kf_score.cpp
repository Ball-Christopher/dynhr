// kf_score.cpp -- C++ Kalman log-likelihood SCORE recursion for dynhr.
//
// Ports the hot per-step loop of R/analytic-gradient.R::.kf_loglik_score_sigma()
// to RcppArmadillo. Differentiates the exact per-step Kalman filter w.r.t. a set
// of parameters whose only effect is through the shock covariance Sigma_e
// (shock-std parameters, by first-order certainty equivalence). The decision
// rule -- TT, RR, ZZ, DD -- is held fixed.
//
// The cheap, one-time setup (QQ/HH/SS and the Lyapunov-initialised P0, dP0_k)
// is done in R and passed in, so this routine needs no Lyapunov solver. The R
// fallback (.kf_loglik_score_sigma with use_cpp = FALSE) is bit-for-bit
// equivalent; see test-kf-score-parity.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// [[Rcpp::export]]
List kf_score_sigma_cpp(const arma::mat& Yd,        // n_obs x T (already Y - d)
                        const arma::mat& TT,
                        const arma::mat& RR,
                        const arma::mat& ZZ,
                        const arma::mat& DD,
                        const arma::mat& Sigma_e,
                        const arma::mat& HH,
                        const arma::mat& SS,
                        const arma::mat& me_diag,
                        const arma::mat& P0,
                        const List& dSigma,         // K matrices (n_exo x n_exo)
                        const List& dH,             // K matrices (n_obs x n_obs)
                        const List& dS,             // K matrices (n_state x n_obs)
                        const List& dP0,            // K matrices (n_state x n_state)
                        double ss_tol = 1e-11) {
  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_T     = Yd.n_cols;
  const arma::uword K       = dSigma.size();
  const double ll_const = -0.5 * static_cast<double>(n_obs) * std::log(2.0 * M_PI);

  // Convert derivative-input lists to std::vector<arma::mat>; allocate running
  // derivative state (ds_k, dP_k).
  std::vector<arma::mat> dSig(K), dHl(K), dSl(K), dP(K);
  std::vector<arma::vec> ds(K);
  for (arma::uword k = 0; k < K; ++k) {
    dSig[k] = Rcpp::as<arma::mat>(dSigma[k]);
    dHl[k]  = Rcpp::as<arma::mat>(dH[k]);
    dSl[k]  = Rcpp::as<arma::mat>(dS[k]);
    dP[k]   = Rcpp::as<arma::mat>(dP0[k]);
    ds[k]   = arma::zeros<arma::vec>(n_state);
  }

  arma::vec dll = arma::zeros<arma::vec>(K);
  arma::mat P = P0;
  arma::vec s = arma::zeros<arma::vec>(n_state);
  double loglik = 0.0;
  bool ok = true;

  // Steady-state lock. Once P and every dP_k stop moving, the innovation
  // covariance F, gain K, and their derivatives dF_k, dK_k are constant, so the
  // O(n^3) Riccati / derivative-Riccati updates can be frozen and the remaining
  // steps run as a cheap constant-gain O(n^2) recursion (only ds_k and the
  // data-driven likelihood terms vary). This is the same trick the numerical
  // Kalman filter uses, and it is what the per-step derivative recursion needs
  // to be competitive.
  bool ss = false;
  arma::mat K_ss, Fi_ss;
  double logdetF_ss = 0.0;
  std::vector<arma::mat> dFss(K), dKss(K);
  arma::vec trconst = arma::zeros<arma::vec>(K);   // trace(Fi_ss dF_k) (constant)

  for (arma::uword t = 0; t < n_T; ++t) {
    arma::vec v = Yd.col(t) - ZZ * s;

    if (ss) {
      // ---- Constant-gain tail ------------------------------------------
      arma::vec Fiv = Fi_ss * v;
      double llt = ll_const - 0.5 * (logdetF_ss + arma::dot(v, Fiv));
      if (!std::isfinite(llt)) { ok = false; break; }
      loglik += llt;
      for (arma::uword k = 0; k < K; ++k) {
        arma::vec dv = -(ZZ * ds[k]);
        dll[k] += -0.5 * (trconst[k]
                          + 2.0 * arma::dot(dv, Fiv)
                          - arma::dot(Fiv, dFss[k] * Fiv));
        ds[k] = TT * ds[k] + dKss[k] * v + K_ss * dv;
      }
      s = TT * s + K_ss * v;
      continue;
    }

    // ---- Full transient step --------------------------------------------
    arma::mat PZ = P * ZZ.t();
    arma::mat F  = ZZ * PZ + HH + me_diag;
    F = 0.5 * (F + F.t());

    arma::mat Fc;
    if (!arma::chol(Fc, F)) { ok = false; break; }   // not PD -> infeasible
    double logdetF = 2.0 * arma::accu(arma::log(Fc.diag()));
    arma::mat Fi = arma::inv_sympd(F);

    arma::vec Fiv = Fi * v;
    double llt = ll_const - 0.5 * (logdetF + arma::dot(v, Fiv));
    if (!std::isfinite(llt)) { ok = false; break; }
    loglik += llt;

    arma::mat Kg   = (TT * PZ + SS) * Fi;
    arma::mat TmKZ = TT - Kg * ZZ;
    arma::mat RmKD = RR - Kg * DD;

    std::vector<arma::mat> dF_t(K), dK_t(K);
    double dP_drift = 0.0;
    for (arma::uword k = 0; k < K; ++k) {
      arma::vec dv  = -(ZZ * ds[k]);
      arma::mat dPZ = dP[k] * ZZ.t();
      arma::mat dF  = ZZ * dPZ + dHl[k];

      dll[k] += -0.5 * (arma::trace(Fi * dF)
                        + 2.0 * arma::dot(dv, Fiv)
                        - arma::dot(Fiv, dF * Fiv));

      arma::mat dKg   = (TT * dPZ + dSl[k]) * Fi - Kg * (dF * Fi);
      arma::mat dTmKZ = -(dKg * ZZ);
      arma::mat dRmKD = -(dKg * DD);

      ds[k] = TT * ds[k] + dKg * v + Kg * dv;

      arma::mat dPn = dTmKZ * P * TmKZ.t()
                    + TmKZ * dP[k] * TmKZ.t()
                    + TmKZ * P * dTmKZ.t()
                    + dRmKD * Sigma_e * RmKD.t()
                    + RmKD * dSig[k] * RmKD.t()
                    + RmKD * Sigma_e * dRmKD.t();
      dPn = 0.5 * (dPn + dPn.t());
      dP_drift = std::max(dP_drift, arma::abs(dPn - dP[k]).max());
      dP[k]  = dPn;
      dF_t[k] = dF; dK_t[k] = dKg;
    }

    s = TT * s + Kg * v;
    arma::mat Pn = TmKZ * P * TmKZ.t() + RmKD * Sigma_e * RmKD.t();
    Pn = 0.5 * (Pn + Pn.t());
    double P_drift = arma::abs(Pn - P).max();
    P = Pn;

    // Lock once both P and all dP_k have converged (after >=2 steps).
    if (t >= 2 && P_drift < ss_tol && dP_drift < ss_tol) {
      ss = true;
      K_ss = Kg; Fi_ss = Fi; logdetF_ss = logdetF;
      for (arma::uword k = 0; k < K; ++k) {
        dFss[k] = dF_t[k]; dKss[k] = dK_t[k];
        trconst[k] = arma::trace(Fi_ss * dFss[k]);
      }
    }
  }

  return List::create(_["loglik"] = loglik,
                      _["score"]  = dll,
                      _["ok"]     = ok,
                      _["ss_locked"] = ss);
}
