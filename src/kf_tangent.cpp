// kf_tangent.cpp -- C++ tangent (forward-sensitivity) Kalman-filter
// log-likelihood + gradient recursion for dynhr.
//
// Ports R/gradient-tangent-kf.R::.kf_loglik_tangent() to RcppArmadillo.
// Mirrors the base filter (.kf_step, Joseph-form update, (X+X')/2
// symmetrisations) together with the per-parameter tangent recursions for
// dQQ/dHH/dSS, the Lyapunov dP0, and the per-step dv/dF/dll/dK/ds/dP. The R
// implementation in R/gradient-tangent-kf.R is the spec; this file follows
// it line-for-line. The R path remains the reference fallback (see
// .HAS_RCPP_KF_TANGENT() and the dispatch shim at the top of
// .kf_loglik_tangent).
//
// Per-parameter derivative stacks are passed as arma::cube (n_par slices;
// a zero slice where the R side had a NULL dXX entry); dd is passed as an
// n_obs x n_par matrix.
//
// Numerical-failure behaviour mirrors the R reference exactly: a Cholesky
// failure on F, a non-finite per-period log-likelihood, or ll_t < ll_min
// sets ok = FALSE and the recursion stops immediately (the R wrapper then
// returns loglik = -Inf, grad = NA). A non-finite Lyapunov solve (P0 or any
// dP0) also sets ok = FALSE before the main loop starts.
//
// Time-varying inputs:
//   shock_scale_mat  -- n_exo x n_T  (pass all-ones for constant case)
//   me_extra_mat     -- n_obs x n_T  (pass all-zeros for constant case)
// When shock_scale_mat is all-ones: Se_t = Sigma_e % (1*1') = Sigma_e (IEEE exact).
// When me_extra_mat is all-zeros: me_diag_t = me_diag + 0 = me_diag (exact).
// P0 / dP0 Lyapunov solves always use baseline Sigma_e (tv is a sample-period
// phenomenon; stationary initialisation uses the baseline model).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

namespace {

inline arma::mat sym(const arma::mat& X) {
  return 0.5 * (X + X.t());
}

}  // namespace

// [[Rcpp::export]]
List kf_tangent_cpp(const arma::mat& Y,
                    const arma::mat& TT,
                    const arma::mat& RR,
                    const arma::mat& ZZ,
                    const arma::mat& DD,
                    const arma::vec& d,
                    const arma::mat& Sigma_e,
                    const arma::cube& dTT_cube,
                    const arma::cube& dRR_cube,
                    const arma::cube& dZZ_cube,
                    const arma::cube& dDD_cube,
                    const arma::mat& dd_mat,
                    const arma::cube& dSigma_cube,
                    double me_variance,
                    double ll_min,
                    const arma::mat& shock_scale_mat,
                    const arma::mat& me_extra_mat) {
  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_T     = Y.n_cols;
  const arma::uword n_par   = dTT_cube.n_slices;

  const arma::mat tZZ = ZZ.t();
  const arma::mat QQ  = sym(RR * Sigma_e * RR.t());
  const arma::mat HH  = sym(DD * Sigma_e * DD.t());
  const arma::mat SS  = RR * Sigma_e * DD.t();
  const arma::mat me_diag = me_variance * arma::eye(n_obs, n_obs);
  const double ll_const = -0.5 * static_cast<double>(n_obs) * std::log(2.0 * M_PI);

  // Pre-compute flags for tv inputs.
  // has_shock_scale: true if any element of shock_scale_mat differs from 1.0
  // has_me_extra:    true if any element of me_extra_mat differs from 0.0
  const bool has_shock_scale = !arma::all(arma::vectorise(
      arma::abs(shock_scale_mat - 1.0) < 1e-15));
  const bool has_me_extra = !arma::all(arma::vectorise(
      arma::abs(me_extra_mat) < 1e-15));

  arma::vec grad = arma::vec(n_par, arma::fill::value(NA_REAL));

  // -- Lyapunov solve: P0 = solve_lyapunov(TT, QQ) and each dP0 -------------
  // .solve_lyapunov(A, Q): n==1 special case Q / (1 - A^2); otherwise solve
  // (I - kron(A,A)) vec(P) = vec(Q), with rcond(M) < eps -> NaN matrix.
  // P0 and every dP0 share the SAME operator (I - kron(TT,TT)); factor once
  // by stacking all RHS columns into a single arma::solve call.
  // IMPORTANT: always uses baseline Sigma_e (tv scaling is a sample-period
  // phenomenon; stationary init uses the baseline model).
  arma::mat P0(n_state, n_state, arma::fill::zeros);
  std::vector<arma::mat> dP0(n_par);
  bool lyap_ok = true;

  if (n_state == 1) {
    double a2 = TT(0, 0) * TT(0, 0);
    double denom = 1.0 - a2;
    P0(0, 0) = QQ(0, 0) / denom;
    if (!std::isfinite(P0(0, 0))) lyap_ok = false;
    for (arma::uword j = 0; j < n_par && lyap_ok; ++j) {
      const arma::mat& dTT = dTT_cube.slice(j);
      const arma::mat& dRR = dRR_cube.slice(j);
      const arma::mat& dSig = dSigma_cube.slice(j);
      arma::mat dQQ = sym(dRR * Sigma_e * RR.t() + RR * dSig * RR.t() +
                           RR * Sigma_e * dRR.t());
      arma::mat rhs = sym(dTT * P0 * TT.t() + TT * P0 * dTT.t() + dQQ);
      double dval = rhs(0, 0) / denom;
      if (!std::isfinite(dval)) { lyap_ok = false; break; }
      dP0[j] = arma::mat(1, 1);
      dP0[j](0, 0) = dval;
    }
  } else {
    const arma::uword n2 = n_state * n_state;
    arma::mat M = arma::eye(n2, n2) - arma::kron(TT, TT);

    double rc = arma::rcond(M);
    if (!(rc >= std::numeric_limits<double>::epsilon())) {
      lyap_ok = false;
    } else {
      // Build dQQ for every parameter first (need P0 itself for the dP0
      // RHS, so solve P0 first, then dP0's).
      arma::vec qq_vec = arma::vectorise(QQ);
      arma::vec p0_vec;
      bool solve1 = arma::solve(p0_vec, M, qq_vec, arma::solve_opts::no_approx);
      if (!solve1) {
        lyap_ok = false;
      } else {
        P0 = arma::reshape(p0_vec, n_state, n_state);
        if (!P0.is_finite()) lyap_ok = false;
      }

      if (lyap_ok && n_par > 0) {
        // Stack the n_par dP0 RHS vectors into one matrix and solve once.
        arma::mat rhs_mat(n2, n_par);
        for (arma::uword j = 0; j < n_par; ++j) {
          const arma::mat& dTT = dTT_cube.slice(j);
          const arma::mat& dRR = dRR_cube.slice(j);
          const arma::mat& dSig = dSigma_cube.slice(j);
          arma::mat dQQ = sym(dRR * Sigma_e * RR.t() + RR * dSig * RR.t() +
                               RR * Sigma_e * dRR.t());
          arma::mat rhs = sym(dTT * P0 * TT.t() + TT * P0 * dTT.t() + dQQ);
          rhs_mat.col(j) = arma::vectorise(rhs);
        }
        arma::mat dP0_mat;
        bool solve2 = arma::solve(dP0_mat, M, rhs_mat, arma::solve_opts::no_approx);
        if (!solve2) {
          lyap_ok = false;
        } else {
          for (arma::uword j = 0; j < n_par; ++j) {
            arma::mat dPj = arma::reshape(dP0_mat.col(j), n_state, n_state);
            if (!dPj.is_finite()) { lyap_ok = false; break; }
            dP0[j] = dPj;
          }
        }
      }
    }
  }

  if (!lyap_ok) {
    return List::create(_["loglik"] = R_NegInf,
                        _["grad"]   = grad,
                        _["ok"]     = false);
  }

  // -- Per-parameter dQQ (dP0 already computed above) -----------------------
  // dHH_l and dSS_l are now per-period when has_shock_scale; we keep a
  // constant-case copy here that is used only when !has_shock_scale (avoids
  // recomputing every period in the constant path).
  std::vector<arma::mat> dQQ_l(n_par), dHH_l(n_par), dSS_l(n_par);
  std::vector<arma::vec> ds_l(n_par);
  for (arma::uword j = 0; j < n_par; ++j) {
    const arma::mat& dRR = dRR_cube.slice(j);
    const arma::mat& dDD = dDD_cube.slice(j);
    const arma::mat& dSig = dSigma_cube.slice(j);

    dQQ_l[j] = sym(dRR * Sigma_e * RR.t() + RR * dSig * RR.t() +
                   RR * Sigma_e * dRR.t());
    dHH_l[j] = sym(dDD * Sigma_e * DD.t() + DD * dSig * DD.t() +
                   DD * Sigma_e * dDD.t());
    dSS_l[j] = dRR * Sigma_e * DD.t() + RR * dSig * DD.t() +
               RR * Sigma_e * dDD.t();

    ds_l[j] = arma::zeros<arma::vec>(n_state);
  }

  arma::vec s = arma::zeros<arma::vec>(n_state);
  arma::mat P = P0;
  double loglik = 0.0;
  arma::vec grad_acc = arma::zeros<arma::vec>(n_par);
  bool ok = true;

  for (arma::uword t = 0; t < n_T; ++t) {
    // -- Per-period tv substitutions (mirrors R reference) -------------------
    // When !has_shock_scale: Se_t = Sigma_e, HH_t = HH, SS_t = SS (constants).
    // When !has_me_extra:    me_diag_t = me_diag (constant).
    arma::mat Se_t, HH_t, SS_t, me_diag_t;
    if (has_shock_scale) {
      const arma::vec sc_t = shock_scale_mat.col(t);
      Se_t = Sigma_e % (sc_t * sc_t.t());
      HH_t = sym(DD * Se_t * DD.t());
      SS_t = RR * Se_t * DD.t();
    } else {
      Se_t = Sigma_e;
      HH_t = HH;
      SS_t = SS;
    }
    if (has_me_extra) {
      me_diag_t = me_diag + arma::diagmat(me_extra_mat.col(t));
    } else {
      me_diag_t = me_diag;
    }

    // -- Base step (mirrors .kf_step exactly) --------------------------------
    arma::mat PZ = P * tZZ;
    arma::mat Ft = sym(ZZ * PZ + HH_t + me_diag_t);

    arma::mat Fc;
    if (!arma::chol(Fc, Ft)) { ok = false; break; }
    arma::mat Fi = arma::inv_sympd(Ft);
    double ldf = 2.0 * arma::accu(arma::log(Fc.diag()));

    arma::vec v = Y.col(t) - d - ZZ * s;
    arma::vec Fiv = Fi * v;

    double ll_t = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
    if (!std::isfinite(ll_t) || ll_t < ll_min) { ok = false; break; }
    loglik += ll_t;

    arma::mat K = (TT * PZ + SS_t) * Fi;
    arma::mat A = TT - K * ZZ;          // TT - K ZZ
    arma::mat B = RR - K * DD;          // RR - K DD

    arma::vec s_n = TT * s + K * v;
    arma::mat P_n_raw = A * P * A.t() + B * Se_t * B.t();
    // Joseph true-noise term for me_extra (mirrors kalman_filter's standard
    // path): P' += K diag(me_extra[, t]) K'. me_variance stays F-only
    // (regularizer convention; no Joseph term). Kme = K diag(me_x_t) is
    // hoisted for the per-parameter tangent recursion below.
    arma::vec me_x_t;
    arma::mat Kme;
    if (has_me_extra) {
      me_x_t = me_extra_mat.col(t);
      Kme    = K * arma::diagmat(me_x_t);
      P_n_raw += Kme * K.t();
    }
    arma::mat P_n = sym(P_n_raw);

    // -- Hoisted parameter-independent pieces for the tangent recursion ----
    arma::mat TtP  = TT * P;
    arma::mat AP   = A * P;
    arma::mat BSig = B * Se_t;

    // Precompute sc_t for the j-loop (only needed when has_shock_scale).
    arma::vec sc_t_j;
    if (has_shock_scale) sc_t_j = shock_scale_mat.col(t);

    for (arma::uword j = 0; j < n_par; ++j) {
      const arma::mat& dTT = dTT_cube.slice(j);
      const arma::mat& dRR = dRR_cube.slice(j);
      const arma::mat& dZZ = dZZ_cube.slice(j);
      const arma::mat& dDD = dDD_cube.slice(j);
      const arma::mat& dSig = dSigma_cube.slice(j);
      arma::mat& dP = dP0[j];
      arma::vec& ds = ds_l[j];

      // Per-period dHH and dSS: if has_shock_scale, use Se_t; else use
      // the constant pre-computed values (bit-identical to old code path).
      arma::mat dHH_t, dSS_t;
      if (has_shock_scale) {
        // dSig_t = dSig % outer(sc_t, sc_t)  (chain rule through Se_t)
        arma::mat dSig_t = dSig % (sc_t_j * sc_t_j.t());
        dHH_t = sym(dDD * Se_t * DD.t() + DD * dSig_t * DD.t() + DD * Se_t * dDD.t());
        dSS_t = dRR * Se_t * DD.t() + RR * dSig_t * DD.t() + RR * Se_t * dDD.t();
      } else {
        dHH_t = dHH_l[j];  // constant-case: exact reference
        dSS_t = dSS_l[j];
      }

      arma::vec dd_j = dd_mat.col(j);
      arma::mat tdZZ = dZZ.t();

      // dv = -dd - dZZ s - ZZ ds
      arma::vec dv = -dd_j - dZZ * s - ZZ * ds;

      // dF = dZZ P ZZ' + ZZ dP ZZ' + ZZ P dZZ' + dHH_t   (symmetrise)
      arma::mat dF = sym(dZZ * P * ZZ.t() + ZZ * dP * ZZ.t() +
                         ZZ * P * tdZZ + dHH_t);

      // dll_t = -0.5 * ( tr(Fi dF) - (Fi v)' dF (Fi v) ) - (Fi v)' dv
      arma::vec Fi_dF_Fiv = dF * Fiv;
      double dll_t = -0.5 * (arma::accu(Fi % dF) - arma::dot(Fiv, Fi_dF_Fiv)) -
                      arma::dot(Fiv, dv);
      grad_acc[j] += dll_t;

      // dK = (dTT P ZZ' + TT dP ZZ' + TT P dZZ' + dSS_t) Fi - K dF Fi
      arma::mat dPZ = dP * tZZ;
      arma::mat dK_num = dTT * P * ZZ.t() + TT * dPZ + TtP * tdZZ + dSS_t;
      arma::mat dK = dK_num * Fi - K * (dF * Fi);

      // ds' = dTT s + TT ds + dK v + K dv
      arma::vec ds_n = dTT * s + TT * ds + dK * v + K * dv;

      // dA = dTT - dK ZZ - K dZZ;  dB = dRR - dK DD - K dDD
      arma::mat dA = dTT - dK * ZZ - K * dZZ;
      arma::mat dB = dRR - dK * DD - K * dDD;

      // dP' = dA P A' + A dP A' + A P dA' + dB Se_t B' + B dSig_t B' + B Se_t dB'
      // (when !has_shock_scale: Se_t = Sigma_e, dSig_t = dSig => same as before)
      arma::mat dSig_eff;
      if (has_shock_scale) {
        dSig_eff = dSig % (sc_t_j * sc_t_j.t());
      } else {
        dSig_eff = dSig;
      }
      arma::mat dP_n_raw = dA * P * A.t() + A * dP * A.t() + AP * dA.t() +
                           dB * Se_t * B.t() + B * dSig_eff * B.t() +
                           BSig * dB.t();
      // Tangent of the me_extra Joseph term P' += K me_x K' (me_extra is
      // data, not differentiated): dP' += dK me_x K' + K me_x dK'.
      if (has_me_extra)
        dP_n_raw += dK * arma::diagmat(me_x_t) * K.t() + Kme * dK.t();
      arma::mat dP_n = sym(dP_n_raw);

      ds = ds_n;
      dP = dP_n;
    }

    s = s_n;
    P = P_n;
  }

  if (!ok) {
    return List::create(_["loglik"] = R_NegInf,
                        _["grad"]   = grad,
                        _["ok"]     = false);
  }

  return List::create(_["loglik"] = loglik,
                      _["grad"]   = grad_acc,
                      _["ok"]     = true);
}
