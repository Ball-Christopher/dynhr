// kf_loglik_dG.cpp -- C++ port of .kf_loglik_dG() from
// R/hessian-adjoint-analytic.R.
//
// Computes the DIRECTIONAL DERIVATIVE of the six adjoint gradient matrices
// (dG_TT, dG_RR, dG_ZZ, dG_DD, dg_d, dG_Sig) of the Kalman-filter
// log-likelihood along a given direction dX = (dTT, dRR, dZZ, dDD, dd, dSig).
//
// Method: "forward-over-reverse" -- the forward pass stores both primal
// quantities AND their tangents; the backward sweep carries both the primal
// bar_* adjoints AND their tangents dbar_* in parallel, exactly mirroring
// the R reference line-for-line.
//
// Two Lyapunov solves via kron+solve:
//   (1) dP0:   (I - TT (x) TT) vec(dP0) = vec(dTT P0 TT' + TT P0 dTT' + dQQ)
//   (2) dbar_QQ: (I - TT' (x) TT') vec(dbar_QQ) = vec(rhs_dbar_QQ)
//
// Faithful scalar (n_state == 1) branch mirrors the R scalar branch exactly.
// Returns ok=false on non-finite Lyapunov (mirroring kalman_adjoint_uni.cpp).
// Stops with error if Y contains NAs.
//
// Exported as: kf_loglik_dG_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
//                               dTT, dRR, dZZ, dDD, dd_dir, dSig,
//                               me_variance, ll_min)
// Returns: List(dG_TT, dG_RR, dG_ZZ, dG_DD, dg_d, dG_Sig, ok)

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

namespace {

inline arma::mat sym(const arma::mat& X) {
  return 0.5 * (X + X.t());
}

// Solve the discrete Lyapunov equation  X - A X A' = Q  for X.
// For n_state==1: X = Q / (1 - a^2).
// For n_state>1:  (I - A (x) A) vec(X) = vec(Q), via kron+solve.
// Returns all-NaN on failure.
arma::mat solve_lyapunov(const arma::mat& A, const arma::mat& Q,
                         arma::uword n) {
  if (n == 1) {
    double denom = 1.0 - A(0, 0) * A(0, 0);
    arma::mat X(1, 1);
    X(0, 0) = Q(0, 0) / denom;
    return X;
  }
  const arma::uword n2 = n * n;
  arma::mat M = arma::eye(n2, n2) - arma::kron(A, A);
  double rc = arma::rcond(M);
  if (!(rc >= std::numeric_limits<double>::epsilon())) {
    return arma::mat(n, n, arma::fill::value(
                       std::numeric_limits<double>::quiet_NaN()));
  }
  arma::vec x_vec;
  bool ok = arma::solve(x_vec, M, arma::vectorise(Q),
                        arma::solve_opts::no_approx);
  if (!ok) {
    return arma::mat(n, n, arma::fill::value(
                       std::numeric_limits<double>::quiet_NaN()));
  }
  return arma::reshape(x_vec, n, n);
}

}  // namespace

// [[Rcpp::export]]
List kf_loglik_dG_cpp(const arma::mat& Y,
                      const arma::mat& TT,
                      const arma::mat& RR,
                      const arma::mat& ZZ,
                      const arma::mat& DD,
                      const arma::vec& d,
                      const arma::mat& Sigma_e,
                      const arma::mat& dTT,
                      const arma::mat& dRR,
                      const arma::mat& dZZ,
                      const arma::mat& dDD,
                      const arma::vec& dd_dir,
                      const arma::mat& dSig,
                      double me_variance,
                      double ll_min) {

  // Guard: Y must have no NAs.
  if (!Y.is_finite()) {
    Rcpp::stop("kf_loglik_dG_cpp: Y must not contain missing values.");
  }

  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_exo   = RR.n_cols;
  const arma::uword n_T     = Y.n_cols;

  const double ll_const = -0.5 * static_cast<double>(n_obs) * std::log(2.0 * M_PI);

  // Failure list (ok = false, all-NA outputs).
  auto fail = [&]() -> List {
    return List::create(
      _["dG_TT"]  = arma::mat(n_state, n_state, arma::fill::value(NA_REAL)),
      _["dG_RR"]  = arma::mat(n_state, n_exo,   arma::fill::value(NA_REAL)),
      _["dG_ZZ"]  = arma::mat(n_obs,   n_state, arma::fill::value(NA_REAL)),
      _["dG_DD"]  = arma::mat(n_obs,   n_exo,   arma::fill::value(NA_REAL)),
      _["dg_d"]   = arma::vec(n_obs,           arma::fill::value(NA_REAL)),
      _["dG_Sig"] = arma::mat(n_exo,   n_exo,   arma::fill::value(NA_REAL)),
      _["ok"]     = false
    );
  };

  // -------------------------------------------------------------------------
  // Pre-compute noise matrices and their tangents.
  // -------------------------------------------------------------------------
  const arma::mat tZZ  = ZZ.t();
  const arma::mat tdZZ = dZZ.t();

  const arma::mat QQ  = sym(RR * Sigma_e * RR.t());
  // HH = DD Sig DD',  SS = RR Sig DD'
  const arma::mat HH  = sym(DD * Sigma_e * DD.t());
  // SS used only in the original R; in the forward recursion we compute inline.

  // dQQ = dRR Sig RR' + RR dSig RR' + RR Sig dRR'   (symmetrised)
  const arma::mat dQQ = sym(dRR * Sigma_e * RR.t() +
                            RR  * dSig    * RR.t() +
                            RR  * Sigma_e * dRR.t());

  // dHH = dDD Sig DD' + DD dSig DD' + DD Sig dDD'   (symmetrised)
  const arma::mat dHH = sym(dDD * Sigma_e * DD.t() +
                            DD  * dSig    * DD.t() +
                            DD  * Sigma_e * dDD.t());

  // dSS = dRR Sig DD' + RR dSig DD' + RR Sig dDD'
  const arma::mat dSS = dRR * Sigma_e * DD.t() +
                        RR  * dSig    * DD.t() +
                        RR  * Sigma_e * dDD.t();

  const arma::mat me_diag = me_variance * arma::eye(n_obs, n_obs);

  // -------------------------------------------------------------------------
  // Lyapunov: P0 and dP0.
  // -------------------------------------------------------------------------
  arma::mat P0 = solve_lyapunov(TT, QQ, n_state);
  if (!P0.is_finite()) return fail();
  P0 = sym(P0);

  // rhs_dP0 = dTT P0 TT' + TT P0 dTT' + dQQ   (symmetrised)
  const arma::mat rhs_dP0 = sym(dTT * P0 * TT.t() + TT * P0 * dTT.t() + dQQ);
  arma::mat dP0 = solve_lyapunov(TT, rhs_dP0, n_state);
  if (!dP0.is_finite()) return fail();
  dP0 = sym(dP0);

  // -------------------------------------------------------------------------
  // Forward pass: store primal + tangent at each step.
  // -------------------------------------------------------------------------
  std::vector<arma::vec>  s_store(n_T),  ds_store(n_T);
  std::vector<arma::mat>  P_store(n_T),  dP_store(n_T);
  std::vector<arma::vec>  v_store(n_T),  dv_store(n_T);
  std::vector<arma::mat>  Fi_store(n_T), dFi_store(n_T);
  std::vector<arma::mat>  K_store(n_T),  dK_store(n_T);
  std::vector<arma::mat>  A_store(n_T),  dA_store(n_T);
  std::vector<arma::mat>  B_store(n_T),  dB_store(n_T);

  arma::vec s  = arma::zeros<arma::vec>(n_state);
  arma::mat P  = P0;
  arma::vec ds = arma::zeros<arma::vec>(n_state);
  arma::mat dP = dP0;

  for (arma::uword t = 0; t < n_T; ++t) {
    s_store[t]  = s;
    P_store[t]  = P;
    ds_store[t] = ds;
    dP_store[t] = dP;

    // Primal innovation covariance: Ft = sym(ZZ P ZZ' + HH + me_diag)
    const arma::mat PZ  = P * tZZ;    // n_state x n_obs
    const arma::mat Ft  = sym(ZZ * PZ + HH + me_diag);

    // Cholesky
    arma::mat Fc;
    if (!arma::chol(Fc, Ft)) return fail();
    const arma::mat Fi  = arma::inv_sympd(Ft);
    const double ldf    = 2.0 * arma::accu(arma::log(Fc.diag()));

    const arma::vec v   = Y.col(t) - d - ZZ * s;
    const arma::vec Fiv = Fi * v;

    const double ll_t = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
    if (!std::isfinite(ll_t) || ll_t < ll_min) return fail();

    // Gain K, A, B (primal)
    const arma::mat SS_t = RR * Sigma_e * DD.t();   // n_state x n_obs
    const arma::mat K    = (TT * PZ + SS_t) * Fi;
    const arma::mat A    = TT - K * ZZ;
    const arma::mat B    = RR - K * DD;

    // Tangent of v: dv = -dd_dir - dZZ s - ZZ ds
    const arma::vec dv_t = -dd_dir - dZZ * s - ZZ * ds;

    // Tangent of F: dF = sym(dZZ P ZZ' + ZZ dP ZZ' + ZZ P dZZ' + dHH)
    const arma::mat dF_t = sym(dZZ * PZ +
                               ZZ * dP * tZZ +
                               ZZ * P  * tdZZ +
                               dHH);

    // Tangent of Fi: dFi = -Fi dF Fi
    const arma::mat dFi_t = -(Fi * dF_t * Fi);

    // Tangent of K:
    //   dMnum = dTT P ZZ' + TT dP ZZ' + TT P dZZ' + dSS
    //   dK = dMnum Fi - K dF Fi
    const arma::mat dPZ   = dP * tZZ;
    const arma::mat dK_num = dTT * P  * tZZ +
                             TT  * dPZ +
                             TT  * P  * tdZZ +
                             dSS;
    const arma::mat dK_t  = dK_num * Fi - K * (dF_t * Fi);

    // Tangent of A: dA = dTT - dK ZZ - K dZZ
    const arma::mat dA_t = dTT - dK_t * ZZ - K * dZZ;

    // Tangent of B: dB = dRR - dK DD - K dDD
    const arma::mat dB_t = dRR - dK_t * DD - K * dDD;

    // Store
    v_store[t]   = v;
    Fi_store[t]  = Fi;
    K_store[t]   = K;
    A_store[t]   = A;
    B_store[t]   = B;
    dv_store[t]  = dv_t;
    dFi_store[t] = dFi_t;
    dK_store[t]  = dK_t;
    dA_store[t]  = dA_t;
    dB_store[t]  = dB_t;

    // Advance primal state
    const arma::vec s_new = TT * s + K * v;
    const arma::mat AP    = A * P;
    const arma::mat BSig  = B * Sigma_e;
    const arma::mat P_new = sym(AP * A.t() + BSig * B.t());

    // Advance tangent state:
    // ds' = dTT s + TT ds + dK v + K dv
    const arma::vec ds_new = dTT * s + TT * ds + dK_t * v + K * dv_t;

    // dP' = sym(dA P A' + A dP A' + A P dA' + dB Sig B' + B dSig B' + B Sig dB')
    const arma::mat dP_new = sym(dA_t * P    * A.t() +
                                 A    * dP   * A.t() +
                                 AP   * dA_t.t() +
                                 dB_t * Sigma_e * B.t() +
                                 B    * dSig    * B.t() +
                                 BSig * dB_t.t());

    s  = s_new;  P  = P_new;
    ds = ds_new; dP = dP_new;
  }

  // -------------------------------------------------------------------------
  // Backward sweep: primal adjoints (bar_*) + their tangents (dbar_*).
  // -------------------------------------------------------------------------
  arma::vec bar_s  = arma::zeros<arma::vec>(n_state);
  arma::mat bar_P  = arma::zeros<arma::mat>(n_state, n_state);
  arma::vec dbar_s = arma::zeros<arma::vec>(n_state);
  arma::mat dbar_P = arma::zeros<arma::mat>(n_state, n_state);

  arma::mat G_TT  = arma::zeros<arma::mat>(n_state, n_state);
  arma::mat G_RR  = arma::zeros<arma::mat>(n_state, n_exo);
  arma::mat G_ZZ  = arma::zeros<arma::mat>(n_obs,   n_state);
  arma::mat G_DD  = arma::zeros<arma::mat>(n_obs,   n_exo);
  arma::vec g_d   = arma::zeros<arma::vec>(n_obs);
  arma::mat G_Sig = arma::zeros<arma::mat>(n_exo,   n_exo);

  arma::mat dG_TT  = arma::zeros<arma::mat>(n_state, n_state);
  arma::mat dG_RR  = arma::zeros<arma::mat>(n_state, n_exo);
  arma::mat dG_ZZ  = arma::zeros<arma::mat>(n_obs,   n_state);
  arma::mat dG_DD  = arma::zeros<arma::mat>(n_obs,   n_exo);
  arma::vec dg_d   = arma::zeros<arma::vec>(n_obs);
  arma::mat dG_Sig = arma::zeros<arma::mat>(n_exo,   n_exo);

  for (arma::sword t_signed = static_cast<arma::sword>(n_T) - 1;
       t_signed >= 0; --t_signed) {
    const arma::uword t = static_cast<arma::uword>(t_signed);

    const arma::vec& s_prev  = s_store[t];
    const arma::mat& P_prev  = P_store[t];
    const arma::vec& ds_prev = ds_store[t];
    const arma::mat& dP_prev = dP_store[t];
    const arma::vec& v       = v_store[t];
    const arma::mat& Fi      = Fi_store[t];
    const arma::mat& K       = K_store[t];
    const arma::mat& A       = A_store[t];
    const arma::mat& B       = B_store[t];
    const arma::vec& dv      = dv_store[t];
    const arma::mat& dFi     = dFi_store[t];
    const arma::mat& dK      = dK_store[t];
    const arma::mat& dA      = dA_store[t];
    const arma::mat& dB      = dB_store[t];

    const arma::vec Fiv  = Fi * v;
    const arma::vec dFiv = dFi * v + Fi * dv;

    // ================================================================
    // Step 1: Adjoint of P_t = sym(A P_prev A' + B Sigma_e B')
    // ================================================================

    bar_P  = sym(bar_P);
    dbar_P = sym(dbar_P);

    // bar_A = 2 bar_P A P_prev
    const arma::mat bar_A  = 2.0 * bar_P * A * P_prev;
    // dbar_A = 2 (dbar_P A P_prev + bar_P dA P_prev + bar_P A dP_prev)
    const arma::mat dbar_A = 2.0 * (dbar_P * A  * P_prev +
                                    bar_P  * dA * P_prev +
                                    bar_P  * A  * dP_prev);

    // bar_B = 2 bar_P B Sigma_e
    const arma::mat bar_B  = 2.0 * bar_P * B * Sigma_e;
    // dbar_B = 2 (dbar_P B Sigma_e + bar_P dB Sigma_e + bar_P B dSig)
    const arma::mat dbar_B = 2.0 * (dbar_P * B  * Sigma_e +
                                    bar_P  * dB * Sigma_e +
                                    bar_P  * B  * dSig);

    // bar_P_prev_from_AP = A' bar_P A
    const arma::mat bar_P_prev_from_AP  = A.t() * bar_P  * A;
    // dbar_P_prev_from_AP = dA' bar_P A + A' dbar_P A + A' bar_P dA
    const arma::mat dbar_P_prev_from_AP = dA.t() * bar_P  * A  +
                                          A.t()  * dbar_P * A  +
                                          A.t()  * bar_P  * dA;

    // G_Sig += B' bar_P B
    G_Sig  += B.t()  * bar_P  * B;
    dG_Sig += dB.t() * bar_P  * B  +
              B.t()  * dbar_P * B  +
              B.t()  * bar_P  * dB;

    // ================================================================
    // Step 2: Adjoint of s_t = TT s_prev + K v
    // ================================================================

    // G_TT += outer(bar_s, s_prev)
    G_TT  += bar_s  * s_prev.t();
    dG_TT += dbar_s * s_prev.t() + bar_s * ds_prev.t();

    // bar_K = outer(bar_s, v)
    arma::mat bar_K  = bar_s  * v.t();
    arma::mat dbar_K = dbar_s * v.t() + bar_s * dv.t();

    // bar_v = K' bar_s
    arma::vec bar_v  = K.t()  * bar_s;
    arma::vec dbar_v = dK.t() * bar_s + K.t() * dbar_s;

    // bar_s_prev = TT' bar_s
    arma::vec bar_s_prev  = TT.t()  * bar_s;
    arma::vec dbar_s_prev = dTT.t() * bar_s + TT.t() * dbar_s;

    // ================================================================
    // Step 3: Adjoint of ll_t = -0.5*(log|F| + v' Fi v)
    // ================================================================

    // bar_F = -0.5 * (Fi - outer(Fiv, Fiv))
    arma::mat bar_F  = -0.5 * (Fi  - Fiv  * Fiv.t());
    // dbar_F = -0.5 * (dFi - outer(dFiv, Fiv) - outer(Fiv, dFiv))
    arma::mat dbar_F = -0.5 * (dFi - dFiv * Fiv.t() - Fiv * dFiv.t());

    // bar_v -= Fiv
    bar_v  -= Fiv;
    dbar_v -= dFiv;

    // ================================================================
    // Step 4: Adjoint of A = TT - K ZZ
    // ================================================================

    G_TT  += bar_A;
    dG_TT += dbar_A;

    G_ZZ  -= K.t()  * bar_A;
    dG_ZZ -= dK.t() * bar_A + K.t() * dbar_A;

    bar_K  -= bar_A  * tZZ;
    dbar_K -= dbar_A * tZZ + bar_A * tdZZ;

    // ================================================================
    // Step 5: Adjoint of B = RR - K DD
    // ================================================================

    G_RR  += bar_B;
    dG_RR += dbar_B;

    G_DD  -= K.t()  * bar_B;
    dG_DD -= dK.t() * bar_B + K.t() * dbar_B;

    bar_K  -= bar_B  * DD.t();
    dbar_K -= dbar_B * DD.t() + bar_B * dDD.t();

    // ================================================================
    // Step 6: Adjoint of K = Mnum Fi  (Mnum = TT P_prev ZZ' + SS)
    // ================================================================

    // bar_F -= K' bar_K Fi   (then symmetrize)
    bar_F  -= K.t()  * bar_K  * Fi;
    dbar_F -= dK.t() * bar_K  * Fi  +
              K.t()  * dbar_K * Fi  +
              K.t()  * bar_K  * dFi;

    bar_F  = sym(bar_F);
    dbar_F = sym(dbar_F);

    // bar_Mnum = bar_K Fi
    const arma::mat bar_Mnum  = bar_K  * Fi;
    const arma::mat dbar_Mnum = dbar_K * Fi + bar_K * dFi;

    // From Mnum = TT P_prev ZZ' in Mnum:
    G_TT  += bar_Mnum  * ZZ  * P_prev;
    dG_TT += dbar_Mnum * ZZ  * P_prev  +
             bar_Mnum  * dZZ * P_prev  +
             bar_Mnum  * ZZ  * dP_prev;

    G_ZZ  += bar_Mnum.t()  * TT  * P_prev;
    dG_ZZ += dbar_Mnum.t() * TT  * P_prev  +
             bar_Mnum.t()  * dTT * P_prev  +
             bar_Mnum.t()  * TT  * dP_prev;

    const arma::mat bar_P_prev_from_Mnum  = TT.t()  * bar_Mnum  * ZZ;
    const arma::mat dbar_P_prev_from_Mnum = dTT.t() * bar_Mnum  * ZZ  +
                                            TT.t()  * dbar_Mnum * ZZ  +
                                            TT.t()  * bar_Mnum  * dZZ;

    // From SS = RR Sig DD' in Mnum:
    G_RR  += bar_Mnum  * DD  * Sigma_e;
    dG_RR += dbar_Mnum * DD  * Sigma_e  +
             bar_Mnum  * dDD * Sigma_e  +
             bar_Mnum  * DD  * dSig;

    G_DD  += bar_Mnum.t()  * RR  * Sigma_e;
    dG_DD += dbar_Mnum.t() * RR  * Sigma_e  +
             bar_Mnum.t()  * dRR * Sigma_e  +
             bar_Mnum.t()  * RR  * dSig;

    G_Sig  += RR.t()  * bar_Mnum  * DD;
    dG_Sig += dRR.t() * bar_Mnum  * DD  +
              RR.t()  * dbar_Mnum * DD  +
              RR.t()  * bar_Mnum  * dDD;

    // ================================================================
    // Step 7: Adjoint of F = sym(ZZ P_prev ZZ' + HH + me_diag)
    // ================================================================

    // bar_F already symmetrized above.
    const arma::mat bar_P_prev_from_F  = tZZ  * bar_F  * ZZ;
    const arma::mat dbar_P_prev_from_F = tdZZ * bar_F  * ZZ  +
                                         tZZ  * dbar_F * ZZ  +
                                         tZZ  * bar_F  * dZZ;

    G_ZZ  += 2.0 * bar_F  * ZZ  * P_prev;
    dG_ZZ += 2.0 * (dbar_F * ZZ  * P_prev  +
                    bar_F  * dZZ * P_prev  +
                    bar_F  * ZZ  * dP_prev);

    // From HH = DD Sig DD':
    G_DD  += 2.0 * bar_F  * DD  * Sigma_e;
    dG_DD += 2.0 * (dbar_F * DD  * Sigma_e  +
                    bar_F  * dDD * Sigma_e  +
                    bar_F  * DD  * dSig);

    G_Sig  += DD.t()  * bar_F  * DD;
    dG_Sig += dDD.t() * bar_F  * DD  +
              DD.t()  * dbar_F * DD  +
              DD.t()  * bar_F  * dDD;

    // ================================================================
    // Step 8: Adjoint of v_t = y_t - d - ZZ s_prev
    // ================================================================

    g_d  -= bar_v;
    dg_d -= dbar_v;

    G_ZZ  -= bar_v  * s_prev.t();
    dG_ZZ -= dbar_v * s_prev.t() + bar_v * ds_prev.t();

    bar_s_prev  -= tZZ  * bar_v;
    dbar_s_prev -= tdZZ * bar_v + tZZ * dbar_v;

    // ================================================================
    // Step 9: collect bar_P_{t-1}
    // ================================================================
    bar_s  = bar_s_prev;
    bar_P  = sym(bar_P_prev_from_AP  + bar_P_prev_from_Mnum  + bar_P_prev_from_F);
    dbar_s = dbar_s_prev;
    dbar_P = sym(dbar_P_prev_from_AP + dbar_P_prev_from_Mnum + dbar_P_prev_from_F);
  }

  // -------------------------------------------------------------------------
  // Lyapunov adjoint: bar_QQ = solve_lyapunov(TT', bar_P0)
  // -------------------------------------------------------------------------
  const arma::mat bar_P0  = sym(bar_P);
  const arma::mat dbar_P0 = sym(dbar_P);

  arma::mat bar_QQ;
  arma::mat dbar_QQ;

  if (n_state == 1) {
    // Scalar branch: P0 = QQ / (1 - TT^2)
    const double tt11  = TT(0, 0);
    const double denom = 1.0 - tt11 * tt11;
    bar_QQ  = arma::mat(1, 1, arma::fill::value(bar_P0(0, 0) / denom));
    dbar_QQ = arma::mat(1, 1, arma::fill::value(dbar_P0(0, 0) / denom));

    // Primal G_TT[1,1] += bar_P0 * 2 TT P0 / denom
    G_TT(0, 0) += bar_P0(0, 0)  * 2.0 * tt11 * P0(0, 0) / denom;

    // Tangent: see R scalar branch derivation
    dG_TT(0, 0) += dbar_P0(0, 0) * 2.0 * tt11        * P0(0, 0) / denom +
                   bar_P0(0, 0)  * 2.0 * dTT(0, 0)   * P0(0, 0) / denom +
                   bar_P0(0, 0)  * 2.0 * tt11         * dP0(0, 0) / denom +
                   bar_P0(0, 0)  * 2.0 * tt11         * P0(0, 0) *
                     2.0 * tt11 * dTT(0, 0) / (denom * denom);
  } else {
    // General: bar_QQ = solve_lyapunov(TT', bar_P0)
    const arma::mat tTT = TT.t();
    bar_QQ = solve_lyapunov(tTT, bar_P0, n_state);
    if (!bar_QQ.is_finite()) return fail();
    bar_QQ = sym(bar_QQ);

    // dbar_QQ = solve_lyapunov(TT', rhs_dbar_QQ)
    // where rhs_dbar_QQ = sym(dbar_P0 + dTT' bar_QQ TT + TT' bar_QQ dTT)
    const arma::mat rhs_dbar_QQ = sym(dbar_P0 +
                                      dTT.t() * bar_QQ * TT +
                                      tTT     * bar_QQ * dTT);
    dbar_QQ = solve_lyapunov(tTT, rhs_dbar_QQ, n_state);
    if (!dbar_QQ.is_finite()) return fail();
    dbar_QQ = sym(dbar_QQ);

    // Primal: G_TT += 2 bar_QQ TT P0
    G_TT  += 2.0 * bar_QQ  * TT  * P0;
    dG_TT += 2.0 * (dbar_QQ * TT  * P0  +
                    bar_QQ  * dTT * P0  +
                    bar_QQ  * TT  * dP0);
  }

  // G_Sig += RR' bar_QQ RR
  G_Sig  += RR.t()  * bar_QQ  * RR;
  dG_Sig += dRR.t() * bar_QQ  * RR  +
            RR.t()  * dbar_QQ * RR  +
            RR.t()  * bar_QQ  * dRR;

  // G_RR += 2 bar_QQ RR Sigma_e
  G_RR  += 2.0 * bar_QQ  * RR  * Sigma_e;
  dG_RR += 2.0 * (dbar_QQ * RR  * Sigma_e  +
                  bar_QQ  * dRR * Sigma_e  +
                  bar_QQ  * RR  * dSig);

  return List::create(
    _["dG_TT"]  = dG_TT,
    _["dG_RR"]  = dG_RR,
    _["dG_ZZ"]  = dG_ZZ,
    _["dG_DD"]  = dG_DD,
    _["dg_d"]   = dg_d,
    _["dG_Sig"] = dG_Sig,
    _["ok"]     = true
  );
}
