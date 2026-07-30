// kalman_adjoint_uni.cpp -- C++ port of the missing-data univariate adjoint
// log-likelihood + gradient recursion for dynhr.
//
// Ports R/gradient-adjoint-uni.R::.kf_loglik_adjoint_uni() to RcppArmadillo.
// The R implementation is the spec, independently validated to 1e-8 vs
// numDeriv (missing data) and 1e-10 (no-NA vs dense). This file follows it
// line-for-line through Steps 1-9 + the Lyapunov adjoint.
//
// Convention (identical to R reference):
//   QQ = RR Sigma_e RR',  HHo = DDo Sigma_e DDo',  SSo = RR Sigma_e DDo'
//   v_t = y[O,t] - d[O] - ZZo s_{t-1}
//   F_t = sym(ZZo P_{t-1} ZZo' + HHo + me * I_q)
//   K_t = (TT P_{t-1} ZZo' + SSo) F_t^{-1}
//   A_t = TT - K_t ZZo,   B_t = RR - K_t DDo
//   s_t = TT s_{t-1} + K_t v_t
//   P_t = sym(A_t P_{t-1} A_t' + B_t Sigma_e B_t')
//   P_0 = solve_lyapunov(TT, QQ)  [or supplied]
//
// Per-period observed-row subsetting:
//   O_t = arma::find_finite(Y.col(t));  empty O => pure-prediction step.
// Scatter G_ZZ / G_DD / g_d contributions back to rows O_t only.
//
// p0_supplied: when true, P0 is treated as a fixed parameter-independent input
// (no Lyapunov adjoint). When false, the Lyapunov adjoint is computed after
// the backward sweep (mirrors R reference's p0_supplied logic exactly).
//
// Contract:
//   kf_adjoint_uni_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
//                      dTT_cube, dRR_cube, dZZ_cube, dDD_cube, dd_mat,
//                      dSigma_cube, me_variance, ll_min, P0, p0_supplied)
//   -> list(ok, loglik, grad)

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
List kf_adjoint_uni_cpp(const arma::mat& Y,
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
                        const arma::mat& P0_in,
                        bool p0_supplied) {

  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_exo   = RR.n_cols;
  const arma::uword n_T     = Y.n_cols;
  const arma::uword n_par   = dTT_cube.n_slices;

  const double ll_2pi = std::log(2.0 * M_PI);

  arma::vec grad_na = arma::vec(n_par, arma::fill::value(NA_REAL));

  // -----------------------------------------------------------------------
  // Initialisation: P0 from Lyapunov or caller-supplied.
  // Mirrors R reference: p0_supplied branch checks all(is.finite(P0)).
  // -----------------------------------------------------------------------
  const arma::mat QQ = sym(RR * Sigma_e * RR.t());
  arma::mat P0(n_state, n_state, arma::fill::zeros);
  bool lyap_ok = true;

  if (p0_supplied) {
    if (!P0_in.is_finite()) {
      return List::create(_["loglik"] = R_NegInf,
                          _["grad"]   = grad_na,
                          _["ok"]     = false);
    }
    P0 = P0_in;
  } else {
    // Lyapunov: P0 = solve(I - kron(TT,TT), vec(QQ))
    if (n_state == 1) {
      double denom = 1.0 - TT(0, 0) * TT(0, 0);
      P0(0, 0) = QQ(0, 0) / denom;
      if (!std::isfinite(P0(0, 0))) lyap_ok = false;
    } else {
      const arma::uword n2 = n_state * n_state;
      arma::mat M = arma::eye(n2, n2) - arma::kron(TT, TT);
      double rc = arma::rcond(M);
      if (!(rc >= std::numeric_limits<double>::epsilon())) {
        lyap_ok = false;
      } else {
        arma::vec p0_vec;
        bool ok1 = arma::solve(p0_vec, M, arma::vectorise(QQ),
                               arma::solve_opts::no_approx);
        if (!ok1) {
          lyap_ok = false;
        } else {
          P0 = arma::reshape(p0_vec, n_state, n_state);
          if (!P0.is_finite()) lyap_ok = false;
        }
      }
    }
    if (!lyap_ok) {
      return List::create(_["loglik"] = R_NegInf,
                          _["grad"]   = grad_na,
                          _["ok"]     = false);
    }
  }

  // -----------------------------------------------------------------------
  // Forward pass: run the filter with per-period observed-row subsetting.
  // Store s_prev, P_prev, v, Fi, K, A, B per period (variable q per period).
  // We use std::vector<arma::mat> since q varies.
  // -----------------------------------------------------------------------
  // Per-period stores (indexed 0..n_T-1):
  std::vector<arma::vec>  s_store(n_T);
  std::vector<arma::mat>  P_store(n_T);
  std::vector<arma::uvec> O_store(n_T);
  // These are only filled when q > 0 (observation step):
  std::vector<arma::vec>  v_store(n_T);
  std::vector<arma::mat>  Fi_store(n_T);
  std::vector<arma::mat>  K_store(n_T);
  std::vector<arma::mat>  A_store(n_T);
  std::vector<arma::mat>  B_store(n_T);

  arma::vec s = arma::zeros<arma::vec>(n_state);
  arma::mat P = P0;
  double loglik = 0.0;
  bool ok = true;

  for (arma::uword t = 0; t < n_T; ++t) {
    s_store[t] = s;
    P_store[t] = P;

    // Observed rows: indices where Y(i,t) is finite.
    arma::uvec O = arma::find_finite(Y.col(t));
    O_store[t] = O;
    const arma::uword q = O.n_elem;

    if (q == 0) {
      // Pure prediction step (K = 0): A = TT, B = RR.
      // s_t = TT s_{t-1},  P_t = sym(TT P TT' + QQ)
      s = TT * s;
      P = sym(TT * P * TT.t() + QQ);
      continue;
    }

    // Per-period observed subsets.
    const arma::mat ZZo = ZZ.rows(O);        // q x n_state
    const arma::mat DDo = DD.rows(O);        // q x n_exo
    const arma::vec do_ = d.elem(O);         // q
    const arma::mat HHo = sym(DDo * Sigma_e * DDo.t());  // q x q
    const arma::mat SSo = RR * Sigma_e * DDo.t();        // n_state x q
    const arma::mat me_o = me_variance * arma::eye(q, q);

    arma::mat PZ = P * ZZo.t();                         // n_state x q
    arma::mat Ft = sym(ZZo * PZ + HHo + me_o);         // q x q

    arma::mat Fc;
    if (!arma::chol(Fc, Ft)) { ok = false; break; }
    // Non-throwing form: chol() success does not imply inv_sympd() success
    // (see kalman_adjoint.cpp) -- degrade gracefully instead of throwing.
    arma::mat Fi;
    if (!arma::inv_sympd(Fi, Ft)) { ok = false; break; }
    double ldf    = 2.0 * arma::accu(arma::log(Fc.diag()));

    arma::vec yt  = arma::vec(Y.col(t));
    arma::vec v   = yt.elem(O) - do_ - ZZo * s;         // q
    arma::vec Fiv = Fi * v;                              // q

    double ll_t = -0.5 * (static_cast<double>(q) * ll_2pi + ldf + arma::dot(v, Fiv));
    if (!std::isfinite(ll_t) || ll_t < ll_min) { ok = false; break; }
    loglik += ll_t;

    arma::mat K = (TT * PZ + SSo) * Fi;    // n_state x q
    arma::mat A = TT - K * ZZo;            // n_state x n_state
    arma::mat B = RR - K * DDo;            // n_state x n_exo

    Fi_store[t] = Fi;
    v_store[t]  = v;
    K_store[t]  = K;
    A_store[t]  = A;
    B_store[t]  = B;

    s = TT * s + K * v;
    P = sym(A * P * A.t() + B * Sigma_e * B.t());
  }

  if (!ok) {
    return List::create(_["loglik"] = R_NegInf,
                        _["grad"]   = grad_na,
                        _["ok"]     = false);
  }

  // -----------------------------------------------------------------------
  // Backward sweep.
  // bar_s and bar_P carry the adjoint of (s_t, P_t) produced by step t.
  // Both initialise to zero.
  // G_ZZ, G_DD, g_d have full n_obs rows; contributions scatter to O_t rows.
  // -----------------------------------------------------------------------
  arma::vec bar_s = arma::zeros<arma::vec>(n_state);
  arma::mat bar_P = arma::zeros<arma::mat>(n_state, n_state);

  arma::mat G_TT  = arma::zeros<arma::mat>(n_state, n_state);
  arma::mat G_RR  = arma::zeros<arma::mat>(n_state, n_exo);
  arma::mat G_ZZ  = arma::zeros<arma::mat>(n_obs,   n_state);
  arma::mat G_DD  = arma::zeros<arma::mat>(n_obs,   n_exo);
  arma::vec g_d   = arma::zeros<arma::vec>(n_obs);
  arma::mat G_Sig = arma::zeros<arma::mat>(n_exo,   n_exo);

  for (arma::sword t_signed = static_cast<arma::sword>(n_T) - 1;
       t_signed >= 0; --t_signed) {
    const arma::uword t = static_cast<arma::uword>(t_signed);

    const arma::vec   s_prev = s_store[t];
    const arma::mat   P_prev = P_store[t];
    const arma::uvec& O      = O_store[t];
    const arma::uword q      = O.n_elem;

    bar_P = sym(bar_P);

    if (q == 0) {
      // Pure-prediction adjoint: A = TT, B = RR.
      // Mirrors R reference lines 128-138:
      //   bar_A = 2 bar_P TT P_prev  (bar_P is bar_P_t, P_prev is P_{t-1})
      //   bar_B = 2 bar_P RR Sigma_e
      //   bar_P_prev = TT' bar_P TT
      //   G_Sig += RR' bar_P RR
      //   G_TT  += outer(bar_s, s_prev)
      //   bar_s_prev = TT' bar_s
      //   G_TT  += bar_A
      //   G_RR  += bar_B
      arma::mat bar_A    = 2.0 * bar_P * TT * P_prev;
      arma::mat bar_B    = 2.0 * bar_P * RR * Sigma_e;
      arma::mat bar_P_prev = TT.t() * bar_P * TT;
      G_Sig += RR.t() * bar_P * RR;
      G_TT  += bar_s * s_prev.t();
      arma::vec bar_s_prev = TT.t() * bar_s;
      G_TT  += bar_A;
      G_RR  += bar_B;
      bar_s = bar_s_prev;
      bar_P = sym(bar_P_prev);
      continue;
    }

    // Retrieve stored quantities.
    const arma::vec& v  = v_store[t];
    const arma::mat& Fi = Fi_store[t];
    const arma::mat& K  = K_store[t];
    const arma::mat& A  = A_store[t];
    const arma::mat& B  = B_store[t];
    const arma::vec  Fiv = Fi * v;

    // Observed subsets of ZZ, DD (needed in backward sweep too).
    const arma::mat ZZo = ZZ.rows(O);   // q x n_state
    const arma::mat DDo = DD.rows(O);   // q x n_exo

    // ---- Step 1: P_t = sym(A P_prev A' + B Sigma_e B') --------------------
    arma::mat bar_A           = 2.0 * bar_P * A * P_prev;
    arma::mat bar_B           = 2.0 * bar_P * B * Sigma_e;
    arma::mat bar_P_prev_AP   = A.t() * bar_P * A;
    G_Sig += B.t() * bar_P * B;

    // ---- Step 2: s_t = TT s_prev + K v ------------------------------------
    G_TT  += bar_s * s_prev.t();
    arma::mat bar_K     = bar_s * v.t();
    arma::vec bar_v     = K.t() * bar_s;
    arma::vec bar_s_prev = TT.t() * bar_s;

    // ---- Step 3: ll_t = -0.5*(q*log(2pi) + log|F| + v' Fi v) -------------
    arma::mat bar_F = -0.5 * (Fi - Fiv * Fiv.t());
    bar_v -= Fiv;

    // ---- Step 4: A = TT - K ZZo -------------------------------------------
    G_TT          += bar_A;
    G_ZZ.rows(O)  -= K.t() * bar_A;          // scatter to rows O
    bar_K         -= bar_A * ZZo.t();

    // ---- Step 5: B = RR - K DDo -------------------------------------------
    G_RR          += bar_B;
    G_DD.rows(O)  -= K.t() * bar_B;          // scatter to rows O
    bar_K         -= bar_B * DDo.t();

    // ---- Step 6: K = (TT P_prev ZZo' + SSo) Fi ---------------------------
    bar_F -= K.t() * bar_K * Fi;
    bar_F  = sym(bar_F);
    arma::mat bar_Mnum = bar_K * Fi;

    // From Mnum = TT P_prev ZZo':
    G_TT          += bar_Mnum * ZZo * P_prev;
    G_ZZ.rows(O)  += bar_Mnum.t() * TT * P_prev;    // scatter to rows O
    arma::mat bar_P_prev_Mnum = TT.t() * bar_Mnum * ZZo;

    // From SSo = RR Sigma_e DDo':
    G_RR          += bar_Mnum * DDo * Sigma_e;
    G_DD.rows(O)  += bar_Mnum.t() * RR * Sigma_e;   // scatter to rows O
    G_Sig         += RR.t() * bar_Mnum * DDo;

    // ---- Step 7: F = sym(ZZo P_prev ZZo' + HHo + me_o) -------------------
    // bar_F already symmetrized above.
    arma::mat bar_P_prev_F = ZZo.t() * bar_F * ZZo;
    G_ZZ.rows(O)  += 2.0 * bar_F * ZZo * P_prev;   // scatter to rows O
    G_DD.rows(O)  += 2.0 * bar_F * DDo * Sigma_e;   // scatter to rows O
    G_Sig         += DDo.t() * bar_F * DDo;

    // ---- Step 8: v = y[O,t] - d[O] - ZZo s_prev --------------------------
    g_d.elem(O)   -= bar_v;                          // scatter to rows O
    G_ZZ.rows(O)  -= bar_v * s_prev.t();             // scatter to rows O
    bar_s_prev    -= ZZo.t() * bar_v;

    // ---- Step 9: collect --------------------------------------------------
    bar_s = bar_s_prev;
    bar_P = sym(bar_P_prev_AP + bar_P_prev_Mnum + bar_P_prev_F);
  }

  // -----------------------------------------------------------------------
  // Lyapunov adjoint (only when P0 was computed internally, not supplied).
  // Mirrors R reference lines 200-213.
  // -----------------------------------------------------------------------
  arma::mat bar_P0 = sym(bar_P);

  if (!p0_supplied) {
    if (n_state == 1) {
      double denom = 1.0 - TT(0, 0) * TT(0, 0);
      double bar_QQ_11 = bar_P0(0, 0) / denom;
      G_TT(0, 0) += bar_P0(0, 0) * 2.0 * TT(0, 0) * P0(0, 0) / denom;
      G_Sig       += RR.t() * arma::mat(1, 1, arma::fill::value(bar_QQ_11)) * RR;
      G_RR        += 2.0 * arma::mat(1, 1, arma::fill::value(bar_QQ_11)) * RR * Sigma_e;
    } else {
      // bar_QQ = solve_lyapunov(TT', bar_P0)
      // i.e. solve (I - kron(TT', TT')) vec(bar_QQ) = vec(bar_P0)
      const arma::uword n2 = n_state * n_state;
      arma::mat tTT = TT.t();
      arma::mat M   = arma::eye(n2, n2) - arma::kron(tTT, tTT);

      double rc = arma::rcond(M);
      if (!(rc >= std::numeric_limits<double>::epsilon())) {
        return List::create(_["loglik"] = R_NegInf,
                            _["grad"]   = grad_na,
                            _["ok"]     = false);
      }

      arma::vec bq_vec;
      bool ok2 = arma::solve(bq_vec, M, arma::vectorise(bar_P0),
                             arma::solve_opts::no_approx);
      if (!ok2) {
        return List::create(_["loglik"] = R_NegInf,
                            _["grad"]   = grad_na,
                            _["ok"]     = false);
      }

      arma::mat bar_QQ = arma::reshape(bq_vec, n_state, n_state);
      if (!bar_QQ.is_finite()) {
        return List::create(_["loglik"] = R_NegInf,
                            _["grad"]   = grad_na,
                            _["ok"]     = false);
      }
      bar_QQ = sym(bar_QQ);

      // G_TT from P_0's dependence on TT: += 2 bar_QQ TT P0
      G_TT  += 2.0 * bar_QQ * TT * P0;
      G_Sig += RR.t() * bar_QQ * RR;
      G_RR  += 2.0 * bar_QQ * RR * Sigma_e;
    }
  }

  // -----------------------------------------------------------------------
  // Final contraction: grad[j] = <G_TT, dTT_j> + ... (Frobenius dot products)
  // -----------------------------------------------------------------------
  arma::vec grad = arma::zeros<arma::vec>(n_par);

  for (arma::uword j = 0; j < n_par; ++j) {
    double gj = 0.0;
    gj += arma::accu(G_TT  % dTT_cube.slice(j));
    gj += arma::accu(G_RR  % dRR_cube.slice(j));
    gj += arma::accu(G_ZZ  % dZZ_cube.slice(j));
    gj += arma::accu(G_DD  % dDD_cube.slice(j));
    gj += arma::dot(g_d,    dd_mat.col(j));
    gj += arma::accu(G_Sig % dSigma_cube.slice(j));
    grad(j) = gj;
  }

  return List::create(_["loglik"] = loglik,
                      _["grad"]   = grad,
                      _["ok"]     = true);
}
