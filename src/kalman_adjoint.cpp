// kalman_adjoint.cpp -- C++ adjoint (reverse-mode) Kalman-filter
// log-likelihood + gradient recursion for dynhr.
//
// Ports R/gradient-adjoint-kf.R::.kf_loglik_adjoint() to RcppArmadillo.
// The R implementation is the spec and is independently validated to 1e-13
// against the tangent filter and finite differences; this file follows it
// line-for-line.  See R/gradient-adjoint-kf.R for a full derivation of every
// adjoint rule.
//
// State-space convention (identical to kf_tangent.cpp):
//   QQ = RR Sigma_e RR',  HH = DD Sigma_e DD',  SS = RR Sigma_e DD'
//   v_t = y_t - d - ZZ s_{t-1}
//   F_t = sym(ZZ P_{t-1} ZZ' + HH_t + me_diag_t)
//   K_t = (TT P_{t-1} ZZ' + SS_t) F_t^{-1}
//   A_t = TT - K_t ZZ,   B_t = RR - K_t DD
//   s_t = TT s_{t-1} + K_t v_t
//   P_t = sym(A_t P_{t-1} A_t' + B_t Se_t B_t')
//   P_0 = solve_lyapunov(TT, QQ),   s_0 = 0
//
// Time-varying inputs (new trailing args):
//   shock_scale_mat  -- n_exo x n_T  (pass all-ones for constant case)
//   me_extra_mat     -- n_obs x n_T  (pass all-zeros for constant case)
// P0 Lyapunov solve always uses baseline Sigma_e (tv is a sample-period
// phenomenon; stationary initialisation uses the baseline model).
//
// Numerical-failure behaviour mirrors the R reference exactly: Cholesky
// failure on F, non-finite ll_t, ll_t < ll_min, or non-finite Lyapunov
// solve all set ok = false and the function returns immediately.
//
// Contract: kf_adjoint_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
//             dTT_cube, dRR_cube, dZZ_cube, dDD_cube, dd_mat, dSigma_cube,
//             me_variance, ll_min, shock_scale_mat, me_extra_mat, return_bars)
//   -> list(ok, loglik, grad[, bars])
// With return_bars = true the raw adjoint (bar) matrices wrt the state-space
// system are appended as bars = list(G_TT, G_RR, G_ZZ, G_DD, g_d, G_Sig) --
// the inputs .solution_adjoint() (Tier 18 A2) contracts against the analytic
// primitive derivatives. Failure returns carry no bars (same as the R kernel).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

namespace {

inline arma::mat sym(const arma::mat& X) {
  return 0.5 * (X + X.t());
}

// Doubling recovery for the two stationary Lyapunov solves below (X = A X A'
// + B). The kron vec-solve's rcond gate fires on HIGHLY NON-NORMAL stable A
// (e.g. Reiter-HANK transition matrices), where rcond(I - A (x) A) underflows
// machine eps while the Lyapunov equation itself is well-posed. Doubling
// (X_{k+1} = X_k + A_k X_k A_k', A_{k+1} = A_k^2) converges for any spectral
// radius < 1 and any symmetric (possibly indefinite) B; for unit/explosive
// roots A_k fails to decay and X overflows to non-finite, so it returns false
// rather than a spurious solution -- the caller's fail contract is preserved.
// Mirrors R solve_lyapunov()'s RELATIVE convergence tolerance (an absolute
// tolerance never converges for near-unit-root systems whose X is huge).
bool lyap_doubling(const arma::mat& A, const arma::mat& B, arma::mat& X) {
  X = B;
  arma::mat Apow = A;
  for (int it = 0; it < 200; ++it) {
    arma::mat Xn = X + Apow * X * Apow.t();
    if (!Xn.is_finite()) return false;
    const double diff  = arma::abs(Xn - X).max();
    const double scale = std::max(1.0, arma::abs(Xn).max());
    if (diff < 1e-14 * scale) { X = Xn; return true; }
    Apow = Apow * Apow;
    if (!Apow.is_finite()) return false;
    X = Xn;
  }
  return false;
}

}  // namespace

// [[Rcpp::export]]
List kf_adjoint_cpp(const arma::mat& Y,
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
                    const arma::mat& me_extra_mat,
                    const bool return_bars = false) {

  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_exo   = RR.n_cols;
  const arma::uword n_T     = Y.n_cols;
  const arma::uword n_par   = dTT_cube.n_slices;

  const arma::mat tZZ     = ZZ.t();
  const arma::mat QQ      = sym(RR * Sigma_e * RR.t());
  const arma::mat HH      = sym(DD * Sigma_e * DD.t());
  const arma::mat SS      = RR * Sigma_e * DD.t();
  const arma::mat me_diag = me_variance * arma::eye(n_obs, n_obs);
  const double ll_const   = -0.5 * static_cast<double>(n_obs) * std::log(2.0 * M_PI);

  // Pre-compute flags for tv inputs.
  const bool has_shock_scale = !arma::all(arma::vectorise(
      arma::abs(shock_scale_mat - 1.0) < 1e-15));
  const bool has_me_extra = !arma::all(arma::vectorise(
      arma::abs(me_extra_mat) < 1e-15));

  arma::vec grad_na = arma::vec(n_par, arma::fill::value(NA_REAL));

  // -----------------------------------------------------------------------
  // Stationary initialisation: P0 = solve_lyapunov(TT, QQ)
  // Always uses baseline Sigma_e (tv scaling is a sample-period phenomenon).
  // Mirrors kf_tangent.cpp: n_state==1 scalar branch; general kron solve.
  // -----------------------------------------------------------------------
  arma::mat P0(n_state, n_state, arma::fill::zeros);
  bool lyap_ok = true;

  if (n_state == 1) {
    double denom = 1.0 - TT(0, 0) * TT(0, 0);
    P0(0, 0) = QQ(0, 0) / denom;
    if (!std::isfinite(P0(0, 0))) lyap_ok = false;
  } else {
    const arma::uword n2 = n_state * n_state;
    arma::mat M = arma::eye(n2, n2) - arma::kron(TT, TT);
    double rc = arma::rcond(M);
    bool kron_ok = (rc >= std::numeric_limits<double>::epsilon());
    if (kron_ok) {
      arma::vec p0_vec;
      kron_ok = arma::solve(p0_vec, M, arma::vectorise(QQ),
                            arma::solve_opts::no_approx);
      if (kron_ok) {
        P0 = arma::reshape(p0_vec, n_state, n_state);
        kron_ok = P0.is_finite();
      }
    }
    // Non-normal-TT recovery: see lyap_doubling() above.
    if (!kron_ok) lyap_ok = lyap_doubling(TT, QQ, P0) && P0.is_finite();
  }

  if (!lyap_ok) {
    return List::create(_["loglik"] = R_NegInf,
                        _["grad"]   = grad_na,
                        _["ok"]     = false);
  }

  // -----------------------------------------------------------------------
  // Forward pass: run the filter, storing per-step quantities.
  // Stored at step t: s_{t-1}, P_{t-1}, v_t, Fi_t, K_t, A_t, B_t.
  // Using arma::cube for P/Fi/K/A/B and arma::mat for s/v stores.
  // For tv: sc_t is recomputed from shock_scale_mat in the backward sweep
  // (no Se cube stored -- recompute from shock_scale_mat.col(t)).
  // -----------------------------------------------------------------------
  // s_store: n_state x n_T  (column t holds s entering step t)
  // v_store: n_obs   x n_T
  arma::mat s_store(n_state, n_T, arma::fill::zeros);
  arma::mat v_store(n_obs,   n_T, arma::fill::zeros);

  // P_store, Fi_store: n_state x n_state x n_T  (or n_obs x n_obs for Fi)
  arma::cube P_store(n_state, n_state, n_T, arma::fill::zeros);
  arma::cube Fi_store(n_obs,  n_obs,   n_T, arma::fill::zeros);
  // K: n_state x n_obs x n_T
  arma::cube K_store(n_state, n_obs,   n_T, arma::fill::zeros);
  // A: n_state x n_state x n_T
  arma::cube A_store(n_state, n_state, n_T, arma::fill::zeros);
  // B: n_state x n_exo x n_T
  arma::cube B_store(n_state, n_exo,   n_T, arma::fill::zeros);

  arma::vec  s = arma::zeros<arma::vec>(n_state);
  arma::mat  P = P0;
  double loglik = 0.0;
  bool ok = true;

  for (arma::uword t = 0; t < n_T; ++t) {
    s_store.col(t)   = s;
    P_store.slice(t) = P;

    // Per-period tv substitutions (forward pass).
    arma::mat Se_t_f, HH_t, SS_t, me_diag_t;
    if (has_shock_scale) {
      const arma::vec sc_t_f = shock_scale_mat.col(t);
      Se_t_f = Sigma_e % (sc_t_f * sc_t_f.t());
      HH_t   = sym(DD * Se_t_f * DD.t());
      SS_t   = RR * Se_t_f * DD.t();
    } else {
      Se_t_f = Sigma_e;
      HH_t   = HH;
      SS_t   = SS;
    }
    if (has_me_extra) {
      me_diag_t = me_diag + arma::diagmat(me_extra_mat.col(t));
    } else {
      me_diag_t = me_diag;
    }

    arma::mat PZ = P * tZZ;
    arma::mat Ft = sym(ZZ * PZ + HH_t + me_diag_t);

    arma::mat Fc;
    if (!arma::chol(Fc, Ft)) { ok = false; break; }
    arma::mat Fi  = arma::inv_sympd(Ft);
    double ldf    = 2.0 * arma::accu(arma::log(Fc.diag()));

    arma::vec v   = Y.col(t) - d - ZZ * s;
    arma::vec Fiv = Fi * v;

    double ll_t = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
    if (!std::isfinite(ll_t) || ll_t < ll_min) { ok = false; break; }
    loglik += ll_t;

    arma::mat K = (TT * PZ + SS_t) * Fi;
    arma::mat A = TT - K * ZZ;   // A_t = TT - K ZZ
    arma::mat B = RR - K * DD;   // B_t = RR - K DD

    Fi_store.slice(t) = Fi;
    v_store.col(t)    = v;
    K_store.slice(t)  = K;
    A_store.slice(t)  = A;
    B_store.slice(t)  = B;

    s = TT * s + K * v;
    arma::mat P_raw = A * P * A.t() + B * Se_t_f * B.t();
    // Joseph true-noise term for me_extra (mirrors kalman_filter's standard
    // path): P' += K diag(me_extra[, t]) K'. me_variance stays F-only
    // (regularizer convention; no Joseph term).
    if (has_me_extra)
      P_raw += (K * arma::diagmat(me_extra_mat.col(t))) * K.t();
    P = sym(P_raw);
  }

  if (!ok) {
    return List::create(_["loglik"] = R_NegInf,
                        _["grad"]   = grad_na,
                        _["ok"]     = false);
  }

  // -----------------------------------------------------------------------
  // Backward sweep.
  // bar_s and bar_P carry the adjoint of (s_t, P_t) produced by step t.
  // Both initialise to zero (loglik does not depend on (s_T, P_T)).
  // -----------------------------------------------------------------------
  arma::vec bar_s = arma::zeros<arma::vec>(n_state);
  arma::mat bar_P = arma::zeros<arma::mat>(n_state, n_state);

  arma::mat G_TT  = arma::zeros<arma::mat>(n_state, n_state);
  arma::mat G_RR  = arma::zeros<arma::mat>(n_state, n_exo);
  arma::mat G_ZZ  = arma::zeros<arma::mat>(n_obs,   n_state);
  arma::mat G_DD  = arma::zeros<arma::mat>(n_obs,   n_exo);
  arma::vec g_d   = arma::zeros<arma::vec>(n_obs);
  arma::mat G_Sig = arma::zeros<arma::mat>(n_exo,   n_exo);

  // Iterate t from n_T-1 down to 0 (0-based, matching the stored slices).
  for (arma::sword t_signed = static_cast<arma::sword>(n_T) - 1;
       t_signed >= 0; --t_signed) {
    const arma::uword t = static_cast<arma::uword>(t_signed);

    const arma::vec s_prev = s_store.col(t);
    const arma::mat P_prev = P_store.slice(t);
    const arma::vec v      = v_store.col(t);
    const arma::mat Fi     = Fi_store.slice(t);
    const arma::mat K      = K_store.slice(t);
    const arma::mat A      = A_store.slice(t);   // TT - K ZZ
    const arma::mat B      = B_store.slice(t);   // RR - K DD

    // Per-period tv substitutions for backward sweep.
    // Recompute Se_t from shock_scale_mat.col(t) (no cube stored).
    arma::mat Se_t;
    arma::vec sc_t;
    if (has_shock_scale) {
      sc_t  = shock_scale_mat.col(t);
      Se_t  = Sigma_e % (sc_t * sc_t.t());
    } else {
      Se_t  = Sigma_e;
    }

    arma::vec Fiv = Fi * v;

    // ---- Step 1: adjoint of P_t = sym(A P_prev A' + B Se_t B') -----------
    // Adjoint through sym(): bar_P <- sym(bar_P).
    bar_P = sym(bar_P);

    // bar_A = 2 bar_P A P_prev
    arma::mat bar_A = 2.0 * bar_P * A * P_prev;                // [n x n]

    // bar_B = 2 bar_P B Se_t  (tv: Se_t replaces Sigma_e)
    arma::mat bar_B = 2.0 * bar_P * B * Se_t;                  // [n x p]

    // bar_P_prev from A P A': += A' bar_P A
    arma::mat bar_P_prev_from_AP = A.t() * bar_P * A;          // [n x n]

    // G_Sig from B Se_t B':
    //   d/dSig_e tr(bar_P B Se_t B') where Se_t = Sigma_e % outer(sc_t,sc_t)
    //   => outer(sc_t,sc_t) .* (B' bar_P B) when has_shock_scale
    //   => B' bar_P B when constant
    if (has_shock_scale) {
      G_Sig += (sc_t * sc_t.t()) % (B.t() * bar_P * B);       // [p x p]
    } else {
      G_Sig += B.t() * bar_P * B;                              // [p x p]
    }

    // ---- Step 2: adjoint of s_t = TT s_prev + K v -------------------------
    // G_TT += outer(bar_s, s_prev)
    G_TT += bar_s * s_prev.t();                                // [n x n]

    // bar_K_from_s: += outer(bar_s, v)
    arma::mat bar_K = bar_s * v.t();                           // [n x q]

    // Adjoint of the me_extra Joseph term in P_t (P_t += K me_x K', with
    // me_x = diag(me_extra_mat.col(t)) being DATA, not differentiated):
    // with bar_P symmetrized in Step 1, d tr(bar_P K me_x K') / dK
    // = 2 bar_P K me_x. Mirrors R/gradient-adjoint-kf.R.
    if (has_me_extra)
      bar_K += 2.0 * bar_P * K * arma::diagmat(me_extra_mat.col(t));

    // bar_v_from_s: K' bar_s
    arma::vec bar_v = K.t() * bar_s;                           // [q]

    // bar_s_prev from TT s_prev: TT' bar_s
    arma::vec bar_s_prev = TT.t() * bar_s;                     // [n]

    // ---- Step 3: adjoint of ll_t = -0.5*(log|F| + v' Fi v) ---------------
    // bar_F += -0.5*(Fi - outer(Fiv, Fiv))
    arma::mat bar_F = -0.5 * (Fi - Fiv * Fiv.t());            // [q x q]

    // bar_v += -Fiv
    bar_v -= Fiv;                                              // [q]

    // ---- Step 4: adjoint of A = TT - K ZZ ---------------------------------
    // G_TT += bar_A
    G_TT += bar_A;

    // G_ZZ -= K' bar_A
    G_ZZ -= K.t() * bar_A;                                    // [q x n]

    // bar_K from A: -= bar_A ZZ'
    bar_K -= bar_A * tZZ;                                      // [n x q]

    // ---- Step 5: adjoint of B = RR - K DD ---------------------------------
    // G_RR += bar_B
    G_RR += bar_B;

    // G_DD -= K' bar_B
    G_DD -= K.t() * bar_B;                                    // [q x p]

    // bar_K from B: -= bar_B DD'
    bar_K -= bar_B * DD.t();                                   // [n x q]

    // ---- Step 6: adjoint of K = Mnum Fi (Mnum = TT P_prev ZZ' + SS_t) ----
    // bar_F += -(K' bar_K Fi)'  [cyclic: -K' bar_K Fi because Fi symmetric]
    bar_F -= K.t() * bar_K * Fi;                               // [q x q]

    // Symmetrize bar_F (F came from sym())
    bar_F = sym(bar_F);

    // bar_Mnum = bar_K Fi
    arma::mat bar_Mnum = bar_K * Fi;                           // [n x q]

    // From Mnum = TT P_prev ZZ' + SS_t:
    // G_TT += bar_Mnum ZZ P_prev
    G_TT += bar_Mnum * ZZ * P_prev;                           // [n x n]
    // G_ZZ += bar_Mnum' TT P_prev
    G_ZZ += bar_Mnum.t() * TT * P_prev;                       // [q x n]
    // bar_P_prev from Mnum: += TT' bar_Mnum ZZ
    arma::mat bar_P_prev_from_Mnum = TT.t() * bar_Mnum * ZZ; // [n x n]

    // From SS_t = RR Se_t DD' in Mnum (bar_SS = bar_Mnum):
    //   G_RR += bar_Mnum DD Se_t
    //   G_DD += bar_Mnum' RR Se_t
    //   G_Sig: outer(sc_t,sc_t) .* (RR' bar_Mnum DD) when tv; else RR' bar_Mnum DD
    G_RR  += bar_Mnum * DD * Se_t;                            // [n x p]
    G_DD  += bar_Mnum.t() * RR * Se_t;                        // [q x p]
    if (has_shock_scale) {
      G_Sig += (sc_t * sc_t.t()) % (RR.t() * bar_Mnum * DD); // [p x p]
    } else {
      G_Sig += RR.t() * bar_Mnum * DD;                        // [p x p]
    }

    // ---- Step 7: adjoint of F = sym(ZZ P_prev ZZ' + HH_t + me_diag_t) ----
    // bar_F already symmetrized above.
    // bar_P_prev from ZZ P ZZ': += ZZ' bar_F ZZ
    arma::mat bar_P_prev_from_F = tZZ * bar_F * ZZ;           // [n x n]

    // G_ZZ from ZZ P ZZ': += 2 bar_F ZZ P_prev
    G_ZZ += 2.0 * bar_F * ZZ * P_prev;                        // [q x n]

    // From HH_t = DD Se_t DD':
    //   G_DD += 2 bar_F DD Se_t
    //   G_Sig: outer(sc_t,sc_t) .* (DD' bar_F DD) when tv; else DD' bar_F DD
    G_DD  += 2.0 * bar_F * DD * Se_t;                         // [q x p]
    if (has_shock_scale) {
      G_Sig += (sc_t * sc_t.t()) % (DD.t() * bar_F * DD);    // [p x p]
    } else {
      G_Sig += DD.t() * bar_F * DD;                           // [p x p]
    }

    // ---- Step 8: adjoint of v_t = y_t - d - ZZ s_prev --------------------
    // g_d -= bar_v
    g_d -= bar_v;

    // G_ZZ -= outer(bar_v, s_prev)
    G_ZZ -= bar_v * s_prev.t();                                // [q x n]

    // bar_s_prev -= ZZ' bar_v
    bar_s_prev -= tZZ * bar_v;                                 // [n]

    // ---- Step 9: collect bar_P_{t-1} --------------------------------------
    arma::mat bar_P_prev_new = sym(bar_P_prev_from_AP +
                                   bar_P_prev_from_Mnum +
                                   bar_P_prev_from_F);

    // ---- Update carry variables for next (earlier) step ------------------
    bar_s = bar_s_prev;
    bar_P = bar_P_prev_new;
  }

  // -----------------------------------------------------------------------
  // Adjoint through P_0 = solve_lyapunov(TT, QQ)
  //
  // bar_P after the full backward loop is bar_P_0.
  // bar_QQ = solve_lyapunov(TT', bar_P_0)  [general]; scalar branch for n==1.
  // -----------------------------------------------------------------------
  arma::mat bar_P0 = sym(bar_P);
  arma::mat bar_QQ(n_state, n_state, arma::fill::zeros);

  if (n_state == 1) {
    // Scalar: P_0 = QQ / (1 - TT^2)
    double denom = 1.0 - TT(0, 0) * TT(0, 0);
    bar_QQ(0, 0) = bar_P0(0, 0) / denom;
    // G_TT from P_0(TT): bar_P_0 * 2 TT P_0 / (1 - TT^2)
    G_TT(0, 0) += bar_P0(0, 0) * 2.0 * TT(0, 0) * P0(0, 0) / denom;
  } else {
    // General: bar_QQ = solve_lyapunov(TT', bar_P_0)
    // i.e. solve (I - kron(TT', TT')) vec(bar_QQ) = vec(bar_P_0)
    const arma::uword n2 = n_state * n_state;
    arma::mat tTT = TT.t();
    arma::mat M   = arma::eye(n2, n2) - arma::kron(tTT, tTT);

    double rc = arma::rcond(M);
    bool kron_ok = (rc >= std::numeric_limits<double>::epsilon());
    if (kron_ok) {
      arma::vec bq_vec;
      kron_ok = arma::solve(bq_vec, M, arma::vectorise(bar_P0),
                            arma::solve_opts::no_approx);
      if (kron_ok) {
        bar_QQ = arma::reshape(bq_vec, n_state, n_state);
        kron_ok = bar_QQ.is_finite();
      }
    }
    // Non-normal-TT recovery: see lyap_doubling() above (doubling handles
    // the indefinite symmetric bar_P0 RHS; spectral radius of TT' == TT).
    if (!kron_ok &&
        (!lyap_doubling(tTT, bar_P0, bar_QQ) || !bar_QQ.is_finite())) {
      return List::create(_["loglik"] = R_NegInf,
                          _["grad"]   = grad_na,
                          _["ok"]     = false);
    }
    if (!bar_QQ.is_finite()) {
      return List::create(_["loglik"] = R_NegInf,
                          _["grad"]   = grad_na,
                          _["ok"]     = false);
    }
    bar_QQ = sym(bar_QQ);

    // G_TT from P_0's dependence on TT (via LHS of Lyapunov):
    // G_TT += 2 bar_QQ TT P_0
    G_TT += 2.0 * bar_QQ * TT * P0;                           // [n x n]
  }

  // G_Sig from QQ = RR Sigma_e RR' (baseline Sigma_e -- no sc_t factor here):
  // += RR' bar_QQ RR
  G_Sig += RR.t() * bar_QQ * RR;                              // [p x p]

  // G_RR from QQ = RR Sigma_e RR': += 2 bar_QQ RR Sigma_e
  G_RR += 2.0 * bar_QQ * RR * Sigma_e;                        // [n x p]

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

  if (return_bars) {
    return List::create(_["loglik"] = loglik,
                        _["grad"]   = grad,
                        _["ok"]     = true,
                        _["bars"]   = List::create(
                            _["G_TT"]  = G_TT,
                            _["G_RR"]  = G_RR,
                            _["G_ZZ"]  = G_ZZ,
                            _["G_DD"]  = G_DD,
                            _["g_d"]   = g_d,
                            _["G_Sig"] = G_Sig));
  }
  return List::create(_["loglik"] = loglik,
                      _["grad"]   = grad,
                      _["ok"]     = true);
}
