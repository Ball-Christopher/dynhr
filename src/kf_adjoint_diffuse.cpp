// kf_adjoint_diffuse.cpp -- C++ port of the exact-diffuse adjoint Kalman
// filter gradient for dynhr.
//
// Ports R/gradient-adjoint-diffuse.R::.kf_loglik_adjoint_diffuse() to
// RcppArmadillo. The R implementation is the authoritative spec; this file
// follows it section-by-section.
//
// Two stages:
//   Stage 1: Forward diffuse Kalman filter (Case A / Case B / stationary tail)
//            + backward adjoint sweep accumulating G_TT, G_RR, G_ZZ, G_DD,
//            g_d, G_Sig.
//   Stage 2: Analytic init-adjoint through the Schur decomposition (Schur-swap
//            to order unit-roots first, Sylvester + Lyapunov column-by-column
//            Jacobian build, O(n^2) one-time cost).
//
// Calling convention (matches what the R wrapper passes):
//   kf_adjoint_diffuse_cpp(Y, TT, RR, ZZ, DD, d_obs, Sigma_e,
//                          me_variance, ur_tol)
//   -> List(loglik, G_TT, G_RR, G_ZZ, G_DD, g_d, G_Sig,
//           stage1_ok, stage2_ok, nunit, d_diffuse,
//           bar_P_inf_0, bar_P_star_0, ok)
//
// Constraints:
//   - No NAs during diffuse phase (documented requirement).
//   - Return ok=false / null-list on non-finite solves.
//
// Reference: R/gradient-adjoint-diffuse.R (601 lines, extensively commented).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// ---- anonymous helpers --------------------------------------------------------
namespace {

inline arma::mat sym(const arma::mat& X) {
  return 0.5 * (X + X.t());
}

// Lyapunov solve: find P such that P - A P A' = Q
// via (I - kron(A,A)) vec(P) = vec(Q).
// Returns false (and leaves P unmodified) if the system is singular or
// non-finite.
static bool solve_lyapunov(const arma::mat& A, const arma::mat& Q,
                            arma::mat& P_out) {
  const arma::uword n = A.n_rows;
  if (n == 1) {
    double denom = 1.0 - A(0, 0) * A(0, 0);
    if (std::abs(denom) < std::numeric_limits<double>::epsilon() * 100)
      return false;
    double val = Q(0, 0) / denom;
    if (!std::isfinite(val)) return false;
    P_out = arma::mat(1, 1, arma::fill::value(val));
    return true;
  }
  const arma::uword n2 = n * n;
  arma::mat M = arma::eye(n2, n2) - arma::kron(A, A);
  double rc = arma::rcond(M);
  if (!(rc >= std::numeric_limits<double>::epsilon())) return false;
  arma::vec p_vec;
  bool ok = arma::solve(p_vec, M, arma::vectorise(Q),
                        arma::solve_opts::no_approx);
  if (!ok) return false;
  P_out = arma::reshape(p_vec, n, n);
  if (!P_out.is_finite()) return false;
  return true;
}

// Sylvester solve: find X such that A X - X B = C
// via (I_nb ⊗ A - B^T ⊗ I_na) vec(X) = vec(C).
// ns = nrows(A)=nrows(X)=nrows(C), nu = ncols(B)=ncols(X)=ncols(C).
static bool solve_sylvester(const arma::mat& A, const arma::mat& B,
                             const arma::mat& C, arma::mat& X_out) {
  const arma::uword ns = A.n_rows;
  const arma::uword nu = B.n_cols;
  arma::mat M = arma::kron(arma::eye(nu, nu), A) -
                arma::kron(B.t(), arma::eye(ns, ns));
  double rc = arma::rcond(M);
  if (!(rc >= std::numeric_limits<double>::epsilon())) return false;
  arma::vec x_vec;
  bool ok = arma::solve(x_vec, M, arma::vectorise(C),
                        arma::solve_opts::no_approx);
  if (!ok) return false;
  X_out = arma::reshape(x_vec, ns, nu);
  if (!X_out.is_finite()) return false;
  return true;
}

// 2×2 Schur swap at position k (0-indexed): swap diagonal blocks k and k+1.
// Mirrors R .swap_schur_11(T, Q, k) where k is 1-indexed; here k is 0-indexed.
// Input:  Ts (n×n upper quasi-triangular), U (n×n orthogonal Schur vectors)
// Output: Ts and U are modified in-place (same as R's list(T,Q) return).
static void swap_schur_11(arma::mat& Ts, arma::mat& U, arma::uword k) {
  // k and k+1 are 0-indexed positions (R k and k+1 are 1-indexed, so R k => C++ k-1)
  double a   = Ts(k, k);
  double b   = Ts(k + 1, k + 1);
  double tol = 1e-15;
  if (std::abs(a - b) < tol) {
    // eigenvalues equal: just swap the diagonal labels (no rotation needed)
    Ts(k, k)         = b;
    Ts(k + 1, k + 1) = a;
    return;
  }
  double t12 = Ts(k, k + 1);
  double r   = std::sqrt(t12 * t12 + (a - b) * (a - b));
  double cc  = t12 / r;
  double ss  = (a - b) / r;

  // Build the full n×n Givens rotation G
  const arma::uword n = Ts.n_rows;
  arma::mat G = arma::eye(n, n);
  G(k,     k)     =  cc;
  G(k,     k + 1) =  ss;
  G(k + 1, k)     = -ss;
  G(k + 1, k + 1) =  cc;

  // Ts <- G' Ts G
  Ts = G.t() * Ts * G;
  // U <- U G
  U = U * G;
  // Force exact zero in sub-diagonal (numerical cleanup)
  Ts(k + 1, k) = 0.0;
}

// Case-tag for the per-step record
enum StepCase { CASE_A, CASE_B, CASE_STAT };

// Per-step store (union of Case A / B / stat fields)
struct StepRecord {
  StepCase  cas;
  arma::vec s_prev;
  arma::mat P_star_prev;
  arma::mat P_inf_prev;   // only for A and B

  // Common innovation + gain (Case A and stat)
  arma::vec v;
  arma::mat Fi;       // F^{-1}
  arma::mat K;
  arma::mat A;        // TT - K ZZ
  arma::mat B_mat;    // RR - K DD
  arma::vec Fiv;

  // Case B specific
  arma::mat F_inf;
  arma::mat F_inf_inv;
  arma::mat F_star;
  arma::mat M_inf;
  arma::mat M_star;
  arma::mat K0;
};

}  // namespace

// [[Rcpp::export]]
List kf_adjoint_diffuse_cpp(const arma::mat& Y,
                             const arma::mat& TT,
                             const arma::mat& RR,
                             const arma::mat& ZZ,
                             const arma::mat& DD,
                             const arma::vec& d_obs,
                             const arma::mat& Sigma_e,
                             double me_variance,
                             double ur_tol) {

  const arma::uword n_state = TT.n_rows;
  const arma::uword n_obs   = ZZ.n_rows;
  const arma::uword n_exo   = RR.n_cols;
  const arma::uword n_T     = Y.n_cols;

  // Fail list (returned on any non-recoverable error).
  auto fail = [&](bool stage1_ok = false, bool stage2_ok = false) -> List {
    return List::create(
      _["loglik"]       = R_NegInf,
      _["G_TT"]         = R_NilValue,
      _["G_RR"]         = R_NilValue,
      _["G_ZZ"]         = R_NilValue,
      _["G_DD"]         = R_NilValue,
      _["g_d"]          = R_NilValue,
      _["G_Sig"]        = R_NilValue,
      _["stage1_ok"]    = stage1_ok,
      _["stage2_ok"]    = stage2_ok,
      _["nunit"]        = NA_INTEGER,
      _["d_diffuse"]    = NA_INTEGER,
      _["bar_P_inf_0"]  = R_NilValue,
      _["bar_P_star_0"] = R_NilValue,
      _["ok"]           = false
    );
  };

  // ---- Derived constants ----------------------------------------------------
  const arma::mat QQ     = sym(RR * Sigma_e * RR.t());
  const arma::mat HH     = sym(DD * Sigma_e * DD.t());
  const arma::mat SS     = RR * Sigma_e * DD.t();
  const arma::mat tZZ    = ZZ.t();
  const arma::mat me_diag = me_variance * arma::eye(n_obs, n_obs);
  const arma::mat HHme   = HH + me_diag;
  const double ll_const  = -0.5 * static_cast<double>(n_obs) * std::log(2.0 * M_PI);
  const double ll_min    = -1e300;

  // ---- Diffuse initialization via Schur decomposition ----------------------
  // We replicate .kf_diffuse_P0(TT, QQ, ur_tol) entirely in C++ using
  // arma::schur + our inline swap_schur_11.
  arma::mat Ts, U_mat;
  arma::schur(U_mat, Ts, TT);   // TT = U_mat * Ts * U_mat'

  // Push unit-root blocks to the leading block (same swap logic as R).
  arma::uword ur_target = 0;
  for (arma::uword i = 0; i < n_state; ++i) {
    if (std::abs(std::abs(Ts(i, i)) - 1.0) < ur_tol) {
      arma::uword j = i;
      while (j > ur_target) {
        // C++ k is 0-indexed; swap positions j-1 and j
        swap_schur_11(Ts, U_mat, j - 1);
        j--;
      }
      ur_target++;
    }
  }
  const arma::uword nunit = ur_target;
  // U_mat is now the ordered Schur basis.

  // P_inf_0 = U_u U_u'
  arma::mat P_inf_0 = arma::zeros<arma::mat>(n_state, n_state);
  if (nunit > 0) {
    arma::mat U_u = U_mat.cols(0, nunit - 1);
    P_inf_0 = sym(U_u * U_u.t());
  }

  // P_star_0 = U_s Pa_ss U_s'  where Pa_ss solves Lyapunov(T_ss, QQ_ss)
  arma::mat P_star_0 = arma::zeros<arma::mat>(n_state, n_state);
  bool has_stable = (nunit < n_state);
  arma::mat Pa_ss;    // for Stage 2
  bool pa_ss_ok = false;
  if (has_stable) {
    arma::mat U_s   = U_mat.cols(nunit, n_state - 1);
    arma::mat T_ss  = Ts.submat(nunit, nunit, n_state - 1, n_state - 1);
    arma::mat QQ_ss = U_s.t() * QQ * U_s;
    arma::mat Pa_tmp;
    if (solve_lyapunov(T_ss, QQ_ss, Pa_tmp)) {
      Pa_ss = sym(Pa_tmp);
      pa_ss_ok = true;
      P_star_0 = sym(U_s * Pa_ss * U_s.t());
    }
  }

  // ---- Forward pass ---------------------------------------------------------
  const double diffuse_tol = 1e-10;
  const double conv_tol    = 1e-8;
  const arma::uword cap    = std::min(n_T, (arma::uword)100);

  std::vector<StepRecord> step_store(n_T);
  std::vector<bool>       step_valid(n_T, false);

  arma::vec s      = arma::zeros<arma::vec>(n_state);
  arma::mat P_inf  = P_inf_0;
  arma::mat P_star = P_star_0;
  double loglik    = 0.0;
  int d_diffuse    = NA_INTEGER;
  bool in_diffuse  = (nunit > 0);

  for (arma::uword t = 0; t < n_T; ++t) {
    arma::vec y_t = Y.col(t);
    // No NAs during diffuse phase (documented requirement).
    if (!y_t.is_finite()) return fail();

    arma::vec v = y_t - d_obs - ZZ * s;

    StepRecord rec;

    if (in_diffuse) {
      arma::mat F_inf  = sym(ZZ * P_inf  * tZZ);
      arma::mat F_star = sym(ZZ * P_star * tZZ + HHme);
      double scale_star = std::max(1.0, arma::abs(F_star).max());

      double max_Finf = arma::abs(F_inf).max();

      if (max_Finf < diffuse_tol * scale_star) {
        // ---- Case A: standard KF step on (s, P_star); propagate P_inf -------
        arma::mat PZ = P_star * tZZ;
        arma::mat Ft = sym(ZZ * PZ + HHme);
        arma::mat Fc;
        if (!arma::chol(Fc, Ft)) return fail();
        arma::mat Fi   = arma::inv_sympd(Ft);
        double ldf     = 2.0 * arma::accu(arma::log(Fc.diag()));
        arma::vec Fiv  = Fi * v;
        double ll_t    = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
        if (!std::isfinite(ll_t) || ll_t < ll_min) return fail();
        loglik += ll_t;

        arma::mat K = (TT * PZ + SS) * Fi;
        arma::mat A = TT - K * ZZ;
        arma::mat B = RR - K * DD;

        arma::mat P_inf_new  = sym(TT * P_inf * TT.t());
        arma::mat P_star_new = sym(A * P_star * A.t() + B * Sigma_e * B.t());

        rec.cas         = CASE_A;
        rec.s_prev      = s;
        rec.P_star_prev = P_star;
        rec.P_inf_prev  = P_inf;
        rec.v           = v;
        rec.Fi          = Fi;
        rec.K           = K;
        rec.A           = A;
        rec.B_mat       = B;
        rec.Fiv         = Fiv;
        step_valid[t]   = true;

        s      = TT * s + K * v;
        P_star = P_star_new;
        P_inf  = P_inf_new;

      } else {
        // ---- Case B: full diffuse step ---------------------------------------
        // Check conditioning of F_inf
        double rc = arma::rcond(F_inf);
        arma::mat Fc_inf;
        if (!arma::chol(Fc_inf, F_inf) || !std::isfinite(rc) ||
            rc <= diffuse_tol) return fail();

        arma::mat F_inf_inv = arma::inv_sympd(F_inf);
        double ll_t = -0.5 * 2.0 * arma::accu(arma::log(Fc_inf.diag()));
        if (!std::isfinite(ll_t)) return fail();
        loglik += ll_t;

        arma::mat M_inf  = TT * P_inf  * tZZ;
        arma::mat M_star = TT * P_star * tZZ + SS;
        arma::mat K0     = M_inf * F_inf_inv;

        arma::mat P_inf_new  = sym(TT * P_inf  * TT.t() - K0 * F_inf * K0.t());
        arma::mat P_star_new = sym(TT * P_star * TT.t() + QQ -
                                   K0 * M_star.t() - M_star * K0.t() +
                                   K0 * F_star * K0.t());

        rec.cas         = CASE_B;
        rec.s_prev      = s;
        rec.P_star_prev = P_star;
        rec.P_inf_prev  = P_inf;
        rec.v           = v;
        rec.F_inf       = F_inf;
        rec.F_inf_inv   = F_inf_inv;
        rec.F_star      = F_star;
        rec.M_inf       = M_inf;
        rec.M_star      = M_star;
        rec.K0          = K0;
        step_valid[t]   = true;

        s      = TT * s + K0 * v;
        P_star = P_star_new;
        P_inf  = P_inf_new;
      }

      // Check for diffuse phase convergence
      double max_Pinf  = arma::abs(P_inf).max();
      double max_Pstar = std::max(1.0, arma::abs(P_star).max());
      if (max_Pinf < conv_tol * max_Pstar) {
        d_diffuse  = static_cast<int>(t + 1);   // 1-indexed like R
        in_diffuse = false;
      } else if (t + 1 >= cap) {
        // Did not converge
        return fail();
      }

    } else {
      // ---- Stationary tail --------------------------------------------------
      arma::mat PZ = P_star * tZZ;
      arma::mat Ft = sym(ZZ * PZ + HHme);
      arma::mat Fc;
      if (!arma::chol(Fc, Ft)) return fail();
      arma::mat Fi   = arma::inv_sympd(Ft);
      double ldf     = 2.0 * arma::accu(arma::log(Fc.diag()));
      arma::vec Fiv  = Fi * v;
      double ll_t    = ll_const - 0.5 * (ldf + arma::dot(v, Fiv));
      if (!std::isfinite(ll_t) || ll_t < ll_min) return fail();
      loglik += ll_t;

      arma::mat K = (TT * PZ + SS) * Fi;
      arma::mat A = TT - K * ZZ;
      arma::mat B = RR - K * DD;

      rec.cas         = CASE_STAT;
      rec.s_prev      = s;
      rec.P_star_prev = P_star;
      rec.v           = v;
      rec.Fi          = Fi;
      rec.K           = K;
      rec.A           = A;
      rec.B_mat       = B;
      rec.Fiv         = Fiv;
      step_valid[t]   = true;

      s      = TT * s + K * v;
      P_star = sym(A * P_star * A.t() + B * Sigma_e * B.t());
    }

    step_store[t] = rec;
  }

  // ---- Backward sweep -------------------------------------------------------
  arma::vec bar_s      = arma::zeros<arma::vec>(n_state);
  arma::mat bar_P_star = arma::zeros<arma::mat>(n_state, n_state);
  arma::mat bar_P_inf  = arma::zeros<arma::mat>(n_state, n_state);

  arma::mat G_TT  = arma::zeros<arma::mat>(n_state, n_state);
  arma::mat G_RR  = arma::zeros<arma::mat>(n_state, n_exo);
  arma::mat G_ZZ  = arma::zeros<arma::mat>(n_obs,   n_state);
  arma::mat G_DD  = arma::zeros<arma::mat>(n_obs,   n_exo);
  arma::vec g_d   = arma::zeros<arma::vec>(n_obs);
  arma::mat G_Sig = arma::zeros<arma::mat>(n_exo,   n_exo);

  for (arma::sword t_signed = static_cast<arma::sword>(n_T) - 1;
       t_signed >= 0; --t_signed) {
    const arma::uword t = static_cast<arma::uword>(t_signed);
    if (!step_valid[t]) continue;

    const StepRecord& st = step_store[t];

    if (st.cas == CASE_STAT || st.cas == CASE_A) {
      // ---- Standard KF step adjoint (operates on bar_P_star) ---------------
      // Mirrors R lines 263-339.
      const arma::vec& s_prev  = st.s_prev;
      const arma::mat& P_prev  = st.P_star_prev;
      const arma::vec& v       = st.v;
      const arma::mat& Fi      = st.Fi;
      const arma::mat& K       = st.K;
      const arma::mat& A       = st.A;
      const arma::mat& B       = st.B_mat;
      const arma::vec& Fiv     = st.Fiv;

      // Adjoint of P_t = sym(A P_prev A' + B Sig B')
      bar_P_star = sym(bar_P_star);
      arma::mat bar_A           = 2.0 * bar_P_star * A * P_prev;
      arma::mat bar_B           = 2.0 * bar_P_star * B * Sigma_e;
      arma::mat bar_Pprev_AP    = A.t() * bar_P_star * A;
      G_Sig += B.t() * bar_P_star * B;

      // Adjoint of s_t = TT s_prev + K v
      G_TT      += bar_s * s_prev.t();
      arma::mat bar_K     = bar_s * v.t();
      arma::vec bar_v     = K.t() * bar_s;
      arma::vec bar_s_prev = TT.t() * bar_s;

      // Adjoint of ll_t = -0.5*(log|F| + v'Fi v)
      arma::mat bar_F = -0.5 * (Fi - Fiv * Fiv.t());
      bar_v -= Fiv;

      // Adjoint of A = TT - K ZZ
      G_TT  += bar_A;
      G_ZZ  -= K.t() * bar_A;
      bar_K -= bar_A * tZZ;        // tZZ = ZZ.t(), so bar_A * tZZ = bar_A * ZZ'

      // Adjoint of B = RR - K DD
      G_RR  += bar_B;
      G_DD  -= K.t() * bar_B;
      bar_K -= bar_B * DD.t();

      // Adjoint of K = (TT P_prev ZZ' + SS) Fi
      bar_F -= K.t() * bar_K * Fi;
      bar_F  = sym(bar_F);
      arma::mat bar_Mnum = bar_K * Fi;

      G_TT         += bar_Mnum * ZZ * P_prev;
      G_ZZ         += bar_Mnum.t() * TT * P_prev;
      arma::mat bar_Pprev_Mnum = TT.t() * bar_Mnum * ZZ;
      G_RR         += bar_Mnum * DD * Sigma_e;
      G_DD         += bar_Mnum.t() * RR * Sigma_e;
      G_Sig        += RR.t() * bar_Mnum * DD;

      // Adjoint of F = sym(ZZ P_prev ZZ' + HH + me_diag)
      arma::mat bar_Pprev_F = tZZ * bar_F * ZZ;
      G_ZZ  += 2.0 * bar_F * ZZ * P_prev;
      G_DD  += 2.0 * bar_F * DD * Sigma_e;
      G_Sig += DD.t() * bar_F * DD;

      // Adjoint of v = y - d - ZZ s_prev
      g_d        -= bar_v;
      G_ZZ       -= bar_v * s_prev.t();
      bar_s_prev -= tZZ * bar_v;

      arma::mat bar_P_star_new = sym(bar_Pprev_AP + bar_Pprev_Mnum + bar_Pprev_F);

      // In Case A: also propagate bar_P_inf through P_inf <- TT P_inf TT'
      if (st.cas == CASE_A) {
        const arma::mat& P_inf_prev = st.P_inf_prev;
        arma::mat bar_P_inf_sym = sym(bar_P_inf);
        // G_TT from TT P_inf_prev TT': 2 * bar_P_inf * TT * P_inf_prev
        G_TT      += 2.0 * bar_P_inf_sym * TT * P_inf_prev;
        bar_P_inf  = TT.t() * bar_P_inf_sym * TT;
      }

      bar_s      = bar_s_prev;
      bar_P_star = bar_P_star_new;

    } else {
      // ---- Case B adjoint ---------------------------------------------------
      // Mirrors R lines 341-466.
      const arma::vec& s_prev      = st.s_prev;
      const arma::mat& P_star_prev = st.P_star_prev;
      const arma::mat& P_inf_prev  = st.P_inf_prev;
      const arma::vec& v           = st.v;
      const arma::mat& F_inf       = st.F_inf;
      const arma::mat& Fi          = st.F_inf_inv;   // F_inf^{-1}
      const arma::mat& F_star      = st.F_star;
      const arma::mat& M_star      = st.M_star;
      const arma::mat& K0          = st.K0;
      arma::mat tK0 = K0.t();

      arma::mat bar_Pstar_t = sym(bar_P_star);
      arma::mat bar_Pinf_t  = sym(bar_P_inf);

      // ----- 1. Adjoint of P_star' = sym(TT P_star TT' + QQ
      //                               - K0 M_star' - M_star K0' + K0 F_star K0')
      G_TT += 2.0 * bar_Pstar_t * TT * P_star_prev;
      arma::mat bar_Pstar_from_TT = TT.t() * bar_Pstar_t * TT;

      // QQ = RR Sig RR'
      G_Sig += RR.t() * bar_Pstar_t * RR;
      G_RR  += 2.0 * bar_Pstar_t * RR * Sigma_e;

      // -sym(K0 M_star')
      arma::mat bar_K0    = -2.0 * bar_Pstar_t * M_star;
      arma::mat bar_Mstar = -2.0 * bar_Pstar_t * K0;

      // +sym(K0 F_star K0')
      bar_K0             += 2.0 * bar_Pstar_t * K0 * F_star;
      arma::mat bar_Fstar = tK0 * bar_Pstar_t * K0;

      // ----- 2. Adjoint of P_inf' = sym(TT P_inf TT' - K0 F_inf K0')
      G_TT += 2.0 * bar_Pinf_t * TT * P_inf_prev;
      arma::mat bar_Pinf_from_TT = TT.t() * bar_Pinf_t * TT;

      // -sym(K0 F_inf K0')
      bar_K0             -= 2.0 * bar_Pinf_t * K0 * F_inf;
      arma::mat bar_F_inf = -tK0 * bar_Pinf_t * K0;

      // ----- 3. Adjoint of ll_t = -0.5 log|F_inf|
      bar_F_inf -= 0.5 * Fi;

      // ----- 4. Adjoint of s' = TT s_prev + K0 v
      G_TT              += bar_s * s_prev.t();
      bar_K0            += bar_s * v.t();
      arma::vec bar_v    = tK0 * bar_s;
      arma::vec bar_s_prev = TT.t() * bar_s;

      // ----- 5. Adjoint of v = y - d - ZZ s_prev
      g_d        -= bar_v;
      G_ZZ       -= bar_v * s_prev.t();
      bar_s_prev -= tZZ * bar_v;

      // ----- 6. Adjoint of K0 = M_inf F_inf^{-1}
      arma::mat bar_Minf = bar_K0 * Fi;
      bar_F_inf         -= tK0 * bar_K0 * Fi;
      bar_F_inf          = sym(bar_F_inf);

      // ----- 7. Adjoint of M_star = TT P_star ZZ' + SS
      G_TT  += bar_Mstar * ZZ * P_star_prev;
      G_ZZ  += bar_Mstar.t() * TT * P_star_prev;
      arma::mat bar_Pstar_from_Mstar = TT.t() * bar_Mstar * ZZ;
      // From SS = RR Sig DD'
      G_RR  += bar_Mstar * DD * Sigma_e;
      G_DD  += bar_Mstar.t() * RR * Sigma_e;
      G_Sig += RR.t() * bar_Mstar * DD;

      // ----- 8. Adjoint of M_inf = TT P_inf ZZ'
      G_TT  += bar_Minf * ZZ * P_inf_prev;
      G_ZZ  += bar_Minf.t() * TT * P_inf_prev;
      arma::mat bar_Pinf_from_Minf = TT.t() * bar_Minf * ZZ;

      // ----- 9. Adjoint of F_star = sym(ZZ P_star ZZ' + HH + me_diag)
      bar_Fstar = sym(bar_Fstar);
      G_ZZ  += 2.0 * bar_Fstar * ZZ * P_star_prev;
      arma::mat bar_Pstar_from_Fstar = tZZ * bar_Fstar * ZZ;
      G_DD  += 2.0 * bar_Fstar * DD * Sigma_e;
      G_Sig += DD.t() * bar_Fstar * DD;

      // ----- 10. Adjoint of F_inf = sym(ZZ P_inf ZZ')
      // bar_F_inf already symmetrized above.
      G_ZZ  += 2.0 * bar_F_inf * ZZ * P_inf_prev;
      arma::mat bar_Pinf_from_Finf = tZZ * bar_F_inf * ZZ;

      // ----- Collect bar_P_{t-1}
      bar_P_inf  = sym(bar_Pinf_from_TT + bar_Pinf_from_Minf + bar_Pinf_from_Finf);
      bar_P_star = sym(bar_Pstar_from_TT + bar_Pstar_from_Mstar + bar_Pstar_from_Fstar);
      bar_s      = bar_s_prev;
    }
  }  // end backward loop

  // At t=0: bar_P_star = bar_P_star_0, bar_P_inf = bar_P_inf_0
  arma::mat bar_P_inf_0  = bar_P_inf;
  arma::mat bar_P_star_0 = bar_P_star;

  // ============================================================
  // Stage 2: Adjoint through P_star_0 (Lyapunov in Schur basis)
  //          and P_inf_0 (spectral projector)
  // ============================================================
  // We use the SAME Ts / U_mat already computed in diffuse init above.
  // (They were already sorted with unit roots in the leading block.)

  bool stage2_ok = true;

  if (nunit > 0) {
    // Unit and stable index partitions (0-indexed)
    arma::uvec idx_u = arma::regspace<arma::uvec>(0, nunit - 1);
    arma::mat  U_u   = U_mat.cols(idx_u);
    arma::mat  T_uu  = Ts.submat(0, 0, nunit - 1, nunit - 1);

    arma::mat U_s, T_ss, QQ_ss_cur;
    arma::uvec idx_s;
    if (has_stable) {
      idx_s   = arma::regspace<arma::uvec>(nunit, n_state - 1);
      U_s     = U_mat.cols(idx_s);
      T_ss    = Ts.submat(nunit, nunit, n_state - 1, n_state - 1);
      QQ_ss_cur = U_s.t() * QQ * U_s;
      if (!pa_ss_ok) stage2_ok = false;
    }

    if (stage2_ok) {
      // bar_TT_init: for each basis direction e_{ab}, compute the forward
      // directional derivative of (P_inf_0, P_star_0) and contract with the
      // backward adjoint bars.

      arma::mat bar_TT_init = arma::zeros<arma::mat>(n_state, n_state);
      arma::mat bar_QQ_init = arma::zeros<arma::mat>(n_state, n_state);

      // Lambda for computing dP_inf, dP_star given (dTT_full, dQQ_full)
      // Mirrors R's fwd() function in the Stage 2 section.
      auto fwd = [&](const arma::mat& dTT,
                     const arma::mat& dQQ,
                     arma::mat& dP_inf_out,
                     arma::mat& dP_star_out,
                     bool use_dTT, bool use_dQQ) -> bool {

        arma::mat Om_su;
        if (has_stable && use_dTT) {
          // G_su = (U' dTT U)[idx_s, idx_u]
          arma::mat tmp_UdTTU = U_mat.t() * dTT * U_mat;
          arma::mat Gsu = tmp_UdTTU.submat(nunit, 0, n_state - 1, nunit - 1);
          // Sylvester: T_ss Om - Om T_uu = -G_su
          if (!solve_sylvester(T_ss, T_uu, -Gsu, Om_su)) return false;
        } else {
          Om_su = arma::zeros<arma::mat>(n_state - nunit, nunit);
        }

        arma::mat dU_u;
        if (has_stable) {
          dU_u = U_s * Om_su;
        } else {
          dU_u = arma::zeros<arma::mat>(n_state, nunit);
        }
        dP_inf_out = sym(dU_u * U_u.t() + U_u * dU_u.t());

        dP_star_out = arma::zeros<arma::mat>(n_state, n_state);
        if (has_stable) {
          arma::mat dU_s = -U_u * Om_su.t();
          arma::mat dTm  = use_dTT ? dTT : arma::zeros<arma::mat>(n_state, n_state);
          arma::mat dQm  = use_dQQ ? dQQ : arma::zeros<arma::mat>(n_state, n_state);
          arma::mat dT_ss = dU_s.t() * TT * U_s +
                            U_s.t() * dTm * U_s +
                            U_s.t() * TT * dU_s;
          arma::mat dQ_ss = dU_s.t() * QQ * U_s +
                            U_s.t() * dQm * U_s +
                            U_s.t() * QQ * dU_s;
          arma::mat rhs = dT_ss * Pa_ss * T_ss.t() +
                          T_ss  * Pa_ss * dT_ss.t() + dQ_ss;
          arma::mat dPa;
          if (!solve_lyapunov(T_ss, sym(rhs), dPa)) return false;
          dPa = sym(dPa);
          dP_star_out = sym(dU_s * Pa_ss * U_s.t() +
                            U_s  * dPa   * U_s.t() +
                            U_s  * Pa_ss * dU_s.t());
        }
        return true;
      };

      // --- bar_TT_init: contract with per-column basis directions of TT ------
      arma::mat E = arma::zeros<arma::mat>(n_state, n_state);
      for (arma::uword a = 0; a < n_state && stage2_ok; ++a) {
        for (arma::uword b = 0; b < n_state && stage2_ok; ++b) {
          E(a, b) = 1.0;
          arma::mat dPinf, dPstar;
          if (!fwd(E, E, dPinf, dPstar, true, false)) {
            stage2_ok = false;
            break;
          }
          bar_TT_init(a, b) = arma::accu(bar_P_inf_0 % dPinf) +
                               arma::accu(bar_P_star_0 % dPstar);
          E(a, b) = 0.0;
        }
      }
      if (stage2_ok) G_TT += bar_TT_init;

      // --- bar_QQ_init: only P_star_0 depends on QQ --------------------------
      if (stage2_ok && has_stable) {
        E.zeros();
        for (arma::uword a = 0; a < n_state && stage2_ok; ++a) {
          for (arma::uword b = 0; b < n_state && stage2_ok; ++b) {
            E(a, b) = 1.0;
            arma::mat dPinf, dPstar;
            // dTT = 0, dQQ = E
            if (!fwd(E, E, dPinf, dPstar, false, true)) {
              stage2_ok = false;
              break;
            }
            bar_QQ_init(a, b) = arma::accu(bar_P_star_0 % dPstar);
            E(a, b) = 0.0;
          }
        }
        if (stage2_ok) {
          bar_QQ_init = sym(bar_QQ_init);
          // QQ = RR Sigma_e RR' -> G_Sig, G_RR
          G_Sig += RR.t() * bar_QQ_init * RR;
          G_RR  += 2.0 * bar_QQ_init * RR * Sigma_e;
        }
      }
    }
  }

  return List::create(
    _["loglik"]       = loglik,
    _["G_TT"]         = G_TT,
    _["G_RR"]         = G_RR,
    _["G_ZZ"]         = G_ZZ,
    _["G_DD"]         = G_DD,
    _["g_d"]          = g_d,
    _["G_Sig"]        = G_Sig,
    _["stage1_ok"]    = true,
    _["stage2_ok"]    = stage2_ok,
    _["nunit"]        = (int)nunit,
    _["d_diffuse"]    = d_diffuse,
    _["bar_P_inf_0"]  = bar_P_inf_0,
    _["bar_P_star_0"] = bar_P_star_0,
    _["ok"]           = true
  );
}
