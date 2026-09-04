// src/tpf_propagate.cpp
// -------------------------------------------------------------------------
// Rcpp/Armadillo kernels for the Tempered Particle Filter (TPF).
//
// Exports two functions consumed by tpf_run_period() via .HAS_RCPP_TPF()
// dispatch in R/tpf-likelihood.R (mirrors .HAS_RCPP_KALMAN_UNI() pattern).
//
// Kronecker convention — CRITICAL (Landmine 4):
//   hxx columns: (state FAST x state SLOW) -> matching vector is (x1 %x% x1)
//   hxu columns: (state FAST x exo SLOW)   -> matching vector is (e  %x% x1)
//   huu columns: (exo  FAST x exo  SLOW)   -> matching vector is (e  %x% e)
// This exactly replicates simulate_model_order2() lines 915-923 in
// R/solve-perturbation-order2.R. Any reversal silently produces wrong second-
// order terms that the linear-parity test (ghxx=0) will NOT catch.
//
// ghss sign/factor (Landmine 5):
//   x2 update:  + 0.5 * hss    (hss = ghss[state_idx], n_s-vector)
//   obs mean:   + 0.5 * ghss[obs_idx]    (applied in R caller)
//
// References:
//   Herbst & Schorfheide (2019) J. Econometrics 210(1):26-44.
//   Andreasen, Fernandez-Villaverde & Rubio-Ramirez (2018) RES 85:1-49.
//   src/kalman_univariate.cpp for the include/depend style.
// -------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using Rcpp::List;
using Rcpp::_;

// [[Rcpp::export]]
arma::mat tpf_propagate_particles(
    const arma::mat& particles,   // (2*n_s) x N
    const arma::mat& shocks,      // n_e x N
    const arma::mat& hx,          // n_s x n_s
    const arma::mat& hu,          // n_s x n_e
    const arma::mat& hxx,         // n_s x (n_s*n_s)
    const arma::mat& hxu,         // n_s x (n_s*n_e)  note: (state FAST, exo SLOW)
    const arma::mat& huu,         // n_s x (n_e*n_e)
    const arma::vec& hss          // n_s
) {
  const arma::uword n_s = hx.n_rows;
  const arma::uword n_e = hu.n_cols;
  const arma::uword N   = particles.n_cols;

  arma::mat x1 = particles.rows(0, n_s - 1);               // n_s x N
  arma::mat x2 = particles.rows(n_s, 2 * n_s - 1);         // n_s x N

  // First-order state update: x1_new = hx x1 + hu e
  arma::mat x1_new = hx * x1 + hu * shocks;

  // Second-order state update (pruned):
  //   x2_new = hx x2_prev
  //           + 0.5 hxx (x1_prev %x% x1_prev)
  //           + hxu     (e %x% x1_prev)
  //           + 0.5 huu (e %x% e)
  //           + 0.5 hss
  arma::mat x2_new = hx * x2;

  // Check whether each second-order matrix is numerically zero to skip
  // expensive Kronecker computations when running on a linear model.
  const bool hxx_nonzero = arma::any(arma::vectorise(arma::abs(hxx)) > 0.0);
  const bool hxu_nonzero = arma::any(arma::vectorise(arma::abs(hxu)) > 0.0);
  const bool huu_nonzero = arma::any(arma::vectorise(arma::abs(huu)) > 0.0);
  const bool hss_nonzero = arma::any(arma::abs(hss) > 0.0);

  if (hxx_nonzero) {
    // kron_xx: (n_s^2) x N where column i = x1[:,i] kron x1[:,i]
    // fast-first: outer loop = x1[j], inner loop = x1[k]
    // => element (j*n_s + k, i) = x1[j,i] * x1[k,i]
    // This is the same ordering as R's kronecker(v, v) which is v %o% v
    // vectorised column-major: element [j + k*n_s] = v[j] * v[k]
    // R %x% is Kronecker with LEFT argument varying FASTEST:
    //   kronecker(a, b)[i] = a[i %% length(a)] * b[i %/% length(a)]
    // So kronecker(x1, x1)[j + k*n_s] = x1[j] * x1[k]  (j fast, k slow)
    arma::mat kron_xx(n_s * n_s, N);
    for (arma::uword i = 0; i < N; ++i) {
      const arma::vec& v = x1.col(i);
      for (arma::uword k = 0; k < n_s; ++k) {
        for (arma::uword j = 0; j < n_s; ++j) {
          kron_xx(j + k * n_s, i) = v(j) * v(k);   // (j fast, k slow)
        }
      }
    }
    x2_new += 0.5 * hxx * kron_xx;
  }

  if (hxu_nonzero) {
    // kron_ex: (n_e * n_s) x N where column i = e[:,i] kron x1[:,i]
    // hxu cols: (state FAST, exo SLOW) -> matching Kronecker is e %x% x1.
    //
    // 2026-08-04 obs-tensor-fix session fix (PRE-EXISTING bug, found while
    // adding the R-vs-C++ obs-tensor parity test on rbc2shock, a 2-shock
    // model): R's kronecker(a,b) has `a` SLOW (outer) and `b` FAST (inner)
    // -- verified empirically, e.g. kronecker(c(10,20), c(1,2,3)) ==
    // c(10,20,30, 20,40,60) == a[1]*b, a[2]*b concatenated. So
    // kronecker(e, x1)[j + k*n_s] = e[k] * x1[j]  (j FAST over x1, k SLOW
    // over e). The code below previously had e FAST / x1 SLOW with an
    // n_e-strided index -- the WRONG layout whenever n_e != n_s (invisible
    // on the single-shock fixtures that exercised this kernel before).
    // R reference: .tpf_propagate_R's `kron_ex[, i] <- kronecker(shocks[, i], x1[, i])`.
    arma::mat kron_ex(n_e * n_s, N);
    for (arma::uword i = 0; i < N; ++i) {
      const arma::vec& ei  = shocks.col(i);
      const arma::vec& xi  = x1.col(i);
      for (arma::uword k = 0; k < n_e; ++k) {
        for (arma::uword j = 0; j < n_s; ++j) {
          kron_ex(j + k * n_s, i) = ei(k) * xi(j);   // x1 FAST, e SLOW
        }
      }
    }
    x2_new += hxu * kron_ex;
  }

  if (huu_nonzero) {
    // kron_ee: (n_e^2) x N where column i = e[:,i] kron e[:,i]
    arma::mat kron_ee(n_e * n_e, N);
    for (arma::uword i = 0; i < N; ++i) {
      const arma::vec& ei = shocks.col(i);
      for (arma::uword k = 0; k < n_e; ++k) {
        for (arma::uword j = 0; j < n_e; ++j) {
          kron_ee(j + k * n_e, i) = ei(j) * ei(k);
        }
      }
    }
    x2_new += 0.5 * huu * kron_ee;
  }

  if (hss_nonzero) {
    // Broadcast hss (n_s column vector) across all N columns
    // x2_new += 0.5 * hss (repeated N times)
    x2_new.each_col() += 0.5 * hss;
  }

  return arma::join_cols(x1_new, x2_new);
}


// -------------------------------------------------------------------------
// tpf_run_period_cpp: C++ per-period tempering loop
//
// Ports the entire within-period adaptive phi-tempering loop from
// tpf_run_period() (R/tpf-likelihood.R:219-333).
//
// HARD INVARIANTS (from commit 2267e8d):
//   A. Always resample at phi=1 regardless of ESS (final-stage weights).
//   B. RWMH mutation only on equally-weighted clouds (after resample).
//
// RNG: uses Rcpp::RNGScope + Rcpp::rnorm / Rcpp::runif so that the same
// seed gives IDENTICAL log-likelihood as the R reference path.
//
// Returns: List with
//   particles       - updated (2*n_s) x N matrix  (period-t particles)
//   log_lik_contrib - scalar double
//   phi_schedule    - arma::vec of phi values visited (for diagnostics)
// -------------------------------------------------------------------------

// Helper: log-sum-exp of an arma::vec
static inline double lse(const arma::vec& x) {
  double mx = arma::max(x);
  if (!std::isfinite(mx)) return -arma::datum::inf;
  return mx + std::log(arma::sum(arma::exp(x - mx)));
}

// Helper: ESS from log weights
static inline double ess_from_logw(const arma::vec& log_w) {
  arma::vec lw = log_w - arma::max(log_w);
  arma::vec w  = arma::exp(lw);
  w /= arma::sum(w);
  return 1.0 / arma::dot(w, w);
}

// Helper: systematic resample -- uses ONE runif(1) draw from R's RNG stream
static arma::uvec systematic_resample_cpp(const arma::vec& w_norm, arma::uword N) {
  arma::vec cw = arma::cumsum(w_norm);
  cw(N - 1) = 1.0;  // guard against rounding

  // Draw ONE uniform via R's RNG (advances R's stream by 1, same as R path)
  double u0 = Rcpp::as<double>(Rcpp::runif(1, 0.0, 1.0));

  arma::uvec idx(N);
  arma::uword j = 0;
  for (arma::uword i = 0; i < N; ++i) {
    double u = ((double)i + u0) / (double)N;
    while (j < N - 1 && cw(j) < u) ++j;
    idx(i) = j;
  }
  return idx;
}

// Helper: sorted systematic resample (Deligiannidis et al. CPM).
//
// Sorts particles by their first state dimension (x1[0], i.e. particles(0,:))
// before applying systematic resampling. Sorting makes the particle-index map
// a continuous function of the driving uniform u0, which is critical for the
// CPM loglik correlation: when u0 is correlated across theta/theta', the
// resampled particle trajectories are also correlated.
//
// Sort choice: first state dimension (x1[0]).
//   - For the AR(1) SBC fixture (n_state=1) this is the only dimension.
//   - For multi-dimensional state it is an approximation; a Hilbert-curve
//     projection sort would be optimal but is much more complex. First-dim
//     sort is acceptable for the small-n_state TPF models in this package.
//
// u0_in: pre-drawn uniform in (0,1). If NaN, draw from R's RNG stream
//        (fallback to non-CPM, bit-identical to systematic_resample_cpp).
//
// Returns: (sort_order, resample_idx) packed as:
//   output[0..N-1]   = sort_order (arma::uvec of pre-sort particle indices)
//   The actual resampled indices (in original particle order) are computed
//   internally and returned as idx.
static arma::uvec systematic_resample_sorted_cpp(
    const arma::mat& particles,   // (2*n_s) x N
    const arma::vec& w_norm,      // N  (must sum to 1)
    arma::uword N,
    double u0_in                  // NaN -> draw from RNG (fallback)
) {
  // Sort particle indices by first state dimension (particles(0, :))
  arma::uvec sort_order = arma::sort_index(particles.row(0).t());

  // Reorder weights according to sort
  arma::vec w_sorted(N);
  for (arma::uword k = 0; k < N; ++k) w_sorted(k) = w_norm(sort_order(k));
  arma::vec cw = arma::cumsum(w_sorted);
  cw(N - 1) = 1.0;

  // Uniform: use supplied value or draw from R's RNG
  double u0;
  if (std::isnan(u0_in)) {
    u0 = Rcpp::as<double>(Rcpp::runif(1, 0.0, 1.0));
  } else {
    u0 = u0_in;
  }

  // Systematic resampling on sorted weights
  arma::uvec idx_sorted(N);
  arma::uword j = 0;
  for (arma::uword i = 0; i < N; ++i) {
    double u = ((double)i + u0) / (double)N;
    while (j < N - 1 && cw(j) < u) ++j;
    idx_sorted(i) = j;
  }

  // Map back to original particle indices via sort_order
  arma::uvec idx(N);
  for (arma::uword i = 0; i < N; ++i) idx(i) = sort_order(idx_sorted(i));
  return idx;
}

// Helper: adaptive bisection to find next phi
//   Returns phi in (phi_curr, 1] s.t. ESS(delta_phi * log_liks) ~ ess_target*N
static double next_phi_cpp(const arma::vec& log_liks, double phi_curr,
                            double ess_target, arma::uword N) {
  double target_ess = ess_target * N;

  // Try phi=1 first
  {
    arma::vec inc = (1.0 - phi_curr) * log_liks;
    if (ess_from_logw(inc) >= target_ess) return 1.0;
  }

  double lo = phi_curr, hi = 1.0;
  for (int iter = 0; iter < 50; ++iter) {
    double mid = 0.5 * (lo + hi);
    arma::vec inc = (mid - phi_curr) * log_liks;
    double e = ess_from_logw(inc);
    if (e > target_ess) lo = mid; else hi = mid;
    if (hi - lo < 1e-8) break;
  }
  return lo;
}

// Helper: per-column Kronecker product, R convention verified empirically
// (kronecker(a,b) with length(a)=2, length(b)=3 gives a[1]*b, a[2]*b
// concatenated -- i.e. A is SLOW (outer), B is FAST (inner)):
//   out[j + i*nb, col] = A[i,col] * B[j,col]   (j in [0,nb), i in [0,na))
// This is the SAME convention .tpf_propagate_R's kron_ex etc. rely on via
// literal kronecker() calls; it matters (produces a different result than
// the A-fast/B-slow layout) whenever nA != nB, e.g. kron(shocks, x1).
static arma::mat kron_cols_2(const arma::mat& A, const arma::mat& B) {
  const arma::uword na = A.n_rows, nb = B.n_rows, N = A.n_cols;
  arma::mat out(na * nb, N);
  for (arma::uword c = 0; c < N; ++c) {
    for (arma::uword i = 0; i < na; ++i) {
      for (arma::uword j = 0; j < nb; ++j) {
        out(j + i * nb, c) = A(i, c) * B(j, c);
      }
    }
  }
  return out;
}

// Helper: inline log-weights kernel (avoids a separate Rcpp round-trip)
//
// 2026-08-04 obs-tensor fix: the observation mean also carries the
// quadratic terms in (x1_prev, e_t) that simulate_model_order2's observable
// reconstruction keeps -- 0.5*hxx_obs(x1 kron x1) + hxu_obs(e kron x1) +
// 0.5*huu_obs(e kron e) -- mirroring R/tpf-likelihood.R's .tpf_log_weights_R.
// hxx_obs/hxu_obs/huu_obs are all-zero matrices (never NULL) when the model
// has no nonlinear obs tensors, so the added terms are exact zeros and the
// LINEAR-model nesting stays bit-identical.
static arma::vec tpf_lw_inline(
    const arma::mat& particles,
    const arma::vec& y_t,
    const arma::mat& ZZ,
    const arma::mat& DD,
    const arma::mat& shocks,
    const arma::vec& d_obs,
    const arma::vec& ghss_obs,
    double me_variance,
    double phi,
    const arma::mat& hxx_obs,
    const arma::mat& hxu_obs,
    const arma::mat& huu_obs
) {
  const arma::uword n_s  = ZZ.n_cols;
  const arma::uword n_obs = ZZ.n_rows;
  const arma::uword N    = particles.n_cols;

  arma::mat x1 = particles.rows(0, n_s - 1);
  arma::mat x2 = particles.rows(n_s, 2 * n_s - 1);

  arma::mat fitted = ZZ * (x1 + x2) + DD * shocks;

  const bool hxx_obs_nz = arma::any(arma::vectorise(arma::abs(hxx_obs)) > 0.0);
  const bool hxu_obs_nz = arma::any(arma::vectorise(arma::abs(hxu_obs)) > 0.0);
  const bool huu_obs_nz = arma::any(arma::vectorise(arma::abs(huu_obs)) > 0.0);
  if (hxx_obs_nz) fitted += 0.5 * hxx_obs * kron_cols_2(x1, x1);
  if (hxu_obs_nz) fitted += hxu_obs * kron_cols_2(shocks, x1);
  if (huu_obs_nz) fitted += 0.5 * huu_obs * kron_cols_2(shocks, shocks);

  arma::vec offset = d_obs + ghss_obs;
  fitted.each_col() += offset;

  const double sd_phi      = std::sqrt(me_variance / phi);
  const double log_sd_phi  = std::log(sd_phi);
  const double log2pi_half = 0.5 * std::log(2.0 * arma::datum::pi);
  const double const_term  = -(log2pi_half + log_sd_phi);
  const double inv_var     = 1.0 / (me_variance / phi);

  arma::vec log_w(N, arma::fill::zeros);
  for (arma::uword j = 0; j < n_obs; ++j) {
    arma::rowvec resid = arma::conv_to<arma::rowvec>::from(fitted.row(j));
    resid -= y_t(j);
    log_w += const_term -
             0.5 * inv_var * arma::conv_to<arma::vec>::from(resid % resid);
  }
  return log_w;
}


// [[Rcpp::export]]
List tpf_run_period_cpp(
    arma::mat       particles,     // (2*n_s) x N  (copied in, returned updated)
    const arma::vec y_t,           // n_obs
    const arma::mat L_e,           // n_e x n_e  lower Cholesky of Sigma_e
    const arma::mat hx,            // n_s x n_s
    const arma::mat hu,            // n_s x n_e
    const arma::mat hxx,           // n_s x (n_s*n_s)
    const arma::mat hxu,           // n_s x (n_s*n_e)
    const arma::mat huu,           // n_s x (n_e*n_e)
    const arma::vec hss,           // n_s
    const arma::mat ZZ,            // n_obs x n_s
    const arma::mat DD,            // n_obs x n_e
    const arma::vec d_obs,         // n_obs
    const arma::vec ghss_obs,      // n_obs
    const arma::mat hxx_obs,       // n_obs x (n_s*n_s), obs-row slice of ghxx
    const arma::mat hxu_obs,       // n_obs x (n_s*n_e), obs-row slice of ghxu
    const arma::mat huu_obs,       // n_obs x (n_e*n_e), obs-row slice of ghuu
    double          me_variance,   // > 0
    double          ess_target,    // default 0.5
    int             n_mh,          // mutation steps (default 1)
    int             max_stages,    // default 200
    Rcpp::Nullable<arma::mat> U_normals  = R_NilValue,
    Rcpp::Nullable<double>    U_resample = R_NilValue,
    Rcpp::Nullable<arma::vec> U_mid      = R_NilValue,
    Rcpp::Nullable<arma::mat> U_mutation = R_NilValue
) {
  // CPM args:
  //   U_normals  = n_e x N standard normals (NULL -> RNG, bit-identical).
  //   U_resample = z ~ N(0,1) for phi=1 resample uniform (NULL -> RNG).
  //   U_mid      = legacy K-vector z_k for mid-stage resamples (never fires; kept for compat).
  //   U_mutation = Option A mutation noise (n_2s+n_e+1) x (max_stages_u*n_mh*N).
  //                Column index: stage*(n_mh*N) + step*N + particle (0-based).
  //                Rows [0,n_2s): z_s, [n_2s,n_2s+n_e): z_e, n_2s+n_e: z_u -> log_u = log(Phi(z_u)).
  //                NULL -> RNG (bit-identical to pre-Option-A code).
  // Activate R's RNG stream so all rnorm/runif calls draw from R's RNG state
  Rcpp::RNGScope rng_scope;

  const arma::uword N    = particles.n_cols;
  const arma::uword n_e  = L_e.n_cols;
  const arma::uword n_2s = particles.n_rows;

  // ---- Step 1: draw shocks and compute full log-likelihoods (phi=1) --------
  // CPM: if U_normals supplied, use them; otherwise draw from R's RNG stream
  // (bit-identical to the pre-CPM path when U_normals = NULL).
  arma::mat z_mat(n_e, N);
  if (U_normals.isNotNull()) {
    z_mat = Rcpp::as<arma::mat>(U_normals);  // use pre-drawn standard normals
  } else {
    Rcpp::NumericVector rn_vec = Rcpp::rnorm((int)(n_e * N), 0.0, 1.0);
    // copy into z_mat (non-owning view is unsafe here since rn_vec goes out
    // of scope; use fill from pointer instead)
    std::copy(rn_vec.begin(), rn_vec.end(), z_mat.memptr());
  }
  arma::mat shocks = L_e * z_mat;  // n_e x N (theta-dependent, materialized)

  arma::vec log_liks = tpf_lw_inline(
      particles, y_t, ZZ, DD, shocks, d_obs, ghss_obs, me_variance, 1.0,
      hxx_obs, hxu_obs, huu_obs);

  // ---- Step 2: adaptive phi tempering loop ---------------------------------
  // CPM: extract pre-drawn z for the phi=1 resampling uniform.
  //   U_resample stores a standard normal z; the uniform is u = Phi(z).
  //   When U_resample = NULL, draw from R's RNG (standard non-CPM path).
  double u_phi1 = std::numeric_limits<double>::quiet_NaN();  // NaN = draw from RNG
  double z_resample_used = std::numeric_limits<double>::quiet_NaN();
  if (U_resample.isNotNull()) {
    double z_r = Rcpp::as<double>(U_resample);
    z_resample_used = z_r;
    // Transform N(0,1) -> Uniform(0,1) via standard normal CDF Phi.
    // R::pnorm(x, 0, 1, lower_tail=true, log_p=false) = Phi(x).
    u_phi1 = R::pnorm(z_r, 0.0, 1.0, 1, 0);
    // Clamp to (0,1) to avoid exact 0 or 1 from extreme z values
    u_phi1 = std::max(1e-15, std::min(1.0 - 1e-15, u_phi1));
  }

  double    phi_curr        = 0.0;
  arma::vec log_w(N, arma::fill::zeros);
  double    log_lik_contrib = 0.0;
  std::vector<double> phi_sched;
  phi_sched.reserve(max_stages + 1);

  // U_mid: legacy pre-allocated mid-stage resample z slots.
  // Kept for backward compatibility; never fires in practice because
  // next_phi_cpp maintains ESS >= target with a single phi=1 jump.
  arma::vec u_mid_vec;
  bool has_u_mid = U_mid.isNotNull();
  if (has_u_mid) u_mid_vec = Rcpp::as<arma::vec>(U_mid);
  int u_mid_idx = 0;   // next slot to consume (0-based)

  // Option A: U_mutation noise buffer for correlated RWMH proposals.
  // Layout: (n_2s+n_e+1) x (max_stages_u * n_mh * N) matrix.
  //   Column c = stage*(n_mh*N) + step*N + particle.
  //   Rows [0, n_2s):       z_s  ~ N(0,1) for state proposal
  //   Rows [n_2s, n_2s+n_e): z_e  ~ N(0,1) for shock proposal
  //   Row  n_2s+n_e:         log_u pre-stored as log(Uniform(0,1)) = log(pnorm(z_u))
  //                          where z_u ~ N(0,1) stored in the buffer.
  // When U_mutation is NULL, draw from R's RNG (bit-identical to pre-Option-A).
  // On exhaustion (col index >= n_cols) fall back to fresh draws and REPORT it
  // to R (see below); the fallback draws are statistically identical
  // N(0,1)/Uniform, so correctness is unaffected either way.
  arma::mat u_mut_mat;
  bool has_u_mut = U_mutation.isNotNull();
  arma::uword u_mut_ncols = 0;
  if (has_u_mut) {
    u_mut_mat   = Rcpp::as<arma::mat>(U_mutation);
    u_mut_ncols = u_mut_mat.n_cols;
  }
  // B5: these two conditions used to call Rcpp::warning() from inside the
  // kernel. Under options(warn = 2) R turns a warning into an ERROR, and R's
  // error mechanism is a longjmp: it unwinds PAST every C++ destructor on the
  // stack, including ~RNGScope, which is what writes the advanced RNG state
  // back to .Random.seed. The result was a silently stale .Random.seed (and
  // leaked Armadillo buffers) on exactly the runs a user asked to be strict.
  // So the kernel only COUNTS/records, and the R wrapper raises the warning
  // after the call returns -- where a longjmp is harmless.
  int u_mid_exhausted = 0;        // times the U_mid fallback fired
  int u_mut_need_col  = -1;       // first column index past the buffer, or -1

  for (int stage = 0; stage < max_stages; ++stage) {

    double phi_next  = next_phi_cpp(log_liks, phi_curr, ess_target, N);
    double delta_phi = phi_next - phi_curr;

    arma::vec log_w_new = log_w + delta_phi * log_liks;

    // Accumulate log normalisation constant
    log_lik_contrib += lse(log_w_new) - lse(log_w);

    log_w    = log_w_new;
    phi_curr = phi_next;
    phi_sched.push_back(phi_curr);

    // Normalise weights
    arma::vec lw_c  = log_w - arma::max(log_w);
    arma::vec w_norm = arma::exp(lw_c);
    w_norm /= arma::sum(w_norm);

    // ---- phi=1: ALWAYS resample (Invariant A) then break -------------------
    // CPM path (U_resample supplied): sorted systematic resampling.
    //   Particles are sorted by first state dimension (x1[0]) before resampling
    //   using the pre-drawn uniform u_phi1 = Phi(z_r). Sorting makes the
    //   particle-index map continuous in u_phi1, maximizing loglik correlation
    //   across CPM iterations (Deligiannidis et al. 2018, Section 2.3).
    //
    // Non-CPM path (U_resample=NULL, u_phi1=NaN): original unsorted resampling
    //   (T1 bit-identity gate). Distribution is unchanged (unbiased regardless
    //   of sort order).
    if (phi_curr >= 1.0 - 1e-10) {
      arma::uvec idx;
      if (!std::isnan(u_phi1)) {
        // CPM path: sorted resampling with injected uniform
        idx = systematic_resample_sorted_cpp(particles, w_norm, N, u_phi1);
      } else {
        // Non-CPM path: original unsorted resampling (bit-identical to pre-CPM)
        idx = systematic_resample_cpp(w_norm, N);
      }
      particles = particles.cols(idx);
      shocks    = shocks.cols(idx);
      break;
    }

    // ---- Mid-stage resample if ESS below threshold -------------------------
    bool resampled = false;
    double ess_val = ess_from_logw(log_w);
    if (ess_val < ess_target * (double)N) {
      arma::uvec idx;
      if (has_u_mid && u_mid_idx < (int)u_mid_vec.n_elem) {
        // Legacy CPM path: sorted systematic resample with pre-drawn z_k
        double z_k = u_mid_vec(u_mid_idx++);
        double u_k = R::pnorm(z_k, 0.0, 1.0, 1, 0);
        u_k = std::max(1e-15, std::min(1.0 - 1e-15, u_k));
        idx = systematic_resample_sorted_cpp(particles, w_norm, N, u_k);
      } else if (has_u_mid) {
        // U_mid supplied but K slots exhausted: fall back to fresh draw and
        // record it for the R wrapper to warn about (B5).
        ++u_mid_exhausted;
        idx = systematic_resample_cpp(w_norm, N);
      } else {
        // Non-CPM path: fresh RNG draw (bit-identical to pre-CPM code)
        idx = systematic_resample_cpp(w_norm, N);
      }
      particles = particles.cols(idx);
      shocks    = shocks.cols(idx);
      log_liks  = log_liks.elem(idx);
      log_w.fill(0.0);
      w_norm.fill(1.0 / (double)N);
      resampled = true;
    }

    // ---- RWMH mutation (Invariant B: only immediately after resample) ------
    // Herbst & Schorfheide (2019): mutate ONLY the period-t shock e_t with
    // the ancestor state s_{t-1} FIXED. An independence proposal e' ~
    // N(0, Sigma_e) (the shock prior, via L_e below) has acceptance ratio
    // exactly phi_curr * (loglik(e') - loglik(e)) -- the formula already in
    // use here. Moving the state s_{t-1} in a random walk (the old
    // behaviour) drops the filtering density p(s_{t-1}|Y_{1:t-1}) that
    // lives only in the resampled cloud, breaking invariance and biasing
    // the likelihood estimate upward. z_s is still drawn/consumed below
    // (RNG/CPM stream compatibility) but is deliberately unused.
    if (resampled && n_mh > 0 && N > 1) {

      // Per-particle RWMH loop
      // Option A: consume mutation noise from U_mutation buffer when supplied.
      // Column layout: stage*(n_mh*N) + step*N + particle (all 0-based).
      // Row layout: [0, n_2s) = z_s, [n_2s, n_2s+n_e) = z_e, n_2s+n_e = log_u.
      // The log_u row stores a pre-drawn N(0,1) z_u; log(u) = log(pnorm(z_u)).
      // Using N(0,1) z's (transformed to log-uniform via log(Phi(z_u))) preserves
      // the AR(1) CRN structure: correlated z_u across theta/theta' gives
      // correlated accept/reject decisions when the log-ratio is similar.
      for (arma::uword i = 0; i < N; ++i) {
        arma::vec s_i(particles.col(i));
        arma::vec e_i(shocks.col(i));

        arma::mat s_mat(s_i.memptr(), n_2s, 1, false);
        arma::mat e_mat(e_i.memptr(), n_e,  1, false);
        double tlp_i = phi_curr *
            tpf_lw_inline(s_mat, y_t, ZZ, DD, e_mat,
                           d_obs, ghss_obs, me_variance, 1.0,
                           hxx_obs, hxu_obs, huu_obs)(0);

        for (int step = 0; step < n_mh; ++step) {
          // Column index into U_mutation for (stage, step, particle i)
          arma::uword col_idx = (arma::uword)stage * (arma::uword)n_mh * N
                                + (arma::uword)step * N + i;
          bool use_u_mut = has_u_mut && (col_idx < u_mut_ncols);
          // Record the FIRST offending column only; the R wrapper decides how
          // loudly to report it (B5).
          if (has_u_mut && col_idx >= u_mut_ncols && u_mut_need_col < 0)
            u_mut_need_col = (int)col_idx;

          arma::vec z_sv(n_2s);
          arma::vec z_ev(n_e);
          double log_u;

          if (use_u_mut) {
            // Consume pre-drawn normals from U_mutation
            const arma::vec& col_v = u_mut_mat.col(col_idx);
            z_sv  = col_v.subvec(0, n_2s - 1);
            z_ev  = col_v.subvec(n_2s, n_2s + n_e - 1);
            // Row n_2s+n_e stores z_u ~ N(0,1); transform to log(Uniform(0,1))
            double z_u = col_v(n_2s + n_e);
            log_u = std::log(R::pnorm(z_u, 0.0, 1.0, 1, 0));
          } else {
            // Fallback: fresh draws from R's RNG (bit-identical to pre-Option-A)
            Rcpp::NumericVector z_s_r = Rcpp::rnorm((int)n_2s, 0.0, 1.0);
            std::copy(z_s_r.begin(), z_s_r.end(), z_sv.memptr());
            Rcpp::NumericVector z_e_r = Rcpp::rnorm((int)n_e, 0.0, 1.0);
            std::copy(z_e_r.begin(), z_e_r.end(), z_ev.memptr());
            log_u = std::log(Rcpp::as<double>(Rcpp::runif(1, 0.0, 1.0)));
          }

          // Ancestor state fixed (H&S 2019); z_sv drawn above but unused --
          // keeps the mutation-buffer layout / RNG stream bit-identical.
          arma::vec s_prop = s_i;
          arma::vec e_prop = L_e * z_ev;

          arma::mat sp_mat(s_prop.memptr(), n_2s, 1, false);
          arma::mat ep_mat(e_prop.memptr(), n_e,  1, false);
          double tlp_p = phi_curr *
              tpf_lw_inline(sp_mat, y_t, ZZ, DD, ep_mat,
                             d_obs, ghss_obs, me_variance, 1.0,
                             hxx_obs, hxu_obs, huu_obs)(0);

          if (std::isfinite(tlp_p) && log_u < tlp_p - tlp_i) {
            s_i   = s_prop;
            e_i   = e_prop;
            tlp_i = tlp_p;
          }
        }  // end mutation steps

        particles.col(i) = s_i;
        shocks.col(i)    = e_i;
      }  // end per-particle loop

      // Recompute full log-likelihoods after mutation
      log_liks = tpf_lw_inline(particles, y_t, ZZ, DD, shocks,
                                d_obs, ghss_obs, me_variance, 1.0,
                                hxx_obs, hxu_obs, huu_obs);
      log_w.fill(0.0);   // reset to uniform: temper from phi_curr next iter
    }
  }  // end phi loop

  // ---- Step 3: propagate (s_{t-1}, e_t) -> s_t ----------------------------
  arma::mat particles_new = tpf_propagate_particles(
      particles, shocks, hx, hu, hxx, hxu, huu, hss);

  // phi_schedule as arma::vec
  arma::vec phi_vec(phi_sched.size());
  for (arma::uword k = 0; k < phi_sched.size(); ++k)
    phi_vec(k) = phi_sched[k];

  // z_resample_used: the standard normal z used for the phi=1 resampling
  // uniform (u = Phi(z)). NaN when U_resample=NULL (non-CPM path) because
  // we do not recover the RNG-drawn uniform from the R stream. The R layer
  // stores NA for the non-CPM path and the actual z for the CPM path.
  Rcpp::NumericVector z_res_ret(1);
  z_res_ret[0] = z_resample_used;  // NaN -> NA in R

  return List::create(
      _["particles"]        = particles_new,
      _["log_lik_contrib"]  = log_lik_contrib,
      _["phi_schedule"]     = phi_vec,
      _["U_used"]           = z_mat,         // CPM: shock normals used (pre-L_e)
      _["z_resample_used"]  = z_res_ret,     // CPM: z for phi=1 uniform (NaN if RNG-drawn)
      _["u_mid_slots_used"] = u_mid_idx,     // CPM: number of legacy U_mid slots consumed
      // B5 warning channel: raised by the R wrapper, never from the kernel.
      _["u_mid_exhausted"]  = u_mid_exhausted,   // times the U_mid fallback fired
      _["u_mid_slots"]      = (int)(has_u_mid ? u_mid_vec.n_elem : 0),
      _["u_mut_need_col"]   = u_mut_need_col,    // first column past the buffer, -1 = none
      _["u_mut_have_cols"]  = (int)u_mut_ncols
  );
}


// [[Rcpp::export]]
arma::vec tpf_log_weights(
    const arma::mat& particles,   // (2*n_s) x N
    const arma::vec& y_t,         // n_obs
    const arma::mat& ZZ,          // n_obs x n_s
    const arma::mat& DD,          // n_obs x n_e (shock-to-obs loading)
    const arma::mat& shocks,      // n_e x N (current-period shock draws)
    const arma::vec& d_obs,       // n_obs (steady-state obs mean)
    const arma::vec& ghss_obs,    // n_obs (= 0.5 * ghss[obs_idx])
    double me_variance,           // scalar > 0 (nominal)
    double phi                    // tempering level in (0, 1]
) {
  const arma::uword n_s  = ZZ.n_cols;
  const arma::uword n_obs = ZZ.n_rows;
  const arma::uword N    = particles.n_cols;

  arma::mat x1 = particles.rows(0, n_s - 1);           // n_s x N
  arma::mat x2 = particles.rows(n_s, 2 * n_s - 1);     // n_s x N

  // Observation means: ZZ (x1 + x2) + DD e_t + d_obs + ghss_obs  (n_obs x N)
  arma::mat fitted = ZZ * (x1 + x2) + DD * shocks;  // n_obs x N
  // Add the constant offset (d_obs + ghss_obs broadcast over columns)
  arma::vec offset = d_obs + ghss_obs;
  fitted.each_col() += offset;

  // Tempered SD: sqrt(me_variance / phi)
  const double sd_phi     = std::sqrt(me_variance / phi);
  const double log_sd_phi = std::log(sd_phi);
  const double log2pi_half = 0.5 * std::log(2.0 * arma::datum::pi);

  // Log-weight for each particle: sum over observables of
  //   -0.5 log(2pi) - log(sd_phi) - 0.5 ((y_tj - mu_ij) / sd_phi)^2
  arma::vec log_w(N, arma::fill::zeros);
  const double const_term = -(log2pi_half + log_sd_phi);  // per observable
  const double inv_var    = 1.0 / (me_variance / phi);

  for (arma::uword j = 0; j < n_obs; ++j) {
    arma::rowvec resid = arma::conv_to<arma::rowvec>::from(fitted.row(j));
    resid -= y_t(j);
    log_w += const_term - 0.5 * inv_var * arma::conv_to<arma::vec>::from(resid % resid);
  }

  return log_w;
}
