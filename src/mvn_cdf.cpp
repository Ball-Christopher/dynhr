// mvn_cdf.cpp -- deterministic log multivariate-normal CDF for the PSKF.
//
// mvn_logcdf_cpp(b, C) returns log Phi_q(b; 0, C) = log P(X <= b),
// X ~ N(0, C) (upper bounds only: the PSKF's CSN normalising constants are
// orthant probabilities). Used by logcdf_ME_r() (R/pskf-likelihood.R, via
// .mvn_logcdf_sov) for 3 <= q <= miwa_qmax when pskf_cdf = "accurate".
// Also here (A3, 2026-09-29; at the end of the file): mvn_logcdf2_cpp, Genz's
// bivariate BVND for p >= 1e-3, mvn_logcdf3_cpp, Genz's exact trivariate
// TVN for p >= 1e-3 (A3b; the q = 3 blocks no longer reach the lattice
// there), and mvn_logcdf_dispatch_cpp, a bit-identical C++ replay of
// logcdf_ME_r's accurate-mode dispatch.
//
// Method (all on the LOG scale, so orthant tails down to log p ~ -1e3 keep
// their RELATIVE accuracy):
//   1. Genz (1992) separation of variables (SOV) with the Gibson-Glasbey-
//      Elston / Genz-Bretz variable reordering: at each Cholesky step the
//      remaining variable with the smallest conditional probability, given
//      the truncated means of the variables before it, goes next.
//   2. Botev (2017, JRSS-B 79:125-148) minimax exponential tilting: the
//      conditional truncated normals are drawn with shifted means mu_k, mu
//      solving the saddle-point equations grad psi(x, mu) = 0 (damped
//      Newton), which makes the SOV integrand nearly constant even in deep
//      tails (without it: errors up to 4.4 nats with 8000 lattice points
//      on the W75 development set). mu is a smooth function of (b, C).
//   3. Rank-1 lattice rules (prime N; generating vectors by a CBC
//      construction for the Korobov alpha = 4 criterion with weights 1/j,
//      computed once offline), Sidi's sin^2 periodising transform
//      psi(t) = t - sin(2 pi t) / (2 pi), and kShifts FIXED shifts. No RNG:
//      identical inputs give bit-identical output. The spread of the kShifts
//      shifted-rule estimates is the error estimate (relative error of p =
//      absolute error of log p).
//   4. Lattice ladder N = 61, 127, ..., 32749 (x kShifts), from a
//      dimension-dependent start level. Level l's estimate v_l is used with
//      weight a(e_l / tau), tau = max(abs_tol, rel_tol |v_l|), a = 1 for
//      e_l <= kRlo tau, 0 for e_l >= kRhi tau, smoothstep in log(e_l)
//      between; the remaining weight goes to the next level:
//        R_l = a v_l + (1 - a) R_{l+1},  R_last = v_last.
//      Changing the lattice size is therefore CONTINUOUS (C1) in (b, C):
//      e_l and v_l are continuous, so no jump appears when a nearby theta
//      needs one more level. The variable reordering (step 1) is a discrete
//      choice: where it switches, the value jumps by the difference of two
//      accepted rules' errors (each <= ~kRhi tau; measured in the W75
//      tests).
//
//   5. A memo of the last kMemo calls keyed by the exact input bits (the
//      PSKF repeats ~40% of its calls within a period); a hit returns
//      exactly what recomputing would.
//
// Accuracy (W75 tests, independent oracles): |log p error| <
// max(1e-5, 1e-7 |log p|) on 200 random problems, q = 3..7, log p down to
// -300 (worst 0.39 of the tolerance). Cost on recorded PSKF calls (-O2):
// ~0.1 ms at q = 3, ~0.3 ms at q = 4, ~1.4 ms at q = 5, ~3 ms at q = 6,
// ~10 ms at q = 7.
//
// Returns c(log p, error estimate (absolute, log scale) of the last level
// evaluated, index of the last level evaluated (0-based), tilt converged
// (1 / 0)). log p is NA when a Cholesky pivot is <= 1e-10 (numerically
// singular C in correlation form) or q > kMaxDim + 1 -- the R caller then
// falls back.

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cfloat>

namespace {

const int kShifts = 4;
const int kMaxDim = 12;     // integrand dimension q - 1  =>  q <= 13
const int kLevels = 10;
const int kLatN[kLevels] = {61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381,
                            32749};
// CBC generating vectors (Korobov space alpha = 4, product weights 1 / j),
// computed offline; first component 1.
const int kLatZ[kLevels][kMaxDim] = {
  {1, 18, 28, 13, 8, 23, 3, 25, 23, 5, 3, 10},
  {1, 29, 53, 36, 48, 21, 34, 43, 60, 12, 55, 5},
  {1, 70, 19, 97, 119, 11, 5, 105, 32, 73, 43, 111},
  {1, 151, 232, 244, 87, 15, 209, 240, 203, 199, 124, 172},
  {1, 374, 156, 285, 342, 141, 260, 325, 332, 250, 369, 22},
  {1, 598, 182, 699, 449, 403, 287, 620, 681, 429, 957, 173},
  {1, 1715, 452, 289, 799, 1457, 1409, 525, 1381, 109, 1513, 1358},
  {1, 2431, 3799, 1141, 1651, 2571, 2362, 557, 568, 528, 978, 3175},
  {1, 6208, 959, 4195, 2660, 6310, 3163, 4567, 1900, 5670, 33, 1209},
  {1, 12764, 15397, 3099, 6829, 4296, 1258, 8427, 3286, 11584, 3609, 1807}
};
// Fixed shifts: set.seed(20260926); matrix(runif(8 * 12), 8, 12)[1:4, ]
// (R's Mersenne-Twister), hard-coded so the rule is the same everywhere.
// Unstructured on purpose: Richtmyer shifts frac(s * sqrt(p_k)) share the
// lattice's error pattern -- their spread under-estimated the error 65x on
// a W75 test problem (q = 4, 1.5e-5 nat error reported as 2.3e-7).
const double kShift[kShifts][kMaxDim] = {
  {0.84941126219928265, 0.19242514367215335, 0.26801304169930518,
   0.98639770364388824, 0.84860375942662358, 0.78368994500488043,
   0.26797405071556568, 0.49743939703330398, 0.92218650551512837,
   0.086382880108430982, 0.022812416544184089, 0.97129672206938267},
  {0.20063198101706803, 0.75948845641687512, 0.97152024996466935,
   0.39623415842652321, 0.82542240247130394, 0.55471664411015809,
   0.69288962730206549, 0.24464827147312462, 0.7811759680043906,
   0.29497016500681639, 0.70350086502730846, 0.95969398971647024},
  {0.78024318837560713, 0.92149605858139694, 0.037294144509360194,
   0.72300614486448467, 0.21912278328090906, 0.90469539328478277,
   0.94430848397314548, 0.13848731992766261, 0.14404741721227765,
   0.91972178081050515, 0.56904939445666969, 0.72569805174134672},
  {0.82569045177660882, 0.77758033177815378, 0.73679690132848918,
   0.018687590723857284, 0.44339623698033392, 0.31245535938069224,
   0.58903712104074657, 0.7342048182617873, 0.91297554899938405,
   0.47039325931109488, 0.55053859716281295, 0.84040472842752934}
};
// First lattice level by q = 3, 4, 5, 6, 7 (q >= 8: 7) -- the level most
// PSKF calls of that dimension end at; only a speed choice (the ladder
// climbs from there as needed).
const int kStartLevel[5] = {0, 2, 3, 5, 6};
// Acceptance band for the error estimate, as a fraction of the tolerance.
// W75 development set (200 problems, q = 3..7): max |error| / tolerance
// 0.46 with (0.1, 0.2).
const double kRlo = 0.1, kRhi = 0.2;

const double kLogSqrt2Pi = 0.918938533204672741780329736406;
const double kPi = 3.141592653589793238462643383280;

inline double lPhi(double x) { return R::pnorm(x, 0.0, 1.0, 1, 1); }
inline double lphi(double x) { return -0.5 * x * x - kLogSqrt2Pi; }
// phi(x) / Phi(x), stable in the lower tail; 0 at +Inf
inline double mills(double x) {
  if (x == R_PosInf) return 0.0;
  return std::exp(lphi(x) - lPhi(x));
}

// Sidi's sin^2 transform of t in [0, 1): log w for w = psi(t) and
// log psi'(t) = log(2 sin^2(pi t)); accurate at both ends (series for small
// t, psi(1 - t) = 1 - psi(t)).
inline double sidi_psi_half(double t) {  // t in [0, 0.5]
  if (t < 0.01) {
    const double x = 2.0 * kPi * t, x2 = x * x;
    return x * x2 / 6.0 *
           (1.0 - x2 / 20.0 * (1.0 - x2 / 42.0 * (1.0 - x2 / 72.0))) /
           (2.0 * kPi);
  }
  return t - std::sin(2.0 * kPi * t) / (2.0 * kPi);
}
inline void sidi(double t, double& logw, double& logjac) {
  const double tm = (t <= 0.5) ? t : 1.0 - t;
  logjac = M_LN2 + 2.0 * std::log(std::sin(kPi * tm));
  const double ps = sidi_psi_half(tm);
  logw = (t <= 0.5) ? std::log(ps) : std::log1p(-ps);
}

// Gaussian elimination with partial pivoting (n x n, row-major); false if
// singular.
bool gsolve(std::vector<double>& A, std::vector<double>& r, int n) {
  for (int k = 0; k < n; ++k) {
    int p = k;
    double mx = std::fabs(A[k * n + k]);
    for (int i = k + 1; i < n; ++i)
      if (std::fabs(A[i * n + k]) > mx) { mx = std::fabs(A[i * n + k]); p = i; }
    if (!(mx > 0.0) || !std::isfinite(mx)) return false;
    if (p != k) {
      for (int j = 0; j < n; ++j) std::swap(A[k * n + j], A[p * n + j]);
      std::swap(r[k], r[p]);
    }
    for (int i = k + 1; i < n; ++i) {
      const double f = A[i * n + k] / A[k * n + k];
      if (f == 0.0) continue;
      for (int j = k; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
      r[i] -= f * r[k];
    }
  }
  for (int i = n - 1; i >= 0; --i) {
    double s = r[i];
    for (int j = i + 1; j < n; ++j) s -= A[i * n + j] * r[j];
    r[i] = s / A[i * n + i];
  }
  return true;
}

// Gradient (and optionally Jacobian) of Botev's psi(x, mu) for upper bounds
// only:  psi = sum_k log Phi(u_k - mu_k - c_k) + mu_k^2 / 2 - x_k mu_k,
// c = Lt x, with Lt strictly lower triangular (rows scaled by the Cholesky
// diagonal; d x d row-major), u the scaled bounds, v = (x_1..x_m,
// mu_1..mu_m), m = d - 1 (x_d = mu_d = 0).
void tilt_grad(const std::vector<double>& Lt, const std::vector<double>& u,
               int d, const std::vector<double>& v, std::vector<double>& g,
               std::vector<double>* J) {
  const int m = d - 1;
  std::vector<double> P(d), dP(d);
  for (int k = 0; k < d; ++k) {
    double c = 0.0;
    for (int j = 0; j < std::min(k, m); ++j) c += Lt[k * d + j] * v[j];
    const double mu = (k < m) ? v[m + k] : 0.0;
    const double ut = u[k] - mu - c;
    const double lam = mills(ut);
    P[k]  = -lam;                                    // d psi_k / d ut ...
    dP[k] = (lam == 0.0) ? 0.0 : -lam * (lam + ut);  // ... and d P_k / d mu_k
  }
  g.assign(2 * m, 0.0);
  for (int j = 0; j < m; ++j) {
    double s = 0.0;
    for (int k = j + 1; k < d; ++k) s += Lt[k * d + j] * P[k];
    g[j]     = -v[m + j] + s;
    g[m + j] = v[m + j] - v[j] + P[j];
  }
  if (J) {
    const int n = 2 * m;
    J->assign(n * n, 0.0);
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < m; ++i) {
        double s = 0.0;
        for (int k = std::max(i, j) + 1; k < d; ++k)
          s += Lt[k * d + j] * dP[k] * Lt[k * d + i];
        (*J)[j * n + i] = s;
        const double mx_ij = (i == j ? -1.0 : 0.0) +
                             (i > j ? dP[i] * Lt[i * d + j] : 0.0);
        const double mx_ji = (i == j ? -1.0 : 0.0) +
                             (j > i ? dP[j] * Lt[j * d + i] : 0.0);
        (*J)[j * n + (m + i)] = mx_ij;
        (*J)[(m + j) * n + i] = mx_ji;
      }
      (*J)[(m + j) * n + (m + j)] = 1.0 + dP[j];
    }
  }
}

// One lattice level: log-mean-exp of the tilted SOV integrand over each of
// the kShifts shifted rules; returns the pooled estimate (log p) and the
// relative standard error of p across the shifts.
void lattice_level(int level, int d, const std::vector<double>& Lt,
                   const std::vector<double>& u, const std::vector<double>& mu,
                   double& est, double& err) {
  const int m = d - 1;
  const int N = kLatN[level];
  const int* z = kLatZ[level];
  std::vector<double> lw(N), ls(kShifts), y(d);
  double tilt_const = 0.0;
  for (int k = 0; k < m; ++k) tilt_const += 0.5 * mu[k] * mu[k];
  for (int s = 0; s < kShifts; ++s) {
    double mx = R_NegInf;
    for (int i = 0; i < N; ++i) {
      double acc = tilt_const;
      for (int k = 0; k < d; ++k) {
        double c = 0.0;
        for (int j = 0; j < k; ++j) c += Lt[k * d + j] * y[j];
        const double ut = u[k] - mu[k] - c;
        const double lp = lPhi(ut);
        acc += lp;
        if (k < m) {
          double t = (double) (((long long) i * z[k]) % N) / N + kShift[s][k];
          t -= std::floor(t);
          double logw, logjac;
          sidi(t, logw, logjac);
          acc += logjac;
          if (!(logw > R_NegInf) || !(acc > R_NegInf)) { acc = R_NegInf; break; }
          // N(mu_k, 1) truncated to (-Inf, u_k - c_k]: mu_k + Phi^{-1}(w
          // Phi(ut)), on the log scale
          y[k] = mu[k] + R::qnorm(logw + lp, 0.0, 1.0, 1, 1);
          acc -= mu[k] * y[k];
        }
      }
      lw[i] = acc;
      if (acc > mx) mx = acc;
    }
    double sm = 0.0;
    if (mx > R_NegInf)
      for (int i = 0; i < N; ++i) sm += std::exp(lw[i] - mx);
    ls[s] = (mx > R_NegInf) ? mx + std::log(sm / N) : R_NegInf;
  }
  const double mx = *std::max_element(ls.begin(), ls.end());
  if (!(mx > R_NegInf) || !std::isfinite(mx)) {
    est = mx;
    err = R_PosInf;
    return;
  }
  double sm = 0.0, sq = 0.0;
  for (int s = 0; s < kShifts; ++s) {
    const double e = std::exp(ls[s] - mx);
    sm += e;
    sq += e * e;
  }
  const double mean = sm / kShifts;
  double var = (sq / kShifts - mean * mean) * kShifts / (kShifts - 1.0);
  if (var < 0.0) var = 0.0;
  est = mx + std::log(mean);
  err = std::sqrt(var / kShifts) / mean;
}

// Memo of the last kMemo calls, keyed by the exact input bits. The PSKF
// repeats ~40% of its calls within a period (e.g. the post-prune
// compensation's normaliser IS the next likelihood denominator); the
// evaluator is deterministic, so a hit returns exactly what recomputing
// would.
const int kMemo = 16;
struct MemoEntry {
  std::vector<double> key;
  double out[4];
};
std::vector<MemoEntry> g_memo(kMemo);
int g_memo_next = 0;

bool memo_find(const std::vector<double>& key, double* out) {
  for (const MemoEntry& e : g_memo) {
    if (e.key.size() == key.size() && !key.empty() &&
        std::memcmp(e.key.data(), key.data(), key.size() * sizeof(double)) == 0) {
      std::copy(e.out, e.out + 4, out);
      return true;
    }
  }
  return false;
}

void memo_store(const std::vector<double>& key, const double* out) {
  MemoEntry& e = g_memo[g_memo_next];
  e.key = key;
  std::copy(out, out + 4, e.out);
  g_memo_next = (g_memo_next + 1) % kMemo;
}

Rcpp::NumericVector mvn_logcdf_impl(Rcpp::NumericVector b_in,
                                    Rcpp::NumericMatrix C_in,
                                    double abs_tol, double rel_tol,
                                    int start_level);

}  // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector mvn_logcdf_cpp(Rcpp::NumericVector b_in,
                                   Rcpp::NumericMatrix C_in,
                                   double abs_tol = 1e-5,
                                   double rel_tol = 1e-7,
                                   int start_level = -1) {
  const int d = b_in.size();
  if (C_in.nrow() != d || C_in.ncol() != d)
    Rcpp::stop("mvn_logcdf_cpp: C must be a %d x %d matrix", d, d);
  std::vector<double> key;
  key.reserve(4 + d + d * d);
  key.push_back(d);
  key.push_back(abs_tol);
  key.push_back(rel_tol);
  key.push_back(start_level);
  key.insert(key.end(), b_in.begin(), b_in.end());
  key.insert(key.end(), C_in.begin(), C_in.end());
  double out[4];
  if (memo_find(key, out))
    return Rcpp::NumericVector(out, out + 4);
  Rcpp::NumericVector res = mvn_logcdf_impl(b_in, C_in, abs_tol, rel_tol,
                                            start_level);
  std::copy(res.begin(), res.end(), out);
  memo_store(key, out);
  return res;
}

namespace {

Rcpp::NumericVector mvn_logcdf_impl(Rcpp::NumericVector b_in,
                                    Rcpp::NumericMatrix C_in,
                                    double abs_tol, double rel_tol,
                                    int start_level) {
  const double NA = NA_REAL;
  const int d = b_in.size();
  if (d == 0) return Rcpp::NumericVector::create(0.0, 0.0, 0.0, 1.0);
  if (d - 1 > kMaxDim) return Rcpp::NumericVector::create(NA, NA, NA, NA);

  // ---- correlation form ---------------------------------------------------
  std::vector<double> b(d), C(d * d), sd(d);
  for (int i = 0; i < d; ++i) {
    if (ISNAN(b_in[i]) || !(C_in(i, i) > 0.0) || !std::isfinite(C_in(i, i)))
      return Rcpp::NumericVector::create(NA, NA, NA, NA);
    sd[i] = std::sqrt(C_in(i, i));
    b[i] = b_in[i] / sd[i];
  }
  for (int i = 0; i < d; ++i)
    if (b[i] == R_NegInf)
      return Rcpp::NumericVector::create(R_NegInf, 0.0, 0.0, 1.0);
  for (int i = 0; i < d; ++i)
    for (int j = 0; j < d; ++j) {
      const double cij = 0.5 * (C_in(i, j) + C_in(j, i));
      if (!std::isfinite(cij))
        return Rcpp::NumericVector::create(NA, NA, NA, NA);
      C[i * d + j] = (i == j) ? 1.0 : cij / (sd[i] * sd[j]);
    }

  // ---- pivoted Cholesky with Genz-Bretz / GGE reordering ------------------
  std::vector<double> L(d * d, 0.0), ymean(d, 0.0);
  for (int i = 0; i < d; ++i) {
    int jbest = i;
    double lbest = R_PosInf;
    for (int j = i; j < d; ++j) {
      double s2 = C[j * d + j], mj = 0.0;
      for (int h = 0; h < i; ++h) {
        s2 -= L[j * d + h] * L[j * d + h];
        mj += L[j * d + h] * ymean[h];
      }
      const double lp = (s2 > 0.0) ? lPhi((b[j] - mj) / std::sqrt(s2))
                                   : ((b[j] - mj >= 0.0) ? 0.0 : R_NegInf);
      if (lp < lbest) { lbest = lp; jbest = j; }
    }
    if (jbest != i) {
      std::swap(b[i], b[jbest]);
      for (int k = 0; k < d; ++k) std::swap(C[i * d + k], C[jbest * d + k]);
      for (int k = 0; k < d; ++k) std::swap(C[k * d + i], C[k * d + jbest]);
      for (int k = 0; k < d; ++k) std::swap(L[i * d + k], L[jbest * d + k]);
    }
    double lii2 = C[i * d + i];
    for (int h = 0; h < i; ++h) lii2 -= L[i * d + h] * L[i * d + h];
    if (!(lii2 > 1e-10)) return Rcpp::NumericVector::create(NA, NA, NA, NA);
    L[i * d + i] = std::sqrt(lii2);
    for (int k = i + 1; k < d; ++k) {
      double s = C[k * d + i];
      for (int h = 0; h < i; ++h) s -= L[k * d + h] * L[i * d + h];
      L[k * d + i] = s / L[i * d + i];
    }
    double mi = 0.0;
    for (int h = 0; h < i; ++h) mi += L[i * d + h] * ymean[h];
    ymean[i] = -mills((b[i] - mi) / L[i * d + i]);   // E[Z | Z <= a]
  }
  // scaled bounds and the strictly lower, row-scaled factor
  std::vector<double> u(d), Lt(d * d, 0.0);
  for (int k = 0; k < d; ++k) {
    u[k] = b[k] / L[k * d + k];
    for (int j = 0; j < k; ++j) Lt[k * d + j] = L[k * d + j] / L[k * d + k];
  }
  if (d == 1) return Rcpp::NumericVector::create(lPhi(u[0]), 0.0, 0.0, 1.0);

  // ---- minimax tilting: damped Newton on grad psi = 0 from 0 --------------
  const int m = d - 1;
  std::vector<double> v(2 * m, 0.0), g, J, g2, vn, dlt(2 * m);
  double converged = 0.0;
  for (int it = 0; it < 100; ++it) {
    tilt_grad(Lt, u, d, v, g, &J);
    double gn = 0.0;
    for (double gi : g) gn += gi * gi;
    if (!std::isfinite(gn)) break;
    if (std::sqrt(gn) < 1e-11) { converged = 1.0; break; }
    for (int i = 0; i < 2 * m; ++i) dlt[i] = -g[i];
    if (!gsolve(J, dlt, 2 * m)) break;
    double t = 1.0;
    bool ok = false;
    while (t >= 1e-10) {
      vn = v;
      for (int i = 0; i < 2 * m; ++i) vn[i] += t * dlt[i];
      tilt_grad(Lt, u, d, vn, g2, nullptr);
      double g2n = 0.0;
      for (double gi : g2) g2n += gi * gi;
      if (std::isfinite(g2n) && g2n <= (1.0 - 1e-4 * t) * gn) { ok = true; break; }
      t *= 0.5;
    }
    if (!ok) break;
    v = vn;
  }
  // no convergence: the untilted Genz SOV (mu = 0) -- still consistent; the
  // lattice ladder absorbs its larger error
  std::vector<double> mu(d, 0.0);
  if (converged == 1.0)
    for (int k = 0; k < m; ++k) mu[k] = v[m + k];

  // ---- lattice ladder with continuous hand-over ---------------------------
  int lev = start_level;
  if (lev < 0) lev = (d <= 7) ? kStartLevel[d - 3 < 0 ? 0 : d - 3] : 7;
  if (lev > kLevels - 1) lev = kLevels - 1;
  const double llo = std::log(kRlo), lhi = std::log(kRhi);
  double result = 0.0, weight_left = 1.0, est = NA, err = NA;
  for (; lev < kLevels; ++lev) {
    lattice_level(lev, d, Lt, u, mu, est, err);
    if (!std::isfinite(est)) {
      if (est == R_NegInf)
        return Rcpp::NumericVector::create(R_NegInf, err, (double) lev, converged);
      return Rcpp::NumericVector::create(NA, NA, (double) lev, converged);
    }
    double a = 1.0;
    if (lev < kLevels - 1) {
      const double tau = std::max(abs_tol, rel_tol * std::fabs(est));
      const double lr = std::log(err / tau);
      if (lr >= lhi) {
        a = 0.0;
      } else if (lr > llo) {
        const double x = (lhi - lr) / (lhi - llo);
        a = x * x * (3.0 - 2.0 * x);
      }
    }
    result += weight_left * a * est;
    weight_left *= (1.0 - a);
    if (weight_left == 0.0) break;
  }
  if (lev > kLevels - 1) lev = kLevels - 1;
  return Rcpp::NumericVector::create(result, err, (double) lev, converged);
}

}  // namespace


// ---------------------------------------------------------------------------
// mvn_logcdf2_cpp(h1, h2, rho): log Phi_2(h1, h2; rho) for STANDARDISED
// bounds by Genz's (2004, Statistics and Computing 14:251-260) BVND -- the
// Drezner-Wesolowsky (1990) Gauss-Legendre rule with Genz's |rho| >= 0.925
// expansion, the algorithm mvtnorm uses for q = 2. Its ABSOLUTE error in p
// is ~1e-16 (A3, 2026-09-29: max |p - mvtnorm::pmvnorm| 5.6e-17 on 200
// random problems), so log p is accurate to ~1e-16 / p RELATIVE. The value
// is therefore returned only when p >= kBvnMinP = 1e-3 (log p error
// <= ~1e-13) and NA otherwise; the R caller (.mvn_logcdf2) then runs its
// log-scale quadrature, which keeps relative accuracy in deep tails.
// Measured against that quadrature on 8000 random problems with p >= 1e-3
// (|rho| up to 0.9999 and down to 1e-12): max |diff in log p| 6.6e-14, the
// quadrature's own rel.tol level. ~0.1 us per call against ~90 us for the
// quadrature (the PSKF makes ~1200 such calls on a T = 200 Reiter-HANK
// likelihood).
namespace {

const double kBvnMinP = 1e-3;
// Gauss-Legendre half-rules (nodes in [-1, 0) with weights), N = 6, 12, 20.
const double kBvnW[3][10] = {
  {0.1713244923791705, 0.3607615730481384, 0.4679139345726904},
  {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
   0.2031674267230659, 0.2334925365383547, 0.2491470458134029},
  {0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
   0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
   0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
   0.1527533871307259}};
const double kBvnX[3][10] = {
  {-0.9324695142031522, -0.6612093864662647, -0.2386191860831970},
  {-0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
   -0.5873179542866171, -0.3678314989981802, -0.1252334085114692},
  {-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
   -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
   -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
   -0.07652652113349733}};

inline double Phi1(double x) { return R::pnorm(x, 0.0, 1.0, 1, 0); }

// P(X > dh, Y > dk), corr(X, Y) = r, |r| <= 1 (Genz's BVND).
double bvnu(double dh, double dk, double r) {
  int ng, lg;
  if (std::fabs(r) < 0.3) { ng = 0; lg = 3; }
  else if (std::fabs(r) < 0.75) { ng = 1; lg = 6; }
  else { ng = 2; lg = 10; }
  double h = dh, k = dk, hk = h * k, bvn = 0.0;
  if (std::fabs(r) < 0.925) {
    const double hs = (h * h + k * k) / 2.0, asr = std::asin(r);
    for (int i = 0; i < lg; ++i) {
      double sn = std::sin(asr * (kBvnX[ng][i] + 1.0) / 2.0);
      bvn += kBvnW[ng][i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
      sn = std::sin(asr * (-kBvnX[ng][i] + 1.0) / 2.0);
      bvn += kBvnW[ng][i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
    }
    bvn = bvn * asr / (4.0 * kPi) + Phi1(-h) * Phi1(-k);
  } else {
    if (r < 0.0) { k = -k; hk = -hk; }
    if (std::fabs(r) < 1.0) {
      const double as = (1.0 - r) * (1.0 + r);
      double a = std::sqrt(as);
      const double bs = (h - k) * (h - k);
      const double c = (4.0 - hk) / 8.0, d = (12.0 - hk) / 16.0;
      bvn = a * std::exp(-(bs / as + hk) / 2.0) *
            (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 +
             c * d * as * as / 5.0);
      if (hk > -160.0) {
        const double b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2.0) * std::sqrt(2.0 * kPi) * Phi1(-b / a) * b *
               (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }
      a /= 2.0;
      for (int i = 0; i < lg; ++i) {
        for (int s = 0; s < 2; ++s) {
          const double xi = (s == 0) ? kBvnX[ng][i] : -kBvnX[ng][i];
          const double xs = (a * (xi + 1.0)) * (a * (xi + 1.0));
          const double rs = std::sqrt(1.0 - xs);
          const double asr = -(bs / xs + hk) / 2.0;
          if (asr > -100.0)
            bvn += a * kBvnW[ng][i] * std::exp(asr) *
                   (std::exp(-hk * xs / (2.0 * (1.0 + rs) * (1.0 + rs))) / rs -
                    (1.0 + c * xs * (1.0 + d * xs)));
        }
      }
      bvn = -bvn / (2.0 * kPi);
    }
    if (r > 0.0) {
      bvn += Phi1(-std::max(h, k));
    } else {
      bvn = -bvn;
      if (k > h) {
        if (h < 0.0) bvn += Phi1(k) - Phi1(h);
        else bvn += Phi1(-h) - Phi1(-k);
      }
    }
  }
  return bvn;
}

}  // namespace

// [[Rcpp::export(rng = false)]]
double mvn_logcdf2_cpp(double h1, double h2, double rho) {
  if (ISNAN(h1) || ISNAN(h2) || ISNAN(rho) || !std::isfinite(h1) ||
      !std::isfinite(h2) || !(std::fabs(rho) <= 1.0))
    return NA_REAL;
  const double p = bvnu(-h1, -h2, rho);
  if (!(p >= kBvnMinP) || !(p <= 1.0 + 1e-12)) return NA_REAL;
  return std::log(std::min(p, 1.0));
}


// ---------------------------------------------------------------------------
// mvn_logcdf3_cpp(b, C): log Phi_3(b; 0, C) = log P(X1 <= b1, X2 <= b2,
// X3 <= b3) by Genz's (2004, Statistics and Computing 14:251-260) TVN
// method: Plackett's (1954) reduction to a 1-D integral of bivariate terms,
// integrated by adaptive Gauss-Kronrod. Written from the paper (not from
// TVPACK's source). Bounds and correlations are standardised here
// (h_i = b_i / sqrt(C_ii), r_ij = C_ij / sqrt(C_ii C_jj)).
//
//   The variables are relabelled so that |r23| is the largest correlation;
//   r12 and r13 are then moved along r12(x) = sin(x asin r12), r13(x) =
//   sin(x asin r13), x in [0, 1], with r23 fixed. At x = 0 the probability
//   factors, Phi(h1) Phi_2(h2, h3; r23) (Genz's BVND above); Plackett's
//   identity dPhi_3 / dr_ij = phi_2(h_i, h_j; r_ij) Phi(conditional bound of
//   the third variable) gives
//     Phi_3 = Phi(h1) Phi_2(h2, h3; r23)
//           + (1 / 2 pi) int_0^1 sum_{j = 2, 3} asin(r1j)
//               exp(-(h1^2 - 2 r h1 hj + hj^2) / (2 cos^2 t)) Phi(B_j(x)) dx,
//   t = x asin r1j, r = sin t (the substitution r = sin t cancels the
//   1 / sqrt(1 - r^2) of phi_2), B_j = the standardised bound of the third
//   variable given X1 = h1, Xj = hj under the correlations at x. In the
//   angle parametrisation r_ij = cos(phi_ij) the path is a straight line in
//   phi from (pi/2, pi/2, phi23) to the target, and the positive-
//   semidefinite set is convex in phi, so every matrix on the path is
//   positive definite.
//   The integral is computed by global adaptive (7, 15) Gauss-Kronrod
//   (bisect the interval with the largest |K15 - G7|) until the summed
//   |K15 - G7| / (2 pi) -- a bound on the error in p -- is <= kTvnTol.
//
// Threshold rule (the SAME function serves the R dispatch, .mvn_logcdf_sov,
// and the C++ dispatch below, so the rule is identical by construction).
// The value is returned only when
//   (a) every input is finite and every C_ii > 0;
//   (b) the correlation matrix is safely nonsingular:
//       det = 1 - r12^2 - r13^2 - r23^2 + 2 r12 r13 r23 >= kTvnMinDet;
//   (c) the adaptive rule reached its tolerance within kTvnMaxIntervals;
//   (d) p >= kTvnMinP = 1e-3 (as for the bivariate rule: an absolute error
//       ~1e-15 in p is <= ~1e-12 in log p there; deeper tails need the
//       lattice evaluator's log-scale RELATIVE accuracy).
// Otherwise NA, and the caller runs the lattice evaluator mvn_logcdf_cpp.
// A3b (2026-09-29): replaces the lattice's ~3.6e-7-per-call error (its
// max(1e-5, 1e-7 |log p|) tolerance) on the ~600 q = 3 calls of a T = 200
// Reiter-HANK PSKF likelihood.
namespace {

const double kTvnMinP = 1e-3;
const double kTvnMinDet = 1e-8;
const double kTvnTol = 1e-15;
const int kTvnMaxIntervals = 100;

// (7, 15) Gauss-Kronrod: Kronrod nodes in (0, 1) descending, then 0, with
// weights; the Gauss-7 weights belong to kGK15X[1], [3], [5] and 0.
const double kGK15X[8] = {
  0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
  0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
  0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
  0.207784955007898467600689403773245, 0.0};
const double kGK15W[8] = {
  0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
  0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
  0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
  0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
const double kG7W[4] = {
  0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
  0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

struct TvnPath { double h1, h2, h3, r23, a12, a13; };

// d Phi_3 / d theta for the pair (a, b) at r_ab = sin(theta), without the
// 1 / (2 pi): exp(-quadratic form / 2) Phi(conditional bound of c), with
// ra = r_ac, rb = r_bc, rr = cos^2(theta) = 1 - r^2.
inline double tvn_pair(double ba, double bb, double bc, double ra, double rb,
                       double r, double rr) {
  // (1 - r^2) det(R)
  const double dt = rr * (rr - (ra - rb) * (ra - rb) - 2.0 * ra * rb * (1.0 - r));
  if (!(dt > 0.0)) return 0.0;
  const double ft = (ba - r * bb) * (ba - r * bb) / rr + bb * bb;
  const double bt = (bc * rr + ba * (r * rb - ra) + bb * (r * ra - rb)) /
                    std::sqrt(dt);
  return std::exp(-0.5 * ft) * Phi1(bt);
}

inline double tvn_integrand(double x, const TvnPath& P) {
  const double t12 = P.a12 * x, t13 = P.a13 * x;
  const double s12 = std::sin(t12), c12 = std::cos(t12);
  const double s13 = std::sin(t13), c13 = std::cos(t13);
  double f = 0.0;
  if (P.a12 != 0.0)
    f += P.a12 * tvn_pair(P.h1, P.h2, P.h3, s13, P.r23, s12, c12 * c12);
  if (P.a13 != 0.0)
    f += P.a13 * tvn_pair(P.h1, P.h3, P.h2, s12, P.r23, s13, c13 * c13);
  return f;
}

// (7, 15) Gauss-Kronrod on [lo, hi]: returns K15, sets err = |K15 - G7|
double tvn_gk15(double lo, double hi, const TvnPath& P, double& err) {
  const double c = 0.5 * (lo + hi), hw = 0.5 * (hi - lo);
  const double fc = tvn_integrand(c, P);
  double k = kGK15W[7] * fc, g = kG7W[3] * fc;
  for (int j = 0; j < 7; ++j) {
    const double f2 = tvn_integrand(c - hw * kGK15X[j], P) +
                      tvn_integrand(c + hw * kGK15X[j], P);
    k += kGK15W[j] * f2;
    if (j % 2 == 1) g += kG7W[j / 2] * f2;
  }
  err = std::fabs(hw * (k - g));
  return hw * k;
}

// p = Phi_3(b; 0, C) (C column-major 3 x 3). Returns false when the rule
// declines (conditions (a)-(c) in the header); p is not thresholded here.
bool tvn_prob(const double* b, const double* C, double& p) {
  double sd[3], h[3];
  for (int i = 0; i < 3; ++i) {
    if (!std::isfinite(b[i]) || !(C[i + 3 * i] > 0.0) ||
        !std::isfinite(C[i + 3 * i]))
      return false;
    sd[i] = std::sqrt(C[i + 3 * i]);
    h[i] = b[i] / sd[i];
    if (!std::isfinite(h[i])) return false;
  }
  double r[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double v = (i == j) ? 1.0
                                : C[i + 3 * j] / std::sqrt(C[i + 3 * i] * C[j + 3 * j]);
      if (!std::isfinite(v) || (i != j && !(std::fabs(v) < 1.0))) return false;
      r[i][j] = v;
    }
  // relabel: (k2, k3) = the pair with the largest |r|, k1 the remaining one
  int k1 = 0, k2 = 1, k3 = 2;
  if (std::fabs(r[0][2]) > std::fabs(r[k2][k3])) { k1 = 1; k2 = 0; k3 = 2; }
  if (std::fabs(r[0][1]) > std::fabs(r[k2][k3])) { k1 = 2; k2 = 0; k3 = 1; }
  const double r12 = r[k1][k2], r13 = r[k1][k3], r23 = r[k2][k3];
  const double det = 1.0 - r12 * r12 - r13 * r13 - r23 * r23 +
                     2.0 * r12 * r13 * r23;
  if (!(det >= kTvnMinDet)) return false;
  const TvnPath P = {h[k1], h[k2], h[k3], r23, std::asin(r12), std::asin(r13)};

  double lo[kTvnMaxIntervals], hi[kTvnMaxIntervals], val[kTvnMaxIntervals],
         er[kTvnMaxIntervals];
  int n = 1;
  lo[0] = 0.0; hi[0] = 1.0;
  val[0] = tvn_gk15(0.0, 1.0, P, er[0]);
  const double tol = kTvnTol * 2.0 * kPi;
  bool converged = false;
  for (;;) {
    double etot = 0.0;
    int iw = 0;
    for (int i = 0; i < n; ++i) {
      etot += er[i];
      if (er[i] > er[iw]) iw = i;
    }
    if (etot <= tol) { converged = true; break; }
    if (n == kTvnMaxIntervals) break;
    const double mid = 0.5 * (lo[iw] + hi[iw]);
    lo[n] = mid; hi[n] = hi[iw]; hi[iw] = mid;
    val[iw] = tvn_gk15(lo[iw], hi[iw], P, er[iw]);
    val[n] = tvn_gk15(lo[n], hi[n], P, er[n]);
    ++n;
  }
  if (!converged) return false;
  double integral = 0.0;
  for (int i = 0; i < n; ++i) integral += val[i];
  p = Phi1(P.h1) * bvnu(-P.h2, -P.h3, r23) + integral / (2.0 * kPi);
  return std::isfinite(p);
}

// log p under the threshold rule (a)-(d); false = decline (lattice).
bool tvn_logcdf(const double* b, const double* C, double& out) {
  double p;
  if (!tvn_prob(b, C, p)) return false;
  if (!(p >= kTvnMinP) || !(p <= 1.0 + 1e-12)) return false;
  out = std::log(std::min(p, 1.0));
  return true;
}

}  // namespace

// [[Rcpp::export(rng = false)]]
double mvn_logcdf3_cpp(Rcpp::NumericVector b, Rcpp::NumericMatrix C) {
  if (b.size() != 3 || C.nrow() != 3 || C.ncol() != 3) return NA_REAL;
  const std::vector<double> bv(b.begin(), b.end()), Cv(C.begin(), C.end());
  double out;
  if (!tvn_logcdf(bv.data(), Cv.data(), out)) return NA_REAL;
  return out;
}


// ---------------------------------------------------------------------------
// mvn_logcdf_dispatch_cpp(x, S, miwa_qmax): the check = TRUE ("accurate")
// path of logcdf_ME_r() (R/pskf-likelihood.R) moved out of the interpreter,
// operation for operation, so it returns BIT-IDENTICAL values: correlation
// form (sdv = sqrt(max(diag S, eps)), Cm = S / (sdv_i sdv_j), x / sdv), snap
// of |off-diagonal| < 1e-12 to 0, the connected components of Cm != 0
// (summed in the same order, each re-entering the recursion), then pnorm
// (q = 1), the log-scale bivariate rule of .mvn_logcdf2 (q = 2: its swap,
// the |rho| <= 0.9999 cap, mvn_logcdf2_cpp), for q = 3 the exact TVN
// (tvn_logcdf, the body of mvn_logcdf3_cpp, which .mvn_logcdf_sov calls on
// the same inputs) and, when that declines or q = 4..miwa_qmax, the lattice
// evaluator mvn_logcdf_cpp with its default tolerances.
// It returns NA -- and the R code then runs its own path from the top --
// whenever that path would leave these evaluators: a non-finite input, a
// bivariate term below mvn_logcdf2_cpp's p >= 1e-3 range (quadrature), a
// numerically singular lattice call (Miwa / Mendell-Elston), or a block
// larger than miwa_qmax (Mendell-Elston). A3 (2026-09-29): the R-level
// dispatch cost ~30 us per call, ~60 ms of a 200-period Reiter-HANK PSKF
// likelihood (~2200 calls).
namespace {

bool logcdf_dispatch(const std::vector<double>& x, const std::vector<double>& S,
                     int q, double miwa_qmax, double& out) {
  if (q == 0) { out = 0.0; return true; }
  if (q == 1) {
    const double b = x[0] / std::sqrt(S[0]);
    if (!std::isfinite(b)) return false;
    out = R::pnorm(b, 0.0, 1.0, 1, 1);
    return true;
  }
  // correlation form + snap (logcdf_ME_r, check = TRUE)
  std::vector<double> sdv(q), Cm(q * q), xs(q);
  for (int i = 0; i < q; ++i)
    sdv[i] = std::sqrt(std::max(S[i + i * q], DBL_EPSILON));
  bool all_adj = true;
  for (int j = 0; j < q; ++j)
    for (int i = 0; i < q; ++i) {
      double c = S[i + j * q] / (sdv[i] * sdv[j]);
      if (i != j && std::fabs(c) < 1e-12) c = 0.0;
      Cm[i + j * q] = c;
      if (c == 0.0) all_adj = false;
    }
  for (int i = 0; i < q; ++i) {
    xs[i] = x[i] / sdv[i];
    if (!std::isfinite(xs[i])) return false;
  }
  if (!all_adj) {
    // connected components of the (row-wise) adjacency Cm != 0, numbered in
    // order of their smallest member, exactly as the R traversal does
    std::vector<int> comp(q, 0);
    int n_comp = 0;
    for (int s0 = 0; s0 < q; ++s0) {
      if (comp[s0] > 0) continue;
      ++n_comp;
      std::vector<int> stack(1, s0);
      while (!stack.empty()) {
        const int v = stack.front();
        stack.erase(stack.begin());
        if (comp[v] > 0) continue;
        comp[v] = n_comp;
        for (int w = 0; w < q; ++w)
          if (Cm[v + w * q] != 0.0 && comp[w] == 0) stack.push_back(w);
      }
    }
    if (n_comp > 1) {
      double tot = 0.0;
      for (int cc = 1; cc <= n_comp; ++cc) {
        std::vector<int> ii;
        for (int i = 0; i < q; ++i) if (comp[i] == cc) ii.push_back(i);
        const int m = ii.size();
        std::vector<double> xb(m), Sb(m * m);
        for (int a = 0; a < m; ++a) {
          xb[a] = xs[ii[a]];
          for (int b = 0; b < m; ++b) Sb[a + b * m] = Cm[ii[a] + ii[b] * q];
        }
        double v;
        if (!logcdf_dispatch(xb, Sb, m, miwa_qmax, v)) return false;
        tot += v;
      }
      out = tot;
      return true;
    }
  }
  if (q >= 4) {
    // weak couplings inside a block that needs the lattice or
    // Mendell-Elston: zero them and split again (logcdf_ME_r, same rule)
    bool weak = false;
    std::vector<double> C2(Cm);
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < q; ++i)
        if (i != j && Cm[i + j * q] != 0.0 && std::fabs(Cm[i + j * q]) < 1e-8) {
          C2[i + j * q] = 0.0;
          weak = true;
        }
    if (weak) return logcdf_dispatch(xs, C2, q, miwa_qmax, out);
  }
  if (q == 2) {
    // .mvn_logcdf2(x[1], x[2], S[1, 2])
    double h1 = xs[0], h2 = xs[1];
    if (h2 < h1) std::swap(h1, h2);
    const double rho = std::max(-0.9999, std::min(0.9999, Cm[0 + 1 * q]));
    const double v = mvn_logcdf2_cpp(h1, h2, rho);
    if (ISNAN(v)) return false;
    out = v;
    return true;
  }
  if ((double) q > miwa_qmax) return false;
  // .mvn_logcdf_sov(x, S) with every bound finite: q = 3 first tries
  // mvn_logcdf3_cpp (Genz TVN, same threshold rule), then the lattice
  if (q == 3) {
    double v;
    if (tvn_logcdf(xs.data(), Cm.data(), v)) { out = v; return true; }
  }
  // mvn_logcdf_cpp(b, C)[1]
  Rcpp::NumericVector b(xs.begin(), xs.end());
  Rcpp::NumericMatrix C(q, q);
  std::copy(Cm.begin(), Cm.end(), C.begin());
  const Rcpp::NumericVector res = mvn_logcdf_cpp(b, C, 1e-5, 1e-7, -1);
  if (ISNAN(res[0])) return false;
  out = res[0];
  return true;
}

}  // namespace

// [[Rcpp::export(rng = false)]]
double mvn_logcdf_dispatch_cpp(Rcpp::NumericVector x, Rcpp::NumericMatrix S,
                               double miwa_qmax) {
  const int q = x.size();
  if (S.nrow() != q || S.ncol() != q) return NA_REAL;
  for (int i = 0; i < q; ++i) if (!std::isfinite(x[i])) return NA_REAL;
  for (int i = 0; i < q * q; ++i) if (!std::isfinite(S[i])) return NA_REAL;
  std::vector<double> xv(x.begin(), x.end()), Sv(S.begin(), S.end());
  double out;
  if (!logcdf_dispatch(xv, Sv, q, miwa_qmax, out)) return NA_REAL;
  return out;
}
