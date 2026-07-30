// hank_egm2.cpp -- compiled two-asset (liquid/illiquid) EGM household-solver
// kernels.
//
// Ports R/hank-egm2.R::.hank_egm2_step() / hank_egm2_solve()'s iteration loop
// to Rcpp. The R implementation is the REFERENCE SPEC and remains available
// via backend = "R" (or options(dynhr.hank_backend = "R")); these kernels
// follow it line-for-line, INCLUDING the tiny = 1e-12 transient-NaN floors on
// Wa/Wb/c (see the header of R/hank-egm2.R, deviation F3) and the deliberate
// absence of any clamp on a' (deviation F2 -- SSJ clamps only b; the Young
// lottery handles overshoot downstream).
//
// ALL input validation, the k_grid direction check, the feasibility pre-check
// and the default Vb_init/Va_init construction stay in the R wrapper, so error
// messages and rejected inputs are identical across backends and the kernels
// never need their own default-init logic to match bit-for-bit.
//
// The upstream reference for the algorithm itself is SSJ's
// hetblocks/hh_twoasset.py; the transcription and its three deviations are
// pinned in briefs/19-twoasset-hank-scope.md section 3. Do not change a
// formula here without changing R/hank-egm2.R and re-running
// test-hank-egm2-cpp-parity.R.
//
// ARRAY CONVENTION: (e, b, a) with e slowest and ILLIQUID a FASTEST, matching
// R's column-major dim = c(n_e, n_b, n_a). Linear offset of (e, j, k) is
//   e + n_e*j + n_e*n_mid*k
// so the (e,mid)-by-a unfolding used by the crossing search is a free reshape,
// exactly as in the R code.

//
// THREADING (2026-07-29). The step is parallelised over its natural disjoint
// outer indices; see the Egm2Ctx / phase comments below for the dependency
// analysis and why there are only TWO parallel phases (hence three barrier
// arrivals) per backward step. Because every write is disjoint and there is no
// reduction anywhere, output is bit-identical at every thread count by
// construction, not by calibration -- test-hank-egm2-threads.R pins that with
// identical(), and the SERIAL path (threads <= 1) is the same code with a
// single worker, so it cannot drift from the reference either.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <atomic>
#include <barrier>
#include <memory>
#include <thread>
#include <vector>
using namespace Rcpp;

// Every power here goes through R's own R_pow (R_ext/Arith.h), NOT std::pow,
// so that `^` can never be a source of R-vs-cpp drift.
//
// Measured over 2e5 random x in (0.001, 500), comparing against R's `^`:
//   exponent    1   1.5    2    2.5    3    0.5   -0.5   -2   1/3
//   std::pow    =    =    278d   =   1093d   =     =     =     =
//   R_pow       =    =     =     =     =     =     =     =     =
// i.e. std::pow disagrees by 1 ULP at INTEGER exponents >= 2 (R evaluates
// those by repeated multiplication rather than a libm call), and agrees
// elsewhere -- including the non-integer 1.5 / 2.5 that chi2 = 2.5 would use.
// R_pow matches R's `^` at every exponent tested, so it is bit-identical by
// construction rather than by luck of calibration, and it is faster where the
// special case applies.
//
// This matters here because core's exponent is chi2 - 1: chi2 = 3 lands
// exactly on the integer-exponent case. (Historical note: this was suspected
// of causing a 1e-12 gap at chi2 = 3, and switching to R_pow did NOT fix it --
// the real cause was an infeasible calibration pinning c to the floor. Both
// facts are worth keeping: the gap had a different cause, and routing through
// R_pow is still correct.)
static inline double rpow(double x, double y) { return R_pow(x, y); }

// Adjustment cost and derivatives -- port of R .hank_psi(). Scalar form.
static inline void psi_and_deriv(double ap, double a, double ra, double chi0,
                                 double chi1, double chi2,
                                 double *Psi, double *Psi1, double *Psi2) {
  double a_with_return = (1.0 + ra) * a;
  double a_change      = ap - a_with_return;
  double abs_change    = std::fabs(a_change);
  double sign_change   = (a_change > 0.0) - (a_change < 0.0);  // R sign()
  double adj_denom     = a_with_return + chi0;
  double core          = rpow(abs_change / adj_denom, chi2 - 1.0);
  double P  = (chi1 / chi2) * abs_change * core;
  double P1 = chi1 * sign_change * core;
  if (Psi)  *Psi  = P;
  if (Psi1) *Psi1 = P1;
  if (Psi2) *Psi2 = -(1.0 + ra) * (P1 + (chi2 - 1.0) * P / adj_denom);
}

// Interpolation coordinates by SSJ's FORWARD SWEEP over increasing query
// points -- port of R .hank2_interp_coord_rows(). Writes 0-indexed lower
// bracket and lower weight for each query; extrapolates linearly out of domain
// (the weight is not clamped, only the index).
//
// Deliberately NOT a binary search. On strictly increasing knots the two agree
// (up to which side of an exact tie is bracketed, which interpolates to the
// same value), but the endogenous liquid grid b_endo is non-monotone in
// TRANSIENT iterates -- notably the very first step from SSJ's hh_init at
// plausible calibrations -- and the sweep degrades gracefully there exactly as
// the reference does, while a binary search would return a silently wrong
// bracket and R's findInterval would error. Keeping the sweep in both backends
// is what makes them agree on non-monotone input as well as monotone.
// See the long comment on .hank2_interp_coord_rows() in R/hank-egm2.R.
static inline void interp_coord_sweep(const double *x, int n_x,
                                      const double *xq, int n_q,
                                      int *iout, double *pout) {
  int xi = 0;
  double x_low = x[0], x_high = x[1];
  for (int q = 0; q < n_q; ++q) {
    double xqc = xq[q];
    while (xi < n_x - 2) {
      if (x_high >= xqc) break;
      ++xi; x_low = x_high; x_high = x[xi + 1];
    }
    iout[q] = xi;
    pout[q] = (x_high - xqc) / (x_high - x_low);
  }
}

// Crossing search -- port of R .hank_lhs_eq_rhs() for ONE lhs row.
// lhs: length n_i (decreasing net of rhs); rhs: n_i x n_j, column-major.
// Writes 0-indexed lower index and lower weight for each j.
// Matches the R per-column "first i with lhs - rhs < 0" form (brief 19, F5),
// including both corners: bottom -> (0, 1) exactly a_grid[0]; no crossing by
// the top -> (n_i-2, p) i.e. linear extrapolation above the grid.
static inline void lhs_eq_rhs_row(const double *lhs, const double *rhs,
                                  int n_i, int n_j, int *iout, double *pout) {
  for (int j = 0; j < n_j; ++j) {
    const double *rc = rhs + (std::size_t)j * n_i;   // column j
    int m = 0;
    for (int i = 0; i < n_i; ++i) if (lhs[i] - rc[i] >= 0.0) ++m;
    int i1 = m;                       // 0-indexed "first i with D < 0"
    if (i1 > n_i - 1) i1 = n_i - 1;   // cap at the top
    if (i1 == 0) { iout[j] = 0; pout[j] = 1.0; continue; }
    int lo = i1 - 1;
    double Dl = lhs[lo] - rc[lo];
    double Du = lhs[i1] - rc[i1];
    iout[j] = lo;
    pout[j] = -Du / (Dl - Du);
  }
}

// ---------------------------------------------------------------------------
// Shared step context.
//
// Everything the backward step touches lives here so that the SOLVE loop can
// allocate it ONCE and keep a worker pool alive across iterations (design A of
// the threading brief). Two-asset per-step work is ~1 ms, so a spawn/join per
// step -- let alone per stage -- would be a large fraction of the runtime; the
// pool is created once per solve and the workers are parked on a std::barrier
// between phases and between iterations.
//
// The buffers are plain std::vector, never Rcpp vectors: a worker thread must
// never touch the R API (allocation, protection, error signalling). Results are
// copied into NumericVectors on the MAIN thread after the last join.
// ---------------------------------------------------------------------------
struct Egm2Ctx {
  int n_e, n_b, n_a, n_k, n_lead_b, n_col;
  std::size_t n_cell, n_ck;
  const double *b_grid, *a_grid, *k_grid, *y;
  double rb, ra, beta, eis, chi0, chi1, chi2;
  double tiny;
  // Step 2 outputs / step 3 inputs
  std::vector<double> rhs, Wb, Wa, W_ratio;
  // Step 3 / 4 workspace
  std::vector<double> a_endo_unc, c_endo_unc, b_endo, a_unc, b_unc;
  // Steps 5-6 workspace
  std::vector<double> a_endo_con, c_endo_con, b_endo_k, a_con;
  // Step 7 outputs
  std::vector<double> Vb_o, Va_o, b_o, a_o, c_o, chi_o;

  Egm2Ctx(int ne, int nb, int na, int nk,
          const double *bg, const double *ag, const double *kg,
          const double *yy, const double *Psi1g,
          double rb_, double ra_, double beta_, double eis_,
          double chi0_, double chi1_, double chi2_)
    : n_e(ne), n_b(nb), n_a(na), n_k(nk), n_lead_b(ne * nb), n_col(nb * na),
      n_cell((std::size_t)ne * nb * na), n_ck((std::size_t)ne * nk * na),
      b_grid(bg), a_grid(ag), k_grid(kg), y(yy),
      rb(rb_), ra(ra_), beta(beta_), eis(eis_),
      chi0(chi0_), chi1(chi1_), chi2(chi2_), tiny(1e-12),
      rhs((std::size_t)na * na), Wb(n_cell), Wa(n_cell), W_ratio(n_cell),
      a_endo_unc(n_cell), c_endo_unc(n_cell), b_endo(n_cell),
      a_unc(n_cell), b_unc(n_cell),
      a_endo_con(n_ck), c_endo_con(n_ck), b_endo_k(n_ck), a_con(n_cell),
      Vb_o(n_cell), Va_o(n_cell), b_o(n_cell), a_o(n_cell),
      c_o(n_cell), chi_o(n_cell) {
    // rhs = 1 + Psi1_grid (n_a x n_a, column-major)
    for (std::size_t t = 0; t < rhs.size(); ++t) rhs[t] = 1.0 + Psi1g[t];
  }
};

// Even split of [0, ntask) over nthr workers, worker w's half-open range.
static inline void egm2_range(std::size_t ntask, int w, int nthr,
                              std::size_t *lo, std::size_t *hi) {
  const std::size_t chunk = (ntask + (std::size_t)nthr - 1) / (std::size_t)nthr;
  *lo = std::min(ntask, (std::size_t)w * chunk);
  *hi = std::min(ntask, *lo + chunk);
}

// --- Step 2: Wb, Wa = beta * Pi %*% V_p, floored at tiny --------------------
// Serial, on the main thread. Pi contracts over TODAY's e; e is the leading
// axis, so the unfolding to n_e x (n_b*n_a) is a free reinterpretation of the
// same buffer. This is a 3 x 3 by 3 x (n_b*n_a) GEMM plus one pass over
// n_cell -- a few percent of the step, and it is a REDUCTION over e, which is
// exactly the thing that could not be split without risking a reordered sum.
// Keeping it serial is what makes bit-identity across thread counts structural.
static void egm2_step2_run(Egm2Ctx &X, double *Vb_p, double *Va_p,
                           double *Pi_ptr) {
  arma::mat Pi(Pi_ptr, X.n_e, X.n_e, false);
  arma::mat Vb_p_m(Vb_p, X.n_e, X.n_col, false);
  arma::mat Va_p_m(Va_p, X.n_e, X.n_col, false);
  arma::mat Wb_m = X.beta * (Pi * Vb_p_m);
  arma::mat Wa_m = X.beta * (Pi * Va_p_m);
  const double tiny = X.tiny;
  for (std::size_t t = 0; t < X.n_cell; ++t) {
    double wb = Wb_m[t] < tiny ? tiny : Wb_m[t];
    double wa = Wa_m[t] < tiny ? tiny : Wa_m[t];
    X.Wb[t] = wb; X.Wa[t] = wa; X.W_ratio[t] = wa / wb;
  }
}

// --- Phase 1 ---------------------------------------------------------------
// Step 3 (unconstrained illiquid FOC, over l = (e, b'), n_e*n_b tasks) and
// Steps 5 (constrained branch's endogenous grids, over (kappa, e), n_k*n_e
// tasks). They are INDEPENDENT of each other -- both read only step 2's
// Wb / W_ratio and the fixed rhs -- so they share one parallel region with no
// barrier between them, and each worker simply does its slice of both.
//
// Writes are disjoint: step 3 writes a_endo_unc / c_endo_unc at offsets
// l + n_lead_b*j (l fixed per task), step 5 writes a_endo_con / c_endo_con /
// b_endo_k at e + n_e*kk + n_e*n_k*j ((kk, e) fixed per task). The crossing
// search's lhs / index / weight buffers are PER-WORKER scratch -- sharing them
// is the classic race in this kernel.
static void egm2_phase1(Egm2Ctx &X, int w, int nthr) {
  const int n_a = X.n_a, n_b = X.n_b, n_e = X.n_e, n_k = X.n_k;
  const int n_lead_b = X.n_lead_b;
  const double tiny = X.tiny, eis = X.eis;
  std::vector<double> lhs_buf(n_a), pw(n_a);
  std::vector<int> iw(n_a);
  std::size_t lo, hi;

  // Step 3
  egm2_range((std::size_t)n_lead_b, w, nthr, &lo, &hi);
  for (std::size_t lt = lo; lt < hi; ++lt) {
    const int l = (int)lt;
    // row l of the (n_e*n_b) x n_a unfolding: element (l, i) at l + n_lead_b*i
    for (int i = 0; i < n_a; ++i)
      lhs_buf[i] = X.W_ratio[l + (std::size_t)n_lead_b * i];
    lhs_eq_rhs_row(lhs_buf.data(), X.rhs.data(), n_a, n_a, iw.data(), pw.data());
    for (int j = 0; j < n_a; ++j) {
      std::size_t o = l + (std::size_t)n_lead_b * j;
      int i0 = iw[j]; double p = pw[j];
      X.a_endo_unc[o] = p * X.a_grid[i0] + (1.0 - p) * X.a_grid[i0 + 1];
      // Floor the INTERPOLATED Wb before the power (F3): at the top corner the
      // crossing search extrapolates (weight outside [0,1]) and, since Wb
      // falls in a', a far extrapolation can drive it negative -- whereupon
      // pow(neg, -eis) is NaN. No-op wherever the extrapolation stays
      // positive. Mirrors R .hank_egm2_step step 3.
      double wbi = p * X.Wb[l + (std::size_t)n_lead_b * i0] +
                   (1.0 - p) * X.Wb[l + (std::size_t)n_lead_b * (i0 + 1)];
      if (wbi < tiny) wbi = tiny;
      X.c_endo_unc[o] = rpow(wbi, -eis);
    }
  }

  // --- Steps 5: liquid-constrained branch (b' = b_grid[0]) -----------------
  // For each (e, kappa, a): W_ratio(e, 0, .)/(1+kappa) == 1 + Psi1.
  egm2_range((std::size_t)n_e * n_k, w, nthr, &lo, &hi);
  for (std::size_t t = lo; t < hi; ++t) {
    const int kk = (int)(t / (std::size_t)n_e);
    const int e  = (int)(t % (std::size_t)n_e);
    double kfac = 1.0 + X.k_grid[kk];
    // lhs = W_ratio(e, b'=0, .)/(1+kappa); Wb slice at b'=0 is stride n_e
    for (int i = 0; i < n_a; ++i)
      lhs_buf[i] = X.W_ratio[e + (std::size_t)n_e * n_b * i] / kfac;
    lhs_eq_rhs_row(lhs_buf.data(), X.rhs.data(), n_a, n_a, iw.data(), pw.data());
    for (int j = 0; j < n_a; ++j) {
      std::size_t o = e + (std::size_t)n_e * kk + (std::size_t)n_e * n_k * j;
      int i0 = iw[j]; double p = pw[j];
      double ae = p * X.a_grid[i0] + (1.0 - p) * X.a_grid[i0 + 1];
      double wbi = p * X.Wb[e + (std::size_t)n_e * n_b * i0] +
                   (1.0 - p) * X.Wb[e + (std::size_t)n_e * n_b * (i0 + 1)];
      if (wbi < tiny) wbi = tiny;            // floored base, as in step 3
      double ce = rpow(kfac, -eis) * rpow(wbi, -eis);
      X.a_endo_con[o] = ae;
      X.c_endo_con[o] = ce;
      double P; psi_and_deriv(ae, X.a_grid[j], X.ra, X.chi0, X.chi1, X.chi2,
                              &P, nullptr, nullptr);
      X.b_endo_k[o] = (ce + ae + X.b_grid[0] - X.y[e]
                       - (1.0 + X.ra) * X.a_grid[j] + P) / (1.0 + X.rb);
    }
  }
}

// --- Phase 2 ---------------------------------------------------------------
// Step 4 (budget inversion + the b_endo -> b' remap), step 6 (the kappa -> b
// remap of the constrained branch) and step 7 (combine + budget residual), ALL
// indexed by the illiquid state k and all writing only at offsets carrying that
// same k. Chunking them identically over k means the worker that produced
// a_unc / b_unc / a_con for a given k is the one that consumes them, so no
// barrier is needed between the three -- they collapse into one parallel phase.
//
// Step 4's b_endo for a given k is likewise read only at that k, so its two
// halves fuse as well.
static void egm2_phase2(Egm2Ctx &X, int w, int nthr) {
  const int n_a = X.n_a, n_b = X.n_b, n_e = X.n_e, n_k = X.n_k;
  const double tiny = X.tiny, eis = X.eis, ra = X.ra, rb = X.rb;
  std::vector<double> knots(n_b), vals(n_b), pco(n_b), knk(n_k), vak(n_k);
  std::vector<int> ico(n_b);
  std::size_t lo, hi;
  egm2_range((std::size_t)n_a, w, nthr, &lo, &hi);

  for (std::size_t kt = lo; kt < hi; ++kt) {
    const int k = (int)kt;

    // --- Step 4: invert the budget -> b'(e, b, a), a'(e, b, a) -------------
    // b_endo(e, b', a) = (c + a' + b' - y - (1+ra)a + Psi(a', a)) / (1+rb)
    for (int j = 0; j < n_b; ++j)
      for (int e = 0; e < n_e; ++e) {
        std::size_t o = e + (std::size_t)n_e * j + (std::size_t)n_e * n_b * k;
        double P; psi_and_deriv(X.a_endo_unc[o], X.a_grid[k], ra,
                                X.chi0, X.chi1, X.chi2, &P, nullptr, nullptr);
        X.b_endo[o] = (X.c_endo_unc[o] + X.a_endo_unc[o] + X.b_grid[j] - X.y[e]
                       - (1.0 + ra) * X.a_grid[k] + P) / (1.0 + rb);
      }
    // Interpolate the b_endo -> b' map at the fixed b_grid, per (e, a).
    for (int e = 0; e < n_e; ++e) {
      for (int j = 0; j < n_b; ++j) {
        std::size_t o = e + (std::size_t)n_e * j + (std::size_t)n_e * n_b * k;
        knots[j] = X.b_endo[o]; vals[j] = X.a_endo_unc[o];
      }
      interp_coord_sweep(knots.data(), n_b, X.b_grid, n_b,
                         ico.data(), pco.data());
      for (int j = 0; j < n_b; ++j) {
        int i0 = ico[j]; double p = pco[j];
        std::size_t o = e + (std::size_t)n_e * j + (std::size_t)n_e * n_b * k;
        X.a_unc[o] = p * vals[i0] + (1.0 - p) * vals[i0 + 1];
        X.b_unc[o] = p * X.b_grid[i0] + (1.0 - p) * X.b_grid[i0 + 1];
      }
    }

    // --- Step 6: map kappa -> b. b_endo_k is DECREASING in kappa and k_grid
    // is decreasing, so b_endo_k is INCREASING along the kappa axis, as
    // required here (F1).
    for (int e = 0; e < n_e; ++e) {
      for (int kk = 0; kk < n_k; ++kk) {
        std::size_t o = e + (std::size_t)n_e * kk + (std::size_t)n_e * n_k * k;
        knk[kk] = X.b_endo_k[o]; vak[kk] = X.a_endo_con[o];
      }
      interp_coord_sweep(knk.data(), n_k, X.b_grid, n_b,
                         ico.data(), pco.data());
      for (int j = 0; j < n_b; ++j) {
        int i0 = ico[j]; double p = pco[j];
        X.a_con[e + (std::size_t)n_e * j + (std::size_t)n_e * n_b * k] =
          p * vak[i0] + (1.0 - p) * vak[i0 + 1];
      }
    }

    // --- Step 7: combine, then c from the budget residual ------------------
    for (int j = 0; j < n_b; ++j)
      for (int e = 0; e < n_e; ++e) {
        std::size_t o = e + (std::size_t)n_e * j + (std::size_t)n_e * n_b * k;
        double bp = X.b_unc[o], ap = X.a_unc[o];
        if (bp <= X.b_grid[0]) { bp = X.b_grid[0]; ap = X.a_con[o]; }
        double P, P2;
        psi_and_deriv(ap, X.a_grid[k], ra, X.chi0, X.chi1, X.chi2,
                      &P, nullptr, &P2);
        double c = X.y[e] + (1.0 + rb) * X.b_grid[j] + (1.0 + ra) * X.a_grid[k]
                   - P - ap - bp;
        if (c < tiny) c = tiny;
        double uc = rpow(c, -1.0 / eis);
        X.Vb_o[o]  = (1.0 + rb) * uc;
        X.Va_o[o]  = (1.0 + ra - P2) * uc;
        X.b_o[o]   = bp;
        X.a_o[o]   = ap;
        X.c_o[o]   = c;
        X.chi_o[o] = P;
      }
  }
}

// ---------------------------------------------------------------------------
// Persistent worker pool: spawned ONCE (per solve, or per exported step call)
// and parked on a barrier between phases and between backward iterations.
//
// Two-asset per-step work is ~1 ms at realistic grid sizes, so spawning threads
// per step (~7 spawns) -- let alone per stage (~28) -- would spend a double-
// digit percentage of the runtime in thread creation. Holding the pool across
// the whole backward iteration amortises the spawn cost to nothing.
//
// The main thread is worker 0 and does the serial step-2 section while the
// others wait at the top barrier. Sequence per step, from the main thread's
// point of view:
//
//   [serial step 2] -> arrive (releases workers) -> phase 1 -> arrive ->
//   phase 2 -> arrive -> [serial convergence check] -> loop
//
// so three barrier arrivals per step and zero thread creations after the first.
// Workers never touch the R API, and never signal errors: this kernel has no
// error condition inside the step at all (the only failure mode, consumption
// pinned at the numerical floor, is diagnosed by the R wrapper AFTER the solve
// returns), so there is nothing to marshal back to the main thread.
class Egm2Pool {
 public:
  Egm2Pool(Egm2Ctx *ctx, int nthr)
    : ctx_(ctx), nthr_(nthr), bar_((std::ptrdiff_t)nthr), done_(false) {
    pool_.reserve(nthr - 1);
    for (int w = 1; w < nthr; ++w)
      pool_.emplace_back([this, w] { this->worker(w); });
  }
  // Run one backward step's parallel phases. Step 2 must already be done.
  void run_step() {
    bar_.arrive_and_wait();          // release the workers into phase 1
    egm2_phase1(*ctx_, 0, nthr_);
    bar_.arrive_and_wait();
    egm2_phase2(*ctx_, 0, nthr_);
    bar_.arrive_and_wait();          // workers park at the top barrier again
  }
  // Releases and joins the workers. Runs on the normal path AND during stack
  // unwinding if the main thread throws between steps, so a worker can never
  // outlive the buffers it points at.
  ~Egm2Pool() {
    done_.store(true, std::memory_order_release);
    bar_.arrive_and_wait();
    for (auto &t : pool_) t.join();
  }
  Egm2Pool(const Egm2Pool &) = delete;
  Egm2Pool &operator=(const Egm2Pool &) = delete;

 private:
  void worker(int w) {
    for (;;) {
      bar_.arrive_and_wait();
      if (done_.load(std::memory_order_acquire)) return;
      egm2_phase1(*ctx_, w, nthr_);
      bar_.arrive_and_wait();
      egm2_phase2(*ctx_, w, nthr_);
      bar_.arrive_and_wait();
    }
  }
  Egm2Ctx *ctx_;
  int nthr_;
  std::barrier<> bar_;
  std::atomic<bool> done_;
  std::vector<std::thread> pool_;
};

// Number of workers actually usable: never more than the smallest task axis
// that has to be split, and never less than 1.
static inline int egm2_nthreads(const Egm2Ctx &X, int threads) {
  const int max_useful = std::max(X.n_lead_b, std::max(X.n_e * X.n_k, X.n_a));
  int n = threads < 1 ? 1 : threads;
  if (n > max_useful) n = max_useful;
  return n < 1 ? 1 : n;
}

// Copy one context buffer out to a dim'd NumericVector (main thread only).
static inline NumericVector egm2_out(const std::vector<double> &v,
                                     const IntegerVector &d) {
  NumericVector o(v.begin(), v.end());
  o.attr("dim") = d;
  return o;
}

// One two-asset EGM backward step (compiled). C++ port of .hank_egm2_step.
// Vb_p_, Va_p_ are length n_e*n_b*n_a vectors in (e, b, a) column-major order.
// Psi1_grid_ is the n_a x n_a matrix Psi1(a'=a_grid[i], a=a_grid[j]).
// Returns a list with Vb, Va, b, a, c, chi (each n_e*n_b*n_a, dim'd in R).
//
// `threads` is resolved in R (see hank_resolve_threads); <= 1 runs the serial
// path, which is the code the parallel path must reproduce bit for bit. A bare
// step call pays ONE pool spawn; the fused solve below pays one per SOLVE.
// [[Rcpp::export]]
List hank_egm2_step_cpp(NumericVector Vb_p_, NumericVector Va_p_,
                        NumericVector b_grid_, NumericVector a_grid_,
                        NumericVector k_grid_, NumericVector y_,
                        double rb, double ra, double beta, double eis,
                        double chi0, double chi1, double chi2,
                        NumericMatrix Pi_, NumericMatrix Psi1_grid_,
                        int threads = 1) {
  Egm2Ctx X(Pi_.nrow(), b_grid_.size(), a_grid_.size(), k_grid_.size(),
            b_grid_.begin(), a_grid_.begin(), k_grid_.begin(), y_.begin(),
            Psi1_grid_.begin(), rb, ra, beta, eis, chi0, chi1, chi2);
  egm2_step2_run(X, Vb_p_.begin(), Va_p_.begin(), Pi_.begin());
  const int nthr = egm2_nthreads(X, threads);
  if (nthr <= 1) {
    egm2_phase1(X, 0, 1);
    egm2_phase2(X, 0, 1);
  } else {
    Egm2Pool pool(&X, nthr);
    pool.run_step();
  }
  IntegerVector d = IntegerVector::create(X.n_e, X.n_b, X.n_a);
  return List::create(_["Vb"] = egm2_out(X.Vb_o, d),
                      _["Va"] = egm2_out(X.Va_o, d),
                      _["b"]  = egm2_out(X.b_o, d),
                      _["a"]  = egm2_out(X.a_o, d),
                      _["c"]  = egm2_out(X.c_o, d),
                      _["chi"] = egm2_out(X.chi_o, d));
}


// Fused backward-iteration loop (compiled). C++ port of hank_egm2_solve()'s
// loop: iterates the compiled step from the caller-supplied initial marginal
// values until max(|db'|, |da'|) < tol. Vb_init/Va_init come from the R
// wrapper so the two backends start from bit-identical guesses.
//
// `threads` spawns the worker pool ONCE here and holds it for the whole
// backward iteration (design A), so per-step spawn overhead is zero.
// [[Rcpp::export]]
List hank_egm2_solve_cpp(NumericVector Vb_init, NumericVector Va_init,
                         NumericVector b_grid_, NumericVector a_grid_,
                         NumericVector k_grid_, NumericVector y_,
                         double rb, double ra, double beta, double eis,
                         double chi0, double chi1, double chi2,
                         NumericMatrix Pi_, NumericMatrix Psi1_grid_,
                         double tol, int maxit, int threads = 1) {
  Egm2Ctx X(Pi_.nrow(), b_grid_.size(), a_grid_.size(), k_grid_.size(),
            b_grid_.begin(), a_grid_.begin(), k_grid_.begin(), y_.begin(),
            Psi1_grid_.begin(), rb, ra, beta, eis, chi0, chi1, chi2);
  const int nthr = egm2_nthreads(X, threads);
  // Constructed only when it will actually be used, so the serial path spawns
  // nothing and cannot be perturbed by the threading machinery at all.
  std::unique_ptr<Egm2Pool> pool;
  if (nthr > 1) pool.reset(new Egm2Pool(&X, nthr));

  std::vector<double> b_old(X.n_cell), a_old(X.n_cell);
  bool converged = false, have_old = false;
  double *vb_src = Vb_init.begin(), *va_src = Va_init.begin();
  int it = 0;
  for (it = 1; it <= maxit; ++it) {
    egm2_step2_run(X, vb_src, va_src, Pi_.begin());
    if (pool) pool->run_step(); else { egm2_phase1(X, 0, 1); egm2_phase2(X, 0, 1); }
    // The step's Vb/Va output IS the next step's input; no copy needed.
    vb_src = X.Vb_o.data(); va_src = X.Va_o.data();
    if (have_old) {
      double gap = 0.0;
      bool bad = false;
      for (std::size_t t = 0; t < X.n_cell; ++t) {
        double db = std::fabs(X.b_o[t] - b_old[t]);
        double da = std::fabs(X.a_o[t] - a_old[t]);
        if (ISNAN(db) || ISNAN(da)) { bad = true; break; }
        if (db > gap) gap = db;
        if (da > gap) gap = da;
      }
      // NaN-aware: a non-finite iterate must never report convergence (the
      // R path stops with an error; here we hand the flag back so the R
      // wrapper raises the identical message).
      if (bad) { converged = false; break; }
      if (gap < tol) { converged = true; break; }
    }
    b_old.assign(X.b_o.begin(), X.b_o.end());
    a_old.assign(X.a_o.begin(), X.a_o.end());
    have_old = true;
  }
  if (it > maxit) it = maxit;
  pool.reset();   // join before anything else touches X

  IntegerVector d = IntegerVector::create(X.n_e, X.n_b, X.n_a);
  return List::create(_["Vb"] = egm2_out(X.Vb_o, d),
                      _["Va"] = egm2_out(X.Va_o, d),
                      _["b"]  = egm2_out(X.b_o, d),
                      _["a"]  = egm2_out(X.a_o, d),
                      _["c"]  = egm2_out(X.c_o, d),
                      _["chi"] = egm2_out(X.chi_o, d),
                      _["iterations"] = it, _["converged"] = converged);
}


// ---------------------------------------------------------------------------
// FUSED BACKWARD SWEEP for the two-asset fake-news Jacobian.
//
// Ports the s = 2 .. T_h anticipation loop of R .hank_curly_sweep2() (see
// R/hank-jacobian2.R) so that the WHOLE sweep runs under ONE worker pool.
//
// WHY. Since 0.9.0.0034 each .hank_block_step2() reaches hank_egm2_step_cpp,
// and a bare step call spawns and joins its own pool. At T_h = 50 that is 98
// spawn/join pairs around a sub-millisecond kernel, so threading bought
// nothing (t = 8 measured a hair SLOWER than t = 1) even though pure
// compilation bought 3.2-3.5x. hank_egm2_solve_cpp never had this problem: it
// spawns one pool and parks the workers on a std::barrier across all backward
// iterations. This routine does exactly the same thing for the sweep, reusing
// Egm2Ctx / Egm2Pool / egm2_phase1 / egm2_phase2 unchanged rather than
// introducing a second threading scheme.
//
// SCOPE. COLLATERAL IS NOT SUPPORTED, because the compiled step is not: there
// is no theta_coll / dtheta_next anywhere in this file. The R caller must keep
// routing theta_coll != 0 blocks AND the "theta_coll" input (whose s = 2 step
// carries a non-zero dtheta_next) to the R path -- silently dropping
// collateral would be a wrong answer, not a slow one.
//
// The s = 1 term, the curly-D pushes and the curlyY aggregation stay in R:
// they are a different algorithm (the forward lottery), they are cheap
// relative to the 2*(T_h-1) backward steps, and keeping them in R keeps the
// C++/R boundary where the risk is lowest.
// ---------------------------------------------------------------------------

// Bit-identical port of R's `V_ss + h * dV`.
//
// The multiply and the add must be rounded SEPARATELY, exactly as R rounds
// them. A fused multiply-add (clang defaults to -ffp-contract=on, which is
// licensed to contract this very expression) would move the perturbed argument
// by 1 ULP; the sweep's 1/(2h) ~ 5e5 amplification turns that into a ~1e-10
// relative wobble in dB, which then COMPOUNDS date after date because dVb/dVa
// feed the next step. `volatile` on the intermediate product pins the rounding
// regardless of the compiler's contraction default. It costs one store per
// cell per step and is invisible next to the EGM step itself.
//
// Passing -h reproduces R's `V_ss - h * dV` exactly too: negation is exact, so
// v + (-(h*d)) and v - (h*d) round identically.
static inline void egm2_axpy_exact(const double *v, double h, const double *d,
                                   double *out, std::size_t n) {
  for (std::size_t t = 0; t < n; ++t) {
    volatile double p = h * d[t];
    out[t] = v[t] + p;
  }
}

// Port of R's `max(1, max(abs(dVb)), max(abs(dVa)))`, NaN-propagating like R's
// max(). max() is exact (pure comparison), so flattening R's three-way max
// into one running maximum cannot change the result by even a ULP -- and it
// MUST not, because h is the single quantity every later date inherits.
static inline double egm2_sweep_scale(const double *dvb, const double *dva,
                                      std::size_t n) {
  double m = 1.0;
  bool nan_seen = false;
  for (std::size_t t = 0; t < n; ++t) {
    double u = std::fabs(dvb[t]);
    if (ISNAN(u)) { nan_seen = true; continue; }
    if (u > m) m = u;
  }
  for (std::size_t t = 0; t < n; ++t) {
    double u = std::fabs(dva[t]);
    if (ISNAN(u)) { nan_seen = true; continue; }
    if (u > m) m = u;
  }
  return nan_seen ? R_NaN : m;
}

// Fused s = 2 .. T_h backward sweep (compiled).
//
// Vb_ss_, Va_ss_ are the steady-state marginals; dVb0_, dVa0_ are the s = 1
// directional derivatives the R caller has already computed. Returns the
// per-date policy derivatives dB, dA, dC, dCHI as n_cell x (T_h - 1) matrices,
// column j holding date s = j + 2. The R side still builds curly-D from
// dB/dA and aggregates curlyY -- both cheap, and both a different algorithm.
//
// `threads` spawns ONE pool for the entire sweep; the two FD evaluations of a
// date and every date share it, so the spawn count is 1 rather than
// 2*(T_h-1). Output is independent of `threads` for the same reason the step
// is (disjoint writes, no reduction in either phase).
// [[Rcpp::export]]
List hank_curly_sweep2_cpp(NumericVector Vb_ss_, NumericVector Va_ss_,
                           NumericVector dVb0_, NumericVector dVa0_,
                           NumericVector b_grid_, NumericVector a_grid_,
                           NumericVector k_grid_, NumericVector y_,
                           double rb, double ra, double beta, double eis,
                           double chi0, double chi1, double chi2,
                           NumericMatrix Pi_, NumericMatrix Psi1_grid_,
                           double delta_va, int T_h, int threads = 1) {
  Egm2Ctx X(Pi_.nrow(), b_grid_.size(), a_grid_.size(), k_grid_.size(),
            b_grid_.begin(), a_grid_.begin(), k_grid_.begin(), y_.begin(),
            Psi1_grid_.begin(), rb, ra, beta, eis, chi0, chi1, chi2);
  const std::size_t n = X.n_cell;
  const int S = T_h > 1 ? T_h - 1 : 0;     // dates s = 2 .. T_h

  NumericMatrix dB((int)n, S), dA((int)n, S), dC((int)n, S), dCHI((int)n, S);
  std::vector<double> dVb_prev(dVb0_.begin(), dVb0_.end());
  std::vector<double> dVa_prev(dVa0_.begin(), dVa0_.end());
  std::vector<double> vbp(n), vap(n);
  // +h outputs, held while the -h step overwrites the context buffers.
  std::vector<double> Vb_p(n), Va_p(n), b_p(n), a_p(n), c_p(n), chi_p(n);

  const int nthr = egm2_nthreads(X, threads);
  std::unique_ptr<Egm2Pool> pool;
  if (nthr > 1 && S > 0) pool.reset(new Egm2Pool(&X, nthr));

  const double *Vb_ss = Vb_ss_.begin(), *Va_ss = Va_ss_.begin();
  for (int s = 0; s < S; ++s) {
    const double h = delta_va / egm2_sweep_scale(dVb_prev.data(),
                                                 dVa_prev.data(), n);
    // --- +h evaluation -----------------------------------------------------
    egm2_axpy_exact(Vb_ss, h, dVb_prev.data(), vbp.data(), n);
    egm2_axpy_exact(Va_ss, h, dVa_prev.data(), vap.data(), n);
    egm2_step2_run(X, vbp.data(), vap.data(), Pi_.begin());
    if (pool) pool->run_step(); else { egm2_phase1(X, 0, 1); egm2_phase2(X, 0, 1); }
    std::copy(X.Vb_o.begin(), X.Vb_o.end(), Vb_p.begin());
    std::copy(X.Va_o.begin(), X.Va_o.end(), Va_p.begin());
    std::copy(X.b_o.begin(),  X.b_o.end(),  b_p.begin());
    std::copy(X.a_o.begin(),  X.a_o.end(),  a_p.begin());
    std::copy(X.c_o.begin(),  X.c_o.end(),  c_p.begin());
    std::copy(X.chi_o.begin(), X.chi_o.end(), chi_p.begin());

    // --- -h evaluation (leaves its outputs in the context) -----------------
    egm2_axpy_exact(Vb_ss, -h, dVb_prev.data(), vbp.data(), n);
    egm2_axpy_exact(Va_ss, -h, dVa_prev.data(), vap.data(), n);
    egm2_step2_run(X, vbp.data(), vap.data(), Pi_.begin());
    if (pool) pool->run_step(); else { egm2_phase1(X, 0, 1); egm2_phase2(X, 0, 1); }

    // --- central difference, denominator formed ONCE as R does -------------
    const double den = 2.0 * h;
    double *cB = &dB(0, s), *cA = &dA(0, s);
    double *cC = &dC(0, s), *cCHI = &dCHI(0, s);
    for (std::size_t t = 0; t < n; ++t) {
      cB[t]   = (b_p[t]   - X.b_o[t])   / den;
      cA[t]   = (a_p[t]   - X.a_o[t])   / den;
      cC[t]   = (c_p[t]   - X.c_o[t])   / den;
      cCHI[t] = (chi_p[t] - X.chi_o[t]) / den;
    }
    for (std::size_t t = 0; t < n; ++t) {
      dVb_prev[t] = (Vb_p[t] - X.Vb_o[t]) / den;
      dVa_prev[t] = (Va_p[t] - X.Va_o[t]) / den;
    }
  }
  pool.reset();   // join before anything else touches X

  return List::create(_["dB"] = dB, _["dA"] = dA, _["dC"] = dC,
                      _["dCHI"] = dCHI);
}
