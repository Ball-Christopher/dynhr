// Matrix-free three-asset transition application.
//
// Policy arrays arrive in R's native (e,d,f,a) column-major layout. State
// vectors use the package distribution order (e slowest, then d, f, a with a
// fastest), matching hank_forward_operator3().

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <array>
#include <functional>
#include <atomic>
#include <barrier>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <string>
#include <vector>
using namespace Rcpp;

static inline std::size_t pidx3(int e, int d, int f, int a,
                                int ne, int nd, int nf) {
  return e + (std::size_t)ne * (d + nd * (f + nf * a));
}

static inline std::size_t sidx3(int e, int d, int f, int a,
                                int nd, int nf, int na) {
  return a + (std::size_t)na * (f + nf * (d + nd * e));
}

struct Lottery3 {
  int lo, hi;
  double plo, phi;
};

// Raw-pointer form so the same lottery serves both the exported entry points
// (which hold NumericVectors) and the fused sweep's plain-memory buffers.
static inline Lottery3 lottery3(const double* grid, int n, double z) {
  if (n == 1) return Lottery3{0, 0, 1.0, 0.0};
  int lo = std::upper_bound(grid, grid + n, z) - grid - 1;
  if (lo < 0) lo = 0;
  if (lo > n - 2) lo = n - 2;
  double p = (grid[lo + 1] - z) / (grid[lo + 1] - grid[lo]);
  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;
  return Lottery3{lo, lo + 1, p, 1.0 - p};
}

// Apply the joint Young operator without materializing it.
// transpose = FALSE: Lambda %*% x (backward expectation application).
// transpose = TRUE:  t(Lambda) %*% x (forward distribution application).
// [[Rcpp::export]]
NumericVector hank_forward_apply3_cpp(NumericVector d_pol,
                                      NumericVector f_pol,
                                      NumericVector a_pol,
                                      NumericVector d_grid,
                                      NumericVector f_grid,
                                      NumericVector a_grid,
                                      NumericMatrix Pi,
                                      NumericVector x,
                                      bool transpose = false) {
  IntegerVector dm = d_pol.attr("dim");
  if (dm.size() != 4) stop("hank_forward_apply3_cpp: policies need four dimensions");
  const int ne = dm[0], nd = dm[1], nf = dm[2], na = dm[3];
  const std::size_t N = (std::size_t)ne * nd * nf * na;
  if ((std::size_t)x.size() != N) stop("hank_forward_apply3_cpp: x has wrong length");
  NumericVector out(N);

  for (int e = 0; e < ne; ++e)
    for (int d = 0; d < nd; ++d)
      for (int f = 0; f < nf; ++f)
        for (int a = 0; a < na; ++a) {
          const std::size_t po = pidx3(e,d,f,a,ne,nd,nf);
          const std::size_t from = sidx3(e,d,f,a,nd,nf,na);
          Lottery3 ld = lottery3(d_grid.begin(), (int)d_grid.size(), d_pol[po]);
          Lottery3 lf = lottery3(f_grid.begin(), (int)f_grid.size(), f_pol[po]);
          Lottery3 la = lottery3(a_grid.begin(), (int)a_grid.size(), a_pol[po]);
          const int di[2] = {ld.lo, ld.hi}; const double dw[2] = {ld.plo, ld.phi};
          const int fi[2] = {lf.lo, lf.hi}; const double fw[2] = {lf.plo, lf.phi};
          const int ai[2] = {la.lo, la.hi}; const double aw[2] = {la.plo, la.phi};
          for (int ep = 0; ep < ne; ++ep)
            for (int jd = 0; jd < 2; ++jd)
              for (int jf = 0; jf < 2; ++jf)
                for (int ja = 0; ja < 2; ++ja) {
                  double pr = Pi(e,ep) * dw[jd] * fw[jf] * aw[ja];
                  if (pr == 0.0) continue;
                  std::size_t to = sidx3(ep,di[jd],fi[jf],ai[ja],nd,nf,na);
                  if (transpose) out[to] += pr * x[from];
                  else out[from] += pr * x[to];
                }
        }
  return out;
}

// Apply the central directional derivative of the transposed Young operator
// to a distribution without constructing either perturbed sparse matrix.
//
// Pi_p / Pi_m (both NULL by default) are the income transition matrices used
// by the PLUS and MINUS legs. They exist for the transition-probability
// Jacobian columns (HANK+SAM f/s): a date-s perturbation of such an input
// moves Pi_s, which enters Lambda_s DIRECTLY as well as through the policy, so
// curly-D at the shock date is the JOINT (policy, Pi) directional derivative.
// Supplying Pi(x + delta) / Pi(x - delta) alongside the policy directions --
// same step, same direction -- is what makes it joint. Left NULL both legs use
// Pi and the result is bit-identical to the pre-existing policy-only form.
// [[Rcpp::export]]
NumericVector hank_forward_direction3_cpp(NumericVector d_pol,
                                          NumericVector f_pol,
                                          NumericVector a_pol,
                                          NumericVector dd,
                                          NumericVector df,
                                          NumericVector da,
                                          NumericVector d_grid,
                                          NumericVector f_grid,
                                          NumericVector a_grid,
                                          NumericMatrix Pi,
                                          NumericVector dist,
                                          double delta,
                                          Nullable<NumericMatrix> Pi_p = R_NilValue,
                                          Nullable<NumericMatrix> Pi_m = R_NilValue) {
  IntegerVector dm = d_pol.attr("dim");
  if (dm.size() != 4) stop("hank_forward_direction3_cpp: policies need four dimensions");
  const int ne = dm[0], nd = dm[1], nf = dm[2], na = dm[3];
  const std::size_t N = (std::size_t)ne * nd * nf * na;
  if ((std::size_t)dist.size() != N || (std::size_t)dd.size() != N ||
      (std::size_t)df.size() != N || (std::size_t)da.size() != N)
    stop("hank_forward_direction3_cpp: state arrays have wrong length");
  if (!R_finite(delta) || delta <= 0.0)
    stop("hank_forward_direction3_cpp: delta must be finite and positive");
  NumericMatrix Pip = Pi_p.isNotNull() ? NumericMatrix(Pi_p.get()) : Pi;
  NumericMatrix Pim = Pi_m.isNotNull() ? NumericMatrix(Pi_m.get()) : Pi;
  if (Pip.nrow() != ne || Pip.ncol() != ne ||
      Pim.nrow() != ne || Pim.ncol() != ne)
    stop("hank_forward_direction3_cpp: Pi_p/Pi_m must be n_e x n_e");
  // Take the joint arithmetic path only when the two legs GENUINELY differ,
  // not merely when the arguments were supplied. An INERT transition input --
  // one its Pi_fn ignores, the transition analogue of px on a singleton-f
  // reduction -- produces Pi_p == Pi_m bitwise, and must yield an exactly
  // zero column, not a contraction-level 1e-12 one. Comparing here is O(ne^2)
  // once against an O(N*ne*8) loop, so it costs nothing and buys the exact
  // zero.
  bool joint = false;
  for (int i = 0; i < ne && !joint; ++i)
    for (int j = 0; j < ne; ++j)
      if (Pip(i, j) != Pim(i, j)) { joint = true; break; }
  NumericVector out(N);
  const double scale = 0.5 / delta;

  for (int e = 0; e < ne; ++e)
    for (int d = 0; d < nd; ++d)
      for (int f = 0; f < nf; ++f)
        for (int a = 0; a < na; ++a) {
          const std::size_t po = pidx3(e,d,f,a,ne,nd,nf);
          const std::size_t from = sidx3(e,d,f,a,nd,nf,na);
          Lottery3 ldp = lottery3(d_grid.begin(), (int)d_grid.size(), d_pol[po] + delta * dd[po]);
          Lottery3 lfp = lottery3(f_grid.begin(), (int)f_grid.size(), f_pol[po] + delta * df[po]);
          Lottery3 lap = lottery3(a_grid.begin(), (int)a_grid.size(), a_pol[po] + delta * da[po]);
          Lottery3 ldm = lottery3(d_grid.begin(), (int)d_grid.size(), d_pol[po] - delta * dd[po]);
          Lottery3 lfm = lottery3(f_grid.begin(), (int)f_grid.size(), f_pol[po] - delta * df[po]);
          Lottery3 lam = lottery3(a_grid.begin(), (int)a_grid.size(), a_pol[po] - delta * da[po]);
          const int dpi[2] = {ldp.lo, ldp.hi}; const double dpw[2] = {ldp.plo, ldp.phi};
          const int fpi[2] = {lfp.lo, lfp.hi}; const double fpw[2] = {lfp.plo, lfp.phi};
          const int api[2] = {lap.lo, lap.hi}; const double apw[2] = {lap.plo, lap.phi};
          const int dmi[2] = {ldm.lo, ldm.hi}; const double dmw[2] = {ldm.plo, ldm.phi};
          const int fmi[2] = {lfm.lo, lfm.hi}; const double fmw[2] = {lfm.plo, lfm.phi};
          const int ami[2] = {lam.lo, lam.hi}; const double amw[2] = {lam.plo, lam.phi};
          const double mass = dist[from] * scale;
          for (int ep = 0; ep < ne; ++ep)
            for (int jd = 0; jd < 2; ++jd)
              for (int jf = 0; jf < 2; ++jf)
                for (int ja = 0; ja < 2; ++ja) {
                  // Per-LEG transition weights. With Pi_p == Pi_m == Pi these
                  // are the same number and everything below reduces to the
                  // original policy-only expression exactly; when they differ
                  // the gap between them IS the direct-Pi term.
                  const double bp = mass * Pip(e,ep);
                  const double bm = mass * Pim(e,ep);
                  const double wp = dpw[jd] * fpw[jf] * apw[ja];
                  const double wm = dmw[jd] * fmw[jf] * amw[ja];
                  const std::size_t top = sidx3(ep,dpi[jd],fpi[jf],api[ja],nd,nf,na);
                  const std::size_t tom = sidx3(ep,dmi[jd],fmi[jf],ami[ja],nd,nf,na);
                  // Combine coincident destinations before accumulation, which
                  // preserves exact structural zeros (notably singleton f).
                  //
                  // The two branches are NOT a micro-optimization. On the
                  // policy-only path the legs share a transition matrix, so
                  // factoring it out cancels wp - wm EXACTLY; writing that as
                  // bp*wp - bm*wm instead lets the compiler contract one
                  // product into an FMA while rounding the other, leaving a
                  // ~1e-17 relative residue that the 0.5/delta scale lifts to
                  // ~1e-12 -- enough to destroy the exact zero this branch
                  // exists to protect, and with it the guarantee that adding
                  // the Pi arguments left every price column bit-identical.
                  //
                  // Once the legs carry DIFFERENT transition matrices the
                  // factoring is not available: it would cancel the direct-Pi
                  // derivative to zero wherever the destinations coincide,
                  // which on a singleton-f reduction (or any unmoved
                  // coordinate) is everywhere -- i.e. exactly the term a
                  // transition-probability column is made of.
                  if (top == tom) {
                    out[top] += joint ? (bp * wp - bm * wm) : (bp * (wp - wm));
                  } else {
                    if (wp != 0.0) out[top] += bp * wp;
                    if (wm != 0.0) out[tom] -= bm * wm;
                  }
                }
        }
  return out;
}

// Shared body of hank_forward_legs3_cpp(), on raw memory so the fused sweep
// (hank_curly_sweep3_cpp) can run the SAME accumulation on its plain buffers.
// `out` must be zero-initialized by the caller (both callers hand freshly
// allocated R memory, which is). The exact-zero contracts live HERE and only
// here: the bitwise Pi_plus == Pi_minus `joint` flag, and the coincident-
// destination `bp * (wp - wm)` branch -- see the comments inside, and do not
// duplicate this loop anywhere else.
static void forward_legs3_core(const double* d_plus, const double* f_plus,
                               const double* a_plus, const double* d_minus,
                               const double* f_minus, const double* a_minus,
                               const double* dg, int ndg,
                               const double* fg, int nfg,
                               const double* ag, int nag,
                               const double* Pip, const double* Pim,
                               int ne, int nd, int nf, int na,
                               const double* dist, double step, double* out) {
  bool joint = false;
  for (int i = 0; i < ne && !joint; ++i)
    for (int j = 0; j < ne; ++j)
      if (Pip[i + (std::size_t)ne * j] != Pim[i + (std::size_t)ne * j]) {
        joint = true; break;
      }

  const double scale = 0.5 / step;
  for (int e = 0; e < ne; ++e)
    for (int d = 0; d < nd; ++d)
      for (int f = 0; f < nf; ++f)
        for (int a = 0; a < na; ++a) {
          const std::size_t po = pidx3(e,d,f,a,ne,nd,nf);
          const std::size_t from = sidx3(e,d,f,a,nd,nf,na);
          Lottery3 ldp = lottery3(dg, ndg, d_plus[po]);
          Lottery3 lfp = lottery3(fg, nfg, f_plus[po]);
          Lottery3 lap = lottery3(ag, nag, a_plus[po]);
          Lottery3 ldm = lottery3(dg, ndg, d_minus[po]);
          Lottery3 lfm = lottery3(fg, nfg, f_minus[po]);
          Lottery3 lam = lottery3(ag, nag, a_minus[po]);
          const int dpi[2] = {ldp.lo, ldp.hi};
          const double dpw[2] = {ldp.plo, ldp.phi};
          const int fpi[2] = {lfp.lo, lfp.hi};
          const double fpw[2] = {lfp.plo, lfp.phi};
          const int api[2] = {lap.lo, lap.hi};
          const double apw[2] = {lap.plo, lap.phi};
          const int dmi[2] = {ldm.lo, ldm.hi};
          const double dmw[2] = {ldm.plo, ldm.phi};
          const int fmi[2] = {lfm.lo, lfm.hi};
          const double fmw[2] = {lfm.plo, lfm.phi};
          const int ami[2] = {lam.lo, lam.hi};
          const double amw[2] = {lam.plo, lam.phi};
          const double mass = dist[from] * scale;
          for (int ep = 0; ep < ne; ++ep)
            for (int jd = 0; jd < 2; ++jd)
              for (int jf = 0; jf < 2; ++jf)
                for (int ja = 0; ja < 2; ++ja) {
                  const double bp = mass * Pip[e + (std::size_t)ne * ep];
                  const double bm = mass * Pim[e + (std::size_t)ne * ep];
                  const double wp = dpw[jd] * fpw[jf] * apw[ja];
                  const double wm = dmw[jd] * fmw[jf] * amw[ja];
                  const std::size_t top =
                    sidx3(ep,dpi[jd],fpi[jf],api[ja],nd,nf,na);
                  const std::size_t tom =
                    sidx3(ep,dmi[jd],fmi[jf],ami[ja],nd,nf,na);
                  // Combine coincident destinations before accumulation, which
                  // preserves exact structural zeros (notably singleton f).
                  // With a SHARED transition matrix, bp * (wp - wm) cancels
                  // exactly; bp*wp - bm*wm would let the compiler contract one
                  // product into an FMA and leave a ~1e-17 residue that the
                  // 0.5/step scale amplifies. With genuinely different Pi legs
                  // the factoring would instead cancel the direct-Pi term to
                  // zero wherever destinations coincide -- exactly the term a
                  // transition-probability column is made of. Both branches
                  // are load-bearing; see hank_forward_direction3_cpp.
                  if (top == tom) {
                    out[top] += joint ? (bp * wp - bm * wm)
                                      : (bp * (wp - wm));
                  } else {
                    if (wp != 0.0) out[top] += bp * wp;
                    if (wm != 0.0) out[tom] -= bm * wm;
                  }
                }
        }
}

// Apply the central difference of two ACTUAL policy/transition legs of the
// transposed Young operator to a distribution. Unlike
// hank_forward_direction3_cpp(), this does not reconstruct a symmetric pair
// around the steady policy. That distinction is binding at active asset
// bounds: reconstructing steady +/- delta*dpolicy can clip one leg after the
// local policy derivative has already averaged an asymmetric active-set
// response, breaking exact first-moment accounting.
// [[Rcpp::export]]
NumericVector hank_forward_legs3_cpp(NumericVector d_plus,
                                     NumericVector f_plus,
                                     NumericVector a_plus,
                                     NumericVector d_minus,
                                     NumericVector f_minus,
                                     NumericVector a_minus,
                                     NumericVector d_grid,
                                     NumericVector f_grid,
                                     NumericVector a_grid,
                                     NumericMatrix Pi_plus,
                                     NumericMatrix Pi_minus,
                                     NumericVector dist,
                                     double step) {
  IntegerVector dm = d_plus.attr("dim");
  if (dm.size() != 4)
    stop("hank_forward_legs3_cpp: policies need four dimensions");
  const int ne = dm[0], nd = dm[1], nf = dm[2], na = dm[3];
  const std::size_t N = (std::size_t)ne * nd * nf * na;
  const NumericVector policies[6] =
    {d_plus, f_plus, a_plus, d_minus, f_minus, a_minus};
  for (int k = 0; k < 6; ++k) {
    IntegerVector dk = policies[k].attr("dim");
    if (dk.size() != 4 || dk[0] != ne || dk[1] != nd ||
        dk[2] != nf || dk[3] != na || (std::size_t)policies[k].size() != N)
      stop("hank_forward_legs3_cpp: policy legs must be conformable");
  }
  if ((std::size_t)dist.size() != N)
    stop("hank_forward_legs3_cpp: dist has wrong length");
  if (!R_finite(step) || step <= 0.0)
    stop("hank_forward_legs3_cpp: step must be finite and positive");
  if (Pi_plus.nrow() != ne || Pi_plus.ncol() != ne ||
      Pi_minus.nrow() != ne || Pi_minus.ncol() != ne)
    stop("hank_forward_legs3_cpp: Pi legs must be n_e x n_e");

  NumericVector out(N);
  forward_legs3_core(d_plus.begin(), f_plus.begin(), a_plus.begin(),
                     d_minus.begin(), f_minus.begin(), a_minus.begin(),
                     d_grid.begin(), (int)d_grid.size(),
                     f_grid.begin(), (int)f_grid.size(),
                     a_grid.begin(), (int)a_grid.size(),
                     Pi_plus.begin(), Pi_minus.begin(),
                     ne, nd, nf, na, dist.begin(), step, out.begin());
  return out;
}

static inline void psi3(double ap, double a, double r, double c0,
                        double c1, double c2, double *P, double *P1,
                        double *P2) {
  const double ar = (1.0 + r) * a;
  const double dx = ap - ar, adx = std::fabs(dx);
  const double den = ar + c0;
  const double core = R_pow(adx / den, c2 - 1.0);
  const double p = (c1 / c2) * adx * core;
  const double p1 = c1 * ((dx > 0.0) - (dx < 0.0)) * core;
  if (P) *P = p;
  if (P1) *P1 = p1;
  if (P2) *P2 = -(1.0 + r) * (p1 + (c2 - 1.0) * p / den);
}

// Takes a raw pointer rather than a std::vector reference so the expectation
// slabs can be read from worker threads without touching any R or container
// machinery.
static inline double interp3_2(const double *Z, int e, int d,
                               double f, double a, const double *fg,
                               const double *ag, int ne, int nd, int nf,
                               int na) {
  int jf = std::upper_bound(fg, fg + nf, f) - fg - 1;
  int ja = std::upper_bound(ag, ag + na, a) - ag - 1;
  jf = std::max(0, std::min(jf, nf - 2));
  ja = std::max(0, std::min(ja, na - 2));
  const double wf = (f - fg[jf]) / (fg[jf + 1] - fg[jf]);
  const double wa = (a - ag[ja]) / (ag[ja + 1] - ag[ja]);
  const std::size_t o00 = pidx3(e,d,jf,ja,ne,nd,nf);
  const std::size_t o10 = pidx3(e,d,jf+1,ja,ne,nd,nf);
  const std::size_t o01 = pidx3(e,d,jf,ja+1,ne,nd,nf);
  const std::size_t o11 = pidx3(e,d,jf+1,ja+1,ne,nd,nf);
  return (1-wf)*(1-wa)*Z[o00] + wf*(1-wa)*Z[o10] +
         (1-wf)*wa*Z[o01] + wf*wa*Z[o11];
}

using Foc3 = std::function<std::array<double,2>(double,double)>;

static inline bool valid_kkt3(const Foc3& foc, const std::array<double,2>& z,
                              const std::array<double,2>& lo,
                              const std::array<double,2>& hi) {
  const auto g = foc(z[0], z[1]);
  for (int j=0;j<2;++j) {
    const bool in = z[j] > lo[j] + 1e-6 && z[j] < hi[j] - 1e-6;
    if (in) { if (std::fabs(g[j]) > 1e-5) return false; }
    else if (z[j] <= lo[j] + 1e-6) { if (g[j] > 1e-5) return false; }
    else if (g[j] < -1e-5) return false;
  }
  return R_finite(g[0]) && R_finite(g[1]);
}

static inline double golden3(const std::function<double(double)>& objective,
                             double lo, double hi) {
  const double gr = 0.6180339887498948482;
  double x1 = hi - gr*(hi-lo), x2 = lo + gr*(hi-lo);
  double f1 = objective(x1), f2 = objective(x2);
  for (int it=0; it<200 && hi-lo > 1e-13; ++it) {
    if (f1 <= f2) { hi=x2; x2=x1; f2=f1; x1=hi-gr*(hi-lo); f1=objective(x1); }
    else { lo=x1; x1=x2; f1=f2; x2=lo+gr*(hi-lo); f2=objective(x2); }
  }
  return f1 <= f2 ? x1 : x2;
}

// `status` replaces what used to be an Rcpp::stop() from inside this routine.
// The active-set search runs on worker threads, and R's error mechanism
// (longjmp out of a C++ throw) is main-thread-only: unwinding from a worker is
// undefined behaviour, not a caught error. Failures are therefore reported as
// a code, propagated up, and re-raised as an R error by the main thread after
// every worker has joined. 0 = ok, 2 = no active-set candidate satisfies KKT.
static std::array<double,2> active_foc3(const Foc3& foc,
                                        std::array<double,2> z,
                                        const std::array<double,2>& lo,
                                        const std::array<double,2>& hi,
                                        int *status) {
  for (int j=0;j<2;++j) z[j]=std::max(lo[j],std::min(hi[j],z[j]));
  // Safeguarded finite-difference Newton, followed by the reference projected
  // ascent iteration. The active-set enumeration below is the binding guard.
  for (int it=0; it<80; ++it) {
    auto g=foc(z[0],z[1]);
    if (!R_finite(g[0]) || !R_finite(g[1])) break;
    if (std::max(std::fabs(g[0]),std::fabs(g[1])) < 1e-10) break;
    double J[2][2];
    for (int j=0;j<2;++j) {
      double h=1e-6*std::max(1.0,std::fabs(z[j]));
      auto zp=z, zm=z; zp[j]=std::min(hi[j],z[j]+h); zm[j]=std::max(lo[j],z[j]-h);
      auto gp=foc(zp[0],zp[1]), gm=foc(zm[0],zm[1]);
      const double den=zp[j]-zm[j];
      J[0][j]=(gp[0]-gm[0])/den; J[1][j]=(gp[1]-gm[1])/den;
    }
    const double det=J[0][0]*J[1][1]-J[0][1]*J[1][0];
    if (!R_finite(det)||std::fabs(det)<1e-14) break;
    const std::array<double,2> dz = {
      (J[1][1]*g[0]-J[0][1]*g[1])/det,
      (-J[1][0]*g[0]+J[0][0]*g[1])/det};
    const double old=std::max(std::fabs(g[0]),std::fabs(g[1]));
    bool moved=false;
    for(double step=1.0;step>=1.0/128.0;step*=.5) {
      std::array<double,2> zn={std::max(lo[0],std::min(hi[0],z[0]-step*dz[0])),
                               std::max(lo[1],std::min(hi[1],z[1]-step*dz[1]))};
      auto gn=foc(zn[0],zn[1]);
      if(R_finite(gn[0])&&R_finite(gn[1])&&std::max(std::fabs(gn[0]),std::fabs(gn[1]))<old){z=zn;moved=true;break;}
    }
    if(!moved) break;
  }
  for(int it=0;it<1000&&!valid_kkt3(foc,z,lo,hi);++it){
    auto g=foc(z[0],z[1]); std::array<double,2> zn;
    for(int j=0;j<2;++j)zn[j]=std::max(lo[j],std::min(hi[j],z[j]+.05*g[j]));
    if(std::max(std::fabs(zn[0]-z[0]),std::fabs(zn[1]-z[1]))<1e-7){z=zn;break;} z=zn;
  }
  std::vector<std::array<double,2>> cand;
  if(valid_kkt3(foc,z,lo,hi)) {
    const bool interior=z[0]>lo[0]+1e-6&&z[0]<hi[0]-1e-6&&
                        z[1]>lo[1]+1e-6&&z[1]<hi[1]-1e-6;
    if(interior) return z;
    cand.push_back(z);
  }
  for(double ff: {lo[0],hi[0]}) {
    double aa=golden3([&](double q){auto g=foc(ff,q);return g[1]*g[1];},lo[1],hi[1]);
    cand.push_back({ff,aa});
  }
  for(double aa: {lo[1],hi[1]}) {
    double ff=golden3([&](double q){auto g=foc(q,aa);return g[0]*g[0];},lo[0],hi[0]);
    cand.push_back({ff,aa});
  }
  for(double ff:{lo[0],hi[0]})for(double aa:{lo[1],hi[1]})cand.push_back({ff,aa});
  bool found=false; double best=R_PosInf; std::array<double,2> ans=z;
  for(auto q:cand)if(valid_kkt3(foc,q,lo,hi)){auto g=foc(q[0],q[1]);double v=std::fabs(g[0])+std::fabs(g[1]);if(v<best){best=v;ans=q;found=true;}}
  if(!found&&status) *status=2;
  return ans;
}

// Everything one (e, f0, a0) triple needs. Raw pointers only: a worker thread
// may not construct, destroy or index an R object, so all Rcpp interaction
// stays on the main thread and the parallel region sees plain memory.
namespace {
struct Egm3Ctx {
  const double *dg,*fg,*ag,*y,*Ed,*Ef,*Ea;
  int ne,nd,nf,na;
  double rd,rf,ra,px,beta,eis,chi0,chi1,chi2,phi0,phi1,phi2;
  double *D,*F,*A,*C,*Chi,*Phi,*Vdn,*Vfn,*Van;
};
}

// One (e, f0, a0) triple: solves the nd endogenous-grid problems, then maps
// back onto the exogenous liquid grid. Writes only to indices
// pidx3(e, ., f0, a0), which are disjoint across triples -- that disjointness,
// plus read-only expectation slabs, is what makes the loop over triples
// parallel with no synchronisation and BIT-IDENTICAL output: each element is
// produced by exactly the same operations in the same order regardless of how
// the triples are distributed over threads.
//
// Scratch (de/fe/ae/ce) is passed in per worker rather than shared: it used to
// be hoisted outside the loop, which is correct serially but is a data race
// the moment two triples run at once.
//
// Returns 0 ok, 1 non-monotone endogenous grid, 2 no valid active-set candidate.
static int egm3_triple(const Egm3Ctx& X, int e, int f0, int a0,
                       std::vector<double>& de, std::vector<double>& fe,
                       std::vector<double>& ae, std::vector<double>& ce) {
  const double *dg=X.dg,*fg=X.fg,*ag=X.ag,*y=X.y;
  const int ne=X.ne,nd=X.nd,nf=X.nf,na=X.na;
  const double rd=X.rd,rf=X.rf,ra=X.ra,px=X.px,beta=X.beta,eis=X.eis;
  const double chi0=X.chi0,chi1=X.chi1,chi2=X.chi2;
  const double phi0=X.phi0,phi1=X.phi1,phi2=X.phi2;
  const std::array<double,2> lo={fg[0],ag[0]},hi={fg[nf-1],ag[na-1]};
  int status=0;
  std::array<double,2> start={fg[f0],ag[a0]};
  for(int d=0;d<nd;++d){
    Foc3 foc=[&](double ff,double aa){
      double wd=interp3_2(X.Ed,e,d,ff,aa,fg,ag,ne,nd,nf,na),p1f,p1a;
      psi3(ff,fg[f0],rf,phi0,phi1,phi2,nullptr,&p1f,nullptr);
      psi3(aa,ag[a0],ra,chi0,chi1,chi2,nullptr,&p1a,nullptr);
      // Foreign claims trade at the world price px, so a unit of f' costs px
      // domestic goods plus the marginal adjustment cost -- hence `-px-p1f`
      // where a unit-priced asset would give `-1-p1f`. px enters the PURCHASE
      // side, not only the revaluation of the existing stock: the paper's
      // Stage-4 contract has foreign market clearing determine px, and a price
      // that only revalued the predetermined position could not clear a market
      // for new purchases. It would also be exactly collinear with rf (both
      // multiplying f alone), making it a redundant Jacobian input.
      return std::array<double,2>{interp3_2(X.Ef,e,d,ff,aa,fg,ag,ne,nd,nf,na)/wd-px-p1f,
        interp3_2(X.Ea,e,d,ff,aa,fg,ag,ne,nd,nf,na)/wd-1-p1a};};
    auto z=active_foc3(foc,start,lo,hi,&status);
    if(status) return status;
    start=z;
    double wd=interp3_2(X.Ed,e,d,z[0],z[1],fg,ag,ne,nd,nf,na);
    double cc=R_pow(beta*wd,-eis),pf,pa;psi3(z[0],fg[f0],rf,phi0,phi1,phi2,&pf,nullptr,nullptr);psi3(z[1],ag[a0],ra,chi0,chi1,chi2,&pa,nullptr,nullptr);
    de[d]=(cc+dg[d]+px*z[0]+z[1]+pf+pa-y[e]-px*(1+rf)*fg[f0]-(1+ra)*ag[a0])/(1+rd);
    fe[d]=z[0];ae[d]=z[1];ce[d]=cc;
  }
  for(int d=1;d<nd;++d)if(de[d]<=de[d-1])return 1;
  start={fg[f0],ag[a0]};
  for(int id=0;id<nd;++id){
    double dp,fp,ap,cp;
    if(dg[id]<de[0]){
      Foc3 foc=[&](double ff,double aa){double pf,p1f,pa,p1a;psi3(ff,fg[f0],rf,phi0,phi1,phi2,&pf,&p1f,nullptr);psi3(aa,ag[a0],ra,chi0,chi1,chi2,&pa,&p1a,nullptr);double cc=y[e]+(1+rd)*dg[id]+px*(1+rf)*fg[f0]+(1+ra)*ag[a0]-dg[0]-px*ff-aa-pf-pa;if(!R_finite(cc)||cc<=0)return std::array<double,2>{1e12,1e12};double uc=R_pow(cc,-1/eis);return std::array<double,2>{beta*interp3_2(X.Ef,e,0,ff,aa,fg,ag,ne,nd,nf,na)/uc-px-p1f,beta*interp3_2(X.Ea,e,0,ff,aa,fg,ag,ne,nd,nf,na)/uc-1-p1a};};
      auto z=active_foc3(foc,start,lo,hi,&status);
      if(status) return status;
      start=z;double pf,pa;psi3(z[0],fg[f0],rf,phi0,phi1,phi2,&pf,nullptr,nullptr);psi3(z[1],ag[a0],ra,chi0,chi1,chi2,&pa,nullptr,nullptr);dp=dg[0];fp=z[0];ap=z[1];cp=y[e]+(1+rd)*dg[id]+px*(1+rf)*fg[f0]+(1+ra)*ag[a0]-dp-px*fp-ap-pf-pa;
    } else {
      if(dg[id]>=de[nd-1]){dp=dg[nd-1];fp=fe[nd-1];ap=ae[nd-1];}
      else {int j=std::upper_bound(de.begin(),de.end(),dg[id])-de.begin()-1;double w=(dg[id]-de[j])/(de[j+1]-de[j]);dp=dg[j]+w*(dg[j+1]-dg[j]);fp=fe[j]+w*(fe[j+1]-fe[j]);ap=ae[j]+w*(ae[j+1]-ae[j]);}
      // Rebuild c from the BUDGET at the interpolated portfolio, rather than
      // interpolating the endogenous-grid c alongside it. The budget is
      // nonlinear in (f', a') through Phi and Psi, so a fourth independent
      // interpolation does not satisfy it off a source knot -- measured 1.98e-4
      // on the 54-state reference, shrinking only as the grid refines, i.e. a
      // real discretisation error rather than solver noise. The top branch had
      // the same defect for a different reason: it paired ce[nd-1] with dg[id]
      // != de[nd-1]. hank_egm2_solve has always done it this way (c_pol is
      // formed from the budget AFTER b_pol/a_pol are fixed); this brings the
      // three-asset kernel onto the same convention, so the household budget
      // holds by construction at every state and the aggregate resource
      // identity is no longer capped at 1e-4.
      double pf2,pa2;
      psi3(fp,fg[f0],rf,phi0,phi1,phi2,&pf2,nullptr,nullptr);
      psi3(ap,ag[a0],ra,chi0,chi1,chi2,&pa2,nullptr,nullptr);
      cp=y[e]+(1+rd)*dg[id]+px*(1+rf)*fg[f0]+(1+ra)*ag[a0]-dp-px*fp-ap-pf2-pa2;
    }
    // The adjustment technologies are REAL RESOURCE costs, so the block has to
    // be able to report them: they are the wedge between what households
    // finance and what the goods market must supply. Recovered here from the
    // same psi3 call that already produces Psi2 for the envelope, at the ONE
    // place where (fp, ap) are final -- recomputing them in the caller would
    // re-introduce exactly the off-policy evaluation the budget fix removed.
    const std::size_t o=pidx3(e,id,f0,a0,ne,nd,nf);X.D[o]=dp;X.F[o]=fp;X.A[o]=ap;X.C[o]=cp;double p2f,p2a,pfc,pac;psi3(fp,fg[f0],rf,phi0,phi1,phi2,&pfc,nullptr,&p2f);psi3(ap,ag[a0],ra,chi0,chi1,chi2,&pac,nullptr,&p2a);X.Phi[o]=pfc;X.Chi[o]=pac;double uc=R_pow(cp,-1/eis);X.Vdn[o]=uc*(1+rd);X.Vfn[o]=uc*(px*(1+rf)-p2f);X.Van[o]=uc*((1+ra)-p2a);
  }
  return 0;
}

// Expectation contraction (main thread only; fills the Ed/Ef/Ea slabs the
// triples read). Extracted verbatim from hank_egm3_step_cpp so the fused
// sweep runs the SAME computation -- Pi is column-major with Pi(e, ep) =
// Pi[e + ne*ep], exactly as NumericMatrix indexes it.
static void egm3_expect(const double* Vd, const double* Vf, const double* Va,
                        const double* Pi, int ne, int nd, int nf, int na,
                        double* Ed, double* Ef, double* Ea) {
  for(int e=0;e<ne;++e)for(int d=0;d<nd;++d)for(int f=0;f<nf;++f)for(int a=0;a<na;++a){
    const std::size_t o=pidx3(e,d,f,a,ne,nd,nf); double vd=0,vf=0,va=0;
    for(int ep=0;ep<ne;++ep){const std::size_t op=pidx3(ep,d,f,a,ne,nd,nf);
      vd+=Pi[e+(std::size_t)ne*ep]*Vd[op];vf+=Pi[e+(std::size_t)ne*ep]*Vf[op];va+=Pi[e+(std::size_t)ne*ep]*Va[op];}
    Ed[o]=vd;Ef[o]=vf;Ea[o]=va;
  }
}

// Shared failure channel for the triple loop. Codes: 0 ok, 1 non-monotone
// endogenous grid, 2 no valid active-set candidate, 3 a C++ exception escaped
// the task body (bad_alloc from the per-worker scratch, an Armadillo
// logic_error, ...). Code 3 exists because an exception that leaves a
// std::thread body calls std::terminate -- the R session dies with no error at
// all -- and because Rcpp::stop() from a worker would be worse still: R's
// error mechanism longjmps, which is undefined behaviour off the main thread.
namespace {
struct Egm3Err {
  std::atomic<int> code;
  std::atomic<std::size_t> task;
  std::mutex mu;
  std::string what;                 // populated only for code 3
  Egm3Err() : code(0), task(0) {}
  void set(int c, std::size_t t) {
    int expect = 0;
    if (code.compare_exchange_strong(expect, c)) task.store(t);
  }
  void set_exception(const char* w, std::size_t t) {
    try {
      std::lock_guard<std::mutex> lk(mu);
      int expect = 0;
      if (code.compare_exchange_strong(expect, 3)) { task.store(t); what = w; }
    } catch (...) {
      set(3, t);                    // even copying the message can fail
    }
  }
};
}

// The triple loop over a contiguous task range. Shared by the per-call
// spawn path (hank_egm3_step_cpp) and the persistent pool (Egm3Pool): each
// task writes a disjoint pidx3(e, ., f0, a0) slice, so output is bit-identical
// for ANY partition of [0, ntask) and any thread count. First failure wins;
// the losers' codes are discarded, which is fine because any one of them is
// enough to abort the step.
//
// NOTHING may propagate out of here: this body IS a std::thread's entry point.
static void egm3_run_tasks(const Egm3Ctx& X, std::size_t lo_t, std::size_t hi_t,
                           Egm3Err& err) {
  std::size_t t = lo_t;
  try {
    const int nd=X.nd,nf=X.nf,na=X.na;
    std::vector<double> de(nd),fe(nd),ae(nd),ce(nd);   // per-worker scratch
    for(;t<hi_t;++t){
      if(err.code.load(std::memory_order_relaxed)) return;  // abandon early
      const int a0=(int)(t%(std::size_t)na);
      const int f0=(int)((t/(std::size_t)na)%(std::size_t)nf);
      const int e=(int)(t/((std::size_t)na*(std::size_t)nf));
      const int s=egm3_triple(X,e,f0,a0,de,fe,ae,ce);
      if(s){err.set(s,t);return;}
    }
  } catch (const std::exception& e) {
    err.set_exception(e.what(), t);
  } catch (...) {
    err.set_exception("unknown C++ exception", t);
  }
}

// Errors are raised on the MAIN thread, after every worker has joined. The
// message text is the step's regardless of which entry point ran it, so the
// fused sweep fails exactly like the per-step R loop it replaces.
static void egm3_raise(const Egm3Err& err, int nf, int na) {
  const int e_code = err.code.load();
  const std::size_t t = err.task.load();
  const int a0=(int)(t%(std::size_t)na);
  const int f0=(int)((t/(std::size_t)na)%(std::size_t)nf);
  const int e=(int)(t/((std::size_t)na*(std::size_t)nf));
  if(e_code==1)
    stop("hank_egm3_step_cpp: non-monotone endogenous liquid grid at "
         "(e=%d, f=%d, a=%d)",e+1,f0+1,a0+1);
  if(e_code==3)
    stop("hank_egm3_step_cpp: worker thread failed at "
         "(e=%d, f=%d, a=%d): %s",e+1,f0+1,a0+1,
         err.what.empty() ? "C++ exception" : err.what.c_str());
  stop("hank_egm3_step_cpp: no active-set candidate satisfies joint KKT "
       "conditions at (e=%d, f=%d, a=%d)",e+1,f0+1,a0+1);
}

// Input-shape validation (B4). `dim(Vd)` used to be read with no check at all:
// a Vd whose dim attribute had been dropped by an R-level drop()/as.vector()
// yielded a zero-length IntegerVector, and dm[0..3] read four ints off the end
// of it -- garbage extents, then out-of-bounds reads over every grid. Every
// extent is checked against the grid lengths here, and the message names the
// offending argument and the shape it must have. Mirrored in R by
// .hank_egm3_step() (R/hank-egm3.R).
static void egm3_check_shapes(const char* fn,
                              const NumericVector& Vd, const NumericVector& Vf,
                              const NumericVector& Va, const NumericVector& dg,
                              const NumericVector& fg, const NumericVector& ag,
                              const NumericVector& y, const NumericMatrix& Pi,
                              int* ne_o, int* nd_o, int* nf_o, int* na_o) {
  RObject dim_attr = Vd.attr("dim");
  if (dim_attr.isNULL())
    stop("%s: `Vd` must be a four-dimensional array "
         "(n_e x n_d x n_f x n_a); it has no `dim` attribute.", fn);
  IntegerVector dm(dim_attr);
  if (dm.size() != 4)
    stop("%s: `Vd` must be a four-dimensional array "
         "(n_e x n_d x n_f x n_a); dim(Vd) has length %d.",
         fn, (int)dm.size());
  const int ne = dm[0], nd = dm[1], nf = dm[2], na = dm[3];
  if (ne < 1 || nd < 1 || nf < 1 || na < 1)
    stop("%s: every dimension of `Vd` must be >= 1; got %d x %d x %d x %d.",
         fn, ne, nd, nf, na);
  if (Pi.nrow() != ne || Pi.ncol() != ne)
    stop("%s: `Pi` must be n_e x n_e = %d x %d (n_e = dim(Vd)[1]); "
         "got %d x %d.", fn, ne, ne, (int)Pi.nrow(), (int)Pi.ncol());
  if (y.size() != ne)
    stop("%s: `y` must have length n_e = %d (= dim(Vd)[1]); got %d.",
         fn, ne, (int)y.size());
  if (dg.size() != nd)
    stop("%s: `dg` must have length n_d = %d (= dim(Vd)[2]); got %d.",
         fn, nd, (int)dg.size());
  if (fg.size() != nf)
    stop("%s: `fg` must have length n_f = %d (= dim(Vd)[3]); got %d.",
         fn, nf, (int)fg.size());
  if (ag.size() != na)
    stop("%s: `ag` must have length n_a = %d (= dim(Vd)[4]); got %d.",
         fn, na, (int)ag.size());
  if (nf < 2 || na < 2)
    stop("%s: the compiled kernel interpolates between grid points, so "
         "`fg` and `ag` need at least 2 entries each; got n_f = %d, "
         "n_a = %d.", fn, nf, na);
  const R_xlen_t N = (R_xlen_t)ne * nd * nf * na;
  if (Vf.size() != N)
    stop("%s: `Vf` must have length n_e*n_d*n_f*n_a = %d (the length of "
         "`Vd`); got %d.", fn, (int)N, (int)Vf.size());
  if (Va.size() != N)
    stop("%s: `Va` must have length n_e*n_d*n_f*n_a = %d (the length of "
         "`Vd`); got %d.", fn, (int)N, (int)Va.size());
  *ne_o = ne; *nd_o = nd; *nf_o = nf; *na_o = na;
}

// Worker count actually usable: never more tasks than exist, never more OS
// threads than the machine can run (the triple loop is compute-bound and the
// output is thread-count-invariant by construction, so this clamp cannot
// change a single bit), never fewer than 1.
static inline int egm3_nthreads(int threads, std::size_t ntask) {
  int n = std::max(1, std::min(threads, (int)ntask));
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc > 0u && n > (int)hc) n = (int)hc;
  return n < 1 ? 1 : n;
}

// One compiled three-asset EGM backward step. The R implementation remains
// the oracle and owns validation/default initialization.
//
// `threads` is resolved in R (see hank_resolve_threads); <= 1 runs the
// serial path, which is the code the parallel path must reproduce bit for bit.
// A bare step call pays ONE thread spawn/join per call; the fused sweep below
// holds one pool across all its steps instead.
// [[Rcpp::export]]
List hank_egm3_step_cpp(NumericVector Vd_, NumericVector Vf_, NumericVector Va_,
                        NumericVector dg_, NumericVector fg_, NumericVector ag_,
                        NumericVector y_, NumericMatrix Pi_, double rd, double rf,
                        double ra, double beta, double eis, double chi0,
                        double chi1, double chi2, double phi0, double phi1,
                        double phi2, double px = 1.0, int threads = 1) {
  int ne, nd, nf, na;
  egm3_check_shapes("hank_egm3_step_cpp", Vd_, Vf_, Va_, dg_, fg_, ag_, y_,
                    Pi_, &ne, &nd, &nf, &na);
  const std::size_t N=(std::size_t)ne*nd*nf*na;
  const double *dg=dg_.begin(),*fg=fg_.begin(),*ag=ag_.begin(),*y=y_.begin();
  // SERIAL: no worker exists yet, so the longjmp-free throw this raises can
  // only unwind main-thread frames.
  Rcpp::checkUserInterrupt();
  std::vector<double> Ed(N),Ef(N),Ea(N);
  egm3_expect(Vd_.begin(),Vf_.begin(),Va_.begin(),Pi_.begin(),
              ne,nd,nf,na,Ed.data(),Ef.data(),Ea.data());
  NumericVector D(N),F(N),A(N),C(N),Chi(N),Phi(N),Vdn(N),Vfn(N),Van(N);
  Egm3Ctx X{dg,fg,ag,y,Ed.data(),Ef.data(),Ea.data(),ne,nd,nf,na,
            rd,rf,ra,px,beta,eis,chi0,chi1,chi2,phi0,phi1,phi2,
            D.begin(),F.begin(),A.begin(),C.begin(),
            Chi.begin(),Phi.begin(),
            Vdn.begin(),Vfn.begin(),Van.begin()};

  const std::size_t ntask=(std::size_t)ne*nf*na;
  const int nthr=egm3_nthreads(threads,ntask);
  Egm3Err err;

  if(nthr<=1){
    egm3_run_tasks(X,0,ntask,err);
  } else {
    std::vector<std::thread> pool; pool.reserve(nthr-1);
    const std::size_t chunk=(ntask+nthr-1)/nthr;
    // Exception-safe spawn: if the k-th thread cannot be created, the ones
    // already running must be told to stop and JOINED before `pool` is
    // destroyed -- ~std::vector<std::thread> on a joinable thread calls
    // std::terminate, i.e. kills the R session.
    try {
      for(int k=1;k<nthr;++k){
        const std::size_t lo_t=std::min(ntask,(std::size_t)k*chunk);
        const std::size_t hi_t=std::min(ntask,lo_t+chunk);
        if(lo_t<hi_t)
          pool.emplace_back([&X,lo_t,hi_t,&err]{
            egm3_run_tasks(X,lo_t,hi_t,err);});
      }
    } catch (...) {
      err.set(3,0);                     // makes the live workers abandon early
      for(auto& th: pool) if(th.joinable()) th.join();
      throw;
    }
    egm3_run_tasks(X,0,std::min(ntask,chunk),err); // main takes chunk 0
    for(auto& th: pool) th.join();
  }

  if(err.code.load()) egm3_raise(err,nf,na);

  IntegerVector dims=IntegerVector::create(ne,nd,nf,na);for(auto v:{D,F,A,C,Chi,Phi,Vdn,Vfn,Van})v.attr("dim")=dims;
  return List::create(_["d"]=D,_["f"]=F,_["a"]=A,_["c"]=C,_["chi"]=Chi,_["phi"]=Phi,
                      _["Vd"]=Vdn,_["Vf"]=Vfn,_["Va"]=Van);
}

// Fused fixed-point iteration around the compiled three-asset step.
// [[Rcpp::export]]
List hank_egm3_solve_cpp(NumericVector Vd, NumericVector Vf, NumericVector Va,
                         NumericVector dg, NumericVector fg, NumericVector ag,
                         NumericVector y, NumericMatrix Pi, double rd, double rf,
                         double ra, double beta, double eis, double chi0,
                         double chi1, double chi2, double phi0, double phi1,
                         double phi2, double tol, int maxit, double relax,
                         double px = 1.0, int threads = 1) {
  Vd=clone(Vd);Vf=clone(Vf);Va=clone(Va);
  NumericVector Dold,Fold,Aold,Cold; bool have_old=false,ok=false;
  double gap=R_PosInf,pol_gap=R_PosInf; int it=0; List step;
  // A step-guard failure (non-monotone endogenous grid / no KKT candidate) is
  // RECOVERABLE, not fatal: it means this particular update overshot into an
  // infeasible region, and a damped retry from the last good marginal values
  // usually walks straight past it. Previously the Rcpp::stop() propagated out
  // of the fused loop and destroyed every completed iterate -- measured on the
  // paper's cross-rung continuation, 196 perfectly good iterations were thrown
  // away by a failure at iteration 197. So catch it, stop iterating, and hand
  // the LAST GOOD state back to R with a status, which lets hank_egm3_solve()
  // lower `relax` and continue instead of starting over.
  //
  // A failure on the FIRST iteration is rethrown: there is no good iterate to
  // return, and the caller's initial values really are unusable.
  std::string step_error; int step_status=0;
  for(it=1;it<=maxit;++it){
    // SERIAL section of the outer iteration: hank_egm3_step_cpp spawns and
    // joins its pool inside the call, so no worker is alive at this point and
    // the interrupt's C++ throw can only unwind main-thread frames. Placed
    // BEFORE the try so the step's own catch cannot swallow it (Rcpp's
    // InterruptedException does not derive from std::exception either).
    Rcpp::checkUserInterrupt();
    try {
      step=hank_egm3_step_cpp(Vd,Vf,Va,dg,fg,ag,y,Pi,rd,rf,ra,beta,eis,
                              chi0,chi1,chi2,phi0,phi1,phi2,px,threads);
    } catch (std::exception &e) {
      if (it==1) throw;
      step_error=e.what(); step_status=1; it=it-1; break;
    }
    NumericVector vn_d=step["Vd"],vn_f=step["Vf"],vn_a=step["Va"];
    NumericVector D=step["d"],F=step["f"],A=step["a"],C=step["c"];
    gap=0.0;pol_gap=have_old?0.0:R_PosInf;
    for(R_xlen_t k=0;k<Vd.size();++k){
      gap=std::max(gap,std::fabs(vn_d[k]-Vd[k]));gap=std::max(gap,std::fabs(vn_f[k]-Vf[k]));gap=std::max(gap,std::fabs(vn_a[k]-Va[k]));
      if(have_old){pol_gap=std::max(pol_gap,std::fabs(D[k]-Dold[k]));pol_gap=std::max(pol_gap,std::fabs(F[k]-Fold[k]));pol_gap=std::max(pol_gap,std::fabs(A[k]-Aold[k]));pol_gap=std::max(pol_gap,std::fabs(C[k]-Cold[k]));}
      Vd[k]=(1-relax)*Vd[k]+relax*vn_d[k];Vf[k]=(1-relax)*Vf[k]+relax*vn_f[k];Va[k]=(1-relax)*Va[k]+relax*vn_a[k];
    }
    Dold=clone(D);Fold=clone(F);Aold=clone(A);Cold=clone(C);have_old=true;
    if(gap<tol){ok=true;break;}
  }
  if(!ok && step_status==0) it=maxit;
  step["Vd"]=Vd;step["Vf"]=Vf;step["Va"]=Va;step["iterations"]=it;
  step["converged"]=ok;step["last_value_gap"]=gap;step["last_policy_gap"]=pol_gap;
  step["step_status"]=step_status;step["step_error"]=step_error;
  return step;
}

// ---------------------------------------------------------------------------
// FUSED BACKWARD SWEEP for the three-asset fake-news Jacobian.
//
// Ports the s = 2 .. T_h anticipation loop of R .hank_curly_sweep3() (see
// R/hank-jacobian3.R) so the WHOLE recursion -- both EGM step legs per date,
// the differencing/aggregation record() does, and the curly-D leg
// accumulation -- runs in one compiled call under ONE worker pool. The
// two-asset analogue is hank_curly_sweep2_cpp() in src/hank_egm2.cpp; this
// follows the same discipline (persistent pool parked on a std::barrier;
// FMA-pinned axpy; exact-max shared scale) and, like it, REUSES the step and
// forward-legs kernels (egm3_expect / egm3_run_tasks / forward_legs3_core)
// rather than re-implementing any accumulation.
// ---------------------------------------------------------------------------

// Persistent worker pool over the triple tasks. Spawned ONCE per fused sweep
// and held across every date and both FD legs; workers park on the barrier
// between steps. Partition and per-task arithmetic are egm3_run_tasks', so
// output is bit-identical to the per-call spawn path at every thread count
// (disjoint writes, no cross-thread reductions).
class Egm3Pool {
 public:
  Egm3Pool(const Egm3Ctx* ctx, std::size_t ntask, int nthr, Egm3Err* err)
    : ctx_(ctx), ntask_(ntask), nthr_(nthr),
      chunk_((ntask + nthr - 1) / nthr), err_(err),
      bar_((std::ptrdiff_t)nthr), done_(false), joined_(false) {
    pool_.reserve(nthr - 1);
    try {
      for (int w = 1; w < nthr; ++w)
        pool_.emplace_back([this, w] { this->worker(w); });
    } catch (...) {
      // The workers already created are parked at the top barrier and
      // ~std::vector<std::thread> would call std::terminate on them. Release
      // and join them first, THEN rethrow -- the destructor never runs for an
      // object whose constructor threw.
      shutdown();
      throw;
    }
  }
  // Run one backward step's triple loop. The caller must have refilled the
  // context's Ed/Ef/Ea slabs first (egm3_expect, main thread).
  void run_step() {
    err_->code.store(0);             // workers are parked; no race with this
    bar_.arrive_and_wait();          // release the workers into the tasks
    run_chunk(0);                    // main thread takes chunk 0
    bar_.arrive_and_wait();          // workers park at the top barrier again
  }
  // Releases and joins the workers. Runs on the normal path AND during stack
  // unwinding if the main thread throws between steps (egm3_raise), so a
  // worker can never outlive the buffers it points at.
  ~Egm3Pool() { shutdown(); }
  Egm3Pool(const Egm3Pool&) = delete;
  Egm3Pool& operator=(const Egm3Pool&) = delete;

 private:
  void run_chunk(int w) {
    const std::size_t lo = std::min(ntask_, (std::size_t)w * chunk_);
    const std::size_t hi = std::min(ntask_, lo + chunk_);
    // egm3_run_tasks swallows every exception into err_ (code 3); nothing can
    // propagate out of a worker body from here.
    if (lo < hi) egm3_run_tasks(*ctx_, lo, hi, *err_);
  }
  // Idempotent: the constructor's failure path and the destructor both call
  // it. Participants that were never created are DROPPED, or the barrier waits
  // forever for arrivals that can never come.
  void shutdown() {
    if (joined_) return;
    joined_ = true;
    done_.store(true, std::memory_order_release);
    for (std::size_t k = pool_.size() + 1; k < (std::size_t)nthr_; ++k)
      bar_.arrive_and_drop();
    bar_.arrive_and_wait();
    for (auto& t : pool_) t.join();
  }
  void worker(int w) {
    for (;;) {
      bar_.arrive_and_wait();
      if (done_.load(std::memory_order_acquire)) return;
      run_chunk(w);
      bar_.arrive_and_wait();
    }
  }
  const Egm3Ctx* ctx_;
  std::size_t ntask_;
  int nthr_;
  std::size_t chunk_;
  Egm3Err* err_;
  std::barrier<> bar_;
  std::atomic<bool> done_;
  bool joined_;
  std::vector<std::thread> pool_;
};

// Bit-identical port of R's `V_ss + h * dV` (see egm2_axpy_exact in
// src/hank_egm2.cpp for the full story). The multiply and the add must round
// SEPARATELY, exactly as R rounds them: clang's default -ffp-contract=on is
// licensed to fuse this very expression, moving the perturbed argument by
// 1 ULP, which the sweep's 1/(2h) amplification then COMPOUNDS date after
// date because dVd/dVf/dVa feed the next step. `volatile` on the intermediate
// product pins the rounding regardless of the contraction default. Passing -h
// reproduces R's `V_ss - h * dV` exactly: negation is exact, so v + (-(h*d))
// and v - (h*d) round identically.
static inline void egm3_axpy_exact(const double* v, double h, const double* d,
                                   double* out, std::size_t n) {
  for (std::size_t t = 0; t < n; ++t) {
    volatile double p = h * d[t];
    out[t] = v[t] + p;
  }
}

// Port of R's `max(1, max(abs(dVd)), max(abs(dVf)), max(abs(dVa)))`,
// NaN-propagating like R's max(). max() is exact (pure comparison), so
// flattening the nested maxima into one running maximum cannot change the
// result by even a ULP -- and it MUST not: this REDUCTION feeds the shared FD
// scale h that every later date inherits, so any reassociation would compound
// down the sweep. The three arrays are scanned in the same order R evaluates
// them.
static inline double egm3_sweep_scale(const double* dvd, const double* dvf,
                                      const double* dva, std::size_t n) {
  double m = 1.0;
  bool nan_seen = false;
  const double* arrs[3] = {dvd, dvf, dva};
  for (int k = 0; k < 3; ++k)
    for (std::size_t t = 0; t < n; ++t) {
      double u = std::fabs(arrs[k][t]);
      if (ISNAN(u)) { nan_seen = true; continue; }
      if (u > m) m = u;
    }
  return nan_seen ? R_NaN : m;
}

// Fused s = 2 .. T_h backward sweep (compiled).
//
// Vd_ss_/Vf_ss_/Va_ss_ are the steady-state marginals; dVd0_/dVf0_/dVa0_ the
// s = 1 directional derivatives the R caller has already computed (the s = 1
// term itself -- the input perturbation, including any Pi legs -- stays in R,
// where the input-type dispatch lives). D_ss_ is the stationary distribution
// in state order. `outputs` is the REQUESTED aggregate set (any subset of
// D/F/A/C/CHI/PHI; the internal "Y" is identically zero at s >= 2 and is left
// to the R side's zero initialization) -- the swept set is variable and must
// never be hardcoded here.
//
// Per date: shared scale h from (dVd, dVf, dVa); two EGM step legs through
// the SAME egm3_expect + egm3_run_tasks bodies as hank_egm3_step_cpp; central
// differences formed with the denominator 2*h computed ONCE as R does;
// aggregation against D_ss in state order with a long-double accumulator
// (matching R's sum()); and the curly-D column through forward_legs3_core on
// the ACTUAL legs, Pi_p = Pi_m = Pi (so the `joint` flag is false and every
// exact-zero contract of that kernel holds by construction).
//
// Returns curlyY (named list over the requested outputs, each length
// T_h - 1, entry j = date s = j + 1... i.e. dates 2..T_h) and curlyD
// (n_cell x (T_h - 1)), which the R driver splices into its own s = 1 column.
// [[Rcpp::export]]
List hank_curly_sweep3_cpp(NumericVector Vd_ss_, NumericVector Vf_ss_,
                           NumericVector Va_ss_, NumericVector dVd0_,
                           NumericVector dVf0_, NumericVector dVa0_,
                           NumericVector dg_, NumericVector fg_,
                           NumericVector ag_, NumericVector y_,
                           NumericMatrix Pi_, double rd, double rf, double ra,
                           double beta, double eis, double chi0, double chi1,
                           double chi2, double phi0, double phi1, double phi2,
                           double px, NumericVector D_ss_,
                           CharacterVector outputs, double delta_v, int T_h,
                           int threads = 1) {
  IntegerVector dm = Vd_ss_.attr("dim");
  if (dm.size() != 4) stop("hank_curly_sweep3_cpp: Vd needs four dimensions");
  const int ne = dm[0], nd = dm[1], nf = dm[2], na = dm[3];
  const std::size_t N = (std::size_t)ne * nd * nf * na;
  if ((std::size_t)Vf_ss_.size() != N || (std::size_t)Va_ss_.size() != N ||
      (std::size_t)dVd0_.size() != N || (std::size_t)dVf0_.size() != N ||
      (std::size_t)dVa0_.size() != N || (std::size_t)D_ss_.size() != N)
    stop("hank_curly_sweep3_cpp: state arrays have wrong length");
  if (Pi_.nrow() != ne || Pi_.ncol() != ne)
    stop("hank_curly_sweep3_cpp: Pi must be n_e x n_e");
  if (T_h < 2) stop("hank_curly_sweep3_cpp: T_h must be >= 2");
  const int S = T_h - 1;                      // dates s = 2 .. T_h

  bool want[6] = {false, false, false, false, false, false};
  static const char* nm[6] = {"D", "F", "A", "C", "CHI", "PHI"};
  for (int i = 0; i < outputs.size(); ++i) {
    const std::string o = as<std::string>(outputs[i]);
    for (int k = 0; k < 6; ++k) if (o == nm[k]) want[k] = true;
  }

  const double *dg = dg_.begin(), *fg = fg_.begin(), *ag = ag_.begin();
  const double *dss = D_ss_.begin();

  // Step working set: perturbed inputs, expectation slabs, step outputs.
  std::vector<double> vdp(N), vfp(N), vap(N), Ed(N), Ef(N), Ea(N);
  std::vector<double> D(N), F(N), A(N), C(N), Chi(N), Phi(N),
                      Vdn(N), Vfn(N), Van(N);
  // +h leg outputs, held while the -h leg overwrites the context buffers.
  std::vector<double> d_p(N), f_p(N), a_p(N), c_p(N), chi_p(N), phi_p(N),
                      Vd_p(N), Vf_p(N), Va_p(N);
  std::vector<double> dVd(dVd0_.begin(), dVd0_.end());
  std::vector<double> dVf(dVf0_.begin(), dVf0_.end());
  std::vector<double> dVa(dVa0_.begin(), dVa0_.end());

  Egm3Ctx X{dg, fg, ag, y_.begin(), Ed.data(), Ef.data(), Ea.data(),
            ne, nd, nf, na,
            rd, rf, ra, px, beta, eis, chi0, chi1, chi2, phi0, phi1, phi2,
            D.data(), F.data(), A.data(), C.data(), Chi.data(), Phi.data(),
            Vdn.data(), Vfn.data(), Van.data()};

  const std::size_t ntask = (std::size_t)ne * nf * na;
  const int nthr = egm3_nthreads(threads, ntask);
  Egm3Err err;
  // Spawned ONCE, held across all dates and both legs -- the whole point.
  std::unique_ptr<Egm3Pool> pool;
  if (nthr > 1) pool.reset(new Egm3Pool(&X, ntask, nthr, &err));
  auto run_leg = [&](const double* vd, const double* vf, const double* va) {
    egm3_expect(vd, vf, va, Pi_.begin(), ne, nd, nf, na,
                Ed.data(), Ef.data(), Ea.data());
    if (pool) { pool->run_step(); }
    else { err.code.store(0); egm3_run_tasks(X, 0, ntask, err); }
    if (err.code.load()) {
      pool.reset();                   // join before longjmp-ing out
      egm3_raise(err, nf, na);
    }
  };

  NumericMatrix curlyD((int)N, S);
  std::vector<NumericVector> cy(6);
  for (int k = 0; k < 6; ++k) if (want[k]) cy[k] = NumericVector(S);

  for (int s = 0; s < S; ++s) {
    // SERIAL section of the date loop: the workers are parked at the top
    // barrier, so this R-API call (and the C++ throw it raises, which unwinds
    // through `pool`'s unique_ptr and JOINS them) is safe here and nowhere
    // inside a leg.
    Rcpp::checkUserInterrupt();
    // One SHARED scale for all three marginal values: they describe a single
    // perturbation direction, and scaling them separately would change it.
    const double h = delta_v / egm3_sweep_scale(dVd.data(), dVf.data(),
                                                dVa.data(), N);
    // --- +h leg ------------------------------------------------------------
    egm3_axpy_exact(Vd_ss_.begin(), h, dVd.data(), vdp.data(), N);
    egm3_axpy_exact(Vf_ss_.begin(), h, dVf.data(), vfp.data(), N);
    egm3_axpy_exact(Va_ss_.begin(), h, dVa.data(), vap.data(), N);
    run_leg(vdp.data(), vfp.data(), vap.data());
    std::copy(D.begin(), D.end(), d_p.begin());
    std::copy(F.begin(), F.end(), f_p.begin());
    std::copy(A.begin(), A.end(), a_p.begin());
    std::copy(C.begin(), C.end(), c_p.begin());
    std::copy(Chi.begin(), Chi.end(), chi_p.begin());
    std::copy(Phi.begin(), Phi.end(), phi_p.begin());
    std::copy(Vdn.begin(), Vdn.end(), Vd_p.begin());
    std::copy(Vfn.begin(), Vfn.end(), Vf_p.begin());
    std::copy(Van.begin(), Van.end(), Va_p.begin());

    // --- -h leg (leaves its outputs in the context buffers) ----------------
    egm3_axpy_exact(Vd_ss_.begin(), -h, dVd.data(), vdp.data(), N);
    egm3_axpy_exact(Vf_ss_.begin(), -h, dVf.data(), vfp.data(), N);
    egm3_axpy_exact(Va_ss_.begin(), -h, dVa.data(), vap.data(), N);
    run_leg(vdp.data(), vfp.data(), vap.data());

    // --- record(): differences + aggregation, denominator formed ONCE -----
    // Aggregation replicates R's sum(D_ss * .hank3_arr_to_vec(dx)) exactly:
    // walk the STATE order (a fastest, then f, d, e -- the arr_to_vec
    // permutation), round each product to double (the named temporary is
    // volatile so no FMA can fuse it into the add), and accumulate in long
    // double, which is what R's sum() does internally (LDOUBLE). CHI/PHI are
    // differenced as the KERNEL reports them -- never re-derived from da/df,
    // which would drop the Psi2 term (see R/hank-jacobian3.R record()).
    const double den = 2.0 * h;
    long double aD = 0, aF = 0, aA = 0, aC = 0, aCHI = 0, aPHI = 0;
    std::size_t from = 0;
    for (int e = 0; e < ne; ++e)
      for (int d = 0; d < nd; ++d)
        for (int f = 0; f < nf; ++f)
          for (int a = 0; a < na; ++a, ++from) {
            const std::size_t po = pidx3(e, d, f, a, ne, nd, nf);
            const double w = dss[from];
            if (want[0]) { volatile double p = w * ((d_p[po]   - D[po])   / den); aD   += p; }
            if (want[1]) { volatile double p = w * ((f_p[po]   - F[po])   / den); aF   += p; }
            if (want[2]) { volatile double p = w * ((a_p[po]   - A[po])   / den); aA   += p; }
            if (want[3]) { volatile double p = w * ((c_p[po]   - C[po])   / den); aC   += p; }
            if (want[4]) { volatile double p = w * ((chi_p[po] - Chi[po]) / den); aCHI += p; }
            if (want[5]) { volatile double p = w * ((phi_p[po] - Phi[po]) / den); aPHI += p; }
          }
    if (want[0]) cy[0][s] = (double)aD;
    if (want[1]) cy[1][s] = (double)aF;
    if (want[2]) cy[2][s] = (double)aA;
    if (want[3]) cy[3][s] = (double)aC;
    if (want[4]) cy[4][s] = (double)aCHI;
    if (want[5]) cy[5][s] = (double)aPHI;

    // curly-D from the ACTUAL legs through the shared forward-legs kernel;
    // Pi_p == Pi_m == Pi, so `joint` is false and the exact-zero branch runs.
    // The NumericMatrix column is freshly zeroed, as the core requires.
    forward_legs3_core(d_p.data(), f_p.data(), a_p.data(),
                       D.data(), F.data(), A.data(),
                       dg, (int)dg_.size(), fg, (int)fg_.size(),
                       ag, (int)ag_.size(),
                       Pi_.begin(), Pi_.begin(),
                       ne, nd, nf, na, dss, h, &curlyD(0, s));

    // --- next date's direction --------------------------------------------
    for (std::size_t t = 0; t < N; ++t) {
      dVd[t] = (Vd_p[t] - Vdn[t]) / den;
      dVf[t] = (Vf_p[t] - Vfn[t]) / den;
      dVa[t] = (Va_p[t] - Van[t]) / den;
    }
  }
  pool.reset();   // join before anything else touches X

  List curlyY;
  for (int k = 0; k < 6; ++k) if (want[k]) curlyY[nm[k]] = cy[k];
  return List::create(_["curlyY"] = curlyY, _["curlyD"] = curlyD);
}
