// hank_egm.cpp -- compiled EGM household-solver + Young's-method
// forward-operator kernels for the one-asset HANK household block.
//
// Ports R/hank-egm.R::.hank_egm_step()/hank_egm_solve() and
// R/hank-distribution.R::hank_stationary_dist() to RcppArmadillo. The R
// implementations are the reference spec and remain available as a fallback
// path (opt out via getOption("dynhr.hank_backend", "cpp") == "R"); as of
// the "cpp" default flip these kernels are the package's DEFAULT execution
// path. They follow the R spec line-for-line, INCLUDING the tiny = 1e-12
// transient-NaN floors
// (see the long comment on those floors in .hank_egm_step / hank_egm_solve).
// All input validation, feasibility (min_coh) pre-checks, and the default
// Va_init construction stay in the R wrapper (R/hank-egm.R) so error
// messages/tests are unaffected by the backend; the R wrapper computes the
// initial Va and passes it in, so hank_egm_solve_cpp never needs its own
// default-Va logic to match the R path bit-for-bit.
//
// State convention: an n_e x n_a matrix (rows = income states, cols = asset
// gridpoints) exactly as in the R code. hank_stationary_dist_cpp additionally
// flattens to the package's distribution order index(e,a) = (e-1)*n_a + a
// (income OUTER/slow, asset INNER/fast; see R/hank-distribution.R), matching
// hank_forward_operator()'s Lambda layout.
//
// hank_stationary_dist_lambda_cpp is a SEPARATE kernel (not in the original
// prototype) backing hank_stationary_dist()'s cpp path for an ARBITRARY
// caller-supplied sparse Lambda (dgCMatrix): it power-iterates
// d_next = t(Lambda) %*% d directly from Lambda's CSC slots (@p, @i, @x)
// without ever forming the transpose, matching hank_stationary_dist()'s R
// semantics exactly. hank_het_block() itself uses the FASTER fused
// hank_stationary_dist_cpp (a_pol/a_grid/Pi -> d, skipping Lambda for the
// power iteration entirely) since it already has the policy in hand; Lambda
// is still built via hank_forward_operator() for the returned object's
// contract (block$Lambda is read by downstream code/tests regardless of
// backend).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Replicates R findInterval(xq, x) semantics (all.inside=FALSE, default),
// but returns the SAME clamped index used by .hank_interp1 / .hank_lottery:
// 1-indexed idx in [1, n-1] such that x[idx] <= xq < x[idx+1] (R indexing).
// We return it 0-indexed here (idx0 in [0, n-2]) for direct array access:
// x[idx0], x[idx0+1].
static inline int clamped_lower_idx0(const double *x, int n, double xq) {
  // pos = # of x[i] <= xq (1-indexed count) == R findInterval(xq, x)
  int lo = 0, hi = n; // binary search over [0,n): find first pos with x[pos] > xq
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (x[mid] <= xq) lo = mid + 1; else hi = mid;
  }
  int idxR = lo;               // R findInterval value (0..n)
  if (idxR < 1) idxR = 1;
  if (idxR > n - 1) idxR = n - 1;
  return idxR - 1;              // 0-indexed access into x[idxR-1], x[idxR]
}

// One EGM backward step (compiled). C++ port of .hank_egm_step: given
// next-period marginal value Va_p (n_e x n_a), returns the updated marginal
// value and this period's savings/consumption policies at fixed prices
// (r, y). amin is already validated by the R wrapper (amin >= a_grid[1]).
// Returns a list with Va, a, c (each n_e x n_a).
// [[Rcpp::export]]
List hank_egm_step_cpp(NumericMatrix Va_p_, NumericVector a_grid_,
                        NumericVector y_, double r, double beta, double eis,
                        NumericMatrix Pi_, double amin,
                        Nullable<NumericMatrix> coh_extra_ = R_NilValue) {
  int n_e = Va_p_.nrow(), n_a = Va_p_.ncol();
  arma::mat Va_p(Va_p_.begin(), n_e, n_a, false);
  arma::mat Pi(Pi_.begin(), n_e, n_e, false);
  arma::vec a_grid(a_grid_.begin(), n_a, false);
  arma::vec y(y_.begin(), n_e, false);

  // Tier 2: optional additive (e,a) incidence matrix on cash-on-hand, at the
  // FIXED (beginning-of-period) grid only -- see the R doc on .hank_egm_step
  // / hank_het_block's Tr_incidence. Empty/NULL is a strict no-op, so the
  // Tier-1 vector-incidence and no-incidence paths (y already carries any
  // e-indexed transfer) are bit-identical to before this parameter existed.
  bool has_extra = coh_extra_.isNotNull();
  arma::mat coh_extra;
  if (has_extra) {
    NumericMatrix ce(coh_extra_);
    coh_extra = arma::mat(ce.begin(), n_e, n_a, false);
  }

  const double tiny = 1e-12;
  double inv_eis = 1.0 / eis;

  arma::mat Wa = beta * (Pi * Va_p);
  Wa.transform([&](double v) { return v < tiny ? tiny : v; });

  arma::mat c_endog(n_e, n_a);
  for (int e = 0; e < n_e; e++)
    for (int a = 0; a < n_a; a++) {
      double v = std::pow(Wa(e, a), -eis);
      c_endog(e, a) = v < tiny ? tiny : v;
    }

  arma::mat coh_endog(n_e, n_a), coh(n_e, n_a);
  for (int e = 0; e < n_e; e++)
    for (int a = 0; a < n_a; a++) {
      coh_endog(e, a) = c_endog(e, a) + a_grid(a);
      coh(e, a)       = (1.0 + r) * a_grid(a) + y(e) +
                         (has_extra ? coh_extra(e, a) : 0.0);
    }

  arma::mat a_pol(n_e, n_a);
  std::vector<double> xk(n_a);
  for (int e = 0; e < n_e; e++) {
    for (int a = 0; a < n_a; a++) xk[a] = coh_endog(e, a);
    for (int a = 0; a < n_a; a++) {
      double xq = coh(e, a);
      int idx0 = clamped_lower_idx0(xk.data(), n_a, xq);
      double x0 = xk[idx0], x1 = xk[idx0 + 1];
      double y0 = a_grid(idx0), y1 = a_grid(idx0 + 1);
      double w = (xq - x0) / (x1 - x0);
      double val = y0 + w * (y1 - y0);
      a_pol(e, a) = val < amin ? amin : val;
    }
  }

  arma::mat c_pol(n_e, n_a), Va(n_e, n_a);
  for (int e = 0; e < n_e; e++)
    for (int a = 0; a < n_a; a++) {
      double cc = coh(e, a) - a_pol(e, a);
      if (cc < tiny) cc = tiny;
      c_pol(e, a) = cc;
      Va(e, a) = (1.0 + r) * std::pow(cc, -inv_eis);
    }

  return List::create(_["Va"] = Va, _["a"] = a_pol, _["c"] = c_pol);
}

// Full EGM solve loop (compiled). C++ port of hank_egm_solve()'s
// backward-iteration loop: runs ALL EGM iterations inside this single call
// (no per-iteration R boundary crossing). Va_init_ is always supplied by the
// R wrapper (which computes the default guess itself, so the default and
// user-supplied cases are handled identically here). tol/maxit: convergence
// tolerance (max abs change in a) and iteration cap. Returns a list with Va,
// a, c, iterations, converged.
// [[Rcpp::export]]
List hank_egm_solve_cpp(NumericVector a_grid_, NumericVector y_, double r,
                         double beta, double eis, NumericMatrix Pi_,
                         double amin, double tol, int maxit,
                         Nullable<NumericMatrix> Va_init_ = R_NilValue,
                         Nullable<NumericMatrix> coh_extra_ = R_NilValue) {
  int n_e = y_.size(), n_a = a_grid_.size();
  arma::mat Pi(Pi_.begin(), n_e, n_e, false);
  arma::vec a_grid(a_grid_.begin(), n_a, false);
  arma::vec y(y_.begin(), n_e, false);
  const double tiny = 1e-12;
  double inv_eis = 1.0 / eis;

  // Tier 2: see hank_egm_step_cpp -- empty/NULL is a strict no-op.
  bool has_extra = coh_extra_.isNotNull();
  arma::mat coh_extra;
  if (has_extra) {
    NumericMatrix ce(coh_extra_);
    coh_extra = arma::mat(ce.begin(), n_e, n_a, false);
  }

  arma::mat Va(n_e, n_a);
  if (Va_init_.isNotNull()) {
    NumericMatrix vi(Va_init_);
    Va = arma::mat(vi.begin(), n_e, n_a);
  } else {
    for (int e = 0; e < n_e; e++)
      for (int a = 0; a < n_a; a++) {
        double coh0 = (1.0 + r) * a_grid(a) + y(e);
        double g = 0.1 * coh0;
        if (g < tiny) g = tiny;
        Va(e, a) = (1.0 + r) * std::pow(g, -inv_eis);
      }
  }

  arma::mat a_old(n_e, n_a, arma::fill::value(-arma::datum::inf));
  arma::mat a_pol(n_e, n_a), c_pol(n_e, n_a);
  bool converged = false;
  int it = 0;

  std::vector<double> xk(n_a);
  arma::mat Wa(n_e, n_a), c_endog(n_e, n_a), coh_endog(n_e, n_a), coh(n_e, n_a);

  for (it = 1; it <= maxit; it++) {
    Wa = beta * (Pi * Va);
    Wa.transform([&](double v) { return v < tiny ? tiny : v; });

    for (int e = 0; e < n_e; e++)
      for (int a = 0; a < n_a; a++) {
        double v = std::pow(Wa(e, a), -eis);
        c_endog(e, a) = v < tiny ? tiny : v;
        coh_endog(e, a) = c_endog(e, a) + a_grid(a);
        coh(e, a) = (1.0 + r) * a_grid(a) + y(e) +
                    (has_extra ? coh_extra(e, a) : 0.0);
      }

    for (int e = 0; e < n_e; e++) {
      for (int a = 0; a < n_a; a++) xk[a] = coh_endog(e, a);
      for (int a = 0; a < n_a; a++) {
        double xq = coh(e, a);
        int idx0 = clamped_lower_idx0(xk.data(), n_a, xq);
        double x0 = xk[idx0], x1 = xk[idx0 + 1];
        double y0 = a_grid(idx0), y1 = a_grid(idx0 + 1);
        double w = (xq - x0) / (x1 - x0);
        double val = y0 + w * (y1 - y0);
        a_pol(e, a) = val < amin ? amin : val;
      }
    }

    double maxdiff = 0.0;
    for (int e = 0; e < n_e; e++)
      for (int a = 0; a < n_a; a++) {
        double cc = coh(e, a) - a_pol(e, a);
        if (cc < tiny) cc = tiny;
        c_pol(e, a) = cc;
        Va(e, a) = (1.0 + r) * std::pow(cc, -inv_eis);
        double d = std::fabs(a_pol(e, a) - a_old(e, a));
        if (d > maxdiff) maxdiff = d;
      }

    a_old = a_pol;
    if (maxdiff < tol) { converged = true; break; }
  }
  if (it > maxit) it = maxit;

  return List::create(_["Va"] = Va, _["a"] = a_pol, _["c"] = c_pol,
                       _["iterations"] = it, _["converged"] = converged);
}

// Fused Young's-method stationary distribution from a savings policy
// (compiled). Direct scatter power iteration: builds the (i, p) lottery once
// from (a_pol, a_grid) and scatters mass through Pi each iteration WITHOUT
// ever materializing the sparse Lambda matrix. Fuses hank_forward_operator()
// + hank_stationary_dist() for the hank_het_block() hot path, where Lambda
// is not needed for the power iteration itself (it is still built
// separately, via the R hank_forward_operator, when the caller's object
// needs it). tol/maxit: convergence tolerance and iteration cap. Returns a
// list with d (length n_e*n_a, distribution order (e-1)*n_a+a), iterations,
// converged.
// [[Rcpp::export]]
List hank_stationary_dist_cpp(NumericMatrix a_pol_, NumericVector a_grid_,
                               NumericMatrix Pi_, double tol, int maxit) {
  int n_e = a_pol_.nrow(), n_a = a_pol_.ncol();
  arma::mat a_pol(a_pol_.begin(), n_e, n_a, false);
  arma::vec a_grid(a_grid_.begin(), n_a, false);
  arma::mat Pi(Pi_.begin(), n_e, n_e, false);

  IntegerMatrix lot_i(n_e, n_a);
  NumericMatrix lot_p(n_e, n_a);
  for (int e = 0; e < n_e; e++)
    for (int a = 0; a < n_a; a++) {
      double ap = a_pol(e, a);
      int idx0 = clamped_lower_idx0(a_grid.memptr(), n_a, ap); // 0-indexed lower
      double aL = a_grid(idx0), aU = a_grid(idx0 + 1);
      double p = (aU - ap) / (aU - aL);
      if (p < 0) p = 0; if (p > 1) p = 1;
      lot_i(e, a) = idx0; // 0-indexed lower bracket
      lot_p(e, a) = p;
    }

  int n = n_e * n_a;
  std::vector<double> d(n, 1.0 / n), d_new(n);
  bool converged = false;
  int it = 0;

  for (it = 1; it <= maxit; it++) {
    std::fill(d_new.begin(), d_new.end(), 0.0);
    for (int e = 0; e < n_e; e++) {
      int base_from = e * n_a;
      for (int a = 0; a < n_a; a++) {
        double mass = d[base_from + a];
        if (mass == 0.0) continue;
        int ii = lot_i(e, a);
        double pp = lot_p(e, a);
        for (int ep = 0; ep < n_e; ep++) {
          double pr = Pi(e, ep);
          if (pr == 0.0) continue;
          int base_to = ep * n_a;
          d_new[base_to + ii]     += pr * pp * mass;
          d_new[base_to + ii + 1] += pr * (1.0 - pp) * mass;
        }
      }
    }
    double maxdiff = 0.0;
    for (int k = 0; k < n; k++) {
      double diff = std::fabs(d_new[k] - d[k]);
      if (diff > maxdiff) maxdiff = diff;
    }
    d.swap(d_new);
    if (maxdiff < tol) { converged = true; break; }
  }
  if (it > maxit) it = maxit;

  double s = 0.0;
  for (int k = 0; k < n; k++) s += d[k];
  NumericVector dout(n);
  for (int k = 0; k < n; k++) dout[k] = d[k] / s;

  return List::create(_["d"] = dout, _["iterations"] = it, _["converged"] = converged);
}

// Stationary distribution power iteration from an explicit sparse Lambda
// (compiled). Backs hank_stationary_dist(Lambda, ..., backend = "cpp") for
// an ARBITRARY caller-supplied dgCMatrix Lambda (not necessarily built by
// hank_forward_operator). Iterates d_next = t(Lambda) %*% d directly from
// Lambda's CSC slots (@p, @i, @x) without ever forming the transpose: column
// j of Lambda (entries p[j]..p[j+1]-1) IS row j of t(Lambda), so
// d_next[j] = sum_k x[k]*d[i[k]] over that column's entries. Matches
// hank_stationary_dist()'s R semantics (same convergence rule, same final
// sum-normalization). p/i/x: Lambda's CSC slots (@p length n+1, @i and @x
// length nnz). n: dimension (nrow(Lambda) == ncol(Lambda)). d0: initial
// distribution (length n), already normalized to sum to 1 by the R wrapper.
// tol/maxit: convergence tolerance and iteration cap. Returns a list with d,
// iterations, converged.
// [[Rcpp::export]]
List hank_stationary_dist_lambda_cpp(IntegerVector p, IntegerVector i,
                                      NumericVector x, int n,
                                      NumericVector d0, double tol, int maxit) {
  // Defensive contracts (adversarial review 2026-07-13, P1): the R wrapper
  // validates first, but a compiled kernel must not rely on that forever.
  // A short d0 previously read past the end of the std::vector (UB); the
  // NaN-swallowing `diff > maxdiff` comparison previously reported
  // converged = TRUE on an all-NaN iterate.
  if ((int)d0.size() != n)
    stop("hank_stationary_dist_lambda_cpp: d0 has length %d, expected n = %d",
         (int)d0.size(), n);
  if ((int)p.size() != n + 1)
    stop("hank_stationary_dist_lambda_cpp: p has length %d, expected n+1 = %d",
         (int)p.size(), n + 1);
  if (i.size() != x.size())
    stop("hank_stationary_dist_lambda_cpp: i and x lengths differ (%d vs %d)",
         (int)i.size(), (int)x.size());
  if (p[n] > (int)i.size() || p[0] != 0)
    stop("hank_stationary_dist_lambda_cpp: invalid CSC pointer slot p");
  for (int k = 0; k < (int)i.size(); k++)
    if (i[k] < 0 || i[k] >= n)
      stop("hank_stationary_dist_lambda_cpp: row index i[%d] = %d out of "
           "range [0, %d)", k, (int)i[k], n);

  std::vector<double> d(d0.begin(), d0.end()), d_new(n);
  bool converged = false;
  int it = 0;

  for (it = 1; it <= maxit; it++) {
    std::fill(d_new.begin(), d_new.end(), 0.0);
    for (int j = 0; j < n; j++) {
      double acc = 0.0;
      for (int k = p[j]; k < p[j + 1]; k++) acc += x[k] * d[i[k]];
      d_new[j] = acc;
    }
    double maxdiff = 0.0;
    bool bad = false;
    for (int k = 0; k < n; k++) {
      double diff = std::fabs(d_new[k] - d[k]);
      if (std::isnan(diff)) { bad = true; break; }
      if (diff > maxdiff) maxdiff = diff;
    }
    d.swap(d_new);
    if (bad) { converged = false; break; }  // NaN iterate: fail loudly, fast
    if (maxdiff < tol) { converged = true; break; }
  }
  if (it > maxit) it = maxit;

  double s = 0.0;
  for (int k = 0; k < n; k++) s += d[k];
  NumericVector dout(n);
  for (int k = 0; k < n; k++) dout[k] = d[k] / s;

  return List::create(_["d"] = dout, _["iterations"] = it, _["converged"] = converged);
}
