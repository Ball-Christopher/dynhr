// kf_guard.h -- cheap per-step risk indicators for the simple covariance
// update P' = T P T' + R Q R' - K F K' (see kalman_ss.cpp / kalman_struct.cpp).
//
// The simple form subtracts K F K' from T P T' + R Q R'; the Joseph form is
// a sum of positive semidefinite terms and cannot lose a digit that way. The
// simple form is therefore exact to round-off only while the subtraction does
// not cancel catastrophically, and these indicators, all O(n) or O(p) per
// step on quantities the step already holds (amp adds one diagonal of
// Z B Z'), flag the calls where it might:
//   piv    min over steps of (min pivot / max pivot)^2 of the Cholesky factor
//          of F (an rcond proxy; scale dependent);
//   r2     min over steps and observables of 1 / (F_ii (F^{-1})_ii), which is
//          1 - R^2 of observable i on the others: a scale-free conditioning
//          measure of F;
//   ret    min over steps and states of P'_ii / (T P T' + R Q R')_ii, the
//          fraction of the prior variance that survives the subtraction;
//   amp    max over steps and observables of (Z B Z')_ii / F'_ii, where B =
//          T P T' + R Q R' is the matrix the subtraction acts on and F' the
//          NEXT step's innovation variance: the factor by which the absolute
//          round-off of the subtraction (about eps * B) is magnified in the
//          next F, relative to F itself;
//   bad    a non-finite diagonal of P'.
// A threshold < 0 only records its indicator; a threshold >= 0 flags the call
// as soon as the indicator crosses it (piv, r2, ret: below; amp: above).
#ifndef DYNHR_KF_GUARD_H
#define DYNHR_KF_GUARD_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

struct KfGuard {
  double thr_piv, thr_r2, thr_ret, thr_amp;
  bool on;
  double piv = 1.0, r2 = 1.0, ret = 1.0, amp = 0.0;
  bool bad = false;
  bool tripped = false;
  arma::vec pending;                 // (Z B Z')_ii from the previous step

  KfGuard(double tp, double tr2, double tret, double tamp)
      : thr_piv(tp), thr_r2(tr2), thr_ret(tret), thr_amp(tamp),
        on(tr2 >= 0.0 || tret >= 0.0 || tp >= 0.0 || tamp >= 0.0) {}

  // Rc: upper Cholesky factor of F; Ft, Fi: F and its inverse; O: the
  // observed rows when F covers a subset (nullptr = all rows in order).
  void note_F(const arma::mat& Rc, const arma::mat& Ft, const arma::mat& Fi,
              const arma::uword* O = nullptr) {
    const arma::uword p = Ft.n_rows;
    double mn = 1e300, mx = 0.0;
    double r2m = 1.0;
    for (arma::uword i = 0; i < p; ++i) {
      const double d = Rc(i, i);
      mn = std::min(mn, d);
      mx = std::max(mx, d);
      const double q = 1.0 / (Ft(i, i) * Fi(i, i));
      if (q < r2m) r2m = q;
    }
    if (mx > 0.0) {
      const double r = mn / mx;
      piv = std::min(piv, r * r);
    }
    r2 = std::min(r2, r2m);
    if (pending.n_elem > 0) {
      for (arma::uword a = 0; a < p; ++a) {
        const arma::uword i = O ? O[a] : a;
        const double r = pending[i] / Ft(a, a);
        if (r > amp) amp = r;
      }
      pending.reset();
    }
    if (on && ((thr_piv >= 0.0 && piv < thr_piv) ||
               (thr_r2 >= 0.0 && r2 < thr_r2) ||
               (thr_amp >= 0.0 && amp > thr_amp))) tripped = true;
  }

  // base: T P T' + R Q R'; nd: the updated P' (only diagonals are read).
  void note_P(const arma::mat& base, const arma::mat& nd) {
    const arma::uword n = nd.n_rows;
    for (arma::uword i = 0; i < n; ++i) {
      const double d = nd(i, i);
      if (!std::isfinite(d)) { bad = true; continue; }
      const double b = base(i, i);
      if (b > 0.0) {
        const double r = d / b;
        if (r < ret) ret = r;
      }
    }
    if (on && (bad || (thr_ret >= 0.0 && ret < thr_ret))) tripped = true;
  }
};

#endif
