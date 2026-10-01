// order2_deriv.cpp -- Kronecker-structured solve for the order-2 solution
// derivative and its adjoint.
//
// The forward derivative d(ghxx)/dtheta and the reverse-mode cotangent both
// solve the generalised Sylvester equation
//
//     A X + F X (H (x) H) = RHS,        X, RHS : n x ns^2,
//
// (forward: A = A_L, F = f_plus, H = hx; transposed / adjoint: A = A_L',
// F = f_plus', H = hx').  Written as the vec-system K_xx vec(X) = vec(RHS)
// with K_xx = I (x) A + (H' (x) H') (x) F this is a dense (n ns^2)^2 matrix
// -- 16000 x 16000 (2 GB) for a 40-equation / 20-state model -- so the dense
// route neither builds nor factors it here.
//
// Method (no inverse of K_xx, no ns^2 x ns^2 matrix):
//   complex generalised Schur  A = Q^H S Z^H,  F = Q^H P Z^H   (S, P upper)
//   complex Schur              H = U T U^H                     (T upper)
//   W = Z^H X (U (x) U)  solves  S W + P W (T (x) T) = Q RHS (U (x) U),
//   and T (x) T is upper triangular, so W is found column by column with one
//   n x n triangular solve each.  The right-multiplication by a Kronecker
//   square is done row by row as  vec(M V) = reshape( V' M_mat V ).
// Iterative refinement against the REAL residual A X + F X (H (x) H) - RHS
// drives the result to the conditioning floor; the achieved relative residual
// is returned so the caller can route to the dense reference if it is not
// small (near-resonant pencil, failed Schur).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::cx_mat;
using arma::cx_vec;
using arma::mat;
using arma::uword;

namespace {

// Y = X (A (x) A) for the row-wise Kronecker square, X: r x (p*p), A: p x p
// (complex). Row vector index is a*p + k (a = block of the left factor), so a
// row reshaped column-major to p x p is M with M(k, a); the product is
// A^T M A (transposed, not conjugated).
cx_mat kron2_right(const cx_mat& X, const cx_mat& A) {
  const uword r = X.n_rows;
  const uword p = A.n_rows;
  cx_mat Y(r, p * p);
  cx_mat M(p, p);
  const cx_mat At = A.st();
  for (uword i = 0; i < r; ++i) {
    for (uword a = 0; a < p; ++a)
      for (uword k = 0; k < p; ++k) M(k, a) = X(i, a * p + k);
    cx_mat R = At * M * A;
    for (uword a = 0; a < p; ++a)
      for (uword k = 0; k < p; ++k) Y(i, a * p + k) = R(k, a);
  }
  return Y;
}

// Real version: Y = X (A (x) A).
mat kron2_right_real(const mat& X, const mat& A) {
  const uword r = X.n_rows;
  const uword p = A.n_rows;
  mat Y(r, p * p);
  mat M(p, p);
  const mat At = A.t();
  for (uword i = 0; i < r; ++i) {
    for (uword a = 0; a < p; ++a)
      for (uword k = 0; k < p; ++k) M(k, a) = X(i, a * p + k);
    mat R = At * M * A;
    for (uword a = 0; a < p; ++a)
      for (uword k = 0; k < p; ++k) Y(i, a * p + k) = R(k, a);
  }
  return Y;
}

struct Kron2Solver {
  uword n = 0, ns = 0, m = 0;
  mat A, F, H;
  cx_mat Q, Z, S, P, U, T;
  cx_mat Uh;  // U^H
  bool ok = false;
  mutable bool singular = false;   // an exactly / numerically singular column pencil

  Kron2Solver(const mat& A_, const mat& F_, const mat& H_)
      : n(A_.n_rows), ns(H_.n_rows), m(H_.n_rows * H_.n_rows),
        A(A_), F(F_), H(H_) {
    cx_mat Ac = arma::conv_to<cx_mat>::from(A);
    cx_mat Fc = arma::conv_to<cx_mat>::from(F);
    cx_mat Hc = arma::conv_to<cx_mat>::from(H);
    cx_mat AA, BB, Qz, Zz;
    if (!arma::qz(AA, BB, Qz, Zz, Ac, Fc)) return;
    // Armadillo: A = Qz^H AA Zz^H  =>  Qz A Zz = AA.
    S = AA; P = BB; Q = Qz; Z = Zz;
    if (!arma::schur(U, T, Hc)) return;
    Uh = U.t();
    ok = true;
  }

  // One Schur-basis solve of A X + F X (H (x) H) = B (B real n x m).
  mat solve_once(const mat& B) const {
    cx_mat Bc = arma::conv_to<cx_mat>::from(B);
    cx_mat G = kron2_right(Q * Bc, U);
    cx_mat W(n, m, arma::fill::zeros);
    cx_vec s(n), rhs(n);
    for (uword b = 0; b < ns; ++b) {
      for (uword l = 0; l < ns; ++l) {
        const uword j = b * ns + l;
        s.zeros();
        // strictly-upper coupling: i = a*ns + k with a <= b, k <= l, i != j.
        for (uword a = 0; a <= b; ++a) {
          const std::complex<double> tab = T(a, b);
          for (uword k = 0; k <= l; ++k) {
            const uword i = a * ns + k;
            if (i == j) continue;
            s += W.col(i) * (tab * T(k, l));
          }
        }
        rhs = G.col(j) - P * s;
        const std::complex<double> tau = T(b, b) * T(l, l);
        cx_mat M = S + tau * P;
        // The column pencil S + tau P is upper triangular; a (numerically) zero
        // diagonal entry means the Kronecker system is singular (resonant
        // pencil). Flag it and let the caller route to the dense reference.
        bool bad = false;
        for (uword d = 0; d < n; ++d) {
          const double mag = std::abs(S(d, d)) + std::abs(tau) * std::abs(P(d, d));
          if (!(std::abs(M(d, d)) > 1e-14 * mag)) { bad = true; break; }
        }
        if (bad) {
          singular = true;
          W.col(j).fill(std::complex<double>(R_NaN, 0.0));
          continue;
        }
        W.col(j) = arma::solve(arma::trimatu(M), rhs);
      }
    }
    cx_mat Xc = Z * kron2_right(W, Uh);
    return arma::real(Xc);
  }

  mat apply(const mat& X) const {
    return A * X + F * kron2_right_real(X, H);
  }
};

}  // namespace

// [[Rcpp::export]]
Rcpp::List o2_sylvester_kron2_cpp(const arma::mat& A, const arma::mat& F,
                                  const arma::mat& H, const arma::mat& RHS,
                                  int max_refine = 4) {
  const uword n = A.n_rows;
  const uword ns = H.n_rows;
  if (A.n_cols != n || F.n_rows != n || F.n_cols != n ||
      H.n_cols != ns || RHS.n_rows != n || RHS.n_cols != ns * ns)
    Rcpp::stop("o2_sylvester_kron2_cpp: dimension mismatch");
  if (ns == 0)
    return Rcpp::List::create(Rcpp::Named("X") = arma::mat(n, 0),
                              Rcpp::Named("rel_resid") = 0.0,
                              Rcpp::Named("ok") = true);
  Kron2Solver sol(A, F, H);
  if (!sol.ok || !RHS.is_finite())
    return Rcpp::List::create(Rcpp::Named("X") = arma::mat(n, ns * ns, arma::fill::zeros),
                              Rcpp::Named("rel_resid") = R_PosInf,
                              Rcpp::Named("ok") = false);
  const double scale = std::max(1.0, arma::abs(RHS).max());
  mat X = sol.solve_once(RHS);
  double rel = R_PosInf;
  if (X.is_finite() && !sol.singular) {
    for (int it = 0; it <= max_refine; ++it) {
      mat R = RHS - sol.apply(X);
      rel = arma::abs(R).max() / scale;
      if (!std::isfinite(rel) || rel < 1e-15 || it == max_refine) break;
      mat dX = sol.solve_once(R);
      if (!dX.is_finite()) { rel = R_PosInf; break; }
      X += dX;
    }
  }
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("rel_resid") = rel,
                            Rcpp::Named("ok") = std::isfinite(rel));
}

// Reverse-mode contraction of the third-cumulant tensor-Lyapunov operator with
// respect to hx: the three mode contractions
//   bar[a,i] += sum_{b,c,j,k} Ma[a,b,c] hx[b,j] hx[c,k] T3[i,j,k]
//   bar[b,j] += sum_{a,c,i,k} Ma[a,b,c] hx[a,i] hx[c,k] T3[i,j,k]
//   bar[c,k] += sum_{a,b,i,j} Ma[a,b,c] hx[a,i] hx[b,j] T3[i,j,k]
// with T3[i,j,k] = c3[i, k*ns + j] and Ma[a,b,c] = M[a, c*ns + b]. Each mode
// reduces to one ns x ns sandwich  hx S hx'  per slice, so the cost is
// O(ns^4) rather than the O(ns^6) of the scalar loops.
// [[Rcpp::export]]
arma::mat o2_third_cross_bar_hx_cpp(const arma::mat& hx, const arma::mat& c3,
                                    const arma::mat& M) {
  const uword ns = hx.n_rows;
  if (hx.n_cols != ns || c3.n_rows != ns || c3.n_cols != ns * ns ||
      M.n_rows != ns || M.n_cols != ns * ns)
    Rcpp::stop("o2_third_cross_bar_hx_cpp: dimension mismatch");
  mat bar(ns, ns, arma::fill::zeros);
  if (ns == 0) return bar;
  const mat hxt = hx.t();
  mat S(ns, ns);
  // mode 1: slice over i, S(j,k) = T3[i,j,k]
  for (uword i = 0; i < ns; ++i) {
    for (uword k = 0; k < ns; ++k)
      for (uword j = 0; j < ns; ++j) S(j, k) = c3(i, k * ns + j);
    const mat P = hx * S * hxt;
    for (uword a = 0; a < ns; ++a) {
      double acc = 0.0;
      for (uword c = 0; c < ns; ++c)
        for (uword b = 0; b < ns; ++b) acc += M(a, c * ns + b) * P(b, c);
      bar(a, i) += acc;
    }
  }
  // mode 2: slice over j, S(i,k) = T3[i,j,k]
  for (uword j = 0; j < ns; ++j) {
    for (uword k = 0; k < ns; ++k)
      for (uword i = 0; i < ns; ++i) S(i, k) = c3(i, k * ns + j);
    const mat D = hx * S * hxt;
    for (uword b = 0; b < ns; ++b) {
      double acc = 0.0;
      for (uword c = 0; c < ns; ++c)
        for (uword a = 0; a < ns; ++a) acc += M(a, c * ns + b) * D(a, c);
      bar(b, j) += acc;
    }
  }
  // mode 3: slice over k, S(i,j) = T3[i,j,k]
  for (uword k = 0; k < ns; ++k) {
    for (uword j = 0; j < ns; ++j)
      for (uword i = 0; i < ns; ++i) S(i, j) = c3(i, k * ns + j);
    const mat E = hx * S * hxt;
    for (uword c = 0; c < ns; ++c) {
      double acc = 0.0;
      for (uword b = 0; b < ns; ++b)
        for (uword a = 0; a < ns; ++a) acc += M(a, c * ns + b) * E(a, b);
      bar(c, k) += acc;
    }
  }
  return bar;
}
