// hank_push2.cpp -- compiled scatter kernel for the TWO-ASSET matrix-free
// forward push, R/hank-distribution2.R::.hank_forward_push2().
//
// SCOPE. This kernel covers exactly the LOTTERIES + BILINEAR SCATTER half of
// that function and returns the raw accumulator `acc`, indexed by the SOURCE
// income block. The final income-mixing crossprod(Pi, .) deliberately stays in
// R: it is one small BLAS matmul against a scatter over n_e*n_b*n_a cells, so
// reimplementing it here could only differ from BLAS at round-off for no
// measurable gain. Keeping the R/C++ boundary there is what makes the compiled
// path bit-identical to the R reference rather than merely close.
//
// WHY. Measured on the installed -O2 build (n_e = 3, n_b = n_a = 50,
// T_h = 50, 3 inputs x 3 outputs), .hank_forward_push2() was ~0.600 s of a
// 0.794 s two-asset Jacobian build -- 76% -- called 2*T_h*n_inputs = 300 times
// at ~2.0 ms each. Everything else on that path was already compiled. The R
// cost is NOT the arithmetic: it is materializing `key` and `w` at 4*n_cell
// elements each and handing them to rowsum(). This kernel accumulates directly
// into `acc` in one pass per corner and allocates nothing of that size.
//
// BIT-IDENTITY (verified with identical(), not assumed). rowsum(w, key,
// reorder = FALSE) sums each key's contributions in INCREASING ELEMENT ORDER
// of the concatenated vectors, and R builds those as
//   key = c(base, base + n_a, base + 1, base + n_a + 1)
//   w   = c(Dv*pb*pa, Dv*(1-pb)*pa, Dv*pb*(1-pa), Dv*(1-pb)*(1-pa)).
// So a destination cell receives ALL corner-0 contributions (in cell order),
// then all corner-1, and so on. This kernel therefore runs FOUR separate
// corner-major passes over the cells; folding them into one pass with four
// adds per cell would interleave the additions and lose the last ULP. The
// per-weight multiplication order (Dv*pb)*pa is matched literally, and each
// product is computed in its own statement so that -ffp-contract cannot fuse
// the final multiply-add into an FMA (contraction is statement-local).
//
// NOT THREADED, deliberately. A single call is ~2 ms even at n = 50, which is
// the regime where per-step worker-pool spawning was already measured to be a
// net loss on this codebase (0.9.0.0038 notes). A scatter is also a reduction
// over shared destinations, so a threaded version would need either per-thread
// accumulators (destroying the summation order above, hence bit-identity) or
// atomics (slower than the serial loop at this size). Serial-but-compiled is
// the design.
//
// CELL ORDER is the package's two-asset convention (R/hank-distribution2.R):
//   cell(e, j, k) = e*n_b*n_a + j*n_a + k   (0-based; income slowest, illiquid
// a fastest), while the POLICY arrays are R's column-major dim c(n_e,n_b,n_a),
// offset e + n_e*j + n_e*n_b*k. The R reference reaches the same pairing via
// aperm() gymnastics on the (n_e*n_a) x n_b and (n_e*n_b) x n_a unfoldings;
// here it is the index arithmetic above. Getting that transposed is the most
// likely bug in this file AND IT STILL CONSERVES MASS, so the gate is the
// sparse-operator cross-oracle in test-hank-forward-push.R, never mass alone.

#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

namespace {

// findInterval(x, grid) with the lottery's clamp, matching
// R/hank-distribution.R::.hank_lottery() exactly: the raw index is the number
// of grid nodes <= x (so 0 below the grid, n at/above the top), clamped to
// [1, n-1] so that i and i+1 are both valid; the returned value is 0-BASED.
// The weight is (U - x) / (U - L), clamped to [0, 1] -- which is what makes an
// off-grid policy put its whole mass on one bracketing node and hand the other
// an EXACT zero.
inline void lottery1(const double x, const double* g, const int n,
                     int& i0, double& p) {
  int i = static_cast<int>(std::upper_bound(g, g + n, x) - g);
  if (i < 1) i = 1;
  if (i > n - 1) i = n - 1;
  const double L = g[i - 1];
  const double U = g[i];
  double w = (U - x) / (U - L);
  if (w < 0) w = 0;
  if (w > 1) w = 1;
  i0 = i - 1;
  p = w;
}

}  // namespace

//' Compiled lotteries + bilinear scatter for the two-asset forward push
//'
//' Kernel behind \code{.hank_forward_push2(backend = "cpp")}.  Returns the
//' scatter accumulator BEFORE the income mixing, i.e. the vector \code{acc}
//' indexed \code{(e-1)*n_b*n_a + j} with \code{e} the SOURCE income state; the
//' caller applies \code{crossprod(Pi, .)}.
//'
//' @param b_pol,a_pol Numeric length-\code{n_e*n_b*n_a} vectors: the liquid and
//'   illiquid policies in R's column-major \code{c(n_e, n_b, n_a)} layout.
//' @param b_grid,a_grid Numeric increasing grids of length \code{n_b} /
//'   \code{n_a} (both must have length >= 2; the R wrapper enforces this).
//' @param Dv Numeric length-\code{n_e*n_b*n_a} distribution in the package's
//'   two-asset CELL order.
//' @param n_e,n_b,n_a Integer dimensions.
//' @return Numeric length-\code{n_e*n_b*n_a} accumulator in cell order.
//' @keywords internal
// [[Rcpp::export]]
NumericVector hank_forward_push2_scatter_cpp(NumericVector b_pol,
                                             NumericVector a_pol,
                                             NumericVector b_grid,
                                             NumericVector a_grid,
                                             NumericVector Dv,
                                             int n_e, int n_b, int n_a) {
  const int n_ba = n_b * n_a;
  const int n_cell = n_e * n_ba;
  if (b_pol.size() != n_cell || a_pol.size() != n_cell ||
      Dv.size() != n_cell || b_grid.size() != n_b || a_grid.size() != n_a)
    stop("hank_forward_push2_scatter_cpp(): dimension mismatch.");
  if (n_b < 2 || n_a < 2)
    stop("hank_forward_push2_scatter_cpp(): both grids need >= 2 nodes.");

  const double* bp = b_pol.begin();
  const double* ap = a_pol.begin();
  const double* bg = b_grid.begin();
  const double* ag = a_grid.begin();
  const double* dv = Dv.begin();

  // Pass 0: per-cell lotteries on both axes, plus the lower-corner
  // destination. `base` already carries the source income block offset, so the
  // scatter stays within the source state (the mixing happens afterwards).
  std::vector<int> base(n_cell);
  std::vector<double> pb(n_cell), pa(n_cell);
  for (int e = 0; e < n_e; ++e) {
    const int eoff = e * n_ba;
    for (int j = 0; j < n_b; ++j) {
      for (int k = 0; k < n_a; ++k) {
        const int c = eoff + j * n_a + k;
        const int src = e + n_e * j + n_e * n_b * k;
        int ib0, ia0;
        double wb, wa;
        lottery1(bp[src], bg, n_b, ib0, wb);
        lottery1(ap[src], ag, n_a, ia0, wa);
        base[c] = eoff + ib0 * n_a + ia0;
        pb[c] = wb;
        pa[c] = wa;
      }
    }
  }

  NumericVector out(n_cell);
  double* acc = out.begin();

  // Four CORNER-MAJOR passes, in the order R concatenates them. Each product
  // is formed in its own statement before the accumulate, so the compiler
  // cannot contract `acc += x*y` into an FMA and drift from R's `*` then `+`.
  for (int c = 0; c < n_cell; ++c) {          // (b_lo, a_lo)
    const double v = dv[c] * pb[c] * pa[c];
    acc[base[c]] += v;
  }
  for (int c = 0; c < n_cell; ++c) {          // (b_hi, a_lo): base + n_a
    const double v = dv[c] * (1 - pb[c]) * pa[c];
    acc[base[c] + n_a] += v;
  }
  for (int c = 0; c < n_cell; ++c) {          // (b_lo, a_hi): base + 1
    const double v = dv[c] * pb[c] * (1 - pa[c]);
    acc[base[c] + 1] += v;
  }
  for (int c = 0; c < n_cell; ++c) {          // (b_hi, a_hi): base + n_a + 1
    const double v = dv[c] * (1 - pb[c]) * (1 - pa[c]);
    acc[base[c] + n_a + 1] += v;
  }

  return out;
}
