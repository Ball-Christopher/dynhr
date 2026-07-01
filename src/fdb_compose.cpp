// fdb_compose.cpp -- C++ inner loop for the folded Faa-di-Bruno composition.
//
// Provides fdb_compose_folded_cpp(), the compiled backend for
// .fdb_compose_folded() in R/solve-perturbation-symbolic.R.  The pure-R
// fallback (.fdb_compose_folded_R) is kept and is value-identical up to
// floating-point summation order; see test-fdb-rcpp-parity.R.
//
// The data-independent fold plan (.fdb_fold_plan, cached in R) is built in R and
// passed in.  This routine evaluates, for each canonical output column, the
// symmetric tensor value
//
//   D[a, c] = sum_shapes sum_{triplet (a,b,val)} val
//               sum_{partition instances} sum_{bp in distinct perms(b)}
//                 prod_i  Hlist[[csizes_i]][ bp_i , inst_i[c] ]
//
// Permutations of the (canonical, sorted-ascending) b-vector are generated with
// std::next_permutation, which yields exactly the distinct multiset permutations.

#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fdb_compose_folded_cpp(int Nc, int n_out,
                                     List shapes,   // plan$shapes
                                     List Glist,    // Glist[[m]] = list(eq, cols, val)
                                     List Hlist) {  // Hlist[[s]] = NumericMatrix
  NumericMatrix D(n_out, Nc);
  std::vector<double> contrib(Nc);
  std::vector<double> pbuf(Nc);

  const int n_shapes = shapes.size();
  for (int sh = 0; sh < n_shapes; ++sh) {
    List shape = shapes[sh];
    const int m = as<int>(shape["m"]);

    // Outer map derivative of order m (number of blocks).
    if (Rf_isNull(Glist[m - 1])) continue;
    List G = Glist[m - 1];
    NumericVector val = G["val"];
    const int ntrip = val.size();
    if (ntrip == 0) continue;
    IntegerVector eq = G["eq"];
    IntegerMatrix cols = G["cols"];

    IntegerVector csizes = shape["csizes"];      // length m, block sizes (desc)
    List insts = shape["insts"];                 // list of instances
    const int n_inst = insts.size();

    // Per-slot inner-map matrix (Hlist[[csizes_i]]) and its row count.
    std::vector<NumericMatrix> Hsz(m);
    std::vector<const double*> Hdat(m);
    std::vector<int> Hnrow(m);
    for (int i = 0; i < m; ++i) {
      NumericMatrix Hi = Hlist[csizes[i] - 1];
      Hsz[i] = Hi;
      Hdat[i] = Hi.begin();
      Hnrow[i] = Hi.nrow();
    }

    // Pre-extract instance index vectors (1-based flat columns) as raw pointers:
    // inst_ptr[inst][i] -> length-Nc int* into Hsz[i]'s columns.
    std::vector< std::vector<const int*> > inst_ptr(n_inst,
                                                    std::vector<const int*>(m));
    std::vector<IntegerVector> inst_keep;          // keep vectors alive
    inst_keep.reserve(n_inst * m);
    for (int q = 0; q < n_inst; ++q) {
      List inst = insts[q];
      for (int i = 0; i < m; ++i) {
        inst_keep.push_back(as<IntegerVector>(inst[i]));
        inst_ptr[q][i] = inst_keep.back().begin();
      }
    }

    for (int t = 0; t < ntrip; ++t) {
      const double vt = val[t];
      if (vt == 0.0) continue;
      const int a = eq[t] - 1;                     // 0-based output row

      // Canonical b-vector (sorted ascending) -> distinct perms via next_permutation.
      std::vector<int> b(m);
      for (int i = 0; i < m; ++i) b[i] = cols(t, i);
      std::sort(b.begin(), b.end());

      std::fill(contrib.begin(), contrib.end(), 0.0);
      do {
        for (int q = 0; q < n_inst; ++q) {
          // Slot-major: build the length-Nc product in pbuf one inner-map matrix
          // at a time (sequential access per Hsz[i], one reused buffer, no alloc),
          // then accumulate into contrib.  Hsz[i] is column-major:
          // element (row, col) = data[row + col * nrow].
          {
            const double* Hd = Hdat[0];
            const int nr = Hnrow[0], row = b[0] - 1;
            const int* ip = inst_ptr[q][0];
            for (int c = 0; c < Nc; ++c)
              pbuf[c] = Hd[row + (std::size_t)(ip[c] - 1) * nr];
          }
          for (int i = 1; i < m; ++i) {
            const double* Hd = Hdat[i];
            const int nr = Hnrow[i], row = b[i] - 1;
            const int* ip = inst_ptr[q][i];
            for (int c = 0; c < Nc; ++c)
              pbuf[c] *= Hd[row + (std::size_t)(ip[c] - 1) * nr];
          }
          for (int c = 0; c < Nc; ++c) contrib[c] += pbuf[c];
        }
      } while (std::next_permutation(b.begin(), b.end()));

      double* drow = &D(a, 0);
      // D is column-major: D(a, c) = D[a + c * n_out].
      for (int c = 0; c < Nc; ++c) drow[(std::size_t)c * n_out] += vt * contrib[c];
    }
  }
  return D;
}
