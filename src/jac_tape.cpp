// jac_tape.cpp -- stack-machine evaluators for the dynamic-derivative "tapes".
//
// The dynamic Jacobian / Hessians are otherwise interpreted R closures (built by
// eval(parse()) in R/compile-dynamic.R) whose per-draw cost is dominated by
// named-vector string lookups. .compile_ast_tape() (R/ast-tape.R) flattens the
// differentiated-residual ASTs into postfix opcode/operand arrays; these
// functions execute them with positional integer indexing in a single call.
//
//   eval_jac_tape_cpp     scatters each expression's result into J(out_row,out_col).
//   eval_triplet_tape_cpp packs a dense value vector v[e] (one per expression),
//                         used by the Hessian2/Hessian3 tapes whose R consumers
//                         do the scatter+symmetry themselves.
//
// One generic evaluator handles every model -- no runtime C++ compilation, no
// toolchain dependency. Bit-parity with the R closures is to machine precision
// (only R `^` vs std::pow can differ in the last bit); asserted by
// test-jac-tape-parity.R / test-hess-tape-parity.R. The opcode switch is shared
// via tape_vm.h (single source of truth for "opcodes MUST match .TAPE_OP").

#include <RcppArmadillo.h>
#include <vector>
#include "tape_vm.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat eval_jac_tape_cpp(const arma::vec& dy,
                            const arma::vec& params,
                            const arma::vec& ss,
                            const arma::ivec& op,
                            const arma::ivec& ia,
                            const arma::vec&  da,
                            const arma::ivec& expr_len,
                            const arma::ivec& out_row,
                            const arma::ivec& out_col,
                            int n_eq,
                            int total_cols) {
  arma::mat J(n_eq, total_cols, arma::fill::zeros);
  const int n_expr = expr_len.n_elem;
  std::vector<double> stk;
  stk.reserve(128);
  int ip = 0;
  for (int e = 0; e < n_expr; ++e) {
    const double val = dynhr_tape::eval_one_expr(ip, expr_len[e], op, ia, da,
                                                 dy, params, ss, stk);
    J(out_row[e], out_col[e]) = val;
  }
  return J;
}

// [[Rcpp::export]]
arma::vec eval_triplet_tape_cpp(const arma::vec& dy,
                                const arma::vec& params,
                                const arma::vec& ss,
                                const arma::ivec& op,
                                const arma::ivec& ia,
                                const arma::vec&  da,
                                const arma::ivec& expr_len) {
  const int n_expr = expr_len.n_elem;
  arma::vec out(n_expr, arma::fill::zeros);
  std::vector<double> stk;
  stk.reserve(128);
  int ip = 0;
  for (int e = 0; e < n_expr; ++e) {
    out[e] = dynhr_tape::eval_one_expr(ip, expr_len[e], op, ia, da,
                                       dy, params, ss, stk);
  }
  return out;
}
