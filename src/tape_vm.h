// tape_vm.h -- shared stack-machine inner loop for the dynamic-derivative
// "tapes" (Jacobian + Hessian2/Hessian3).
//
// Both eval_jac_tape_cpp (scatters each expression into J(out_row,out_col)) and
// eval_triplet_tape_cpp (packs a dense value vector v[e]) execute the SAME
// postfix opcode stream produced by R/ast-tape.R (.compile_ast_tape). To make
// the "opcodes MUST match .TAPE_OP in R/ast-tape.R" invariant impossible to
// violate per-kernel, the per-expression evaluation lives here as a single
// inline function -- one source of truth for the opcode switch.
//
// Opcodes MUST match .TAPE_OP in R/ast-tape.R.

#ifndef DYNHR_TAPE_VM_H
#define DYNHR_TAPE_VM_H

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

namespace dynhr_tape {

enum {
  OP_CONST = 0, OP_DY = 1, OP_PARAM = 2, OP_SS = 3,
  OP_ADD = 10, OP_SUB = 11, OP_MUL = 12, OP_DIV = 13, OP_POW = 14, OP_NEG = 20,
  OP_EXP = 30, OP_LOG = 31, OP_SQRT = 32, OP_ABS = 33, OP_SIGN = 34,
  OP_SIN = 35, OP_COS = 36, OP_TAN = 37, OP_ASIN = 38, OP_ACOS = 39, OP_ATAN = 40,
  OP_SINH = 41, OP_COSH = 42, OP_TANH = 43,
  OP_NORMCDF = 44, OP_NORMPDF = 45
};

// Evaluate ONE postfix expression of `len` instructions starting at op/ia/da
// index `ip` (advanced in place). `stk` is cleared on entry and the single
// remaining stack value (the expression result) is returned.
static inline double eval_one_expr(int& ip, const int len,
                                   const arma::ivec& op,
                                   const arma::ivec& ia,
                                   const arma::vec&  da,
                                   const arma::vec&  dy,
                                   const arma::vec&  params,
                                   const arma::vec&  ss,
                                   std::vector<double>& stk) {
  stk.clear();
  for (int k = 0; k < len; ++k, ++ip) {
    switch (op[ip]) {
      case OP_CONST: stk.push_back(da[ip]); break;
      case OP_DY:    stk.push_back(dy[ia[ip]]); break;
      case OP_PARAM: stk.push_back(params[ia[ip]]); break;
      case OP_SS:    stk.push_back(ss[ia[ip]]); break;
      case OP_ADD: { double b=stk.back(); stk.pop_back(); stk.back()+=b; break; }
      case OP_SUB: { double b=stk.back(); stk.pop_back(); stk.back()-=b; break; }
      case OP_MUL: { double b=stk.back(); stk.pop_back(); stk.back()*=b; break; }
      case OP_DIV: { double b=stk.back(); stk.pop_back(); stk.back()/=b; break; }
      case OP_POW: { double b=stk.back(); stk.pop_back();
                     stk.back()=std::pow(stk.back(), b); break; }
      case OP_NEG:  stk.back() = -stk.back();          break;
      case OP_EXP:  stk.back() = std::exp(stk.back());  break;
      case OP_LOG:  stk.back() = std::log(stk.back());  break;
      case OP_SQRT: stk.back() = std::sqrt(stk.back()); break;
      case OP_ABS:  stk.back() = std::fabs(stk.back()); break;
      case OP_SIGN: stk.back() = (stk.back()>0) - (stk.back()<0); break;
      case OP_SIN:  stk.back() = std::sin(stk.back());  break;
      case OP_COS:  stk.back() = std::cos(stk.back());  break;
      case OP_TAN:  stk.back() = std::tan(stk.back());  break;
      case OP_ASIN: stk.back() = std::asin(stk.back()); break;
      case OP_ACOS: stk.back() = std::acos(stk.back()); break;
      case OP_ATAN: stk.back() = std::atan(stk.back()); break;
      case OP_SINH: stk.back() = std::sinh(stk.back()); break;
      case OP_COSH: stk.back() = std::cosh(stk.back()); break;
      case OP_TANH: stk.back() = std::tanh(stk.back()); break;
      // normcdf(x) = standard normal CDF = pnorm(x,0,1); normpdf(x) = dnorm(x,0,1)
      case OP_NORMCDF: stk.back() = R::pnorm(stk.back(), 0.0, 1.0, 1, 0); break;
      case OP_NORMPDF: stk.back() = R::dnorm(stk.back(), 0.0, 1.0, 0);    break;
      default: Rcpp::stop("tape_vm: unknown opcode %d", (int)op[ip]);
    }
  }
  return stk.back();
}

} // namespace dynhr_tape

#endif // DYNHR_TAPE_VM_H
