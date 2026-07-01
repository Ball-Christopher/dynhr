## R/ast-tape.R
## --------------------------------------------------------------------------
## Compile differentiated-residual ASTs into a flat stack-machine "tape" that
## is evaluated in C++ (src/jac_tape.cpp, eval_jac_tape_cpp). This replaces the
## per-draw eval() of the interpreted jacobian_fn closure -- whose cost is
## dominated by named-vector string lookups (dy["k__m1"], params["alp"]) -- with
## positional integer indexing in one Rcpp call. No runtime C++ compilation and
## no toolchain dependency: a single generic evaluator runs any model's tape.
##
## Bit-parity with the R closure is to machine precision (the only difference is
## R `^` vs std::pow in the last bit); asserted by test-jac-tape-parity.R.
##
## Opcodes MUST stay in sync with src/jac_tape.cpp.
## --------------------------------------------------------------------------

.TAPE_OP <- c(
  CONST = 0L, DY = 1L, PARAM = 2L, SS = 3L,
  ADD = 10L, SUB = 11L, MUL = 12L, DIV = 13L, POW = 14L, NEG = 20L,
  EXP = 30L, LOG = 31L, SQRT = 32L, ABS = 33L, SIGN = 34L,
  SIN = 35L, COS = 36L, TAN = 37L, ASIN = 38L, ACOS = 39L, ATAN = 40L,
  SINH = 41L, COSH = 42L, TANH = 43L,
  NORMCDF = 44L, NORMPDF = 45L
)

.TAPE_FUNC_OP <- c(
  exp = 30L, log = 31L, ln = 31L, sqrt = 32L, abs = 33L, sign = 34L,
  sin = 35L, cos = 36L, tan = 37L, asin = 38L, acos = 39L, atan = 40L,
  sinh = 41L, cosh = 42L, tanh = 43L,
  normcdf = 44L, normpdf = 45L
)

.TAPE_BINOP <- c("+" = 10L, "-" = 11L, "*" = 12L, "/" = 13L, "^" = 14L)

.tape_ll_suffix <- function(ll)
  if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else paste0("__m", abs(ll))

## Recursively emit postfix instructions for one AST node via the `emit` closure
## (supplied by .compile_ast_tape, which writes into its own pre-sized buffers in
## place). Throws on any construct the tape VM does not support (caller catches
## -> R-closure fallback).
.tape_compile_node <- function(node, emit, dy_idx, par_idx, ss_idx) {
  switch(node$type,
    "number" = emit(0L, -1L, node$value),
    "variable" = {
      key <- paste0(node$name, .tape_ll_suffix(node$lead_lag))
      i <- dy_idx[[key]]
      if (is.null(i)) stop("tape: dy key not found: ", key)
      emit(1L, i, 0)
    },
    "parameter" = {
      i <- par_idx[[node$name]]
      if (is.null(i)) stop("tape: parameter not found: ", node$name)
      emit(2L, i, 0)
    },
    "binop" = {
      .tape_compile_node(node$left,  emit, dy_idx, par_idx, ss_idx)
      .tape_compile_node(node$right, emit, dy_idx, par_idx, ss_idx)
      o <- .TAPE_BINOP[[node$op]]
      if (is.null(o)) stop("tape: unsupported binop ", node$op)
      emit(o)
    },
    "unaryop" = {
      .tape_compile_node(node$operand, emit, dy_idx, par_idx, ss_idx)
      if (node$op == "-") emit(20L)            # unary + is a no-op
      else if (node$op != "+") stop("tape: unsupported unaryop ", node$op)
    },
    "funcall" = {
      if (node$name %in% c("STEADY_STATE", "steady_state")) {
        arg <- node$args[[1]]
        if (arg$type == "variable" && !is.null(ss_idx[[arg$name]])) {
          emit(3L, ss_idx[[arg$name]], 0); return(invisible())
        }
        return(.tape_compile_node(arg, emit, dy_idx, par_idx, ss_idx))
      }
      o <- .TAPE_FUNC_OP[[node$name]]
      if (length(node$args) != 1L || is.null(o))
        stop("tape: unsupported funcall ", node$name, "/", length(node$args))
      .tape_compile_node(node$args[[1]], emit, dy_idx, par_idx, ss_idx)
      emit(o)
    },
    "local_variable" = stop("tape: local variables must be substituted first"),
    stop("tape: unknown node type ", node$type)
  )
  invisible()
}

#' Compile a list of differentiated-residual ASTs into a flat stack-machine tape.
#'
#' Generic core shared by the Jacobian tape (which adds out_row/out_col so the
#' C++ kernel can scatter each expression into J) and the Hessian tapes (which
#' just pack a dense value vector v[e], one per AST). Returns NULL on any
#' construct the tape VM does not support (caught from .tape_compile_node).
#'
#' @param asts        list of AST nodes (one per output expression).
#' @param dy_keys     character vector of timing-suffixed variable keys, in the
#'                    positional order the dy vector is built.
#' @param param_names character vector of parameter names (positional order).
#' @param endo_names  character vector of endogenous names (ss positional order).
#' @return list(op, ia, da, expr_len) or NULL on unsupported construct.
#' @noRd
.compile_ast_tape <- function(asts, dy_keys, param_names, endo_names) {
  if (length(asts) == 0L)
    return(list(op = integer(0), ia = integer(0), da = numeric(0),
                expr_len = integer(0)))
  dy_idx  <- setNames(as.list(seq_along(dy_keys) - 1L), dy_keys)
  par_idx <- setNames(as.list(seq_along(param_names) - 1L), param_names)
  ss_idx  <- setNames(as.list(seq_along(endo_names) - 1L), endo_names)

  tryCatch({
    # Growable instruction buffers held as LOCALS and written IN PLACE.
    #
    # The buffers `op`/`ia`/`da` are ordinary locals of this frame; emit() writes
    # them with `<<-`, which modifies the enclosing-frame binding in place
    # (refcount stays 1) -- O(1) per instruction. An earlier version stored the
    # buffers in an environment and wrote `tp$op[n] <- v`; reading `tp$op` for the
    # `[<-` elevates its refcount so R COPIES the whole vector on every write ->
    # O(n^2) over the tape (this was ~50-70% of order-3 compile, super-linear per
    # emit). Element writes via the enclosing-frame local are ~25x faster (verified
    # micro-bench). Growth is geometric (amortized O(1)); emit() is allocated once.
    cap <- 1024L
    n   <- 0L
    op  <- integer(cap); ia <- integer(cap); da <- numeric(cap)
    emit <- function(opc, ia_arg = -1L, da_arg = 0) {
      n2 <- n + 1L
      if (n2 > cap) {
        nc <- cap * 2L
        length(op) <<- nc; length(ia) <<- nc; length(da) <<- nc
        cap <<- nc
      }
      op[n2] <<- opc; ia[n2] <<- ia_arg; da[n2] <<- da_arg
      n <<- n2
    }
    nt <- length(asts)
    expr_len <- integer(nt)
    for (k in seq_len(nt)) {
      before <- n
      .tape_compile_node(asts[[k]], emit, dy_idx, par_idx, ss_idx)
      expr_len[k] <- n - before
    }
    list(op = op[seq_len(n)], ia = ia[seq_len(n)],
         da = da[seq_len(n)], expr_len = expr_len)
  }, error = function(e) NULL)
}

#' Compile a list of Jacobian triplets (each with $ast, $row, $col) into a tape.
#'
#' @param jac_triplets list of list(row, col, ast, ...).
#' @param dy_keys      character vector of timing-suffixed variable keys, in the
#'                     positional order extract_system_matrices_fast builds dy.
#' @param param_names  character vector of parameter names (positional order).
#' @param endo_names   character vector of endogenous names (ss positional order).
#' @return list(op, ia, da, expr_len, out_row, out_col) with 0-based out indices,
#'         or NULL if any expression uses an unsupported construct.
#' @noRd
compile_jacobian_tape <- function(jac_triplets, dy_keys, param_names, endo_names) {
  if (length(jac_triplets) == 0L)
    return(list(op = integer(0), ia = integer(0), da = numeric(0),
                expr_len = integer(0), out_row = integer(0), out_col = integer(0)))
  base <- .compile_ast_tape(lapply(jac_triplets, `[[`, "ast"),
                            dy_keys, param_names, endo_names)
  if (is.null(base)) return(NULL)
  base$out_row <- vapply(jac_triplets, function(t) t$row - 1L, integer(1))
  base$out_col <- vapply(jac_triplets, function(t) t$col - 1L, integer(1))
  base
}

#' Compile a list of Hessian triplets (each with $ast) into a value-vector tape.
#'
#' Unlike the Jacobian tape there are no out_row/out_col: the C++ kernel
#' (eval_triplet_tape_cpp) returns a dense value vector v[e], one entry per
#' triplet, in input order. The consumers (.compute_model_hessian_symbolic /
#' .compute_model_hessian3_symbolic) do the scatter+symmetry in R from the
#' triplet index fields. NULL on any unsupported construct.
#'
#' @param triplets    list of list(ast, ...).
#' @param dy_keys     timing-suffixed variable keys (dy positional order).
#' @param param_names parameter names (positional order).
#' @param endo_names  endogenous names (ss positional order).
#' @return list(op, ia, da, expr_len) or NULL on unsupported construct.
#' @noRd
compile_hessian_tape <- function(triplets, dy_keys, param_names, endo_names) {
  .compile_ast_tape(lapply(triplets, function(t) t$ast),
                    dy_keys, param_names, endo_names)
}
