## R/parse-macro.R
## --------------------------------------------------------------------------
## Dynare macro-language preprocessor (Dynare 6/7 language).
##
## A line-oriented expander for Dynare's `@#` macro directives:
##
##   @#define NAME = EXPR           macro variable (any macro value)
##   @#define NAME                  macro variable set to `true`
##   @#define f(x, y) = EXPR        user-defined macro function
##   @#for V in EXPR [when COND]    loop over an array (or a range A:B, A:S:B)
##   @#for (a, b) in EXPR ...       loop over an array of tuples
##   @#endfor
##   @#if COND / @#ifdef NAME / @#ifndef NAME
##     [@#elseif COND]... [@#else] @#endif
##   @#include "path" | EXPR        textual splice, searched in the including
##                                  file's directory, then @#includepath dirs
##   @#includepath "dir" | EXPR     add a directory to the include search path
##   @#error "msg"                  abort (classed `dynhr_error_macro_error`)
##                                  when reached in a taken branch
##   @#echo ..., @#echomacrovars    informational, no effect on the model
##
## Directive lines may end in a `//` or `/* */` comment.
##
## MACRO EXPRESSIONS are parsed and evaluated by a small interpreter written
## here -- no R code is ever evaluated, so a .mod file cannot reach anything
## outside the macro language (the sandbox contract: a call
## to a function that is neither built in nor `@#define`d is a
## `dynhr_error_unsafe_mod_expression`).  The language follows the Dynare 7.1
## macro processor, checked against its `savemacro` output:
##   values     real, string ("..."; '...' also accepted), boolean
##              (true/false), array [a, b], tuple (a, b), user function
##   operators  (lowest to highest) || ; && ; == != ; < > <= >= ; in ;
##              A:B and A:S:B ranges ; | (union) ; & (intersection) ;
##              + - ; * / ; unary - + ! ; ^ ; casts (bool) (real) (string)
##              (tuple) (array) ; indexing v[i], v[[i, j]], v[i:j]
##              `+` concatenates strings and arrays; `-` on arrays is set
##              difference; `*` on arrays is the cartesian product and
##              `A^n` the cartesian power; `==` compares any two values
##              (values of different types are unequal).
##   comprehensions  [EXPR for V in A], [EXPR for V in A when C],
##                   [V in A when C]
##   functions  length isempty isboolean isreal isstring istuple isarray
##              defined exp ln log log10 sin cos tan asin acos atan sqrt cbrt
##              sign floor ceil trunc round erf erfc gamma lgamma sum mod max
##              min normpdf normcdf
## Interpolation `@{EXPR}` renders reals without a decimal point when they
## are integers (and otherwise with enough digits to read back as the same
## double), booleans as true/false, strings without their quotes, arrays as
## `[a, b]` and tuples as `(a, b)`.
##
## `A:B` with B < A is an EMPTY range, as in Dynare.
##
## @#include splices the named file textually (relative to the including
## file's directory) BEFORE macro expansion, exactly as Dynare does.  Nested
## includes are supported; a depth guard (max 32) prevents infinite cycles.
##
## The pass is a strict no-op when the source contains no `@#` directive and
## no `@{...}` interpolation: `expand_macros()` returns the input unchanged
## (byte-for-byte) in that case, so macro-free models are unaffected.
##
## Design: FAIL LOUD on any directive or construct we do not understand rather
## than silently dropping it (a dropped @#for would silently produce a wrong,
## under-specified model).  The error message names the offending directive.
## --------------------------------------------------------------------------

#' Does this .mod source contain any Dynare macro construct?
#'
#' Detects bare `@#` directive lines (ignoring those inside // or /* */
#' comments) and `@{...}` interpolation.  Used to short-circuit the
#' preprocessor to an exact no-op when there is nothing to expand.
#' @noRd
.macro_has_directives <- function(txt) {
  # Strip comments for *detection* only (we never want a commented-out
  # directive to trigger the expander; Dynare does expand macros inside
  # comments, but the replication models never rely on that, and treating
  # them as inert is the conservative choice).
  scan <- gsub("(?s)/\\*.*?\\*/", "\n", txt, perl = TRUE)
  scan <- gsub("//[^\n]*", "", scan, perl = TRUE)
  scan <- gsub("%[^\n]*", "", scan, perl = TRUE)
  # (?m) so `^` matches at the start of every line, not just the whole string.
  grepl("(?m)^[ \t]*@#", scan, perl = TRUE) ||
    grepl("@\\{", scan, perl = TRUE)
}

## Abort from the macro processor (every macro error carries this class).
.macro_abort <- function(..., class = NULL) {
  .dynhr_abort("dynhr macro preprocessor: ", ...,
               class = c(class, "dynhr_error_macro"))
}

# ---------------------------------------------------------------------------
# Macro values
#   real -> numeric(1); string -> character(1); boolean -> logical(1);
#   array -> unclassed list; tuple -> list of class dynhr_macro_tuple;
#   function -> list(args, body, name) of class dynhr_macro_fn.
# ---------------------------------------------------------------------------

.mv_tuple <- function(elems) structure(elems, class = "dynhr_macro_tuple")

.mv_type <- function(x) {
  if (inherits(x, "dynhr_macro_fn")) return("function")
  if (inherits(x, "dynhr_macro_tuple")) return("tuple")
  if (is.list(x)) return("array")
  if (is.logical(x)) return("boolean")
  if (is.numeric(x)) return("real")
  if (is.character(x)) return("string")
  "unknown"
}

## Elements of an array or tuple as a plain list.
.mv_elems <- function(x) {
  x <- unclass(x)
  attributes(x) <- NULL
  x
}

.mv_equal <- function(a, b) {
  ta <- .mv_type(a)
  if (!identical(ta, .mv_type(b))) return(FALSE)
  if (ta %in% c("array", "tuple")) {
    a <- .mv_elems(a); b <- .mv_elems(b)
    if (length(a) != length(b)) return(FALSE)
    for (k in seq_along(a)) if (!.mv_equal(a[[k]], b[[k]])) return(FALSE)
    return(TRUE)
  }
  if (ta == "function") return(identical(a, b))
  isTRUE(a == b)
}

.mv_contains <- function(arr, x) {
  for (el in .mv_elems(arr)) if (.mv_equal(el, x)) return(TRUE)
  FALSE
}

## Boolean value of a condition: a boolean, or a real (nonzero is true).
.mv_bool <- function(x, what) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  if (is.numeric(x) && length(x) == 1L && !is.na(x)) return(x != 0)
  .macro_abort(what, " needs a boolean or a real, not a ", .mv_type(x), ".",
               class = "dynhr_error_macro_type")
}

.mv_real <- function(x, what) {
  if (is.numeric(x) && length(x) == 1L) return(x)
  .macro_abort(what, " needs a real, not a ", .mv_type(x), ".",
               class = "dynhr_error_macro_type")
}

#' Render a macro value for splicing into model text.
#'
#' Reals render without a decimal point when integer-valued (so `@{j}` with
#' j=2 yields `2`, not `2.0` -- Dynare builds identifiers like `ln_p2`), and
#' otherwise with the fewest digits that read back as the same double.
#' Booleans render as `true`/`false`, strings verbatim (no quotes), arrays as
#' `[a, b]` and tuples as `(a, b)`, as Dynare's macro processor prints them.
#' @noRd
.macro_render <- function(val) {
  switch(.mv_type(val),
    real = {
      if (is.finite(val) && val == round(val) &&
          abs(val) < .Machine$integer.max)
        return(format(as.integer(round(val)), trim = TRUE))
      # Shortest text that reads back as the SAME double (up to 17 significant
      # digits). format()'s default of 7 digits silently rounded e.g. 1/3 to
      # 0.3333333.
      if (is.finite(val)) return(.wm_num(val, "macro value"))
      if (is.nan(val)) return("nan")
      if (val > 0) "inf" else "-inf"
    },
    string  = val,
    boolean = if (isTRUE(val)) "true" else "false",
    array   = paste0("[", paste(vapply(.mv_elems(val), .macro_render,
                                       character(1)), collapse = ", "), "]"),
    tuple   = paste0("(", paste(vapply(.mv_elems(val), .macro_render,
                                       character(1)), collapse = ", "), ")"),
    .macro_abort("a ", .mv_type(val), " cannot be rendered into the model ",
                 "text.", class = "dynhr_error_macro_type"))
}

# ---------------------------------------------------------------------------
# Tokenizer and parser for macro expressions
# ---------------------------------------------------------------------------

.MACRO_TOKEN_RES <- c(
  ws  = "^\\s+",
  num = "^(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?",
  str = "^(?:\"[^\"]*\"|'[^']*')",
  ## `.` is accepted inside a name only so that an R-looking call such as
  ## `file.create(...)` is reported as an unknown (refused) function.
  id  = "^[A-Za-z_][A-Za-z0-9_.]*",
  op  = "^(?:&&|\\|\\||==|!=|<=|>=|[-+*/^<>!:|&()\\[\\],])"
)

.macro_tokenize <- function(expr, label) {
  types <- character(0)
  vals  <- character(0)
  rest  <- expr
  while (nzchar(rest)) {
    hit <- FALSE
    for (ty in names(.MACRO_TOKEN_RES)) {
      m <- regexpr(.MACRO_TOKEN_RES[[ty]], rest, perl = TRUE)
      if (m == 1L) {
        len <- attr(m, "match.length")
        if (ty != "ws") {
          types <- c(types, ty)
          vals  <- c(vals, substr(rest, 1L, len))
        }
        rest <- substr(rest, len + 1L, nchar(rest))
        hit <- TRUE
        break
      }
    }
    if (!hit)
      .macro_abort("unexpected character `", substr(rest, 1L, 1L),
                   "` in the macro expression `", expr, "` (", label, ").",
                   class = "dynhr_error_macro_syntax")
  }
  list(type = c(types, "end"), val = c(vals, ""))
}

## Binary operator precedences (higher binds tighter).  Unary - + ! bind at
## 11, below `^` (so -2^2 is -4, as in Dynare) and above everything else.
.MACRO_BINARY_PREC <- c("||" = 1, "&&" = 2, "==" = 3, "!=" = 3,
                        "<" = 4, ">" = 4, "<=" = 4, ">=" = 4, "in" = 5,
                        ":" = 6, "|" = 7, "&" = 8, "+" = 9, "-" = 9,
                        "*" = 10, "/" = 10, "^" = 12)
.MACRO_UNARY_PREC <- 11
.MACRO_CASTS <- c("bool", "real", "string", "tuple", "array")

#' Parse macro-expression text into an AST
#'
#' @param expr  Expression text.
#' @param label Directive label for error messages.
#' @param what  "expr" for a plain expression, "for" for a `@#for` header
#'   (`V in EXPR [when COND]` or `(a, b) in EXPR [when COND]`).
#' @return For "expr" the AST; for "for" list(vars, iter, cond).
#' @noRd
.macro_parse <- function(expr, label, what = "expr") {
  tk <- .macro_tokenize(expr, label)
  ## Cursor kept in an environment so the nested parse functions can advance
  ## it without `<<-`.
  st <- new.env(parent = emptyenv())
  st$pos <- 1L
  peek_t <- function(k = 0L) tk$type[min(st$pos + k, length(tk$type))]
  peek_v <- function(k = 0L) tk$val[min(st$pos + k, length(tk$val))]
  is_op  <- function(v, k = 0L) peek_t(k) == "op" && peek_v(k) == v
  is_kw  <- function(v, k = 0L) peek_t(k) == "id" && peek_v(k) == v
  fail   <- function(msg)
    .macro_abort(msg, " in the macro expression `", expr, "` (", label, ").",
                 class = "dynhr_error_macro_syntax")
  expect_op <- function(v) {
    if (!is_op(v))
      fail(paste0("expected `", v, "` but found `",
                  if (peek_t() == "end") "end of expression" else peek_v(),
                  "`"))
    st$pos <- st$pos + 1L
  }

  parse_expr <- function(min_prec = 1) {
    left <- parse_unary()
    repeat {
      op <- if (peek_t() == "op" && peek_v() %in% names(.MACRO_BINARY_PREC))
        peek_v() else if (is_kw("in")) "in" else NA_character_
      if (is.na(op)) break
      prec <- .MACRO_BINARY_PREC[[op]]
      if (prec < min_prec) break
      st$pos <- st$pos + 1L
      if (op == ":") {
        mid <- parse_expr(prec + 1)
        if (is_op(":")) {
          st$pos <- st$pos + 1L
          hi <- parse_expr(prec + 1)
          left <- list(k = "range", lo = left, by = mid, hi = hi)
        } else {
          left <- list(k = "range", lo = left, by = NULL, hi = mid)
        }
        next
      }
      right <- parse_expr(prec + 1)
      left <- list(k = "bin", op = op, l = left, r = right)
    }
    left
  }

  parse_unary <- function() {
    if (peek_t() == "op" && peek_v() %in% c("-", "+", "!")) {
      op <- peek_v()
      st$pos <- st$pos + 1L
      return(list(k = "un", op = op, x = parse_expr(.MACRO_UNARY_PREC)))
    }
    parse_postfix(parse_primary())
  }

  parse_postfix <- function(node) {
    force(node)   # parse the primary BEFORE looking at the token after it
    while (is_op("[")) {
      st$pos <- st$pos + 1L
      idx <- list(parse_expr())
      while (is_op(",")) { st$pos <- st$pos + 1L; idx <- c(idx, list(parse_expr())) }
      expect_op("]")
      node <- list(k = "index", x = node,
                   i = if (length(idx) == 1L) idx[[1L]]
                       else list(k = "array", el = idx))
    }
    node
  }

  parse_loop_vars <- function() {
    if (peek_t() == "id") {
      v <- peek_v(); st$pos <- st$pos + 1L
      return(v)
    }
    if (is_op("(")) {
      st$pos <- st$pos + 1L
      vs <- character(0)
      repeat {
        if (peek_t() != "id") fail("expected a loop-variable name")
        vs <- c(vs, peek_v()); st$pos <- st$pos + 1L
        if (is_op(",")) { st$pos <- st$pos + 1L; next }
        break
      }
      expect_op(")")
      return(vs)
    }
    fail("expected a loop variable or a `(a, b)` tuple of loop variables")
  }

  ## Loop variables named by the left operand of `in` in `[V in A when C]`.
  loop_vars_of <- function(node) {
    if (node$k == "var") return(node$name)
    if (node$k == "tuple" &&
        all(vapply(node$el, function(e) e$k == "var", logical(1))))
      return(vapply(node$el, function(e) e$name, character(1)))
    NULL
  }

  parse_primary <- function() {
    ty <- peek_t(); v <- peek_v()
    if (ty == "num") {
      st$pos <- st$pos + 1L
      return(list(k = "const", v = as.numeric(v)))
    }
    if (ty == "str") {
      st$pos <- st$pos + 1L
      return(list(k = "const", v = substr(v, 2L, nchar(v) - 1L)))
    }
    if (ty == "id") {
      st$pos <- st$pos + 1L
      if (v %in% c("true", "false")) return(list(k = "const", v = v == "true"))
      if (is_op("(")) {
        st$pos <- st$pos + 1L
        args <- list()
        if (!is_op(")")) {
          args <- list(parse_expr())
          while (is_op(",")) { st$pos <- st$pos + 1L; args <- c(args, list(parse_expr())) }
        }
        expect_op(")")
        return(list(k = "call", name = v, args = args))
      }
      return(list(k = "var", name = v))
    }
    if (ty == "op" && v == "(") {
      if (peek_t(1L) == "id" && peek_v(1L) %in% .MACRO_CASTS && is_op(")", 2L)) {
        to <- peek_v(1L)
        st$pos <- st$pos + 3L
        return(list(k = "cast", to = to, x = parse_expr(.MACRO_UNARY_PREC)))
      }
      st$pos <- st$pos + 1L
      first <- parse_expr()
      if (is_op(",")) {
        el <- list(first)
        while (is_op(",")) { st$pos <- st$pos + 1L; el <- c(el, list(parse_expr())) }
        expect_op(")")
        return(list(k = "tuple", el = el))
      }
      expect_op(")")
      return(first)
    }
    if (ty == "op" && v == "[") {
      st$pos <- st$pos + 1L
      if (is_op("]")) { st$pos <- st$pos + 1L; return(list(k = "array", el = list())) }
      first <- parse_expr()
      if (is_kw("for")) {
        st$pos <- st$pos + 1L
        vars <- parse_loop_vars()
        if (!is_kw("in")) fail("expected `in` in a comprehension")
        st$pos <- st$pos + 1L
        iter <- parse_expr(.MACRO_BINARY_PREC[["in"]] + 1)
        cond <- NULL
        if (is_kw("when")) { st$pos <- st$pos + 1L; cond <- parse_expr() }
        expect_op("]")
        return(list(k = "comp", map = first, vars = vars, iter = iter,
                    cond = cond))
      }
      if (is_kw("when") && first$k == "bin" && first$op == "in" &&
          !is.null(loop_vars_of(first$l))) {
        st$pos <- st$pos + 1L
        cond <- parse_expr()
        expect_op("]")
        return(list(k = "comp", map = NULL, vars = loop_vars_of(first$l),
                    iter = first$r, cond = cond))
      }
      el <- list(first)
      while (is_op(",")) { st$pos <- st$pos + 1L; el <- c(el, list(parse_expr())) }
      expect_op("]")
      return(list(k = "array", el = el))
    }
    fail(paste0("unexpected `",
                if (ty == "end") "end of expression" else v, "`"))
  }

  if (identical(what, "for")) {
    vars <- parse_loop_vars()
    if (!is_kw("in")) fail("expected `in` after the loop variable")
    st$pos <- st$pos + 1L
    iter <- parse_expr(.MACRO_BINARY_PREC[["in"]] + 1)
    cond <- NULL
    if (is_kw("when")) { st$pos <- st$pos + 1L; cond <- parse_expr() }
    if (peek_t() != "end") fail(paste0("unexpected `", peek_v(), "`"))
    return(list(vars = vars, iter = iter, cond = cond))
  }
  out <- parse_expr()
  if (peek_t() != "end") fail(paste0("unexpected `", peek_v(), "`"))
  out
}

# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------

## Look a macro variable up in the environment chain (ends at emptyenv()).
.macro_lookup <- function(name, env) {
  if (!exists(name, envir = env, inherits = TRUE))
    .macro_abort("unknown macro variable `", name, "`.",
                 class = "dynhr_error_macro_undefined")
  get(name, envir = env, inherits = TRUE)
}

## Bind loop variable(s) to one iteration value in `env`.
.macro_bind_loop <- function(vars, value, env, label) {
  if (length(vars) == 1L) {
    assign(vars, value, envir = env)
    return(invisible(NULL))
  }
  if (!(.mv_type(value) %in% c("tuple", "array")) ||
      length(.mv_elems(value)) != length(vars))
    .macro_abort("in ", label, ", each element must be a tuple of ",
                 length(vars), " values to bind (",
                 paste(vars, collapse = ", "), ").",
                 class = "dynhr_error_macro_type")
  el <- .mv_elems(value)
  for (k in seq_along(vars)) assign(vars[[k]], el[[k]], envir = env)
  invisible(NULL)
}

## Elements an iterable yields (`@#for` and comprehensions).  A scalar yields
## itself once (a `@#for` over a scalar macro variable runs one iteration).
.macro_iter_values <- function(val) {
  if (.mv_type(val) %in% c("array", "tuple")) return(.mv_elems(val))
  list(val)
}

.macro_builtin <- function(name, a, label) {
  n <- length(a)
  need <- function(k) if (!(n %in% k))
    .macro_abort("`", name, "()` takes ", paste(k, collapse = " or "),
                 " argument(s), not ", n, " (", label, ").",
                 class = "dynhr_error_macro_syntax")
  r1 <- function() { need(1L); .mv_real(a[[1L]], paste0("`", name, "()`")) }
  switch(name,
    length = {
      need(1L)
      x <- a[[1L]]
      if (.mv_type(x) %in% c("array", "tuple")) length(.mv_elems(x))
      else if (is.character(x)) nchar(x)
      else .macro_abort("`length()` needs an array, a tuple or a string.",
                        class = "dynhr_error_macro_type")
    },
    isempty = {
      need(1L)
      x <- a[[1L]]
      if (.mv_type(x) %in% c("array", "tuple")) length(.mv_elems(x)) == 0L
      else if (is.character(x)) !nzchar(x)
      else .macro_abort("`isempty()` needs an array, a tuple or a string.",
                        class = "dynhr_error_macro_type")
    },
    isboolean = { need(1L); .mv_type(a[[1L]]) == "boolean" },
    isreal    = { need(1L); .mv_type(a[[1L]]) == "real" },
    isstring  = { need(1L); .mv_type(a[[1L]]) == "string" },
    istuple   = { need(1L); .mv_type(a[[1L]]) == "tuple" },
    isarray   = { need(1L); .mv_type(a[[1L]]) == "array" },
    exp = exp(r1()), ln = log(r1()), log = log(r1()), log10 = log10(r1()),
    sin = sin(r1()), cos = cos(r1()), tan = tan(r1()),
    asin = asin(r1()), acos = acos(r1()), atan = atan(r1()),
    sqrt = sqrt(r1()),
    cbrt = { x <- r1(); sign(x) * abs(x)^(1 / 3) },
    sign = sign(r1()), floor = floor(r1()), ceil = ceiling(r1()),
    trunc = trunc(r1()),
    ## std::round: halves away from zero (R's round() rounds half to even).
    round = { x <- r1(); sign(x) * floor(abs(x) + 0.5) },
    erf  = { x <- r1(); 2 * stats::pnorm(x * sqrt(2)) - 1 },
    erfc = { x <- r1(); 2 * stats::pnorm(-x * sqrt(2)) },
    gamma = gamma(r1()), lgamma = lgamma(r1()),
    sum = {
      need(1L)
      if (.mv_type(a[[1L]]) != "array")
        .macro_abort("`sum()` needs an array of reals.",
                     class = "dynhr_error_macro_type")
      s <- 0
      for (el in .mv_elems(a[[1L]])) s <- s + .mv_real(el, "`sum()`")
      s
    },
    ## std::fmod: the result has the sign of the dividend (mod(-7, 3) = -1).
    mod = {
      need(2L)
      x <- .mv_real(a[[1L]], "`mod()`"); y <- .mv_real(a[[2L]], "`mod()`")
      x - y * trunc(x / y)
    },
    max = { need(2L); max(.mv_real(a[[1L]], "`max()`"), .mv_real(a[[2L]], "`max()`")) },
    min = { need(2L); min(.mv_real(a[[1L]], "`min()`"), .mv_real(a[[2L]], "`min()`")) },
    normpdf = , normcdf = {
      need(c(1L, 3L))
      x  <- .mv_real(a[[1L]], paste0("`", name, "()`"))
      mu <- if (n == 3L) .mv_real(a[[2L]], paste0("`", name, "()`")) else 0
      sg <- if (n == 3L) .mv_real(a[[3L]], paste0("`", name, "()`")) else 1
      if (name == "normpdf") stats::dnorm(x, mu, sg) else stats::pnorm(x, mu, sg)
    },
    NULL)
}

.MACRO_BUILTINS <- c("length", "isempty", "isboolean", "isreal", "isstring",
                     "istuple", "isarray", "exp", "ln", "log", "log10", "sin",
                     "cos", "tan", "asin", "acos", "atan", "sqrt", "cbrt",
                     "sign", "floor", "ceil", "trunc", "round", "erf", "erfc",
                     "gamma", "lgamma", "sum", "mod", "max", "min", "normpdf",
                     "normcdf")

.macro_binop <- function(op, l, r) {
  tl <- .mv_type(l); tr <- .mv_type(r)
  mismatch <- function()
    .macro_abort("type mismatch for the `", op, "` operator (", tl, " ", op,
                 " ", tr, ").", class = "dynhr_error_macro_type")
  both <- function(t) tl == t && tr == t
  switch(op,
    "==" = .mv_equal(l, r),
    "!=" = !.mv_equal(l, r),
    "<" = , ">" = , "<=" = , ">=" = {
      if (!(both("real") || both("string"))) mismatch()
      switch(op, "<" = l < r, ">" = l > r, "<=" = l <= r, ">=" = l >= r)
    },
    "in" = {
      if (!(tr %in% c("array", "tuple"))) mismatch()
      .mv_contains(r, l)
    },
    "+" = {
      if (both("real")) l + r
      else if (both("string")) paste0(l, r)
      else if (both("array")) c(.mv_elems(l), .mv_elems(r))
      else mismatch()
    },
    "-" = {
      if (both("real")) l - r
      else if (both("array"))
        Filter(function(el) !.mv_contains(r, el), .mv_elems(l))
      else mismatch()
    },
    "*" = {
      if (both("real")) return(l * r)
      if (!(tl %in% c("array", "tuple") && tr == "array")) mismatch()
      out <- list()
      for (x in .mv_elems(l)) for (y in .mv_elems(r)) {
        xe <- if (.mv_type(x) == "tuple") .mv_elems(x) else list(x)
        out[[length(out) + 1L]] <- .mv_tuple(c(xe, list(y)))
      }
      out
    },
    "/" = { if (!both("real")) mismatch(); l / r },
    "^" = {
      if (both("real")) return(l^r)
      if (!(tl == "array" && tr == "real") || r < 1 || r != round(r)) mismatch()
      out <- lapply(.mv_elems(l), function(x) .mv_tuple(list(x)))
      if (r >= 2) for (k in 2:r) {
        nxt <- list()
        for (x in out) for (y in .mv_elems(l))
          nxt[[length(nxt) + 1L]] <- .mv_tuple(c(.mv_elems(x), list(y)))
        out <- nxt
      }
      out
    },
    "|" = {
      if (!both("array")) mismatch()
      out <- .mv_elems(l)
      for (el in .mv_elems(r)) if (!.mv_contains(out, el))
        out[[length(out) + 1L]] <- el
      out
    },
    "&" = {
      if (!both("array")) mismatch()
      Filter(function(el) .mv_contains(r, el), .mv_elems(l))
    },
    mismatch())
}

.macro_cast <- function(to, x) {
  tx <- .mv_type(x)
  bad <- function()
    .macro_abort("cannot convert a ", tx, " to ", to, ".",
                 class = "dynhr_error_macro_type")
  switch(to,
    bool = {
      if (tx == "boolean") x
      else if (tx == "real") x != 0
      else if (tx == "string" && x %in% c("true", "false")) x == "true"
      else if (tx %in% c("array", "tuple") && length(.mv_elems(x)) == 1L)
        .macro_cast("bool", .mv_elems(x)[[1L]])
      else bad()
    },
    real = {
      if (tx == "real") x
      else if (tx == "boolean") as.numeric(x)
      else if (tx == "string") {
        s <- trimws(x)
        if (!grepl("^[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?$", s,
                   perl = TRUE)) bad()
        as.numeric(s)
      } else if (tx %in% c("array", "tuple") && length(.mv_elems(x)) == 1L)
        .macro_cast("real", .mv_elems(x)[[1L]])
      else bad()
    },
    string = .macro_render(x),
    tuple  = if (tx %in% c("array", "tuple")) .mv_tuple(.mv_elems(x))
             else .mv_tuple(list(x)),
    array  = if (tx %in% c("array", "tuple")) .mv_elems(x) else list(x))
}

.macro_index <- function(x, i) {
  tx <- .mv_type(x)
  if (!(tx %in% c("array", "tuple", "string")))
    .macro_abort("only an array, a tuple or a string can be indexed, not a ",
                 tx, ".", class = "dynhr_error_macro_type")
  idx <- if (.mv_type(i) == "array")
    vapply(.mv_elems(i), .mv_real, numeric(1), what = "an index")
  else .mv_real(i, "an index")
  len <- if (tx == "string") nchar(x) else length(.mv_elems(x))
  if (any(idx < 1 | idx > len | idx != round(idx)))
    .macro_abort("index out of range (", paste(idx, collapse = ", "),
                 " of a ", tx, " of length ", len, ").",
                 class = "dynhr_error_macro_type")
  if (tx == "string")
    return(paste(substring(x, idx, idx), collapse = ""))
  el <- .mv_elems(x)
  if (.mv_type(i) == "array") el[idx] else el[[idx]]
}

#' Evaluate a macro-expression AST in the macro environment.
#' @noRd
.macro_ev <- function(node, env, label) {
  switch(node$k,
    const = node$v,
    var   = .macro_lookup(node$name, env),
    array = lapply(node$el, .macro_ev, env = env, label = label),
    tuple = .mv_tuple(lapply(node$el, .macro_ev, env = env, label = label)),
    un = {
      x <- .macro_ev(node$x, env, label)
      switch(node$op,
        "-" = -.mv_real(x, "unary `-`"),
        "+" = .mv_real(x, "unary `+`"),
        "!" = !.mv_bool(x, "`!`"))
    },
    bin = {
      if (node$op %in% c("&&", "||")) {
        l <- .mv_bool(.macro_ev(node$l, env, label), paste0("`", node$op, "`"))
        if (node$op == "&&" && !l) return(FALSE)
        if (node$op == "||" && l) return(TRUE)
        return(.mv_bool(.macro_ev(node$r, env, label),
                        paste0("`", node$op, "`")))
      }
      .macro_binop(node$op, .macro_ev(node$l, env, label),
                   .macro_ev(node$r, env, label))
    },
    range = {
      lo <- .mv_real(.macro_ev(node$lo, env, label), "a range bound")
      hi <- .mv_real(.macro_ev(node$hi, env, label), "a range bound")
      by <- if (is.null(node$by)) 1
            else .mv_real(.macro_ev(node$by, env, label), "a range increment")
      if (by == 0)
        .macro_abort("a range increment of 0 (", label, ").",
                     class = "dynhr_error_macro_type")
      # Dynare ranges are EMPTY when they run the wrong way: `1:0` is
      # zero iterations (seq() would give 1, 0).
      if ((by > 0 && hi < lo) || (by < 0 && hi > lo)) return(list())
      as.list(as.numeric(seq(lo, hi, by = by)))
    },
    index = .macro_index(.macro_ev(node$x, env, label),
                         .macro_ev(node$i, env, label)),
    cast = .macro_cast(node$to, .macro_ev(node$x, env, label)),
    comp = {
      out <- list()
      for (v in .macro_iter_values(.macro_ev(node$iter, env, label))) {
        it_env <- new.env(parent = env)
        .macro_bind_loop(node$vars, v, it_env, label)
        if (!is.null(node$cond) &&
            !.mv_bool(.macro_ev(node$cond, it_env, label), "a `when` filter"))
          next
        out[[length(out) + 1L]] <- if (is.null(node$map)) v
                                   else .macro_ev(node$map, it_env, label)
      }
      out
    },
    call = {
      nm <- node$name
      if (nm == "defined") {
        if (length(node$args) != 1L || node$args[[1L]]$k != "var")
          .macro_abort("`defined()` takes one macro-variable name (", label,
                       ").", class = "dynhr_error_macro_syntax")
        return(exists(node$args[[1L]]$name, envir = env, inherits = TRUE))
      }
      args <- lapply(node$args, .macro_ev, env = env, label = label)
      if (nm %in% .MACRO_BUILTINS) return(.macro_builtin(nm, args, label))
      fn <- if (exists(nm, envir = env, inherits = TRUE))
        get(nm, envir = env, inherits = TRUE) else NULL
      if (!inherits(fn, "dynhr_macro_fn"))
        .macro_abort(
          "refuses to evaluate `", nm, "(...)` in ", label, ": `", nm, "` is ",
          "neither a Dynare macro function nor one defined with ",
          "`@#define ", nm, "(...) = ...`. Parsing a .mod file must never ",
          "execute arbitrary code.",
          class = c("dynhr_error_unsafe_mod_expression",
                    "dynhr_error_macro_undefined"))
      if (length(args) != length(fn$args))
        .macro_abort("macro function `", nm, "` takes ", length(fn$args),
                     " argument(s), not ", length(args), " (", label, ").",
                     class = "dynhr_error_macro_syntax")
      call_env <- new.env(parent = env)
      for (k in seq_along(args)) assign(fn$args[[k]], args[[k]], envir = call_env)
      .macro_ev(fn$body, call_env, paste0("macro function `", nm, "`"))
    },
    .macro_abort("internal: unknown macro AST node `", node$k, "`."))
}

#' Evaluate macro-expression text in the macro environment.
#' @noRd
.macro_eval <- function(expr, env, directive_label) {
  .macro_ev(.macro_parse(trimws(expr), directive_label), env, directive_label)
}

## Directive text without a trailing `//` or `/* */` comment (Dynare 6+
## allows inline comments on directive lines).  Quoted strings are skipped.
.macro_strip_inline_comment <- function(s) {
  out <- ""
  rest <- s
  repeat {
    m <- regexpr("\"[^\"]*\"|'[^']*'|//|/\\*", rest, perl = TRUE)
    if (m < 0L) return(trimws(paste0(out, rest)))
    tok <- regmatches(rest, m)
    if (tok == "//") return(trimws(paste0(out, substr(rest, 1L, m - 1L))))
    if (tok == "/*") {
      after <- substr(rest, m + 2L, nchar(rest))
      close <- regexpr("*/", after, fixed = TRUE)
      out <- paste0(out, substr(rest, 1L, m - 1L), " ")
      if (close < 0L) return(trimws(out))
      rest <- substr(after, close + 2L, nchar(after))
      next
    }
    out <- paste0(out, substr(rest, 1L, m + nchar(tok) - 1L))
    rest <- substr(rest, m + nchar(tok), nchar(rest))
  }
}

#' Interpolate every `@{expr}` occurrence in a single line.
#' @noRd
.macro_interp_line <- function(line, env) {
  if (!grepl("@\\{", line, perl = TRUE)) return(line)
  out <- ""
  rest <- line
  repeat {
    m <- regexpr("@\\{", rest, perl = TRUE)
    if (m[1] < 0) {
      out <- paste0(out, rest)
      break
    }
    start <- m[1]
    # The closing `}` is the first one outside a quoted string.
    after <- substring(rest, start + 2L)
    ch <- strsplit(after, "", fixed = TRUE)[[1]]
    quote <- ""
    close <- NA_integer_
    for (k in seq_along(ch)) {
      if (nzchar(quote)) { if (ch[k] == quote) quote <- ""; next }
      if (ch[k] %in% c("\"", "'")) { quote <- ch[k]; next }
      if (ch[k] == "}") { close <- k; break }
    }
    if (is.na(close)) {
      stop("dynhr macro preprocessor: unterminated @{ interpolation in line:\n  ",
           line, call. = FALSE)
    }
    expr <- substring(after, 1L, close - 1L)
    val  <- .macro_eval(expr, env, paste0("@{", expr, "}"))
    out  <- paste0(out, substring(rest, 1L, start - 1L), .macro_render(val))
    rest <- substring(after, close + 1L)
  }
  out
}

## Resolve an @#include target: absolute as is; otherwise the including
## file's directory first, then every @#includepath directory in order.
.macro_resolve_include <- function(inc_rel, env) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", inc_rel)) return(inc_rel)
  inc_dir <- get(".macro_include_dir", envir = env, inherits = TRUE)
  paths <- get(".macro_include_paths", envir = env, inherits = TRUE)
  cands <- c(if (!is.na(inc_dir)) file.path(inc_dir, inc_rel) else inc_rel,
             if (length(paths)) file.path(paths, inc_rel))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[[1L]] else cands[[1L]]
}

## The value of an @#include / @#includepath argument: a quoted string, or a
## macro expression that evaluates to a string.
.macro_path_arg <- function(arg, env, label) {
  arg <- trimws(arg)
  if (grepl('^"[^"]*"$', arg) || grepl("^'[^']*'$", arg))
    return(substring(arg, 2L, nchar(arg) - 1L))
  val <- .macro_eval(arg, env, label)
  if (!is.character(val) || length(val) != 1L)
    stop("dynhr macro preprocessor: ", label, " must name a file with a ",
         "quoted string or a string-valued macro expression; got: ", arg,
         call. = FALSE)
  val
}

#' Recursively expand a block of macro lines.
#'
#' @param lines      Character vector of source lines for this block.
#' @param env        Environment holding @#define vars + active loop indices +
#'                   `.macro_include_dir` (character(1), may be NA),
#'                   `.macro_include_paths` (character) and
#'                   `.macro_include_depth` (integer(1)).
#' @return Character vector of fully expanded lines (no `@#` directives, all
#'         `@{}` interpolated).
#' @noRd
.macro_expand_lines <- function(lines, env) {
  out <- character(0)
  i <- 1L
  n <- length(lines)
  while (i <= n) {
    line <- lines[i]
    trimmed <- trimws(line)
    if (grepl("^@#", trimmed, perl = TRUE)) {
      directive <- .macro_strip_inline_comment(sub("^@#\\s*", "", trimmed))

      # ---- @#define NAME [= EXPR]  |  @#define f(args) = EXPR --------
      if (grepl("^define\\b", directive)) {
        m <- regmatches(directive, regexec(paste0(
          "^define\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(\\(([^)]*)\\))?\\s*",
          "(?:=(?!=)\\s*(.*))?$"), directive, perl = TRUE))[[1]]
        if (length(m) != 5L) {
          stop("dynhr macro preprocessor: @#define must have the form ",
               "`@#define NAME = VALUE`, `@#define NAME` or ",
               "`@#define f(x) = EXPR`; got: ", trimmed, call. = FALSE)
        }
        nm <- m[2]
        has_rhs <- grepl("=", sub("^define\\s+[A-Za-z_][A-Za-z0-9_]*\\s*(\\([^)]*\\))?",
                                  "", directive, perl = TRUE), fixed = TRUE)
        if (nzchar(m[3])) {
          if (!has_rhs || !nzchar(trimws(m[5])))
            stop("dynhr macro preprocessor: the macro function `", nm,
                 "` needs a body: `@#define ", nm, "(...) = EXPR`.",
                 call. = FALSE)
          fargs <- trimws(strsplit(m[4], ",", fixed = TRUE)[[1]])
          fargs <- fargs[nzchar(fargs)]
          if (!all(grepl("^[A-Za-z_][A-Za-z0-9_]*$", fargs)))
            stop("dynhr macro preprocessor: malformed argument list in ",
                 trimmed, call. = FALSE)
          val <- structure(
            list(args = fargs,
                 body = .macro_parse(trimws(m[5]), paste0("@#define ", nm)),
                 name = nm),
            class = "dynhr_macro_fn")
        } else if (!has_rhs) {
          # Dynare 5+: a variable defined without a value is `true`.
          val <- TRUE
        } else {
          val <- .macro_eval(m[5], env, paste0("@#define ", nm))
        }
        assign(nm, val, envir = env)
        i <- i + 1L
        next
      }

      # ---- @#for V in <iter> [when COND] ... @#endfor ------------------
      if (grepl("^for\\b", directive)) {
        header <- trimws(sub("^for\\b", "", directive))
        label  <- paste0("@#for ", header)
        spec   <- .macro_parse(header, label, what = "for")
        # Collect the body up to the matching @#endfor (respecting nesting).
        depth <- 1L
        body_lines <- character(0)
        j <- i + 1L
        while (j <= n) {
          tj <- trimws(lines[j])
          if (grepl("^@#\\s*for\\b", tj, perl = TRUE)) depth <- depth + 1L
          else if (grepl("^@#\\s*endfor\\b", tj, perl = TRUE)) {
            depth <- depth - 1L
            if (depth == 0L) break
          }
          body_lines <- c(body_lines, lines[j])
          j <- j + 1L
        }
        if (depth != 0L) {
          stop("dynhr macro preprocessor: @#for without matching @#endfor ",
               "(starting at: ", trimmed, ").", call. = FALSE)
        }
        values <- .macro_iter_values(.macro_ev(spec$iter, env, label))
        for (v in values) {
          # Child env so the loop index does not clobber an outer define and
          # is scoped to this iteration; nested loops/defines stack here.
          iter_env <- new.env(parent = env)
          .macro_bind_loop(spec$vars, v, iter_env, label)
          if (!is.null(spec$cond) &&
              !.mv_bool(.macro_ev(spec$cond, iter_env, label),
                        "a `when` filter"))
            next
          out <- c(out, .macro_expand_lines(body_lines, iter_env))
        }
        i <- j + 1L
        next
      }

      # ---- @#if / @#ifdef / @#ifndef [@#elseif]... [@#else] @#endif ---
      if (grepl("^(if|ifdef|ifndef)\\b", directive)) {
        # Collect the branches up to @#endif, respecting nested @#if.  Each
        # branch is list(kind, arg, lines); the @#else branch has kind "else".
        branches <- list(list(kind = regmatches(directive, regexpr("^(ifdef|ifndef|if)", directive)),
                              arg = trimws(sub("^(ifdef|ifndef|if)", "", directive)),
                              lines = character(0)))
        depth <- 1L
        seen_else <- FALSE
        j <- i + 1L
        while (j <= n) {
          tj <- trimws(lines[j])
          if (grepl("^@#\\s*(if|ifdef|ifndef)\\b", tj, perl = TRUE)) {
            depth <- depth + 1L
          } else if (grepl("^@#\\s*endif\\b", tj, perl = TRUE)) {
            depth <- depth - 1L
            if (depth == 0L) break
          } else if (depth == 1L &&
                     grepl("^@#\\s*(else|elseif)\\b", tj, perl = TRUE)) {
            d <- .macro_strip_inline_comment(sub("^@#\\s*", "", tj))
            kind <- regmatches(d, regexpr("^(elseif|else)", d))
            if (seen_else)
              stop("dynhr macro preprocessor: @#", kind, " after @#else ",
                   "(in the @#if starting at: ", trimmed, ").", call. = FALSE)
            seen_else <- kind == "else"
            branches[[length(branches) + 1L]] <-
              list(kind = kind, arg = trimws(sub("^(elseif|else)", "", d)),
                   lines = character(0))
            j <- j + 1L
            next
          }
          b <- length(branches)
          branches[[b]]$lines <- c(branches[[b]]$lines, lines[j])
          j <- j + 1L
        }
        if (depth != 0L) {
          stop("dynhr macro preprocessor: @#if without matching @#endif ",
               "(starting at: ", trimmed, ").", call. = FALSE)
        }
        # Conditions are evaluated in order, and only until one is taken.
        for (b in branches) {
          taken <- switch(b$kind,
            ifdef  = exists(b$arg, envir = env, inherits = TRUE),
            ifndef = !exists(b$arg, envir = env, inherits = TRUE),
            "else" = TRUE,
            .mv_bool(.macro_eval(b$arg, env, paste0("@#", b$kind, " ", b$arg)),
                     paste0("the condition of @#", b$kind)))
          if (taken) {
            out <- c(out, .macro_expand_lines(b$lines, env))
            break
          }
        }
        i <- j + 1L
        next
      }

      # ---- @#includepath "dir" ----------------------------------------
      if (grepl("^includepath\\b", directive)) {
        d <- .macro_path_arg(sub("^includepath\\b", "", directive), env,
                             "@#includepath")
        inc_dir <- get(".macro_include_dir", envir = env, inherits = TRUE)
        if (!grepl("^(/|[A-Za-z]:[/\\\\])", d) && !is.na(inc_dir))
          d <- file.path(inc_dir, d)
        root <- get(".macro_root_env", envir = env, inherits = TRUE)
        assign(".macro_include_paths",
               c(get(".macro_include_paths", envir = root), d), envir = root)
        i <- i + 1L
        next
      }

      # ---- @#include "path" | EXPR -------------------------------------
      if (grepl("^include\\b", directive)) {
        inc_rel  <- .macro_path_arg(sub("^include\\b", "", directive), env,
                                    "@#include")
        inc_path <- .macro_resolve_include(inc_rel, env)
        if (!file.exists(inc_path)) {
          stop("dynhr macro preprocessor: @#include file not found: ",
               inc_path, call. = FALSE)
        }
        # Depth guard against infinite include cycles.
        depth_now <- get(".macro_include_depth", envir = env, inherits = TRUE)
        if (depth_now >= 32L) {
          stop("dynhr macro preprocessor: @#include nesting depth exceeded 32 ",
               "(possible cycle); offending file: ", inc_path, call. = FALSE)
        }
        # Read the included file and splice its lines (textual include, just
        # like Dynare -- macro expansion of the included content happens in
        # the current env context, exactly as if the lines were inline).
        inc_txt   <- paste(readLines(inc_path, warn = FALSE), collapse = "\n")
        inc_txt   <- iconv(inc_txt, from = "LATIN1", to = "UTF-8", sub = "?")
        inc_lines <- strsplit(inc_txt, "\n", fixed = TRUE)[[1]]
        # Child env that updates the include dir + depth for nested includes.
        inc_env <- new.env(parent = env)
        assign(".macro_include_dir",   normalizePath(dirname(inc_path), mustWork = FALSE),
               envir = inc_env)
        assign(".macro_include_depth", depth_now + 1L, envir = inc_env)
        out <- c(out, .macro_expand_lines(inc_lines, inc_env))
        # Propagate any @#define assignments made inside the include back to
        # the calling env (mirror Dynare: defines in an included file are
        # visible after the @#include).
        inc_names <- ls(envir = inc_env, all.names = FALSE)
        for (nm in inc_names) {
          assign(nm, get(nm, envir = inc_env), envir = env)
        }
        i <- i + 1L
        next
      }

      # ---- @#error "msg" : abort (only reached in a TAKEN branch) ------
      # Untaken @#if branches and zero-iteration @#for bodies are never
      # expanded, so reaching this line means Dynare would stop here too.
      if (grepl("^error\\b", directive)) {
        msg <- trimws(sub("^error\\b", "", directive))
        if (grepl('^".*"$', msg) || grepl("^'.*'$", msg)) {
          msg <- substring(msg, 2L, nchar(msg) - 1L)
        } else if (nzchar(msg)) {
          msg <- .macro_render(.macro_eval(msg, env, "@#error"))
        }
        .dynhr_abort("dynhr macro preprocessor: @#error: ", msg,
                     class = "dynhr_error_macro_error")
      }

      # ---- @#echo / @#echomacrovars : informational, no model effect ---
      if (grepl("^(echo|echomacrovars)\\b", directive)) {
        i <- i + 1L
        next
      }

      # ---- Anything else: FAIL LOUD ----------------------------------
      stop("dynhr macro preprocessor: unsupported macro directive `@#",
           sub("[\\s(].*$", "", directive, perl = TRUE), "` in line:\n  ",
           trimmed,
           "\nSupported: @#define, @#for/@#endfor, @#if/@#ifdef/@#ifndef/",
           "@#elseif/@#else/@#endif, @#include, @#includepath, @#error, ",
           "@#echo, @#echomacrovars.", call. = FALSE)
    }

    # Plain model line: interpolate any @{...} and emit.
    out <- c(out, .macro_interp_line(line, env))
    i <- i + 1L
  }
  out
}

#' Expand Dynare `@#` macro directives in .mod source text.
#'
#' Runs before the lexer / declaration parsing in `parse_mod()`.  When the
#' source contains no macro directive and no `@{}` interpolation this is an
#' exact byte-for-byte no-op.
#'
#' @param txt     Character string (full .mod file content).
#' @param mod_dir Character(1) or NULL.  Directory used to resolve
#'   `@#include` paths.  `parse_mod()` passes `dirname(source_file)` when
#'   reading from a file; inline-text callers pass NULL (includes disabled
#'   unless an absolute path is given).
#' @return The macro-expanded source text.
#' @noRd
expand_macros <- function(txt, mod_dir = NULL) {
  if (!.macro_has_directives(txt)) return(txt)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  ## Macro expressions are evaluated by the interpreter above, never by
  ## R's eval(); the macro environment only holds macro VALUES and is parented
  ## at emptyenv(), so name lookup cannot leave it.
  env <- new.env(parent = emptyenv())
  if (!is.null(mod_dir) && nzchar(mod_dir)) {
    assign(".macro_include_dir",   normalizePath(mod_dir, mustWork = FALSE),
           envir = env)
  } else {
    assign(".macro_include_dir",   NA_character_, envir = env)
  }
  assign(".macro_include_depth", 0L, envir = env)
  assign(".macro_include_paths", character(0), envir = env)
  assign(".macro_root_env", env, envir = env)
  expanded <- .macro_expand_lines(lines, env)
  paste(expanded, collapse = "\n")
}
