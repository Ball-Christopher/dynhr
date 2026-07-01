## R/ast-codegen.R
## --------------------------------------------------------------------------
## AST <-> R-source-string conversion utilities, plus equation -> residual
## transformation, local-variable substitution, lead/lag collapsing for
## the static model, and direct numerical evaluation of an AST.
##
## Phase-1 split from jacobian-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Replace all local_variable AST nodes with their defining expressions
#'
#' @param node       AST node.
#' @param local_vars Named list of (name -> AST expression).
#' @return AST with local variables expanded.
#' @noRd
## TRUE iff NO #-local definition references an endogenous/exogenous VARIABLE
## (directly or via the funcall/binop tree). When TRUE, a local's derivative
## w.r.t. any model variable is 0, so the I20 let-binding fast path (which keeps
## LOCAL_x un-substituted and treats it as a leaf during variable-direction
## differentiation) is CORRECT. When FALSE (a local like `#wedge = q - tot` over
## endo q/tot, e.g. nz_base_nonlinear), the variable Jacobian must chain-rule
## through the local — so the caller MUST substitute locals before differentiating
## w.r.t. variables instead of using the fast path. (local_variable nodes are NOT
## counted: a local referencing a param-only prior local stays param-only; if a
## referenced local DOES touch a variable, that local trips this check itself.)
.locals_are_param_only <- function(local_vars) {
    if (length(local_vars) == 0L) return(TRUE)
    has_var <- function(node) {
        if (is.null(node)) return(FALSE)
        switch(node$type,
            "variable" = TRUE,
            "exo"      = TRUE,
            "number"   = FALSE,
            "parameter" = FALSE,
            "local_variable" = FALSE,
            "binop"   = has_var(node$left) || has_var(node$right),
            "unaryop" = has_var(node$operand),
            "funcall" = any(vapply(node$args, has_var, logical(1))),
            FALSE)
    }
    !any(vapply(local_vars, has_var, logical(1)))
}

ast_substitute_locals <- function(node, local_vars) {
    if (is.null(node) || length(local_vars) == 0) return(node)

    switch(node$type,
        "number"    = node,
        "variable"  = node,
        "parameter" = node,
        "local_variable" = {
            if (node$name %in% names(local_vars))
                ast_substitute_locals(local_vars[[node$name]], local_vars)
            else node
        },
        "binop" = {
            ast_binop(node$op,
                ast_substitute_locals(node$left, local_vars),
                ast_substitute_locals(node$right, local_vars))
        },
        "unaryop" = {
            ast_unaryop(node$op,
                ast_substitute_locals(node$operand, local_vars))
        },
        "funcall" = {
            ast_funcall(node$name,
                lapply(node$args, ast_substitute_locals, local_vars = local_vars))
        },
        node
    )
}


#' Convert an equation to residual form: lhs - rhs
#'
#' @param equation An equation object with $lhs and $rhs.
#' @return AST node representing the residual (should be zero).
#' @noRd
equation_to_residual <- function(equation) {
    ast_simplify(ast_binop("-", equation$lhs, equation$rhs))
}


#' Evaluate an AST numerically at a given point
#'
#' @param node         AST node.
#' @param var_values   Named numeric: names are "varname__0", "varname__m1",
#'                     "varname__p1", etc.
#' @param param_values Named numeric vector of parameter values.
#' @param ss_values    Named numeric vector of steady state values (optional).
#' @return Numeric scalar.
#' @noRd
ast_eval <- function(node, var_values = numeric(0),
                     param_values = numeric(0),
                     ss_values = numeric(0)) {
    if (is.null(node)) return(0)

    switch(node$type,
        "number" = node$value,
        "variable" = {
            suffix <- if (node$lead_lag == 0L) "__0"
                      else if (node$lead_lag > 0L) paste0("__p", node$lead_lag)
                      else paste0("__m", abs(node$lead_lag))
            key <- paste0(node$name, suffix)
            if (key %in% names(var_values)) var_values[[key]]
            else stop("ast_eval: variable '", key, "' not found")
        },
        "parameter" = {
            if (node$name %in% names(param_values)) param_values[[node$name]]
            else stop("ast_eval: parameter '", node$name, "' not found")
        },
        "local_variable" = {
            stop("ast_eval: local variables must be substituted before evaluation")
        },
        "binop" = {
            l <- ast_eval(node$left, var_values, param_values, ss_values)
            r <- ast_eval(node$right, var_values, param_values, ss_values)
            switch(node$op,
                "+" = l + r,  "-" = l - r,
                "*" = l * r,  "/" = l / r,
                "^" = l ^ r,
                # Relational operators: return numeric 1.0/0.0, not logical.
                "<=" = as.numeric(l <= r),
                ">=" = as.numeric(l >= r),
                "<"  = as.numeric(l <  r),
                ">"  = as.numeric(l >  r),
                stop("unknown op: ", node$op)
            )
        },
        "unaryop" = {
            inner <- ast_eval(node$operand, var_values, param_values, ss_values)
            if (node$op == "-") -inner else inner
        },
        "funcall" = {
            if (node$name %in% c("STEADY_STATE", "steady_state")) {
                arg <- node$args[[1]]
                if (arg$type == "variable" && arg$name %in% names(ss_values))
                    return(ss_values[[arg$name]])
                return(ast_eval(arg, var_values, param_values, ss_values))
            }
            vals <- vapply(node$args, ast_eval, numeric(1),
                           var_values = var_values,
                           param_values = param_values,
                           ss_values = ss_values)
            fn <- switch(node$name,
                "exp" = exp, "log" = log, "ln" = log,
                "sqrt" = sqrt, "abs" = abs, "sign" = sign,
                "sin" = sin, "cos" = cos, "tan" = tan,
                "asin" = asin, "acos" = acos, "atan" = atan,
                "sinh" = sinh, "cosh" = cosh, "tanh" = tanh,
                "normcdf" = pnorm, "normpdf" = dnorm,
                "erf" = function(x) 2*pnorm(x*sqrt(2)) - 1,
                "cbrt" = function(x) sign(x) * abs(x)^(1/3),
                "max" = max, "min" = min,
                # Active-regime indicator helpers (used by max/min derivatives).
                # Return 1.0/0.0 so the product with the derivative is numeric.
                ".ind_ge" = function(x, y) as.numeric(x >= y),
                ".ind_gt" = function(x, y) as.numeric(x >  y),
                ".ind_le" = function(x, y) as.numeric(x <= y),
                ".ind_lt" = function(x, y) as.numeric(x <  y),
                stop("ast_eval: unknown function '", node$name, "'")
            )
            do.call(fn, as.list(vals))
        },
        stop("ast_eval: unknown node type '", node$type, "'")
    )
}


#' Replace all variable lead/lags with 0 (current period)
#'
#' Used to construct the static version of the model.
#'
#' @param node AST node.
#' @return AST with all variable timings set to 0.
#' @noRd
ast_collapse_timing <- function(node) {
    if (is.null(node)) return(node)
    switch(node$type,
        "number"    = node,
        "parameter" = node,
        "local_variable" = node,
        "variable"  = ast_variable(node$name, 0L),
        "binop"     = ast_binop(node$op,
                          ast_collapse_timing(node$left),
                          ast_collapse_timing(node$right)),
        "unaryop"   = ast_unaryop(node$op,
                          ast_collapse_timing(node$operand)),
        "funcall"   = {
            if (node$name %in% c("STEADY_STATE", "steady_state")) return(node)  # don't touch
            ast_funcall(node$name,
                lapply(node$args, ast_collapse_timing))
        },
        node
    )
}


#' Shift every variable lead/lag in an AST by a fixed integer offset
#'
#' Used by the Ramsey FOC derivation (\code{ramsey_augment_mod}) to re-time a
#' derivative expression.  The Lagrangian FOC for \eqn{y_i(t)} collects a
#' contribution from the period-(t-1) constraint when \eqn{y_i} appears at a
#' lead (and from period-(t+1) when it appears at a lag).  The symbolic
#' derivative \eqn{\partial g_j / \partial y_i(\pm 1)} is computed in the
#' constraint's own time frame, so to place it in the FOC at time t every
#' variable timing must be shifted by \code{-target_ll} (e.g. a lead derivative,
#' \code{target_ll = +1}, shifts by \code{-1}).  Steady-state references are
#' time-invariant and left untouched.
#'
#' @param node  AST node.
#' @param shift Integer offset added to every variable's \code{lead_lag}.
#' @return AST with all variable lead/lags shifted by \code{shift}.
#' @noRd
ast_shift_timing <- function(node, shift) {
    if (is.null(node) || shift == 0L) return(node)
    shift <- as.integer(shift)
    switch(node$type,
        "number"    = node,
        "parameter" = node,
        "local_variable" = node,
        "variable"  = ast_variable(node$name, node$lead_lag + shift),
        "binop"     = ast_binop(node$op,
                          ast_shift_timing(node$left,  shift),
                          ast_shift_timing(node$right, shift)),
        "unaryop"   = ast_unaryop(node$op,
                          ast_shift_timing(node$operand, shift)),
        "funcall"   = {
            if (node$name %in% c("STEADY_STATE", "steady_state")) return(node)  # don't touch
            ast_funcall(node$name,
                lapply(node$args, ast_shift_timing, shift = shift))
        },
        node
    )
}


#' Shift the lead/lag of NAMED variables in an AST by a fixed offset
#'
#' Selective counterpart to \code{\link{ast_shift_timing}}: only variables whose
#' name is in \code{names} have their \code{lead_lag} adjusted; every other node
#' is left untouched.  Used to implement Dynare's \code{predetermined_variables}
#' convention, where a beginning-of-period stock variable \code{k} written in the
#' model means standard-timing \code{k(-1)}.  EVERY occurrence of such a variable
#' must be re-timed by \code{-1} (so \code{k(+1)} -> \code{k}, \code{k} -> \code{k(-1)},
#' \code{k(+2)} -> \code{k(+1)}, ...), which makes the variable appear as a plain
#' lagged state to the LLI, dynamic Jacobian, classification, and QZ solver.
#'
#' Steady-state references (\code{STEADY_STATE(...)}) are time-invariant and are
#' never descended into, matching \code{ast_shift_timing}.
#'
#' @param node  AST node.
#' @param names Character vector of variable names to shift.
#' @param shift Integer offset added to the matched variables' \code{lead_lag}
#'   (default \code{-1L} for the predetermined convention).
#' @return AST with the named variables' lead/lags shifted; all other nodes
#'   structurally unchanged.
#' @noRd
ast_shift_named_timing <- function(node, names, shift = -1L) {
    if (is.null(node) || length(names) == 0L || shift == 0L) return(node)
    shift <- as.integer(shift)
    switch(node$type,
        "number"    = node,
        "parameter" = node,
        "local_variable" = node,
        "variable"  = if (node$name %in% names)
                          ast_variable(node$name, node$lead_lag + shift)
                      else node,
        "binop"     = ast_binop(node$op,
                          ast_shift_named_timing(node$left,  names, shift),
                          ast_shift_named_timing(node$right, names, shift)),
        "unaryop"   = ast_unaryop(node$op,
                          ast_shift_named_timing(node$operand, names, shift)),
        "funcall"   = {
            if (node$name %in% c("STEADY_STATE", "steady_state")) return(node)
            ast_funcall(node$name,
                lapply(node$args, ast_shift_named_timing,
                       names = names, shift = shift))
        },
        node
    )
}


#' Convert an AST to an R expression string with all variables replaced by SS
#'
#' Used to lower \code{steady_state(EXPR)} when EXPR is a compound expression:
#' every variable reference (at any lead/lag, since SS is time-invariant) is
#' emitted as \code{ss["varname"]} and every parameter as
#' \code{params["paramname"]}, producing a plain arithmetic expression whose
#' value is the steady-state evaluation of EXPR.
#'
#' @param node  AST node (the argument of a steady_state() call).
#' @return Character string of valid R code.
#' @noRd
ast_ss_subst_str <- function(node) {
    if (is.null(node)) return("0")
    switch(node$type,
        "number" = {
            if (node$value == as.integer(node$value) && abs(node$value) < 1e15)
                as.character(as.integer(node$value))
            else
                format(node$value, digits = 15, scientific = FALSE)
        },
        "variable"  = paste0("ss[\"", node$name, "\"]"),
        "parameter" = paste0("params[\"", node$name, "\"]"),
        "local_variable" = paste0("LOCAL_", node$name),
        "binop" = {
            l <- ast_ss_subst_str(node$left)
            r <- ast_ss_subst_str(node$right)
            if (node$op %in% c("<=", ">=", "<", ">"))
                paste0("as.numeric(", l, " ", node$op, " ", r, ")")
            else
                paste0("(", l, " ", node$op, " ", r, ")")
        },
        "unaryop" = {
            inner <- ast_ss_subst_str(node$operand)
            paste0("(", node$op, inner, ")")
        },
        "funcall" = {
            rname <- switch(node$name,
                "ln"       = "log",
                "normcdf"  = "pnorm",
                "normpdf"  = "dnorm",
                "cbrt"     = "((function(x) sign(x)*abs(x)^(1/3)))",
                # Active-regime indicator helpers.
                ".ind_ge" = "((function(x,y) as.numeric(x>=y)))",
                ".ind_gt" = "((function(x,y) as.numeric(x>y)))",
                ".ind_le" = "((function(x,y) as.numeric(x<=y)))",
                ".ind_lt" = "((function(x,y) as.numeric(x<y)))",
                # nested steady_state(...) inside: distribute again
                "STEADY_STATE" =,
                "steady_state" = {
                    if (length(node$args) == 1 &&
                        node$args[[1]]$type == "variable")
                        return(paste0("ss[\"", node$args[[1]]$name, "\"]"))
                    # compound nested: recurse into argument
                    return(ast_ss_subst_str(node$args[[1]]))
                },
                node$name
            )
            args_str <- paste(vapply(node$args, ast_ss_subst_str, character(1)),
                              collapse = ", ")
            paste0(rname, "(", args_str, ")")
        },
        stop("ast_ss_subst_str: unknown node type '", node$type, "'")
    )
}


#' Convert an AST to an R expression string using named vector lookups
#'
#' Variables are referenced as y["name"], exogenous as x["name"],
#' parameters as params["name"], steady state as ss["name"].
#'
#' For the dynamic model with use_timing=TRUE, variable names include
#' timing suffixes: dy["name__0"], dy["name__m1"], dy["name__p1"].
#'
#' @param node         AST node.
#' @param endo_names   Character vector of endogenous variable names.
#' @param exo_names    Character vector of exogenous variable names.
#' @param param_names  Character vector of parameter names.
#' @param use_timing   Logical: if TRUE, include timing suffixes.
#' @param y_prefix     Prefix for endogenous vector (default "y").
#' @param x_prefix     Prefix for exogenous vector (default "x").
#' @return Character string of valid R code.
#' @noRd
ast_to_fn_body <- function(node, endo_names = character(0),
                           exo_names = character(0),
                           param_names = character(0),
                           use_timing = FALSE,
                           y_prefix = "y", x_prefix = "x") {
    if (is.null(node)) return("0")

    switch(node$type,
        "number" = {
            if (node$value == as.integer(node$value) && abs(node$value) < 1e15)
                as.character(as.integer(node$value))
            else
                format(node$value, digits = 15, scientific = FALSE)
        },
        "variable" = {
            if (use_timing) {
                suffix <- if (node$lead_lag == 0L) "__0"
                          else if (node$lead_lag > 0L) paste0("__p", node$lead_lag)
                          else paste0("__m", abs(node$lead_lag))
                key <- paste0(node$name, suffix)
                paste0("dy[\"", key, "\"]")
            } else {
                if (node$name %in% endo_names)
                    paste0(y_prefix, "[\"", node$name, "\"]")
                else if (node$name %in% exo_names)
                    paste0(x_prefix, "[\"", node$name, "\"]")
                else
                    paste0(y_prefix, "[\"", node$name, "\"]")  # fallback
            }
        },
        "parameter" = {
            paste0("params[\"", node$name, "\"]")
        },
        "local_variable" = {
            paste0("LOCAL_", node$name)
        },
        "binop" = {
            l <- ast_to_fn_body(node$left, endo_names, exo_names,
                                param_names, use_timing, y_prefix, x_prefix)
            r <- ast_to_fn_body(node$right, endo_names, exo_names,
                                param_names, use_timing, y_prefix, x_prefix)
            if (node$op %in% c("<=", ">=", "<", ">"))
                paste0("as.numeric(", l, " ", node$op, " ", r, ")")
            else
                paste0("(", l, " ", node$op, " ", r, ")")
        },
        "unaryop" = {
            inner <- ast_to_fn_body(node$operand, endo_names, exo_names,
                                    param_names, use_timing, y_prefix, x_prefix)
            paste0("(", node$op, inner, ")")
        },
        "funcall" = {
            rname <- switch(node$name,
                "ln"           = "log",
                "normcdf"      = "pnorm",
                "normpdf"      = "dnorm",
                "cbrt"         = "((function(x) sign(x)*abs(x)^(1/3)))",
                # Active-regime indicator helpers: emit as.numeric(x >= y) etc.
                ".ind_ge" = "((function(x,y) as.numeric(x>=y)))",
                ".ind_gt" = "((function(x,y) as.numeric(x>y)))",
                ".ind_le" = "((function(x,y) as.numeric(x<=y)))",
                ".ind_lt" = "((function(x,y) as.numeric(x<y)))",
                "STEADY_STATE" =,
                "steady_state" = {
                    # Distribute SS substitution over the argument expression:
                    # bare var -> ss["var"], compound -> recurse replacing each
                    # variable node with its SS value (timing-collapsed).
                    return(ast_ss_subst_str(node$args[[1]]))
                },
                node$name
            )
            args_str <- paste(
                vapply(node$args, ast_to_fn_body, character(1),
                       endo_names = endo_names, exo_names = exo_names,
                       param_names = param_names, use_timing = use_timing,
                       y_prefix = y_prefix, x_prefix = x_prefix),
                collapse = ", ")
            paste0(rname, "(", args_str, ")")
        },
        stop("ast_to_fn_body: unknown node type '", node$type, "'")
    )
}


#' Common-subexpression-eliminated code generation for a list of ASTs
#'
#' Emits a flat list of temporary-variable declarations plus one shallow
#' reference per input AST. Every interior node (binop/unaryop/funcall) is
#' computed once into a temporary \code{.t<i>}; identical subexpressions
#' (hash-consed across ALL the input ASTs) share a temporary. Leaves
#' (numbers / dy[...] / params[...] / ss[...]) are inlined.
#'
#' This serves two purposes that a single flat expression per output cannot:
#'   1. No generated statement is deeply nested, so \code{parse()} never hits
#'      "contextstack overflow" on large higher-order models (e.g. the 134-eq
#'      Andreasen model, whose 3rd derivatives are thousands of chars deep).
#'   2. Common subexpressions are evaluated once (Dynare-style), shrinking the
#'      code and speeding evaluation.
#'
#' The numeric result is identical to evaluating each AST's flat expression:
#' the exact binary structure of every operation is preserved verbatim in the
#' temp RHS, and CSE only de-duplicates structurally-identical subtrees.
#'
#' @param asts list of AST nodes (NULL entries are allowed and skipped).
#' @return list(decls = character vector of "\\.t# <- rhs" lines,
#'              refs  = character vector, same length as `asts`; refs[[k]] is the
#'                      R expression (a temp name or inlined leaf) for asts[[k]],
#'                      or NA for NULL entries).
#' @noRd
ast_cse_emit <- function(asts, endo_names = character(0),
                         exo_names = character(0),
                         param_names = character(0),
                         use_timing = FALSE,
                         y_prefix = "y", x_prefix = "x") {
    memo <- new.env(parent = emptyenv())          # key -> temp name
    # Growable declaration buffer held as LOCALS, written IN PLACE via `<<-`.
    # Storing `v` in an environment and writing `db$v[[i]] <- line` copies the
    # whole list spine on every push (reading db$v for `[[<-` elevates its
    # refcount) -> O(n^2). Writing an enclosing-frame local with `<<-` keeps
    # refcount 1 and modifies in place -- ~1000x faster on large lists (same
    # fix as the tape builder in ast-tape.R). Growth is geometric.
    v  <- vector("list", 1024L)
    n  <- 0L
    tc <- 0L                                       # temporary-variable counter
    push <- function(line) {
        n <<- n + 1L
        if (n > length(v)) length(v) <<- 2L * length(v)
        v[[n]] <<- line
    }
    leaf <- function(node) {
        switch(node$type,
            "number" = {
                if (node$value == as.integer(node$value) && abs(node$value) < 1e15)
                    as.character(as.integer(node$value))
                else
                    format(node$value, digits = 15, scientific = FALSE)
            },
            "variable" = {
                if (use_timing) {
                    suffix <- if (node$lead_lag == 0L) "__0"
                              else if (node$lead_lag > 0L) paste0("__p", node$lead_lag)
                              else paste0("__m", abs(node$lead_lag))
                    paste0("dy[\"", node$name, suffix, "\"]")
                } else if (node$name %in% endo_names) {
                    paste0(y_prefix, "[\"", node$name, "\"]")
                } else if (node$name %in% exo_names) {
                    paste0(x_prefix, "[\"", node$name, "\"]")
                } else {
                    paste0(y_prefix, "[\"", node$name, "\"]")
                }
            },
            "parameter" = paste0("params[\"", node$name, "\"]"),
            "local_variable" = paste0("LOCAL_", node$name),
            NULL  # not a leaf
        )
    }
    go <- function(node) {
        lf <- leaf(node)
        if (!is.null(lf)) return(lf)
        if (node$type == "binop") {
            l <- go(node$left); r <- go(node$right)
            key <- paste0("b", node$op, "", l, "", r)
            m <- memo[[key]]; if (!is.null(m)) return(m)
            if (node$op %in% c("<=", ">=", "<", ">"))
                rhs <- paste0("as.numeric(", l, " ", node$op, " ", r, ")")
            else
                rhs <- paste0("(", l, " ", node$op, " ", r, ")")
        } else if (node$type == "unaryop") {
            x <- go(node$operand)
            key <- paste0("u", node$op, "", x)
            m <- memo[[key]]; if (!is.null(m)) return(m)
            rhs <- paste0("(", node$op, x, ")")
        } else if (node$type == "funcall") {
            rname <- switch(node$name,
                "ln"           = "log",
                "normcdf"      = "pnorm",
                "normpdf"      = "dnorm",
                "cbrt"         = "((function(x) sign(x)*abs(x)^(1/3)))",
                # Active-regime indicator helpers (from max/min derivatives).
                ".ind_ge" = "((function(x,y) as.numeric(x>=y)))",
                ".ind_gt" = "((function(x,y) as.numeric(x>y)))",
                ".ind_le" = "((function(x,y) as.numeric(x<=y)))",
                ".ind_lt" = "((function(x,y) as.numeric(x<y)))",
                "STEADY_STATE" =,
                "steady_state" = {
                    # Distribute SS substitution: replace every variable in the
                    # argument expression with ss["varname"] and emit a temp.
                    expr <- ast_ss_subst_str(node$args[[1]])
                    key  <- paste0("ss_subst:", expr)
                    m    <- memo[[key]]; if (!is.null(m)) return(m)
                    tc <<- tc + 1L
                    nm  <- paste0(".t", tc)
                    push(paste0(nm, " <- ", expr))
                    memo[[key]] <- nm
                    return(nm)
                },
                node$name
            )
            args <- vapply(node$args, go, character(1))
            key <- paste0("f", rname, "", paste(args, collapse = ""))
            m <- memo[[key]]; if (!is.null(m)) return(m)
            rhs <- paste0(rname, "(", paste(args, collapse = ", "), ")")
        } else {
            stop("ast_cse_emit: unknown node type '", node$type, "'")
        }
        tc <<- tc + 1L
        nm <- paste0(".t", tc)
        push(paste0(nm, " <- ", rhs))
        memo[[key]] <- nm
        nm
    }
    refs <- vapply(asts, function(a) if (is.null(a)) NA_character_ else go(a),
                   character(1))
    decls <- if (n > 0L) unlist(v[seq_len(n)], use.names = FALSE)
             else character(0)
    list(decls = decls, refs = refs)
}
