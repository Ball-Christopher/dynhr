## R/ast-differentiate.R
## --------------------------------------------------------------------------
## Symbolic differentiation of AST nodes: chain rule, power rule, product
## rule, elementary functions (exp, log, sin, ...), with special handling
## for STEADY_STATE() expressions.
##
## Phase-1 split from jacobian-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Compute the partial derivative of an AST w.r.t. a variable or parameter
#'
#' With \code{wrt = "variable"} (the default) this differentiates w.r.t. an
#' endogenous/exogenous variable at a given lead/lag -- the behaviour used to
#' build the model Jacobian/Hessian. With \code{wrt = "parameter"} it instead
#' differentiates w.r.t. a structural parameter (the \code{var_ll} argument is
#' ignored), holding all variables AND \code{STEADY_STATE()} references fixed.
#' This yields the *explicit* parameter partial used by the symbolic
#' solution-derivative layer (Tier 11 #3); the steady-state chain rule is
#' applied separately by the caller.
#'
#' @param node     AST node to differentiate.
#' @param var_name Name of the variable (or parameter) to differentiate wrt.
#' @param var_ll   Lead/lag of the variable (ignored when \code{wrt="parameter"}).
#' @param wrt      Either \code{"variable"} (default), \code{"parameter"}, or
#'   \code{"local"} (treat a named local_variable as an independent variable for
#'   chain-rule differentiation through #-locals).
#' @return Simplified AST node of the derivative.
#' @noRd
ast_differentiate <- function(node, var_name, var_ll = 0L, wrt = "variable") {
    if (is.null(node)) return(ast_number(0))

    d <- switch(node$type,
        # Constants
        "number"    = ast_number(0),
        "parameter" = {
            if (wrt == "parameter" && node$name == var_name)
                ast_number(1)
            else
                ast_number(0)
        },
        "local_variable" = {
            ## With wrt = "local", treat LOCAL_var_name as an independent variable
            ## (used for the chain-rule pre-computation path in compile-static/dynamic).
            if (wrt == "local" && node$name == var_name)
                ast_number(1)
            else
                ast_number(0)
        },

        # Variable
        "variable" = {
            if (wrt == "variable" && node$name == var_name &&
                node$lead_lag == var_ll)
                ast_number(1)
            else
                ast_number(0)
        },

        # Binary operations
        "binop" = {
            dL <- ast_differentiate(node$left,  var_name, var_ll, wrt)
            dR <- ast_differentiate(node$right, var_name, var_ll, wrt)

            switch(node$op,
                # Relational operators (indicator functions): derivative is 0 a.e.
                # No warning — user-written <= >= < > are legitimate indicators.
                "<="  = ast_number(0),
                ">="  = ast_number(0),
                "<"   = ast_number(0),
                ">"   = ast_number(0),
                "+" = .ds_binop("+", dL, dR),
                "-" = .ds_binop("-", dL, dR),
                # Product rule: d(f*g) = f'*g + f*g'
                "*" = .ds_binop("+",
                    .ds_binop("*", dL, node$right),
                    .ds_binop("*", node$left, dR)
                ),
                # Quotient rule: d(f/g) = (f'*g - f*g') / g^2
                "/" = .ds_binop("/",
                    .ds_binop("-",
                        .ds_binop("*", dL, node$right),
                        .ds_binop("*", node$left, dR)
                    ),
                    .ds_binop("^", node$right, ast_number(2))
                ),
                # Power rule variants
                "^" = {
                    f <- node$left
                    g <- node$right
                    df <- dL
                    dg <- dR
                    f_const <- !ast_contains_variable(f) && ast_is_zero(df)
                    g_const <- !ast_contains_variable(g) && ast_is_zero(dg)

                    if (f_const && g_const) {
                        ast_number(0)
                    } else if (g_const) {
                        # g * f^(g-1) * f'
                        .ds_binop("*",
                            .ds_binop("*", g,
                                .ds_binop("^", f,
                                    .ds_binop("-", g, ast_number(1)))),
                            df)
                    } else if (f_const) {
                        # f^g * log(f) * g'
                        .ds_binop("*",
                            .ds_binop("*",
                                .ds_binop("^", f, g),
                                ast_funcall("log", list(f))),
                            dg)
                    } else {
                        # General: f^g * (g' * log(f) + g * f'/f)
                        .ds_binop("*",
                            .ds_binop("^", f, g),
                            .ds_binop("+",
                                .ds_binop("*", dg, ast_funcall("log", list(f))),
                                .ds_binop("*", g,
                                    .ds_binop("/", df, f))))
                    }
                },
                stop("ast_differentiate: unknown binop '", node$op, "'")
            )
        },

        # Unary operations
        "unaryop" = {
            d_inner <- ast_differentiate(node$operand, var_name, var_ll, wrt)
            if (node$op == "-") .ds_unaryop("-", d_inner)
            else d_inner  # unary + is identity
        },

        # Function calls (chain rule)
        "funcall" = {
            fname <- node$name
            args  <- node$args

            # STEADY_STATE / steady_state are treated as constants here: their
            # parameter dependence (x̄(θ)) is handled by the caller's chain rule.
            if (fname %in% c("STEADY_STATE", "steady_state")) return(ast_number(0))

            # For single-argument functions: d/dx[f(g(x))] = f'(g) * g'
            if (length(args) == 1) {
                g  <- args[[1]]
                dg <- ast_differentiate(g, var_name, var_ll, wrt)

                outer_deriv <- switch(fname,
                    "exp"  = ast_funcall("exp", list(g)),   # exp(g)
                    "log"  = .ds_binop("/", ast_number(1), g),  # 1/g
                    "ln"   = .ds_binop("/", ast_number(1), g),
                    "sqrt" = .ds_binop("/", ast_number(1),
                                 .ds_binop("*", ast_number(2),
                                     ast_funcall("sqrt", list(g)))),
                    "abs"  = ast_funcall("sign", list(g)),
                    "sign" = ast_number(0),
                    "sin"  = ast_funcall("cos", list(g)),
                    "cos"  = .ds_unaryop("-", ast_funcall("sin", list(g))),
                    "tan"  = .ds_binop("/", ast_number(1),
                                 .ds_binop("^",
                                     ast_funcall("cos", list(g)),
                                     ast_number(2))),
                    "asin" = .ds_binop("/", ast_number(1),
                                 ast_funcall("sqrt", list(
                                     .ds_binop("-", ast_number(1),
                                         .ds_binop("^", g, ast_number(2)))))),
                    "acos" = .ds_unaryop("-",
                                 .ds_binop("/", ast_number(1),
                                     ast_funcall("sqrt", list(
                                         .ds_binop("-", ast_number(1),
                                             .ds_binop("^", g, ast_number(2))))))),
                    "atan" = .ds_binop("/", ast_number(1),
                                 .ds_binop("+", ast_number(1),
                                     .ds_binop("^", g, ast_number(2)))),
                    "sinh" = ast_funcall("cosh", list(g)),
                    "cosh" = ast_funcall("sinh", list(g)),
                    "tanh" = .ds_binop("/", ast_number(1),
                                 .ds_binop("^",
                                     ast_funcall("cosh", list(g)),
                                     ast_number(2))),
                    "normcdf" = ast_funcall("normpdf", list(g)),
                    "normpdf" = .ds_binop("*",
                                    .ds_unaryop("-", g),
                                    ast_funcall("normpdf", list(g))),
                    "erf"  = .ds_binop("*",
                                 .ds_binop("/", ast_number(2),
                                     ast_funcall("sqrt", list(
                                         ast_parameter("pi")))),
                                 ast_funcall("exp", list(
                                     .ds_unaryop("-",
                                         .ds_binop("^", g, ast_number(2)))))),
                    "cbrt" = .ds_binop("/", ast_number(1),
                                 .ds_binop("*", ast_number(3),
                                     .ds_binop("^", g,
                                         .ds_binop("/", ast_number(2),
                                             ast_number(3))))),
                    stop("ast_differentiate: unknown function '", fname, "'")
                )

                # Chain rule: f'(g) * g'  (children already in normal form, so a
                # single shallow simplification at the product node suffices).
                return(.ds_binop("*", outer_deriv, dg))
            }

            # max / min: active-regime indicator derivative
            # d/dx max(a,b,...) = sum_i da_i * [a_i >= max of the rest]
            # d/dx min(a,b,...) = sum_i da_i * [a_i <= min of the rest]
            # For the 2-arg case (by far the most common):
            #   d/dx max(a,b) = da*.ind_ge(a,b) + db*.ind_lt(a,b)
            #   d/dx min(a,b) = da*.ind_le(a,b) + db*.ind_gt(a,b)
            # .ind_ge/.ind_gt/.ind_le/.ind_lt are internal helper funcalls that
            # evaluate to 1.0/0.0 -- registered in ast_eval, ast_to_fn_body, and
            # ast_cse_emit so the indicator is numeric in every codegen path.
            if (fname %in% c("max", "min") && length(args) >= 2L) {
                dargs <- lapply(args, ast_differentiate, var_name = var_name,
                                var_ll = var_ll, wrt = wrt)
                # For each argument i, the indicator is 1 iff arg_i is the
                # active branch.  For 2-arg max: ind_i = .ind_ge(a_i, a_{3-i}).
                # For N-arg, approximate by pairwise dominance vs the funcall
                # of the remaining args -- exact at a non-tie evaluation point.
                if (length(args) == 2L) {
                    a <- args[[1]]; b <- args[[2]]
                    da <- dargs[[1]]; db <- dargs[[2]]
                    if (fname == "max") {
                        # indicator for a: a >= b  -> .ind_ge(a,b)
                        # indicator for b: b >  a  -> .ind_gt(b,a)  [tie: a wins]
                        ind_a <- ast_funcall(".ind_ge", list(a, b))
                        ind_b <- ast_funcall(".ind_gt", list(b, a))
                    } else {
                        # min: indicator for a: a <= b -> .ind_le(a,b)
                        # indicator for b: b <  a -> .ind_lt(b,a)  [tie: a wins]
                        ind_a <- ast_funcall(".ind_le", list(a, b))
                        ind_b <- ast_funcall(".ind_lt", list(b, a))
                    }
                    # da*ind_a + db*ind_b  (smart constructors simplify 0*ind=0)
                    return(.ds_binop("+",
                        .ds_binop("*", da, ind_a),
                        .ds_binop("*", db, ind_b)))
                } else {
                    # N > 2 args: for each arg, compare against the max/min of
                    # the remaining N-1 args (exact at a non-tie eval point).
                    terms <- vector("list", length(args))
                    for (i in seq_along(args)) {
                        rest <- args[-i]
                        rest_node <- if (length(rest) == 1L) rest[[1L]]
                                     else ast_funcall(fname, rest)
                        ind <- if (fname == "max")
                            ast_funcall(".ind_ge", list(args[[i]], rest_node))
                        else
                            ast_funcall(".ind_le", list(args[[i]], rest_node))
                        terms[[i]] <- .ds_binop("*", dargs[[i]], ind)
                    }
                    result <- terms[[1L]]
                    for (i in seq(2L, length(terms)))
                        result <- .ds_binop("+", result, terms[[i]])
                    return(result)
                }
            }

            # Other multi-argument functions -- not differentiable
            warning("ast_differentiate: multi-arg function '", fname,
                    "' treated as non-differentiable -- returning 0")
            ast_number(0)
        },

        stop("ast_differentiate: unknown node type '", node$type, "'")
    )

    # `d` is already in normal form: every branch builds its result bottom-up
    # with the smart constructors (.ds_binop/.ds_unaryop), which apply the same
    # algebraic identities as ast_simplify but only to freshly-created nodes
    # (children are already simplified). This replaces the previous per-level
    # ast_simplify() fixpoint -- which re-walked the entire subtree at every
    # recursion level (O(size x depth x iterations)) and was the dominant cost
    # of higher-order symbolic compilation -- with O(result size) work.
    # NOTE: callers must pass an already-simplified `node` (all internal callers
    # differentiate ast_simplify'd residuals or prior normal-form derivatives).
    d
}


#' Explicit partial derivative of an AST w.r.t. a structural parameter
#'
#' Convenience wrapper around \code{ast_differentiate(..., wrt = "parameter")}.
#' Returns the *explicit* partial (variables and \code{STEADY_STATE()} held
#' fixed); the steady-state chain rule is the caller's responsibility.
#'
#' @param node       AST node to differentiate.
#' @param param_name Name of the parameter to differentiate with respect to.
#' @return Simplified AST of \eqn{\partial node / \partial param}.
#' @noRd
ast_differentiate_param <- function(node, param_name) {
    ast_differentiate(node, param_name, 0L, wrt = "parameter")
}


#' Is an AST safe for explicit symbolic parameter differentiation?
#'
#' The explicit-partial path (Tier 11 #3) is exact only when the residual
#' contains no \code{STEADY_STATE()} reference (whose parameter dependence is
#' not captured by the explicit partial) and no non-differentiable construct
#' (\code{max}, \code{min}, or any other multi-argument / unknown function for
#' which \code{ast_differentiate} would silently return 0). When this returns
#' FALSE the caller must fall back to finite differences.
#'
#' @param node AST node (a residual or sub-expression).
#' @return TRUE if every node is smoothly differentiable wrt parameters.
#' @noRd
.ast_param_deriv_safe <- function(node) {
    if (is.null(node)) return(TRUE)
    switch(node$type,
        "number"         = TRUE,
        "parameter"      = TRUE,
        "variable"       = TRUE,
        "local_variable" = TRUE,
        "binop"   = .ast_param_deriv_safe(node$left) &&
                    .ast_param_deriv_safe(node$right),
        "unaryop" = .ast_param_deriv_safe(node$operand),
        "funcall" = {
            # STEADY_STATE() carries an uncaptured parameter dependence; any
            # multi-argument call (max/min/...) has no smooth symbolic rule.
            if (node$name %in% c("STEADY_STATE", "steady_state")) return(FALSE)
            if (length(node$args) != 1L) return(FALSE)
            # Single-arg call: differentiable only if ast_differentiate knows it.
            if (!node$name %in% .AST_DIFF_KNOWN_FUNS) return(FALSE)
            .ast_param_deriv_safe(node$args[[1]])
        },
        FALSE
    )
}

## Single-argument functions for which ast_differentiate has a smooth rule
## (mirrors the switch() in the funcall branch above).
.AST_DIFF_KNOWN_FUNS <- c(
    "exp", "log", "ln", "sqrt", "abs", "sign", "sin", "cos", "tan",
    "asin", "acos", "atan", "sinh", "cosh", "tanh",
    "normcdf", "normpdf", "erf", "cbrt"
)


#' Compute the second partial derivative of an AST w.r.t. two variables
#'
#' Applies ast_differentiate twice: first w.r.t. (var1_name, var1_ll),
#' then w.r.t. (var2_name, var2_ll).  Implements Faà di Bruno / Schwarz
#' symmetry: d²f/(dx dy) = d²f/(dy dx) for smooth f, so callers may pass
#' (var1, var2) in either order.
#'
#' @param node      AST node to differentiate (already simplified).
#' @param var1_name Name of the first variable.
#' @param var1_ll   Lead/lag of the first variable.
#' @param var2_name Name of the second variable.
#' @param var2_ll   Lead/lag of the second variable.
#' @return Simplified AST of d²f/(dvar1 dvar2).
#' @noRd
ast_differentiate2 <- function(node, var1_name, var1_ll = 0L,
                                var2_name, var2_ll = 0L) {
    d1 <- ast_differentiate(node, var1_name, var1_ll)
    ast_differentiate(d1, var2_name, var2_ll)  # already simplified inside
}
