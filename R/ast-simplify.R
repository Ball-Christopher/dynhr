## R/ast-simplify.R
## --------------------------------------------------------------------------
## Algebraic simplification of AST nodes: constant folding, identity rules,
## one-pass and fixpoint variants. Consumed by ast-differentiate.R and the
## compile-static / compile-dynamic builders.
##
## Phase-1 split from jacobian-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Check if an AST node represents zero
#' @noRd
ast_is_zero <- function(node) {
    !is.null(node) && node$type == "number" && node$value == 0
}

#' Check if an AST node represents one
#' @noRd
ast_is_one <- function(node) {
    !is.null(node) && node$type == "number" && node$value == 1
}

#' Check structural equality of two AST nodes
#' @param a,b   AST nodes to compare.
#' @param depth Current recursion depth (internal, for stack overflow guard).
#' @noRd
ast_equal <- function(a, b, depth = 0L) {
    if (depth > 2000L) return(TRUE)  # bail out: assume equal to avoid overflow
    if (is.null(a) && is.null(b)) return(TRUE)
    if (is.null(a) || is.null(b)) return(FALSE)
    if (a$type != b$type) return(FALSE)
    switch(a$type,
        "number"    = a$value == b$value,
        "variable"  = a$name == b$name && a$lead_lag == b$lead_lag,
        "parameter" = a$name == b$name,
        "local_variable" = a$name == b$name,
        "binop"     = a$op == b$op && ast_equal(a$left, b$left, depth + 1L) &&
                      ast_equal(a$right, b$right, depth + 1L),
        "unaryop"   = a$op == b$op && ast_equal(a$operand, b$operand, depth + 1L),
        "funcall"   = {
            a$name == b$name && length(a$args) == length(b$args) &&
            all(mapply(function(x, y) ast_equal(x, y, depth + 1L), a$args, b$args))
        },
        FALSE
    )
}

#' Check if an AST contains any variable nodes (non-constant)
#' @noRd
ast_contains_variable <- function(node) {
    if (is.null(node)) return(FALSE)
    switch(node$type,
        "number"    = FALSE,
        "parameter" = FALSE,
        "local_variable" = FALSE,
        "variable"  = TRUE,
        "binop"     = ast_contains_variable(node$left) ||
                      ast_contains_variable(node$right),
        "unaryop"   = ast_contains_variable(node$operand),
        "funcall"   = {
            if (node$name == "STEADY_STATE") return(FALSE)
            any(vapply(node$args, ast_contains_variable, logical(1)))
        },
        FALSE
    )
}


#' Simplify an AST algebraically (one pass)
#'
#' @param node  AST node to simplify.
#' @param depth Current recursion depth (internal, for stack overflow guard).
#' @noRd
ast_simplify_once <- function(node, depth = 0L) {
    if (is.null(node)) return(ast_number(0))
    # Guard against node stack overflow: R's default recursion limit is ~5000.
    # Deeply nested ASTs from repeated differentiation can hit this limit.
    # If we exceed threshold, return the node unsimplified rather than crashing.
    if (depth > 2000L) return(node)

    switch(node$type,
        "number"    = node,
        "variable"  = node,
        "parameter" = node,
        "local_variable" = node,
        "binop" = {
            L <- ast_simplify_once(node$left, depth + 1L)
            R <- ast_simplify_once(node$right, depth + 1L)

            # Constant folding: number OP number
            if (L$type == "number" && R$type == "number") {
                val <- switch(node$op,
                    "+" = L$value + R$value,
                    "-" = L$value - R$value,
                    "*" = L$value * R$value,
                    "/" = if (R$value != 0) L$value / R$value else NULL,
                    "^" = L$value ^ R$value,
                    NULL
                )
                if (!is.null(val) && is.finite(val))
                    return(ast_number(val))
            }

            switch(node$op,
                "+" = {
                    if (ast_is_zero(L)) return(R)
                    if (ast_is_zero(R)) return(L)
                    ast_binop("+", L, R)
                },
                "-" = {
                    if (ast_is_zero(R)) return(L)
                    if (ast_is_zero(L)) return(ast_unaryop("-", R))
                    if (ast_equal(L, R)) return(ast_number(0))
                    ast_binop("-", L, R)
                },
                "*" = {
                    if (ast_is_zero(L) || ast_is_zero(R)) return(ast_number(0))
                    if (ast_is_one(L)) return(R)
                    if (ast_is_one(R)) return(L)
                    if (L$type == "number" && L$value == -1)
                        return(ast_unaryop("-", R))
                    if (R$type == "number" && R$value == -1)
                        return(ast_unaryop("-", L))
                    ast_binop("*", L, R)
                },
                "/" = {
                    if (ast_is_zero(L)) return(ast_number(0))
                    if (ast_is_one(R)) return(L)
                    if (ast_equal(L, R)) return(ast_number(1))
                    ast_binop("/", L, R)
                },
                "^" = {
                    if (ast_is_zero(R)) return(ast_number(1))
                    if (ast_is_one(R))  return(L)
                    if (ast_is_zero(L)) return(ast_number(0))
                    if (ast_is_one(L))  return(ast_number(1))
                    ast_binop("^", L, R)
                },
                ast_binop(node$op, L, R)
            )
        },
        "unaryop" = {
            inner <- ast_simplify_once(node$operand, depth + 1L)
            if (node$op == "-") {
                if (ast_is_zero(inner)) return(ast_number(0))
                if (inner$type == "number") return(ast_number(-inner$value))
                if (inner$type == "unaryop" && inner$op == "-")
                    return(inner$operand)
            }
            ast_unaryop(node$op, inner)
        },
        "funcall" = {
            args <- lapply(node$args, function(a) ast_simplify_once(a, depth + 1L))
            ast_funcall(node$name, args)
        },
        node
    )
}


#' Shallow (single-level) simplification assuming children are already simplified
#'
#' Applies exactly the same algebraic identities as \code{ast_simplify_once} to
#' ONE node, but WITHOUT recursing into the children -- the caller guarantees
#' \code{node}'s immediate children are already in normal form. Because
#' \code{ast_simplify_once} is itself bottom-up (it normalises children before
#' applying any rule), simplifying-on-construction with this helper produces the
#' same normal form as a full \code{ast_simplify} fixpoint, but in O(1) per
#' constructed node instead of re-walking the whole subtree. This is what makes
#' \code{ast_differentiate} O(result size) rather than O(result size x depth x
#' fixpoint-iterations): the differentiator builds its output bottom-up with the
#' smart constructors below, so no node is ever re-simplified.
#'
#' @param node AST node whose children are already simplified.
#' @noRd
ast_simplify_shallow <- function(node) {
    switch(node$type,
        "binop" = {
            L <- node$left; R <- node$right
            if (L$type == "number" && R$type == "number") {
                val <- switch(node$op,
                    "+" = L$value + R$value,
                    "-" = L$value - R$value,
                    "*" = L$value * R$value,
                    "/" = if (R$value != 0) L$value / R$value else NULL,
                    "^" = L$value ^ R$value,
                    NULL
                )
                if (!is.null(val) && is.finite(val)) return(ast_number(val))
            }
            switch(node$op,
                "+" = {
                    if (ast_is_zero(L)) return(R)
                    if (ast_is_zero(R)) return(L)
                    node
                },
                "-" = {
                    if (ast_is_zero(R)) return(L)
                    if (ast_is_zero(L)) return(ast_unaryop("-", R))
                    if (ast_equal(L, R)) return(ast_number(0))
                    node
                },
                "*" = {
                    if (ast_is_zero(L) || ast_is_zero(R)) return(ast_number(0))
                    if (ast_is_one(L)) return(R)
                    if (ast_is_one(R)) return(L)
                    if (L$type == "number" && L$value == -1)
                        return(ast_unaryop("-", R))
                    if (R$type == "number" && R$value == -1)
                        return(ast_unaryop("-", L))
                    node
                },
                "/" = {
                    if (ast_is_zero(L)) return(ast_number(0))
                    if (ast_is_one(R)) return(L)
                    if (ast_equal(L, R)) return(ast_number(1))
                    node
                },
                "^" = {
                    if (ast_is_zero(R)) return(ast_number(1))
                    if (ast_is_one(R))  return(L)
                    if (ast_is_zero(L)) return(ast_number(0))
                    if (ast_is_one(L))  return(ast_number(1))
                    node
                },
                node
            )
        },
        "unaryop" = {
            inner <- node$operand
            if (node$op == "-") {
                if (ast_is_zero(inner)) return(ast_number(0))
                if (inner$type == "number") return(ast_number(-inner$value))
                if (inner$type == "unaryop" && inner$op == "-")
                    return(inner$operand)
            }
            node
        },
        node
    )
}

## Smart constructors: build a node from already-simplified children and apply
## the shallow simplification rules in one step. Used by ast_differentiate to
## emit normal-form output without any post-hoc ast_simplify pass.
.ds_binop <- function(op, L, R) ast_simplify_shallow(ast_binop(op, L, R))
.ds_unaryop <- function(op, x)  ast_simplify_shallow(ast_unaryop(op, x))

#' Simplify an AST to fixpoint (repeated application)
#'
#' @param node AST node.
#' @param max_iter Maximum iterations to prevent infinite loops.
#' @return Simplified AST node.
#' @noRd
ast_simplify <- function(node, max_iter = 20L) {
    for (i in seq_len(max_iter)) {
        simplified <- ast_simplify_once(node)
        if (ast_equal(simplified, node)) break
        node <- simplified
    }
    node
}
