## utils-pipe.R
## --------------------------------------------------------------------------
## The single canonical definition of `%||%` for the package.
##
## During phase 0, the monolith files have their local `%||%` definitions
## commented out (with a TODO marker) so this one wins. Once everything is
## split apart in phase 1+, this file becomes the only place `%||%` lives.
##
## Internal utility consolidation.
## --------------------------------------------------------------------------

#' Null-coalescing operator
#'
#' Returns `a` if it is not `NULL`, otherwise `b`. Useful for supplying
#' defaults when a value may be absent.
#'
#' @param a Primary value.
#' @param b Fallback value if `a` is `NULL`.
#'
#' @return `a` if `!is.null(a)`, else `b`.
#'
#' @noRd
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
