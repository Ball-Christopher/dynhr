## R/hank-grids.R
## --------------------------------------------------------------------------
## Grid construction for heterogeneous-agent (HANK) household blocks.
##
## The asset grid is non-uniform: points cluster near the borrowing constraint
## amin (where policy functions are most curved) and thin out toward amax.  The
## construction is equidistant in logs of a shifted variable, matching the
## `sequence-jacobian` reference package (utilities.discretize.agrid) so its
## published goldens are reproducible.
## --------------------------------------------------------------------------


#' Non-uniform asset grid (log-spaced, clustered near the borrowing constraint)
#'
#' Builds an \code{n}-point grid on \eqn{[a_{min}, a_{max}]} that is equidistant
#' in the logs of \eqn{a + \text{pivot}}, so gridpoints cluster near
#' \eqn{a_{min}}.  The first point is set exactly to \code{amin}.
#'
#' @param amax Numeric: upper bound of the asset grid.
#' @param n Integer >= 2: number of gridpoints.
#' @param amin Numeric: borrowing constraint / lower bound (default 0).
#' @param pivot Numeric > 0: curvature shift; smaller = more clustering near
#'   \code{amin}.  Default \code{abs(amin) + 0.25} matches sequence-jacobian.
#'
#' @return A strictly increasing numeric vector of length \code{n}, with
#'   \code{[1] == amin} and \code{[n] == amax}.
#'
#' @examples
#' a <- hank_asset_grid(amax = 200, n = 500, amin = 0)
#' a[1]      # 0
#' length(a) # 500
#' @export
hank_asset_grid <- function(amax, n, amin = 0, pivot = abs(amin) + 0.25) {
  n <- as.integer(n)
  if (n < 2L) stop("asset grid requires n >= 2")
  if (!is.finite(amax) || !is.finite(amin) || amax <= amin)
    stop("require amax > amin, both finite")
  if (!is.finite(pivot) || pivot <= 0) stop("pivot must be > 0")

  ## geomspace(amin+pivot, amax+pivot, n): equidistant in logs.
  lo <- amin + pivot
  hi <- amax + pivot
  a  <- exp(seq(log(lo), log(hi), length.out = n)) - pivot
  a[1L] <- amin   # enforce exact lower bound (undo the -pivot rounding)
  a[n]  <- amax   # enforce exact upper bound
  a
}
