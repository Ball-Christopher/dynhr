## R/hank-debt-blocks.R
## --------------------------------------------------------------------------
## Debt-side primitives D2 and D4 (briefs/19 section 9.7): staggered rate
## repricing and the Fisher channel. Both are exact-algebra simple blocks in
## the DAG plus, for Fisher, a date-1 stock-revaluation helper.
##
## The intended NZ composition (with D3's r_minus household input):
##     rb (policy-linked) --> hank_repricing_block --> rb_eff
##                            rb_eff + wedge --> r_minus --> household
## so the household's BORROWING cost reprices with the ~2-year NZ mortgage-fix
## lag while its deposit rate tracks the market -- the cash-flow channel of
## the criticised literature, done properly and separately named.
## --------------------------------------------------------------------------


#' Staggered rate-repricing block (the NZ mortgage-fix structure, D2)
#'
#' The effective rate on an outstanding debt stock when a fraction
#' \code{phi_r} of loans reprices to the market rate each period:
#' \deqn{\tilde r_t = (1 - \phi_r)\,\tilde r_{t-1} + \phi_r\, r_t,}
#' with the pre-shock stock rate \eqn{\tilde r_0 = r_{ss}} predetermined.
#' The average repricing lag is \eqn{1/\phi_r} periods (NZ mortgage fixes:
#' roughly two years, \code{phi_r} about 0.125 quarterly).
#'
#' The Jacobian is analytic and exactly geometric:
#' \eqn{\partial \tilde r_t / \partial r_s = \phi_r (1-\phi_r)^{t-s}} for
#' \eqn{s \le t}, zero otherwise.  \code{phi_r = 1} is instant repricing
#' (\eqn{\tilde r = r} identically) and \code{phi_r -> 0} freezes the stock
#' rate -- the two limits pin the block from both ends.
#'
#' @param name Character block name.
#' @param phi_r Repricing fraction per period, in \code{(0, 1]}.
#' @param r_ss Steady-state rate (the predetermined date-0 stock rate).
#' @param r_input Character: name of the market-rate input variable.
#' @param r_eff_name Character: name of the effective-rate output variable.
#' @return A \code{\link{hank_simple_block}} with an analytic Jacobian.
#' @seealso \code{\link{hank_bond_block}} (the duration primitive),
#'   \code{\link{hank_fisher_block}}
#' @export
hank_repricing_block <- function(name = "repricing", phi_r, r_ss,
                                 r_input = "rb", r_eff_name = "rb_eff") {
  if (!is.numeric(phi_r) || length(phi_r) != 1L || !is.finite(phi_r) ||
      phi_r <= 0 || phi_r > 1)
    stop("hank_repricing_block: 'phi_r' must be a scalar in (0, 1].")
  if (!is.numeric(r_ss) || length(r_ss) != 1L || !is.finite(r_ss))
    stop("hank_repricing_block: 'r_ss' must be a finite scalar.")
  fn <- function(paths, ss) {
    r <- paths[[r_input]]
    r_eff <- numeric(length(r))
    prev <- r_ss                                # predetermined stock rate
    for (t in seq_along(r)) {
      r_eff[t] <- (1 - phi_r) * prev + phi_r * r[t]
      prev <- r_eff[t]
    }
    setNames(list(r_eff), r_eff_name)
  }
  jac <- function(ss, T_h) {
    Jm <- matrix(0, T_h, T_h)
    for (t in seq_len(T_h)) {
      s <- seq_len(t)
      Jm[t, s] <- phi_r * (1 - phi_r)^(t - s)
    }
    setNames(list(setNames(list(Jm), r_input)), r_eff_name)
  }
  hank_simple_block(name, inputs = r_input, outputs = r_eff_name,
                    fn = fn, jac = jac)
}


#' Fisher-equation block: the real rate from nominal rate and inflation (D4)
#'
#' \deqn{r_t = (1 + i_t)/(1 + \pi_t) - 1.}
#' Under perfect foresight this is ALL of the Fisher channel for dates
#' \eqn{t \ge 2}: anticipated inflation is just a real-rate movement.  The
#' unanticipated date-1 part -- the revaluation of the PREDETERMINED nominal
#' stock by the new price level -- is a distribution shock, handled by
#' \code{\link{hank_liquid_reval_d0}} and fed through the transition's
#' \code{D0} argument.  Keeping the two pieces separate is deliberate: they
#' are different channels with different incidence (the anticipated part
#' works through choices; the surprise part is a pure balance-sheet
#' redistribution from lenders to borrowers).
#'
#' @param name Character block name.
#' @param i_ss,pi_ss Steady-state nominal rate and inflation.
#' @param i_input,pi_input,r_name Variable names in the DAG.
#' @return A \code{\link{hank_simple_block}} with an analytic (diagonal)
#'   Jacobian.
#' @seealso \code{\link{hank_liquid_reval_d0}}, \code{\link{hank_repricing_block}}
#' @export
hank_fisher_block <- function(name = "fisher", i_ss, pi_ss,
                              i_input = "i", pi_input = "pi_infl",
                              r_name = "rb") {
  for (x in list(i_ss, pi_ss))
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
      stop("hank_fisher_block: 'i_ss' and 'pi_ss' must be finite scalars.")
  if (1 + pi_ss <= 0)
    stop("hank_fisher_block: need 1 + pi_ss > 0.")
  fn <- function(paths, ss) {
    setNames(list((1 + paths[[i_input]]) / (1 + paths[[pi_input]]) - 1),
             r_name)
  }
  jac <- function(ss, T_h) {
    I <- diag(T_h)
    setNames(list(setNames(list(I / (1 + pi_ss),
                                -I * (1 + i_ss) / (1 + pi_ss)^2),
                           c(i_input, pi_input))), r_name)
  }
  hank_simple_block(name, inputs = c(i_input, pi_input), outputs = r_name,
                    fn = fn, jac = jac)
}


#' Date-1 revaluation of nominal liquid positions (the Fisher surprise, D4)
#'
#' An unanticipated price-level jump deflates every NOMINAL position on the
#' liquid axis by the same factor: real \eqn{b \to \lambda b} with
#' \eqn{\lambda = (1+\pi_{ss})/(1+\pi_1)}.  Savers lose, borrowers gain (debt
#' shrinks toward zero) -- the Fisher redistribution.  \strong{The AGGREGATE
#' consumption sign is not unconditional}: it is governed by the net nominal
#' position and the covariance of MPCs with positions (Auclert's
#' redistribution channel).  Measured on a net-SAVER test block, surprise
#' inflation LOWERS aggregate consumption (the wealth effect on positive net
#' positions dominates the borrower gain); for a net-DEBTOR household sector
#' -- the NZ configuration, where mortgages exceed deposits -- the sign
#' flips.  Condition any aggregate claim on the calibrated net position; the
#' per-group incidence is what is robust.  This helper remaps the block's
#' stationary distribution accordingly (Young's lottery on the scaled
#' positions, per income state and -- for two-asset blocks -- per illiquid
#' slice, which stays put: housing is real), for use as
#' \code{hank_td_nonlinear(..., D0 = )} / \code{hank_td2_nonlinear(..., D0 = )}.
#'
#' Positions scaled beyond the grid are clamped by the lottery (mass cannot
#' leave the grid); with the F13 grid discipline the clamped mass is
#' negligible.
#'
#' @param block A \code{\link{hank_het_block}}, \code{\link{hank_het2_block}}
#'   or \code{\link{hank_het2d_block}}.
#' @param factor The revaluation factor \eqn{\lambda > 0} applied to the
#'   liquid axis (\code{< 1} for surprise inflation).
#' @return A distribution vector in the block's cell order, suitable for the
#'   corresponding transition's \code{D0}.
#' @seealso \code{\link{hank_fisher_block}}
#' @export
hank_liquid_reval_d0 <- function(block, factor) {
  if (!is.numeric(factor) || length(factor) != 1L || !is.finite(factor) ||
      factor <= 0)
    stop("hank_liquid_reval_d0: 'factor' must be a finite scalar > 0.")
  if (inherits(block, "hank_het_block")) {
    ## one-asset: the whole asset axis is the (nominal) liquid position
    n_e <- block$n_e; n_a <- block$n_a
    ## block$D is in cell order (asset fastest): reshape to n_e x n_a
    Dm <- t(matrix(block$D, n_a, n_e))
    pos <- matrix(block$a_grid, n_e, n_a, byrow = TRUE) * factor
    lot <- .hank_lottery(pos, block$a_grid)
    out <- matrix(0, n_e, n_a)
    for (e in seq_len(n_e)) for (j in seq_len(n_a)) {
      i0 <- lot$i[e, j]; p0 <- lot$p[e, j]
      out[e, i0]      <- out[e, i0]      + p0 * Dm[e, j]
      out[e, i0 + 1L] <- out[e, i0 + 1L] + (1 - p0) * Dm[e, j]
    }
    return(as.numeric(t(out)))
  }
  if (inherits(block, "hank_het2_block") || inherits(block, "hank_het2d_block")) {
    ## two-asset: scale the LIQUID (b) axis only; the illiquid axis is real
    n_e <- block$n_e; n_b <- block$n_b; n_a <- block$n_a
    Darr <- .hank2_vec_to_arr(block$D, n_e, n_b, n_a)
    pos <- block$b_grid * factor
    lot <- .hank_lottery(matrix(pos, 1L), block$b_grid)
    out <- array(0, c(n_e, n_b, n_a))
    for (j in seq_len(n_b)) {
      i0 <- lot$i[1L, j]; p0 <- lot$p[1L, j]
      out[, i0, ]      <- out[, i0, ]      + p0 * Darr[, j, ]
      out[, i0 + 1L, ] <- out[, i0 + 1L, ] + (1 - p0) * Darr[, j, ]
    }
    return(.hank2_arr_to_vec(out))
  }
  stop("hank_liquid_reval_d0: 'block' must be a hank_het_block, ",
       "hank_het2_block or hank_het2d_block.")
}
