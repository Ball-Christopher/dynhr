## R/hank-bond-block.R
## --------------------------------------------------------------------------
## The geometric (delta-coupon) bond: the package's DURATION PRIMITIVE and the
## cleanest asset-price revaluation instrument (briefs/19 sections 9.6-9.7).
##
## A geometric bond pays 1 this period, delta next period, delta^2 after, ...
## Under perfect foresight and a required-return path r_t, no-arbitrage prices
## it by the backward recursion
##     q_t = (1 + delta * q_{t+1}) / (1 + r_{t+1}),
## with steady state q_ss = 1 / (1 + r_ss - delta), and the REALIZED return on
## holdings carried into period t is
##     ra_t = (1 + delta * q_t) / q_{t-1} - 1        (q_0 = q_ss predetermined).
##
## Why this is the right first pricing block:
##   - the coupon stream is FIXED, so every movement in ra_t is pure
##     REVALUATION -- monetary policy moves the price with zero cash-flow
##     ambiguity, which is exactly the instrument the wealth-vs-cash-flow
##     decomposition needs;
##   - dq/dr has a CLOSED FORM (a geometric sequence in the horizon gap), so
##     the block's analytic Jacobian is testable against algebra, not just
##     finite differences;
##   - delta = 0 collapses to the one-period bond with ra_t = r_t exactly:
##     the revaluation channel is KILLED by a single parameter, giving a
##     built-in no-revaluation counterfactual;
##   - a delta-coupon LIABILITY is a mortgage-like long-duration debt, so the
##     same primitive carries the debt-side items D2/D6 in briefs/19 9.7.
## --------------------------------------------------------------------------


#' Geometric (delta-coupon) bond pricing block for the sequence-space DAG
#'
#' Builds a \code{\link{hank_simple_block}} that maps a required-return path
#' (input \code{r_input}) to the bond price path \code{q_name} and the
#' REALIZED holding return \code{ra_name}, via the perfect-foresight
#' no-arbitrage recursion \eqn{q_t = (1 + \delta q_{t+1})/(1 + r_{t+1})}
#' (terminal condition: return to steady state beyond the horizon) and
#' \eqn{ra_t = (1 + \delta q_t)/q_{t-1} - 1} with \eqn{q_0 = q_{ss}}
#' predetermined -- so \code{ra} at date 1 carries the IMPACT REVALUATION of
#' an unanticipated shock, which is what delivers stock-scaled incidence to
#' bond holders.
#'
#' The block's Jacobian is ANALYTIC. At steady state,
#' \deqn{\partial q_t / \partial r_s = -\frac{q_{ss}}{1+r_{ss}}
#'       \left(\frac{\delta}{1+r_{ss}}\right)^{s-t-1}, \quad s > t,}
#' zero otherwise, and \code{ra} chains through
#' \eqn{\partial ra_t/\partial q_t = \delta/q_{ss}},
#' \eqn{\partial ra_t/\partial q_{t-1} = -(1+r_{ss})/q_{ss}}.  The geometric
#' decay is asserted against the block in \code{test-hank-bond-block.R}, and
#' the \eqn{\delta = 0} limit must reproduce \eqn{ra_t = r_t} exactly (a
#' one-period bond has no revaluation) -- the built-in no-revaluation
#' counterfactual.
#'
#' Steady-state values the model's \code{ss} list must carry:
#' \code{q_name = 1/(1 + r_ss - delta)} and \code{ra_name = r_ss} (the
#' realized return equals the required return at steady state);
#' \code{\link{hank_bond_ss}} computes them.
#'
#' @param name Character block name.
#' @param delta Coupon decay in \code{[0, 1]}; effective (Macaulay-like)
#'   duration is \eqn{(1+r_{ss})/(1+r_{ss}-\delta)}.  Must satisfy
#'   \code{delta < 1 + r_ss} for a finite price.
#' @param r_ss Steady-state required return.
#' @param r_input Character: name of the required-return input variable in the
#'   DAG (e.g. \code{"rb"} to price off the liquid rate, or a spread-adjusted
#'   rate from an upstream block).
#' @param q_name,ra_name Character: output variable names.
#'
#' @return A \code{\link{hank_simple_block}} with inputs \code{r_input} and
#'   outputs \code{c(q_name, ra_name)}, carrying an analytic Jacobian.
#' @seealso \code{\link{hank_bond_ss}}, \code{\link{hank_td2_reval_decompose}},
#'   \code{\link{hank_simple_block}}
#' @examples
#' bb <- hank_bond_block(delta = 0.95, r_ss = 0.01)
#' hank_bond_ss(delta = 0.95, r_ss = 0.01)
#' @export
hank_bond_block <- function(name = "bond", delta, r_ss, r_input = "rb",
                            q_name = "q", ra_name = "ra") {
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta < 0 || delta > 1)
    stop("hank_bond_block: 'delta' must be a scalar in [0, 1].")
  if (!is.numeric(r_ss) || length(r_ss) != 1L || !is.finite(r_ss))
    stop("hank_bond_block: 'r_ss' must be a finite scalar.")
  if (1 + r_ss - delta <= 0)
    stop("hank_bond_block: need delta < 1 + r_ss for a finite price ",
         "(q_ss = 1/(1 + r_ss - delta)).")
  q_ss <- 1 / (1 + r_ss - delta)

  fn <- function(paths, ss) {
    r <- paths[[r_input]]
    T_h <- length(r)
    q <- numeric(T_h)
    ## backward from the return-to-steady-state terminal condition; q_t uses
    ## the NEXT period's required return (the buyer's holding period)
    q_next <- q_ss
    for (t in T_h:1L) {
      r_next <- if (t < T_h) r[t + 1L] else r_ss
      q[t] <- (1 + delta * q_next) / (1 + r_next)
      q_next <- q[t]
    }
    q_lag <- c(q_ss, q[-T_h])                # q_0 = q_ss predetermined
    ra <- (1 + delta * q) / q_lag - 1
    setNames(list(q, ra), c(q_name, ra_name))
  }

  jac <- function(ss, T_h) {
    ## dq_t/dr_s: geometric in the horizon gap, zero for s <= t
    dq <- matrix(0, T_h, T_h)
    base <- -q_ss / (1 + r_ss)
    fac  <- delta / (1 + r_ss)
    for (t in seq_len(T_h)) {
      s <- seq_len(T_h)
      k <- s - t - 1L
      dq[t, s > t] <- base * fac^k[s > t]
    }
    ## ra chains through q_t and q_{t-1} (row 1 has no q_{t-1} term: q_0 is
    ## predetermined at q_ss and does not respond)
    dra <- (delta / q_ss) * dq
    if (T_h >= 2L)
      dra[2:T_h, ] <- dra[2:T_h, ] -
        ((1 + r_ss) / q_ss) * dq[1:(T_h - 1L), , drop = FALSE]
    setNames(list(setNames(list(dq),  r_input),
                  setNames(list(dra), r_input)),
             c(q_name, ra_name))
  }

  hank_simple_block(name, inputs = r_input, outputs = c(q_name, ra_name),
                    fn = fn, jac = jac)
}


#' Steady-state values for a geometric bond
#'
#' @inheritParams hank_bond_block
#' @return List with \code{q} (\code{= 1/(1 + r_ss - delta)}), \code{ra}
#'   (\code{= r_ss}: realized equals required at steady state), and
#'   \code{duration} (\code{= (1 + r_ss)/(1 + r_ss - delta)}).
#' @seealso \code{\link{hank_bond_block}}
#' @export
hank_bond_ss <- function(delta, r_ss) {
  if (1 + r_ss - delta <= 0)
    stop("hank_bond_ss: need delta < 1 + r_ss.")
  list(q = 1 / (1 + r_ss - delta), ra = r_ss,
       duration = (1 + r_ss) / (1 + r_ss - delta))
}
