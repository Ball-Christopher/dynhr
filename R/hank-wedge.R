## R/hank-wedge.R
## --------------------------------------------------------------------------
## D3 (briefs/19 section 9.7): the BORROWING WEDGE on the one-asset household.
## The liquid return becomes state-dependent -- r_plus on savings (b >= 0),
## r_minus on debt (b < 0) -- and r_minus is a SECOND AGGREGATE INPUT, not a
## parameter, so the DAG can compose the NZ chain
##     rb --> hank_repricing_block --> rb_eff (+ wedge) --> r_minus.
##
## PARALLEL MACHINERY, VALIDATED PATHS UNTOUCHED: the symmetric
## .hank_egm_step / hank_egm_solve (R + cpp) are not modified; blocks built
## without r_minus are byte-identical to before. The wedge variant is R-only
## (cpp kernel a documented follow-up) and is dispatched from the block layer
## when r_minus is supplied.
##
## THE KINK. coh(b) = (1 + r(b)) b + y is continuous with a slope drop at
## b = 0, and the envelope Vb(b) = (1 + r(b)) u'(c(b)) JUMPS at b = 0, so the
## Euler equation has a gap there: households whose Wb falls inside the gap
## optimally park at EXACTLY b' = 0 ("neither borrow nor save"), producing a
## point mass at zero -- the classic wedge signature, and a TESTED qualitative
## oracle. On the grid the jump is smeared across the cells adjacent to zero
## (one grid gap of smoothing, shrinking with refinement); include 0 as a
## grid point so the smearing is one-sided. The endogenous grid stays
## MONOTONE under the wedge (both branches of the piecewise inversion are
## increasing and meet continuously at x = y), so the symmetric solver's
## interpolation machinery carries over unchanged.
## --------------------------------------------------------------------------


#' One EGM backward step with a borrowing wedge (one-asset)
#'
#' As \code{\link{.hank_egm_step}}, with the liquid return \code{r_plus} on
#' \code{b >= 0} and \code{r_minus} on \code{b < 0}.  See the file header for
#' the kink treatment.
#' @keywords internal
.hank_egm_step_wedge <- function(Va_p, a_grid, y, r_plus, r_minus,
                                 beta, eis, Pi, amin = a_grid[1L],
                                 coh_extra = NULL) {
  n_e <- length(y); n_a <- length(a_grid)
  tiny <- 1e-12
  Wb <- pmax(beta * (Pi %*% Va_p), tiny)         # over choices a'
  c_endo <- Wb^(-eis)
  x <- c_endo + matrix(a_grid, n_e, n_a, byrow = TRUE)   # resources needed
  ## piecewise inversion of the kinked coh: b >= 0 iff x - y >= 0
  ## (coh_extra is a beginning-of-period (e,a) quantity -- like the symmetric
  ## kernel, it enters only the ACTUAL coh below, never this endogenous-grid
  ## inversion, which is indexed by the CHOICE a', not the current state.)
  gap <- x - y
  a_endo <- ifelse(gap >= 0, gap / (1 + r_plus), gap / (1 + r_minus))
  r_state <- ifelse(a_grid < 0, r_minus, r_plus)          # today's rate by state
  coh <- matrix((1 + r_state) * a_grid, n_e, n_a, byrow = TRUE) + y
  if (!is.null(coh_extra)) coh <- coh + coh_extra
  a_pol <- matrix(0, n_e, n_a); c_pol <- a_pol
  for (e in seq_len(n_e)) {
    ap <- .hank_interp1(a_endo[e, ], a_grid, a_grid)
    ap <- pmax(ap, amin)
    a_pol[e, ] <- ap
    c_pol[e, ] <- coh[e, ] - ap                            # exact budget
  }
  c_pol_raw <- c_pol
  if (!is.null(coh_extra) && any(c_pol_raw <= tiny))
    .hank_egm_warn_once(
      "coh_extra_floor_wedge",
      "hank_egm (wedge): a transition-path evaluation with a matrix ",
      "'Tr_incidence' (Tier 2) drove cash-on-hand to the EGM tiny-floor ",
      "(1e-12) region at some (e, a) point. This warning fires at most ",
      "once per session.")
  c_pol <- pmax(c_pol_raw, tiny)
  uc <- c_pol^(-1 / eis)
  Va <- matrix(1 + r_state, n_e, n_a, byrow = TRUE) * uc   # envelope, state rate
  list(Va = Va, a = a_pol, c = c_pol)
}


#' Solve the one-asset household with a borrowing wedge to a stationary policy
#'
#' Fixed point of \code{\link{.hank_egm_step_wedge}}.  Validation and the
#' starting point come from a SYMMETRIC \code{\link{hank_egm_solve}} at
#' \code{r = r_plus} (which also makes the \code{r_minus = r_plus} case an
#' exact-reduction oracle: the wedge solver must then reproduce the validated
#' symmetric solution).
#' @keywords internal
.hank_egm_solve_wedge <- function(a_grid, y, r_plus, r_minus, beta, eis, Pi,
                                  amin = a_grid[1L], tol = 1e-11,
                                  maxit = 5000L, Va_init = NULL,
                                  coh_extra = NULL) {
  if (!is.numeric(r_minus) || length(r_minus) != 1L || !is.finite(r_minus))
    stop(".hank_egm_solve_wedge: 'r_minus' must be a finite scalar.")
  base <- hank_egm_solve(a_grid, y = y, r = r_plus, beta = beta, eis = eis,
                         Pi = Pi, amin = amin, tol = tol, maxit = maxit,
                         Va_init = Va_init, coh_extra = coh_extra)
  Va <- base$Va
  converged <- FALSE; it <- 0L; a_old <- NULL; step <- NULL
  for (it in seq_len(as.integer(maxit))) {
    step <- .hank_egm_step_wedge(Va, a_grid, y, r_plus, r_minus, beta, eis,
                                 Pi, amin = amin, coh_extra = coh_extra)
    Va <- step$Va
    if (!is.null(a_old)) {
      d <- max(abs(step$a - a_old))
      if (!is.finite(d)) break
      if (d < tol) { converged <- TRUE; break }
    }
    a_old <- step$a
  }
  c(step, list(iterations = it, converged = converged, a_grid = a_grid,
               y = y, r = r_plus, r_minus = r_minus, beta = beta, eis = eis,
               Pi = Pi, amin = amin))
}


#' Reject a wedge block at a wedge-unaware entry point, loudly
#'
#' Several one-asset routines difference or invert the household problem with
#' a SINGLE rate (\code{block$r}): \code{hank_impc}, \code{hank_mpc},
#' \code{hank_euler_residual}, \code{hank_sam_reiter_linearize},
#' \code{hank_het_dist_jacobian}. On a borrowing-wedge block their answers
#' would be wrong precisely on the debt side -- the side the wedge exists to
#' price -- so they refuse rather than compute plausible nonsense (the F19
#' discipline).
#' @keywords internal
.hank_reject_wedge <- function(block, caller) {
  if (is.null(block$r_minus)) return(invisible(NULL))
  stop(caller, "(): this block carries a borrowing wedge (r_minus = ",
       format(block$r_minus), " != r = ", format(block$r), "), and this ",
       "routine is not wedge-aware -- it would price the DEBT side at the ",
       "saving rate. Use the wedge-aware paths (hank_td_nonlinear, ",
       "hank_het_jacobian) or rebuild the block without 'r_minus'.",
       call. = FALSE)
}
