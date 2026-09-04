## R/hank-finite-ramsey-ss.R
## --------------------------------------------------------------------------
## Global-EGM (LeGrand-Ragot-style) level-Ramsey steady state for the finite
## HANK capital-tax instrument, and the re-anchor/re-emit driver around it.
##
## Motivation (see R/hank-finite.R roxygen VALIDITY CAVEAT and
## tests/testthat/test-hank-finite-ramsey.R): the frozen-bracket .mod emission
## (`hank_finite_mod(instrument = "captax")`) is only a LOCAL model around the
## COMPETITIVE steady state. Unrestricted level-Ramsey via the augmented FOC
## on that emission is INVALID -- dW/dtau < 0 at tau = 0 (golden-rule
## under-accumulation) drives Newton to tau ~ -0.37, outside the frozen-
## bracket neighborhood, where Blanchard-Kahn fails.
##
## THE FIX implemented here: solve the Ramsey steady state directly in the
## full nonlinear EGM household problem (no interpolation-bracket freezing),
## using the KEY EQUIVALENCE that a linear capital tax `tau` with a
## balanced-budget lump-sum rebate `tr = tau*r*K` is EXACTLY a standard EGM
## household problem with a NET return `r_net = (1-tau)*r` and an income
## vector `y = w*e + tr` (tr is a scalar added to every income state, since
## it is lump-sum). This lets `hank_egm_solve()` be reused UNMODIFIED --
## R/hank-egm.R and R/hank-het-block.R are not touched by this file.
##
## Once the Ramsey (r, w, K, tau) triple is found, `hank_ks_coarse_anchored()`
## re-anchors a coarse grid at that point (still reusing the EXISTING
## anchoring machinery in R/hank-truncation.R -- see the small additive hook
## added there, `hank_ks_taxed_coarse_anchored()`, for the taxed-household
## case) and `hank_finite_mod(instrument = "captax", tau_rule = "tau=<tau*>;")`
## re-emits the frozen-bracket local model with brackets frozen AT THE RAMSEY
## SS instead of the competitive SS.
## --------------------------------------------------------------------------


#' Taxed household block at fixed (tau, r, w)
#'
#' Solves the household EGM problem under a linear capital tax \code{tau} with
#' balanced-budget lump-sum rebate \code{tr = tau*r*K}, using the KEY
#' EQUIVALENCE documented in the file header: this is EXACTLY
#' \code{\link{hank_egm_solve}} with net return \code{r_net = (1-tau)*r} and
#' income \code{y = w*e + tr}. Mirrors \code{\link{hank_het_block}}'s
#' packaging (policies, distribution, aggregates) but for the taxed problem;
#' does not modify \code{hank_het_block} itself.
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels.
#' @param beta,eis Discount factor and EIS.
#' @param tau Capital-income tax rate.
#' @param r,w Pre-tax return and wage.
#' @param K Aggregate capital (needed only to compute the lump-sum rebate
#'   \code{tr = tau*r*K}; pass the candidate market-clearing \code{K(r)}).
#' @param tol,maxit Passed to \code{\link{hank_egm_solve}}.
#'
#' @return A list with the same shape as \code{\link{hank_het_block}}
#'   (\code{a}, \code{c}, \code{Va}, \code{Lambda}, \code{D}, \code{A},
#'   \code{C}), plus \code{tau}, \code{tr}, \code{r} (pre-tax), \code{r_net}.
#' @keywords internal
.hank_taxed_het_block <- function(a_grid, Pi, e, beta, eis, tau, r, w, K,
                                  tol = 1e-11, maxit = 5000L) {
  tr <- tau * r * K
  r_net <- (1 - tau) * r
  y <- w * e + tr
  hh <- hank_egm_solve(a_grid, y = y, r = r_net, beta = beta, eis = eis,
                       Pi = Pi, tol = tol, maxit = maxit)
  if (!hh$converged)
    warning(".hank_taxed_het_block: household EGM did not converge")
  Lam <- hank_forward_operator(hh$a, a_grid, Pi)
  sd  <- hank_stationary_dist(Lam)
  D   <- sd$d
  list(a_grid = a_grid, Pi = Pi, e = e, beta = beta, eis = eis,
       tau = tau, tr = tr, r = r, r_net = r_net, w = w,
       a = hh$a, c = hh$c, Va = hh$Va,
       Lambda = Lam, D = D,
       A = hank_aggregate(D, hh$a),
       C = hank_aggregate(D, hh$c),
       n_e = length(e), n_a = length(a_grid),
       dist_converged = sd$converged)
}


#' Market-clearing steady state of the finite HANK at a fixed capital tax
#'
#' Given a candidate \code{tau}, clears the CAPITAL market of the taxed
#' economy: outer \code{uniroot} over the pre-tax return \code{r}; at each
#' candidate \code{r}, firm FOCs give \code{K(r)} and \code{w(r)}, the
#' balanced-budget transfer is \code{tr = tau*r*K(r)}, and the household
#' block is \code{\link{.hank_taxed_het_block}}. This mirrors
#' \code{\link{hank_ks_steady}} exactly except for the tax/transfer wedge in
#' the household problem; at \code{tau = 0} it reproduces
#' \code{hank_ks_steady} to machine precision (verified in the test file).
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels.
#' @param beta,eis Discount factor and EIS.
#' @param alpha,delta Capital share and depreciation.
#' @param Z Steady-state TFP (default 1).
#' @param tau Capital-income tax rate (scalar).
#' @param r_bracket Optional length-2 search bracket for \eqn{r^*}. Defaults
#'   to \code{c(-delta + 1e-4, 1/beta - 1 - 1e-4)} at \code{tau = 0} (same as
#'   \code{hank_ks_steady}), but is adjusted when \code{tau != 0} (see
#'   Details): the lower end is shrunk away from \code{-delta}, and the upper
#'   end is widened for \code{tau > 0} since it is the household's NET return
#'   \code{r_net = (1-tau)*r} (not the pre-tax \code{r}) that must stay below
#'   the patience threshold \code{1/beta}. When the default bracket fails to
#'   bracket a root (\code{uniroot}'s "not of opposite sign" or the taxed
#'   household's cash-on-hand grid going non-monotonic from an extreme
#'   \code{tr}), it is automatically widened/narrowed -- see Details.
#'
#' @details
#' The lower endpoint of the default bracket is \code{-delta/2} when
#' \code{tau != 0} (vs. \code{-delta + 1e-4} at \code{tau = 0}) -- still far
#' below any plausible equilibrium \code{r}, just not so close to
#' \code{-delta} that \code{K_of_r} explodes and drags the balanced-budget
#' rebate \code{tr = tau*r*K(r)} with it (the untaxed problem has
#' \code{tr = 0} identically and is immune to this). This economy can be
#' DYNAMICALLY INEFFICIENT at high risk aversion (competitive \code{r < 0},
#' i.e. capital above the golden rule) so the root can sit close to this
#' lower bound; if the initial default bracket does not straddle a root (or
#' hits the non-monotone-coh failure at an endpoint), the search retries with
#' the lower end moved halfway back toward \code{-delta} (up to 6 times, each
#' a genuinely WIDER interval) before giving up. The upper endpoint is
#' widened to \code{(1/beta - 1e-4)/(1-tau) - 1} for \code{tau > 0}: there the
#' pre-tax \code{r} can exceed \code{1/beta - 1} while \code{r_net} stays
#' sub-threshold, so the untaxed upper bound would otherwise exclude the true
#' root (observed empirically at \code{tau = 0.2}).
#'
#' @return A list of class \code{hank_ks_taxed} with \code{r}, \code{w},
#'   \code{K}, \code{tau}, \code{tr}, \code{Z}, calibration, the taxed
#'   \code{block} (see \code{.hank_taxed_het_block}), and \code{mkt_residual}.
#' @export
hank_ks_taxed_steady <- function(a_grid, Pi, e, beta, eis, alpha, delta,
                                 tau = 0, Z = 1, r_bracket = NULL) {
  K_of_r <- function(r) ((r + delta) / (alpha * Z))^(1 / (alpha - 1))
  w_of_r <- function(r) (1 - alpha) * Z * K_of_r(r)^alpha
  A_of_r <- function(r) {
    K <- K_of_r(r); w <- w_of_r(r)
    blk <- .hank_taxed_het_block(a_grid, Pi, e, beta = beta, eis = eis,
                                 tau = tau, r = r, w = w, K = K)
    blk$A
  }
  f <- function(r) A_of_r(r) - K_of_r(r)

  user_bracket <- !is.null(r_bracket)
  if (is.null(r_bracket)) {
    r_lo <- if (tau == 0) -delta + 1e-4 else -delta / 2
    r_hi <- if (tau > 0) (1 / beta - 1e-4) / (1 - tau) - 1
            else 1 / beta - 1 - 1e-4
    r_bracket <- c(r_lo, r_hi)
  }
  sol <- tryCatch(stats::uniroot(f, interval = r_bracket, tol = 1e-10),
                  error = function(e) e)
  if (inherits(sol, "error")) {
    if (user_bracket) stop(sol)
    ## Retry with a SEQUENCE of candidate lower endpoints between -delta and
    ## r_bracket[1] (both directions: the pathology is the extreme-|r| end
    ## of the DEFAULT bracket, either because K_of_r(r_lo) -> Inf drags the
    ## balanced-budget rebate tr = tau*r*K(r) with it (need r_lo CLOSER to 0),
    ## or because the default r_lo is simply not extreme enough to bracket a
    ## root that sits close to -delta (need r_lo CLOSER to -delta). Each
    ## candidate here is evaluated independently (not nested widening), so
    ## whichever direction is needed is found; only genuine sign-changing
    ## brackets are accepted.
    r_hi <- r_bracket[2L]
    cands <- unique(c(
      ## closer to 0 (narrower): fixes the extreme-tr blowup failure mode
      r_bracket[1L] * c(0.5, 0.25, 0.1, 0.05, 0.02, 0.01),
      ## closer to -delta (wider): fixes the no-sign-change failure mode
      -delta + (1 - c(0.75, 0.9, 0.97, 0.99, 0.999, 0.9999)) *
        (r_bracket[1L] - (-delta))))
    ok <- FALSE
    for (r_lo_try in cands) {
      sol <- tryCatch(stats::uniroot(f, interval = c(r_lo_try, r_hi),
                                     tol = 1e-10),
                      error = function(e) NULL)
      if (!is.null(sol)) { ok <- TRUE; break }
    }
    if (!ok)
      stop(sprintf(paste0(
        "hank_ks_taxed_steady(): uniroot() failed to bracket a root of the ",
        "taxed asset-market clearing condition at tau = %.6g after trying ",
        "%d candidate lower endpoints. Provide an explicit r_bracket."),
        tau, length(cands)))
  }
  r <- sol$root; w <- w_of_r(r); K <- K_of_r(r)
  blk <- .hank_taxed_het_block(a_grid, Pi, e, beta = beta, eis = eis,
                               tau = tau, r = r, w = w, K = K)
  structure(list(r = r, w = w, K = K, tau = tau, tr = blk$tr, Z = Z,
                 alpha = alpha, delta = delta, beta = beta, eis = eis,
                 block = blk, mkt_residual = blk$A - K),
            class = c("hank_ks_taxed", "hank_block"))
}


#' Steady-state utilitarian welfare of a taxed (or untaxed) equilibrium
#'
#' \eqn{W = \sum_i D_i\, u(c_i)} with CRRA utility (\code{gamma = 1/eis};
#' log utility if \code{eis} implies \code{gamma} within \code{1e-8} of 1),
#' evaluated at the cleared household consumption policy of a
#' \code{\link{hank_ks_taxed_steady}} (or \code{\link{hank_ks_steady}})
#' result.
#'
#' @param ks_taxed A \code{hank_ks_taxed} (or \code{hank_ks}) object.
#' @return Scalar steady-state utilitarian welfare.
#' @export
hank_ss_welfare <- function(ks_taxed) {
  blk <- ks_taxed$block
  eis <- blk$eis
  gam <- 1 / eis
  cc  <- .hank_mat_to_vec(blk$c)
  D   <- as.numeric(blk$D)
  u <- if (abs(gam - 1) < 1e-8) log(cc) else (cc^(1 - gam)) / (1 - gam)
  sum(D * u)
}


#' Global-EGM Ramsey capital-tax steady state
#'
#' Solves for the utilitarian-optimal LEVEL capital tax \code{tau*} directly
#' in the full nonlinear (global EGM) household problem -- i.e. WITHOUT the
#' frozen-interpolation-bracket approximation of \code{\link{hank_finite_mod}}
#' -- by root-finding the steady-state welfare gradient
#' \code{dW/dtau = 0} via central finite differences of
#' \code{\link{hank_ss_welfare}} evaluated through
#' \code{\link{hank_ks_taxed_steady}}.
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels.
#' @param beta,eis Discount factor and EIS.
#' @param alpha,delta Capital share and depreciation.
#' @param Z Steady-state TFP (default 1).
#' @param tau_bracket Search bracket for \code{tau*} (default \code{c(-0.2,
#'   0.2)}; widened automatically -- see Details -- if the sign of
#'   \code{dW/dtau} does not change across it).
#' @param fd_step Absolute finite-difference step used for the central
#'   difference of \code{dW/dtau} (default \code{2e-3}).
#' @param tol \code{uniroot} tolerance on \code{tau}.
#'
#' @details
#' If \code{dW/dtau} does not change sign on \code{tau_bracket}, the bracket
#' is widened geometrically (doubling the interval, up to 4 times) before
#' giving up -- this only ever produces a WIDER bracket, so any root found is
#' still the genuine welfare stationary point.
#'
#' @return A list of class \code{hank_ramsey_ss} with:
#'   \code{tau_ramsey} (the Ramsey capital tax), \code{ks_ramsey} (the
#'   \code{hank_ks_taxed_steady} result at \code{tau_ramsey}),
#'   \code{ks_competitive} (the \code{hank_ks_taxed_steady} result at
#'   \code{tau = 0}, i.e. the competitive benchmark -- exactly reproduces
#'   \code{\link{hank_ks_steady}}), \code{W_ramsey}, \code{W_competitive}
#'   (steady-state welfare at each), \code{dW_dtau_ramsey} (the FD derivative
#'   at \code{tau_ramsey}, should be near zero).
#' @export
hank_ramsey_ss_solve <- function(a_grid, Pi, e, beta, eis, alpha, delta,
                                 Z = 1, tau_bracket = c(-0.2, 0.2),
                                 fd_step = 2e-3, tol = 1e-8) {
  W_of_tau <- function(tau) {
    ks <- hank_ks_taxed_steady(a_grid, Pi, e, beta = beta, eis = eis,
                               alpha = alpha, delta = delta, tau = tau, Z = Z)
    hank_ss_welfare(ks)
  }
  dW_dtau <- function(tau) {
    (W_of_tau(tau + fd_step) - W_of_tau(tau - fd_step)) / (2 * fd_step)
  }

  lo <- tau_bracket[1L]; hi <- tau_bracket[2L]
  f_lo <- dW_dtau(lo); f_hi <- dW_dtau(hi)
  n_widen <- 0L
  while (sign(f_lo) == sign(f_hi) && n_widen < 4L) {
    width <- hi - lo
    lo <- lo - width / 2
    hi <- hi + width / 2
    f_lo <- dW_dtau(lo); f_hi <- dW_dtau(hi)
    n_widen <- n_widen + 1L
  }
  if (sign(f_lo) == sign(f_hi))
    stop(sprintf(paste0(
      "hank_ramsey_ss_solve(): dW/dtau does not change sign on the widened ",
      "bracket [%.4f, %.4f] (f_lo = %.4g, f_hi = %.4g). Provide a wider ",
      "tau_bracket."), lo, hi, f_lo, f_hi))

  sol <- stats::uniroot(dW_dtau, interval = c(lo, hi), f.lower = f_lo,
                        f.upper = f_hi, tol = tol)
  tau_ramsey <- sol$root

  ks_ramsey <- hank_ks_taxed_steady(a_grid, Pi, e, beta = beta, eis = eis,
                                    alpha = alpha, delta = delta,
                                    tau = tau_ramsey, Z = Z)
  ks_comp <- hank_ks_taxed_steady(a_grid, Pi, e, beta = beta, eis = eis,
                                  alpha = alpha, delta = delta,
                                  tau = 0, Z = Z)
  W_ramsey <- hank_ss_welfare(ks_ramsey)
  W_comp   <- hank_ss_welfare(ks_comp)

  structure(
    list(tau_ramsey = tau_ramsey, ks_ramsey = ks_ramsey,
         ks_competitive = ks_comp, W_ramsey = W_ramsey,
         W_competitive = W_comp, dW_dtau_ramsey = dW_dtau(tau_ramsey),
         a_grid = a_grid, Pi = Pi, e = e, beta = beta, eis = eis,
         alpha = alpha, delta = delta, Z = Z),
    class = c("hank_ramsey_ss", "hank_block"))
}


#' Re-anchor and re-emit the finite HANK captax model at the Ramsey SS
#'
#' Given a \code{\link{hank_ramsey_ss_solve}} result, builds the coarse-grid
#' anchored economy AT THE RAMSEY STEADY STATE (via
#' \code{\link{hank_ks_taxed_coarse_anchored}}, R/hank-truncation.R) and
#' re-emits the frozen-bracket finite HANK captax model
#' (\code{\link{hank_finite_mod}}) with \code{tau} PINNED at
#' \code{tau_ramsey} -- i.e. the frozen interpolation brackets/constraint set
#' now describe the model LOCALLY AROUND THE RAMSEY SS, not the competitive
#' SS, so a Newton/perturbation solve of this emission stays in-bracket where
#' the competitive-anchored emission would walk to the BK boundary (see
#' R/hank-finite.R's VALIDITY CAVEAT).
#'
#' @param ramsey A \code{\link{hank_ramsey_ss_solve}} result.
#' @param n_a Coarse grid size (default 8, matching
#'   \code{tests/testthat/test-hank-finite-ramsey.R}).
#' @param rho_z,sigma_z TFP calibration passed to \code{\link{hank_finite_mod}}.
#' @param order,path,verbose,sparse Passed to \code{\link{hank_finite_solve}}.
#'
#' @return A list with \code{ks_coarse} (the Ramsey-anchored coarse
#'   \code{hank_ks_taxed}/\code{hank_ks}), \code{solve} (the
#'   \code{\link{hank_finite_solve}} result, pinned at
#'   \code{tau = tau_ramsey}), and \code{tau_ramsey}.
#' @export
hank_ramsey_reanchor_emit <- function(ramsey, n_a = 8L, rho_z = 0.9,
                                      sigma_z = 0.01, order = 1L,
                                      path = NULL, verbose = FALSE,
                                      sparse = NULL) {
  if (!inherits(ramsey, "hank_ramsey_ss"))
    stop("hank_ramsey_reanchor_emit(): `ramsey` must be a hank_ramsey_ss object.")

  ks_coarse <- hank_ks_taxed_coarse_anchored(ramsey$ks_ramsey, n_a = n_a)

  tau_rule <- sprintf("tau=%.17g;", ramsey$tau_ramsey)
  sol <- hank_finite_solve(ks_coarse, order = order, rho_z = rho_z,
                           sigma_z = sigma_z, path = path, verbose = verbose,
                           instrument = "captax", tau_rule = tau_rule,
                           sparse = sparse)

  list(ks_coarse = ks_coarse, solve = sol, tau_ramsey = ramsey$tau_ramsey)
}
