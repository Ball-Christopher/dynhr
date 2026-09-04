## R/hank-nk.R
## --------------------------------------------------------------------------
## One-asset New Keynesian HANK, assembled through the general block-DAG engine
## (R/hank-model.R).  A canonical demand-side HANK (in the spirit of McKay-
## Nakamura-Steinsson / the sequence-jacobian one-asset-HANK notebook, in a
## deliberately simplified reduced-form): heterogeneous households save in a
## single bond, monetary policy follows a Taylor rule, and inflation follows a
## reduced-form (output-gap) New Keynesian Phillips curve.
##
## Closure (levels; L = E[e] = 1, all output accrues to households as income):
##   Taylor:   i_t   = r_ss + phi * pi_t + rstar_t
##   Fisher:   r_t   = i_{t-1} - pi_t                 (real return on nominal bonds)
##   income:   w_t   = Y_t - r_t B                    (output net of the tax that
##                                                     funds bond interest r_t B)
##   household:(r, w) -> (A, C)                        (het block; income y_e = w e)
##   NKPC:     pi_t  = beta_p * pi_{t+1} + kappa (Y_t - Y_ss)
##   asset:    A_t   = B                               (bond-market clearing)
## Unknowns {Y, pi}; targets {nkpc_res, asset_mkt}; exogenous {rstar}. The tax
## r_t B funds the government's bond-interest bill, so the household budget gives
## C_t = Y_t and goods market (Y = C) clears by Walras given bond clearing. The
## steady-state bond supply solves the fixed point B = A(r_ss, Y_ss - r_ss B).
##
## This is a REPRESENTATIVE EXAMPLE model: its blocks are standard textbook
## relations and the whole object is validated for internal consistency (GE
## residuals ~0), linearization consistency (linear IRF is the O(scale) limit of
## the nonlinear transition), determinacy, and economic sanity -- not bit-matched
## to an external reference (published one-asset-HANK IRFs are plot-only).
## --------------------------------------------------------------------------


#' Build a one-asset New Keynesian HANK model
#'
#' Constructs the steady state (calibrating the bond supply to household asset
#' demand at the target real rate) and assembles the NK-HANK as a
#' \code{\link{hank_model}}.
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels.
#' @param beta,eis Household discount factor and elasticity of intertemporal
#'   substitution.
#' @param r_ss Target steady-state real rate (also the Taylor-rule intercept).
#' @param phi Taylor-rule inflation coefficient (\code{> 1} for determinacy).
#' @param kappa NKPC slope on the output gap.
#' @param beta_p NKPC discount factor (default \code{1/(1+r_ss)}).
#' @param T_h Integer horizon.
#'
#' @return A list of class \code{hank_nk} with the assembled \code{model}
#'   (a \code{\link{hank_model}}), the steady-state household \code{block}, the
#'   bond supply \code{B}, and the calibration.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(60, 60, 0)
#' nk  <- hank_nk_hank(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                     r_ss = 0.005, phi = 1.5, kappa = 0.1, T_h = 100)
#' irf <- hank_model_irf(nk$model, dZ = list(rstar = 0.01 * 0.6^(0:99)))
#' @export
hank_nk_hank <- function(a_grid, Pi, e, beta, eis, r_ss = 0.005,
                         phi = 1.5, kappa = 0.1, beta_p = NULL, T_h = 200L) {
  if (is.null(beta_p)) beta_p <- 1 / (1 + r_ss)
  Y_ss <- 1
  ## Steady-state bond supply solves B = A(r_ss, Y_ss - r_ss B): household asset
  ## demand at income net of the tax that funds the interest bill r_ss B.
  A_of_B <- function(B)
    hank_het_block(a_grid, Pi, e, beta = beta, eis = eis,
                   r = r_ss, w = Y_ss - r_ss * B)$A
  B <- stats::uniroot(function(B) A_of_B(B) - B,
                      interval = c(1e-6, Y_ss / r_ss - 1e-6), tol = 1e-10)$root
  w_ss  <- Y_ss - r_ss * B
  block <- hank_het_block(a_grid, Pi, e, beta = beta, eis = eis,
                          r = r_ss, w = w_ss)

  taylor <- hank_simple_block(
    "taylor", inputs = c("pi", "rstar"), outputs = "i",
    fn = function(paths, ss) list(i = r_ss + phi * paths$pi + paths$rstar))

  fisher <- hank_simple_block(
    "fisher", inputs = c("i", "pi"), outputs = "r",
    fn = function(paths, ss) {
      i_lag <- c(ss$i, paths$i[-length(paths$i)])
      list(r = i_lag - paths$pi)
    })

  income <- hank_simple_block(
    "income", inputs = c("Y", "r"), outputs = "w",
    fn = function(paths, ss) list(w = paths$Y - paths$r * B))

  household <- hank_het_block_spec("household", block,
                                  inputs = c("r", "w"), outputs = c("A", "C"))

  nkpc <- hank_simple_block(
    "nkpc", inputs = c("pi", "Y"), outputs = "nkpc_res",
    fn = function(paths, ss) {
      pi_lead <- c(paths$pi[-1], ss$pi)
      list(nkpc_res = paths$pi - beta_p * pi_lead - kappa * (paths$Y - ss$Y))
    })

  asset <- hank_simple_block(
    "asset", inputs = "A", outputs = "asset_mkt",
    fn = function(paths, ss) list(asset_mkt = paths$A - B))

  ss <- list(Y = Y_ss, pi = 0, i = r_ss, r = r_ss, w = w_ss,
             A = B, C = block$C, nkpc_res = 0, asset_mkt = 0, rstar = 0)
  model <- hank_model(list(taylor, fisher, income, household, nkpc, asset),
                      unknowns = c("Y", "pi"),
                      targets = c("asset_mkt", "nkpc_res"),
                      exogenous = "rstar", ss = ss, T_h = T_h)

  structure(list(model = model, block = block, B = B, Y_ss = Y_ss,
                 r_ss = r_ss, phi = phi, kappa = kappa, beta_p = beta_p),
            class = c("hank_nk", "hank_block"))
}


#' Build a one-asset New Keynesian HANK with a discount-factor MIXTURE household
#'
#' Mirrors \code{\link{hank_nk_hank}} exactly (same Taylor rule, Fisher
#' equation, income/tax closure, NKPC, and asset-market target), replacing the
#' single-\code{beta} household with a \code{K}-type discount-factor MIXTURE
#' household (see \code{\link{hank_mixture_blocks}} /
#' \code{\link{hank_mixture_block_spec}}): every type shares the asset grid,
#' income process, and EIS, differing only in \code{beta}, and interacts with
#' the rest of the economy only through the common aggregate prices
#' \code{(r, w)}. Because of that, the household block's steady-state
#' aggregates and sequence-space Jacobian are EXACT omega-weighted sums of the
#' per-type objects -- no cross term between types -- so this constructor lets
#' the level-vs-response identification hierarchy (see
#' \code{\link{hank_partial_id_level_response}}) be tested inside a full sticky-price GE
#' model instead of only in partial equilibrium.
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels (shared
#'   by every type).
#' @param betas Numeric length-K vector of discount factors, one per type.
#' @param omega Numeric length-K mixture weights, non-negative, summing to 1.
#' @param eis Elasticity of intertemporal substitution, shared (scalar; only
#'   \code{beta} varies across types here).
#' @param r_ss Target steady-state real rate (also the Taylor-rule intercept).
#' @param phi Taylor-rule inflation coefficient (\code{> 1} for determinacy).
#' @param kappa NKPC slope on the output gap.
#' @param beta_p NKPC discount factor (default \code{1/(1+r_ss)}).
#' @param T_h Integer horizon.
#'
#' @return A list of class \code{c("hank_nk_mixture", "hank_nk")} with the
#'   assembled \code{model} (a \code{\link{hank_model}}), the per-type
#'   steady-state \code{blocks} (list of \code{K}
#'   \code{\link{hank_het_block}} objects, all solved at the clearing
#'   \code{(r_ss, w_ss)}), the mixture weights \code{omega}, the discount
#'   factors \code{betas}, the bond supply \code{B}, and the calibration.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(60, 60, 0)
#' nkm <- hank_nk_hank_mixture(ag, inc$Pi, inc$e,
#'                             betas = c(0.94, 0.97), omega = c(0.5, 0.5),
#'                             eis = 1, r_ss = 0.005, phi = 1.5, kappa = 0.1,
#'                             T_h = 100)
#' ## Steady-state consistency: household budget + bond clearing => goods
#' ## clearing by Walras, so C_ss = Y_ss exactly (as in hank_nk_hank).
#' irf <- hank_model_irf(nkm$model, dZ = list(rstar = 0.01 * 0.6^(0:99)))
#' @export
hank_nk_hank_mixture <- function(a_grid, Pi, e, betas, omega, eis,
                                 r_ss = 0.005, phi = 1.5, kappa = 0.1,
                                 beta_p = NULL, T_h = 200L) {
  if (length(omega) != length(betas))
    stop(sprintf(
      "hank_nk_hank_mixture(): length(omega) (%d) must equal length(betas) (%d).",
      length(omega), length(betas)))
  if (any(!is.finite(omega)) || any(omega < 0))
    stop("hank_nk_hank_mixture(): 'omega' must be finite and non-negative.")
  if (abs(sum(omega) - 1) > 1e-8)
    stop(sprintf(
      "hank_nk_hank_mixture(): 'omega' must sum to 1 (sum = %.6f).", sum(omega)))
  if (max(betas) * (1 + r_ss) >= 1)
    stop(sprintf(paste0(
      "hank_nk_hank_mixture(): max(betas)*(1+r_ss) = %.6f >= 1 -- the most-",
      "patient type would have unbounded asset demand at r_ss = %.6f. Lower ",
      "max(betas) or r_ss."), max(betas) * (1 + r_ss), r_ss))

  if (is.null(beta_p)) beta_p <- 1 / (1 + r_ss)
  Y_ss <- 1
  ## Steady-state bond supply solves the MIXTURE fixed point
  ## B = Sum_k omega_k * A_k(r_ss, Y_ss - r_ss B): each type's asset demand is
  ## evaluated at the SAME (r_ss, w) since types share prices; only the
  ## aggregation (the omega-weighted sum) differs from hank_nk_hank's single
  ## type.
  A_of_B <- function(B) {
    w <- Y_ss - r_ss * B
    acc <- 0
    for (k in seq_along(betas))
      acc <- acc + omega[k] * hank_het_block(a_grid, Pi, e, beta = betas[k],
                                             eis = eis, r = r_ss, w = w)$A
    acc
  }
  B <- stats::uniroot(function(B) A_of_B(B) - B,
                      interval = c(1e-6, Y_ss / r_ss - 1e-6), tol = 1e-10)$root
  w_ss   <- Y_ss - r_ss * B
  blocks <- hank_mixture_blocks(a_grid, Pi, e, betas = betas, eis = eis,
                                r = r_ss, w = w_ss)

  taylor <- hank_simple_block(
    "taylor", inputs = c("pi", "rstar"), outputs = "i",
    fn = function(paths, ss) list(i = r_ss + phi * paths$pi + paths$rstar))

  fisher <- hank_simple_block(
    "fisher", inputs = c("i", "pi"), outputs = "r",
    fn = function(paths, ss) {
      i_lag <- c(ss$i, paths$i[-length(paths$i)])
      list(r = i_lag - paths$pi)
    })

  income <- hank_simple_block(
    "income", inputs = c("Y", "r"), outputs = "w",
    fn = function(paths, ss) list(w = paths$Y - paths$r * B))

  household <- hank_mixture_block_spec("household", blocks, omega,
                                       inputs = c("r", "w"), outputs = c("A", "C"))

  nkpc <- hank_simple_block(
    "nkpc", inputs = c("pi", "Y"), outputs = "nkpc_res",
    fn = function(paths, ss) {
      pi_lead <- c(paths$pi[-1], ss$pi)
      list(nkpc_res = paths$pi - beta_p * pi_lead - kappa * (paths$Y - ss$Y))
    })

  asset <- hank_simple_block(
    "asset", inputs = "A", outputs = "asset_mkt",
    fn = function(paths, ss) list(asset_mkt = paths$A - B))

  C_ss <- sum(vapply(seq_along(blocks), function(k) omega[k] * blocks[[k]]$C,
                     numeric(1)))
  ss <- list(Y = Y_ss, pi = 0, i = r_ss, r = r_ss, w = w_ss,
             A = B, C = C_ss, nkpc_res = 0, asset_mkt = 0, rstar = 0)
  model <- hank_model(list(taylor, fisher, income, household, nkpc, asset),
                      unknowns = c("Y", "pi"),
                      targets = c("asset_mkt", "nkpc_res"),
                      exogenous = "rstar", ss = ss, T_h = T_h)

  structure(list(model = model, blocks = blocks, omega = omega, betas = betas,
                 B = B, Y_ss = Y_ss, r_ss = r_ss, phi = phi, kappa = kappa,
                 beta_p = beta_p),
            class = c("hank_nk_mixture", "hank_nk", "hank_block"))
}
