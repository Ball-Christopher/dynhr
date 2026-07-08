## R/hank-welfare.R
## --------------------------------------------------------------------------
## Welfare machinery (W0) for a one-asset HANK het block: the LEVEL value
## function and consumption-equivalent variation (CEV).
##
## Het blocks (see R/hank-het-block.R, R/hank-egm.R) store the steady-state
## policies and the MARGINAL value Va = u'(c) (needed for the EGM Euler-
## equation inversion), but not the level value V itself.  Welfare comparisons
## (e.g. across calibrations, or of a policy counterfactual) need the level:
##
##   V(x) = u(c(x)) + beta * E_x[V(x')]
##
## which at the steady-state policy is a LINEAR fixed point in V, because the
## policy (and hence the transition) is fixed and the aggregate-price inputs
## (r, w) do not appear directly except through the already-solved c and
## Lambda:
##
##   V = u(c) + beta * (Lambda %*% V)
##
## Objects (all in the block's CELL / distribution order, i.e. index
## (e-1)*n_a + a, asset-fast; see R/hank-distribution.R):
##   - c      : n_e x n_a steady-state consumption policy matrix (block$c);
##              flattened to cell order via .hank_mat_to_vec() before use.
##   - Lambda : (n_e*n_a) x (n_e*n_a) sparse EXPECTATION operator (no
##              transpose): E_x[f(x')] = (Lambda %*% f)(x).  This is the same
##              convention hank_het_jacobian() uses for its expectation
##              vectors (R/hank-jacobian.R): "E_s = Lam %*% E_{s-1} with NO
##              transpose", as opposed to distributions, which push FORWARD
##              via t(Lambda) %*% d.
##   - V      : length-(n_e*n_a) LEVEL value vector, cell order matching
##              block$D.  Solved here; not stored on the block by
##              hank_het_block().
##
## CEV convention: hank_cev() answers "what constant proportional scaling of
## the consumption path makes baseline lifetime utility V0 equal to
## alternative lifetime utility V1?".  This is implemented for GENERAL CRRA
## (see hank_cev() below) via the closed-form mapping
##   V0 -> (1+lambda)^rho (V0 + const) - const,   const = 1/(rho (1-beta)),
## with rho = 1 - 1/eis, inverted for lambda.  The eis == 1 (log-utility) case
## is handled as a special additive-shift branch (the rho -> 0 limit of the
## above), since it requires dividing by rho = 0.
## --------------------------------------------------------------------------


#' CRRA flow utility (vectorized)
#'
#' \eqn{u(c) = \log(c)} if \code{eis == 1}, else the "-1"-normalized CRRA form
#' \eqn{u(c) = (c^{1 - 1/eis} - 1) / (1 - 1/eis)}.  The normalization (the
#' \code{-1} in the numerator) is chosen so the \code{eis != 1} expression has
#' a finite limit as \code{eis -> 1} (it converges to \eqn{\log(c)} by
#' L'Hopital, since both numerator and denominator vanish there), keeping the
#' family continuous in \code{eis} even though only the \code{eis == 1} branch
#' is evaluated exactly by this function.
#'
#' @param c Numeric vector: consumption levels.  Must be strictly positive.
#' @param eis Numeric scalar: elasticity of intertemporal substitution.
#' @return Numeric vector, same length as \code{c}.
#' @export
hank_utility <- function(c, eis) {
  if (any(!is.finite(c)) || any(c <= 0))
    stop("hank_utility(): c must be strictly positive and finite")
  if (isTRUE(all.equal(eis, 1))) {
    log(c)
  } else {
    (c^(1 - 1 / eis) - 1) / (1 - 1 / eis)
  }
}


#' Level value function of a het block's steady-state policy
#'
#' Solves the linear fixed point \eqn{V = u(c) + \beta (\Lambda V)} for the
#' LEVEL value \code{V} (as opposed to the marginal value \code{block$Va}
#' already stored on the block), by fixed-point iteration from
#' \eqn{V_0 = u(c) / (1 - \beta)} (the value of consuming \code{c} forever).
#' \code{Lambda} is row-stochastic (spectral radius 1) and \code{beta < 1}, so
#' the iteration is a contraction and converges geometrically at rate
#' \code{beta}.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param tol Numeric: convergence tolerance on \code{max|V_{k+1} - V_k|}.
#' @param maxit Integer: maximum iterations.
#'
#' @return Numeric length-\code{n_e*n_a} vector: the level value \code{V},
#'   in the same cell order as \code{block$D}.
#' @export
hank_value_function <- function(block, tol = 1e-11, maxit = 100000L) {
  c_vec <- .hank_mat_to_vec(block$c)
  u_c   <- hank_utility(c_vec, block$eis)
  beta  <- block$beta
  Lam   <- block$Lambda

  V <- u_c / (1 - beta)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    V_new <- u_c + beta * as.numeric(Lam %*% V)
    if (max(abs(V_new - V)) < tol) { V <- V_new; converged <- TRUE; break }
    V <- V_new
  }
  if (!converged)
    warning("hank_value_function: Bellman fixed-point iteration did not ",
            "converge in ", maxit, " iterations")
  V
}


#' Per-cell consumption-equivalent variation between two value functions
#'
#' The constant proportional consumption scaling \eqn{\lambda} such that
#' scaling the baseline consumption path by \eqn{(1+\lambda)} in every period
#' and state raises baseline lifetime welfare \code{V0} to the alternative
#' lifetime welfare \code{V1}.
#'
#' Implemented for GENERAL CRRA utility.  With the "-1"-normalized CRRA flow
#' utility \eqn{u(c) = (c^\rho - 1)/\rho}, \eqn{\rho = 1 - 1/eis}, scaling the
#' entire consumption path by \eqn{(1+\lambda)} maps the value function as
#' \deqn{V_0 \to (1+\lambda)^\rho (V_0 + k) - k, \qquad k = \frac{1}{\rho(1-\beta)},}
#' because \eqn{u((1+\lambda)c) = (1+\lambda)^\rho u(c) + ((1+\lambda)^\rho - 1)/\rho}
#' and the extra (state- and period-independent) term discount-sums to
#' \eqn{((1+\lambda)^\rho - 1) k}.  Inverting \eqn{V_1 = (1+\lambda)^\rho (V_0 + k) - k}
#' gives
#' \deqn{\lambda = \left(\frac{V_1 + k}{V_0 + k}\right)^{1/\rho} - 1.}
#' \eqn{(V_0 + k) = (1/\rho) E[\sum_t \beta^t c_t^\rho]} shares the sign of
#' \eqn{1/\rho} with \eqn{(V_1 + k)}, so the ratio is positive and the power is
#' well-defined.
#'
#' LOG utility (\code{eis == 1}) is handled as a special branch, since the
#' general formula divides by \eqn{\rho = 0} there: scaling consumption by
#' \eqn{(1+\lambda)} in every period instead adds
#' \eqn{\log(1+\lambda) \sum_{t \ge 0} \beta^t = \log(1+\lambda)/(1-\beta)} to
#' the value, and inverting
#' \eqn{V_1 = V_0 + \log(1+\lambda)/(1-\beta)} gives
#' \eqn{\lambda = \exp((1-\beta)(V_1 - V_0)) - 1} -- the \eqn{eis \to 1}
#' (\eqn{\rho \to 0}) limit of the general-CRRA formula above.
#'
#' @param V0 Numeric: baseline per-cell level value (from
#'   \code{\link{hank_value_function}}).
#' @param V1 Numeric: alternative per-cell level value, same length/cell order
#'   as \code{V0}.
#' @param block A \code{\link{hank_het_block}} (used for \code{beta} and
#'   \code{eis}).
#'
#' @return Numeric vector, same length as \code{V0}/\code{V1}: the per-cell
#'   CEV \code{lambda}.
#' @export
hank_cev <- function(V0, V1, block) {
  beta <- block$beta
  eis  <- block$eis
  if (isTRUE(all.equal(eis, 1))) {
    ## Log utility: scaling consumption by (1+lambda) adds a state-independent
    ## log(1+lambda)/(1-beta) to the value; invert directly.
    return(exp((1 - beta) * (V1 - V0)) - 1)
  }
  ## General CRRA ("-1"-normalized u(c) = (c^rho - 1)/rho, rho = 1 - 1/eis):
  ## scaling the whole consumption path by (1+lambda) maps
  ##   V0 -> (1+lambda)^rho (V0 + const) - const,   const = 1/(rho (1-beta)),
  ## since u((1+lambda)c) = (1+lambda)^rho u(c) + ((1+lambda)^rho - 1)/rho and the
  ## extra term discount-sums to ((1+lambda)^rho - 1)*const. Invert for lambda.
  ## (V0 + const) = (1/rho) E[sum beta^t c^rho] shares the sign of 1/rho with
  ## (V1 + const), so the ratio is positive and the power well-defined; the
  ## eis -> 1 limit recovers the log branch above.
  rho   <- 1 - 1 / eis
  const <- 1 / (rho * (1 - beta))
  ((V1 + const) / (V0 + const))^(1 / rho) - 1
}


#' First-order per-cell welfare response to a transitory price path
#'
#' The FIRST-ORDER change in lifetime value \eqn{dV_0} induced by a transitory
#' path of price DEVIATIONS \code{r_path - block$r}, \code{w_path - block$w}
#' around the steady state, holding the household's policy and distribution
#' fixed at their steady-state values.  By the envelope theorem this
#' first-order approximation is exact to first order in the price path even
#' though it ignores the household's re-optimization: the household's policy
#' is itself chosen to maximize value, so its response to the price change is
#' second-order (Auclert 2019's interest-rate/labour-income "exposure"
#' representation of MPCs and welfare incidence). Concretely, only the DIRECT
#' budget windfall each period -- the extra resources \code{a*dr + e*dw} a
#' household would receive at its steady-state asset holdings \code{a} and
#' labour supply \code{e} -- matters, marginal-utility-weighted and
#' compounded forward through the steady-state transition operator:
#' \deqn{dV_0 = \sum_{p=0}^{T_h-1} \beta^p \Lambda^p g_{p+1}, \qquad
#'   g_t = u'(c_{ss}) \odot (a \, dr_t + e \, dw_t),}
#' with \eqn{dr_t = r_t - r_{ss}}, \eqn{dw_t = w_t - w_{ss}}, \code{a}/\code{e}
#' the steady-state asset/labour-efficiency grids broadcast to cells (asset
#' fast, income slow -- see \code{\link{hank_forward_operator}}), and
#' \eqn{\Lambda} the block's steady-state expectation operator (no
#' transpose): \eqn{\Lambda^p g} is the expected exposure \code{p} periods
#' ahead of each starting cell. \code{dV_0} is computed by backward Horner
#' recursion in \code{p} (no explicit matrix powers), assuming the price path
#' returns to steady state by \code{T_h} so the post-horizon continuation
#' value is unaffected (tail contribution 0).
#'
#' The accompanying \code{lambda} rescales \code{dV_0} into
#' consumption-equivalent units by dividing by the first-order value effect
#' of a constant 1-unit proportional scaling of the steady-state consumption
#' path forever: writing \code{up_ss} for the marginal utility
#' \eqn{u'(c_{ss})}, \code{dV_scale = (I - beta*Lambda)^{-1} (up_ss * c_ss)}
#' (solved by the same fixed-point recursion used by
#' \code{\link{hank_value_function}}, since \code{dV_scale} is itself a value
#' function for the constant flow payoff \code{up_ss * c_ss}). For log utility
#' (\code{eis == 1}), \code{up_ss * c_ss} equals 1 for every cell, so
#' \code{dV_scale == 1/(1-beta)} uniformly and \code{lambda == (1-beta)*dV0}
#' exactly -- a convenient internal check echoing the log branch of
#' \code{\link{hank_cev}}.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param r_path,w_path Numeric length-\code{T_h} paths of the aggregate
#'   return/wage LEVEL (not deviation).  Missing entries default to the
#'   steady-state constant \code{block$r} / \code{block$w}.
#' @param T_h Integer horizon (default \code{max(length(r_path),
#'   length(w_path))}; if both paths are \code{NULL}, defaults to 1).
#'
#' @return A list with:
#'   \item{dV0}{Numeric length-\code{n_e*n_a} vector: the first-order
#'     period-0 lifetime-value change per cell, in the block's cell order.}
#'   \item{lambda}{Numeric length-\code{n_e*n_a} vector: \code{dV0} rescaled
#'     into consumption-equivalent (CEV) units, \code{dV0 / dV_scale}.}
#'   \item{dV_scale}{Numeric length-\code{n_e*n_a} vector: the first-order
#'     value effect of a 1-unit constant proportional consumption-path
#'     scaling (the denominator used for \code{lambda}).}
#' @seealso \code{\link{hank_value_function}}, \code{\link{hank_cev}}
#' @export
hank_welfare_response <- function(block, r_path = NULL, w_path = NULL,
                                  T_h = NULL) {
  if (is.null(T_h))
    T_h <- max(length(r_path), length(w_path),
               if (is.null(r_path) && is.null(w_path)) 1L else 0L)
  if (is.null(r_path)) r_path <- rep(block$r, T_h)
  if (is.null(w_path)) w_path <- rep(block$w, T_h)
  stopifnot(length(r_path) == T_h, length(w_path) == T_h)

  beta <- block$beta
  Lam  <- block$Lambda
  n_e  <- block$n_e; n_a <- block$n_a
  n_cell <- n_e * n_a

  c_ss  <- .hank_mat_to_vec(block$c)
  up_ss <- c_ss^(-1 / block$eis)
  a_cell <- rep(block$a_grid, times = n_e)
  e_cell <- rep(block$e,      each  = n_a)

  dr <- r_path - block$r
  dw <- w_path - block$w

  ## Per-period direct budget windfall, marginal-utility weighted.
  g <- vector("list", T_h)
  for (t in seq_len(T_h))
    g[[t]] <- up_ss * (a_cell * dr[t] + e_cell * dw[t])

  ## dV0 = sum_{p=0}^{T_h-1} beta^p Lambda^p g_{p+1}, by backward Horner
  ## recursion (tail = 0 beyond T_h, since the path is back at steady state).
  acc <- numeric(n_cell)
  for (t in T_h:1L)
    acc <- g[[t]] + beta * as.numeric(Lam %*% acc)
  dV0 <- acc

  ## dV_scale: value effect of scaling c_ss by a constant 1 unit forever,
  ## i.e. the value function for constant flow payoff gscale = u'(c_ss)*c_ss,
  ## solved by the identical fixed-point recursion as hank_value_function().
  gscale <- up_ss * c_ss
  x <- gscale / (1 - beta)
  repeat {
    x_new <- gscale + beta * as.numeric(Lam %*% x)
    if (max(abs(x_new - x)) < 1e-11) { x <- x_new; break }
    x <- x_new
  }
  dV_scale <- x

  list(dV0 = dV0, lambda = dV0 / dV_scale, dV_scale = dV_scale)
}


#' Interest-rate vs labour-income channel decomposition of a welfare response
#'
#' Splits the FIRST-ORDER per-cell welfare response of
#' \code{\link{hank_welfare_response}} into its two additive
#' partial-equilibrium EXPOSURE channels: the interest-rate/asset channel and
#' the labour-income channel. Since the per-period direct budget windfall is
#' itself additive in the two price deviations,
#' \deqn{g_t = u'(c_{ss}) \odot (a \, dr_t + e \, dw_t) = g^r_t + g^w_t,
#'   \qquad g^r_t = u'(c_{ss}) \odot a \, dr_t, \quad
#'         g^w_t = u'(c_{ss}) \odot e \, dw_t,}
#' each channel is accumulated through the identical forward recursion used by
#' \code{\link{hank_welfare_response}},
#' \deqn{dV_0^r = \sum_{p=0}^{T_h-1} \beta^p \Lambda^p g^r_{p+1}, \qquad
#'       dV_0^w = \sum_{p=0}^{T_h-1} \beta^p \Lambda^p g^w_{p+1},}
#' so that, by linearity of the (backward Horner) recursion in its per-period
#' forcing term, \code{dV0_interest + dV0_labour} equals
#' \code{\link{hank_welfare_response}}'s \code{dV0} EXACTLY (to machine
#' precision, since both are computed by the same recursion, merely applied to
#' the two additive pieces of \code{g_t} instead of their sum). The
#' \code{lambda} channels use the SAME \code{dV_scale} denominator as
#' \code{\link{hank_welfare_response}} (the value effect of a constant 1-unit
#' proportional consumption-path scaling), so \code{lambda_interest +
#' lambda_labour} also equals \code{hank_welfare_response}'s \code{lambda}
#' exactly.
#'
#' This is a partial-equilibrium EXPOSURE decomposition at the household-block
#' level, in the spirit of Auclert (2019)'s unhedged-interest-rate-exposure
#' (URE) and earnings-exposure channels: it holds the aggregate prices'
#' PATHS as given (exogenous) and asks how much of a given household's
#' welfare response comes from its asset exposure to \code{dr} versus its
#' labour-income exposure to \code{dw}. It is NOT the full general-equilibrium
#' multi-channel decomposition (e.g. Auclert 2019's five-channel Fisher/
#' interest-rate/income/human-wealth/... split, or a sequence-space
#' price-feedback attribution): those require propagating the GE
#' price responses themselves back through the economy's Jacobians, which is
#' out of scope here.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param r_path,w_path Numeric length-\code{T_h} paths of the aggregate
#'   return/wage LEVEL (not deviation).  Missing entries default to the
#'   steady-state constant \code{block$r} / \code{block$w}.
#' @param T_h Integer horizon (default \code{max(length(r_path),
#'   length(w_path))}; if both paths are \code{NULL}, defaults to 1).
#'
#' @return A list with:
#'   \item{dV0_interest}{Numeric length-\code{n_e*n_a} vector: the
#'     interest-rate/asset-channel period-0 lifetime-value change per cell.}
#'   \item{dV0_labour}{Numeric length-\code{n_e*n_a} vector: the
#'     labour-income-channel period-0 lifetime-value change per cell.}
#'   \item{lambda_interest,lambda_labour}{The two channels rescaled into
#'     consumption-equivalent (CEV) units by the same \code{dV_scale} used by
#'     \code{\link{hank_welfare_response}}.}
#'   \item{dV0}{\code{dV0_interest + dV0_labour}, equal to
#'     \code{\link{hank_welfare_response}}'s \code{dV0}.}
#'   \item{lambda}{\code{lambda_interest + lambda_labour}, equal to
#'     \code{\link{hank_welfare_response}}'s \code{lambda}.}
#' @seealso \code{\link{hank_welfare_response}}, \code{\link{hank_value_function}}
#' @export
hank_welfare_channels <- function(block, r_path = NULL, w_path = NULL,
                                  T_h = NULL) {
  if (is.null(T_h))
    T_h <- max(length(r_path), length(w_path),
               if (is.null(r_path) && is.null(w_path)) 1L else 0L)
  if (is.null(r_path)) r_path <- rep(block$r, T_h)
  if (is.null(w_path)) w_path <- rep(block$w, T_h)
  stopifnot(length(r_path) == T_h, length(w_path) == T_h)

  beta <- block$beta
  Lam  <- block$Lambda
  n_e  <- block$n_e; n_a <- block$n_a
  n_cell <- n_e * n_a

  c_ss  <- .hank_mat_to_vec(block$c)
  up_ss <- c_ss^(-1 / block$eis)
  a_cell <- rep(block$a_grid, times = n_e)
  e_cell <- rep(block$e,      each  = n_a)

  dr <- r_path - block$r
  dw <- w_path - block$w

  ## Per-period direct budget windfall, split into its two additive pieces.
  g_r <- vector("list", T_h)
  g_w <- vector("list", T_h)
  for (t in seq_len(T_h)) {
    g_r[[t]] <- up_ss * (a_cell * dr[t])
    g_w[[t]] <- up_ss * (e_cell * dw[t])
  }

  ## Same backward-Horner recursion as hank_welfare_response(), applied
  ## separately to each additive piece of g_t (tail = 0 beyond T_h).
  acc_r <- numeric(n_cell)
  acc_w <- numeric(n_cell)
  for (t in T_h:1L) {
    acc_r <- g_r[[t]] + beta * as.numeric(Lam %*% acc_r)
    acc_w <- g_w[[t]] + beta * as.numeric(Lam %*% acc_w)
  }
  dV0_interest <- acc_r
  dV0_labour   <- acc_w

  ## dV_scale: identical computation to hank_welfare_response() (the value
  ## effect of scaling c_ss by a constant 1 unit forever), so the channel
  ## lambdas sum EXACTLY to hank_welfare_response()'s lambda.
  gscale <- up_ss * c_ss
  x <- gscale / (1 - beta)
  repeat {
    x_new <- gscale + beta * as.numeric(Lam %*% x)
    if (max(abs(x_new - x)) < 1e-11) { x <- x_new; break }
    x <- x_new
  }
  dV_scale <- x

  list(dV0_interest = dV0_interest, dV0_labour = dV0_labour,
       lambda_interest = dV0_interest / dV_scale,
       lambda_labour   = dV0_labour   / dV_scale,
       dV0 = dV0_interest + dV0_labour,
       lambda = dV0_interest / dV_scale + dV0_labour / dV_scale)
}


#' Efficiency/redistribution decomposition of a cross-sectional welfare change
#'
#' The population-measure, first-order form of the Dávila-Schaab (2024)
#' efficiency/redistribution split (their Proposition 1 / eq. 10; see
#' \code{references/HANK_WELFARE_SURVEY.md} for the verified statement quoted
#' from the paper). Given a per-cell consumption-equivalent welfare gain
#' \code{lambda} (e.g. \code{\link{hank_welfare_response}}'s \code{lambda}, or a
#' per-cell \code{\link{hank_cev}}), a population mass \code{D} over the same
#' cells, and a welfarist planner's per-cell Pareto weights \code{weights}, it
#' splits the planner's aggregate welfare change
#' \deqn{\Delta W = \sum_i D_i\, \tilde\omega_i\, \lambda_i}
#' (with \eqn{\tilde\omega} the weights renormalised to population-mass-mean 1)
#' into an efficiency and a redistribution component:
#' \deqn{\Xi^E = \sum_i D_i \lambda_i \;\;(\text{aggregate/Kaldor-Hicks efficiency}),}
#' \deqn{\Xi^{RD} = \sum_i D_i(\tilde\omega_i - 1)\lambda_i = \mathrm{Cov}_D(\tilde\omega,\lambda)\;\;(\text{redistribution}).}
#'
#' \eqn{\Xi^E} (the mass-weighted mean gain) is the compensation-principle
#' aggregate and is INVARIANT to the planner's weights (Dávila-Schaab Prop. 2a):
#' all planner disagreement lives in \eqn{\Xi^{RD}}, the mass-weighted covariance
#' of the normalised planner weight with the individual gain. With the default
#' uniform weights, \eqn{\Xi^{RD}=0} and \eqn{\Delta W = \Xi^E}.
#'
#' SCOPE / caveats (see the survey, section 5-6): (i) this is the FIRST-ORDER
#' (marginal) decomposition -- exact for the marginal welfare gains of
#' \code{\link{hank_welfare_response}}; applying it to a large discrete
#' steady-state \code{\link{hank_cev}} is a linearisation, not a theorem. (ii)
#' The finer aggregate-efficiency / risk-sharing / intertemporal-sharing split
#' (Dávila-Schaab Prop. 3) needs date/history-indexed welfare gains and its
#' specialisation to a heterogeneous-discount-factor steady state is NOT verified
#' -- deliberately not implemented here. (iii) The planner weights are the user's
#' normative choice; \code{dispersion} is a planner-FREE summary of how much
#' redistribution is \emph{available} to any non-uniform planner.
#'
#' @param lambda Numeric per-cell consumption-equivalent welfare gain.
#' @param D Numeric per-cell population mass (non-negative; used as-is as the
#'   population measure and renormalised to sum 1 internally -- for a mixture
#'   population pass the pooled mass \code{omega_k * block_k$D}).
#' @param weights Numeric per-cell welfarist Pareto weights (any positive scale;
#'   renormalised to population-mass-mean 1). \code{NULL} (default) = uniform
#'   (\eqn{\Xi^{RD}=0}).
#' @return A list with \code{efficiency} (\eqn{\Xi^E}), \code{redistribution}
#'   (\eqn{\Xi^{RD}}), \code{total} (\eqn{\Delta W}), and \code{dispersion} (the
#'   population-mass-weighted standard deviation of \code{lambda}).
#' @seealso \code{\link{hank_welfare_response}}, \code{\link{hank_cev}}
#' @export
hank_welfare_decompose <- function(lambda, D, weights = NULL) {
  stopifnot(length(lambda) == length(D), all(D >= 0), sum(D) > 0)
  Dn <- D / sum(D)                              # population measure, mass 1
  efficiency <- sum(Dn * lambda)                # Xi^E: mass-weighted mean gain
  if (is.null(weights)) {
    redistribution <- 0
  } else {
    stopifnot(length(weights) == length(lambda), all(weights >= 0))
    w_norm <- weights / sum(Dn * weights)       # renormalise to mass-mean 1
    redistribution <- sum(Dn * (w_norm - 1) * lambda)   # Cov_D(w_norm, lambda)
  }
  list(efficiency = efficiency, redistribution = redistribution,
       total = efficiency + redistribution,
       dispersion = sqrt(sum(Dn * (lambda - efficiency)^2)))
}
