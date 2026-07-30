## R/hank-twoasset-model.R
## --------------------------------------------------------------------------
## Wiring the two-asset household into the general block DAG, plus a
## demonstration two-asset GE model.
##
## The DAG engine itself is asset-count agnostic: topological ordering, the
## chain rule and the H_U inversion are all name-based. Only the two dispatch
## points in R/hank-model.R (.hank_block_jacobian and .hank_model_eval) need to
## know that kind = "het2" routes to hank_het2_jacobian / hank_td2_nonlinear.
##
## ROADMAP -- deliberately deferred variants of the household side (this file
## is the natural home for them; see briefs/19-twoasset-hank-scope.md section 2):
##   * KMV kinked cost (linear + convex, an inaction region). hank_egm2_solve
##     requires chi2 > 1 precisely because its illiquid-FOC crossing search
##     assumes a strictly convex, hence strictly monotone, marginal cost; the
##     kinked cost breaks that and needs different machinery.
##   * discrete adjust / don't-adjust choice with taste shocks (discrete-choice
##     EGM). The trigger for this is whether the convex spec can hit a target
##     hand-to-mouth share at defensible parameters -- a calibration question,
##     not a numerical one.
##   * two-asset Reiter / finite-state emission (see the W5 checkpoint).
## --------------------------------------------------------------------------


#' Wrap a two-asset household as a sequence-space block
#'
#' The two-asset counterpart of \code{\link{hank_het_block_spec}}.
#'
#' @param name Character: block name.
#' @param block A \code{\link{hank_het2_block}}.
#' @param inputs Character: aggregate inputs, a subset of
#'   \code{c("rb", "ra", "w", "Tr")} plus, for a block built with
#'   \code{Pi_fn}/\code{Pi_inputs}, the block's named transition-probability
#'   inputs.  Validated here so a bad wiring fails at spec time rather than
#'   inside the Jacobian dispatch.
#' @param outputs Character: a subset of \code{c("B", "A", "C", "CHI")}
#'   (aggregate liquid, illiquid, consumption, adjustment cost).  \code{CHI} is
#'   a real resource cost, so a GE resource constraint that nets it out must
#'   request it.
#'
#' @return An object of class \code{hank_block} (kind \code{"het2"}).
#' @seealso \code{\link{hank_het_block_spec}} (one-asset),
#'   \code{\link{hank_het2_block}}, \code{\link{hank_twoasset_model}}
#' @export
hank_het2_block_spec <- function(name, block, inputs = c("rb", "ra", "w"),
                                 outputs = c("B", "A", "C")) {
  if (inherits(block, "hank_het_block"))
    stop("hank_het2_block_spec(): this is a ONE-asset block ",
         "(hank_het_block); use hank_het_block_spec(), whose inputs are ",
         "('r', 'w') and outputs ('A', 'C').")
  .hank_het2_check_inputs(block, inputs)
  .hank_het2_check_outputs(outputs)
  structure(list(name = name, kind = "het2", inputs = inputs,
                 outputs = outputs, block = block),
            class = "hank_block")
}


#' Steady state of the demonstration two-asset GE economy
#'
#' A small-open-economy Krusell-Smith-style economy with TWO assets: a liquid
#' foreign bond position in fixed net supply \code{Bg}, and illiquid capital.
#' A representative firm rents
#' capital, so the illiquid return and the wage come from its marginal
#' products; the liquid return sits a \code{spread} below the illiquid one
#' (the reduced-form intermediation wedge that makes the liquid asset the
#' inferior store of value, which is what gives the adjustment cost something
#' to trade off against).
#'
#' \code{beta} is calibrated so the household's illiquid demand clears the
#' capital market at the supplied \code{K}, and \code{Bg} is then set to the
#' household's liquid demand, so the steady state clears both asset markets by
#' construction. The foreign sector pays net factor income \code{rb * Bg};
#' \code{hank_twoasset_model()} includes it explicitly in the aggregate resource
#' identity rather than treating bond interest as an unfinanced transfer.
#'
#' @param K Numeric: target capital stock (the illiquid market clears here).
#' @param Z Numeric: TFP level.
#' @param alpha,delta Capital share and depreciation rate.
#' @param spread Numeric >= 0: illiquid-minus-liquid return wedge.
#' @param eis,chi0,chi1,chi2 Household preference / adjustment-cost parameters.
#' @param b_grid,a_grid Liquid and illiquid grids.  Size them so essentially no
#'   stationary mass reaches either top: where a policy overshoots its grid the
#'   Young lottery clamps it, which breaks the aggregate resource identity and
#'   hence market clearing (see \code{briefs/19-twoasset-hank-scope.md} F13).
#'   Widen the span and refine the grid together as a remedy, then require
#'   \code{\link{hank_twoasset_grid_check}} to clear; there is no sufficient
#'   span-to-points formula because failures are non-monotone in the geometric
#'   knot placement. Widening
#'   \code{a_grid} alone stretches its top gaps, and the budget residual that
#'   \code{\link{.hank_egm2_step}} computes from separately-interpolated
#'   \eqn{a'} and \eqn{b'} then accumulates enough interpolation error at the
#'   richest cells to drive consumption negative -- which
#'   \code{\link{hank_egm2_solve}} rejects outright.  The defaults here
#'   (\code{amax = 100}, \code{n_a = 20}) leave ~5e-7 stationary mass at the
#'   illiquid top and satisfy the resource identity to ~2e-7 (F16).
#' @param Pi,e Income transition matrix and levels.
#' @param n_k Multiplier-grid size for the constrained branch.
#' @param beta_bracket Numeric length-2: search bracket for the calibrated
#'   \code{beta}.
#' @param tol Numeric: market-clearing tolerance for the \code{beta} solve.
#'
#' @return A list with the calibrated \code{beta}, the solved
#'   \code{\link{hank_het2_block}}, and the steady-state aggregates/prices.
#' @seealso \code{\link{hank_twoasset_model}}
#' @export
hank_twoasset_steady <- function(K = 3, Z = 1, alpha = 0.11, delta = 0.02,
                                 spread = 0.01, eis = 0.5,
                                 chi0 = 0.25, chi1 = 6.5, chi2 = 2,
                                 b_grid = hank_asset_grid(40, 12L, 0),
                                 a_grid = hank_asset_grid(100, 20L, 0),
                                 Pi = NULL, e = NULL, n_k = 10L,
                                 beta_bracket = c(0.90, 0.995), tol = 1e-8) {
  if (is.null(Pi) || is.null(e)) {
    inc <- hank_income_rouwenhorst(0.9, 0.7, 3L)
    if (is.null(Pi)) Pi <- inc$Pi
    if (is.null(e))  e  <- inc$e
  }
  ra <- alpha * Z * K^(alpha - 1) - delta
  w  <- (1 - alpha) * Z * K^alpha
  rb <- ra - spread
  if (rb <= -1)
    stop("hank_twoasset_steady: implied rb <= -1; reduce 'spread'.")

  mk <- function(beta)
    hank_het2_block(b_grid, a_grid, Pi, e, beta = beta, eis = eis,
                    rb = rb, ra = ra, w = w, chi0 = chi0, chi1 = chi1,
                    chi2 = chi2, n_k = n_k)
  ## Calibrate beta so illiquid demand clears the capital market at K.
  f <- function(beta) mk(beta)$A - K
  lo <- f(beta_bracket[1L]); hi <- f(beta_bracket[2L])
  if (lo > 0 || hi < 0)
    stop("hank_twoasset_steady: illiquid demand does not bracket K = ", K,
         " over beta in [", beta_bracket[1L], ", ", beta_bracket[2L],
         "] (A ranges ", format(lo + K), " to ", format(hi + K),
         "). Widen 'beta_bracket' or move 'K'.")
  beta <- stats::uniroot(f, beta_bracket, tol = tol)$root
  block <- mk(beta)
  grid_check <- hank_twoasset_grid_check(block)
  ## A failing diagnostic must not pass silently: every calibrated quantity
  ## below is contaminated by it. Bg is DEFINED as the household's liquid
  ## demand, and where a policy overshoots its grid the Young lottery clamps,
  ## so B (hence Bg, NFI, and the beta that uniroot just solved for) carries
  ## the clamping error. Warning rather than erroring because the tolerances
  ## are judgment calls the caller can retune; the returned $grid_check has
  ## the measured numbers either way.
  if (!isTRUE(grid_check$ok))
    warning("hank_twoasset_steady: the grid-adequacy diagnostic FAILED ",
            "(top-boundary policy mass ", format(grid_check$top_policy_mass),
            ", household resource residual ",
            format(grid_check$resource_residual), "). The calibrated beta, ",
            "Bg and NFI inherit this discretization error. ",
            grid_check$action, " -- see ?hank_twoasset_grid_check.",
            call. = FALSE)
  list(beta = beta, block = block, K = K, Z = Z, alpha = alpha, delta = delta,
       spread = spread, ra = ra, rb = rb, w = w,
       Bg = block$B, A = block$A, C = block$C, CHI = block$CHI,
       NFI = rb * block$B, grid_check = grid_check,
       b_grid = b_grid, a_grid = a_grid, Pi = Pi, e = e)
}


#' Diagnose two-asset steady-state grid adequacy
#'
#' Reports the two observable conditions behind the F13/F16 grid guidance.
#' Essentially no stationary mass should choose a policy at either upper grid
#' boundary, where Young's lottery clamps instead of preserving the policy
#' mean, and the stationary household resource identity
#' \deqn{C + CHI = Y + r_b B + r_a A}
#' should hold using beginning-of-period asset stocks. These are acceptance
#' diagnostics, not a claim that a fixed span-to-points formula is sufficient:
#' geometric-knot placement makes failures non-monotone in the span, so callers
#' should re-run this check after every grid change.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param top_tol Maximum stationary mass whose liquid or illiquid policy is at
#'   the corresponding upper grid boundary.
#' @param resource_tol Maximum absolute aggregate household-resource residual.
#' @return A list with \code{ok}, \code{top_policy_mass},
#'   \code{resource_residual}, beginning/end asset gaps, tolerances, and an
#'   actionable \code{action} string. A numerical-floor failure occurs before a
#'   block can be returned and is reported directly by \code{hank_egm2_solve()}.
#' @seealso \code{\link{hank_twoasset_steady}},
#'   \code{\link{hank_twoasset_htm_stats}}
#' @export
hank_twoasset_grid_check <- function(block, top_tol = 1e-5,
                                     resource_tol = 1e-5) {
  if (!inherits(block, "hank_het2_block"))
    stop("hank_twoasset_grid_check: 'block' must be a hank_het2_block.")
  for (x in list(top_tol, resource_tol))
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0)
      stop("hank_twoasset_grid_check: tolerances must be finite, non-negative scalars.")

  n_e <- block$n_e; n_b <- block$n_b; n_a <- block$n_a
  b_state <- .hank2_bcast_mid(block$b_grid, n_e, n_b, n_a)
  a_state <- .hank2_bcast_a(block$a_grid, n_e, n_b, n_a)
  income  <- .hank2_bcast_e(block$w * block$e, n_e, n_b, n_a)
  B_begin <- hank_aggregate2(block$D, b_state)
  A_begin <- hank_aggregate2(block$D, a_state)
  Y       <- hank_aggregate2(block$D, income)
  b_top <- .hank2_arr_to_vec(block$b) >= max(block$b_grid) - 1e-9
  a_top <- .hank2_arr_to_vec(block$a) >= max(block$a_grid) - 1e-9
  top_mass <- max(sum(block$D[b_top]), sum(block$D[a_top]))
  resid <- block$C + block$CHI - Y - block$rb * B_begin - block$ra * A_begin
  ok <- top_mass <= top_tol && abs(resid) <= resource_tol
  list(ok = ok, top_policy_mass = top_mass, resource_residual = resid,
       B_begin_minus_end = B_begin - block$B,
       A_begin_minus_end = A_begin - block$A,
       top_tol = top_tol, resource_tol = resource_tol,
       action = if (ok) "grid diagnostics clear" else paste0(
         "increase the relevant grid span until upper-bound policy mass is negligible, ",
         "and refine point counts with the span; accept only after this diagnostic clears"))
}


#' Wealthy hand-to-mouth statistics for a two-asset block
#'
#' Separates two statistics that are related at stationarity but not identical:
#' \code{state_at_floor} is the beginning-of-period mass located at the liquid
#' grid floor, while \code{policy_constrained} is the mass whose liquid policy
#' chooses that floor. Each reports the group's mean beginning-of-period
#' illiquid wealth, making release and paper claims reproducible without
#' silently switching definitions.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param tol Numerical tolerance around the liquid floor.
#' @return A list containing \code{population_mean_illiquid},
#'   \code{state_at_floor}, and \code{policy_constrained}; the latter two each
#'   contain \code{mass} and \code{mean_illiquid}.
#' @seealso \code{\link{hank_twoasset_grid_check}}
#' @export
hank_twoasset_htm_stats <- function(block, tol = 1e-10) {
  ## Discrete-adjustment blocks: "policy chooses the floor" is probabilistic
  ## (P mixes two branch policies), so the policy-constrained indicator is the
  ## P-weighted mixture of the branch floor indicators. State-at-floor is
  ## unchanged (a property of the distribution, not the policy).
  if (inherits(block, "hank_het2d_block")) {
    b_state <- .hank2_arr_to_vec(
      .hank2_bcast_mid(block$b_grid, block$n_e, block$n_b, block$n_a))
    a_state <- .hank2_arr_to_vec(
      .hank2_bcast_a(block$a_grid, block$n_e, block$n_b, block$n_a))
    Pv <- .hank2_arr_to_vec(block$P)
    fl_pol <- Pv * .hank2_arr_to_vec(block$b_A <= block$b_grid[1L] + tol) +
      (1 - Pv) * .hank2_arr_to_vec(block$b_N <= block$b_grid[1L] + tol)
    wt_state <- block$D * (b_state <= block$b_grid[1L] + tol)
    wt_pol   <- block$D * fl_pol
    summ <- function(wt) list(
      mass = sum(wt),
      mean_illiquid = if (sum(wt) > 0) sum(wt * a_state) / sum(wt) else NA_real_)
    return(list(population_mean_illiquid = sum(block$D * a_state),
                state_at_floor = summ(wt_state),
                policy_constrained = summ(wt_pol)))
  }
  if (!inherits(block, "hank_het2_block"))
    stop("hank_twoasset_htm_stats: 'block' must be a hank_het2_block or ",
         "hank_het2d_block.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("hank_twoasset_htm_stats: 'tol' must be a finite non-negative scalar.")
  b_state <- .hank2_arr_to_vec(
    .hank2_bcast_mid(block$b_grid, block$n_e, block$n_b, block$n_a))
  a_state <- .hank2_arr_to_vec(
    .hank2_bcast_a(block$a_grid, block$n_e, block$n_b, block$n_a))
  b_policy <- .hank2_arr_to_vec(block$b)
  summarize <- function(ix) {
    mass <- sum(block$D[ix])
    list(mass = mass,
         mean_illiquid = if (mass > 0)
           sum(block$D[ix] * a_state[ix]) / mass else NA_real_)
  }
  list(population_mean_illiquid = sum(block$D * a_state),
       state_at_floor = summarize(b_state <= block$b_grid[1L] + tol),
       policy_constrained = summarize(b_policy <= block$b_grid[1L] + tol))
}


#' Build the demonstration two-asset GE model through the block DAG
#'
#' Composes \code{\link{hank_twoasset_steady}} into a
#' \code{\link{hank_model}}: a firm block producing \code{(ra, w)} from lagged
#' capital and TFP, a foreign-income block producing \code{NFI = rb * Bg}, the
#' two-asset household, two asset-market targets, and an explicit aggregate
#' resource diagnostic.
#'
#' The system is deliberately \strong{2x2} -- unknowns \code{c("K", "rb")}
#' against targets \code{c("asset_mkt", "bond_mkt")} -- which exercises the
#' multi-unknown \code{H_U} packing with a heterogeneous-agent block, a
#' combination no one-asset model in the package covers (they are all 1x1 in
#' \code{K}).  \code{rb} is an unknown rather than a pass-through so the liquid
#' market clears a fixed net foreign bond position, which pins the liquid rate.
#' The transition spread \code{ra - rb} is therefore endogenous; there is
#' deliberately no spread block after the steady-state calibration.
#'
#' @param ts A \code{\link{hank_twoasset_steady}} steady state.
#' @param T_h Integer horizon.
#'
#' @return A \code{\link{hank_model}} with unknowns \code{c("K", "rb")},
#'   targets \code{c("asset_mkt", "bond_mkt")}, exogenous \code{Z}.
#' @seealso \code{\link{hank_twoasset_steady}}, \code{\link{hank_ks_model}}
#'   (the one-asset analogue)
#' @examples
#' \donttest{
#' ts <- hank_twoasset_steady()
#' m  <- hank_twoasset_model(ts, T_h = 96)
#' irf <- hank_model_irf(m, list(Z = 0.01 * 0.8^(0:95)))
#' head(irf$C)
#' }
#' @export
hank_twoasset_model <- function(ts, T_h) {
  alpha <- ts$alpha; delta <- ts$delta; Bg <- ts$Bg

  firm <- hank_simple_block(
    "firm", inputs = c("K", "Z"), outputs = c("ra", "w"),
    fn = function(paths, ss) {
      Klag <- c(ss$K, paths$K[-length(paths$K)])
      list(ra = alpha * paths$Z * Klag^(alpha - 1) - delta,
           w  = (1 - alpha) * paths$Z * Klag^(alpha))
    },
    jac = function(ss, T_h) {
      K <- ss$K; Z <- ss$Z
      lag <- rbind(0, cbind(diag(T_h - 1L), 0))
      list(ra = list(K = alpha * (alpha - 1) * Z * K^(alpha - 2) * lag,
                     Z = alpha * K^(alpha - 1) * diag(T_h)),
           w  = list(K = (1 - alpha) * alpha * Z * K^(alpha - 1) * lag,
                     Z = (1 - alpha) * K^(alpha) * diag(T_h)))
    })

  household <- hank_het2_block_spec("household", ts$block,
                                    inputs = c("rb", "ra", "w"),
                                    outputs = c("B", "A", "C", "CHI"))

  ## The liquid asset is a net foreign bond position. Its interest income is a
  ## real external resource, not an unfinanced domestic transfer.
  foreign <- hank_simple_block(
    "foreign", inputs = "rb", outputs = "NFI",
    fn = function(paths, ss) list(NFI = paths$rb * Bg),
    jac = function(ss, T_h) list(NFI = list(rb = Bg * diag(T_h))))

  ## Illiquid market: household illiquid demand == capital.
  asset <- hank_simple_block(
    "asset", inputs = c("A", "K"), outputs = "asset_mkt",
    fn = function(paths, ss) list(asset_mkt = paths$A - paths$K),
    jac = function(ss, T_h)
      list(asset_mkt = list(A = diag(T_h), K = -diag(T_h))))

  ## Liquid market: household liquid demand == fixed bond supply.
  bond <- hank_simple_block(
    "bond", inputs = "B", outputs = "bond_mkt",
    fn = function(paths, ss) list(bond_mkt = paths$B - Bg),
    jac = function(ss, T_h) list(bond_mkt = list(B = diag(T_h))))

  ## Walras-law diagnostic (not a third target): output plus foreign factor
  ## income finances consumption, adjustment costs, and net investment.
  resource <- hank_simple_block(
    "resource", inputs = c("K", "Z", "C", "CHI", "NFI"),
    outputs = "resource_mkt",
    fn = function(paths, ss) {
      Klag <- c(ss$K, paths$K[-length(paths$K)])
      Y <- paths$Z * Klag^alpha
      investment <- paths$K - (1 - delta) * Klag
      list(resource_mkt = Y + paths$NFI - paths$C - paths$CHI - investment)
    },
    jac = function(ss, T_h) {
      lag <- rbind(0, cbind(diag(T_h - 1L), 0))
      I <- diag(T_h)
      list(resource_mkt = list(
        K = (alpha * ss$Z * ss$K^(alpha - 1) + 1 - delta) * lag - I,
        Z = ss$K^alpha * I,
        C = -I, CHI = -I, NFI = I))
    })

  resource_ss <- ts$Z * ts$K^alpha + ts$NFI - ts$C - ts$CHI - ts$delta * ts$K
  ss <- list(K = ts$K, Z = ts$Z, ra = ts$ra, rb = ts$rb, w = ts$w,
             A = ts$A, B = ts$Bg, C = ts$C, CHI = ts$CHI, NFI = ts$NFI,
             asset_mkt = 0, bond_mkt = 0, resource_mkt = resource_ss)
  ## 2x2: (K, rb) against (asset_mkt, bond_mkt). K clears the illiquid market
  ## against capital; rb clears the liquid market against a FIXED net foreign
  ## bond position,
  ## which is what makes rb an unknown rather than a pass-through from ra --
  ## the steady-state `spread` is then just ra - rb, and moves endogenously
  ## along the transition rather than being imposed.
  hank_model(list(firm, foreign, household, asset, bond, resource),
             unknowns = c("K", "rb"),
             targets = c("asset_mkt", "bond_mkt"),
             exogenous = "Z", ss = ss, T_h = T_h)
}


#' Wrap a discrete-adjustment two-asset household as a sequence-space block
#'
#' The \code{\link{hank_het2_block_spec}} analogue for the fixed-cost
#' household (kind \code{"het2d"}): inputs a subset of
#' \code{c("rb", "ra", "w", "Tr")}, outputs a subset of
#' \code{c("B", "A", "C", "CHI", "ADJ")} (\code{ADJ} is the aggregate
#' adjustment probability -- observable against microdata on rebalancing
#' frequency).
#'
#' @param name Character block name.
#' @param block A \code{\link{hank_het2d_block}}.
#' @param inputs,outputs Character vectors, validated at spec time.
#' @return An object of class \code{hank_block} (kind \code{"het2d"}).
#' @seealso \code{\link{hank_het2d_block}}, \code{\link{hank_het2d_jacobian}}
#' @export
hank_het2d_block_spec <- function(name, block, inputs = c("rb", "ra", "w"),
                                  outputs = c("B", "A", "C")) {
  if (inherits(block, "hank_het2_block"))
    stop("hank_het2d_block_spec(): this is the SMOOTH two-asset block ",
         "(hank_het2_block); use hank_het2_block_spec().")
  if (inherits(block, "hank_het_block"))
    stop("hank_het2d_block_spec(): this is a ONE-asset block ",
         "(hank_het_block); use hank_het_block_spec().")
  .hank_het2d_check_inputs(block, inputs)
  .hank_het2d_check_outputs(outputs)
  structure(list(name = name, kind = "het2d", inputs = inputs,
                 outputs = outputs, block = block),
            class = "hank_block")
}
