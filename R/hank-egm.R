## R/hank-egm.R
## --------------------------------------------------------------------------
## Endogenous gridpoints method (EGM, Carroll 2006) household solver for a
## one-asset consumption-savings problem with CRRA utility and a discretized
## idiosyncratic income Markov chain.
##
## Household problem (recursive):
##   V(a, e) = max_{a' >= amin} u(c) + beta E[ V(a', e') | e ]
##   s.t.  c + a' = (1 + r) a + y(e),   u(c) = c^{1-1/eis} / (1 - 1/eis)
##
## EGM avoids root-finding: given next-period marginal value Va'(a', e'), it
## inverts the Euler equation on a fixed grid of end-of-period assets a' to
## recover consumption, then reads off the endogenous beginning-of-period asset
## that rationalizes each a'.  Policies are then interpolated back onto the
## fixed asset grid and the borrowing constraint is imposed.
##
## Matches the sequence-jacobian reference EGM (hetblocks/hh_sim.py) so its
## household-Jacobian goldens are reproducible: Va = (1+r) c^{-1/eis}, and the
## savings policy is interpolated in cash-on-hand space.
## --------------------------------------------------------------------------


#' Vectorized monotone linear interpolation with linear extrapolation
#'
#' Interpolates \code{y} (defined at strictly increasing knots \code{x}) at
#' query points \code{xq}.  Outside \code{[x[1], x[n]]} it extrapolates using
#' the nearest interval's slope (matching sequence-jacobian's \code{interpolate_y}).
#'
#' @param x Strictly increasing numeric knots.
#' @param y Numeric values at the knots (same length as \code{x}).
#' @param xq Numeric query points.
#' @return Numeric vector, same length as \code{xq}.
#' @keywords internal
.hank_interp1 <- function(x, y, xq) {
  n <- length(x)
  ## Interval index for each query: clamp to [1, n-1] so ends extrapolate.
  idx <- findInterval(xq, x)
  idx[idx < 1L]      <- 1L
  idx[idx > n - 1L]  <- n - 1L
  x0 <- x[idx]; x1 <- x[idx + 1L]
  y0 <- y[idx]; y1 <- y[idx + 1L]
  w  <- (xq - x0) / (x1 - x0)
  y0 + w * (y1 - y0)
}


#' One EGM backward step
#'
#' Given next-period marginal value \code{Va_p} on the income-by-asset grid,
#' returns updated marginal value and this period's savings/consumption policies.
#'
#' @param Va_p Numeric \code{n_e x n_a} matrix: next-period marginal value of
#'   assets, \eqn{V_a(a', e')}, rows = income states, cols = asset gridpoints.
#' @param a_grid Numeric length-\code{n_a}: fixed asset grid (increasing).
#' @param y Numeric length-\code{n_e}: income by state this period.
#' @param r Numeric: real return on assets this period.
#' @param beta Numeric: discount factor.
#' @param eis Numeric: elasticity of intertemporal substitution (\eqn{1/\gamma}).
#' @param Pi Numeric \code{n_e x n_e}: income transition matrix (row-stochastic).
#' @param amin Numeric: borrowing constraint (minimum end-of-period assets).
#'   Must satisfy \code{amin >= a_grid[1]} (the grid cannot represent assets
#'   below its floor). Defaults to \code{a_grid[1L]} (\code{NULL} resolves to
#'   the grid floor), matching the historical hardcoded behavior.
#' @param coh_extra Optional \code{n_e x n_a} numeric matrix, ADDITIVE on cash
#'   on hand at the FIXED (beginning-of-period) grid, i.e. \code{coh(e,a) <-
#'   coh(e,a) + coh_extra(e,a)}, entering only the actual-coh side of the
#'   interpolation (never the endogenous-grid side, which is a function of
#'   the CHOICE \eqn{a'}, not of the current state). \code{NULL} (default)
#'   is a pure no-op -- bit-identical to the pre-Tier-2 kernel -- and is how a
#'   matrix (asset-dependent) \code{Tr_incidence} on
#'   \code{\link{hank_het_block}} enters the household problem; see that
#'   constructor's Tier-2 documentation. When supplied, a warning fires
#'   ONCE PER SESSION if it drives \eqn{coh - a'} to the EGM tiny-floor
#'   region (see the transient-NaN guard below).
#'
#' @return A list with \code{Va} (updated marginal value), \code{a} (savings
#'   policy \eqn{a'(a,e)}), and \code{c} (consumption policy), each
#'   \code{n_e x n_a}.
#' @keywords internal
.hank_egm_step <- function(Va_p, a_grid, y, r, beta, eis, Pi, amin = NULL,
                           coh_extra = NULL) {
  n_e <- nrow(Va_p); n_a <- ncol(Va_p)
  if (is.null(amin)) amin <- a_grid[1L]
  if (amin < a_grid[1L])
    stop("hank_egm: amin (", amin, ") is below the asset grid floor (",
         a_grid[1L], "); the grid cannot represent assets below its floor.")
  ## Feasibility: the lowest cash-on-hand on the grid is
  ## (1+r)*a_grid[1] + min(y).  If amin exceeds it, a household at the grid
  ## floor with the lowest income cannot afford to sit at the constraint --
  ## the forced choice a'=amin gives c = coh - amin <= 0 (=> a NaN value that
  ## otherwise surfaces only as a cryptic downstream interpolation error).
  ## (This is exactly the natural-borrowing-limit condition r*a_grid[1] +
  ## min(y) >= amin - a_grid[1].)
  min_coh <- (1 + r) * a_grid[1L] + min(y)
  if (amin > min_coh)
    stop("hank_egm: amin (", format(amin), ") exceeds the lowest feasible ",
         "cash-on-hand on the grid ((1+r)*a_grid[1] + min(y) = ",
         format(min_coh), "); a household at the grid floor with the lowest ",
         "income could not afford to save amin, forcing non-positive ",
         "consumption. Reduce amin or raise the income/grid floor.")

  ## Transient-NaN guard (robustness for negative amin near the natural limit).
  ## When the lowest grid cash-on-hand (1+r)a_grid[1]+min(y) is negative (i.e.
  ## amin < -min(y)/(1+r)), transient EGM iterates can drive the discounted
  ## marginal value Wa or consumption non-positive, so Wa^(-eis) / c^(-1/eis)
  ## become NaN and the endogenous coh grid turns non-monotone -- crashing the
  ## findInterval inside .hank_interp1. Flooring at a tiny positive epsilon
  ## keeps every intermediate finite and coh_endog strictly increasing WITHOUT
  ## changing any converged policy: at the fixed point (and for any amin >= 0)
  ## these quantities are bounded well away from `tiny`, so the pmax() calls are
  ## exact no-ops and amin >= 0 results are byte-identical.
  tiny <- 1e-12
  ## Expected next-period marginal value, discounted (rows indexed by TODAY's e).
  Wa <- pmax(beta * (Pi %*% Va_p), tiny)          # n_e x n_a, over end-of-period a'
  ## Invert marginal utility: u'(c) = c^{-1/eis}  =>  c = Wa^{-eis}.
  c_endog   <- pmax(Wa^(-eis), tiny)               # consumption at each (e, a')
  ## Endogenous cash-on-hand that rationalizes choosing a': coh = c + a'.
  coh_endog <- c_endog + matrix(a_grid, n_e, n_a, byrow = TRUE)

  ## Actual cash-on-hand on the fixed grid: (1+r) a + y(e) [+ coh_extra(e,a)].
  coh <- (1 + r) * matrix(a_grid, n_e, n_a, byrow = TRUE) + y
  if (!is.null(coh_extra)) coh <- coh + coh_extra   ## Tier 2: (e,a) incidence

  a_pol <- matrix(0, n_e, n_a)
  for (e in seq_len(n_e)) {
    ## Interpolate a'(coh): knots = endogenous coh, values = a' grid.
    a_pol[e, ] <- .hank_interp1(coh_endog[e, ], a_grid, coh[e, ])
  }
  ## Impose the borrowing constraint.
  a_pol[a_pol < amin] <- amin

  ## Floor consumption before the envelope power (same transient-NaN guard as
  ## above; a no-op wherever coh - a_pol > 0, i.e. always at the fixed point and
  ## for amin >= 0, so converged / non-negative-amin results are byte-identical).
  c_pol_raw <- coh - a_pol
  ## Only checked when coh_extra is present (Tier 2), so the NULL/vector
  ## incidence path is a strict no-op here -- no new work, no new warning.
  if (!is.null(coh_extra) && any(c_pol_raw <= tiny))
    .hank_egm_warn_once(
      "coh_extra_floor",
      "hank_egm: a transition-path evaluation with a matrix 'Tr_incidence' ",
      "(Tier 2, asset-dependent incidence) drove cash-on-hand to the EGM ",
      "tiny-floor (1e-12) region at some (e, a) point -- consumption there ",
      "is being floored rather than reflecting the true budget. This ",
      "warning fires at most once per session; if unexpected, check the ",
      "sign/magnitude of 'Tr' against the incidence weight in that corner.")
  c_pol <- pmax(c_pol_raw, tiny)
  Va    <- (1 + r) * c_pol^(-1 / eis)              # envelope condition
  list(Va = Va, a = a_pol, c = c_pol)
}


## Package-private store for one-time warnings (avoids transition/Jacobian
## sweep warning spam -- the same pattern as .cumulant_warn_once).
.hank_egm_warn_env <- new.env(parent = emptyenv())

#' Emit a warning at most once per session, keyed by \code{key}.
#' @keywords internal
.hank_egm_warn_once <- function(key, ...) {
  if (isTRUE(.hank_egm_warn_env[[key]])) return(invisible(NULL))
  .hank_egm_warn_env[[key]] <- TRUE
  warning(paste0(...), call. = FALSE)
}


#' Solve the household problem to a steady-state policy (fixed prices)
#'
#' Iterates the EGM backward step at constant prices \code{(r, y)} until the
#' savings policy converges.
#'
#' @inheritParams .hank_egm_step
#' @param tol Numeric: convergence tolerance on the max absolute change in the
#'   savings policy between backward iterations.
#' @param maxit Integer: maximum backward iterations.
#' @param Va_init Optional \code{n_e x n_a} initial marginal value; if
#'   \code{NULL}, a standard consume-10\%-of-cash-on-hand guess is used.
#' @param amin Numeric: borrowing constraint (minimum end-of-period assets).
#'   Must satisfy \code{amin >= a_grid[1]}. Defaults to \code{NULL}, which
#'   resolves to \code{a_grid[1L]} (the historical hardcoded behavior, so
#'   existing calls are byte-identical). A tighter (higher) \code{amin}
#'   decouples the borrowing limit from the grid floor, letting different
#'   household types face different constraints on a shared grid.
#' @param backend Character: \code{"cpp"} (default) or \code{"R"}. \code{"cpp"}
#'   dispatches the whole backward-iteration loop to a compiled
#'   (RcppArmadillo) kernel -- no per-iteration R boundary crossing -- for
#'   the SAME numerics (same tiny = 1e-12 floors, same interpolation rule).
#'   Defaults to \code{getOption("dynhr.hank_backend", "cpp")}, so the
#'   package default is the compiled path; the R reference implementation
#'   remains available by passing \code{backend = "R"} or setting
#'   \code{options(dynhr.hank_backend = "R")} globally. Input validation and
#'   the feasibility (\code{min_coh}) pre-check run in this R wrapper
#'   regardless of backend, so error messages are unaffected by the choice.
#'
#' @return A list with converged \code{Va}, \code{a} (savings policy), \code{c}
#'   (consumption policy), the calibration (\code{a_grid}, \code{y}, \code{r},
#'   \code{beta}, \code{eis}, \code{Pi}, and the resolved \code{amin}), and the
#'   run metadata \code{\link{hank_het_manifest}} reads:
#'   \describe{
#'     \item{\code{iterations}}{Number of backward iterations run.}
#'     \item{\code{converged}}{Logical: whether \code{tol} was met.}
#'     \item{\code{last_policy_gap}}{The final \eqn{\max|a'_k - a'_{k-1}|},
#'       i.e. the convergence criterion's own value on the last iteration.
#'       \code{NA_real_} on \code{backend = "cpp"}: the compiled kernel tests
#'       the gap internally and does not return it, and re-deriving it in R
#'       would mean re-running a step, so it is reported as unknown rather
#'       than invented.}
#'     \item{\code{last_value_gap}}{The final \eqn{\max|V_a^k - V_a^{k-1}|}.
#'       \code{NA_real_} on \code{backend = "cpp"} for the same reason.}
#'     \item{\code{backend}}{Character: the backend that actually ran.}
#'     \item{\code{threads}}{Always \code{1L}. The one-asset kernel is
#'       deliberately serial (a whole solve is milliseconds, with only
#'       \code{n_e} parallel tasks -- a measured no-go), so \code{1} is the
#'       worker count that genuinely ran, not a placeholder.}
#'     \item{\code{elapsed}}{Wall-clock seconds
#'       (\code{proc.time()[["elapsed"]]} difference) spent in the backward
#'       iteration, excluding this wrapper's argument validation.}
#'   }
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 5)
#' a   <- hank_asset_grid(50, 100, 0)
#' hh  <- hank_egm_solve(a, y = inc$e, r = 0.01, beta = 0.96,
#'                       eis = 1, Pi = inc$Pi)
#' hh$converged
#' @param coh_extra Optional \code{n_e x n_a} numeric matrix, additive on cash
#'   on hand; see \code{\link{.hank_egm_step}}. \code{NULL} (default) is a
#'   pure no-op on both backends.
#' @export
hank_egm_solve <- function(a_grid, y, r, beta, eis, Pi,
                           tol = 1e-11, maxit = 5000L, Va_init = NULL,
                           amin = NULL,
                           backend = getOption("dynhr.hank_backend", "cpp"),
                           coh_extra = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  n_e <- length(y); n_a <- length(a_grid)
  ## Shared Markov contract (same validator as hank_forward_operator, so all
  ## public HANK routes agree on what a transition matrix is).
  .hank_check_markov(Pi, n_e, caller = "hank_egm_solve")
  if (is.null(amin)) amin <- a_grid[1L]
  if (amin < a_grid[1L])
    stop("hank_egm_solve: amin (", amin, ") is below the asset grid floor (",
         a_grid[1L], "); the grid cannot represent assets below its floor.")
  ## Feasibility pre-check, hoisted from .hank_egm_step so it fires
  ## IDENTICALLY for both backends before any compiled dispatch (the cpp
  ## backend bypasses .hank_egm_step entirely, so this can't rely on the
  ## per-iteration guard inside it). min_coh depends only on (a_grid, r, y),
  ## not on the iterate, so hoisting it here is a no-op for the R backend
  ## (same check, same message, just evaluated once up front instead of on
  ## the first iteration) -- existing error-message tests are unaffected.
  min_coh <- (1 + r) * a_grid[1L] + min(y)
  if (amin > min_coh)
    stop("hank_egm: amin (", format(amin), ") exceeds the lowest feasible ",
         "cash-on-hand on the grid ((1+r)*a_grid[1] + min(y) = ",
         format(min_coh), "); a household at the grid floor with the lowest ",
         "income could not afford to save amin, forcing non-positive ",
         "consumption. Reduce amin or raise the income/grid floor.")

  if (is.null(Va_init)) {
    coh0 <- (1 + r) * matrix(a_grid, n_e, n_a, byrow = TRUE) + y
    ## Floor the initial cash-on-hand before the negative power: with a negative
    ## borrowing limit the lowest grid nodes can have coh0 < 0, making the
    ## default guess (0.1*coh0)^(-1/eis) NaN. pmax(., tiny) keeps the guess
    ## finite; it is a no-op wherever coh0 > 0 (all amin >= 0 calibrations), so
    ## existing results are byte-identical. (See the transient-NaN guard in
    ## .hank_egm_step.)
    Va   <- (1 + r) * pmax(0.1 * coh0, 1e-12)^(-1 / eis)
  } else {
    Va <- Va_init
  }

  ## Timed from HERE, not from the top of the function: the manifest's
  ## elapsed_solve is meant to be the cost of the backward iteration itself,
  ## and argument validation is not backend-dependent (so including it would
  ## make the two backends' timings non-comparable). Same convention as
  ## hank_egm3_solve().
  t0 <- proc.time()[["elapsed"]]

  if (backend == "cpp") {
    ## The whole backward-iteration loop runs inside this one .Call -- no
    ## per-iteration R boundary crossing. Va (either the default guess just
    ## computed above, or the caller's Va_init) is passed through explicitly
    ## so the cpp path never needs its own default-Va formula to keep in
    ## sync with the R one.
    res <- hank_egm_solve_cpp(a_grid, y, r, beta, eis, Pi, amin, tol,
                              as.integer(maxit), Va, coh_extra)
    return(list(Va = res$Va, a = res$a, c = res$c,
                iterations = as.integer(res$iterations),
                converged = as.logical(res$converged),
                ## The compiled kernel computes the policy gap only to TEST
                ## it and returns Va/a/c alone, so neither gap survives the
                ## .Call. Recovering one would mean running an extra backward
                ## step, which is a DIFFERENT quantity (the gap of one further
                ## iterate, not the gap between the last two) -- reported as
                ## unknown instead of quietly substituted.
                last_value_gap = NA_real_, last_policy_gap = NA_real_,
                backend = backend, threads = 1L,
                elapsed = proc.time()[["elapsed"]] - t0,
                a_grid = a_grid, y = y, r = r, beta = beta, eis = eis, Pi = Pi,
                amin = amin))
  }

  a_old <- matrix(-Inf, n_e, n_a)
  converged <- FALSE
  it <- 0L
  ## The R loop has both iterates in hand, so both gaps are MEASURED here.
  ## policy_gap is exactly the quantity the convergence test reads, so the
  ## recorded number and the accept/reject decision can never disagree.
  policy_gap <- NA_real_; value_gap <- NA_real_
  for (it in seq_len(maxit)) {
    Va_prev <- Va
    step <- .hank_egm_step(Va, a_grid, y, r, beta, eis, Pi, amin = amin,
                           coh_extra = coh_extra)
    Va   <- step$Va
    value_gap  <- max(abs(Va - Va_prev))
    policy_gap <- max(abs(step$a - a_old))
    if (policy_gap < tol) { converged <- TRUE; a_old <- step$a; break }
    a_old <- step$a
  }
  list(Va = Va, a = step$a, c = step$c,
       iterations = it, converged = converged,
       last_value_gap = value_gap, last_policy_gap = policy_gap,
       backend = backend, threads = 1L,
       elapsed = proc.time()[["elapsed"]] - t0,
       a_grid = a_grid, y = y, r = r, beta = beta, eis = eis, Pi = Pi,
       amin = amin)
}


#' Euler-equation residual of a converged household policy
#'
#' Self-contained correctness oracle for \code{\link{hank_egm_solve}}: at every
#' UNCONSTRAINED gridpoint the intertemporal Euler equation
#' \eqn{u'(c) = \beta (1+r) E[u'(c')]} must hold to numerical precision.
#'
#' @param hh A solved household object from \code{\link{hank_egm_solve}}. The
#'   borrowing constraint \code{amin} is read from \code{hh$amin} (falling
#'   back to \code{hh$a_grid[1L]} for older objects that predate the
#'   \code{amin} field).
#' @param constraint_tol Points whose savings policy is within
#'   \code{constraint_tol} of \code{amin} are treated as borrowing-constrained
#'   and excluded (the Euler equation holds with inequality there).
#'
#' @return A list with \code{max_abs} (max absolute residual over unconstrained
#'   points) and \code{residual} (the full \code{n_e x n_a} residual matrix,
#'   \code{NA} at constrained points).
#' @export
hank_euler_residual <- function(hh, constraint_tol = 1e-8) {
  .hank_reject_het2(hh, "hank_euler_residual", use = NULL)
  .hank_reject_wedge(hh, "hank_euler_residual")
  a_grid <- hh$a_grid; Pi <- hh$Pi; r <- hh$r; beta <- hh$beta; eis <- hh$eis
  n_e <- length(hh$y); n_a <- length(a_grid)
  amin <- if (!is.null(hh$amin)) hh$amin else a_grid[1L]

  ## Next-period consumption at chosen a', for each future state e'.
  ## c'(a'(a,e), e') : interpolate the consumption policy at the chosen savings.
  uc      <- hh$c^(-1 / eis)                       # u'(c) today, n_e x n_a
  exp_ucp <- matrix(0, n_e, n_a)
  for (e in seq_len(n_e)) {
    ap_e <- hh$a[e, ]
    ucp_next <- matrix(0, n_e, n_a)               # over e', at query a' = ap_e
    for (ep in seq_len(n_e)) {
      c_ep <- .hank_interp1(a_grid, hh$c[ep, ], ap_e)
      ucp_next[ep, ] <- c_ep^(-1 / eis)
    }
    exp_ucp[e, ] <- as.numeric(Pi[e, ] %*% ucp_next)
  }
  resid <- uc - beta * (1 + r) * exp_ucp

  constrained <- hh$a <= (amin + constraint_tol)
  resid[constrained] <- NA_real_
  list(max_abs = max(abs(resid), na.rm = TRUE), residual = resid)
}
