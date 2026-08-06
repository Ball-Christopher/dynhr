## R/hank-het-block.R
## --------------------------------------------------------------------------
## Heterogeneous-agent "het block": packages a steady-state household solution
## (EGM policies + Young's-method distribution) with the operations needed to
## compute its sequence-space Jacobian (see R/hank-jacobian.R) and nonlinear
## perfect-foresight transitions.
##
## The reference household is a one-asset Krusell-Smith / one-asset-HANK
## household: CRRA utility, a single asset with return r, labour income
## y(e) = w * e over the idiosyncratic productivity states e.  Aggregate block
## OUTPUTS are:
##   A_t = sum_x D_t(x) a'(x)   (aggregate end-of-period assets / savings)
##   C_t = sum_x D_t(x) c(x)    (aggregate consumption)
## aggregated with the BEGINNING-of-period-t distribution D_t.  Block INPUTS are
## the aggregate paths {r_t, w_t}, plus -- when the block is built with a
## Pi_fn/Pi_inputs pair (HANK+SAM: hank_employment_income) -- named
## transition-probability inputs (e.g. the job-finding rate f_t and separation
## rate s_t) that rebuild the income transition matrix Pi_t period by period.
##
## TIMING CONVENTION for a time-varying Pi: Pi_t is the transition applied
## BETWEEN periods t and t+1.  It enters period t twice, consistently:
##   - backward step at t: expectations over date-t+1 idiosyncratic states use
##     Pi_t (Wa = beta * Pi_t %*% Va_{t+1});
##   - forward step at t: the distribution pushes D_{t+1} = Lambda_t' D_t with
##     Lambda_t built from (date-t savings policy, Pi_t).
## So a date-s perturbation of a transition input moves policies at all t <= s
## (anticipation via the value function) and the distribution from date s+1 on.
## --------------------------------------------------------------------------


#' Construct a one-asset heterogeneous-agent household block at steady state
#'
#' Solves the household EGM problem and its stationary distribution at fixed
#' aggregate prices \code{(r, w)}, and stores everything the sequence-space
#' Jacobian and nonlinear-transition routines need.
#'
#' @param a_grid Numeric: asset grid (see \code{\link{hank_asset_grid}}).
#' @param Pi Numeric \code{n_e x n_e}: income transition matrix.
#' @param e Numeric length-\code{n_e}: income levels (see
#'   \code{\link{hank_income_rouwenhorst}}).
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param r,w Steady-state real return and wage.
#' @param tol,maxit Passed to \code{\link{hank_egm_solve}}.
#' @param amin Numeric: borrowing constraint (minimum end-of-period assets).
#'   Must satisfy \code{amin >= a_grid[1]}. Defaults to \code{NULL}, which
#'   resolves to \code{a_grid[1L]} (the historical hardcoded behavior, so
#'   existing calls are byte-identical). Lets different household types on a
#'   shared \code{a_grid} face distinct borrowing limits (the wealth axis).
#' @param backend Character: \code{"cpp"} (default) or \code{"R"}, passed to
#'   \code{\link{hank_egm_solve}}. When \code{"cpp"}, the stationary
#'   distribution is also computed by the fused compiled kernel (savings
#'   policy -> distribution directly, skipping Lambda for the power
#'   iteration); \code{Lambda} itself is still always built via
#'   \code{\link{hank_forward_operator}} for the returned object, since it is
#'   part of this function's contract regardless of backend. Defaults to
#'   \code{getOption("dynhr.hank_backend", "cpp")}; the R reference path
#'   remains available via \code{backend = "R"} or
#'   \code{options(dynhr.hank_backend = "R")}.
#' @param Tr Lump-sum transfer added to income (\code{y = w*e + Tr*omega},
#'   with \code{omega} the incidence weight of \code{Tr_incidence}); default
#'   \code{0} reproduces the transfer-free household byte-for-byte. This is
#'   the aggregate input a fiscal block's profit-rebate or transfer closure
#'   feeds (name it \code{"Tr"} in the DAG; \code{T} itself is avoided because
#'   it is \code{TRUE} in R).
#' @param r_minus Borrowing rate applied on \code{a < 0}. \code{NULL} (the
#'   default) is the SYMMETRIC household, in which borrowers and savers both
#'   face \code{r}; otherwise a finite scalar, typically \code{r} plus a
#'   wedge. Reserved as a \code{Pi_fn} input name, so a transition function
#'   may respond to it.
#' @param Tr_incidence Incidence weight \eqn{\omega} distributing \code{Tr}
#'   across households: household in idiosyncratic state \code{(e, a)}
#'   receives \code{Tr * omega(e, a)}. Two forms are accepted:
#'   \describe{
#'     \item{\code{NULL} or a length-\code{length(e)} vector (Tier 1)}{The
#'       uniform rule (\code{NULL}, \code{omega == 1}) or any income-state-only
#'       weight, byte-for-byte identical to the pre-incidence household when
#'       uniform. Normalised here to \eqn{\sum_e \bar\pi_e \omega_e = 1}
#'       against the \code{Pi}-invariant distribution \eqn{\bar\pi}, so
#'       \code{Tr} always means the per-capita transfer and the aggregate
#'       outlay is exactly \code{Tr} at every date of a transition (the
#'       \code{e}-marginal of the distribution is invariant under the forward
#'       operator: \code{e} is an exogenous Markov chain and the steady-state
#'       marginal is already \code{Pi}-invariant, so
#'       \eqn{\sum_i D_{t,i} \omega_i \equiv 1} for all \code{t}). See
#'       \code{\link{hank_incidence_earnings}} for the earnings-proportional
#'       rule.}
#'     \item{A finite numeric \code{length(e) x length(a_grid)} matrix
#'       (Tier 2)}{An \code{(e, a)}-varying incidence rule, e.g.
#'       wealth-proportional. \code{omega} MUST be a function of the
#'       BEGINNING-of-period state only (indexed by the current grid, never
#'       \eqn{a'}): \code{Tr * omega(e, a)} then enters cash on hand
#'       additively without touching the marginal return to saving, so the
#'       Euler equation is undisturbed. Unlike the vector form, a matrix
#'       \code{omega} is taken AS SUPPLIED -- it is NOT normalised, because
#'       that would be circular (\code{D_ss} depends on the household policy,
#'       which depends on \code{omega}); the realised aggregate outlay is
#'       instead reported on the returned block as \code{Omega_ss} (see
#'       Value), and its date-\code{t} counterpart is available as the
#'       \code{"Omega"} het-block output (\code{\link{hank_het_jacobian}}).
#'       Because a non-uniform \code{omega} can push \code{w*e + Tr*omega(e,a)}
#'       toward or through zero in a low-\code{e}/low-\code{a} corner, the
#'       steady-state solve here ASSERTS strictly positive cash on hand at
#'       every \code{(e, a)} grid point (hard error, naming the offending
#'       corner and the sign of \code{Tr}, if violated); a transition-path
#'       evaluation that later drives cash on hand into the EGM tiny-floor
#'       region instead WARNS, once per session (see
#'       \code{\link{.hank_egm_step}}).}
#'   }
#' @param Va_init Optional \code{n_e x n_a} initial marginal value, passed
#'   straight through to \code{\link{hank_egm_solve}} (and to the wedge
#'   solver's symmetric starting solve). \code{NULL} (default) uses that
#'   solver's own cash-on-hand guess.
#' @param Pi_fn Optional function rebuilding the income transition matrix from
#'   named transition-probability inputs (e.g. the \code{Pi_fn(f, s)} returned
#'   by \code{\link{hank_employment_income}}). Supplying it makes those inputs
#'   perturbable aggregate inputs of the block alongside \code{(r, w)}: the
#'   nonlinear transition (\code{\link{hank_td_nonlinear}}, via
#'   \code{pi_input_paths}) and the fake-news Jacobian
#'   (\code{\link{hank_het_jacobian}}) then accept them by name. Must satisfy
#'   \code{do.call(Pi_fn, Pi_inputs) == Pi} at steady state (checked here).
#' @param Pi_inputs Named list of the steady-state values of \code{Pi_fn}'s
#'   inputs (e.g. \code{list(f = 0.7, s = 0.05)}). Required together with
#'   \code{Pi_fn} (and only then); names must not collide with
#'   \code{"r"}/\code{"w"}.
#'
#' @return An object of class \code{hank_het_block} with the steady-state
#'   policies (\code{a}, \code{c}), marginal value \code{Va}, distribution
#'   \code{D} (vector) and forward operator \code{Lambda}, aggregate
#'   steady-state outputs \code{A}, \code{C}, the borrowing constraint
#'   \code{amin}, the calibration, (when supplied) \code{Pi_fn} /
#'   \code{Pi_inputs}, \code{Omega_ss} (the realised steady-state transfer
#'   outlay per unit of \code{Tr}, \code{sum(D * Tr_incidence)}; exactly
#'   \code{1} to numerical precision for the Tier-1 vector/uniform form, and
#'   the genuine -- possibly drifting-along-a-transition -- outlay for a
#'   Tier-2 matrix \code{Tr_incidence}), and the run metadata
#'   \code{\link{hank_het_manifest}} reads:
#'   \describe{
#'     \item{\code{backend}}{The backend that ACTUALLY ran: the
#'       \code{backend} argument on a symmetric block, and \code{"R"} on a
#'       borrowing-wedge block (\code{r_minus} non-\code{NULL}), whose fixed
#'       point has no compiled path.}
#'     \item{\code{threads}}{Always \code{1L} -- the one-asset kernel is
#'       deliberately serial; see \code{\link{hank_egm_solve}}.}
#'     \item{\code{iterations}, \code{converged}}{From the household solve.}
#'     \item{\code{last_value_gap}, \code{last_policy_gap}}{The solve's final
#'       marginal-value and savings-policy gaps; \code{NA_real_} on the
#'       compiled backend and on the wedge path, neither of which returns
#'       them (see \code{\link{hank_egm_solve}}).}
#'     \item{\code{elapsed_solve}, \code{elapsed_dist}}{Wall-clock seconds
#'       (\code{proc.time()[["elapsed"]]} differences) spent in the household
#'       solve and in the stationary-distribution iteration, timed
#'       separately.}
#'     \item{\code{dist_converged}}{Logical: whether the stationary
#'       distribution's power iteration converged.}
#'   }
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 5)
#' ag  <- hank_asset_grid(50, 100, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                       r = 0.01, w = 1.0)
#' blk$A   # aggregate assets
#' @export
hank_het_block <- function(a_grid, Pi, e, beta, eis, r, w,
                           tol = 1e-11, maxit = 5000L, amin = NULL,
                           backend = getOption("dynhr.hank_backend", "cpp"),
                           Pi_fn = NULL, Pi_inputs = NULL, Tr = 0,
                           Tr_incidence = NULL, r_minus = NULL,
                           Va_init = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  if (!is.null(r_minus) &&
      (!is.numeric(r_minus) || length(r_minus) != 1L || !is.finite(r_minus)))
    stop("hank_het_block: 'r_minus' must be NULL (symmetric) or a finite ",
         "scalar (the borrowing rate on b < 0; typically r + a wedge).")
  if (!is.numeric(Tr) || length(Tr) != 1L || !is.finite(Tr))
    stop("hank_het_block: 'Tr' must be a finite numeric scalar (the ",
         "lump-sum transfer; 0 restores the transfer-free household).")
  omega <- .hank_normalize_incidence(Tr_incidence, e, Pi, "hank_het_block",
                                     n_a = length(a_grid))
  if (is.null(amin)) amin <- a_grid[1L]
  .hank_check_pi_fn(Pi_fn, Pi_inputs, Pi,
                    reserved = c("r", "w", "Tr", "r_minus"),
                    caller = "hank_het_block")
  inc       <- .hank_income_extra(w, e, Tr, omega)
  y         <- inc$y
  coh_extra <- inc$coh_extra
  if (!is.null(coh_extra)) {
    ## Tier 2 economic guard: strictly positive cash on hand at the SS solve
    ## (spec "Two economic points" -- see the Tr_incidence doc). A negative or
    ## zero corner would otherwise be silently absorbed by the EGM tiny-floor
    ## (see .hank_egm_step), so this fails loudly INSTEAD, naming the corner
    ## and the sign of Tr.
    coh_ss <- (1 + r) * matrix(a_grid, length(e), length(a_grid),
                               byrow = TRUE) + y + coh_extra
    if (any(coh_ss <= 0)) {
      bad <- which(coh_ss <= 0, arr.ind = TRUE)[1L, ]
      stop("hank_het_block: cash-on-hand is non-positive (",
           format(coh_ss[bad[1L], bad[2L]]), ") at (e = ", bad[1L],
           ", a_grid[", bad[2L], "] = ", format(a_grid[bad[2L]]),
           ") under the supplied matrix 'Tr_incidence' with Tr = ",
           format(Tr), " (", if (Tr < 0) "negative" else "non-negative",
           "). A wealth-/state-proportional incidence rule combined with ",
           "this Tr drives that household's beginning-of-period income to ",
           "zero or below. Reduce |Tr| or adjust the incidence weight in ",
           "the offending (e, a) corner.", call. = FALSE)
    }
  }
  t0_solve <- proc.time()[["elapsed"]]
  hh <- if (is.null(r_minus))
    hank_egm_solve(a_grid, y = y, r = r, beta = beta, eis = eis, Pi = Pi,
                   tol = tol, maxit = maxit, amin = amin, backend = backend,
                   Va_init = Va_init, coh_extra = coh_extra)
  else
    .hank_egm_solve_wedge(a_grid, y = y, r_plus = r, r_minus = r_minus,
                          beta = beta, eis = eis, Pi = Pi, amin = amin,
                          tol = tol, maxit = maxit, Va_init = Va_init,
                          coh_extra = coh_extra)
  elapsed_solve <- proc.time()[["elapsed"]] - t0_solve
  if (!hh$converged)
    warning("hank_het_block: household EGM did not converge at steady state")
  Lam <- hank_forward_operator(hh$a, a_grid, Pi)
  ## Timed separately from the solve, matching hank_het3_block(): the two
  ## stages have very different scaling, and a single combined figure cannot
  ## tell a slow household problem from a slow power iteration. Lambda is
  ## built outside the window because the cpp path does not use it for the
  ## distribution at all (fused kernel), so including it would make the two
  ## backends' elapsed_dist non-comparable.
  t0_dist <- proc.time()[["elapsed"]]
  sd  <- if (backend == "cpp")
    hank_stationary_dist_cpp(hh$a, a_grid, Pi, 1e-13, 200000L)
  else
    hank_stationary_dist(Lam)
  elapsed_dist <- proc.time()[["elapsed"]] - t0_dist
  D   <- sd$d
  ## What ACTUALLY ran, not what was asked for. The borrowing-wedge solver is
  ## a pure-R fixed point (.hank_egm_solve_wedge) that has no compiled path,
  ## so a wedge block ran the R backend whatever `backend` said -- recording
  ## the argument here would misreport the run.
  backend_used <- if (is.null(r_minus)) hh$backend else "R"
  ## Realised steady-state outlay, sum(D * omega): exactly 1 (to numerical
  ## precision) for the normalised Tier-1 vector/uniform form -- an
  ## independent sanity echo of the "outlay invariant" the normalisation is
  ## FOR -- and the genuine (possibly < 1 or > 1) Tier-2 outlay otherwise
  ## (spec item 2: reported rather than forced to 1, since forcing it would
  ## require normalising against the circular D_ss).
  omega_full <- if (is.matrix(omega)) omega else matrix(omega, length(e),
                                                         length(a_grid))
  structure(
    list(a_grid = a_grid, Pi = Pi, e = e, beta = beta, eis = eis,
         r = r, w = w, Tr = Tr, Tr_incidence = omega,
         Omega_ss = hank_aggregate(D, omega_full),
         r_minus = r_minus, amin = amin,
         a = hh$a, c = hh$c, Va = hh$Va,
         Lambda = Lam, D = D,
         A = hank_aggregate(D, hh$a),
         C = hank_aggregate(D, hh$c),
         n_e = length(e), n_a = length(a_grid),
         ## Run metadata, read by hank_het_manifest(). Exact indexing on the
         ## gaps: the wedge solver predates them and returns neither, and a
         ## missing field must surface as NA rather than as a partial match.
         backend = backend_used, threads = 1L,
         iterations = hh$iterations, converged = hh$converged,
         last_value_gap = .hank_na_real(hh[["last_value_gap", exact = TRUE]]),
         last_policy_gap = .hank_na_real(hh[["last_policy_gap", exact = TRUE]]),
         elapsed_solve = elapsed_solve, elapsed_dist = elapsed_dist,
         dist_converged = sd$converged,
         Pi_fn = Pi_fn, Pi_inputs = Pi_inputs),
    class = "hank_het_block")
}


#' A missing run-metadata scalar, as NA_real_
#'
#' Solve routines that predate a metadata field return \code{NULL} for it.
#' Recording \code{NULL} would DROP the field from the block's list (and so
#' silently shift every later field), while recording a made-up number would
#' be worse still, so the block stores \code{NA_real_}: present, and honest
#' about being unknown.
#'
#' @param v The value read off the solve, possibly \code{NULL}.
#' @return \code{v} when it is a length-1 numeric, \code{NA_real_} otherwise.
#' @keywords internal
.hank_na_real <- function(v) {
  if (is.null(v) || !is.numeric(v) || length(v) != 1L) NA_real_ else v
}


#' Validate an endogenous-transition (Pi_fn, Pi_inputs) pair
#'
#' The contract every block constructor that accepts endogenous transition
#' probabilities must enforce, in ONE place: supplied together or not at all;
#' \code{Pi_inputs} a non-empty, uniquely named list; names disjoint from the
#' block's own aggregate inputs (otherwise a Jacobian column name would be
#' ambiguous); and \code{Pi_fn} reproducing the steady-state \code{Pi} at
#' \code{Pi_inputs}.
#'
#' That last check is the one that matters. Without it a block can carry a
#' transition rule inconsistent with the \code{Pi} its own steady state was
#' solved at, and every derivative taken from it -- fake-news column, ND
#' oracle, Reiter emission -- is then differentiating around a point the block
#' is not actually sitting at, silently.
#'
#' @param Pi_fn,Pi_inputs The pair to validate (both \code{NULL} is valid).
#' @param Pi The block's steady-state transition matrix.
#' @param reserved Character vector of the caller's own aggregate input names.
#' @param caller Name of the calling constructor, for error messages.
#' @return \code{invisible(NULL)}; called for its errors.
#' @keywords internal
.hank_check_pi_fn <- function(Pi_fn, Pi_inputs, Pi, reserved, caller) {
  if (xor(is.null(Pi_fn), is.null(Pi_inputs)))
    stop(caller, ": supply Pi_fn and Pi_inputs together (or neither).")
  if (is.null(Pi_fn)) return(invisible(NULL))
  if (!is.function(Pi_fn))
    stop(caller, ": Pi_fn must be a function.")
  nm <- names(Pi_inputs)
  if (!is.list(Pi_inputs) || length(Pi_inputs) == 0L ||
      is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop(caller, ": Pi_inputs must be a non-empty, fully and ",
         "uniquely named list of steady-state transition-input values.")
  if (any(nm %in% reserved))
    stop(caller, ": Pi_inputs names must not collide with the ",
         "aggregate inputs (", paste0("'", reserved, "'", collapse = ", "),
         ").")
  Pi_check <- do.call(Pi_fn, Pi_inputs)
  if (!isTRUE(all.equal(unname(as.matrix(Pi_check)), unname(as.matrix(Pi)),
                        tolerance = 1e-10)))
    stop(caller, ": Pi_fn evaluated at Pi_inputs does not reproduce ",
         "the steady-state Pi (the transition inputs and Pi are ",
         "inconsistent).")
  invisible(NULL)
}


#' One nonlinear EGM backward step at given aggregate prices
#'
#' Thin wrapper mapping block inputs \code{(r, w)} to the household EGM step.
#' @param block A \code{\link{hank_het_block}}.
#' @param Va_p Next-period marginal value (\code{n_e x n_a}).
#' @param r,w Aggregate return and wage this period.
#' @param Pi Optional income transition matrix overriding the steady-state
#'   \code{block$Pi} for this period (time-varying transition probabilities;
#'   see the timing convention in the file header). Default \code{NULL} keeps
#'   the steady-state \code{Pi}, so existing calls are byte-identical.
#' @return List with \code{Va}, \code{a}, \code{c} (see \code{.hank_egm_step}).
#' @keywords internal
.hank_block_step <- function(block, Va_p, r, w, Pi = NULL, Tr = NULL,
                             r_minus = NULL) {
  amin <- if (!is.null(block$amin)) block$amin else block$a_grid[1L]
  if (is.null(Pi)) Pi <- block$Pi
  if (is.null(Tr)) Tr <- .hank_block_tr(block)
  if (is.null(r_minus)) r_minus <- block$r_minus     # NULL for symmetric blocks
  inc <- .hank_income_extra(w, block$e, Tr, .hank_block_omega(block))
  if (is.null(r_minus))
    .hank_egm_step(Va_p, block$a_grid, y = inc$y, r = r,
                   beta = block$beta, eis = block$eis, Pi = Pi,
                   amin = amin, coh_extra = inc$coh_extra)
  else
    .hank_egm_step_wedge(Va_p, block$a_grid, y = inc$y,
                         r_plus = r, r_minus = r_minus, beta = block$beta,
                         eis = block$eis, Pi = Pi, amin = amin,
                         coh_extra = inc$coh_extra)
}


#' Split incidence-weighted income into a per-\code{e} vector and an optional
#' \code{coh_extra} matrix
#'
#' Shared by \code{\link{hank_het_block}} (steady state) and
#' \code{\link{.hank_block_step}} (nonlinear transitions / the fake-news
#' sweep): the Tier-1 vector incidence rule folds \code{Tr * omega} into the
#' per-\code{e} income vector \code{y} exactly as before this function
#' existed (so that path is bit-identical); the Tier-2 matrix rule instead
#' keeps \code{y = w*e} and returns \code{Tr * omega} as a separate
#' \code{n_e x n_a} \code{coh_extra} matrix, additive on cash on hand (see
#' \code{\link{.hank_egm_step}}).
#'
#' @param w Numeric scalar wage.
#' @param e Numeric length-\code{n_e} income-state levels.
#' @param Tr Numeric scalar transfer.
#' @param omega Length-\code{n_e} vector or \code{n_e x n_a} matrix incidence
#'   weight (as stored on the block / returned by
#'   \code{\link{.hank_normalize_incidence}}).
#' @return List with \code{y} (length-\code{n_e}) and \code{coh_extra}
#'   (\code{n_e x n_a} matrix, or \code{NULL} for the Tier-1 form).
#' @keywords internal
.hank_income_extra <- function(w, e, Tr, omega) {
  if (is.matrix(omega)) list(y = w * e, coh_extra = Tr * omega)
  else list(y = w * e + Tr * omega, coh_extra = NULL)
}


#' Steady-state transfer of a household block, NULL-safe
#'
#' Blocks saved before the \code{Tr} field existed lack it; treat them as
#' transfer-free rather than erroring, so old objects keep working unchanged
#' (and note \code{y + 0} is bit-identical to \code{y} in IEEE arithmetic,
#' so the default path is byte-for-byte the pre-Tr behaviour).
#' Note the EXACT indexing: \code{block$Tr} would partial-match the sibling
#' field \code{Tr_incidence} on a block whose \code{Tr} has been dropped,
#' silently returning a length-\code{n_e} vector where a scalar is required.
#' @keywords internal
.hank_block_tr <- function(block) {
  tr <- block[["Tr", exact = TRUE]]
  if (is.null(tr)) 0 else tr
}


#' Transfer incidence weight of a household block, NULL-safe
#'
#' Blocks saved before the \code{Tr_incidence} field existed lack it; treat
#' them as uniform-incidence rather than erroring. Returns a length-\code{n_e}
#' vector of ones in that case, and since \code{Tr * 1} is bit-identical to
#' \code{Tr} in IEEE arithmetic the default path is byte-for-byte the
#' pre-incidence behaviour (the same argument \code{\link{.hank_block_tr}}
#' relies on for \code{y + 0}). On a Tier-2 block this returns the stored
#' \code{n_e x n_a} MATRIX unchanged (see \code{\link{hank_het_block}}'s
#' \code{Tr_incidence} doc); callers that need a length-\code{n_e} vector
#' either check \code{is.matrix()} first or go through
#' \code{\link{.hank_income_extra}}, which dispatches on it.
#' @keywords internal
.hank_block_omega <- function(block) {
  om <- block[["Tr_incidence", exact = TRUE]]
  if (is.null(om)) rep(1, length(block$e)) else om
}


#' Validate and normalise a transfer incidence weight
#'
#' Shared by \code{\link{hank_het_block}} and any caller needing the same
#' contract. \code{NULL} yields the uniform rule. A vector weight (Tier 1) is
#' scaled so that \eqn{\sum_e \bar\pi_e \omega_e = 1}, where \eqn{\bar\pi} is
#' the invariant distribution of \code{Pi}. A matrix weight (Tier 2,
#' \code{n_e x n_a}, when \code{n_a} is supplied) is validated but returned
#' AS-IS, deliberately not normalised (see the \code{Tr_incidence}
#' documentation on \code{\link{hank_het_block}}: normalising an
#' asset-dependent weight against \code{D_ss} would be circular).
#'
#' Normalising the VECTOR form against \eqn{\bar\pi} rather than against the
#' steady-state joint distribution \code{D} is what keeps that path
#' NON-CIRCULAR: \eqn{\bar\pi} is a property of \code{Pi} alone, available
#' before the household problem is solved, whereas \code{D} depends on the
#' policy which depends on \code{omega}. The two agree exactly, because the
#' \code{e}-marginal of \code{D} IS \eqn{\bar\pi} at any stationary
#' distribution. No such non-circular anchor exists for an \code{a}-dependent
#' weight, which is exactly why the matrix form is left unnormalised instead.
#'
#' @param omega \code{NULL} (uniform), a finite numeric length-\code{n_e}
#'   vector weight, or (when \code{n_a} is supplied) a finite numeric
#'   \code{n_e x n_a} matrix weight.
#' @param e Income-state grid (its length sets \code{n_e}).
#' @param Pi Income transition matrix.
#' @param who Calling function name, for error messages.
#' @param n_a Asset-grid length. \code{NULL} (default) means "no matrix form
#'   available here"; a matrix \code{omega} is then always a validation error
#'   (a caller that cannot make sense of a matrix incidence should not have
#'   to special-case that itself).
#' @return Numeric length-\code{n_e} normalised weight (Tier 1), or the
#'   validated \code{n_e x n_a} matrix unchanged (Tier 2).
#' @keywords internal
.hank_normalize_incidence <- function(omega, e, Pi, who = "hank_het_block",
                                      n_a = NULL) {
  n_e <- length(e)
  if (is.null(omega)) return(rep(1, n_e))
  if (is.matrix(omega)) {
    if (is.null(n_a))
      stop(who, ": a matrix 'Tr_incidence' is not supported here (no ",
           "asset-grid length available to validate against).")
    if (!is.numeric(omega) || nrow(omega) != n_e || ncol(omega) != n_a ||
        !all(is.finite(omega)))
      stop(who, ": a matrix 'Tr_incidence' must be a finite numeric ",
           n_e, " x ", n_a, " (n_e x n_a) matrix of beginning-of-period ",
           "(e, a) incidence weights (got ",
           if (is.numeric(omega))
             paste0(nrow(omega), " x ", ncol(omega)) else class(omega)[1L],
           ").")
    ## Tier 2: NOT normalised (see the roxygen note above) -- returned as
    ## supplied; the realised outlay is reported as Omega_ss on the block.
    return(omega)
  }
  if (!is.numeric(omega) || length(omega) != n_e || !all(is.finite(omega)))
    stop(who, ": 'Tr_incidence' must be NULL (uniform) or a finite numeric ",
         "vector of length length(e) = ", n_e, " (got ",
         if (is.numeric(omega)) paste0("length ", length(omega)) else
           class(omega)[1L], "). Incidence is indexed by the INCOME state ",
         "only -- an asset-indexed weight is not supported, because its ",
         "aggregate outlay would drift along a transition. Pass an ", n_e,
         " x n_a MATRIX (Tier 2) if you deliberately want an (e, a)-varying ",
         "rule; it is not normalised and its realised outlay is reported as ",
         "Omega_ss, not forced to 1.")
  pi_bar <- .hank_stationary(Pi)
  scale  <- sum(pi_bar * omega)
  if (!is.finite(scale) || abs(scale) < 1e-12)
    stop(who, ": 'Tr_incidence' has (near-)zero mass under the Pi-invariant ",
         "distribution (sum(pi_bar * omega) = ", format(scale), "), so it ",
         "cannot be normalised to a per-capita transfer. Supply a weight ",
         "with nonzero mean incidence.")
  omega / scale
}


#' Earnings-proportional transfer incidence weight
#'
#' The incidence rule under which each household's share of a transfer is
#' proportional to its labour earnings, \eqn{\omega_e \propto e}. Pass the
#' result as \code{Tr_incidence} to \code{\link{hank_het_block}}.
#'
#' This is the incidence counterpart of scaling the wage: it distributes
#' \code{Tr} exactly as wage income is distributed, while leaving \code{w}
#' itself free to be a separate aggregate input. That separation is the point
#' -- routing an accounting residual through \code{w} changes BOTH the
#' incidence of the residual and the household's exposure to the wage, and
#' those two effects are not otherwise separable.
#'
#' Requires \eqn{\sum_e \bar\pi_e e > 0} (true whenever \code{e} is a
#' nonnegative income grid that is not identically zero).
#'
#' @param e Income-state grid, as passed to \code{\link{hank_het_block}}.
#' @param Pi Income transition matrix, as passed to
#'   \code{\link{hank_het_block}}.
#' @return Numeric length-\code{length(e)} weight, normalised to
#'   \eqn{\sum_e \bar\pi_e \omega_e = 1}.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 5)
#' om  <- hank_incidence_earnings(inc$e, inc$Pi)
#' sum(inc$pi * om)   # 1
#' @export
hank_incidence_earnings <- function(e, Pi) {
  .hank_normalize_incidence(as.numeric(e), e, Pi, "hank_incidence_earnings")
}


#' Per-period transition-matrix path from transition-input paths
#'
#' Resolves \code{pi_input_paths} (named list of length-\code{T_h} LEVEL paths
#' for a subset of \code{names(block$Pi_inputs)}) into a length-\code{T_h}
#' list of transition matrices via \code{block$Pi_fn}, with missing inputs
#' held at their steady-state values. Returns \code{NULL} when no
#' transition-input path is supplied (caller then uses the steady-state
#' \code{block$Pi} everywhere -- the zero-allocation fast path).
#'
#' @param block Any block carrying a \code{Pi_fn}/\code{Pi_inputs} pair --
#'   \code{\link{hank_het_block}} or \code{\link{hank_het3_block}}; this reads
#'   only those two fields, so it is shared by both transition routes.
#' @param pi_input_paths Named list of length-\code{T_h} transition-input
#'   LEVEL paths, or \code{NULL}.
#' @param T_h Integer horizon.
#' @return \code{NULL}, or a length-\code{T_h} list of \code{n_e x n_e}
#'   transition matrices (\code{Pi_t}, applied between periods \code{t} and
#'   \code{t+1}).
#' @keywords internal
.hank_pi_path <- function(block, pi_input_paths, T_h) {
  if (is.null(pi_input_paths) || length(pi_input_paths) == 0L) return(NULL)
  if (is.null(block$Pi_fn))
    stop(".hank_pi_path: transition-input paths supplied, but the block has ",
         "no Pi_fn/Pi_inputs (rebuild it with the Pi_fn/Pi_inputs arguments ",
         "of its constructor).")
  nm  <- names(pi_input_paths)
  bad <- setdiff(nm, names(block$Pi_inputs))
  if (length(bad))
    stop(".hank_pi_path: unknown transition input(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; the block's Pi_inputs are ",
         paste0("'", names(block$Pi_inputs), "'", collapse = ", "), ".")
  ok <- vapply(pi_input_paths, function(p) length(p) == T_h, logical(1))
  if (!all(ok))
    stop(".hank_pi_path: every transition-input path must have length T_h (",
         T_h, ").")
  lapply(seq_len(T_h), function(t) {
    args <- block$Pi_inputs
    for (k in nm) args[[k]] <- pi_input_paths[[k]][t]
    do.call(block$Pi_fn, args)
  })
}


#' Nonlinear perfect-foresight transition of a het block
#'
#' Given aggregate input PATHS \code{r_path}, \code{w_path} over horizon
#' \code{T} (with terminal conditions returning to steady state), computes the
#' aggregate output paths by a backward household solve followed by a forward
#' distribution simulation.  This is the block's nonlinear map from input paths
#' to output paths, and is the brute-force reference used to validate the
#' sequence-space Jacobian.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param r_path,w_path Numeric length-\code{T} input paths (levels).  Missing
#'   entries default to the steady-state value.
#' @param T_h Integer horizon (default \code{length(r_path)}).
#' @param r_minus_path Optional length-\code{T} LEVEL path of the BORROWING
#'   rate (the rate applied on \code{a < 0}). \code{NULL} (default) holds it
#'   at the block's own \code{r_minus}, and if the block is symmetric
#'   (\code{r_minus = NULL}) the borrowing rate tracks \code{r_path}, which is
#'   the pre-wedge household exactly.
#' @param Tr_path Optional length-\code{T} LEVEL path of the lump-sum
#'   transfer; default holds it at the block's steady-state \code{Tr}. The
#'   block's own \code{Tr_incidence} weight distributes it across income
#'   states in every period (the path scales the aggregate, not the rule).
#' @param pi_input_paths Optional named list of length-\code{T} LEVEL paths
#'   for (a subset of) the block's transition-probability inputs
#'   (\code{names(block$Pi_inputs)}, e.g. \code{list(f = f_path)}); missing
#'   inputs stay at their steady-state values. Requires a block built with
#'   \code{Pi_fn}/\code{Pi_inputs}. \code{Pi_t = Pi_fn(inputs_t)} is applied
#'   between periods \code{t} and \code{t+1}: it enters the date-\code{t}
#'   backward expectation AND the date-\code{t} forward push (see the timing
#'   convention in the file header). Default \code{NULL} (constant
#'   steady-state \code{Pi}; existing calls are byte-identical).
#' @param keep_policies logical (default \code{FALSE}). When \code{TRUE},
#'   also return the per-period objects the solve builds and otherwise
#'   discards: \code{c_pol}/\code{a_pol} (length-\code{T_h} lists of
#'   \code{n_e x n_a} policy matrices from the backward pass) and
#'   \code{Lambda} (length-\code{T_h} list of sparse per-period forward
#'   operators from the forward pass; same no-transpose convention as
#'   \code{block$Lambda} -- distributions push forward via
#'   \code{t(Lambda[[t]]) \%*\% d}). Inputs for date-indexed distributional/
#'   welfare analysis along the transition. Off by default (\code{Lambda}
#'   costs \code{O(T_h)} sparse \code{n_cell x n_cell} matrices of memory).
#' @param D0 Optional initial (beginning-of-period-1) distribution, a
#'   length-\code{n_e*n_a} nonnegative vector summing to 1 in the package
#'   distribution order (\code{.hank_mat_to_vec}: asset index fastest).
#'   Default \code{NULL} starts from the steady-state \code{block$D}. Use for
#'   STATE-DEPENDENCE experiments (e.g. a high-debt vs low-debt initial
#'   cross-section facing the same shock). The backward pass is unchanged
#'   (terminal \code{Va = Va_ss}); only the forward simulation reweights.
#'
#' @return A list with numeric length-\code{T} paths \code{A} and \code{C}, and
#'   \code{Dpath} (\code{(n_e*n_a) x T} matrix): column \code{t} is the
#'   beginning-of-period-\code{t} distribution over which \code{A[t]},
#'   \code{C[t]} are aggregated. When \code{keep_policies = TRUE}, also
#'   \code{c_pol}, \code{a_pol}, and \code{Lambda} (see above); the exact
#'   aggregation identity \code{C[t] == hank_aggregate(Dpath[, t],
#'   c_pol[[t]])} holds by construction.
#' @export
hank_td_nonlinear <- function(block, r_path = NULL, w_path = NULL, T_h = NULL,
                              keep_policies = FALSE, pi_input_paths = NULL,
                              D0 = NULL, Tr_path = NULL,
                              r_minus_path = NULL) {
  .hank_reject_het2(block, "hank_td_nonlinear", use = "hank_td2_nonlinear")
  if (!is.null(D0)) {
    D0 <- as.numeric(D0)
    if (length(D0) != length(block$D) || any(D0 < -1e-12) ||
        abs(sum(D0) - 1) > 1e-8)
      stop("hank_td_nonlinear: D0 must be a nonnegative length-",
           length(block$D), " distribution summing to 1 (package order: ",
           "asset index fastest).", call. = FALSE)
  }
  if (is.null(T_h)) {
    lens <- c(length(r_path), length(w_path), length(Tr_path),
              length(r_minus_path),
              if (!is.null(pi_input_paths))
                vapply(pi_input_paths, length, integer(1)))
    T_h <- max(lens, if (all(lens == 0L)) 1L else 0L)
  }
  if (is.null(r_path))  r_path  <- rep(block$r, T_h)
  if (is.null(w_path))  w_path  <- rep(block$w, T_h)
  if (is.null(Tr_path)) Tr_path <- rep(.hank_block_tr(block), T_h)
  if (is.null(r_minus_path) && !is.null(block$r_minus))
    r_minus_path <- rep(block$r_minus, T_h)
  stopifnot(length(r_path) == T_h, length(w_path) == T_h,
            length(Tr_path) == T_h,
            is.null(r_minus_path) || length(r_minus_path) == T_h)
  ## Per-period transition matrices (NULL = steady-state Pi everywhere).
  Pi_path <- .hank_pi_path(block, pi_input_paths, T_h)
  Pi_at   <- function(t) if (is.null(Pi_path)) block$Pi else Pi_path[[t]]

  ## Backward: terminal Va_{T+1} = Va_ss.
  a_pol <- vector("list", T_h)
  c_pol <- vector("list", T_h)
  Va <- block$Va
  for (t in T_h:1L) {
    step <- .hank_block_step(block, Va, r_path[t], w_path[t], Pi = Pi_at(t),
                             Tr = Tr_path[t],
                             r_minus = if (is.null(r_minus_path)) NULL
                                       else r_minus_path[t])
    a_pol[[t]] <- step$a
    c_pol[[t]] <- step$c
    Va <- step$Va
  }

  ## Forward: beginning distribution D_1 = D_ss (or the caller's D0).
  D <- if (is.null(D0)) block$D else D0
  A <- numeric(T_h); C <- numeric(T_h)
  Dpath <- matrix(0, length(D), T_h)
  Lam_keep <- if (keep_policies) vector("list", T_h) else NULL
  for (t in seq_len(T_h)) {
    Dpath[, t] <- D
    A[t] <- hank_aggregate(D, a_pol[[t]])
    C[t] <- hank_aggregate(D, c_pol[[t]])
    Lam  <- hank_forward_operator(a_pol[[t]], block$a_grid, Pi_at(t))
    if (keep_policies) Lam_keep[[t]] <- Lam
    D    <- as.numeric(Matrix::t(Lam) %*% D)
  }
  out <- list(A = A, C = C, Dpath = Dpath)
  if (keep_policies) {
    out$c_pol  <- c_pol
    out$a_pol  <- a_pol
    out$Lambda <- Lam_keep
  }
  out
}
