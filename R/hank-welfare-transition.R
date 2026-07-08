## R/hank-welfare-transition.R
## --------------------------------------------------------------------------
## EXACT (finite, non-marginal) transition-path consumption-equivalent
## variation (CEV) for a HANK het block, as opposed to the FIRST-ORDER
## envelope welfare response of hank_welfare_response() (R/hank-welfare.R).
##
## hank_welfare_response() answers "what is the first-order change in
## lifetime value from a marginal price-path perturbation, holding the
## household's policy fixed at its steady-state optimum" (an envelope-theorem
## shortcut that is exact only to first order in the path). This file answers
## the harder, EXACT question: given an arbitrary (possibly large) transitory
## price path, what constant proportional consumption-equivalent scaling of
## the household's STEADY-STATE consumption path would make the household
## indifferent to actually living through the transition, accounting for its
## FULL nonlinear re-optimization along the path (not just the direct budget
## windfall)? The welfare survey (references/HANK_WELFARE_SURVEY.md) flags
## this exact finite-horizon transition CEV as missing from the literature's
## marginal/envelope framework -- most treatments stop at the first-order
## response because the exact object requires a full nonlinear backward solve
## of both the POLICY (EGM) and the VALUE (Bellman) recursions along the
## transition, which is what hank_value_transition() below does.
##
## Two backward passes over the SAME per-period distribution/policy objects
## (cell order, asset-fast: cell = (e-1)*n_a + a; see R/hank-distribution.R):
##   1. POLICY pass (identical to hank_td_nonlinear(), R/hank-het-block.R):
##      backward EGM from terminal Va_{T_h+1} = block$Va, storing each
##      period's (a_pol, c_pol).
##   2. VALUE pass: backward Bellman recursion using each period's OWN
##      (already-solved, fully re-optimized) policy and forward operator,
##      from a terminal continuation V_{T_h+1} = V_ss (the STATIONARY level
##      value from hank_value_function()) -- this is the object
##      test-hank-welfare-response.R's dV_full() oracle already computes
##      to validate the first-order response; here it IS the deliverable,
##      not just a cross-check.
##
## TERMINAL CONDITION: the price path is assumed to have returned to (and
## stayed at) its steady-state value by period T_h, so the household's
## continuation value beyond T_h is exactly the stationary value V_ss. The
## caller must choose T_h long enough for the path to have settled -- if the
## path has not converged back to steady state by T_h, the terminal
## condition is wrong and V_1 (hence the CEV) will be biased. There is no
## internal check for this (mirroring hank_td_nonlinear()'s identical
## assumption); see the first-order-agreement test in
## test-hank-welfare-transition.R for a worked example of a path chosen to
## have decayed to machine-zero well before its horizon.
##
## EXACT CEV FROM THE TRANSITION VALUE: hank_cev(V0, V1, block) answers "what
## constant proportional consumption scaling maps baseline lifetime value V0
## to alternative lifetime value V1". Applying it here with V0 = V_ss (the
## value of the household's STEADY-STATE consumption path forever) and
## V1 = V_1 (the value of actually living through the transition starting
## from period 1) is valid for exactly the same reason hank_cev() is valid in
## the steady-state comparison it was written for: it is a purely algebraic
## fact about the CRRA utility family that scaling ANY consumption path
## (steady state or transition path, held fixed) by a constant (1+lambda)
## forever maps its value via the closed form in hank_cev()'s roxygen. The
## BASELINE whose consumption path is being scaled is the steady-state path
## (value V_ss); V_1 is the alternative. No approximation is introduced by
## reusing hank_cev() here -- it is the same exact inversion, just applied to
## a transition-path V1 instead of a counterfactual-steady-state V1.
## --------------------------------------------------------------------------


#' Exact (finite, non-marginal) transition-path CEV for a HANK het block
#'
#' Computes the EXACT (not first-order) per-cell consumption-equivalent
#' variation between a household's stationary steady-state consumption path
#' and actually living through a transitory aggregate price path
#' \code{(r_path, w_path)}, fully accounting for the household's nonlinear
#' re-optimization (policy response) along the transition -- as opposed to
#' \code{\link{hank_welfare_response}}'s first-order envelope approximation,
#' which holds the policy fixed at its steady-state optimum.
#'
#' Two backward passes, sharing the exact recursion patterns of
#' \code{\link{hank_td_nonlinear}} (policy) and \code{\link{hank_value_function}}
#' (value fixed point):
#' \enumerate{
#'   \item \strong{Policy pass}: backward EGM household solve from terminal
#'     \code{Va_{T_h+1} = block$Va} down to period 1 (identical to
#'     \code{\link{hank_td_nonlinear}}), storing each period's asset and
#'     consumption policy (\code{a_pol[[t]]}, \code{c_pol[[t]]}).
#'   \item \strong{Value pass}: backward Bellman recursion
#'     \deqn{V_t = u(c_t) + \beta\, \Lambda_t V_{t+1}, \qquad t = T_h, \dots, 1,}
#'     using each period's OWN forward operator \eqn{\Lambda_t} (built from
#'     that period's re-optimized savings policy via
#'     \code{\link{hank_forward_operator}}) and consumption policy, from
#'     terminal continuation \eqn{V_{T_h+1} = V_{ss}} -- the STATIONARY level
#'     value from \code{\link{hank_value_function}}. This terminal condition
#'     assumes the price path has returned to (and stays at) steady state by
#'     \code{T_h}, so the true continuation beyond the horizon is exactly the
#'     stationary value; \strong{the caller must choose \code{T_h} long
#'     enough for the path to have settled}, or this terminal condition (and
#'     hence \code{cev}) will be biased. No internal check is made for this.
#' }
#'
#' The exact per-cell CEV is then \code{\link{hank_cev}(V0 = V_ss, V1 = V_1,
#' block)}: the constant proportional scaling of the STEADY-STATE consumption
#' path (value \code{V_ss}) that would make the household indifferent to
#' actually experiencing the transition (value \code{V_1}). This reuses
#' \code{hank_cev}'s closed-form CRRA inversion exactly (no approximation is
#' introduced): the algebraic fact that scaling any fixed consumption path by
#' \eqn{(1+\lambda)} forever maps its value through \code{hank_cev}'s closed
#' form holds regardless of whether the baseline path is the steady state or
#' a transition path -- only the BASELINE here (steady state, value
#' \code{V_ss}) needs to be a well-defined discounted-utility value, which it
#' is by construction.
#'
#' This exact finite-horizon transition CEV is the object the welfare survey
#' (\code{references/HANK_WELFARE_SURVEY.md}) flags as absent from the
#' literature's marginal/envelope-theorem framework
#' (\code{\link{hank_welfare_response}}): it agrees with the first-order
#' envelope response to \eqn{O(\epsilon)} for a small price-path perturbation
#' of scale \eqn{\epsilon}, but captures the full \eqn{O(\epsilon^2)} and
#' higher re-optimization terms that the envelope approximation discards.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param r_path,w_path Numeric length-\code{T_h} paths of the aggregate
#'   return/wage LEVEL (not deviation). Missing entries default to the
#'   steady-state constant \code{block$r} / \code{block$w} (in which case
#'   \code{cev} is exactly 0).
#' @param T_h Integer horizon (default \code{max(length(r_path),
#'   length(w_path))}; if both paths are \code{NULL}, defaults to 1). Must be
#'   long enough for the paths to have settled back to steady state (see
#'   Details).
#' @param tol,maxit Passed to \code{\link{hank_value_function}} for computing
#'   the terminal stationary value \code{V_ss}.
#' @param keep_policies logical (default \code{FALSE}). When \code{TRUE},
#'   the per-period objects the two backward passes build (and otherwise
#'   discard) are returned as well: \code{c_pol}/\code{a_pol} (length-
#'   \code{T_h} lists of \code{n_e x n_a} consumption/savings policy
#'   matrices) and \code{Lambda} (length-\code{T_h} list of sparse
#'   per-period expectation operators, same no-transpose convention as
#'   \code{block$Lambda}: \eqn{E_x[f(x')] = \Lambda_t f}; distributions push
#'   forward via \code{t(Lambda[[t]]) \%*\% d}). These are the inputs for
#'   date-indexed welfare analysis and goods-side (consumption) valuation
#'   along the transition. Off by default: \code{Lambda} costs
#'   \code{O(T_h)} sparse \code{n_cell x n_cell} matrices of memory.
#'
#' @return A list with:
#'   \item{V0}{Numeric length-\code{n_e*n_a} vector: the period-1 transition
#'     lifetime value \code{V_1} (named \code{V0} for the period-0/period-1
#'     convention shared with \code{\link{hank_welfare_response}}'s
#'     \code{dV0}: the value as seen from the start of the transition).}
#'   \item{V_path}{Numeric \code{(n_e*n_a) x T_h} matrix: column \code{t} is
#'     the transition value \code{V_t}.}
#'   \item{V_ss}{Numeric length-\code{n_e*n_a} vector: the terminal/baseline
#'     stationary value from \code{\link{hank_value_function}}.}
#'   \item{dV_path}{\code{V_path - V_ss} (recycled column-wise): the EXACT
#'     (non-marginal) per-period deviation from steady-state value.}
#'   \item{cev}{Numeric length-\code{n_e*n_a} vector: the exact per-cell
#'     consumption-equivalent variation, \code{\link{hank_cev}(V_ss, V0,
#'     block)}.}
#'   \item{c_pol, a_pol, Lambda}{Only when \code{keep_policies = TRUE}: the
#'     per-period policies and expectation operators (see the parameter
#'     description).}
#' @seealso \code{\link{hank_welfare_response}}, \code{\link{hank_cev}},
#'   \code{\link{hank_value_function}}, \code{\link{hank_td_nonlinear}}
#' @export
hank_value_transition <- function(block, r_path = NULL, w_path = NULL,
                                  T_h = NULL, tol = 1e-11, maxit = 100000L,
                                  keep_policies = FALSE) {
  if (is.null(T_h))
    T_h <- max(length(r_path), length(w_path),
               if (is.null(r_path) && is.null(w_path)) 1L else 0L)
  if (is.null(r_path)) r_path <- rep(block$r, T_h)
  if (is.null(w_path)) w_path <- rep(block$w, T_h)
  stopifnot(length(r_path) == T_h, length(w_path) == T_h)

  beta <- block$beta
  n_e <- block$n_e; n_a <- block$n_a
  n_cell <- n_e * n_a

  ## 1. Backward EGM policy pass (identical pattern to hank_td_nonlinear()):
  ## terminal Va_{T_h+1} = Va_ss, march backward storing each period's policy.
  a_pol <- vector("list", T_h)
  c_pol <- vector("list", T_h)
  Va <- block$Va
  for (t in T_h:1L) {
    step <- .hank_block_step(block, Va, r_path[t], w_path[t])
    a_pol[[t]] <- step$a
    c_pol[[t]] <- step$c
    Va <- step$Va
  }

  ## 2. Terminal/baseline stationary value.
  V_ss <- hank_value_function(block, tol = tol, maxit = maxit)

  ## 3. Backward Bellman value pass using each period's OWN re-optimized
  ## policy and forward operator, from terminal continuation V_{T_h+1} = V_ss.
  V_path <- matrix(0, n_cell, T_h)
  V_next <- V_ss
  Lam_keep <- if (keep_policies) vector("list", T_h) else NULL
  for (t in T_h:1L) {
    Lam_t <- hank_forward_operator(a_pol[[t]], block$a_grid, block$Pi)
    if (keep_policies) Lam_keep[[t]] <- Lam_t
    u_t   <- hank_utility(.hank_mat_to_vec(c_pol[[t]]), block$eis)
    V_t   <- u_t + beta * as.numeric(Lam_t %*% V_next)
    V_path[, t] <- V_t
    V_next <- V_t
  }
  V_1 <- V_path[, 1]

  ## 4. Exact CEV: baseline = the stationary steady-state consumption path
  ## (value V_ss), alternative = actually living through the transition
  ## (value V_1). Reuses hank_cev()'s closed-form CRRA inversion exactly.
  cev <- hank_cev(V0 = V_ss, V1 = V_1, block = block)

  out <- list(V0 = V_1, V_path = V_path, V_ss = V_ss,
              dV_path = V_path - V_ss, cev = cev)
  if (keep_policies) {
    out$c_pol  <- c_pol
    out$a_pol  <- a_pol
    out$Lambda <- Lam_keep
  }
  out
}
