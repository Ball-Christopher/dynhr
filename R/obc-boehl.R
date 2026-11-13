## R/obc-boehl.R
## --------------------------------------------------------------------------
## Boehl (2022) / OccBin piecewise-linear OBC path solver.
##
## Provides:
##   boehl_simulate()            -- forward simulation under a fixed regime path
##   boehl_fb_residual()         -- Fischer-Burmeister complementarity residuals
##   boehl_end_of_spell_gap()    -- end-of-spell consistency check (scalar)
##   boehl_find_spell_duration() -- bisection root-finder: scalar k_star
##   boehl_solve_regime_path()   -- main solver; dispatches to spell-duration
##                                  fast path (single_spell=TRUE) or iterative
##   compute_irfs_obc()          -- OBC-aware IRF computation (obc_solver arg)
##
## THE PIECEWISE-LINEAR SOLUTION
##
##   A regime path is simulated with TIME-VARYING rules from the backward
##   recursion of .obc_pwl_rules() (R/obc-binding.R): the rule of a binding
##   period depends on the regimes ahead of it, and a slack period that
##   anticipates a later binding period does not follow the slack rule.  The
##   per-regime policies of obc_ensure_policy() assume the NEXT period is
##   slack, so simulating a spell of 2+ periods with them solved every period
##   but the last with the wrong expectation.
##
##   A regime path is accepted when it reproduces itself on the CONSTRAINED
##   path (.obc_pwl_check()): slack periods satisfy every bound, binding
##   periods have a multiplier of the right sign.  The old search flagged a
##   constraint only when the SLACK policy violated it, so a bound violated
##   only because another constraint binds was never imposed.
##
##   Shock sequences are surprises (Dynare occbin_solver): at each period with
##   a non-zero shock the expected regime path is re-solved from the current
##   state.  With one shock in period 1 (an IRF) the result is the exact
##   perfect-foresight piecewise-linear path.
##
## TWO SOLVER PATHS
##
##   Iterative (single_spell = FALSE, default)
##     OccBin guess-and-verify on the constrained path.  Handles arbitrary
##     regime patterns: multi-spell, staggered and interacting constraints.
##
##   Spell-duration root-finder (single_spell = TRUE) — Boehl 2022 §3
##     Parameterises the regime path by a single scalar k (spell length) and
##     bisects for the shortest k whose end-of-spell gap is non-negative.
##     Assumption: exactly one contiguous binding block [1, k_star] followed
##     by permanent slack.  When the resulting path is not an equilibrium
##     (the assumption fails) the solver falls through, with a classed
##     message (dynhr_obc_spell_fallthrough), to the iterative path.
##
## Reference: Boehl G. (2022) "Efficient solution and computation of models
##   with occasionally binding constraints." J. Econ. Dyn. Control 143, 104495.
##   Guerrieri L. & Iacoviello M. (2015) "OccBin". J. Monet. Econ. 70, 22-38.
## --------------------------------------------------------------------------


# =============================================================================
# Forward simulation
# =============================================================================

#' Forward simulation of an OBC model under a fixed regime path
#'
#' Simulates the model from an initial state (default: zero) for T periods
#' under the piecewise-linear rules of \code{regime_path}: period t follows
#' the rule obtained by the backward recursion over the regimes of periods
#' t, t+1, ..., slack after the last binding period (see
#' \code{.obc_pwl_rules()}).  Shocks are surprises.
#'
#' Endo equation:  y_t = ghx_t * s_{t-1} + ghu_t * eps_t + c_t
#'
#' @param shock_seq    n_exo x T numeric matrix of structural shocks
#' @param dr_slack     Slack-regime DecisionRules (unused; kept for the
#'                     interface -- the cache carries the system)
#' @param regime_cache R environment seeded by obc_ensure_policy()
#' @param regime_path  Integer vector (length T): bitfield regime per period
#'                     (0 = all slack; bit j = 1 means spec j binds)
#' @param state_init   Numeric vector (length n_state) initial state;
#'                     NULL or omitted → zero initial state
#' @return List with:
#'   $paths  -- n_endo x T matrix (rownames = endo variable names)
#'   $states -- n_state x T matrix: s_t AFTER applying the period-t rule
#' @noRd
boehl_simulate <- function(shock_seq, dr_slack, regime_cache, regime_path,
                            state_init = NULL) {
  ctx   <- .obc_pwl_cache_context(regime_cache)
  rules <- .obc_pwl_rules(ctx, regime_path)
  .obc_pwl_forward(ctx, rules, shock_seq, state_init)
}


# =============================================================================
# Fischer-Burmeister complementarity residual
# =============================================================================

#' Fischer-Burmeister complementarity residuals of a regime path
#'
#' For each period t and constraint j evaluates
#'   phi(a_jt, b_jt) = a_jt + b_jt - sqrt(a_jt^2 + b_jt^2)
#' with b_jt = sgn_j (x_j(t) - bound_j) the slack of the bound on the
#' simulated path and a_jt = sgn_j F_j(t) the multiplier (residual of the
#' tagged equation on the path, with E_t y_{t+1}); sgn_j = +1 for a lower
#' bound, -1 for an upper bound.  phi = 0 iff a >= 0, b >= 0 and ab = 0: at
#' an equilibrium regime path every entry is zero to rounding.
#'
#' @param sim          List from boehl_simulate() under regime_path
#' @param shock_seq    n_exo x T matrix of structural shocks
#' @param regime_cache R environment seeded by obc_ensure_policy()
#' @param regime_path  Integer vector (length T): the simulated regime path
#' @param state_init   Initial state (NULL = zero)
#' @return n_spec x T numeric matrix of FB residuals (zero at a solution)
#' @noRd
boehl_fb_residual <- function(sim, shock_seq, regime_cache, regime_path,
                              state_init = NULL) {
  ctx   <- .obc_pwl_cache_context(regime_cache)
  rules <- .obc_pwl_rules(ctx, regime_path)
  .obc_pwl_check(ctx, rules, sim, shock_seq, regime_path, state_init)$phi
}


# =============================================================================
# Spell-duration root-finder (Boehl 2022 §3)
# =============================================================================

#' End-of-spell consistency gap for a candidate spell length k
#'
#' The candidate regime path binds \code{spell_regime} in periods 1..k and is
#' slack afterwards; it is simulated with its time-varying rules.  The gap is
#'
#'   min_j sign_j * (x_j(k+1) - bound_j)
#'
#' on that path: positive iff every constraint is slack in period k+1, so the
#' spell can end at k.  For k = T, x(T+1) is the slack rule at s_T.
#'
#' @param k            Integer spell length (0 <= k <= T)
#' @param ctx          Context from .obc_pwl_context()
#' @param spell_regime Integer bitfield binding during the spell
#' @param shock_seq    n_exo x T shock matrix
#' @param state_init   Initial state (NULL = zero)
#' @return Scalar minimum gap (positive = all slack at k+1)
#' @noRd
boehl_end_of_spell_gap <- function(k, ctx, spell_regime, shock_seq,
                                   state_init = NULL) {
  n_T <- ncol(shock_seq)
  rp  <- c(rep(as.integer(spell_regime), k), integer(n_T - k))
  sim <- .obc_pwl_forward(ctx, .obc_pwl_rules(ctx, rp), shock_seq, state_init)
  x_next <- if (k < n_T) sim$paths[, k + 1L]
            else drop(ctx$ghx %*% sim$states[, n_T])
  min(ctx$sgn * (x_next[ctx$var] - ctx$bnd))
}


#' Find the spell duration k_star via bisection (Boehl 2022 §3)
#'
#' For an IRF with a single binding spell [1, k_star] the end-of-spell gap
#' f(k) = boehl_end_of_spell_gap(k, ...) is taken to be monotone
#' non-decreasing in k.  Bisection finds k_star = min{k : f(k) >= 0} in
#' ceiling(log2(T)) evaluations; each evaluation simulates the candidate path
#' with its own time-varying rules.
#'
#' Special cases:
#'   f(0) >= 0: no binding at all → k_star = 0 (all slack)
#'   f(T) < 0: the bisection bracket is invalid. f need NOT be monotone: with
#'     two or more constraints, holding \code{spell_regime} for long can push
#'     a constraint OUTSIDE it through its bound (nk_2obc: y pinned at its
#'     floor for 4+ periods drives r below the ZLB, so f(1..3) > 0 but
#'     f(4..T) < 0). The smallest k with f(k) >= 0 is then found by a linear
#'     scan (k_star = T only if there is none). Before 0.9.4.77 this case
#'     returned k_star = T outright, and the caller reported non-convergence
#'     for a shock whose path IS a single spell (k_star = 1).
#'
#' @param shock_seq    n_exo x T shock matrix (unit shock at t=1 for IRF)
#' @param ctx          Context from .obc_pwl_context()
#' @param spell_regime Integer bitfield binding during the spell
#' @param state_init   Numeric vector (length n_state) initial state; NULL = zero
#' @return Named list:
#'   $k_star       -- integer spell length (0 = all slack, T = binding throughout)
#'   $gap_at_kstar -- scalar gap value at k_star (should be >= 0)
#'   $n_bisect     -- integer: bisection steps taken
#' @noRd
boehl_find_spell_duration <- function(shock_seq, ctx, spell_regime,
                                       state_init = NULL) {
  n_T <- ncol(shock_seq)
  gap <- function(k)
    boehl_end_of_spell_gap(k, ctx, spell_regime, shock_seq, state_init)

  g0 <- gap(0L)
  if (g0 >= 0)
    return(list(k_star = 0L, gap_at_kstar = g0, n_bisect = 0L))
  gap_T <- gap(n_T)
  if (gap_T < 0) {
    ## Invalid bracket (non-monotone f): first crossing by linear scan.
    for (k in seq_len(n_T - 1L)) {
      g_k <- gap(k)
      if (g_k >= 0) return(list(k_star = k, gap_at_kstar = g_k, n_bisect = k))
    }
    return(list(k_star = n_T, gap_at_kstar = gap_T, n_bisect = n_T - 1L))
  }

  # Invariant: f(lo) < 0, f(hi) >= 0.
  lo <- 0L
  hi <- n_T
  n_bisect <- 0L
  g_hi <- gap_T
  while (hi - lo > 1L) {
    n_bisect <- n_bisect + 1L
    mid <- (lo + hi) %/% 2L
    g_mid <- gap(mid)
    if (g_mid < 0) {
      lo <- mid
    } else {
      hi <- mid
      g_hi <- g_mid
    }
  }
  list(k_star = hi, gap_at_kstar = g_hi, n_bisect = n_bisect)
}


#' Validate user-supplied regime inputs of the OBC path solvers
#'
#' A regime is a bitfield over the \code{n_spec} constraints (bit j-1 set =
#' constraint j binds), so a code lies in 0..2^n_spec - 1; a spell regime must
#' bind at least one constraint. The internals silently REPLACED a
#' wrong-length \code{state_init} by the zero state and dropped a wrong-length
#' \code{regime_path_init}; an \code{NA} warm start came back as an all-NA
#' regime path with \code{converged = TRUE}, and an out-of-range or NA
#' \code{spell_regime} (a constraint that does not exist) was simulated.
#'
#' @param what Caller name for the messages.
#' @param n_spec Number of constraints.
#' @param n_T Horizon (NULL: not checked).
#' @param n_state Number of state variables (NULL: not checked).
#' @param spell_regime,regime_path_init,state_init The arguments (NULL = unset).
#' @return \code{invisible(TRUE)}; aborts with class
#'   \code{dynhr_error_bad_argument} otherwise.
#' @noRd
.obc_check_regime_inputs <- function(what, n_spec, n_T = NULL, n_state = NULL,
                                     spell_regime = NULL,
                                     regime_path_init = NULL,
                                     state_init = NULL) {
  max_code <- 2^n_spec - 1
  is_code <- function(v, lo)
    (is.numeric(v) || is.logical(v)) && length(v) > 0L && !anyNA(v) &&
      all(is.finite(v)) &&
      all(v == round(v)) && all(v >= lo) && all(v <= max_code)
  if (!is.null(spell_regime) &&
      (length(spell_regime) != 1L || !is_code(spell_regime, 1)))
    .dynhr_abort(sprintf(paste0(
      "%s: `spell_regime` must be ONE integer bitfield in 1..%d (bit j-1 set = ",
      "constraint j binds; %d constraint(s))."), what, as.integer(max_code),
      as.integer(n_spec)), class = "dynhr_error_bad_argument")
  if (!is.null(regime_path_init) &&
      (!is_code(regime_path_init, 0) ||
       (!is.null(n_T) && length(regime_path_init) != n_T)))
    .dynhr_abort(sprintf(paste0(
      "%s: `regime_path_init` must be an integer bitfield vector of length ",
      "%s (one code per period) with codes in 0..%d and no NA."), what,
      if (is.null(n_T)) "T" else as.character(n_T), as.integer(max_code)),
      class = "dynhr_error_bad_argument")
  if (!is.null(state_init) &&
      (!is.numeric(state_init) || any(!is.finite(state_init)) ||
       (!is.null(n_state) && length(state_init) != n_state)))
    .dynhr_abort(sprintf(paste0(
      "%s: `state_init` must be a finite numeric vector of length %s ",
      "(the number of state variables)."), what,
      if (is.null(n_state)) "n_state" else as.character(n_state)),
      class = "dynhr_error_bad_argument")
  invisible(TRUE)
}


# =============================================================================
# Main solver (iterative path + spell-duration fast path)
# =============================================================================

#' Boehl (2022) / OccBin OBC path solver
#'
#' Finds the binding/slack regime path and the piecewise-linear path for a
#' deterministic shock sequence.  A regime path is simulated with the
#' time-varying rules of the OccBin backward recursion (the rule of a period
#' depends on the regimes ahead of it) and accepted when it reproduces itself
#' on the constrained path: every slack period satisfies every bound and
#' every binding period has a multiplier of the right sign (the residual F =
#' lhs - rhs of the constraint's tagged equation, F >= 0 at a lower bound, F
#' <= 0 at an upper bound, as in Dynare's MCP convention and
#' \code{\link{pf_newton_solve}}).
#'
#' Shocks are surprises: at period 1 and at every period with a non-zero
#' shock the regime path expected from then on is re-solved from the current
#' state (Dynare's \code{occbin_solver}).  With a single shock in period 1
#' the result is the exact perfect-foresight piecewise-linear path (slack
#' after the horizon \code{nrow(shock_seq)}).
#'
#' Two dispatch paths:
#'
#' \strong{single_spell = FALSE} (default): OccBin guess-and-verify.  Starting
#'   from \code{regime_path_init} (all slack by default), every period is
#'   re-checked on the constrained path and the regime path updated until it
#'   reproduces itself.  Handles multi-spell, staggered and interacting
#'   constraints (a bound violated only because another one binds).
#'
#' \strong{single_spell = TRUE}: spell-duration root-finder (Boehl 2022 §3).
#'   Parameterises the regime path as a single contiguous block [1, k_star]
#'   of \code{spell_regime} followed by permanent slack and bisects for
#'   k_star.  Intended for single-spell ZLB / ELB IRFs.  When no spell length
#'   gives an equilibrium (the constraints bind in different periods, a
#'   constraint outside \code{spell_regime} binds, or the regime re-binds) it
#'   falls through to the iterative path, announced by a message of class
#'   \code{dynhr_obc_spell_fallthrough}; \code{attr(result,
#'   "obc_path_solver")} is \code{"spell"}, \code{"spell_fallthrough"} or
#'   \code{"iterative"}.
#'
#' @param shock_seq        T x n_exo numeric matrix of structural shocks
#'                         (rows = periods, columns = shocks; name the
#'                         columns by \code{dr_slack$exo_names} -- a matrix
#'                         named by shock on its ROWS, the old n_exo x T
#'                         orientation, is transposed with a message, an
#'                         unnamed n_exo x T or square one is an error)
#' @param dr_slack         Slack-regime DecisionRules (from solve_perturbation)
#' @param sys              System matrices (from extract_system_matrices_fast)
#' @param specs            OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param obs_idx          Unused (kept for interface compatibility).
#' @param single_spell     Logical: use the spell-duration fast path
#'                         (default FALSE → iterative path)
#' @param spell_regime     Integer bitfield active during the binding spell;
#'                         used only when single_spell = TRUE.
#'                         NULL → 2^n_spec - 1 (all specs bind simultaneously)
#' @param max_iter         Maximum guess-and-verify iterations per re-solve
#'                         (ignored when single_spell = TRUE; default 50L)
#' @param tol              Relative complementarity tolerance (default 1e-8):
#'                         a bound is violated when the slack
#'                         \eqn{s (x - b)} is below \eqn{-tol (|x| + |b|)}, a
#'                         multiplier (the residual of the tagged equation)
#'                         has the wrong sign when below -tol times the sum of
#'                         the absolute terms of that residual.  Scale-free:
#'                         rescaling the model's units leaves the regimes
#'                         unchanged (an absolute tolerance previously,
#'                         2026-09).
#' @param regime_path_init Integer vector (length T) warm start for the
#'                         iterative path; ignored when single_spell = TRUE
#' @param state_init       Numeric vector (length n_state) initial model state
#'                         (deviations); NULL → zero initial state
#' @return Named list with:
#'   $regime_path -- integer vector (length T): final bitfield regime per period
#'   $paths       -- T x n_endo matrix of simulated paths (deviations;
#'                   rows = periods)
#'   $fb_norm     -- scalar: Fischer-Burmeister residual norm of the returned
#'                   path (zero to rounding at an equilibrium)
#'   $converged   -- logical
#'   $n_iter      -- integer: guess-and-verify iterations (the largest
#'                   over the re-solves) or bisection steps (spell)
#' @export
boehl_solve_regime_path <- function(shock_seq, dr_slack, sys, specs,
                                     obs_idx          = NULL,
                                     single_spell     = FALSE,
                                     spell_regime     = NULL,
                                     max_iter         = 50L,
                                     tol              = 1e-8,
                                     regime_path_init = NULL,
                                     state_init       = NULL) {
  ## Rows = periods (R/time-orientation.R); the solver internals run on the
  ## n_exo x T transpose and the paths are returned T x n_endo.
  shock_seq <- t(.orient_time_matrix(shock_seq, dr_slack$exo_names,
                                     "shock_seq", "boehl_solve_regime_path"))
  n_T    <- ncol(shock_seq)
  n_spec <- length(specs)
  ctx    <- .obc_pwl_context(sys, dr_slack, specs)
  .obc_check_regime_inputs("boehl_solve_regime_path", n_spec, n_T,
                           ctx$n_state, spell_regime = spell_regime,
                           regime_path_init = regime_path_init,
                           state_init = state_init)

  # ------------------------------------------------------------------
  # FAST PATH: spell-duration root-finder (Boehl 2022 §3)
  # ------------------------------------------------------------------
  if (single_spell) {
    if (is.null(spell_regime))
      spell_regime <- 2L^n_spec - 1L
    spell_regime <- as.integer(spell_regime)

    spell  <- boehl_find_spell_duration(shock_seq, ctx, spell_regime,
                                        state_init = state_init)
    k_star <- spell$k_star
    regime_path <- c(rep(spell_regime, k_star), integer(n_T - k_star))

    rules <- .obc_pwl_rules(ctx, regime_path)
    sim   <- .obc_pwl_forward(ctx, rules, shock_seq, state_init)
    chk   <- .obc_pwl_check(ctx, rules, sim, shock_seq, regime_path,
                            state_init, tol)

    if (identical(chk$regime_path, regime_path)) {
      out <- list(
        regime_path = regime_path,
        paths       = t(sim$paths),
        fb_norm     = sqrt(sum(chk$phi^2)),
        converged   = TRUE,
        n_iter      = spell$n_bisect
      )
      attr(out, "obc_path_solver") <- "spell"
      return(out)
    }
    ## The single-spell assumption does not hold for this shock (no k makes
    ## "spell_regime on [1, k], slack after" an equilibrium: the constraints
    ## bind in different periods, a constraint outside spell_regime binds, or
    ## the regime re-binds). Fall through, announced, to the general solver
    ## below rather than returning a non-equilibrium path.
    .dynhr_inform(sprintf(paste0(
      "boehl_solve_regime_path: the single-spell fast path (regime %d on ",
      "periods 1..k, slack after) has no equilibrium for this shock sequence ",
      "(best k = %d); falling through to the general OccBin guess-and-verify ",
      "solver."), spell_regime, k_star),
      class = "dynhr_obc_spell_fallthrough")
  }

  # ------------------------------------------------------------------
  # ITERATIVE PATH: OccBin guess-and-verify, re-solved at each surprise
  # ------------------------------------------------------------------
  inner <- function(ctx, e, s, ini)
    .obc_pwl_solve(ctx, e, s, ini, max_iter = max_iter, tol = tol)
  res <- .obc_pwl_solve_surprise(ctx, shock_seq, state_init,
                                 regime_path_init, inner)
  out <- list(
    regime_path = res$regime_path,
    paths       = t(res$paths),
    fb_norm     = sqrt(sum(res$phi^2)),
    converged   = res$converged,
    n_iter      = res$n_iter
  )
  attr(out, "obc_path_solver") <- if (single_spell) "spell_fallthrough"
                                  else "iterative"
  out
}


# =============================================================================
# OBC-aware IRF computation
# =============================================================================

#' Compute impulse response functions for an OBC model
#'
#' Extends the standard linear compute_irfs() with OBC-regime switching.
#' For each shock, constructs an impulse at t=1 (zero thereafter), runs the
#' chosen OBC solver to identify the binding regime path, simulates the full
#' n_endo x T path, and packages it as an IRFCollection with the same
#' structure as compute_irfs().
#'
#' The impulse of shock k is the k-th column of the lower Cholesky factor of
#' the shock covariance times \code{shock_size} -- the same rule as
#' \code{\link{compute_irfs}()}, so correlated shocks propagate and, when no
#' constraint binds, the OBC IRF equals the linear IRF for every solver.
#'
#' Four solver options:
#' \describe{
#'   \item{\code{"boehl"}}{Spell-duration fast path (single_spell = TRUE).
#'     Bisection over the spell length.  Valid when the binding spell is a
#'     single contiguous block starting at t=1 in which the constraints of
#'     \code{spell_regime} bind; otherwise it falls through (announced) to
#'     the \code{"occbin"} guess-and-verify path.}
#'   \item{\code{"occbin"}}{OccBin guess-and-verify (single_spell = FALSE).
#'     Handles multi-spell, multi-constraint, and escape-and-rebind
#'     scenarios.}
#'   \item{\code{"lcp"}}{LCP solver (\code{\link{solve_obc_lcp}}, Lemke on
#'     the stacked multiplier LCP).  Same equilibrium as \code{"occbin"}.}
#'   \item{\code{"mcp"}}{MCP semi-smooth Newton solver (\code{\link{mcp_solve_path}}).
#'     Uses the Fischer-Burmeister FB reformulation on the full nonlinear model.
#'     Works for both linear and nonlinear models.  On linear models, matches
#'     the Boehl and LCP solvers to within numerical tolerance.  The
#'     perfect-foresight problem is solved on an internal horizon of
#'     \code{max(2 * n_periods, n_periods + 40)} periods (return to steady
#'     state imposed at its end) and truncated to \code{n_periods}, so the
#'     returned periods do not depend on \code{n_periods}.}
#' }
#'
#' \strong{Default solver.}  With a single constraint (or an explicit
#' \code{spell_regime}) the default is \code{"boehl"}.  With two or more
#' constraints and \code{spell_regime = NULL}, the boehl fast path would
#' assume every constraint binds in the same contiguous block, which is wrong
#' whenever the constraints bind in different periods (or one binds and the
#' other does not); the default is then \code{"occbin"}.  \code{"boehl"} stays
#' available when requested explicitly.
#'
#' \strong{Verification.}  Every returned path is checked against every
#' constraint; when a bound is violated beyond tolerance, or the solver did
#' not converge, a warning of class
#' \code{dynhr_warning_obc_path_infeasible} is emitted.
#'
#' @param dr_slack    Slack-regime DecisionRules (from solve_perturbation)
#' @param model       dynhr_mod (for shock standard deviations)
#' @param sys         System matrices (from extract_system_matrices_fast)
#' @param specs       OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param obs_idx     Integer vector: observable positions in endo vector;
#'                    NULL → all endo indices
#' @param n_periods   Number of IRF periods (default 40L)
#' @param shock_size  Shock size multiplier applied to each shock's impulse
#'                    (default 1: one-standard-deviation impulse)
#' @param params      Named numeric parameter vector used ONLY to select the
#'                    shock covariance of the impulse, exactly as in
#'                    \code{\link{compute_irfs}()}; the model is NOT
#'                    re-solved at \code{params} by any solver (the
#'                    dynamics are those of \code{dr_slack} / \code{sys};
#'                    the \code{"mcp"} path uses \code{model$param_values}).
#'                    NULL → \code{dr_slack$Sigma_e}, else the model's
#'                    shocks block.
#' @param obc_solver  Character: \code{"boehl"}, \code{"occbin"}, \code{"lcp"}
#'                    or \code{"mcp"}.  Default: \code{"boehl"} for one
#'                    constraint, \code{"occbin"} for two or more
#'                    constraints with \code{spell_regime = NULL} (see
#'                    Details).
#' @param compiled    Optional compiled model (needed by the \code{"mcp"} solver path)
#' @param spell_regime Integer bitfield for the binding block when
#'                    obc_solver = "boehl"; NULL → all specs simultaneously
#' @return Object of class "IRFCollection": a named list (one entry per shock)
#'   of n_periods x n_endo matrices.  Attributes: n_periods, endo_names,
#'   exo_names (same structure as compute_irfs()), plus
#'   \code{attr(., "converged")} (named logical vector, one per shock: the
#'   solver converged AND the path satisfies every constraint),
#'   \code{attr(., "regime_paths")} (named list of the solver's regime /
#'   active-set paths) and \code{attr(., "obc_solver")} (the solver used) --
#'   the convergence attributes of \code{\link{occbin_compute_irfs}()}.
#' @export
compute_irfs_obc <- function(dr_slack, model, sys, specs,
                              obs_idx      = NULL,
                              n_periods    = 40L,
                              shock_size   = 1,
                              params       = NULL,
                              obc_solver   = c("boehl", "occbin", "lcp", "mcp"),
                              spell_regime = NULL,
                              compiled     = NULL) {
  ## Multi-constraint default: the boehl fast path's "all constraints bind in
  ## one block" assumption is not a safe default (adversarial review
  ## 2026-11-11 F3), so default to the exact OccBin guess-and-verify.
  if (missing(obc_solver) && length(specs) > 1L && is.null(spell_regime))
    obc_solver <- "occbin"
  obc_solver <- match.arg(obc_solver)
  n_periods <- .check_irf_n_periods(n_periods, "compute_irfs_obc")
  .obc_check_regime_inputs("compute_irfs_obc", length(specs),
                           spell_regime = spell_regime)

  endo      <- dr_slack$endo_names
  exo       <- dr_slack$exo_names
  n_exo     <- length(exo)

  ## `params` selects the shock covariance only (the compute_irfs() rule).
  E_imp <- .irf_impulse_matrix(dr_slack, model, params, shock_size)

  irfs <- vector("list", n_exo)
  names(irfs) <- exo
  converged_v  <- stats::setNames(logical(n_exo), exo)
  regime_paths <- stats::setNames(vector("list", n_exo), exo)

  ## Internal perfect-foresight horizon for the MCP path: the terminal
  ## return-to-steady-state condition otherwise bites inside the window.
  n_solve <- if (obc_solver == "mcp") .pf_irf_horizon(n_periods) else n_periods

  # For MCP solver, build compiled model once and reuse across shocks.
  # Structural parameters are the model's own (params = shock covariance only).
  mcp_compiled <- NULL
  if (obc_solver == "mcp") {
    bridge <- .mcp_build_from_obc(model, specs, model$param_values,
                                  verbose = FALSE)
    mcp_compiled <- bridge$compiled
    mcp_specs_local <- bridge$mcp_specs
  }

  ctx <- .obc_pwl_context(sys, dr_slack, specs)
  tol <- 1e-6
  bad <- character(0)

  for (k in seq_along(exo)) {
    # Impulse at t=1 only (T x n_exo, rows = periods)
    shock_seq <- matrix(0, nrow = n_solve, ncol = n_exo,
                        dimnames = list(NULL, exo))
    shock_seq[1L, ] <- E_imp[, k]

    res <- if (obc_solver == "lcp") {
      solve_obc_lcp(shock_seq, dr_slack, sys, specs, obs_idx = obs_idx)
    } else if (obc_solver == "mcp") {
      # MCP solver — operates on the full nonlinear compiled model
      y0_num <- rep(0, length(endo))
      names(y0_num) <- endo
      mcp_res <- mcp_solve_path(
        compiled   = mcp_compiled,
        y0         = y0_num,
        y_ss       = y0_num,
        shock_path = shock_seq,
        params     = model$param_values,
        mcp_specs  = mcp_specs_local
      )
      list(paths = mcp_res$Y, regime_path = mcp_res$active_set,
           converged = mcp_res$converged)
    } else {
      boehl_solve_regime_path(
        shock_seq, dr_slack, sys, specs,
        obs_idx      = obs_idx,
        single_spell = obc_solver == "boehl",
        spell_regime = spell_regime
      )
    }

    irf_mat           <- res$paths[seq_len(n_periods), , drop = FALSE]
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))
    irfs[[k]]         <- irf_mat

    ## Verify the returned path against every constraint.
    x    <- irf_mat[, ctx$var, drop = FALSE]
    bnd  <- matrix(ctx$bnd, nrow(x), ncol(x), byrow = TRUE)
    sgn  <- matrix(ctx$sgn, nrow(x), ncol(x), byrow = TRUE)
    viol <- sgn * (x - bnd) < -tol * (abs(x) + abs(bnd))
    feasible <- !any(viol)
    converged_v[[k]]  <- isTRUE(res$converged) && feasible
    regime_paths[[k]] <- res$regime_path
    if (!converged_v[[k]]) bad <- c(bad, exo[k])
  }

  if (length(bad) > 0L)
    .dynhr_warn("compute_irfs_obc(obc_solver = \"", obc_solver, "\"): the ",
                "path for shock(s) ", paste(bad, collapse = ", "),
                " did not converge or violates a constraint beyond tolerance, ",
                "so it is not an equilibrium of the constrained model.",
                if (obc_solver == "boehl")
                  " The boehl fast path assumes one contiguous binding block of spell_regime starting at t = 1; use obc_solver = \"occbin\" or \"lcp\"."
                else "",
                " See attr(., \"converged\").",
                class = "dynhr_warning_obc_path_infeasible")

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")    <- n_periods
  attr(irfs, "endo_names")   <- endo
  attr(irfs, "exo_names")    <- exo
  attr(irfs, "converged")    <- converged_v
  attr(irfs, "regime_paths") <- regime_paths
  attr(irfs, "obc_solver")   <- obc_solver
  irfs
}
