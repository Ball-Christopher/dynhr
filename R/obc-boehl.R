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
##     by permanent slack.  converged = FALSE when the resulting path is not
##     an equilibrium (the assumption fails).
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
#'   f(T) < 0: binding through the end of sample → k_star = T
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
  if (gap_T < 0)
    return(list(k_star = n_T, gap_at_kstar = gap_T, n_bisect = 0L))

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
#' after the horizon \code{ncol(shock_seq)}).
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
#'   k_star.  Intended for single-spell ZLB / ELB IRFs; \code{converged} is
#'   FALSE when the resulting path is not an equilibrium.
#'
#' @param shock_seq        n_exo x T numeric matrix of structural shocks
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
#'   $paths       -- n_endo x T matrix of simulated paths (deviations)
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
  n_T    <- ncol(shock_seq)
  n_spec <- length(specs)
  ctx    <- .obc_pwl_context(sys, dr_slack, specs)

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

    return(list(
      regime_path = regime_path,
      paths       = sim$paths,
      fb_norm     = sqrt(sum(chk$phi^2)),
      converged   = identical(chk$regime_path, regime_path),
      n_iter      = spell$n_bisect
    ))
  }

  # ------------------------------------------------------------------
  # ITERATIVE PATH: OccBin guess-and-verify, re-solved at each surprise
  # ------------------------------------------------------------------
  inner <- function(ctx, e, s, ini)
    .obc_pwl_solve(ctx, e, s, ini, max_iter = max_iter, tol = tol)
  res <- .obc_pwl_solve_surprise(ctx, shock_seq, state_init,
                                 regime_path_init, inner)
  list(
    regime_path = res$regime_path,
    paths       = res$paths,
    fb_norm     = sqrt(sum(res$phi^2)),
    converged   = res$converged,
    n_iter      = res$n_iter
  )
}


# =============================================================================
# OBC-aware IRF computation
# =============================================================================

#' Compute impulse response functions for an OBC model
#'
#' Extends the standard linear compute_irfs() with OBC-regime switching.
#' For each shock, constructs a unit-impulse shock sequence (shock hits at
#' t=1, zero thereafter), runs the chosen OBC solver to identify the binding
#' regime path, simulates the full n_endo x T path, and packages it as an
#' IRFCollection with the same structure as compute_irfs().
#'
#' Four solver options:
#' \describe{
#'   \item{\code{"boehl"}}{Spell-duration fast path (single_spell = TRUE).
#'     Bisection over the spell length.  Valid when the binding spell is a
#'     single contiguous block starting at t=1; \code{converged} is FALSE
#'     otherwise.}
#'   \item{\code{"occbin"}}{OccBin guess-and-verify (single_spell = FALSE).
#'     Handles multi-spell, multi-constraint, and escape-and-rebind
#'     scenarios.}
#'   \item{\code{"lcp"}}{LCP solver (\code{\link{solve_obc_lcp}}, Lemke on
#'     the stacked multiplier LCP).  Same equilibrium as \code{"occbin"}.}
#'   \item{\code{"mcp"}}{MCP semi-smooth Newton solver (\code{\link{mcp_solve_path}}).
#'     Uses the Fischer-Burmeister FB reformulation on the full nonlinear model.
#'     Works for both linear and nonlinear models.  On linear models, matches
#'     the Boehl and LCP solvers to within numerical tolerance.}
#' }
#'
#' @param dr_slack    Slack-regime DecisionRules (from solve_perturbation)
#' @param model       dynhr_mod (for shock standard deviations)
#' @param sys         System matrices (from extract_system_matrices_fast)
#' @param specs       OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param obs_idx     Integer vector: observable positions in endo vector;
#'                    NULL → all endo indices
#' @param n_periods   Number of IRF periods (default 40L)
#' @param shock_size  Shock size multiplier applied to each shock's stderr
#'                    (default 1: one-standard-deviation impulse)
#' @param params      Named numeric parameter vector; NULL → model$param_values
#' @param obc_solver  Character: \code{"boehl"} (default), \code{"occbin"}, or \code{"lcp"}
#' @param compiled    Optional compiled model (needed by the \code{"mcp"} solver path)
#' @param spell_regime Integer bitfield for the binding block when
#'                    obc_solver = "boehl"; NULL → all specs simultaneously
#' @return Object of class "IRFCollection": a named list (one entry per shock)
#'   of n_periods x n_endo matrices.  Attributes: n_periods, endo_names,
#'   exo_names.  Same structure as compute_irfs().
#' @export
compute_irfs_obc <- function(dr_slack, model, sys, specs,
                              obs_idx      = NULL,
                              n_periods    = 40L,
                              shock_size   = 1,
                              params       = NULL,
                              obc_solver   = c("boehl", "occbin", "lcp", "mcp"),
                              spell_regime = NULL,
                              compiled     = NULL) {
  obc_solver <- match.arg(obc_solver)
  if (is.null(params)) params <- model$param_values

  endo      <- dr_slack$endo_names
  exo       <- dr_slack$exo_names
  n_exo     <- length(exo)

  shock_stderr <- .get_shock_stderr(model, exo, params)

  irfs <- vector("list", n_exo)
  names(irfs) <- exo

  # For MCP solver, build compiled model once and reuse across shocks
  mcp_compiled <- NULL
  if (obc_solver == "mcp") {
    bridge <- .mcp_build_from_obc(model, specs, params, verbose = FALSE)
    mcp_compiled <- bridge$compiled
    mcp_specs_local <- bridge$mcp_specs
    mcp_y_ss <- bridge$y_ss
  }

  for (k in seq_along(exo)) {
    # Unit impulse at t=1 only
    shock_seq <- matrix(0, nrow = n_exo, ncol = n_periods)
    shock_seq[k, 1L] <- shock_stderr[exo[k]] * shock_size

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
        shock_path = t(shock_seq),
        params     = params,
        mcp_specs  = mcp_specs_local
      )
      list(paths = t(mcp_res$Y), regime_path = mcp_res$active_set)
    } else {
      boehl_solve_regime_path(
        shock_seq, dr_slack, sys, specs,
        obs_idx      = obs_idx,
        single_spell = obc_solver == "boehl",
        spell_regime = spell_regime
      )
    }

    irf_mat           <- t(res$paths)   # n_periods x n_endo
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))
    irfs[[k]]         <- irf_mat
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")  <- n_periods
  attr(irfs, "endo_names") <- endo
  attr(irfs, "exo_names")  <- exo
  irfs
}
