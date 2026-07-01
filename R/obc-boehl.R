## R/obc-boehl.R
## --------------------------------------------------------------------------
## Boehl (2022) complementarity-based OBC simulation solver.
##
## Provides:
##   boehl_simulate()            -- forward simulation under a fixed regime path
##   boehl_fb_residual()         -- Fischer-Burmeister complementarity residuals
##   boehl_spell_trajectory()    -- state trajectory under full binding spell
##   boehl_end_of_spell_gap()    -- end-of-spell consistency check (scalar)
##   boehl_find_spell_duration() -- bisection root-finder: scalar k_star
##   boehl_solve_regime_path()   -- main solver; dispatches to spell-duration
##                                  fast path (single_spell=TRUE) or iterative
##   compute_irfs_obc()          -- OBC-aware IRF computation (obc_solver arg)
##
## TWO SOLVER PATHS
##
##   Iterative (single_spell = FALSE, default)
##     Replaces guess-and-verify with complementarity iteration on simulation.
##     O(n_iter × T) cost.  Handles arbitrary regime patterns: multi-spell,
##     staggered multi-constraint.  Typically 2–5 iterations.
##
##   Spell-duration root-finder (single_spell = TRUE) — Boehl 2022 §3
##     Parameterises the regime path by a single scalar k (spell length).
##     Pre-computes the full binding-regime state trajectory in one pass (O(T)),
##     then finds k_star via bisection (O(log T) evaluations, each O(1)).
##     Total cost: O(T + log T) ≈ O(T) vs O(n_iter × T) for the iterative path.
##
##     Assumption: exactly one contiguous binding block [1, k_star] followed
##     by permanent slack.  Valid for standard ZLB / ELB IRFs where the shock
##     hits at t=1 and the constraint relaxes monotonically.  Falls back to the
##     iterative solver if single_spell = FALSE.
##
## Reference: Boehl G. (2022) "Efficient solution and computation of models
##   with occasionally binding constraints." J. Econ. Dyn. Control 143, 104495.
## --------------------------------------------------------------------------


# =============================================================================
# Forward simulation
# =============================================================================

#' Forward simulation of an OBC model under a fixed regime path
#'
#' Simulates the model from an initial state (default: zero) for T periods,
#' using per-regime policy matrices from regime_cache.  Returns the full
#' endogenous variable paths and the state trajectory needed by the
#' boehl_fb_residual() binding check.
#'
#' State equation:  s_t = TT_r * s_{t-1} + RR_r * eps_t + c_state_r
#' Endo equation:   y_t = ghx_r * s_{t-1} + ghu_r * eps_t + c_full_r
#'
#' @param shock_seq    n_exo x T numeric matrix of structural shocks
#' @param dr_slack     Slack-regime DecisionRules (provides endo_names,
#'                     state_idx, and dimensions)
#' @param regime_cache R environment of per-regime policies built by
#'                     obc_ensure_policy() (keys: character regime index)
#' @param regime_path  Integer vector (length T): bitfield regime per period
#'                     (0 = all slack; bit j = 1 means spec j binds)
#' @param state_init   Numeric vector (length n_state) initial state;
#'                     NULL or omitted → zero initial state
#' @return List with:
#'   $paths  -- n_endo x T matrix (rownames = endo variable names)
#'   $states -- n_state x T matrix: s_t AFTER applying the period-t policy
#' @noRd
boehl_simulate <- function(shock_seq, dr_slack, regime_cache, regime_path,
                            state_init = NULL) {
  n_T     <- ncol(shock_seq)
  n_endo  <- nrow(dr_slack$ghx)
  si      <- dr_slack$state_idx
  n_state <- length(si)

  paths  <- matrix(NA_real_, n_endo, n_T)
  states <- matrix(NA_real_, n_state, n_T)
  rownames(paths) <- dr_slack$endo_names

  s <- if (!is.null(state_init) && length(state_init) == n_state)
         as.numeric(state_init)
       else numeric(n_state)

  for (t in seq_len(n_T)) {
    pol <- get(as.character(regime_path[t]), envir = regime_cache, inherits = FALSE)
    eps <- shock_seq[, t]

    paths[, t]  <- drop(pol$dr$ghx %*% s) + drop(pol$dr$ghu %*% eps) + pol$c_full
    s           <- drop(pol$TT  %*% s) + drop(pol$RR  %*% eps) + pol$c_state
    states[, t] <- s
  }

  list(paths = paths, states = states)
}


# =============================================================================
# Fischer-Burmeister complementarity residual
# =============================================================================

#' Compute Fischer-Burmeister complementarity residuals
#'
#' For each period t and constraint j evaluates
#'   φ(a_jt, b_jt) = a_jt + b_jt - sqrt(a_jt^2 + b_jt^2)
#' where:
#'   a_jt = max(0, ±(bound_j - x_j^slack(t)))  -- shadow-multiplier proxy
#'   b_jt = max(0, ±(x_j(t)  - bound_j))       -- constraint-slack proxy
#' with sign convention adjusted for lower ">=" vs upper "<=" bounds.
#'
#' φ = 0 iff the complementarity condition holds: (a >= 0 AND b >= 0 AND ab = 0).
#' At a converged regime path, |φ_jt| is at most O(machine epsilon) for all
#' (t, j) because the simulation exactly pins binding variables at the bound.
#'
#' @param states    n_state x T matrix: state trajectory from boehl_simulate()
#' @param shock_seq n_exo x T matrix of structural shocks
#' @param paths     n_endo x T matrix: simulated endo paths from boehl_simulate()
#' @param specs     OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param dr_slack  Slack-regime DecisionRules
#' @return n_spec x T numeric matrix of FB residuals (zero at convergence)
#' @noRd
boehl_fb_residual <- function(states, shock_seq, paths, specs, dr_slack) {
  n_T     <- ncol(shock_seq)
  n_spec  <- length(specs)
  n_state <- nrow(states)
  phi     <- matrix(NA_real_, n_spec, n_T)
  s_prev  <- numeric(n_state)

  for (t in seq_len(n_T)) {
    eps <- shock_seq[, t]

    for (j in seq_len(n_spec)) {
      s <- specs[[j]]

      # Slack-policy prediction of x_j at period t
      x_slack <- sum(dr_slack$ghx[s$var_idx, ] * s_prev) +
                 sum(dr_slack$ghu[s$var_idx, ] * eps)

      # Current simulated value (= bound when binding, free when slack)
      x_cur <- paths[s$var_idx, t]

      if (s$op == ">") {
        # Lower bound: want x_j >= bound.
        #   a = max(0, bound - x_slack)  > 0 when slack policy would violate
        #   b = max(0, x_cur  - bound)   > 0 when not at bound (slack regime)
        a_jt <- max(0, s$bound - x_slack)
        b_jt <- max(0, x_cur  - s$bound)
      } else {
        # Upper bound: want x_j <= bound.
        #   a = max(0, x_slack - bound)  > 0 when slack policy would violate
        #   b = max(0, bound   - x_cur)  > 0 when not at bound (slack regime)
        a_jt <- max(0, x_slack - s$bound)
        b_jt <- max(0, s$bound - x_cur)
      }

      phi[j, t] <- a_jt + b_jt - sqrt(a_jt^2 + b_jt^2)
    }
    s_prev <- states[, t]
  }

  phi
}


# =============================================================================
# Spell-duration root-finder (Boehl 2022 §3)
# =============================================================================

#' Accumulate the model state under a fixed binding regime for T periods
#'
#' Computes s_t for t = 0, 1, ..., T by applying the binding-regime policy
#' recursively (one forward pass):
#'   s_t = TT_b * s_{t-1} + RR_b * eps_t + c_state_b
#'
#' The result is an n_state x (T+1) matrix: column 1 holds s_0 (= state_init),
#' column t+1 holds s_t.  All T+1 states are returned so that
#' boehl_end_of_spell_gap() can evaluate any k ∈ {0, ..., T} in O(1).
#'
#' @param shock_seq   n_exo x T numeric matrix of structural shocks
#' @param pol_bind    Binding-regime policy entry from obc_ensure_policy()
#'                    (list with $TT, $RR, $c_state)
#' @param state_init  Numeric vector (length n_state) initial state;
#'                    NULL → zero initial state
#' @return n_state x (T+1) numeric matrix; column k+1 = s_k
#' @noRd
boehl_spell_trajectory <- function(shock_seq, pol_bind, state_init = NULL) {
  n_T     <- ncol(shock_seq)
  n_state <- nrow(pol_bind$TT)

  states        <- matrix(NA_real_, n_state, n_T + 1L)
  s             <- if (!is.null(state_init) && length(state_init) == n_state)
                     as.numeric(state_init)
                   else numeric(n_state)
  states[, 1L]  <- s   # s_0

  for (t in seq_len(n_T)) {
    s            <- drop(pol_bind$TT %*% s) +
                    drop(pol_bind$RR %*% shock_seq[, t]) +
                    pol_bind$c_state
    states[, t + 1L] <- s
  }
  states
}


#' Evaluate the end-of-spell consistency gap at period k
#'
#' After k periods of binding the model state is s_k (column k+1 of
#' states_bind).  The gap measures whether the slack policy at period k+1
#' satisfies ALL constraints simultaneously:
#'
#'   gap_j(k) = sign_j * (x_j^slack(k+1; s_k, eps_{k+1}) - bound_j)
#'
#' where sign_j = +1 for lower bounds (op ">") and −1 for upper bounds (op "<").
#' gap_j > 0 means constraint j would be slack under the slack policy at k+1.
#'
#' Returns min_j gap_j(k): positive iff ALL constraints are slack at k+1
#' (spell can end at k); negative iff at least one still binds (spell continues).
#'
#' For k = T, eps_{k+1} is taken as zero (terminal check beyond the sample).
#'
#' @param k            Integer spell length to evaluate (0 ≤ k ≤ T)
#' @param states_bind  n_state x (T+1) matrix from boehl_spell_trajectory()
#' @param shock_seq    n_exo x T shock matrix
#' @param specs        OBC spec list
#' @param dr_slack     Slack-regime DecisionRules
#' @return Scalar minimum gap (positive = all slack at k+1)
#' @noRd
boehl_end_of_spell_gap <- function(k, states_bind, shock_seq, specs, dr_slack) {
  n_T <- ncol(shock_seq)
  s_k <- states_bind[, k + 1L]                  # state after k binding periods

  # Shock at k+1; zero beyond the sample horizon
  eps_next <- if (k < n_T) shock_seq[, k + 1L] else numeric(nrow(shock_seq))

  min(vapply(specs, function(s) {
    x_slack <- sum(dr_slack$ghx[s$var_idx, ] * s_k) +
               sum(dr_slack$ghu[s$var_idx, ] * eps_next)
    if (s$op == ">") x_slack - s$bound else s$bound - x_slack
  }, numeric(1)))
}


#' Find the spell duration k_star via bisection (Boehl 2022 §3)
#'
#' For an IRF with a single binding spell [1, k_star] the end-of-spell gap
#' f(k) = boehl_end_of_spell_gap(k, ...) is monotone non-decreasing in k
#' (for well-behaved DSGE models the constraint relaxes over time after a shock).
#' Bisection finds k_star = min{k : f(k) ≥ 0} in ⌈log₂(T)⌉ evaluations.
#'
#' Special cases:
#'   f(0) ≥ 0: no binding at all → k_star = 0 (all slack)
#'   f(T) < 0: binding through the end of sample → k_star = T
#'
#' Cost: one O(T) trajectory pre-computation, then O(log T) O(1) gap queries.
#'
#' @param shock_seq   n_exo x T shock matrix (unit shock at t=1 for IRF)
#' @param pol_bind    Binding-regime policy entry from obc_ensure_policy()
#' @param specs       OBC spec list
#' @param dr_slack    Slack-regime DecisionRules
#' @param state_init  Numeric vector (length n_state) initial state; NULL = zero
#' @return Named list:
#'   $k_star         -- integer spell length (0 = all slack, T = binding throughout)
#'   $states_bind    -- n_state x (T+1) trajectory (reusable for boehl_simulate)
#'   $gap_at_kstar   -- scalar gap value at k_star (should be ≥ 0)
#' @noRd
boehl_find_spell_duration <- function(shock_seq, pol_bind, specs, dr_slack,
                                       state_init = NULL) {
  n_T         <- ncol(shock_seq)
  states_bind <- boehl_spell_trajectory(shock_seq, pol_bind, state_init)

  # k=0: no binding periods at all
  if (boehl_end_of_spell_gap(0L, states_bind, shock_seq, specs, dr_slack) >= 0) {
    return(list(k_star = 0L, states_bind = states_bind, gap_at_kstar = 0))
  }

  # k=T: still binding at the sample horizon
  gap_T <- boehl_end_of_spell_gap(n_T, states_bind, shock_seq, specs, dr_slack)
  if (gap_T < 0) {
    return(list(k_star = n_T, states_bind = states_bind, gap_at_kstar = gap_T))
  }

  # Bisection: find smallest k in [1, T] such that f(k) >= 0.
  # Invariant: f(lo) < 0, f(hi) >= 0.
  lo <- 0L
  hi <- n_T

  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (boehl_end_of_spell_gap(mid, states_bind, shock_seq, specs, dr_slack) < 0)
      lo <- mid
    else
      hi <- mid
  }

  list(
    k_star      = hi,
    states_bind = states_bind,
    gap_at_kstar = boehl_end_of_spell_gap(hi, states_bind, shock_seq, specs, dr_slack)
  )
}


# =============================================================================
# Main solver (iterative path + spell-duration fast path)
# =============================================================================

#' Boehl (2022) OBC simulation solver
#'
#' Finds the binding/slack regime path for a given deterministic shock sequence.
#' Two dispatch paths:
#'
#' \strong{single_spell = FALSE} (default): iterative complementarity path.
#'   At each iteration the model is forward-simulated, the slack-policy
#'   prediction at every period determines which constraints bind, and the
#'   regime path is updated.  Repeats until the path stabilises.  Handles
#'   multi-spell, staggered, and multi-constraint scenarios.
#'
#' \strong{single_spell = TRUE}: spell-duration root-finder (Boehl 2022 §3).
#'   Parameterises the regime path as a single contiguous block [1, k_star]
#'   followed by permanent slack.  Pre-computes the full binding-regime state
#'   trajectory in one pass, then bisects for k_star in O(log T) steps.
#'   Cost: O(T + log T) vs O(n_iter × T) for the iterative path.
#'   Assumption: the shock hits at t=1 and the constraint relaxes
#'   monotonically — the standard IRF convention.  Intended for single-spell
#'   ZLB / ELB IRF analysis.  spell_regime controls which constraint bitfield
#'   is active during the binding block (default: all specs simultaneously).
#'
#' @param shock_seq        n_exo x T numeric matrix of structural shocks
#' @param dr_slack         Slack-regime DecisionRules (from solve_perturbation)
#' @param sys              System matrices (from extract_system_matrices_fast)
#' @param specs            OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param obs_idx          Integer vector: observable positions in endo vector.
#'                         NULL → all endo indices (safe for simulation-only use)
#' @param single_spell     Logical: use the O(T + log T) spell-duration fast path
#'                         (default FALSE → iterative path)
#' @param spell_regime     Integer bitfield active during the binding spell;
#'                         used only when single_spell = TRUE.
#'                         NULL → 2^n_spec - 1 (all specs bind simultaneously)
#' @param max_iter         Maximum outer iterations for the iterative path (ignored
#'                         when single_spell = TRUE; default 50L)
#' @param tol              FB-norm tolerance; also used as the end-of-spell gap
#'                         tolerance when single_spell = TRUE (default 1e-8)
#' @param regime_path_init Integer vector (length T) for warm-starting the
#'                         iterative path; ignored when single_spell = TRUE
#' @param state_init       Numeric vector (length n_state) initial model state;
#'                         NULL → zero initial state (standard IRF convention)
#' @return Named list with:
#'   $regime_path -- integer vector (length T): final bitfield regime per period
#'   $paths       -- n_endo x T matrix of simulated endogenous variable paths
#'   $fb_norm     -- scalar: FB residual norm (iterative) or end-of-spell gap
#'                   magnitude (spell-duration)
#'   $converged   -- logical
#'   $n_iter      -- integer: iterations (iterative) or bisection steps (spell)
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

  if (is.null(obs_idx)) obs_idx <- seq_len(nrow(dr_slack$ghx))

  # Lazy policy cache — always seed with the slack policy
  regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)

  # ------------------------------------------------------------------
  # FAST PATH: spell-duration root-finder (Boehl 2022 §3)
  # ------------------------------------------------------------------
  if (single_spell) {
    # Determine which regime bitfield is active during the binding spell.
    # Default: all n_spec specs bind simultaneously = 2^n_spec - 1.
    if (is.null(spell_regime))
      spell_regime <- 2L^n_spec - 1L

    # Build the binding policy for spell_regime (lazy)
    obc_ensure_policy(as.integer(spell_regime), regime_cache,
                      sys, dr_slack, specs, obs_idx)
    pol_bind <- get(as.character(spell_regime), envir = regime_cache,
                    inherits = FALSE)

    # Find k_star via bisection
    spell <- boehl_find_spell_duration(shock_seq, pol_bind, specs, dr_slack,
                                        state_init = state_init)
    k_star <- spell$k_star

    # Construct regime path: binding for periods 1..k_star, slack thereafter
    regime_path <- c(rep(as.integer(spell_regime), k_star),
                     integer(n_T - k_star))

    # Forward simulate under this regime path
    sim <- boehl_simulate(shock_seq, dr_slack, regime_cache, regime_path,
                          state_init = state_init)

    return(list(
      regime_path = regime_path,
      paths       = sim$paths,
      fb_norm     = abs(spell$gap_at_kstar),
      converged   = spell$gap_at_kstar >= -tol,
      n_iter      = ceiling(log2(max(n_T, 1L)))  # bisection steps taken
    ))
  }

  # ------------------------------------------------------------------
  # ITERATIVE PATH: complementarity iteration
  # ------------------------------------------------------------------
  regime_path <- if (!is.null(regime_path_init) &&
                     length(regime_path_init) == n_T)
                   as.integer(regime_path_init)
                 else integer(n_T)

  for (r in unique(regime_path[regime_path != 0L])) {
    if (!exists(as.character(r), envir = regime_cache, inherits = FALSE))
      obc_ensure_policy(r, regime_cache, sys, dr_slack, specs, obs_idx)
  }

  n_state <- length(dr_slack$state_idx)
  sim     <- NULL

  for (iter in seq_len(max_iter)) {
    sim <- boehl_simulate(shock_seq, dr_slack, regime_cache, regime_path,
                          state_init = state_init)

    new_regime <- integer(n_T)
    s_prev     <- if (!is.null(state_init) && length(state_init) == n_state)
                    as.numeric(state_init)
                  else numeric(n_state)

    for (t in seq_len(n_T)) {
      eps <- shock_seq[, t]
      bind_flags <- vapply(specs, function(s) {
        x_slack <- sum(dr_slack$ghx[s$var_idx, ] * s_prev) +
                   sum(dr_slack$ghu[s$var_idx, ] * eps)
        if (s$op == ">") x_slack < s$bound else x_slack > s$bound
      }, logical(1))
      new_regime[t] <- obc_regime_idx(bind_flags)
      s_prev        <- sim$states[, t]
    }

    for (r in unique(new_regime[new_regime != 0L])) {
      if (!exists(as.character(r), envir = regime_cache, inherits = FALSE))
        obc_ensure_policy(r, regime_cache, sys, dr_slack, specs, obs_idx)
    }

    if (identical(new_regime, regime_path)) {
      phi     <- boehl_fb_residual(sim$states, shock_seq, sim$paths, specs, dr_slack)
      fb_norm <- sqrt(sum(phi^2))
      return(list(
        regime_path = regime_path,
        paths       = sim$paths,
        fb_norm     = fb_norm,
        converged   = TRUE,
        n_iter      = iter
      ))
    }

    regime_path <- new_regime
  }

  phi     <- boehl_fb_residual(sim$states, shock_seq, sim$paths, specs, dr_slack)
  fb_norm <- sqrt(sum(phi^2))
  list(
    regime_path = regime_path,
    paths       = sim$paths,
    fb_norm     = fb_norm,
    converged   = FALSE,
    n_iter      = max_iter
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
#'     O(T + log T) cost per shock.  Valid when the binding spell is a single
#'     contiguous block starting at t=1 — the standard IRF assumption.}
#'   \item{\code{"occbin"}}{Iterative complementarity solver (single_spell = FALSE).
#'     O(n_iter × T) cost per shock.  Handles multi-spell, multi-constraint,
#'     and escape-and-rebind scenarios.}
#'   \item{\code{"lcp"}}{LCP-Newton solver (solve_obc_lcp).
#'     Builds a T*n_spec × T*n_spec finite-difference Jacobian and takes Newton
#'     steps on the stacked FB system.  Handles multi-spell and coupled-regime
#'     scenarios where the iterative path may cycle.}
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
