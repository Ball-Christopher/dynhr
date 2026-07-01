## R/occbin-nonlinear-solve.R
## --------------------------------------------------------------------------
## Nonlinear OccBin path solver with regime-switching outer loop.
##
## Provides:
##   occbin_build_regime_fns()   -- precompute residual/jacobian selectors per regime
##   occbin_solve_path()         -- main solver: stacked Newton + regime outer loop
##   occbin_compute_irfs()       -- OBC IRF computation using the nonlinear solver
##
## ALGORITHM
##
##   For a model with k bind/relax constraint pairs (n_eq > n_endo), the solver:
##
##   1. Given a regime_path (integer bitfield per period), select the correct
##      equation variant for each endogenous variable per period.
##   2. Stack T × n_endo equations and solve via Newton with a block-tridiagonal
##      sparse linear solve (Matrix::sparseQR / Eigen via RcppArmadillo).
##   3. After Newton convergence, check complementarity: for each constraint j
##      and period t, does the current path violate the bound?  Flip periods that
##      violate, update regime_path, and re-solve (outer loop).
##   4. Repeat until no regime flips occur.
##
## RELATION TO EXISTING SOLVERS
##
##   pf_newton_solve()  — nonlinear, but requires n_eq == n_endo.  Handles OBCs
##                        via active-set modification of specific equation rows.
##   solve_obc_lcp()    — linear LCP Newton on regime assignment (binary Newton).
##   boehl_solve_regime_path() — iterative complementarity on linear models.
##
##   This solver fills the gap: nonlinear models with bind/relax pairs where
##   n_eq > n_endo and the regime determines which equation variant to use.
##
## SPARSE BLOCK-TRIDIAGONAL STRUCTURE
##
##   The stacked Jacobian J (T*n_endo × T*n_endo) has block-tridiagonal form:
##
##     [ A_1  B_1   0    0   ... ]
##     [ C_2  A_2  B_2   0   ... ]
##     [  0   C_3  A_3  B_3  ... ]
##     [  :    :    :    :    \   ]
##
##   where at period t:
##     C_t = ∂F_t/∂y_{t-1}  (n_endo × n_endo, lag Jacobian cols)
##     A_t = ∂F_t/∂y_t      (n_endo × n_endo, current Jacobian cols)
##     B_t = ∂F_t/∂y_{t+1}  (n_endo × n_endo, lead Jacobian cols)
##
##   The full Jacobian is assembled as a sparse dgCMatrix and solved via
##   Matrix::sparseQR or Matrix::solve().
##
## References:
##   Guerrieri & Iacoviello (2015) — OccBin
##   Adjemian & Juillard (2025) — Stochastic Extended Path
##   Holden & Paetz (2012) — OccBin stacked Newton
## --------------------------------------------------------------------------


# =============================================================================
# Precompute column mapping for the block-tridiagonal Jacobian
# =============================================================================

#' Precompute column metadata for block-tridiagonal Jacobian assembly
#'
#' For each equation, identifies which columns of the full Jacobian correspond
#' to lagged, current, and lead endogenous variables.  Used by the stacked
#' Newton solver to route Jacobian entries into the correct block of the
#' block-tridiagonal system.
#'
#' @param dyn compiled$dynamic (from build_dynamic_model)
#' @return List with:
#'   $lag_cols  — integer vector: column indices for y_{t-1} in the full Jacobian
#'   $cur_cols  — integer vector: column indices for y_t in the full Jacobian
#'   $lead_cols — integer vector: column indices for y_{t+1} in the full Jacobian
#'   $exo_cols  — integer vector: column indices for exogenous variables (ignored)
#'   $lag_var   — integer vector: endo variable index for each lag column
#'   $cur_var   — integer vector: endo variable index for each cur column
#'   $lead_var  — integer vector: endo variable index for each lead column
#'   $n_endo    — integer: number of endogenous variables
#'   $total_cols — integer: total number of columns in the full Jacobian
#' @noRd
.occbin_col_meta <- function(dyn) {
  cmap <- dyn$dyn_col_map
  n    <- nrow(cmap)

  lag_cols  <- integer(0)
  cur_cols  <- integer(0)
  lead_cols <- integer(0)
  lag_var   <- integer(0)
  cur_var   <- integer(0)
  lead_var  <- integer(0)
  exo_cols  <- integer(0)

  for (k in seq_len(n)) {
    nm <- cmap$name[k]
    ll <- cmap$lead_lag[k]
    cl <- cmap$col[k]

    if (nm %in% dyn$exo_names) {
      exo_cols <- c(exo_cols, cl)
      next
    }

    vi <- match(nm, dyn$endo_names)
    if (is.na(vi)) next

    if (ll < 0L) {
      lag_cols  <- c(lag_cols, cl)
      lag_var   <- c(lag_var, vi)
    } else if (ll == 0L) {
      cur_cols  <- c(cur_cols, cl)
      cur_var   <- c(cur_var, vi)
    } else {
      lead_cols <- c(lead_cols, cl)
      lead_var  <- c(lead_var, vi)
    }
  }

  list(
    lag_cols   = lag_cols,
    cur_cols   = cur_cols,
    lead_cols  = lead_cols,
    exo_cols   = exo_cols,
    lag_var    = lag_var,
    cur_var    = cur_var,
    lead_var   = lead_var,
    n_endo     = length(dyn$endo_names),
    total_cols = dyn$total_cols
  )
}


# =============================================================================
# Build regime-specific residual and Jacobian function selectors
# =============================================================================

#' Build residual/Jacobian selector functions for all regimes
#'
#' For each regime (0 .. 2^k-1), computes which equations are active and wraps
#' the compiled dynamic model's residuals_fn/jacobian_fn to select only those
#' equation rows.  This avoids recompiling the model for each regime.
#'
#' The selectors are simple wrappers: given (dy, params, y_ss), they call the
#' full residuals_fn/jacobian_fn and subset rows by the active equation indices.
#'
#' @param compiled     dynhr_compiled (from compile_model)
#' @param regime_map   Regime-to-equation mapping from occbin_build_regime_map()
#' @return List of length 2^k, each entry with:
#'   $eq_indices  — integer vector of active equation indices
#'   $res_fn      — function(dy, params, y_ss) -> numeric vector (n_endo)
#'   $jac_fn      — function(dy, params, y_ss) -> numeric matrix (n_endo × total_cols)
#' @noRd
.make_regime_fn <- function(eq_idx, res_fn, jac_fn) {
  # Force evaluation: eq_idx is captured by value in this function's scope
  force(eq_idx)
  force(res_fn)
  force(jac_fn)
  list(
    eq_indices = eq_idx,
    res_fn     = function(dy, params, y_ss) res_fn(dy, params, y_ss)[eq_idx],
    jac_fn     = function(dy, params, y_ss) jac_fn(dy, params, y_ss)[eq_idx, , drop = FALSE]
  )
}

#' Build residual/Jacobian selector functions for all regimes
#'
#' For each regime (0 .. 2^k-1), computes which equations are active and wraps
#' the compiled dynamic model's residuals_fn/jacobian_fn to select only those
#' equation rows.  This avoids recompiling the model for each regime.
#'
#' The selectors are simple wrappers: given (dy, params, y_ss), they call the
#' full residuals_fn/jacobian_fn and subset rows by the active equation indices.
#'
#' @param compiled     dynhr_compiled (from compile_model)
#' @param regime_map   Regime-to-equation mapping from occbin_build_regime_map()
#' @return List of length 2^k, each entry with:
#'   $eq_indices  — integer vector of active equation indices
#'   $res_fn      — function(dy, params, y_ss) -> numeric vector (n_endo)
#'   $jac_fn      — function(dy, params, y_ss) -> numeric matrix (n_endo × total_cols)
#' @noRd
occbin_build_regime_fns <- function(compiled, regime_map) {
  dyn     <- compiled$dynamic
  res_fn  <- dyn$residuals_fn
  jac_fn  <- dyn$jacobian_fn
  n_total <- length(regime_map)

  regime_fns <- vector("list", n_total)

  for (r in seq_len(n_total)) {
    eq_idx <- regime_map[[r]]$eq_indices

    # Factory function forces capture of eq_idx by value, not by reference.
    # In R, for-loop bodies share a single scope, so direct closures would
    # all reference the LAST value of eq_idx.
    regime_fns[[r]] <- .make_regime_fn(eq_idx, res_fn, jac_fn)
  }

  regime_fns
}


# =============================================================================
# Build the stacked sparse Jacobian from a path + regime
# =============================================================================

#' Build the stacked block-tridiagonal sparse Jacobian
#'
#' Assembles a (T*n_endo) × (T*n_endo) sparse dgCMatrix from the per-period
#' regime-specific Jacobians.  The matrix has block-tridiagonal structure:
#'   diagonal blocks A_t = ∂F_t/∂y_t
#'   sub-diagonal blocks C_t = ∂F_t/∂y_{t-1}
#'   super-diagonal blocks B_t = ∂F_t/∂y_{t+1}
#'
#' @param Y            T × n_endo numeric matrix: current path
#' @param y0_num       Length-n_endo initial state at t=0
#' @param y_ss_num     Length-n_endo steady state (terminal condition)
#' @param eps_mat      T × n_exo shock matrix
#' @param pf_meta      Column metadata from .pf_col_meta() (for .pf_make_dy)
#' @param col_meta     Column metadata from .occbin_col_meta() (for Jacobian routing)
#' @param regime_fns   Regime function list from occbin_build_regime_fns()
#' @param regime_path  Integer vector (length T): per-period regime bitmask
#' @param params       Named numeric parameter vector
#' @param T            Integer horizon
#' @param n_endo       Integer number of endogenous variables
#' @return List with:
#'   $J   — sparse dgCMatrix (T*n_endo × T*n_endo)
#'   $R   — numeric vector (T*n_endo): stacked residual
#' @noRd
.occbin_build_stacked_system <- function(Y, y0_num, y_ss_num, eps_mat,
                                          pf_meta, col_meta,
                                          regime_fns, regime_path,
                                          params, T, n_endo) {
  n_total <- T * n_endo

  # Pre-allocate residual vector
  R <- numeric(n_total)

  # Build the sparse Jacobian incrementally
  # We'll use the Matrix package to create a sparse matrix
  # For efficiency, collect triplets (i, j, value)
  i_triplet <- integer(0)
  j_triplet <- integer(0)
  v_triplet <- numeric(0)

  # Cache for dy vectors (reuse when regime is same across periods)
  # We still need to evaluate each period individually for residuals

  for (t in seq_len(T)) {
    row_off <- (t - 1L) * n_endo
    col_off <- (t - 1L) * n_endo

    # Build dy vector for this period (same as pf_newton_solve)
    dy <- .pf_make_dy(pf_meta, Y, y0_num, y_ss_num, eps_mat[t, ], t, T)

    # Get regime functions for this period
    r_idx <- regime_path[t] + 1L  # 1-based indexing
    rf    <- regime_fns[[r_idx]]

    # Compute residual and Jacobian.
    # Pass y_ss_num as the 'ss' argument so that steady_state(X) references
    # in bind/relax equations (e.g. iv = PHI*steady_state(iv)) resolve correctly.
    Rt <- rf$res_fn(dy, params, y_ss_num)
    Jt <- rf$jac_fn(dy, params, y_ss_num)

    # Store residual
    R[row_off + seq_len(n_endo)] <- Rt

    # Route Jacobian columns to global sparse entries
    # Lag columns -> block (t, t-1)
    for (k in seq_along(col_meta$lag_cols)) {
      dc <- col_meta$lag_cols[k]
      vi <- col_meta$lag_var[k]
      gcol <- if (t > 1L) (t - 2L) * n_endo + vi else NA_integer_
      if (is.na(gcol)) next
      for (eq in seq_len(n_endo)) {
        val <- Jt[eq, dc]
        if (val != 0) {
          i_triplet <- c(i_triplet, row_off + eq)
          j_triplet <- c(j_triplet, gcol)
          v_triplet <- c(v_triplet, val)
        }
      }
    }

    # Current columns -> block (t, t)
    for (k in seq_along(col_meta$cur_cols)) {
      dc <- col_meta$cur_cols[k]
      vi <- col_meta$cur_var[k]
      gcol <- col_off + vi
      for (eq in seq_len(n_endo)) {
        val <- Jt[eq, dc]
        if (val != 0) {
          i_triplet <- c(i_triplet, row_off + eq)
          j_triplet <- c(j_triplet, gcol)
          v_triplet <- c(v_triplet, val)
        }
      }
    }

    # Lead columns -> block (t, t+1)
    for (k in seq_along(col_meta$lead_cols)) {
      dc <- col_meta$lead_cols[k]
      vi <- col_meta$lead_var[k]
      gcol <- if (t < T) t * n_endo + vi else NA_integer_
      if (is.na(gcol)) next
      for (eq in seq_len(n_endo)) {
        val <- Jt[eq, dc]
        if (val != 0) {
          i_triplet <- c(i_triplet, row_off + eq)
          j_triplet <- c(j_triplet, gcol)
          v_triplet <- c(v_triplet, val)
        }
      }
    }
  }

  # Build sparse matrix from triplets
  J <- Matrix::sparseMatrix(
    i = i_triplet, j = j_triplet, x = v_triplet,
    dims = c(n_total, n_total)
  )

  list(J = J, R = R)
}


# =============================================================================
# Bind-condition evaluator
# =============================================================================

#' Evaluate OccBin bind conditions against a simulated path
#'
#' For each constraint j and period t, determines whether the constraint
#' SHOULD be binding given the current solution path Y.
#'
#' The bind direction is encoded in \code{cn$condition}:
#' \itemize{
#'   \item \code{"less"}    — constraint binds when \code{Y[t, var] < bound}
#'   \item \code{"greater"} — constraint binds when \code{Y[t, var] > bound}
#' }
#'
#' The bound is resolved from \code{cn$bound_expr} by:
#' \enumerate{
#'   \item Replacing all occurrences of \code{steady_state(X)} with the
#'         numeric steady-state value of variable \code{X} from
#'         \code{y_ss_num}.
#'   \item Evaluating the resulting expression in an environment that
#'         contains all model parameters (from \code{params}).
#' }
#'
#' Y must be in LEVELS (not deviations from steady state).
#'
#' @param Y           T x n_endo numeric matrix: simulated endogenous path in levels.
#' @param dyn         compiled$dynamic: contains endo_names for column look-up.
#' @param y_ss_num    Named numeric vector (length n_endo): steady-state values,
#'                    used to resolve steady_state(X) calls in the bound expression.
#' @param constraints Named list of constraint specs (each element from
#'                    occbin_parse_bind_relax()); must contain:
#'                    \code{var_name}, \code{condition}, \code{bound_expr}.
#' @param params      Named numeric vector of model parameters.
#' @return Logical matrix (n_constraints x T): TRUE where constraint j should
#'         bind at period t.
#' @noRd
.occbin_eval_bind_condition <- function(Y, dyn, y_ss_num, constraints, params) {
  n_constraints <- length(constraints)
  T             <- nrow(Y)

  should_bind <- matrix(FALSE, nrow = n_constraints, ncol = T)
  if (n_constraints == 0L) return(should_bind)

  # Build evaluation environment: params + endo SS values
  # (endo SS values are needed for steady_state(X) resolution)
  eval_env <- as.list(params)

  for (j in seq_len(n_constraints)) {
    cn <- if (is.list(constraints) && !is.null(names(constraints))) {
      constraints[[j]]
    } else {
      constraints[[j]]
    }

    # Column index of the constrained variable
    vi <- match(cn$var_name, dyn$endo_names)
    if (is.na(vi)) {
      warning(sprintf(
        ".occbin_eval_bind_condition: variable '%s' not found in endo_names; skipping.",
        cn$var_name))
      next
    }

    # Resolve the bound expression
    # Step 1: substitute steady_state(X) → numeric value
    bound_text <- if (!is.null(cn$bound_expr) && !is.na(cn$bound_expr)) {
      expr <- cn$bound_expr
      # Replace every steady_state(VAR) with the SS numeric value
      # Use gsub with a callback via Reduce over matched variable names
      m_all <- gregexpr("steady_state\\s*\\(\\s*(\\w+)\\s*\\)", expr, perl = TRUE)[[1]]
      if (m_all[1] != -1L) {
        ml   <- attr(m_all, "match.length")
        caps <- regmatches(expr,
                           gregexpr("steady_state\\s*\\(\\s*(\\w+)\\s*\\)",
                                    expr, perl = TRUE))[[1]]
        for (cap in caps) {
          vname <- sub("steady_state\\s*\\(\\s*(\\w+)\\s*\\)", "\\1", cap, perl = TRUE)
          vi_ss <- match(vname, dyn$endo_names)
          if (is.na(vi_ss)) {
            warning(sprintf(
              ".occbin_eval_bind_condition: steady_state(%s) not found; replacing with NA.",
              vname))
            ss_val <- NA_real_
          } else {
            ss_val <- y_ss_num[vi_ss]
          }
          expr <- gsub(cap, as.character(ss_val), expr, fixed = TRUE)
        }
      }
      expr
    } else {
      NA_character_
    }

    # Step 2: evaluate the resulting expression to get a numeric bound
    bound_val <- if (is.na(bound_text)) {
      NA_real_
    } else {
      tryCatch(
        eval(parse(text = bound_text), envir = eval_env),
        error = function(e) {
          warning(sprintf(
            ".occbin_eval_bind_condition: bound_expr '%s' evaluation failed: %s",
            bound_text, conditionMessage(e)))
          NA_real_
        }
      )
    }

    if (is.na(bound_val)) {
      warning(sprintf(
        ".occbin_eval_bind_condition: constraint '%s' has no evaluable bound; all slack.",
        cn$name %||% paste0("j=", j)))
      next
    }

    # Step 3: compare Y[t, var] to the bound using the bind direction
    var_path <- Y[, vi]

    should_bind[j, ] <- if (identical(cn$condition, "less")) {
      var_path < bound_val
    } else {
      var_path > bound_val
    }
  }

  should_bind
}


# =============================================================================
# Piecewise-linear solver (method = "pwlinear")
# =============================================================================

#' Build binding-regime decision rules for piecewise-linear OccBin
#'
#' Computes the slack DR and the binding DR (via OccBin terminal substitution)
#' for each constraint.  For single-constraint models (the common case), the
#' binding DR covers the fully-binding regime; for multi-constraint models each
#' constraint is solved independently under the assumption that only one binds
#' at a time (Guerrieri-Iacoviello 2015 approximation).
#'
#' The critical fix for dual-variable constraints (relax var ≠ bind var):
#'   eq_idx = match(cn$eq_var_name, endo_names)
#' where eq_var_name is the LHS variable of the RELAX equation.  The
#' extract_system_matrices() reorders Jacobian rows into declaration order,
#' so the row for the relax equation is at its declaration position.
#'
#' @param compiled dynhr_compiled
#' @param y_ss_num Named numeric (n_endo): steady-state levels
#' @param params   Named numeric: parameter vector
#' @return List with $dr_slack, $dr_bind, $c_bind, $state_idx, or NULL on error
#' @noRd
.occbin_build_pwlinear_drs <- function(compiled, y_ss_num, params) {
  dyn        <- compiled$dynamic
  endo_names <- dyn$endo_names
  n_endo     <- length(endo_names)

  pr        <- compiled$occbin$parse_result
  n_con     <- pr$n_constraints
  if (n_con == 0L) return(NULL)

  model <- compiled$model

  # Slack DR: standard first-order perturbation in the relax regime
  dr_slack <- tryCatch(
    solve_perturbation(model, compiled, y_ss_num, params),
    error = function(e) {
      warning(".occbin_build_pwlinear_drs: slack DR failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dr_slack)) return(NULL)

  state_idx <- dr_slack$state_idx

  # System matrices (Jacobians at SS, in declaration order)
  sys <- tryCatch(
    extract_system_matrices(compiled, y_ss_num, params),
    error = function(e) {
      warning(".occbin_build_pwlinear_drs: extract_system_matrices failed: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(sys)) return(NULL)

  # Binding DR: OccBin terminal-substitution for the fully-binding regime
  # (all constraints bind simultaneously, using the same slack DR as terminal)
  #
  # Build specs for obc_solve_binding.  CRITICAL: eq_idx is the declaration-
  # order row of the RELAX equation, which is match(eq_var_name, endo_names)
  # because extract_system_matrices reorders rows to declaration order.
  # Using cn$eq_relax (equation number in .mod file) is WRONG here.
  obc_specs <- list()
  for (j in seq_len(n_con)) {
    cn <- pr$constraints[[j]]

    # Row of the relax equation in declaration-order f_zero
    eq_var_nm  <- cn$eq_var_name   # LHS of relax equation, e.g. "lam"
    eq_idx_j   <- if (nzchar(eq_var_nm)) {
      match(eq_var_nm, endo_names)
    } else {
      # Fallback: single-variable case (bind var == relax LHS var)
      match(cn$var_name, endo_names)
    }

    if (is.na(eq_idx_j)) {
      warning(sprintf(
        ".occbin_build_pwlinear_drs: cannot locate relax equation for '%s'; aborting DR build.",
        cn$name))
      return(NULL)
    }

    var_idx_j <- match(cn$var_name, endo_names)
    if (is.na(var_idx_j)) {
      warning(sprintf(
        ".occbin_build_pwlinear_drs: bind var '%s' not in endo_names; aborting.",
        cn$var_name))
      return(NULL)
    }

    # Resolve the bound value using SS + params (handles PHI*steady_state(iv))
    rs <- model$occbin_constraints[[j]]
    bound_val <- if (!is.null(rs) && !is.na(rs$bound)) {
      as.numeric(rs$bound)
    } else if (!is.null(rs)) {
      tmp <- tryCatch(
        obc_resolve_block_specs(list(rs), model, ss = y_ss_num, params = params),
        error = function(e) list()
      )
      if (length(tmp) > 0L) tmp[[1L]]$bound else NA_real_
    } else {
      NA_real_
    }

    if (is.na(bound_val)) {
      warning(sprintf(
        ".occbin_build_pwlinear_drs: cannot resolve bound for constraint '%s'; aborting.",
        cn$name))
      return(NULL)
    }

    obc_specs[[j]] <- list(
      eq_idx  = eq_idx_j,
      var_idx = var_idx_j,
      bound   = bound_val
    )
  }

  obs_idx <- seq_len(n_endo)
  bind_res <- tryCatch(
    obc_solve_binding(sys, dr_slack, obc_specs, obs_idx),
    error = function(e) {
      warning(".occbin_build_pwlinear_drs: obc_solve_binding failed: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(bind_res)) return(NULL)

  # Build the binding system matrices once (for backward recursion reuse)
  built_b <- obc_build_binding_sys(sys, obc_specs)
  const_b  <- built_b$const_b  # zero before SS correction (done in solve loop)

  list(
    dr_slack  = dr_slack,
    dr_bind   = bind_res$dr,
    c_bind    = bind_res$c_full,
    state_idx = state_idx,
    sys_b     = built_b$sys_b,
    obc_specs = obc_specs,
    const_b   = const_b
  )
}


#' Piecewise-linear backward recursion for OccBin (Dynare-exact)
#'
#' Replicates Dynare's OccBin piecewise-linear algorithm exactly.
#'
#' For a fixed regime path, Dynare builds PERIOD-SPECIFIC decision rules
#' backward from t=T:
#'   - At each slack period: use ghx_slack (infinite-horizon DR)
#'   - At each binding period t: compute a binding DR with the NEXT period's
#'     forward DR as terminal condition (obc_solve_binding with ghx_fwd(t+1))
#'
#' For consecutive binding periods, each has a different ghx because the
#' terminal condition changes.  Slack periods all share ghx_slack.
#'
#' Once period-specific DRs are computed backward, propagate forward from t=1.
#'
#' For single-period binding episodes, this collapses to the standard
#' obc_solve_binding approach.  For multi-period episodes, this correctly
#' handles the backward telescoping.
#'
#' @param drs       List from .occbin_build_pwlinear_drs() with $dr_slack,
#'                  $sys_b, $obc_specs, $state_idx (sys_b and specs for the
#'                  backward recursion)
#' @param y0_num    Named numeric (n_endo): initial condition in levels
#' @param y_ss_num  Named numeric (n_endo): steady-state levels
#' @param eps_mat   T x n_exo matrix: structural shocks
#' @param regime_path Integer (T): bitmask (0=slack, nonzero=binding)
#' @param T         Integer: horizon
#' @param n_endo    Integer: number of endogenous variables
#' @param endo_names Character (n_endo)
#' @return List with $Y (T x n_endo levels matrix) and $converged (TRUE)
#' @noRd
.occbin_pwlinear_solve <- function(drs, y0_num, y_ss_num, eps_mat,
                                    regime_path, T, n_endo, endo_names) {
  dr_slack  <- drs$dr_slack
  state_idx <- drs$state_idx
  n_state   <- length(state_idx)
  n_exo     <- ncol(eps_mat)

  # Selection matrix E_s: E_s[i, state_idx[i]] = 1 (n_state x n_endo)
  E_s <- matrix(0, n_state, n_endo)
  for (k in seq_len(n_state)) E_s[k, state_idx[k]] <- 1

  # ==========================================================================
  # Backward pass: compute period-specific DRs
  # ==========================================================================
  # period_dr[[t]] = list(ghx, ghu, c) for the forward DR at period t.
  # For slack periods: period_dr = (ghx_slack, ghu_slack, 0).
  # For binding periods: computed via backward substitution.
  #
  # ghx_fwd[t+1] is the DR used as terminal condition when solving period t.
  # ghx_fwd[T+1] = ghx_slack (beyond horizon, return to slack).
  period_ghx <- vector("list", T)
  period_ghu <- vector("list", T)
  period_c   <- vector("list", T)

  # Forward DR at t+1 (start with slack for t+1 = T+1 = beyond horizon)
  ghx_fwd <- dr_slack$ghx
  ghu_fwd <- dr_slack$ghu

  has_bind_sys <- !is.null(drs$sys_b) && !is.null(drs$obc_specs)

  for (t in seq(T, 1L)) {
    is_binding <- (regime_path[t] != 0L)

    if (!is_binding) {
      # Slack period: use slack DR; next-period terminal is also slack DR
      period_ghx[[t]] <- dr_slack$ghx
      period_ghu[[t]] <- dr_slack$ghu
      period_c[[t]]   <- numeric(n_endo)
      # Forward DR remains slack for next (earlier) period
      ghx_fwd <- dr_slack$ghx
      ghu_fwd <- dr_slack$ghu
    } else if (!has_bind_sys) {
      # No system matrices available: fall back to precomputed single bind DR
      period_ghx[[t]] <- drs$dr_bind$ghx
      period_ghu[[t]] <- drs$dr_bind$ghu
      period_c[[t]]   <- drs$c_bind
      ghx_fwd <- drs$dr_bind$ghx
      ghu_fwd <- drs$dr_bind$ghu
    } else {
      # Binding period: solve bind DR with current ghx_fwd as terminal condition.
      # A_b(t) = f_zero_b + f_plus_b * ghx_fwd * E_s
      sys_b   <- drs$sys_b
      A_b     <- sys_b$f_zero + sys_b$f_plus %*% ghx_fwd %*% E_s
      const_b <- drs$const_b

      # Adjust const_b for SS deviations (same as obc_solve_binding)
      ys_vals <- y_ss_num
      for (s in drs$obc_specs) const_b[s$eq_idx] <- s$bound - ys_vals[s$var_idx]

      rcond_Ab <- rcond(A_b)
      if (!is.finite(rcond_Ab) || rcond_Ab <= .Machine$double.eps) {
        # Singular: fall back to single bind DR
        period_ghx[[t]] <- drs$dr_bind$ghx
        period_ghu[[t]] <- drs$dr_bind$ghu
        period_c[[t]]   <- drs$c_bind
        ghx_fwd <- drs$dr_bind$ghx
        ghu_fwd <- drs$dr_bind$ghu
      } else {
        rhs_ghx_t <- -sys_b$f_minus[, state_idx, drop = FALSE]
        rhs_ghu_t <- -sys_b$f_exo
        rhs_c_t   <- const_b

        ghx_t <- tryCatch(solve(A_b, rhs_ghx_t), error = function(e) NULL)
        ghu_t <- tryCatch(solve(A_b, rhs_ghu_t), error = function(e) NULL)
        c_t   <- tryCatch(solve(A_b, rhs_c_t),   error = function(e) NULL)

        if (is.null(ghx_t) || is.null(ghu_t) || is.null(c_t)) {
          period_ghx[[t]] <- drs$dr_bind$ghx
          period_ghu[[t]] <- drs$dr_bind$ghu
          period_c[[t]]   <- drs$c_bind
          ghx_fwd <- drs$dr_bind$ghx
        } else {
          period_ghx[[t]] <- ghx_t
          period_ghu[[t]] <- ghu_t
          period_c[[t]]   <- as.numeric(c_t)
          ghx_fwd <- ghx_t
          ghu_fwd <- ghu_t
        }
      }
    }
  }

  # ==========================================================================
  # Forward pass: propagate using period-specific DRs
  # ==========================================================================
  Y      <- matrix(0, nrow = T, ncol = n_endo)
  colnames(Y) <- endo_names
  s_prev <- (y0_num - y_ss_num)[state_idx]

  for (t in seq_len(T)) {
    eps_t  <- eps_mat[t, ]
    delta  <- as.numeric(period_ghx[[t]] %*% s_prev +
                         period_ghu[[t]] %*% eps_t +
                         period_c[[t]])
    Y[t, ] <- y_ss_num + delta
    s_prev <- delta[state_idx]
  }

  list(Y = Y, converged = TRUE)
}


# =============================================================================
# Main solver: stacked Newton + regime-switching outer loop
# =============================================================================

#' Stacked Newton path solver for nonlinear OccBin models
#'
#' Solves the T-period deterministic path for a nonlinear DSGE model with
#' occasionally binding constraints declared via bind/relax equation pairs.
#'
#' The solver has two nested loops:
#' \enumerate{
#'   \item \strong{Inner loop (Newton or linear propagation depending on method)}:
#'     For a fixed regime path, solve the T-period system.
#'   \item \strong{Outer loop (regime switching)}: After propagation, check
#'     complementarity for each constraint and period.  Flip periods where the
#'     constraint is violated, update the regime path, and re-solve.
#' }
#'
#' When \code{method = "pwlinear"}, the inner loop is replaced by a simple
#' O(T) linear forward pass using the first-order decision rules for the slack
#' and binding regimes (Dynare OccBin piecewise-linear method).  This gives
#' exact parity with Dynare's OccBin output.  When
#' \code{method = "nonlinear"} (default), the original Newton solver is used.
#'
#' Unlike pf_newton_solve(), this solver handles models where n_eq > n_endo
#' (i.e., models with bind/relax equation pairs).  The regime selects which
#' equation variant is active per endogenous variable per period.
#'
#' @param compiled       dynhr_compiled (from compile_model)
#' @param y0             Named numeric vector: endogenous state at t=0
#' @param y_ss           Named numeric vector: steady state (terminal condition)
#' @param shock_path     T × n_exo numeric matrix of structural shocks.
#'                       Column names must match varexo_names.
#' @param params         Named numeric parameter vector
#' @param constraints    Constraint list from occbin_parse_bind_relax(); if NULL
#'                       or zero-length, solves the unrestricted path (all slack)
#' @param regime_path_init Integer vector (length T) for warm-starting the
#'                       regime path; NULL → all-slack initial guess
#' @param method         Character: "nonlinear" (default) or "pwlinear".
#'                       "pwlinear" uses Dynare's piecewise-linear approach
#'                       (first-order DRs per regime, O(T) forward pass) for
#'                       exact Dynare parity.  "nonlinear" uses the stacked
#'                       Newton solver for higher accuracy.
#' @param max_iter       Maximum Newton iterations per regime (default 50; ignored for pwlinear)
#' @param tol            Newton convergence tolerance on max|R| (default 1e-8; ignored for pwlinear)
#' @param max_regime_iter Maximum outer regime-switching iterations (default 30)
#' @param step_size      Newton step-length (default 1.0; ignored for pwlinear)
#' @param line_search    Logical: perform backtracking line search (default TRUE; ignored for pwlinear)
#' @return List with:
#'   $Y           T × n_endo solution matrix
#'   $regime_path Integer vector (length T): final per-period regime bitmask
#'   $irf         T × n_endo deviation from steady state
#'   $converged   Logical: TRUE if converged
#'   $outer_iter  Integer: outer regime-switching iterations used
#'   $n_iter      Integer: Newton iterations on final regime solve (0 for pwlinear)
#'   $max_res     Numeric: final max|R| (0 for pwlinear)
#'   $endo_names  Character: variable ordering of Y columns
#'
#' @references
#'   Guerrieri, L., & Iacoviello, M. (2015). OccBin: A toolkit for solving
#'     dynamic models with occasionally binding constraints easily.
#'     \emph{Journal of Monetary Economics}, 70, 22-38.
#' @export
occbin_solve_path <- function(compiled,
                               y0,
                               y_ss,
                               shock_path,
                               params,
                               constraints       = NULL,
                               regime_path_init  = NULL,
                               method            = c("nonlinear", "pwlinear"),
                               max_iter          = 50L,
                               tol               = 1e-8,
                               max_regime_iter   = 30L,
                               step_size         = 1.0,
                               line_search       = TRUE) {
  method <- match.arg(method)

  dyn    <- compiled$dynamic
  n_endo <- length(dyn$endo_names)
  n_eq   <- dyn$n_eq

  # Normalise shock_path
  if (!is.matrix(shock_path)) shock_path <- matrix(shock_path, nrow = 1L)
  T <- nrow(shock_path)
  eps_mat <- matrix(0, nrow = T, ncol = length(dyn$exo_names))
  colnames(eps_mat) <- dyn$exo_names
  if (!is.null(colnames(shock_path))) {
    for (nm in intersect(colnames(shock_path), dyn$exo_names))
      eps_mat[, nm] <- shock_path[, nm]
  } else {
    nc <- min(ncol(shock_path), length(dyn$exo_names))
    eps_mat[, seq_len(nc)] <- shock_path[, seq_len(nc)]
  }

  # Align y0 and y_ss.
  # Keep names on y_ss_num so steady_state(X) → ss["X"] look-ups in compiled
  # residual/Jacobian functions resolve correctly for equations like
  # iv = PHI*steady_state(iv) or ivhat = 100*(iv/steady_state(iv)-1).
  y0_num   <- as.numeric(y0[dyn$endo_names])
  y_ss_num <- y_ss[dyn$endo_names]
  if (!is.numeric(y_ss_num)) y_ss_num <- as.numeric(y_ss_num)
  names(y_ss_num) <- dyn$endo_names

  # Precompute column metadata
  cmap_meta <- .occbin_col_meta(dyn)
  # Also need .pf_col_meta for .pf_make_dy
  meta <- .pf_col_meta(dyn)

  # NEW-OCC1: default to the COMPILED OccBin constraints (which carry the
  # eq_bind/eq_relax regime-equation indices the solver needs) when the caller
  # passes none but the model is an OccBin model. The raw model$occbin_constraints
  # returned by parse_mod() do NOT carry these indices.
  if ((is.null(constraints) || length(constraints) == 0L) &&
      !is.null(compiled$occbin) &&
      length(compiled$occbin$parse_result$constraints) > 0L) {
    constraints <- compiled$occbin$parse_result$constraints
  }
  # NEW-OCC1: fail loud if supplied constraints lack eq_bind/eq_relax (e.g. the
  # raw model$occbin_constraints) — NEVER silently solve the unconstrained path
  # (the old pwlinear footgun: conv=TRUE, no warning, constraint ignored).
  if (length(constraints) > 0L) {
    bad <- vapply(constraints, function(cn)
      is.null(cn$eq_bind) || is.null(cn$eq_relax) ||
        anyNA(c(cn$eq_bind, cn$eq_relax)), logical(1))
    if (any(bad)) {
      stop("occbin_solve_path: supplied `constraints` lack eq_bind/eq_relax ",
           "regime-equation indices, so the regime machinery cannot select the ",
           "bind/relax equations. Pass the COMPILED constraints ",
           "`compile_model(model)$occbin$parse_result$constraints` (or leave ",
           "`constraints = NULL` to use them automatically) -- NOT the raw ",
           "`model$occbin_constraints` from parse_mod(), which carry no indices ",
           "and would otherwise silently solve the UNCONSTRAINED path.",
           call. = FALSE)
    }
  }

  # Parse constraints if provided
  if (is.null(constraints) || length(constraints) == 0L) {
    # No constraints — use all equations (must have n_eq == n_endo)
    if (n_eq != n_endo) {
      stop(sprintf(
        "occbin_solve_path: n_eq (%d) != n_endo (%d) and no usable constraints. ",
        n_eq, n_endo),
        "This is an OccBin model: pass `constraints = NULL` with a model compiled ",
        "by compile_model() (so compiled$occbin is populated), or supply ",
        "compiled$occbin$parse_result$constraints. Or use pf_newton_solve() for a ",
        "single-regime path.")
    }
    # Use all equations for all periods (single regime)
    regime_fns <- list(list(
      eq_indices = seq_len(n_eq),
      res_fn     = dyn$residuals_fn,
      jac_fn     = dyn$jacobian_fn
    ))
    n_regimes <- 1L
    constraints <- list()
  } else {
    # Build regime map and regime functions
    parse_result <- list(
      constraints   = constraints,
      neutral_eqs   = .find_neutral_eqs(compiled, constraints),
      n_constraints = length(constraints)
    )
    regime_map    <- occbin_build_regime_map(parse_result)
    regime_fns    <- occbin_build_regime_fns(compiled, regime_map)
    n_regimes     <- length(regime_fns)
  }

  n_constraints <- length(constraints)

  # Initialise path at steady state
  Y <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
  colnames(Y) <- dyn$endo_names

  # Initialise regime path
  if (!is.null(regime_path_init) && length(regime_path_init) == T) {
    regime_path <- as.integer(regime_path_init)
    # Validate: no regime index out of range
    regime_path <- pmin(pmax(regime_path, 0L), as.integer(2L^n_constraints - 1L))
  } else {
    regime_path <- integer(T)  # all slack (regime 0)
  }

  converged      <- FALSE
  n_outer_iter   <- 0L
  n_newton_iter  <- 0L
  max_res_final  <- Inf

  # ==========================================================================
  # Piecewise-linear path: build DRs before outer loop
  # ==========================================================================
  pwl_drs <- NULL
  if (method == "pwlinear" && n_constraints > 0L) {
    pwl_drs <- tryCatch(
      .occbin_build_pwlinear_drs(compiled, y_ss_num, params),
      error = function(e) {
        warning("occbin_solve_path: pwlinear DR build failed: ",
                conditionMessage(e), "; falling back to nonlinear.")
        NULL
      }
    )
    if (is.null(pwl_drs)) {
      method <- "nonlinear"
      warning("occbin_solve_path: pwlinear setup failed, using nonlinear.")
    }
  }

  # ==========================================================================
  # Outer regime-switching loop
  # ==========================================================================
  for (outer_iter in seq_len(max_regime_iter)) {
    n_outer_iter <- outer_iter

    if (method == "pwlinear" && !is.null(pwl_drs)) {
      # ------------------------------------------------------------------
      # Piecewise-linear DR-recursive solve (Dynare OccBin method)
      # Uses precomputed slack + binding DRs; propagates period by period.
      # ------------------------------------------------------------------
      pwl_res <- .occbin_pwlinear_solve(
        pwl_drs, y0_num, y_ss_num, eps_mat, regime_path,
        T, n_endo, dyn$endo_names
      )
      Y                <- pwl_res$Y
      newton_converged <- pwl_res$converged
      n_newton_iter    <- 0L
      max_res_final    <- 0
    } else {
      # ------------------------------------------------------------------
      # Inner Newton loop for the current regime path
      # ------------------------------------------------------------------
      newton_converged <- FALSE

      for (newton_iter in seq_len(max_iter)) {
        n_newton_iter <- newton_iter

        # Build stacked system
        sys <- .occbin_build_stacked_system(
          Y, y0_num, y_ss_num, eps_mat, meta, cmap_meta,
          regime_fns, regime_path, params, T, n_endo
        )

        max_res <- max(abs(sys$R))
        max_res_final <- max_res

        if (max_res < tol) {
          newton_converged <- TRUE
          break
        }

        # Solve J * delta = -R via sparse QR
        delta_vec <- {
          # Use Matrix's sparse solve
          as.numeric(Matrix::solve(sys$J, -sys$R))
        }

        if (anyNA(delta_vec)) {
          warning("occbin_solve_path: singular Jacobian at Newton iter ", newton_iter,
                  "; aborting Newton for current regime.")
          break
        }

        # Line search (backtracking)
        alpha <- step_size
        if (line_search) {
          # Evaluate residual norm at current point
          norm_R0 <- max_res
          for (ls_iter in seq_len(10L)) {
            # Trial point
            Y_trial <- Y
            for (t in seq_len(T)) {
              idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
              Y_trial[t, ] <- Y[t, ] + alpha * delta_vec[idx_t]
            }
            # Evaluate residual at trial point
            trial_res <- .occbin_eval_residual(
              Y_trial, y0_num, y_ss_num, eps_mat, meta,
              regime_fns, regime_path, params, T, n_endo
            )
            norm_R1 <- max(abs(trial_res))
            # Armijo condition: norm_R1 <= (1 - alpha/2) * norm_R0
            if (norm_R1 <= (1 - alpha * 0.5) * norm_R0 || norm_R1 < tol) break
            alpha <- alpha * 0.5
          }
        }

        # Update path
        for (t in seq_len(T)) {
          idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
          Y[t, ] <- Y[t, ] + alpha * delta_vec[idx_t]
        }
      }
    }

    if (n_constraints == 0L) {
      # No constraints: done after first solve
      converged <- newton_converged
      break
    }

    # ------------------------------------------------------------------
    # Complementarity check (regime flipping)
    # ------------------------------------------------------------------
    # For each constraint and period, check if the current regime is consistent
    # with the solved path.  Use the parsed bind condition from the
    # occbin_constraints block (via .occbin_eval_bind_condition) instead of
    # the previous hardcoded heuristic.
    #
    # Y is in LEVELS; the evaluator handles steady_state(X) resolution and
    # parameter substitution in the bound expression.

    # Evaluate which constraint/period pairs should be binding.
    #
    # The complementarity check must use the NOTIONAL value of each constrained
    # variable (what it would be under the relax/slack equation), NOT the realized
    # Y. For a SLACK period the realized var already IS the notional (relax eq
    # active), so Y is fine. For a BINDING period the realized var is pinned to the
    # bound (e.g. r = r_lb), so `Y[t,var] < bound` is ALWAYS false → the period
    # would immediately relax, re-violate, and the outer loop oscillates forever
    # (the Phase-2 bug: nk_zlb hit max_regime_iter with 0 binding). The notional is
    # recovered from the slack-regime residual of the relax equation:
    #   relax_eq residual = var - notional  (equation_to_residual is lhs - rhs,
    #   relax eq is `var = <expr>`), so  notional = Y[t,var] - relax_resid.
    # Evaluating the all-slack residual at the current path gives relax_resid for
    # every period; for slack periods relax_resid ~ 0 (eq satisfied) so the notional
    # equals Y (no-op). We then run the SAME bind-condition evaluator on this
    # notional path, so the bound/direction logic is shared and sign-robust.
    #
    # For pwlinear: compute the notional using the slack DR propagation path.
    # The notional at binding periods = y_ss + ghx_slack * s_{t-1}^{actual} + ghu * eps_t
    # where s_{t-1}^{actual} is the state from the actual (regime-aware) path.
    Y_notional <- Y
    any_binding <- any(regime_path != 0L)
    if (any_binding) {
      # Compute the notional (slack-regime) path to check complementarity.
      # The notional is the all-slack-regime linear solve at the current Y.
      # For slack periods: Y_notional = Y (relax eq is satisfied).
      # For binding periods: Y_notional[t, bind_var] = what the bind var
      # would be if not constrained (the slack-regime prediction).
      # We get this from the all-slack residual of the RELAX equation:
      #   relax_eq residual = var - notional, so notional = var - residual.
      if (method == "pwlinear" && !is.null(pwl_drs)) {
        # For pwlinear: compute the notional path using the all-slack DR
        # (propagate with dr_slack for every period).
        slack_path_res <- .occbin_pwlinear_solve(
          pwl_drs, y0_num, y_ss_num, eps_mat, integer(T),
          T, n_endo, dyn$endo_names
        )
        if (slack_path_res$converged) {
          Y_notional_slack <- slack_path_res$Y
          # Override notional only for binding periods and bind vars
          for (j in seq_len(n_constraints)) {
            bind_bit <- 2L^(j - 1L)
            vj       <- match(constraints[[j]]$var_name, dyn$endo_names)
            if (is.na(vj)) next
            for (t in seq_len(T)) {
              if (bitwAnd(regime_path[t], bind_bit) != 0L) {
                Y_notional[t, vj] <- Y_notional_slack[t, vj]
              }
            }
          }
        }
      } else {
        R_slack <- .occbin_eval_residual(
          Y, y0_num, y_ss_num, eps_mat, meta,
          regime_fns, integer(T), params, T, n_endo
        )
        slack_eq_indices <- regime_fns[[1L]]$eq_indices   # regime 0 = all slack
        for (j in seq_len(n_constraints)) {
          bind_bit <- 2L^(j - 1L)
          vj   <- match(constraints[[j]]$var_name, dyn$endo_names)
          posj <- match(constraints[[j]]$eq_relax, slack_eq_indices)
          if (is.na(vj) || is.na(posj)) next
          for (t in seq_len(T)) {
            if (bitwAnd(regime_path[t], bind_bit) != 0L) {
              relax_resid <- R_slack[(t - 1L) * n_endo + posj]
              Y_notional[t, vj] <- Y[t, vj] - relax_resid
            }
          }
        }
      }
    }

    bind_mat <- .occbin_eval_bind_condition(Y_notional, dyn, y_ss_num, constraints, params)

    new_regime <- regime_path
    any_flip   <- FALSE

    for (j in seq_len(n_constraints)) {
      bind_bit <- 2L^(j - 1L)

      for (t in seq_len(T)) {
        is_binding  <- bitwAnd(regime_path[t], bind_bit) != 0L
        should_bind <- bind_mat[j, t]

        if (should_bind && !is_binding) {
          # Should bind but currently slack → flip to binding
          new_regime[t] <- bitwOr(new_regime[t], bind_bit)
          any_flip <- TRUE
        } else if (!should_bind && is_binding) {
          # Currently binding but should be slack → flip to slack
          new_regime[t] <- bitwAnd(new_regime[t], bitwNot(bind_bit))
          any_flip <- TRUE
        }
      }
    }

    if (!any_flip) {
      # Regime path converged
      converged <- newton_converged
      break
    }

    regime_path <- new_regime
  }

  # Build output
  y_ss_mat <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)

  list(
    Y           = Y,
    regime_path = regime_path,
    irf         = Y - y_ss_mat,
    converged   = converged,
    outer_iter  = n_outer_iter,
    n_iter      = n_newton_iter,
    max_res     = max_res_final,
    endo_names  = dyn$endo_names
  )
}


# =============================================================================
# Helper: evaluate stacked residual (no Jacobian)
# =============================================================================

#' Evaluate the stacked residual for a given path and regime
#'
#' @inheritParams .occbin_build_stacked_system
#' @return Numeric vector (T*n_endo)
#' @noRd
.occbin_eval_residual <- function(Y, y0_num, y_ss_num, eps_mat,
                                   meta, regime_fns, regime_path,
                                   params, T, n_endo) {
  R <- numeric(T * n_endo)

  for (t in seq_len(T)) {
    row_off <- (t - 1L) * n_endo
    dy <- .pf_make_dy(meta, Y, y0_num, y_ss_num, eps_mat[t, ], t, T)
    r_idx <- regime_path[t] + 1L
    rf <- regime_fns[[r_idx]]
    R[row_off + seq_len(n_endo)] <- rf$res_fn(dy, params, y_ss_num)
  }

  R
}


# =============================================================================
# Helper: find neutral equations
# =============================================================================

#' Find neutral equations (not referenced by any bind/relax constraint)
#'
#' @param compiled    dynhr_compiled
#' @param constraints List of constraint specs from occbin_parse_bind_relax()
#' @return Integer vector of equation indices that are not bind/relax variants
#' @noRd
.find_neutral_eqs <- function(compiled, constraints) {
  all_eqs <- seq_len(compiled$dynamic$n_eq)
  annotated <- c()
  for (cn in constraints) {
    annotated <- c(annotated, cn$eq_bind, cn$eq_relax)
  }
  setdiff(all_eqs, unique(annotated))
}


# =============================================================================
# OBC IRF computation using nonlinear OccBin solver
# =============================================================================

#' Compute IRFs using the OccBin path solver (nonlinear or piecewise-linear)
#'
#' Wraps occbin_solve_path() into the standard IRFCollection interface, similar
#' to compute_irfs_obc() but using the regime-aware solver.  Supports both the
#' nonlinear stacked-Newton solver (default, more accurate) and the
#' piecewise-linear Dynare-style solver (method="pwlinear", exact Dynare parity).
#'
#' @param compiled      dynhr_compiled
#' @param y_ss          Named numeric vector: steady state
#' @param model         dynhr_mod (for shock standard deviations)
#' @param params        Named numeric parameter vector; NULL -> model$param_values
#' @param constraints   Constraint list from occbin_parse_bind_relax(); NULL ->
#'   auto-resolved from compiled$occbin$parse_result$constraints when available
#' @param n_periods     Number of IRF periods (default 40L)
#' @param shock_size    Shock size multiplier (default 1)
#' @param method        Character: "nonlinear" (default) or "pwlinear".
#'   "pwlinear" uses Dynare's piecewise-linear approach for exact Dynare parity.
#'   "nonlinear" uses the stacked Newton solver for higher accuracy.
#' @param ...           Additional arguments passed to occbin_solve_path()
#' @return IRFCollection object (same structure as compute_irfs()); each shock's
#'   entry is the DEVIATION from steady state (T x n_endo).  The object also
#'   carries \code{attr(., "regime_paths")} (a named list, one integer vector
#'   per shock) and \code{attr(., "converged")} (a named logical vector).
#' @export
occbin_compute_irfs <- function(compiled, y_ss, model, params = NULL,
                                 constraints = NULL,
                                 n_periods   = 40L,
                                 shock_size  = 1,
                                 method      = c("nonlinear", "pwlinear"),
                                 ...) {
  method <- match.arg(method)
  if (is.null(params)) params <- model$param_values

  # Auto-fill constraints from compiled$occbin when not explicitly provided.
  # This is the normal path for a model parsed with an occbin_constraints block.
  if (is.null(constraints) && !is.null(compiled$occbin)) {
    constraints <- compiled$occbin$parse_result$constraints
  }

  endo  <- compiled$dynamic$endo_names
  exo   <- compiled$dynamic$exo_names
  n_exo <- length(exo)

  shock_stderr <- .get_shock_stderr(model, exo, params)

  irfs         <- vector("list", n_exo)
  regime_paths <- vector("list", n_exo)
  converged_v  <- logical(n_exo)
  names(irfs)         <- exo
  names(regime_paths) <- exo
  names(converged_v)  <- exo

  y0_num <- as.numeric(y_ss[endo])
  names(y0_num) <- endo

  for (k in seq_along(exo)) {
    shock_seq <- matrix(0, nrow = n_exo, ncol = n_periods)
    shock_seq[k, 1L] <- shock_stderr[exo[k]] * shock_size

    res <- occbin_solve_path(
      compiled     = compiled,
      y0           = y0_num,
      y_ss         = y0_num,
      shock_path   = t(shock_seq),
      params       = params,
      constraints  = constraints,
      method       = method,
      ...
    )

    # Use res$irf (deviation from SS) not res$Y (absolute levels).
    # res$irf = Y - y_ss_mat is computed inside occbin_solve_path.
    irf_mat <- res$irf
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))
    irfs[[k]]         <- irf_mat
    regime_paths[[k]] <- res$regime_path
    converged_v[[k]]  <- isTRUE(res$converged)
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")    <- n_periods
  attr(irfs, "endo_names")   <- endo
  attr(irfs, "exo_names")    <- exo
  attr(irfs, "regime_paths") <- regime_paths
  attr(irfs, "converged")    <- converged_v
  irfs
}


# =============================================================================
# Extract shock standard errors (internal helper)
# =============================================================================
