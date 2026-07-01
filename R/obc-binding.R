## R/obc-binding.R
## --------------------------------------------------------------------------
## OBC binding-regime algebra.
##
## Provides:
##   obc_build_binding_sys()   -- substitute constraint equations into Jacobians
##   obc_binding_constant()    -- compute bound-induced constant offsets
##   obc_solve_binding()       -- OccBin terminal-substitution policy solve
##
## The OccBin approach (Guerrieri-Iacoviello 2015) constructs the binding-regime
## policy by replacing each tagged equation with its constraint equality and
## substituting the SLACK policy as the next-period terminal condition.  This
## gives a determinate linear system even when the Taylor principle is lost
## (e.g. at the ZLB).
##
## Design note: a near-zero lag coefficient (1e-10) is retained on the
## constrained variable so its classification as a state variable is preserved
## and the Kalman filter state dimension is unchanged across regimes.
## --------------------------------------------------------------------------


# =============================================================================
# Binding system construction
# =============================================================================

#' Build the binding-regime system matrices
#'
#' For each OBC spec, replaces row eq_idx of the system Jacobians with the
#' constraint equation: var_idx_t = bound. To preserve the state-space
#' dimension (critical for the Kalman filter), a near-zero lag coefficient
#' (1e-10) is retained for the constrained variable so the classifier keeps
#' it as a state variable. The bound itself is captured in const_b.
#'
#' @param sys   System matrices list (from extract_system_matrices_fast)
#' @param specs List of OBC specs from obc_parse_tags
#' @return List:
#'   $sys_b   -- modified system matrices (same structure as sys)
#'   $const_b -- numeric vector (length n_endo), nonzero at bound equations
#' @noRd
obc_build_binding_sys <- function(sys, specs) {
  f_minus_b <- sys$f_minus
  f_zero_b  <- sys$f_zero
  f_plus_b  <- sys$f_plus
  f_exo_b   <- sys$f_exo
  const_b   <- numeric(sys$n_endo)

  for (s in specs) {
    i <- s$eq_idx
    j <- s$var_idx
    b <- s$bound

    # Replace row i with the binding constraint: var_j_t = bound
    # f_zero[i, j] = 1 (contemporaneous var_j)
    # const_b[i]   = b (the bound, which becomes the RHS forcing)
    #
    # Near-zero lag on j (1e-10 << any economic quantity) retains var_j in
    # the state-variable classification, keeping T_b and T_s the same size.
    f_minus_b[i, ] <- 0
    f_zero_b[i, ]  <- 0
    f_plus_b[i, ]  <- 0
    f_exo_b[i, ]   <- 0

    f_zero_b[i, j]  <- 1
    f_minus_b[i, j] <- 1e-10   # keeps var_j classified as state variable
    const_b[i]      <- b
  }

  sys_b         <- sys
  sys_b$f_minus <- f_minus_b
  sys_b$f_zero  <- f_zero_b
  sys_b$f_plus  <- f_plus_b
  sys_b$f_exo   <- f_exo_b

  list(sys_b = sys_b, const_b = const_b)
}


#' Compute the binding-regime constant offsets
#'
#' The bound value (e.g. -0.00401606 for ZLB) is non-zero in deviation-from-SS
#' form, so pinning var_t = bound creates a non-homogeneous system. The
#' "binding-regime steady state" (where all variables settle if the OBC keeps
#' binding, with shocks at zero) satisfies:
#'
#'   (f_minus_b + f_zero_b + f_plus_b) * c_full = const_b
#'
#' Splits c_full into a state-level constant (c_state, for the transition
#' equation) and an observable-level constant (c_obs, for the measurement
#' equation).
#'
#' @param sys_b     Binding-regime system matrices (from obc_build_binding_sys)
#' @param const_b   Constant vector from obc_build_binding_sys
#' @param state_idx Integer vector: state variable indices (from dr$state_idx)
#' @param obs_idx   Integer vector: observable variable indices in endo vector
#' @return List with $c_full (n_endo), $c_state (n_state), $c_obs (n_obs)
#' @noRd
obc_binding_constant <- function(sys_b, const_b, state_idx, obs_idx) {
  A <- sys_b$f_minus + sys_b$f_zero + sys_b$f_plus
  c_full <- solve(A, const_b)
  list(
    c_full  = c_full,
    c_state = c_full[state_idx],
    c_obs   = c_full[obs_idx]
  )
}


# =============================================================================
# Binding policy function
# =============================================================================

#' Compute the binding-regime policy function via the OccBin terminal substitution
#'
#' A permanently pegged rate (or any binding OBC) makes the linearised model
#' indeterminate: the Taylor principle is lost and there is no nominal anchor.
#' The OccBin solution (Guerrieri-Iacoviello 2015) avoids this by substituting
#' the SLACK policy function as the terminal condition for the next period:
#'
#'   E\[y_{t+1}\] ~= dr_slack$ghx * s_t    (next period is expected to be slack)
#'
#' Substituting into the binding-regime system:
#'   (f_zero_b + f_plus_b * ghx_slack * E_s) * y_t = const_b - f_minus_b\[:,state_idx\] * s_{t-1}
#'
#' where E_s is the n_state x n_endo selection matrix that picks state variables.
#' This linear system is generically non-singular (slack policy provides the anchor)
#' and gives a well-defined, determinate binding-regime policy function.
#'
#' Returns a pseudo-DecisionRules object with the same interface as dr_slack
#' (ghx, ghu, state_idx, etc.) so kalman_filter_obc can use it directly.
#'
#' @param sys       System matrices (from extract_system_matrices_fast)
#' @param dr_slack  Slack-regime DecisionRules (provides terminal condition)
#' @param specs     OBC spec list from obc_parse_tags
#' @param obs_idx   Integer vector: observable positions in endo vector
#' @return List($dr, $c_state, $c_obs, $c_full), or NULL if linear solve fails
#' @noRd
obc_solve_binding <- function(sys, dr_slack, specs, obs_idx) {
  built   <- obc_build_binding_sys(sys, specs)
  sys_b   <- built$sys_b
  const_b <- built$const_b

  # Correct const_b to deviation form.
  # The MCP bound annotation (e.g. [mcp = 'i > 1']) expresses the constraint
  # bound in the SAME UNITS as the model variable.  For model(linear) the
  # steady state is 0 for all vars, so const_b = bound already.  For model
  # (nonlinear), the perturbation Jacobians are in DEVIATION-FROM-SS form, so
  # the binding constant must be (bound - ss_level):
  #   v_t = bound  =>  (v_t - v_ss) = bound - v_ss
  # dr_slack$ys holds the SS vector; ys[j] = 0 for any linear model, so this
  # correction is safe and correct for both model types.
  # dr$ys is the full dynhr_steady object (from solve_steady); extract the
  # numeric steady-state vector.  For model(linear) this is always zero.
  ys_vals <- if (inherits(dr_slack$ys, "dynhr_steady")) dr_slack$ys$values else dr_slack$ys
  for (s in specs) {
    const_b[s$eq_idx] <- s$bound - ys_vals[s$var_idx]
  }

  n_endo    <- sys$n_endo
  n_exo     <- sys$n_exo
  state_idx <- dr_slack$state_idx
  n_state   <- length(state_idx)

  # Selection matrix E_s (n_state x n_endo): picks state variables from y_t
  # so that y_t[state_idx] = E_s %*% y_t
  E_s <- matrix(0, n_state, n_endo)
  for (i in seq_len(n_state)) E_s[i, state_idx[i]] <- 1

  # Binding-regime coefficient matrix with slack terminal condition substituted:
  #   A_b = f_zero_b + f_plus_b * ghx_slack * E_s   (n_endo x n_endo)
  A_b <- sys_b$f_zero + sys_b$f_plus %*% dr_slack$ghx %*% E_s

  # Solve for the three components of the binding policy function:
  #   y_t = ghx_b * s_{t-1} + ghu_b * eps_t + c_b
  # where:
  #   ghx_b (n_endo x n_state) = -A_b^{-1} * f_minus_b[:, state_idx]
  #   ghu_b (n_endo x n_exo)   = -A_b^{-1} * f_exo_b
  #   c_b   (n_endo)           =  A_b^{-1} * const_b
  rhs_ghx <- -sys_b$f_minus[, state_idx, drop = FALSE]
  rhs_ghu <- -sys_b$f_exo
  rhs_c   <- const_b

  # Check condition number before solve; fall back to SVD pseudoinverse if singular
  rcond_Ab <- rcond(A_b)
  if (is.finite(rcond_Ab) && rcond_Ab > .Machine$double.eps) {
    soln <- list(
      ghx = solve(A_b, rhs_ghx),
      ghu = solve(A_b, rhs_ghu),
      c   = solve(A_b, rhs_c)
    )
  } else {
    # Fallback: pseudoinverse via SVD
    sv    <- svd(A_b)
    d_inv <- ifelse(sv$d > 1e-10, 1 / sv$d, 0)
    A_inv <- sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
    soln <- list(ghx = A_inv %*% rhs_ghx, ghu = A_inv %*% rhs_ghu, c = A_inv %*% rhs_c)
  }

  rownames(soln$ghx) <- dr_slack$endo_names
  colnames(soln$ghx) <- dr_slack$endo_names[state_idx]
  rownames(soln$ghu) <- dr_slack$endo_names
  if (length(dr_slack$exo_names) > 0) colnames(soln$ghu) <- dr_slack$exo_names

  # Construct a pseudo-DecisionRules object for the binding regime.
  # Uses the same state_idx as dr_slack so kalman_filter_obc dimensions match.
  dr_b <- structure(list(
    ghx          = soln$ghx,
    ghu          = soln$ghu,
    ys           = dr_slack$ys,
    endo_names   = dr_slack$endo_names,
    exo_names    = dr_slack$exo_names,
    state_vars   = dr_slack$state_vars,
    state_idx    = state_idx,
    n_state      = n_state,
    n_stable     = n_state,
    n_exo        = n_exo,
    eigenvalues  = NULL,
    n_unstable   = 0L,
    bk_satisfied = TRUE
  ), class = "DecisionRules")

  c_full  <- soln$c
  c_state <- c_full[state_idx]
  c_obs   <- c_full[obs_idx]

  list(dr = dr_b, c_state = c_state, c_obs = c_obs, c_full = c_full)
}
