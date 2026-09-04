## R/mcp-utils.R
## --------------------------------------------------------------------------
## MCP solver utility functions: Fischer-Burmeister complementarity function,
## merit function, column metadata, active-set extraction, and other helpers.
##
## Provides:
##   mcp_fb()              -- Fischer-Burmeister function φ(a, b)
##   mcp_fb_deriv()        -- FB partial derivatives ∂φ/∂a, ∂φ/∂b
##   mcp_compute_merit()   -- squared FB residual merit function θ = ½||R||²
##   mcp_assemble_active_set() -- extract binding pattern from solved path
##   mcp_col_meta()        -- precompute column metadata for sparse assembly
##   mcp_natural_residual() -- min(a, b) complementarity residual
##
## REFERENCES
##   Fischer, A. (1992). "A special Newton-type optimization method."
##     Optimization, 24(3-4), 269-284.
##   De Luca, T., Facchinei, F., & Kanzow, C. (1996). "A semismooth equation
##     approach to the solution of nonlinear complementarity problems."
##     Mathematical Programming, 75(3), 407-439.
## --------------------------------------------------------------------------


# =============================================================================
# Fischer-Burmeister function
# =============================================================================

#' Fischer-Burmeister complementarity function
#'
#' Computes φ(a, b) = a + b − √(a² + b²).  This function has the property:
#'   φ(a, b) = 0  ⇔  a ≥ 0, b ≥ 0, a·b = 0
#' i.e., a and b form a complementary pair.
#'
#' For vector inputs, returns element-wise φ values.
#'
#' @param a,b Numeric vectors (same length, or one scalar).  Typically:
#'   For lower bound: a = x − bound, b = F(x) (equation residual)
#'   For upper bound: a = bound − x, b = −F(x)
#' @return Numeric vector of φ values.  Zero means complementarity holds.
#' @noRd
#' @examples
#' mcp_fb(0, 0)    # = 0 (complementarity satisfied at the kink)
#' mcp_fb(1, 0)    # = 0 (complementarity: a > 0, b = 0)
#' mcp_fb(0, 1)    # = 0 (complementarity: a = 0, b > 0)
#' mcp_fb(1, 1)    # > 0 (violation: both a and b positive)
#' mcp_fb(-1, 0)   # > 0 (violation: a negative)
mcp_fb <- function(a, b) {
  a + b - sqrt(a * a + b * b)
}


#' FB function partial derivatives
#'
#' Computes ∂φ/∂a and ∂φ/∂b for φ(a, b) = a + b − √(a² + b²).
#' At (a, b) = (0, 0), the function is not differentiable; we return
#' the subgradient (0, 0) which corresponds to "no correction needed"
#' (the Newton direction for this row is undefined, so we zero it out).
#'
#' @param a,b Numeric vectors (same length)
#' @return List with $da and $db: partial derivatives (same length as a, b)
#' @noRd
mcp_fb_deriv <- function(a, b) {
  nrm <- sqrt(a * a + b * b)
  zero <- nrm == 0
  da <- ifelse(zero, 0, 1 - a / nrm)
  db <- ifelse(zero, 0, 1 - b / nrm)
  list(da = da, db = db)
}


# =============================================================================
# Complementarity residual for MCP
# =============================================================================


# =============================================================================
# Merit function
# =============================================================================

#' Squared FB residual merit function
#'
#' θ(Y) = ½ Σ_t Σ_j φ(a_{j,t}, F_j(t))²
#'
#' The squared L2 norm of all FB residuals.  Used as the merit function
#' for backtracking line search.  Zero iff all complementarity conditions
#' are satisfied.
#'
#' @param fb_mat  n_spec × T matrix of FB residuals
#' @return Scalar: θ = ½ Σ φ²
#' @noRd
mcp_compute_merit <- function(fb_mat) {
  0.5 * sum(fb_mat * fb_mat, na.rm = TRUE)
}


#' Natural (min-max) complementarity residual
#'
#' An alternative to FB: ψ(a, b) = min(a, b).  Simpler but has weaker
#' differentiability properties.  Used for diagnostic comparison only.
#'
#' ψ(a, b) = 0  ⇔  a ≥ 0, b ≥ 0, a·b = 0  (same as FB)
#'
#' @param a,b Numeric vectors
#' @return Numeric vector of ψ values
#' @noRd
mcp_natural_residual <- function(a, b) {
  pmin(a, b)
}


# =============================================================================
# Active-set extraction
# =============================================================================

#' Extract binding pattern from a solved MCP path
#'
#' Given a converged MCP path, determines which constraints are binding
#' at each period.  A constraint is binding at period t if the constrained
#' variable is at (or very near) its bound.
#'
#' The result is a regime bitfield vector compatible with the existing
#' \code{obc_regime_idx} format: bit j is set if constraint j binds at period t.
#'
#' @param Y         T × n_endo solution matrix
#' @param mcp_specs List of MCP specs (from mcp_parse_tags)
#' @param tol       Numeric tolerance for bound proximity (default 1e-8)
#' @return Integer vector (length T): bitfield per period
#' @noRd
mcp_assemble_active_set <- function(Y, mcp_specs, tol = 1e-8) {
  n_spec <- length(mcp_specs)
  T <- nrow(Y)

  if (n_spec == 0L) return(integer(T))

  regime_path <- integer(T)

  for (t in seq_len(T)) {
    bits <- 0L
    for (j in seq_len(n_spec)) {
      sp  <- mcp_specs[[j]]
      val <- Y[t, sp$var_idx]

      if (sp$op == ">") {
        # Lower bound: binding if x is close to bound
        if (abs(val - sp$bound) < tol) bits <- bits + 2L^(j - 1L)
      } else {
        # Upper bound: binding if x is close to bound
        if (abs(val - sp$bound) < tol) bits <- bits + 2L^(j - 1L)
      }
    }
    regime_path[t] <- bits
  }

  regime_path
}


# =============================================================================
# Column metadata for sparse block-tridiagonal assembly
# =============================================================================

#' Precompute column metadata for MCP stacked Jacobian assembly
#'
#' Extends \code{.occbin_col_meta()} to also track which Jacobian columns
#' correspond to MCP-constrained variables.  Used by \code{mcp_solve_path()}
#' to correctly route FB derivative contributions into the sparse matrix.
#'
#' @param dyn       compiled$dynamic (from build_dynamic_model)
#' @param mcp_specs List of MCP specs (from mcp_parse_tags)
#' @return List with all fields from \code{.occbin_col_meta()} plus:
#'   $cur_var_to_spec — integer vector: for each current-period endo variable,
#'     which spec index modifies it (0 = unconstrained)
#'   $spec_cur_dc     — integer vector: dyn column index for each spec's
#'     variable in the current period
#' @noRd
mcp_col_meta <- function(dyn, mcp_specs) {
  base <- .occbin_col_meta(dyn)
  n_spec <- length(mcp_specs)
  n_endo <- base$n_endo

  # Build per-variable spec mapping
  cur_var_to_spec <- integer(n_endo)
  spec_cur_dc     <- integer(n_spec)

  if (n_spec > 0L) {
    cmap <- dyn$dyn_col_map
    for (j in seq_len(n_spec)) {
      vi <- mcp_specs[[j]]$var_idx
      cur_var_to_spec[vi] <- j

      # Find dyn column for var__0 (current period)
      k <- which(cmap$name == mcp_specs[[j]]$var_name & cmap$lead_lag == 0L)
      if (length(k) > 0L) {
        spec_cur_dc[j] <- cmap$col[k[1L]]
      }
    }
  }

  c(base, list(
    cur_var_to_spec = cur_var_to_spec,
    spec_cur_dc     = spec_cur_dc
  ))
}


# =============================================================================
# Convergence check
# =============================================================================

#' Check MCP solver convergence
#'
#' Evaluates whether the MCP solver has converged by checking both the
#' stacked Newton residual and the FB complementarity residual.
#'
#' @param R_stack  Numeric vector: stacked Newton residual (T × n_endo)
#' @param fb_mat   n_spec × T matrix of FB residuals
#' @param tol      Numeric tolerance (default 1e-8)
#' @return List with:
#'   $converged    — Logical
#'   $max_res      — max|R_stack|
#'   $max_fb       — max|fb_mat|
#'   $fb_zero      — Logical: TRUE if FB residuals are all near zero
#' @noRd
mcp_check_convergence <- function(R_stack, fb_mat, tol = 1e-8) {
  max_res <- max(abs(R_stack))
  max_fb  <- if (length(fb_mat) > 0L) max(abs(fb_mat), na.rm = TRUE) else 0

  list(
    converged = max_res < tol && max_fb < tol,
    max_res   = max_res,
    max_fb    = max_fb,
    fb_zero   = max_fb < tol
  )
}
