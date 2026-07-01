## R/sequence-jacobian.R
## --------------------------------------------------------------------------
## Representative-agent sequence-space Jacobian (RA-SSJ) spike.
##
## For a model's equilibrium map H(X, eps) = 0 over horizon T (X stacks
## endogenous paths, eps is the exogenous shock path), the sequence-space
## Jacobian is:
##
##   J = dH/dX  (T*n_endo x T*n_endo, block-tridiagonal)
##
## evaluated at the deterministic steady state (all periods at ss, eps = 0).
## This is IDENTICAL to the stacked Jacobian already assembled by the
## perfect-foresight Newton solver at convergence.
##
## The general-equilibrium impulse response to a unit shock at t=1 is:
##
##   dX = -J^{-1} * (dH/deps_1) * scale
##
## where dH/deps_1 is a T*n_endo vector: in period 1 it is the partial of
## the dynamic residuals w.r.t. eps_a (extracted from jacobian_fn's exo
## columns); in all other periods it is zero (AR shock enters through the
## state, not directly as an exogenous perturbation in t > 1).
##
## Validation oracle: the SSJ IRF must match the order-1 perturbation IRF
## (ghx/ghu recursion) to ~1e-8 (exact linear identity for a linear model).
##
## Reuses WITHOUT MODIFICATION:
##   .pf_col_meta()             (R/pf-newton.R)
##   .pf_build_sparse_system()  (R/pf-newton.R)
##   .pf_make_dy_exo()          (R/perfect-foresight-solve.R)
##   pf_newton_solve()          (R/pf-newton.R)
##   compute_irfs()             (R/stochsimul-monolith.R)
## --------------------------------------------------------------------------


# =============================================================================
# Constructor: build the steady-state sequence Jacobian
# =============================================================================

#' Build a representative-agent sequence-space Jacobian (RA-SSJ) object
#'
#' Assembles the T*n_endo x T*n_endo stacked Jacobian \eqn{J = \partial H /
#' \partial X} evaluated at the deterministic steady state and the shock-impact
#' vector \eqn{\partial H / \partial \varepsilon_1} at t = 1.  Both are
#' extracted from the same compiled-derivative machinery used by
#' \code{\link{pf_newton_solve}}.
#'
#' The resulting object supports one operation, \code{\link{sj_irf}}, which
#' computes the general-equilibrium impulse response
#' \eqn{dX = -J^{-1} (\partial H/\partial\varepsilon_1) \cdot \text{scale}}.
#'
#' @param model     A parsed dynhr model (\code{\link{parse_mod}}).
#' @param compiled  A compiled model (\code{\link{compile_model}}).
#' @param ss        Named numeric vector: deterministic steady state
#'   (from \code{\link{solve_steady}}).
#' @param params    Named numeric parameter vector (default: \code{model$param_values}).
#' @param horizon   Integer: sequence horizon T (default 40).
#'
#' @return An object of class \code{ra_ssj} with components:
#'   \describe{
#'     \item{\code{J}}{Sparse dgCMatrix (T*n_endo x T*n_endo): dH/dX at ss.}
#'     \item{\code{dH_deps}}{Named list, one entry per shock: each a numeric
#'       vector of length T*n_endo giving dH/deps at t=1.}
#'     \item{\code{J_lu}}{LU factorisation of J (via \code{Matrix::lu}).}
#'     \item{\code{T}}{Integer: horizon.}
#'     \item{\code{n_endo}}{Integer: number of endogenous variables.}
#'     \item{\code{endo_names}}{Character: endogenous variable names.}
#'     \item{\code{exo_names}}{Character: exogenous variable names.}
#'   }
#'
#' @details
#' \strong{Sign/transpose/ordering conventions}
#'
#' The PF Newton system solves \eqn{H(X) = 0} where the stacked residual is
#' indexed as \eqn{H_{(t-1) n + i}} = residual of equation \eqn{i} at period
#' \eqn{t}.  The stacked state vector \eqn{X} is ordered identically.
#'
#' The steady-state Jacobian is built by evaluating
#' \code{.pf_build_sparse_system} at the all-ss path with zero shocks.  At the
#' steady state \eqn{H(X_{ss}) = 0} (by definition), so the assembled residual
#' is zero; only the Jacobian is used.
#'
#' The shock-impact vector \code{dH_deps} is extracted by calling
#' \code{jacobian_fn} at t = 1 (at the ss path) and reading the columns
#' corresponding to exogenous variables.  For periods t > 1 the direct shock
#' impact is zero (the shock enters future periods only through lagged state
#' transitions, which are already captured in \eqn{J}).
#'
#' The GE impulse response formula \eqn{dX = -J^{-1} \partial H/\partial\varepsilon}
#' follows from total differentiation of \eqn{H(X, \varepsilon) = 0}:
#' \eqn{J \, dX + (\partial H/\partial\varepsilon) \, d\varepsilon = 0}.
#'
#' @seealso \code{\link{sj_irf}}
#' @export
ra_ssj <- function(model, compiled, ss, params = NULL, horizon = 40L) {
  if (is.null(params)) params <- model$param_values

  dyn    <- compiled$dynamic
  n_endo <- length(dyn$endo_names)
  n_exo  <- length(dyn$exo_names)
  T      <- as.integer(horizon)

  # Align ss to endo_names order (same convention as pf_newton_solve)
  y_ss_num <- as.numeric(ss[dyn$endo_names])

  # Build a T x n_endo path at the steady state (all periods = ss)
  Y_ss <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
  colnames(Y_ss) <- dyn$endo_names

  # Zero shock path
  eps_mat <- matrix(0, nrow = T, ncol = n_exo)
  colnames(eps_mat) <- dyn$exo_names

  # Precompute column metadata (shared with pf_newton_solve)
  meta <- .pf_col_meta(dyn)

  # Preallocate triplet buffer
  trip_buf <- .pf_preallocate_triplets(dyn, T)

  # Build the stacked Jacobian J = dH/dX at the ss path with zero shocks.
  # The named steady-state vector (for STEADY_STATE() references)
  y_ss_named <- setNames(y_ss_num, dyn$endo_names)

  sys <- .pf_build_sparse_system(
    Y        = Y_ss,
    y0_num   = y_ss_num,
    y_ss_num = y_ss_num,
    eps_mat  = eps_mat,
    meta     = meta,
    dyn      = dyn,
    params   = params,
    y_ss     = y_ss_named,
    T        = T,
    n_eq     = n_endo,
    n_endo   = n_endo,
    obc_specs    = list(),
    regime       = matrix(FALSE, nrow = 1L, ncol = T),
    spec_cur_dc  = integer(0),
    trip_buf     = trip_buf
  )

  J <- sys$J  # T*n_endo x T*n_endo sparse Jacobian

  # -- Build dH/deps for each shock: T*n_endo vector --
  # At t=1, evaluate jacobian_fn at the ss dy vector with zero shocks.
  # The exo columns of Jt[, dc] give dH_t / deps (rows = equations at t=1;
  # col = dynamic column for the exo variable).
  # For t > 1, the direct shock impact is zero.
  dy1 <- .pf_make_dy_exo(meta, Y_ss, y_ss_num, y_ss_num, eps_mat, 1L, T, n_exo)
  Jt1 <- dyn$jacobian_fn(dy1, params, y_ss_named)

  # dH_deps: for each shock, a T*n_endo vector (nonzero only at rows 1..n_endo)
  dH_deps <- vector("list", n_exo)
  names(dH_deps) <- dyn$exo_names

  for (k in seq_len(n_exo)) {
    exo_nm  <- dyn$exo_names[k]
    vec     <- numeric(T * n_endo)  # zero everywhere

    # Find all dynamic columns for this exo variable (at any lead/lag)
    # In the standard linear RBC, eps_a enters at lag 0 only.
    # We sum over all exo dynamic columns for this variable.
    for (j in seq_along(meta$kind)) {
      if (meta$kind[j] != "exo") next
      if (meta$var_idx[j] != k)  next
      dc <- meta$dyn_col[j]
      # Add partial of equations at t=1 w.r.t. this exo column
      vec[seq_len(n_endo)] <- vec[seq_len(n_endo)] + Jt1[, dc]
    }

    dH_deps[[exo_nm]] <- vec
  }

  # LU factorisation of J (reused across sj_irf calls)
  J_lu <- Matrix::lu(J)

  structure(
    list(
      J          = J,
      dH_deps    = dH_deps,
      J_lu       = J_lu,
      T          = T,
      n_endo     = n_endo,
      endo_names = dyn$endo_names,
      exo_names  = dyn$exo_names
    ),
    class = "ra_ssj"
  )
}


# =============================================================================
# sj_irf: compute the GE impulse response via SSJ inversion
# =============================================================================

#' Compute the general-equilibrium impulse response via sequence-space Jacobian
#'
#' Given an \code{\link{ra_ssj}} object, computes the linear impulse response
#' \deqn{dX = -J^{-1} \frac{\partial H}{\partial \varepsilon_1} \cdot \text{scale}}
#' and returns it as a T x n_endo matrix in the same format as
#' \code{\link{compute_irfs}} (rows = periods, cols = endogenous variables).
#'
#' @param ssj    An \code{ra_ssj} object from \code{\link{ra_ssj}}.
#' @param shock  Character: name of the shock (must be in \code{ssj$exo_names}).
#' @param scale  Numeric scalar: shock size (default 1.0).  Set to
#'   \code{model$param_values["sigma_a"]} to match the one-standard-deviation
#'   IRF convention used by \code{compute_irfs}.
#'
#' @return A T x n_endo numeric matrix (rows = periods, cols = endo vars),
#'   with column names equal to \code{ssj$endo_names}.  Compatible with
#'   the slot format of an \code{IRFCollection} object returned by
#'   \code{compute_irfs}.
#'
#' @details
#' The GE response formula follows from total differentiation of
#' \eqn{H(X, \varepsilon) = 0} at the steady state:
#' \deqn{J \, dX + \frac{\partial H}{\partial \varepsilon} \, d\varepsilon = 0
#'       \;\Rightarrow\; dX = -J^{-1} \frac{\partial H}{\partial \varepsilon} \, d\varepsilon.}
#'
#' @seealso \code{\link{ra_ssj}}, \code{\link{compute_irfs}}
#' @export
sj_irf <- function(ssj, shock, scale = 1.0) {
  stopifnot(inherits(ssj, "ra_ssj"))
  if (!shock %in% ssj$exo_names)
    stop(sprintf("sj_irf: shock '%s' not in exo_names (%s)",
                 shock, paste(ssj$exo_names, collapse = ", ")))

  dH <- ssj$dH_deps[[shock]] * scale  # T*n_endo vector

  # Solve J * dX = -dH  => dX = -J^{-1} dH
  dX_vec <- as.numeric(Matrix::solve(ssj$J_lu, -dH))

  # Reshape to T x n_endo (row t = period t deviation from ss).
  # dX_vec is indexed: dX_vec[(t-1)*n_endo + i] = deviation of var i at period t.
  # matrix(..., byrow=TRUE) fills row-major: row t = dX_vec[(t-1)*n_endo + 1..n_endo].
  T      <- ssj$T
  n_endo <- ssj$n_endo
  dX <- matrix(dX_vec, nrow = T, ncol = n_endo, byrow = TRUE)
  colnames(dX) <- ssj$endo_names
  dX
}


# =============================================================================
# print method
# =============================================================================

#' @export
print.ra_ssj <- function(x, ...) {
  cat(sprintf(
    "<ra_ssj>  horizon T=%d  n_endo=%d  n_exo=%d\n",
    x$T, x$n_endo, length(x$exo_names)
  ))
  cat("  endo:", paste(x$endo_names, collapse = ", "), "\n")
  cat("  exo: ", paste(x$exo_names,  collapse = ", "), "\n")
  invisible(x)
}
