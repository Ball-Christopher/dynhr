## R/mcp-solve.R
## --------------------------------------------------------------------------
## MCP (Mixed Complementarity Problem) path solver via semi-smooth Newton
## with Fischer-Burmeister reformulation and sparse block-tridiagonal linear
## algebra.
##
## Provides:
##   mcp_solve_path()  -- main MCP stacked Newton path solver
##
## ALGORITHM
##
##   Reformulate the complementarity conditions using the Fischer-Burmeister
##   function φ(a, b) = a + b − √(a² + b²).  For each period t and MCP
##   constraint j on variable x with bound b:
##
##     Lower bound (op = ">"):  φ(x_t − b,  F_j(y_{t-1}, y_t, y_{t+1}, ε_t)) = 0
##     Upper bound (op = "<"):  φ(b − x_t, −F_j(y_{t-1}, y_t, y_{t+1}, ε_t)) = 0
##
##   The stacked system is T·n_endo equations in T·n_endo unknowns, solved via
##   semi-smooth Newton with a sparse block-tridiagonal Jacobian.
##
##   The Jacobian has the standard block-tridiagonal structure, with MCP rows
##   modified by the FB chain rule:
##     J_new[i, :] = ∂φ/∂b · J_orig[i, :] + δ_{i, var} · ∂φ/∂a
##   where δ_{i, var} is 1 in the column of the constrained variable.
##
## REFERENCES
##   Fischer, A. (1992). "A special Newton-type optimization method."
##     Optimization, 24(3-4), 269-284.
##   De Luca, T., Facchinei, F., & Kanzow, C. (1996). "A semismooth equation
##     approach to the solution of nonlinear complementarity problems."
##     Mathematical Programming, 75(3), 407-439.
##   Ferris, M. C., & Munson, T. S. (2000). "Complementarity problems in
##     GAMS and the PATH solver." Journal of Economic Dynamics and Control.
##   Guerrieri, L., & Iacoviello, M. (2015). "OccBin: A toolkit for solving
##     dynamic models with occasionally binding constraints easily."
##     Journal of Monetary Economics, 70, 22-38.
## --------------------------------------------------------------------------


# =============================================================================
# Build the stacked system with FB-modified residuals and Jacobian
# =============================================================================

#' Build the stacked MCP system with FB-modified residuals and Jacobian
#'
#' For each period t:
#'   1. Evaluate the full dynamic residual F_t and Jacobian J_t.
#'   2. For MCP-constrained equations, replace the residual with the FB
#'      complementarity function and adjust the Jacobian row via the
#'      FB chain rule.
#'   3. Assemble into a sparse block-tridiagonal dgCMatrix.
#'
#' @param Y            T × n_endo numeric matrix: current path
#' @param y0_num       Length-n_endo initial state at t=0
#' @param y_ss_num     Length-n_endo steady state (terminal condition)
#' @param eps_mat      T × n_exo shock matrix
#' @param pf_meta      Column metadata from .pf_col_meta()
#' @param mcp_meta     Column metadata from mcp_col_meta() (with spec info)
#' @param mcp_specs    List of MCP specs (from mcp_parse_tags)
#' @param compiled     dynhr_compiled (for residuals_fn, jacobian_fn)
#' @param params       Named numeric parameter vector
#' @param T, n_endo    Integer: horizon and number of endogenous variables
#' @return List with:
#'   $J   — sparse dgCMatrix (T*n_endo × T*n_endo)
#'   $R   — numeric vector (T*n_endo): stacked residual with FB modifications
#'   $fb  — n_spec × T matrix of FB residuals (for convergence check)
#' @noRd
.mcp_build_stacked_system <- function(Y, y0_num, y_ss_num, eps_mat,
                                       pf_meta, mcp_meta,
                                       mcp_specs, compiled,
                                       params, T, n_endo) {
  dyn      <- compiled$dynamic
  n_spec   <- length(mcp_specs)
  n_total  <- T * n_endo
  n_eq     <- dyn$n_eq

  # Preallocate
  R <- numeric(n_total)

  # Sparse triplet storage
  i_triplet <- integer(0)
  j_triplet <- integer(0)
  v_triplet <- numeric(0)

  # FB residual matrix
  fb_mat <- if (n_spec > 0L) matrix(0, n_spec, T) else matrix(0, 0, T)

  # Determine which equation rows to use: we need exactly n_endo rows per period.
  # If n_eq == n_endo: use all rows (standard case).
  # If n_eq > n_endo: select n_endo rows — MCP-specified equations come first,
  # then neutral equations fill the rest.  This handles bind/relax models
  # where n_eq > n_endo due to equation variants.
  # If n_eq < n_endo: pad with zero rows (shouldn't happen in practice).
  eq_selector <- if (n_eq == n_endo) {
    seq_len(n_endo)
  } else if (n_eq > n_endo) {
    # Collect MCP-specified eq indices
    mcp_eqs <- sort(unique(vapply(mcp_specs, `[[`, integer(1), "eq_idx")))
    # Neutral equations: indices in 1..n_eq that are not MCP equations
    neutral_eqs <- setdiff(seq_len(n_eq), mcp_eqs)
    # Select: all MCP equations + enough neutral equations to reach n_endo
    needed <- n_endo - length(mcp_eqs)
    selected <- c(mcp_eqs, neutral_eqs[seq_len(min(needed, length(neutral_eqs)))])
    # If we still don't have n_endo, pad
    if (length(selected) < n_endo) {
      selected <- c(selected, setdiff(seq_len(n_eq), selected)[seq_len(n_endo - length(selected))])
    }
    selected[seq_len(n_endo)]
  } else {
    # n_eq < n_endo: pad with zeros
    seq_len(n_eq)
  }

  for (t in seq_len(T)) {
    row_off <- (t - 1L) * n_endo
    col_off <- (t - 1L) * n_endo

    # Build dy vector for this period
    dy <- .pf_make_dy(pf_meta, Y, y0_num, y_ss_num, eps_mat[t, ], t, T)

    # Compute full residual and Jacobian
    Rt <- dyn$residuals_fn(dy, params, NULL)
    Jt <- dyn$jacobian_fn(dy, params, NULL)

    # Select the n_endo rows we actually use
    if (length(eq_selector) != n_endo || any(eq_selector != seq_len(n_endo))) {
      Rt <- Rt[eq_selector]
      Jt <- Jt[eq_selector, , drop = FALSE]
    }

    # Remap: MCP spec eq_idx -> row in the selected system
    # Find which selected row corresponds to each spec's eq_idx
    spec_to_selected_row <- integer(n_spec)
    for (j in seq_len(n_spec)) {
      spec_to_selected_row[j] <- match(mcp_specs[[j]]$eq_idx, eq_selector)
    }

    # Apply FB modifications to MCP-constrained equations
    if (n_spec > 0L) {
      for (j in seq_len(n_spec)) {
        sp <- mcp_specs[[j]]
        sr <- spec_to_selected_row[j]   # selected row index (in 1..n_endo)
        vi <- sp$var_idx

        # Current value of constrained variable
        x_cur <- Y[t, vi]

        if (sp$op == ">") {
          # Lower bound: a = x - bound, b = F(x)
          a <- x_cur - sp$bound
          b <- Rt[sr]
        } else {
          # Upper bound: a = bound - x, b = -F(x)
          a <- sp$bound - x_cur
          b <- -Rt[sr]
        }

        # Store FB residual
        fb <- mcp_fb(a, b)
        fb_mat[j, t] <- fb

        # Replace original residual with FB residual
        Rt[sr] <- fb

        # Compute FB derivatives
        deriv <- mcp_fb_deriv(a, b)
        dphi_da <- deriv$da
        dphi_db <- deriv$db

        # If both derivatives are ~0 (at the kink), skip row modification.
        # isTRUE guards NaN derivatives from invalid (e.g. C<0) trial steps;
        # such steps are penalised by the merit function and discarded.
        if (isTRUE(abs(dphi_da) < .Machine$double.eps &&
                   abs(dphi_db) < .Machine$double.eps)) next

        # Modify Jacobian row sr:
        #   J_new[sr, :] = dphi_db * J_orig[sr, :]
        #   J_new[sr, vi_col] += dphi_da * 1
        Jt[sr, ] <- dphi_db * Jt[sr, ]

        # Add dphi_da to the current-period column of the constrained variable
        dc_vi <- mcp_meta$spec_cur_dc[j]
        if (dc_vi > 0L && dc_vi <= ncol(Jt)) {
          Jt[sr, dc_vi] <- Jt[sr, dc_vi] + dphi_da
        }
      }
    }

    # Store residual and route Jacobian columns into sparse triplets
    R[row_off + seq_len(n_endo)] <- Rt

    # Route Jacobian columns to global sparse entries
    # Lag columns -> block (t, t-1)
    for (k in seq_along(mcp_meta$lag_cols)) {
      dc  <- mcp_meta$lag_cols[k]
      vi  <- mcp_meta$lag_var[k]
      gcol <- if (t > 1L) (t - 2L) * n_endo + vi else NA_integer_
      if (is.na(gcol)) next
      for (eq in seq_len(n_endo)) {
        val <- Jt[eq, dc]
        if (isTRUE(val != 0) && is.finite(val)) {
          i_triplet <- c(i_triplet, row_off + eq)
          j_triplet <- c(j_triplet, gcol)
          v_triplet <- c(v_triplet, val)
        }
      }
    }

    # Current columns -> block (t, t)
    for (k in seq_along(mcp_meta$cur_cols)) {
      dc  <- mcp_meta$cur_cols[k]
      vi  <- mcp_meta$cur_var[k]
      gcol <- col_off + vi
      for (eq in seq_len(n_endo)) {
        val <- Jt[eq, dc]
        if (isTRUE(val != 0) && is.finite(val)) {
          i_triplet <- c(i_triplet, row_off + eq)
          j_triplet <- c(j_triplet, gcol)
          v_triplet <- c(v_triplet, val)
        }
      }
    }

    # Lead columns -> block (t, t+1)
    for (k in seq_along(mcp_meta$lead_cols)) {
      dc  <- mcp_meta$lead_cols[k]
      vi  <- mcp_meta$lead_var[k]
      gcol <- if (t < T) t * n_endo + vi else NA_integer_
      if (is.na(gcol)) next
      for (eq in seq_len(n_endo)) {
        val <- Jt[eq, dc]
        if (isTRUE(val != 0) && is.finite(val)) {
          i_triplet <- c(i_triplet, row_off + eq)
          j_triplet <- c(j_triplet, gcol)
          v_triplet <- c(v_triplet, val)
        }
      }
    }
  }

  # Build sparse matrix from triplets — repr = "C" replaces deprecated
  # giveCsparse = TRUE (Landmine 5).
  if (length(i_triplet) == 0L) {
    J <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 0,
                              dims = c(n_total, n_total),
                              repr = "C")
  } else {
    J <- Matrix::sparseMatrix(
      i = i_triplet, j = j_triplet, x = v_triplet,
      dims = c(n_total, n_total),
      index1 = TRUE,  # R-style 1-based indexing
      repr = "C"
    )
  }

  list(J = J, R = R, fb = fb_mat)
}


# =============================================================================
# Backtracking line search
# =============================================================================

#' Armijo backtracking line search for MCP merit function
#'
#' Finds step length α ∈ (0, 1] such that:
#'   θ(Y + α·Δ) ≤ θ(Y) + σ·α·∇θ(Y)·Δ
#' where θ(Y) = ½ Σ φ(a, F)² is the squared FB merit function.
#'
#' The gradient ∇θ = J_stack' · R_stack is computed from the current
#' Jacobian and residual.
#'
#' @param Y           T × n_endo matrix: current path
#' @param delta_vec   Numeric vector (length T*n_endo): Newton direction
#' @param R_stack     Numeric vector (length T*n_endo): current residual
#' @param J_stack     dgCMatrix: current Jacobian
#' @param theta_cur   Scalar: current merit function value
#' @param fn_merit    Function to evaluate θ at a new Y
#' @param sigma       Armijo parameter (default 1e-4)
#' @param max_ls      Maximum line search iterations (default 20)
#' @param n_endo      Integer: number of endogenous variables
#' @param T           Integer: horizon
#' @return List with:
#'   $alpha      — step length
#'   $theta_new  — merit value at Y + alpha*delta
#'   $Y_new      — updated path matrix (or NULL if no tried step reduced merit)
#'   $ls_iter    — iterations used
#'   $accepted   — logical: TRUE if the Armijo condition was met
#'
#' @details
#' On Armijo failure the search returns the smallest-merit trial point actually
#' evaluated (when it strictly improves on the incumbent), rather than \code{NULL}
#' with a discarded step. The previous behaviour — return NULL, caller applies a
#' blind half Newton step — could blow the path up and singularise the next
#' Jacobian on cold starts (issue M15).
#' @noRd
.mcp_line_search <- function(Y, delta_vec, R_stack, J_stack,
                              theta_cur, fn_merit,
                              sigma = 1e-4, max_ls = 20L,
                              n_endo, T) {
  # Compute gradient: ∇θ = J' · R
  grad <- as.numeric(Matrix::crossprod(J_stack, R_stack))
  directional_deriv <- sum(grad * delta_vec)

  # If directional derivative is positive, Newton direction is not a descent
  # direction.  Fall back to steepest descent: Δ = -∇θ
  if (directional_deriv >= 0) {
    delta_vec <- -grad
    directional_deriv <- -sum(grad * grad)
  }

  alpha <- 1.0

  # Track the best merit-decreasing trial point seen during backtracking (M15).
  best_theta <- theta_cur
  best_Y     <- NULL
  best_alpha <- 0

  for (ls_iter in seq_len(max_ls)) {
    # Trial point
    Y_trial <- Y
    for (t in seq_len(T)) {
      idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
      Y_trial[t, ] <- Y[t, ] + alpha * delta_vec[idx_t]
    }

    theta_new <- fn_merit(Y_trial)

    if (is.finite(theta_new) && theta_new < best_theta) {
      best_theta <- theta_new
      best_Y     <- Y_trial
      best_alpha <- alpha
    }

    # Armijo condition (isTRUE guards NaN merit values from invalid steps).
    if (isTRUE(theta_new <= theta_cur + sigma * alpha * directional_deriv)) {
      return(list(alpha = alpha, theta_new = theta_new,
                  Y_new = Y_trial, ls_iter = ls_iter, accepted = TRUE))
    }

    alpha <- alpha * 0.5
  }

  # Armijo never satisfied: return the best vetted (merit-decreasing) point.
  list(alpha = best_alpha, theta_new = best_theta,
       Y_new = best_Y, ls_iter = max_ls, accepted = FALSE)
}


# =============================================================================
# Main MCP solver
# =============================================================================

#' Semi-smooth Newton path solver for MCP-constrained DSGE models
#'
#' Solves the T-period deterministic path for a DSGE model with occasionally
#' binding constraints expressed as Mixed Complementarity Problems (MCPs).
#'
#' The solver uses a **semi-smooth Newton method** with the **Fischer-Burmeister
#' complementarity function** to handle lower and upper bounds on endogenous
#' variables.  The stacked system (T × n_endo equations) is assembled as a
#' sparse block-tridiagonal dgCMatrix and solved via \code{Matrix::sparseQR}.
#'
#' This is the MCP analogue of Dynare's \code{perfect_foresight_solver} with
#' \code{stack_solve_algo = 7} (PATH solver), but using sparse linear algebra
#' and the FB reformulation instead of PATH's dense banded LU.
#'
#' @param compiled  dynhr_compiled (from \code{\link{compile_model}})
#' @param y0        Named numeric vector: endogenous state at t=0.
#' @param y_ss      Named numeric vector: steady state (terminal condition).
#' @param shock_path  T × n_exo numeric matrix of structural shocks.  Column
#'   names must match \code{varexo_names}; missing columns are zero.
#' @param params    Named numeric parameter vector.
#' @param mcp_specs List of MCP specs from \code{\link{mcp_parse_tags}} or
#'   \code{obc_parse_tags}.  Empty list → no constraints (standard Newton).
#' @param max_iter       Maximum Newton iterations (default 50).
#' @param tol            Newton convergence tolerance on max|R| (default 1e-8).
#' @param step_size      Initial Newton step length (default 1.0).
#' @param line_search    Logical: perform Armijo backtracking (default TRUE).
#' @param sigma_line     Armijo parameter for line search (default 1e-4).
#' @param sparse_fallback Logical: use dense \code{solve()} if sparse fails
#'   (default TRUE).
#' @param homotopy       Continuation strategy for cold starts. \code{"auto"}
#'   (default) tries a direct solve first and, only if it fails, escalates
#'   through a ladder of continuation runs (4, 16, 64 stages) that ramp the
#'   shock path from zero (steady state) to its full value, warm-starting each
#'   stage. So easy problems are unchanged and hard cold starts are rescued
#'   (issue M15). An integer \code{>= 1} forces exactly that many stages;
#'   \code{1} is a single direct solve.
#' @param method         Sparse-solve method: \code{"sparse"} (uses
#'   \code{Matrix::sparseQR}), \code{"dense"} (uses \code{solve}), or
#'   \code{"auto"} (chooses dense for \code{T * n_endo < 100}, sparse otherwise).
#' @param backend        Sparse solve backend: \code{"Matrix"} (default, uses
#'   \code{Matrix::solve} on dgCMatrix) or \code{"Rcpp"} (uses Armadillo
#'   \code{spsolve} via SuperLU, typically 2-5x faster for systems > 1000×1000).
#' @param Y_init         Optional initial path (T × n_endo matrix) for
#'   warm-starting.  NULL → initialize at steady state (default).
#' @return List with:
#'   \item{Y}{T × n_endo solution matrix (rows = periods, cols = variables).}
#'   \item{irf}{T × n_endo deviation from steady state.}
#'   \item{converged}{Logical: TRUE if Newton converged.}
#'   \item{n_iter}{Integer: Newton iterations used.}
#'   \item{max_res}{Numeric: final max|R_stack|.}
#'   \item{fb_final}{n_spec × T matrix of final FB residuals.}
#'   \item{active_set}{Integer vector (length T): per-period binding bitfield.}
#'   \item{endo_names}{Character: variable ordering of Y columns.}
#'   \item{merit_history}{Numeric vector: merit function value per iteration.}
#' @export
#'
#' @examples
#' \dontrun{
#' model <- parse_mod("nk_zlb_dynare.mod")
#' compiled <- compile_model(model)
#' ss <- solve_steady_state(model, compiled, model$param_values)
#' specs <- mcp_parse_tags(model)
#' shock_path <- matrix(0, nrow = 40, ncol = length(model$varexo_names))
#' shock_path[1, 1] <- -0.01  # negative demand shock
#' result <- mcp_solve_path(compiled, ss$values, ss$values,
#'                           shock_path, model$param_values, specs)
#' }
mcp_solve_path <- function(compiled,
                            y0,
                            y_ss,
                            shock_path,
                            params,
                            mcp_specs       = list(),
                            max_iter        = 50L,
                            tol             = 1e-8,
                            step_size       = 1.0,
                            line_search     = TRUE,
                            sigma_line      = 1e-4,
                            sparse_fallback = TRUE,
                            homotopy        = "auto",
                            method          = c("sparse", "dense", "auto"),
                            backend         = c("Matrix", "Rcpp"),
                            Y_init          = NULL) {

  dyn    <- compiled$dynamic
  n_endo <- length(dyn$endo_names)
  n_eq   <- dyn$n_eq
  n_spec <- length(mcp_specs)

  # Validate: with MCP constraints, we need n_eq >= n_endo - n_distinct_specs
  # (some equations get replaced by FB constraints)
  if (n_spec > 0L && n_eq < n_endo) {
    warning(sprintf(
      "mcp_solve_path: n_eq (%d) < n_endo (%d). MCP constraints replace equations,",
      n_eq, n_endo),
      " but fewer equations than variables suggests an underdetermined system.")
  }

  # Normalize shock_path
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

  # Align y0 and y_ss
  y0_num   <- as.numeric(y0[dyn$endo_names])
  y_ss_num <- as.numeric(y_ss[dyn$endo_names])

  # Precompute column metadata
  pf_meta  <- .pf_col_meta(dyn)
  mcp_meta <- mcp_col_meta(dyn, mcp_specs)

  # Select linear algebra method and backend
  method  <- match.arg(method)
  backend <- match.arg(backend)
  if (method == "auto") {
    method <- if (T * n_endo < 100L) "dense" else "sparse"
  }

  # Initialize path: warm-start or steady state
  if (!is.null(Y_init)) {
    if (nrow(Y_init) != T || ncol(Y_init) != n_endo) {
      stop(sprintf("Y_init must be %d x %d (T x n_endo), got %d x %d",
                   T, n_endo, nrow(Y_init), ncol(Y_init)))
    }
    Y <- Y_init
    colnames(Y) <- dyn$endo_names
  } else {
    Y <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
    colnames(Y) <- dyn$endo_names
  }

  # ==========================================================================
  # Inner semi-smooth Newton solve, parameterised by the shock path + warm
  # start. Reused by the homotopy/continuation driver below. Globalised with an
  # Armijo line search, a best-vetted-point fallback on line-search failure,
  # and a Levenberg-Marquardt step when the FB Jacobian singularises (M15).
  # ==========================================================================
  .mcp_run_newton <- function(Y_start, eps_use) {
    Y_loc <- Y_start
    colnames(Y_loc) <- dyn$endo_names

    merit_fn <- function(Y_try) {
      ss <- .mcp_build_stacked_system(
        Y_try, y0_num, y_ss_num, eps_use,
        pf_meta, mcp_meta, mcp_specs, compiled,
        params, T, n_endo
      )
      if (anyNA(ss$R) || any(!is.finite(ss$R))) return(1e30)  # invalid step
      0.5 * sum(ss$R * ss$R)
    }

    conv_loc    <- FALSE
    n_iter_loc  <- 0L
    max_res_loc <- Inf
    fb_loc      <- matrix(0, n_spec, T)
    merit_loc   <- numeric(0)

    for (iter in seq_len(max_iter)) {
      n_iter_loc <- iter

      sys <- .mcp_build_stacked_system(
        Y_loc, y0_num, y_ss_num, eps_use,
        pf_meta, mcp_meta, mcp_specs, compiled,
        params, T, n_endo
      )

      # Non-finite residuals (e.g. C^(-sigma) with C<0 after an over-large
      # cold-start step): abort this run cleanly so the homotopy driver can
      # retry with finer staging. Guards against `if (NA)` downstream.
      if (anyNA(sys$R) || any(!is.finite(sys$R))) {
        max_res_loc <- Inf
        break
      }

      max_res     <- max(abs(sys$R))
      max_res_loc <- max_res
      fb_loc      <- sys$fb

      theta_cur <- 0.5 * sum(sys$R * sys$R)
      merit_loc <- c(merit_loc, theta_cur)

      conv <- mcp_check_convergence(sys$R, sys$fb, tol = tol)
      if (isTRUE(conv$converged)) { conv_loc <- TRUE; break }

      # ---- Solve J * delta = -R (with singular-Jacobian LM fallback) ----
      delta_vec <- if (method == "dense") {
        tryCatch(solve(as.matrix(sys$J), -sys$R), error = function(e) NULL)
      } else if (backend == "Rcpp") {
        tryCatch({
          J_sp <- sys$J
          i <- J_sp@i + 1L
          j <- rep(seq_len(ncol(J_sp)), diff(J_sp@p))
          x <- J_sp@x
          n_dim <- nrow(J_sp)
          n_elems <- length(i)
          if (n_elems != length(j) || n_elems != length(x))
            stop("mcp_solve_path: sparse matrix triplet extraction mismatch")
          if (n_dim < 1L) stop("mcp_solve_path: Jacobian has zero rows")
          if (length(sys$R) != n_dim)
            stop("mcp_solve_path: RHS length ", length(sys$R), " != n_dim ", n_dim)
          mcp_sparse_solve_cpp(i, j, x, -sys$R, n_dim)
        }, error = function(e) {
          if (sparse_fallback)
            tryCatch(as.numeric(Matrix::solve(sys$J, -sys$R)),
                     error = function(e2) NULL)
          else NULL
        })
      } else {
        tryCatch(as.numeric(Matrix::solve(sys$J, -sys$R)),
                 error = function(e) NULL)
      }

      # Singular / non-finite plain solve → Levenberg-Marquardt damped step.
      if ((is.null(delta_vec) || anyNA(delta_vec) || any(!is.finite(delta_vec)))) {
        if (sparse_fallback) {
          sol <- .pf_robust_solve(sys$J, sys$R, sparse_fallback = TRUE)
          delta_vec <- sol$delta
        } else if (is.null(delta_vec)) {
          stop("mcp_solve_path: singular Jacobian at iter ", iter,
               ". Sys.J: ", deparse(dim(sys$J)))
        }
      }

      if (is.null(delta_vec) || anyNA(delta_vec) || any(!is.finite(delta_vec))) {
        warning("mcp_solve_path: non-finite Newton step at iter ", iter,
                "; aborting.")
        break
      }

      # ---- Line search or direct step ----
      if (line_search) {
        ls <- .mcp_line_search(
          Y_loc, delta_vec, sys$R, sys$J, theta_cur, merit_fn,
          sigma = sigma_line, max_ls = 20L, n_endo = n_endo, T = T
        )

        if (is.null(ls$Y_new)) {
          # No tried step reduced the merit — do NOT take a blind half step
          # (the old behaviour blew the path up; M15). Stop this run; the
          # homotopy driver can retry with continuation.
          break
        }
        Y_loc <- ls$Y_new
      } else {
        alpha <- step_size
        for (t in seq_len(T)) {
          idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
          Y_loc[t, ] <- Y_loc[t, ] + alpha * delta_vec[idx_t]
        }
      }
    }

    list(Y = Y_loc, converged = conv_loc, n_iter = n_iter_loc,
         max_res = max_res_loc, fb = fb_loc, merit_history = merit_loc)
  }

  # ==========================================================================
  # Homotopy / continuation driver. Ramp the shock path from zero (steady
  # state) to its full value over `n_stage` stages, warm-starting each stage.
  # ==========================================================================
  run_homotopy <- function(n_stage) {
    n_stage <- max(1L, as.integer(n_stage))
    Y_cur   <- Y
    res     <- NULL
    for (st in seq_len(n_stage)) {
      frac   <- if (n_stage == 1L) 1 else st / n_stage
      eps_st <- frac * eps_mat
      res    <- .mcp_run_newton(Y_cur, eps_st)
      Y_cur  <- res$Y
      if (!res$converged && st < n_stage && !all(is.finite(res$Y))) break
    }
    res
  }

  auto_mode <- identical(homotopy, "auto")
  n_stage   <- if (auto_mode) 1L else max(1L, as.integer(homotopy))

  result <- run_homotopy(n_stage)

  if (auto_mode && !result$converged) {
    for (ns in c(4L, 16L, 64L)) {
      result <- run_homotopy(ns)
      if (result$converged) break
    }
  }

  Y             <- result$Y
  converged     <- result$converged
  n_iter        <- result$n_iter
  max_res_final <- result$max_res
  fb_final      <- result$fb
  merit_hist    <- result$merit_history

  # Extract active set
  active_set <- mcp_assemble_active_set(Y, mcp_specs, tol = 1e-6)

  # Build output
  y_ss_mat <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
  irf      <- Y - y_ss_mat
  colnames(irf) <- dyn$endo_names
  rownames(irf) <- paste0("t", seq_len(T))
  colnames(Y) <- dyn$endo_names
  rownames(Y) <- paste0("t", seq_len(T))

  list(
    Y             = Y,
    irf           = irf,
    converged     = converged,
    n_iter        = n_iter,
    max_res       = max_res_final,
    fb_final      = fb_final,
    active_set    = active_set,
    endo_names    = dyn$endo_names,
    merit_history = merit_hist
  )
}


# =============================================================================
# Steady-state MCP solver
# =============================================================================

#' Solve for steady state with MCP (bound) constraints
#'
#' Extends the standard steady-state solver to handle models with occasionally
#' binding constraints at the steady state.  Uses the same Fischer-Burmeister
#' semi-smooth Newton approach as \code{mcp_solve_path()}, but on the STATIC
#' model (no time dimension).
#'
#' For each MCP-constrained variable x with bound b:
#'   φ(x − b, f_static(x, params)) = 0
#' where f_static is the steady-state residual function and φ is the
#' Fischer-Burmeister function.
#'
#' This is useful when a constraint binds in the steady state (e.g., a
#' permanently binding ZLB or a capacity constraint that is always active).
#'
#' @param compiled  dynhr_compiled (from \code{\link{compile_model}})
#' @param params    Named numeric parameter vector
#' @param mcp_specs List of MCP specs from \code{\link{mcp_parse_tags}}.
#'   Empty list → standard steady-state solve (no bounds).
#' @param y0        Optional initial guess for steady state (named numeric).
#' @param max_iter  Maximum Newton iterations (default 50)
#' @param tol       Convergence tolerance on max|R| (default 1e-10)
#' @param verbose   Logical: print progress messages (default FALSE)
#' @return List with:
#'   \item{values}{Named numeric: steady-state values}
#'   \item{converged}{Logical}
#'   \item{n_iter}{Integer: Newton iterations used}
#'   \item{max_res}{Numeric: final max|residual|}
#' @export
mcp_solve_steady <- function(compiled, params, mcp_specs = list(),
                              y0 = NULL, max_iter = 50L, tol = 1e-10,
                              verbose = FALSE) {
  static <- compiled$static
  n_endo <- length(static$endo_names)
  # Exogenous names come from the dynamic model (static model may not store them)
  exo_names <- compiled$dynamic$exo_names %||% compiled$model$varexo_names %||% character(0)
  n_exo  <- length(exo_names)
  n_spec <- length(mcp_specs)

  # Steady-state residual and Jacobian
  # Signature: res_fn(y, x, params, y_ss)
  #   y      = named endogenous vector
  #   x      = named exogenous vector (all zero at SS)
  #   params = parameter vector
  #   y_ss   = steady state
  res_fn <- static$residuals_fn
  jac_fn <- static$jacobian_fn

  # Exogenous variables at steady state (all zero)
  x_ss <- setNames(rep(0, n_exo), exo_names)

  # Ensure params is a plain numeric vector
  if (is.list(params)) params <- unlist(params)

  # Initial guess for endogenous variables
  y <- if (!is.null(y0)) setNames(as.numeric(y0[static$endo_names]), static$endo_names)
       else setNames(rep(0, n_endo), static$endo_names)
  y[is.na(y)] <- 0

  converged <- FALSE
  n_iter <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter <- iter
    R <- res_fn(y, x_ss, params, y)
    J <- jac_fn(y, x_ss, params, y)

    # Apply FB modifications for MCP-constrained equations
    if (n_spec > 0L) {
      for (j in seq_len(n_spec)) {
        sp <- mcp_specs[[j]]
        ei <- sp$eq_idx
        vi <- sp$var_idx
        x_cur <- y[vi]

        a <- if (sp$op == ">") x_cur - sp$bound else sp$bound - x_cur
        b <- if (sp$op == ">") R[ei] else -R[ei]

        fb <- mcp_fb(a, b)
        deriv <- mcp_fb_deriv(a, b)

        # Replace equation residual with FB function
        R[ei] <- fb
        # Modify Jacobian row
        J[ei, ] <- deriv$db * J[ei, ]
        J[ei, vi] <- J[ei, vi] + deriv$da
      }
    }

    max_res <- max(abs(R))
    if (max_res < tol) {
      converged <- TRUE
      break
    }

    # Newton step
    delta <- solve(J, -R)
    if (anyNA(delta)) break

    # Line search
    alpha <- 1.0
    for (ls in seq_len(20L)) {
      y_try <- y + alpha * delta
      names(y_try) <- static$endo_names
      R_try <- res_fn(y_try, x_ss, params, y_try)

      # Apply FB to MCP rows for merit
      if (n_spec > 0L) {
        for (j in seq_len(n_spec)) {
          sp <- mcp_specs[[j]]
          ei <- sp$eq_idx
          vi <- sp$var_idx
          a <- if (sp$op == ">") y_try[vi] - sp$bound else sp$bound - y_try[vi]
          b <- if (sp$op == ">") R_try[ei] else -R_try[ei]
          R_try[ei] <- mcp_fb(a, b)
        }
      }

      if (max(abs(R_try)) <= (1 - alpha * 0.5) * max_res || max(abs(R_try)) < tol)
        break
      alpha <- alpha * 0.5
    }

    y <- y + alpha * delta
  }

  # Build named output
  values <- rep(0, n_endo)
  names(values) <- static$endo_names
  values[static$endo_names] <- y

  list(
    values    = values,
    converged = converged,
    n_iter    = n_iter,
    max_res   = if (converged) max_res else Inf
  )
}
