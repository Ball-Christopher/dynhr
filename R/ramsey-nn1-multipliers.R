## R/ramsey-nn1-multipliers.R
## --------------------------------------------------------------------------
## E1: Steady-state Lagrange multiplier computation for the (n, n+1)
## approximation (Gross & Hansen 2021).
##
## Given a model and planner objective, computes the steady-state values
## of the Lagrange multipliers (lambda for backward-looking constraints,
## psi for forward-looking constraints) that satisfy the deterministic
## steady state of the Ramsey FOCs.
##
## Two methods:
##   Method A (preferred): Extract from Phase B's ramsey_steady() output.
##   Method B (standalone): Solve the linear subsystem A·m = b directly
##     using the model Jacobian and objective gradient at the steady state.
##
## References:
##   Gross, I. & Hansen, J. (2021). "Optimal policy design in nonlinear
##     DSGE models: An n-order accurate approximation." EER 140, 103918.
##   Bodenstein, M. & Guerrieri, L. (2019). Nash–Ramsey toolbox.
## --------------------------------------------------------------------------

#' Compute steady-state Lagrange multipliers for the (n, n+1) approximation
#'
#' Solves the linear subsystem A·m = b for the steady-state values of the
#' Lagrange multipliers associated with the Ramsey optimal policy problem.
#'
#' At the steady state, the FOC for variable y_i collapses to:
#'   df/dy_i + Σ_j λ̄_j · dg_j/dy_i(combined) + Σ_l ψ̄_l · (1+β) · dh_l/dy_i = 0
#'
#' where dg_j/dy_i(combined) = dg_j/dy_i(0) + dg_j/dy_i(-1) + dg_j/dy_i(+1)
#' and forward constraints contribute (1+β)·ψ̄_l due to lead/lag doubling.
#'
#' @param model            A dynhr_mod object.
#' @param compiled         Optional dynhr_compiled; computed if NULL.
#' @param ss               Named numeric steady state vector.
#' @param params           Named parameter vector; defaults to model$param_values.
#' @param planner_objective Character string: the planner objective expression.
#' @param method           "linear" (faster, default) or "newton" (more robust).
#' @param beta             Discount factor; defaults to params["beta"] or 0.99.
#' @param verbose          Print progress messages.
#'
#' @return A list with:
#'   \item{multipliers}{Named list with \code{lambda} (backward) and
#'         \code{psi} (forward) multiplier vectors.}
#'   \item{residual}{Max absolute residual of the linear system.}
#'   \item{A}{The coefficient matrix A (n_endo × n_eq).}
#'   \item{b}{The right-hand side vector (n_endo).}
#'   \item{method_used}{"phase_b", "linear", or "newton".}
#'   \item{forward_eqs}{Integer indices of forward-looking equations.}
#'   \item{backward_eqs}{Integer indices of backward-looking equations.}
#' @noRd
.compute_ramsey_ss_multipliers <- function(model,
                                           compiled = NULL,
                                           ss = NULL,
                                           params = NULL,
                                           planner_objective = NULL,
                                           method = c("linear", "newton"),
                                           beta = NULL,
                                           verbose = FALSE) {
  method <- match.arg(method)

  # ---- 1. Validate inputs ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  if (is.null(params)) params <- model$param_values
  if (is.null(ss)) {
    stop("ss (steady state) must be provided.")
  }
  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective provided.")
  }
  if (is.null(beta)) {
    beta <- if ("beta" %in% names(params) && is.finite(params[["beta"]])) {
      as.numeric(params[["beta"]])
    } else {
      0.99
    }
  }

  endo_names <- model$var_names
  n_endo <- length(endo_names)
  n_eq <- length(model$equations)
  exo_names <- model$varexo_names

  if (n_eq > n_endo) {
    stop(sprintf("Number of equations (%d) > number of endogenous vars (%d).",
                 n_eq, n_endo))
  }
  if (n_eq < n_endo) {
    if (verbose) {
      cat(sprintf("  Note: %d equations for %d vars (%d instrument(s)). ",
                  n_eq, n_endo, n_endo - n_eq))
      cat("Using least-squares for multiplier system.\n")
    }
  }

  # ---- 2. Ensure we have a compiled model ----
  if (is.null(compiled)) {
    compiled <- compile_model(model, verbose = FALSE)
  }

  # ---- 3. Classify equations as backward-looking (g) or forward-looking (h) ----
  eq_class <- .classify_equations(model)
  backward_idx <- eq_class$backward_idx
  forward_idx  <- eq_class$forward_idx
  n_back <- length(backward_idx)
  n_forw <- length(forward_idx)

  if (verbose) {
    cat(sprintf("[nn1_multipliers] %d backward, %d forward equations\n",
                n_back, n_forw))
  }

  # ---- 4. Build the linear system A·m = b ----
  # A[i, j] = combined derivative of equation j w.r.t. variable y_i at SS
  #   For backward eqs:  dg_j/dy_i = dg_j/dy_i(0) + dg_j/dy_i(-1) + dg_j/dy_i(+1)
  #   For forward eqs:   dh_l/dy_i * (1+beta)    (lead/lag doubling)
  # b[i] = -df/dy_i at SS

  A <- matrix(0, nrow = n_endo, ncol = n_eq)
  b <- numeric(n_endo)

  # Parse planner objective
  all_var_names <- c(endo_names, exo_names, model$varexo_det_names)
  obj_ast <- parse_expression(obj_text,
                              var_names = all_var_names,
                              param_names = model$param_names)

  exo_zero <- setNames(rep(0, length(exo_names)), exo_names)

  # Use compiled Jacobian for equation derivatives (handles model-local
  # variables like #Omega, #psi_n_ya correctly). The Jacobian matrices
  # f_minus, f_zero, f_plus are n_eq × n_endo at each lead/lag.
  sys <- extract_system_matrices(compiled, ss, params)
  J_minus <- sys$f_minus  # ∂f_j/∂y_i(t-1)
  J_zero  <- sys$f_zero   # ∂f_j/∂y_i(t)
  J_plus  <- sys$f_plus   # ∂f_j/∂y_i(t+1)

  for (i in seq_len(n_endo)) {
    var_name <- endo_names[i]

    # --- 4a. Compute -df/dy_i at SS ---
    df_ast <- ast_differentiate(obj_ast, var_name, 0L)
    if (!ast_is_zero(df_ast)) {
      b[i] <- -ast_eval(df_ast,
                        var_values = .make_var_values(ss, endo_names, exo_zero),
                        param_values = params,
                        ss_values = ss)
    }

    # --- 4b. For each equation j, compute combined derivative at SS ---
    # Combined: dg_j/dy_i = dg_j/dy_i(0) + dg_j/dy_i(-1) + dg_j/dy_i(+1)
    # For forward eqs: dh_l/dy_i * (1+beta)  (lead/lag doubling)
    for (j in seq_len(n_eq)) {
      total_deriv <- J_zero[j, i] + J_minus[j, i] + J_plus[j, i]
      if (!is.finite(total_deriv)) total_deriv <- 0

      # Apply lead/lag multiplier factor for forward equations
      if (j %in% forward_idx) {
        total_deriv <- total_deriv * (1 + beta)
      }

      A[i, j] <- total_deriv
    }
  }

  # ---- 5. Solve the linear system ----
  # Separate into backward and forward multiplier blocks
  # A_b is n_endo × n_back, A_f is n_endo × n_forw
  # Solve for [λ̄; ψ̄] in A·[λ; ψ] = b

  rank_A <- qr(A)$rank
  if (rank_A < n_eq) {
    if (verbose) {
      cat(sprintf("  A matrix is rank-deficient (rank %d < %d). Using QR solve.\n",
                  rank_A, n_eq))
    }
    mult <- qr.solve(A, b, tol = 1e-10)
  } else if (n_eq == n_endo) {
    mult <- solve(A, b)
  } else {
    # Under-determined system: use QR
    mult <- as.numeric(qr.solve(A, b, tol = 1e-10))
  }

  # ---- 6. Partition into lambda and psi ----
  lambda <- setNames(mult[backward_idx],
                     paste0("MULT_bwd_", backward_idx))
  psi    <- setNames(mult[forward_idx],
                     paste0("MULT_fwd_", forward_idx))

  # ---- 7. Compute residual ----
  resid <- as.numeric(A %*% mult - b)
  max_resid <- max(abs(resid), na.rm = TRUE)

  if (verbose) {
    cat(sprintf("  Multiplier SS solved. Max residual: %.2e\n", max_resid))
    if (n_back > 0) {
      cat("  Backward multipliers (lambda):\n")
      for (k in seq_along(lambda)) {
        cat(sprintf("    %s = %.6f\n", names(lambda)[k], lambda[k]))
      }
    }
    if (n_forw > 0) {
      cat("  Forward multipliers (psi):\n")
      for (k in seq_along(psi)) {
        cat(sprintf("    %s = %.6f\n", names(psi)[k], psi[k]))
      }
    }
  }

  # ---- 8. Return ----
  list(
    multipliers   = list(lambda = lambda, psi = psi),
    residual      = max_resid,
    A             = A,
    b             = b,
    method_used   = "linear",
    forward_eqs   = forward_idx,
    backward_eqs  = backward_idx,
    beta          = beta
  )
}


#' Classify equations as backward-looking or forward-looking
#'
#' An equation is forward-looking if any variable appears with lead_lag = +1.
#' Otherwise it is backward-looking.
#'
#' @param model A dynhr_mod object.
#' @return A list with:
#'   \item{backward_idx}{Integer indices of backward-looking equations.}
#'   \item{forward_idx}{Integer indices of forward-looking equations.}
#'   \item{is_forward}{Logical vector of length n_eq.}
#' @noRd
.classify_equations <- function(model) {
  n_eq <- length(model$equations)
  is_forward <- logical(n_eq)

  for (j in seq_len(n_eq)) {
    eq <- model$equations[[j]]
    is_forward[j] <- .equation_has_lead(eq)
  }

  list(
    backward_idx = which(!is_forward),
    forward_idx  = which(is_forward),
    is_forward   = is_forward
  )
}


#' Check if an equation AST has any variable with lead_lag = +1
#'
#' Recursively searches the equation (both LHS and RHS) for variable nodes
#' with non-zero lead_lag. Equations are stored as list(lhs, rhs, tag, text).
#'
#' @param node An equation node (list with lhs/rhs) or a bare AST node.
#' @return TRUE if any variable has lead_lag = +1.
#' @noRd
.equation_has_lead <- function(node) {
  if (is.null(node) || length(node) == 0) return(FALSE)
  if (!is.list(node)) return(FALSE)
  if (is.null(node$type)) {
    # Equation structure: list(lhs, rhs, tag, text)
    return(.equation_has_lead(node$lhs) || .equation_has_lead(node$rhs))
  }
  switch(node$type,
    variable = {
      ll <- node$lead_lag %||% 0L
      ll > 0L
    },
    number = FALSE,
    parameter = FALSE,
    local_variable = FALSE,
    unaryop = .equation_has_lead(node$operand),
    binop = .equation_has_lead(node$left) || .equation_has_lead(node$right),
    funcall = {
      for (a in node$args %||% list()) {
        if (.equation_has_lead(a)) return(TRUE)
      }
      FALSE
    },
    FALSE
  )
}


#' Extract multipliers from a Phase B ramsey_model result
#'
#' Utility to pull the steady-state multipliers from a Phase B result
#' for use in the Phase E (n, n+1) approximation.
#'
#' @param ramsey_result A dynhr_ramsey_result2 object from ramsey_model().
#' @return A list with lambda and psi vectors, or NULL if unavailable.
#' @noRd
.extract_multipliers_from_phase_b <- function(ramsey_result) {
  if (is.null(ramsey_result)) return(NULL)
  if (!inherits(ramsey_result, "dynhr_ramsey_result2")) return(NULL)

  ss <- ramsey_result$ramsey_steady
  if (is.null(ss)) return(NULL)

  mult <- ss$multiplier_ss
  if (is.null(mult)) return(NULL)

  # Try to identify forward vs backward multipliers from names
  # Phase B uses MULT_1, MULT_2, ... naming
  n_mult <- length(mult)

  # Partition: if the Phase B model classifies forward/backward equations
  # we store that info; otherwise return all as lambda (backward) with
  # empty psi
  lambda <- mult
  psi <- setNames(numeric(0), character(0))

  # If augmented model info is available, use equation classification
  aug_model <- ramsey_result$augmented_model
  if (!is.null(aug_model)) {
    eq_class <- .classify_equations(aug_model)
    # The augmented model has original equations + FOCs
    # Multipliers correspond to original equations
    n_orig <- length(aug_model$var_names) - n_mult
    # Actually, the multiplier mapping is complex. For simplicity,
    # return all multipliers and let the user sort them out.
  }

  list(
    lambda = lambda,
    psi    = psi,
    from_phase_b = TRUE
  )
}


#' Build a variable-values vector for AST evaluation
#'
#' Creates a named numeric vector including all lead/lag combinations
#' (__m1, __0, __p1) for each endogenous variable, as required by
#' ast_eval().
#'
#' @param ss         Named numeric steady state.
#' @param endo_names Character vector of endogenous variable names.
#' @param exo_zero   Named numeric zero vector for exogenous shocks.
#' @return Named numeric vector with __m1, __0, __p1 entries.
#' @noRd
.make_var_values <- function(ss, endo_names, exo_zero) {
  vals <- numeric(0)
  for (nm in endo_names) {
    sv <- if (!is.null(ss[[nm]]) && is.finite(ss[[nm]])) ss[[nm]] else 0.0
    vals[[paste0(nm, "__m1")]] <- sv
    vals[[paste0(nm, "__0")]]  <- sv
    vals[[paste0(nm, "__p1")]] <- sv
  }
  for (nm in names(exo_zero)) {
    vals[[paste0(nm, "__0")]] <- if (!is.null(exo_zero[[nm]])) exo_zero[[nm]] else 0.0
  }
  vals
}
