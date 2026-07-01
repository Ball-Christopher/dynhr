## R/ramsey-steady.R
## --------------------------------------------------------------------------
## Augmented steady-state solver for Ramsey optimal policy.
##
## Given a competitive-equilibrium steady state and the augmented model
## (original variables + Lagrange multipliers), solves for the steady state
## of the full augmented system.
##
## Approach (Bodenstein & Guerrieri 2019):
##   1. The original variables' SS values are known from the competitive-
##      equilibrium steady state.
##   2. At the steady state, leads/lags collapse. The FOCs become a linear
##      system in the multiplier values:
##        A · λ = b
##      where A_{i,j} = ∂g_j/∂y_i(combined) evaluated at SS, and
##      b_i = -∂f/∂y_i evaluated at SS.
##   3. Solve the linear system for λ_ss.
##   4. If refinement is needed, run one Newton step on the full augmented
##      system.
##
## References:
##   Bodenstein, M. & Guerrieri, L. (2019). Nash–Ramsey toolbox.
##   Schmitt-Grohé, S. & Uribe, M. (2004). Optimal fiscal and monetary
##     policy under imperfect competition.
## --------------------------------------------------------------------------

#' Solve the steady state of the augmented Ramsey system
#'
#' Takes the competitive-equilibrium steady state and augments it with
#' Lagrange multiplier steady-state values by solving the linear subsystem
#' implied by the FOCs at the steady state.
#'
#' @param model             Original dynhr_mod object.
#' @param compiled          Compiled version of the ORIGINAL model.
#' @param aug_model         Augmented dynhr_mod object (from ramsey_parse_augmented).
#' @param aug_compiled      Compiled version of the AUGMENTED model.
#' @param orig_ss           Competitive-equilibrium steady state from the original model.
#' @param params            Named parameter vector.
#' @param planner_objective Optional objective expression. If NULL, uses the
#'   parsed \code{planner_objective(...)} from \code{model}.
#' @param prefix            Multiplier variable name prefix (default "MULT").
#' @param refine            If TRUE, run one Newton step for refinement.
#' @param verbose           Print progress messages.
#'
#' @return A list with:
#'   \item{values}{Named numeric vector: full augmented steady state.}
#'   \item{converged}{Logical.}
#'   \item{multiplier_ss}{Named numeric vector of multiplier SS values.}
#'   \item{orig_ss}{Named numeric vector of original variable SS values.}
#'   \item{max_residual}{Max absolute residual at the solution.}
#'   \item{method_used}{String describing the method used.}
#' @export
ramsey_steady <- function(model,
                           compiled = NULL,
                           aug_model = NULL,
                           aug_compiled = NULL,
                           orig_ss = NULL,
                           params = NULL,
                           planner_objective = NULL,
                           prefix = "MULT",
                           refine = TRUE,
                           verbose = FALSE) {
  # ---- 1. Validate inputs ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }

  if (is.null(params)) params <- model$param_values
  if (is.null(orig_ss)) {
    stop("orig_ss (competitive-equilibrium steady state) must be provided.")
  }

  endo_names <- model$var_names
  n_endo <- length(endo_names)
  n_eq <- length(model$equations)
  multiplier_names <- paste0(prefix, "_", seq_len(n_eq))

  # Ensure orig_ss has the right shape
  orig_ss_vec <- setNames(numeric(n_endo), endo_names)
  for (nm in endo_names) {
    orig_ss_vec[nm] <- if (!is.null(orig_ss[[nm]]) && is.finite(orig_ss[[nm]])) {
      as.numeric(orig_ss[[nm]])
    } else {
      0.0
    }
  }

  if (verbose) {
    cat(sprintf("[ramsey_steady] Solving for %d multipliers from FOC linear system\n",
                n_eq))
  }

  # ---- 2. Build the Jacobian matrices symbolically and evaluate at SS ----
  # At steady state, the FOC for variable y_i is:
  #   df/dy_i + Σ_j λ_j * [dg_j/dy_i(0) + dg_j/dy_i(-1) + dg_j/dy_i(+1)] = 0
  #
  # This is a linear system: A @ λ = b
  #   A[i, j] = dg_j/dy_i(combined) evaluated at SS
  #   b[i]    = -df/dy_i evaluated at SS

  # Parse planner objective
  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective found in model.")
  }
  all_var_names <- c(endo_names, model$varexo_names, model$varexo_det_names)
  obj_ast <- parse_expression(obj_text,
                              var_names = all_var_names,
                              param_names = model$param_names)

  # Build matrix A and vector b
  A <- matrix(0, nrow = n_endo, ncol = n_eq)
  b <- numeric(n_endo)

  # Evaluate df/dy_i at steady state
  exo_zero <- setNames(rep(0, length(model$varexo_names)), model$varexo_names)

  for (i in seq_len(n_endo)) {
    var_name <- endo_names[i]

    # Compute -df/dy_i at SS
    df_ast <- ast_differentiate(obj_ast, var_name, 0L)
    if (!ast_is_zero(df_ast)) {
      b[i] <- -ast_eval(df_ast,
                        var_values = .make_var_values(orig_ss_vec, endo_names, exo_zero),
                        param_values = params,
                        ss_values = orig_ss_vec)
    }

    # For each equation j, compute combined derivative at SS
    for (j in seq_len(n_eq)) {
      eq <- model$equations[[j]]
      # Substitute #local variables before evaluation
      eq_resid <- equation_to_residual(eq)
      eq_resid <- ast_substitute_locals(eq_resid, model$local_variables)

      total_deriv <- ast_number(0)
      for (ll in c(-1L, 0L, 1L)) {
        deriv_ast <- ast_differentiate(eq_resid, var_name, ll)
        if (!ast_is_zero(deriv_ast)) {
          val <- ast_eval(deriv_ast,
                          var_values = .make_var_values(orig_ss_vec, endo_names, exo_zero),
                          param_values = params,
                          ss_values = orig_ss_vec)
          total_deriv_val <- if (!is.null(total_deriv$value)) total_deriv$value else 0
          total_deriv <- ast_number(total_deriv_val + val)
        }
      }
      A[i, j] <- if (total_deriv$type == "number") total_deriv$value else 0.0
    }
  }

  # ---- 3. Solve the linear system A @ λ = b ----
  # When n_eq < n_endo (instruments), the system is overdetermined.
  # Use ginv for non-square, or qr.solve for rank-deficient.
  rank_A <- qr(A)$rank
  if (rank_A < min(dim(A)) || n_endo > n_eq) {
    if (verbose) {
      cat(sprintf("  Solving %d x %d system via ginv (rank %d)\n",
                  n_endo, n_eq, rank_A))
    }
    lambda_ss <- MASS::ginv(A) %*% b
  } else {
    lambda_ss <- solve(A, b)
  }

  # Name the multipliers
  names(lambda_ss) <- multiplier_names

  if (verbose) {
    cat(sprintf("  Multiplier SS values:\n"))
    for (k in seq_along(lambda_ss)) {
      cat(sprintf("    %s = %.6f\n", multiplier_names[k], lambda_ss[k]))
    }
  }

  # ---- 4. Assemble full augmented steady state ----
  aug_values <- c(orig_ss_vec, lambda_ss)

  # ---- 5. Compute residuals of the augmented system ----
  max_resid <- NA_real_
  converged <- TRUE

  # If we have the compiled augmented model, verify residuals
  if (!is.null(aug_compiled)) {
    exo_zero_aug <- setNames(rep(0, length(aug_model$varexo_names)),
                              aug_model$varexo_names)
    r <- aug_compiled$static$residuals_fn(aug_values, exo_zero_aug,
                                           params, aug_values)
    max_resid <- max(abs(r), na.rm = TRUE)
    converged <- max_resid < 1e-6
    if (verbose) {
      cat(sprintf("  Augmented SS max residual: %.3e\n", max_resid))
    }
  }

  # ---- 6. Optional Newton refinement ----
  if (refine && !is.null(aug_compiled) && !converged) {
    if (verbose) cat("  Running Newton refinement on augmented system...\n")
    # Use the linear-subsystem solution as initial guess for Newton
    refined <- .refine_augmented_ss(aug_compiled, aug_values, params,
                                     aug_model$var_names, aug_model$varexo_names,
                                     verbose = verbose)
    if (!is.null(refined)) {
      aug_values <- refined$values
      max_resid <- refined$max_residual
      converged <- refined$converged
      # Update multiplier values
      for (k in seq_along(multiplier_names)) {
        mn <- multiplier_names[k]
        if (mn %in% names(aug_values)) {
          lambda_ss[mn] <- aug_values[mn]
        }
      }
    }
  }

  # ---- 7. Return ----
  result <- list(
    values          = aug_values,
    converged       = converged,
    multiplier_ss   = lambda_ss,
    orig_ss         = orig_ss_vec,
    max_residual    = max_resid,
    method_used     = if (refine && !is.null(aug_compiled) && !converged) "linear+newton" else "linear"
  )
  class(result) <- "dynhr_ramsey_steady"
  result
}


#' Refine augmented steady state with Newton's method
#'
#' Starting from the linear-subsystem solution, run Newton iterations
#' on the full augmented system to refine the steady state.
#'
#' @param compiled   Compiled augmented model.
#' @param y0         Initial guess for augmented steady state.
#' @param params     Parameter vector.
#' @param endo_names Endogenous variable names.
#' @param exo_names  Exogenous variable names.
#' @param max_iter   Maximum Newton iterations.
#' @param tol        Convergence tolerance.
#' @param verbose    Print progress.
#' @return List with $values, $converged, $max_residual.
#' @noRd
.refine_augmented_ss <- function(compiled, y0, params,
                                  endo_names, exo_names,
                                  max_iter = 50L, tol = 1e-10,
                                  verbose = FALSE) {
  res_fn <- compiled$static$residuals_fn
  jac_fn <- compiled$static$jacobian_fn
  n <- length(endo_names)
  x <- setNames(rep(0, length(exo_names)), exo_names)
  y <- y0[endo_names]
  y[is.na(y) | !is.finite(y)] <- 0.0

  converged <- FALSE
  iter <- 0L

  for (k in seq_len(max_iter)) {
    iter <- k
    r <- res_fn(y, x, params, y)
    if (any(!is.finite(r))) break

    max_r <- max(abs(r))
    if (verbose && (k <= 3 || k %% 10 == 0 || max_r < tol)) {
      cat(sprintf("    Newton iter %3d: max|r| = %.3e\n", k, max_r))
    }
    if (max_r < tol) {
      converged <- TRUE
      break
    }

    J <- jac_fn(y, x, params, y)
    if (any(!is.finite(J))) break

    dy <- solve(J, -r)
    if (any(!is.finite(dy))) break

    # Line search
    step <- 1.0
    y_new <- y + step * dy
    r_new <- res_fn(y_new, x, params, y_new)
    norm_r <- sum(r^2, na.rm = TRUE)
    for (ls in seq_len(10)) {
      norm_new <- if (all(is.finite(r_new))) sum(r_new^2, na.rm = TRUE) else Inf
      if (norm_new < norm_r) break
      step <- step * 0.5
      y_new <- y + step * dy
      r_new <- res_fn(y_new, x, params, y_new)
    }
    y <- y_new
  }

  r_final <- res_fn(y, x, params, y)

  list(
    values       = y,
    converged    = converged,
    iterations   = iter,
    max_residual = if (all(is.finite(r_final))) max(abs(r_final)) else Inf
  )
}


# .make_var_values is defined in ramsey-nn1-multipliers.R (the canonical version).
# The definition was removed here to prevent namespace collision.


#' @export
print.dynhr_ramsey_steady <- function(x, ...) {
  conv <- if (isTRUE(x$converged)) "CONVERGED" else "NOT CONVERGED"
  cat(sprintf("<dynhr_ramsey_steady>  [%s, method = %s]\n", conv, x$method_used))
  cat(sprintf("  Original vars:     %d\n", length(x$orig_ss)))
  cat(sprintf("  Multipliers:       %d\n", length(x$multiplier_ss)))
  cat(sprintf("  Max residual:      %.3e\n", x$max_residual))
  cat("  Multiplier values:\n")
  for (nm in names(x$multiplier_ss)) {
    cat(sprintf("    %s = %.6f\n", nm, x$multiplier_ss[nm]))
  }
  invisible(x)
}
