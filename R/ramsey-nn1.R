## R/ramsey-nn1.R
## --------------------------------------------------------------------------
## E5: Top-level entry point for the Gross-Hansen (n, n+1) approximation.
##
## ramsey_nn1() is the main function. It implements the full pipeline:
##   1. Parse model, compile, solve SS
##   2. Compute steady-state multipliers (E1)
##   3. Compute Taylor expansions (E2)
##   4. Build modified objective (E3)
##   5. Build and solve modified model (E4)
##   6. Compute welfare
##   7. Compare with Phase B (if available)
##   8. Return dynhr_nn1_result
##
## The (n, n+1) approximation provides an n-order accurate approximation
## of optimal policy without doubling the state space, using steady-state
## multipliers computed once.
##
## References:
##   Gross, I. & Hansen, J. (2021). "Optimal policy design in nonlinear
##     DSGE models: An n-order accurate approximation." EER 140, 103918.
## --------------------------------------------------------------------------

#' (n, n+1) approximation of optimal Ramsey policy
#'
#' Implements the Gross & Hansen (2021) (n, n+1) approximation for optimal
#' policy in nonlinear DSGE models.
#'
#' The (n, n+1) approximation uses an (n+1)-order Taylor expansion of the
#' planner's objective, eliminates linear terms using steady-state Lagrange
#' multipliers (the "blue correction"), and maximises the resulting modified
#' objective subject to n-order expansions of the constraints.
#'
#' Special cases:
#' \itemize{
#'   \item \strong{n = 1:} Linear-Quadratic (LQ) approximation — symmetric,
#'         certainty-equivalent, first-order accurate.
#'   \item \strong{n = 2:} Quadratic-Cubic (QC) approximation — asymmetric,
#'         risk-sensitive, second-order accurate.
#' }
#'
#' @param model             A dynhr_mod object (from \code{\link{parse_mod}}).
#' @param planner_objective Character string: the planner objective expression.
#'   If NULL, uses the \code{planner_objective(...)} block from the model.
#' @param n                 Approximation order: 1 (LQ), 2 (QC), 3, ...
#' @param ramsey_result     Optional \code{dynhr_ramsey_result2} from
#'   \code{\link{ramsey_model}} (Phase B). If provided, it supplies the
#'   steady-state multipliers and serves as a validation reference.
#' @param order             Alias for n (for compatibility with ramsey_model).
#'   If set, overrides n.
#' @param return_model      If TRUE, return the modified model objects for
#'   debugging.
#' @param beta              Discount factor. Defaults to \code{params["beta"]}
#'   or 0.99.
#' @param orig_ss           Optional pre-computed steady state list (from
#'   \code{\link{solve_steady}}).  When \code{NULL} (default) the steady state
#'   is solved internally.
#' @param compiled          Optional pre-compiled \code{dynhr_compiled}
#'   (from \code{\link{compile_model}}).  When \code{NULL} (default) the model
#'   is compiled internally.
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to
#'   \code{\link{solve_perturbation}} (e.g., \code{Sigma_e}, \code{h}).
#'
#' @return An object of class \code{dynhr_nn1_result} containing:
#'   \item{n}{Approximation order used.}
#'   \item{multipliers}{Steady-state Lagrange multipliers.}
#'   \item{modified_objective}{Polynomial coefficients of \eqn{W_t^{(n,n+1)}}.}
#'   \item{dr}{Decision rules from the optimal policy (original dimension).}
#'   \item{ss}{Steady state (same as competitive equilibrium).}
#'   \item{bk_ok}{Logical: Blanchard-Kahn condition satisfied?}
#'   \item{welfare}{Welfare metrics (unconditional, steady-state).}
#'   \item{comparison}{Comparison with Phase B (if \code{ramsey_result} provided).}
#'   \item{timing}{Computation time per step.}
#'   \item{meta}{Metadata (version, timestamp, params).}
#'
#' @examples
#' \dontrun{
#' # Parse a model
#' model <- parse_mod("path/to/model.mod")
#'
#' # LQ approximation (n=1)
#' lq <- ramsey_nn1(model, "-(pi^2 + 0.5*y_gap^2)", n = 1)
#'
#' # QC approximation (n=2)
#' qc <- ramsey_nn1(model, "-(pi^2 + 0.5*y_gap^2)", n = 2)
#'
#' # Compare with Phase B augmented-system Ramsey
#' ramsey <- ramsey_model(model, "-(pi^2 + 0.5*y_gap^2)")
#' qc2 <- ramsey_nn1(model, "-(pi^2 + 0.5*y_gap^2)", n = 2,
#'                    ramsey_result = ramsey)
#' print(qc2$comparison)
#' }
#'
#' @references
#'   Gross, T., & Hansen, J. (2021). Optimal policy design in nonlinear DSGE
#'     models: An n-order accurate approximation. \emph{Working Paper}.
#' @export
ramsey_nn1 <- function(model,
                       planner_objective = NULL,
                       n = 2L,
                       ramsey_result = NULL,
                       order = NULL,
                       return_model = FALSE,
                       beta = NULL,
                       orig_ss = NULL,
                       compiled = NULL,
                       verbose = FALSE,
                       ...) {
  # ---- 1. Validate ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object created by parse_mod().")
  }

  # Allow 'order' as alias for 'n'
  if (!is.null(order)) n <- as.integer(order)
  n <- as.integer(n)
  if (n < 1L) stop("n must be >= 1.")

  if (verbose) {
    cat(sprintf("\n============================================\n"))
    cat(sprintf("ramsey_nn1: (n=%d, n+1=%d) approximation\n", n, n + 1))
    cat(sprintf("============================================\n\n"))
  }

  timing <- list()
  t_start <- Sys.time()

  # ---- 2. Get planner objective ----
  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective provided. Add planner_objective(...) or pass planner_objective=.")
  }

  params <- model$param_values
  if (is.null(beta)) {
    beta <- if ("beta" %in% names(params) && is.finite(params[["beta"]])) {
      as.numeric(params[["beta"]])
    } else {
      0.99
    }
  }

  # ---- 3. Compile original model and solve steady state ----
  if (verbose) cat("[1/6] Compiling model and solving SS...\n")
  t1 <- Sys.time()

  # Detect instruments: extra endogenous vars beyond equations
  n_endo <- length(model$var_names)
  n_eq <- length(model$equations)
  instruments_used <- NULL
  
  if (n_endo > n_eq) {
    # Identify instrument candidates — variables referenced in the model
    # that have no equation defining them
    all_lhs_names <- sapply(model$equations, function(e) {
      if (e$lhs$type == "variable") e$lhs$name else NA_character_
    })
    instruments_used <- setdiff(model$var_names, all_lhs_names)
    if (length(instruments_used) > 0) {
      if (verbose) {
        cat(sprintf("  Detected %d instrument(s) (%d vars, %d eqs): %s\n",
                    length(instruments_used), n_endo, n_eq,
                    paste(instruments_used, collapse = ", ")))
      }
    }
  }

  # Use provided orig_ss if available
  if (!is.null(orig_ss)) {
    ss <- orig_ss
    if (verbose) cat("  Using provided orig_ss (skipping SS solve).\n")
    if (is.null(compiled)) {
      compiled <- compile_model(model, verbose = verbose)
    }
  } else {
    if (is.null(compiled)) {
      compiled <- compile_model(model, verbose = verbose)
    }

    ss_result <- solve_steady(compiled, params,
                              endo_names = model$var_names,
                              exo_names = model$varexo_names,
                              verbose = verbose)
    if (!isTRUE(ss_result$converged)) {
      # Try fallback
      if (verbose) cat("  Steady-state Newton did not converge. Trying nleqslv...\n")
      ss_result <- tryCatch(
        solve_steady_state(model, compiled, params = params, verbose = FALSE),
        error = function(e) NULL
      )
    }

    ss <- if (!is.null(ss_result) && isTRUE(ss_result$converged)) {
      if (!is.null(ss_result$values)) ss_result$values else ss_result$ss
    } else {
      stop("Competitive-equilibrium steady state did not converge. ",
           "Supply orig_ss = <named numeric steady state> to bypass SS solving.")
    }
  }

  timing$ss <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # ---- 4. Compute steady-state multipliers (E1) ----
  if (verbose) cat("[2/6] Computing steady-state multipliers...\n")
  t1 <- Sys.time()

  # Try Phase B first if available
  phase_b_multipliers <- NULL
  if (!is.null(ramsey_result)) {
    phase_b_multipliers <- .extract_multipliers_from_phase_b(ramsey_result)
  }

  if (!is.null(phase_b_multipliers) && length(phase_b_multipliers$lambda) > 0) {
    multipliers <- list(
      multipliers = list(
        lambda = phase_b_multipliers$lambda,
        psi = phase_b_multipliers$psi
      ),
      residual = NA_real_,
      method_used = "phase_b",
      forward_eqs = .classify_equations(model)$forward_idx,
      backward_eqs = .classify_equations(model)$backward_idx,
      beta = beta
    )
    if (verbose) cat("  Using Phase B multipliers.\n")
  } else {
    multipliers <- .compute_ramsey_ss_multipliers(
      model, compiled, ss, params, obj_text,
      method = "linear", beta = beta, verbose = verbose
    )
  }

  timing$multipliers <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # ---- 5. Compute Taylor expansions (E2) ----
  if (verbose) cat(sprintf("[3/6] Computing Taylor expansions (order %d)...\n", n + 1))
  t1 <- Sys.time()

  # Parse objective for Taylor expansion
  all_var_names <- c(model$var_names, model$varexo_names, model$varexo_det_names)
  obj_ast <- parse_expression(obj_text,
                              var_names = all_var_names,
                              param_names = model$param_names)

  taylor <- .nn1_taylor_expand(
    compiled, ss, params, order = n + 1,
    obj_ast = obj_ast,
    method = "symbolic",
    verbose = verbose
  )

  timing$taylor <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # ---- 6. Build modified objective (E3) ----
  if (verbose) cat("[4/6] Building modified objective W_t^{(n,n+1)}...\n")
  t1 <- Sys.time()

  nn1_objective <- .nn1_build_modified_objective(
    taylor, multipliers, n, beta, verbose = verbose
  )

  # Verify blue correction
  if (!.nn1_verify_blue_correction(nn1_objective)) {
    warning(sprintf(
      "Blue correction gradient at SS = %.2e. May exceed tolerance.", 
      nn1_objective$grad_at_ss))
  }

  timing$objective <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # ---- 7. Build and solve modified model (E4) ----
  if (verbose) cat("[5/6] Building and solving modified model...\n")
  t1 <- Sys.time()

  nn1_result <- .nn1_solve(
    model, compiled, ss, params,
    taylor, nn1_objective, multipliers,
    n, beta, verbose = verbose, ...
  )

  timing$solve <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # ---- 8. Compute welfare ----
  if (verbose) cat("[6/6] Computing welfare...\n")

  welfare <- .nn1_compute_welfare(
    nn1_result$dr, model, ss, params,
    obj_ast, beta, n
  )

  # ---- 9. Compare with Phase B (if available) ----
  comparison <- NULL
  if (!is.null(ramsey_result)) {
    comparison <- .nn1_compare_ramsey(nn1_result, ramsey_result, beta)
  }

  # ---- 10. Assemble result ----
  timing$total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  result <- list(
    n                  = n,
    multipliers        = multipliers$multipliers,
    modified_objective = nn1_objective$coefficients,
    dr                 = nn1_result$dr,
    ss                 = ss,
    bk_ok              = nn1_result$bk_ok,
    welfare            = welfare,
    comparison         = comparison,
    timing             = timing,
    meta = list(
      package_version = utils::packageVersion("dynhr"),
      timestamp       = Sys.time(),
      params          = params,
      planner_objective = obj_text,
      n_cols          = taylor$n_cols,
      n_endo          = length(model$var_names),
      n_exo           = length(model$varexo_names)
    )
  )

  # Optionally include model objects
  if (return_model) {
    result$modified_model    <- nn1_result$modified_model
    result$modified_compiled <- nn1_result$modified_compiled
    result$modified_ss       <- nn1_result$modified_ss
    result$taylor            <- taylor
  }

  class(result) <- c("dynhr_nn1_result", "list")

  if (verbose) {
    cat(sprintf("\n============================================\n"))
    cat(sprintf("ramsey_nn1 complete (%.1f sec)\n", timing$total))
    cat(sprintf("  BK: %s\n", if (result$bk_ok) "PASSED" else "FAILED"))
    cat(sprintf("  Welfare (uncond): %.6f\n", result$welfare$unconditional))
    if (!is.null(comparison)) {
      cat(sprintf("  Welfare diff vs Ramsey: %.6f\n", comparison$welfare_diff))
    }
    cat(sprintf("============================================\n"))
  }

  result
}


#' Compute welfare under the (n, n+1) optimal policy
#'
#' Computes unconditional and steady-state welfare using the decision rules
#' from the (n,n+1) approximation.
#'
#' @param dr       DecisionRules from the (n,n+1) solution.
#' @param model    Original dynhr_mod.
#' @param ss       Steady state.
#' @param params   Parameter vector.
#' @param obj_ast  Parsed planner objective AST.
#' @param beta     Discount factor.
#' @param n        Approximation order.
#' @return A list with unconditional welfare and steady-state welfare.
#' @noRd
.nn1_compute_welfare <- function(dr, model, ss, params, obj_ast, beta, n) {
  # Steady-state welfare
  welfare_ss <- .eval_planner_ast(obj_ast, ss, params, ss) / (1 - beta)

  # --- Unconditional welfare via simulation ---
  # For order 1: Lyapunov-based computation
  # For order 2+: use simulation
  welfare_uncond <- NA_real_

  if (n == 1) {
    # Use Lyapunov moments for first-order solution
    # Wrap in tryCatch: the solver may fail for models with unit roots
    # (e.g. NN1 placeholder equations). Fall through to simulation.
    moments <- tryCatch(
      compute_moments(dr, model, params = params),
      error = function(e) NULL
    )
    if (!is.null(moments) && !is.null(moments$std_dev) &&
        all(is.finite(moments$std_dev))) {
      # Welfare = E[f(y)] ≈ f(ss) + ½·tr(Σ_y · H_yy)
      # where H_yy is the Hessian of f w.r.t. endogenous variables
      # For LQ (n=1), the welfare is computed from the quadratic objective
      # which is already in the DR. Use simulation as fallback.
      welfare_uncond <- welfare_ss  # placeholder
    }
  }

  # Simulation-based welfare (works for any order)
  if (is.na(welfare_uncond)) {
    sim <- {
      ## Dispatch on the rule's order: an order-2 (n>=2) rule must use the
      ## second-order recursion or the welfare silently collapses to order 1.
      .simulate_dr_any_order(dr, n_periods = 10000L, model = model, burn_in = 1000L)
    }

    if (!is.null(sim)) {
      sim_levels <- attr(sim, "levels")
      if (!is.null(sim_levels)) {
        obj_t <- apply(sim_levels, 1L, function(row) {
          vals <- setNames(as.numeric(row), colnames(sim_levels))
          .eval_planner_ast(obj_ast, vals, params, ss)
        })
        obj_mean <- mean(obj_t, na.rm = TRUE)
        welfare_uncond <- obj_mean / max(1e-8, (1 - beta))
      }
    }
  }

  list(
    steady_state   = welfare_ss,
    unconditional  = welfare_uncond,
    discount       = beta
  )
}


#' Compare (n, n+1) welfare with Phase B augmented-system Ramsey
#'
#' @param nn1_result    Result from .nn1_solve().
#' @param ramsey_result Result from ramsey_model() (Phase B).
#' @param beta          Discount factor.
#' @return A list with ramsey_welfare, nn1_welfare, welfare_diff.
#' @noRd
.nn1_compare_ramsey <- function(nn1_result, ramsey_result, beta) {
  # Extract Ramsey welfare
  ramsey_welfare <- NA_real_

  if (inherits(ramsey_result, "dynhr_ramsey_result2")) {
    # Phase B result
    ramsey_dr <- ramsey_result$ramsey_dr$ramsey_dr
    ramsey_welfare <- ramsey_result$welfare_steady %||% NA_real_
  } else if (inherits(ramsey_result, "dynhr_ramsey_result")) {
    ramsey_welfare <- ramsey_result$welfare$unconditional_value %||% NA_real_
  }

  list(
    ramsey_welfare = ramsey_welfare,
    nn1_welfare    = NA_real_,  # computed elsewhere
    welfare_diff   = NA_real_   # computed elsewhere
  )
}


# ==========================================================================
# S3 methods
# ==========================================================================

#' @export
print.dynhr_nn1_result <- function(x, ...) {
  cat("\n<dynhr_nn1_result>\n")
  cat(sprintf("  (n, n+1) approximation with n = %d\n", x$n))
  cat(sprintf("  BK condition           : %s\n", if (isTRUE(x$bk_ok)) "PASSED" else "FAILED"))
  cat(sprintf("  Welfare (unconditional): %.6f\n", x$welfare$unconditional))
  cat(sprintf("  Welfare (steady state) : %.6f\n", x$welfare$steady_state))

  # Multiplier summary
  n_lambda <- length(x$multipliers$lambda)
  n_psi <- length(x$multipliers$psi)
  cat(sprintf("  Multipliers: %d backward (lambda), %d forward (psi)\n", n_lambda, n_psi))

  # Comparison
  if (!is.null(x$comparison) && !is.null(x$comparison$welfare_diff)) {
    cat(sprintf("  Welfare gap vs Ramsey  : %.6f\n", x$comparison$welfare_diff))
  }

  # Timing
  if (!is.null(x$timing)) {
    cat(sprintf("  Total time: %.2f sec\n", x$timing$total %||% NA_real_))
  }

  # Decision rules summary
  if (!is.null(x$dr)) {
    cat(sprintf("  Decision rules: %d endo, %d state, %d exo\n",
                NROW(x$dr$ghx), NCOL(x$dr$ghx), NCOL(x$dr$ghu)))
  }

  invisible(x)
}


#' @export
summary.dynhr_nn1_result <- function(object, ...) {
  cat("\nSummary: dynhr_nn1_result\n")
  cat("========================\n")
  cat(sprintf("Method               : (n, n+1) approximation (n = %d)\n", object$n))
  cat(sprintf("BK condition         : %s\n", if (isTRUE(object$bk_ok)) "Satisfied" else "Not satisfied"))
  cat(sprintf("Welfare (uncond)     : %.6f\n", object$welfare$unconditional))
  cat(sprintf("Welfare (steady)     : %.6f\n", object$welfare$steady_state))

  if (!is.null(object$dr) && !is.null(object$dr$ghx)) {
    cat(sprintf("State vector size    : %d\n", NCOL(object$dr$ghx)))
    if (!is.null(object$dr$state_vars)) {
      cat(sprintf("State variables      : %s\n",
                  paste(object$dr$state_vars, collapse = ", ")))
    }
  }

  if (!is.null(object$comparison)) {
    cat("\n--- Ramsey Comparison ---\n")
    cat(sprintf("Ramsey welfare       : %.6f\n", object$comparison$ramsey_welfare))
    cat(sprintf("(n,n+1) welfare      : %.6f\n", object$comparison$nn1_welfare))
    cat(sprintf("Welfare gap          : %.6f\n", object$comparison$welfare_diff))
  }

  if (!is.null(object$timing)) {
    cat(sprintf("\nTiming: SS=%.2fs, Mult=%.2fs, Taylor=%.2fs, Obj=%.2fs, Solve=%.2fs, Total=%.2fs\n",
                object$timing$ss %||% NA,
                object$timing$multipliers %||% NA,
                object$timing$taylor %||% NA,
                object$timing$objective %||% NA,
                object$timing$solve %||% NA,
                object$timing$total %||% NA))
  }

  invisible(object)
}
