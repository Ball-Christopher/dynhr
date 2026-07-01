## R/osr.R
## --------------------------------------------------------------------------
## Phase D — Optimal Simple Rules (OSR)
##
## Numerically minimise a quadratic loss over policy-rule parameters.
## For each candidate parameter vector:
##   1. Update model parameters with the candidate coefficients
##   2. Solve first-order perturbation
##   3. Compute unconditional variances via Lyapunov equation
##   4. Evaluate the quadratic loss L = Σ w_i · var(var_i)
##   5. CMA-ES drives the search
##
## Design decisions:
##   - Free parameters are existing model parameters (e.g. phi_pi, phi_y in
##     the NK Taylor rule). The user provides initial values and bounds.
##   - The loss is constructed from variable names + weights.
##   - Solver failures (BK violations, non-convergence) → large penalty.
## --------------------------------------------------------------------------

#' Optimal Simple Rules (OSR)
#'
#' Minimises a quadratic loss function over policy-rule parameters using
#' CMA-ES. For each candidate parameter vector, the model is solved at
#' first order, unconditional variances are computed via the Lyapunov
#' equation, and the loss \eqn{L = \sum_i w_i \cdot \mathrm{var}(y_i)} is
#' evaluated.
#'
#' @param model             A \code{dynhr_mod} object.
#' @param compiled          Optional pre-compiled model (will be recompiled
#'   each iteration if \code{recompile = TRUE}).
#' @param params            Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param free_params       Named numeric vector of initial values for the
#'   policy-rule coefficients to optimise. Each name must be a parameter
#'   name in the model.
#' @param loss_vars         Character vector of endogenous variable names
#'   that enter the loss function.
#' @param loss_weights      Named numeric vector of weights for each loss
#'   variable. Names must match entries in \code{loss_vars}. If \code{NULL},
#'   all weights are set to 1.
#' @param lower,upper       Lower/upper bounds for the free parameters.
#'   Recycled to length of \code{free_params}. Default \code{lower = 0.1},
#'   \code{upper = 10}.
#' @param recompile         If \code{TRUE} (default), recompiles the model
#'   at each candidate vector. Set to \code{FALSE} if the free parameters
#'   do not affect the Jacobian structure (only parameter values change).
#' @param max_iter          CMA-ES generation budget.
#' @param sigma0            Initial CMA-ES step size. \code{NULL} = auto.
#' @param ramsey_result     Optional \code{dynhr_ramsey_result2} for welfare
#'   comparison against Ramsey-optimal commitment.
#' @param order             Perturbation order at which the rule is solved
#'   (\code{1} or \code{2}). A variance loss is order-invariant, so \code{order
#'   = 2} is only meaningful together with \code{planner_objective}.
#' @param planner_objective Optional planner objective expression (string). When
#'   supplied, the loss is \code{-E[welfare]} (welfare from
#'   \code{\link{welfare_compute}}) instead of the weighted variance loss; at
#'   \code{order = 2} this captures the second-order volatility correction.
#' @param welfare_n_periods,welfare_burn_in Simulation length / burn-in for the
#'   welfare objective (used only when \code{planner_objective} is supplied).
#' @param welfare_seed      Fixed RNG seed making the simulated welfare a
#'   deterministic function of the policy coefficients (so CMA-ES optimises a
#'   reproducible objective); the caller's RNG state is restored.
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to
#'   \code{cmaes_optimize()}.
#'
#' @return An object of class \code{dynhr_osr_result} with components:
#'   \describe{
#'     \item{optimal_par}{Named numeric: optimal parameter values.}
#'     \item{optimal_loss}{Numeric: loss at the optimum.}
#'     \item{optimal_welfare}{Numeric: unconditional welfare at the optimum
#'       (\code{NA} unless \code{planner_objective} is supplied).}
#'     \item{order}{Perturbation order used.}
#'     \item{objective}{\code{"welfare"} or \code{"variance"}.}
#'     \item{free_params}{Character: names of optimised parameters.}
#'     \item{loss_vars, loss_weights}{Input specification.}
#'     \item{dr}{DecisionRules object at the optimum.}
#'     \item{moments}{Unconditional moments at the optimum (from
#'       \code{compute_moments()}).}
#'     \item{convergence}{CMA-ES convergence code.}
#'     \item{iterations}{Number of CMA-ES evaluations.}
#'
#'     \item{ramsey_comparison}{If \code{ramsey_result} provided, a list with
#'       Ramsey loss, OSR loss, and welfare gap.}
#'     \item{optimiser_trace}{Data frame of iteration history (if
#'       \code{verbose = TRUE}).}
#'   }
#'
#' @references
#'   Söderlind, P. (1999). Solution and estimation of RE macromodels with
#'     optimal policy. \emph{European Economic Review}, 43(4-6), 813-823.
#'   Hansen, N. (2016). The CMA evolution strategy: A tutorial.
#'     \emph{arXiv:1604.00772}.
#' @export
osr <- function(model,
                compiled = NULL,
                params = NULL,
                free_params,
                loss_vars,
                loss_weights = NULL,
                lower = NULL,
                upper = NULL,
                recompile = TRUE,
                max_iter = 10000L,
                sigma0 = NULL,
                ramsey_result = NULL,
                order = 1L,
                planner_objective = NULL,
                welfare_n_periods = 10000L,
                welfare_burn_in = 1000L,
                welfare_seed = 1L,
                verbose = FALSE,
                ...) {

  # ---- 1. Validate inputs ----
  if (!inherits(model, "dynhr_mod")) {
    stop("'model' must be a dynhr_mod object.")
  }
  order <- as.integer(order)
  if (!order %in% c(1L, 2L))
    stop("osr(): 'order' must be 1 or 2.")
  use_welfare <- !is.null(planner_objective)
  ## A variance loss L = sum w_i var(y_i) is invariant to the perturbation
  ## order (second-order terms shift means / add O(sigma^2) variance terms that
  ## the Lyapunov variance does not see), so order = 2 is only meaningful with a
  ## welfare objective (planner_objective), which captures the volatility
  ## correction. Warn rather than stop so an explicit order = 2 still solves.
  if (order >= 2L && !use_welfare)
    warning("osr(): order = 2 without 'planner_objective' has no effect on the ",
            "variance loss (it is order-invariant). Supply 'planner_objective' ",
            "for a second-order welfare objective.", call. = FALSE)
  if (is.null(params)) params <- model$param_values
  if (is.null(names(free_params)) || any(!nzchar(names(free_params)))) {
    stop("'free_params' must be a named vector.")
  }

  # Check free params exist in model params
  missing <- setdiff(names(free_params), names(params))
  if (length(missing) > 0) {
    stop("Free parameter(s) not found in model parameters: ",
         paste(missing, collapse = ", "))
  }

  # Check loss vars exist in model
  missing_vars <- setdiff(loss_vars, model$var_names)
  if (length(missing_vars) > 0) {
    stop("Loss variable(s) not found in model: ",
         paste(missing_vars, collapse = ", "))
  }

  # Set default weights
  if (is.null(loss_weights)) {
    loss_weights <- setNames(rep(1, length(loss_vars)), loss_vars)
  } else {
    if (is.null(names(loss_weights))) {
      # Positional weights
      if (length(loss_weights) != length(loss_vars)) {
        stop("Length of loss_weights must match length of loss_vars.")
      }
      loss_weights <- setNames(as.numeric(loss_weights), loss_vars)
    } else {
      # Named weights — subset to loss_vars
      loss_weights <- loss_weights[loss_vars]
      loss_weights[is.na(loss_weights)] <- 0
    }
  }

  # Default bounds
  n_free <- length(free_params)
  if (is.null(lower)) lower <- rep(0.1, n_free)
  if (is.null(upper)) upper <- rep(10, n_free)
  if (length(lower) == 1) lower <- rep(lower, n_free)
  if (length(upper) == 1) upper <- rep(upper, n_free)
  names(lower) <- names(free_params)
  names(upper) <- names(free_params)

  # ---- 2. Compile model (once, if not recompiling each iteration) ----
  ## Order-2 OSR needs the second-order derivatives; ensure the cached compiled
  ## object carries them (one-time; the loss uses it when recompile = FALSE).
  max_ord <- if (order >= 2L) 2L else 1L
  if (is.null(compiled)) {
    compiled <- compile_model(model, max_order = max_ord, verbose = FALSE)
  } else if (order >= 2L) {
    compiled <- compile_model(model, max_order = 2L, verbose = FALSE)
  }

  # ---- 3. Build the loss function for CMA-ES ----
  loss_fn <- .make_osr_loss_fn(
    model         = model,
    compiled      = compiled,
    base_params   = params,
    free_names    = names(free_params),
    loss_vars     = loss_vars,
    loss_weights  = loss_weights,
    recompile     = recompile,
    order         = order,
    planner_objective = planner_objective,
    welfare_n_periods = welfare_n_periods,
    welfare_burn_in   = welfare_burn_in,
    welfare_seed      = welfare_seed,
    verbose       = verbose
  )

  # ---- 4. Run CMA-ES optimisation ----
  if (verbose) {
    cat(sprintf("\n[osr] Optimising %d free parameter(s) via CMA-ES\n", n_free))
    cat(sprintf("      Loss vars: %s\n", paste(loss_vars, collapse = ", ")))
    cat(sprintf("      Init:      %s\n",
                paste(sprintf("%s = %.4g", names(free_params), free_params),
                      collapse = ", ")))
  }

  opt_result <- cmaes_optimize(
    fn       = loss_fn,
    par      = free_params,
    lower    = lower,
    upper    = upper,
    max_iter = max_iter,
    sigma0   = sigma0,
    verbose  = verbose,
    ...
  )

  # ---- 5. Compute the optimal model solution ----
  opt_params <- params
  opt_params[names(opt_result$par)] <- opt_result$par

  # Re-solve at the optimum
  opt_dr <- .osr_solve_model(model, compiled, opt_params, order = order,
                              recompile = recompile, verbose = verbose)

  opt_welfare <- NA_real_
  if (is.null(opt_dr)) {
    warning("Could not solve model at optimal parameters.")
    opt_moments <- NULL
  } else {
    opt_moments <- compute_moments(opt_dr, model, params = opt_params)
    if (use_welfare) {
      wf <- tryCatch(
        welfare_compute(opt_dr, model, opt_params, planner_objective,
                        n_periods = welfare_n_periods, burn_in = welfare_burn_in,
                        seed = welfare_seed, verbose = FALSE),
        error = function(e) NULL)
      opt_welfare <- if (!is.null(wf)) wf$unconditional %||% NA_real_ else NA_real_
    }
  }

  # ---- 6. Ramsey comparison (optional) ----
  ramsey_comp <- NULL
  if (!is.null(ramsey_result)) {
    ramsey_comp <- .osr_compare_ramsey(
      osr_loss      = opt_result$value,
      osr_moments   = opt_moments,
      ramsey_result = ramsey_result,
      loss_vars     = loss_vars,
      loss_weights  = loss_weights
    )
  }

  # ---- 7. Assemble result ----
  result <- list(
    optimal_par   = opt_result$par,
    optimal_loss  = opt_result$value,
    optimal_welfare = opt_welfare,
    order         = order,
    objective     = if (use_welfare) "welfare" else "variance",
    free_params   = names(free_params),
    loss_vars     = loss_vars,
    loss_weights  = loss_weights,
    dr            = opt_dr,
    moments       = opt_moments,
    convergence   = opt_result$convergence,
    iterations    = opt_result$iterations,
    ramsey_comparison = ramsey_comp,
    meta = list(
      lower    = lower,
      upper    = upper,
      max_iter = max_iter,
      params   = opt_params,
      timestamp = Sys.time()
    )
  )
  class(result) <- c("dynhr_osr_result", "list")

  if (verbose) {
    cat(sprintf("\n[osr] Done. Optimal loss = %.6f\n", result$optimal_loss))
    cat(sprintf("      Optimal params: %s\n",
                paste(sprintf("%s = %.4g", names(result$optimal_par),
                              result$optimal_par), collapse = ", ")))
  }

  result
}


# ==========================================================================
# Internal helpers
# ==========================================================================

#' Build the CMA-ES loss function closure
#'
#' Creates a function of the free-parameter vector that:
#'   1. Updates model parameters
#'   2. (Recompiles if needed)
#'   3. Solves first-order perturbation
#'   4. Computes unconditional variances via Lyapunov
#'   5. Returns the weighted loss (or a large penalty on failure)
#'
#' @noRd
.make_osr_loss_fn <- function(model, compiled, base_params,
                               free_names, loss_vars, loss_weights,
                               recompile, order = 1L,
                               planner_objective = NULL,
                               welfare_n_periods = 10000L,
                               welfare_burn_in = 1000L,
                               welfare_seed = 1L,
                               verbose) {

  use_welfare <- !is.null(planner_objective)

  function(x) {
    # Set names for the parameter vector
    names(x) <- free_names

    # Update params
    trial_params <- base_params
    trial_params[free_names] <- x

    # Solve the model (at the requested perturbation order)
    dr <- .osr_solve_model(model, compiled, trial_params, order = order,
                            recompile = recompile, verbose = FALSE)

    # If solution failed, return large penalty
    if (is.null(dr) || !isTRUE(dr$bk_satisfied)) {
      return(1e20)
    }

    if (use_welfare) {
      ## Welfare objective: minimise -E[welfare]. welfare_compute is made a
      ## deterministic function of x via a fixed simulation seed (RNG restored
      ## inside, so CMA-ES's own search randomness is untouched). At order 2 it
      ## simulates the pruned second-order rule, capturing the volatility
      ## correction the order-1 variance loss cannot see.
      wf <- tryCatch(
        welfare_compute(dr, model, trial_params, planner_objective,
                        n_periods = welfare_n_periods, burn_in = welfare_burn_in,
                        seed = welfare_seed, verbose = FALSE),
        error = function(e) NULL)
      if (is.null(wf) || !is.finite(wf$unconditional)) return(1e20)
      return(-wf$unconditional)
    }

    # Variance loss: weighted unconditional variances (Lyapunov; order-invariant)
    moments <- compute_moments(dr, model, params = trial_params)
    if (is.null(moments)) return(1e20)

    vars <- moments$std_dev^2
    loss <- sum(loss_weights * vars[loss_vars], na.rm = TRUE)

    if (!is.finite(loss)) return(1e20)
    loss
  }
}


#' Solve the model at given parameters
#'
#' Handles both recompile and cached-compiled modes.
#' Returns NULL on failure.
#'
#' @noRd
.osr_solve_model <- function(model, compiled, params, order = 1L,
                              recompile = TRUE, verbose = FALSE) {
  need_o2 <- order >= 2L
  if (recompile) {
    # Recompile: the Jacobian may change with new parameter values
    comp <- compile_model(model, max_order = if (need_o2) 2L else 1L,
                          verbose = FALSE)
  } else {
    comp <- compiled
  }

  # Steady state (linear models: ss = 0, handled automatically)
  ss_result <- solve_steady_state(model, comp, params = params,
                                   verbose = FALSE)

  if (!isTRUE(ss_result$converged)) return(NULL)

  # First-order perturbation.
  # BK violations are handled by the caller (loss function returns 1e20
  # penalty when !isTRUE(dr$bk_satisfied)), so no tryCatch needed here.
  dr1 <- solve_perturbation(model, comp, ss_result$ss, params, verbose = FALSE)
  if (!need_o2) return(dr1)
  if (is.null(dr1) || !isTRUE(dr1$bk_satisfied)) return(dr1)  # caller penalises

  # Lift to second order for the welfare objective.
  Sigma_e <- .get_shock_cov(model, model$varexo_names, params)
  dr2 <- tryCatch(
    solve_perturbation_order2(model, comp, ss_result$ss, params, dr1 = dr1,
                              Sigma_e = Sigma_e, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(dr2)) dr1 else dr2
}


#' Compare OSR loss against Ramsey-optimal commitment
#'
#' Computes the loss under the Ramsey policy using the same loss
#' function specification, then reports the welfare gap.
#'
#' @noRd
.osr_compare_ramsey <- function(osr_loss, osr_moments,
                                 ramsey_result, loss_vars, loss_weights) {
  # Extract Ramsey decision rules (restricted to original variables)
  ramsey_dr <- ramsey_result$ramsey_dr$ramsey_dr

  if (is.null(ramsey_dr) || is.null(ramsey_dr$ghx)) {
    return(list(
      ramsey_loss   = NA_real_,
      osr_loss      = osr_loss,
      welfare_gap   = NA_real_,
      message       = "Ramsey DR not available for comparison"
    ))
  }

  # Compute Ramsey unconditional moments
  ramsey_moments <- compute_moments(ramsey_dr, ramsey_result$augmented_model,
                                      params = ramsey_result$meta$params %||%
                                        ramsey_result$competitive_ss)

  if (is.null(ramsey_moments)) {
    ramsey_loss <- NA_real_
  } else {
    ramsey_vars <- ramsey_moments$std_dev^2
    ramsey_loss <- sum(loss_weights * ramsey_vars[loss_vars], na.rm = TRUE)
  }

  list(
    ramsey_loss = ramsey_loss,
    osr_loss    = osr_loss,
    welfare_gap = osr_loss - ramsey_loss  # positive = Ramsey better
  )
}


# ==========================================================================
# S3 methods
# ==========================================================================

#' @export
print.dynhr_osr_result <- function(x, ...) {
  cat(sprintf("\n<dynhr_osr_result>\n"))
  cat(sprintf("  Optimal loss:     %.6f\n", x$optimal_loss))

  cat("  Optimal params:\n")
  for (nm in names(x$optimal_par)) {
    cat(sprintf("    %-15s = %.6f\n", nm, x$optimal_par[nm]))
  }

  cat(sprintf("  Loss variables:   %s\n",
              paste(x$loss_vars, collapse = ", ")))
  cat(sprintf("  Loss weights:     %s\n",
              paste(sprintf("%s=%.2f", names(x$loss_weights), x$loss_weights),
                    collapse = ", ")))
  cat(sprintf("  CMA-ES evals:     %d\n", x$iterations))
  cat(sprintf("  Converged:        %s\n",
              if (isTRUE(x$convergence == 0)) "yes" else "no"))

  if (!is.null(x$ramsey_comparison)) {
    cat(sprintf("\n  Ramsey comparison:\n"))
    cat(sprintf("    Ramsey loss:    %.6f\n", x$ramsey_comparison$ramsey_loss))
    cat(sprintf("    OSR loss:       %.6f\n", x$ramsey_comparison$osr_loss))
    cat(sprintf("    Welfare gap:    %.6f (positive = Ramsey better)\n",
                x$ramsey_comparison$welfare_gap))
  }

  if (!is.null(x$dr)) {
    cat(sprintf("\n  DR:  %d x %d (ghx), %d x %d (ghu)\n",
                nrow(x$dr$ghx), ncol(x$dr$ghx),
                nrow(x$dr$ghu), ncol(x$dr$ghu)))
  }
  invisible(x)
}
