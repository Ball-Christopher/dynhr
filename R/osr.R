## R/osr.R
## --------------------------------------------------------------------------
## Phase D — Optimal Simple Rules (OSR)
##
## Numerically minimise a quadratic loss over policy-rule parameters.
## For each candidate parameter vector:
##   1. Update model parameters with the candidate coefficients
##   2. Solve first-order perturbation
##   3. Compute the unconditional covariance via the Lyapunov equation
##   4. Evaluate the quadratic loss L = sum_ij W_ij Cov(y_i, y_j)
##      (Dynare's osr.objective: loss = W(:)' * vx(:))
##   5. CMA-ES drives the search, inside the box [lower, upper]
##
## Design decisions:
##   - Free parameters are existing model parameters (e.g. phi_pi, phi_y in
##     the NK Taylor rule).
##   - The loss is a full symmetric weight matrix W over the loss variables:
##     a named weight vector gives a diagonal W; a matrix or the .mod's
##     optim_weights block (with `y, pie w;` covariance entries) gives cross
##     terms.
##   - A .mod carrying osr_params / optim_weights / osr_params_bounds drives
##     osr() directly (model$osr), as Dynare's `osr;` command does.
##   - Solver failures (BK violations, non-convergence) and a loss that is not
##     finite (a loss variable with a unit root) -> large penalty.
## --------------------------------------------------------------------------

#' Optimal Simple Rules (OSR)
#'
#' Minimises a quadratic loss function over policy-rule parameters using
#' CMA-ES. For each candidate parameter vector, the model is solved at
#' first order, the unconditional covariance matrix \eqn{V} is computed via
#' the Lyapunov equation, and the loss
#' \deqn{L = \sum_{i,j} W_{ij} \, \mathrm{Cov}(y_i, y_j)}{L = sum_ij W_ij Cov(y_i, y_j)}
#' is evaluated (Dynare's OSR objective). With a diagonal \eqn{W} this is
#' \eqn{\sum_i w_i \mathrm{var}(y_i)}.
#'
#' @section Driving osr() from a .mod file:
#' When the model was parsed from a \code{.mod} file carrying Dynare's OSR
#' statements (stored in \code{model$osr} by \code{\link{parse_mod}}), the
#' arguments can be omitted and \code{osr()} runs the problem Dynare's
#' \code{osr;} command runs:
#' \itemize{
#'   \item \code{osr_params a b;} supplies \code{free_params} (started at their
#'     values in \code{params});
#'   \item \code{optim_weights} supplies \code{loss_vars} and the weight
#'     matrix. A variance entry \code{y w;} sets \eqn{W_{yy} = w}; a
#'     covariance entry \code{y, pie w;} adds \eqn{w \, \mathrm{Cov}(y, pie)}
#'     to the loss once, exactly as Dynare 7.1 does (it writes \code{w} into
#'     the single cell (y, pie) of \code{M_.osr.variable_weights}); the
#'     symmetric matrix used here carries \eqn{w/2} in each off-diagonal cell,
#'     which gives the same loss;
#'   \item \code{osr_params_bounds; a, lo, hi; end;} supplies the bounds.
#'     A free parameter without an entry is unbounded, as in Dynare.
#' }
#' Weights and bounds may be parameter expressions; they are evaluated at
#' \code{params}. Explicit arguments override the parsed specification.
#'
#' @section planner_objective without optim_weights (Dynare 7):
#' A \code{.mod} with a \code{planner_objective} and no \code{optim_weights}
#' block (and no explicit \code{loss_vars} / \code{loss_weights} /
#' \code{planner_objective} argument) runs Dynare 7's welfare-based OSR: the
#' loss is \eqn{-E[W]}, the unconditional welfare of the second-order (LQ)
#' approximation,
#' \deqn{E[W] = \frac{U + \tfrac12 \sum_{ij} U_{ij} \mathrm{Cov}(y_i, y_j)}{1 - \beta},}{E[W] = (U + 0.5 sum_ij U_ij Cov(y_i, y_j)) / (1 - beta),}
#' where \eqn{U} and \eqn{U_{ij}} are the objective and its Hessian at the
#' steady state (taken at the initial rule coefficients) and \eqn{\beta} is
#' \code{discount}, by default the \code{planner_discount} option of the
#' \code{.mod}'s \code{osr} command (required: Dynare's default of 1 makes the
#' objective infinite). \code{optimal_loss} is then Dynare's
#' \code{oo_.osr.objective_function}. This is exact for a linear model, and for
#' any model whose objective has a zero gradient at the steady state; for a
#' nonlinear model with a nonzero gradient a warning of class
#' \code{dynhr_warning_osr_lq_approx} says that the gradient terms of the
#' second-order welfare are omitted.
#'
#' @param model             A \code{dynhr_mod} object.
#' @param compiled          Optional pre-compiled model (will be recompiled
#'   each iteration if \code{recompile = TRUE}).
#' @param params            Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param free_params       Named numeric vector of initial values for the
#'   policy-rule coefficients to optimise (each name must be a model
#'   parameter), or a character vector of parameter names (started at their
#'   values in \code{params}). \code{NULL} (default) takes the \code{.mod}'s
#'   \code{osr_params}.
#' @param loss_vars         Character vector of endogenous variable names
#'   that enter the loss function. \code{NULL} (default) takes the names of
#'   \code{loss_weights}, or the variables of the \code{.mod}'s
#'   \code{optim_weights} block when \code{loss_weights} is also \code{NULL}.
#' @param loss_weights      Either a numeric vector of variance weights (named
#'   by variable, or positional along \code{loss_vars}; a diagonal weight
#'   matrix), or a square weight matrix whose row and column names are the
#'   loss variables (cross terms allowed; it is symmetrised, which leaves the
#'   loss unchanged). If both \code{loss_vars} and \code{loss_weights} are
#'   \code{NULL} the \code{.mod}'s \code{optim_weights} are used; if only
#'   \code{loss_weights} is \code{NULL}, all variance weights are 1.
#' @param lower,upper       Lower/upper bounds for the free parameters: a
#'   scalar (recycled), a vector along \code{free_params}, or a named vector
#'   (overriding only the named parameters). Parameters not bounded here take
#'   the \code{.mod}'s \code{osr_params_bounds}, and are otherwise unbounded
#'   (\code{-Inf}, \code{Inf}), as in Dynare.
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
#' @param discount          Planner discount factor for the
#'   \code{planner_objective}-only loss (see the section above); \code{NULL}
#'   (default) takes the \code{.mod}'s \code{osr(planner_discount = ...)}.
#'   Ignored otherwise.
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to
#'   \code{cmaes_optimize()}.
#'
#' @return An object of class \code{dynhr_osr_result} with components:
#'   \describe{
#'     \item{optimal_par}{Named numeric: optimal parameter values.}
#'     \item{optimal_loss}{Numeric: loss at the optimum.}
#'     \item{optimal_welfare}{Numeric: unconditional welfare at the optimum
#'       (\code{NA} for the variance loss).}
#'     \item{order}{Perturbation order used.}
#'     \item{objective}{\code{"welfare"} (simulated, \code{planner_objective}
#'       argument), \code{"planner_objective"} (the \code{.mod}'s objective, LQ
#'       welfare) or \code{"variance"}.}
#'     \item{discount}{Planner discount of the \code{"planner_objective"} loss
#'       (\code{NA} otherwise).}
#'     \item{free_params}{Character: names of optimised parameters.}
#'     \item{loss_vars}{Loss variables.}
#'     \item{loss_weights}{The symmetric weight matrix \eqn{W} over
#'       \code{loss_vars} (the loss is \code{sum(W * V)}).}
#'     \item{loss_offset}{Constant added to \code{sum(W * V)} (nonzero only
#'       for the \code{"planner_objective"} loss: \eqn{-U/(1-\beta)}).}
#'     \item{dr}{DecisionRules object at the optimum.}
#'     \item{moments}{Unconditional moments at the optimum (from
#'       \code{compute_moments()}).}
#'     \item{convergence}{CMA-ES convergence code.}
#'     \item{iterations}{Number of CMA-ES evaluations.}
#'
#'     \item{ramsey_comparison}{If \code{ramsey_result} provided, a list with
#'       Ramsey loss, OSR loss, and welfare gap.}
#'     \item{meta}{List with the bounds used (\code{lower}, \code{upper}),
#'       \code{max_iter}, the optimal \code{params} and a timestamp.}
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
                free_params = NULL,
                loss_vars = NULL,
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
                discount = NULL,
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
    .dynhr_warn("osr(): order = 2 without 'planner_objective' has no effect on the ",
            "variance loss (it is order-invariant). Supply 'planner_objective' ",
            "for a second-order welfare objective.", call. = FALSE)
  if (is.null(params)) params <- model$param_values
  spec <- model$osr

  # Free parameters: explicit, else the .mod's osr_params.
  if (is.null(free_params)) {
    if (length(spec$params) == 0L)
      .dynhr_abort("osr(): no 'free_params' supplied and the model has no ",
                   "osr_params statement.", class = "dynhr_error_osr_spec")
    free_params <- spec$params
  }
  if (is.character(free_params)) {
    missing <- setdiff(free_params, names(params))
    if (length(missing) > 0)
      stop("Free parameter(s) not found in model parameters: ",
           paste(missing, collapse = ", "))
    free_params <- params[free_params]
    if (any(!is.finite(free_params)))
      .dynhr_abort("osr(): the initial value of free parameter(s) ",
                   paste(names(free_params)[!is.finite(free_params)],
                         collapse = ", "), " is not finite.",
                   class = "dynhr_error_osr_spec")
  }
  if (is.null(names(free_params)) || any(!nzchar(names(free_params)))) {
    stop("'free_params' must be a named vector.")
  }

  # Check free params exist in model params
  missing <- setdiff(names(free_params), names(params))
  if (length(missing) > 0) {
    stop("Free parameter(s) not found in model parameters: ",
         paste(missing, collapse = ", "))
  }

  # Loss specification -> symmetric weight matrix W over loss_vars.
  has_parsed_w <- is.data.frame(spec$weights) && nrow(spec$weights) > 0L
  mod_objective <- model$planner_objective$text %||% ""
  use_lq_planner <- FALSE
  loss_offset <- 0
  if (is.null(loss_vars) && is.null(loss_weights)) {
    if (has_parsed_w) {
      loss_W <- .osr_parsed_weight_matrix(spec$weights, params)
      loss_vars <- rownames(loss_W)
    } else if (use_welfare) {
      loss_vars <- character(0)
      loss_W <- matrix(0, 0L, 0L)
    } else if (nzchar(mod_objective)) {
      ## Dynare 7: `osr` with a planner_objective and no optim_weights
      ## minimises -E[W]; the weight matrix is built below, once the model is
      ## compiled (it needs the steady state).
      use_lq_planner <- TRUE
    } else {
      .dynhr_abort("osr(): no 'loss_vars' / 'loss_weights' supplied and the ",
                   "model has neither an optim_weights block nor a ",
                   "planner_objective.",
                   class = "dynhr_error_osr_spec")
    }
  } else {
    if (is.null(loss_vars)) {
      loss_vars <- if (is.matrix(loss_weights)) rownames(loss_weights)
                   else names(loss_weights)
      if (is.null(loss_vars))
        stop("osr(): 'loss_vars' must be supplied when 'loss_weights' is ",
             "unnamed.")
    }
    loss_W <- .osr_weight_matrix(loss_vars, loss_weights)
  }
  # Check loss vars exist in model
  missing_vars <- setdiff(loss_vars, model$var_names)
  if (length(missing_vars) > 0) {
    stop("Loss variable(s) not found in model: ",
         paste(missing_vars, collapse = ", "))
  }

  # Bounds: explicit lower/upper, else the .mod's osr_params_bounds, else
  # unbounded (Dynare: M_.osr.param_bounds defaults to [-Inf, Inf]).
  n_free <- length(free_params)
  bnd <- .osr_resolve_bounds(names(free_params), lower, upper, spec$bounds,
                             params)
  lower <- bnd$lower
  upper <- bnd$upper
  # A calibrated start outside its bounds (e.g. a .mod whose calibration
  # violates osr_params_bounds) starts the search at the nearest bound.
  free_params <- pmin(pmax(free_params, lower), upper)

  # ---- 2. Compile model (once, if not recompiling each iteration) ----
  ## Order-2 OSR needs the second-order derivatives; ensure the cached compiled
  ## object carries them (one-time; the loss uses it when recompile = FALSE).
  max_ord <- if (order >= 2L) 2L else 1L
  if (is.null(compiled)) {
    compiled <- compile_model(model, max_order = max_ord, verbose = FALSE)
  } else if (order >= 2L) {
    compiled <- compile_model(model, max_order = 2L, verbose = FALSE)
  }

  # ---- 2b. planner_objective-only OSR (Dynare 7) -> quadratic loss ----
  if (use_lq_planner) {
    if (is.null(discount)) {
      pd <- .mod_command_options(model, "osr")$planner_discount
      if (is.null(pd))
        .dynhr_abort(
          "osr(): the loss comes from the model's planner_objective, which ",
          "needs the planner's discount factor: pass discount = <value> or ",
          "write osr(planner_discount = ...) in the .mod (Dynare's default, ",
          "1, makes the objective infinite).", class = "dynhr_error_osr_spec")
      discount <- .eval_mod_scalar_option(pd, params, "osr(planner_discount)")
    }
    start_params <- params
    start_params[names(free_params)] <- free_params
    lq <- .osr_planner_lq_loss(model, compiled, start_params, mod_objective,
                               discount)
    loss_W <- lq$W
    loss_vars <- rownames(loss_W)
    loss_offset <- lq$offset
  }

  # ---- 3. Build the loss function for CMA-ES ----
  loss_fn <- .make_osr_loss_fn(
    model         = model,
    compiled      = compiled,
    base_params   = params,
    free_names    = names(free_params),
    loss_W        = loss_W,
    loss_offset   = loss_offset,
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
    .dynhr_cat(sprintf("\n[osr] Optimising %d free parameter(s) via CMA-ES\n", n_free))
    .dynhr_cat(sprintf("      Loss vars: %s\n", paste(loss_vars, collapse = ", ")))
    .dynhr_cat(sprintf("      Init:      %s\n",
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
    .dynhr_warn("Could not solve model at optimal parameters.")
    opt_moments <- NULL
  } else {
    opt_moments <- compute_moments(opt_dr, model, params = opt_params)
    if (use_welfare) {
      wf <- tryCatch(
        welfare_compute(opt_dr, model, opt_params, planner_objective,
                        n_periods = welfare_n_periods, burn_in = welfare_burn_in,
                        seed = welfare_seed, verbose = FALSE),
        error = function(e) .dynhr_reraise_bug(e, NULL))
      opt_welfare <- if (!is.null(wf)) wf$unconditional %||% NA_real_ else NA_real_
    } else if (use_lq_planner) {
      ## The loss IS -E[W] (Dynare's objective_function).
      opt_welfare <- -opt_result$value
    }
  }

  # ---- 6. Ramsey comparison (optional) ----
  ramsey_comp <- NULL
  if (!is.null(ramsey_result)) {
    ramsey_comp <- .osr_compare_ramsey(
      osr_loss      = opt_result$value,
      osr_moments   = opt_moments,
      ramsey_result = ramsey_result,
      loss_W        = loss_W,
      loss_offset   = loss_offset
    )
  }

  # ---- 7. Assemble result ----
  result <- list(
    optimal_par   = opt_result$par,
    optimal_loss  = opt_result$value,
    optimal_welfare = opt_welfare,
    order         = order,
    objective     = if (use_welfare) "welfare"
                    else if (use_lq_planner) "planner_objective"
                    else "variance",
    discount      = if (use_lq_planner) discount else NA_real_,
    free_params   = names(free_params),
    loss_vars     = loss_vars,
    loss_weights  = loss_W,
    loss_offset   = loss_offset,
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
    .dynhr_cat(sprintf("\n[osr] Done. Optimal loss = %.6f\n", result$optimal_loss))
    .dynhr_cat(sprintf("      Optimal params: %s\n",
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
#'   4. Computes the unconditional covariance via Lyapunov
#'   5. Returns the quadratic loss sum(loss_W * V) (or a large penalty on
#'      failure)
#'
#' @noRd
.make_osr_loss_fn <- function(model, compiled, base_params,
                               free_names, loss_W,
                               recompile, order = 1L,
                               loss_offset = 0,
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
        error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(wf) || !is.finite(wf$unconditional)) return(1e20)
      return(-wf$unconditional)
    }

    # Quadratic loss sum(W * V) over the unconditional covariance (Lyapunov;
    # order-invariant). A loss variable with a unit root has a NaN (infinite)
    # variance: the loss is then not finite and the candidate is penalised,
    # never scored as if that variable were absent.
    moments <- compute_moments(dr, model, params = trial_params)
    if (is.null(moments)) return(1e20)

    loss <- .osr_quadratic_loss(moments$var_cov, loss_W) + loss_offset

    if (!is.finite(loss)) return(1e20)
    loss
  }
}


#' Quadratic OSR loss sum_ij W_ij V_ij (Dynare: W(:)' * vx(:))
#'
#' @param V Unconditional covariance matrix with variable dimnames.
#' @param W Symmetric weight matrix over the loss variables (dimnames).
#' @return Scalar; \code{NA} when a loss variable is missing from \code{V},
#'   \code{NaN} when a weighted (co)variance is not finite.
#' @noRd
.osr_quadratic_loss <- function(V, W) {
  if (length(W) == 0L) return(0)
  if (is.null(V) || !all(rownames(W) %in% rownames(V))) return(NA_real_)
  Vs <- V[rownames(W), colnames(W), drop = FALSE]
  nz <- W != 0
  sum(W[nz] * Vs[nz])
}


#' Quadratic OSR loss from the model's planner_objective (Dynare 7)
#'
#' Dynare 7's `osr` without optim_weights minimises -E[W], the unconditional
#' welfare of the second-order (LQ) approximation
#' (evaluate_planner_objective.m):
#'   E[W] = (U + 0.5 * sum_ij Uyy_ij Cov(y_i, y_j)) / (1 - beta),
#' with U and Uyy the planner objective and its Hessian at the steady state.
#' So the loss is sum(W * V) + offset with W = -Uyy / (2 (1 - beta)) and
#' offset = -U / (1 - beta).  This is exact for a linear model (Dynare's LQ
#' branch) and for any model whose objective has a zero gradient at the
#' steady state; otherwise Dynare's full second-order welfare adds
#' Uy * (second-order policy terms), which this loss omits (warned).
#'
#' @param params Parameters at which the steady state and Hessian are taken.
#' @return list(W, offset): W over the variables the objective involves.
#' @noRd
.osr_planner_lq_loss <- function(model, compiled, params, obj_text, discount) {
  if (!is.numeric(discount) || length(discount) != 1L || !is.finite(discount) ||
      discount <= 0 || discount >= 1)
    .dynhr_abort("osr(): the planner discount factor must lie in (0, 1) (got ",
                 paste(discount, collapse = " "), "); the unconditional ",
                 "welfare E[U]/(1 - beta) is not finite otherwise.",
                 class = "dynhr_error_osr_spec")
  ss_res <- solve_steady_state(model, compiled, params = params, verbose = FALSE)
  if (!isTRUE(ss_res$converged))
    .dynhr_abort("osr(): the steady state at the initial osr_params did not ",
                 "converge, so the planner_objective cannot be expanded ",
                 "around it.", class = "dynhr_error_osr_spec")
  ss   <- ss_res$ss
  endo <- model$var_names
  x0   <- as.numeric(ss[endo])
  f    <- .planner_objective_fn(obj_text, model, params, ss, endo,
                                model$varexo_names)
  U0   <- f(x0)
  Uyy  <- .planner_objective_hessian(obj_text, model, params, ss, endo,
                                     model$varexo_names, at = x0)
  if (!is.finite(U0) || any(!is.finite(Uyy)))
    .dynhr_abort("osr(): the planner_objective is not finite at the steady ",
                 "state.", class = "dynhr_error_osr_spec")
  keep <- which(rowSums(abs(Uyy)) > 0)
  if (length(keep) == 0L)
    .dynhr_abort("osr(): the planner_objective has a zero Hessian at the ",
                 "steady state, so its second-order welfare does not depend ",
                 "on the rule (every rule gives the same loss).",
                 class = "dynhr_error_osr_spec")
  if (!isTRUE(model$model_options$linear)) {
    h  <- 1e-6
    Uy <- vapply(seq_along(endo), function(i) {
      xp <- x0; xp[i] <- xp[i] + h
      xm <- x0; xm[i] <- xm[i] - h
      (f(xp) - f(xm)) / (2 * h)
    }, numeric(1))
    if (any(abs(Uy) > 1e-8 * max(1, abs(U0))))
      .dynhr_warn(
        "osr(): the planner_objective has a nonzero gradient at the steady ",
        "state of this nonlinear model; the quadratic loss omits the ",
        "gradient-times-second-order-policy terms of Dynare's order-2 welfare. ",
        "Use osr(planner_objective = ..., order = 2) for the simulated ",
        "second-order welfare.", class = "dynhr_warning_osr_lq_approx")
  }
  W <- -Uyy[keep, keep, drop = FALSE] / (2 * (1 - discount))
  list(W = W, offset = -U0 / (1 - discount))
}


#' Symmetric OSR weight matrix from a vector or matrix specification
#'
#' @param loss_vars Loss variables (row/column order of the result).
#' @param loss_weights NULL (unit variance weights), a named or positional
#'   numeric vector (diagonal), or a square matrix with the loss variables as
#'   row and column names (symmetrised; the loss sum(W * V) is unchanged).
#' @noRd
.osr_weight_matrix <- function(loss_vars, loss_weights) {
  n <- length(loss_vars)
  if (is.null(loss_weights)) {
    W <- diag(1, n)
  } else if (is.matrix(loss_weights)) {
    rn <- rownames(loss_weights); cn <- colnames(loss_weights)
    if (nrow(loss_weights) != ncol(loss_weights) || is.null(rn) ||
        is.null(cn) || !setequal(rn, cn) || anyDuplicated(rn) > 0L)
      stop("osr(): a 'loss_weights' matrix must be square, with the loss ",
           "variables as (identical sets of) row and column names.")
    if (!setequal(rn, loss_vars))
      stop("osr(): the names of the 'loss_weights' matrix (",
           paste(rn, collapse = ", "), ") must be the loss variables (",
           paste(loss_vars, collapse = ", "), ").")
    W <- loss_weights[loss_vars, loss_vars, drop = FALSE]
  } else if (is.null(names(loss_weights))) {
    # Positional weights
    if (length(loss_weights) != n)
      stop("Length of loss_weights must match length of loss_vars.")
    W <- diag(as.numeric(loss_weights), n)
  } else {
    # Named weights -- subset to loss_vars; unnamed loss vars get weight 0
    w <- as.numeric(loss_weights[loss_vars])
    w[is.na(w)] <- 0
    W <- diag(w, n)
  }
  storage.mode(W) <- "double"
  W <- (W + t(W)) / 2
  dimnames(W) <- list(loss_vars, loss_vars)
  if (any(!is.finite(W)))
    .dynhr_abort("osr(): every loss weight must be finite.",
                 class = "dynhr_error_osr_spec")
  W
}


#' Symmetric OSR weight matrix from a parsed optim_weights block
#'
#' Each entry adds its value to cell (var1, var2) of Dynare's (possibly
#' asymmetric) M_.osr.variable_weights; the symmetrised matrix gives the same
#' loss sum(W * V).
#'
#' @param wdf data.frame(var1, var2, expr) from parse_optim_weights_block().
#' @param params Named parameter vector the expressions are evaluated at.
#' @noRd
.osr_parsed_weight_matrix <- function(wdf, params) {
  vars <- unique(as.vector(rbind(wdf$var1, wdf$var2)))
  W <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  env <- .dynhr_param_eval_env(params)
  for (k in seq_len(nrow(wdf))) {
    v <- .dynhr_sandbox_eval(wdf$expr[k], env,
                             context = "the optim_weights expression")
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      .dynhr_abort("osr(): the optim_weights entry for ", wdf$var1[k],
                   if (wdf$var2[k] != wdf$var1[k]) paste0(", ", wdf$var2[k]),
                   " (`", wdf$expr[k], "`) does not evaluate to a finite ",
                   "number at the supplied parameters.",
                   class = "dynhr_error_osr_spec")
    W[wdf$var1[k], wdf$var2[k]] <- W[wdf$var1[k], wdf$var2[k]] + v
  }
  (W + t(W)) / 2
}


#' Resolve OSR parameter bounds
#'
#' Precedence per parameter: explicit lower/upper argument, then the parsed
#' osr_params_bounds entry, then unbounded (-Inf / Inf, Dynare's default).
#'
#' @param free_names Names of the free parameters.
#' @param lower,upper NULL, a scalar, a vector along free_names, or a named
#'   vector (overriding only the named parameters).
#' @param bounds_df data.frame(name, lower, upper) of expression text, or NULL.
#' @param params Named parameter vector the bound expressions are evaluated at.
#' @return list(lower, upper), named numeric vectors along free_names.
#' @noRd
.osr_resolve_bounds <- function(free_names, lower, upper, bounds_df, params) {
  n  <- length(free_names)
  lo <- stats::setNames(rep(-Inf, n), free_names)
  hi <- stats::setNames(rep(Inf, n), free_names)

  if (is.data.frame(bounds_df) && nrow(bounds_df) > 0L) {
    env <- .dynhr_param_eval_env(params)
    ev <- function(txt, nm) {
      ## Dynare spells infinity `inf` or `Inf`.
      txt <- gsub("(?<![A-Za-z0-9_])inf(?![A-Za-z0-9_])", "Inf", txt,
                  perl = TRUE)
      v <- .dynhr_sandbox_eval(txt, env,
                               context = "the osr_params_bounds expression")
      if (!is.numeric(v) || length(v) != 1L || is.na(v))
        .dynhr_abort("osr(): the osr_params_bounds entry for ", nm,
                     " (`", txt, "`) does not evaluate to a number at the ",
                     "supplied parameters.", class = "dynhr_error_osr_spec")
      v
    }
    for (k in seq_len(nrow(bounds_df))) {
      nm <- bounds_df$name[k]
      if (!(nm %in% free_names)) next
      lo[[nm]] <- ev(bounds_df$lower[k], nm)
      hi[[nm]] <- ev(bounds_df$upper[k], nm)
    }
  }

  overlay <- function(base, arg, what) {
    if (is.null(arg)) return(base)
    arg <- unlist(arg)
    if (!is.numeric(arg) || anyNA(arg))
      stop("osr(): '", what, "' must be numeric.")
    if (!is.null(names(arg)) && all(nzchar(names(arg)))) {
      bad <- setdiff(names(arg), free_names)
      if (length(bad) > 0L)
        stop("osr(): '", what, "' names parameter(s) that are not free: ",
             paste(bad, collapse = ", "))
      base[names(arg)] <- arg
    } else if (length(arg) == 1L) {
      base[] <- arg
    } else if (length(arg) == n) {
      base[] <- arg
    } else {
      stop("osr(): '", what, "' must have length 1 or ", n, ".")
    }
    base
  }
  lo <- overlay(lo, lower, "lower")
  hi <- overlay(hi, upper, "upper")
  if (any(lo > hi))
    .dynhr_abort("osr(): lower bound above upper bound for ",
                 paste(free_names[lo > hi], collapse = ", "), ".",
                 class = "dynhr_error_osr_spec")
  list(lower = lo, upper = hi)
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
    error = function(e) .dynhr_reraise_bug(e, NULL))
  if (is.null(dr2)) dr1 else dr2
}


#' Compare OSR loss against Ramsey-optimal commitment
#'
#' Computes the loss under the Ramsey policy using the same loss
#' function specification, then reports the welfare gap.
#'
#' @noRd
.osr_compare_ramsey <- function(osr_loss, osr_moments,
                                 ramsey_result, loss_W, loss_offset = 0) {
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

  ramsey_loss <- if (is.null(ramsey_moments)) NA_real_
                 else .osr_quadratic_loss(ramsey_moments$var_cov, loss_W) +
                        loss_offset

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
  W <- x$loss_weights
  if (is.matrix(W) && length(W) > 0L) {
    ## Diagonal entries, then each cross term as its total loss weight 2*W_ij
    ## (the coefficient on Cov(i, j) in sum(W * V)).
    lab <- sprintf("%s=%.4g", rownames(W), diag(W))
    ij <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
    if (nrow(ij) > 0L)
      lab <- c(lab, sprintf("Cov(%s,%s)=%.4g", rownames(W)[ij[, 1]],
                            colnames(W)[ij[, 2]], 2 * W[ij]))
    cat(sprintf("  Loss weights:     %s\n", paste(lab, collapse = ", ")))
  }
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
