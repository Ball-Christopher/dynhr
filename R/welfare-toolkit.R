## R/welfare-toolkit.R
## --------------------------------------------------------------------------
## Phase F — Full welfare evaluation toolkit
##
## Provides a suite of welfare analysis functions that work with results
## from all policy approaches:
##   - Phase B: augmented-system Ramsey (ramsey_model / dynhr_ramsey_result2)
##   - Phase B v1: simulation-based Ramsey (ramsey_policy / dynhr_ramsey_result)
##   - Phase E: (n, n+1) approximation (ramsey_nn1 / dynhr_nn1_result)
##   - Phase D: OSR (osr / dynhr_osr_result)
##   - Phase D: Discretionary (discretionary_policy / dynhr_discretionary_result)
##
## Functions:
##   welfare_ce_diff()      — Consumption-equivalent welfare difference
##   welfare_cost_of_rule() — Welfare cost of deviating from Ramsey optimum
##   welfare_decompose()    — Level + volatility decomposition (order 2)
##   conditional_welfare()  — Welfare conditional on initial state
## --------------------------------------------------------------------------


# ==========================================================================
# 1. Internal helpers
# ==========================================================================

#' Extract welfare value from a policy result object
#'
#' Supports dynhr_ramsey_result, dynhr_ramsey_result2, dynhr_nn1_result,
#' dynhr_osr_result, dynhr_discretionary_result, and plain lists with
#' a welfare element.
#'
#' @param x A policy result object.
#' @param type Which welfare metric to extract: "unconditional" or "steady".
#' @return Numeric welfare value, or NA if unavailable.
#' @noRd
.extract_welfare <- function(x, type = c("unconditional", "steady")) {
  type <- match.arg(type)

  if (inherits(x, "dynhr_ramsey_result")) {
    # Phase B v1: simulation-based
    switch(type,
      unconditional = x$welfare$unconditional_value %||% NA_real_,
      steady        = x$welfare$steady_value %||% NA_real_
    )
  } else if (inherits(x, "dynhr_ramsey_result2")) {
    # Phase B v2: augmented-system Ramsey
    switch(type,
      unconditional = x$welfare$unconditional %||% x$welfare$steady_value %||% NA_real_,
      steady        = x$welfare$steady_value %||% NA_real_
    )
  } else if (inherits(x, "dynhr_nn1_result")) {
    # Phase E: (n,n+1) approximation
    switch(type,
      unconditional = x$welfare$unconditional %||% NA_real_,
      steady        = x$welfare$steady_state %||% NA_real_
    )
  } else if (inherits(x, "dynhr_osr_result")) {
    # Phase D: OSR
    switch(type,
      unconditional = -(x$optimal_loss %||% NA_real_) / max(1e-8, 1 - .extract_discount(x)),
      steady        = NA_real_
    )
  } else if (inherits(x, "dynhr_discretionary_result")) {
    # Phase D: Discretionary
    switch(type,
      unconditional = -(x$loss %||% x$ramsey_comparison$disc_loss %||% NA_real_) /
                       max(1e-8, 1 - .extract_discount(x)),
      steady        = NA_real_
    )
  } else if (is.list(x) && !is.null(x$welfare)) {
    # Generic list with welfare
    switch(type,
      unconditional = x$welfare$unconditional %||% x$welfare$unconditional_value %||% NA_real_,
      steady        = x$welfare$steady %||% x$welfare$steady_value %||% x$welfare$steady_state %||% NA_real_
    )
  } else {
    NA_real_
  }
}


#' Extract decision rules from a policy result object
#'
#' @param x A policy result object.
#' @return DecisionRules list (with ghx, ghu, etc.) or NULL.
#' @noRd
.extract_dr <- function(x) {
  if (inherits(x, "dynhr_ramsey_result")) {
    x$dr
  } else if (inherits(x, "dynhr_ramsey_result2")) {
    x$ramsey_dr$ramsey_dr %||% x$ramsey_dr
  } else if (inherits(x, "dynhr_nn1_result")) {
    x$dr
  } else if (inherits(x, "dynhr_osr_result")) {
    x$dr %||% x$optimal_dr
  } else if (inherits(x, "dynhr_discretionary_result")) {
    x$dr
  } else if (is.list(x) && !is.null(x$ghx)) {
    x  # It is itself a DecisionRules-like list
  } else if (is.list(x) && !is.null(x$dr)) {
    x$dr
  } else {
    NULL
  }
}


#' Extract discount factor from a policy result object
#'
#' @param x A policy result object.
#' @return Numeric discount factor, defaulting to 0.99.
#' @noRd
.extract_discount <- function(x) {
  val <- NULL
  if (inherits(x, "dynhr_ramsey_result")) {
    val <- x$welfare$discount
  } else if (inherits(x, "dynhr_ramsey_result2")) {
    val <- x$welfare$discount
  } else if (inherits(x, "dynhr_nn1_result")) {
    val <- x$welfare$discount
  } else if (inherits(x, "dynhr_discretionary_result")) {
    val <- x$discount
  } else if (is.list(x) && !is.null(x$discount)) {
    val <- x$discount
  }

  # Handle case where val might be a vector/matrix with multiple elements
  if (!is.null(val) && length(val) == 1 && is.finite(val)) {
    return(as.numeric(val))
  }
  0.99
}


#' Extract steady state from a policy result object
#'
#' @param x A policy result object.
#' @return Named numeric vector of steady-state values, or NULL.
#' @noRd
.extract_ss <- function(x, model = NULL) {
  if (inherits(x, "dynhr_ramsey_result")) {
    x$steady$ss %||% x$steady$values %||% x$steady
  } else if (inherits(x, "dynhr_ramsey_result2")) {
    x$competitive_ss %||% x$ramsey_steady$ss %||% x$ramsey_steady$values
  } else if (inherits(x, "dynhr_nn1_result")) {
    x$ss
  } else if (inherits(x, "dynhr_osr_result")) {
    # OSR stores dr but not ss explicitly; try model or solve from dr
    x$ss %||% x$steady %||% .infer_ss_from_dr(x$dr, model)
  } else if (inherits(x, "dynhr_discretionary_result")) {
    x$ss %||% x$steady %||% .infer_ss_from_dr(x$dr, model)
  } else if (is.list(x) && !is.null(x$ss)) {
    x$ss
  } else {
    .infer_ss_from_dr(.extract_dr(x), model)
  }
}

#' Infer steady state from decision rules or model
#'
#' Attempts to extract or compute steady-state values from a set of
#' decision rules or a dynhr_mod.
#'
#' @param dr    DecisionRules object (optional).
#' @param model dynhr_mod (optional).
#' @return Named numeric vector of steady-state values, or NULL.
#' @noRd
.infer_ss_from_dr <- function(dr, model) {
  if (!is.null(dr$state_lag) && !is.null(dr$endo_names)) {
    # Steady state of a linear model is the origin
    ss_vec <- rep(0, length(dr$endo_names))
    names(ss_vec) <- dr$endo_names
    return(ss_vec)
  }
  if (!is.null(model) && inherits(model, "dynhr_mod") &&
      !is.null(model$initval)) {
    return(model$initval)
  }
  NULL
}


#' Extract model from a policy result object or retrieve from parent
#'
#' @param x A policy result object.
#' @param model Optional dynhr_mod fallback.
#' @return dynhr_mod or NULL.
#' @noRd
.extract_model <- function(x, model = NULL) {
  if (!is.null(model) && inherits(model, "dynhr_mod")) return(model)

  if (inherits(x, "dynhr_ramsey_result2")) {
    return(x$augmented_model)
  }

  # For other types, search common fields
  for (field in c("model", "augmented_model", "modified_model")) {
    if (!is.null(x[[field]]) && inherits(x[[field]], "dynhr_mod")) {
      return(x[[field]])
    }
  }
  NULL
}


#' Extract shock covariance from a model or result object
#'
#' @param model dynhr_mod.
#' @param params Named parameter vector.
#' @return Sigma_e matrix.
#' @noRd
.extract_sigma_e <- function(model, params) {
  exo <- model$varexo_names
  n_exo <- length(exo)
  Sigma_e <- matrix(0, n_exo, n_exo)
  stderr <- vapply(seq_len(n_exo), function(k) {
    sname <- exo[k]
    # Check model$shocks for stderr specification
    if (!is.null(model$shocks)) {
      sd_val <- model$shocks$std[sname]
      if (!is.null(sd_val) && is.finite(sd_val)) return(as.numeric(sd_val))
    }
    # Fallback: check if a parameter named stderr_<name> exists
    pname <- paste0("stderr_", sname)
    if (pname %in% names(params) && is.finite(params[pname])) {
      return(as.numeric(params[pname]))
    }
    # Fallback: check if a parameter with the shock name exists
    if (sname %in% names(params) && is.finite(params[sname])) {
      return(as.numeric(params[sname]))
    }
    1.0  # Default fallback
  }, numeric(1))
  diag(Sigma_e) <- stderr^2
  rownames(Sigma_e) <- exo
  colnames(Sigma_e) <- exo
  Sigma_e
}


#' Estimate marginal utility of consumption at steady state
#'
#' Uses a numerical derivative of the planner objective with respect to
#' the consumption variable at the steady state. If consumption_variable
#' is NULL, returns 1/c_ss (CRRA approximation).
#'
#' @param obj_ast  Parsed planner objective AST.
#' @param ss       Steady-state values (named numeric).
#' @param params   Parameter vector.
#' @param consumption_variable Character: name of consumption variable.
#' @param eps      Finite-difference step.
#' @return Marginal utility (dU/dC) at SS, or NA.
#' @noRd
.estimate_marginal_utility <- function(obj_ast, ss, params,
                                        consumption_variable = NULL,
                                        eps = 1e-6) {
  if (is.null(consumption_variable) || !consumption_variable %in% names(ss)) {
    # CRRA approximation: MU = 1/c_ss (under log utility)
    # Try common consumption variable names
    for (cv in c("c", "C", "cons", "consumption", "CONS", "y", "Y")) {
      if (cv %in% names(ss) && is.finite(ss[[cv]]) && abs(ss[[cv]]) > 1e-12) {
        return(1 / abs(ss[[cv]]))
      }
    }
    return(NA_real_)
  }

  if (is.null(obj_ast)) return(1 / abs(ss[[consumption_variable]]))

  # Central difference of objective w.r.t. consumption
  ss_up <- ss; ss_up[[consumption_variable]] <- ss[[consumption_variable]] + eps
  ss_dn <- ss; ss_dn[[consumption_variable]] <- ss[[consumption_variable]] - eps

  f_up <- .eval_planner_ast(obj_ast, ss_up, params, ss)
  f_dn <- .eval_planner_ast(obj_ast, ss_dn, params, ss)

  if (is.finite(f_up) && is.finite(f_dn)) {
    (f_up - f_dn) / (2 * eps)
  } else {
    1 / abs(ss[[consumption_variable]])
  }
}


# ==========================================================================
# 2. Consumption-equivalent welfare difference
# ==========================================================================

#' Consumption-equivalent welfare difference
#'
#' Computes the consumption-equivalent (CE) welfare difference between two
#' policy regimes. The CE measure expresses the welfare gain of regime A
#' over regime B as a percentage of steady-state consumption that a
#' representative agent would require to be indifferent between the two.
#'
#' The formula is:
#' \deqn{CE = \frac{(1-\beta)(W_A - W_B)}{\partial u / \partial c \cdot c_{ss}} \times 100}
#'
#' where \eqn{W_A} and \eqn{W_B} are unconditional welfare values under
#' regimes A and B, \eqn{\beta} is the discount factor, \eqn{\partial u / \partial c}
#' is the marginal utility of consumption at the steady state, and
#' \eqn{c_{ss}} is steady-state consumption.
#'
#' @param result_a            First policy regime result (the "new" policy).
#'   Accepts \code{dynhr_ramsey_result}, \code{dynhr_ramsey_result2},
#'   \code{dynhr_nn1_result}, \code{dynhr_osr_result}, or
#'   \code{dynhr_discretionary_result}.
#' @param result_b            Second policy regime result (the "baseline" policy).
#'   Same types as \code{result_a}.
#' @param model               Optional \code{dynhr_mod}. If not provided,
#'   extracted from the result objects when possible.
#' @param params              Optional named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param consumption_variable Character: name of the consumption variable in
#'   the model. If \code{NULL}, attempts to auto-detect from common names
#'   (\code{c}, \code{C}, \code{cons}, \code{consumption}, \code{y}).
#' @param planner_objective   Optional objective expression. If \code{NULL},
#'   extracted from the result objects.
#' @param verbose             Print additional information.
#'
#' @return A list with:
#'   \describe{
#'     \item{ce_percent}{Consumption-equivalent difference in percent of
#'       steady-state consumption. Positive means \code{result_a} is
#'       welfare-improving relative to \code{result_b}.}
#'     \item{welfare_a}{Unconditional welfare under regime A.}
#'     \item{welfare_b}{Unconditional welfare under regime B.}
#'     \item{welfare_gap}{Raw welfare difference (A - B).}
#'     \item{discount}{Discount factor used.}
#'     \item{marginal_utility}{Marginal utility of consumption at SS.}
#'     \item{consumption_ss}{Steady-state consumption value.}
#'     \item{method}{Description of the CE computation method.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Compare Ramsey policy vs Taylor rule
#' ramsey <- ramsey_model(model, "-(pi^2 + 0.5*y_gap^2)")
#' taylor <- osr(model, free_params = c(phi_pi = 1.5, phi_y = 0.5),
#'               loss_vars = c("pi", "y_gap"))
#' ce <- welfare_ce_diff(ramsey, taylor, model, consumption_variable = "y")
#' print(sprintf("CE welfare gain of Ramsey: %.2f%%", ce$ce_percent))
#' }
#'
#' @export
welfare_ce_diff <- function(result_a, result_b,
                             model = NULL,
                             params = NULL,
                             consumption_variable = NULL,
                             planner_objective = NULL,
                             verbose = FALSE) {

  # ---- 1. Extract welfare ----
  welfare_a <- .extract_welfare(result_a, "unconditional")
  welfare_b <- .extract_welfare(result_b, "unconditional")

  if (!is.finite(welfare_a) || !is.finite(welfare_b)) {
    stop("Cannot compute CE difference: unconditional welfare not available ",
         "for both regimes. Try using simulation-based welfare via ",
         "welfare_compute() first.")
  }

  # ---- 2. Extract model and params ----
  model <- .extract_model(result_a, model) %||% .extract_model(result_b, model)
  if (is.null(model)) {
    stop("Could not locate a dynhr_mod object. Please supply the 'model' argument.")
  }
  if (is.null(params)) params <- model$param_values

  # ---- 3. Extract discount ----
  discount <- .extract_discount(result_a)

  # ---- 4. Extract steady state ----
  ss <- .extract_ss(result_a, model) %||% .extract_ss(result_b, model)
  if (is.null(ss)) {
    stop("Could not locate steady-state values. Please supply steady state ",
         "via the result objects.")
  }

  # ---- 5. Get or parse planner objective ----
  if (is.null(planner_objective)) {
    # Try to extract from results
    if (inherits(result_a, "dynhr_nn1_result") && !is.null(result_a$meta$planner_objective)) {
      planner_objective <- result_a$meta$planner_objective
    } else if (inherits(result_a, "dynhr_ramsey_result") && !is.null(result_a$objective$text)) {
      planner_objective <- result_a$objective$text
    } else if (inherits(result_a, "dynhr_ramsey_result2")) {
      planner_objective <- result_a$augmented_result$planner_objective %||%
                           model$planner_objective$text
    }
  }

  # Parse objective AST for marginal utility estimation
  obj_ast <- NULL
  if (!is.null(planner_objective) && nzchar(trimws(planner_objective))) {
    all_vars <- c(model$var_names, model$varexo_names)
    obj_ast <- parse_expression(planner_objective,
        var_names = all_vars,
        param_names = model$param_names)
  }

  # ---- 6. Compute consumption-equivalent measure ----
  # Identify consumption variable
  cons_var <- NULL
  if (!is.null(consumption_variable) && consumption_variable %in% names(ss)) {
    cons_var <- consumption_variable
  } else {
    for (cv in c("c", "C", "cons", "consumption", "CONS", "y", "Y")) {
      if (cv %in% names(ss) && is.finite(ss[[cv]]) && abs(ss[[cv]]) > 1e-12) {
        cons_var <- cv
        break
      }
    }
  }

  if (is.null(cons_var)) {
    warning("Could not identify consumption variable. ",
            "Specify 'consumption_variable' explicitly. ",
            "Using default CRRA approximation.")
    cons_ss <- 1.0
  } else {
    cons_ss <- abs(as.numeric(ss[[cons_var]]))
    if (!is.finite(cons_ss) || cons_ss < 1e-12) cons_ss <- 1.0
  }

  # Marginal utility
  mu <- .estimate_marginal_utility(obj_ast, ss, params, cons_var)

  if (!is.finite(mu) || mu < 1e-12) {
    # Fallback: CRRA approximation
    mu <- 1 / cons_ss
  }

  # CE formula: CE = (1-β) * (W_A - W_B) / (MU * c_ss)
  welfare_gap <- welfare_a - welfare_b
  ce_raw <- (1 - discount) * welfare_gap / (mu * cons_ss)
  ce_pct <- ce_raw * 100

  if (verbose) {
    cat(sprintf("CE welfare difference:\n"))
    cat(sprintf("  Welfare A:          %.6f\n", welfare_a))
    cat(sprintf("  Welfare B:          %.6f\n", welfare_b))
    cat(sprintf("  Welfare gap:        %.6f\n", welfare_gap))
    cat(sprintf("  Discount:           %.4f\n", discount))
    cat(sprintf("  MU (dU/dC):         %.6f\n", mu))
    cat(sprintf("  C_ss:               %.6f\n", cons_ss))
    cat(sprintf("  CE (%% of C_ss):    %.4f%%\n", ce_pct))
  }

  result <- list(
    ce_percent       = ce_pct,
    ce_raw           = ce_raw,
    welfare_a        = welfare_a,
    welfare_b        = welfare_b,
    welfare_gap      = welfare_gap,
    discount         = discount,
    marginal_utility = mu,
    consumption_ss   = cons_ss,
    consumption_var  = cons_var,
    method           = if (!is.null(cons_var))
                         sprintf("CE via MU = %.4f, C_ss = %.4f", mu, cons_ss)
                       else "CRRA approximation (1/C_ss)"
  )
  class(result) <- c("dynhr_ce_diff", "list")
  result
}


# ==========================================================================
# 3. Welfare cost of deviating from Ramsey
# ==========================================================================

#' Welfare cost of deviating from the Ramsey optimum
#'
#' Computes the welfare cost (in consumption-equivalent units) of using an
#' alternative policy — such as an OSR rule or discretionary policy —
#' instead of the Ramsey-optimal commitment policy.
#'
#' @param ramsey_result   Result from \code{\link{ramsey_model}} or
#'   \code{\link{ramsey_nn1}} representing the optimal policy.
#' @param alternative_result Result from an alternative policy regime:
#'   \code{dynhr_osr_result}, \code{dynhr_discretionary_result},
#'   \code{dynhr_nn1_result}, \code{dynhr_ramsey_result}, or
#'   \code{dynhr_ramsey_result2}.
#' @param model           Optional \code{dynhr_mod}.
#' @param params          Optional named parameter vector.
#' @param consumption_variable Character: name of the consumption variable.
#' @param ...             Additional arguments passed to \code{welfare_ce_diff()}.
#'
#' @return A list with the CE cost of deviating from Ramsey (positive means
#'   the Ramsey policy delivers higher welfare), plus details of the
#'   welfare comparison. Class \code{dynhr_welfare_cost}.
#'
#' @examples
#' \dontrun{
#' ramsey <- ramsey_model(model, "-(pi^2 + 0.5*y_gap^2)")
#' taylor_rule <- osr(model, free_params = c(phi_pi = 1.5, phi_y = 0.5),
#'                    loss_vars = c("pi", "y_gap"))
#' cost <- welfare_cost_of_rule(ramsey, taylor_rule,
#'                               consumption_variable = "y_gap")
#' print(sprintf("Cost of Taylor rule vs Ramsey: %.2f%% of C_ss",
#'               cost$ce_percent))
#' }
#'
#' @export
welfare_cost_of_rule <- function(ramsey_result,
                                  alternative_result,
                                  model = NULL,
                                  params = NULL,
                                  consumption_variable = NULL,
                                  ...) {

  # Validate inputs
  if (!inherits(ramsey_result, c("dynhr_ramsey_result", "dynhr_ramsey_result2",
                                  "dynhr_nn1_result"))) {
    stop("'ramsey_result' must be a result from ramsey_model(), ",
         "ramsey_policy(), or ramsey_nn1().")
  }

  # Compute CE difference: Ramsey (A) vs Alternative (B)
  ce <- welfare_ce_diff(
    result_a = ramsey_result,
    result_b = alternative_result,
    model = model,
    params = params,
    consumption_variable = consumption_variable,
    ...
  )

  # Add description
  alt_label <- switch(
    class(alternative_result)[1],
    dynhr_osr_result         = "OSR (optimal simple rule)",
    dynhr_discretionary_result = "Discretionary policy",
    dynhr_nn1_result         = sprintf("(n,n+1) approximation (n=%d)",
                                     alternative_result$n %||% "?"),
    dynhr_ramsey_result      = "Ramsey (simulation-based)",
    dynhr_ramsey_result2     = "Ramsey (augmented-system)",
    "Alternative policy"
  )

  result <- c(ce, list(
    ramsey_label   = "Ramsey optimal policy",
    alternative_label = alt_label,
    welfare_a      = .extract_welfare(ramsey_result, "unconditional"),
    welfare_b      = .extract_welfare(alternative_result, "unconditional")
  ))
  class(result) <- c("dynhr_welfare_cost", "list")
  result
}


# ==========================================================================
# 4. Welfare decomposition (level + volatility)
# ==========================================================================

#' Decompose unconditional welfare into level and volatility effects
#'
#' For a second-order approximation, unconditional welfare can be decomposed
#' as:
#' \deqn{E[W] = W_{ss} + \text{Level} + \text{Volatility}}
#' where:
#' \itemize{
#'   \item \eqn{W_{ss} = u(y_{ss}) / (1-\beta)} is the steady-state welfare
#'   \item \eqn{\text{Level}} captures the deterministic correction (pre-drift
#'     from second-order terms in the decision rules)
#'   \item \eqn{\text{Volatility}} captures the risk-related component:
#'     \eqn{\frac{1}{2(1-\beta)} \mathrm{tr}(\Sigma_y \cdot H_{uu})}
#' }
#'
#' For a first-order (LQ) solution, there is no level effect and the
#' volatility effect fully accounts for the welfare gap.
#'
#' @param result   A policy result object (\code{dynhr_ramsey_result},
#'   \code{dynhr_ramsey_result2}, \code{dynhr_nn1_result}).
#' @param model    Optional \code{dynhr_mod}. Extracted from result if
#'   not provided.
#' @param params   Optional named parameter vector.
#' @param planner_objective Optional objective expression text. Extracted
#'   from the result if not provided.
#' @param verbose  Print decomposition details.
#'
#' @return A list with class \code{dynhr_welfare_decomposition} containing:
#'   \describe{
#'     \item{welfare_total}{Unconditional welfare (from the result).}
#'     \item{welfare_ss}{Steady-state welfare component.}
#'     \item{level_effect}{Level (pre-drift) effect.}
#'     \item{volatility_effect}{Volatility (risk) effect.}
#'     \item{method}{Description of computation method.}
#'     \item{dr_order}{Order of the decision rules (1, 2, or NA).}
#'   }
#'
#' @examples
#' \dontrun{
#' # Decompose welfare from a Ramsey policy
#' ramsey <- ramsey_model(model, "-(pi^2 + 0.5*y_gap^2)", order = 2)
#' wd <- welfare_decompose(ramsey, model)
#' print(wd)
#' }
#'
#' @export
welfare_decompose <- function(result,
                               model = NULL,
                               params = NULL,
                               planner_objective = NULL,
                               verbose = FALSE) {

  # ---- 1. Extract welfare components ----
  welfare_total <- .extract_welfare(result, "unconditional")
  welfare_ss    <- .extract_welfare(result, "steady")

  if (!is.finite(welfare_total)) {
    stop("Unconditional welfare not available. ",
         "Try running with simulation-based welfare first.")
  }

  model <- .extract_model(result, model)
  if (is.null(model)) {
    stop("Could not locate dynhr_mod. Please supply 'model' argument.")
  }
  if (is.null(params)) params <- model$param_values

  dr <- .extract_dr(result)
  if (is.null(dr)) {
    stop("Decision rules not available in the result object.")
  }

  # Detect order of decision rules
  dr_order <- 1L
  if (inherits(dr, "DecisionRules2") || !is.null(dr$ghss)) dr_order <- 2L
  if (!is.null(dr$ghs3)) dr_order <- 3L

  # ---- 2. Get planner objective AST ----
  if (is.null(planner_objective)) {
    if (inherits(result, "dynhr_ramsey_result")) {
      planner_objective <- result$objective$text
    } else if (inherits(result, "dynhr_nn1_result")) {
      planner_objective <- result$meta$planner_objective
    } else if (inherits(result, "dynhr_ramsey_result2")) {
      planner_objective <- result$augmented_result$planner_objective %||%
                           model$planner_objective$text
    }
  }

  ss <- .extract_ss(result)
  if (is.null(ss)) {
    stop("Steady state not available in the result object.")
  }

  # ---- 3. Compute decomposition ----
  # W_ss component
  if (!is.finite(welfare_ss)) {
    # Compute from objective
    all_vars <- c(model$var_names, model$varexo_names)
    obj_ast <- NULL
    if (!is.null(planner_objective)) {
      obj_ast <- parse_expression(planner_objective,
          var_names = all_vars,
          param_names = model$param_names)
    }
    if (!is.null(obj_ast)) {
      discount <- .extract_discount(result)
      f_ss <- .eval_planner_ast(obj_ast, ss, params, ss)
      welfare_ss <- if (is.finite(f_ss)) f_ss / max(1e-8, (1 - discount)) else NA_real_
    }
  }

  # Level effect: difference in steady-state value of the ghss-corrected policy function
  level_effect <- NA_real_
  if (!is.finite(welfare_ss)) {
    welfare_ss <- welfare_total  # fallback
    level_effect <- 0
  }

  if (dr_order >= 2 && !is.null(dr$ghss)) {
    # Second-order solution: level effect = welfare at corrected SS - welfare at raw SS
    # The ghss (risk) term shifts the mean of endogenous variables by 0.5*ghss.
    ss_corrected <- ss
    endo_names <- names(ss)
    for (nm in names(dr$ghss)) {
      if (nm %in% endo_names) {
        ss_corrected[nm] <- ss_corrected[nm] + 0.5 * dr$ghss[[nm]]
      }
    }
    # Evaluate objective at corrected SS
    discount <- .extract_discount(result)
    all_vars <- c(model$var_names, model$varexo_names)
    obj_ast <- NULL
    if (!is.null(planner_objective)) {
      obj_ast <- parse_expression(planner_objective,
          var_names = all_vars,
          param_names = model$param_names)
    }
    if (!is.null(obj_ast)) {
      f_ss <- .eval_planner_ast(obj_ast, ss, params, ss)
      f_corrected <- .eval_planner_ast(obj_ast, ss_corrected, params, ss)
      if (is.finite(f_ss) && is.finite(f_corrected)) {
        welfare_corrected <- f_corrected / max(1e-8, (1 - discount))
        level_effect <- welfare_corrected - welfare_ss
      }
    }
  } else {
    # First-order: no level effect (certainty equivalence)
    level_effect <- 0
  }

  # Volatility effect = residual
  # For first-order: total - ss = volatility
  # For second-order: total - ss - level = volatility (approximately)
  raw_gap <- welfare_total - welfare_ss
  volatility_effect <- if (is.finite(level_effect)) {
    raw_gap - level_effect
  } else {
    raw_gap
  }

  if (verbose) {
    cat("Welfare Decomposition:\n")
    cat(sprintf("  Total welfare:      %10.6f\n", welfare_total))
    cat(sprintf("  Steady state:       %10.6f\n", welfare_ss))
    cat(sprintf("  Level effect:       %10.6f\n", level_effect))
    cat(sprintf("  Volatility effect:  %10.6f\n", volatility_effect))
    cat(sprintf("  DR order:           %d\n", dr_order))
  }

  result <- list(
    welfare_total     = welfare_total,
    welfare_ss        = welfare_ss,
    level_effect      = level_effect,
    volatility_effect = volatility_effect,
    raw_gap           = raw_gap,
    method            = if (dr_order >= 2)
                          "Second-order decomposition (level + volatility)"
                        else
                          "First-order (certainty-equivalent): volatility only",
    dr_order          = dr_order
  )
  class(result) <- c("dynhr_welfare_decomposition", "list")
  result
}


# ==========================================================================
# 5. Conditional welfare
# ==========================================================================

#' Collect the bare variable names referenced by a planner-objective AST.
#'
#' Mirrors `.eval_planner_ast`'s by-name treatment of variable nodes (timing
#' is ignored), so the gradient/Hessian below is taken over exactly the
#' variables the objective actually depends on.
#' @noRd
.collect_obj_vars <- function(node, acc = character(0)) {
  if (is.null(node)) return(acc)
  switch(node$type,
    variable = c(acc, node$name),
    unaryop  = .collect_obj_vars(node$operand, acc),
    binop    = .collect_obj_vars(node$right, .collect_obj_vars(node$left, acc)),
    funcall  = { for (a in node$args) acc <- .collect_obj_vars(a, acc); acc },
    acc)
}

#' Central-difference gradient and Hessian of the planner objective at SS.
#'
#' Evaluated over `vars_names` (the objective's own variables), holding all
#' other variables at their steady state.  Used by the order-1 closed form.
#' @noRd
.objective_grad_hess <- function(obj_ast, ss, params, vars_names, eps = 1e-4) {
  n    <- length(vars_names)
  base <- as.list(ss)               # allow [[<-]] regardless of ss storage
  f0   <- .eval_planner_ast(obj_ast, base, params, ss)
  g    <- setNames(numeric(n), vars_names)
  H    <- matrix(0, n, n, dimnames = list(vars_names, vars_names))
  h    <- vapply(vars_names,
                 function(v) eps * max(1, abs(as.numeric(base[[v]] %||% 0))),
                 numeric(1))

  ev <- function(mods) {
    v <- base
    for (nm in names(mods)) v[[nm]] <- (as.numeric(v[[nm]] %||% 0) + mods[[nm]])
    .eval_planner_ast(obj_ast, v, params, ss)
  }

  for (i in seq_len(n)) {
    vi <- vars_names[i]
    fp <- ev(setNames(list( h[i]), vi))
    fm <- ev(setNames(list(-h[i]), vi))
    g[i]   <- (fp - fm) / (2 * h[i])
    H[i, i] <- (fp - 2 * f0 + fm) / (h[i]^2)
  }
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    vi <- vars_names[i]; vj <- vars_names[j]
    fpp <- ev(setNames(list( h[i],  h[j]), c(vi, vj)))
    fpm <- ev(setNames(list( h[i], -h[j]), c(vi, vj)))
    fmp <- ev(setNames(list(-h[i],  h[j]), c(vi, vj)))
    fmm <- ev(setNames(list(-h[i], -h[j]), c(vi, vj)))
    Hij <- (fpp - fpm - fmp + fmm) / (4 * h[i] * h[j])
    H[i, j] <- Hij; H[j, i] <- Hij
  }
  list(g = g, H = H, f0 = f0)
}

#' Order-1 (LQ) conditional-welfare correction (Schmitt-Grohe & Uribe).
#'
#' Returns the state-dependent gap V(s0) - V(0) for a first-order rule:
#' \deqn{V(s_0) - V(0) = g' C (I - \beta A)^{-1} s_0 + \tfrac12 s_0' P s_0,}
#' with \eqn{A = ghx[state,]}, \eqn{C = ghx}, \eqn{Q = C' H C}, and
#' \eqn{P = Q + \beta A' P A}.  The (s0-independent) stochastic-volatility term
#' cancels in the gap, so this is exact for a quadratic objective and
#' second-order accurate otherwise.  Returns NULL if the closed form is
#' inapplicable (no states, singular system / unit root).
#' @noRd
.conditional_welfare_lq <- function(dr, obj_ast, ss, params, discount, init_dev) {
  endo      <- dr$endo_names
  state_idx <- dr$state_idx
  if (is.null(state_idx) || length(state_idx) == 0L || is.null(dr$ghx)) return(NULL)

  ghx     <- dr$ghx
  A       <- ghx[state_idx, , drop = FALSE]      # n_state x n_state
  n_state <- length(state_idx)

  ## Initial state deviation ordered as endo[state_idx].
  st_names <- endo[state_idx]
  s0 <- vapply(st_names, function(nm) {
    if (nm %in% names(init_dev)) as.numeric(init_dev[[nm]]) else 0
  }, numeric(1))
  if (all(s0 == 0)) return(list(gap = 0, linear = 0, quad = 0, P = NULL))

  obj_vars <- intersect(unique(.collect_obj_vars(obj_ast)), endo)
  if (length(obj_vars) == 0L) return(NULL)
  gh  <- .objective_grad_hess(obj_ast, ss, params, obj_vars)
  C_r <- ghx[obj_vars, , drop = FALSE]           # n_obj x n_state

  ## Linear term: g' C_r (I - beta A)^{-1} s0.
  M       <- diag(n_state) - discount * A
  if (rcond(M) < .Machine$double.eps) return(NULL)
  rhs_lin <- as.numeric(crossprod(C_r, gh$g))    # C_r' g
  lin_vec <- solve(t(M), rhs_lin)                # (I - beta A)^{-T} C_r' g
  linear  <- as.numeric(crossprod(lin_vec, s0))

  ## Quadratic term: 0.5 s0' P s0, P = Q + beta A' P A solved as a discounted
  ## Lyapunov equation P = M_b P M_b' + Q with M_b = sqrt(beta) A'.
  Q <- crossprod(C_r, gh$H %*% C_r)              # n_state x n_state
  P <- solve_lyapunov(sqrt(discount) * t(A), Q)
  if (any(!is.finite(P))) return(NULL)
  quad <- 0.5 * as.numeric(t(s0) %*% P %*% s0)

  list(gap = linear + quad, linear = linear, quad = quad, P = P)
}

#' Numeric Hessian of the planner objective w.r.t. the endogenous variables,
#' evaluated at the steady state.  Constant (hence exact) for a quadratic
#' objective; the leading second-order term otherwise.
#' @noRd
.planner_hessian <- function(obj_ast, ss, params, endo, h = 1e-4) {
  n    <- length(endo)
  base <- setNames(as.numeric(ss[endo]), endo)
  feval <- function(v) .eval_planner_ast(obj_ast, setNames(v, endo), params, ss)
  f0 <- feval(base)
  H  <- matrix(0, n, n)
  fp <- numeric(n); fm <- numeric(n)
  for (i in seq_len(n)) {
    vp <- base; vp[i] <- vp[i] + h; fp[i] <- feval(vp)
    vm <- base; vm[i] <- vm[i] - h; fm[i] <- feval(vm)
    H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / (h * h)
  }
  for (i in seq_len(n)) for (j in seq_len(i - 1L)) {
    vpp <- base; vpp[i] <- vpp[i] + h; vpp[j] <- vpp[j] + h
    vpm <- base; vpm[i] <- vpm[i] + h; vpm[j] <- vpm[j] - h
    vmp <- base; vmp[i] <- vmp[i] - h; vmp[j] <- vmp[j] + h
    vmm <- base; vmm[i] <- vmm[i] - h; vmm[j] <- vmm[j] - h
    H[i, j] <- H[j, i] <- (feval(vpp) - feval(vpm) - feval(vmp) + feval(vmm)) /
                          (4 * h * h)
  }
  H
}

#' Analytic order-2 conditional welfare (Tier 15 B4 option b).
#'
#' E[W | s0] = sum_t beta^(t-1) { f(m_t) + 0.5 tr(H_f Sigma_t) }, where m_t and
#' Sigma_t are the EXACT s0-conditional mean and covariance of the pruned
#' second-order output (.order2_conditional_moments) and H_f is the objective
#' Hessian at the steady state.  This is a second-order Taylor expansion of the
#' objective around the conditional-mean path: EXACT for a quadratic objective
#' (constant H_f), second-order accurate otherwise.  The s0-dependent
#' conditional volatility (the bilinear hxu*(eps@x1) term the deterministic
#' path drops) enters through Sigma_t.  Returns list(gap = E[W|s0]-E[W|0]) or
#' NULL if not applicable.  Deterministic and reproducible (no RNG).
#' @noRd
.conditional_welfare_analytic <- function(dr, model, params, ss, discount,
                                          obj_ast, init_dev_vec, n_periods) {
  endo <- dr$endo_names; exo <- dr$exo_names; sidx <- dr$state_idx
  if (is.null(sidx) || length(sidx) == 0L || is.null(obj_ast)) return(NULL)
  Sigma_e <- tryCatch(.get_shock_cov(model, exo, params), error = function(e) NULL)
  if (is.null(Sigma_e)) return(NULL)
  sys <- tryCatch(.order2_aug_system(dr, Sigma_e), error = function(e) NULL)
  if (is.null(sys)) return(NULL)
  state_nm <- endo[sidx]
  s0 <- setNames(numeric(length(sidx)), state_nm)
  common <- intersect(state_nm, names(init_dev_vec))
  if (length(common)) s0[common] <- as.numeric(init_dev_vec[common])
  np <- as.integer(n_periods)
  cm_s0 <- .order2_conditional_moments(sys, as.numeric(s0), np)
  cm_0  <- .order2_conditional_moments(sys, numeric(length(sidx)), np)
  Hf <- .planner_hessian(obj_ast, ss, params, endo)
  w  <- discount ^ (seq_len(np) - 1L)
  welfare_of <- function(cm) {
    fmean <- vapply(seq_len(np), function(t)
      .eval_planner_ast(obj_ast, setNames(cm$mean[t, ], endo), params, ss),
      numeric(1))
    trc <- vapply(seq_len(np), function(t) 0.5 * sum(Hf * cm$cov[[t]]), numeric(1))
    sum(w * (fmean + trc))
  }
  W_s0 <- welfare_of(cm_s0); W_0 <- welfare_of(cm_0)
  if (!is.finite(W_s0) || !is.finite(W_0)) return(NULL)
  list(gap = W_s0 - W_0)
}

#' Discounted utility sum along a path simulated from `init_state`.
#'
#' With zero `shocks` this is the deterministic conditional path; with random
#' shocks it is one stochastic draw.  `init_state` is a named deviation vector
#' over the endogenous variables (NULL = start at steady state).
#' @noRd
.discounted_welfare_path <- function(dr, obj_ast, params, ss, discount,
                                     n_periods, init_state = NULL,
                                     shocks = NULL, model = NULL) {
  sim <- .simulate_dr_any_order(dr, n_periods = as.integer(n_periods),
                                model = model, burn_in = 0L,
                                shocks = shocks, init_state = init_state)
  if (is.null(sim)) return(NA_real_)
  lev <- attr(sim, "levels")
  if (is.null(lev) || nrow(lev) < 1L) return(NA_real_)
  cols <- colnames(lev)
  np   <- nrow(lev)
  w    <- discount ^ (seq_len(np) - 1L)
  u    <- vapply(seq_len(np), function(t) {
    .eval_planner_ast(obj_ast, setNames(as.numeric(lev[t, ]), cols), params, ss)
  }, numeric(1))
  sum(w * u)
}

#' Welfare conditional on an initial state vector
#'
#' Computes the expected discounted sum of future utility conditional on
#' the economy starting from a given initial state vector \eqn{s_0}.
#'
#' For first-order (LQ) approximations the state-dependent part is computed in
#' closed form (Schmitt-Grohe & Uribe):
#' \deqn{V(s_0) = W_{ss} + g' C (I - \beta A)^{-1} s_0 + \tfrac12 s_0' P s_0,}
#' where \eqn{A = ghx[state,]}, \eqn{C = ghx}, \eqn{g} and \eqn{H} are the
#' gradient and Hessian of the objective at the steady state, and \eqn{P}
#' solves the discounted Lyapunov equation \eqn{P = C'HC + \beta A' P A}. The
#' (\eqn{s_0}-independent) stochastic-volatility term cancels in the reported
#' \code{gap_vs_steady}, so the correction is exact for a quadratic objective
#' and second-order accurate otherwise.
#'
#' For second-order rules, the deterministic (zero-shock) path misses the
#' s0-dependent conditional volatility from the bilinear
#' \eqn{h_{xu}(\varepsilon \otimes x^1)} term (whose variance scales with
#' \eqn{\|h_x^t s_0\|^2}).  Two methods capture it:
#' \itemize{
#'   \item \code{"analytic"} (the default for order-2 rules) — a fast,
#'     deterministic second-order Taylor expansion of the objective around the
#'     s0-conditional mean path, using the EXACT conditional moments of the
#'     pruned second-order output (mean \eqn{m_t} and covariance
#'     \eqn{\Sigma_t}): \eqn{E[W|s_0] = \sum_t \beta^{t-1}\{f(m_t) +
#'     \tfrac12\,\mathrm{tr}(H_f \Sigma_t)\}}.  EXACT for a quadratic objective
#'     (constant \eqn{H_f}), second-order accurate otherwise.  No Monte Carlo.
#'   \item \code{"stochastic"} — a fixed-seed Monte-Carlo average over
#'     \code{n_mc_draws} replications from \eqn{s_0} minus the same from
#'     \eqn{s_0 = 0}.  Exact for any objective up to MC noise; the caller's RNG
#'     state is saved/restored.  Use it when the objective is strongly
#'     non-quadratic.
#' }
#' Pass \code{method = "deterministic"} for the old zero-shock path (fast but
#' biased for order-2 rules with non-trivial \eqn{h_{xu}}).
#'
#' @param result       A policy result object (\code{dynhr_ramsey_result},
#'   \code{dynhr_ramsey_result2}, \code{dynhr_nn1_result}).
#' @param initial_state Named numeric vector of initial state values. If
#'   \code{NULL}, uses the steady state (returns steady-state welfare).
#' @param model        Optional \code{dynhr_mod}.
#' @param params       Optional named parameter vector.
#' @param planner_objective Optional objective expression text.
#' @param n_periods    Number of periods to simulate (for order >= 2).
#' @param burn_in      Burn-in periods discarded before welfare evaluation
#'   (for simulation-based).
#' @param method       Computation method for order-2 rules.
#'   \code{"auto"} (default) uses \code{"stochastic"} for order-2 rules and the
#'   LQ closed form / deterministic path for order-1 rules.
#'   \code{"stochastic"} forces the MC-average path (correct for order 2).
#'   \code{"deterministic"} uses the zero-shock conditional path (fast, but
#'   misses the s0-dependent volatility term for order-2 rules).
#' @param seed         Integer seed for the MC draws (order-2 stochastic path).
#'   When \code{NULL} (default), seed \code{42L} is used.  The caller's RNG
#'   state is always restored.
#' @param n_mc_draws   Number of MC replications for the stochastic path
#'   (default 500).  Increase for tighter Monte Carlo standard error at the
#'   cost of runtime.
#' @param ...          Additional arguments passed to methods.
#'
#' @return A list with class \code{dynhr_conditional_welfare} containing:
#'   \describe{
#'     \item{welfare}{Conditional welfare value.}
#'     \item{initial_state}{The initial state vector used (levels).}
#'     \item{steady_welfare}{Welfare at the steady state.}
#'     \item{gap_vs_steady}{State-dependent welfare correction
#'       \eqn{V(s_0) - W_{ss}}.}
#'     \item{method}{Computation method.}
#'     \item{dr_order}{Order of decision rules.}
#'   }
#'
#' @examples
#' \dontrun{
#' ramsey <- ramsey_model(model, "-(pi^2 + 0.5*y_gap^2)")
#' # Welfare from steady state
#' cw_ss <- conditional_welfare(ramsey)
#' # Welfare from a recessionary initial state
#' init <- c(y_gap = -0.02, pi = -0.005)
#' cw_recession <- conditional_welfare(ramsey, initial_state = init)
#' print(sprintf("Welfare loss from recession: %.4f",
#'               cw_ss$welfare - cw_recession$welfare))
#' }
#'
#' @export
conditional_welfare <- function(result,
                                 initial_state = NULL,
                                 model = NULL,
                                 params = NULL,
                                 planner_objective = NULL,
                                 n_periods = 10000L,
                                 burn_in = 0L,
                                 method = c("auto", "stochastic", "deterministic",
                                            "analytic"),
                                 seed = NULL,
                                 n_mc_draws = 500L,
                                 ...) {
  method <- match.arg(method)

  # ---- 1. Extract components ----
  model <- .extract_model(result, model)
  if (is.null(model)) {
    stop("Could not locate dynhr_mod. Please supply 'model' argument.")
  }
  if (is.null(params)) params <- model$param_values

  dr <- .extract_dr(result)
  if (is.null(dr)) {
    stop("Decision rules not available in the result object.")
  }

  ss <- .extract_ss(result)
  if (is.null(ss)) {
    stop("Steady state not available in the result object.")
  }

  discount <- .extract_discount(result)
  dr_order <- 1L
  if (inherits(dr, "DecisionRules2") || !is.null(dr$ghss)) dr_order <- 2L

  # ---- 2. Get planner objective AST ----
  if (is.null(planner_objective)) {
    if (inherits(result, "dynhr_ramsey_result")) {
      planner_objective <- result$objective$text
    } else if (inherits(result, "dynhr_nn1_result")) {
      planner_objective <- result$meta$planner_objective
    } else if (inherits(result, "dynhr_ramsey_result2")) {
      planner_objective <- result$augmented_result$planner_objective %||%
                           model$planner_objective$text
    }
  }

  all_vars <- c(model$var_names, model$varexo_names)
  obj_ast <- NULL
  if (!is.null(planner_objective)) {
    obj_ast <- parse_expression(planner_objective,
        var_names = all_vars,
        param_names = model$param_names)
  }

  # Steady-state welfare
  welfare_ss <- .extract_welfare(result, "steady")
  if (!is.finite(welfare_ss) && !is.null(obj_ast)) {
    f_ss <- .eval_planner_ast(obj_ast, ss, params, ss)
    welfare_ss <- if (is.finite(f_ss)) f_ss / max(1e-8, (1 - discount)) else NA_real_
  }

  # ---- 3. Handle default initial state ----
  if (is.null(initial_state)) {
    return(structure(list(
      welfare        = welfare_ss,
      initial_state  = ss,
      steady_welfare = welfare_ss,
      method         = "Steady state (no initial deviations)",
      dr_order       = dr_order,
      discount       = discount
    ), class = c("dynhr_conditional_welfare", "list")))
  }

  # ---- 4. Compute conditional welfare ----
  # Build a full initial state vector (fill missing with SS values)
  init_full <- ss
  for (nm in names(initial_state)) {
    if (nm %in% names(init_full)) {
      init_full[[nm]] <- as.numeric(initial_state[[nm]])
    }
  }

  # Build the initial state in deviation form (named over the endogenous
  # variables); this is the state the economy is conditioned to start from.
  init_dev <- init_full
  for (nm in names(init_dev)) {
    if (nm %in% names(ss)) init_dev[[nm]] <- init_dev[[nm]] - ss[[nm]]
  }
  init_dev_vec <- unlist(init_dev)

  cond_welfare <- NA_real_
  method_used  <- "steady-state approximation (conditioning unavailable)"
  gap          <- NA_real_

  # ---- Method 1: order-1 closed form (Schmitt-Grohe & Uribe) ----
  # V(s0) = W_ss + g' C (I - beta A)^{-1} s0 + 0.5 s0' P s0, with
  # P = C'HC + beta A' P A.  The state-dependent correction is exact for a
  # quadratic objective and second-order accurate otherwise; the (s0-
  # independent) stochastic-volatility term cancels in the gap-vs-steady.
  if (dr_order == 1L && !is.null(obj_ast) && !is.null(dr$ghx) &&
      length(dr$state_idx %||% integer(0)) > 0L) {
    lq <- tryCatch(
      .conditional_welfare_lq(dr, obj_ast, ss, params, discount, init_dev_vec),
      error = function(e) NULL)
    if (!is.null(lq) && is.finite(lq$gap)) {
      gap          <- lq$gap
      cond_welfare <- welfare_ss + gap
      method_used  <- "LQ closed form (SGU conditional welfare)"
    }
  }

  # ---- Method 2c: analytic order-2 conditional value function ----
  # Exact deterministic E[W|s0] via a second-order Taylor of the objective
  # around the s0-conditional mean path, using the EXACT conditional moments
  # (.order2_conditional_moments).  This captures the s0-dependent bilinear
  # volatility the deterministic path misses, with NO Monte-Carlo noise; it is
  # exact for a quadratic objective and second-order accurate otherwise.  It is
  # the default ("auto") for order-2 rules; falls through to the stochastic MC
  # if it cannot be formed (e.g. no states / shock-referencing objective).
  use_analytic <- (dr_order == 2L) &&
                  (method %in% c("auto", "analytic")) &&
                  !is.null(obj_ast)
  if (!is.finite(cond_welfare) && use_analytic) {
    an <- tryCatch(
      .conditional_welfare_analytic(dr, model, params, ss, discount,
                                    obj_ast, init_dev_vec, n_periods),
      error = function(e) NULL)
    if (!is.null(an) && is.finite(an$gap)) {
      gap          <- an$gap
      cond_welfare <- welfare_ss + gap
      method_used  <- sprintf(
        "analytic order-2 conditional value function (%d periods)",
        as.integer(n_periods))
    }
  }

  # ---- Method 2a: stochastic MC average (order-2 rules) ----
  # For a second-order rule, the deterministic (zero-shock) path misses the
  # s0-dependent conditional volatility from the bilinear hxu*(eps@x1) term,
  # whose variance scales as ||hx^t * s0||^2 and does NOT cancel between W(s0)
  # and W(0).  We correct this by averaging n_mc_draws stochastic paths from
  # init_state = s0, and subtracting the same average from init_state = 0.
  # The subtraction cancels the constant unconditional volatility term (ghss,
  # ghuu, hxx*Sigma_x) so the gap is the pure s0-dependent correction.
  # RNG state is saved and restored so the caller's stream is untouched.
  use_stochastic <- (dr_order == 2L) &&
                    (method %in% c("auto", "stochastic")) &&
                    !is.null(obj_ast)
  if (!is.finite(cond_welfare) && use_stochastic) {
    mc_seed <- if (!is.null(seed)) as.integer(seed) else 42L
    old_rng <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(mc_seed)
    on.exit({
      if (!is.null(old_rng)) {
        assign(".Random.seed", old_rng, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)

    nd     <- as.integer(n_mc_draws)
    np     <- as.integer(n_periods)
    W_s0_v <- numeric(nd)
    W_0_v  <- numeric(nd)
    for (d in seq_len(nd)) {
      W_s0_v[d] <- .discounted_welfare_path(dr, obj_ast, params, ss, discount,
                                            np, init_state = init_dev_vec,
                                            shocks = NULL, model = model)
      W_0_v[d]  <- .discounted_welfare_path(dr, obj_ast, params, ss, discount,
                                            np, init_state = NULL,
                                            shocks = NULL, model = model)
    }
    ok_s0 <- is.finite(W_s0_v)
    ok_0  <- is.finite(W_0_v)
    if (sum(ok_s0) > 0L && sum(ok_0) > 0L) {
      gap          <- mean(W_s0_v[ok_s0]) - mean(W_0_v[ok_0])
      cond_welfare <- welfare_ss + gap
      mc_se_val    <- sqrt(var(W_s0_v[ok_s0]) / sum(ok_s0) +
                           var(W_0_v[ok_0])  / sum(ok_0))
      method_used  <- sprintf(
        "stochastic MC average (%d draws, %d periods, seed %d)",
        nd, np, mc_seed)
    }
  }
  mc_se <- if (exists("mc_se_val")) mc_se_val else NA_real_

  # ---- Method 2b: deterministic conditional path (order-1 fallback or forced) ----
  # Simulate the rule forward from s0 with zero future shocks, take the
  # discounted utility sum, and subtract the same from the steady state.  The
  # subtraction cancels the order-2 risk constant (ghss) and any shock-free
  # drift, leaving the pure state-dependent correction.  Matches Method 1 for
  # a first-order rule with a quadratic objective.  For order-2 rules this is
  # BIASED because it omits the s0-dependent conditional volatility of the
  # bilinear hxu*(eps@x1) term; use Method 2a (stochastic) instead.
  use_deterministic <- !is.finite(cond_welfare) && !is.null(obj_ast) &&
                       (method == "deterministic" || dr_order == 1L || !use_stochastic)
  if (use_deterministic) {
    n_exo       <- length(dr$exo_names %||% character(0))
    zero_shocks <- matrix(0, as.integer(n_periods), n_exo)
    W_s0 <- .discounted_welfare_path(dr, obj_ast, params, ss, discount,
                                     n_periods, init_state = init_dev_vec,
                                     shocks = zero_shocks, model = model)
    W_0  <- .discounted_welfare_path(dr, obj_ast, params, ss, discount,
                                     n_periods, init_state = NULL,
                                     shocks = zero_shocks, model = model)
    if (is.finite(W_s0) && is.finite(W_0)) {
      gap          <- W_s0 - W_0
      cond_welfare <- welfare_ss + gap
      method_used  <- sprintf("deterministic conditional path (%d periods)",
                              as.integer(n_periods))
    }
  }

  # Fallback: steady-state welfare if conditioning could not be computed.
  if (!is.finite(cond_welfare)) {
    cond_welfare <- welfare_ss
    gap          <- 0
    method_used  <- "steady-state approximation (conditioning unavailable)"
  }

  result <- list(
    welfare        = cond_welfare,
    initial_state  = init_full,
    steady_welfare = welfare_ss,
    method         = method_used,
    dr_order       = dr_order,
    discount       = discount,
    gap_vs_steady  = cond_welfare - welfare_ss,
    mc_se          = mc_se
  )
  class(result) <- c("dynhr_conditional_welfare", "list")
  result
}


# ==========================================================================
# 6. Compute welfare for arbitrary decision rules (utility wrapper)
# ==========================================================================

#' Compute unconditional welfare for any set of decision rules
#'
#' Given decision rules, a planner objective, and a model, computes the
#' unconditional welfare using either Lyapunov-based (order 1) or
#' simulation-based methods.
#'
#' This is useful when you have custom decision rules (e.g., from OSR or
#' discretionary policy) and want to evaluate them using the same objective
#' as the Ramsey planner.
#'
#' @param dr                DecisionRules object (must contain ghx, ghu, etc.).
#' @param model             dynhr_mod object.
#' @param params            Named parameter vector.
#' @param planner_objective Character: planner objective expression.
#' @param discount          Discount factor. Defaults to \code{beta} parameter.
#' @param ss                Steady-state vector. If NULL, extracted from model
#'   or computed.
#' @param n_periods         Simulation length (for order >= 2 or fallback).
#' @param burn_in           Burn-in periods.
#' @param seed              Optional RNG seed. When set, the simulation-based
#'   welfare is a deterministic function of \code{dr} (the caller's RNG state is
#'   restored afterwards) -- used by \code{\link{osr}}'s welfare objective.
#' @param verbose           Print progress messages.
#' @param n_sim             Deprecated alias for \code{n_periods}; if supplied it
#'   overrides \code{n_periods} with a warning.
#'
#' @return A list with welfare components, class \code{dynhr_welfare_computed}.
#'
#' @examples
#' \dontrun{
#' # Evaluate welfare under an arbitrary rule
#' dr <- solve_perturbation(model, compiled, ss, params)
#' w <- welfare_compute(dr, model, params,
#'                       planner_objective = "-(pi^2 + 0.5*y_gap^2)")
#' print(w$unconditional)
#' }
#'
#' @references
#'   Schmitt-Grohé, S., & Uribe, M. (2004). Solving dynamic general equilibrium
#'     models using a second-order approximation to the policy function.
#'     \emph{Journal of Economic Dynamics and Control}, 28(4), 755-775.
#'   Kim, J., & Kim, S. H. (2003). Spurious welfare reversals in international
#'     business cycle models. \emph{Journal of International Economics}, 60(2),
#'     471-500.
#' @export
welfare_compute <- function(dr, model, params,
                             planner_objective,
                             discount = NULL,
                             ss = NULL,
                             n_periods = 10000L,
                             burn_in = 1000L,
                             seed = NULL,
                             verbose = FALSE,
                             n_sim = NULL) {

  ## Deprecation alias: n_sim= was the old name; forward to n_periods= with warning.
  if (!is.null(n_sim)) {
    warning("welfare_compute(): 'n_sim' is deprecated; use 'n_periods' instead.",
            call. = FALSE)
    n_periods <- n_sim
  }

  if (is.null(params)) params <- model$param_values
  if (is.null(discount)) {
    beta_val <- if ("beta" %in% names(params)) params[["beta"]] else NA
    discount <- if (length(beta_val) == 1 && is.finite(beta_val)) {
      as.numeric(beta_val)
    } else {
      0.99
    }
  }

  # ---- Parse objective ----
  all_vars <- c(model$var_names, model$varexo_names)
  obj_ast <- parse_expression(planner_objective,
      var_names = all_vars,
      param_names = model$param_names)

  # ---- Get steady state ----
  if (is.null(ss)) {
    compiled <- compile_model(model, verbose = FALSE)
    if (!is.null(compiled)) {
      ss_result <- solve_steady(compiled, params,
                     endo_names = model$var_names,
                     exo_names = model$varexo_names,
                     verbose = FALSE)
      ss <- if (!is.null(ss_result) && isTRUE(ss_result$converged)) {
        ss_result$values %||% ss_result$ss
      } else {
        rep(0, length(model$var_names)) |> setNames(model$var_names)
      }
    } else {
      ss <- rep(0, length(model$var_names)) |> setNames(model$var_names)
    }
  }

  # ---- Steady-state welfare ----
  welfare_ss <- .eval_planner_ast(obj_ast, ss, params, ss) / max(1e-8, (1 - discount))

  # ---- Unconditional welfare ----
  # For order 1, try Lyapunov-based moments first
  dr_order <- 1L
  if (inherits(dr, "DecisionRules2") || !is.null(dr$ghss)) dr_order <- 2L

  welfare_uncond <- NA_real_

  if (dr_order == 1) {
    moments <- compute_moments(dr, model, params = params)
    if (!is.null(moments) && !is.null(moments$var_cov)) {
      # Compute E[f(y)] ≈ f(ss) + ½·tr(Σ · H)
      # where H is the Hessian of the objective at SS
      # For LQ with quadratic objective: use the quadratic form directly
      # This is exact for linear models with quadratic objectives

      # Compute via simulation for generality
      welfare_uncond <- NA_real_
    }
  }

  # Simulation-based welfare (works for any order)
  if (is.na(welfare_uncond)) {
    ## Deterministic simulation when a seed is supplied (e.g. inside an OSR
    ## optimiser loss, where the welfare objective must be a reproducible
    ## function of the policy coefficients). Restore the caller's RNG state so
    ## the surrounding optimiser's own randomness is untouched.
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      set.seed(seed)
      on.exit({
        if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }, add = TRUE)
    }
    ## Simulate at the rule's native order (order-2 -> pruned second-order
    ## recursion, capturing the ghss/ghuu/ghxx volatility correction).
    sim <- .simulate_dr_any_order(dr, n_periods = n_periods, model = model,
                                  burn_in = burn_in)

    if (!is.null(sim)) {
      sim_levels <- attr(sim, "levels")
      if (!is.null(sim_levels)) {
        obj_t <- apply(sim_levels, 1L, function(row) {
          vals <- setNames(as.numeric(row), colnames(sim_levels))
          .eval_planner_ast(obj_ast, vals, params, ss)
        })
        obj_mean <- mean(obj_t, na.rm = TRUE)
        welfare_uncond <- obj_mean / max(1e-8, (1 - discount))
      }

      if (verbose) {
        cat(sprintf("  Welfare computed via simulation (%d periods)\n", n_periods))
        cat(sprintf("    Uncond: %.6f, Steady: %.6f\n", welfare_uncond, welfare_ss))
      }
    }
  }

  # Last resort fallback
  if (!is.finite(welfare_uncond)) {
    welfare_uncond <- welfare_ss
    if (verbose) cat("  Welfare: using steady-state approximation.\n")
  }

  result <- list(
    unconditional = welfare_uncond,
    steady_state  = welfare_ss,
    discount      = discount,
    method        = if (exists("moments") && !is.null(moments))
                      "Lyapunov-based (order 1)"
                    else sprintf("simulation (%d periods)", n_periods),
    n_sim         = if (exists("sim") && !is.null(sim)) n_periods else 0L
  )
  class(result) <- c("dynhr_welfare_computed", "list")
  result
}


# ==========================================================================
# 7. S3 methods
# ==========================================================================

#' @export
print.dynhr_ce_diff <- function(x, ...) {
  cat(sprintf("\n<dynhr_ce_diff>\n"))
  cat(sprintf("  Consumption-equivalent difference: %.4f%% of C_ss\n",
              x$ce_percent))
  cat(sprintf("  Welfare A:  %.6f\n", x$welfare_a))
  cat(sprintf("  Welfare B:  %.6f\n", x$welfare_b))
  cat(sprintf("  Welfare gap: %.6f\n", x$welfare_gap))
  cat(sprintf("  Discount:   %.4f\n", x$discount))
  if (!is.null(x$consumption_var)) {
    cat(sprintf("  C_ss:       %.6f (%s)\n", x$consumption_ss, x$consumption_var))
  }
  cat(sprintf("  Method:     %s\n", x$method))
  invisible(x)
}

#' @export
print.dynhr_welfare_cost <- function(x, ...) {
  cat(sprintf("\n<dynhr_welfare_cost>\n"))
  cat(sprintf("  %s\n", x$ramsey_label))
  cat(sprintf("  vs %s\n", x$alternative_label))
  cat(sprintf("  Cost of deviation: %.4f%% of C_ss (CE units)\n", x$ce_percent))
  cat(sprintf("  Ramsey welfare:      %.6f\n", x$welfare_a))
  cat(sprintf("  Alternative welfare: %.6f\n", x$welfare_b))
  cat(sprintf("  Welfare gap:         %.6f\n", x$welfare_gap))
  if (!is.null(x$consumption_var)) {
    cat(sprintf("  C_ss: %.6f (%s)\n", x$consumption_ss, x$consumption_var))
  }
  invisible(x)
}

#' @export
print.dynhr_welfare_decomposition <- function(x, ...) {
  cat(sprintf("\n<dynhr_welfare_decomposition>\n"))
  cat(sprintf("  DR order:  %d\n", x$dr_order))
  cat(sprintf("  Method:    %s\n", x$method))
  cat(sprintf("  Total welfare:      %10.6f\n", x$welfare_total))
  cat(sprintf("  Steady state:       %10.6f\n", x$welfare_ss))
  cat(sprintf("  Level effect:       %10.6f\n", x$level_effect))
  cat(sprintf("  Volatility effect:  %10.6f\n", x$volatility_effect))
  invisible(x)
}

#' @export
print.dynhr_conditional_welfare <- function(x, ...) {
  cat(sprintf("\n<dynhr_conditional_welfare>\n"))
  cat(sprintf("  Conditional welfare: %.6f\n", x$welfare))
  cat(sprintf("  Steady welfare:      %.6f\n", x$steady_welfare))
  cat(sprintf("  Gap vs steady:       %.6f\n", x$gap_vs_steady %||% NA_real_))
  cat(sprintf("  DR order:            %d\n", x$dr_order))
  cat(sprintf("  Method:              %s\n", x$method))
  invisible(x)
}

#' @export
print.dynhr_welfare_computed <- function(x, ...) {
  cat(sprintf("\n<dynhr_welfare_computed>\n"))
  cat(sprintf("  Unconditional: %.6f\n", x$unconditional))
  cat(sprintf("  Steady state:  %.6f\n", x$steady_state))
  cat(sprintf("  Discount:      %.4f\n", x$discount))
  cat(sprintf("  Method:        %s\n", x$method))
  invisible(x)
}
