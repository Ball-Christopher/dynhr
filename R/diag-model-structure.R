## R/diag-model-structure.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## model_structure_summary() and extraction helpers
## --------------------------------------------------------------------------

# -- Object accessor helpers -----------------------------------------------

.get_decision_rules <- function(x) {
  if (inherits(x, "DecisionRules")) return(x)
  if (!is.null(x$dr) && is.list(x$dr))     return(x$dr)
  if (!is.null(x$decision_rules))           return(x$decision_rules)
  if (!is.null(x$ghx) && !is.null(x$ghu))  return(x)
  stop("Cannot find DecisionRules in the supplied object")
}

.get_ss_vector <- function(ss) {
  if (is.numeric(ss) && !is.list(ss)) return(ss)
  if (is.list(ss) && !is.null(ss$values)) return(ss$values)
  if (is.list(ss) && !is.null(ss$ss))     return(ss$ss)
  if (is.list(ss) && !is.null(ss$ys))     return(ss$ys)
  stop("Cannot extract numeric steady-state vector")
}

.get_model_object <- function(x) {
  if (inherits(x, "dynhr_mod")) return(x)
  if (!is.null(x$model) && inherits(x$model, "dynhr_mod")) return(x$model)
  if (!is.null(x$model) && is.list(x$model) && !is.null(x$model$var_names)) return(x$model)
  if (!is.null(x$var_names)) return(x)
  stop("Cannot find dynhr_mod in the supplied object")
}

.get_endo_names <- function(x) {
  m <- .get_model_object(x)
  for (f in c("var_names", "endo_names", "variables", "endo")) {
    if (!is.null(m[[f]])) return(m[[f]])
  }
  NULL
}

.get_exo_names <- function(x) {
  m <- .get_model_object(x)
  for (f in c("varexo_names", "exo_names", "shocks", "exo")) {
    v <- m[[f]]
    if (!is.null(v) && is.character(v)) return(v)
  }
  NULL
}


#' Model structure summary
#'
#' Extracts and reports the structural dimensions of a DSGE model: number of
#' endogenous variables, exogenous shocks, parameters, and equations;
#' variable classification (static, predetermined, forward-looking, mixed);
#' and, when decision rules are provided, the BK condition rank check,
#' state dimension, and norm of the policy matrices.
#'
#' This function is called automatically by the diagnostic orchestrator and
#' its output appears as the first entry in any diagnostic report.
#'
#' @param model  A \code{dynhr_mod} object, or any object from which one
#'   can be extracted (e.g. the output of \code{run_full_estimation()}).
#' @param dr     Optional \code{DecisionRules} object.  When provided,
#'   augments the summary with policy-matrix norms and state dimension.
#' @param ss     Optional named numeric steady-state vector.  When provided,
#'   the maximum absolute steady-state value is reported for scale context.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{Named list of structural dimensions and, when \code{dr}
#'     is supplied, BK-condition indicators.}
#'   \item{pass}{Logical -- TRUE if all BK conditions are satisfied; NA when
#'     \code{dr} is not supplied.}
#'   \item{plots}{Empty list (no plots produced by this diagnostic).}
#'   \item{summary}{Human-readable multi-line summary of model dimensions.}
#'
#' @export
model_structure_summary <- function(model, dr = NULL, ss = NULL) {

  m <- .get_model_object(model)

  endo        <- .get_endo_names(m)
  exo         <- .get_exo_names(m)
  n_endo      <- length(endo)
  n_exo       <- length(exo)
  param_names <- m$param_names
  param_vals  <- m$param_values
  n_params    <- length(param_names)
  n_equations <- length(m$equations)

  n_static        <- m$n_static        %||% NA
  n_predetermined <- m$n_predetermined %||% NA
  n_forward       <- m$n_forward       %||% NA
  n_mixed         <- m$n_mixed         %||% NA
  var_class       <- m$variable_classification
  is_linear       <- isTRUE(m$model_options$linear)

  # Decision rules
  eigenvalues  <- NULL; bk_satisfied <- NULL
  n_unstable   <- NULL; n_state <- NULL
  ghx_dim      <- NULL; ghu_dim <- NULL
  state_vars   <- NULL

  dr_inner <- tryCatch(if (!is.null(dr)) .get_decision_rules(dr) else NULL,
                       error = function(e) NULL)
  if (!is.null(dr_inner)) {
    eigenvalues  <- dr_inner$eigenvalues
    bk_satisfied <- dr_inner$bk_satisfied
    n_unstable   <- dr_inner$n_unstable
    n_state      <- dr_inner$n_state
    state_vars   <- dr_inner$state_vars
    if (!is.null(dr_inner$ghx)) ghx_dim <- dim(dr_inner$ghx)
    if (!is.null(dr_inner$ghu)) ghu_dim <- dim(dr_inner$ghu)
  }

  # Steady state - try to extract numeric vector.  Tolerant: a missing or
  # unrecognised `ss` must not crash the structural summary (L32) -- the rest
  # of the report is still useful without a steady-state vector.
  ss_values <- tryCatch(if (!is.null(ss)) .get_ss_vector(ss) else NULL,
                        error = function(e) NULL)
  if (is.null(ss_values) && !is.null(dr) && !is.null(dr$ss)) {
    ss_values <- tryCatch(.get_ss_vector(dr$ss), error = function(e) NULL)
  }
  # Ensure ss_values is a plain numeric vector, not a list
  if (is.list(ss_values)) ss_values <- unlist(ss_values)

  # Eigenvalue table
  eig_table <- NULL
  if (!is.null(eigenvalues)) {
    mods <- Mod(eigenvalues)
    eig_table <- data.frame(
      index    = seq_along(eigenvalues),
      real     = Re(eigenvalues),
      imag     = Im(eigenvalues),
      modulus  = mods,
      stable   = mods < 1,
      stringsAsFactors = FALSE
    )
  }

  # Parameter table.  `param_vals` may be SHORTER than `param_names` when some
  # parameters are computed in an external *_steadystate.m and never assigned a
  # top-level value (M17 family).  Align by name and fill NA for the unvalued
  # ones so the data.frame never trips "differing number of rows" (L32).
  param_table <- NULL
  if (!is.null(param_names)) {
    pv <- rep(NA_real_, length(param_names))
    names(pv) <- param_names
    if (!is.null(param_vals)) {
      nm <- names(param_vals)
      if (!is.null(nm)) {
        common <- intersect(param_names, nm)
        if (length(common)) pv[common] <- as.numeric(param_vals[common])
      } else if (length(param_vals) == length(param_names)) {
        pv[] <- as.numeric(param_vals)
      }
    }
    param_table <- data.frame(
      parameter = param_names,
      value     = unname(pv),
      stringsAsFactors = FALSE
    )
  }

  result <- list(
    n_endo          = n_endo,
    n_exo           = n_exo,
    n_params        = n_params,
    n_equations     = n_equations,
    endo_names      = endo,
    exo_names       = exo,
    param_names     = param_names,
    param_values    = param_vals,
    param_table     = param_table,
    n_static        = n_static,
    n_predetermined = n_predetermined,
    n_forward       = n_forward,
    n_mixed         = n_mixed,
    variable_classification = var_class,
    steady_state    = ss_values,
    eigenvalues     = eigenvalues,
    eigenvalue_table = eig_table,
    bk_satisfied    = bk_satisfied,
    n_unstable      = n_unstable,
    n_state         = n_state,
    state_vars      = state_vars,
    ghx_dim         = ghx_dim,
    ghu_dim         = ghu_dim,
    is_linear       = is_linear
  )

  summary_lines <- c(
    sprintf("Model: %d endo, %d exo, %d params, %d equations",
            n_endo, n_exo, n_params, n_equations),
    sprintf("Classification: %d static, %d predetermined, %d forward, %d mixed",
            n_static, n_predetermined, n_forward, n_mixed),
    sprintf("States: %d  |  Shocks: %d",
            n_state %||% n_predetermined, n_exo),
    sprintf("BK satisfied: %s  |  Unstable: %s / %s eigenvalues",
            as.character(bk_satisfied %||% "unknown"),
            as.character(n_unstable %||% "?"),
            as.character(length(eigenvalues) %||% "?")),
    sprintf("Linear: %s", is_linear),
    if (!is.null(ghx_dim))
      sprintf("ghx: %d x %d  |  ghu: %d x %d",
              ghx_dim[1], ghx_dim[2], ghu_dim[1], ghu_dim[2]),
    if (!is.null(ss_values))
      sprintf("Steady state: %d values (all zero = %s)",
              length(ss_values), all(ss_values == 0))
  )

  structure(
    list(
      result  = result,
      pass    = if (is.null(bk_satisfied)) NA else isTRUE(bk_satisfied),
      plots   = list(),
      summary = paste(Filter(Negate(is.null), summary_lines), collapse = "\n")
    ),
    class = "dynhr_diagnostic"
  )
}
