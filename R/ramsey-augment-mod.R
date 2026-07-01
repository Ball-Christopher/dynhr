## R/ramsey-augment-mod.R
## --------------------------------------------------------------------------
## Augmented .mod file generator for Ramsey optimal policy.
##
## Given a parsed dynhr_mod and a planner_objective expression, this module
## builds the symbolic Lagrangian, derives first-order conditions (FOCs) for
## each endogenous variable, and writes an augmented .mod file where:
##   - Lagrange multipliers MULT_1, ..., MULT_n are added as new variables
##   - Original model equations are preserved
##   - FOCs are appended as additional equations
##
## The augmented .mod is then re-parseable by parse_mod() and feedable into
## the existing steady-state solver and perturbation pipeline.
##
## References:
##   Bodenstein, M. & Guerrieri, L. (2019). Nash–Ramsey toolbox.
##   Schmitt-Grohé, S. & Uribe, M. (2004). Optimal fiscal and monetary
##     policy under imperfect competition.
## --------------------------------------------------------------------------

#' Generate an augmented .mod file for Ramsey optimal policy
#'
#' Takes a parsed dynhr_mod and a planner objective expression, constructs
#' the symbolic Lagrangian, derives first-order conditions via symbolic
#' differentiation of the AST, and writes an augmented .mod file.
#'
#' The Lagrangian is:
#'   \deqn{L = E_0 \sum_{t=0}^{\infty} \beta^t \left[ f(\cdot) + \sum_j \lambda_{j,t} g_j(\cdot) \right]}
#'
#' The FOC for each endogenous variable \eqn{y_{i,t}}, after dividing by
#' \eqn{\beta^t}, is:
#'   \deqn{\frac{\partial f}{\partial y_{i,t}} + \sum_j \left[ \lambda_{j,t} \frac{\partial g_j}{\partial y_{i,t}}(t) + \beta^{-1}\,\lambda_{j,t-1} \frac{\partial g_j}{\partial y_{i,t}}(t-1) + \beta\,\lambda_{j,t+1} \frac{\partial g_j}{\partial y_{i,t}}(t+1) \right] = 0}
#'
#' where the discount factors \eqn{\beta^{-1}} (on the lag-of-the-multiplier,
#' i.e. lead-of-the-variable contribution) and \eqn{\beta} (on the lead-of-the-
#' multiplier, i.e. lag-of-the-variable contribution) arise from the
#' \eqn{\beta^t} weighting in the Lagrangian.  Omitting them yields decision
#' rules that are wrong on every model whose constraints carry an intertemporal
#' term (e.g. the wrong commitment targeting rule on the NK Phillips curve).
#'
#' @param model            A dynhr_mod object.
#' @param planner_objective Character string: the planner objective expression.
#' @param prefix            Prefix for multiplier variable names (default "MULT").
#' @param discount          Discount factor for the Lagrangian.  May be a
#'   character string naming a model parameter (applied symbolically/exactly,
#'   e.g. \code{"beta"}), a numeric value (applied as a literal), or \code{NULL}
#'   (default) to auto-detect a parameter named \code{"beta"} or \code{"betta"}.
#'   If no discount can be resolved, a value of 1 is used (no discounting) and a
#'   warning is issued.
#' @param verbose           Print progress messages.
#'
#' @return A list with:
#'   \item{augmented_text}{Character string of the augmented .mod file.}
#'   \item{multiplier_names}{Character vector of multiplier variable names.}
#'   \item{foc_equations}{List of FOC equation ASTs (for introspection).}
#'   \item{original_var_names}{Original endogenous variable names.}
#'   \item{multiplier_map}{Named integer vector mapping multiplier name -> original equation index.}
#' @export
ramsey_augment_mod <- function(model,
                                planner_objective = NULL,
                                prefix = "MULT",
                                discount = NULL,
                                verbose = FALSE) {
  # ---- 1. Validate inputs ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object created by parse_mod().")
  }

  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective provided. Supply planner_objective= or ",
         "include planner_objective(...); in the .mod file.")
  }

  endo_names <- model$var_names
  n_endo <- length(endo_names)
  equations <- model$equations
  n_eq <- length(equations)
  all_var_names <- c(endo_names, model$varexo_names, model$varexo_det_names)
  param_names <- model$param_names

  if (n_eq == 0) {
    stop("Model has no equations parsed.")
  }
  # Allow n_eq < n_endo for Ramsey models where the planner instrument(s)
  # are not pinned down by an equation — the planner objective provides
  # the extra FOC(s).  The augmented system will have n_endo + n_eq
  # equations (original + FOCs) for n_endo + n_eq augmented vars.
  if (n_eq > n_endo) {
    stop(sprintf("Number of equations (%d) > number of endogenous vars (%d).",
                 n_eq, n_endo))
  }
  if (n_eq < n_endo) {
    if (verbose) {
      cat(sprintf("  Note: %d equations for %d vars (%d instrument(s)).\n",
                  n_eq, n_endo, n_endo - n_eq))
    }
  }

  if (verbose) {
    cat(sprintf("[ramsey_augment_mod] Parsing objective: %s\n", obj_text))
    cat(sprintf("  Endogenous vars: %d, Equations: %d\n", n_endo, n_eq))
  }

  # ---- 2. Parse planner objective ----
  obj_ast <- parse_expression(obj_text,
                              var_names = all_var_names,
                              param_names = param_names)

  # ---- 2b. Resolve the discount factor (H4 BUG 2) ----
  # The Lagrangian carries a beta^t weight, so the FOC for y_i(t) -- divided
  # by beta^t -- discounts the multiplier contributions originating from
  # neighbouring periods:
  #   * the y_i-at-LEAD contribution comes from g_j(t-1) with MULT_j(t-1), and
  #     beta^{t-1}/beta^t = 1/beta  -> multiply by beta^{-1}
  #   * the y_i-at-LAG  contribution comes from g_j(t+1) with MULT_j(t+1), and
  #     beta^{t+1}/beta^t = beta    -> multiply by beta
  # Apply the discount SYMBOLICALLY by parameter name when possible (exact);
  # fall back to a numeric literal, else 1 (no discounting) with a warning.
  disc_info <- .resolve_discount_ast(discount, model)
  discount_ast     <- disc_info$ast       # AST for beta
  discount_inv_ast <- disc_info$inv_ast   # AST for 1/beta
  if (verbose) {
    cat(sprintf("[ramsey_augment_mod] Discount: %s\n", disc_info$label))
  }

  # ---- 3. Derive FOCs ----
  multiplier_names <- paste0(prefix, "_", seq_len(n_eq))
  foc_equations <- vector("list", n_endo)

  for (i in seq_len(n_endo)) {
    var_name <- endo_names[i]

    if (verbose) {
      cat(sprintf("  Deriving FOC for %s ...\n", var_name))
    }

    # --- 3a. Derivative of planner objective w.r.t. y_i (current period) ---
    df <- ast_number(0)
    if (!is.null(obj_ast) && is.list(obj_ast) && !is.null(obj_ast$type)) {
      df <- ast_differentiate(obj_ast, var_name, 0L)
    }

    # --- 3b. Sum over equations: contributions from multipliers ----
    # For each equation j:
    #   - If y_i appears with lead_lag=0 in g_j:  MULT_j * dg_j/dy_i(0)
    #   - If y_i appears with lead_lag=-1 in g_j: MULT_j(+1) * dg_j/dy_i(-1)
    #   - If y_i appears with lead_lag=+1 in g_j: MULT_j(-1) * dg_j/dy_i(+1)
    #
    # The multiplier timing shift:
    #   y_i at lag in g_j(t) → y_i(t-1) → multiplier at t is MULT_j associated
    #     with g_j at time t. But from the Lagrangian, d(L)/d(y_i(t)) gets
    #     contributions from g_j(t+1) where y_i(t) is the lag. So MULT_j(t+1).
    #     In the augmented mod: MULT_j(+1)
    #   y_i at lead in g_j(t) → y_i(t+1) → MULT_j(t-1). In augmented mod: MULT_j(-1)
    #
    # So we check each lead_lag in the equation and create the appropriate
    # multiplier timing.

    # Build list of (ll, mult_ll) pairs for each equation
    # ll = the lead_lag at which y_i appears in the equation
    # mult_ll = the lead_lag of the multiplier for that contribution

    contributions <- list()

    for (j in seq_len(n_eq)) {
      eq <- equations[[j]]
      # Get the residual: lhs - rhs (should be 0)
      eq_residual <- equation_to_residual(eq)
      mult_name <- multiplier_names[j]

      # Collect where y_i appears in this equation (check both lhs and rhs)
      # We differentiate the residual w.r.t. y_i at each lead/lag

      for (target_ll in c(-1L, 0L, 1L)) {
        deriv <- {
          ast_differentiate(eq_residual, var_name, target_ll)
        }

        if (is.null(deriv)) next

        # Check if derivative is non-zero (not just the number 0)
        if (ast_is_zero(deriv)) next

        # Determine multiplier timing
        # target_ll = 0  → multiplier is current: MULT_j
        # target_ll = -1 → multiplier is lead:    MULT_j(+1)
        # target_ll = +1 → multiplier is lag:     MULT_j(-1)
        mult_ll <- 0L
        if (target_ll == -1L) mult_ll <- 1L    # MULT_j(+1)
        else if (target_ll == 1L) mult_ll <- -1L  # MULT_j(-1)

        # Time-shift the derivative into the FOC's time frame (H4).
        # When y_i appears at lead/lag in g_j, the contribution originates from
        # g_j(t-target_ll), so the derivative -- computed in g_j's own frame --
        # must have every variable lead_lag shifted by -target_ll before being
        # placed in the FOC at time t.  For target_ll == 0 this is a no-op.
        # Omitting this drops backward-looking dependencies on the original
        # variables for NONLINEAR constraints, under-dimensioning the state
        # vector (linear constraints have constant derivatives -> shift is a
        # no-op, which is why this bug was invisible to linear-model tests).
        deriv <- ast_shift_timing(deriv, -target_ll)

        # Create multiplier AST node with proper timing
        mult_node <- list(
          type = "variable",
          name = mult_name,
          lead_lag = mult_ll
        )

        # Term: MULT_j(t+mult_ll) * shifted derivative
        term <- ast_binop("*", mult_node, deriv)

        # Apply the Lagrangian discount factor (H4 BUG 2):
        #   target_ll == +1 (MULT_j(-1), variable at LEAD) -> multiply by 1/beta
        #   target_ll == -1 (MULT_j(+1), variable at LAG)  -> multiply by beta
        #   target_ll ==  0 (MULT_j current)               -> no factor
        # In the canonical NK commitment case the beta^{-1} exactly cancels the
        # structural beta on beta*pi(+1), leaving the textbook coefficient -1 on
        # MULT_1(-1) (and hence the targeting rule pi = -(lambda/kappa)(y-y(-1))).
        if (target_ll == 1L) {
          term <- ast_binop("*", discount_inv_ast, term)
        } else if (target_ll == -1L) {
          term <- ast_binop("*", discount_ast, term)
        }

        contributions <- c(contributions, list(term))
      }
    }

    # --- 3c. Assemble FOC: df/dy_i + sum(contributions) = 0 ---
    if (length(contributions) == 0) {
      # No contributions from any equation — FOC is just df/dy_i = 0
      # This is the "static" or "irrelevant" case; still include it
      foc_expr <- df
    } else {
      # Sum all contributions
      sum_contrib <- contributions[[1]]
      if (length(contributions) > 1) {
        for (k in 2:length(contributions)) {
          sum_contrib <- ast_binop("+", sum_contrib, contributions[[k]])
        }
      }
      foc_expr <- ast_binop("+", df, sum_contrib)
    }

    # Simplify and store
    foc_expr <- ast_simplify(foc_expr)
    foc_equations[[i]] <- foc_expr

    if (verbose && ast_is_zero(foc_expr)) {
      cat(sprintf("    Warning: FOC for %s is identically zero.\n", var_name))
    }
  }

  # ---- 4. Generate augmented .mod text ----
  augmented_text <- .build_augmented_mod_text(
    model = model,
    endo_names = endo_names,
    multiplier_names = multiplier_names,
    foc_equations = foc_equations,
    verbose = verbose
  )

  # ---- 5. Return ----
  multiplier_map <- setNames(seq_len(n_eq), multiplier_names)

  result <- list(
    augmented_text    = augmented_text,
    multiplier_names  = multiplier_names,
    foc_equations     = foc_equations,
    original_var_names = endo_names,
    multiplier_map    = multiplier_map
  )
  class(result) <- "dynhr_ramsey_augmented"
  result
}


#' Resolve the Ramsey discount factor into beta and 1/beta AST nodes
#'
#' Used by \code{ramsey_augment_mod()} to attach the Lagrangian discount factor
#' to the multiplier contributions in the FOCs.  When \code{discount} names a
#' model parameter, the factor is applied symbolically (a \code{parameter} AST
#' node), keeping the augmented .mod exact and re-parseable; the parameter's
#' calibrated value flows through automatically.  A numeric \code{discount}, or
#' a numeric auto-detected value, is applied as a literal.
#'
#' @param discount Character (parameter name), numeric value, or NULL (auto).
#' @param model    dynhr_mod (for parameter names/values when auto-detecting).
#' @return List with \code{ast} (beta), \code{inv_ast} (1/beta), and a
#'   human-readable \code{label}.
#' @noRd
.resolve_discount_ast <- function(discount, model) {
  param_names  <- model$param_names
  param_values <- model$param_values

  # 1. Character: a parameter name -> apply symbolically (exact).
  if (is.character(discount) && length(discount) == 1L && nzchar(discount)) {
    if (discount %in% param_names) {
      p <- ast_parameter(discount)
      return(list(
        ast     = p,
        inv_ast = ast_binop("/", ast_number(1), p),
        label   = sprintf("parameter '%s' (symbolic)", discount)
      ))
    }
    # Named but not a parameter: treat as numeric-string fallback if possible.
    num <- suppressWarnings(as.numeric(discount))
    if (!is.na(num)) {
      return(.discount_numeric_ast(num,
        label = sprintf("%g (from string '%s')", num, discount)))
    }
    stop("discount='", discount, "' is neither a model parameter nor numeric.")
  }

  # 2. Numeric: literal.
  if (is.numeric(discount) && length(discount) == 1L && is.finite(discount)) {
    return(.discount_numeric_ast(discount,
      label = sprintf("%g (numeric literal)", discount)))
  }

  # 3. NULL: auto-detect a parameter named 'beta' (or 'betta') and apply it
  #    symbolically.  Falls back to its numeric value, then to no discounting.
  if (is.null(discount)) {
    for (cand in c("beta", "betta")) {
      if (cand %in% param_names) {
        p <- ast_parameter(cand)
        return(list(
          ast     = p,
          inv_ast = ast_binop("/", ast_number(1), p),
          label   = sprintf("auto-detected parameter '%s' (symbolic)", cand)
        ))
      }
    }
    # No 'beta' parameter: try a calibrated value, else default to 1 (warn).
    for (cand in c("beta", "betta")) {
      if (cand %in% names(param_values) &&
          is.finite(param_values[[cand]])) {
        return(.discount_numeric_ast(as.numeric(param_values[[cand]]),
          label = sprintf("auto-detected value of '%s'", cand)))
      }
    }
    warning("ramsey_augment_mod: no discount factor found (no 'beta'/'betta' ",
            "parameter and none supplied via discount=). Using beta = 1 ",
            "(no discounting); commitment FOCs that depend on beta will be ",
            "wrong. Pass discount= explicitly.")
    return(.discount_numeric_ast(1, label = "1 (DEFAULT, no discounting)"))
  }

  stop("discount must be a parameter name (character), a numeric value, or NULL.")
}

#' Build beta / 1-beta AST nodes from a numeric discount value
#' @noRd
.discount_numeric_ast <- function(value, label) {
  list(
    ast     = ast_number(value),
    inv_ast = ast_number(1 / value),
    label   = label
  )
}


#' Build the augmented .mod file text
#'
#' Assembles the augmented .mod file from the original model and the derived
#' FOC equations. Preserves declarations, calibrations, shocks, and adds
#' multiplier variables and FOC equations.
#'
#' @param model             Original dynhr_mod.
#' @param endo_names        Original endogenous variable names.
#' @param multiplier_names  Multiplier variable names.
#' @param foc_equations     List of FOC ASTs.
#' @param verbose           Print progress.
#' @return Character string of the augmented .mod file.
#' @noRd
.build_augmented_mod_text <- function(model, endo_names,
                                       multiplier_names, foc_equations,
                                       verbose = FALSE) {
  lines <- character(0)

  # ---- Header ----
  lines <- c(lines, "// Augmented .mod file generated by dynhr ramsey_augment_mod()")
  lines <- c(lines, "// Adapted from: Bodenstein & Guerrieri (2019) Nash-Ramsey toolbox")
  lines <- c(lines, "")

  # ---- var block ----
  all_vars <- c(endo_names, multiplier_names)
  lines <- c(lines, "var")
  # Write variables in columns (up to 6 per line for readability)
  var_line <- "  "
  for (k in seq_along(all_vars)) {
    comma <- if (k < length(all_vars)) "," else ";"
    if (nchar(var_line) + nchar(all_vars[k]) + 1 > 72) {
      lines <- c(lines, var_line)
      var_line <- "  "
    }
    var_line <- paste0(var_line, all_vars[k], comma, " ")
  }
  if (nchar(trimws(var_line)) > 0) {
    lines <- c(lines, var_line)
  }
  lines <- c(lines, "")

  # ---- varexo block ----
  exo_names <- model$varexo_names
  if (length(exo_names) > 0) {
    lines <- c(lines, "varexo")
    exo_line <- "  "
    for (k in seq_along(exo_names)) {
      comma <- if (k < length(exo_names)) "," else ";"
      if (nchar(exo_line) + nchar(exo_names[k]) + 1 > 72) {
        lines <- c(lines, exo_line)
        exo_line <- "  "
      }
      exo_line <- paste0(exo_line, exo_names[k], comma, " ")
    }
    if (nchar(trimws(exo_line)) > 0) {
      lines <- c(lines, exo_line)
    }
    lines <- c(lines, "")
  }

  # ---- varexo_det block ----
  exo_det_names <- model$varexo_det_names
  if (length(exo_det_names) > 0) {
    lines <- c(lines, "varexo_det")
    det_line <- "  "
    for (k in seq_along(exo_det_names)) {
      comma <- if (k < length(exo_det_names)) "," else ";"
      if (nchar(det_line) + nchar(exo_det_names[k]) + 1 > 72) {
        lines <- c(lines, det_line)
        det_line <- "  "
      }
      det_line <- paste0(det_line, exo_det_names[k], comma, " ")
    }
    if (nchar(trimws(det_line)) > 0) {
      lines <- c(lines, det_line)
    }
    lines <- c(lines, "")
  }

  # ---- parameters block ----
  param_names <- model$param_names
  if (length(param_names) > 0) {
    lines <- c(lines, "parameters")
    param_line <- "  "
    for (k in seq_along(param_names)) {
      comma <- if (k < length(param_names)) "," else ";"
      if (nchar(param_line) + nchar(param_names[k]) + 1 > 72) {
        lines <- c(lines, param_line)
        param_line <- "  "
      }
      param_line <- paste0(param_line, param_names[k], comma, " ")
    }
    if (nchar(trimws(param_line)) > 0) {
      lines <- c(lines, param_line)
    }
    lines <- c(lines, "")
  }

  # ---- Calibration (parameter values) ----
  param_vals <- model$param_values
  if (length(param_vals) > 0) {
    lines <- c(lines, "// Parameter calibration")
    for (nm in names(param_vals)) {
      val <- param_vals[[nm]]
      if (is.finite(val)) {
        val_str <- format(val, scientific = FALSE, trim = TRUE)
        lines <- c(lines, sprintf("%s = %s;", nm, val_str))
      }
    }
    # Evaluate steady_state_model to fill in derived parameters (e.g. kappa,
    # Omega, lambda, vartheta) that are declared as parameters but computed
    # analytically in steady_state_model rather than calibrated directly.
    ssm <- model$steady_state_model
    if (length(ssm) > 0) {
      env <- new.env(parent = baseenv())
      for (nm in names(param_vals)) assign(nm, param_vals[[nm]], envir = env)
      for (assignment in ssm) {
        val <- tryCatch(eval(parse(text = assignment$text), envir = env),
                        error = function(e) NA)
        if (is.numeric(val) && length(val) == 1 && is.finite(val)) {
          assign(assignment$name, val, envir = env)
          # Emit calibration if not already present in param_vals
          if (!assignment$name %in% names(param_vals)) {
            val_str <- format(val, scientific = FALSE, trim = TRUE)
            lines <- c(lines, sprintf("%s = %s; // from steady_state_model",
                                      assignment$name, val_str))
          }
        }
      }
    }
    lines <- c(lines, "")
  }

  # ---- model block ----
  lines <- c(lines, "model;")
  if (isTRUE(model$model_options$linear)) {
    # Need to insert linear option carefully
  }

  # Original equations (with #local variables substituted)
  local_vars <- model$local_variables
  for (j in seq_along(model$equations)) {
    eq <- model$equations[[j]]
    lhs_ast <- ast_substitute_locals(eq$lhs, local_vars)
    rhs_ast <- ast_substitute_locals(eq$rhs, local_vars)
    lhs_str <- ast_to_string(lhs_ast)
    rhs_str <- ast_to_string(rhs_ast)
    # Preserve tags if any
    tag_str <- if (!is.na(eq$tag) && nzchar(eq$tag)) {
      sprintf(" [name='%s']", eq$tag)
    } else ""
    eq_text <- sprintf("  %s%s = %s;", lhs_str, tag_str, rhs_str)
    lines <- c(lines, eq_text)
  }

  # FOC equations
  lines <- c(lines, "")
  lines <- c(lines, "  // ---- Ramsey FOCs (derived analytically) ----")

  for (i in seq_along(foc_equations)) {
    foc <- foc_equations[[i]]
    if (ast_is_zero(foc)) {
      lines <- c(lines,
        sprintf("  // FOC for %s is identically zero (variable is static in objective)", endo_names[i]))
      next
    }
    foc <- ast_substitute_locals(foc, local_vars)
    foc_str <- ast_to_string(foc)
    lines <- c(lines, sprintf("  %s = 0; // FOC for %s", foc_str, endo_names[i]))
  }

  lines <- c(lines, "end;")
  lines <- c(lines, "")

  # ---- steady_state_model block (if present) ----
  # We deliberately omit the steady_state_model block for the augmented model
  # because the multipliers have no direct analytical SS representation.
  # The augmented SS will be solved numerically.

  # ---- initval block (if present) ----
  if (length(model$initval) > 0) {
    lines <- c(lines, "initval;")
    for (nm in names(model$initval)) {
      val <- model$initval[[nm]]
      if (is.finite(val)) {
        lines <- c(lines, sprintf("  %s = %s;", nm,
                                  format(val, scientific = FALSE, trim = TRUE)))
      }
    }
    # Add initial guesses for multipliers (small positive)
    for (mn in multiplier_names) {
      lines <- c(lines, sprintf("  %s = 0.01;", mn))
    }
    lines <- c(lines, "end;")
    lines <- c(lines, "")
  }

  # ---- shocks block ----
  shocks <- model$shocks
  has_shocks <- nrow(shocks$variances) > 0 || nrow(shocks$correlations) > 0
  if (has_shocks) {
    lines <- c(lines, "shocks;")
    if (nrow(shocks$variances) > 0) {
      for (k in seq_len(nrow(shocks$variances))) {
        row <- shocks$variances[k, ]
        if (!is.na(row$stderr) && is.finite(row$stderr)) {
          lines <- c(lines, sprintf("  var %s = %s;", row$name,
                                    format(row$stderr^2, scientific = FALSE, trim = TRUE)))
        } else if (!is.na(row$variance) && is.finite(row$variance)) {
          lines <- c(lines, sprintf("  var %s = %s;", row$name,
                                    format(row$variance, scientific = FALSE, trim = TRUE)))
        }
      }
    }
    if (nrow(shocks$correlations) > 0) {
      for (k in seq_len(nrow(shocks$correlations))) {
        row <- shocks$correlations[k, ]
        lines <- c(lines, sprintf("  corr %s, %s = %s;", row$var1, row$var2,
                                  format(row$corr, scientific = FALSE, trim = TRUE)))
      }
    }
    lines <- c(lines, "end;")
    lines <- c(lines, "")
  }

  # ---- Estimated params block ----
  if (nrow(model$estimated_params) > 0) {
    lines <- c(lines, "estimated_params;")
    ep <- model$estimated_params
    for (k in seq_len(nrow(ep))) {
      # Format varies by type
      if (ep$type[k] == "stderr") {
        lines <- c(lines, sprintf("  stderr %s, %s, %s, %s, %s, %s, %s;",
          ep$name[k], ep$prior[k], ep$p1[k], ep$p2[k], ep$p3[k], ep$p4[k], ""))
      } else {
        lines <- c(lines, sprintf("  %s, %s, %s, %s, %s, %s, %s;",
          ep$name[k], ep$prior[k], ep$p1[k], ep$p2[k], ep$p3[k], ep$p4[k], ""))
      }
    }
    lines <- c(lines, "end;")
    lines <- c(lines, "")
  }

  # ---- Model-local variables (# definitions) ----
  # Skipped: the original model's #-local-variable definitions are already
  # substituted into the equations during parsing.  Including them as bare
  # `#name = expr;` lines outside the model block is invalid Dynare syntax.
  # The FOCs and original equations are written with fully expanded forms.

  result <- paste(lines, collapse = "\n")
  result
}


#' Print method for dynhr_ramsey_augmented
#' @param x   A \code{dynhr_ramsey_augmented} object.
#' @param ... Unused; included for S3 compatibility.
#' @export
print.dynhr_ramsey_augmented <- function(x, ...) {
  cat(sprintf("<dynhr_ramsey_augmented>\n"))
  cat(sprintf("  Original vars:     %d\n", length(x$original_var_names)))
  cat(sprintf("  Multipliers:       %d\n", length(x$multiplier_names)))
  cat(sprintf("  FOC equations:     %d\n", length(x$foc_equations)))
  cat(sprintf("  Augmented .mod:\n"))
  lines <- strsplit(x$augmented_text, "\n")[[1]]
  for (l in head(lines, 20)) {
    cat(sprintf("    | %s\n", l))
  }
  if (length(lines) > 20) {
    cat(sprintf("    | ... (%d more lines)\n", length(lines) - 20))
  }
  invisible(x)
}


#' Convert augmented model text to a temporary file and parse it
#'
#' Convenience function: writes the augmented .mod text to a temporary file
#' and re-parses it with parse_mod().
#'
#' @param aug_result Result from ramsey_augment_mod().
#' @param verbose    Print progress.
#' @return A dynhr_mod object for the augmented system.
#' @export
ramsey_parse_augmented <- function(aug_result, verbose = FALSE) {
  if (!inherits(aug_result, "dynhr_ramsey_augmented")) {
    stop("aug_result must be from ramsey_augment_mod().")
  }

  ## parse_mod() accepts inline .mod text directly, so feed augmented_text
  ## straight in -- avoids a tempfile write+read round-trip per call (L10 perf).
  aug_model <- parse_mod(aug_result$augmented_text, verbose = verbose)
  aug_model
}
