## R/obc-spec.R
## --------------------------------------------------------------------------
## OBC constraint specification: model validation and MCP tag parsing.
##
## Provides:
##   obc_assert_linear()       -- guard: OBC only valid for model(linear)
##   obc_parse_tags()          -- extract OBC specs from [mcp = '...'] tags
##   obc_resolve_block_specs() -- resolve raw occbin_constraints block specs
##   obc_collect_specs()       -- primary entry point: merge tag + block specs
##
## An OBC spec is a list with fields:
##   $eq_idx   -- equation index in model$equations (1-based)
##   $var_idx  -- index of the constrained variable in model$var_names (1-based)
##   $var_name -- name of the constrained variable (for diagnostics)
##   $op       -- ">" (lower bound) or "<" (upper bound)
##   $bound    -- numeric bound value (in model deviation units)
## --------------------------------------------------------------------------


# =============================================================================
# Validation
# =============================================================================

#' Validate OBC model linearity
#'
#' The piecewise-linear Kalman filter (PKF) operates on the **first-order**
#' decision rules regardless of the perturbation order used for the main
#' solution.  The first-order DR is a valid linear approximation around the
#' deterministic steady state for both \code{model(linear)} and fully
#' nonlinear \code{model;} blocks.
#'
#' For nonlinear models we therefore issue a \emph{warning} rather than an
#' error: the PKF is mathematically well-defined against the order-1 DR, but
#' the user should be aware that the constraint bounds are expressed in
#' deviation units at the linearisation point, not in levels.
#'
#' @param model dynhr_mod
#' @noRd
obc_assert_linear <- function(model) {
  if (!isTRUE(model$model_options$linear)) {
    warning(
      "OBC/PKF is being applied to a nonlinear model (model; block detected). ",
      "The piecewise-linear Kalman filter will use the first-order (linear) ",
      "decision rules as its approximation.  OBC bounds should be expressed in ",
      "deviation units at the linearisation point.  For a fully global solution ",
      "with large constraint violations, consider a sparse-grid or projection method.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# =============================================================================
# Tag parsing
# =============================================================================

#' Parse MCP equation tags from a dynhr_mod object
#'
#' Reads `model$equations[[i]]$tag`, which the parser already populates from
#' \[ mcp = '...' \] annotations in the .mod/.txt file. Extracts the constrained
#' variable name, direction, and numeric bound from each tag.
#'
#' Tag format (from Dynare): mcp = 'VAR OP BOUND'
#'   e.g. "mcp = 'r > -0.00401606'"  (ZLB lower bound on r)
#'        "mcp = 'y < 0.01'"         (capacity upper bound on y)
#'
#' Note: the constrained variable (VAR) is the one that gets the box constraint
#' in PATH. It is NOT necessarily the LHS variable of the tagged equation.
#' E.g., \[ mcp = 'y < 0.01' \] on "tau_y = 0" pins y, not tau_y.
#'
#' @param model dynhr_mod (already parsed by parse_mod)
#' @return List of OBC spec lists, each with:
#'   $eq_idx   -- equation index in model$equations (1-based)
#'   $var_idx  -- index of the constrained variable in model$var_names (1-based)
#'   $var_name -- name of the constrained variable (for diagnostics)
#'   $op       -- ">" (lower bound) or "<" (upper bound)
#'   $bound    -- numeric bound value (in model deviation units)
#' Stops if the model has no MCP tags, or if a tag references an unknown variable.
#' @export
obc_parse_tags <- function(model) {
  eqs   <- model$equations
  specs <- list()

  for (i in seq_along(eqs)) {
    # Prefer tag_raw (full bracket content: name='...', mcp='...') over tag
    # (which only stores the name= attribute).  Fall back to tag for models
    # whose parser does not populate tag_raw.
    tag <- eqs[[i]]$tag_raw %||% eqs[[i]]$tag
    if (is.null(tag) || is.na(tag) || !nzchar(trimws(tag))) next

    # Match:  mcp = 'VAR OP BOUND'  (single or double quotes, flexible spacing)
    # VAR:    word characters (alphanumeric + underscore)
    # OP:     > >= < <=
    # BOUND:  numeric literal  (e.g. -0.004)  OR  parameter name  (e.g. zlb_floor)
    m <- regmatches(tag, regexec(
      paste0("mcp\\s*=\\s*['\"]\\s*(\\w+)\\s*([<>]=?)\\s*",
             "([+-]?(?:[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?|\\w+))\\s*['\"]"),
      tag, perl = TRUE
    ))[[1]]

    if (length(m) < 4) next   # tag present but not an mcp tag

    var_name  <- m[2]
    op_raw    <- m[3]
    bound_str <- m[4]

    # Evaluate bound: literal number or parameter name
    bound <- suppressWarnings(as.numeric(bound_str))
    if (is.na(bound)) {
      # Try to look up as a parameter name
      pv <- model$param_values[[bound_str]]
      if (!is.null(pv) && is.finite(pv)) {
        bound <- as.numeric(pv)
      } else {
        stop(sprintf(
          "MCP bound '%s' on equation %d is neither a numeric literal nor a known parameter.\n",
          bound_str, i
        ))
      }
    }

    op <- if (startsWith(op_raw, ">")) ">" else "<"

    var_idx <- match(var_name, model$var_names)
    if (is.na(var_idx)) {
      stop(sprintf(
        "MCP tag on equation %d references variable '%s', which is not in var_names.\n",
        i, var_name
      ))
    }

    specs <- c(specs, list(list(
      eq_idx   = i,
      var_idx  = var_idx,
      var_name = var_name,
      op       = op,
      bound    = bound
    )))
  }

  if (length(specs) == 0L) {
    stop(
      "No MCP equation tags found in the model. ",
      "Add [ mcp = 'var > bound' ] or [ mcp = 'var < bound' ] tags to the ",
      "relevant equations in the .mod/.txt file, then re-parse."
    )
  }

  specs
}


#' Resolve raw occbin_constraints block specs into the internal OBC spec form
#'
#' Takes the list produced by parse_occbin_constraints_block() (which contains
#' var_name, op, bound, bound_expr, name, relax_str, bind_eqs) and resolves
#' eq_idx, var_idx, and any parameter-reference bound values.
#'
#' eq_idx resolution: finds the equation whose LHS (as an AST string) matches
#' var_name. For linear models the LHS is always a bare variable name.
#'
#' When \code{ss} and \code{params} are supplied (both non-NULL), expression
#' bounds such as \code{PHI*steady_state(iv)} are evaluated: each
#' \code{steady_state(X)} is replaced with the numeric value \code{ss[X]},
#' then the resulting expression is eval'd in an environment containing all
#' model parameters.  This mirrors the resolution logic in
#' \code{.occbin_eval_bind_condition} and is required for the
#' \code{method="pwlinear"} OccBin path.
#'
#' @param raw_specs  List from parse_occbin_constraints_block()
#' @param model      dynhr_mod (already parsed, with equations and var_names)
#' @param ss         Named numeric SS vector (optional; needed for expression bounds)
#' @param params     Named numeric parameter vector (optional; needed for expression bounds)
#' @return List of resolved OBC specs, same shape as obc_parse_tags() output
#'   plus extra fields: $name, $bind_eqs, $relax_str, $source = "occbin_block"
#' @noRd
obc_resolve_block_specs <- function(raw_specs, model, ss = NULL, params = NULL) {
  if (length(raw_specs) == 0L) return(list())

  eqs       <- model$equations
  var_names <- model$var_names
  specs     <- list()

  for (rs in raw_specs) {
    # var_idx: look up in model$var_names
    var_idx <- match(rs$var_name, var_names)
    if (is.na(var_idx)) {
      stop(sprintf(
        "occbin_constraints: bind clause references variable '%s', which is not in var_names.",
        rs$var_name
      ))
    }

    # eq_idx: find equation whose LHS is the constrained variable.
    # ast_to_string() renders the LHS AST; for linear models it is just the
    # bare variable name (no lag suffix).
    eq_idx <- NA_integer_
    for (i in seq_along(eqs)) {
      lhs_str <- ast_to_string(eqs[[i]]$lhs)
      if (trimws(lhs_str) == rs$var_name) {
        eq_idx <- i
        break
      }
    }

    if (is.na(eq_idx)) {
      stop(sprintf(
        paste0("occbin_constraints: cannot find equation with LHS '%s' for ",
               "constraint '%s'. Ensure the model has an equation whose ",
               "left-hand side is exactly '%s'."),
        rs$var_name, rs$name, rs$var_name
      ))
    }

    # Resolve bound value: numeric literal, parameter reference, or full expression
    bound_val <- rs$bound
    if (is.na(bound_val) && is.character(rs$bound_expr) && nzchar(rs$bound_expr)) {
      expr <- rs$bound_expr

      # Try 1: bare parameter name (fast path, no eval needed)
      param_val <- model$param_values[expr]
      if (!is.na(param_val)) {
        bound_val <- as.numeric(param_val)
      } else if (!is.null(model$params) && expr %in% names(model$params)) {
        bound_val <- as.numeric(model$params[[expr]])
      }

      # Try 2: full expression with steady_state(X) calls (needs ss + params)
      if (is.na(bound_val) && !is.null(ss) && !is.null(params)) {
        # Build evaluation environment: all model parameters
        eval_env <- as.list(params)

        # Substitute steady_state(X) → numeric value from ss
        caps <- regmatches(expr,
                           gregexpr("steady_state\\s*\\(\\s*(\\w+)\\s*\\)",
                                    expr, perl = TRUE))[[1]]
        for (cap in caps) {
          vname  <- sub("steady_state\\s*\\(\\s*(\\w+)\\s*\\)", "\\1", cap, perl = TRUE)
          vi_ss  <- match(vname, names(ss))
          if (is.na(vi_ss)) vi_ss <- match(vname, var_names)
          ss_val <- if (!is.na(vi_ss)) ss[vi_ss] else NA_real_
          if (is.na(ss_val)) {
            warning(sprintf(
              "obc_resolve_block_specs: steady_state(%s) not found in ss; ",
              vname), "bound expression may be NA.")
          }
          expr <- gsub(cap, as.character(ss_val), expr, fixed = TRUE)
        }

        # Evaluate the expression
        bound_val <- tryCatch(
          as.numeric(eval(parse(text = expr), envir = eval_env)),
          error = function(e) NA_real_
        )
      }

      if (is.na(bound_val)) {
        stop(sprintf(
          paste0("occbin_constraints: cannot resolve bound expression '%s' for ",
                 "constraint '%s'. Bound must be a numeric literal, a parameter name, ",
                 "or an expression involving parameters and steady_state(X) calls. ",
                 "For expression bounds, supply ss and params to obc_resolve_block_specs."),
          rs$bound_expr, rs$name))
      }
    }

    specs <- c(specs, list(list(
      eq_idx    = eq_idx,
      var_idx   = var_idx,
      var_name  = rs$var_name,
      op        = rs$op,
      bound     = bound_val,
      name      = rs$name,
      bind_eqs  = rs$bind_eqs,
      relax_str = rs$relax_str,
      source    = "occbin_block"
    )))
  }

  specs
}


#' Collect all OBC specs from a model (MCP tags + occbin_constraints block)
#'
#' This is the primary entry point for the OBC pipeline. It merges specs from
#' both [mcp] equation tags (via obc_parse_tags) and the occbin_constraints
#' block (via obc_resolve_block_specs), deduplicating by eq_idx if both
#' sources declare the same equation.
#'
#' Backward-compatible: models with only [mcp] tags work exactly as before.
#' Models with only an occbin_constraints block also work. Mixed models get
#' the union, with [mcp] tags taking precedence on duplicates.
#'
#' @param model dynhr_mod
#' @return List of OBC spec lists (same shape as obc_parse_tags() output)
#' @noRd
obc_collect_specs <- function(model) {
  has_tags <- any(vapply(model$equations,
                         function(e) !is.na(e$tag) && nzchar(trimws(e$tag)),
                         logical(1)))
  has_block <- length(model$occbin_constraints) > 0L

  if (!has_tags && !has_block) {
    stop(
      "No OBC constraints found. Add [ mcp = 'var > bound' ] tags to model ",
      "equations, or add an occbin_constraints block to the .mod file."
    )
  }

  tag_specs   <- if (has_tags)  obc_parse_tags(model) else list()
  block_specs <- if (has_block) obc_resolve_block_specs(model$occbin_constraints, model) else list()

  # Dedup: if both sources reference the same eq_idx, tag wins.
  if (length(tag_specs) > 0L && length(block_specs) > 0L) {
    tag_eq_idxs <- vapply(tag_specs, `[[`, integer(1), "eq_idx")
    block_specs <- Filter(function(s) !(s$eq_idx %in% tag_eq_idxs), block_specs)
  }

  c(tag_specs, block_specs)
}
