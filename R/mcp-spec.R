## R/mcp-spec.R
## --------------------------------------------------------------------------
## MCP (Mixed Complementarity Problem) constraint specification: parsing,
## validation, and constraint-to-variable mapping.
##
## Provides:
##   mcp_parse_tags()        -- parse standard Dynare [mcp = 'var OP bound'] tags
##   mcp_validate_specs()    -- cross-validate specs against model structure
##   mcp_constraint_map()    -- build per-variable constraint index map
##   mcp_specs_from_occbin() -- convert occbin bind/relax constraints to MCP specs
##
## MCP SPEC FIELDS
##   Each spec is a list with:
##     $var_name  -- constrained variable name
##     $var_idx   -- index in model$var_names / endo_names (1-based)
##     $eq_idx    -- original equation index (1-based)
##     $op        -- ">" (lower bound) or "<" (upper bound)
##     $bound     -- numeric bound value (NA_real_ when bound is an expression)
##     $bound_ast -- AST of the bound expression (NULL when bound is a literal)
##     $tag_type  -- "mcp" or "bind_relax" (source format)
##     $name      -- optional constraint name (from bind/relax tags)
##
## REFERENCES
##   Dynare MCP tag format: [mcp = 'VAR OP BOUND']
##     e.g. [mcp = 'r > -0.00401606']   (ZLB lower bound on r, numeric)
##          [mcp = 'y < 0.01']          (capacity upper bound on y, numeric)
##          [mcp = 'gam2 > i+1/beta-1'] (ZLB on gam2, expression bound)
##
##   The constrained variable (VAR) is the one that gets the box constraint
##   in PATH. It is NOT necessarily the LHS variable of the tagged equation.
##   E.g., [mcp = 'y < 0.01'] on "tau_y = 0" pins y, not tau_y.
##
##   Ferris, M. C., & Munson, T. S. (2000). "Complementarity problems in
##     GAMS and the PATH solver." Journal of Economic Dynamics and Control.
## --------------------------------------------------------------------------


# =============================================================================
# Parse MCP tags from dynhr_mod
# =============================================================================

#' Parse MCP equation tags from a dynhr_mod object
#'
#' Reads \code{model$equations[[i]]$tag} entries and extracts MCP constraint
#' specifications from \code{[mcp = 'VAR OP BOUND']} annotations.  Returns
#' a uniform list of MCP spec objects that can be used by \code{mcp_solve_path()}.
#'
#' The tag parser handles both single and double quotes, optional spacing,
#' and scientific notation in numeric bound values.  When the bound RHS is
#' an expression (e.g. \code{i+1/beta-1}), it is parsed into an AST and
#' stored in the \code{bound_ast} field; \code{bound} is \code{NA_real_}.
#'
#' If no MCP tags are found, attempts to fall back to OccBin bind/relax
#' tags (via \code{mcp_specs_from_occbin()}) for backward compatibility
#' with existing nonlinear model annotations.
#'
#' Models containing bare \code{0 = min(VAR, EXPR)} or \code{0 = max(VAR, EXPR)}
#' equations raise a clear error directing the user to annotate with
#' \code{[mcp = '...']} instead.
#'
#' @param model dynhr_mod (from \code{\link{parse_mod}})
#' @param verbose Logical: print diagnostic messages (default FALSE)
#' @return List of MCP spec lists, each with:
#'   \item{var_name}{Character: constrained variable name}
#'   \item{var_idx}{Integer: index in \code{model$var_names} (1-based)}
#'   \item{eq_idx}{Integer: equation index in \code{model$equations} (1-based)}
#'   \item{op}{Character: \code{">"} for lower bound, \code{"<"} for upper bound}
#'   \item{bound}{Numeric: bound literal, or \code{NA_real_} for expression bounds}
#'   \item{bound_ast}{AST list or \code{NULL}: parsed expression when bound is not
#'     a numeric literal}
#'   \item{tag_type}{Character: \code{"mcp"} or \code{"bind_relax"}}
#'   \item{name}{Character or NA: optional constraint name (from bind/relax)}
#' @export
#'
#' @examples
#' \dontrun{
#' model <- parse_mod("nk_zlb_dynare.mod")
#' specs <- mcp_parse_tags(model)
#' # specs[[1]]$var_name == "r", specs[[1]]$op == ">", specs[[1]]$bound == -0.01
#'
#' # Expression bound: [mcp = 'gam2 > i+1/beta-1']
#' # specs[[1]]$var_name == "gam2", specs[[1]]$op == ">",
#' # is.na(specs[[1]]$bound), !is.null(specs[[1]]$bound_ast)
#' }
mcp_parse_tags <- function(model, verbose = FALSE) {
  eqs   <- model[["equations"]]
  specs <- list()
  n_mcp <- 0L

  # Collect var/param names for expression parsing
  var_names_all   <- c(model$var_names, model$exo_names)
  param_names_all <- names(model$param_values)

  # I4: detect bare 0 = min(VAR, EXPR) / 0 = max(VAR, EXPR) equations that
  # are missing an [mcp=...] tag.  These are silently ignored by the solver,
  # which is wrong.  Emit a clear error so the user knows to annotate.
  for (i in seq_along(eqs)) {
    eq      <- eqs[[i]]
    tag_raw <- eq[["tag_raw"]]
    has_mcp_tag <- !is.null(tag_raw) && !is.na(tag_raw) &&
                   grepl("mcp\\s*=", tag_raw, perl = TRUE)
    has_bind_tag <- !is.null(tag_raw) && !is.na(tag_raw) &&
                    grepl("\\bbind\\s*=", tag_raw, perl = TRUE)
    if (has_mcp_tag || has_bind_tag) next  # already annotated

    lhs <- eq[["lhs"]]
    rhs <- eq[["rhs"]]
    # Pattern: 0 = min(VAR, EXPR) or 0 = max(VAR, EXPR)
    # LHS is a number node with value 0, RHS is a funcall of min/max
    lhs_is_zero <- !is.null(lhs) && identical(lhs$type, "number") &&
                   isTRUE(lhs$value == 0)
    rhs_is_minmax <- !is.null(rhs) && identical(rhs$type, "funcall") &&
                     rhs$name %in% c("min", "max")
    if (lhs_is_zero && rhs_is_minmax) {
      fn <- rhs$name
      stop(sprintf(paste0(
        "Equation %d contains `0 = %s(...)` without an [mcp=...] tag.\n",
        "  This equation is silently ignored by the OBC solver.\n",
        "  To register it as an OBC constraint, annotate it with an MCP tag, e.g.:\n",
        "    [mcp = 'VAR > BOUND']  0 = %s(VAR, EXPR)\n",
        "  or use an occbin_constraints block with bind/relax equation pairs.\n",
        "  See ?mcp_parse_tags for details."),
        i, fn, fn))
    }
  }

  for (i in seq_along(eqs)) {
    # Try both tag and tag_raw.  The 'tag' field stores the name= value
    # (e.g. "LOM susceptible"), while 'tag_raw' stores the full bracket
    # content (e.g. "mcp='S>0',name='LOM susceptible'") which may contain
    # the mcp= spec.
    tag     <- eqs[[i]][["tag"]]
    tag_raw <- eqs[[i]][["tag_raw"]]
    if ((is.null(tag) || is.na(tag) || !nzchar(trimws(tag))) &&
        (is.null(tag_raw) || is.na(tag_raw) || !nzchar(trimws(tag_raw)))) next

    # Try matching mcp= in tag_raw first (it has the full spec), then tag
    search_str <- if (!is.null(tag_raw) && !is.na(tag_raw) && nzchar(trimws(tag_raw))) tag_raw else tag

    # Two-pass match:
    #   Pass 1: try numeric bound  mcp = 'VAR OP NUMBER'
    #   Pass 2: try expression bound  mcp = 'VAR OP EXPR'
    # In both cases capture: (1) VAR, (2) OP, (3) RHS text (everything inside quotes after OP)
    m <- regmatches(search_str, regexec(
      "mcp\\s*=\\s*(['\"])\\s*(\\w+)\\s*([<>]=?)\\s*(.+?)\\s*\\1",
      search_str, perl = TRUE
    ))[[1]]

    if (length(m) < 5) next

    var_name  <- m[3]
    op_raw    <- m[4]
    bound_str <- trimws(m[5])

    op <- if (startsWith(op_raw, ">")) ">" else "<"

    var_idx <- match(var_name, model$var_names)
    if (is.na(var_idx)) {
      stop(sprintf(
        "MCP tag on equation %d references variable '%s', which is not in var_names.",
        i, var_name))
    }

    # Try to parse bound_str as a numeric literal first (backward compat)
    bound     <- suppressWarnings(as.numeric(bound_str))
    bound_ast <- NULL

    if (is.na(bound)) {
      # Not a plain number — parse as an expression AST
      bound_ast <- tryCatch(
        parse_expression(bound_str,
                         var_names   = var_names_all,
                         param_names = param_names_all),
        error = function(e) {
          stop(sprintf(
            "MCP tag on equation %d: cannot parse bound expression '%s': %s",
            i, bound_str, conditionMessage(e)))
        }
      )
    }

    n_mcp <- n_mcp + 1L
    specs[[n_mcp]] <- list(
      var_name  = var_name,
      var_idx   = var_idx,
      eq_idx    = i,
      op        = op,
      bound     = bound,
      bound_ast = bound_ast,
      tag_type  = "mcp",
      name      = NA_character_
    )

    if (verbose) {
      bound_repr <- if (is.finite(bound)) sprintf("%.6g", bound) else
                    sprintf("<expr: %s>", ast_to_string(bound_ast))
      message(sprintf("  [mcp] eq %d: %s %s %s", i, var_name, op, bound_repr))
    }
  }

  # Fall back to OccBin bind/relax if no MCP tags found
  if (n_mcp == 0L) {
    if (verbose) message("  No MCP tags found; trying OccBin bind/relax fallback...")
    specs <- mcp_specs_from_occbin(model, verbose = verbose)
  }

  if (length(specs) == 0L) {
    if (verbose) message("  No MCP constraints found in model.")
  }

  specs
}


# =============================================================================
# Convert OccBin bind/relax constraints to MCP specs
# =============================================================================

#' Convert OccBin bind/relax constraints to MCP spec format
#'
#' Parses Dynare's bind/relax equation tags (used by OccBin for nonlinear
#' models) and converts them to the standard MCP spec format.  This enables
#' \code{mcp_solve_path()} to handle models annotated with \code{bind='...'}
#' and \code{relax='...'} attributes alongside standard MCP tags.
#'
#' The conversion identifies the constrained variable from the equation LHS,
#' and extracts the bound condition from the \code{occbin_constraints} block
#' (if present) or from the equation name.
#'
#' @param model dynhr_mod (from \code{\link{parse_mod}})
#' @param verbose Logical: print diagnostic messages (default FALSE)
#' @return List of MCP spec lists (same format as \code{mcp_parse_tags()});
#'   empty list if no bind/relax constraints found.
#' @noRd
mcp_specs_from_occbin <- function(model, verbose = FALSE) {
  eqs <- model$equations
  n_eq <- length(eqs)

  # Collect bind annotations: name -> list(eq_idx, var_name)
  bind_eqs   <- list()
  relax_eqs  <- list()

  for (i in seq_len(n_eq)) {
    tag_raw <- eqs[[i]]$tag_raw
    if (is.na(tag_raw) || !nzchar(trimws(tag_raw))) next

    # Extract bind='NAME'
    bm <- regmatches(tag_raw,
      regexec("\\bbind\\s*=\\s*'([^']*)'", tag_raw, perl = TRUE))[[1]]
    if (length(bm) > 1) {
      cn <- bm[2]
      lhs_str <- ast_to_string(eqs[[i]]$lhs)
      bind_eqs[[cn]] <- list(eq_idx = i, var_name = lhs_str)
    }

    # Extract relax='NAME'
    rm <- regmatches(tag_raw,
      regexec("\\brelax\\s*=\\s*'([^']*)'", tag_raw, perl = TRUE))[[1]]
    if (length(rm) > 1) {
      cn <- rm[2]
      lhs_str <- ast_to_string(eqs[[i]]$lhs)
      relax_eqs[[cn]] <- list(eq_idx = i, var_name = lhs_str)
    }
  }

  if (length(bind_eqs) == 0L) return(list())

  # Build specs from bind constraints
  constraint_names <- names(bind_eqs)
  specs <- vector("list", length(constraint_names))

  for (j in seq_along(constraint_names)) {
    cn <- constraint_names[[j]]
    be <- bind_eqs[[cn]]

    # Determine bound from occbin_constraints block if available
    bound <- NA_real_
    op    <- ">"  # default: lower bound
    bound_expr <- NA_character_

    if (length(model$occbin_constraints) > 0L) {
      for (oc in model$occbin_constraints) {
        if (identical(oc$name, cn)) {
          # Derive op from condition field
          # condition = "less" (binds when variable is low) -> lower bound ">"
          # condition = "greater" (binds when variable is high) -> upper bound "<"
          if (identical(oc$condition, "less")) {
            op <- ">"
          } else if (identical(oc$condition, "greater")) {
            op <- "<"
          } else {
            op <- if (identical(oc$op, ">")) ">" else "<"
          }
          bound_expr <- oc$bound_expr

          # Try to evaluate bound_expr as a numeric
          if (!is.na(bound_expr) && nzchar(bound_expr)) {
            bound_val <- tryCatch(
              eval(parse(text = bound_expr),
                   envir = as.list(model$param_values)),
              error = function(e) NA_real_
            )
            if (is.finite(bound_val)) bound <- bound_val
          }
          break
        }
      }
    }

    var_idx <- match(be$var_name, model$var_names)
    if (is.na(var_idx)) next

    # Use the RELAX equation index as eq_idx for the FB residual.
    #
    # The Fischer-Burmeister complementarity is φ(a, b) = 0 where
    #   a = x - bound  (slack for the lower bound)
    #   b = F_{eq}(y)  (residual of the model equation at eq_idx)
    #
    # For a lower bound on x (e.g. r >= r_lb):
    #   - When the constraint is SLACK (r > r_lb): the RELAX equation should
    #     hold (e.g. the standard Taylor rule F_relax = 0).  That makes
    #     a > 0, b ≈ 0  →  φ ≈ 0 ✓
    #   - When the constraint BINDS (r = r_lb): a = 0, b ≠ 0  →  φ = 0 ✓
    #
    # If we instead used the BIND equation (e.g. r = r_lb → F_bind = r - r_lb),
    # then b = r - r_lb = a always, so φ(a, a) = (2−√2)·a = 0 forces a = 0,
    # i.e. r = r_lb for all periods — the explosive all-bound path (NEW-MCP1).
    re <- relax_eqs[[cn]]
    eq_idx_use <- if (!is.null(re)) re$eq_idx else be$eq_idx

    specs[[j]] <- list(
      var_name = be$var_name,
      var_idx  = var_idx,
      eq_idx   = eq_idx_use,
      op       = op,
      bound    = bound,
      tag_type = "bind_relax",
      name     = cn
    )

    if (verbose) {
      message(sprintf("  [bind/relax] '%s': %s %s %s",
                      cn, be$var_name, op,
                      if (is.finite(bound)) sprintf("%.6g", bound) else "?"))
    }
  }

  # Remove NULL entries (from failed var_idx lookups)
  specs[!vapply(specs, is.null, logical(1))]
}


# =============================================================================
# Validate MCP specs against model
# =============================================================================

#' Validate MCP constraint specifications
#'
#' Cross-checks a list of MCP specs against the model structure:
#' - All referenced variables exist in endo_names / var_names
#' - All equation indices are valid
#' - Bound values are finite and sensible
#' - Detects conflicting constraints (same variable with both > and <)
#'
#' @param model dynhr_mod (from \code{\link{parse_mod}})
#' @param mcp_specs List of MCP spec lists (from \code{mcp_parse_tags()})
#' @param verbose Logical: print validation messages (default FALSE)
#' @return Invisibly returns TRUE if all specs are valid; stops with an
#'   error message describing the first invalid spec.
#' @export
mcp_validate_specs <- function(model, mcp_specs, verbose = FALSE) {
  if (length(mcp_specs) == 0L) {
    if (verbose) message("  [mcp_validate] No specs to validate.")
    return(invisible(TRUE))
  }

  n_endo <- length(model$var_names)
  n_eq   <- length(model$equations)

  # Track which variables have lower/upper bounds for conflict detection
  var_lower <- character(0)
  var_upper <- character(0)

  for (i in seq_along(mcp_specs)) {
    sp <- mcp_specs[[i]]

    # Check required fields
    if (is.null(sp$var_name) || is.null(sp$var_idx)) {
      stop(sprintf("mcp_specs[[%d]]: missing var_name or var_idx.", i))
    }
    if (is.null(sp$eq_idx)) {
      stop(sprintf("mcp_specs[[%d]]: missing eq_idx for var '%s'.", i, sp$var_name))
    }
    if (is.null(sp$op) || !sp$op %in% c(">", "<")) {
      stop(sprintf("mcp_specs[[%d]]: op must be '>' (lower) or '<' (upper), got '%s'.",
                   i, sp$op %||% "NULL"))
    }
    if (is.null(sp$bound) || (!is.finite(sp$bound) && is.null(sp$bound_ast))) {
      stop(sprintf(
        "mcp_specs[[%d]]: bound must be finite numeric or an expression AST, got %s.",
        i, deparse(sp$bound)))
    }

    # Check variable index
    if (sp$var_idx < 1 || sp$var_idx > n_endo) {
      stop(sprintf("mcp_specs[[%d]]: var_idx %d out of range [1, %d].",
                   i, sp$var_idx, n_endo))
    }

    # Check variable name consistency
    expected_name <- model$var_names[sp$var_idx]
    if (sp$var_name != expected_name) {
      warning(sprintf(
        "mcp_specs[[%d]]: var_name '%s' != var_names[%d] = '%s'. Using var_idx.",
        i, sp$var_name, sp$var_idx, expected_name))
      sp$var_name <- expected_name
    }

    # Check equation index
    if (sp$eq_idx < 1 || sp$eq_idx > n_eq) {
      stop(sprintf("mcp_specs[[%d]]: eq_idx %d out of range [1, %d] for var '%s'.",
                   i, sp$eq_idx, n_eq, sp$var_name))
    }

    # Detect conflicting bounds
    if (sp$op == ">") {
      if (sp$var_name %in% var_lower) {
        bound_repr <- if (is.finite(sp$bound)) sprintf("%.6g", sp$bound) else "<expr>"
        warning(sprintf(
          "mcp_specs[[%d]]: var '%s' has multiple lower bounds (%s).",
          i, sp$var_name, bound_repr))
      }
      var_lower <- c(var_lower, sp$var_name)
    } else {
      if (sp$var_name %in% var_upper) {
        warning(sprintf(
          "mcp_specs[[%d]]: var '%s' has multiple upper bounds.", i, sp$var_name))
      }
      var_upper <- c(var_upper, sp$var_name)
    }

    # Both bounds on same variable
    if (sp$var_name %in% var_lower && sp$var_name %in% var_upper) {
      if (verbose) message(sprintf(
        "  [mcp_validate] var '%s' has both lower and upper bounds.", sp$var_name))
    }
  }

  if (verbose) {
    message(sprintf("  [mcp_validate] %d specs validated (%d vars with bounds).",
                    length(mcp_specs), length(unique(c(var_lower, var_upper)))))
  }

  invisible(TRUE)
}


# =============================================================================
# Build constraint-to-variable mapping
# =============================================================================

#' Build per-variable constraint index map
#'
#' For each constrained variable, identifies which MCP specs apply and whether
#' it has a lower bound, upper bound, or both.  Used by the solver to route
#' FB modifications to the correct Jacobian rows.
#'
#' @param mcp_specs List of MCP spec lists (from \code{mcp_parse_tags()})
#' @return A list with:
#'   \item{n_spec}{Integer: total number of MCP specs}
#'   \item{var_to_spec}{Named list: for each constrained variable, the list of
#'     spec indices that apply}
#'   \item{var_to_bound_type}{Named character: "lower", "upper", or "both"}
#'   \item{var_to_spec_idx}{Named integer: for each constrained variable, the
#'     FIRST spec index (for convenience when there's only one)}
#'   \item{eq_to_spec}{Integer vector (length = max eq_idx): maps equation index
#'     to spec index; 0 = unconstrained equation}
#'   \item{constrained_vars}{Character: names of all constrained variables}
#' @noRd
mcp_constraint_map <- function(mcp_specs) {
  n_spec <- length(mcp_specs)

  if (n_spec == 0L) {
    return(list(
      n_spec              = 0L,
      var_to_spec         = list(),
      var_to_bound_type   = character(0),
      var_to_spec_idx     = integer(0),
      eq_to_spec          = integer(0),
      constrained_vars    = character(0)
    ))
  }

  var_to_spec       <- list()
  var_to_bound_type <- character(0)

  for (j in seq_len(n_spec)) {
    sp  <- mcp_specs[[j]]
    vn  <- sp$var_name

    if (is.null(var_to_spec[[vn]])) {
      var_to_spec[[vn]] <- integer(0)
    }
    var_to_spec[[vn]] <- c(var_to_spec[[vn]], j)

    current_type <- var_to_bound_type[vn]
    if (is.na(current_type)) {
      var_to_bound_type[vn] <- sp$op
    } else {
      if (current_type != sp$op) {
        var_to_bound_type[vn] <- "both"
      }
    }
  }

  # Convenience: first spec index per variable
  var_to_spec_idx <- setNames(
    vapply(var_to_spec, `[`, integer(1), 1L),
    names(var_to_spec)
  )

  # Equation-to-spec mapping
  max_eq <- max(vapply(mcp_specs, `[[`, integer(1), "eq_idx"))
  eq_to_spec <- integer(max_eq)
  for (j in seq_len(n_spec)) {
    eq_to_spec[mcp_specs[[j]]$eq_idx] <- j
  }

  list(
    n_spec            = n_spec,
    var_to_spec       = var_to_spec,
    var_to_bound_type = var_to_bound_type,
    var_to_spec_idx   = var_to_spec_idx,
    eq_to_spec        = eq_to_spec,
    constrained_vars  = names(var_to_spec)
  )
}


# =============================================================================
# Bound value resolution
# =============================================================================

#' Resolve MCP bound values that may be expressed as parameter names
#'
#' Some models express bounds as parameter references (e.g., \code{r_lb})
#' rather than numeric literals.  This resolves them using the model's
#' parameter values and the OccBin constraints block.
#'
#' @param mcp_specs List of MCP spec lists (from \code{mcp_parse_tags()} or
#'   \code{mcp_specs_from_occbin()})
#' @param model dynhr_mod (for \code{occbin_constraints} block)
#' @param params Named numeric parameter vector
#' @return Modified MCP spec list with resolved \code{bound} fields.
#' @noRd
mcp_resolve_bounds <- function(mcp_specs, model, params) {
  if (length(mcp_specs) == 0L) return(mcp_specs)

  for (j in seq_along(mcp_specs)) {
    sp <- mcp_specs[[j]]

    # If bound is already numeric and finite, skip
    if (is.finite(sp$bound)) next

    # Try to resolve from occbin_constraints block
    if (sp$tag_type == "bind_relax" && !is.na(sp$name) &&
        length(model$occbin_constraints) > 0L) {
      for (oc in model$occbin_constraints) {
        if (identical(oc$name, sp$name) && !is.na(oc$bound_expr)) {
          bound_val <- tryCatch(
            eval(parse(text = oc$bound_expr), envir = as.list(params)),
            error = function(e) NA_real_
          )
          if (is.finite(bound_val)) {
            sp$bound <- bound_val
            mcp_specs[[j]] <- sp
          }
          break
        }
      }
    }
  }

  # Check for any still-unresolved bounds
  for (j in seq_along(mcp_specs)) {
    if (!is.finite(mcp_specs[[j]]$bound)) {
      warning(sprintf(
        "mcp_specs[[%d]] ('%s'): bound could not be resolved from model.",
        j, mcp_specs[[j]]$var_name))
    }
  }

  mcp_specs
}
