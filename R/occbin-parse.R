## R/occbin-parse.R
## --------------------------------------------------------------------------
## OccBin bind/relax tag parsing for nonlinear DSGE models.
##
## Provides:
##   occbin_parse_bind_relax()   -- extract bind/relax constraints from eq tags
##   occbin_build_regime_map()   -- build equation-to-regime mapping table
##
## Dynare's OccBin for nonlinear models uses equation pairs annotated with
## bind='NAME' and relax='NAME' attributes in the equation tags.  E.g.:
##
##   model;
##     [name='Euler', bind='INEG']  C^(-sigma) = beta * R * ... (with adj cost)
##     [name='Euler', relax='INEG'] C^(-sigma) = beta * R * ... (standard Euler)
##     ...
##   end;
##   occbin_constraints;
##     name 'INEG'; bind log_Invest - log(steady_state(Invest)) < -1e-6;
##   end;
##
## Each constraint defines a "regime bit".  A model with k constraints has
## 2^k regimes.  For each regime, exactly n_endo equations are selected:
## one per endogenous variable.
##
## The occbin_constraints block (when present) provides the bound condition.
## --------------------------------------------------------------------------


# =============================================================================
# Parse bind/relax tags from equations
# =============================================================================

#' Parse bind/relax equation tags from a dynhr_mod
#'
#' Scans model equations for bind='...' and relax='...' attributes in their
#' equation tags (stored in eq$tag_raw).  Groups equations into bind/relax
#' variants for each constraint and identifies neutral equations (no regime
#' annotation).
#'
#' The detection identifies regime tag attributes using the tag_raw field:
#'   tag_raw = "name='Euler', bind='INEG'"
#'   tag_raw = "name='Euler', relax='INEG'"
#'
#' Equations without bind/relax attributes are "neutral" — used in every regime.
#'
#' @param model dynhr_mod (from parse_mod)
#' @return List with:
#'   $constraints  -- list of RegimeConstraint objects, each with:
#'     $name       - constraint name (e.g., "ZLB", "INEG")
#'     $eq_bind    - integer vector of equation indices for bind variant
#'     $eq_relax   - integer vector of equation indices for relax variant
#'     $base_name  - the shared equation name (from name='...' attribute)
#'     $var_name   - constrained variable name (LHS of the equation)
#'     $condition  - "less" or "greater" (from occbin_constraints block)
#'     $bound_expr - raw bound expression string
#'   $neutral_eqs  -- integer vector of equation indices with no regime tag
#'   $n_constraints -- integer: number of distinct constraints (k)
#' @noRd
occbin_parse_bind_relax <- function(model) {
  eqs <- model$equations
  n_eq <- length(eqs)

  # Collect all bind/relax annotations
  # Each entry: list(name, kind ("bind"|"relax"), eq_idx, base_name, var_name)
  annotations <- list()

  for (i in seq_len(n_eq)) {
    tag_raw <- eqs[[i]]$tag_raw
    if (is.na(tag_raw) || !nzchar(trimws(tag_raw))) next

    # Extract base name (from name='...')
    base_name <- NA_character_
    nm <- regmatches(tag_raw,
                     regexec("name\\s*=\\s*'([^']*)'",
                             tag_raw, perl = TRUE))[[1]]
    if (length(nm) > 1) base_name <- nm[2]

    # Extract bind='...' or relax='...' attributes
    # Format: bind='CONSTRAINT_NAME' or relax='CONSTRAINT_NAME'
    # The constraint name allows single or double quotes
    bind_match <- regmatches(tag_raw,
                             regexec("\\bbind\\s*=\\s*'([^']*)'",
                                     tag_raw, perl = TRUE))[[1]]
    relax_match <- regmatches(tag_raw,
                              regexec("\\brelax\\s*=\\s*'([^']*)'",
                                      tag_raw, perl = TRUE))[[1]]

    # Determine LHS variable name for grouping
    lhs_str <- ast_to_string(eqs[[i]]$lhs)

    if (length(bind_match) > 1) {
      annotations <- c(annotations, list(list(
        name      = bind_match[2],
        kind      = "bind",
        eq_idx    = i,
        base_name = base_name,
        var_name  = lhs_str
      )))
    }
    if (length(relax_match) > 1) {
      annotations <- c(annotations, list(list(
        name      = relax_match[2],
        kind      = "relax",
        eq_idx    = i,
        base_name = base_name,
        var_name  = lhs_str
      )))
    }
  }

  if (length(annotations) == 0L) {
    return(list(
      constraints   = list(),
      neutral_eqs   = seq_len(n_eq),
      n_constraints = 0L
    ))
  }

  # Group annotations by constraint name
  constraint_names <- unique(vapply(annotations, `[[`, character(1), "name"))
  constraints <- vector("list", length(constraint_names))
  names(constraints) <- constraint_names

  for (cn in constraint_names) {
    ann_cn <- annotations[vapply(annotations, function(a) a$name == cn, logical(1))]
    bind_eqs   <- sort(vapply(ann_cn[vapply(ann_cn, `[[`, character(1), "kind") == "bind"],
                             `[[`, integer(1), "eq_idx"))
    relax_eqs  <- sort(vapply(ann_cn[vapply(ann_cn, `[[`, character(1), "kind") == "relax"],
                              `[[`, integer(1), "eq_idx"))

    # Determine LHS variable name (should be same for both variants)
    var_name <- ann_cn[[1]]$var_name
    base_name <- ann_cn[[1]]$base_name

    # Extract condition and bound from occbin_constraints block if present
    # Derive condition from the constraint op:
    #   op = ">" (lower bound, e.g. r > r_lb) → binds when r < bound → "less"
    #   op = "<" (upper bound, e.g. r < r_lb) → binds when r > bound → "greater"
    condition    <- "greater"
    bound_expr   <- NA_character_
    # bind_var_name: variable whose VALUE is compared to the bound.
    # This is the variable named in the occbin_constraints block bind clause
    # (e.g. 'iv' in 'bind iv < PHI*steady_state(iv)').  It may DIFFER from
    # the LHS of the bind equation (e.g. the bind equation is 'iv = PHI*...'
    # which has LHS 'iv', but the relax equation has LHS 'lam = 0').
    # We override with the occbin_constraints block variable when available.
    bind_var_name <- var_name   # fallback: first annotation's LHS

    if (length(model$occbin_constraints) > 0L) {
      for (oc in model$occbin_constraints) {
        if (identical(oc$name, cn)) {
          condition     <- if (identical(oc$op, ">")) "less" else "greater"
          bound_expr    <- oc$bound_expr
          # oc$var_name is the authoritative bind variable from the
          # occbin_constraints block (e.g. 'iv'), which is correct even when
          # the bind/relax equation pair uses a different variable for its LHS.
          if (!is.null(oc$var_name) && nzchar(oc$var_name))
            bind_var_name <- oc$var_name
          break
        }
      }
    }

    constraints[[cn]] <- list(
      name          = cn,
      eq_bind       = bind_eqs,
      eq_relax      = relax_eqs,
      base_name     = base_name,
      var_name      = bind_var_name,  # bind-condition variable (from occbin_constraints)
      eq_var_name   = var_name,       # first equation LHS (kept for backward compat)
      condition     = condition,
      bound_expr    = bound_expr
    )
  }

  # Neutral equations: those not referenced in any bind/relax annotation
  annotated_eqs <- unique(vapply(annotations, `[[`, integer(1), "eq_idx"))
  neutral_eqs <- setdiff(seq_len(n_eq), annotated_eqs)

  list(
    constraints    = constraints,
    neutral_eqs    = neutral_eqs,
    n_constraints  = length(constraint_names)
  )
}


# =============================================================================
# Build regime-to-equation mapping
# =============================================================================

#' Build regime-to-equation selection map
#'
#' For each possible regime (bitmask over constraints), determines which
#' equations are active.  A regime bitmask of 0 means all constraints are
#' slack (relax variants used).  Bit j = 1 means constraint j binds.
#'
#' @param parse_result Output of occbin_parse_bind_relax()
#' @return List of length 2^n_constraints, each entry a list with:
#'   $regime     -- integer bitmask
#'   $eq_indices -- integer vector of active equation indices (length n_endo)
#' @noRd
occbin_build_regime_map <- function(parse_result) {
  constraints   <- parse_result$constraints
  neutral_eqs   <- parse_result$neutral_eqs
  n_constraints <- parse_result$n_constraints
  n_regimes     <- 2L^n_constraints

  # Map constraint name to index
  if (n_constraints > 0L) {
    cn_to_idx <- setNames(seq_len(n_constraints), names(constraints))
  } else {
    cn_to_idx <- integer(0)
  }

  regime_map <- vector("list", n_regimes)

  for (r in seq_len(n_regimes) - 1L) {
    eq_set <- neutral_eqs  # start with neutral equations

    if (n_constraints > 0L) {
      for (cn in names(constraints)) {
        j <- cn_to_idx[cn]
        bind_bit <- bitwAnd(r, 2L^(j - 1L)) != 0L

        if (bind_bit) {
          eq_set <- c(eq_set, constraints[[cn]]$eq_bind)
        } else {
          eq_set <- c(eq_set, constraints[[cn]]$eq_relax)
        }
      }
    }

    regime_map[[r + 1L]] <- list(
      regime     = r,
      eq_indices = sort(unique(eq_set))
    )
  }

  regime_map
}


# =============================================================================
# Constraint condition evaluation
# =============================================================================


# =============================================================================
# Regime path from current constraint violations
# =============================================================================
