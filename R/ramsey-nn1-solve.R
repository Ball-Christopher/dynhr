
## R/ramsey-nn1-solve.R
## E4: (n, n+1) system solver. Modifies planner objective in-memory, solves via perturbation.

.nn1_solve <- function(model, compiled, ss, params,
                        taylor, nn1_objective, multipliers,
                        n, beta, verbose = FALSE, ...) {
  endo_names <- model$var_names
  n_endo <- length(endo_names)
  n_eq <- length(model$equations)

  if (verbose) cat("[nn1_solve] Building modified objective expression...
")
  obj_expr <- .nn1_build_objective_expr(nn1_objective[["coefficients"]], endo_names, ss, n)
  if (verbose) cat(sprintf("  Objective: %s ...
", substr(obj_expr, 1, 150)))

  modified_model <- model
  modified_model$planner_objective <- list(text = obj_expr)

  # Handle instruments: if n_endo > n_eq, add placeholder equations
  # to make the system square. The modified objective ensures optimality.
  instruments_added <- character(0)
  if (n_endo > n_eq) {
    # Identify instrument variables: endogenous vars that don't have a
    # dedicated equation. We use the compiled Jacobian: a variable that
    # appears in f_zero with a |coeff| >= 0.5 in exactly ONE row is
    # solved by that equation. Variables that DON'T match this pattern
    # but still have no LHS equation are instruments.
    sys <- extract_system_matrices(compiled, ss, params)
    J_zero <- sys$f_zero
    all_lhs_names <- sapply(model$equations, function(e) {
      if (e$lhs$type == "variable") e$lhs$name else NA_character_
    })
    lhs_candidates <- setdiff(model$var_names, all_lhs_names)
    # Among LHS candidates, find ones that appear in f_zero with a
    # near-unity coefficient in exactly one equation row.
    has_dedicated_eq <- logical(length(lhs_candidates))
    for (k in seq_along(lhs_candidates)) {
      v <- lhs_candidates[k]
      j <- which(model$var_names == v)
      if (length(j) == 0) next
      # Count rows where |coeff| >= 0.5
      strong_rows <- sum(abs(J_zero[, j]) >= 0.5)
      has_dedicated_eq[k] <- (strong_rows >= 1)
    }
    n_instruments <- n_endo - n_eq
    # Pick the instruments: LHS candidates WITHOUT a dedicated equation
    inst_candidates <- lhs_candidates[!has_dedicated_eq]
    # If too many or too few, fall back: pick variable that appears in
    # the MOST equations (least likely to have a single dedicated eq)
    if (length(inst_candidates) > n_instruments) {
      col_appearances <- sapply(inst_candidates, function(v) {
        j <- which(model$var_names == v)
        if (length(j) == 0) return(0L)
        sum(abs(J_zero[, j]) > 1e-10)
      })
      inst_candidates <- inst_candidates[order(col_appearances, decreasing = TRUE)[seq_len(n_instruments)]]
    } else if (length(inst_candidates) < n_instruments) {
      inst_candidates <- lhs_candidates[seq_len(n_instruments)]
    }
    instruments_added <- inst_candidates
    if (length(instruments_added) > 0) {
      if (verbose) {
        cat(sprintf("  Adding %d placeholder equation(s) for instrument(s): %s
",
                    length(instruments_added),
                    paste(instruments_added, collapse = ", ")))
      }
      # Build simple AST: instrument_name = 0
      # Also add a marker equation to signal to solve_perturbation that this
      # is a Ramsey/NN1 context (unit-root eigenvalues are forward-looking).
      for (inst_name in instruments_added) {
        placeholder_eq <- list(
          lhs = list(type = "variable", name = inst_name, lead_lag = 0L),
          rhs = list(type = "number", value = 0),
          type = "equation"
        )
        class(placeholder_eq) <- "ast_node"
        modified_model$equations[[length(modified_model$equations) + 1]] <- placeholder_eq
      }
      # Add a marker variable so the BK check can detect NN1 context
      modified_model$ramsey_context <- TRUE
    }
  }

  modified_compiled <- compile_model(modified_model, verbose = verbose)
  if (verbose) cat(sprintf("[nn1_solve] Running perturbation(order=%d)...
", n))
  dr <- solve_perturbation(modified_model, modified_compiled, ss, params,
                            order = as.integer(n), verbose = verbose, ...)
  bk_ok <- isTRUE(dr$bk_satisfied)
  if (verbose) cat(sprintf("  BK: %s
", if (bk_ok) "PASSED" else "FAILED"))
  list(dr = dr, modified_model = modified_model,
       modified_compiled = modified_compiled,
       modified_ss = ss, bk_ok = bk_ok, method = "in_memory",
       instruments_added = instruments_added)
}

.nn1_build_objective_expr <- function(coefficients, endo_names, ss, n) {
  terms <- character(0)
  n_endo <- length(endo_names)
  dev_vars <- character(n_endo)
  for (i in seq_len(n_endo)) {
    nm <- endo_names[i]
    ssv <- ss[[nm]] %||% 0
    dev_vars[i] <- if (abs(ssv) < 1e-15) nm else sprintf("(%s-%.10g)", nm, ssv)
  }
  quad <- coefficients[["quad"]]
  if (!is.null(quad)) {
    n_q <- min(n_endo, ncol(quad))
    for (i in seq_len(n_q)) {
      for (j in i:n_q) {
        coeff <- quad[i, j]
        if (abs(coeff) < 1e-14) next
        cv <- coeff * (if (i == j) 0.5 else 1.0)
        if (abs(cv) < 1e-14) next
        if (i == j) {
          terms <- c(terms, sprintf("%.10g*%s^2", cv, dev_vars[i]))
        } else {
          terms <- c(terms, sprintf("%.10g*%s*%s", cv, dev_vars[i], dev_vars[j]))
        }
      }
    }
  }
  cubic <- coefficients[["cubic"]]
  if (!is.null(cubic) && n >= 2) {
    n_c <- min(n_endo, dim(cubic)[1])
    for (i in seq_len(n_c)) {
      for (j in i:n_c) {
        for (k in j:n_c) {
          coeff <- cubic[i, j, k]
          if (abs(coeff) < 1e-14) next
          cv <- coeff / 6
          if (abs(cv) < 1e-14) next
          vars <- c(dev_vars[i], dev_vars[j], dev_vars[k])
          vc <- table(vars)
          factors <- c()
          for (nm in names(vc)) {
            factors <- c(factors, if (vc[[nm]] == 1) nm else sprintf("%s^%d", nm, vc[[nm]]))
          }
          terms <- c(terms, sprintf("%.10g*%s", cv, paste(factors, collapse = "*")))
        }
      }
    }
  }
  quartic <- coefficients[["quartic"]]
  if (!is.null(quartic) && n >= 3) {
    n_q4 <- min(n_endo, dim(quartic)[1])
    for (i in seq_len(n_q4)) {
      for (j in i:n_q4) {
        for (k in j:n_q4) {
          for (l in k:n_q4) {
            coeff <- quartic[i, j, k, l]
            if (abs(coeff) < 1e-14) next
            cv <- coeff / 24
            if (abs(cv) < 1e-14) next
            vars <- c(dev_vars[i], dev_vars[j], dev_vars[k], dev_vars[l])
            vc <- table(vars)
            factors <- c()
            for (nm in names(vc)) {
              factors <- c(factors, if (vc[[nm]] == 1) nm else sprintf("%s^%d", nm, vc[[nm]]))
            }
            terms <- c(terms, sprintf("%.10g*%s", cv, paste(factors, collapse = "*")))
          }
        }
      }
    }
  }
  if (length(terms) == 0) return("0")
  paste0("(", paste(terms, collapse = " + "), ")")
}

.parse_col_names <- function(col_names, endo_names, exo_names) {
  info <- data.frame(name = character(length(col_names)),
                     lead_lag = integer(length(col_names)),
                     is_exo = logical(length(col_names)),
                     stringsAsFactors = FALSE)
  for (i in seq_along(col_names)) {
    cn <- col_names[i]
    parts <- strsplit(cn, "__")[[1]]
    info$name[i] <- parts[1]
    info$is_exo[i] <- parts[1] %in% exo_names
    if (length(parts) >= 2) {
      ll_str <- parts[2]
      info$lead_lag[i] <- switch(substr(ll_str, 1, 1),
        "m" = -as.integer(substr(ll_str, 2, nchar(ll_str))),
        "p" = as.integer(substr(ll_str, 2, nchar(ll_str))), 0L)
    }
  }
  info
}

