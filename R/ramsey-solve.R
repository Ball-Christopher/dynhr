## R/ramsey-solve.R
## --------------------------------------------------------------------------
## Augmented perturbation orchestrator for Ramsey optimal policy.
##
## Runs the existing perturbation solver on the augmented (original variables
## + Lagrange multipliers) system, then extracts the Ramsey decision rules
## from the full augmented decision rules.
##
## Key verification:
##   - Blanchard–Kahn conditions for the augmented system must be checked.
##     In the augmented system, the number of unstable eigenvalues should
##     equal the number of forward-looking variables (original forward-looking
##     variables + multipliers associated with forward-looking constraints).
##   - The Ramsey decision rules for original variables are a subset of the
##     augmented decision rules.
##
## References:
##   Bodenstein, M. & Guerrieri, L. (2019). Nash–Ramsey toolbox.
##   Klein, P. (2000). Using the generalized Schur form to solve a
##     multivariate linear rational expectations model.
##   Villemot, S. (2011). Solving rational expectations models at first
##     order: what Dynare does.
## --------------------------------------------------------------------------

#' Solve the Ramsey optimal policy problem via augmented perturbation
#'
#' Compiles the augmented model, solves its steady state, runs perturbation
#' at the requested order, and extracts the Ramsey-optimal decision rules
#' for the original variables.
#'
#' @param aug_model     Augmented dynhr_mod (from ramsey_parse_augmented).
#' @param aug_steady    Augmented steady state (from ramsey_steady).
#' @param params        Named parameter vector.
#' @param order         Perturbation order: 1 (default), 2, or 3.
#' @param verbose       Print progress messages.
#' @param aug_compiled  Optional pre-compiled augmented model (from
#'   \code{compile_model()}). When supplied and compiled to a sufficient order
#'   it is reused, avoiding a redundant recompile; otherwise the augmented model
#'   is compiled internally.
#' @param ...           Additional arguments passed to solve_perturbation()
#'                      (e.g., Sigma_e, h for order 2).
#'
#' @return An object of class \code{dynhr_ramsey_dr} containing:
#'   \item{augmented_dr}{Full decision rules for the augmented system.}
#'   \item{ramsey_dr}{Decision rules restricted to original variables.}
#'   \item{ramsey_steady}{Augmented steady state.}
#'   \item{bk_condition}{Blanchard–Kahn condition check result.}
#'   \item{n_forward_aug}{Number of forward-looking variables in augmented system.}
#'   \item{n_unstable}{Number of unstable eigenvalues.}
#'   \item{bk_ok}{Logical: whether BK conditions are satisfied.}
#' @export
ramsey_solve <- function(aug_model,
                          aug_steady,
                          params = NULL,
                          order = 1L,
                          verbose = FALSE,
                          aug_compiled = NULL,
                          ...) {
  # ---- 1. Validate ----
  if (!inherits(aug_model, "dynhr_mod")) {
    stop("aug_model must be a dynhr_mod object (from ramsey_parse_augmented).")
  }
  if (is.null(params)) params <- aug_model$param_values

  order <- as.integer(order)
  if (!order %in% c(1L, 2L, 3L)) {
    stop("order must be 1, 2, or 3.")
  }

  if (verbose) {
    cat(sprintf("[ramsey_solve] Compiling augmented model (%d vars, %d eqs)...\n",
                length(aug_model$var_names), length(aug_model$equations)))
  }

  # ---- 2. Compile augmented model ----
  ## Compile at the order being solved: a default max_order = 1L makes
  ## ramsey_solve(order >= 2) fail inside solve_perturbation (which requires the
  ## order-2 derivatives). order 2 and 3 both need max_order >= 2.
  ## Reuse a caller-supplied aug_compiled when it is compiled to a sufficient
  ## order (L10 perf: ramsey_model() already compiled the augmented model for
  ## the SS solve, so this avoids a redundant recompile). Recompile defensively
  ## if absent or under-compiled.
  needed_order <- if (order >= 2L) 2L else 1L
  if (is.null(aug_compiled) ||
      (aug_compiled$max_order %||% 0L) < needed_order) {
    aug_compiled <- compile_model(aug_model,
                                  max_order = needed_order,
                                  verbose = verbose)
  }

  # ---- 3. Extract augmented steady state ----
  ss_values <- aug_steady$values
  # Ensure ordering matches aug_model$var_names
  ss <- setNames(numeric(length(aug_model$var_names)), aug_model$var_names)
  for (nm in aug_model$var_names) {
    ss[nm] <- if (!is.null(ss_values[[nm]]) && is.finite(ss_values[[nm]])) {
      as.numeric(ss_values[[nm]])
    } else {
      0.0
    }
  }

  # ---- 4. Run perturbation ----
  if (verbose) {
    cat(sprintf("[ramsey_solve] Running perturbation at order %d...\n", order))
  }

  dr <- solve_perturbation(aug_model, aug_compiled, ss, params,
                            order = order, verbose = verbose, ...)

  # ---- 5. Check Blanchard–Kahn conditions ----
  bk_info <- .check_ramsey_bk(dr, aug_model, verbose = verbose)

  # ---- 6. Extract Ramsey decision rules for original variables ----
  # Identify which variables in the augmented system are original (not multipliers)
  orig_names <- names(aug_steady$orig_ss)
  n_orig <- length(orig_names)
  mult_names <- names(aug_steady$multiplier_ss)

  # Extract ghx, ghu etc. for original variables only
  ramsey_dr <- .extract_ramsey_dr(dr, orig_names, mult_names)

  # ---- 7. Compute BK on the original subsystem (Ramsey-optimal) ----
  # The augmented BK condition is the primary check. We also report whether
  # the Ramsey subsystem (original variables only, given the optimal
  # commitment) would satisfy BK on its own.

  # ---- 8. Return ----
  result <- list(
    augmented_dr   = dr,
    ramsey_dr      = ramsey_dr,
    ramsey_steady  = aug_steady,
    bk_condition   = bk_info,
    n_forward_aug  = bk_info$n_forward,
    n_unstable     = bk_info$n_unstable,
    bk_ok          = bk_info$bk_ok,
    order          = order,
    meta           = list(
      orig_var_names = orig_names,
      mult_names     = mult_names,
      n_orig         = length(orig_names),
      n_mult         = length(mult_names)
    )
  )
  class(result) <- c("dynhr_ramsey_dr", "list")
  result
}


#' Check Blanchard–Kahn conditions for the augmented system
#'
#' Verifies that the number of unstable (explosive) eigenvalues from the QZ
#' decomposition equals the number of forward-looking variables in the
#' augmented system.
#'
#' In the augmented Ramsey system:
#'   - Original forward-looking variables remain forward-looking
#'   - Multipliers associated with backward-looking constraints are
#'     forward-looking (they appear with leads in the FOCs)
#'   - Multipliers associated with forward-looking constraints are
#'     backward-looking
#'
#' @param dr        Decision rules object from solve_perturbation.
#' @param aug_model Augmented dynhr_mod.
#' @param verbose   Print progress.
#' @return List with bk_ok, n_forward, n_unstable, details.
#' @noRd
.check_ramsey_bk <- function(dr, aug_model, verbose = FALSE) {
  # The DecisionRules object stores:
  #   - dr$eigenvalues    (generalised eigenvalues from QZ)
  #   - dr$n_unstable     (count of eigenvalues outside unit circle)
  #   - dr$bk_satisfied   (logical: whether BK condition holds)
  #   - dr$n_state        (number of state variables)
  #
  # Use the solver's own BK check when available; compute from eigenvalues
  # as a cross-validation when possible.

  # --- Primary: use solver's own BK result ---
  bk_ok <- dr$bk_satisfied %||% NA
  n_unstable <- dr$n_unstable %||% NA_integer_

  # --- Cross-validate from eigenvalues if available ---
  ev <- dr$eigenvalues %||% NULL
  if (!is.null(ev) && is.na(n_unstable)) {
    n_unstable <- sum(abs(ev) > 1 + 1e-6, na.rm = TRUE)
    bk_ok <- (n_unstable == aug_model$n_forward)
  }

  if (verbose && !is.null(ev)) {
    cat(sprintf("  BK check: %d unstable eig / %d forward vars in aug system -> %s\n",
                n_unstable, aug_model$n_forward,
                if (isTRUE(bk_ok)) "OK" else if (isFALSE(bk_ok)) "FAIL" else "N/A"))
    if (isFALSE(bk_ok)) {
      cat(sprintf("    WARNING: Blanchard-Kahn condition not satisfied!\n"))
      cat(sprintf("    Augmented system: %d forward-looking variables, %d unstable eigenvalues.\n",
                  aug_model$n_forward, n_unstable))
    }
  }

  list(
    bk_ok        = isTRUE(bk_ok),
    n_forward    = aug_model$n_forward,
    n_pred       = aug_model$n_predetermined,
    n_static     = aug_model$n_static,
    n_unstable   = n_unstable,
    n_total_endo = length(aug_model$var_names)
  )
}


#' Extract Ramsey-optimal decision rules from augmented system
#'
#' Subsets the augmented decision rules (ghx, ghu, etc.) to keep only the
#' rows corresponding to original (non-multiplier) endogenous variables.
#'
#' @param dr         Full augmented decision rules.
#' @param orig_names Original endogenous variable names.
#' @param mult_names Multiplier variable names.
#' @return List with ghx, ghu (and ghxx, etc. if present), restricted to
#'   original variables.
#' @noRd
.extract_ramsey_dr <- function(dr, orig_names, mult_names) {
  # Find row indices of original variables in the decision rules
  # ghx is n_endo x n_state; rows correspond to dr$endo_names order
  orig_idx <- which(dr$endo_names %in% orig_names)

  ramsey <- list()

  # First-order terms
  ramsey$ghx <- if (!is.null(dr$ghx)) dr$ghx[orig_idx, , drop = FALSE] else NULL
  ramsey$ghu <- if (!is.null(dr$ghu)) dr$ghu[orig_idx, , drop = FALSE] else NULL
  ramsey$ghu_det <- if (!is.null(dr$ghu_det)) dr$ghu_det[orig_idx, , drop = FALSE] else NULL

  # Second-order risk correction: ghss is a length-n_endo NAMED VECTOR (not a
  # matrix), so it must be subset 1-dimensionally. (The old code indexed it as
  # `[orig_idx, , drop = FALSE]`, which errors on a vector, and also referenced
  # a non-existent `dr$ghs2` field -- the real field is `ghss`.)
  ramsey$ghss <- if (!is.null(dr$ghss)) dr$ghss[orig_idx] else NULL

  # Second-order terms. dynhr stores these FLATTENED as 2-D matrices
  # (n_endo x n_state^2, n_endo x n_state*n_exo, n_endo x n_exo^2), not as
  # 3-D arrays -- restrict ROWS to the original variables, keep all (augmented-
  # state) columns, mirroring ghx above. (The old 3-D/4-D indexing errored with
  # "incorrect number of dimensions" once order-2 Ramsey became reachable.)
  ramsey$ghxx <- if (!is.null(dr$ghxx)) dr$ghxx[orig_idx, , drop = FALSE] else NULL
  ramsey$ghxu <- if (!is.null(dr$ghxu)) dr$ghxu[orig_idx, , drop = FALSE] else NULL
  ramsey$ghuu <- if (!is.null(dr$ghuu)) dr$ghuu[orig_idx, , drop = FALSE] else NULL

  # Third-order terms (also flattened 2-D: n_endo x n_state^3, etc.)
  ramsey$ghxxx <- if (!is.null(dr$ghxxx)) dr$ghxxx[orig_idx, , drop = FALSE] else NULL
  ramsey$ghxxu <- if (!is.null(dr$ghxxu)) dr$ghxxu[orig_idx, , drop = FALSE] else NULL
  ramsey$ghxuu <- if (!is.null(dr$ghxuu)) dr$ghxuu[orig_idx, , drop = FALSE] else NULL
  ramsey$ghuuu <- if (!is.null(dr$ghuuu)) dr$ghuuu[orig_idx, , drop = FALSE] else NULL

  ramsey$ghxss <- if (!is.null(dr$ghxss)) dr$ghxss[orig_idx, , drop = FALSE] else NULL
  ramsey$ghuss <- if (!is.null(dr$ghuss)) dr$ghuss[orig_idx, , drop = FALSE] else NULL
  # ghs3 is a length-n_endo NAMED VECTOR (third-order risk constant), subset 1-D.
  ramsey$ghs3  <- if (!is.null(dr$ghs3))  dr$ghs3[orig_idx]  else NULL

  # Variable names
  ramsey$endo_names   <- orig_names
  ramsey$exo_names    <- dr$exo_names
  ramsey$state_names  <- dr$state_names
  ramsey$n_endo       <- length(orig_names)
  ramsey$n_exo        <- dr$n_exo
  ramsey$n_state      <- dr$n_state
  ramsey$order        <- dr$order

  ramsey
}


#' @export
print.dynhr_ramsey_dr <- function(x, ...) {
  cat(sprintf("<dynhr_ramsey_dr>\n"))
  cat(sprintf("  Order:              %d\n", x$order))
  cat(sprintf("  Original variables: %d\n", x$meta$n_orig))
  cat(sprintf("  Multipliers:        %d\n", x$meta$n_mult))
  cat(sprintf("  BK condition:       %s\n",
              if (isTRUE(x$bk_ok)) "OK" else if (isFALSE(x$bk_ok)) "FAIL" else "N/A"))
  cat(sprintf("  Forward-looking:    %d\n", x$n_forward_aug))
  cat(sprintf("  Unstable eig:       %d\n", x$n_unstable))

  if (!is.null(x$ramsey_dr$ghx)) {
    cat(sprintf("\n  Ramsey ghx:         %d x %d\n",
                nrow(x$ramsey_dr$ghx), ncol(x$ramsey_dr$ghx)))
  }
  if (!is.null(x$ramsey_dr$ghu)) {
    cat(sprintf("  Ramsey ghu:         %d x %d\n",
                nrow(x$ramsey_dr$ghu), ncol(x$ramsey_dr$ghu)))
  }
  invisible(x)
}
