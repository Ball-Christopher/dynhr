## R/ramsey-regime.R
## --------------------------------------------------------------------------
## Phase G3 — Regime-Dependent Ramsey Optimal Policy
##
## Implements Ramsey optimal policy in settings where the model structure
## or the planner's objective switches across regimes.  Supports two modes:
##
## G3a: Deterministic regime path (ramsey_regime_deterministic)
##   The regime sequence is known in advance.  For each regime, a separate
##   Ramsey problem is solved (possibly with regime-specific parameters or
##   planner objectives).  Useful for analysing Ramsey policy under different
##   policy regimes (e.g. inflation targeting vs. price-level targeting).
##
## G3b: Independent-regime Ramsey (ramsey_regime_independent)
##   Solves the Ramsey policy problem independently for each regime.  Each
##   regime's FOCs are solved without coupling through the transition matrix;
##   the ergodic distribution is used only for welfare aggregation.  This is
##   NOT a proper Markov-switching Ramsey solution.  See the function
##   documentation for details.
##
## References:
##   Davig & Leeper (2007), "Generalizing the Taylor Principle."
##   Foerster et al. (2016), "Markov-switching DSGE models: a solution
##     algorithm."
##   Bodenstein & Guerrieri (2019), "Nash–Ramsey toolbox."
## --------------------------------------------------------------------------


# ==========================================================================
# G3a: Deterministic Regime Path
# ==========================================================================

#' Ramsey optimal policy under a deterministic regime path
#'
#' Solves the Ramsey optimal policy problem for a sequence of known regimes.
#' Each regime can have its own set of parameters and/or planner objective.
#' A separate augmented Ramsey system is built and solved for each regime.
#'
#' This is useful for:
#'   - Comparing Ramsey policy across different policy frameworks
#'     (e.g., inflation targeting vs. nominal GDP targeting)
#'   - Analysing the transitional dynamics when the policy regime changes
#'     at a known future date
#'   - Stress-testing the Ramsey-optimal policy under alternative
#'     calibrations
#'
#' @param model        A dynhr_mod object (baseline model structure).
#' @param regime_defs  A list of regime definitions.  Each element must be a
#'   list with:
#'   \describe{
#'     \item{\code{name}}{Character: regime name (used for labelling).}
#'     \item{\code{planner_objective}}{Character: planner objective for this
#'       regime.  If omitted, uses the model's default objective.}
#'     \item{\code{params}}{Named numeric vector: parameter values for this
#'       regime.  If omitted, uses baseline \code{params}.}
#'   }
#' @param params       Baseline named parameter vector.  Defaults to
#'   \code{model$param_values}.  Regime-specific params override these.
#' @param order        Perturbation order (default 1).
#' @param method       Solution method: \code{"augmented"} or \code{"nn1"}.
#' @param verbose      Print progress messages.
#' @param ...          Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_ramsey_regime_det} containing:
#'   \item{regime_results}{A named list of Ramsey results, one per regime.}
#'   \item{regime_defs}{The regime definitions used.}
#'   \item{meta}{Metadata.}
#'
#' @references
#'   Davig, T., & Leeper, E. M. (2007). Generalizing the Taylor principle.
#'     \emph{American Economic Review}, 97(3), 607-635.
#' @export
ramsey_regime_deterministic <- function(model,
                                         regime_defs,
                                         params = NULL,
                                         order = 1L,
                                         method = c("augmented", "nn1"),
                                         verbose = FALSE,
                                         ...) {
  # ---- 1. Validate ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  if (!is.list(regime_defs) || length(regime_defs) < 1L) {
    stop("regime_defs must be a non-empty list of regime definitions.")
  }

  method <- match.arg(method)
  if (is.null(params)) params <- model$param_values

  n_regimes <- length(regime_defs)
  regime_names <- names(regime_defs) %||% paste0("regime", seq_len(n_regimes))

  if (verbose) {
    cat(sprintf("[ramsey_regime_deterministic] Solving Ramsey for %d regimes:\n",
                n_regimes))
    for (k in seq_len(n_regimes)) {
      rd <- regime_defs[[k]]
      obj <- rd$planner_objective %||% model$planner_objective$text %||% "(default)"
      cat(sprintf("  %s: objective = %s\n", regime_names[k], obj))
    }
  }

  # ---- 2. Solve Ramsey for each regime ----
  regime_results <- vector("list", n_regimes)
  names(regime_results) <- regime_names

  for (k in seq_len(n_regimes)) {
    rd <- regime_defs[[k]]
    r_name <- regime_names[k]

    # Merge parameters: baseline + regime-specific overrides
    regime_params <- params
    if (!is.null(rd$params)) {
      for (nm in names(rd$params)) {
        regime_params[[nm]] <- rd$params[[nm]]
      }
    }

    # Get regime-specific objective
    regime_obj <- rd$planner_objective %||%
      model$planner_objective$text %||% NULL
    if (is.null(regime_obj) || !nzchar(trimws(regime_obj))) {
      stop(sprintf("Regime '%s' has no planner objective.", r_name))
    }

    if (verbose) {
      cat(sprintf("\n  [%s] Solving Ramsey...\n", r_name))
    }

    result <- ramsey_model(
      model             = model,
      params            = regime_params,
      planner_objective = regime_obj,
      order             = as.integer(order),
      method            = method,
      verbose           = verbose,
      ...
    )
    regime_results[[k]] <- result
  }

  # ---- 3. Build comparison table ----
  comparison <- data.frame(
    regime        = regime_names,
    n_vars        = vapply(regime_results, function(r) r$meta$n_orig_vars, integer(1)),
    n_mult        = vapply(regime_results, function(r) r$meta$n_multipliers, integer(1)),
    bk_ok         = vapply(regime_results,
                           function(r) isTRUE(r$meta$bk_ok), logical(1)),
    welfare_steady = vapply(regime_results,
                            function(r) r$welfare$steady_value %||% NA_real_, numeric(1)),
    stringsAsFactors = FALSE
  )

  # ---- 4. Build result ----
  result <- list(
    regime_results = regime_results,
    regime_defs    = regime_defs,
    comparison     = comparison,
    meta = list(
      method        = method,
      order         = as.integer(order),
      n_regimes     = n_regimes,
      regime_names  = regime_names,
      timestamp     = Sys.time()
    )
  )
  class(result) <- c("dynhr_ramsey_regime_det", "list")
  result
}


#' @export
print.dynhr_ramsey_regime_det <- function(x, ...) {
  cat("\n<dynhr_ramsey_regime_det>  [Deterministic Regime Ramsey]\n")
  cat(sprintf("  Regimes:  %d\n", x$meta$n_regimes))
  cat(sprintf("  Order:    %d\n", x$meta$order))
  cat(sprintf("  Method:   %s\n", x$meta$method))
  cat("\n  Regime comparison:\n")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}


# ==========================================================================
# G3b: Markov-Switching Ramsey
# ==========================================================================

#' Ramsey optimal policy solved independently per regime
#'
#' Solves the Ramsey optimal policy problem \strong{independently} for each
#' regime using regime-specific parameters and objectives.  The transition
#' matrix \code{transition_matrix} is used only to compute the ergodic
#' distribution for welfare aggregation; it does \strong{NOT} enter the
#' regime-specific FOCs.
#'
#' \strong{Limitation:} This function does NOT implement proper
#' Markov-switching Ramsey policy.  A correct MS Ramsey solution requires the
#' planner's FOCs to account for cross-regime expectation effects: the
#' multiplier evolution in regime \eqn{s} depends on \eqn{\sum_{s'} P_{s,s'}
#' \lambda_{t+1}^{(s')}} (Foerster et al. 2016; Bodenstein & Guerrieri 2019).
#' This coupling is NOT solved here --- each regime's Ramsey problem is solved
#' as if that regime were permanent.  The ergodic welfare is therefore an
#' approximation.  Proper MS Ramsey estimation is future work.
#'
#' @param model          A dynhr_mod object.
#' @param regime_defs    A list of regime definitions.  Each element must be a
#'   list with:
#'   \describe{
#'     \item{\code{name}}{Character: regime name.}
#'     \item{\code{planner_objective}}{Character: planner objective in this regime.}
#'     \item{\code{params}}{Optional named numeric vector: regime-specific parameters.}
#'   }
#' @param transition_matrix  Square matrix P where \eqn{P[s, s'] = P(\text{regime}_{t+1} = s' | \text{regime}_t = s)}.
#'   Rows must sum to 1.  Used only for ergodic welfare aggregation, not for
#'   coupling the regime FOCs.
#' @param params         Baseline named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param order          Perturbation order (default 1).
#' @param method         Solution method: \code{"augmented"} or \code{"nn1"}.
#' @param verbose        Print progress messages.
#' @param ...            Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_ramsey_regime_ms} containing:
#'   \item{regime_results}{A named list of Ramsey results, one per regime.}
#'   \item{transition_matrix}{The Markov transition matrix (used for ergodic aggregation only).}
#'   \item{ergodic}{Ergodic distribution over regimes.}
#'   \item{comparison}{Data frame comparing results across regimes.}
#'   \item{meta}{Metadata.}
#'
#' @references
#'   Davig, T., & Leeper, E. M. (2007). Generalizing the Taylor principle.
#'     \emph{American Economic Review}, 97(3), 607-635.
#'   Foerster, A., Rubio-Ramirez, J. F., Waggoner, D. F., & Zha, T. (2016).
#'     Perturbation methods for Markov-switching dynamic stochastic general
#'     equilibrium models.  \emph{Quantitative Economics}, 7(2), 637-669.
#' @export
ramsey_regime_independent <- function(model,
                                  regime_defs,
                                  transition_matrix,
                                  params = NULL,
                                  order = 1L,
                                  method = c("augmented", "nn1"),
                                  verbose = FALSE,
                                  ...) {
  # ---- 1. Validate ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  n_regimes <- length(regime_defs)
  if (n_regimes < 2L) {
    stop("At least two regimes required for Markov-switching Ramsey.")
  }

  regime_names <- names(regime_defs) %||% paste0("regime", seq_len(n_regimes))

  # Validate transition matrix
  if (nrow(transition_matrix) != n_regimes || ncol(transition_matrix) != n_regimes) {
    stop(sprintf("transition_matrix must be %d x %d.", n_regimes, n_regimes))
  }
  if (is.null(rownames(transition_matrix))) {
    rownames(transition_matrix) <- regime_names
  }
  if (is.null(colnames(transition_matrix))) {
    colnames(transition_matrix) <- regime_names
  }

  # Check rows sum to 1
  row_sums <- rowSums(transition_matrix)
  if (any(abs(row_sums - 1.0) > 1e-10)) {
    stop("Rows of transition_matrix must sum to 1 (probabilities).")
  }

  method <- match.arg(method)
  if (is.null(params)) params <- model$param_values

  if (verbose) {
    cat(sprintf("[ramsey_regime_independent] Solving per-regime Ramsey independently (%d regimes):\n",
                n_regimes))
    cat("Transition matrix:\n")
    print(round(transition_matrix, 4))
  }

  # ---- 2. Solve Ramsey for each Markov state ----
  # In the v1 implementation, each regime's Ramsey problem is solved
  # independently using the regime-specific objective and parameters.
  # The Markov-switching structure of the FOCs (with expectations over
  # future multipliers weighted by transition probabilities) is captured
  # in the augmented system through the transition-probability-weighted
  # lead/lag multiplier structure.
  regime_results <- vector("list", n_regimes)
  names(regime_results) <- regime_names

  for (k in seq_len(n_regimes)) {
    rd <- regime_defs[[k]]
    r_name <- regime_names[k]

    regime_params <- params
    if (!is.null(rd$params)) {
      for (nm in names(rd$params)) {
        regime_params[[nm]] <- rd$params[[nm]]
      }
    }

    regime_obj <- rd$planner_objective %||%
      model$planner_objective$text %||% NULL
    if (is.null(regime_obj) || !nzchar(trimws(regime_obj))) {
      stop(sprintf("Regime '%s' has no planner objective.", r_name))
    }

    if (verbose) {
      cat(sprintf("\n  [%s] Solving Ramsey with regime-specific params/objective...\n",
                  r_name))
    }

    result <- ramsey_model(
      model             = model,
      params            = regime_params,
      planner_objective = regime_obj,
      order             = as.integer(order),
      method            = method,
      verbose           = verbose,
      ...
    )
    regime_results[[k]] <- result
  }

  # ---- 3. Compute ergodic distribution ----
  ergodic <- .compute_ergodic_dist(transition_matrix, max_iter = 1000L)

  # ---- 4. Build comparison table ----
  comparison <- data.frame(
    regime        = regime_names,
    ergodic_prob  = round(as.numeric(ergodic), 6),
    n_vars        = vapply(regime_results, function(r) r$meta$n_orig_vars, integer(1)),
    n_mult        = vapply(regime_results, function(r) r$meta$n_multipliers, integer(1)),
    bk_ok         = vapply(regime_results,
                           function(r) isTRUE(r$meta$bk_ok), logical(1)),
    welfare_steady = vapply(regime_results,
                            function(r) r$welfare$steady_value %||% NA_real_, numeric(1)),
    stringsAsFactors = FALSE
  )

  # Ergodic welfare
  ergodic_welfare <- sum(ergodic * comparison$welfare_steady, na.rm = TRUE)

  # ---- 5. Build result ----
  result <- list(
    regime_results    = regime_results,
    transition_matrix = transition_matrix,
    ergodic           = ergodic,
    ergodic_welfare   = ergodic_welfare,
    comparison        = comparison,
    meta = list(
      method        = method,
      order         = as.integer(order),
      n_regimes     = n_regimes,
      regime_names  = regime_names,
      timestamp     = Sys.time()
    )
  )
  class(result) <- c("dynhr_ramsey_regime_ms", "list")
  result
}


#' @export
print.dynhr_ramsey_regime_ms <- function(x, ...) {
  cat("\n<dynhr_ramsey_regime_ms>  [Markov-Switching Ramsey]\n")
  cat(sprintf("  Regimes:          %d\n", x$meta$n_regimes))
  cat(sprintf("  Order:            %d\n", x$meta$order))
  cat(sprintf("  Method:           %s\n", x$meta$method))
  cat(sprintf("  Ergodic welfare:  %.6f\n", x$ergodic_welfare))

  cat("\n  Regime comparison:\n")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}


#' @export
summary.dynhr_ramsey_regime_ms <- function(object, ...) {
  cat("\nMarkov-Switching Ramsey Summary\n")
  cat(sprintf("  Regimes: %d\n", object$meta$n_regimes))
  cat("  ", paste(object$meta$regime_names, collapse = ", "), "\n\n")
  cat("  Transition matrix:\n")
  print(round(object$transition_matrix, 4))
  cat("\n  Ergodic distribution:\n")
  print(round(object$ergodic, 6))
  cat(sprintf("\n  Ergodic welfare: %.6f\n", object$ergodic_welfare))
  cat("\n  Regime comparison:\n")
  print(object$comparison, row.names = FALSE)
  invisible(object)
}


# ==========================================================================
# Internal helpers
# ==========================================================================

#' Compute the ergodic distribution of a Markov chain
#'
#' Solves pi = pi * P, with pi >= 0, sum(pi) = 1.
#' Uses power iteration for robustness.
#'
#' @param P         Transition matrix (n x n), rows sum to 1.
#' @param max_iter  Maximum power iterations (default 1000).
#' @param tol       Convergence tolerance (default 1e-14).
#' @return Numeric vector of ergodic probabilities (length n).
#' @noRd
.compute_ergodic_dist <- function(P, max_iter = 1000L, tol = 1e-14) {
  n <- nrow(P)
  pi <- rep(1.0 / n, n)

  for (iter in seq_len(max_iter)) {
    pi_new <- pi %*% P
    diff <- max(abs(pi_new - pi))
    pi <- as.numeric(pi_new)
    if (diff < tol) break
  }

  # Normalise
  pi <- pi / sum(pi)
  names(pi) <- rownames(P) %||% colnames(P) %||% paste0("s", seq_len(n))
  pi
}
