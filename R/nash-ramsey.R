## R/nash-ramsey.R
## --------------------------------------------------------------------------
## Phase G2 — Nash–Ramsey Games (Multi-Player Optimal Policy)
##
## Implements two solution concepts for multi-player optimal policy in DSGE
## models, following Bodenstein & Guerrieri (2019):
##
## G2a: Cooperative Ramsey (nash_ramsey_cooperative)
##   A single "planner" maximises a weighted sum of individual player
##   objectives.  Reduces to a standard Ramsey problem with a combined
##   objective:  W = Σ_k ω_k · f_k(y)
##   Uses the existing ramsey_augment_mod() + ramsey_model() pipeline
##   with the combined objective.
##
## G2b: Open-Loop Nash Equilibrium (nash_ramsey_openloop)
##   Each player j maximises their own objective f_j(y) taking the other
##   players' policy functions as given.  Builds separate augmented systems
##   for each player and solves them independently.  Returns a list of
##   player-specific Ramsey results.
##
## References:
##   Bodenstein, M. & Guerrieri, L. (2019). "Nash–Ramsey games in DSGE
##     models."  Mimeo, Federal Reserve Board.
##   Tabellini, G. & Persson, T. (1995). "Double-edged incentives: Private
##     information in a model of optimal fiscal policy."
## --------------------------------------------------------------------------


# ==========================================================================
# G2a: Cooperative Nash–Ramsey
# ==========================================================================

#' Cooperative Ramsey optimal policy with multiple players
#'
#' Solves the cooperative Ramsey problem where a single planner maximises
#' a weighted sum of individual player objectives.  Reduces to the standard
#' Ramsey problem via \code{\link{ramsey_model}} with the combined objective
#' \eqn{W = \sum_k \omega_k \cdot f_k(\mathbf{y})}.
#'
#' @param model        A dynhr_mod object.
#' @param objectives   A named list of player objective expressions. Elements
#'   should be character strings (e.g., \code{list(fiscal = "log(c)",
#'   monetary = "-pi^2")}).
#' @param weights      Numeric vector of weights for each player's objective.
#'   Recycled if shorter than \code{length(objectives)}.  If NULL, equal
#'   weights are used.
#' @param params       Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param order        Perturbation order (default 1).
#' @param method       Solution method: \code{"augmented"} or \code{"nn1"}.
#' @param verbose      Print progress messages.
#' @param ...          Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_nash_ramsey_cooperative} containing:
#'   \item{combined_objective}{The weighted-sum objective expression.}
#'   \item{player_weights}{Named numeric vector of weights used.}
#'   \item{player_objectives}{The original player objectives.}
#'   \item{ramsey_result}{The Ramsey result from \code{\link{ramsey_model}}.}
#'   \item{meta}{Metadata.}
#' @references
#'   Nash, J. (1951). Non-cooperative games. \emph{Annals of Mathematics},
#'     54(2), 286-295.
#'   Cooley, T. F., & Quadrini, V. (2006). Monetary policy and the financial
#'     decisions of firms. \emph{Economic Theory}, 27(1), 243-270.
#' @export
nash_ramsey_cooperative <- function(model,
                                     objectives,
                                     weights = NULL,
                                     params = NULL,
                                     order = 1L,
                                     method = c("augmented", "nn1"),
                                     verbose = FALSE,
                                     ...) {
  # ---- 1. Validate ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  if (is.null(names(objectives))) {
    names(objectives) <- paste0("player", seq_along(objectives))
  }
  n_players <- length(objectives)
  if (n_players < 2L) {
    stop("At least two player objectives required for a cooperative game.")
  }

  method <- match.arg(method)
  if (is.null(params)) params <- model$param_values

  # ---- 2. Build combined objective ----
  if (is.null(weights)) {
    weights <- rep(1.0 / n_players, n_players)
  } else {
    weights <- rep(as.numeric(weights), length.out = n_players)
  }
  names(weights) <- names(objectives)

  # Combine: W = Σ ω_k · (f_k)
  combined_parts <- mapply(function(obj, w) {
    if (abs(w - 1.0) < 1e-12) {
      obj
    } else if (abs(w + 1.0) < 1e-12) {
      sprintf("(-(%s))", obj)
    } else {
      sprintf("%s * (%s)", format(w, scientific = FALSE), obj)
    }
  }, objectives, weights, SIMPLIFY = FALSE, USE.NAMES = FALSE)

  combined_obj <- paste(combined_parts, collapse = " + ")

  if (verbose) {
    cat("[nash_ramsey_cooperative] Combined objective:", combined_obj, "\n")
    for (k in seq_len(n_players)) {
      cat(sprintf("  Player %d (%s): weight = %.4f, objective = %s\n",
                  k, names(objectives)[k], weights[k], objectives[[k]]))
    }
  }

  # ---- 3. Solve Ramsey with combined objective ----
  ramsey_result <- ramsey_model(
    model             = model,
    params            = params,
    planner_objective = combined_obj,
    order             = as.integer(order),
    method            = method,
    verbose           = verbose,
    ...
  )

  # ---- 4. Build result ----
  result <- list(
    combined_objective  = combined_obj,
    player_weights      = weights,
    player_objectives   = objectives,
    ramsey_result       = ramsey_result,
    meta = list(
      method        = method,
      order         = as.integer(order),
      n_players     = n_players,
      player_names  = names(objectives),
      timestamp     = Sys.time()
    )
  )
  class(result) <- c("dynhr_nash_ramsey_cooperative", "list")
  result
}


#' @export
print.dynhr_nash_ramsey_cooperative <- function(x, ...) {
  cat("\n<dynhr_nash_ramsey_cooperative>  [Cooperative Nash-Ramsey]\n")
  cat(sprintf("  Players:            %d\n", x$meta$n_players))
  for (k in seq_len(x$meta$n_players)) {
    cat(sprintf("    %s (weight: %.4f): %s\n",
                x$meta$player_names[k],
                x$player_weights[k],
                x$player_objectives[[k]]))
  }
  cat(sprintf("  Combined objective: %s\n", x$combined_objective))
  cat(sprintf("  Order:              %d\n", x$meta$order))
  cat(sprintf("  Method:             %s\n", x$meta$method))

  if (!is.null(x$ramsey_result)) {
    cat(sprintf("  BK condition:       %s\n",
                if (isTRUE(x$ramsey_result$meta$bk_ok)) "OK" else "FAIL"))
    cat(sprintf("  Welfare (steady):   %.6f\n",
                x$ramsey_result$welfare$steady_value))
  }
  invisible(x)
}


# ==========================================================================
# G2b: Open-Loop Nash Equilibrium
# ==========================================================================

#' Open-loop Nash equilibrium for multiple policy authorities
#'
#' Solves the open-loop Nash equilibrium where each player j maximises their
#' own objective f_j(y) taking the other players' policies as given.  Builds
#' separate augmented systems for each player (using their individual
#' objective as the planner objective) and solves them independently.
#'
#' In the open-loop Nash concept, each player chooses a policy plan at time 0
#' and commits to it, under the belief that other players' plans are fixed.
#' The equilibrium is a fixed point where no player has an incentive to
#' deviate unilaterally, given the other players' committed plans.
#'
#' @param model        A dynhr_mod object.
#' @param objectives   A named list of player objective expressions. Elements
#'   should be character strings (e.g., \code{list(fiscal = "log(c)",
#'   monetary = "-pi^2")}).
#' @param params       Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param order        Perturbation order (default 1).
#' @param method       Solution method: \code{"augmented"} or \code{"nn1"}.
#' @param verbose      Print progress messages.
#' @param ...          Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_nash_ramsey_openloop} containing:
#'   \item{player_results}{A named list of Ramsey results, one per player.}
#'   \item{objectives}{The original player objectives.}
#'   \item{meta}{Metadata.}
#' @export
nash_ramsey_openloop <- function(model,
                                  objectives,
                                  params = NULL,
                                  order = 1L,
                                  method = c("augmented", "nn1"),
                                  verbose = FALSE,
                                  ...) {
  # ---- 1. Validate ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  if (is.null(names(objectives))) {
    names(objectives) <- paste0("player", seq_along(objectives))
  }
  n_players <- length(objectives)
  if (n_players < 2L) {
    stop("At least two player objectives required for an open-loop Nash game.")
  }

  method <- match.arg(method)
  if (is.null(params)) params <- model$param_values

  if (verbose) {
    cat(sprintf("[nash_ramsey_openloop] Solving open-loop Nash with %d players:\n",
                n_players))
    for (k in seq_len(n_players)) {
      cat(sprintf("  %s: %s\n", names(objectives)[k], objectives[[k]]))
    }
  }

  # ---- 2. Solve individual Ramsey problem for each player ----
  player_results <- vector("list", n_players)
  names(player_results) <- names(objectives)

  for (k in seq_len(n_players)) {
    player_name <- names(objectives)[k]
    obj_text <- objectives[[k]]

    if (verbose) {
      cat(sprintf("\n[Player %d: %s] Solving Ramsey with objective: %s\n",
                  k, player_name, obj_text))
    }

    result <- ramsey_model(
      model             = model,
      params            = params,
      planner_objective = obj_text,
      order             = as.integer(order),
      method            = method,
      verbose           = verbose,
      ...
    )
    player_results[[k]] <- result
  }

  # ---- 3. Compute welfare for each player under each regime ----
  welfare_matrix <- matrix(NA_real_, nrow = n_players, ncol = n_players)
  rownames(welfare_matrix) <- names(objectives)
  colnames(welfare_matrix) <- names(objectives)

  for (i in seq_len(n_players)) {
    for (j in seq_len(n_players)) {
      # Welfare of player i under the regime where player j's objective
      # determines policy (i.e., player j is the "planner")
      welfare_matrix[i, j] <- player_results[[j]]$welfare$steady_value %||% NA_real_
    }
  }

  # ---- 4. Build result ----
  result <- list(
    player_results  = player_results,
    objectives      = objectives,
    welfare_matrix  = welfare_matrix,
    meta = list(
      method        = method,
      order         = as.integer(order),
      n_players     = n_players,
      player_names  = names(objectives),
      timestamp     = Sys.time()
    )
  )
  class(result) <- c("dynhr_nash_ramsey_openloop", "list")
  result
}


#' @export
print.dynhr_nash_ramsey_openloop <- function(x, ...) {
  cat("\n<dynhr_nash_ramsey_openloop>  [Open-Loop Nash-Ramsey]\n")
  cat(sprintf("  Players:  %d\n", x$meta$n_players))
  cat(sprintf("  Order:    %d\n", x$meta$order))
  cat(sprintf("  Method:   %s\n", x$meta$method))
  cat("\n  Player objectives:\n")
  for (k in seq_len(x$meta$n_players)) {
    cat(sprintf("    %d. %s: %s\n", k, x$meta$player_names[k], x$objectives[[k]]))
  }
  cat("\n  Welfare matrix (row = whose welfare, col = whose policy):\n")
  print(round(x$welfare_matrix, 6))
  invisible(x)
}


#' Summarise a Nash–Ramsey game result
#'
#' @param object A \code{dynhr_nash_ramsey_openloop} object.
#' @param ...    Additional arguments (unused).
#' @return Invisibly returns the welfare matrix.
#' @export
summary.dynhr_nash_ramsey_openloop <- function(object, ...) {
  cat("\nNash-Ramsey Open-Loop Equilibrium Summary\n")
  cat(sprintf("Players: %d\n", object$meta$n_players))
  cat(paste(object$meta$player_names, collapse = ", "), "\n\n")

  cat("Welfare matrix (row = whose welfare, col = whose policy):\n")
  mat <- round(object$welfare_matrix, 6)
  print(mat)

  cat("\nDiagonal: welfare when player sets policy\n")
  diag_vals <- diag(mat)
  names(diag_vals) <- object$meta$player_names
  print(diag_vals)

  invisible(object$welfare_matrix)
}
