## R/discretionary-policy.R
## --------------------------------------------------------------------------
## Phase D — Discretionary Policy (Markov-perfect equilibrium)
##
## Solves the LQ discretionary policy problem via fixed-point iteration on
## the Riccati equation:
##
##   F = (R + β B' P B)^{-1} β B' P A
##   P = Q + F' R F + β (A - B F)' P (A - B F)
##
## where s' = A s + B u is the state-space representation of the model
## (without the policy rule), and the loss is:
##
##   L = E[Σ β^t (s' Q s + u' R u)]
##
## The algorithm:
##   1. Extract A, B from the model's structural Jacobian
##      (removing the policy-instrument equation)
##   2. Build Q (state penalty) and R (control penalty) from user loss spec
##   3. Iterate P ← Q + β A' P A - β² A' P B (R + β B' P B)^{-1} B' P A
##   4. Extract optimal policy rule u = -F s
##   5. Compute resulting dynamics, IRFs, and unconditional moments
##   6. Compare against Ramsey-optimal commitment (if provided)
## --------------------------------------------------------------------------

#' Discretionary (Markov-perfect) optimal policy
#'
#' Solves for the Markov-perfect equilibrium under discretion using
#' LQ fixed-point iteration. The model's structural Jacobian is used
#' to extract the state-space matrices A and B (without the policy
#' instrument equation). The optimal policy takes the form
#' \eqn{u_t = -F s_t} where \eqn{s_t} are the predetermined
#' (state) variables.
#'
#' @param model             A \code{dynhr_mod} object.
#' @param compiled          Optional pre-compiled model.
#' @param params            Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param policy_instrument Character: name of the policy instrument
#'   variable (e.g., \code{"i"}).
#' @param loss_vars         Character vector of variable names whose
#'   variances enter the loss function (typically a subset of endogenous
#'   variables, possibly including the instrument).
#' @param loss_weights      Named numeric vector of weights for each loss
#'   variable. Names must match \code{loss_vars}. If \code{NULL}, equal
#'   weights are used.
#' @param control_penalty   Optional extra penalty weight on the policy
#'   instrument, added to the instrument's loss weight. If \code{NULL}
#'   (default) no penalty is added: the Dennis (2007) recursion already
#'   regularises the instrument block through the value function, and adding
#'   a penalty biases the rule away from the true optimum. Set a small
#'   positive value only to regularise a genuinely singular instrument block.
#' @param discount          Discount factor. Defaults to \code{beta}
#'   parameter in the model, or 0.99.
#' @param max_iter          Maximum fixed-point iterations.
#' @param tol               Convergence tolerance on ||P_new - P||.
#' @param ramsey_result     Optional \code{dynhr_ramsey_result2} for
#'   welfare comparison against Ramsey-optimal commitment.
#' @param verbose           Print progress messages.
#'
#' @return An object of class \code{dynhr_discretionary_result} with:
#'   \describe{
#'     \item{F}{Policy rule matrix: u_t = -F s_t.}
#'     \item{P}{Value function matrix: V(s) = s' P s.}
#'     \item{A, B}{State-space matrices of the uncontrolled system.}
#'     \item{Q, R}{Loss matrices.}
#'     \item{policy_instrument}{Name of the policy instrument.}
#'     \item{state_vars}{Character vector of state variable names.}
#'     \item{discount}{Discount factor used.}
#'     \item{converged}{Logical: did the fixed-point iteration converge?}
#'     \item{iterations}{Number of iterations.}
#'     \item{dr}{DecisionRules under the optimal discretionary policy.}
#'     \item{moments}{Unconditional moments under discretion.}
#'     \item{irfs}{IRFs under discretion (20 periods).}
#'     \item{ramsey_comparison}{Comparison with Ramsey policy (if
#'       \code{ramsey_result} provided).}
#'   }
#'
#' @references
#'   Dennis, R. (2007). Optimal policy in rational expectations models: New
#'     solution algorithms. \emph{Macroeconomic Dynamics}, 11(1), 31-55.
#'   Söderlind, P. (1999). Solution and estimation of RE macromodels with
#'     optimal policy. \emph{European Economic Review}, 43(4-6), 813-823.
#' @export
discretionary_policy <- function(model,
                                  compiled = NULL,
                                  params = NULL,
                                  policy_instrument,
                                  loss_vars,
                                  loss_weights = NULL,
                                  control_penalty = NULL,
                                  discount = NULL,
                                  max_iter = 1000L,
                                  tol = 1e-12,
                                  ramsey_result = NULL,
                                  verbose = FALSE) {

  # ---- 1. Validate inputs ----
  if (!inherits(model, "dynhr_mod")) {
    stop("'model' must be a dynhr_mod object.")
  }
  if (!policy_instrument %in% model$var_names) {
    stop("'policy_instrument' must be an endogenous variable name.")
  }
  if (is.null(params)) params <- model$param_values

  missing_vars <- setdiff(loss_vars, model$var_names)
  if (length(missing_vars) > 0) {
    stop("Loss variable(s) not found in model: ",
         paste(missing_vars, collapse = ", "))
  }

  # Default discount factor
  if (is.null(discount)) {
    discount <- if ("beta" %in% names(params) && is.finite(params[["beta"]])) {
      as.numeric(params[["beta"]])
    } else {
      0.99
    }
  }

  # Default loss weights
  if (is.null(loss_weights)) {
    loss_weights <- setNames(rep(1, length(loss_vars)), loss_vars)
  } else if (is.null(names(loss_weights))) {
    if (length(loss_weights) != length(loss_vars)) {
      stop("Length of loss_weights must match length of loss_vars.")
    }
    loss_weights <- setNames(as.numeric(loss_weights), loss_vars)
  } else {
    loss_weights <- loss_weights[loss_vars]
    loss_weights[is.na(loss_weights)] <- 0
  }

  # Default control penalty.
  # Dennis (2007) inverts (Q + A3' D^-T P D^-1 A3), where the second term is
  # the value-function cost of the instrument; this is generically positive
  # definite even when the raw instrument weight Q is zero.  So NO control
  # penalty is needed for a determinate problem, and adding one biases the
  # rule away from the true optimum (and away from Dynare, whose objective
  # carries no instrument penalty unless written into planner_objective).
  # Default to 0; the user can still pass a positive value to regularise a
  # near-singular instrument block.
  if (is.null(control_penalty)) {
    control_penalty <- 0
  }

  # ---- 2. Compute steady state ----
  if (is.null(compiled)) {
    compiled <- compile_model(model, verbose = FALSE)
  }

  n_endo_check <- length(model$var_names)
  n_eq_check   <- length(model$equations)
  if (n_eq_check < n_endo_check) {
    # Non-square model (free instrument): solve_steady_state() expects a square
    # system and would error. For linear discretionary-policy models the SS is
    # zero everywhere; use initval or zeros.
    if (verbose) {
      cat(sprintf(
        "  Non-square model (%d eqs, %d vars): using initval as SS.\n",
        n_eq_check, n_endo_check))
    }
    ss <- setNames(rep(0, n_endo_check), model$var_names)
    if (length(model$initval) > 0) {
      for (nm in names(model$initval)) {
        if (nm %in% model$var_names) ss[nm] <- as.numeric(model$initval[[nm]])
      }
    }
  } else {
    ss_result <- solve_steady_state(model, compiled, params = params,
                                     verbose = FALSE)
    if (!isTRUE(ss_result$converged)) {
      stop("Steady state did not converge.")
    }
    ss <- ss_result$ss
  }

  # ---- 3. Build the LQ discretion inputs ----
  # Cast the model into the Dennis (2007) form
  #   AAlag y_{t-1} + AA0 y_t + AAlead y_{t+1} + BB e_t = 0
  # over the FULL endogenous vector y_t (predetermined states AND forward-
  # looking jumps), with the planner loss y_t' W y_t + x_t' Q x_t where x_t is
  # the instrument.  The weight matrix W is the Hessian of the planner
  # objective over ALL variables -- this is what correctly carries the weights
  # of forward-looking loss variables (pi, y_gap) that the previous state-only
  # Q-build dropped (issue H1).
  inp <- .build_discretion_inputs(
    model, compiled, ss, params,
    policy_instrument = policy_instrument,
    loss_vars = loss_vars, loss_weights = loss_weights,
    control_penalty = control_penalty, verbose = verbose
  )

  if (verbose) {
    cat(sprintf("[discretionary_policy] endo=%d eqs=%d instrument=%s\n",
                inp$n_endo, inp$n_eq, policy_instrument))
    wnz <- which(abs(diag(inp$bigw)) > 0)
    cat(sprintf("  Loss weights (W diag): %s\n",
                paste(sprintf("%s=%.4g", inp$endo[wnz], diag(inp$bigw)[wnz]),
                      collapse = ", ")))
  }

  # ---- 4. Solve the Dennis (2007) discretion recursion ----
  eng <- .discretion_dennis_engine(
    inp$AAlag, inp$AA0, inp$AAlead, inp$BB, inp$bigw, inp$instr_id,
    beta = discount, max_iter = max_iter, tol = tol, verbose = verbose
  )

  if (!isTRUE(eng$converged)) {
    stop(
      "Discretionary-policy fixed point did not converge (retcode=", eng$retcode,
      ", |Delta|=", signif(eng$diff, 3), " after ", eng$iterations,
      " iterations). The model may be indeterminate under discretion, or the ",
      "loss / discount specification may be ill-posed."
    )
  }

  H <- eng$H   # n_endo x n_endo : y_t = H y_{t-1}
  G <- eng$G   # n_endo x n_exo  : impact y_t = G e_t
  endo <- inp$endo
  exo  <- inp$exo

  # Feedback rule on the instrument: i_t = (H[i,] over lagged endo) y_{t-1}
  #                                       + G[i,] e_t
  inst_row <- which(endo == policy_instrument)
  # State variables = endogenous columns of H with any nonzero feedback.
  state_idx  <- which(apply(abs(H), 2, max) > 1e-12)
  if (length(state_idx) == 0L) state_idx <- integer(0)
  state_vars <- endo[state_idx]
  # F_policy follows the documented API convention u_t = -F s_{t-1}, so the
  # closed-loop instrument feedback i_t = H[i, state] s_{t-1} corresponds to
  # F = -H[i, state] (a positive F_j is a stabilising response to state j).
  F_policy <- if (length(state_idx) > 0L) -H[inst_row, state_idx] else numeric(0)
  names(F_policy) <- state_vars

  if (verbose) {
    cat(sprintf("  Converged in %d iterations (|Delta|=%.3e)\n",
                eng$iterations, eng$diff))
    cat(sprintf("  %d state variable(s) feed the instrument rule.\n",
                length(state_idx)))
  }

  # ---- 5. Build the closed-loop DecisionRules ----
  dr_disc <- .discretion_build_dr(model, ss, endo, exo, H, G, state_idx)

  # Compute unconditional moments
  moments_disc <- tryCatch(
    compute_moments(dr_disc, model, params = params),
    error = function(e) NULL
  )

  # Compute IRFs (20 periods, matching Dynare's default for these models)
  irfs_disc <- tryCatch(
    compute_irfs(dr_disc, model, n_periods = 20L, params = params),
    error = function(e) NULL
  )

  # ---- 6. Ramsey comparison (optional) ----
  ramsey_comp <- NULL
  if (!is.null(ramsey_result)) {
    ramsey_comp <- .disc_compare_ramsey(
      dr_disc, moments_disc, loss_vars, loss_weights,
      ramsey_result, discount
    )
  }

  # ---- 7. Assemble result ----
  result <- list(
    F                 = F_policy,
    P                 = eng$P,
    H                 = H,
    G                 = G,
    A                 = eng$A_state,
    B                 = eng$B_state,
    Q                 = inp$bigw,
    R                 = inp$bigw[inp$instr_id, inp$instr_id, drop = FALSE],
    W                 = inp$bigw,
    policy_instrument = policy_instrument,
    state_vars        = state_vars,
    state_idx         = state_idx,
    discount          = discount,
    converged         = eng$converged,
    iterations        = eng$iterations,
    diff              = eng$diff,
    dr                = dr_disc,
    moments           = moments_disc,
    irfs              = irfs_disc,
    ramsey_comparison = ramsey_comp,
    meta = list(
      params    = params,
      loss_vars = loss_vars,
      loss_weights = loss_weights,
      control_penalty = control_penalty,
      timestamp = Sys.time()
    )
  )
  class(result) <- c("dynhr_discretionary_result", "list")

  if (verbose) {
    cat(sprintf("\n[discretionary_policy] Done.\n"))
    cat(sprintf("  Converged: %s (%d iters, |Delta| = %.3e)\n",
                eng$converged, eng$iterations, eng$diff))
    if (length(state_idx) > 0L) {
      cat(sprintf("  Policy rule: %s_t = ", policy_instrument))
      # F follows u = -F s, so the actual loading is -F.
      cat(paste(sprintf("%.4f*%s(-1)", -F_policy, state_vars), collapse = " + "))
      cat("\n")
    }
  }

  result
}


# ==========================================================================
# Internal helpers
# ==========================================================================

# --------------------------------------------------------------------------
# Söderlind (1999) / Dennis (2007) discretion solver
# --------------------------------------------------------------------------
# These three helpers implement the proper LQ discretion solution that
# handles forward-looking / non-state loss variables (issue H1).  They
# replace the previous "reduce to A/B then DARE on the predetermined states"
# approach, which dropped the contribution of jump loss variables (pi,
# y_gap) and produced an all-zero feedback rule.
#
# The algorithm mirrors Dynare 7.1's discretionary_policy_engine.m
# (validated to machine precision against it on Gali (2008) Ch.5).

#' Build the Dennis (2007) discretion inputs from a model
#'
#' Casts the model into the form
#'   AAlag y_{t-1} + AA0 y_t + AAlead y_{t+1} + BB e_t = 0
#' over the FULL endogenous vector (states and jumps), and builds the loss
#' weight matrix bigw (n_endo x n_endo) as the Hessian of the planner
#' objective.  The instrument equation is absent (non-square model) -- the
#' planner is free to choose the instrument optimally.
#'
#' bigw is built from (in priority order):
#'   1. the model's planner_objective text (Hessian over all variables) when
#'      available -- the Dynare-faithful path; or
#'   2. the supplied loss_vars / loss_weights: a diagonal with 2*weight on
#'      each loss variable (matching the Hessian of sum-of-squares).
#' A control penalty (2*control_penalty) is added on the instrument diagonal
#' when the instrument carries no own loss weight, to keep Q (= bigw on the
#' instrument block) invertible.
#'
#' @noRd
.build_discretion_inputs <- function(model, compiled, ss, params,
                                     policy_instrument,
                                     loss_vars, loss_weights,
                                     control_penalty, verbose = FALSE) {
  dyn  <- compiled$dynamic
  endo <- model$var_names
  exo  <- model$varexo_names
  n_endo <- length(endo)
  n_exo  <- length(exo)
  n_eq   <- dyn$n_eq
  dcm    <- dyn$dyn_col_map

  # ---- Dynamic Jacobian at steady state ----
  dy   <- numeric(nrow(dcm))
  keys <- character(nrow(dcm))
  for (k in seq_len(nrow(dcm))) {
    nm <- dcm$name[k]; ll <- dcm$lead_lag[k]
    sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll)
           else paste0("__m", abs(ll))
    keys[k] <- paste0(nm, sfx)
    dy[k]   <- if (nm %in% names(ss)) ss[[nm]] else 0
  }
  names(dy) <- keys
  J <- dyn$jacobian_fn(dy, params, ss)
  J[!is.finite(J)] <- 0

  # ---- Partition columns into lag / contemporaneous / lead / shock ----
  AAlag  <- matrix(0, n_eq, n_endo)
  AA0    <- matrix(0, n_eq, n_endo)
  AAlead <- matrix(0, n_eq, n_endo)
  BB     <- matrix(0, n_eq, n_exo)
  for (col in seq_len(nrow(dcm))) {
    nm <- dcm$name[col]; ll <- dcm$lead_lag[col]
    vi <- match(nm, endo)
    if (!is.na(vi)) {
      if (ll == -1L) AAlag[, vi]  <- J[, col]
      else if (ll == 0L) AA0[, vi]   <- J[, col]
      else if (ll ==  1L) AAlead[, vi] <- J[, col]
      # higher-order leads/lags are not supported by the LQ discretion path
    } else {
      ei <- match(nm, exo)
      if (!is.na(ei)) BB[, ei] <- J[, col]
    }
  }
  colnames(AAlag) <- colnames(AA0) <- colnames(AAlead) <- endo
  if (n_exo > 0) colnames(BB) <- exo

  # ---- Instrument index ----
  instr_id <- match(policy_instrument, endo)
  if (is.na(instr_id)) {
    stop("policy_instrument '", policy_instrument,
         "' not found in endogenous variables.")
  }

  # ---- Make the system non-square (free instrument) ----
  # The Dennis (2007) recursion needs the reduced system: one equation fewer
  # than variables, with NO equation pinning the instrument (the planner sets
  # it optimally).  A free-instrument model already has n_eq = n_endo -
  # n_instr; a "square" model carries an explicit policy equation (e.g. a
  # Taylor rule i = phi_pi*pi + ...) that must be dropped before solving for
  # the OPTIMAL discretionary rule.
  n_instr <- length(instr_id)
  if (n_eq == n_endo - n_instr) {
    # Already non-square: nothing to drop.
  } else if (n_eq == n_endo) {
    # Square model: find and remove the equation that pins the instrument.
    # The policy equation is the one whose contemporaneous coefficient on the
    # instrument is largest relative to its coefficients on the other vars
    # (typically the row where solving the .mod isolates i on the LHS).
    inst_coef  <- abs(AA0[, instr_id])
    other_load <- rowSums(abs(AA0[, -instr_id, drop = FALSE])) +
                  rowSums(abs(AAlag)) + rowSums(abs(AAlead)) +
                  rowSums(abs(BB))
    # Prefer an equation of the form i = f(other vars): instrument coefficient
    # nonzero, and the instrument does NOT appear at lead/lag in that row.
    cand <- which(inst_coef > 1e-10 &
                  abs(AAlag[, instr_id]) < 1e-10 &
                  abs(AAlead[, instr_id]) < 1e-10)
    if (length(cand) == 0L) {
      stop(
        "Cannot identify the policy equation to drop for instrument '",
        policy_instrument, "' in this square model. For optimal discretionary ",
        "policy the instrument must be free; either write the model with the ",
        "instrument equation omitted (n_eq = n_endo - 1), or ensure exactly ",
        "one equation pins the instrument contemporaneously."
      )
    }
    # Among candidates, drop the one most dominated by the instrument
    # (largest instrument coefficient relative to the rest of the row).
    score    <- inst_coef[cand] / (other_load[cand] + 1e-12)
    drop_eq  <- cand[which.max(score)]
    if (verbose) {
      cat(sprintf("  Square model: dropping policy equation row %d (pins %s).\n",
                  drop_eq, policy_instrument))
    }
    keep <- setdiff(seq_len(n_eq), drop_eq)
    AAlag  <- AAlag[keep, , drop = FALSE]
    AA0    <- AA0[keep, , drop = FALSE]
    AAlead <- AAlead[keep, , drop = FALSE]
    BB     <- BB[keep, , drop = FALSE]
    n_eq   <- length(keep)
  } else {
    stop(
      "Unexpected equation count for discretionary policy: n_eq=", n_eq,
      ", n_endo=", n_endo, ", n_instruments=", n_instr,
      ". Expected either a free-instrument model (n_eq = n_endo - n_instr) ",
      "or a square model with an explicit policy equation (n_eq = n_endo)."
    )
  }

  # ---- Loss weight matrix bigw (Hessian of the planner objective) ----
  bigw <- .build_loss_weight_matrix(
    model, params, ss, endo, exo,
    loss_vars = loss_vars, loss_weights = loss_weights,
    instr_id = instr_id, control_penalty = control_penalty,
    verbose = verbose
  )

  # ---- Guard: degenerate (identically-zero) loss => fail loud ----
  # Discretion with W = 0 and Q = 0 gives F = 0 (do-nothing): a silent wrong
  # answer.  Keep H1's fail-loud behaviour for genuinely degenerate problems.
  if (all(abs(bigw) < 1e-300)) {
    stop(
      "Loss weight matrix is identically zero (loss_vars = (",
      paste(loss_vars, collapse = ", "), "), and the model carries no usable ",
      "planner_objective). The discretion problem is degenerate: every policy ",
      "yields the same (zero) loss, so the optimal rule is undefined. Supply ",
      "non-zero loss_weights or a planner_objective."
    )
  }

  list(AAlag = AAlag, AA0 = AA0, AAlead = AAlead, BB = BB, bigw = bigw,
       instr_id = instr_id, endo = endo, exo = exo,
       n_endo = n_endo, n_exo = n_exo, n_eq = n_eq)
}


#' Build the loss-weight (Hessian) matrix over all endogenous variables
#'
#' @return n_endo x n_endo symmetric matrix.
#' @noRd
.build_loss_weight_matrix <- function(model, params, ss, endo, exo,
                                      loss_vars, loss_weights,
                                      instr_id, control_penalty,
                                      verbose = FALSE) {
  n_endo <- length(endo)
  bigw <- matrix(0, n_endo, n_endo)
  rownames(bigw) <- colnames(bigw) <- endo

  obj_text <- model$planner_objective$text
  used_objective <- FALSE

  if (!is.null(obj_text) && nzchar(obj_text)) {
    # Hessian of the planner objective w.r.t. the contemporaneous endogenous
    # vector at SS, mirroring Dynare's W = build_two_dim_hessian(objective).
    Wobj <- tryCatch(
      .planner_objective_hessian(obj_text, model, params, ss, endo, exo),
      error = function(e) {
        if (verbose) {
          cat("  planner_objective Hessian failed (", conditionMessage(e),
              "); falling back to loss_vars.\n", sep = "")
        }
        NULL
      }
    )
    if (!is.null(Wobj) && any(abs(Wobj) > 0)) {
      bigw <- Wobj
      used_objective <- TRUE
      if (verbose) cat("  Loss matrix W built from planner_objective.\n")
    }
  }

  if (!used_objective) {
    # Build from loss_vars / loss_weights: weight w on v^2 => Hessian 2*w.
    for (v in loss_vars) {
      w <- loss_weights[v]
      if (is.na(w) || w == 0) next
      idx <- which(endo == v)
      if (length(idx) == 1L) bigw[idx, idx] <- bigw[idx, idx] + 2 * w
    }
    if (verbose) cat("  Loss matrix W built from loss_vars/loss_weights.\n")
  }

  # Control penalty on the instrument (only if it carries no own loss weight).
  if (!is.null(control_penalty) && control_penalty > 0 &&
      abs(bigw[instr_id, instr_id]) < 1e-300) {
    bigw[instr_id, instr_id] <- 2 * control_penalty
  }

  bigw
}


#' Hessian of the planner objective over current-period endogenous variables
#'
#' @return n_endo x n_endo symmetric matrix.
#' @noRd
.planner_objective_hessian <- function(obj_text, model, params, ss, endo, exo) {
  all_vn  <- c(endo, exo)
  obj_ast <- parse_expression(obj_text, var_names = all_vn,
                              param_names = model$param_names)
  n_endo  <- length(endo)

  # Evaluate the objective at a vector of current-period endo deviations.
  obj_fn <- function(v) {
    var_values <- numeric(0)
    for (i in seq_along(endo)) {
      nm <- endo[i]
      var_values[[paste0(nm, "__m1")]] <- v[i]
      var_values[[paste0(nm, "__0")]]  <- v[i]
      var_values[[paste0(nm, "__p1")]] <- v[i]
    }
    for (nm in exo) var_values[[paste0(nm, "__0")]] <- 0
    ast_eval(obj_ast, var_values = var_values, param_values = params,
             ss_values = ss)
  }

  # Central second differences.  Loss functions in LQ optimal-policy models
  # are quadratic, for which central 2nd differences are analytically exact;
  # the only error is floating-point round-off, which shrinks as the step h
  # grows, so a relatively large step (1e-3) gives the cleanest result.
  h  <- 1e-3
  W  <- matrix(0, n_endo, n_endo)
  x0 <- rep(0, n_endo)
  for (i in seq_len(n_endo)) {
    for (j in i:n_endo) {
      xpp <- x0; xpp[i] <- xpp[i] + h; xpp[j] <- xpp[j] + h
      xpm <- x0; xpm[i] <- xpm[i] + h; xpm[j] <- xpm[j] - h
      xmp <- x0; xmp[i] <- xmp[i] - h; xmp[j] <- xmp[j] + h
      xmm <- x0; xmm[i] <- xmm[i] - h; xmm[j] <- xmm[j] - h
      hij <- (obj_fn(xpp) - obj_fn(xpm) - obj_fn(xmp) + obj_fn(xmm)) /
             (4 * h * h)
      W[i, j] <- hij
      W[j, i] <- hij
    }
  }
  W[abs(W) < 1e-12] <- 0
  rownames(W) <- colnames(W) <- endo
  W
}


#' Dennis (2007) discretionary-policy fixed-point solver
#'
#' Solves the Markov-perfect (time-consistent, no-commitment) LQ problem
#'   min E_0 sum beta^t (y_t' W y_t + x_t' Q x_t)
#'   s.t. AAlag y_{t-1} + AA0 y_t + AAlead y_{t+1} + BB e_t = 0
#' where y_t is the full endogenous vector and x_t = y_t[instr_id] are the
#' instruments.  Returns the solution y_t = H y_{t-1} + G e_t (with the
#' instrument rows of H/G giving the optimal feedback / impact rule).
#'
#' Mirrors Dynare 7.1's discretionary_policy_engine.m (Dennis 2007,
#' Macroeconomic Dynamics 11, 31-55), including the auxiliary-variable
#' handling for instruments that appear with a lag.
#'
#' @noRd
.discretion_dennis_engine <- function(AAlag, AA0, AAlead, BB, bigw, instr_id,
                                      beta, max_iter = 3000L, tol = 1e-12,
                                      verbose = FALSE) {
  eq_nbr  <- nrow(AAlag)
  endo_nbr <- ncol(AAlag)
  exo_nbr <- ncol(BB)
  yidx    <- setdiff(seq_len(endo_nbr), instr_id)
  instr_nbr <- length(instr_id)

  # ---- Dennis matrices (A0 y = A1 y(-1) + A2 y(+1) + A3 x + A4 x(+1) + A5 e) ----
  A0 <- AA0[, yidx, drop = FALSE]
  A1 <- -AAlag[, yidx, drop = FALSE]
  A2 <- -AAlead[, yidx, drop = FALSE]
  A3 <- -AA0[, instr_id, drop = FALSE]
  A4 <- -AAlead[, instr_id, drop = FALSE]
  A5 <- -BB
  W  <- bigw[yidx, yidx, drop = FALSE]
  Q  <- bigw[instr_id, instr_id, drop = FALSE]

  # Auxiliary equations for instruments that appear with a lag.
  A6  <- -AAlag[, instr_id, drop = FALSE]
  aux <- apply(A6, 2, function(col) any(abs(col) > 0))
  n_aux <- sum(aux)

  ny <- eq_nbr
  m  <- eq_nbr + n_aux
  A00 <- matrix(0, m, m); A00[1:ny, 1:ny] <- A0
  if (n_aux > 0) A00[(ny + 1):m, (ny + 1):m] <- diag(n_aux)
  A11 <- matrix(0, m, m); A11[1:ny, 1:ny] <- A1
  if (n_aux > 0) A11[1:ny, (ny + 1):m] <- A6[, aux, drop = FALSE]
  A22 <- matrix(0, m, m); A22[1:ny, 1:ny] <- A2
  A33 <- matrix(0, m, instr_nbr); A33[1:ny, ] <- A3
  if (n_aux > 0) A33[(ny + 1):m, aux] <- diag(n_aux)
  A44 <- matrix(0, m, instr_nbr); A44[1:ny, ] <- A4
  A55 <- matrix(0, m, exo_nbr);   A55[1:ny, ] <- A5
  WW  <- matrix(0, m, m);         WW[1:ny, 1:ny] <- W
  endo_augm_id <- setdiff(seq_len(endo_nbr + n_aux), instr_id)

  # ---- Fixed-point iteration ----
  H1 <- matrix(0, m, m)
  F1 <- matrix(0, instr_nbr, m)
  H10 <- H1; F10 <- F1
  converged <- FALSE; retcode <- 0L; diff <- Inf
  A3DPD <- NULL; Dinv <- NULL; iter <- 0L

  for (iter in seq_len(max_iter)) {
    # P solves the Sylvester eq:  P = (WW + beta F1'Q F1) + beta H1' P H1
    P <- .sylvester_doubling(WW + beta * t(F1) %*% Q %*% F1,
                             beta * t(H1), H1, tol, max_iter)
    if (any(!is.finite(P))) { retcode <- 2L; break }

    D    <- A00 - A22 %*% H1 - A44 %*% F1
    Dinv <- tryCatch(solve(D), error = function(e) NULL)
    if (is.null(Dinv)) { retcode <- 2L; break }
    A3DPD <- t(A33) %*% t(Dinv) %*% P %*% Dinv
    F1 <- -solve(Q + A3DPD %*% A33) %*% (A3DPD %*% A11)
    H1 <- Dinv %*% (A11 + A33 %*% F1)

    diff <- max(abs(rbind(H1, F1) - rbind(H10, F10)))
    if (!is.finite(diff)) { retcode <- 3L; break }
    if (diff < tol) { converged <- TRUE; break }
    H10 <- H1; F10 <- F1
  }

  if (!converged) {
    return(list(H = NULL, G = NULL, P = P, converged = FALSE,
                retcode = if (retcode == 0L) 1L else retcode,
                iterations = iter, diff = diff,
                A_state = NULL, B_state = NULL))
  }

  # ---- Impact (shock) response ----
  F2 <- -solve(Q + A3DPD %*% A33) %*% (A3DPD %*% A55)
  H2 <- Dinv %*% (A55 + A33 %*% F2)

  # ---- Re-assemble full H (endo x endo) and G (endo x exo) ----
  H <- matrix(0, endo_nbr + n_aux, endo_nbr + n_aux)
  G <- matrix(0, endo_nbr + n_aux, exo_nbr)
  H[endo_augm_id, endo_augm_id] <- H1
  H[instr_id, endo_augm_id]     <- F1
  G[endo_augm_id, ]             <- H2
  G[instr_id, ]                 <- F2

  if (n_aux > 0) {
    # Map auxiliary-variable columns back onto the lagged instrument.
    aux_cols <- (endo_nbr + n_aux) - (rev(seq_len(n_aux)) - 1L)
    H[, instr_id[aux]] <- H[, aux_cols, drop = FALSE]
  }
  H <- H[1:endo_nbr, 1:endo_nbr, drop = FALSE]
  G <- G[1:endo_nbr, , drop = FALSE]

  # State-space (A,B) over the predetermined columns of H, for back-compat.
  state_cols <- which(apply(abs(H), 2, max) > 1e-12)
  A_state <- H[state_cols, state_cols, drop = FALSE]
  B_state <- G[state_cols, , drop = FALSE]

  list(H = H, G = G, P = P, converged = TRUE, retcode = 0L,
       iterations = iter, diff = diff,
       A_state = A_state, B_state = B_state)
}


#' Sylvester doubling:  v = d + g v h  (Dynare's SylvesterDoubling)
#' @noRd
.sylvester_doubling <- function(d, g, h, tol, max_iter) {
  v <- d
  for (i in seq_len(max_iter)) {
    vadd <- g %*% v %*% h
    v <- v + vadd
    nv <- norm(v, "1")
    if (norm(vadd, "1") <= tol * (if (nv > 0) nv else 1)) break
    g <- g %*% g
    h <- h %*% h
  }
  v
}


#' Build a closed-loop DecisionRules object from the discretion solution
#'
#' @param H n_endo x n_endo state-transition matrix (y_t = H y_{t-1}).
#' @param G n_endo x n_exo impact matrix (y_t = G e_t).
#' @param state_idx Integer indices of the endogenous columns of H that carry
#'   feedback (the effective state variables).
#' @noRd
.discretion_build_dr <- function(model, ss, endo, exo, H, G, state_idx) {
  n_endo <- length(endo)
  n_exo  <- length(exo)
  if (length(state_idx) == 0L) state_idx <- integer(0)

  ghx <- H[, state_idx, drop = FALSE]   # n_endo x n_state ; y_t = ghx y_{t-1}[state]
  ghu <- G                              # n_endo x n_exo

  ys <- if (length(ss) == n_endo) ss else setNames(rep(0, n_endo), endo)

  ghx_state <- ghx[state_idx, , drop = FALSE]
  ev <- if (length(state_idx) > 0L) {
    tryCatch(eigen(ghx_state, only.values = TRUE)$values,
             error = function(e) rep(NA_complex_, length(state_idx)))
  } else complex(0)
  n_unstable <- if (length(ev) > 0L) sum(Mod(ev) > 1 + 1e-8) else 0L

  dr <- list(
    ghx          = ghx,
    ghu          = ghu,
    ys           = ys,
    endo_names   = endo,
    exo_names    = exo,
    state_vars   = endo[state_idx],
    state_idx    = state_idx,
    n_state      = length(state_idx),
    n_stable     = length(state_idx) - n_unstable,
    n_exo        = n_exo,
    eigenvalues  = ev,
    n_unstable   = n_unstable,
    bk_satisfied = (n_unstable == 0L)
  )
  class(dr) <- "DecisionRules"
  dr
}

#' Solve the discrete algebraic Riccati equation via fixed-point iteration
#'
#' P = Q + β A' P A - β² A' P B (R + β B' P B)^{-1} B' P A
#' F = (R + β B' P B)^{-1} β B' P A
#'
#' @noRd
.solve_dare_fixed_point <- function(A, B, Q, R, beta,
                                     max_iter = 1000L, tol = 1e-12,
                                     verbose = FALSE) {
  n <- nrow(A)
  P <- Q  # Initial guess

  converged <- FALSE
  diff <- Inf

  for (iter in seq_len(max_iter)) {
    # Compute B' P B and B' P A
    BPB <- t(B) %*% P %*% B
    BPA <- t(B) %*% P %*% A

    # F = (R + β B' P B)^{-1} β B' P A
    R_inv <- solve(R + beta * BPB)
    F_new <- R_inv %*% (beta * BPA)  # (1 × n)

    # P_new = Q + β A' P A - β² A' P B (R + β B' P B)^{-1} B' P A
    #       = Q + β A' P A - β A' P B * F_new
    #       = Q + β (A - B F_new)' P (A - B F_new) + F_new' R F_new
    Acl <- A - B %*% F_new
    P_new <- Q + t(F_new) %*% R %*% F_new + beta * t(Acl) %*% P %*% Acl

    # Check convergence
    diff <- max(abs(P_new - P))

    if (diff < tol) {
      converged <- TRUE
      if (verbose && iter > 1) {
        cat(sprintf("    DARE converged in %d iterations, |DeltaP| = %.3e\n",
                    iter, diff))
      }
      P <- P_new
      break
    }

    P <- P_new
  }

  if (!converged && verbose) {
    cat(sprintf("    DARE did NOT converge in %d iterations, |DeltaP| = %.3e\n",
                max_iter, diff))
  }

  # Final F
  BPB <- t(B) %*% P %*% B
  BPA <- t(B) %*% P %*% A
  R_inv <- solve(R + beta * BPB + 1e-10 * diag(nrow(BPB)))
  F_final <- R_inv %*% (beta * BPA)

  list(F = F_final, P = P, converged = converged,
       iterations = iter, diff = diff)
}


#' Compare discretionary policy against Ramsey-optimal commitment
#'
#' @noRd
.disc_compare_ramsey <- function(dr_disc, moments_disc,
                                  loss_vars, loss_weights,
                                  ramsey_result, discount) {

  # Extract Ramsey DR
  ramsey_dr <- ramsey_result$ramsey_dr$ramsey_dr

  if (is.null(ramsey_dr) || is.null(ramsey_dr$ghx)) {
    return(list(
      ramsey_loss   = NA_real_,
      disc_loss     = NA_real_,
      welfare_gap   = NA_real_,
      message       = "Ramsey DR not available"
    ))
  }

  # Ramsey moments
  ramsey_moments <- compute_moments(ramsey_dr, ramsey_result$augmented_model)

  ramsey_loss <- if (!is.null(ramsey_moments)) {
    rv <- ramsey_moments$std_dev^2
    sum(loss_weights * rv[loss_vars], na.rm = TRUE)
  } else NA_real_

  # Discretionary loss
  disc_loss <- if (!is.null(moments_disc)) {
    dv <- moments_disc$std_dev^2
    sum(loss_weights * dv[loss_vars], na.rm = TRUE)
  } else NA_real_

  list(
    ramsey_loss = ramsey_loss,
    disc_loss   = disc_loss,
    welfare_gap = disc_loss - ramsey_loss  # positive = Ramsey better
  )
}


# ==========================================================================
# S3 methods
# ==========================================================================

#' @export
print.dynhr_discretionary_result <- function(x, ...) {
  cat(sprintf("\n<dynhr_discretionary_result>\n"))
  cat(sprintf("  Policy instrument: %s\n", x$policy_instrument))
  cat(sprintf("  State variables:   %d\n", length(x$state_vars)))
  cat(sprintf("  Discount factor:   %.4f\n", x$discount))
  cat(sprintf("  Converged:         %s (%d iters, |DeltaP| = %.3e)\n",
              x$converged, x$iterations, x$diff))

  cat("  Optimal policy rule:\n")
  cat(sprintf("    %s = ", x$policy_instrument))
  terms <- sapply(seq_along(x$state_vars), function(j) {
    coef <- -x$F[j]
    if (abs(coef) < 1e-14) return(NULL)
    sprintf("%.4f * %s", coef, x$state_vars[j])
  })
  terms <- terms[!sapply(terms, is.null)]
  cat(paste(terms, collapse = " + "), "\n")

  if (!is.null(x$moments)) {
    cat("\n  Unconditional std devs:\n")
    top_vars <- head(names(sort(x$moments$std_dev, decreasing = TRUE)), 6)
    for (v in top_vars) {
      cat(sprintf("    %-15s = %.6f\n", v, x$moments$std_dev[v]))
    }
  }

  if (!is.null(x$ramsey_comparison)) {
    cat(sprintf("\n  Ramsey comparison:\n"))
    cat(sprintf("    Ramsey loss:   %.6f\n", x$ramsey_comparison$ramsey_loss))
    cat(sprintf("    Disc. loss:    %.6f\n", x$ramsey_comparison$disc_loss))
    cat(sprintf("    Welfare gap:   %.6f (positive = Ramsey better)\n",
                x$ramsey_comparison$welfare_gap))
  }

  invisible(x)
}
