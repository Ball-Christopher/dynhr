## R/ramsey-obc.R
## --------------------------------------------------------------------------
## Phase G1 — Ramsey Optimal Policy with Occasionally Binding Constraints
##
## Provides two complementary approaches:
##
## G1a: Perfect-foresight Ramsey with OBC (ramsey_obc_pf)
##   Combines the Ramsey augmented system with the perfect-foresight Newton
##   path solver (pf_newton_solve).  Solves the deterministic transition path
##   under Ramsey commitment when an OBC (e.g. ZLB) may bind.  Uses the
##   compiled augmented model (original vars + Lagrange multipliers + FOCs)
##   and applies the OBC constraint to the relevant original equation.
##
## G1b: Piecewise-linear optimal commitment (ramsey_obc_pwlinear)
##   Extends the OccBin / Harrison-Waldron (2021) piecewise-linear algorithm
##   to the augmented Ramsey system.  Builds slack and binding regime policy
##   matrices for the full augmented system, then runs the iterative
##   guess-and-verify regime search.  Enables OBC-constrained Ramsey policy
##   within the linear(-quadratic) framework.
##
## References:
##   Eggertsson & Woodford (2003), "The Zero Bound on Interest Rates and
##     Optimal Monetary Policy", Brookings Papers on Economic Activity.
##   Harrison & Waldron (2021), "Optimal Monetary Policy with Occasionally
##     Binding Constraints", J. Econ. Dyn. Control.
##   Guerrieri & Iacoviello (2015), "OccBin: A toolkit for solving dynamic
##     models with occasionally binding constraints easily."
##   Bodenstein & Guerrieri (2019), "Nash–Ramsey toolbox."
## --------------------------------------------------------------------------


# ==========================================================================
# G1a: Perfect-foresight Ramsey with OBC
# ==========================================================================

#' Perfect-foresight Ramsey optimal policy with occasionally binding constraints
#'
#' Solves a deterministic transition path under Ramsey commitment subject to
#' occasionally binding constraints (e.g. ZLB).  Uses the augmented-system
#' Ramsey model (via \code{\link{ramsey_model}}) and the perfect-foresight
#' Newton path solver (\code{\link{pf_newton_solve}}) on the augmented system.
#'
#' The OBC constraint is applied to the specific equation in the AUGMENTED
#' system.  Since the augmented system preserves the original equations in
#' their original positions (1..n) and appends the FOCs as equations (n+1..2n),
#' OBC specs that target an original equation (e.g. the Taylor rule for ZLB)
#' operate on the same equation index in the augmented model.
#'
#' @param model             A dynhr_mod object.
#' @param planner_objective Character string: the planner objective expression.
#'   If NULL, uses \code{model$planner_objective$text}.
#' @param shock_path        T x n_exo matrix of structural shock values for each
#'   period. Column names must match \code{model$varexo_names}.
#' @param obc_specs         List of OBC specs as returned by
#'   \code{\link{obc_parse_tags}} or \code{obc_collect_specs}.  Each spec must
#'   have \code{$eq_idx}, \code{$var_idx}, \code{$var_name}, \code{$op}, \code{$bound}.
#'   If NULL, tries to parse MCP tags from the model.
#' @param params            Named numeric parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param order             Perturbation order for Ramsey solution (default 1).
#' @param T_horizon         Integer: number of periods for the perfect-foresight
#'   path (default 40).
#' @param max_iter          Maximum Newton iterations per regime (default 50).
#' @param tol               Newton convergence tolerance on max|R| (default 1e-8).
#' @param max_regime_iter   Maximum regime switches in active-set iteration (default 30).
#' @param step_size         Newton step-length (default 1.0).
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_ramsey_obc_pf} containing:
#'   \item{Y}{T x n_endo_orig matrix: paths of original endogenous variables.}
#'   \item{Y_full}{T x (n_endo_orig + n_mult) matrix: paths including multipliers.}
#'   \item{regime}{n_spec x T logical matrix (TRUE = binding).}
#'   \item{ramsey_result}{The underlying \code{dynhr_ramsey_result2} object.}
#'   \item{converged}{Logical: whether Newton converged.}
#'   \item{welfare}{List with steady and path-based welfare values.}
#'   \item{meta}{Metadata (horizon, order, timing).}
#'
#' @references
#'   Bodenstein, M., & Guerrieri, L. (2011). The welfare-optimal degree of
#'     central bank transparency. \emph{Journal of Money, Credit and Banking},
#'     43(5), 857-889.
#'   Harrison, R., & Waldron, M. (2018). Optimal policy with occasionally
#'     binding constraints: Piecewise-linear commitment. \emph{Working Paper}.
#' @export
ramsey_obc_pf <- function(model,
                           planner_objective = NULL,
                           shock_path = NULL,
                           obc_specs = NULL,
                           params = NULL,
                           order = 1L,
                           T_horizon = 40L,
                           max_iter = 50L,
                           tol = 1e-8,
                           max_regime_iter = 30L,
                           step_size = 1.0,
                           verbose = FALSE,
                           ...) {
  # ---- 1. Validate and defaults ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  if (is.null(params)) params <- model$param_values
  T_horizon <- as.integer(T_horizon)

  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective provided.")
  }

  # ---- 2. Get or create OBC specs ----
  if (is.null(obc_specs)) {
    obc_specs <- obc_parse_tags(model)
  }

  # ---- 3. Build the Ramsey augmented model ----
  if (verbose) cat("[ramsey_obc_pf] Building augmented Ramsey model (order=1)...\n")

  ramsey_result <- ramsey_model(
    model             = model,
    params            = params,
    planner_objective = obj_text,
    order             = as.integer(order),
    method            = "augmented",
    verbose           = verbose,
    ...
  )

  aug_model <- ramsey_result$augmented_model
  aug_ss    <- ramsey_result$ramsey_steady$values
  n_orig    <- length(model$var_names)
  n_aug     <- length(aug_model$var_names)
  n_mult    <- n_aug - n_orig

  # Compile the augmented model if not already compiled
  aug_compiled <- compile_model(aug_model, verbose = FALSE)

  # ---- 4. Build shock path ----
  if (is.null(shock_path)) {
    # Default: one-time shock at period 1, zero thereafter
    shock_path <- matrix(0, nrow = T_horizon, ncol = length(model$varexo_names))
    colnames(shock_path) <- model$varexo_names
  } else {
    if (!is.matrix(shock_path)) shock_path <- as.matrix(shock_path)
    if (nrow(shock_path) < T_horizon) {
      # Extend with zeros
      ext <- matrix(0, nrow = T_horizon - nrow(shock_path), ncol = ncol(shock_path))
      colnames(ext) <- colnames(shock_path)
      shock_path <- rbind(shock_path, ext)
    }
    # Ensure all exo names present
    full_shocks <- matrix(0, nrow = T_horizon, ncol = length(model$varexo_names))
    colnames(full_shocks) <- model$varexo_names
    for (nm in intersect(colnames(shock_path), model$varexo_names)) {
      full_shocks[, nm] <- shock_path[, nm]
    }
    shock_path <- full_shocks
  }

  # ---- 5. Map OBC specs to augmented model ----
  # OBC specs reference original eq/var indices. In the augmented model:
  #   - Original equations 1..n_orig are unchanged
  #   - Original variable indices 1..n_orig are unchanged
  # So the specs are already valid for the augmented system.
  aug_specs <- obc_specs  # Same eq_idx and var_idx for original equations/vars

  # ---- 6. Run perfect-foresight solver on augmented system ----
  if (verbose) cat("[ramsey_obc_pf] Running perfect-foresight Newton solver",
                   sprintf("(%d periods, %d OBC specs)...\n", T_horizon, length(aug_specs)))

  pf_result <- pf_newton_solve(
    compiled        = aug_compiled,
    y0              = aug_ss,
    y_ss            = aug_ss,
    shock_path      = shock_path,
    params          = params,
    obc_specs       = aug_specs,
    max_iter        = max_iter,
    tol             = tol,
    max_regime_iter = max_regime_iter,
    step_size       = step_size
  )

  # ---- 7. Extract results ----
  Y_full <- pf_result$Y   # T x n_aug matrix (augmented variables)
  regime <- pf_result$regime

  # Extract original variables
  orig_idx <- seq_len(n_orig)
  Y_orig <- Y_full[, orig_idx, drop = FALSE]
  colnames(Y_orig) <- model$var_names

  # ---- 8. Compute welfare along the path ----
  discount <- .get_discount(params) %||% 0.99
  welfare_steady <- .eval_planner_ast(
    parse_expression(obj_text,
      var_names = c(model$var_names, model$varexo_names),
      param_names = model$param_names),
    setNames(as.numeric(model$var_names), model$var_names),
    params, aug_ss[seq_len(n_orig)]
  ) / max(1e-8, 1 - discount)

  # Path welfare: discounted sum of objective along the path
  obj_ast <- parse_expression(obj_text,
    var_names = c(model$var_names, model$varexo_names),
    param_names = model$param_names)

  obj_t <- apply(Y_orig, 1L, function(row) {
    vals <- setNames(as.numeric(row), model$var_names)
    .eval_planner_ast(obj_ast, vals, params, aug_ss[seq_len(n_orig)])
  })

  disc_weights <- discount ^ (seq_len(T_horizon) - 1)
  welfare_path <- sum(obj_t * disc_weights, na.rm = TRUE)

  # ---- 9. Build result ----
  result <- list(
    Y             = Y_orig,
    Y_full        = Y_full,
    regime        = regime,
    ramsey_result = ramsey_result,
    converged     = isTRUE(pf_result$converged),
    welfare       = list(
      steady_value = welfare_steady,
      path_value   = welfare_path,
      discount     = discount
    ),
    meta = list(
      method        = "perfect_foresight",
      order         = as.integer(order),
      T_horizon     = T_horizon,
      n_orig_vars   = n_orig,
      n_multipliers = n_mult,
      n_obc_specs   = length(obc_specs),
      converged     = isTRUE(pf_result$converged),
      n_newton_iter = pf_result$n_iter %||% NA_integer_,
      n_regime_iter = NA_integer_,
      timestamp     = Sys.time()
    )
  )
  class(result) <- c("dynhr_ramsey_obc_pf", "list")
  result
}


#' @export
print.dynhr_ramsey_obc_pf <- function(x, ...) {
  cat("\n<dynhr_ramsey_obc_pf>  [Perfect-foresight Ramsey + OBC]\n")
  cat(sprintf("  Converged:          %s\n", if (isTRUE(x$converged)) "YES" else "NO"))
  cat(sprintf("  Horizon:            %d periods\n", x$meta$T_horizon))
  cat(sprintf("  Original vars:      %d\n", x$meta$n_orig_vars))
  cat(sprintf("  Multipliers:        %d\n", x$meta$n_multipliers))
  cat(sprintf("  OBC constraints:    %d\n", x$meta$n_obc_specs))
  cat(sprintf("  Path welfare:       %.6f\n", x$welfare$path_value))
  cat(sprintf("  Steady welfare:     %.6f\n", x$welfare$steady_value))
  cat(sprintf("  Discount factor:    %.4f\n", x$welfare$discount))

  if (!is.null(x$regime) && ncol(x$regime) > 0) {
    n_bind <- sum(x$regime, na.rm = TRUE)
    cat(sprintf("  Binding periods:    %d / %d\n", n_bind, length(x$regime)))
  }
  invisible(x)
}


# ==========================================================================
# G1b: Piecewise-linear Ramsey optimal commitment (Harrison-Waldron)
# ==========================================================================

#' Piecewise-linear Ramsey optimal policy with OBC (OccBin-style)
#'
#' Extends the OccBin piecewise-linear algorithm to the augmented Ramsey
#' system.  Builds slack and binding regime policy matrices for the full
#' augmented system (original variables + Lagrange multipliers), then runs
#' the iterative guess-and-verify regime search to find the consistent
#' OBC-constrained Ramsey equilibrium.
#'
#' The algorithm:
#'   1. Build the augmented Ramsey system via \code{\link{ramsey_model}}.
#'   2. Extract the linear system matrices from the augmented compiled model.
#'   3. For each OBC spec, construct the binding-regime system by replacing
#'      the constrained equation with the binding constraint.
#'   4. Solve the binding-regime system for each candidate regime combination.
#'   5. Run the OccBin guess-and-verify (or Boehl complementarity) iteration
#'      to find the consistent regime path.
#'   6. Return slack and binding decision rules for the augmented system.
#'
#' @param model             A dynhr_mod object.
#' @param planner_objective Character string: the planner objective expression.
#' @param obc_specs         List of OBC specs from \code{\link{obc_parse_tags}}.
#' @param params            Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param method            Regime search method: \code{"guess_verify"} (OccBin,
#'   default) or \code{"boehl"} (complementarity iteration).
#' @param shock_std         Standard deviation for the shock used in regime path
#'   computation (default 1.0, i.e. one-standard-deviation impulse).
#' @param max_iter          Maximum regime search iterations (default 100).
#' @param tol               Tolerance for regime convergence (default 1e-6).
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_ramsey_obc_pwl} containing:
#'   \item{ramsey_result}{The underlying \code{dynhr_ramsey_result2} object.}
#'   \item{dr_slack}{Slack-regime DecisionRules for the augmented system.}
#'   \item{dr_bind}{Binding-regime DecisionRules for each unique regime.}
#'   \item{regime_path}{Integer vector: optimal regime path (bitfield encoding).}
#'   \item{irf}{IRF matrix (T x n_aug) under the Ramsey OBC policy.}
#'   \item{meta}{Metadata.}
#' @export
ramsey_obc_pwlinear <- function(model,
                                 planner_objective = NULL,
                                 obc_specs = NULL,
                                 params = NULL,
                                 method = c("guess_verify", "boehl"),
                                 shock_std = 1.0,
                                 max_iter = 100L,
                                 tol = 1e-6,
                                 verbose = FALSE,
                                 ...) {
  # ---- 1. Validate ----
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  method <- match.arg(method)

  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective provided.")
  }
  if (is.null(params)) params <- model$param_values

  # ---- 2. Get OBC specs ----
  if (is.null(obc_specs)) {
    obc_specs <- obc_parse_tags(model)
  }
  n_spec <- length(obc_specs)

  # ---- 3. Build augmented Ramsey system (order 1) ----
  if (verbose) cat("[ramsey_obc_pwlinear] Building augmented Ramsey system...\n")

  ramsey_result <- ramsey_model(
    model             = model,
    params            = params,
    planner_objective = obj_text,
    order             = 1L,
    method            = "augmented",
    verbose           = verbose,
    ...
  )

  aug_model    <- ramsey_result$augmented_model
  aug_compiled <- compile_model(aug_model, verbose = FALSE)
  aug_ss       <- ramsey_result$ramsey_steady$values
  n_aug        <- length(aug_model$var_names)
  n_orig       <- length(model$var_names)
  n_mult       <- n_aug - n_orig

  # ---- 4. Extract augmented system matrices ----
  if (verbose) cat("[ramsey_obc_pwlinear] Extracting augmented system matrices...\n")

  sys <- extract_system_matrices(aug_compiled, aug_ss, params)
  # Augment sys with n_endo, endo_names, etc. needed by obc helpers
  sys$n_endo     <- n_aug
  sys$endo_names <- aug_model$var_names
  sys$exo_names  <- aug_model$varexo_names
  sys$n_exo      <- length(sys$exo_names)

  # ---- 5. Solve slack-regime perturbation (order 1) ----
  if (verbose) cat("[ramsey_obc_pwlinear] Solving slack-regime perturbation...\n")

  dr_slack <- solve_perturbation(
    aug_model, aug_compiled, aug_ss, params,
    order = 1L, verbose = FALSE
  )

  # ---- 6. Build binding-regime system ----
  # Map OBC specs to augmented model indices
  aug_specs <- lapply(obc_specs, function(s) {
    list(
      eq_idx   = s$eq_idx,    # Same equation index in augmented system
      var_idx  = s$var_idx,   # Same variable index in augmented system
      var_name = s$var_name,
      op       = s$op,
      bound    = s$bound
    )
  })

  # Build the binding-regime policy using the OccBin terminal-substitution
  # approach: obc_solve_binding replaces the constrained equations with the
  # binding constraints and substitutes the slack policy as the terminal
  # condition for the next period (Guerrieri-Iacoviello 2015).
  obs_idx <- seq_len(n_aug)
  state_idx <- dr_slack$state_idx
  n_state   <- length(state_idx)

  dr_bind <- obc_solve_binding(sys, dr_slack, aug_specs, obs_idx)

  # ---- 7. Run regime search ----
  if (verbose) cat("[ramsey_obc_pwlinear] Running regime search (method = '", method, "')...\n", sep = "")

  # Build shock sequence: one std dev impulse in the first period
  n_T <- 40L  # Default horizon
  shock_seq <- matrix(0, nrow = n_T, ncol = length(aug_model$varexo_names))
  colnames(shock_seq) <- aug_model$varexo_names
  if (ncol(shock_seq) >= 1) {
    shock_seq[1, 1] <- shock_std
  }

  if (method == "boehl") {
    # Use Boehl complementarity iteration
    regime_search <- boehl_solve_regime_path(
      shock_seq    = shock_seq,
      dr_slack     = dr_slack,
      sys          = sys,
      specs        = aug_specs,
      obs_idx      = obs_idx,
      max_iter     = max_iter,
      tol          = tol
    )
    regime_path <- regime_search$regime_path
    irf_paths   <- regime_search$paths
  } else {
    # Use OccBin guess-and-verify Kalman-filter approach
    # Generate Y by simulating under all-slack assumption first
    Y_init <- boehl_simulate(
      shock_seq    = shock_seq,
      dr_slack     = dr_slack,
      regime_cache = .build_regime_cache(dr_slack, n_aug, state_idx, obs_idx),
      regime_path  = integer(n_T)
    )$paths

    regime_search <- obc_guess_verify(
      Y         = Y_init,
      dr_slack  = dr_slack,
      sys       = sys,
      obs_idx   = obs_idx,
      model     = aug_model,
      params    = params,
      obs_vars  = aug_model$var_names,
      specs     = aug_specs,
      max_iter  = max_iter
    )
    regime_path <- regime_search$regime_path

    # Compute IRF by simulating under the found regime path
    sim <- boehl_simulate(
      shock_seq    = shock_seq,
      dr_slack     = dr_slack,
      regime_cache = regime_search$regime_cache,
      regime_path  = regime_path
    )
    irf_paths <- sim$paths
  }

  # ---- 8. Compute welfare ----
  discount <- .get_discount(params) %||% 0.99
  welfare_steady <- ramsey_result$welfare$steady_value

  # ---- 9. Build result ----
  result <- list(
    ramsey_result = ramsey_result,
    dr_slack      = dr_slack,
    dr_bind       = dr_bind,
    regime_path   = regime_path,
    irf           = irf_paths,
    welfare       = list(
      steady_value = welfare_steady,
      discount     = discount
    ),
    meta = list(
      method        = method,
      n_orig_vars   = n_orig,
      n_multipliers = n_mult,
      n_aug_vars    = n_aug,
      n_obc_specs   = n_spec,
      n_regime_iter = regime_search$n_iter %||% NA_integer_,
      timestamp     = Sys.time()
    )
  )
  class(result) <- c("dynhr_ramsey_obc_pwl", "list")
  result
}


#' @export
print.dynhr_ramsey_obc_pwl <- function(x, ...) {
  cat("\n<dynhr_ramsey_obc_pwl>  [Piecewise-linear Ramsey + OBC]\n")
  cat(sprintf("  Method:             %s\n", x$meta$method))
  cat(sprintf("  Original vars:      %d\n", x$meta$n_orig_vars))
  cat(sprintf("  Multipliers:        %d\n", x$meta$n_multipliers))
  cat(sprintf("  Augmented vars:     %d\n", x$meta$n_aug_vars))
  cat(sprintf("  OBC constraints:    %d\n", x$meta$n_obc_specs))
  cat(sprintf("  Welfare (steady):   %.6f\n", x$welfare$steady_value))

  if (!is.null(x$regime_path) && length(x$regime_path) > 0) {
    n_bind <- sum(x$regime_path != 0)
    cat(sprintf("  Binding periods:    %d / %d\n", n_bind, length(x$regime_path)))
  }
  if (!is.null(x$irf)) {
    cat(sprintf("  IRF matrix:         %d x %d\n", nrow(x$irf), ncol(x$irf)))
  }
  invisible(x)
}


#' Compare Ramsey policy welfare across OBC and unconstrained cases
#'
#' Runs both unconstrained Ramsey and OBC-constrained Ramsey, then reports
#' the welfare cost of the OBC constraint.
#'
#' @param model             dynhr_mod object.
#' @param planner_objective Planner objective expression.
#' @param obc_specs         OBC specs (or NULL to parse from model).
#' @param shock_std         Shock standard deviation for IRF (default 1.0).
#' @param params            Parameter vector.
#' @param verbose           Print progress.
#' @return List with welfare comparison and both result objects.
#' @export
ramsey_obc_welfare_cost <- function(model,
                                    planner_objective = NULL,
                                    obc_specs = NULL,
                                    shock_std = 1.0,
                                    params = NULL,
                                    verbose = FALSE) {
  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text)) stop("No planner objective provided.")
  if (is.null(params)) params <- model$param_values
  discount <- .get_discount(params) %||% 0.99

  # Unconstrained Ramsey
  if (verbose) cat("[ramsey_obc_welfare_cost] Solving unconstrained Ramsey...\n")
  ramsey_uncon <- ramsey_model(model, params = params,
    planner_objective = obj_text, order = 1L, verbose = verbose)

  # OBC-constrained Ramsey (perfect-foresight)
  if (verbose) cat("[ramsey_obc_welfare_cost] Solving OBC-constrained Ramsey (PF)...\n")
  ramsey_obc <- ramsey_obc_pf(model,
    planner_objective = obj_text, obc_specs = obc_specs,
    params = params, verbose = verbose)

  welfare_uncon <- .extract_welfare(ramsey_uncon, "steady")
  welfare_obc   <- ramsey_obc$welfare$steady_value

  consumption_equiv <- if (is.finite(welfare_uncon - welfare_obc) &&
                           abs(welfare_uncon) > 1e-12) {
    # Approximate CE difference using log-linear approximation
    (1 - discount) * (welfare_uncon - welfare_obc)
  } else NA_real_

  result <- list(
    ramsey_unconstrained = ramsey_uncon,
    ramsey_obc           = ramsey_obc,
    welfare = list(
      unconstrained = welfare_uncon,
      obc_constrained = welfare_obc,
      gap = welfare_uncon - welfare_obc,
      consumption_equivalent = consumption_equiv,
      discount = discount
    )
  )
  class(result) <- c("dynhr_ramsey_obc_cost", "list")
  result
}


#' @export
print.dynhr_ramsey_obc_cost <- function(x, ...) {
  cat("\n<dynhr_ramsey_obc_cost>  [Welfare cost of OBC constraint]\n")
  cat(sprintf("  Welfare unconstrained:  %.6f\n", x$welfare$unconstrained))
  cat(sprintf("  Welfare OBC-constrained: %.6f\n", x$welfare$obc_constrained))
  cat(sprintf("  Welfare gap:            %.6f\n", x$welfare$gap))
  cat(sprintf("  Consumption-equiv diff: %.6f\n", x$welfare$consumption_equivalent))
  cat(sprintf("  Discount factor:        %.4f\n", x$welfare$discount))
  invisible(x)
}


#' Build a minimal regime cache containing only the slack policy
#'
#' Internal helper for piecewise-linear Ramsey.  Creates an R environment
#' with the slack (regime 0) policy pre-seeded, suitable for initial
#' simulation passes.
#'
#' @param dr_slack  Slack-regime DecisionRules.
#' @param n_aug     Number of augmented endogenous variables.
#' @param state_idx Integer vector: state variable indices.
#' @param obs_idx   Integer vector: observable variable indices.
#' @return An R environment (hash map) with key "0".
#' @noRd
.build_regime_cache <- function(dr_slack, n_aug, state_idx, obs_idx) {
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  n_state <- length(state_idx)
  assign("0", list(
    dr      = dr_slack,
    c_full  = numeric(n_aug),
    c_state = numeric(n_state),
    c_obs   = numeric(length(obs_idx)),
    TT      = dr_slack$ghx[state_idx, , drop = FALSE],
    RR      = dr_slack$ghu[state_idx, , drop = FALSE],
    ZZ      = dr_slack$ghx,
    DD      = dr_slack$ghu
  ), envir = cache)
  cache
}
