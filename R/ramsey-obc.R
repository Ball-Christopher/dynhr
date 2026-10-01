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
#'   If NULL, the specs come from the model: its MCP tags / complementarity
#'   conditions (\code{\link{obc_parse_tags}}) and its
#'   \code{ramsey_constraints} block.  A \code{ramsey_constraints} bound on
#'   \code{x} (in levels) is complementary to the Ramsey first-order
#'   condition with respect to \code{x}, as in Dynare 7's
#'   \code{perfect_foresight_solver(lmmcp)}: in a binding period the FOC is
#'   replaced by \code{x = bound}.
#' @param params            Named numeric parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param order             Perturbation order for Ramsey solution (default 1).
#' @param T_horizon         Integer: number of periods for the perfect-foresight
#'   path (default 40).
#' @param max_iter          Maximum Newton iterations per regime (default 50).
#' @param tol               Tolerance (default 1e-8) passed to
#'   \code{\link{pf_newton_solve}}: its Newton stop (every stacked
#'   equation's residual within \code{tol} times the magnitude of its terms,
#'   and \eqn{\max|R| <} \code{tol}) and the relative round-off band of its
#'   regime decisions (bound violations relative to \eqn{|x| + |b|}, the
#'   release test relative to the relaxed equation's term magnitudes), so the
#'   regimes do not depend on the units of the model.
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
  # With obc_specs = NULL the specs come from the .mod: MCP tags /
  # complementarity conditions on the model equations, plus the
  # ramsey_constraints block (mapped onto the augmented system below).
  use_rc <- is.null(obc_specs) && is.data.frame(model$ramsey_constraints) &&
    nrow(model$ramsey_constraints) > 0L
  if (is.null(obc_specs)) {
    obc_specs <- if (!use_rc || .model_has_obc_tags(model))
      obc_parse_tags(model) else list()
  }

  # ---- 3. Build the Ramsey augmented model ----
  if (verbose) .dynhr_cat("[ramsey_obc_pf] Building augmented Ramsey model (order=1)...\n")

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
  if (use_rc) {
    aug_specs <- c(aug_specs,
                   .ramsey_constraint_specs(model, ramsey_result, params))
    obc_specs <- aug_specs
  }

  # ---- 6. Run perfect-foresight solver on augmented system ----
  if (verbose) .dynhr_cat("[ramsey_obc_pf] Running perfect-foresight Newton solver",
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
  ## The discount ramsey_model() resolved (discount= in ..., the .mod's
  ## planner_discount, else beta/betta) -- not beta with a silent 0.99.
  discount <- ramsey_result$welfare$discount
  ## Evaluated AT the steady state (this used to pass
  ## as.numeric(model$var_names) -- the variable NAMES coerced to NA -- which
  ## warned on every call and returned NA).
  ss_orig <- stats::setNames(as.numeric(aug_ss[model$var_names]),
                             model$var_names)
  welfare_steady <- .eval_planner_ast(
    parse_expression(obj_text,
      var_names = c(model$var_names, model$varexo_names),
      param_names = model$param_names),
    ss_orig, params, ss_orig
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


#' Does any model equation carry an MCP tag or a complementarity condition?
#' @noRd
.model_has_obc_tags <- function(model) {
  any(vapply(model$equations, function(eq) {
    if (!is.null(eq$complementarity)) return(TRUE)
    tag <- eq$tag_raw %||% eq$tag
    !is.null(tag) && length(tag) == 1L && !is.na(tag) &&
      grepl("mcp\\s*=", tag)
  }, logical(1)))
}

#' OBC specs on the augmented Ramsey system from a ramsey_constraints block
#'
#' Dynare 7 pairs a `ramsey_constraints` bound on variable x with the Ramsey
#' first-order condition with respect to x
#' (`+<fname>/dynamic_complementarity_conditions.m`: lb/ub on x, the FOC
#' "Ramsey FOC w.r.t. x" as its complementary equation).  In the augmented
#' model of ramsey_augment_mod() the FOCs follow the original equations in
#' variable order, skipping identically-zero ones, so the FOC of x is equation
#' n_eq + (number of nonzero FOCs up to x).  Bounds are LEVELS (the path
#' solver works in levels), evaluated at `params`.
#' @return list of OBC specs (eq_idx, var_idx, var_name, op, bound).
#' @noRd
.ramsey_constraint_specs <- function(model, ramsey_result, params) {
  rc <- model$ramsey_constraints
  focs <- ramsey_result$augmented_result$foc_equations
  aug_model <- ramsey_result$augmented_model
  if (is.null(focs) || length(focs) != length(model$var_names))
    .dynhr_abort("ramsey_constraints need the augmented Ramsey system ",
                 "(ramsey_model(method = \"augmented\")).",
                 class = "dynhr_error_ramsey_constraints")
  nonzero <- !vapply(focs, ast_is_zero, logical(1))
  env <- .dynhr_param_eval_env(params)
  lapply(seq_len(nrow(rc)), function(k) {
    v <- rc$var[k]
    j <- match(v, model$var_names)
    if (!nonzero[j])
      .dynhr_abort("ramsey_constraints: the Ramsey first-order condition ",
                   "with respect to ", v, " is identically zero, so there is ",
                   "no equation to pair with its bound.",
                   class = "dynhr_error_ramsey_constraints")
    bx <- gsub("(?<![A-Za-z0-9_])inf(?![A-Za-z0-9_])", "Inf", rc$bound[k],
               perl = TRUE)
    b <- .dynhr_sandbox_eval(bx, env, context = "the ramsey_constraints bound")
    if (!is.numeric(b) || length(b) != 1L || is.na(b))
      .dynhr_abort("ramsey_constraints: the bound `", rc$bound[k], "` of ", v,
                   " does not evaluate to a number at the parameters.",
                   class = "dynhr_error_ramsey_constraints")
    eq_idx <- length(model$equations) + sum(nonzero[seq_len(j)])
    if (eq_idx > length(aug_model$equations))
      .dynhr_abort("ramsey_constraints: cannot locate the Ramsey first-order ",
                   "condition with respect to ", v, " in the augmented model.",
                   class = "dynhr_error_ramsey_constraints")
    list(eq_idx   = eq_idx,
         var_idx  = match(v, aug_model$var_names),
         var_name = v,
         op       = rc$op[k],
         bound    = as.numeric(b))
  })
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
#'   4. Run the OccBin guess-and-verify iteration of
#'      \code{\link{boehl_solve_regime_path}} on the impulse path: each
#'      candidate regime path is solved with its time-varying rules (the
#'      backward recursion over the regimes ahead) and verified on the
#'      constrained path.
#'   5. Return slack and binding decision rules for the augmented system.
#'
#' On a linear-quadratic Ramsey problem the result is the perfect-foresight
#' path of \code{\link{ramsey_obc_pf}} (Dynare's \code{lmmcp}), including
#' multi-period spells and a bound that binds only because another one does.
#'
#' @param model             A dynhr_mod object.
#' @param planner_objective Character string: the planner objective expression.
#' @param obc_specs         List of OBC specs from \code{\link{obc_parse_tags}}.
#'   If NULL, the specs come from the model: its MCP tags / complementarity
#'   conditions and its \code{ramsey_constraints} block (a bound on \code{x}
#'   is complementary to the Ramsey first-order condition with respect to
#'   \code{x}, as in \code{\link{ramsey_obc_pf}}).
#' @param params            Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param method            Regime search on the impulse path.  Both run the
#'   deterministic complementarity iteration of
#'   \code{\link{boehl_solve_regime_path}} (every period is re-checked at
#'   every iteration, so a period can bind and later be released):
#'   \code{"guess_verify"} (default, OccBin) starts from binding exactly the
#'   periods where the all-slack path violates a bound, \code{"boehl"} starts
#'   from the all-slack path.
#' @param shock_std         Standard deviation for the shock used in regime path
#'   computation (default 1.0, i.e. one-standard-deviation impulse).
#' @param max_iter          Maximum regime search iterations (default 100).
#' @param tol               Relative round-off band of the regime decisions
#'   (default 1e-6): passed to \code{\link{boehl_solve_regime_path}}, whose
#'   bound and multiplier tests are relative to \eqn{|x| + |b|} and to the
#'   tagged equation's term magnitudes, so the regimes do not depend on the
#'   units of the model.
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to \code{\link{ramsey_model}}.
#'
#' @return An object of class \code{dynhr_ramsey_obc_pwl} containing:
#'   \item{ramsey_result}{The underlying \code{dynhr_ramsey_result2} object.}
#'   \item{dr_slack}{Slack-regime DecisionRules for the augmented system.}
#'   \item{dr_bind}{Binding-regime rule of a period in which every spec binds
#'     and the next period is slack (the last period of a spell).}
#'   \item{regime_path}{Integer vector: optimal regime path (bitfield encoding).}
#'   \item{irf}{n_aug x T matrix of deviations from the Ramsey steady state
#'     under the Ramsey OBC policy.}
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
  # With obc_specs = NULL the specs are the model's MCP tags plus its
  # ramsey_constraints block (mapped onto the augmented system in step 6),
  # exactly as in ramsey_obc_pf().  Until 2026-09-25 a
  # ramsey_constraints block was refused here: the per-regime binding policy
  # assumed the NEXT period slack (a spell of 2+ periods was solved with the
  # wrong expectation) and the regime search flagged a bound only when the
  # SLACK policy violated it (a bound violated only because another one
  # binds was never imposed).  The regime search now uses the time-varying
  # OccBin rules and verifies on the constrained path, and matches
  # ramsey_obc_pf() / Dynare lmmcp.
  use_rc <- is.null(obc_specs) && is.data.frame(model$ramsey_constraints) &&
    nrow(model$ramsey_constraints) > 0L
  if (is.null(obc_specs)) {
    obc_specs <- if (!use_rc || .model_has_obc_tags(model))
      obc_parse_tags(model) else list()
  }

  # ---- 3. Build augmented Ramsey system (order 1) ----
  if (verbose) .dynhr_cat("[ramsey_obc_pwlinear] Building augmented Ramsey system...\n")

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
  if (verbose) .dynhr_cat("[ramsey_obc_pwlinear] Extracting augmented system matrices...\n")

  sys <- extract_system_matrices(aug_compiled, aug_ss, params)
  # Augment sys with n_endo, endo_names, etc. needed by obc helpers
  sys$n_endo     <- n_aug
  sys$endo_names <- aug_model$var_names
  sys$exo_names  <- aug_model$varexo_names
  sys$n_exo      <- length(sys$exo_names)

  # ---- 5. Solve slack-regime perturbation (order 1) ----
  if (verbose) .dynhr_cat("[ramsey_obc_pwlinear] Solving slack-regime perturbation...\n")

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
  if (use_rc) {
    obc_specs <- c(obc_specs,
                   .ramsey_constraint_specs(model, ramsey_result, params))
    aug_specs <- obc_specs
  }
  n_spec <- length(aug_specs)

  # The rule of a period in which every spec binds and the NEXT period is
  # slack (obc_solve_binding(), Guerrieri-Iacoviello 2015 terminal
  # substitution): the last period of a spell.  Earlier spell periods follow
  # the time-varying rules of the regime search below.
  obs_idx <- seq_len(n_aug)

  dr_bind <- obc_solve_binding(sys, dr_slack, aug_specs, obs_idx)

  # ---- 7. Run regime search ----
  if (verbose) .dynhr_cat("[ramsey_obc_pwlinear] Running regime search (method = '", method, "')...\n", sep = "")

  # Build shock sequence: one std dev impulse to the first shock in period 1.
  # boehl_simulate()/boehl_solve_regime_path() take an n_exo x T matrix
  # (shocks in ROWS, periods in COLUMNS); the previous T x n_exo layout was
  # non-conformable and made this function fail on every call.
  n_T <- 40L  # Default horizon
  shock_seq <- matrix(0, nrow = length(aug_model$varexo_names), ncol = n_T)
  rownames(shock_seq) <- aug_model$varexo_names
  if (nrow(shock_seq) >= 1) {
    shock_seq[1, 1] <- shock_std
  }

  # The shock path is KNOWN (an impulse), so the regime search is the
  # deterministic complementarity iteration on the simulated path (each
  # iteration re-checks EVERY period, so a period can bind and later be
  # released).  "boehl" starts it from the all-slack path; "guess_verify"
  # starts from OccBin's initial guess, binding exactly where the all-slack
  # path violates a bound.  This used to call obc_guess_verify(), the
  # Kalman-FILTER regime search, on the all-slack simulation: its observation
  # pre-pass locked every period where the simulated DEVIATION was at or
  # below the LEVEL bound as binding for good, so a one-period ZLB episode
  # came back binding for the whole horizon.
  init <- NULL
  if (method == "guess_verify") {
    cache <- new.env(parent = emptyenv(), hash = TRUE)
    obc_ensure_policy(0L, cache, sys, dr_slack, aug_specs, obs_idx)
    slack_paths <- boehl_simulate(shock_seq, dr_slack, cache, integer(n_T))$paths
    bnd_dev <- .obc_bound_dev(aug_specs, dr_slack)
    init <- integer(n_T)
    for (j in seq_along(aug_specs)) {
      x <- slack_paths[aug_specs[[j]]$var_idx, ]
      viol <- if (aug_specs[[j]]$op == ">") x < bnd_dev[[j]] else x > bnd_dev[[j]]
      init[viol] <- bitwOr(init[viol], 2L^(j - 1L))
    }
  }
  regime_search <- boehl_solve_regime_path(
    shock_seq        = shock_seq,
    dr_slack         = dr_slack,
    sys              = sys,
    specs            = aug_specs,
    obs_idx          = obs_idx,
    max_iter         = max_iter,
    tol              = tol,
    regime_path_init = init
  )
  regime_path <- regime_search$regime_path
  irf_paths   <- regime_search$paths

  # ---- 8. Compute welfare ----
  discount <- ramsey_result$welfare$discount
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
#' @param shock_std         Standard deviation of the one-time impulse applied
#'   to the first exogenous shock at period 1 of the perfect-foresight path
#'   used to compute the OBC-constrained welfare (default 1.0). Larger values
#'   trace out a larger occasionally-binding episode and a larger welfare gap.
#'   NOTE the comparison is steady-state unconstrained welfare vs the
#'   shocked-path OBC welfare, so the reported gap bundles the shock's own
#'   transition cost together with the OBC distortion (an unconstrained
#'   path-welfare leg under the same shock is not currently available).
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
  discount <- .ramsey_discount(NULL, model, params,
                               "ramsey_obc_welfare_cost")$value

  # Unconstrained Ramsey
  if (verbose) .dynhr_cat("[ramsey_obc_welfare_cost] Solving unconstrained Ramsey...\n")
  ramsey_uncon <- ramsey_model(model, params = params,
    planner_objective = obj_text, order = 1L, verbose = verbose)

  # OBC-constrained Ramsey (perfect-foresight)
  # Build shock path: one std-dev impulse to the first exogenous shock at
  # period 1, matching the convention used by ramsey_obc_pwlinear() (see
  # `shock_seq[1, 1] <- shock_std` above).
  if (verbose) .dynhr_cat("[ramsey_obc_welfare_cost] Solving OBC-constrained Ramsey (PF)...\n")
  T_horizon <- 40L
  shock_path <- matrix(0, nrow = T_horizon, ncol = length(model$varexo_names))
  colnames(shock_path) <- model$varexo_names
  if (ncol(shock_path) >= 1) {
    shock_path[1, 1] <- shock_std
  }
  ramsey_obc <- ramsey_obc_pf(model,
    planner_objective = obj_text, obc_specs = obc_specs,
    shock_path = shock_path, T_horizon = T_horizon,
    params = params, verbose = verbose)

  welfare_uncon <- .extract_welfare(ramsey_uncon, "steady")
  # Use the shock-path welfare (not steady_value, which is shock-independent
  # by construction) so that the comparison is sensitive to shock_std.
  welfare_obc   <- ramsey_obc$welfare$path_value

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

