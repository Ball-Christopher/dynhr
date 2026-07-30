## R/solve-perturbation.R
## --------------------------------------------------------------------------
## First-order perturbation solver for DSGE models (Klein 2000 / Villemot
## 2011). Top-level solve_perturbation() orchestrator + the cached fast
## variant solve_perturbation_fast(), the .solve_from_system() core that
## classifies variables, applies QR static-elimination, builds the reduced
## pencil, dispatches to .solve_qz(), and recovers the decision rules
## (ghx, ghu) plus the DecisionRules S3 object.
##
## Phase-1 split from perturbation-monolith.R (no logic changes).
## --------------------------------------------------------------------------

## =====================================================================
## dynhr_perturbation.R
##
## First-order perturbation solver for DSGE models using the
## Villemot (2011) / Klein (2000) approach with static variable
## elimination via QR decomposition.
##
## Key steps:
##   1. Classify variables: static, backward-only, mixed, forward-only
##   2. Reorder columns to [static | backward | mixed | forward]
##   3. QR-eliminate static variables from the system
##   4. Build reduced companion-form pencil (D, E) of dimension p x p
##      where p = n_minus + n_plus
##   5. QZ decomposition -> check Blanchard-Kahn conditions
##   6. Extract decision rules ghx, ghu
##   7. Recover static variable responses
##
## References:
##   Klein, P. (2000). Using the generalized Schur form to solve a
##     multivariate linear rational expectations model.
##   Villemot, S. (2011). Solving rational expectations models at first
##     order: what Dynare does.
## =====================================================================

.reconcile_endo_exo <- function(sys, model) {
  # -----------------------------------------------------------------
  # Case 1: varexo names leaked into endo list -> remove them
  # -----------------------------------------------------------------
  ## match()-based membership is cheaper than intersect() (no sort/unique) and
  ## the common case is no overlap, so we only build the name list when needed.
  hit <- sys$exo_names %in% sys$endo_names
  if (any(hit)) {
    exo_in_endo <- sys$exo_names[hit]
    message(sprintf(
      "Removing %d varexo names found in endo list: %s",
      length(exo_in_endo), paste(exo_in_endo, collapse = ", ")))

    keep <- !(sys$endo_names %in% exo_in_endo)
    keep_idx <- which(keep)

    sys$endo_names <- sys$endo_names[keep_idx]
    sys$n_endo     <- length(sys$endo_names)

    # Drop corresponding columns from Jacobians
    sys$f_zero  <- sys$f_zero[,  keep_idx, drop = FALSE]
    sys$f_minus <- sys$f_minus[, keep_idx, drop = FALSE]
    sys$f_plus  <- sys$f_plus[,  keep_idx, drop = FALSE]

    # Update variable classification vectors
    sys$is_static <- sys$is_static[keep_idx]
    sys$is_pred   <- sys$is_pred[keep_idx]
    sys$is_fwd    <- sys$is_fwd[keep_idx]
    sys$is_mixed  <- sys$is_mixed[keep_idx]
  }

  # -----------------------------------------------------------------
  # Case 2: genuine auxiliary lead/lag equations missing
  # (keep the earlier logic as fallback for future models)
  # -----------------------------------------------------------------
  n_eq   <- nrow(sys$f_zero)
  n_endo <- sys$n_endo
  if (n_eq != n_endo) {
    stop(sprintf(
      "After reconciliation, Jacobian rows (%d) != n_endo (%d). Check parser.",
      n_eq, n_endo))
  }

  sys
}

# =====================================================================
# Loglinear (log-deviation) transform for first-order decision rules
# =====================================================================

#' Apply the loglinear (log-deviation) transform to a first-order DR.
#'
#' Converts ghx and ghu from level deviations to log-deviations for all
#' variables with strictly positive steady state. Variables with non-positive
#' SS are kept in levels (matching Dynare's loglinear option).
#'
#' Derivation: if y_t = y_ss + delta_y_t (levels) then
#'   yhat_t = delta_y_t / y_ss  (log-deviation for small deviations)
#' So:
#'   yhat_t = (ghx_lev * delta_x_{t-1} + ghu_lev * e_t) / ys[i]
#'   yhat_t = ghx_lev * ys[j] * xhat_{t-1} / ys[i] + ghu_lev * e_t / ys[i]
#' i.e.  ghx_log[i,j] = ghx_lev[i,j] * ys[j] / ys[i]   (both log-transformed)
#'       ghu_log[i,j] = ghu_lev[i,j] / ys[i]            (row i log-transformed)
#' Column scaling by ys[j] only applies when the state variable j has positive SS.
#'
#' @param dr   DecisionRules object (order=1) from .solve_from_system.
#' @param ss   Named numeric steady state (all endogenous variables).
#' @return dr with ghx and ghu replaced by log-deviation versions; adds
#'         field dr$loglinear_vars naming the log-transformed variables.
#' @noRd
.apply_loglinear_transform <- function(dr, ss) {
  endo  <- dr$endo_names   # n_endo names (rows of ghx/ghu)
  state <- dr$state_vars   # n_state names (cols of ghx)
  exo   <- dr$exo_names    # n_exo names  (cols of ghu)

  ## ys values aligned to endo and state (by name).
  ys_endo  <- ss[endo]          # n_endo
  ys_state <- ss[state]         # n_state

  ## Flag which variables are log-transformed (strictly positive SS).
  ## Variables with SS <= 0 stay in levels (Dynare convention).
  log_endo  <- is.finite(ys_endo)  & (ys_endo  > 0)
  log_state <- is.finite(ys_state) & (ys_state > 0)

  ghx <- dr$ghx  # n_endo x n_state
  ghu <- dr$ghu  # n_endo x n_exo

  ## --- Scale ghx ---
  ## ghx_log[i,j] = ghx_lev[i,j] * ys[j] / ys[i]
  ## Only scale row i by 1/ys[i] when log_endo[i] is TRUE.
  ## Only scale col j by ys[j]  when log_state[j] is TRUE.
  if (ncol(ghx) > 0L) {
    ## Column scaling: multiply each col j by ys_state[j] when log_state[j].
    ## (sweep over cols)
    col_scale <- ifelse(log_state, ys_state, 1.0)
    ghx <- sweep(ghx, 2L, col_scale, `*`)

    ## Row scaling: divide each row i by ys_endo[i] when log_endo[i].
    row_scale <- ifelse(log_endo, ys_endo, 1.0)
    ghx <- sweep(ghx, 1L, row_scale, `/`)
  }

  ## --- Scale ghu ---
  ## ghu_log[i,j] = ghu_lev[i,j] / ys[i]
  ## Only scale row i by 1/ys[i] when log_endo[i] is TRUE.
  row_scale <- ifelse(log_endo, ys_endo, 1.0)
  ghu <- sweep(ghu, 1L, row_scale, `/`)

  dr$ghx            <- ghx
  dr$ghu            <- ghu
  dr$loglinear      <- TRUE
  dr$loglinear_vars <- endo[log_endo]
  dr
}

# =====================================================================
# Main entry points
# =====================================================================

#' Solve the perturbation of a DSGE model (first to fifth order)
#'
#' Uses the Villemot (2011) method for first order: classify variables,
#' eliminate statics via QR, form reduced companion pencil, apply QZ
#' decomposition, check Blanchard-Kahn conditions, extract decision rules
#' ghx and ghu.
#'
#' Higher orders extend via successive approximation:
#'   Order 2: Schmitt-Grohé & Uribe (2004) Kronecker-product method
#'   Order 3: Binning (2013) recursive Sylvester with Faà di Bruno
#'   Orders 4-5: Levintal (2017) compact tensor Sylvester recursion
#'
#' @param model    dynhr_mod object
#' @param compiled dynhr_compiled
#' @param ss       Named numeric steady state
#' @param params   Named numeric parameters
#' @param verbose  Print progress
#' @param order    Perturbation order: 1 (default), 2, 3, 4, or 5
#'
#' \strong{Dynare name mapping (order >= 2):} the field \code{dr$ghss} holds
#' the order-2 uncertainty-correction vector (the sigma^2 correction to the
#' stochastic steady state).  In Dynare this is \code{oo_.dr.ghs2}; dynhr
#' uses \code{ghss} to follow standard mathematical notation (\eqn{g_{\sigma\sigma}}).
#' There is no \code{ghs2} field in dynhr.
#' @param Sigma_e  (order >= 2) n_exo x n_exo shock covariance. NULL uses
#'                 model shocks block.
#' @param h        Finite-difference step for numerical derivatives.
#'                 Default 1e-4 for order 2, 1e-2 for order 4, 5e-2 for order 5.
#' @param sigma3   (order=3 only) Optional n_u^3 vector of third-order shock
#'                 product moments E[u_i u_j u_k]. Default NULL leaves the
#'                 third-cumulant correction at zero (Gaussian shocks).
#' @param loglinear Logical (default FALSE). When TRUE, returns the first-order
#'                 decision rule expressed in log-deviations, matching Dynare's
#'                 \code{loglinear} option. Each variable with a strictly positive
#'                 steady state is log-transformed (y -> exp(yhat) where
#'                 yhat = log(y) - log(y_ss)), so ghx and ghu give percentage
#'                 deviations from steady state. Variables with non-positive
#'                 steady state are kept in levels. The transform is a simple
#'                 Jacobian scaling of the first-order rule:
#'                 ghx_log[i,j] = ghx_lev[i,j] * ys[j] / ys[i],
#'                 ghu_log[i,j] = ghu_lev[i,j] / ys[i],
#'                 where the scaling only applies to rows/columns corresponding
#'                 to log-transformed variables (positive SS). Currently scoped
#'                 to order=1. The FALSE default path is byte-identical to the
#'                 previous behaviour.
#' @param center  Optional named numeric vector giving the point at which the
#'                 model equations are linearised. Defaults to \code{NULL}
#'                 (the steady state \code{ss}); supply a non-steady centre to
#'                 expand around an alternative path.
#' @return For order=1: DecisionRules object.
#'   For order=2: DecisionRules2 object.
#'   For order=3: DecisionRules3 object.
#'   For order=4: DecisionRules4 object (Levintal compact).
#'   For order=5: DecisionRules5 object (Levintal compact).
#'
#' @references
#'   Klein, P. (2000). Using the generalized Schur form to solve a multivariate
#'     linear rational expectations model. \emph{Journal of Economic Dynamics and
#'     Control}, 24(10), 1405-1423.
#'   Villemot, S. (2011). Solving rational expectations models by first-order
#'     perturbation. \emph{Dynare Working Paper Series}, 6.
#'   Schmitt-Grohé, S., & Uribe, M. (2004). Solving dynamic general equilibrium
#'     models using a second-order approximation to the policy function.
#'     \emph{Journal of Economic Dynamics and Control}, 28(4), 755-775.
#'   Binning, A. (2013). Underidentified SVAR models: A framework for combining
#'     short-run and long-run restrictions. \emph{Norges Bank Working Paper}.
#'   Levintal, O. (2017). Fifth-order perturbation solution of DSGE models.
#'     \emph{Journal of Economic Dynamics and Control}, 80, 1-16.
#' @export
solve_perturbation <- function(model, compiled, ss, params, verbose = FALSE,
                               order = 1L, Sigma_e = NULL, h = NULL,
                               sigma3 = NULL, loglinear = FALSE,
                               center = NULL) {
  order <- as.integer(order)
  if (!order %in% c(1L, 2L, 3L, 4L, 5L)) stop("order must be 1, 2, 3, 4, or 5")

  ## center= validation: only supported for order=1.
  if (!is.null(center)) {
    if (order > 1L) {
      stop(paste0(
        "solve_perturbation: `center` (non-steady-state linearization point) ",
        "is only supported for order = 1. For order >= 2 the Hessian/higher ",
        "tensors are evaluated at the steady state; off-SS higher-order rules ",
        "are not yet implemented."))
    }
    if (!is.numeric(center) || is.null(names(center))) {
      stop("solve_perturbation: `center` must be a named numeric vector.")
    }
  }

  # Guard: the compiled model must carry symbolic derivatives deep enough for
  # `order`. compile_model(max_order=) builds derivatives one past its value
  # for the 2/3 case (max_order>=2 yields BOTH Hessian and 3rd derivatives;
  # only 4th/5th order need max_order>=4). Without enough derivatives the
  # higher-order solvers silently consume zero tensors and return all-zero
  # ghxx/ghuu/etc. instead of failing loudly.
  compiled_order <- compiled$dynamic$max_order
  required_co <- if (order <= 1L) 1L else if (order <= 3L) 2L else 4L
  if (!is.null(compiled_order) && compiled_order < required_co) {
    stop(sprintf(
      paste0("order = %d requires the model compiled with max_order >= %d ",
             "(got %d). Re-run compile_model(model, max_order = %d)."),
      order, required_co, compiled_order, required_co))
  }

  # Default step sizes
  if (is.null(h)) {
    h <- switch(as.character(order),
      "2" = 1e-4, "4" = 1e-2, "5" = 5e-2, 1e-4)
  }

  # ------------------------------------------------------------------
  # Auto-SS: if ss is NULL, solve steady state first.
  # Accepts:
  #   NULL           → solve_steady() internally
  #   dynhr_steady   → extract $values
  #   numeric vector → use directly
  # ------------------------------------------------------------------
  if (is.null(ss)) {
    if (verbose) cat("  Solving steady state (ss=NULL)...\n")
    ss_result <- solve_steady(compiled, params,
                               y0 = setNames(rep(0, length(compiled$model$var_names)),
                                             compiled$model$var_names),
                               verbose = verbose)
    if (isTRUE(ss_result$converged)) {
      ss <- ss_result$values
      if (verbose) cat(sprintf("  Steady state solved (max|ss| = %.6e)\n",
                               max(abs(ss))))
    } else {
      stop("Auto steady-state solve failed. ",
           "Provide an explicit `ss` vector or check the model calibration.")
    }
  } else if (inherits(ss, "dynhr_steady")) {
    # If ss_result has updated params (from analytical SS filling in NAs),
    # merge them before extracting system matrices
    if (!is.null(ss[["params"]])) {
      for (nm in names(ss[["params"]])) {
        val <- ss[["params"]][[nm]]
        if (is.finite(val)) params[[nm]] <- val
      }
    }
    ss <- ss$values
  }

  ## Name an UNNAMED numeric ss positionally in declaration order. Every
  ## downstream consumer (extract_system_matrices_fast, the jac-tape path)
  ## looks the steady state up BY NAME: `ss[dyn$endo_names]`. An unnamed ss
  ## therefore silently misses every lookup and is zero-filled, yielding a
  ## Jacobian evaluated at ss = 0 (exp(0)=1 for log-level models) -- a
  ## plausible-but-wrong system that can still satisfy Blanchard-Kahn. This
  ## bit the Born_Pfeifer_2014 review harness, which pinned to Dynare's
  ## (unnamed) ys and saw O(1) ghx/ghu errors. The only sane convention for
  ## an unnamed ss is declaration order, which is exactly what solve_steady()
  ## returns names for; assign them, and fail loudly on a length mismatch.
  if (is.numeric(ss) && is.null(names(ss))) {
    vn <- compiled$model$var_names
    if (length(ss) != length(vn)) {
      stop(sprintf(
        paste0("solve_perturbation: `ss` is unnamed with length %d but the ",
               "model has %d endogenous variables. Supply a named steady ",
               "state, or an unnamed vector in declaration order."),
        length(ss), length(vn)))
    }
    names(ss) <- vn
  }

  ## M25 fail-fast: a non-finite steady state (solve_steady did not converge, or
  ## undeclared/missing parameters -- e.g. params computed in an external
  ## steadystate.m -- left the SS degenerate) makes the dynamic Jacobian
  ## non-finite. The higher-order Kronecker/Sylvester solves then either error
  ## confusingly ("system is exactly singular") or run effectively unbounded
  ## (the Basu_Bundick_2017 order-3 "hang"). Fail loudly and early with an
  ## actionable message instead of grinding into a doomed solve.
  if (is.numeric(ss) && !is.null(names(ss))) {
    ss_endo <- ss[intersect(names(ss), compiled$model$var_names)]
    bad_ss  <- names(ss_endo)[!is.finite(ss_endo)]
    if (length(bad_ss) > 0L) {
      stop(sprintf(
        paste0("solve_perturbation: the steady state is non-finite for %d ",
               "variable(s): %s. The steady state did not converge -- check the ",
               "calibration and for parameters that are referenced in equations ",
               "but not supplied (e.g. computed in an external steadystate.m; ",
               "supply them with inject_params() for an external/superset vector, ",
               "or set_param_values() for the exact model-param set). Cannot solve the ",
               "perturbation at order %d."),
        length(bad_ss), paste(head(bad_ss, 10L), collapse = ", "), order))
    }
  }

  ## I2-extended fail-loud: declared parameters that are absent or NA at solve
  ## time will silently corrupt the dynamic Jacobian.  The most common causes:
  ##   (a) computed only in an external *_steadystate.m (dynhr cannot run it);
  ##   (b) computed in a steady_state_model block but that block did not run
  ##       (e.g. model(linear) with a degenerate SSM path).
  ## This check fires AFTER the ss$params merge above (line ~190) so that
  ## legitimate SSM-computed params that were absent from the initial param_values
  ## but returned by solve_steady() via ss$params are present by this point.
  ## Scope to declared params that actually feed the Jacobian (referenced in an
  ## equation / planner objective / inlined #-define).  Falls back to all
  ## declared params for older model objects without the field.  This excludes
  ## shock-only params (e.g. a `stderr sigma_e` std), which are absent from
  ## `params` by design and never enter the decision rule.
  .decl <- model$equation_param_names
  if (is.null(.decl)) .decl <- model$param_names
  if (!is.null(.decl) && length(.decl) > 0L) {
    .absent  <- setdiff(.decl, names(params))
    .na_vals <- .decl[.decl %in% names(params) & is.na(params[.decl])]
    .bad_params <- union(.absent, .na_vals)
    if (length(.bad_params) > 0L) {
      stop(sprintf(
        paste0("solve_perturbation: %d declared parameter(s) are NA or absent ",
               "at solve time: %s. These are typically computed in an external ",
               "*_steadystate.m (which dynhr cannot run) or in a ",
               "steady_state_model block that did not propagate. ",
               "Use inject_params(model, named_vec) to fill them from an external ",
               "(possibly superset) vector without errors, or set_param_values() ",
               "for the exact model-param set, or pass the correct params argument. ",
               "Cannot build a valid decision rule."),
        length(.bad_params),
        paste(head(.bad_params, 10L), collapse = ", ")),
        call. = FALSE)
    }
  }

  ## OccBin relax-regime dynamic Jacobian row selection.
  ## When the model has OccBin constraints, dyn$n_eq > n_endo: the dynamic
  ## Jacobian returns all n_eq rows (both bind- and relax-equations).  The row
  ## selection and declaration-order mapping are handled inside
  ## extract_system_matrices (slow path), which already applies eq_to_decl row
  ## reordering.  Use the slow path for OccBin models to get both row selection
  ## and correct row ordering; the fast path skips eq_to_decl reordering and
  ## would need mirrored logic.
  if (inherits(compiled, "dynhr_compiled") &&
      !is.null(compiled$occbin) &&
      compiled$dynamic$n_eq > length(model$var_names)) {
    if (verbose)
      message(sprintf(
        "OccBin perturbation: using slow extract path, selecting relax-regime rows from %d equations",
        compiled$dynamic$n_eq))
    sys <- extract_system_matrices(compiled, ss, params, center = center)
  } else {
    ## Performance: cache the model-structure parsing keyed on the compiled
    ## object's invariant pieces. The cache depends only on `compiled` (lli +
    ## dynamic column map), not on ss/params, so it is safe to reuse. The
    ## cache lives in a package-private env to survive across calls (R copies
    ## `compiled` by value, so attaching it as an attribute would not
    ## persist). R-only; no Rcpp dependency.
    sys_cache <- .get_sys_cache(compiled)
    sys <- extract_system_matrices_fast(sys_cache, ss, params, center = center)
  }

  if (verbose){
    cat("\n=== Full variable classification check ===\n")
    cat(sprintf("  %-15s  %12s  %12s  %-10s  %s\n",
                "variable", "max|f_minus|", "max|f_plus|", "class", "flag"))
    for (i in seq_along(sys$endo_names)) {
      fm <- max(abs(sys$f_minus[, i]))
      fp <- max(abs(sys$f_plus[, i]))
      cls <- if (sys$is_mixed[i]) "mixed"
      else if (sys$is_pred[i]) "backward"
      else if (sys$is_fwd[i]) "forward"
      else "static"
      flag <- ""
      # Classification says no lag, but Jacobian has lag
      if (!sys$is_pred[i] && !sys$is_mixed[i] && fm > 1e-10)
        flag <- "*** HAS LAG BUT NOT CLASSIFIED AS PRED/MIXED"
      # Classification says no lead, but Jacobian has lead
      if (!sys$is_fwd[i] && !sys$is_mixed[i] && fp > 1e-10)
        flag <- "*** HAS LEAD BUT NOT CLASSIFIED AS FWD/MIXED"
      # Classification says lag, but Jacobian has no lag
      if ((sys$is_pred[i] || sys$is_mixed[i]) && fm < 1e-10)
        flag <- "*** CLASSIFIED AS PRED/MIXED BUT NO LAG IN JACOBIAN"
      # Classification says lead, but Jacobian has no lead
      if ((sys$is_fwd[i] || sys$is_mixed[i]) && fp < 1e-10)
        flag <- "*** CLASSIFIED AS FWD/MIXED BUT NO LEAD IN JACOBIAN"
      cat(sprintf("  %-15s  %12.2e  %12.2e  %-10s  %s\n",
                  sys$endo_names[i], fm, fp, cls, flag))
    }
  }

  ## M28: populate dr$Sigma_e so downstream compute_irfs/compute_moments always
  ## use the correct shock covariance (incl. off-diagonal correlations from the
  ## shocks block) rather than re-deriving it on each call or falling back to a
  ## diagonal-only approximation.
  ##
  ## When the caller supplies Sigma_e, use it as-is (no change to prior behaviour).
  ## When Sigma_e is NULL (the common default), derive it from the model's shocks
  ## block via .get_shock_cov(), which honours both stderr and corr/cov entries.
  ## The exo ordering on dr$Sigma_e is dr$exo_names (== sys$exo_names), matching
  ## what compute_irfs / compute_moments expect.
  .attach_Se <- function(dr) {
    if (!is.null(Sigma_e)) {
      dr$Sigma_e <- Sigma_e
    } else {
      dr$Sigma_e <- .get_shock_cov(model, dr$exo_names, params)
    }
    dr
  }

  dr1 <- .solve_from_system(sys, model, compiled, ss, params, verbose)

  ## When the user supplied a non-SS linearization center, attach it so
  ## callers can tell where the rule was evaluated.
  if (!is.null(center)) {
    dr1$center <- center
  }

  ## ------------------------------------------------------------------
  ## Loglinear transform (order=1 only, FALSE by default).
  ## When loglinear=TRUE, the decision rule is expressed in log-deviations
  ## from steady state for each variable with a strictly positive SS value.
  ## The levels rule dr1$ghx / dr1$ghu is converted via:
  ##   ghx_log[i,j] = ghx_lev[i,j] * ys[j] / ys[i]   (both i,j log-transformed)
  ##   ghu_log[i,j] = ghu_lev[i,j] / ys[i]            (row i log-transformed)
  ## Columns (state vars) that are NOT log-transformed are NOT scaled by ys[j].
  ## Rows (output vars) that are NOT log-transformed are NOT scaled by 1/ys[i].
  ## This replicates Dynare's loglinear DR. Currently scoped to order=1.
  ## ------------------------------------------------------------------
  if (isTRUE(loglinear)) {
    dr1 <- .apply_loglinear_transform(dr1, ss)
  }

  if (order >= 2L) {
    dr2 <- solve_perturbation_order2(model, compiled, ss, params,
                                     dr1 = dr1, Sigma_e = Sigma_e,
                                     h = h, verbose = verbose)
    if (order == 2L) return(.attach_Se(dr2))

    dr3 <- solve_perturbation_order3(model, compiled, ss, params,
                                     dr2 = dr2, verbose = verbose)
    dr3 <- solve_sigma_cross(dr3, compiled, ss, params, Sigma_e = Sigma_e)
    if (!is.null(sigma3)) {
      dr3 <- solve_third_cumulant(dr3, compiled, ss, params, sigma3 = sigma3)
    }
    if (order == 3L) return(.attach_Se(dr3))

    dr4 <- solve_perturbation_order4(model, compiled, ss, params,
                                     dr3 = dr3, h = h, verbose = verbose)
    if (order >= 4L) {
      dr4 <- solve_sigma_order4(dr4, compiled, ss, params,
                                 Sigma_e = Sigma_e, verbose = verbose)
    }
    if (order == 4L) return(.attach_Se(dr4))

    dr5 <- solve_perturbation_order5(model, compiled, ss, params,
                                     dr4 = dr4, h = h, verbose = verbose)
    if (order >= 5L) {
      dr5 <- solve_sigma_order5(dr5, compiled, ss, params,
                                 Sigma_e = Sigma_e, verbose = verbose)
    }
    return(.attach_Se(dr5))
  }
  .attach_Se(dr1)
}

#' Wrapper: solve_perturbation with optional cached structure
#'
#' If sys_cache is provided, uses extract_system_matrices_fast instead of
#' the original extract_system_matrices. Everything else is unchanged.
#' @noRd
solve_perturbation_fast <- function(model, compiled, ss, params,
                                   verbose = FALSE, sys_cache = NULL) {
  if (!is.null(sys_cache)) {
    sys <- extract_system_matrices_fast(sys_cache, ss, params)
  } else {
    sys <- extract_system_matrices(compiled, ss, params)
  }
  .solve_from_system(sys, model, compiled, ss, params, verbose)
}

# =====================================================================
# Core solver: Villemot (2011) approach
# =====================================================================

#' Core first-order perturbation solver
#'
#' Implements the full Villemot (2011) algorithm:
#'   - Variable classification and reordering
#'   - Static elimination via QR
#'   - Reduced pencil construction
#'   - QZ decomposition and Blanchard-Kahn check
#'   - Decision rule extraction (ghx, ghu)
#'   - Static variable recovery
#'
#' @param sys     System matrices from extract_system_matrices
#' @param model   dynhr_mod object
#' @param compiled dynhr_compiled
#' @param ss      Named numeric steady state
#' @param params  Named numeric parameters
#' @param verbose Print progress
#' @return Decision rules object (list of class "DecisionRules")
#' @noRd
## TRUE when the compiled C++ static-elimination QR transform is available and
## the Rcpp backend has not been disabled (options(dynhr.use_rcpp = FALSE)).
.HAS_RCPP_STATIC_ELIM <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("qr_static_transform_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

.solve_from_system <- function(sys, model, compiled, ss, params, verbose,
                               pencil_only = FALSE) {
  sys <- .reconcile_endo_exo(sys, model)
  n     <- sys$n_endo
  n_exo <- sys$n_exo

  n_eq <- nrow(sys$f_zero)
  if (n_eq != n) {
    stop(sprintf(
      "Jacobian row count (%d) != n_endo (%d). %d auxiliary equations may be missing.",
      n_eq, n, n - n_eq))
  }

  # ------------------------------------------------------------------
  # Validate system matrices for non-finite values (defense-in-depth)
  # ------------------------------------------------------------------
  # Non-finite entries in f_zero, f_minus, or f_plus will crash QR.
  # The primary fix is the numerical Jacobian fallback in

  # extract_system_matrices; this is the safety net.
  has_bad <- FALSE
  bad_vars_all <- character(0)
  bad_eqs_all  <- integer(0)
  for (mat_name in c("f_zero", "f_minus", "f_plus", "f_exo")) {
    mat <- sys[[mat_name]]
    ## Fast path: clean matrices are the overwhelming common case across MCMC
    ## draws, so skip building the arr.ind index matrix (a per-draw hot spot)
    ## unless a non-finite entry is actually present.
    if (all(is.finite(mat))) next
    bad_cells <- which(!is.finite(mat), arr.ind = TRUE)
    if (nrow(bad_cells) > 0) {
      has_bad <- TRUE
      bad_eqs_all <- c(bad_eqs_all, bad_cells[, 1])
      if (mat_name != "f_exo") {
        bad_vars_all <- c(bad_vars_all, sys$endo_names[bad_cells[, 2]])
      } else {
        bad_vars_all <- c(bad_vars_all, sys$exo_names[bad_cells[, 2]])
      }
      # Replace non-finite with 0 so we can continue
      sys[[mat_name]][!is.finite(sys[[mat_name]])] <- 0
    }
  }
  if (has_bad) {
    bad_vars_all <- unique(bad_vars_all)
    bad_eqs_all  <- sort(unique(bad_eqs_all))
    warning(sprintf(
      paste0("Non-finite values (NA/NaN/Inf) in system Jacobian matrices ",
             "after extract_system_matrices. Affected variables: %s. ",
             "Affected equations: %s. ",
             "Replaced with 0; results may be inaccurate. ",
             "Check the steady state and model equations."),
      paste(bad_vars_all, collapse = ", "),
      paste(head(bad_eqs_all, 10), collapse = ", ")))
  }

  # ------------------------------------------------------------------
  # Step 1: Classify variables
  # ------------------------------------------------------------------
  # static:        no lag, no lead
  # backward-only: lag, no lead
  # mixed:         lag AND lead
  # forward-only:  no lag, lead

  is_static <- sys$is_static
  is_bkw    <- sys$is_pred & !sys$is_mixed   # backward-only
  is_mix    <- sys$is_mixed                   # mixed
  is_fonly  <- sys$is_fwd & !sys$is_mixed     # forward-only

  n_s     <- sum(is_static)
  n_bkw   <- sum(is_bkw)
  n_mix   <- sum(is_mix)
  n_fonly <- sum(is_fonly)
  n_d     <- n_bkw + n_mix + n_fonly          # dynamic variables
  n_minus <- n_bkw + n_mix                    # state (backward) variables
  n_plus  <- n_mix + n_fonly                  # forward variables (pencil dim)
  p       <- n_minus + n_plus                 # companion pencil dimension

  if (verbose) {
    cat("Variable classification:\n")
    cat("  static:", n_s, " backward-only:", n_bkw,
        " mixed:", n_mix, " forward-only:", n_fonly, "\n")
    cat("  n_minus (state):", n_minus, " n_plus (forward):", n_plus,
        " p (pencil dim):", p, "\n")
  }

  # ------------------------------------------------------------------
  # Step 2: Reorder variables to [static | backward | mixed | forward]
  # ------------------------------------------------------------------
  endo <- sys$endo_names
  idx_static <- which(is_static)
  idx_bkw    <- which(is_bkw)
  idx_mix    <- which(is_mix)
  idx_fonly  <- which(is_fonly)
  perm       <- c(idx_static, idx_bkw, idx_mix, idx_fonly)
  inv_perm   <- order(perm)   # inverse permutation to restore original order

  endo_local <- sys$endo_names
  n_plus_bk  <- n_mix + n_fonly

  f_minus_r <- sys$f_minus[, perm, drop = FALSE]
  f_zero_r  <- sys$f_zero[, perm, drop = FALSE]
  f_plus_r  <- sys$f_plus[, perm, drop = FALSE]
  f_exo_r   <- sys$f_exo  # exogenous columns unchanged

  # Column ranges in reordered system
  cols_static <- if (n_s > 0) seq_len(n_s) else integer(0)
  cols_dyn    <- if (n_d > 0) (n_s + 1):n else integer(0)
  cols_minus  <- if (n_minus > 0) (n_s + 1):(n_s + n_minus) else integer(0)
  cols_plus   <- if (n_plus > 0) (n_s + n_bkw + 1):n else integer(0)

  # ------------------------------------------------------------------
  # Step 3: QR elimination of static variables
  # ------------------------------------------------------------------
  if (n_s > 0 && n_d > 0) {
    f_static <- f_zero_r[, cols_static, drop = FALSE]
    if (.HAS_RCPP_STATIC_ELIM()) {
      ## C++ uses non-pivoted QR (matches R's dqrdc2 in the full-rank case),
      ## so the pivot is the identity and inv_piv is trivial. Householder
      ## sign differences cancel downstream; see test-static-elim-parity.R.
      qt <- qr_static_transform_cpp(f_static, f_minus_r, f_zero_r,
                                    f_plus_r, f_exo_r)
      Qf_minus <- qt$Qf_minus; Qf_zero <- qt$Qf_zero
      Qf_plus  <- qt$Qf_plus;  Qf_exo  <- qt$Qf_exo
      inv_piv  <- seq_len(n_s)
    } else {
      qr_obj <- qr(f_static)
      Q_full <- qr.Q(qr_obj, complete = TRUE)
      piv <- qr_obj$pivot
      inv_piv <- order(piv)

      # Transform ALL matrices by Q'
      Qf_minus <- t(Q_full) %*% f_minus_r
      Qf_zero  <- t(Q_full) %*% f_zero_r
      Qf_plus  <- t(Q_full) %*% f_plus_r
      Qf_exo   <- t(Q_full) %*% f_exo_r
    }

    if (verbose) {
      cat("--- QZ diagnostic ---\n")
      cat("  n =", n, " n_s =", n_s, " n_minus =", n_minus, " n_plus =", n_plus, "\n")
      cat("  dim(Qf_minus) =", paste(dim(Qf_minus), collapse=" x "), "\n")
      cat("  row index: (n_s+1):n =", (n_s+1), ":", n, "\n")
      cat("  cols_minus =", paste(cols_minus, collapse=","), "\n")
      cat("  # stable eigenvalues =", n_s, " (BK requires", n_minus, ")\n")
      cat("------------------------\n")
    }

    # Static recovery block: top n_s rows (undo column pivoting)
    R_block <- Qf_zero[1:n_s, cols_static, drop = FALSE]
    R_block <- R_block[, inv_piv, drop = FALSE]

    Ap_s   <- Qf_plus[1:n_s, cols_plus, drop = FALSE]
    A0s_d  <- Qf_zero[1:n_s, cols_dyn, drop = FALSE]
    Am_s   <- Qf_minus[1:n_s, cols_minus, drop = FALSE]
    fe_s   <- Qf_exo[1:n_s, , drop = FALSE]

    # Dynamic block: bottom n_d rows
    Am   <- Qf_minus[(n_s + 1):n, cols_minus, drop = FALSE]
    A0d  <- Qf_zero[(n_s + 1):n, cols_dyn, drop = FALSE]
    Ap   <- Qf_plus[(n_s + 1):n, cols_plus, drop = FALSE]
    fe_d <- Qf_exo[(n_s + 1):n, , drop = FALSE]

  } else if (n_s == 0) {
    Am   <- f_minus_r[, cols_minus, drop = FALSE]
    A0d  <- f_zero_r[, cols_dyn, drop = FALSE]
    Ap   <- f_plus_r[, cols_plus, drop = FALSE]
    fe_d <- f_exo_r

    R_block <- NULL
    Ap_s <- NULL; A0s_d <- NULL; Am_s <- NULL; fe_s <- NULL

  } else {
    # n_d == 0: all variables are static (degenerate case)
    ghx <- matrix(0, nrow = n, ncol = 0)
    M <- sys$f_zero
    ghu <- tryCatch(
      solve(M, -sys$f_exo),
      error = function(e) .safe_inv(M) %*% (-sys$f_exo)
    )
    rownames(ghx) <- endo
    rownames(ghu) <- endo
    if (n_exo > 0) colnames(ghu) <- sys$exo_names
    dr <- list(
      ghx = ghx, ghu = ghu, ys = ss,
      endo_names = endo, exo_names = sys$exo_names,
      state_vars = character(0), state_idx = integer(0),
      n_state = 0L, n_exo = n_exo,
      eigenvalues = complex(0), n_stable = 0L, n_unstable = 0L,
      bk_satisfied = TRUE
    )
    class(dr) <- "DecisionRules"
    return(dr)
  }

  # ------------------------------------------------------------------
  # Step 4: Build reduced companion-form pencil (D, E)
  # ------------------------------------------------------------------
  # State vector: z_t = [y^minus_t ; y^plus_{t+1}]
  #   y^minus = [y^bkw ; y^mix]       (n_minus components)
  #   y^plus  = [y^mix ; y^fonly]      (n_plus components)
  #
  # Pencil: D * z_t = E * z_{t-1} + shock terms
  #
  # D (p x p):
  #   Row 1..n_d:     [A0d_left | Ap]
  #   Row n_d+1..p:   I for mixed in minus block: D[n_d+i, n_bkw+i] = 1
  #
  # E (p x p):
  #   Row 1..n_d:     [-Am | 0_{n_d x n_mix} | -A0d_right]
  #   Row n_d+1..p:   I for mixed in plus block: E[n_d+i, n_minus+i] = 1

  if (p == 0) {
    stop("Internal error: p=0 but n_d>0")
  }

  # Split A0d into left (bkw+mix) and right (fonly)
  if (n_minus > 0 && n_fonly > 0) {
    A0d_left  <- A0d[, 1:n_minus, drop = FALSE]
    A0d_right <- A0d[, (n_minus + 1):n_d, drop = FALSE]
  } else if (n_minus > 0) {
    A0d_left  <- A0d[, 1:n_minus, drop = FALSE]
    A0d_right <- matrix(0, nrow = n_d, ncol = 0)
  } else {
    A0d_left  <- matrix(0, nrow = n_d, ncol = 0)
    A0d_right <- A0d
  }

  # Build D matrix (p x p)
  D_mat <- matrix(0, nrow = p, ncol = p)
  if (n_minus > 0) D_mat[1:n_d, 1:n_minus] <- A0d_left
  if (n_plus > 0)  D_mat[1:n_d, (n_minus + 1):p] <- Ap
  if (n_mix > 0) {
    for (i in seq_len(n_mix)) {
      D_mat[n_d + i, n_bkw + i] <- 1
    }
  }

  # Build E matrix (p x p)
  E_mat <- matrix(0, nrow = p, ncol = p)
  if (n_minus > 0) E_mat[1:n_d, 1:n_minus] <- -Am
  if (n_fonly > 0) {
    E_mat[1:n_d, (n_minus + n_mix + 1):p] <- -A0d_right
  }
  if (n_mix > 0) {
    for (i in seq_len(n_mix)) {
      E_mat[n_d + i, n_minus + i] <- 1
    }
  }

  if (verbose) {
    cat("Pencil dimensions:", p, "x", p, "\n")
    cat("Solving generalized eigenvalue problem...\n")
  }

  ## Early exit for callers that only need the reduced companion-form pencil
  ## (e.g. bk_distance()'s finite-difference pencil derivatives, which must
  ## not pay for -- or depend on the eigenvalue ordering of -- a QZ solve).
  ## The generalized eigenvalues lambda solve E_mat x = lambda * D_mat x
  ## (companion dynamics z_t = D_mat^{-1} E_mat z_{t-1}); these are the same
  ## lambda reported in dr$eigenvalues, and |lambda| = 1 is the BK wall.
  if (isTRUE(pencil_only)) {
    return(list(D = D_mat, E = E_mat, n_minus = n_minus, p = p))
  }

  # ------------------------------------------------------------------
  # Step 5: QZ decomposition and Blanchard-Kahn check
  # ------------------------------------------------------------------
  qz_result <- .solve_qz(D_mat, E_mat, n_minus, verbose)

  if (is.null(qz_result)) {
    stop("QZ decomposition failed. Install 'QZ' or 'geigen' package for better support.")
  }

  n_unstable_base <- qz_result$n_unstable_finite
  n_on_unit <- qz_result$n_on_unit_circle %||% 0L

  # Compute the adjusted forward-looking variable count first (needed
  # for unit-root classification heuristic below)
  n_fwd_model <- model$n_forward %||% n_fonly
  n_mix_model <- model$n_mixed %||% n_mix
  n_pred_mix_model <- sum(model$variable_classification$mixed %in%
                          (model$predetermined_vars %||% character(0)))
  n_plus_bk_model <- n_fwd_model + n_mix_model - n_pred_mix_model

  n_aux_fonly <- sum(grepl("^AUX_", endo_local[idx_fonly]))
  n_plus_bk_alt <- n_plus_bk - n_aux_fonly
  # Use the model-based forward count as the primary target, but fall back
  # to n_plus_bk_alt if it matches n_unstable better.  n_plus_bk_alt excludes
  # AUX_LEAD forward-only variables (created for lead > 1 on predetermined
  # vars like k(+2)), which don't represent genuine jump decisions.
  n_plus_bk_adj <- n_plus_bk_model

  # Decide whether to count unit-root eigenvalues as unstable.
  # Strategy: prefer the count that matches n_plus_bk_adj.
  n_unstable_no_unit <- n_unstable_base + qz_result$n_unstable_infinite
  n_unstable_with_unit <- n_unstable_base + n_on_unit + qz_result$n_unstable_infinite

  # Ramsey models (MULT_* or nn1): always count unit roots as unstable
  is_ramsey <- any(grepl("^MULT_", endo_local)) ||
    isTRUE(model$ramsey_context)

  if (is_ramsey && n_on_unit > 0) {
    n_unstable <- n_unstable_with_unit
  } else if (n_on_unit > 0) {
    # Try multiple unstable counts; pick the one matching n_plus_bk_adj
    if (n_unstable_with_unit == n_plus_bk_adj) {
      n_unstable <- n_unstable_with_unit
    } else if (n_unstable_no_unit == n_plus_bk_adj) {
      n_unstable <- n_unstable_no_unit
    } else if (n_unstable_base == n_plus_bk_adj) {
      # Finite-only eigenvalues may match when QZ Inf eigenvalues are
      # artifacts of static equations (matching Dynare's GEP reduction).
      n_unstable <- n_unstable_base
    } else {
      # Neither matches; prefer the one closest to n_plus_bk_adj
      candidates <- c(n_unstable_with_unit, n_unstable_no_unit, n_unstable_base)
      n_unstable <- candidates[which.min(abs(candidates - n_plus_bk_adj))]
    }
  } else {
    # No unit roots: also consider finite-only count as a candidate
    if (n_unstable_no_unit == n_plus_bk_adj) {
      n_unstable <- n_unstable_no_unit
    } else if (n_unstable_base == n_plus_bk_adj) {
      n_unstable <- n_unstable_base
    } else {
      n_unstable <- n_unstable_no_unit
    }
  }

  if (verbose) {
    cat("  Eigenvalues:", n_unstable_base, "unstable_finite,",
        n_on_unit, "on_unit_circle,", n_unstable, "total_unstable\n")
    cat("  Required stable:", n_minus, "\n")
  }

  if (verbose) {
    # -- Eigenvalue diagnostic --
    ev <- qz_result$eigenvalues
    ev_mod <- Mod(ev)
    ev_ord <- order(ev_mod)
    tol_tag <- 1e-6
    cat("  Eigenvalue moduli (sorted):\n")
    for (j in ev_ord) {
      tag <- if (ev_mod[j] < 1 - tol_tag) "S"
             else if (abs(ev_mod[j] - 1.0) < tol_tag) "~1"
             else "U"
      cat(sprintf("    [%2d] |??| = %.10f  %s\n", j, ev_mod[j], tag))
    }
    cat("  n_unstable reported by .solve_qz:", n_unstable, "\n")
    # -----------------------------
  }

  if (verbose) {
    cat("  n_plus (pencil):", n_plus, " n_plus_bk (Jacobian):", n_plus_bk,
        " n_plus_bk_adj (model-based):", n_plus_bk_adj,
        " n_unstable:", n_unstable, "\n")
    cat("    model: n_fwd=", n_fwd_model, " n_mix=", n_mix_model,
        " n_pred_mix=", n_pred_mix_model, "\n")
  }

  bk_ok <- (n_unstable == n_plus_bk_adj)

  # If the model-based count doesn't match, try the alternative count
  # that excludes AUX_LEAD forward-only variables.  AUX_LEAD vars are
  # created for lead > 1 (e.g. k(+2) => AUX_LEAD_k_1), and their
  # forward-only classification is a computational artifact, not a
  # genuine jump decision.
  if (!bk_ok && n_unstable == n_plus_bk_alt) {
    bk_ok <- TRUE
    n_plus_bk_adj <- n_plus_bk_alt
    if (verbose) {
      cat("  BK satisfied using n_plus_bk_alt (excl. AUX forward-only):",
          n_plus_bk_alt, "\n")
    }
  }

  if (!bk_ok) {
    # Use n_plus_bk_adj for the error message (either the model-based
    # count that was tried first, or the alt that was tried second).
    msg <- if (n_unstable > n_plus_bk_adj) {
      sprintf(
        "Blanchard-Kahn violation: %d unstable roots but %d forward variables (NO SOLUTION).",
        n_unstable, n_plus_bk_adj
      )
    } else {
      sprintf(
        "Blanchard-Kahn violation: %d unstable roots but %d forward variables (INDETERMINACY).",
        n_unstable, n_plus_bk_adj
      )
    }
    warning(msg, " Returning DR with bk_satisfied=FALSE.")
    bk_ok <- FALSE
    # Continue to compute decision rules so callers can inspect the result
    # and decide how to handle (e.g. OSR loss function returns a penalty).
  }

  # ------------------------------------------------------------------
  # Step 6: Extract decision rules from QZ decomposition
  # ------------------------------------------------------------------
  Z <- qz_result$Z
  T_mat <- qz_result$T
  S_mat <- qz_result$S

  if (n_minus > 0) {
    Z11 <- Z[1:n_minus, 1:n_minus, drop = FALSE]
    T11 <- T_mat[1:n_minus, 1:n_minus, drop = FALSE]
    S11 <- S_mat[1:n_minus, 1:n_minus, drop = FALSE]

    Z11_inv <- .safe_inv(Z11)
    T11_inv <- .safe_inv(T11)

    g_minus_y <- Re(Z11 %*% T11_inv %*% S11 %*% Z11_inv)

    if (n_plus > 0) {
      Z21 <- Z[(n_minus + 1):p, 1:n_minus, drop = FALSE]
      g_plus_y <- Re(Z21 %*% Z11_inv)
    } else {
      g_plus_y <- matrix(0, nrow = 0, ncol = n_minus)
    }
  } else {
    g_minus_y <- matrix(0, nrow = 0, ncol = 0)
    g_plus_y  <- matrix(0, nrow = 0, ncol = 0)
  }

  # ------------------------------------------------------------------
  # Step 6a (post-processing): Correct QZ deflating-subspace mixing
  # for repeated eigenvalues.
  #
  # When the pencil has repeated eigenvalues (e.g. 0.75 from both the
  # eps_pref AR(1) and structural Calvo dynamics), the QZ deflating
  # subspace is non-unique.  Different LAPACK implementations may
  # select different bases for the same 2D eigenspace, causing
  # cross-contamination of ghx columns.  Here we restore the known
  # pure-AR(1) structure for exogenous state variables whose columns
  # have been contaminated.
  #
  # Only the eps_pref column is known to be affected in this model
  # (its rho=0.75 coincides with the Calvo structural eigenvalue).
  # Additional entries can be added as they are identified by
  # golden-comparison tests.  If the model does not have the
  # specified variable or parameter, the entry is silently skipped.
  # ------------------------------------------------------------------
  # NOTE (2026-05-31): a hardcoded `ar1_contam` patch used to live here. It
  # zeroed the entire ghx column of `eps_pref` (keeping only the AR(1) diagonal
  # `rho_pref`) whenever it detected off-diagonal entries, on the theory that QZ
  # deflating-subspace basis ambiguity for repeated eigenvalues "contaminated"
  # the column. That is mathematically unfounded: the stable deflating subspace
  # is unique, so g_minus_y = Z11 T11^{-1} S11 Z11^{-1} and g_plus_y =
  # Z21 Z11^{-1} are basis-invariant. For models where the shock process
  # genuinely feeds other equations (e.g. a habit Euler `... + eps_pref`), the
  # off-diagonals are REAL economic transmission; zeroing them produced an
  # INVALID decision rule (a ghx that violates the model equations). On
  # 03j_no_omega_cp_est this dropped the entire preference-shock transmission,
  # giving a stochastically-singular Kalman innovation covariance and -Inf
  # log-likelihood. The patch is removed; the QZ solution is already correct.

  # ------------------------------------------------------------------
  # Step 6b: Build ghx_dyn (n_d x n_minus) for dynamic variables
  # ------------------------------------------------------------------
  if (n_minus > 0) {
    ghx_dyn <- matrix(0, nrow = n_d, ncol = n_minus)

    row_ptr <- 0
    if (n_bkw > 0) {
      ghx_dyn[1:n_bkw, ] <- g_minus_y[1:n_bkw, , drop = FALSE]
      row_ptr <- n_bkw
    }
    if (n_mix > 0) {
      ghx_dyn[(row_ptr + 1):(row_ptr + n_mix), ] <-
        g_minus_y[(n_bkw + 1):n_minus, , drop = FALSE]
      row_ptr <- row_ptr + n_mix
    }
    if (n_fonly > 0) {
      ghx_dyn[(row_ptr + 1):n_d, ] <-
        g_plus_y[(n_mix + 1):n_plus, , drop = FALSE]
    }
  } else {
    ghx_dyn <- matrix(0, nrow = n_d, ncol = 0)
  }

  # ------------------------------------------------------------------
  # Step 7: Recover static variable responses
  # ------------------------------------------------------------------
  if (n_s > 0 && n_minus > 0) {
    R_inv <- .safe_inv(R_block)
    if (n_plus > 0) {
      gpy_gmy <- g_plus_y %*% g_minus_y
    } else {
      gpy_gmy <- matrix(0, nrow = 0, ncol = n_minus)
    }
    rhs <- matrix(0, nrow = n_s, ncol = n_minus)
    if (n_plus > 0) rhs <- rhs + Ap_s %*% gpy_gmy
    rhs <- rhs + A0s_d %*% ghx_dyn
    if (n_minus > 0) rhs <- rhs + Am_s
    ghx_s <- -R_inv %*% rhs
  } else if (n_s > 0) {
    ghx_s <- matrix(0, nrow = n_s, ncol = 0)
  } else {
    ghx_s <- matrix(0, nrow = 0, ncol = n_minus)
  }

  # ------------------------------------------------------------------
  # Step 8: Assemble full ghx and reorder to original variable order
  # ------------------------------------------------------------------
  ghx_reordered <- rbind(ghx_s, ghx_dyn)  # n x n_minus

  # Permute rows back to original variable order
  ghx_full <- ghx_reordered[inv_perm, , drop = FALSE]

  # Permute columns: state vars in reordered order are [idx_bkw, idx_mix]
  state_vars_reordered <- c(idx_bkw, idx_mix)
  state_idx <- which(sys$is_pred | sys$is_mixed)
  col_perm <- match(state_idx, state_vars_reordered)
  if (n_minus > 0 && !any(is.na(col_perm))) {
    ghx <- ghx_full[, col_perm, drop = FALSE]
  } else {
    ghx <- ghx_full
  }

  # ------------------------------------------------------------------
  # Step 9: Compute ghu (shock impact matrix)
  # ------------------------------------------------------------------
  P <- matrix(0, nrow = n, ncol = n)
  if (n_minus > 0) {
    P[, state_idx] <- ghx
  }

  M <- sys$f_zero + sys$f_plus %*% P
  ghu <- tryCatch(
    solve(M, -sys$f_exo),
    error = function(e) .safe_inv(M) %*% (-sys$f_exo)
  )

  # NOTE (2026-05-31): a hardcoded `shock_ar1_map` patch used to live here. It
  # zeroed the ghu column of `eps_pref_` (keeping only the unit impact on
  # `eps_pref`), on the false premise that a shock entering a pure-AR(1)
  # equation "should only" hit that AR(1) variable. But when the AR(1) variable
  # feeds other equations contemporaneously (e.g. the preference shock enters
  # the eps_pref state which enters the Euler), the shock DOES have a
  # contemporaneous impact on those variables, which `ghu = solve(M, -f_exo)`
  # computes correctly. Zeroing it dropped the preference shock's
  # contemporaneous transmission and corrupted the likelihood. Removed (see the
  # companion ghx `ar1_contam` note above).

  # ------------------------------------------------------------------
  # Step 9b: Zero ghx columns for trivially-constant state variables
  #
  # A state variable v is trivially constant if its row in both ghx
  # and ghu is numerically zero: v does not move in response to any
  # lagged state or shock.  In that case, v is identically 0 for all
  # t >= 0 (given a zero steady state), so the ghx column for v
  # (how v's lagged value propagates to other variables) is spurious
  # and should be zeroed.
  #
  # Canonical case: ADPM2007 PIE_BAR = 0.  PIE_BAR(-1) appears in
  # the Robs equation, giving PIE_BAR a structural lag and a nonzero
  # f_minus column.  The perturbation solver computes a spurious
  # ghx[:, PIE_BAR] ≈ 2.16; Dynare constant-folds PIE_BAR away and
  # gets 0.  Zeroing the column here matches Dynare's output without
  # changing H or IRFs (PIE_BAR is always 0, so its lagged value
  # never contributes to any forecast).
  #
  # Safety: we do NOT prune the state from the vector (state_idx is
  # unchanged) to avoid re-dimensioning the companion matrix.  We only
  # zero the column, which is always the correct value for a constant.
  #
  # Iid shocks (e.g. ireland_2004 zt with rho_z=0) have ghx[v,] = 0
  # but ghu[v,] != 0, so they are NOT zeroed -- matching Dynare which
  # keeps iid shock states (harmless zero-root state).
  # ------------------------------------------------------------------
  if (n_minus > 0) {
    const_tol <- 1e-10
    for (ci in seq_len(n_minus)) {
      vi <- state_idx[ci]
      if (max(abs(ghx[vi, ])) < const_tol && max(abs(ghu[vi, ])) < const_tol) {
        ghx[, ci] <- 0
      }
    }
  }

  # ------------------------------------------------------------------
  # Step 10: Set names and return DecisionRules object
  # ------------------------------------------------------------------
  if (n_minus > 0) {
    rownames(ghx) <- endo
    colnames(ghx) <- endo[state_idx]
  } else {
    ghx <- matrix(0, nrow = n, ncol = 0)
    rownames(ghx) <- endo
  }
  rownames(ghu) <- endo
  if (n_exo > 0) colnames(ghu) <- sys$exo_names

  n_state <- length(state_idx)

  ## Post-solve stability guard.
  ## The BK count (n_unstable == n_plus_bk_adj) can be satisfied while the QZ
  ## ordering at a near-boundary draw yields a numerically inconsistent decision
  ## rule whose STATE TRANSITION ghx[state_idx, ] is actually EXPLOSIVE (spectral
  ## radius >> 1). Such a "solution" is invalid -- the model has no bounded
  ## equilibrium there -- yet it would be marked bk_satisfied = TRUE and fed to
  ## the Kalman filter, where it yields a garbage (non-stationary) likelihood
  ## after a slow diffuse/univariate fallback (the dominant NZSIM mode-finding
  ## cost: ~0.5 s/eval on such draws). Verify the realized state transition is
  ## non-explosive and flip to bk_satisfied = FALSE otherwise, so callers (e.g.
  ## make_log_posterior) reject it as -Inf. The 1 + 1e-6 margin admits genuine
  ## unit-root / trend models (state radius ~ 1).
  if (bk_ok && n_state > 0L) {
    state_sr <- tryCatch(
      max(Mod(eigen(ghx[state_idx, , drop = FALSE],
                    symmetric = FALSE, only.values = TRUE)$values)),
      error = function(e) NA_real_)
    if (is.finite(state_sr) && state_sr > 1 + 1e-6) {
      if (verbose)
        cat(sprintf(paste0("  Post-solve guard: realized state transition is ",
                           "explosive (radius %.4g > 1); bk_satisfied = FALSE.\n"),
                    state_sr))
      bk_ok <- FALSE
    }
  }

  dr <- list(
    ghx          = ghx,
    ghu          = ghu,
    ys           = ss,
    endo_names   = endo,
    exo_names    = sys$exo_names,
    state_vars   = sys$state_vars,
    state_idx    = state_idx,
    n_state      = n_state,
    n_stable     = n_state,
    n_exo        = n_exo,
    eigenvalues  = qz_result$eigenvalues,
    n_unstable   = n_unstable,
    bk_satisfied = bk_ok
  )
  class(dr) <- "DecisionRules"

  if (verbose) {
    cat("Decision rules computed.\n")
    cat("  ghx:", nrow(dr$ghx), "x", ncol(dr$ghx), "\n")
    cat("  ghu:", nrow(dr$ghu), "x", ncol(dr$ghu), "\n")
  }

  dr
}
