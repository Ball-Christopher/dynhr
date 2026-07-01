## R/solve-model.R
## --------------------------------------------------------------------------
## solve_model() -- Unified entry point for solving a DSGE model.
##
## Orchestrates the full solution pipeline:
##   parse_mod -> compile_model -> solve_steady -> solve_perturbation
##
## Also supports:
##   - Higher-order perturbation (order 1-5)
##   - Ramsey optimal policy workflow
##   - Stoch_simul at solution (IRFs + moments)
## --------------------------------------------------------------------------


#' Solve a DSGE model from a .mod file
#'
#' Parses, compiles, finds steady state, and solves the perturbation of a
#' DSGE model in a single call.  This is the recommended entry point for
#' all downstream analysis (mode-finding, estimation, diagnostics).
#'
#' @param mod_file  Path to a Dynare \code{.mod} file, or a pre-parsed
#'   \code{dynhr_mod} object.
#' @param params   Named numeric parameter vector.  Defaults to the
#'   calibrated values from the \code{.mod} file.
#' @param order    Perturbation order: 1 (default), 2, 3, 4, or 5.
#' @param Sigma_e  (order >= 2) \eqn{n_{exo} \times n_{exo}} shock
#'   covariance matrix.  \code{NULL} uses the model's \code{shocks} block.
#' @param steady_options  List of options passed to \code{solve_steady()},
#'   e.g. \code{list(y0 = ..., tol = 1e-10)}.
#' @param perturbation_options  List of options passed to
#'   \code{solve_perturbation()}, e.g. \code{list(h = 1e-4, sigma3 = NULL)}.
#' @param stoch_simul  If \code{TRUE} (default), compute IRFs and
#'   theoretical moments at the solution.  Set to \code{FALSE} to skip.
#' @param irf_horizons  Number of IRF periods when \code{stoch_simul = TRUE}
#'   (default 40).
#' @param ramsey  If \code{TRUE}, also run the Ramsey optimal policy
#'   workflow at the solved parameters.  Requires a
#'   \code{planner_objective(...)} in the \code{.mod} file.
#' @param ramsey_order  Perturbation order for the Ramsey solution
#'   (1 or 2, default 1).
#' @param ramsey_options  List of additional options passed to
#'   \code{ramsey_policy()}, e.g. \code{list(n_periods = 400, burn_in = 100)}.
#' @param verbose  Print progress messages.
#' @return An object of class \code{"dynhr_solved"} containing:
#'   \describe{
#'     \item{\code{model}}{Parsed \code{dynhr_mod}}
#'     \item{\code{compiled}}{Compiled model}
#'     \item{\code{ss}}{Steady-state result (\code{dynhr_steady})}
#'     \item{\code{dr}}{Decision rules object (\code{DecisionRules},
#'       \code{DecisionRules2}, etc.)}
#'     \item{\code{params}}{Parameter vector used}
#'     \item{\code{irfs}}{IRF list (if \code{stoch_simul = TRUE})}
#'     \item{\code{moments}}{Theoretical moments (if \code{stoch_simul = TRUE})}
#'     \item{\code{ramsey}}{Ramsey result (if \code{ramsey = TRUE})}
#'     \item{\code{meta}}{Run metadata}
#'   }
#'
#' @examples
#' \dontrun{
#' mod <- solve_model("my_model.mod")
#' mod <- solve_model("my_model.mod", order = 2)
#' mod <- solve_model("my_model.mod", ramsey = TRUE)
#' }
#' @param mcp        Logical: if TRUE, parse MCP tags and compute MCP-constrained
#'   IRFs via \code{\link{mcp_solve_path}}.  Default FALSE.
#' @param mcp_options List of additional arguments passed to \code{\link{mcp_solve_path}}
#'   (e.g. \code{list(method = "sparse", max_iter = 100L)}).
#' @param ss_mcp     Logical: if TRUE and the model has MCP tags, solve the steady
#'   state with bound constraints via the semi-smooth MCP solver.
#'   Default FALSE.
#' @seealso \code{\link{run_mode_finding}}, \code{\link{run_posterior_estimation}},
#'   \code{run_all_diagnostics}
#' @export
solve_model <- function(mod_file,
                        params           = NULL,
                        order            = 1L,
                        Sigma_e          = NULL,
                        steady_options   = list(),
                        perturbation_options = list(),
                        stoch_simul      = TRUE,
                        irf_horizons     = 40L,
                        ramsey           = FALSE,
                        ramsey_order     = 1L,
                        ramsey_options   = list(),
                        mcp              = FALSE,
                        mcp_options      = list(),
                        ss_mcp           = FALSE,
                        verbose          = TRUE) {

  .vcat <- function(...) if (verbose) cat(...)

  t_start <- Sys.time()

  .vcat("\n================================================================\n")
  .vcat("  solve_model\n")
  .vcat("================================================================\n\n")

  # -------------------------------------------------------------------
  # Step 1: Parse
  # -------------------------------------------------------------------
  .vcat("-- Step 1: Parse model --\n")
  if (is.character(mod_file) && length(mod_file) == 1) {
    model <- parse_mod(mod_file, verbose = verbose)
  } else if (inherits(mod_file, "dynhr_mod")) {
    model <- mod_file
  } else {
    stop("'mod_file' must be a path to a .mod file or a dynhr_mod object.")
  }
  .vcat(sprintf("  %d endo, %d exo, %d params, %d eqs\n",
                length(model$var_names), length(model$varexo_names),
                length(model$param_names), length(model$equations)))

  # -------------------------------------------------------------------
  # Step 2: Compile
  # -------------------------------------------------------------------
  .vcat("-- Step 2: Compile model --\n")
  max_order <- max(order, if (isTRUE(ramsey)) ramsey_order else 1L)
  compiled <- compile_model(model, verbose = verbose, max_order = max_order)
  .vcat("  compiled OK\n")

  # -------------------------------------------------------------------
  # Step 3: Parameters
  # -------------------------------------------------------------------
  if (is.null(params)) params <- model$param_values

  # -------------------------------------------------------------------
  # Step 4: Steady state
  # -------------------------------------------------------------------
  .vcat("-- Step 3: Steady state --\n")
  so <- modifyList(list(y0 = NULL, tol = 1e-10, max_iter = 1000L),
                   steady_options)
  ss <- solve_steady(compiled, params,
                     y0        = so$y0,
                     endo_names = model$var_names,
                     exo_names  = model$varexo_names,
                     max_iter  = so$max_iter,
                     tol       = so$tol,
                     verbose   = verbose)
  if (!isTRUE(ss$converged)) {
    stop("Steady state did not converge. Check parameter values or initial guess.")
  }
  .vcat(sprintf("  converged=%s, iter=%d, max_res=%.2e\n",
                ss$converged, ss$iterations, ss$max_residual))

  # Parameters assigned INSIDE a steady_state_model block (e.g.
  # `theta = lambda*w/n^xi;`) are computed during steady-state solving and
  # returned in ss$params.  The dynamic Jacobian must see those updated values,
  # not the stale calibrated defaults (e.g. `theta = 0; // set in SS`), or the
  # perturbation will silently use wrong derivatives.  Adopt the SS-updated
  # parameter vector when the steady-state solver supplied one; this is a no-op
  # for models with no steady_state_model param assignments (ss$params is then
  # identical to params, or NULL via the Newton/nleqslv path -> we keep params).
  if (!is.null(ss$params))
    params <- ss$params

  # -------------------------------------------------------------------
  # Step 5: Perturbation
  # -------------------------------------------------------------------
  .vcat(sprintf("-- Step 4: Perturbation (order = %d) --\n", order))
  po <- perturbation_options
  dr <- solve_perturbation(
    model    = model,
    compiled = compiled,
    ss       = ss$values,
    params   = params,
    order    = order,
    Sigma_e  = Sigma_e %||% po$Sigma_e,
    h        = po$h,
    sigma3   = po$sigma3,
    verbose  = verbose
  )
  .vcat(sprintf("  BK satisfied=%s, n_stable=%d\n",
                dr$bk_satisfied, dr$n_stable))

  # -------------------------------------------------------------------
  # Step 5b: Stoch_simul (IRFs + moments)
  # -------------------------------------------------------------------
  irfs    <- NULL
  moments <- NULL
  if (isTRUE(stoch_simul)) {
    .vcat("-- Step 5: stoch_simul (IRFs + moments) --\n")
    irfs <- compute_irfs(dr, model, n_periods = irf_horizons, params = params)
    moments <- compute_moments(dr, model, params = params)
    if (!is.null(irfs))
      .vcat(sprintf("  IRFs: %d shocks x %d periods\n",
                    length(irfs), irf_horizons))
  }

  # -------------------------------------------------------------------
  # Step 5c: MCP-constrained IRFs (optional)
  # -------------------------------------------------------------------
  irfs_mcp <- NULL
  if (isTRUE(mcp)) {
    .vcat("-- Step 5c: MCP-constrained IRFs --\n")
    mcp_specs <- mcp_parse_tags(model, verbose = verbose)
    if (length(mcp_specs) > 0L) {
      mcp_specs <- mcp_resolve_bounds(mcp_specs, model, params)
      mcp_validate_specs(model, mcp_specs, verbose = verbose)
      mo <- modifyList(list(), mcp_options)
      irfs_mcp <- do.call(compute_irfs_mcp, c(list(
        compiled   = compiled,
        y_ss       = ss$values,
        model      = model,
        params     = params,
        mcp_specs  = mcp_specs,
        n_periods  = irf_horizons
      ), mo))
      .vcat(sprintf("  MCP IRFs: %d shocks x %d periods\n",
                    length(irfs_mcp), irf_horizons))
    } else {
      .vcat("  No MCP tags found -- skipping MCP IRFs\n")
    }
  }

  # -------------------------------------------------------------------
  # Step 6: Ramsey optimal policy (optional)
  # -------------------------------------------------------------------
  ramsey_res <- NULL
  if (isTRUE(ramsey)) {
    .vcat("-- Step 6: Ramsey optimal policy --\n")
    ro <- modifyList(list(n_periods = 400L, burn_in = 100L,
                          discount = NULL, planner_objective = NULL),
                     ramsey_options)
    ramsey_res <- ramsey_policy(model, compiled = compiled, params = params,
                                 planner_objective = ro$planner_objective,
                                 order    = ramsey_order,
                                 n_periods = ro$n_periods,
                                 burn_in   = ro$burn_in,
                                 discount  = ro$discount,
                                 verbose  = verbose)
    if (!is.null(ramsey_res))
      .vcat("  Ramsey policy solved.\n")
  }

  # -------------------------------------------------------------------
  # Wrap up
  # -------------------------------------------------------------------
  t_end   <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))

  .vcat(sprintf("\n  Done in %.1f sec\n\n", elapsed))

  result <- list(
    model     = model,
    compiled  = compiled,
    ss        = ss,
    dr        = dr,
    params    = params,
    irfs      = irfs,
    moments   = moments,
    irfs_mcp  = irfs_mcp,
    ramsey    = ramsey_res,
    meta      = list(
      order      = order,
      stoch_simul = isTRUE(stoch_simul),
      mcp        = isTRUE(mcp),
      ramsey     = isTRUE(ramsey),
      elapsed_sec = elapsed
    )
  )
  class(result) <- c("dynhr_solved", "list")
  invisible(result)
}


#' Print method for dynhr_solved
#' @noRd
#' @export
print.dynhr_solved <- function(x, ...) {
  cat("\n<dynhr_solved>\n")
  cat(sprintf("  Model     : %s\n",
              if (!is.null(x$model$source_file)) x$model$source_file else "(in-memory)"))
  cat(sprintf("  Variables : %d endo, %d exo\n",
              length(x$model$var_names), length(x$model$varexo_names)))
  cat(sprintf("  Parameters: %d\n", length(x$params)))
  conv <- if (isTRUE(x$ss$converged)) "converged" else "NOT converged"
  cat(sprintf("  Steady    : %s  (max|r| = %.2e)\n", conv, x$ss$max_residual))
  if (!is.null(x$dr)) {
    cat(sprintf("  Perturbation: order=%d, BK=%s, n_stable=%d\n",
                x$meta$order %||% 1L,
                if (isTRUE(x$dr$bk_satisfied)) "satisfied" else "VIOLATED",
                x$dr$n_stable %||% 0L))
  }
  if (!is.null(x$irfs))
    cat(sprintf("  IRFs      : %d shocks\n", length(x$irfs)))
  if (!is.null(x$moments))
    cat(sprintf("  Moments   : computed\n"))
  if (!is.null(x$ramsey))
    cat("  Ramsey    : solved\n")
  cat(sprintf("  Elapsed   : %.1f sec\n", x$meta$elapsed_sec %||% NA))
  invisible(x)
}
