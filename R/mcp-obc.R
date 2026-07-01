## R/mcp-obc.R
## --------------------------------------------------------------------------
## Phase 4: Bridge between MCP semi-smooth Newton solver and the existing
## OBC (Occasionally Binding Constraints) linear-model infrastructure.
##
## Provides:
##   mcp_compute_irfs_obc() — IRF computation using the MCP path solver,
##     callable as obc_solver = "mcp" in compute_irfs_obc().
##
## The MCP solver uses the full nonlinear compiled model, so it works for
## both linear and nonlinear models.  On linear models, it should produce
## paths that match boehl_solve_regime_path() and solve_obc_lcp() to within
## numerical tolerance (1e-8), serving as an independent cross-check.
##
## REFERENCES
##   Same as mcp-solve.R — Fischer-Burmeister semi-smooth Newton.
## --------------------------------------------------------------------------


# =============================================================================
# Build a compiled model + MCP specs from a linear OBC model setup
# =============================================================================

#' Build compiled model and MCP specs from OBC model components
#'
#' Given the inputs typically available in \code{compute_irfs_obc()} (a model,
#' system matrices, and OBC specs), builds the compiled dynamic model and
#' MCP specs needed by \code{mcp_solve_path()}.
#'
#' This bridges the linear-model OBC infrastructure (which uses DecisionRules
#' and system matrices) with the nonlinear MCP solver (which needs a
#' dynhr_compiled with residual/jacobian functions).
#'
#' @param model   dynhr_mod (from \code{\link{parse_mod}})
#' @param specs   OBC spec list (from \code{obc_parse_tags} or equivalent)
#' @param params  Named numeric parameter vector
#' @param verbose Logical: print compilation messages (default FALSE)
#' @return List with:
#'   \item{compiled}{dynhr_compiled}
#'   \item{mcp_specs}{MCP spec list (may be converted from OBC specs)}
#'   \item{y_ss}{Named numeric: steady state (zero for linear models)}
#' @noRd
.mcp_build_from_obc <- function(model, specs, params, verbose = FALSE) {
  # Compile the model (needed for mcp_solve_path)
  compiled <- compile_model(model, verbose = verbose)

  # Convert OBC specs to MCP specs if needed
  # (OBC specs have $eq_idx, $var_idx, $var_name, $op, $bound — same as MCP specs)
  if (length(specs) > 0L && is.null(specs[[1]]$tag_type)) {
    # Tag as MCP specs for uniform handling
    for (j in seq_along(specs)) {
      specs[[j]]$tag_type <- "obc_compat"
    }
  }

  # Steady state (zero for linear(deviation) models)
  y_ss <- rep(0, length(model$var_names))
  names(y_ss) <- model$var_names

  list(compiled = compiled, mcp_specs = specs, y_ss = y_ss)
}


# =============================================================================
# IRF computation using MCP solver, callable from compute_irfs_obc
# =============================================================================

#' Compute IRFs using the MCP path solver (bridge function)
#'
#' Wraps \code{mcp_solve_path()} for use as \code{obc_solver = "mcp"} in
#' \code{\link{compute_irfs_obc}}.  Unlike the Boehl and LCP solvers, this
#' works directly on the compiled nonlinear model and does NOT require
#' \code{model(linear)}.
#'
#' For linear models, the results should match \code{boehl_solve_regime_path}
#' and \code{solve_obc_lcp} to within numerical tolerance (1e-8).  This is
#' verified in the parity tests (\code{test-mcp-path.R}).
#'
#' @param dr_slack    Slack-regime DecisionRules (used for endo/exo names)
#' @param model       dynhr_mod (from \code{\link{parse_mod}})
#' @param sys         System matrices (from \code{extract_system_matrices_fast})
#'   — not used by the MCP solver, included for API compatibility.
#' @param specs       OBC spec list (from \code{obc_parse_tags})
#' @param obs_idx     Integer vector: observable positions (not used by MCP solver)
#' @param n_periods   Number of IRF periods (default 40)
#' @param shock_size  Shock size multiplier (default 1)
#' @param params      Named numeric parameter vector; NULL → model$param_values
#' @param compiled    Optional pre-compiled model.  If NULL, compiles from model.
#' @param ...         Additional arguments passed to \code{\link{mcp_solve_path}}
#' @return List with:
#'   \item{paths}{n_endo × n_periods matrix of endogenous variable paths}
#'   \item{regime_path}{Integer vector (length n_periods): per-period binding}
#'   \item{converged}{Logical: TRUE if MCP solver converged}
#' @noRd
mcp_compute_irfs_obc <- function(dr_slack, model, sys, specs,
                                  obs_idx      = NULL,
                                  n_periods    = 40L,
                                  shock_size   = 1,
                                  params       = NULL,
                                  compiled     = NULL,
                                  ...) {
  if (is.null(params)) params <- model$param_values

  endo <- dr_slack$endo_names
  exo  <- dr_slack$exo_names

  # Build or use provided compiled model
  if (is.null(compiled)) {
    bridge <- .mcp_build_from_obc(model, specs, params, verbose = FALSE)
    compiled <- bridge$compiled
    mcp_specs <- bridge$mcp_specs
    y_ss <- bridge$y_ss
  } else {
    # Use the MCP specs as-is (already available)
    mcp_specs <- specs
    y_ss <- rep(0, length(endo))
    names(y_ss) <- endo
  }

  shock_stderr <- .get_shock_stderr(model, exo, params)

  # Run MCP solver for each shock
  res_list <- vector("list", length(exo))
  names(res_list) <- exo

  for (k in seq_along(exo)) {
    shock_seq <- matrix(0, nrow = n_periods, ncol = length(exo))
    colnames(shock_seq) <- exo
    shock_seq[1L, exo[k]] <- shock_stderr[exo[k]] * shock_size

    y0_num <- rep(0, length(endo))
    names(y0_num) <- endo

    res <- mcp_solve_path(
      compiled   = compiled,
      y0         = y0_num,
      y_ss       = y0_num,
      shock_path = shock_seq,
      params     = params,
      mcp_specs  = mcp_specs,
      ...
    )

    # Extract path in the format expected by compute_irfs_obc
    paths <- t(res$Y)   # n_endo × n_periods
    rownames(paths) <- endo

    res_list[[k]] <- list(
      paths       = paths,
      regime_path = res$active_set,
      converged   = res$converged
    )
  }

  # Build combined result (for first shock, mainly)
  list(
    paths       = res_list[[1]]$paths,
    regime_path = res_list[[1]]$regime_path,
    converged   = all(vapply(res_list, `[[`, logical(1), "converged")),
    all_results = res_list
  )
}


# =============================================================================
# MCP-aware decision rule construction for linear models
# =============================================================================

#' Build regime-specific decision rules from MCP solver active set
#'
#' For a linear model, runs the MCP path solver to identify which constraints
#' bind at each period, then constructs regime-specific policy matrices
#' (ghx, ghu, c_state, c_full, TT, RR) for each unique regime.
#'
#' This mirrors what \code{obc_ensure_policy()} does in the Boehl/LCP
#' framework, but uses the MCP solver's active-set output instead of
#' OccBin binding-system construction.
#'
#' @param dr_slack Slack-regime DecisionRules (from \code{solve_perturbation})
#' @param model    dynhr_mod (from \code{\link{parse_mod}})
#' @param sys      System matrices (from \code{extract_system_matrices_fast})
#' @param specs    OBC spec list (from \code{obc_parse_tags})
#' @param params   Named numeric parameter vector; NULL → model$param_values
#' @param compiled Pre-compiled model (optional).  If NULL, compiles from model.
#' @param ...      Additional arguments passed to \code{mcp_solve_path}
#' @return List (same structure as \code{obc_ensure_policy} cache entries):
#'   Each element is a regime-specific list with:
#'   \item{dr}{DecisionRules-like list with ghx, ghu}
#'   \item{c_state}{Constant term for state equation}
#'   \item{c_full}{Constant term for full endogenous equation}
#'   \item{TT}{State transition matrix}
#'   \item{RR}{Shock impact matrix for states}
#' @export
#'
#' @examples
#' \dontrun{
#' model <- parse_mod("nk_2obc.mod")
#' compiled <- compile_model(model)
#' specs <- obc_parse_tags(model)
#' sys <- extract_system_matrices_fast(cache_system_structure(compiled), ss, params)
#' dr_slack <- solve_perturbation(model, compiled, ss, params)
#' dr_map <- mcp_solve_dr(dr_slack, model, sys, specs, params)
#' }
mcp_solve_dr <- function(dr_slack, model, sys, specs,
                          params = NULL, compiled = NULL, ...) {
  if (is.null(params)) params <- model$param_values
  if (is.null(compiled)) {
    bridge <- .mcp_build_from_obc(model, specs, params, verbose = FALSE)
    compiled <- bridge$compiled
  }

  endo <- dr_slack$endo_names
  exo  <- dr_slack$exo_names
  n_endo <- length(endo)
  n_spec <- length(specs)

  # Use a moderate shock to identify the active set
  # (Demand shock that triggers binding)
  n_periods <- 40L
  shock_seq <- matrix(0, nrow = n_periods, ncol = length(exo))
  colnames(shock_seq) <- exo

  # Use whichever exo shock first — if there's a demand/monetary shock, use it
  shock_idx <- 1L
  if (length(exo) >= 1L) {
    shock_seq[shock_idx, 1L] <- -0.10  # moderate demand shock
  }

  y0_num <- rep(0, n_endo)
  names(y0_num) <- endo

  res <- mcp_solve_path(
    compiled   = compiled,
    y0         = y0_num,
    y_ss       = y0_num,
    shock_path = shock_seq,
    params     = params,
    mcp_specs  = specs,
    ...
  )

  # Extract unique regimes from active set
  unique_regimes <- unique(res$active_set)
  n_regimes <- length(unique_regimes)

  dr_map <- vector("list", n_regimes)
  names(dr_map) <- as.character(unique_regimes)

  for (r in seq_len(n_regimes)) {
    regime_bits <- unique_regimes[r]
    flags <- obc_regime_flags(regime_bits, n_spec)

    # Build binding system for this regime
    # Active specs: those with flag = TRUE
    active_specs <- specs[flags]
    slack_specs  <- specs[!flags]

    # Build binding-system matrix using the existing infrastructure
    if (length(active_specs) > 0L) {
      binding_sys <- obc_build_binding_sys(sys, active_specs)
      dr_bind <- .solve_from_system(binding_sys$sys_b, model, compiled,
                                     rep(0, n_endo), params, verbose = FALSE)
    } else {
      dr_bind <- dr_slack
    }

    # Extract policy matrices
    si <- dr_bind$state_idx
    n_state <- length(si)

    dr_map[[r]] <- list(
      dr      = dr_bind,
      c_state = if (!is.null(dr_bind$c_state)) dr_bind$c_state else numeric(n_state),
      c_full  = if (!is.null(dr_bind$c_full)) dr_bind$c_full else numeric(n_endo),
      TT      = dr_bind$ghx[si, , drop = FALSE],
      RR      = dr_bind$ghu[si, , drop = FALSE],
      ZZ      = if (!is.null(dr_bind$ghx)) dr_bind$ghx else NULL,
      DD      = if (!is.null(dr_bind$ghu)) dr_bind$ghu else NULL
    )
  }

  attr(dr_map, "active_set") <- res$active_set
  attr(dr_map, "unique_regimes") <- unique_regimes
  dr_map
}
