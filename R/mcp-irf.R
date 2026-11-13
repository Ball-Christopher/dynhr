## R/mcp-irf.R
## --------------------------------------------------------------------------
## IRF computation using the MCP semi-smooth Newton path solver.
##
## Provides:
##   compute_irfs_mcp()  -- MCP-aware IRF computation (IRFCollection output)
##
## Unlike compute_irfs_obc(), this does NOT require model(linear) or
## separate slack/binding decision rules.  It operates directly on the
## compiled nonlinear dynamic model with MCP constraints.
##
## REFERENCES
##   Guerrieri, L., & Iacoviello, M. (2015). "OccBin: A toolkit for solving
##     dynamic models with occasionally binding constraints easily."
##     Journal of Monetary Economics, 70, 22-38.
## --------------------------------------------------------------------------


# =============================================================================
# IRF computation using MCP path solver
# =============================================================================

## Internal perfect-foresight horizon for an IRF of `n_periods` periods. A
## perfect-foresight path solved over exactly `n_periods` imposes the return to
## steady state AT n_periods, so the terminal condition bites inside the
## reported window and the first k periods change with n_periods. Solving on
## max(2 * n_periods, n_periods + 40) periods and truncating keeps the reported
## periods independent of n_periods. Shared by compute_irfs_mcp() and
## compute_irfs_obc(obc_solver = "mcp").
.pf_irf_horizon <- function(n_periods)
  as.integer(max(2L * n_periods, n_periods + 40L))

#' Compute IRFs using the MCP semi-smooth Newton path solver
#'
#' For each shock in the model, constructs a unit-impulse shock sequence
#' (shock hits at t=1, zero thereafter), runs the MCP path solver to
#' find the nonlinear path with complementarity constraints, and packages
#' the result as an IRFCollection (same structure as \code{\link{compute_irfs}}).
#'
#' This is the MCP analogue of \code{\link{compute_irfs_obc}} for nonlinear
#' models.  Unlike the OBC solvers, it does NOT require \code{model(linear)}
#' and handles both standard MCP tags (\code{[mcp = 'var OP bound']}) and
#' OccBin bind/relax tags.
#'
#' **Dynare parity:** For a model with \code{[mcp = '...']} tags, this should
#' produce the same IRFs as Dynare's \code{perfect_foresight_solver} with
#' \code{stack_solve_algo = 7} (PATH solver) to within 1e-6 (linear models)
#' or 1e-4 (nonlinear models).
#'
#' @param compiled    dynhr_compiled (from \code{\link{compile_model}})
#' @param y_ss        Named numeric vector: steady state (used as initial
#'   condition and terminal condition)
#' @param model       dynhr_mod (for shock standard deviations and variable names)
#' @param params      Named numeric parameter vector; NULL -> model$param_values
#' @param mcp_specs   List of MCP specs (from \code{\link{mcp_parse_tags}}).
#'   If NULL, attempts to parse MCP tags from the model automatically.
#' @param n_periods   Number of IRF periods (default 40). The
#'   perfect-foresight problem is solved on an internal horizon of
#'   \code{max(2 * n_periods, n_periods + 40)} periods (return to steady state
#'   imposed at its end) and truncated to \code{n_periods}, so the returned
#'   periods do not depend on \code{n_periods} (the rule
#'   \code{\link{compute_irfs_obc}} uses for \code{obc_solver = "mcp"}).
#' @param shock_size  Shock size multiplier applied to each shock's standard
#'   deviation (default 1: one-standard-deviation impulse)
#' @param ...         Additional arguments passed to \code{\link{mcp_solve_path}}
#' @return An object of class \code{IRFCollection}: a named list (one entry
#'   per exogenous shock).  Each entry is an \code{n_periods x n_endo} matrix.
#'   Attributes: \code{n_periods}, \code{endo_names}, \code{exo_names},
#'   \code{mcp_specs}.
#' @export
#'
#' @examples
#' \dontrun{
#' model <- parse_mod("nk_zlb_dynare.mod")
#' compiled <- compile_model(model)
#' ss <- solve_steady_state(model, compiled, model$param_values)
#' specs <- mcp_parse_tags(model)
#' irfs <- compute_irfs_mcp(compiled, ss$values, model,
#'                           model$param_values, specs)
#' # Compare with Dynare's PATH solver:
#' #   > dynare nk_zlb_dynare.mod
#' #   > [irfs_MCP, irfs_PATH] = compare_irfs(irfs, dynare_irfs)
#' #   > max(abs(irfs_MCP - irfs_PATH))  # should be < 1e-4
#' }
compute_irfs_mcp <- function(compiled,
                              y_ss,
                              model,
                              params      = NULL,
                              mcp_specs   = NULL,
                              n_periods   = 40L,
                              shock_size  = 1,
                              ...) {
  n_periods <- .check_irf_n_periods(n_periods, "compute_irfs_mcp")
  if (is.null(params)) params <- model$param_values

  endo  <- compiled$dynamic$endo_names
  exo   <- compiled$dynamic$exo_names
  n_exo <- length(exo)

  # Parse MCP specs if not provided
  if (is.null(mcp_specs)) {
    mcp_specs <- mcp_parse_tags(model, verbose = FALSE)
    if (length(mcp_specs) > 0L) {
      mcp_specs <- mcp_resolve_bounds(mcp_specs, model, params)
      mcp_validate_specs(model, mcp_specs, verbose = FALSE)
    }
  }

  n_spec <- length(mcp_specs)

  # Impulses: lower Cholesky of the full Sigma_e (shared IRF rule)
  E_imp <- .irf_impulse_matrix(NULL, model, params, shock_size, exo = exo)

  irfs <- vector("list", n_exo)
  names(irfs) <- exo
  n_solve <- .pf_irf_horizon(n_periods)

  y0_num <- as.numeric(y_ss[endo])
  names(y0_num) <- endo

  for (k in seq_along(exo)) {
    # Build unit impulse: shock hits at t=1, zero thereafter
    shock_seq <- matrix(0, nrow = n_solve, ncol = n_exo)
    colnames(shock_seq) <- exo
    shock_seq[1L, ] <- E_imp[, k]

    # Run MCP solver
    res <- mcp_solve_path(
      compiled   = compiled,
      y0         = y0_num,
      y_ss       = y0_num,
      shock_path = shock_seq,
      params     = params,
      mcp_specs  = mcp_specs,
      ...
    )

    # Extract IRF matrix
    irf_mat <- res$Y[seq_len(n_periods), , drop = FALSE]
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))
    irfs[[k]] <- irf_mat
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")  <- n_periods
  attr(irfs, "endo_names") <- endo
  attr(irfs, "exo_names")  <- exo
  attr(irfs, "mcp_specs")  <- mcp_specs
  irfs
}
