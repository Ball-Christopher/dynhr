## R/obc-simulate.R
## --------------------------------------------------------------------------
## OBC deterministic forward simulation.
##
## Provides:
##   obc_simulate() -- forward simulation under a given regime path
##
## Used primarily for post-estimation scenario analysis and OBC-aware IRF
## comparison against Dynare / Julia PATHSolver output.
## --------------------------------------------------------------------------


#' OccBin forward simulation under a given regime path
#'
#' Simulates the OBC model forward from an initial state (default zero) under
#' the piecewise-linear OccBin rules of \code{regime_path}: the rule of period
#' t comes from the backward recursion over the regimes of periods t, t+1,
#' ... (slack after the last binding period), so a binding period that is
#' followed by another binding period is NOT solved with the one-period
#' binding policy of \code{obc_solve_binding()} (which assumes the next
#' period slack).  Shocks are surprises.
#'
#' y_t = ghx_t * s_{t-1} + ghu_t * eps_t + c_t
#'
#' @param shock_seq   n_exo x T numeric matrix of shock values at each period
#' @param dr_slack    Slack-regime DecisionRules (from solve_perturbation)
#' @param sys         System matrices (from extract_system_matrices)
#' @param specs       OBC spec list (from obc_parse_tags)
#' @param regime_path Integer vector (length T): bitfield regime per period
#'                    (0 = slack; bit j set = spec j binds)
#' @param state_init  Initial state (length n_state); NULL = zero
#' @return n_endo x T matrix of endogenous variable paths (deviations);
#'         rownames = endo names
#' @noRd
obc_simulate <- function(shock_seq, dr_slack, sys, specs, regime_path,
                         state_init = NULL) {
  ctx   <- .obc_pwl_context(sys, dr_slack, specs)
  rules <- .obc_pwl_rules(ctx, regime_path)
  .obc_pwl_forward(ctx, rules, shock_seq, state_init)$paths
}
