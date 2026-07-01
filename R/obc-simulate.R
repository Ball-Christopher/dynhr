## R/obc-simulate.R
## --------------------------------------------------------------------------
## OBC deterministic forward simulation.
##
## Provides:
##   obc_simulate() -- forward simulation switching between slack/binding policy
##
## Used primarily for post-estimation scenario analysis and OBC-aware IRF
## comparison against Dynare / Julia PATHSolver output.
## --------------------------------------------------------------------------


#' OccBin forward simulation
#'
#' Simulates the OBC model forward from a zero initial state, switching
#' between slack and binding policy matrices at each period according to
#' regime_path. Useful for post-estimation scenario analysis and OBC-aware
#' IRF comparison against Julia PF output.
#'
#' y_t = ghx_regime * s_{t-1} + ghu_regime * eps_t + c_full_regime
#' s_t = T_regime   * s_{t-1} + R_regime   * eps_t + c_state_regime
#'
#' @param shock_seq   n_exo x T numeric matrix of shock values at each period
#' @param dr_slack    Slack-regime DecisionRules (from solve_perturbation)
#' @param dr_bind     Binding-regime DecisionRules (from obc_solve_binding)
#' @param c_state     Binding-regime state constant (length n_state)
#' @param c_full      Binding-regime full-endo constant (length n_endo)
#' @param regime_path Integer vector (length T): 0 = slack, nonzero = binding
#' @return n_endo x T matrix of endogenous variable paths; rownames = endo names
#' @noRd
obc_simulate <- function(shock_seq, dr_slack, dr_bind, c_state, c_full, regime_path) {
  si    <- dr_slack$state_idx
  n_endo <- nrow(dr_slack$ghx)
  n_T    <- ncol(shock_seq)

  gs_st <- dr_slack$ghx[si, , drop = FALSE]
  gu_st <- dr_slack$ghu[si, , drop = FALSE]
  gb_st <- dr_bind$ghx[si, , drop = FALSE]
  ub_st <- dr_bind$ghu[si, , drop = FALSE]

  s     <- numeric(length(si))
  paths <- matrix(NA_real_, n_endo, n_T)
  rownames(paths) <- dr_slack$endo_names

  for (t in seq_len(n_T)) {
    eps <- shock_seq[, t]
    if (regime_path[t] != 0L) {
      paths[, t] <- drop(dr_bind$ghx %*% s) + drop(dr_bind$ghu %*% eps) + c_full
      s          <- drop(gb_st %*% s)        + drop(ub_st %*% eps)       + c_state
    } else {
      paths[, t] <- drop(dr_slack$ghx %*% s) + drop(dr_slack$ghu %*% eps)
      s          <- drop(gs_st %*% s)         + drop(gu_st %*% eps)
    }
  }
  paths
}
