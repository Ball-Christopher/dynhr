## R/obc-decomposition.R
## --------------------------------------------------------------------------
## OBC historical decomposition.
##
## Provides:
##   historical_decomposition_obc() -- per-shock + constraint contributions
##
## Decomposes the smoothed endogenous path (from pkf_smoother_obc) into
## contributions from each structural shock plus a "constraint" component
## capturing the cumulative effect of binding-regime intercepts.
##
## Full pipeline:
##   kalman_filter_obc_pkf(return_store=TRUE)
##     -> pkf_smoother_obc()
##     -> historical_decomposition_obc()
## --------------------------------------------------------------------------


#' Historical decomposition for OBC piecewise-linear Kalman filter
#'
#' Decomposes the path of ALL endogenous variables into contributions from
#' each structural shock plus a "constraint" component capturing the
#' cumulative effect of binding-regime intercepts.
#'
#' For each shock j the forward recursion is run with regime-specific matrices:
#'
#'   \eqn{s_t^{(j)} = TT_t \cdot s_{t-1}^{(j)} + RR_t[:,j] \cdot \varepsilon_{j,t|T}}
#'   \eqn{x_t^{(j)} = ghx_t \cdot s_{t-1}^{(j)} + ghu_t[:,j] \cdot \varepsilon_{j,t|T}}
#'
#' The "constraint" contribution carries the accumulated binding constant:
#'
#'   \eqn{s_t^{(c)} = TT_t \cdot s_{t-1}^{(c)} + c\_state_t}
#'   \eqn{x_t^{(c)} = ghx_t \cdot s_{t-1}^{(c)} + c\_full_t}
#'
#' Summing all contributions reconstructs the smoothed endogenous path
#' (up to the initial-condition contribution from \eqn{s_{0|T} \neq 0}):
#'
#'   \eqn{\sum_j x_t^{(j)} + x_t^{(c)} = ghx_t \cdot s_{t-1} + ghu_t \cdot \varepsilon_t + c\_full_t}
#'
#' Initial conditions: the per-shock recursions start from s = 0, so on their
#' own they miss the effect of a smoothed initial state \eqn{s_{0|T} \neq 0}.
#' The \code{"initial"} component supplies it -- the shock-free, intercept-free
#' trajectory started at \code{s0} and propagated with the SAME regime-specific
#' transition matrices along \code{regime_path}:
#'
#'   \eqn{s_t^{(0)} = TT_t \cdot s_{t-1}^{(0)}, \quad
#'        x_t^{(0)} = ghx_t \cdot s_{t-1}^{(0)}}
#'
#' With \code{s0} supplied the adding-up is exact at every t, not just once
#' the initial condition has decayed. \code{s0 = NULL} (the default) leaves it
#' identically zero, i.e. the pre-0.9.2 output.
#'
#' @param smoothed_shocks  n_exo x T matrix from pkf_smoother_obc()
#'                         (rows = shocks, columns = time periods)
#' @param regime_path      Integer vector length T (accepted regime per period)
#' @param regime_cache     R environment with per-regime policies; must have
#'                         c_full stored (populated by obc_ensure_policy)
#' @param dr_slack         Slack-regime DecisionRules
#' @param s0               Numeric length-n_state smoothed initial state, or
#'                         NULL (zero). State ordering is
#'                         \code{dr_slack$state_idx}.
#' @param shock_groups     Optional named list grouping shocks (see
#'                         \code{\link{historical_decomposition}}). The
#'                         \code{constraint} and \code{initial} components are
#'                         never grouped.
#' @param model            Optional parsed model, used only as the fallback
#'                         source of \code{shock_groups}.
#' @return List with:
#'   $contributions  Named list of n_endo x T matrices, one per shock (or
#'                   group) plus $constraint for the binding-constant
#'                   contribution and $initial for the initial condition
#'   $total          n_endo x T matrix: sum of all contributions
#'   $endo_names     Character vector
#'   $exo_names      Character vector (with "constraint", "initial" last)
#' @export
historical_decomposition_obc <- function(smoothed_shocks, regime_path,
                                          regime_cache, dr_slack,
                                          s0 = NULL, shock_groups = NULL,
                                          model = NULL) {

  n_exo   <- nrow(smoothed_shocks)
  n_T     <- ncol(smoothed_shocks)
  n_state <- length(dr_slack$state_idx)
  n_endo  <- nrow(dr_slack$ghx)
  exo_names  <- if (!is.null(rownames(smoothed_shocks))) rownames(smoothed_shocks)
                else dr_slack$exo_names
  endo_names <- dr_slack$endo_names

  if (length(regime_path) != n_T)
    stop(sprintf(
      "historical_decomposition_obc: regime_path length (%d) != n_T (%d)",
      length(regime_path), n_T
    ))

  # Initialize contributions: one n_endo x T per shock + "constraint"
  contributions <- c(
    setNames(lapply(seq_len(n_exo), function(j) matrix(0, n_endo, n_T)),
             exo_names),
    list(constraint = matrix(0, n_endo, n_T))
  )

  # ---- Per-shock contributions ---------------------------------------------
  for (j in seq_len(n_exo)) {
    s_j <- numeric(n_state)
    for (t in seq_len(n_T)) {
      pol   <- get(as.character(regime_path[t]), envir = regime_cache,
                   inherits = FALSE)
      ghx_t <- pol$dr$ghx    # n_endo x n_state
      ghu_t <- pol$dr$ghu    # n_endo x n_exo
      TT_t  <- pol$TT        # n_state x n_state
      RR_t  <- pol$RR        # n_state x n_exo
      eps_j <- smoothed_shocks[j, t]

      # Contribution to all endo vars from shock j at period t
      contributions[[j]][, t] <- drop(ghx_t %*% s_j) + ghu_t[, j] * eps_j

      # Advance state attributable to shock j
      s_j <- drop(TT_t %*% s_j) + RR_t[, j] * eps_j
    }
    rownames(contributions[[j]]) <- endo_names
  }

  # ---- Constraint contribution (accumulated binding-regime intercepts) -----
  s_c <- numeric(n_state)
  for (t in seq_len(n_T)) {
    pol     <- get(as.character(regime_path[t]), envir = regime_cache,
                   inherits = FALSE)
    ghx_t   <- pol$dr$ghx
    TT_t    <- pol$TT
    c_full  <- pol$c_full    # n_endo (stored by obc_ensure_policy; zero for slack)
    c_state <- pol$c_state   # n_state

    # Contribution: accumulated constant in state-space + direct intercept
    contributions$constraint[, t] <- drop(ghx_t %*% s_c) + c_full

    # Advance constant-driven state accumulator
    s_c <- drop(TT_t %*% s_c) + c_state
  }
  rownames(contributions$constraint) <- endo_names

  # ---- Initial-condition contribution (regime-specific propagation) --------
  # Shock-free AND intercept-free: only the initial state is propagated, with
  # the same TT_t / ghx_t the accepted regime_path selects at each period.
  if (is.null(s0)) s0 <- numeric(n_state)
  s0 <- as.numeric(s0)
  if (length(s0) != n_state)
    stop(sprintf(
      "historical_decomposition_obc: s0 has length %d, expected %d.",
      length(s0), n_state), call. = FALSE)

  initial <- matrix(0, n_endo, n_T)
  s_i     <- s0
  for (t in seq_len(n_T)) {
    pol   <- get(as.character(regime_path[t]), envir = regime_cache,
                 inherits = FALSE)
    initial[, t] <- drop(pol$dr$ghx %*% s_i)
    s_i          <- drop(pol$TT %*% s_i)
  }
  rownames(initial) <- endo_names

  # ---- Optional shock groups (shock columns only) --------------------------
  groups <- .resolve_shock_groups(shock_groups, model, exo_names)
  if (!is.null(groups)) {
    shock_part    <- contributions[exo_names]
    grouped       <- .group_contributions(shock_part, groups)
    contributions <- c(grouped, contributions["constraint"])
    exo_names     <- names(grouped)
  }

  contributions$initial <- initial

  # ---- Sum to total --------------------------------------------------------
  total <- Reduce(`+`, contributions)
  rownames(total) <- endo_names

  structure(
    list(
      contributions = contributions,
      total         = total,
      initial       = initial,
      s0            = s0,
      endo_names    = endo_names,
      exo_names     = c(exo_names, "constraint", "initial"),
      components    = names(contributions),
      shock_groups  = groups,
      orientation   = "endo_time"
    ),
    class = c("dynhr_shock_decomposition", "list")
  )
}
