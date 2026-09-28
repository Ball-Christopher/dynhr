## R/obc-regime.R
## --------------------------------------------------------------------------
## OBC regime encoding, lazy policy cache, and OccBin guess-and-verify.
##
## Provides:
##   obc_should_bind()    -- zero-shock regime check (state-based, no KF info)
##   obc_regime_idx()     -- logical bind-flags -> integer bitfield regime index
##   obc_regime_flags()   -- integer regime index -> logical bind-flags
##   obc_ensure_policy()  -- lazily compute and cache per-regime policy matrices
##   obc_guess_verify()   -- regime path of the data via the OccBin PKF
##   .obc_warn_me_floor_lock() -- me-floor hazard guard using the slack regime
##                                as the stationary-F proxy (see below)
##
## Regime encoding: integer bitfield where bit j (0-based) = 1 means spec j
## binds.  Regime 0 = all slack.  For k=1 this degenerates to the legacy 0/1
## encoding used in the original kalman_filter_obc.
##
## Policy cache: an R environment keyed by character(regime_idx).  Each entry
## stores pre-sliced state-space matrices {dr, c_state, c_obs, c_full, TT, RR,
## ZZ, DD} so the Kalman filter can look them up in O(1).
## --------------------------------------------------------------------------

# =============================================================================
# me-floor hazard guard (OBC family)
# =============================================================================

#' Warn-only me_variance-floor hazard guard for the OBC filter family
#'
#' The main Gaussian \code{kalman_filter()} guards its default \code{me_variance}
#' floor against a near-degenerate model-implied innovation covariance (see
#' \code{.pruned_me_floor_ratio()} / \code{.warn_me_floor_lock()} in
#' R/pruned-state-space.R). The OBC filters (\code{kalman_filter_obc()},
#' \code{kalman_filter_obc_pkf()}) are regime-switching, so there is no single
#' stationary \code{F}; instead this uses the SLACK regime (no constraint
#' binding) as the stationary-F proxy -- it is the reference regime whose
#' \code{TT}/\code{ZZ}/\code{QQ}/\code{HH} matrices are exactly the plain
#' first-order linear solution, so \code{.pruned_me_floor_ratio()} applies to
#' it directly (same detector, same threshold, same message).
#'
#' Warn-only: never changes control flow or any numeric result. Guarded by
#' \code{getOption("dynhr.me_floor_check", TRUE)} and the caller-supplied
#' \code{check} flag (used for the once-per-closure latch in the OBC posterior
#' factories, mirroring \code{make_log_posterior()}).
#'
#' @param dr_slack    Slack-regime DecisionRules (from \code{.solve_from_system}
#'                    or \code{solve_perturbation}).
#' @param model       dynhr_mod
#' @param params      Named numeric parameter vector
#' @param obs_vars    Character vector of observed variable names
#' @param obs_idx     Integer vector: observable positions in endo vector
#' @param me_variance Scalar measurement-error floor actually in force
#' @param check       Logical: run the guard at all (already combines the
#'                     caller's once-per-closure latch and the global option)
#' @return \code{invisible(NULL)}; called for the warning side-effect only
#' @noRd
.obc_warn_me_floor_lock <- function(dr_slack, model, params, obs_vars, obs_idx,
                                     me_variance, check) {
  if (!isTRUE(check) || !is.finite(me_variance) || me_variance <= 0)
    return(invisible(NULL))

  si  <- dr_slack$state_idx
  exo <- dr_slack$exo_names
  TT  <- dr_slack$ghx[si,      , drop = FALSE]
  RR  <- dr_slack$ghu[si,      , drop = FALSE]
  ZZ  <- dr_slack$ghx[obs_idx, , drop = FALSE]
  DD  <- dr_slack$ghu[obs_idx, , drop = FALSE]

  Sigma_e <- .get_shock_cov(model, exo, params)
  QQ <- tcrossprod(RR %*% Sigma_e, RR)
  HH <- tcrossprod(DD %*% Sigma_e, DD)
  SS <- RR %*% Sigma_e %*% t(DD)

  Sxi0 <- tryCatch(solve_lyapunov(TT, QQ),
                   error = function(e) .dynhr_reraise_bug(e, NULL))
  if (is.null(Sxi0) || !all(is.finite(Sxi0))) return(invisible(NULL))

  .warn_me_floor_lock(
    .pruned_me_floor_ratio(TT, ZZ, QQ, HH, SS, Sxi0, me_variance),
    obs_vars, me_variance)
  invisible(NULL)
}


# =============================================================================
# Regime detection (zero-shock, state-based)
# =============================================================================

#' Check whether each OBC would bind given the previous filtered state
#'
#' Predicts the constrained variable's value at time t by applying the slack
#' decision rule to state_{t-1} (with zero shock expectation). Returns TRUE
#' for each spec where the predicted value violates the bound.
#'
#' For a lower bound (>): violation when predicted value < bound.
#' For an upper bound (<): violation when predicted value > bound.
#'
#' Note: this zero-shock, one-period, slack-rule check is NOT the OccBin
#' regime check (it ignores the shock and every later period); since W49
#' (0.9.3.93) neither obc_guess_verify() nor the PKF uses it.
#'
#' @param state_prev Numeric vector: filtered state at t-1 (length n_state)
#' @param dr_slack   Slack-regime DecisionRules
#' @param specs      OBC spec list
#' @return Logical vector (length = length(specs))
#' @noRd
obc_should_bind <- function(state_prev, dr_slack, specs) {
  ## `ghx %*% state_prev` is a DEVIATION from steady state while `s$bound` is a
  ## LEVEL (the package convention, ledger A6) -- compare against the
  ## deviation-form bound.  Before 0.9.4 this used the raw level and disagreed
  ## with obc_solve_binding()/boehl whenever the steady state was non-zero.
  bnd <- .obc_bound_dev(specs, dr_slack)
  vapply(seq_along(specs), function(j) {
    s <- specs[[j]]
    var_pred <- sum(dr_slack$ghx[s$var_idx, ] * state_prev)
    if (s$op == ">") var_pred < bnd[[j]] else var_pred > bnd[[j]]
  }, logical(1))
}


# =============================================================================
# Regime encoding
# =============================================================================

#' Encode a logical binding vector as an integer regime index (0-based bitfield)
#'
#' Regime 0 = all constraints slack.  Bit j (0-based) = 1 means spec j binds.
#' For k=1 this degenerates to 0 (slack) / 1 (binding), matching the legacy path.
#'
#' @param bind_flags Logical vector (length k) — TRUE where constraint j binds
#' @return Integer regime index in [0, 2^k - 1]
#' @noRd
obc_regime_idx <- function(bind_flags) {
  sum(as.integer(bind_flags) * (2L ^ (seq_along(bind_flags) - 1L)))
}

#' Decode an integer regime index to a logical binding vector
#'
#' @param regime_idx Integer in [0, 2^k - 1]
#' @param n_specs    Number of OBC constraints (k)
#' @return Logical vector (length k) — TRUE where constraint j binds
#' @noRd
obc_regime_flags <- function(regime_idx, n_specs) {
  as.logical(bitwAnd(regime_idx, 2L ^ (seq_len(n_specs) - 1L)))
}


# =============================================================================
# Lazy policy cache
# =============================================================================

#' Lazily compute and cache the policy matrices for a given regime index
#'
#' regime_cache is an R environment used as a mutable hash map keyed by the
#' character form of regime_idx.  Each entry stores policy matrices and
#' pre-sliced system arrays so the Kalman filter can look them up in O(1).
#'
#' Regime 0 (all slack) is the dr_slack policy with zero constant offsets.
#' Non-zero regimes activate the binding subset of specs via obc_solve_binding.
#'
#' Cache entry fields: dr, c_state, c_obs, c_full, TT, RR, ZZ, DD.
#' c_full is the full n_endo constant (zero for slack, bound offset for binding).
#'
#' When copf_args is non-NULL (a list with Sigma_e, Sigma_e_inv, me_variance),
#' the entry additionally stores COPF-specific per-regime quantities:
#'   Omega     -- n_exo x n_exo posterior covariance of eps | y, s, r
#'   L_Omega   -- lower-triangular Cholesky of Omega (for drawing)
#'   F         -- n_obs x n_obs innovation covariance DD Sigma_e DD' + me*I
#'   F_inv     -- inverse of F
#'   log_det_F -- log|det F|
#'
#' @param regime_idx  Integer regime index (0 = all slack)
#' @param regime_cache R environment (new.env) accumulating computed policies
#' @param sys         System matrices (from extract_system_matrices_fast)
#' @param dr_slack    Slack-regime DecisionRules
#' @param specs       Full OBC spec list from obc_parse_tags
#' @param obs_idx     Integer vector of observable positions in endo vector
#' @param copf_args   NULL (bootstrap) or list(Sigma_e, Sigma_e_inv, me_variance)
#'                    When non-NULL, COPF matrices are computed and stored.
#' @return Named list with $dr, $c_state, $c_obs, $c_full, $TT, $RR, $ZZ, $DD
#'         (plus $Omega, $L_Omega, $F, $F_inv, $log_det_F when copf_args != NULL);
#'         NULL if obc_solve_binding fails (singular binding system)
#' @noRd
obc_ensure_policy <- function(regime_idx, regime_cache, sys, dr_slack, specs, obs_idx,
                               copf_args = NULL) {
  ## The system the cache was built for: boehl_simulate() and the other
  ## deterministic-path solvers derive the TIME-VARYING rules of a regime
  ## path from it (.obc_pwl_rules(); a per-regime entry below is the rule of
  ## a binding period whose NEXT period is slack, i.e. the last period of a
  ## spell only).
  if (!exists(".pwl_src", envir = regime_cache, inherits = FALSE))
    assign(".pwl_src", list(sys = sys, dr_slack = dr_slack, specs = specs),
           envir = regime_cache)
  key <- as.character(regime_idx)
  if (exists(key, envir = regime_cache, inherits = FALSE))
    return(get(key, envir = regime_cache, inherits = FALSE))

  si      <- dr_slack$state_idx
  n_endo  <- nrow(dr_slack$ghx)

  if (regime_idx == 0L) {
    entry <- list(
      dr      = dr_slack,
      c_state = numeric(length(si)),
      c_obs   = numeric(length(obs_idx)),
      c_full  = numeric(n_endo),              # zero constant for slack regime
      TT      = dr_slack$ghx[si,      , drop = FALSE],
      RR      = dr_slack$ghu[si,      , drop = FALSE],
      ZZ      = dr_slack$ghx[obs_idx, , drop = FALSE],
      DD      = dr_slack$ghu[obs_idx, , drop = FALSE]
    )
  } else {
    flags        <- obc_regime_flags(regime_idx, length(specs))
    active_specs <- specs[flags]
    result       <- obc_solve_binding(sys, dr_slack, active_specs, obs_idx)
    if (is.null(result)) return(NULL)
    dr_b <- result$dr
    si_b <- dr_b$state_idx
    entry <- list(
      dr      = dr_b,
      c_state = result$c_state,
      c_obs   = result$c_obs,
      c_full  = result$c_full,                # full-endo constant (n_endo)
      TT      = dr_b$ghx[si_b,    , drop = FALSE],
      RR      = dr_b$ghu[si_b,    , drop = FALSE],
      ZZ      = dr_b$ghx[obs_idx, , drop = FALSE],
      DD      = dr_b$ghu[obs_idx, , drop = FALSE]
    )
  }

  ## -- COPF per-regime precomputation (only when copf_args supplied) ----------
  if (!is.null(copf_args)) {
    Sigma_e     <- copf_args$Sigma_e
    Sigma_e_inv <- copf_args$Sigma_e_inv
    me_variance <- copf_args$me_variance
    DD_r        <- entry$DD
    n_obs_r     <- nrow(DD_r)

    ## Omega_r = (DD_r' DD_r / me + Sigma_e^{-1})^{-1}
    Omega_r_inv <- crossprod(DD_r) / me_variance + Sigma_e_inv
    ch_Oi  <- tryCatch(chol(Omega_r_inv), error = function(e2) NULL)
    if (!is.null(ch_Oi)) {
      Omega_r   <- chol2inv(ch_Oi)
      ch_O      <- tryCatch(chol(Omega_r), error = function(e2) NULL)
      L_Omega_r <- if (!is.null(ch_O)) t(ch_O) else NULL
    } else {
      Omega_r   <- NULL
      L_Omega_r <- NULL
    }

    ## F_r = DD_r Sigma_e DD_r' + me * I
    F_r      <- DD_r %*% Sigma_e %*% t(DD_r) + me_variance * diag(n_obs_r)
    ch_F     <- tryCatch(chol(F_r), error = function(e2) NULL)
    if (!is.null(ch_F)) {
      F_r_inv     <- chol2inv(ch_F)
      log_det_F_r <- 2 * sum(log(diag(ch_F)))
    } else {
      F_r_inv     <- NULL
      log_det_F_r <- NA_real_
    }

    entry$Omega     <- Omega_r
    entry$L_Omega   <- L_Omega_r
    entry$F         <- F_r
    entry$F_inv     <- F_r_inv
    entry$log_det_F <- log_det_F_r
  }

  assign(key, entry, envir = regime_cache)
  entry
}


# =============================================================================
# OccBin guess-and-verify
# =============================================================================

#' OccBin regime path of the data (piecewise-linear Kalman filter)
#'
#' Runs \code{kalman_filter_obc_pkf()} -- Dynare's OccBin PKF: in every
#' period the regime sequence expected at t is solved (guess-and-verify on
#' the constrained path from the one-step smoothed state and the filtered
#' shock) and the period's time-varying rule is used -- and returns its
#' regime path together with the regime cache holding the per-period rules,
#' so a following \code{kalman_filter_obc(..., regime_path)} evaluates the
#' same likelihood.
#'
#' Changed in 0.9.3.93 (W49): this was an outer loop that filtered with one
#' policy per regime (next period slack), re-checked every period against
#' the SLACK rule applied to the lagged filtered state, and locked periods in
#' which an observed constrained variable sat at its bound (an
#' observation pre-pass Dynare does not have).  Both checks are wrong for
#' spells of two or more periods.
#'
#' @param Y           Observation matrix (n_obs x T)
#' @param dr_slack    Slack-regime DecisionRules
#' @param sys         System matrices (from extract_system_matrices_fast)
#' @param obs_idx     Integer vector of observable positions in endo vector
#' @param model       dynhr_mod
#' @param params      Named numeric parameter vector
#' @param obs_vars    Character vector of observed variable names
#' @param specs       OBC spec list
#' @param me_variance Measurement error variance
#' @param max_iter    Maximum re-updates per period (Dynare's
#'                    likelihood.max_number_of_iterations, default 10)
#' @return List with:
#'   $regime_path     -- integer vector (length T): regime in period t
#'   $regime_expected -- list: regime sequence expected in period t
#'   $regime_cache    -- R environment: slack policy plus the PKF's
#'                       per-period rules
#'   $loglik          -- the PKF log-likelihood (-Inf on failure)
#'   $failed_period   -- NA or the period in which no regime converged
#' @noRd
obc_guess_verify <- function(Y, dr_slack, sys, obs_idx,
                              model, params, obs_vars, specs,
                              me_variance = 1e-8, max_iter = 10L) {
  regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)
  kf <- kalman_filter_obc_pkf(
    Y, dr_slack, regime_cache, sys, model, params, obs_vars, specs,
    obs_idx = obs_idx, me_variance = me_variance, max_inner = max_iter)
  list(regime_path     = kf$regime_path,
       regime_expected = kf$regime_expected,
       regime_cache    = regime_cache,
       loglik          = kf$loglik,
       failed_period   = kf$failed_period)
}
