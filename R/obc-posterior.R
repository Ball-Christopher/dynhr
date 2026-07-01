## R/obc-posterior.R
## --------------------------------------------------------------------------
## OBC log-posterior factory functions for MCMC estimation.
##
## Provides:
##   make_log_posterior_obc()      -- legacy: OccBin outer-loop + kalman_filter_obc
##   make_log_posterior_obc_pkf()  -- preferred: Pfeiffer-Ratto PKF (exported)
##
## Both return Function(theta) -> list(logpost, loglik, logprior[, regime_path])
## and are drop-in replacements for make_log_posterior() when the model has OBC
## constraints.
##
## make_log_posterior_obc_pkf() is the preferred estimator: the per-period inner
## convergence loop discovers regimes from the extracted shock eps_{t|t} and
## backward-smoothed state s_{t-1|t}, requiring no separate outer pass.
## Supports warm-starting via attr(theta, "regime_hint").
## --------------------------------------------------------------------------


#' Create an OBC-aware log-posterior evaluator (legacy OccBin outer loop)
#'
#' Mirrors make_log_posterior (dynhr_estimation.R) but adds two steps inside
#' the closure at each parameter draw:
#'   1. Run the OccBin guess-and-verify loop (lazily building per-regime policies).
#'   2. Pass the regime path and policy cache to kalman_filter_obc.
#'
#' Validation guards (obc_assert_linear, obc_parse_tags) run ONCE at factory
#' creation time, not inside the per-draw closure.
#'
#' @param model       dynhr_mod
#' @param data        Observation matrix (T x n_obs or n_obs x T)
#' @param prior_spec  Prior specification data.frame
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled
#' @param specs       OBC spec list (from obc_parse_tags), or NULL to parse
#' @param me_variance Measurement error variance (default 1e-8)
#' @return Function(theta) -> list(logpost, loglik, logprior)
#' @noRd
make_log_posterior_obc <- function(model, data, prior_spec, obs_vars,
                                    compiled, specs = NULL,
                                    me_variance = 1e-8) {

  # ---- Guards: run once, fail fast ----------------------------------------
  obc_assert_linear(model)
  if (is.null(specs)) specs <- obc_parse_tags(model)

  # ---- One-time setup (mirrors make_log_posterior) -------------------------
  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence)) {
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  }
  sys_cache <- cache_system_structure(compiled)

  # Observable indices: fixed by model structure, not by parameter draw
  endo    <- model$var_names
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx))) {
    stop("obs_vars contains names not found in model$var_names: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  }

  # Ensure data is n_obs x T (columns = time periods)
  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  # ---- Closure: evaluated at each MCMC draw --------------------------------
  function(theta) {
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    params <- .apply_theta_to_params(model, theta)

    ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
    if (is.null(ss_result) || !ss_result$converged)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)

    # Slack policy function
    dr_slack <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr_slack) || !dr_slack$bk_satisfied)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    max_eig <- max(Mod(eigen(
      dr_slack$ghx[dr_slack$state_idx, , drop = FALSE],
      only.values = TRUE
    )$values))
    if (max_eig >= 1)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    # OccBin regime path + lazy per-regime policy cache
    gv <- obc_guess_verify(
      Y, dr_slack, sys, obs_idx,
      model, params, obs_vars, specs,
      me_variance = me_variance
    )
    if (is.null(gv))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    # Regime-switching Kalman filter
    kf <- kalman_filter_obc(
      Y, dr_slack, gv$regime_cache,
      model, params, obs_vars, gv$regime_path,
      me_variance     = me_variance,
      return_filtered = FALSE
    )
    if (is.null(kf) || !is.finite(kf$loglik))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    list(logpost = kf$loglik + lp, loglik = kf$loglik, logprior = lp)
  }
}


#' Create a PKF-based log-posterior evaluator (Pfeiffer-Ratto inversion filter)
#'
#' Drop-in replacement for make_log_posterior_obc() that calls
#' kalman_filter_obc_pkf() directly.  Unlike the outer obc_guess_verify()
#' approach, the per-period inner loop discovers the OBC regime from the
#' extracted shock and backward-smoothed state, so no separate pre-pass is
#' needed.
#'
#' Warm-starting: the returned closure optionally accepts a $regime_hint
#' attribute on the theta vector (integer vector length T) to initialise the
#' inner loop from a previous MCMC draw's accepted regime path.
#'
#' @param model       dynhr_mod
#' @param data        Observation matrix (T x n_obs or n_obs x T)
#' @param prior_spec  Prior specification data.frame
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled
#' @param specs       OBC spec list (from obc_parse_tags), or NULL to parse
#' @param me_variance Measurement error variance (default 1e-8)
#' @param max_inner   Max inner iterations per period (default 10)
#' @return Function(theta) -> list(logpost, loglik, logprior, regime_path)
#' @export
make_log_posterior_obc_pkf <- function(model, data, prior_spec, obs_vars,
                                        compiled, specs = NULL,
                                        me_variance = 1e-8,
                                        max_inner   = 10L) {

  ## prior_spec/me_variance/max_inner are only referenced inside the returned
  ## closure, so without forcing they remain unevaluated promises pointing at
  ## the caller's frame. A mirai daemon that ships this closure before it has
  ## been called once would fail to resolve them ("object ... not found").
  force(prior_spec); force(me_variance); force(max_inner)

  if (is.null(specs)) specs <- obc_parse_tags(model)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence)) {
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  }
  sys_cache <- cache_system_structure(compiled)

  endo    <- model$var_names
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("obs_vars contains names not found in model$var_names: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  function(theta) {
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp,
                  regime_path = NULL))

    params <- .apply_theta_to_params(model, theta)

    ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
    if (is.null(ss_result) || !ss_result$converged)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp,
                  regime_path = NULL))

    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)

    dr_slack <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr_slack) || !dr_slack$bk_satisfied)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp,
                  regime_path = NULL))

    max_eig <- max(Mod(eigen(
      dr_slack$ghx[dr_slack$state_idx, , drop = FALSE],
      only.values = TRUE
    )$values))
    if (max_eig >= 1)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp,
                  regime_path = NULL))

    # Seed the regime cache with the slack policy
    regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
    obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)

    # Warm-start from a previously accepted regime path if provided
    regime_hint <- attr(theta, "regime_hint")

    kf <- kalman_filter_obc_pkf(
      Y, dr_slack, regime_cache, sys,
      model, params, obs_vars, specs,
      obs_idx          = obs_idx,
      regime_path_init = regime_hint,
      me_variance      = me_variance,
      max_inner        = max_inner,
      return_filtered  = FALSE,
      return_shocks    = FALSE
    )
    if (is.null(kf) || !is.finite(kf$loglik))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp,
                  regime_path = NULL))

    list(
      logpost     = kf$loglik + lp,
      loglik      = kf$loglik,
      logprior    = lp,
      regime_path = kf$regime_path
    )
  }
}
