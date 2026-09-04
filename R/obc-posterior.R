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
#' @param me_variance Measurement error variance (default 1e-8). Guarded
#'   (warn-only, once per closure) against a near-degenerate slack-regime
#'   innovation covariance via \code{getOption("dynhr.me_floor_check", TRUE)};
#'   see \code{.obc_warn_me_floor_lock()} in R/obc-regime.R.
#' @param power       Power-posterior (generalised-Bayes) tempering exponent
#'   \eqn{\zeta}: \code{$logpost} becomes
#'   \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)} while \code{$loglik}
#'   keeps the RAW (untempered) value. \code{NULL} (default) resolves the
#'   \code{power_posterior} package option ONCE, at factory time.
#' @return Function(theta) -> list(logpost, loglik, logprior)
#' @noRd
make_log_posterior_obc <- function(model, data, prior_spec, obs_vars,
                                    compiled, specs = NULL,
                                    me_variance = 1e-8,
                                    power = NULL) {

  ## Resolve zeta ONCE here, not per draw (see .resolve_power_posterior).
  power <- .resolve_power_posterior(power, "make_log_posterior_obc")

  # ---- Guards: run once, fail fast ----------------------------------------
  obc_assert_linear(model)
  if (is.null(specs)) specs <- obc_parse_tags(model)

  # ---- One-time setup (mirrors make_log_posterior) -------------------------
  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence)) {
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  }
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

  ## Adapter over the shared closure builder (R/posterior-closure.R). Specific
  ## to this branch: the cold (never warm-started) steady-state solve, the
  ## always-eigen() stationarity guard, and the OccBin guess-and-verify pass
  ## feeding kalman_filter_obc(). No system priors on the OBC paths.
  .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    stationarity = "eigen",
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      dr_slack <- sol$dr
      .obc_warn_me_floor_lock(
        dr_slack, model, params, obs_vars, obs_idx, me_variance,
        check = me_floor_check)

      # OccBin regime path + lazy per-regime policy cache
      gv <- obc_guess_verify(
        Y, dr_slack, sol$sys, obs_idx,
        model, params, obs_vars, specs,
        me_variance = me_variance
      )
      if (is.null(gv)) return(NULL)

      # Regime-switching Kalman filter
      kf <- kalman_filter_obc(
        Y, dr_slack, gv$regime_cache,
        model, params, obs_vars, gv$regime_path,
        me_variance     = me_variance,
        return_filtered = FALSE
      )
      if (is.null(kf) || !is.finite(kf$loglik)) return(NULL)
      list(loglik = kf$loglik)
    },
    power      = power,
    warm_start = FALSE)
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
#' @param me_variance Measurement error variance (default 1e-8). Guarded
#'   (warn-only, once per closure) against a near-degenerate slack-regime
#'   innovation covariance via \code{getOption("dynhr.me_floor_check", TRUE)};
#'   see \code{.obc_warn_me_floor_lock()} in R/obc-regime.R.
#' @param max_inner   Max inner iterations per period (default 10)
#' @param power       Power-posterior (generalised-Bayes) tempering exponent
#'   \eqn{\zeta}: \code{$logpost} becomes
#'   \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)} while \code{$loglik}
#'   keeps the RAW (untempered) value. \code{NULL} (default) resolves the
#'   \code{power_posterior} package option ONCE, at factory time.
#' @return Function(theta) -> list(logpost, loglik, logprior, regime_path)
#' @export
make_log_posterior_obc_pkf <- function(model, data, prior_spec, obs_vars,
                                        compiled, specs = NULL,
                                        me_variance = 1e-8,
                                        max_inner   = 10L,
                                        power       = NULL) {

  ## prior_spec/me_variance/max_inner are only referenced inside the returned
  ## closure, so without forcing they remain unevaluated promises pointing at
  ## the caller's frame. A mirai daemon that ships this closure before it has
  ## been called once would fail to resolve them ("object ... not found").
  force(prior_spec); force(me_variance); force(max_inner)
  ## Resolve zeta ONCE here, not per draw (see .resolve_power_posterior).
  power <- .resolve_power_posterior(power, "make_log_posterior_obc_pkf")

  if (is.null(specs)) specs <- obc_parse_tags(model)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence)) {
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  }
  endo    <- model$var_names
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("obs_vars contains names not found in model$var_names: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  ## Adapter over the shared closure builder (R/posterior-closure.R). Two
  ## things are specific here beyond the PKF call itself: the cold
  ## steady-state solve with the always-eigen() stationarity guard, and the
  ## RETURN SHAPE -- this is the only factory that carries a `regime_path`
  ## field, NULL on every rejected draw, which `reject_fields` supplies.
  .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    stationarity = "eigen",
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      dr_slack <- sol$dr
      .obc_warn_me_floor_lock(
        dr_slack, model, params, obs_vars, obs_idx, me_variance,
        check = me_floor_check)

      # Seed the regime cache with the slack policy
      regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
      obc_ensure_policy(0L, regime_cache, sol$sys, dr_slack, specs, obs_idx)

      # Warm-start from a previously accepted regime path if provided
      regime_hint <- attr(theta, "regime_hint")

      kf <- kalman_filter_obc_pkf(
        Y, dr_slack, regime_cache, sol$sys,
        model, params, obs_vars, specs,
        obs_idx          = obs_idx,
        regime_path_init = regime_hint,
        me_variance      = me_variance,
        max_inner        = max_inner,
        return_filtered  = FALSE,
        return_shocks    = FALSE
      )
      if (is.null(kf) || !is.finite(kf$loglik)) return(NULL)
      list(loglik = kf$loglik, extra = list(regime_path = kf$regime_path))
    },
    power         = power,
    warm_start    = FALSE,
    reject_fields = list(regime_path = NULL))
}
