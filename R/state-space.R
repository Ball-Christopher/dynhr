## R/state-space.R
## --------------------------------------------------------------------------
## Tagged state-space object and timing-convention converter.
##
## dynhr canonical state-space convention ("lagged-state", Convention A):
##   s_t = T s_{t-1} + R eps_t,   eps_t ~ N(0, Sigma_e)
##   y_t = Z s_{t-1} + D eps_t + d
##
## Alternative "current-state" convention (Convention B):
##   s_t = T s_{t-1} + R eps_t
##   y_t = C s_t     + D eps_t
##
## Algebraic conversion B -> A (ss_convert_timing):
##   ZZ_new = C T
##   DD_new = C R + D
##   RR_new = R  (unchanged)
##
## When D != 0, the converted DD_new captures the mixed shock term; the
## cross-covariance SS = RR_new Sigma_e DD_new' is non-zero even when the
## original B system has D = 0.  The downstream Kalman filter handles this
## correctly via its SS term.
##
## Consumers of the lagged-state object:
##   kalman_filter, kalman_smoother, gradient-tangent-kf, gradient-adjoint-kf,
##   whittle-likelihood, tpf-likelihood, conditional_forecast
## --------------------------------------------------------------------------


#' Constructor for a tagged DSGE state-space object
#'
#' Creates a \code{dsge_ss} object encoding either the lagged-state or
#' current-state timing convention.  The lagged-state form is dynhr's
#' canonical convention and is the form every filter, smoother and
#' frequency-domain likelihood in the package works in.
#'
#' @param T_mat  n_state x n_state state-transition matrix.
#' @param R_mat  n_state x n_shock shock-impact matrix (unit-shock responses;
#'   \code{Sigma_e} is applied by the filter, not baked in).
#' @param Z_mat  n_obs x n_state observation-loading matrix.  Under the
#'   lagged convention this multiplies \code{s_{t-1}}; under the current
#'   convention it multiplies \code{s_t}.
#' @param D_mat  n_obs x n_shock direct-shock matrix in the observation
#'   equation.
#' @param Sigma_e n_shock x n_shock shock covariance matrix.
#' @param d      n_obs-length constant term in the observation equation
#'   (default zero vector) -- the observation intercept, subtracted from the
#'   data by everything that consumes the object, which is what makes the
#'   data convention LEVELS. \code{\link{build_dsge_state_space}()} fills it
#'   with \code{dr$ys[obs_vars]}; a hand-built state space has no steady
#'   state behind it, so its default zero means its data is in deviations by
#'   construction.
#' @param state_names Character vector of state-variable names (length n_state).
#' @param obs_vars   Character vector of observable names (length n_obs).
#' @param shock_names Character vector of shock names (length n_shock).
#' @param timing  \code{"lagged"} (default, dynhr canonical) or
#'   \code{"current"}.  See \emph{Details}.
#' @param ...    Additional fields stored verbatim in the returned list.
#'
#' @details
#' \strong{Lagged-state (timing = "lagged"):}
#' \deqn{s_t = T s_{t-1} + R \epsilon_t}
#' \deqn{y_t = Z s_{t-1} + D \epsilon_t + d}
#'
#' \strong{Current-state (timing = "current"):}
#' \deqn{s_t = T s_{t-1} + R \epsilon_t}
#' \deqn{y_t = Z s_t + D \epsilon_t + d}
#'
#' The naming convention in the object always uses \code{Z_mat} and
#' \code{D_mat} regardless of the timing field; the \code{timing} tag
#' records which convention the matrices encode.
#'
#' Use \code{\link{ss_convert_timing}} to convert a current-state object to
#' the lagged-state form required by the Kalman filter and smoother.
#'
#' @return A list with class \code{"dsge_ss"} containing all input matrices
#'   plus derived fields \code{n_state}, \code{n_obs}, \code{n_shock}.
#'
#' @seealso \code{\link{ss_convert_timing}}, \code{\link{build_dsge_state_space}}
#' @export
new_dsge_ss <- function(T_mat, R_mat, Z_mat, D_mat, Sigma_e,
                        d = NULL,
                        state_names = NULL,
                        obs_vars   = NULL,
                        shock_names = NULL,
                        timing      = c("lagged", "current"),
                        ...) {
  timing <- match.arg(timing)

  n_state <- nrow(T_mat)
  n_obs   <- nrow(Z_mat)
  n_shock <- ncol(R_mat)

  if (is.null(d)) d <- rep(0, n_obs)

  structure(
    c(
      list(
        T_mat       = T_mat,
        R_mat       = R_mat,
        Z_mat       = Z_mat,
        D_mat       = D_mat,
        Sigma_e     = Sigma_e,
        d           = d,
        state_names = state_names,
        obs_names   = obs_vars,
        shock_names = shock_names,
        n_state     = n_state,
        n_obs       = n_obs,
        n_shock     = n_shock,
        timing      = timing
      ),
      list(...)
    ),
    class = "dsge_ss"
  )
}


## ---- Deterministic observation trends (Dynare `observation_trends`) -------
##
## Dynare's measurement equation with an observation_trends block is
##     y_t = ys + trend_coeff * (first_obs + t - 1) + Z s_{t-1} + D eps_t,
## t = 1..T over the estimation sample (dsge_likelihood.m builds
## `trend = constant + trend_coeff * (first_obs:first_obs+T-1)` and filters
## Y - trend; compute_trend_coefficients.m evaluates each slope expression at
## the current parameters). The trend is therefore a time-varying part of the
## observation INTERCEPT: every filter/smoother here subtracts it from the data
## exactly where it subtracts d = ys[obs], and nothing else changes.

#' Observation-trend slopes at a parameter vector
#'
#' @param model  Parsed model; its `observation_trends` field (from
#'   parse_mod) holds `trends` (observable -> slope expression text) and
#'   `first_obs`.
#' @param params Named parameter vector the slope expressions are evaluated
#'   at (`NULL` = the model's calibration).
#' @param obs_vars Observables, in data-column order.
#' @return Named numeric slope per observable (0 where none), or NULL when no
#'   observable in `obs_vars` carries a trend.
#' @noRd
.has_obs_trends <- function(model, dr = NULL) {
  ot <- model$observation_trends
  (!is.null(ot) && length(ot$trends) > 0L) ||
    ## A balanced-growth model's observable growth is the same kind of
    ## deterministic intercept trend (R/balanced-growth.R).
    isTRUE(model$balanced_growth) || !is.null(.dr_growth(dr))
}

#' Abort unless the caller's filter honours observation_trends.
#'
#' Only kalman_filter() / the Gaussian posterior and the Kalman smoother
#' subtract the deterministic trend. Every other likelihood, gradient and
#' diagnostic builds its own measurement intercept from `dr$ys`, so on a
#' trended model it would silently score untrended data.
#' @noRd
.refuse_obs_trends <- function(model, what, dr = NULL) {
  if (.has_obs_trends(model, dr)) {
    growth <- isTRUE(model$balanced_growth) || !is.null(.dr_growth(dr))
    .dynhr_abort(
      what, " does not support ",
      if (growth) paste0("a balanced-growth model (solved with ",
                         "steady_options = list(growth = TRUE); its growth ",
                         "is a deterministic observation trend)")
      else "observation_trends",
      ": only the Gaussian Kalman ",
      "filter (kalman_filter(), make_log_posterior(likelihood = \"gaussian\")) ",
      "and the Kalman smoother subtract the deterministic trend. Use the ",
      "Gaussian likelihood, or detrend the data and drop the ",
      if (growth) "drift." else "block.",
      class = "dynhr_error_observation_trends_unsupported")
  }
  invisible(TRUE)
}

## `dr` (optional): a decision rule whose `$growth` (balanced growth, see
## R/balanced-growth.R) adds g[obs] to the slope of each observable -- the
## observation intercept of period t is then ys[obs] + g[obs] * (first_obs +
## t - 1), ADDED to any observation_trends slope. The model and dr must agree
## about balanced growth (.check_growth_pair).
.obs_trend_slopes <- function(model, params, obs_vars, dr = NULL) {
  .check_growth_pair(model, dr)
  g <- .dr_growth(dr)
  g_obs <- NULL
  if (!is.null(g)) {
    g_obs <- stats::setNames(numeric(length(obs_vars)), obs_vars)
    hit_g <- intersect(obs_vars, names(g))
    g_obs[hit_g] <- g[hit_g]
    if (!any(g_obs != 0)) g_obs <- NULL
  }
  ot <- model$observation_trends
  hit <- if (is.null(ot) || length(ot$trends) == 0L) character(0)
         else intersect(obs_vars, names(ot$trends))
  if (length(hit) == 0L) return(g_obs)
  if (is.null(params)) params <- model$param_values
  env <- .dynhr_param_eval_env(params)
  slopes <- stats::setNames(numeric(length(obs_vars)), obs_vars)
  if (!is.null(g_obs)) slopes <- slopes + g_obs
  for (nm in hit) {
    v <- .dynhr_sandbox_eval(ot$trends[[nm]], env,
                             context = "the observation_trends expression")
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      .dynhr_abort(
        "observation_trends: the trend slope of ", nm, " (`", ot$trends[[nm]],
        "`) does not evaluate to a finite number at the supplied parameters",
        if (is.null(v)) " (a parameter it uses has no value)" else "", ".",
        class = "dynhr_error_observation_trends")
    slopes[[nm]] <- slopes[[nm]] + v
  }
  slopes
}

#' Observation-trend path, n_obs x n_T (column t is period first_obs + t - 1)
#'
#' @inheritParams .obs_trend_slopes
#' @param n_T Number of sample periods (data rows).
#' @return Numeric matrix, or NULL when there is no trend to subtract.
#' @noRd
.obs_trend_path <- function(model, params, obs_vars, n_T, dr = NULL) {
  slopes <- .obs_trend_slopes(model, params, obs_vars, dr = dr)
  if (is.null(slopes)) return(NULL)
  first_obs <- model$observation_trends$first_obs %||% 1L
  outer(slopes, first_obs - 1 + seq_len(n_T))
}


#' Convert a state-space object from current-state to lagged-state timing
#'
#' Converts a \code{dsge_ss} object tagged \code{timing = "current"} to
#' the lagged-state convention every filter, smoother and frequency-domain
#' likelihood in the package requires.  The observation intercept \code{d} is
#' timing-invariant and carries over unchanged.  If
#' \code{ss$timing} is already \code{"lagged"} the object is returned
#' unchanged.
#'
#' @details
#' The algebraic conversion from current-state to lagged-state form is:
#' \deqn{ZZ_{\text{new}} = C T}
#' \deqn{DD_{\text{new}} = C R + D}
#' \deqn{RR_{\text{new}} = R  \quad \text{(unchanged)}}
#'
#' When \eqn{D \ne 0}, the converted \eqn{DD_{\text{new}}} captures the
#' mixed observation-noise term \eqn{C R \epsilon_t + D \epsilon_t}.  The
#' Kalman filter's cross-covariance \eqn{SS = R \Sigma_e DD_{\text{new}}'}
#' then accounts for the correlation between state and observation noise
#' that arises from the shared \eqn{\epsilon_t}.
#'
#' @param ss A \code{dsge_ss} object (from \code{new_dsge_ss} or
#'   \code{build_dsge_state_space}).
#' @return A \code{dsge_ss} object with \code{timing = "lagged"}.
#'
#' @seealso \code{\link{new_dsge_ss}}, \code{\link{kalman_smoother}}
#' @export
ss_convert_timing <- function(ss) {
  if (!inherits(ss, "dsge_ss"))
    stop("ss_convert_timing: argument must be a dsge_ss object.", call. = FALSE)
  if (ss$timing == "lagged") return(ss)

  T_mat <- ss$T_mat
  R_mat <- ss$R_mat
  C_mat <- ss$Z_mat   # current-state naming: Z = C
  D_mat <- ss$D_mat

  ## Algebraic conversion: ZZ_new = C T,  DD_new = C R + D
  ZZ_new <- C_mat %*% T_mat
  DD_new <- C_mat %*% R_mat + D_mat

  new_dsge_ss(
    T_mat       = T_mat,
    R_mat       = R_mat,
    Z_mat       = ZZ_new,
    D_mat       = DD_new,
    Sigma_e     = ss$Sigma_e,
    d           = ss$d,
    state_names = ss$state_names,
    obs_vars   = ss$obs_names,
    shock_names = ss$shock_names,
    timing      = "lagged"
  )
}
