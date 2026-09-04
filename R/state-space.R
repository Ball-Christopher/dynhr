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
