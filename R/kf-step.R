## R/kf-step.R
## --------------------------------------------------------------------------
## Warm-startable single-period Kalman filter step.
##
## `kalman_filter()` runs the full recursion in one call; it exposes no entry
## point for advancing ONE period from a supplied (s, P) with a per-period
## shock scale.  The Rao-Blackwellized particle filter for stochastic
## volatility (see R/sv-rbpf.R) needs exactly that: for each volatility
## particle it advances the analytic Kalman recursion one step with
## `shock_scale = exp(h_t / 2)`, then resamples.
##
## `kf_step()` is that primitive.  It is the prediction-form recursion lifted
## verbatim from the `has_shock_scale` branch of the "dare"/"standard" paths
## in R/kalman-filter.R (lines ~1207-1262, ~1600-1640):
##
##   Se_t  = diag(scale) Sigma_e diag(scale)
##   v     = y - ZZ s - d
##   F     = ZZ P ZZ' + DD Se_t DD' + me_diag_t
##   ll    = -0.5 n_obs log(2pi) - 0.5 (log|F| + v' F^-1 v)
##   K     = (TT P ZZ' + RR Se_t DD') F^-1
##   s'    = TT s + K v
##   P'    = (TT-K ZZ) P (TT-K ZZ)' + (RR-K DD) Se_t (RR-K DD)'   [+ Joseph me_extra]
##
## The state (s, P) is in the PREDICTION representation: on entry s = E[x_t |
## y_{1:t-1}], P = Var; on exit s' = E[x_{t+1} | y_{1:t}].  Iterating from
## `kf_stationary_init()` reproduces `kalman_filter(..., method = "dare",
## lik_init = "stationary", shock_scale = S, return_ll_contrib = TRUE)$loglik`
## to machine precision (verified in test-kf-step.R).
## --------------------------------------------------------------------------


#' Stationary Kalman initialisation P0
#'
#' Returns the stationary prediction covariance
#' \eqn{P_0 = TT\,P_0\,TT' + RR\,\Sigma_e\,RR'} (the discrete Lyapunov solution),
#' the same \code{P0} that \code{\link{kalman_filter}(lik_init = "stationary")}
#' uses.  Baseline (unscaled) \code{Sigma_e} is used, matching the shipped
#' filter's initialisation even when a per-period \code{shock_scale} is applied
#' from period 1 onward.
#'
#' @param TT      State transition matrix (\code{dr$ghx[state_idx, ]}), n_s x n_s.
#' @param RR      Shock impact matrix (\code{dr$ghu[state_idx, ]}), n_s x n_e.
#' @param Sigma_e Baseline shock covariance, n_e x n_e.
#' @return The n_s x n_s stationary prediction covariance.
#' @seealso \code{\link{kf_step}}
#' @export
kf_stationary_init <- function(TT, RR, Sigma_e) {
  QQ <- tcrossprod(RR %*% Sigma_e, RR)
  P0 <- solve_lyapunov(TT, QQ)
  (P0 + t(P0)) * 0.5
}


#' One warm-started Kalman filter step
#'
#' Advances the prediction-form Kalman recursion one period, optionally with a
#' per-period multiplicative shock-standard-deviation \code{scale} (the SV /
#' heteroskedastic input) and per-period extra measurement-error variances.
#' This is the exact single-period map used inside the SV Rao-Blackwellized
#' particle filter (\code{\link{make_log_posterior_sv_rbpf}}).
#'
#' @param s        Prediction state mean \eqn{E[x_t \mid y_{1:t-1}]}, length n_s.
#' @param P        Prediction state covariance, n_s x n_s.
#' @param y        Observation at period t, length n_obs.
#' @param TT       State transition, n_s x n_s.
#' @param ZZ       State-to-observation loading, n_obs x n_s.
#' @param RR       Shock-to-state impact, n_s x n_e.
#' @param DD       Shock-to-observation impact, n_obs x n_e.
#' @param Sigma_e  Baseline shock covariance, n_e x n_e.
#' @param scale    Optional length-n_e multiplicative factors on shock standard
#'   deviation for THIS period (\code{exp(h_t/2)} in the SV filter). \code{NULL}
#'   (default) = baseline (all ones).
#' @param d        Optional observation intercept (SS level of the observables),
#'   length n_obs. \code{NULL} = zero.
#' @param me_diag  Optional length-n_obs baseline measurement-error variances
#'   (added to the diagonal of F only; documented F-regulariser convention).
#'   \code{NULL} = none.
#' @param me_extra Optional length-n_obs per-period EXTRA measurement-error
#'   variances (added to F AND propagated through the Joseph term, matching the
#'   shipped filter's \code{me_extra} treatment). \code{NULL} = none.
#' @return \code{list(ll, s, P)} — the period log-likelihood increment and the
#'   next-period prediction mean/covariance — or \code{NULL} if the forecast
#'   covariance is not positive definite (Cholesky failure), so a caller can
#'   treat the draw as infeasible.
#' @seealso \code{\link{kf_stationary_init}}, \code{\link{make_log_posterior_sv_rbpf}}
#' @export
kf_step <- function(s, P, y, TT, ZZ, RR, DD, Sigma_e,
                    scale = NULL, d = NULL, me_diag = NULL, me_extra = NULL) {
  n_obs <- length(y)

  ## Per-period scaled shock covariance Se_t = diag(scale) Sigma_e diag(scale)
  ## and the derived observation/cross terms rebuilt from the SAME Se_t so a
  ## shock hitting both state and observation is scaled consistently.
  if (is.null(scale)) {
    Se_t <- Sigma_e
  } else {
    Se_t <- Sigma_e * outer(scale, scale)
  }
  HH_t <- tcrossprod(DD %*% Se_t, DD)
  SS_t <- RR %*% Se_t %*% t(DD)

  ## Innovation.
  v <- y - as.numeric(ZZ %*% s)
  if (!is.null(d)) v <- v - d

  ## Forecast covariance F with measurement-error diagonal.
  PZ <- P %*% t(ZZ)
  Ft <- ZZ %*% PZ + HH_t
  if (!is.null(me_diag))  Ft <- Ft + diag(me_diag,  nrow = n_obs)
  if (!is.null(me_extra)) Ft <- Ft + diag(me_extra, nrow = n_obs)
  Ft <- (Ft + t(Ft)) * 0.5

  Fc <- tryCatch(chol(Ft), error = function(e) NULL)
  if (is.null(Fc)) return(NULL)
  Fi  <- chol2inv(Fc)
  ldf <- 2 * sum(log(diag(Fc)))
  ll  <- -0.5 * n_obs * log(2 * pi) - 0.5 * (ldf + drop(crossprod(v, Fi %*% v)))
  if (!is.finite(ll)) return(NULL)

  ## Kalman gain and prediction-form update.
  K_t  <- (TT %*% PZ + SS_t) %*% Fi
  s_n  <- as.numeric(TT %*% s) + drop(K_t %*% v)
  TmKZ <- TT - K_t %*% ZZ
  RmKD <- RR - K_t %*% DD
  P_n  <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Se_t, RmKD)
  ## Joseph true-noise term for the me_extra part only (matches the shipped
  ## filter: base me_diag stays on the F-only-regulariser convention).
  if (!is.null(me_extra)) P_n <- P_n + K_t %*% (me_extra * t(K_t))
  P_n <- (P_n + t(P_n)) * 0.5

  list(ll = ll, s = s_n, P = P_n)
}
