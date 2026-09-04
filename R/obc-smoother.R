## R/obc-smoother.R
## --------------------------------------------------------------------------
## OBC PKF fixed-interval smoother (Durbin-Koopman backward recursion).
##
## Provides:
##   pkf_smoother_obc() -- DK smoother using kf_store from kalman_filter_obc_pkf
##
## Uses per-period quantities {v, F_inv, L=TT-K*ZZ, P_in, s_in, RR, ZZ, DD}
## stored by kalman_filter_obc_pkf(..., return_store = TRUE).
##
## kf_store invariant: v[[t]] is NULL (not missing) for no-obs periods.
## The smoother uses is.null(v) to distinguish obs vs no-obs periods.
## --------------------------------------------------------------------------


#' Fixed-interval smoother for the OBC piecewise Kalman filter
#'
#' Implements the Durbin-Koopman (2003) backward recursion for the lag-1
#' observation state-space used by dynhr, using the per-period quantities
#' stored by kalman_filter_obc_pkf(..., return_store = TRUE).
#'
#' For each period t = T, T-1, ..., 1 the recursion is:
#' \preformatted{
#'   Step 0 (Kalman gain, recomputed from the stored quantities):
#'     K_t = (TT_t * P_{t-1|t-1} * ZZ_t' + RR_t * Sigma_e * DD_t') * F_t^{-1}
#'
#'   Step 1 (smoothed shock):
#'     eps_{t|T} = Sigma_e * (DD_t' * F_t^{-1} * v_t
#'                            + (RR_t - K_t * DD_t)' * r_t)
#'
#'   Step 2 (backward update of smoother residual):
#'     r_{t-1} = ZZ_t' * F_t^{-1} * v_t + L_t' * r_t
#'     where L_t = TT_t - K_t * ZZ_t  (= TmKZ stored by the filter)
#'
#'   Step 3 (smoothed lagged state):
#'     s_{t-1|T} = s_{t-1|t-1} + P_{t-1|t-1} * r_{t-1}
#' }
#' Initialisation: r_T = 0 (Koopman & Durbin 2003, p.89).
#'
#' \strong{Correlated state/observation noise.}  In dynhr's lag-1 observation
#' timing the SAME shock drives the transition and the observation
#' (\eqn{s_t = TT s_{t-1} + RR \epsilon_t}, \eqn{y_t = ZZ s_{t-1} + DD \epsilon_t}), so
#' \eqn{Cov(RR \epsilon_t, DD \epsilon_t) = RR \Sigma_e DD' \neq 0} and the Kalman gain
#' carries that cross-term.  The future-information loading in step 1 is
#' therefore \eqn{(RR_t - K_t DD_t)'}, not \eqn{RR_t'}: for \eqn{j > t} the innovation
#' \eqn{v_j} sees \eqn{\epsilon_t} only through the prediction error
#' \eqn{\alpha_{t+1} - a_{t+1}}, whose \eqn{\epsilon_t} loading is \eqn{RR_t - K_t DD_t}
#' because \eqn{a_{t+1}} has already absorbed \eqn{K_t v_t}.  Dropping the
#' \eqn{-K_t DD_t} term (the pre-0.9.2.0003 form) makes the smoothed shocks
#' inconsistent with the smoothed states: the transition residual
#' \eqn{s_{t|T} - (TT_t s_{t-1|T} + RR_t \epsilon_{t|T} + c_t)} was ~4e-4 on a
#' 1e-2-scale state, and is ~1e-18 with the cross-term restored.  This is the
#' same \eqn{eps_{t|T} = Q (R' r_t + D' u_t)}, \eqn{u_t = F^{-1} v_t - K_t' r_t}
#' form used by \code{\link{kalman_smoother}} (Durbin & Koopman 2012 §4.5).
#' \eqn{K_t} is not stored by the filter, so it is recomputed in step 0 from
#' the stored \eqn{TT, P_{in}, ZZ, RR, DD, F^{-1}} -- exactly the expression the
#' filter used.
#'
#' When observations are missing at period t (v_t = NULL), the step reduces
#' to pure backward propagation: \eqn{r_{t-1} = TT_t' \cdot r_t} (\eqn{L_t = TT_t} in that
#' case), and the smoothed shock uses only the \eqn{RR_t' \cdot r_t} future-info term
#' (there is no gain to correct for: \eqn{K_t = 0} when nothing is observed).
#'
#' Note on indexing: smoothed_states[, t] is \eqn{s_{t-1|T}}, the smoothed state
#' that appeared *lagged* in period t's observation equation.  This differs
#' from filtered_states[, t] = \eqn{s_{t|t}} returned by kalman_filter_obc_pkf
#' (the current-period updated state).
#'
#' @param kf_store  List returned by kalman_filter_obc_pkf(..., return_store = TRUE)
#' @return List with:
#'   \item{smoothed_states}{n_state x T: \eqn{s_{t-1|T}} (lagged state smoothed over all T obs)}
#'   \item{smoothed_shocks}{n_exo x T: \eqn{\epsilon_{t|T}} (shocks smoothed over all T obs)}
#' @export
pkf_smoother_obc <- function(kf_store) {
  n_T     <- kf_store$n_T
  n_state <- kf_store$n_state
  n_exo   <- kf_store$n_exo
  Sigma_e <- kf_store$Sigma_e

  smoothed_states <- matrix(0, n_state, n_T)
  smoothed_shocks <- matrix(0, n_exo,   n_T)

  r <- numeric(n_state)   # r_T = 0 (Koopman-Durbin initialisation)

  for (t in rev(seq_len(n_T))) {
    v     <- kf_store$v[[t]]
    L     <- kf_store$L[[t]]     # TT_t - K_t*ZZ_t  (or  TT_t  if no obs)
    P_in  <- kf_store$P_in[[t]]  # P_{t-1|t-1}
    s_in  <- kf_store$s_in[, t, drop = TRUE]   # s_{t-1|t-1}
    RR    <- kf_store$RR[[t]]    # n_state x n_exo

    if (is.null(v)) {
      # No observations at period t:
      #   eps_{t|T} = Sigma_e * RR_t' * r_t  (future-info only, no obs term)
      #   r_{t-1}   = L_t' * r_t  = TT_t' * r_t
      smoothed_shocks[, t] <- drop(Sigma_e %*% (t(RR) %*% r))
      r                    <- drop(t(L) %*% r)
      smoothed_states[, t] <- drop(s_in + P_in %*% r)
    } else {
      F_inv <- kf_store$F_inv[[t]]
      ZZ    <- kf_store$ZZ[[t]]   # n_obs_t x n_state (possibly subsetted)
      DD    <- kf_store$DD[[t]]   # n_obs_t x n_exo
      TT    <- kf_store$TT[[t]]   # n_state x n_state

      Finv_v <- F_inv %*% v       # precompute once

      # Step 0: rebuild the filter's Kalman gain
      #   K_t = (TT P_{t-1|t-1} ZZ' + RR Sigma_e DD') F^{-1}
      # (the filter does not store K; this is the same expression it used,
      #  with the same P_in — P is only updated after K is formed).
      K    <- (TT %*% P_in %*% t(ZZ) + RR %*% Sigma_e %*% t(DD)) %*% F_inv
      RmKD <- RR - K %*% DD       # eps_t loading of the t+1 prediction error

      # Step 1: smoothed shock (uses r_t BEFORE the backward update).
      # The future-information term loads RR - K*DD, not RR: see the
      # correlated-noise note in the roxygen block above.
      smoothed_shocks[, t] <- drop(Sigma_e %*% (t(DD) %*% Finv_v + t(RmKD) %*% r))

      # Step 2: backward update of smoother residual
      r <- drop(t(ZZ) %*% Finv_v) + drop(t(L) %*% r)

      # Step 3: smoothed lagged state
      smoothed_states[, t] <- drop(s_in + P_in %*% r)
    }
  }

  if (!is.null(kf_store$endo_state_names))
    rownames(smoothed_states) <- kf_store$endo_state_names
  if (!is.null(kf_store$exo_names))
    rownames(smoothed_shocks) <- kf_store$exo_names

  list(smoothed_states = smoothed_states, smoothed_shocks = smoothed_shocks)
}
