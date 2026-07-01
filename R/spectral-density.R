## R/spectral-density.R
## --------------------------------------------------------------------------
## Spectral density module for dynhr DSGE models.
##
## Public entry point: spectral_density(ss, omega, me_variance = 0)
## Internal kernel:    .spectral_density_core(omega, TT, RR, ZZ, DD, Sigma_e,
##                                             me_variance = 0)
##
## The four prior implementations in whittle-likelihood.R, diag-pre-d23-spectral.R
## (two copies), and diag-pre-d24-global-kl.R are replaced by one-line shims
## that delegate here.  No numeric change: the kernel body is identical to the
## former .whittle_spectral_density.
##
## Lagged-state (dynhr canonical) transfer function:
##   H(e^{iω}) = ZZ * z * (I − TT * z)^{-1} * RR + DD,   z = e^{-iω}
##   S_yy(ω)   = H * Sigma_e * H^* + me_variance * I
##
## Current-state transfer function (used by D24):
##   H(e^{iω}) = obs_mat * (I − T * z)^{-1} * R
##   S_yy(ω)   = H * Sigma_e * H^*
## This is handled by .spectral_density_core_current_no_d (no z factor, no D).
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## Internal kernel — lagged-state convention
## --------------------------------------------------------------------------

## Transfer function for the lagged-state (dynhr canonical) convention.
##   H = ZZ * z * (I - TT z)^{-1} * RR + DD
## Includes a near-unit-root geometric-series fallback.
##
## @param omega     Single frequency in (0, pi].
## @param TT        n_state x n_state state-transition matrix.
## @param RR        n_state x n_shock shock-impact matrix.
## @param ZZ        n_obs   x n_state observation matrix (lagged convention).
## @param DD        n_obs   x n_shock direct-impact matrix.
## @param Sigma_e   n_shock x n_shock shock covariance.
## @param me_variance Scalar >= 0; added to diagonal of S.
## @return n_obs x n_obs complex Hermitian spectral density matrix.
## @noRd
.spectral_density_core <- function(omega, TT, RR, ZZ, DD, Sigma_e,
                                    me_variance = 0) {
  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)

  z     <- exp(-1i * omega)           # e^{-i*omega}
  ## (I - TT z)
  A     <- diag(n_state) - TT * z
  rcond_A <- rcond(A)
  if (is.finite(rcond_A) && rcond_A > .Machine$double.eps * 1e4) {
    ## H = ZZ * z * (I - TT z)^{-1} * RR + DD
    H <- ZZ %*% (z * solve(A, RR)) + DD
  } else {
    ## Fallback: geometric series (valid when spectral radius of TT < 1)
    ## H = ZZ * z * sum_{k>=0} (TT z)^k * RR + DD
    ##   = ZZ * sum_{k>=0} TT^k z^{k+1} * RR + DD
    H_ser <- matrix(0 + 0i, nrow = n_obs, ncol = ncol(RR))
    Ak    <- diag(n_state)
    for (k in 0:500) {
      H_ser <- H_ser + ZZ %*% Ak %*% (RR * z^(k + 1L))
      Ak    <- Ak %*% TT
      if (max(abs(Ak)) < 1e-15) break
    }
    H <- H_ser + DD
  }

  S <- H %*% Sigma_e %*% Conj(t(H))
  if (me_variance > 0) S <- S + me_variance * diag(n_obs)
  S
}


## --------------------------------------------------------------------------
## Internal kernel — current-state convention without D (used by D24)
## --------------------------------------------------------------------------

## Transfer function for the current-state convention with no direct-impact D.
##   H = obs_mat * (I − T * z)^{-1} * R,   z = e^{-iω}
##   S_yy(ω) = H * Sigma_e * H^*
##
## This is numerically identical to the original D24 `.spectral_density_at_theta`
## computation.  No z factor on the state; no D matrix.
##
## @param omega     Single frequency in (0, pi].
## @param T_mat     n_state x n_state state-transition matrix.
## @param R_mat     n_state x n_shock shock-impact matrix.
## @param obs_mat   n_obs   x n_state observation matrix.
## @param Sigma_e   n_shock x n_shock shock covariance.
## @return n_obs x n_obs complex Hermitian spectral density matrix.
## @noRd
.spectral_density_core_current_no_d <- function(omega, T_mat, R_mat, obs_mat,
                                                  Sigma_e) {
  n_state <- nrow(T_mat)
  n_obs   <- nrow(obs_mat)
  n_shock <- ncol(R_mat)

  z   <- exp(-1i * omega)
  I_n <- diag(n_state)
  A   <- I_n - T_mat * z

  rcond_A <- rcond(A)
  if (is.finite(rcond_A) && rcond_A > .Machine$double.eps) {
    temp <- solve(A, R_mat)
  } else {
    ## Series fallback: sum_{k>=0} (T z)^k * R = (I - T z)^{-1} R
    H_series <- matrix(0 + 0i, nrow = n_state, ncol = n_shock)
    Ak <- diag(n_state)
    for (k2 in 0:200) {
      H_series <- H_series + Ak %*% (R_mat * z^k2)
      Ak <- Ak %*% T_mat
      if (max(abs(Ak)) < 1e-15) break
    }
    temp <- H_series
  }

  H <- obs_mat %*% temp
  H %*% Sigma_e %*% Conj(t(H))
}


## --------------------------------------------------------------------------
## Public entry point
## --------------------------------------------------------------------------

#' Spectral density of a DSGE state-space model at a single frequency
#'
#' Computes the one-sided spectral density matrix
#' \deqn{S_{yy}(\omega) = H(e^{i\omega})\,\Sigma_e\,H(e^{i\omega})^* + \sigma^2_{\rm me} I}
#' at a single angular frequency \eqn{\omega}.
#'
#' @details
#' The input \code{ss} may carry either the "lagged" or "current" timing
#' convention.  A "current"-tagged object is automatically converted to the
#' lagged-state form via \code{\link{ss_convert_timing}} before computing
#' the transfer function, so both conventions yield the same spectral density.
#'
#' The computation delegates to \code{.spectral_density_core}, the same
#' kernel used internally by the Whittle likelihood.
#'
#' @param ss          A \code{dsge_ss} object (from \code{\link{new_dsge_ss}} or
#'                    \code{\link{build_dsge_state_space}}).
#' @param omega       A single angular frequency in \eqn{(0, \pi]}.
#' @param me_variance Non-negative scalar; added to the diagonal of \eqn{S}
#'                    as measurement-error variance.
#'
#' @return An \eqn{n_{\rm obs} \times n_{\rm obs}} complex Hermitian matrix.
#'
#' @seealso \code{\link{new_dsge_ss}}, \code{\link{ss_convert_timing}}
#' @export
spectral_density <- function(ss, omega, me_variance = 0) {
  if (!inherits(ss, "dsge_ss"))
    stop("spectral_density: 'ss' must be a dsge_ss object.", call. = FALSE)
  ## Normalise to lagged-state convention
  if (ss$timing != "lagged")
    ss <- ss_convert_timing(ss)
  .spectral_density_core(omega,
                          TT        = ss$T_mat,
                          RR        = ss$R_mat,
                          ZZ        = ss$Z_mat,
                          DD        = ss$D_mat,
                          Sigma_e   = ss$Sigma_e,
                          me_variance = me_variance)
}
