## R/hank-estimation.R
## --------------------------------------------------------------------------
## Estimation bridge for the linearized HANK/GE aggregate dynamics.
##
## The linear GE solve implies a moving-average representation of the aggregate
## observables in the structural shock innovations: the impulse responses ARE
## the MA coefficients (Theta_s).  Following ABRS (2021, eq. 37) the exact
## Gaussian likelihood of aggregate data is built directly from the implied
## autocovariance function
##
##   Gamma(k) = sum_s Theta_{s+k} Sigma_eps Theta_s'
##
## via the block-Toeplitz covariance of the stacked sample -- no state-space or
## Kalman filter required (though the same aggregate law could be packed into
## dynhr's dsge_ss for the Kalman/Whittle/SBC tooling as a follow-up).
##
## Only a scalar shock (KS TFP) is wired here; the MA/autocovariance machinery
## is written for a general n_shock Sigma_eps.  With more observables than
## shocks the model is stochastically singular, so measurement error must be
## supplied (me_var) to keep the covariance positive definite.
## --------------------------------------------------------------------------


#' Moving-average coefficients of aggregate observables (TFP innovation)
#'
#' Computes the impulse responses of the chosen aggregate observables to a unit
#' TFP innovation, which are the MA coefficients \eqn{\Theta_s} of the
#' linearized HANK/GE aggregate law.
#'
#' @param ks A \code{\link{hank_ks_steady}} steady state.
#' @param rho_z Numeric in [0,1): AR(1) persistence of (log) TFP.
#' @param T_h Integer horizon (number of MA terms).
#' @param observables Character subset of \code{c("Y","C","K","r","w")}.
#' @param ge Optional precomputed \code{\link{hank_ks_ge_jacobian}}.
#'
#' @return A \code{T_h x n_obs} matrix of MA coefficients (columns named by
#'   \code{observables}).
#' @export
hank_ma_coefficients <- function(ks, rho_z, T_h,
                                 observables = c("Y", "C"), ge = NULL) {
  if (is.null(ge)) ge <- hank_ks_ge_jacobian(ks, T_h)
  ## Deterministic TFP response to a unit innovation at t=0: dlogZ_t = rho^t.
  dZ  <- ks$Z * rho_z^(seq_len(T_h) - 1L)
  irf <- hank_ks_linear_irf(ks, dZ, ge = ge)
  fld <- c(Y = "dY", C = "dC", K = "dK", r = "dr", w = "dw")
  Theta <- sapply(observables, function(o) irf[[fld[[o]]]])
  matrix(Theta, T_h, length(observables),
         dimnames = list(NULL, observables))
}


#' Autocovariance function from MA coefficients
#'
#' @param Theta \code{T_h x n_obs} MA coefficients (see
#'   \code{\link{hank_ma_coefficients}}).
#' @param sigma_eps Numeric: standard deviation of the (scalar) shock innovation.
#' @param n_lags Integer: maximum lag \eqn{k} to return.
#'
#' @return An \code{n_obs x n_obs x (n_lags+1)} array with slice \code{k+1}
#'   giving \eqn{\Gamma(k)}.
#' @export
hank_autocov <- function(Theta, sigma_eps, n_lags) {
  T_h <- nrow(Theta); n_obs <- ncol(Theta)
  v <- sigma_eps^2
  G <- array(0, dim = c(n_obs, n_obs, n_lags + 1L))
  for (k in 0:n_lags) {
    s_idx <- seq_len(T_h - k)
    ## sum_s Theta_{s+k} Theta_s'
    G[, , k + 1L] <- v * crossprod(Theta[s_idx + k, , drop = FALSE],
                                   Theta[s_idx, , drop = FALSE])
  }
  G
}


#' Block-Toeplitz covariance of a stacked aggregate sample
#'
#' @param G Autocovariance array from \code{\link{hank_autocov}} (max lag must
#'   be >= \code{T_data - 1}).
#' @param T_data Integer: number of time periods in the sample.
#' @param me_var Numeric >= 0: measurement-error variance added on the diagonal
#'   (required when n_obs > n_shocks to avoid singularity).
#'
#' @return A \code{(T_data*n_obs) x (T_data*n_obs)} covariance matrix, stacked
#'   time-major (row \code{(t-1)*n_obs + j}).
#' @keywords internal
.hank_stacked_cov <- function(G, T_data, me_var = 0) {
  n_obs <- dim(G)[1]; n_lags <- dim(G)[3] - 1L
  if (n_lags < T_data - 1L)
    stop("autocovariance max lag < T_data - 1; increase T_h/n_lags")
  N <- T_data * n_obs
  S <- matrix(0, N, N)
  for (t in seq_len(T_data)) for (s in seq_len(T_data)) {
    k <- t - s
    blk <- if (k >= 0) G[, , k + 1L] else t(G[, , -k + 1L])
    S[((t - 1L) * n_obs + 1L):(t * n_obs),
      ((s - 1L) * n_obs + 1L):(s * n_obs)] <- blk
  }
  if (me_var > 0) diag(S) <- diag(S) + me_var
  (S + t(S)) / 2   # symmetrize against round-off
}


#' Gaussian log-likelihood of aggregate HANK data (ABRS autocovariance form)
#'
#' Exact time-domain Gaussian likelihood of the stacked demeaned aggregate
#' sample under the MA representation, via the block-Toeplitz covariance and a
#' Cholesky solve.
#'
#' @param Y \code{T_data x n_obs} matrix of demeaned aggregate observations
#'   (deviations from steady state), columns in the same order as \code{Theta}.
#' @param Theta MA coefficients (see \code{\link{hank_ma_coefficients}}).
#' @param sigma_eps Shock innovation standard deviation.
#' @param me_var Measurement-error variance (default 0).
#'
#' @return The scalar Gaussian log-likelihood.
#' @export
hank_loglik_aggregate <- function(Y, Theta, sigma_eps, me_var = 0) {
  Y <- as.matrix(Y)
  T_data <- nrow(Y); n_obs <- ncol(Y)
  G <- hank_autocov(Theta, sigma_eps, n_lags = T_data - 1L)
  S <- .hank_stacked_cov(G, T_data, me_var = me_var)
  yv <- as.numeric(t(Y))                       # time-major stacking
  ch <- chol(S)                                # upper triangular
  z  <- backsolve(ch, yv, transpose = TRUE)    # solve t(ch) z = yv
  logdet <- 2 * sum(log(diag(ch)))
  N <- T_data * n_obs
  -0.5 * (N * log(2 * pi) + logdet + sum(z^2))
}


#' State-space (companion) form of the aggregate MA representation
#'
#' Packs the moving-average coefficients into a finite lagged-state-space
#' \code{dsge_ss} (companion / shift-register form): the state stacks the last
#' \code{q} shock innovations, so \code{y_t = sum_{s<q} Theta_s eps_{t-s}} is an
#' exact MA(\code{q-1}).  This is the bridge that lets dynhr's Kalman
#' filter/smoother and frequency-domain tooling estimate a linearized HANK model
#' (the state-transition is nilpotent, so the process is stationary and the
#' initial state covariance is \code{sigma_eps^2 I}).
#'
#' @param Theta \code{T_h x n_obs} MA coefficients (see
#'   \code{\link{hank_ma_coefficients}}).
#' @param sigma_eps Innovation standard deviation (scalar shock).
#' @param q Integer state length (number of MA terms retained); default all
#'   \code{nrow(Theta)}.
#'
#' @return A lagged-timing \code{dsge_ss} object (see \code{\link{new_dsge_ss}}).
#' @export
hank_ma_state_space <- function(Theta, sigma_eps, q = NULL) {
  Theta <- as.matrix(Theta)
  T_h <- nrow(Theta); n_obs <- ncol(Theta)
  if (is.null(q)) q <- T_h
  q <- min(q, T_h)

  ## Shift register: s_t = TT s_{t-1} + RR eps_t, s_t = (eps_t,...,eps_{t-q+1}).
  TT <- matrix(0, q, q)
  if (q >= 2L) TT[cbind(2:q, 1:(q - 1L))] <- 1
  RR <- matrix(0, q, 1); RR[1, 1] <- 1

  ## Lagged observation:  y_t = Z_lag s_{t-1} + D_lag eps_t, with
  ## Z_lag = [Theta_1, ..., Theta_{q-1}, 0], D_lag = Theta_0.
  Z_lag <- matrix(0, n_obs, q)
  if (q >= 2L) for (k in 1:(q - 1L)) Z_lag[, k] <- Theta[k + 1L, ]
  D_lag <- matrix(Theta[1L, ], n_obs, 1)

  new_dsge_ss(T_mat = TT, R_mat = RR, Z_mat = Z_lag, D_mat = D_lag,
              Sigma_e = matrix(sigma_eps^2, 1, 1),
              obs_names = colnames(Theta), timing = "lagged")
}


#' Kalman log-likelihood of aggregate HANK data (state-space bridge)
#'
#' Runs dynhr's Kalman filter on the MA companion state space
#' (\code{\link{hank_ma_state_space}}).  Equivalent to the direct autocovariance
#' likelihood \code{\link{hank_loglik_aggregate}} (both are the exact Gaussian
#' likelihood of the same MA process) but routed through the Kalman engine, so it
#' composes with dynhr's missing-data handling, smoother, priors, samplers and
#' SBC tooling.
#'
#' @param Y \code{T_data x n_obs} matrix of demeaned aggregate observations.
#' @param Theta MA coefficients (see \code{\link{hank_ma_coefficients}}).
#' @param sigma_eps Innovation standard deviation.
#' @param me_var Measurement-error variance (needed when n_obs > n_shocks).
#' @param q Optional state length (default all MA terms).
#'
#' @return The scalar Kalman log-likelihood.
#' @export
hank_loglik_ss <- function(Y, Theta, sigma_eps, me_var = 0, q = NULL) {
  Y  <- as.matrix(Y)
  ss <- hank_ma_state_space(Theta, sigma_eps, q = q)
  n_state <- ss$n_state
  out <- .kf_univariate_dispatch(
    Y_minus_d   = t(Y),                       # core expects n_obs x T
    ZZ = ss$Z_mat, TT = ss$T_mat, RR = ss$R_mat, DD = ss$D_mat,
    Sigma_e = ss$Sigma_e,
    s0 = rep(0, n_state),
    P_state = sigma_eps^2 * diag(n_state),    # nilpotent T => stationary cov
    P_inf_state = NULL,
    me_variance = me_var)
  out$loglik
}


#' Simulate an aggregate sample from the MA representation
#'
#' @param Theta MA coefficients.
#' @param sigma_eps Innovation standard deviation.
#' @param T_data Sample length.
#' @param me_sd Measurement-error standard deviation (default 0).
#' @param eps Optional pre-drawn innovations (length \code{T_data}); if
#'   \code{NULL}, drawn N(0, sigma_eps^2).  Provide for reproducibility.
#'
#' @return A \code{T_data x n_obs} matrix of simulated demeaned observations.
#' @export
hank_simulate_aggregate <- function(Theta, sigma_eps, T_data, me_sd = 0,
                                    eps = NULL) {
  T_h <- nrow(Theta); n_obs <- ncol(Theta)
  if (is.null(eps)) eps <- stats::rnorm(T_data, 0, sigma_eps)
  Y <- matrix(0, T_data, n_obs)
  for (t in seq_len(T_data)) {
    s <- 0:min(t - 1L, T_h - 1L)               # Theta_s eps_{t-s}
    Y[t, ] <- colSums(Theta[s + 1L, , drop = FALSE] * eps[t - s])
  }
  if (me_sd > 0) Y <- Y + matrix(stats::rnorm(T_data * n_obs, 0, me_sd),
                                 T_data, n_obs)
  Y
}
