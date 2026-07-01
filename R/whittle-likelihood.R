## R/whittle-likelihood.R
## --------------------------------------------------------------------------
## Whittle (frequency-domain) likelihood for dynhr DSGE models.
##
## The Whittle likelihood approximates the Gaussian log-likelihood in the
## frequency domain (Whittle 1951, 1953). For large T it is asymptotically
## equivalent to the exact likelihood. It is faster than the Kalman filter
## for large T (O(T n_obs^2 n_state^2) → O(J n_obs^3) where J ≪ T after
## discarding many frequencies), supports band-restricted estimation (fit
## only on business-cycle frequencies), and pairs naturally with the D23
## spectral identification diagnostics.
##
## References:
##   Whittle, P. (1951). Hypothesis testing in time series analysis.
##   Whittle, P. (1953). The analysis of multiple stationary time series.
##     Journal of the Royal Statistical Society, Series B, 15(1), 125-139.
##   Christiano, L. J., Eichenbaum, M., & Vigfusson, R. (2003). What happens
##     after a technology shock? NBER WP 9819. (band-restricted estimation)
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## 1.  DATA PERIODOGRAM
## --------------------------------------------------------------------------

## Compute the periodogram matrix at Fourier frequencies omega_j = 2*pi*j/T,
## j = 1, ..., floor((T-1)/2).  Frequency zero (the mean) is excluded because
## we demean first.  For even T the Nyquist frequency (j = T/2) is included
## as a real-valued ordinate.
##
## Returns a list with:
##   $omega   -- numeric vector of frequencies in (0, pi]
##   $I       -- n_obs x n_obs x J complex array, I[,,j] = w_j w_j^H / (2*pi*T)
##
## Y must already be DEMEANED (columns = observables, rows = time; T x n_obs).

.whittle_periodogram <- function(Y) {
  ## Y: T x n_obs matrix, already demeaned
  stopifnot(is.matrix(Y))
  T_len  <- nrow(Y)
  n_obs  <- ncol(Y)

  ## Fourier frequencies omega_j = 2*pi*j/T, j = 1..floor((T-1)/2)
  ## Plus Nyquist (j = T/2) if T is even.
  j_max  <- T_len %/% 2L          # floor(T/2)
  j_seq  <- seq_len(j_max)        # j = 1, ..., floor(T/2)
  omega  <- 2 * pi * j_seq / T_len

  ## DFT of each observable column via stats::fft.
  ## stats::fft returns X_k = sum_{t=0}^{T-1} x_t * e^{-2*pi*i*k*t/T}
  ## So w_j = X_j = DFT output at index j+1 (1-based) for j = 0,...,T-1.
  ## The periodogram ordinate is I(omega_j) = w_j w_j^H / (2*pi*T).
  W <- matrix(0 + 0i, nrow = j_max, ncol = n_obs)
  for (k in seq_len(n_obs)) {
    fft_k  <- stats::fft(Y[, k])          # length T, index 1 = frequency 0
    W[, k] <- fft_k[j_seq + 1L]           # j = 1..j_max  => indices 2..j_max+1
  }

  ## Periodogram matrices: I[,,j] = w_j w_j^H / (2*pi*T)
  ## We return as a list of n_obs x n_obs Hermitian matrices.
  scale   <- 1.0 / (2.0 * pi * T_len)
  I_list  <- vector("list", j_max)
  for (j in seq_len(j_max)) {
    wj        <- W[j, , drop = FALSE]     # 1 x n_obs
    I_list[[j]] <- scale * Conj(t(wj)) %*% wj  # n_obs x n_obs
  }

  list(omega = omega, I = I_list, n_obs = n_obs, T_len = T_len)
}


## --------------------------------------------------------------------------
## 2.  MODEL-IMPLIED SPECTRAL DENSITY
## --------------------------------------------------------------------------

## Compute S_yy(omega) for the dynhr state-space model at a single frequency.
##
## dynhr state-space convention (identical to kalman_filter in R/kalman-filter.R):
##   s_t = TT s_{t-1} + RR eps_t          (state transition)
##   y_t = ZZ s_{t-1} + DD eps_t          (observation; uses LAGGED state)
##
## In the z-domain (z = e^{iw}):
##   Y(z) = [ZZ z^{-1} (I - TT z^{-1})^{-1} RR + DD] E(z)
##
## Transfer function: H(e^{iw}) = ZZ e^{-iw} (I - TT e^{-iw})^{-1} RR + DD
## Spectral density:  S_yy(omega) = H Sigma_e H^* + me_variance * I_{n_obs}
##
## Note: this differs from D23's convention (which has no z^{-1} multiplier
## on the state part). D23 uses a present-state observation equation
## y_t = C s_t + D eps_t, whereas kalman_filter uses y_t = ZZ s_{t-1} + DD eps_t.

.whittle_spectral_density <- function(omega, TT, RR, ZZ, DD, Sigma_e,
                                       me_variance = 0)
  .spectral_density_core(omega, TT, RR, ZZ, DD, Sigma_e, me_variance)


## --------------------------------------------------------------------------
## 3.  DEBIASED WHITTLE: AUTOCOVARIANCE SEQUENCE AND EXPECTED PERIODOGRAM
## --------------------------------------------------------------------------
##
## Reference: Sykulski et al. (2019) "The debiased Whittle likelihood."
## Biometrika, 106(2), 251–266.
##
## For a sample of length T, the expected periodogram at Fourier frequency
## omega_j is not S(omega_j) but a Fejer-kernel convolution:
##
##   E[I(omega_j)] = (1/2pi) sum_{|tau|<T} (1 - |tau|/T) c(tau) exp(-i omega_j tau)
##               = (1/2pi) [c(0) + 2 sum_{tau=1}^{T-1} (1-tau/T) c(tau) cos(omega_j tau)]
##
## where c(tau) = Cov(y_t, y_{t+tau}) is the autocovariance function of the
## observed process.  The O(1/T) bias of the standard Whittle likelihood arises
## because S(omega) != E[I(omega)] at finite T.
##
## dynhr lagged-state convention:
##   s_t = TT s_{t-1} + RR eps_t
##   y_t = ZZ s_{t-1} + DD eps_t
##
## Autocovariance (verified in proto; see brief §2.2–2.3):
##   c(0)   = ZZ P0 ZZ' + DD Sigma_e DD' [+ me_variance * I]
##   c(tau) = ZZ TT^{tau-1} K   for tau >= 1
##
## where K = TT P0 ZZ' + RR Sigma_e DD'  (n_state x n_obs cross-covariance kernel)
##
## P0 = solve_lyapunov(TT, RR Sigma_e RR')  is the stationary state covariance.

## .whittle_compute_ctau(TT, RR, ZZ, DD, Sigma_e, T_len, me_variance)
##   -> array[n_obs, n_obs, T_len] of real autocovariances c(0), c(1), ..., c(T-1)
## All entries are REAL (Im part is zero for a stationary real Gaussian process).
.whittle_compute_ctau <- function(TT, RR, ZZ, DD, Sigma_e, T_len, me_variance = 0) {
  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)

  ## Stationary state covariance P0 (Lyapunov solve)
  Q_rr   <- RR %*% Sigma_e %*% t(RR)
  P0     <- solve_lyapunov(TT, Q_rr)

  ## c(0): ZZ P0 ZZ' + DD Sigma_e DD' [+ me*I]
  c0 <- Re(ZZ %*% P0 %*% t(ZZ) + DD %*% Sigma_e %*% t(DD))
  if (me_variance > 0) c0 <- c0 + me_variance * diag(n_obs)

  if (T_len == 1L) {
    arr <- array(0, dim = c(n_obs, n_obs, 1L))
    arr[,,1L] <- c0
    return(arr)
  }

  ## Cross-covariance kernel K = TT P0 ZZ' + RR Sigma_e DD'  (n_state x n_obs)
  K_mat  <- TT %*% P0 %*% t(ZZ) + RR %*% Sigma_e %*% t(DD)

  ## Fill array: index 1 = c(0), index k = c(k-1)
  arr       <- array(0, dim = c(n_obs, n_obs, T_len))
  arr[,,1L] <- c0

  TT_pow <- diag(n_state)           ## TT^0 = I
  for (tau in seq_len(T_len - 1L)) {
    arr[,, tau + 1L] <- Re(ZZ %*% TT_pow %*% K_mat)
    TT_pow <- TT_pow %*% TT
  }

  arr
}


## .whittle_compute_EI(c_arr, omega, T_len)
##   c_arr: array[n_obs, n_obs, T_len] from .whittle_compute_ctau
##   omega: vector of Fourier frequencies (length J)
##   T_len: sample length (integer)
##   -> list of J real symmetric n_obs x n_obs matrices E[I(omega_j)]
##
## Formula (Fejer kernel):
##   EI_j = (1/2pi) * [c(0) + 2 sum_{tau=1}^{T-1} (1-tau/T) c(tau) cos(omega_j tau)]
##
## Eigenvalue clamp: same as in .whittle_loglik (pmax(ev, eps * max(ev)))
## ensures EI is numerically PD even for near-cancellation at high frequencies.
.whittle_compute_EI <- function(c_arr, omega, T_len) {
  J     <- length(omega)
  n_obs <- dim(c_arr)[1L]

  ## Fejer weights for tau = 1, ..., T_len-1
  tau_seq <- seq_len(T_len - 1L)
  w_tau   <- 1.0 - tau_seq / T_len

  ## Precompute cosine table: cos_mat[j, tau] = cos(omega_j * tau)
  ## dims: J x (T_len-1)
  cos_mat <- outer(omega, tau_seq, function(w, t) cos(w * t))

  EI_list <- vector("list", J)
  inv2pi  <- 1.0 / (2.0 * pi)

  for (j in seq_len(J)) {
    EI_j <- c_arr[,,1L]          ## start with c(0); Fejer weight (1 - 0/T) = 1
    if (T_len > 1L) {
      ## For a stationary multivariate real process, the contribution of lags tau
      ## and -tau to the Fejer sum is (1-tau/T) * (c(tau) + c(-tau)) * cos(omega*tau).
      ## Since c(-tau) = c(tau)' (transpose, not c(tau) itself for multivariate!),
      ## the correct coefficient is (c(tau) + c(tau)') per lag, not 2*c(tau).
      ## For the univariate case c(tau) is scalar so c(tau)' = c(tau) and the two
      ## expressions coincide; for n_obs > 1 they differ when c(tau) is asymmetric.
      for (tau in tau_seq) {
        c_tau <- c_arr[,, tau + 1L]
        EI_j  <- EI_j + w_tau[tau] * (c_tau + t(c_tau)) * cos_mat[j, tau]
      }
    }
    EI_j <- EI_j * inv2pi

    ## Clamp: EI should be PD but floating-point cancellation can yield tiny negatives
    ## at high frequencies for persistent processes.  Same rule as eigenvalue clamp
    ## in .whittle_loglik.
    if (n_obs == 1L) {
      ## Scalar path: direct clamp
      EI_j <- pmax(EI_j, .Machine$double.eps * max(abs(EI_j)))
    } else {
      ## Matrix path: symmetrize (should already be symmetric by construction above)
      EI_j <- 0.5 * (EI_j + t(EI_j))   ## force exact symmetry from roundoff
      ## (Eigenvalue clamp happens inside .whittle_loglik when eigen() is called)
    }

    EI_list[[j]] <- EI_j
  }

  EI_list
}


## --------------------------------------------------------------------------
## 3b. TANGENT LYAPUNOV + dc(tau)/dtheta_k FOR THE DEBIASED GRADIENT
## --------------------------------------------------------------------------
##
## For a single parameter k with ss derivatives (dTT, dRR, dZZ, dDD, dSigma_e):
##
##   dP0_k solves the tangent Lyapunov:
##     dP0 = TT dP0 TT' + dTT P0 TT' + TT P0 dTT' + dRR Se RR' + RR dSe RR' + RR Se dRR'
##
##   dc(0)/dtheta_k = ZZ dP0 ZZ' + dZZ P0 ZZ' + ZZ P0 dZZ'
##                    + dDD Se DD' + DD dSe DD' + DD Se dDD'
##
##   dK/dtheta_k   = dTT P0 ZZ' + TT dP0 ZZ' + TT P0 dZZ'
##                    + dRR Se DD' + RR dSe DD' + RR Se dDD'
##
##   dc(tau)/dtheta_k = dZZ TT^{tau-1} K           (observation-matrix term)
##                     + ZZ TT^{tau-1} dK          (cross-covariance kernel term)
##                     + sum_{s=0}^{tau-2} ZZ TT^s dTT TT^{tau-2-s} K   (tau >= 2)
##
## This function returns dcArr [n_obs x n_obs x T_len] of dc(tau)/dtheta_k.
.whittle_compute_dctau <- function(TT, RR, ZZ, DD, Sigma_e, P0, K_mat,
                                    dTT, dRR, dZZ, dDD, dSe, T_len) {
  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)

  ## Tangent Lyapunov: dP0 solves TT dP0 TT' + rhs_k = dP0
  rhs_dP0 <- dTT %*% P0 %*% t(TT) + TT %*% P0 %*% t(dTT) +
              dRR %*% Sigma_e %*% t(RR) + RR %*% dSe %*% t(RR) + RR %*% Sigma_e %*% t(dRR)
  dP0 <- tryCatch(solve_lyapunov(TT, rhs_dP0), error = function(e) NULL)
  if (is.null(dP0)) {
    ## Tangent Lyapunov failed (near unit root): return zero derivatives
    return(array(0, dim = c(n_obs, n_obs, T_len)))
  }

  ## dc(0)/dtheta_k
  dc0 <- Re(ZZ %*% dP0 %*% t(ZZ) +
            dZZ %*% P0 %*% t(ZZ) +
            ZZ  %*% P0 %*% t(dZZ) +
            dDD %*% Sigma_e %*% t(DD) +
            DD  %*% dSe    %*% t(DD) +
            DD  %*% Sigma_e %*% t(dDD))

  if (T_len == 1L) {
    dcArr <- array(0, dim = c(n_obs, n_obs, 1L))
    dcArr[,,1L] <- dc0
    return(dcArr)
  }

  ## dK/dtheta_k (n_state x n_obs)
  dK <- dTT %*% P0 %*% t(ZZ) + TT %*% dP0 %*% t(ZZ) + TT %*% P0 %*% t(dZZ) +
        dRR %*% Sigma_e %*% t(DD) + RR %*% dSe %*% t(DD) + RR %*% Sigma_e %*% t(dDD)

  ## dc(tau)/dtheta_k for tau >= 1 has three contributions:
  ##   A) dZZ * TT^{tau-1} * K  (observation-matrix sensitivity)
  ##   B) ZZ * TT^{tau-1} * dK  (cross-covariance kernel sensitivity)
  ##   C) ZZ * d(TT^{tau-1})/dTT * K  (state-transition power-series deriv)
  ##
  ## Define M_tau = d(TT^{tau-1} K)/dTT in direction dTT (n_state x n_obs):
  ##   M_1 = 0  (d(TT^0 K)/dTT = 0)
  ##   M_{tau+1} = dTT (TT^{tau-1} K) + TT M_tau
  ##             = dTT * TT_pow_K + TT * M_tau
  ## where TT_pow_K = TT^{tau-1} K is maintained across iterations.
  ##
  ## dc(tau)/dtheta_k = dZZ TT^{tau-1} K + ZZ (TT^{tau-1} dK + M_tau)

  dcArr <- array(0, dim = c(n_obs, n_obs, T_len))
  dcArr[,,1L] <- dc0

  TT_pow   <- diag(n_state)          ## TT^0 = I (for tau=1, TT^{tau-1} = I)
  M_tau    <- matrix(0, n_state, n_obs)  ## M_1 = 0
  TT_pow_K <- K_mat                  ## TT^0 K = K  (updated at end of each iteration)

  for (tau in seq_len(T_len - 1L)) {
    ## dc(tau) = dZZ TT^{tau-1} K + ZZ (TT^{tau-1} dK + M_tau)
    dcArr[,, tau + 1L] <- Re(dZZ %*% TT_pow_K + ZZ %*% (TT_pow %*% dK + M_tau))

    ## Update for next tau:
    ##   M_{tau+1} = dTT * TT^{tau-1} K + TT * M_tau = dTT * TT_pow_K + TT * M_tau
    ##   TT_pow_{tau+1} = TT_pow * TT
    ##   TT_pow_K_{tau+1} = TT * TT_pow_K
    M_tau    <- dTT %*% TT_pow_K + TT %*% M_tau
    TT_pow   <- TT_pow %*% TT
    TT_pow_K <- TT %*% TT_pow_K
  }

  dcArr
}


## .whittle_compute_dEI(dcArr, omega, T_len)
##   dcArr: array[n_obs, n_obs, T_len] from .whittle_compute_dctau
##   omega: Fourier frequencies (length J)
##   T_len: sample length
##   -> list of J real symmetric n_obs x n_obs matrices dEI(omega_j)/dtheta_k
.whittle_compute_dEI <- function(dcArr, omega, T_len) {
  J     <- length(omega)
  n_obs <- dim(dcArr)[1L]

  tau_seq <- seq_len(T_len - 1L)
  w_tau   <- 1.0 - tau_seq / T_len
  cos_mat <- outer(omega, tau_seq, function(w, t) cos(w * t))

  dEI_list <- vector("list", J)
  inv2pi   <- 1.0 / (2.0 * pi)

  for (j in seq_len(J)) {
    dEI_j <- dcArr[,,1L]
    if (T_len > 1L) {
      ## Same symmetry fix as .whittle_compute_EI: use (dc(tau) + dc(tau)') not 2*dc(tau).
      for (tau in tau_seq) {
        dc_tau <- dcArr[,, tau + 1L]
        dEI_j  <- dEI_j + w_tau[tau] * (dc_tau + t(dc_tau)) * cos_mat[j, tau]
      }
    }
    dEI_list[[j]] <- dEI_j * inv2pi
  }

  dEI_list
}


## --------------------------------------------------------------------------
## 4.  WHITTLE LOG-LIKELIHOOD
## --------------------------------------------------------------------------

## Evaluate the Whittle log-likelihood given:
##   pdgm   -- periodogram list from .whittle_periodogram()
##   S_fn   -- function(omega) -> n_obs x n_obs complex Hermitian spectral
##             density matrix (can include me_variance on the diagonal)
##   freq_band -- numeric(2): [lo, hi] in radians; only sum over omega in band
##
## Returns a scalar (the Whittle log-likelihood, a real number).
##
## Formula (per Whittle 1953; multivariate form from Dunsmuir & Hannan 1976):
##   ll = sum_j [ -log det S(omega_j) - tr(S(omega_j)^{-1} I(omega_j)) ]
##          - J * n_obs * log(2*pi)   (normalising constant)
##
## where J = number of included frequencies, and we use the REAL part of the
## trace (the imaginary part is zero for consistent spectral density estimates).

## debias:  logical — if TRUE and EI_list is provided, substitute EI_list[[j]]
##          for S_fn(omega[j]).  When FALSE, behaviour is BIT-IDENTICAL to the
##          original code (S_fn path; EI_list is ignored).
## EI_list: optional precomputed expected periodogram (list of J real symmetric
##          matrices from .whittle_compute_EI).  Ignored when debias=FALSE.
.whittle_loglik <- function(pdgm, S_fn, freq_band = c(0, pi),
                             debias = FALSE, EI_list = NULL) {
  omega   <- pdgm$omega
  I_list  <- pdgm$I
  n_obs   <- pdgm$n_obs
  J_total <- length(omega)

  ## When debiasing, EI_list must cover the full frequency grid (all J)
  ## so that subsetting by in_band is correct.  Checked here once.
  use_debias <- isTRUE(debias) && !is.null(EI_list)

  ## Select frequencies inside the band
  lo <- freq_band[1L]; hi <- freq_band[2L]
  in_band <- which(omega > lo & omega <= hi)

  if (length(in_band) == 0L)
    stop("whittle_loglik: no Fourier frequencies fall inside freq_band = [",
         lo, ", ", hi, "]. Widen the band or increase T.")

  ll      <- 0
  ll_const <- -0.5 * n_obs * log(2 * pi)   # per frequency

  for (j in in_band) {
    ## Spectral density matrix at this frequency:
    ##   debias=FALSE  -> model S(omega_j)  [complex Hermitian]
    ##   debias=TRUE   -> E[I(omega_j)]     [real symmetric]
    S <- if (use_debias) EI_list[[j]] else S_fn(omega[j])

    ## Exact complex Hermitian path — works in base R for any n_obs.
    ## eigen(S, symmetric=TRUE) returns real eigenvalues for Hermitian S;
    ## For the debiased path S is real symmetric so symmetric=TRUE is exact.
    ## solve(S) works on complex matrices in base R.
    ev <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
    if (is.null(ev)) return(-Inf)
    ## Clamp near-zero eigenvalues from roundoff: reject only when S is
    ## identically zero/negative (max(ev) <= 0).  A genuinely PD S whose
    ## smallest eigenvalue is < eps*||S|| due to floating-point error is
    ## clamped and accepted; this matches the tolerance a chol-with-jitter
    ## would allow.  Same rule applied to EI_j (PD up to roundoff).
    ev_max <- max(ev$values)
    if (ev_max <= 0) return(-Inf)
    ev_clamped <- pmax(ev$values, .Machine$double.eps * ev_max)
    log_det_S <- sum(log(ev_clamped))       # log|det S| = sum(log(lambda_k))
    ## Rebuild S_inv from clamped eigenvalues so loglik and inverse are
    ## consistent (V %*% diag(1/ev_clamped) %*% t(V)).
    V     <- ev$vectors
    S_inv <- V %*% diag(1 / ev_clamped, nrow = length(ev_clamped)) %*% Conj(t(V))

    I_j   <- I_list[[j]]                   # n_obs x n_obs Hermitian periodogram
    tr_SI <- Re(sum(diag(S_inv %*% I_j)))  # Im part is O(eps) by Hermitian product

    ll_j  <- ll_const - 0.5 * log_det_S - 0.5 * tr_SI
    if (!is.finite(ll_j)) return(-Inf)
    ll <- ll + ll_j
  }

  ll
}


## --------------------------------------------------------------------------
## 4.  LOG-POSTERIOR FACTORY  (parallel to make_log_posterior in posterior.R)
## --------------------------------------------------------------------------

#' Create a Whittle-likelihood log-posterior evaluator
#'
#' Returns a closure \code{function(theta) -> list(logpost, loglik, logprior)}
#' that uses the Whittle (frequency-domain) likelihood instead of the Kalman
#' filter.
#'
#' @param model       dynhr_mod (from \code{parse_mod()})
#' @param data        Observation matrix (T x n_obs), columns = obs_vars
#' @param prior_spec  Prior specification data.frame
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled (from \code{compile_model()})
#' @param me_variance Scalar measurement-error variance added to diag(S(omega))
#'   at every frequency. Whittle assumes stationarity, so only a
#'   time-invariant (scalar or per-observable diagonal) ME is supported.
#'   Passed as a scalar; applied uniformly across observables.
#' @param freq_band   Numeric(2): \code{c(lo, hi)} in radians. Only Fourier
#'   frequencies omega_j with \code{lo < omega_j <= hi} contribute to the
#'   Whittle log-likelihood. Default \code{c(0, pi)} uses the full spectrum.
#'   Use \code{whittle_business_cycle_band()} for the standard 6-32 quarter
#'   business-cycle band.
#' @return function(theta) -> list(logpost, loglik, logprior)
#' @noRd
make_log_posterior_whittle <- function(model, data, prior_spec, obs_vars,
                                        compiled, me_variance = 0,
                                        freq_band = c(0, pi),
                                        system_priors = NULL,
                                        debias = TRUE) {
  ## Guard: filter_tunes are a time-domain feature (they insert NA data at
  ## specific periods); the Whittle likelihood assumes a complete, stationary
  ## panel. Stop early with a clear message rather than silently wrong results.
  tunes <- model$filter_tunes$tunes
  if (!is.null(tunes) && is.data.frame(tunes) && nrow(tunes) > 0L)
    stop("make_log_posterior_whittle: model has filter_tunes. ",
         "filter_tunes are a time-domain feature incompatible with the ",
         "Whittle (frequency-domain) likelihood. ",
         "Use likelihood = \"gaussian\" instead.",
         call. = FALSE)

  ## Validate freq_band
  if (!is.numeric(freq_band) || length(freq_band) != 2L ||
      !all(is.finite(freq_band)) || freq_band[1L] < 0 ||
      freq_band[2L] > pi || freq_band[1L] >= freq_band[2L])
    stop("make_log_posterior_whittle: freq_band must be c(lo, hi) with ",
         "0 <= lo < hi <= pi.", call. = FALSE)

  ## Prepare data: T x n_obs matrix
  if (is.null(dim(data))) data <- matrix(data, ncol = 1L)
  if (nrow(data) < ncol(data)) data <- t(data)        # ensure T x n_obs
  if (is.null(colnames(data))) colnames(data) <- obs_vars
  T_len <- nrow(data)
  n_obs <- length(obs_vars)

  ## Demean each observable (Whittle assumes zero mean / stationary data)
  data_col_idx <- match(obs_vars, colnames(data))
  if (any(is.na(data_col_idx)))
    stop("make_log_posterior_whittle: obs_vars not found in data columns.")
  Y_raw    <- data[, data_col_idx, drop = FALSE]
  col_means <- colMeans(Y_raw, na.rm = TRUE)
  Y_demeaned <- sweep(Y_raw, 2, col_means, "-")
  if (anyNA(Y_demeaned))
    stop("make_log_posterior_whittle: data contains NA after demeaning. ",
         "Whittle likelihood requires a complete, balanced panel. ",
         "Use likelihood = \"gaussian\" for missing-data support.",
         call. = FALSE)

  ## Pre-compute the periodogram once for this dataset (it does not depend on
  ## model parameters, only on the data).
  pdgm <- .whittle_periodogram(Y_demeaned)

  ## Cache system structure for fast repeated evaluation
  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  sys_cache <- cache_system_structure(compiled)

  ## Warm-start cache for steady-state solve (same trick as gaussian path)
  ss_warm <- NULL

  ## Observable indices (used to slice ghx, ghu identically to kalman_filter)
  exo_names <- model$varexo_names

  function(theta) {
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    params <- .apply_theta_to_params(model, theta)

    ss_result <- solve_steady_state(model, compiled, params,
                                    y0 = ss_warm, verbose = FALSE)
    if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
      if (!is.null(ss_warm))
        ss_result <- solve_steady_state(model, compiled, params,
                                        verbose = FALSE)
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
        ss_warm <<- NULL
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
      }
    }
    ss_warm <<- ss_result$ss

    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
    dr  <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Stationarity guard: Whittle spectral density is undefined at unit roots.
    ## Same logic as the gaussian path in make_log_posterior().
    ns <- length(dr$state_idx)
    ev <- dr$eigenvalues
    spectral_radius <- if (!is.null(ev) && length(ev) >= ns)
      max(Mod(ev[seq_len(ns)]))
    else
      max(Mod(eigen(dr$ghx[dr$state_idx, , drop = FALSE],
                    only.values = TRUE)$values))
    if (spectral_radius >= 1)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Extract state-space matrices (same slicing as kalman_filter)
    state_idx <- dr$state_idx
    endo      <- dr$endo_names
    ghx       <- dr$ghx
    ghu       <- dr$ghu
    obs_idx   <- match(obs_vars, endo)
    if (any(is.na(obs_idx)))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Construct dsge_ss at the Whittle boundary (timing = "lagged").
    ## Unpack to locals immediately — .whittle_spectral_density is a hot kernel
    ## that must not incur S3 dispatch per frequency.
    Sigma_e  <- .get_shock_cov(model, exo_names, params)
    ss <- new_dsge_ss(
      T_mat   = ghx[state_idx, , drop = FALSE],
      R_mat   = ghu[state_idx, , drop = FALSE],
      Z_mat   = ghx[obs_idx,   , drop = FALSE],
      D_mat   = ghu[obs_idx,   , drop = FALSE],
      Sigma_e = Sigma_e,
      timing  = "lagged"
    )
    TT <- ss$T_mat; RR <- ss$R_mat; ZZ <- ss$Z_mat; DD <- ss$D_mat

    ## Build spectral-density function for this parameter draw
    me_var   <- me_variance    # capture in closure
    S_fn     <- function(omega)
      .whittle_spectral_density(omega, TT, RR, ZZ, DD, Sigma_e, me_var)

    ## Debiased Whittle: precompute E[I(omega_j)] for all J frequencies.
    ## Re-computed every draw because TT, ZZ, Sigma_e (and hence P0) vary
    ## with theta.  Cost: O(T * n_state^2) for c(tau) + O(T * J * n_obs^2)
    ## for the Fejer sum — both negligible vs the steady-state solve.
    EI_list_draw <- if (isTRUE(debias)) {
      tryCatch({
        c_arr <- .whittle_compute_ctau(TT, RR, ZZ, DD, Sigma_e, T_len, me_var)
        .whittle_compute_EI(c_arr, pdgm$omega, T_len)
      }, error = function(e) NULL)
    } else {
      NULL
    }

    ll <- tryCatch(
      .whittle_loglik(pdgm, S_fn, freq_band,
                      debias = isTRUE(debias) && !is.null(EI_list_draw),
                      EI_list = EI_list_draw),
      error = function(e) {
        ## Re-throw configuration errors (empty band) immediately; only swallow
        ## numerical errors that arise during spectral density evaluation.
        if (grepl("no Fourier frequencies", conditionMessage(e), fixed = TRUE))
          stop(e)
        -Inf
      }
    )
    if (!is.finite(ll))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## System priors: penalty on model features evaluated from the solved dr.
    if (!is.null(system_priors)) {
      sp_lp <- .eval_system_priors(
        system_priors,
        list(theta   = theta,
             model   = model,
             dr      = dr,
             Sigma_e = .get_shock_cov(model, exo_names, params),
             params  = params))
      if (!is.finite(sp_lp))
        return(list(logpost = -Inf, loglik = ll, logprior = lp))
      lp <- lp + sp_lp
    }

    list(logpost = .dynhr_opt("power_posterior", default = 1) * ll + lp,
         loglik = ll, logprior = lp)
  }
}


## --------------------------------------------------------------------------
## 5.  WHITTLE ANALYTIC GRADIENT
## --------------------------------------------------------------------------

## Compute d log-lik / d theta_j for the Whittle likelihood, given
## per-parameter solution derivatives and the periodogram + frequency grid.
##
## Arguments:
##   pdgm       -- periodogram list from .whittle_periodogram()
##   TT, RR, ZZ, DD, Sigma_e -- state-space matrices at the current theta
##   d_ss_list  -- named list; each element is a list with (some of):
##                   dTT, dRR, dZZ, dDD (solution-mover derivatives)
##                   dSigma_e            (shock-cov derivative, or NULL)
##                 NULL element signals "use FD" (skip analytic for that param).
##   freq_band  -- numeric(2): [lo, hi] in radians; must match the loglik band
##   me_variance -- closure constant (scalar); gradient w.r.t. it is zero
##
## Returns a named numeric vector of length == length(d_ss_list), one entry
## per parameter.  NULL slots return NA_real_ (caller should fall back to FD).
##
## Inner loop: per-frequency precompute (z, A, B = z*A^{-1}*RR, G=ZZ*B,
## H, S, eigen(S), S^{-1}) then per-parameter accumulate via
##   d ll / d theta_j = sum_{k in band} Re[ -0.5 tr(Si dS) + 0.5 tr(Si dS Si I) ]
## where dS = dH Sigma_e H^H + H dSigma_e H^H + H Sigma_e dH^H.
## Si is the full complex inverse of the Hermitian spectral density matrix;
## I_j is the full complex Hermitian periodogram. Re() applied only at the
## final scalar — the imaginary part of tr(Si dS Si I) is O(eps) by Hermitian
## structure and serves as a numerical guard only.
##
## z-factor convention: z = exp(-i*omega), A = I - TT*z (matches loglik).
## freq_band edge semantics: omega > lo & omega <= hi (matches loglik line 151).
## debias:  logical (default FALSE for bit-identical to historical behaviour).
##          When TRUE, EI_list and dEI_list must be provided; the gradient is
##          then computed as d ll_debias / d theta_k using EI_j / dEI_j instead
##          of S_j / dS_j.  The formula is structurally identical (pass-through
##          theorem, verified to 1e-14 in the debiased-Whittle scout).
##
## EI_list:  list of J real symmetric n_obs x n_obs matrices (from
##           .whittle_compute_EI), indexed by full frequency grid j=1..J.
##           Ignored when debias=FALSE.
##
## dEI_per_param: named list (same names as d_ss_list), each element is a
##           list of J n_obs x n_obs dEI matrices from .whittle_compute_dEI.
##           Ignored when debias=FALSE.
.whittle_loglik_grad <- function(pdgm, TT, RR, ZZ, DD, Sigma_e,
                                  d_ss_list, freq_band = c(0, pi),
                                  me_variance = 0,
                                  debias = FALSE,
                                  EI_list = NULL,
                                  dEI_per_param = NULL) {
  omega   <- pdgm$omega
  I_list  <- pdgm$I
  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)
  np      <- length(d_ss_list)
  nm_all  <- names(d_ss_list)

  use_debias <- isTRUE(debias) && !is.null(EI_list) && !is.null(dEI_per_param)

  ## Same band mask as .whittle_loglik (strict-left-open, closed-right)
  lo <- freq_band[1L]; hi <- freq_band[2L]
  in_band <- which(omega > lo & omega <= hi)

  ## Zero-matrices for missing derivative slots
  zero_TT <- matrix(0, n_state, n_state)
  zero_RR <- matrix(0, n_state, n_exo)
  zero_ZZ <- matrix(0, n_obs, n_state)
  zero_DD <- matrix(0, n_obs, n_exo)
  zero_Se <- matrix(0, n_exo, n_exo)

  grad     <- setNames(rep(0, np), nm_all)
  na_slots <- vapply(d_ss_list, is.null, logical(1))
  grad[na_slots] <- NA_real_

  ## Skip entirely if all slots are NULL (no analytic contribution)
  active <- which(!na_slots)
  if (length(active) == 0L || length(in_band) == 0L) return(grad)

  ## Extract derivative arrays once for all parameters
  dTT_arr <- lapply(d_ss_list, function(d) if (!is.null(d)) d$dTT  %||% zero_TT else NULL)
  dRR_arr <- lapply(d_ss_list, function(d) if (!is.null(d)) d$dRR  %||% zero_RR else NULL)
  dZZ_arr <- lapply(d_ss_list, function(d) if (!is.null(d)) d$dZZ  %||% zero_ZZ else NULL)
  dDD_arr <- lapply(d_ss_list, function(d) if (!is.null(d)) d$dDD  %||% zero_DD else NULL)
  dSe_arr <- lapply(d_ss_list, function(d) if (!is.null(d)) d$dSigma_e %||% zero_Se else NULL)

  n_skipped <- 0L

  ## -----------------------------------------------------------------------
  ## DEBIASED PATH: use precomputed EI_list and dEI_per_param
  ## -----------------------------------------------------------------------
  if (use_debias) {
    for (j in in_band) {
      EI_j <- EI_list[[j]]                  ## real symmetric n_obs x n_obs
      I_j  <- I_list[[j]]

      ev_j <- tryCatch(eigen(EI_j, symmetric = TRUE), error = function(e) NULL)
      if (is.null(ev_j)) next
      ev_max_j <- max(ev_j$values)
      if (ev_max_j <= 0) next
      ev_clamped_j <- pmax(ev_j$values, .Machine$double.eps * ev_max_j)
      V_j   <- ev_j$vectors
      EIi_j <- V_j %*% diag(1 / ev_clamped_j, nrow = length(ev_clamped_j)) %*% t(V_j)
      ## Note: EI is real symmetric so t(V) not Conj(t(V))

      ## Per-parameter: grad += Re[-0.5 tr(EIi dEI) + 0.5 tr(EIi dEI EIi I)]
      for (k in active) {
        dEI_k_j <- dEI_per_param[[k]][[j]]   ## real symmetric n_obs x n_obs
        EIi_dEI <- EIi_j %*% dEI_k_j
        grad[k]  <- grad[k] +
          Re(-0.5 * sum(diag(EIi_dEI)) +
              0.5 * sum(diag(EIi_dEI %*% (EIi_j %*% I_j))))
      }
    }

    if (n_skipped > 0L)
      warning(sprintf(
        ".whittle_loglik_grad: skipped %d near-singular-EI frequency(ies).",
        n_skipped), call. = FALSE)

    return(grad)
  }

  ## -----------------------------------------------------------------------
  ## STANDARD PATH (debias=FALSE): BIT-IDENTICAL to historical behaviour
  ## -----------------------------------------------------------------------
  for (j in in_band) {
    om_j <- omega[j]
    z_j  <- exp(-1i * om_j)

    ## A = I - TT*z (same as .whittle_spectral_density line 99)
    A_j     <- diag(n_state) - TT * z_j
    rcond_A <- rcond(A_j)

    if (!is.finite(rcond_A) || rcond_A <= .Machine$double.eps * 1e4) {
      ## Near-unit-root at this frequency: skip gradient contribution.
      ## The spectral_radius >= 1 guard in the log-posterior already rejects
      ## unit-root draws, so this branch is self-limiting in live MCMC.
      n_skipped <- n_skipped + 1L
      next
    }

    ## B_j = z * A^{-1} * RR  (n_state x n_exo; the H computation uses z*B_raw)
    B_j <- z_j * solve(A_j, RR)
    ## G_j = ZZ * B_j  (n_obs x n_exo), H_j = G_j + DD
    G_j <- ZZ %*% B_j
    H_j <- G_j + DD

    ## S_j = H_j Sigma_e H_j^H [+ me_variance * I]
    S_j <- H_j %*% Sigma_e %*% Conj(t(H_j))
    if (me_variance > 0) S_j <- S_j + me_variance * diag(n_obs)

    ## Exact complex Hermitian path — consistent with .whittle_loglik.
    ## Clamp near-zero eigenvalues (same rule as in .whittle_loglik) so that
    ## loglik and gradient are built from the same clamped spectral density.
    ev_j <- tryCatch(eigen(S_j, symmetric = TRUE), error = function(e) NULL)
    if (is.null(ev_j)) next
    ev_max_j <- max(ev_j$values)
    if (ev_max_j <= 0) next   # identically zero / negative — skip frequency
    ev_clamped_j <- pmax(ev_j$values, .Machine$double.eps * ev_max_j)
    V_j  <- ev_j$vectors
    Si_j <- V_j %*% diag(1 / ev_clamped_j, nrow = length(ev_clamped_j)) %*% Conj(t(V_j))

    I_j <- I_list[[j]]                 # n_obs x n_obs full complex Hermitian periodogram

    ## Per-parameter inner loop
    for (k in active) {
      dTT_k <- dTT_arr[[k]]
      dRR_k <- dRR_arr[[k]]
      dZZ_k <- dZZ_arr[[k]]
      dDD_k <- dDD_arr[[k]]
      dSe_k <- dSe_arr[[k]]

      ## dH_k from the brief:
      ##   dH = dZZ * (z*B_raw) + ZZ * z * A^{-1} * (dTT * z * B_raw + dRR) + dDD
      ## With B_j = z * B_raw, so z * B_raw = B_j and B_raw = B_j / z:
      ##   inner = dTT * B_j + dRR_k   [note: dTT*z*B_raw = dTT*B_j]
      inner_k <- dTT_k %*% B_j + dRR_k   # n_state x n_exo
      dG_k    <- dZZ_k %*% B_j + ZZ %*% (z_j * solve(A_j, inner_k))
      dH_k    <- dG_k + dDD_k

      ## dS_k = dH Sigma_e H^H + H dSigma_e H^H + H Sigma_e dH^H
      dS_k <- dH_k %*% Sigma_e  %*% Conj(t(H_j)) +
              H_j  %*% dSe_k    %*% Conj(t(H_j)) +
              H_j  %*% Sigma_e  %*% Conj(t(dH_k))

      ## Exact gradient formula (from ll = -0.5 log det S - 0.5 tr(S^{-1} I)):
      ##   d ll / d theta_k = Re[ -0.5 tr(S^{-1} dS_k) + 0.5 tr(S^{-1} dS_k S^{-1} I_j) ]
      ## Si_j is complex (from solve(S_j)); dS_k is complex; I_j is full complex.
      ## Re() applied only at the final scalar — the imaginary part of the trace
      ## is O(eps) for Hermitian-structured matrices, so Re() is a numerical guard.
      SidS    <- Si_j %*% dS_k          # complex n x n
      grad[k] <- grad[k] +
        Re(-0.5 * sum(diag(SidS)) +
            0.5 * sum(diag(SidS %*% (Si_j %*% I_j))))
    }
  }

  if (n_skipped > 0L)
    warning(sprintf(
      ".whittle_loglik_grad: skipped %d near-singular-A frequency(ies); ",
      "gradient contribution set to zero for those frequencies.",
      n_skipped), call. = FALSE)

  grad
}


## --------------------------------------------------------------------------
## 6.  CONVENIENCE: BUSINESS-CYCLE BAND
## --------------------------------------------------------------------------

#' Business-cycle frequency band for the Whittle likelihood
#'
#' Returns \code{c(lo, hi)} in radians corresponding to periodicities between
#' \code{period_lo} and \code{period_hi} periods per unit time.  For quarterly
#' data the conventional business-cycle band is 6 to 32 quarters.
#'
#' @param period_lo Minimum periodicity (shortest cycles), default 6.
#' @param period_hi Maximum periodicity (longest cycles), default 32.
#' @return Numeric vector \code{c(lo, hi)} in radians, suitable for the
#'   \code{freq_band} argument of \code{make_log_posterior()} and
#'   \code{run_mode_finding()}.
#' @examples
#' ## Business cycle (6--32 quarters)
#' whittle_business_cycle_band()          # c(0.196, 1.047)
#' ## Annual-to-medium cycles (4--20 years = 16--80 quarters)
#' whittle_business_cycle_band(16, 80)
#' @export
whittle_business_cycle_band <- function(period_lo = 6, period_hi = 32) {
  if (period_lo <= 0 || period_hi <= period_lo)
    stop("whittle_business_cycle_band: need 0 < period_lo < period_hi.")
  ## omega = 2*pi / period; shorter period = higher frequency
  c(lo = 2 * pi / period_hi, hi = 2 * pi / period_lo)
}
