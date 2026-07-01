## R/pskf-smoother.R
## --------------------------------------------------------------------------
## PSKF smoother: CSN backward pass on the CSN forward-pass moments.
##
## NOT exact-on-paper: the nu-shift backward update below keeps Gamma and
## Delta fixed at their filtered values (the full GMTrede smoothed CSN
## enlarges the skew dimension each backward step).  Measured residual vs
## the dense-grid oracle: 4.7e-4 mean gap at alpha=2 (T=3, no pruning) —
## ~10x tighter than the 4.7e-3 two-moment Gaussian-RTS gap, but a
## first-order approximation, not the exact recursion.
##
## Reference: Guljanov, Mutschler & Trede (2026), "Pruned Skewed Kalman Filter
##   and Smoother with Application to DSGE Models," JEDC Vol. 187 (Dynare WP
##   #78). Reference implementation: github.com/gguljanov/pruned-skewed-kalman.
##
## IMPLEMENTATION (method = "csn"):
##   Step 1 (forward): run .pskf_filter(store_path=TRUE) to store the full
##     per-period CSN state distributions
##     (mu_{t|t}, Sigma_{t|t}, Gamma_{t|t}, nu_{t|t}, Delta_{t|t})
##     and predicted distributions (mu_{t+1|t}, Sigma_{t+1|t}, ...).
##
##   Step 2 (backward): Rauch-Tung-Striebel (RTS) backward pass.
##     Gaussian part: standard RTS gives smoothed (mu_{t|T}, Sigma_{t|T}).
##     CSN part: backward update of nu (Guljanov et al. 2026, Section 3.3):
##       nu_{t|T} = nu_{t|t} - Gamma_{t|t} G_t (mu_{t+1|T} - mu_{t+1|t})
##       Gamma_{t|T} = Gamma_{t|t}   (unchanged -- only nu shifts backward)
##       Delta_{t|T} = Delta_{t|t}   (unchanged)
##     Initialise at t=T: (Gamma_{T|T}, nu_{T|T}, Delta_{T|T}) from the filter.
##
##   Step 3 (CSN mean correction): the true smoothed mean is
##       E[x_{t|T}] = mu_{t|T} + Sigma_{t|T} Gamma_{t|T}' h
##     where h is the q-dimensional normal hazard rate evaluated at -nu_{t|T}
##     w.r.t. the smoothed CSN covariance D_{t|T} = Delta_{t|T} + Gamma_{t|T} Sigma_{t|T} Gamma_{t|T}'.
##     The hazard rate h satisfies h_i = [D^{-1} phi_q(-nu; 0, D)]_i / Phi_q(-nu; 0, D).
##     For q=1: h = phi(nu/sqrt(D)) / (sqrt(D) * Phi(-nu/sqrt(D))) (Mills ratio).
##     For q>1: we use logcdf_ME_r (same as the filter) to compute Phi_q, and
##              finite differences on nu to approximate h_i numerically.
##
##   When Gamma is all-zero (pure Gaussian), the hazard rate is zero and the CSN
##   correction vanishes -- exact Gaussian RTS (method = "gaussian" also does this).
##
## METHOD = "gaussian" (legacy):
##   The v1 Option-B smoother: Gaussian RTS backward pass on the CSN
##   forward-pass moments WITHOUT propagating or applying the CSN skewness
##   correction. Exact at alpha=0; approximate at alpha != 0 (mean gap ~4.7e-3
##   for alpha=2). Kept for backward compatibility and as a cheap baseline.
##
## Entry point: pskf_smoother()
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## .csn_hazard_rate
##
## Compute the q-dimensional normal hazard rate vector:
##   h(-nu; 0, D) = D^{-1} phi_q(-nu; 0, D) / Phi_q(-nu; 0, D)
##
## This is the gradient of log Phi_q with respect to nu, evaluated at -nu.
## For q = 1: exact via Mills ratio.
## For q >= 2: numerical gradient via centered finite differences on nu,
##   using logcdf_ME_r for Phi_q evaluation.  Step size h = 1e-5.
##
## Returns: numeric(q) hazard rate vector.
#' @noRd
.csn_hazard_rate <- function(nu, D) {
  q <- length(nu)
  if (q == 0L) return(numeric(0))

  D <- as.matrix(D)
  nu <- as.numeric(nu)

  ## q = 1: exact Mills ratio
  if (q == 1L) {
    d_sd <- sqrt(max(D[1L, 1L], .Machine$double.eps))
    z    <- -nu[1L] / d_sd             # standardised argument
    lp   <- pnorm(z, log.p = TRUE)    # log Phi(-nu / sqrt(D))
    lph  <- dnorm(z, log = TRUE) - log(d_sd)  # log phi(-nu; 0, D) = log N(-nu; 0, D)
    ## h = phi / Phi  (both evaluated at -nu)
    return(exp(lph - lp))
  }

  ## q >= 2: numerical gradient of log Phi_q(-nu; 0, D) w.r.t. nu_i
  ## h_i = d/d(nu_i) log Phi_q(-nu; 0, D) = -d/d(x_i) log Phi_q(x; 0, D)|_{x=-nu}
  ## Use centered finite difference with h = 1e-5 * sqrt(D_ii)
  h_vec <- numeric(q)
  lp0   <- logcdf_ME_r(-nu, D)   # log Phi_q(-nu; 0, D)

  fd_step <- 1e-5
  for (i in seq_len(q)) {
    nu_p <- nu;  nu_p[i] <- nu_p[i] + fd_step
    nu_m <- nu;  nu_m[i] <- nu_m[i] - fd_step
    lp_p <- logcdf_ME_r(-nu_p, D)
    lp_m <- logcdf_ME_r(-nu_m, D)
    ## h_i = d/d(nu_i) log Phi_q(-nu; 0, D) = [lp_m - lp_p] / (2 * fd_step)
    ## because d/d(nu_i) (-nu) = -1 => d/dnu lp = -d/dx lp evaluated at x=-nu
    h_vec[i] <- (lp_m - lp_p) / (2 * fd_step)
  }
  h_vec
}


## ---------------------------------------------------------------------------
## pskf_smoother
##
## Run the PSKF forward pass (store_path=TRUE) and then apply the backward
## pass to produce smoothed state estimates.
##
## Two methods:
##   method = "csn"     -- CSN backward recursion (default); propagates the
##                         CSN nu parameter backward (first-order: Gamma/Delta
##                         held at filtered values) and applies the hazard-rate
##                         mean correction. Reduces to Gaussian RTS at alpha=0
##                         (exact, tolerance 1e-10); ~4.7e-4 mean gap vs the
##                         grid oracle at alpha=2.
##   method = "gaussian" -- Gaussian RTS backward pass only; no CSN skewness
##                         correction. Legacy v1 method; preserved for backward
##                         compatibility and as a cheap baseline.
##
## Arguments:
##   Y         n_obs x T matrix of observations (NA = missing).
##   TT        n_state x n_state transition matrix.
##   ZZ        n_obs x n_state observation matrix.
##   mu_eta    length-n_state mean-correction for state noise (E[eta]=0 after).
##   Sigma_eta n_state x n_state state noise covariance.
##   Gamma_eta q_eta x n_state skewness loading matrix for state noise.
##   nu_eta    length-q_eta skewness location for state noise.
##   Delta_eta q_eta x q_eta skewness covariance for state noise.
##   mu_eps    length-n_obs measurement noise mean.
##   Sigma_eps n_obs x n_obs measurement noise covariance.
##   cut_tol   Pruning tolerance passed to .pskf_filter (default 0.01).
##   method    "csn" (default) or "gaussian". See above.
##
## Returns: list(
##   smoothed_means  T x n_state matrix of smoothed state means.
##   smoothed_covs   n_state x n_state x T array of smoothed covariances.
##   filtered_means  T x n_state matrix of filtered state means.
##   filtered_covs   n_state x n_state x T array of filtered covariances.
##   loglik          scalar log p(y_{1:T}) from the PSKF forward pass.
## )
##
## Note: at T=1 there is no backward pass; smoothed == filtered.
##
## Backward compatibility: all existing positional/named calls with 11 arguments
##   (Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, mu_eps,
##    Sigma_eps, cut_tol) work unchanged and return the same fields plus
##   the CSN correction in smoothed_means.  Callers that assume the Gaussian
##   baseline can pass method = "gaussian".
## @noRd
pskf_smoother <- function(Y, TT, ZZ, mu_eta, Sigma_eta, Gamma_eta, nu_eta,
                           Delta_eta, mu_eps, Sigma_eps, cut_tol = 0.01,
                           method = c("csn", "gaussian")) {
  method <- match.arg(method)

  ## Ensure Y is a matrix (n_obs x T)
  if (!is.matrix(Y)) Y <- matrix(Y, nrow = 1L)
  n_obs   <- nrow(Y)
  n_T     <- ncol(Y)
  n_state <- nrow(TT)

  ## ---- Forward pass (with path storage) ------------------------------------
  fwd <- .pskf_filter(
    Y         = Y,
    TT        = TT,
    ZZ        = ZZ,
    mu_eta    = mu_eta,
    Sigma_eta = Sigma_eta,
    Gamma_eta = Gamma_eta,
    nu_eta    = nu_eta,
    Delta_eta = Delta_eta,
    mu_eps    = mu_eps,
    Sigma_eps = Sigma_eps,
    cut_tol   = cut_tol,
    store_path = TRUE
  )

  loglik          <- fwd$ll
  mu_pred_path    <- fwd$mu_pred_path
  Sigma_pred_path <- fwd$Sigma_pred_path
  mu_filt_path    <- fwd$mu_filt_path
  Sigma_filt_path <- fwd$Sigma_filt_path
  Gamma_filt_path <- fwd$Gamma_filt_path
  nu_filt_path    <- fwd$nu_filt_path
  Delta_filt_path <- fwd$Delta_filt_path

  ## ---- Gaussian RTS backward pass ------------------------------------------
  ## Allocate output arrays
  s_smooth <- matrix(0, n_T, n_state)   # smoothed means (T x n_state)
  P_smooth <- array(0, dim = c(n_state, n_state, n_T))

  ## Allocate CSN path arrays for smoothed skewness (method="csn")
  Gamma_smooth_path <- vector("list", n_T)
  nu_smooth_path    <- vector("list", n_T)
  Delta_smooth_path <- vector("list", n_T)

  ## Initialise smoother at t = T: smoothed == filtered
  s_smooth[n_T, ]  <- mu_filt_path[[n_T]]
  P_smooth[, , n_T] <- Sigma_filt_path[[n_T]]
  Gamma_smooth_path[[n_T]] <- Gamma_filt_path[[n_T]]
  nu_smooth_path[[n_T]]    <- nu_filt_path[[n_T]]
  Delta_smooth_path[[n_T]] <- Delta_filt_path[[n_T]]

  ## Backward sweep: t = T-1 down to 1 (empty when n_T == 1)
  for (step in seq_len(n_T - 1L)) {
    t <- n_T - step   # reverse order: T-1, T-2, ..., 1

    mu_f    <- mu_filt_path[[t]]
    Sigma_f <- Sigma_filt_path[[t]]
    mu_p    <- mu_pred_path[[t + 1L]]
    Sigma_p <- Sigma_pred_path[[t + 1L]]

    ## RTS gain: G_t = Sigma_f TT' Sigma_p^{-1}
    ## Solve for G' instead of inverting Sigma_p directly (more stable)
    Sigma_p_sym <- 0.5 * (Sigma_p + t(Sigma_p))   # enforce symmetry
    Sigma_p_reg <- Sigma_p_sym + 1e-10 * diag(n_state)  # small regularisation
    TT_Sf <- TT %*% Sigma_f                        # TT Sigma_f
    G_t <- tryCatch(
      t(solve(Sigma_p_reg, TT_Sf)),                # = Sigma_f TT' Sigma_p^{-1}
      error = function(e) {
        ## Fallback: pseudo-inverse via SVD
        sv <- svd(Sigma_p_sym)
        tol <- max(sv$d) * .Machine$double.eps * n_state
        d_inv <- ifelse(sv$d > tol, 1 / sv$d, 0)
        ## G_t = Sigma_f TT' (Sigma_p^{-1})
        Sigma_f %*% t(TT) %*% (sv$u %*% diag(d_inv, n_state) %*% t(sv$v))
      }
    )

    ## Gaussian smoother update (mean and covariance)
    s_smooth[t, ] <- mu_f + as.numeric(G_t %*% (s_smooth[t + 1L, ] - mu_p))
    P_s_next      <- P_smooth[, , t + 1L]
    diff_P        <- P_s_next - Sigma_p_sym
    P_smooth[, , t] <- Sigma_f + G_t %*% diff_P %*% t(G_t)

    ## Enforce symmetry of smoothed covariance
    P_smooth[, , t] <- 0.5 * (P_smooth[, , t] + t(P_smooth[, , t]))

    ## CSN backward update of nu (method = "csn"):
    ## nu_{t|T} = nu_{t|t} - Gamma_{t|t} G_t (mu_{t+1|T} - mu_{t+1|t})
    ## Gamma and Delta are unchanged (the shift is a mean-shift, not a new
    ## skewness direction; the CSN correction is via the new nu in the
    ## hazard rate calculation below).
    ##
    ## Reference: Guljanov, Mutschler & Trede (2026), Section 3.3, eqs (3.14)-(3.19).
    ## Physical intuition: the backward correction mu_{t+1|T} - mu_{t+1|t} adds
    ## information about x_{t+1} that shifts the effective nu just as the
    ## forward innovation does in the filter update (nu_{t|t} = nu_pred - G_skew v_t).
    Gamma_filt_t  <- Gamma_filt_path[[t]]
    nu_filt_t     <- nu_filt_path[[t]]
    Delta_filt_t  <- Delta_filt_path[[t]]

    q_filt <- nrow(Gamma_filt_t)
    if (method == "csn" && q_filt > 0L) {
      ## Backward mean correction vector: d_t = mu_{t+1|T}^{RTS} - mu_{t+1|t}
      d_t    <- s_smooth[t + 1L, ] - as.numeric(mu_p)
      ## nu shift: Gamma_{t|t} G_t d_t  (q_filt-vector)
      nu_shift <- as.numeric(Gamma_filt_t %*% (G_t %*% d_t))
      nu_smooth_path[[t]]    <- nu_filt_t - nu_shift
      Gamma_smooth_path[[t]] <- Gamma_filt_t
      Delta_smooth_path[[t]] <- Delta_filt_t
    } else {
      ## Pure Gaussian (alpha=0) or Gaussian method: no CSN correction
      nu_smooth_path[[t]]    <- nu_filt_t
      Gamma_smooth_path[[t]] <- Gamma_filt_t
      Delta_smooth_path[[t]] <- Delta_filt_t
    }
  }

  ## ---- Apply CSN mean correction (method = "csn") --------------------------
  ## The Gaussian RTS smoothed mean mu_{t|T} is just the Gaussian part.
  ## The true smoothed mean of the CSN distribution is:
  ##   E[x_{t|T}] = mu_{t|T} + Sigma_{t|T} Gamma_{t|T}' h(-nu_{t|T}; 0, D_{t|T})
  ## where D_{t|T} = Delta_{t|T} + Gamma_{t|T} Sigma_{t|T} Gamma_{t|T}'  and
  ## h is the q-dimensional normal hazard rate vector.

  if (method == "csn") {
    for (t in seq_len(n_T)) {
      Gamma_s <- Gamma_smooth_path[[t]]
      q_s     <- nrow(Gamma_s)

      if (q_s == 0L) next   # pure Gaussian: no correction

      nu_s    <- nu_smooth_path[[t]]
      Delta_s <- Delta_smooth_path[[t]]
      P_s     <- P_smooth[, , t]
      P_s_sym <- 0.5 * (P_s + t(P_s))

      ## D_{t|T} = Delta_{t|T} + Gamma_{t|T} Sigma_{t|T} Gamma_{t|T}'
      D_s     <- Delta_s + Gamma_s %*% P_s_sym %*% t(Gamma_s)
      D_s     <- 0.5 * (D_s + t(D_s))   # enforce symmetry

      ## Hazard rate vector h = h(-nu_s; 0, D_s)
      h_vec <- tryCatch(
        .csn_hazard_rate(nu_s, D_s),
        error = function(e) rep(0, q_s)
      )

      ## CSN mean correction: Sigma_{t|T} Gamma_{t|T}' h
      correction <- as.numeric(P_s_sym %*% t(Gamma_s) %*% h_vec)
      s_smooth[t, ] <- s_smooth[t, ] + correction
    }
  }

  ## Collect filtered means / covariances into tidy arrays for the caller.
  ## For method="csn": apply the CSN hazard-rate correction to the filtered means
  ## too, so filtered_means[t, ] = E[x_t | Y_{1:t}] (true posterior mean).
  ## This ensures smoothed_means == filtered_means at T=1 (invariant: the
  ## T=1 smoothed distribution IS the filtered distribution).
  s_filt <- matrix(0, n_T, n_state)
  P_filt <- array(0, dim = c(n_state, n_state, n_T))
  for (t in seq_len(n_T)) {
    s_filt[t, ]   <- mu_filt_path[[t]]
    P_filt[, , t] <- Sigma_filt_path[[t]]
  }
  if (method == "csn") {
    for (t in seq_len(n_T)) {
      Gamma_f <- Gamma_filt_path[[t]]
      q_f     <- nrow(Gamma_f)
      if (q_f == 0L) next
      nu_f    <- nu_filt_path[[t]]
      Delta_f <- Delta_filt_path[[t]]
      P_f     <- P_filt[, , t]
      P_f_sym <- 0.5 * (P_f + t(P_f))
      D_f     <- Delta_f + Gamma_f %*% P_f_sym %*% t(Gamma_f)
      D_f     <- 0.5 * (D_f + t(D_f))
      h_vec_f <- tryCatch(
        .csn_hazard_rate(nu_f, D_f),
        error = function(e) rep(0, q_f)
      )
      s_filt[t, ] <- s_filt[t, ] +
        as.numeric(P_f_sym %*% t(Gamma_f) %*% h_vec_f)
    }
  }

  list(
    smoothed_means  = s_smooth,
    smoothed_covs   = P_smooth,
    filtered_means  = s_filt,
    filtered_covs   = P_filt,
    loglik          = loglik
  )
}
