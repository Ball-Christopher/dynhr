### ===========================================================================
### dynhr_smoother.R -- Kalman smoother + historical decomposition for dynhr
### ===========================================================================
### Provides:
###   build_dsge_state_space()  -- extract compact state-space from m + dr
###   kalman_smoother()         -- RTS backward smoother for DSGE form
###   historical_decomposition() -- per-shock contribution to all endo vars
###
### DSGE state-space form (first-order perturbation):
###   s_t = T s_{t-1} + R eps_t          (n_state x 1)
###   y_t = Z s_{t-1} + D eps_t          (n_obs x 1)
###
### where s_t = [backward; mixed] variables, y_t = observables.
### Note: y_t depends on s_{t-1} and eps_t shares between both equations
### (correlated state/observation noise).
###
### For historical decomposition, all endo vars are recovered via:
###   yall_t = ghx * s_{t-1} + ghu * eps_t
### ===========================================================================

# ---------------------------------------------------------------------------
#' Build compact DSGE state-space from model and decision rules
#'
#' Extracts T, R, Z, D matrices from ghx/ghu using lead_lag_incidence
#' to determine the correct ghx column ordering.
#'
#' State-space form:
#'   \eqn{s_t   = T s_{t-1} + R \varepsilon_t}  (n_state x 1)
#'   \eqn{yall_t = ghx \cdot s_{t-1} + ghu \cdot \varepsilon_t}  (n_endo x 1, all vars)
#'   \eqn{y_t   = Z s_{t-1} + D \varepsilon_t}  (n_obs x 1, observables only)
#'
#' @param m         Parsed model object (from parse_mod)
#' @param dr        Decision rules (from stoch_simul()$dr)
#' @param obs_names Character vector of observable names
#' @param verbose   Print diagnostic info (default TRUE)
#' @param params    Named parameter vector used to evaluate the shock
#'   covariance \code{Sigma_e} from the \code{shocks;} block (default:
#'   \code{m$param_values}). Pass the draw-specific vector when \code{dr}
#'   was solved at non-default parameters.
#' @return List with T_mat, R_mat, Z_mat, D_mat, Sigma_e, ghx, ghu,
#'   indices, names
#' @export
# ---------------------------------------------------------------------------
build_dsge_state_space <- function(m, dr, obs_names, verbose = TRUE,
                                   params = m$param_values) {
  
  endo_names <- m$var_names
  n_endo     <- length(endo_names)
  n_state    <- ncol(dr$ghx)
  n_shock    <- ncol(dr$ghu)
  
  stopifnot(nrow(dr$ghx) == n_endo,
            nrow(dr$ghu) == n_endo)
  
  ## ---- Step 1: Determine ghx column ordering from lead_lag_incidence ----
  ## lli[1, j] > 0 means variable j appears at t-1 (is a state variable)
  ## The VALUE of lli[1, j] gives the Jacobian column index, which
  ## determines the ordering of ghx columns.
  
  lli <- m$lead_lag_incidence
  ghx_col_to_endo <- NULL
  
  if (!is.null(lli) && is.matrix(lli) && nrow(lli) >= 1) {
    lag_row <- lli[1, ]
    has_lag <- which(lag_row > 0)
    
    if (length(has_lag) == n_state) {
      ## Sort by Jacobian column index to get ghx column ordering
      ghx_col_to_endo <- has_lag[order(lag_row[has_lag])]
      
      if (verbose) {
        state_names_ordered <- endo_names[ghx_col_to_endo]
        cat(sprintf("  ghx column order (from lli): %s\n",
                    paste(state_names_ordered, collapse = ", ")))
      }
    } else {
      warning(sprintf(
        "lead_lag_incidence lag count (%d) != ghx columns (%d). Falling back.",
        length(has_lag), n_state
      ))
    }
  }
  
  ## ---- Fallback: try declaration order of state variables ----
  if (is.null(ghx_col_to_endo)) {
    vc <- m$variable_classification
    back_vars  <- character(0)
    mixed_vars <- character(0)
    
    if (is.list(vc) && !is.data.frame(vc)) {
      if (!is.null(vc$predetermined)) back_vars  <- vc$predetermined
      if (!is.null(vc$backward))      back_vars  <- c(back_vars, vc$backward)
      if (!is.null(vc$mixed))         mixed_vars <- vc$mixed
    }
    
    state_vars <- c(back_vars, mixed_vars)
    state_idx_unsorted <- match(state_vars, endo_names)
    
    ## Try declaration order first (most common for dynhr)
    decl_order <- sort(state_idx_unsorted)
    pred_mixed_order <- state_idx_unsorted   # [pred, mixed] Dynare convention
    
    ## Test both: whichever gives stable T eigenvalues is correct
    T_decl <- dr$ghx[decl_order, , drop = FALSE]
    T_pm   <- dr$ghx[pred_mixed_order, , drop = FALSE]
    
    eig_decl <- max(abs(eigen(T_decl, only.values = TRUE)$values))
    eig_pm   <- max(abs(eigen(T_pm,   only.values = TRUE)$values))
    
    if (eig_decl < 1.0) {
      ghx_col_to_endo <- decl_order
      if (verbose) cat("  ghx column order: declaration order (verified by eigenvalues)\n")
    } else if (eig_pm < 1.0) {
      ghx_col_to_endo <- pred_mixed_order
      if (verbose) cat("  ghx column order: [predetermined, mixed] (verified by eigenvalues)\n")
    } else {
      warning(sprintf(
        "Neither ordering gives stable T. max|eig|: decl=%.4f, pred_mixed=%.4f",
        eig_decl, eig_pm
      ))
      ghx_col_to_endo <- decl_order  # best guess
    }
  }
  
  ## ---- Step 2: Build compact state-space matrices ----
  ## T_mat: state transition (rows = state vars in ghx column order)
  ## R_mat: shock impact on states
  T_mat <- dr$ghx[ghx_col_to_endo, , drop = FALSE]   # n_state x n_state
  R_mat <- dr$ghu[ghx_col_to_endo, , drop = FALSE]    # n_state x n_shock
  
  ## Z_mat: observation equation (obs rows of ghx)
  ## D_mat: direct shock -> obs
  obs_idx <- match(obs_names, endo_names)
  if (any(is.na(obs_idx))) {
    stop("Observables not found: ",
         paste(obs_names[is.na(obs_idx)], collapse = ", "))
  }
  Z_mat <- dr$ghx[obs_idx, , drop = FALSE]             # n_obs x n_state
  D_mat <- dr$ghu[obs_idx, , drop = FALSE]              # n_obs x n_shock
  
  ## ---- Step 3: Verify eigenvalue stability ----
  eig_mod <- abs(eigen(T_mat, only.values = TRUE)$values)
  state_names_ordered <- endo_names[ghx_col_to_endo]
  
  if (verbose) {
    cat(sprintf("  State-space: %d states, %d obs, %d shocks\n",
                n_state, length(obs_names), n_shock))
    cat(sprintf("  State vars: %s\n", paste(state_names_ordered, collapse = ", ")))
    cat(sprintf("  T eigenvalues: [%.4f, %.4f]",
                min(eig_mod), max(eig_mod)))
    if (all(eig_mod < 1.0)) cat(" -- all stable [OK]\n")
    else cat(sprintf(" -- UNSTABLE (max=%.4f) [X]\n", max(eig_mod)))
  }
  
  if (any(eig_mod >= 1.0)) {
    warning(sprintf(
      "T matrix has unstable eigenvalues (max |lambda| = %.4f).",
      max(eig_mod)
    ))
  }
  
  ## Shock covariance from the shocks; block. ghx/ghu are unit-shock
  ## responses (Sigma_e is NOT baked into them), so downstream filters
  ## must use this as Q -- kalman_smoother() defaults to it.
  Sigma_e <- .get_shock_cov(m, m$varexo_names, params)

  structure(
    list(
      T_mat             = T_mat,
      R_mat             = R_mat,
      Z_mat             = Z_mat,
      D_mat             = D_mat,
      Sigma_e           = Sigma_e,
      ghx               = dr$ghx,
      ghu               = dr$ghu,
      ghx_col_to_endo   = ghx_col_to_endo,     # maps ghx column i -> endo index
      obs_idx           = obs_idx,
      state_names       = state_names_ordered,  # in ghx column order
      obs_names         = obs_names,            # observable names (character)
      endo_names        = endo_names,
      shock_names       = m$varexo_names,
      n_state           = n_state,
      n_obs             = length(obs_names),
      n_shock           = n_shock,
      n_endo            = n_endo,
      timing            = "lagged"              # Convention A: y_t = Z s_{t-1} + D eps_t
    ),
    class = "dsge_ss"
  )
}

# ---------------------------------------------------------------------------
#' Kalman smoother for DSGE state-space
#'
#' DSGE form (note timing: \eqn{y_t} depends on \eqn{s_{t-1}}):
#' \preformatted{
#'   s_t = T s_{t-1} + R eps_t
#'   y_t = Z s_{t-1} + D eps_t
#' }
#'
#' State and observation noise are correlated (shared eps_t):
#'   Var(R eps) = R Q R'
#'   Var(D eps) = D Q D'
#'   Cov(R eps, D eps) = R Q D'
#'
#' Forward pass: standard KF with correlated noise.
#' Backward pass: Rauch-Tung-Striebel smoother.
#'
#' @param Y      T x n_obs data matrix
#' @param ss     State-space list from build_dsge_state_space()
#' @param Q      n_shock x n_shock shock covariance. Default \code{NULL}:
#'   use \code{ss$Sigma_e} (the covariance from the \code{shocks;} block),
#'   which makes the forward-pass \code{loglik} identical to
#'   \code{kalman_filter()} on the same data. \code{ghu} holds unit-shock
#'   responses, so the identity is only correct when every shock has
#'   \code{stderr 1}; falls back to the identity (with a warning) only for
#'   hand-built \code{ss} lists lacking \code{Sigma_e}.
#' @param me_extra  \code{n_obs x T} matrix of per-period additive
#'   measurement-error variances (from a \code{filter_tunes} block), or
#'   \code{NULL} (no filter tunes).
#' @param shock_scale  \code{n_exo x T} matrix of per-period shock
#'   standard-deviation scale factors (from a \code{heteroskedastic_shocks}
#'   block), or \code{NULL} (constant shock variances).
#' @return List with smoothed_states, smoothed_shocks, filtered_states, loglik
#' @export
# ---------------------------------------------------------------------------
kalman_smoother <- function(Y, ss, Q = NULL, me_extra = NULL,
                            shock_scale = NULL) {

  ## Convert current-state dsge_ss to lagged-state before extracting matrices.
  ## ss_convert_timing() is a no-op when ss$timing == "lagged".
  if (inherits(ss, "dsge_ss") && !is.null(ss$timing) && ss$timing != "lagged")
    ss <- ss_convert_timing(ss)

  TT      <- nrow(Y)
  n_s     <- ss$n_state
  n_obs   <- ss$n_obs
  n_shk   <- ss$n_shock

  TT_mat  <- ss$T_mat
  R_mat   <- ss$R_mat
  Z_mat   <- ss$Z_mat
  D_mat   <- ss$D_mat

  if (is.null(Q)) {
    if (!is.null(ss$Sigma_e)) {
      Q <- ss$Sigma_e
    } else {
      warning("kalman_smoother: ss has no Sigma_e (hand-built list?); ",
              "using Q = identity, which assumes all shocks have stderr 1.",
              call. = FALSE)
      Q <- diag(n_shk)
    }
  }

  ## ---- me_extra validation ------------------------------------------------
  ## me_extra (n_obs x T) holds per-period per-observable extra ME variances
  ## (filter_tunes soft tunes: stderr^2 at tune periods, 0 elsewhere).
  if (!is.null(me_extra)) {
    if (!is.matrix(me_extra) || nrow(me_extra) != n_obs || ncol(me_extra) != TT)
      stop(sprintf(
        "kalman_smoother: me_extra must be n_obs x T (%d x %d); got %s.",
        n_obs, TT,
        if (is.matrix(me_extra)) paste0(nrow(me_extra), " x ", ncol(me_extra))
        else "non-matrix"), call. = FALSE)
  }
  has_me_extra <- !is.null(me_extra) && any(me_extra != 0)

  ## ---- shock_scale validation ----------------------------------------------
  ## shock_scale (n_shk x T) holds per-period shock std scale factors.
  if (!is.null(shock_scale)) {
    if (!is.matrix(shock_scale) || nrow(shock_scale) != n_shk || ncol(shock_scale) != TT)
      stop(sprintf(
        "kalman_smoother: shock_scale must be n_shk x T (%d x %d); got %s.",
        n_shk, TT,
        if (is.matrix(shock_scale)) paste0(nrow(shock_scale), " x ", ncol(shock_scale))
        else "non-matrix"), call. = FALSE)
  }
  has_shock_scale <- !is.null(shock_scale) && !all(shock_scale == 1)

  ## Pre-compute noise covariances (baseline, used for P0 and when not scaling)
  RQR <- R_mat %*% Q %*% t(R_mat)       # state noise cov
  DQD <- D_mat %*% Q %*% t(D_mat)       # obs noise cov
  RQD <- R_mat %*% Q %*% t(D_mat)       # cross-covariance

  ## ---- Initialisation: unconditional state covariance (Lyapunov) ----
  ## Use the canonical solve_lyapunov() from stochsimul-monolith.R;
  ## fall back to large diagonal for near-unit-root / nonstationary models
  ## (solve_lyapunov returns a NaN matrix when the doubling algorithm diverges
  ## and the vec-Lyapunov system is singular).
  P_ss <- solve_lyapunov(TT_mat, RQR)
  if (anyNA(P_ss)) {
    ## Unit roots detected: use a diffuse (large-diagonal) prior so the
    ## smoother does not crash.  The smoothed states will be valid but the
    ## loglik has a kappa-dependent additive offset (not suitable for
    ## cross-method comparison; use kalman_filter(lik_init="diffuse") for
    ## exact diffuse likelihood evaluation).
    warning("kalman_smoother: unit root(s) detected in TT -- ",
            "solve_lyapunov() returned NaN. ",
            "Falling back to diffuse prior P0 = ", .DIFFUSE_SCALE,
            " * I(", n_s, "). ",
            "Smoothed states are valid; loglik has a kappa-dependent offset. ",
            "Use kalman_filter(lik_init=\"diffuse\") for exact diffuse loglik.",
            call. = FALSE)
    P_ss <- .DIFFUSE_SCALE * diag(n_s)
  }

  ## ---- Forward pass (Kalman filter) ----
  ## Per-period NA observations are handled by dropping the NA rows from
  ## Z_mat / DQD / RQD / the innovation vector.  An all-NA period is treated
  ## as predict-only (no update).  This mirrors the standard filter's NA
  ## branch in R/kalman-filter.R:1193-1218 so the forward quantities fed into
  ## the RTS backward pass are always consistent.
  ##
  ## DK disturbance smoother requires: v_t (innovations), F_inv_t (inverse
  ## innovation covariance), K_t (Kalman gain), and Zt/Dt (obs/shock matrices
  ## subsetting to non-NA rows).  These are stored in lists indexed by period.
  s_filt <- matrix(0, TT, n_s)         # s_{t|t}
  s_pred <- matrix(0, TT, n_s)         # s_{t|t-1}
  P_filt <- array(0, dim = c(n_s, n_s, TT))
  P_pred <- array(0, dim = c(n_s, n_s, TT))
  ## Number of non-NA observations per period (for the loglik correction term).
  n_obs_t_vec <- integer(TT)

  ## Storage for DK disturbance smoother backward pass.
  dk_v     <- vector("list", TT)   # v_t (n_ok x 1)
  dk_Finv  <- vector("list", TT)   # F_t^{-1} (n_ok x n_ok)
  dk_K     <- vector("list", TT)   # K_t (n_s x n_ok)
  dk_Z     <- vector("list", TT)   # Zt  (n_ok x n_s)
  dk_D     <- vector("list", TT)   # Dt  (n_ok x n_shk)
  dk_Q     <- vector("list", TT)   # Q_t (n_shk x n_shk) -- needed when shock_scale active

  s_tt <- rep(0, n_s)
  P_tt <- P_ss
  loglik <- 0

  for (t in seq_len(TT)) {
    ## Per-period scaled shock covariance (when shock_scale is active).
    Q_t   <- if (has_shock_scale) {
      sc_t <- shock_scale[, t]
      Q * outer(sc_t, sc_t)
    } else Q
    RQR_t <- if (has_shock_scale) R_mat %*% Q_t %*% t(R_mat) else RQR
    DQD_t <- if (has_shock_scale) D_mat %*% Q_t %*% t(D_mat) else DQD
    RQD_t <- if (has_shock_scale) R_mat %*% Q_t %*% t(D_mat) else RQD

    ## ---- Predict ----
    s_tp <- as.numeric(TT_mat %*% s_tt)
    P_tp <- TT_mat %*% P_tt %*% t(TT_mat) + RQR_t

    s_pred[t, ]   <- s_tp
    P_pred[, , t] <- P_tp

    ## ---- Innovation (with per-period NA handling) ----
    y_pred <- as.numeric(Z_mat %*% s_tt)
    v_t    <- Y[t, ] - y_pred

    obs_ok <- which(!is.na(v_t))
    n_ok   <- length(obs_ok)
    n_obs_t_vec[t] <- n_ok

    if (n_ok == 0L) {
      ## All observables missing: predict-only, no update.
      s_tt <- s_tp
      P_tt <- P_tp
      s_filt[t, ]   <- s_tt
      P_filt[, , t] <- P_tt
      ## DK: no observation this period -- v, Finv, K left NULL; Q stored.
      dk_Q[[t]] <- Q_t
      next
    }

    ## Subset to non-NA observables.
    if (n_ok < n_obs) {
      v_t   <- v_t[obs_ok]
      Zt    <- Z_mat[obs_ok, , drop = FALSE]
      Dt    <- D_mat[obs_ok, , drop = FALSE]
      DQDt  <- tcrossprod(Dt %*% Q_t, Dt)
      RQDt  <- RQD_t[, obs_ok, drop = FALSE]
    } else {
      Zt   <- Z_mat
      Dt   <- D_mat
      DQDt <- DQD_t
      RQDt <- RQD_t
    }

    ## Innovation covariance: F_t = Z_t P_{t-1|t-1} Z_t' + DQD_t [+ me_extra_t]
    F_t  <- Zt %*% P_tt %*% t(Zt) + DQDt
    if (has_me_extra) diag(F_t) <- diag(F_t) + me_extra[obs_ok, t]
    F_t  <- (F_t + t(F_t)) * 0.5

    ## base chol() THROWS on a non-PD matrix (it never returns NULL), so the
    ## jitter fallback below must catch the error for the regularisation to run.
    F_ch <- tryCatch(chol(F_t), error = function(e) NULL)
    if (is.null(F_ch)) {
      ## F_t not PD (e.g. stochastic singularity): regularise with increasing jitter.
      for (jit in c(1e-8, 1e-6, 1e-4, 1e-2)) {
        F_ch <- tryCatch(chol(F_t + jit * diag(n_ok)), error = function(e) NULL)
        if (!is.null(F_ch)) { F_t <- F_t + jit * diag(n_ok); break }
      }
      if (is.null(F_ch)) F_ch <- chol(F_t + 0.1 * diag(n_ok))
    }
    F_inv     <- chol2inv(F_ch)
    log_det_F <- 2 * sum(log(diag(F_ch)))

    ## Kalman gain: K = (T P_{t-1|t-1} Z_t' + RQD_t) F_t^{-1}
    K_t <- (TT_mat %*% P_tt %*% t(Zt) + RQDt) %*% F_inv

    ## Store forward-pass quantities for DK disturbance smoother.
    dk_v[[t]]    <- v_t
    dk_Finv[[t]] <- F_inv
    dk_K[[t]]    <- K_t
    dk_Z[[t]]    <- Zt
    dk_D[[t]]    <- Dt
    dk_Q[[t]]    <- Q_t

    ## Updated state: s_{t|t} = s_{t|t-1} + K (y_t - Z_t s_{t-1|t-1})
    s_tt <- s_tp + as.numeric(K_t %*% v_t)
    P_tt <- P_tp - K_t %*% F_t %*% t(K_t)
    P_tt <- 0.5 * (P_tt + t(P_tt))

    s_filt[t, ]   <- s_tt
    P_filt[, , t] <- P_tt

    ## Log-likelihood contribution (adjusts constant for n_ok != n_obs).
    loglik <- loglik - 0.5 * (n_ok * log(2 * pi) + log_det_F +
                                as.numeric(t(v_t) %*% F_inv %*% v_t))
  }

  ## ---- Backward pass (RTS smoother) ----
  ## Means: s_{t|T} = s_{t|t} + J_t (s_{t+1|T} - s_{t+1|t}).
  ## Covariances: V_{t|T} = P_{t|t} + J_t (V_{t+1|T} - P_{t+1|t}) J_t'
  ## (Rauch-Tung-Striebel; J_t = P_{t|t} T' P_{t+1|t}^{-1} is the same gain used
  ## for the mean). V_{T|T} = P_{T|T}. Smoothing never increases uncertainty, so
  ## diag(V_{t|T}) <= diag(P_{t|t}) by construction.
  s_smooth <- matrix(0, TT, n_s)
  s_smooth[TT, ] <- s_filt[TT, ]
  V_smooth <- array(0, dim = c(n_s, n_s, TT))
  V_smooth[, , TT] <- P_filt[, , TT]

  for (t in (TT - 1L):1L) {
    Pp1     <- P_pred[, , t + 1L]
    Pp1_inv <- MASS::ginv(Pp1)
    J_t <- P_filt[, , t] %*% t(TT_mat) %*% Pp1_inv
    s_smooth[t, ] <- s_filt[t, ] +
      as.numeric(J_t %*% (s_smooth[t + 1L, ] - s_pred[t + 1L, ]))
    Vt <- P_filt[, , t] + J_t %*% (V_smooth[, , t + 1L] - Pp1) %*% t(J_t)
    V_smooth[, , t] <- 0.5 * (Vt + t(Vt))
  }

  ## ---- DK disturbance smoother: recover smoothed structural shocks ----
  ## Uses the Durbin-Koopman (2012) adjoint backward recursion to compute
  ## eps_{t|T} directly from the backward adjoint, matching Dynare's
  ## calib_smoother convention at ALL t including t=1.
  ##
  ## State-space (lagged form, DK §4.4):
  ##   s_t   = T s_{t-1} + R eps_t       -> alpha_t = s_{t-1}; eta_t = R eps_t
  ##   y_t   = Z s_{t-1} + D eps_t       -> obs noise = D eps_t; cross-cov = RQD'
  ##
  ## DK adjoint backward recursion (initialised r_T = 0):
  ##   r_{t-1} = Zt' Finv_t v_t + (T - K_t Zt)' r_t
  ##   (for periods with no obs: r_{t-1} = T' r_t)
  ##
  ## Smoothed structural shock (DK §4.4 eqns 4.44-4.47, correlated-noise form):
  ##
  ##   eps_{t|T} = Q_t (R' r_t + Dt' u_t)
  ##
  ## where r_t is the adjoint BEFORE the backward step at period t (i.e., the
  ## value coming in from period t+1), and u_t = Finv_t v_t - K_t' r_{t-1}.
  ##
  ## NOTE on timing: in DK's alpha_t = s_{t-1} indexing, the shock eps_t
  ## drives alpha_{t+1} = T alpha_t + R eps_t and enters y_t = Z alpha_t + D eps_t.
  ## The smoothed shock is hat_eta_t = Q_eta r_t (DK 4.44) with Q_eta = RQR'.
  ## Since hat_eta_t = R eps_{t|T}, this gives eps_{t|T} = Q R' r_t (D=0 case).
  ## For D != 0, the correlated-noise correction adds Q Dt' u_t (DK §4.4).
  ##
  ## For t>=2 this is numerically identical to the previous RTS+pseudo-inverse
  ## approach (both are minimum-MSE). At t=1 the DK formula uses r_1 (the
  ## backward adjoint from t=1 forward step) to get eps_{1|T} = Q (R' r_1 + D' u_1),
  ## whereas the previous code used J_0 to compute s_{0|T} residually; these
  ## differ only in the initial-condition treatment (both internally consistent).
  ## This matches Dynare's etahat(:,1) = Q*R'*r(:,1) convention.
  eps_smooth <- matrix(0, TT, n_shk)
  r_t <- rep(0, n_s)   # r_{T} = 0 (terminal condition)

  for (t in TT:1L) {
    Q_t <- dk_Q[[t]]
    ## eps_{t|T} = Q_t (R' r_t + Dt' u_t) -- uses r_t (BEFORE backward step at t)
    if (is.null(dk_v[[t]])) {
      ## All observables missing at period t.
      ## u_t doesn't exist; eps_{t|T} = Q_t R' r_t (only state adjoint, no D term)
      eps_smooth[t, ] <- as.numeric(Q_t %*% (t(R_mat) %*% r_t))
      ## Backward step: r_{t-1} = T' r_t
      r_t <- as.numeric(t(TT_mat) %*% r_t)
    } else {
      v_t    <- dk_v[[t]]
      Finv   <- dk_Finv[[t]]
      K_t_dk <- dk_K[[t]]
      Zt     <- dk_Z[[t]]
      Dt     <- dk_D[[t]]
      ## u_t = Finv v_t - K_t' r_t  (uses r_t BEFORE the backward step)
      ## NOTE: for dynhr's lagged-state form y_t = Z s_{t-1} + D eps_t, the
      ## correct u_t uses r_t (not r_{t-1} as in the standard DK §4.3 form for
      ## y_t = Z alpha_t). See derivation: for t=T, r_T=0 gives u_T = F_T^{-1}v_T,
      ## and eps_{T|T} = Q D' F_T^{-1} v_T -- the terminal shock is identified
      ## purely from the observation (correct when D != 0).
      u_t     <- as.numeric(Finv %*% v_t - t(K_t_dk) %*% r_t)
      ## eps_{t|T} = Q_t (R' r_t + Dt' u_t)  -- both use r_t (BEFORE backward step)
      eps_smooth[t, ] <- as.numeric(Q_t %*% (t(R_mat) %*% r_t + t(Dt) %*% u_t))
      ## Backward step: r_{t-1} = Zt' Finv v_t + (T - K_t Zt)' r_t
      Lt    <- TT_mat - K_t_dk %*% Zt           # L_t = T - K_t Z_t  (n_s x n_s)
      r_t   <- as.numeric(t(Zt) %*% (Finv %*% v_t) + t(Lt) %*% r_t)
    }
  }

  colnames(s_smooth)   <- ss$state_names
  colnames(s_filt)     <- ss$state_names
  colnames(eps_smooth) <- ss$shock_names

  dimnames(P_filt)   <- list(ss$state_names, ss$state_names, NULL)
  dimnames(P_pred)   <- list(ss$state_names, ss$state_names, NULL)
  dimnames(V_smooth) <- list(ss$state_names, ss$state_names, NULL)

  list(
    smoothed_states = s_smooth,
    smoothed_shocks = eps_smooth,
    filtered_states = s_filt,
    ## Per-period state covariances (n_state x n_state x T):
    filtered_cov    = P_filt,            # P_{t|t}
    predicted_cov   = P_pred,            # P_{t|t-1}
    smoothed_cov    = V_smooth,          # P_{t|T} (RTS)
    P_filt_last     = P_filt[, , TT],    # P_{T|T} (kept for back-compat)
    loglik          = loglik
  )
}


# ---------------------------------------------------------------------------
#' Historical decomposition: per-shock contributions to ALL endo variables
#'
#' Uses the full ghx/ghu (not just the compact state-space) to recover
#' contributions to all n_endo variables, not just states.
#'
#' For each shock j:
#'   \deqn{s_t^{(j)} = T s_{t-1}^{(j)} + R[:,j] \epsilon_{j,t}}
#'   \deqn{y_t^{(j)} = ghx \cdot s_{t-1}^{(j)} + ghu[:,j] \epsilon_{j,t}}
#'
#' @param smoothed_shocks  T x n_shock matrix
#' @param ss               State-space list from build_dsge_state_space()
#' @return List with $contributions (named list of T x n_endo matrices)
#'         and $total (T x n_endo)
#' @export
# ---------------------------------------------------------------------------
historical_decomposition <- function(smoothed_shocks, ss) {
  
  TT    <- nrow(smoothed_shocks)
  n_shk <- ss$n_shock
  n_s   <- ss$n_state
  n_end <- ss$n_endo
  
  TT_mat <- ss$T_mat       # n_state x n_state
  R_mat  <- ss$R_mat        # n_state x n_shock
  ghx    <- ss$ghx           # n_endo x n_state
  ghu    <- ss$ghu           # n_endo x n_shock
  
  contributions <- setNames(
    lapply(seq_len(n_shk), function(j) matrix(0, TT, n_end)),
    ss$shock_names
  )
  
  for (j in seq_len(n_shk)) {
    s_j <- rep(0, n_s)        # state attributable to shock j
    r_j <- R_mat[, j]          # state impact column
    g_j <- ghu[, j]            # full endo impact column
    
    for (t in seq_len(TT)) {
      eps_jt <- smoothed_shocks[t, j]
      if (is.na(eps_jt)) eps_jt <- 0
      
      ## All endo vars at t from shock j:
      ## y_t^(j) = ghx * s_{t-1}^(j) + ghu[:,j] * eps_jt
      contributions[[j]][t, ] <- as.numeric(ghx %*% s_j + g_j * eps_jt)
      
      ## State transition for shock j:
      s_j <- as.numeric(TT_mat %*% s_j + r_j * eps_jt)
    }
    colnames(contributions[[j]]) <- ss$endo_names
  }
  
  total <- Reduce(`+`, contributions)
  colnames(total) <- ss$endo_names
  
  list(contributions = contributions, total = total)
}
