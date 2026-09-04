## R/obc-filter.R
## --------------------------------------------------------------------------
## OBC regime-switching Kalman filter implementations.
##
## Provides:
##   kalman_filter_obc()     -- legacy filter consuming a pre-computed regime_path
##   pkf_extract_shock()     -- extract eps_{t|t} from KF quantities (GPR eq. 6)
##   pkf_backward_one_step() -- one-step backward smoother (GPR eq. 7)
##   pkf_check_binding()     -- shock-aware regime check replacing obc_should_bind
##   kalman_filter_obc_pkf() -- PKF with per-period inner convergence (exported)
##
## REGIME-SWITCHING STATE SPACE
##   State eq:  s_t = TT_r * s_{t-1} + RR_r * eps_t + c_state_r
##   Obs eq:    y_t = ZZ_r * s_{t-1} + DD_r * eps_t + d + c_obs_r
## where r = regime_path[t] is an integer bitfield (0 = all slack).
##
## kalman_filter_obc(): legacy outer-loop approach. Consumes a regime_path
##   built by obc_guess_verify() and calls the filter once.
##
## kalman_filter_obc_pkf(): PKF approach (Giovannini, Pfeiffer, Ratto 2021).
##   At each period t, runs an inner guess-and-verify loop that extracts the
##   shock eps_{t|t} (eq. 6) and backward-smooths the state s_{t-1|t} (eq. 7)
##   to check OBC regime consistency — no separate outer pass needed.
##
## kf_store completeness invariant (enforced by two fixes in 5fd8b64):
##   For every t in 1:n_T, kf_store$L[[t]] and kf_store$RR[[t]] are matrices.
##   - Bug A (max_inner exhaustion): if the inner for(j) loop uses all max_inner
##     iterations without convergence, a no-obs prediction fallback is stored
##     rather than leaving kf_store[[t]] as NULL.
##   - Bug B (list-element removal): NULL assignments use kf_store$v[t] <-
##     list(NULL)  (single-bracket), NOT kf_store$v[[t]] <- NULL  (which would
##     remove the element and shift subsequent indices).
## --------------------------------------------------------------------------


# =============================================================================
# Legacy regime-switching Kalman filter (outer-loop approach)
# =============================================================================

#' Regime-switching Kalman filter for OBC models
#'
#' State equation (regime_path[t] == 0, slack):
#'   s_t = TT_s * s_{t-1} + RR_s * eps_t
#'
#' State equation (regime_path[t] != 0, any subset binding):
#'   s_t = TT_r * s_{t-1} + RR_r * eps_t + c_state_r
#'
#' Observation equation (slack):
#'   y_t = ZZ_s * s_{t-1} + DD_s * eps_t + d
#'
#' Observation equation (non-slack):
#'   y_t = ZZ_r * s_{t-1} + DD_r * eps_t + d + c_obs_r
#'
#' P is always updated via the Joseph form; no steady-state shortcut is
#' applied because the gain changes at every regime transition.
#'
#' @param Y             Observation matrix (n_obs x T); columns = time periods
#' @param dr_slack      Slack-regime DecisionRules (from solve_perturbation)
#' @param regime_cache  R environment of per-regime policies built by
#'                      obc_ensure_policy() / obc_guess_verify()
#' @param model         dynhr_mod
#' @param params        Named numeric parameter vector
#' @param obs_vars      Character vector of observed variable names
#' @param regime_path   Integer vector (length T): bitfield regime per period
#'                      (0 = all slack; bit j = 1 means spec j binds)
#' @param me_variance   Scalar measurement error variance added to F (default 1e-8)
#' @param return_filtered Logical: if TRUE, return the n_state x T filtered
#'                        state matrix (needed by obc_guess_verify)
#' @return List with:
#'   $loglik          -- total log-likelihood (scalar)
#'   $filtered_states -- n_state x T matrix if return_filtered, else NULL
#'   $n_obs, $n_T     -- dimensions
#' @noRd
kalman_filter_obc <- function(Y, dr_slack, regime_cache,
                               model, params, obs_vars, regime_path,
                               me_variance = 1e-8, return_filtered = FALSE) {

  endo    <- dr_slack$endo_names
  exo     <- dr_slack$exo_names
  n_state <- length(dr_slack$state_idx)
  n_exo   <- length(exo)
  n_obs   <- length(obs_vars)
  obs_idx <- match(obs_vars, endo)

  if (any(is.na(obs_idx)))
    stop("kalman_filter_obc: observed variables not in model: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  # ---- Slack policy (pre-sliced for fast access) ---------------------------
  pol_s <- get("0", envir = regime_cache, inherits = FALSE)
  TT_s  <- pol_s$TT
  RR_s  <- pol_s$RR
  ZZ_s  <- pol_s$ZZ
  DD_s  <- pol_s$DD

  # SS observable means (deviation-form: nominally zero, but carry for safety)
  d <- dr_slack$ys[obs_vars]

  Sigma_e <- .get_shock_cov(model, exo, params)

  # Noise covariances for slack regime
  QQ_s <- tcrossprod(RR_s %*% Sigma_e, RR_s)
  HH_s <- tcrossprod(DD_s %*% Sigma_e, DD_s)
  SS_s <- RR_s %*% Sigma_e %*% t(DD_s)   # state-obs cross-covariance

  # Local noise cache: regime_idx (character) → list(QQ, HH, SS)
  # Keyed separately from regime_cache because noise depends on Sigma_e (params)
  # but regime_cache (dr, matrices) does not.
  noise_cache <- new.env(parent = emptyenv(), hash = TRUE)
  assign("0", list(QQ = QQ_s, HH = HH_s, SS = SS_s), envir = noise_cache)

  # ---- Data setup ----------------------------------------------------------
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  if (length(regime_path) != n_T)
    stop(sprintf(
      "kalman_filter_obc: regime_path length (%d) != number of time periods (%d).",
      length(regime_path), n_T
    ))

  # ---- Initialisation ------------------------------------------------------
  # Covariance from the slack-regime Lyapunov equation — valid starting point
  # for any first-period regime; P converges within a few periods.
  s <- numeric(n_state)
  P <- solve_lyapunov(TT_s, QQ_s)
  if (any(!is.finite(P))) P <- diag(1e6, n_state)

  loglik   <- 0
  ll_const <- -0.5 * n_obs * log(2 * pi)

  filtered <- if (return_filtered) matrix(0, n_state, n_T) else NULL

  # ---- Main filter loop ----------------------------------------------------
  for (t in seq_len(n_T)) {
    regime_idx <- regime_path[t]
    key        <- as.character(regime_idx)

    if (regime_idx == 0L) {
      TT <- TT_s; RR <- RR_s; ZZ <- ZZ_s; DD <- DD_s
      QQ <- QQ_s; HH <- HH_s; SS <- SS_s
      c_state_t <- pol_s$c_state   # zero vector for slack
      d_eff     <- d
    } else {
      pol <- get(key, envir = regime_cache, inherits = FALSE)
      TT  <- pol$TT; RR <- pol$RR; ZZ <- pol$ZZ; DD <- pol$DD
      c_state_t <- pol$c_state

      if (!exists(key, envir = noise_cache, inherits = FALSE)) {
        QQ_r <- tcrossprod(RR %*% Sigma_e, RR)
        HH_r <- tcrossprod(DD %*% Sigma_e, DD)
        SS_r <- RR %*% Sigma_e %*% t(DD)
        assign(key, list(QQ = QQ_r, HH = HH_r, SS = SS_r), envir = noise_cache)
      }
      nc <- get(key, envir = noise_cache, inherits = FALSE)
      QQ <- nc$QQ; HH <- nc$HH; SS <- nc$SS
      d_eff <- d + pol$c_obs
    }

    # Innovation
    v <- Y[, t] - drop(ZZ %*% s) - d_eff

    # ---- Missing observations (partial obs at this period) -----------------
    if (any(is.na(v))) {
      obs_ok <- which(!is.na(v))

      # Advance state with no update if all obs missing
      if (length(obs_ok) == 0L) {
        s_pred <- drop(TT %*% s) + c_state_t
        P <- tcrossprod(TT %*% P, TT) + QQ
        P <- (P + t(P)) * 0.5
        s <- s_pred
        if (return_filtered) filtered[, t] <- s
        next
      }

      # Partial update: subset to observed rows
      ZZ_t    <- ZZ[obs_ok, , drop = FALSE]
      DD_t    <- DD[obs_ok, , drop = FALSE]
      HH_t    <- tcrossprod(DD_t %*% Sigma_e, DD_t)
      SS_t    <- RR %*% Sigma_e %*% t(DD_t)
      v_t     <- v[obs_ok]
      n_obs_t <- length(obs_ok)

      F_t    <- ZZ_t %*% P %*% t(ZZ_t) + HH_t + me_variance * diag(n_obs_t)
      F_chol <- chol(F_t)

      F_inv     <- chol2inv(F_chol)
      log_det_F <- 2 * sum(log(diag(F_chol)))
      loglik    <- loglik -
        0.5 * n_obs_t * log(2 * pi) -
        0.5 * (log_det_F + drop(crossprod(v_t, F_inv %*% v_t)))

      K      <- (TT %*% P %*% t(ZZ_t) + SS_t) %*% F_inv
      s_pred <- drop(TT %*% s) + drop(K %*% v_t) + c_state_t

      TmKZ <- TT - K %*% ZZ_t
      RmKD <- RR - K %*% DD_t
      P    <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD) +
              me_variance * tcrossprod(K)
      P    <- (P + t(P)) * 0.5

      s <- s_pred
      if (return_filtered) filtered[, t] <- s
      next
    }

    # ---- Full observation update ------------------------------------------
    F      <- ZZ %*% P %*% t(ZZ) + HH + me_variance * diag(n_obs)
    # Check positive definiteness via eigenvalues before Cholesky
    F_eig <- eigen(F, symmetric = TRUE, only.values = TRUE)$values
    if (min(F_eig) <= .Machine$double.eps * max(F_eig) * length(F_eig)) {
      s_pred <- drop(TT %*% s) + c_state_t
      P <- tcrossprod(TT %*% P, TT) + QQ
      P <- (P + t(P)) * 0.5
      s <- s_pred
      if (return_filtered) filtered[, t] <- s
      next
    }

    F_chol    <- chol(F)
    F_inv     <- chol2inv(F_chol)
    log_det_F <- 2 * sum(log(diag(F_chol)))
    ll_t      <- ll_const - 0.5 * (log_det_F + drop(crossprod(v, F_inv %*% v)))

    if (!is.finite(ll_t) || ll_t < -1e8) {
      return(list(loglik = -Inf, filtered_states = NULL, n_obs = n_obs, n_T = n_T))
    }
    loglik <- loglik + ll_t

    K      <- (TT %*% P %*% t(ZZ) + SS) %*% F_inv
    s_pred <- drop(TT %*% s) + drop(K %*% v) + c_state_t

    TmKZ <- TT - K %*% ZZ
    RmKD <- RR - K %*% DD
    P    <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD) +
            me_variance * tcrossprod(K)
    P    <- (P + t(P)) * 0.5

    s <- s_pred
    if (return_filtered) filtered[, t] <- s
  }

  if (return_filtered && !is.null(filtered))
    rownames(filtered) <- endo[dr_slack$state_idx]

  list(loglik = loglik, filtered_states = filtered, n_obs = n_obs, n_T = n_T)
}


# =============================================================================
# PKF helpers — Giovannini, Pfeiffer, Ratto (2021) §2.2
#
# Three helper functions implement eqs. 6-7 (shock extraction and one-step
# backward smoother).  kalman_filter_obc_pkf() wraps them into the full
# per-period inner convergence loop, replacing obc_guess_verify().
#
# The "inversion" label: eq. 6 inverts the innovation to recover the
# underlying shock eps_{t|t}, which is then used (together with the
# backward-smoothed state) to verify the OBC regime — rather than checking
# the slack prediction with zero shocks as obc_should_bind() does.
# =============================================================================

#' Extract the period-t shock estimate from Kalman filter quantities
#'
#' Implements GPR (2021) eq. 6:
#'   eps_{t|t} = Sigma_e * t(DD) * F_inv * v
#'
#' With dynhr's correlated-noise timing (y_t = ZZ*s_{t-1} + DD*eps_t),
#' DD encodes H*R in paper notation, so Cov(eps_t, v_t) = Sigma_e * t(DD).
#' Handles partial-observation subsetting: pass DD[obs_ok, , drop=FALSE].
#'
#' @param Sigma_e n_exo x n_exo shock covariance matrix
#' @param DD      n_obs_t x n_exo observation-shock loading (ghu[obs_ok, ])
#' @param F_inv   n_obs_t x n_obs_t inverse innovation covariance
#' @param v       length-n_obs_t innovation vector
#' @return length-n_exo vector eps_{t|t}
#' @noRd
pkf_extract_shock <- function(Sigma_e, DD, F_inv, v) {
  drop(Sigma_e %*% t(DD) %*% (F_inv %*% v))
}


#' One-step backward smoother: update s_{t-1|t-1} to s_{t-1|t}
#'
#' Implements the first step of the backward recursion in GPR (2021) eq. 7,
#' initialised with r_{t+1} = 0 so only period-t innovations contribute:
#'
#'   r_t       = t(ZZ) * F_inv * v          (r_{t+1} = 0 term vanishes)
#'   s_{t-1|t} = s_{t-1|t-1} + P_{t-1|t-1} * r_t
#'
#' Equivalently: s_{t-1|t} = s + Cov(s_{t-1}, v_t) * F^{-1} * v
#' where Cov(s_{t-1}, v_t) = P * t(ZZ) (the cross-covariance under the
#' lag-1 observation convention y_t = ZZ*s_{t-1} + DD*eps_t).
#'
#' @param s_prev  length-n_state filtered state s_{t-1|t-1}
#' @param P_prev  n_state x n_state filtered covariance P_{t-1|t-1}
#' @param ZZ      n_obs_t x n_state observation-state loading (ghx[obs_ok, ])
#' @param F_inv   n_obs_t x n_obs_t inverse innovation covariance
#' @param v       length-n_obs_t innovation vector
#' @return length-n_state backward-smoothed state s_{t-1|t}
#' @noRd
pkf_backward_one_step <- function(s_prev, P_prev, ZZ, F_inv, v) {
  drop(s_prev + P_prev %*% (t(ZZ) %*% (F_inv %*% v)))
}


#' Check which OBC constraints bind given backward-smoothed state and shock
#'
#' Predicts the constrained variable at t using the SLACK policy applied to
#' (s_backward, eps_hat):
#'
#'   var_t = ghx_slack[var_idx, ] * s_backward + ghu_slack[var_idx, ] * eps_hat
#'
#' This is the shock-aware replacement for obc_should_bind(), which uses
#' the lagged filtered state with zero shock instead.  Because eps_hat is
#' extracted from the actual observations (via pkf_extract_shock), the regime
#' check correctly reflects the driving force behind any constraint violation.
#'
#' @param s_backward length-n_state backward-smoothed state s_{t-1|t}
#' @param eps_hat    length-n_exo extracted shock eps_{t|t}
#' @param specs      OBC spec list from obc_parse_tags
#' @param dr_slack   Slack-regime DecisionRules
#' @return logical vector (length k): TRUE where constraint j binds
#' @noRd
pkf_check_binding <- function(s_backward, eps_hat, specs, dr_slack) {
  vapply(specs, function(s) {
    # var_pred is the predicted deviation of the constrained variable from SS.
    # The MCP bound is expressed in the SAME UNITS as the model variable (level
    # form for nonlinear models; deviation form = level form for linear models
    # where SS = 0).  We convert to level form before comparing:
    #   var_level = var_dev + SS_level
    # For linear models SS_level = 0, so this is a no-op.
    var_dev   <- sum(dr_slack$ghx[s$var_idx, ] * s_backward) +
                 sum(dr_slack$ghu[s$var_idx, ] * eps_hat)
    var_level <- var_dev + dr_slack$ys[s$var_idx]
    if (s$op == ">") var_level < s$bound else var_level > s$bound
  }, logical(1))
}


# =============================================================================
# PKF main filter
# =============================================================================

#' Piecewise Kalman Filter with per-period inner convergence (GPR 2021)
#'
#' Implements the full PKF algorithm from Giovannini, Pfeiffer, Ratto (2021)
#' S2.2.  At each period t the filter runs an inner guess-and-verify loop:
#'
#'   Step 1: Predict using current regime guess.
#'   Step 2: Update state and covariance (standard KF).
#'   Step 3: Extract shock \eqn{eps_{t|t}} (eq. 6) and backward-smooth the state
#'           to get \eqn{s_{t-1|t}} (eq. 7, one step, \eqn{r_{t+1}=0}).
#'   Step 4: Check OBC regime consistency via pkf_check_binding().
#'           If regime changed: cache new policy and restart inner loop.
#'           If regime matches: accept, record loglik, advance to t+1.
#'
#' Unlike kalman_filter_obc() + obc_guess_verify(), there is no separate
#' outer loop over the full sample; regime discovery is fully local to each t.
#'
#' @param data               n_obs x T observation matrix
#' @param dr_slack        Slack-regime DecisionRules
#' @param regime_cache    R environment of per-regime policies (modified
#'                        in-place; will be populated lazily as needed)
#' @param sys             System matrices (from extract_system_matrices_fast);
#'                        needed to lazily build new regime policies
#' @param model           dynhr_mod
#' @param params          Named numeric parameter vector
#' @param obs_vars        Character vector of observed variable names
#' @param specs           OBC spec list from obc_parse_tags
#' @param obs_idx         Integer vector: observable positions in endo vector
#' @param regime_path_init Integer vector (length T) for warm-starting the
#'                        inner loop (e.g. from a previous MCMC draw).
#'                        NULL or wrong length: all-slack initial guess.
#' @param me_variance     Measurement error variance (default 1e-8)
#' @param max_inner       Max inner iterations per period (default 10)
#' @param return_filtered Logical: return n_state x T filtered state matrix
#' @param return_shocks   Logical: return n_exo x T extracted shock matrix
#' @param return_store    Logical: return per-period KF quantities needed by
#'                        pkf_smoother_obc() for fixed-interval smoothing
#' @param return_P_last   Logical: return the final filtered covariance matrix
#'                        \eqn{P_{T|T}} in \code{$P_last} (used by sequential
#'                        estimation routines).
#' @return List with:
#'   $loglik          total log-likelihood
#'   $regime_path     integer vector (length T): accepted regime per period
#'   $filtered_states n_state x T if return_filtered, else NULL
#'   $filtered_shocks n_exo   x T if return_shocks,   else NULL
#'   $kf_store        smoother input list if return_store, else NULL
#'   $n_obs, $n_T
#'
#' @references
#'   Boehl, G. (2020). Efficient solution and computation of models with
#'     occasionally binding constraints. \emph{Deutsche Bundesbank Discussion
#'     Paper}, 38/2020.
#' @export
kalman_filter_obc_pkf <- function(data, dr_slack, regime_cache, sys,
                                   model, params, obs_vars, specs,
                                   obs_idx          = NULL,
                                   regime_path_init = NULL,
                                   me_variance      = 1e-8,
                                   max_inner        = 10L,
                                   return_filtered  = FALSE,
                                   return_shocks    = FALSE,
                                   return_store     = FALSE,
                                   return_P_last    = FALSE) {

  endo    <- dr_slack$endo_names
  exo     <- dr_slack$exo_names
  n_state <- length(dr_slack$state_idx)
  n_exo   <- length(exo)
  n_obs   <- length(obs_vars)

  if (is.null(obs_idx)) obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("kalman_filter_obc_pkf: observed variables not in model: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  # ---- Slack policy (pre-fetched) ------------------------------------------
  pol_s <- get("0", envir = regime_cache, inherits = FALSE)
  TT_s  <- pol_s$TT
  RR_s  <- pol_s$RR

  d       <- dr_slack$ys[obs_vars]
  Sigma_e <- .get_shock_cov(model, exo, params)

  QQ_s <- tcrossprod(RR_s %*% Sigma_e, RR_s)
  HH_s <- tcrossprod(pol_s$DD %*% Sigma_e, pol_s$DD)
  SS_s <- RR_s %*% Sigma_e %*% t(pol_s$DD)

  noise_cache <- new.env(parent = emptyenv(), hash = TRUE)
  assign("0", list(QQ = QQ_s, HH = HH_s, SS = SS_s), envir = noise_cache)

  # ---- Data setup ----------------------------------------------------------
  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  if (nrow(data) != n_obs) data <- t(data)
  n_T <- ncol(data)

  # Warm-start regime path (or all-slack)
  regime_path <- if (!is.null(regime_path_init) &&
                     length(regime_path_init) == n_T)
                   as.integer(regime_path_init)
                 else integer(n_T)

  # ---- Initialisation ------------------------------------------------------
  s <- numeric(n_state)
  P <- solve_lyapunov(TT_s, QQ_s)
  if (any(!is.finite(P))) P <- diag(1e6, n_state)

  loglik   <- 0
  filtered <- if (return_filtered) matrix(0, n_state, n_T) else NULL
  shocks   <- if (return_shocks)   matrix(0, n_exo,   n_T) else NULL

  # Smoother storage (per-period KF quantities for the backward pass).
  # INVARIANT: for every t, kf_store$L[[t]] and kf_store$RR[[t]] are matrices.
  # NULL assignments use single-bracket list(NULL) to preserve list length.
  if (return_store) {
    kf_store <- list(
      n_T   = n_T, n_state = n_state, n_exo = n_exo,
      Sigma_e = Sigma_e,
      endo_state_names = endo[dr_slack$state_idx],
      exo_names        = exo,
      v      = vector("list", n_T),
      F_inv  = vector("list", n_T),
      L      = vector("list", n_T),
      P_in   = vector("list", n_T),
      s_in   = matrix(0, n_state, n_T),
      TT     = vector("list", n_T),
      RR     = vector("list", n_T),
      ZZ     = vector("list", n_T),
      DD     = vector("list", n_T)
    )
  } else {
    kf_store <- NULL
  }

  # ---- Main filter loop ----------------------------------------------------
  for (t in seq_len(n_T)) {
    regime_idx <- regime_path[t]

    # Save incoming state/covariance before inner loop (for smoother storage)
    s_in_t <- s
    P_in_t <- P

    for (j in seq_len(max_inner)) {
      key <- as.character(regime_idx)

      # Lazily build policy for this regime if not yet cached
      if (!exists(key, envir = regime_cache, inherits = FALSE))
        obc_ensure_policy(regime_idx, regime_cache, sys, dr_slack, specs, obs_idx)

      pol       <- get(key, envir = regime_cache, inherits = FALSE)
      TT        <- pol$TT; RR <- pol$RR; ZZ <- pol$ZZ; DD <- pol$DD
      c_state_t <- pol$c_state
      d_eff     <- d + pol$c_obs

      # Noise covariances for this regime (cached)
      if (!exists(key, envir = noise_cache, inherits = FALSE)) {
        assign(key, list(
          QQ = tcrossprod(RR %*% Sigma_e, RR),
          HH = tcrossprod(DD %*% Sigma_e, DD),
          SS = RR %*% Sigma_e %*% t(DD)
        ), envir = noise_cache)
      }
      nc <- get(key, envir = noise_cache, inherits = FALSE)
      QQ <- nc$QQ; HH <- nc$HH; SS <- nc$SS

      # ---- Identify available observations ---------------------------------
      v_full  <- data[, t] - drop(ZZ %*% s) - d_eff
      obs_ok  <- which(!is.na(v_full))
      n_obs_t <- length(obs_ok)

      if (n_obs_t == 0L) {
        # No observations: advance state and covariance, accept current regime.
        # For the smoother: L = TT (pure transition, no observation update).
        # Use single-bracket list(NULL) to store NULL without removing the element.
        if (return_store) {
          kf_store$v[t]      <- list(NULL);  kf_store$F_inv[t] <- list(NULL)
          kf_store$L[[t]]    <- TT;    kf_store$P_in[[t]]  <- P_in_t
          kf_store$s_in[, t] <- s_in_t
          kf_store$TT[[t]]   <- TT;    kf_store$RR[[t]]    <- RR
          kf_store$ZZ[t]     <- list(NULL);  kf_store$DD[t]  <- list(NULL)
        }
        s <- drop(TT %*% s) + c_state_t
        P <- tcrossprod(TT %*% P, TT) + QQ
        P <- (P + t(P)) * 0.5
        regime_path[t] <- regime_idx
        if (return_filtered) filtered[, t] <- s
        break
      }

      # Subset to available observations
      v    <- v_full[obs_ok]
      ZZ_t <- ZZ[obs_ok, , drop = FALSE]
      DD_t <- DD[obs_ok, , drop = FALSE]
      HH_t <- tcrossprod(DD_t %*% Sigma_e, DD_t)
      SS_t <- RR %*% Sigma_e %*% t(DD_t)

      F_t    <- ZZ_t %*% P %*% t(ZZ_t) + HH_t + me_variance * diag(n_obs_t)
      # Check positive definiteness via eigenvalues before Cholesky
      F_eig <- eigen(F_t, symmetric = TRUE, only.values = TRUE)$values
      if (min(F_eig) <= .Machine$double.eps * max(F_eig) * length(F_eig)) {
        # Singular F: advance without update, accept current regime.
        if (return_store) {
          kf_store$v[t]      <- list(NULL);  kf_store$F_inv[t] <- list(NULL)
          kf_store$L[[t]]    <- TT;    kf_store$P_in[[t]]  <- P_in_t
          kf_store$s_in[, t] <- s_in_t
          kf_store$TT[[t]]   <- TT;    kf_store$RR[[t]]    <- RR
          kf_store$ZZ[t]     <- list(NULL);  kf_store$DD[t]  <- list(NULL)
        }
        s <- drop(TT %*% s) + c_state_t
        P <- tcrossprod(TT %*% P, TT) + QQ
        P <- (P + t(P)) * 0.5
        regime_path[t] <- regime_idx
        if (return_filtered) filtered[, t] <- s
        break
      }
      F_chol <- chol(F_t)
      F_inv <- chol2inv(F_chol)

      # ---- PKF inversion step (GPR 2021 eqs 6-7) ---------------------------
      eps_hat    <- pkf_extract_shock(Sigma_e, DD_t, F_inv, v)
      s_backward <- pkf_backward_one_step(s, P, ZZ_t, F_inv, v)

      # ---- Regime consistency check ----------------------------------------
      bind_flags <- pkf_check_binding(s_backward, eps_hat, specs, dr_slack)
      new_regime <- obc_regime_idx(bind_flags)

      if (new_regime != regime_idx) {
        # Regime changed: cache the new policy and re-run this period.
        obc_ensure_policy(new_regime, regime_cache, sys, dr_slack, specs, obs_idx)
        regime_idx <- new_regime

        if (j < max_inner) {
          next
        } else {
          # max_inner exhausted without convergence: fall back to a pure
          # state-prediction step (no observation update) using the last
          # accepted policy.  This avoids leaving kf_store[[t]] as NULL
          # which would crash the backward smoother.
          pol_fb     <- get(as.character(regime_idx), envir = regime_cache,
                            inherits = FALSE)
          TT_fb      <- pol_fb$TT
          QQ_fb      <- tcrossprod(pol_fb$RR %*% Sigma_e, pol_fb$RR)
          c_state_fb <- pol_fb$c_state
          if (return_store) {
            kf_store$v[t]      <- list(NULL);  kf_store$F_inv[t] <- list(NULL)
            kf_store$L[[t]]    <- TT_fb; kf_store$P_in[[t]]  <- P_in_t
            kf_store$s_in[, t] <- s_in_t
            kf_store$TT[[t]]   <- TT_fb; kf_store$RR[[t]]   <- pol_fb$RR
            kf_store$ZZ[t]     <- list(NULL);  kf_store$DD[t]  <- list(NULL)
          }
          s <- drop(TT_fb %*% s) + c_state_fb
          P <- tcrossprod(TT_fb %*% P, TT_fb) + QQ_fb
          P <- (P + t(P)) * 0.5
          regime_path[t] <- regime_idx
          if (return_filtered) filtered[, t] <- s
          break
        }
      }

      # ---- Regime verified: log-likelihood contribution --------------------
      log_det_F <- 2 * sum(log(diag(F_chol)))
      ll_t <- -0.5 * n_obs_t * log(2 * pi) -
              0.5 * (log_det_F + drop(crossprod(v, F_inv %*% v)))

      if (!is.finite(ll_t) || ll_t < -1e8) {
        return(list(
          loglik          = -Inf,
          regime_path     = regime_path,
          filtered_states = NULL,
          filtered_shocks = NULL,
          kf_store        = NULL,
          n_obs = n_obs, n_T = n_T
        ))
      }
      loglik <- loglik + ll_t

      # ---- Standard KF update (Joseph form) --------------------------------
      K      <- (TT %*% P %*% t(ZZ_t) + SS_t) %*% F_inv
      s_new  <- drop(TT %*% s) + drop(K %*% v) + c_state_t
      TmKZ   <- TT - K %*% ZZ_t
      RmKD   <- RR - K %*% DD_t
      P      <- tcrossprod(TmKZ %*% P, TmKZ) +
                tcrossprod(RmKD %*% Sigma_e, RmKD) +
                me_variance * tcrossprod(K)
      P      <- (P + t(P)) * 0.5

      # ---- Smoother storage (uses P_in_t / s_in_t saved before inner loop) --
      if (return_store) {
        kf_store$v[[t]]    <- v;      kf_store$F_inv[[t]] <- F_inv
        kf_store$L[[t]]    <- TmKZ;  kf_store$P_in[[t]]  <- P_in_t
        kf_store$s_in[, t] <- s_in_t
        kf_store$TT[[t]]   <- TT;    kf_store$RR[[t]]    <- RR
        kf_store$ZZ[[t]]   <- ZZ_t;  kf_store$DD[[t]]    <- DD_t
      }

      regime_path[t] <- regime_idx
      s <- s_new
      if (return_filtered) filtered[, t] <- s
      if (return_shocks)   shocks[, t]   <- eps_hat
      break
    }
  }

  if (return_filtered && !is.null(filtered))
    rownames(filtered) <- endo[dr_slack$state_idx]
  if (return_shocks && !is.null(shocks))
    rownames(shocks) <- exo

  ## After the main loop, P holds P_{T|T} (the post-update covariance at the
  ## last time period T).  Capture it when requested by the caller.
  P_last <- if (return_P_last) P else NULL

  list(
    loglik          = loglik,
    regime_path     = regime_path,
    filtered_states = filtered,
    filtered_shocks = shocks,
    kf_store        = kf_store,
    P_last          = P_last,
    n_obs = n_obs, n_T = n_T
  )
}
