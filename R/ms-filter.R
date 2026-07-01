## R/ms-filter.R
## --------------------------------------------------------------------------
## ms_kim_filter() -- Kim-Nelson (GPB(2)) filter for Markov-switching DSGE.
##
## Implements the Kim (1994) / Kim-Nelson (1999) filter for a state-space
## model where ONLY shock variances switch across regimes (structural
## parameters and TT, ZZ are common across regimes).
##
## Reference: Kim, C.-J. (1994), "Dynamic linear models with Markov-switching",
##   Journal of Econometrics 60(1-2), 1-22.
##   Kim, C.-J. & Nelson, C. R. (1999), "State-Space Models with Regime
##   Switching", MIT Press.
##
## dynhr state-space convention (lagged-state form):
##   s_t  = TT * s_{t-1} + RR * eps_t,   eps_t ~ N(0, I)   [ghu excl. Sigma_e]
##   y_t  = ZZ * s_{t-1} + DD * eps_t + d
##
## KEY INSIGHT: Because y_t depends on s_{t-1} (the PREVIOUS state), the
## innovation covariance uses the PRIOR covariance P_i, NOT the one-step-ahead
## predicted P_pred = TT*P_i*TT'+QQ. This is the fundamental difference from
## the standard (contemporaneous) state-space convention.
##
## Correct formulas for dynhr's lagged-state convention:
##   v_ij   = y_t - ZZ * b_i                      (innovation, uses prior b_i)
##   F_ij   = ZZ * P_i * ZZ' + HH_j               (uses PRIOR covariance P_i)
##   K_ij   = (TT * P_i * ZZ' + SS_j) * F_ij^{-1} (uses PRIOR P_i)
##   b_hat  = TT * b_i + K_ij * v_ij              (combined predict+update)
##   P_hat  = Joseph-form using PRIOR P_i
##
## Per-regime shock covariance:
##   Sigma_e^(j) = diag(scale_j) %*% Sigma_e %*% diag(scale_j)
##   QQ^(j)      = RR %*% Sigma_e^(j) %*% t(RR)
##   HH^(j)      = DD %*% Sigma_e^(j) %*% t(DD)
##   SS^(j)      = RR %*% Sigma_e^(j) %*% t(DD)
## --------------------------------------------------------------------------


#' Kim-Nelson filter for Markov-switching DSGE (shock-variance switching)
#'
#' Evaluates the log-likelihood of an MS-DSGE model where only shock variances
#' switch across regimes, using the Kim (1994) / GPB(2) filter.  The
#' structural parameters and decision rules (\code{ghx}, \code{ghu}) are
#' identical across regimes; only the per-regime shock covariance differs.
#'
#' @param Y  Observation matrix (\code{n_obs x T}).  May contain \code{NA}s;
#'   missing observations at period \code{t} are handled by skipping the
#'   measurement update (propagating the state through transition only).
#' @param dr  Decision rule (output of \code{\link{solve_perturbation}}).
#' @param model  Compiled model object (output of \code{\link{compile_model}}).
#' @param params  Named numeric vector of parameter values.
#' @param obs_vars  Character vector of observed variable names.
#' @param ms_spec  An \code{\link{ms_dsge_spec}} object.
#' @param me_variance  Scalar measurement-error variance added to the
#'   innovation covariance \code{F_{ij}} at every period (default \code{0}).
#' @param return_regime_probs  Logical; if \code{TRUE} return a
#'   \code{n_regimes x T} matrix of filtered regime probabilities
#'   \code{Pr[s_t = j | y_{1:t}]}.  Default \code{FALSE}.
#' @param lik_init  Character; state covariance initialisation: \code{"auto"},
#'   \code{"stationary"}, or \code{"kappa"}.  The exact-diffuse path is not
#'   yet supported for the MS filter; models with true unit roots should use
#'   \code{"kappa"}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{loglik}}{Total log-likelihood (scalar).}
#'     \item{\code{regime_probs}}{\code{n_regimes x T} matrix if
#'       \code{return_regime_probs = TRUE}, otherwise \code{NULL}.}
#'     \item{\code{n_obs}}{Number of observed variables.}
#'     \item{\code{n_T}}{Number of time periods.}
#'   }
#' @export
ms_kim_filter <- function(Y, dr, model, params, obs_vars, ms_spec,
                           me_variance = 0,
                           return_regime_probs = FALSE,
                           lik_init = c("auto", "stationary", "kappa")) {

  lik_init <- match.arg(lik_init)

  ## ---- input checks -------------------------------------------------------
  if (!inherits(ms_spec, "ms_dsge_spec"))
    stop("ms_kim_filter: ms_spec must be an ms_dsge_spec object.", call. = FALSE)

  h <- ms_spec$n_regimes
  P <- ms_spec$transition    # h x h, rows sum to 1

  ## ---- extract state-space matrices (common across regimes) ---------------
  state_idx <- dr$state_idx
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)
  n_obs     <- length(obs_vars)

  if (n_obs > n_exo)
    warning(sprintf("Stochastic singularity: %d obs but only %d shocks.", n_obs, n_exo))

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("ms_kim_filter: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  ## Validate that ms_spec shock-scale names match exo names (if names present)
  sc1 <- ms_spec$shock_scales[[1L]]
  if (!is.null(names(sc1))) {
    if (length(sc1) != n_exo ||
        !identical(sort(names(sc1)), sort(exo)))
      stop(sprintf(
        "ms_kim_filter: shock_scales names (%s) do not match model exo names (%s).",
        paste(names(sc1), collapse = ","), paste(exo, collapse = ",")),
        call. = FALSE)
    ## Reorder each regime's scale vector to match exo order
    ms_spec$shock_scales <- lapply(ms_spec$shock_scales, function(v) v[exo])
  } else if (length(sc1) != n_exo) {
    stop(sprintf(
      "ms_kim_filter: shock_scales length (%d) != n_exo (%d).",
      length(sc1), n_exo), call. = FALSE)
  }

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d   <- dr$ys[obs_vars]
  tZZ <- t(ZZ)

  Sigma_e <- .get_shock_cov(model, exo, params)

  ## ---- per-regime covariance matrices (QQ, HH, SS, Sigma_e per regime) ---
  regime_covs <- .ms_build_regime_covs(RR, DD, Sigma_e, ms_spec$shock_scales)

  ## ---- prepare observations -----------------------------------------------
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)
  Y_minus_d <- Y - d

  ## ---- initialise state distributions (one per regime) --------------------
  ## Resolve lik_init = "auto": try stationary P0 via Lyapunov, fall back kappa.
  if (lik_init == "auto") {
    ## Use the baseline (regime-1) QQ for the Lyapunov solve
    QQ0 <- regime_covs[[1L]]$QQ
    P0_try <- tryCatch(solve_lyapunov(TT, QQ0), error = function(e) NULL)
    ok_stat <- !is.null(P0_try) && all(is.finite(P0_try)) &&
      min(Re(eigen((P0_try + t(P0_try)) / 2, symmetric = TRUE,
                   only.values = TRUE)$values)) > -1e-8
    lik_init <- if (ok_stat) "stationary" else "kappa"
    P0 <- if (ok_stat) P0_try else .build_P0(TT, QQ0)
  } else if (lik_init == "stationary") {
    QQ0 <- regime_covs[[1L]]$QQ
    P0 <- tryCatch(solve_lyapunov(TT, QQ0),
                   error = function(e) .build_P0(TT, QQ0))
    if (anyNA(P0)) P0 <- .build_P0(TT, QQ0)
  } else {
    ## kappa: large-variance diffuse init
    QQ0 <- regime_covs[[1L]]$QQ
    P0  <- .build_P0(TT, QQ0)
  }

  ## State storage:
  ##   Beta[[j]]  = s_{t-1|t-1} for paths arriving in regime j
  ##               (the FILTERED state after measurement update at t-1)
  ##   Pvar[[j]]  = P_{t-1|t-1} for paths arriving in regime j
  ##               (the FILTERED covariance -- NOT the predicted P_{t-1|t})
  ## Initialised with s_0 = 0, P_0 = P0 for all regimes.
  s0   <- numeric(n_state)
  Beta <- replicate(h, s0, simplify = FALSE)
  Pvar <- replicate(h, P0, simplify = FALSE)

  ## Regime probability vector: Pr[s_{t-1} = j | y_{1:t-1}]
  regime_prob <- ms_spec$pi0

  ## ---- constant for Gaussian likelihood -----------------------------------
  ll_const <- -0.5 * n_obs * log(2 * pi)
  me_diag  <- me_variance * diag(n_obs)

  loglik  <- 0
  ll_floor <- -1e300

  ## Optional: store filtered regime probs (h x T)
  if (return_regime_probs)
    reg_prob_out <- matrix(0, h, n_T)

  ## ---- main Kim-Nelson loop -----------------------------------------------
  ##
  ## At each period t:
  ##   Beta[[i]], Pvar[[i]] = s_{t-1|t-1} and P_{t-1|t-1} conditional on
  ##                          being in regime i at t-1 (after collapsing).
  ##   regime_prob[i] = Pr[s_{t-1} = i | y_{1:t-1}]
  ##
  ## For each (i,j) path:
  ##   1. Predict: beta_pred = TT * b_i,  P_pred = TT * P_i * TT' + QQ_j
  ##      (one-step forecast of state, using regime-j shock variance)
  ##   2. Innovation (LAGGED-STATE): v_ij = y_t - ZZ * b_i
  ##      (ZZ acts on s_{t-1}, so b_i is the correct base, NOT beta_pred)
  ##   3. Innovation covariance (LAGGED-STATE): F_ij = ZZ * P_i * ZZ' + HH_j
  ##      (uses PRIOR P_i, not the predicted P_pred)
  ##   4. Kalman gain: K_ij = (TT * P_i * ZZ' + SS_j) * F_ij^{-1}
  ##   5. Updated state: b_hat_ij = TT * b_i + K_ij * v_ij
  ##   6. Updated covariance: P_hat_ij = Joseph-form
  ##   7. Likelihood weight: N(v_ij; 0, F_ij) * P[i,j] * regime_prob[i]
  ##
  ## Kim collapsing (h^2 -> h after Hamilton filter):
  ##   b_j = sum_i wts[i] * b_hat_ij
  ##   P_j = sum_i wts[i] * (P_hat_ij + (b_hat_ij - b_j)(b_hat_ij - b_j)')
  ##   The outer-product cross-term is CRITICAL; omitting it is the most common
  ##   implementation bug (causes underestimated variance and -Inf loglik at t>1).

  for (t in seq_len(n_T)) {
    y_t   <- Y_minus_d[, t]
    obs_ok <- is.finite(y_t)
    all_na <- !any(obs_ok)

    ## Storage for this period's h^2 paths
    ## LOG joint density of each (from_i -> to_j) path: lp_ij + log P[i,j] +
    ## log Pr[s_{t-1}=i]. Accumulated in log space (log-sum-exp below) so a
    ## tight regime whose Gaussian density underflows to 0 does NOT collapse
    ## f_y to 0 and spuriously return -Inf -- the failure mode when regimes
    ## have very different shock scales (e.g. a 0.5x volatility regime).
    log_lik_joint <- matrix(-Inf, h, h)   # [from_i, to_j]
    beta_hat  <- array(0, c(n_state, h, h))
    p_hat     <- array(0, c(n_state, n_state, h, h))

    for (i in seq_len(h)) {
      b_i <- Beta[[i]]
      P_i <- Pvar[[i]]   # P_{t-1|t-1} for regime i -- PRIOR covariance
      ## j-invariant blocks (depend only on the FROM-regime i) hoisted out of the
      ## j-loop: ZPZt_i is the state term of F_ij; TPZt_i is the gain numerator.
      ## Only HH_j / SS_j / QQ_j carry the TO-regime j dependence.
      ZPZt_i <- ZZ %*% P_i %*% tZZ
      TPZt_i <- TT %*% P_i %*% tZZ

      for (j in seq_len(h)) {
        cov_j <- regime_covs[[j]]
        QQ_j  <- cov_j$QQ
        HH_j  <- cov_j$HH
        SS_j  <- cov_j$SS

        ## -- Prediction step ------------------------------------------------
        ## beta_pred = TT * b_i (one-step-ahead state forecast)
        ## P_pred    = TT * P_i * TT' + QQ_j (one-step-ahead cov forecast)
        ## P_pred is used for the UPDATED covariance, not for F_ij.
        beta_pred <- drop(TT %*% b_i)

        if (all_na) {
          ## No observations: skip measurement update; propagate state only.
          ## P_pred is referenced only here and on the singular-F fallback, so it
          ## is formed lazily -- the common (observed, non-singular) path, which
          ## uses the Joseph form on P_i directly, never computes it.
          P_pred <- tcrossprod(TT %*% P_i, TT) + QQ_j
          P_pred <- (P_pred + t(P_pred)) * 0.5
          beta_hat[, i, j]  <- beta_pred
          p_hat[, , i, j]   <- P_pred
          log_lik_joint[i, j] <- log(P[i, j] * regime_prob[i])  # density 1
        } else {
          ## -- Innovation (LAGGED-STATE convention) --------------------------
          ## v_ij = y_t - ZZ * b_i   (b_i = s_{t-1|t-1}, not the predicted state)
          v_full <- y_t - as.numeric(ZZ %*% b_i)

          ## Handle partial NAs: zero out missing obs entries
          v_obs   <- v_full
          obs_idx_t <- obs_ok
          if (any(!obs_ok)) v_obs[!obs_ok] <- 0

          ## -- Innovation covariance F_ij (LAGGED-STATE) --------------------
          ## F_ij = ZZ * P_i * ZZ' + HH_j   (PRIOR P_i, not P_pred!)
          ## This is correct because y_t depends on s_{t-1}, so F is the
          ## variance of y_t - ZZ*s_{t-1|t-1} which uses P_{t-1|t-1} = P_i.
          F_ij <- ZPZt_i + HH_j + me_diag
          F_ij <- (F_ij + t(F_ij)) * 0.5

          ## Handle partial missing: project onto observed block
          if (any(!obs_ok)) {
            F_ij_obs <- F_ij[obs_ok, obs_ok, drop = FALSE]
            v_use    <- v_obs[obs_ok]
            n_obs_t  <- sum(obs_ok)
          } else {
            F_ij_obs <- F_ij
            v_use    <- v_obs
            n_obs_t  <- n_obs
          }

          Fc_ij <- tryCatch(chol(F_ij_obs), error = function(e) NULL)
          if (is.null(Fc_ij)) {
            ## Singular F for this path: zero likelihood contribution; keep prior
            ## (lazy P_pred -- see the all_na branch above).
            P_pred <- tcrossprod(TT %*% P_i, TT) + QQ_j
            P_pred <- (P_pred + t(P_pred)) * 0.5
            log_lik_joint[i, j] <- -Inf
            beta_hat[, i, j]  <- beta_pred
            p_hat[, , i, j]   <- P_pred
            next
          }

          Fi_ij     <- chol2inv(Fc_ij)
          log_det_F <- 2 * sum(log(diag(Fc_ij)))
          ll_const_t <- -0.5 * n_obs_t * log(2 * pi)
          quad       <- drop(crossprod(v_use, Fi_ij %*% v_use))
          lp_ij      <- ll_const_t - 0.5 * (log_det_F + quad)

          log_lik_joint[i, j] <- lp_ij + log(P[i, j] * regime_prob[i])

          ## -- Kalman gain (LAGGED-STATE) ------------------------------------
          ## K_ij = (TT * P_i * ZZ' + SS_j) * F_ij^{-1}   (uses PRIOR P_i)
          if (any(!obs_ok)) {
            ZZ_obs  <- ZZ[obs_ok, , drop = FALSE]
            DD_obs  <- DD[obs_ok, , drop = FALSE]
            SS_j_obs <- RR %*% cov_j$Sigma_e %*% t(DD_obs)
            ## TT %*% P_i %*% t(ZZ_obs) is the observed-column subset of TPZt_i.
            K_ij    <- (TPZt_i[, obs_ok, drop = FALSE] + SS_j_obs) %*% Fi_ij
            b_upd   <- beta_pred + drop(K_ij %*% v_use)
            IKZ     <- TT - K_ij %*% ZZ_obs
            RmKD    <- RR - K_ij %*% DD_obs
          } else {
            K_ij <- (TPZt_i + SS_j) %*% Fi_ij
            b_upd <- beta_pred + drop(K_ij %*% v_use)
            IKZ   <- TT - K_ij %*% ZZ
            RmKD  <- RR - K_ij %*% DD
          }

          ## -- Updated covariance (Joseph-form) ---------------------------
          ## P_hat = (TT - K*ZZ) * P_i * (TT - K*ZZ)' + (RR - K*DD) * Se * (RR - K*DD)'
          P_upd <- tcrossprod(IKZ %*% P_i, IKZ) +
                   tcrossprod(RmKD %*% cov_j$Sigma_e, RmKD)
          P_upd <- (P_upd + t(P_upd)) * 0.5

          beta_hat[, i, j]  <- b_upd
          p_hat[, , i, j]   <- P_upd
        }
      }
    }

    ## -- Loglik contribution for period t (log-sum-exp, underflow-safe) ------
    mx <- max(log_lik_joint)
    if (!is.finite(mx)) {            # every path has -Inf log-density
      loglik <- -Inf
      break
    }
    log_f_y <- mx + log(sum(exp(log_lik_joint - mx)))
    loglik  <- loglik + log_f_y
    if (loglik < ll_floor) { loglik <- -Inf; break }

    ## -- Hamilton filter: posterior joint regime probabilities ---------------
    prob_joint <- exp(log_lik_joint - log_f_y)   # h x h, sums to 1

    ## Marginal Pr[s_t = j | y_{1:t}]
    regime_prob_new <- colSums(prob_joint)    # length h

    ## Guard: floor at 1e-300 to avoid division by zero in collapsing
    regime_prob_new <- pmax(regime_prob_new, 1e-300)
    regime_prob_new <- regime_prob_new / sum(regime_prob_new)

    ## -- Kim collapsing: h^2 -> h -------------------------------------------
    ## Collapsed mean (Kim 1994, eq. 4):
    ##   b_j = sum_i Pr[s_{t-1}=i|s_t=j, y_{1:t}] * b_hat_ij
    ##        = sum_i (prob_joint[i,j] / regime_prob_new[j]) * b_hat_ij
    ##
    ## Collapsed covariance with CROSS-TERM (Kim 1994, eq. 5):
    ##   P_j = sum_i wts[i] * (P_hat_ij + (b_hat_ij - b_j)(b_hat_ij - b_j)')
    ##
    ## The outer-product term (b_hat_ij - b_j)(b_hat_ij - b_j)' accounts for
    ## the variance INCREASE from collapsing h separate paths into one weighted
    ## mean.  Omitting it systematically understates the covariance, causing
    ## F to shrink below the true innovation variance in subsequent periods,
    ## which makes the filter assign overly high likelihoods to later
    ## innovations and eventually diverge to -Inf.
    Beta_new <- vector("list", h)
    Pvar_new <- vector("list", h)

    for (j in seq_len(h)) {
      if (regime_prob_new[j] < 1e-300) {
        ## Near-zero probability: carry forward the prior state for this regime
        Beta_new[[j]] <- Beta[[j]]
        Pvar_new[[j]] <- Pvar[[j]]
        next
      }

      wts <- prob_joint[, j] / regime_prob_new[j]   # length h, sums to 1

      ## Collapsed mean
      b_j <- numeric(n_state)
      for (i in seq_len(h)) b_j <- b_j + wts[i] * beta_hat[, i, j]
      Beta_new[[j]] <- b_j

      ## Collapsed covariance WITH CROSS-TERM
      P_j <- matrix(0, n_state, n_state)
      for (i in seq_len(h)) {
        diff_ij <- beta_hat[, i, j] - b_j
        ## Critical: p_hat[,,i,j] + outer-product cross-term
        P_j <- P_j + wts[i] * (p_hat[, , i, j] + tcrossprod(diff_ij))
      }
      Pvar_new[[j]] <- (P_j + t(P_j)) * 0.5
    }

    Beta        <- Beta_new
    Pvar        <- Pvar_new
    regime_prob <- regime_prob_new

    if (return_regime_probs)
      reg_prob_out[, t] <- regime_prob
  }

  list(
    loglik       = loglik,
    regime_probs = if (return_regime_probs) reg_prob_out else NULL,
    n_obs        = n_obs,
    n_T          = n_T
  )
}


## ============================================================================
## Structural MS Kim-Nelson filter (regime-specific TT/ZZ/RR/DD)
## ============================================================================

#' Kim-Nelson filter for structural MS-DSGE (regime-specific decision rules)
#'
#' Evaluates the log-likelihood of a structural MS-DSGE model where BOTH
#' structural parameters AND shock variances can differ across regimes.
#' Each regime has its own state-space matrices \code{TT_s, RR_s, ZZ_s, DD_s}
#' extracted from an \code{MsDecisionRules} object (output of
#' \code{\link{solve_ms_perturbation}}).
#'
#' The shock covariance \eqn{\Sigma_e} is per-regime from the decision rules;
#' for shared-\eqn{\Sigma_e} with only structural params switching, supply the
#' same \code{Sigma_e} for all regimes (the common case).
#'
#' @param Y          Observation matrix (n_obs x T). May contain \code{NA}s.
#' @param ms_dr      An \code{MsDecisionRules} object from
#'                   \code{\link{solve_ms_perturbation}}.
#' @param model      Compiled model object.
#' @param params     Named numeric parameter vector (used only to compute
#'                   \eqn{\Sigma_e} if \code{Sigma_e_by_regime} is \code{NULL}).
#' @param obs_vars   Character vector of observed variable names.
#' @param Sigma_e_by_regime  Optional list of length h giving per-regime
#'                   shock covariance matrices. If \code{NULL}, a common
#'                   \eqn{\Sigma_e} from \code{params} is used for all regimes.
#' @param me_variance  Scalar measurement-error variance (default 0).
#' @param return_regime_probs  Logical; return \code{n_regimes x T} filtered
#'                   regime probability matrix (default FALSE).
#' @param lik_init   State covariance initialisation: \code{"auto"},
#'                   \code{"stationary"}, or \code{"kappa"}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{loglik}}{Total log-likelihood (scalar).}
#'     \item{\code{regime_probs}}{\code{n_regimes x T} matrix if
#'       \code{return_regime_probs = TRUE}, otherwise \code{NULL}.}
#'     \item{\code{n_obs}}{Number of observed variables.}
#'     \item{\code{n_T}}{Number of time periods.}
#'   }
#' @export
ms_kim_filter_struct <- function(Y, ms_dr, model, params, obs_vars,
                                  Sigma_e_by_regime = NULL,
                                  me_variance = 0,
                                  return_regime_probs = FALSE,
                                  lik_init = c("auto", "stationary", "kappa")) {

  lik_init <- match.arg(lik_init)

  if (!inherits(ms_dr, "MsDecisionRules"))
    stop("ms_kim_filter_struct: ms_dr must be an MsDecisionRules object.", call. = FALSE)

  h <- length(ms_dr$dr)
  P <- ms_dr$P

  ## ---- extract per-regime state-space matrices ----------------------------
  ## Use regime 1 to get shared structural info (state_idx, endo, exo).
  dr1       <- ms_dr$dr[[1L]]
  state_idx <- dr1$state_idx
  endo      <- dr1$endo_names
  exo       <- dr1$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)
  n_obs     <- length(obs_vars)

  if (n_obs > n_exo)
    warning(sprintf("Stochastic singularity: %d obs but only %d shocks.", n_obs, n_exo))

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("ms_kim_filter_struct: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  ## Per-regime TT_s, RR_s, ZZ_s, DD_s, d_s (observable SS mean)
  TT_list <- vector("list", h)
  RR_list <- vector("list", h)
  ZZ_list <- vector("list", h)
  DD_list <- vector("list", h)
  d_list  <- vector("list", h)

  for (s in seq_len(h)) {
    dr_s       <- ms_dr$dr[[s]]
    TT_list[[s]] <- dr_s$ghx[state_idx, , drop = FALSE]
    RR_list[[s]] <- dr_s$ghu[state_idx, , drop = FALSE]
    ZZ_list[[s]] <- dr_s$ghx[obs_idx,   , drop = FALSE]
    DD_list[[s]] <- dr_s$ghu[obs_idx,   , drop = FALSE]
    d_list[[s]]  <- dr_s$ys[obs_vars]
  }

  ## ---- per-regime shock covariance ----------------------------------------
  if (is.null(Sigma_e_by_regime)) {
    Sigma_e_common <- .get_shock_cov(model, exo, params)
    Sigma_e_by_regime <- replicate(h, Sigma_e_common, simplify = FALSE)
  }

  ## Build (QQ_s, HH_s, SS_s) per regime
  regime_covs <- vector("list", h)
  for (s in seq_len(h)) {
    RR_s  <- RR_list[[s]]
    DD_s  <- DD_list[[s]]
    Se_s  <- Sigma_e_by_regime[[s]]
    regime_covs[[s]] <- list(
      QQ      = tcrossprod(RR_s %*% Se_s, RR_s),
      HH      = tcrossprod(DD_s %*% Se_s, DD_s),
      SS      = RR_s %*% Se_s %*% t(DD_s),
      Sigma_e = Se_s
    )
  }

  ## ---- prepare observations -----------------------------------------------
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## Observable SS mean for regime 1 (used as reference; each regime has its own d_s)
  ## For the innovation, we use the FROM-regime's d (= dr_s$ys[obs_vars]).

  ## ---- initialise state distributions (one per regime) --------------------
  TT1 <- TT_list[[1L]]
  QQ1 <- regime_covs[[1L]]$QQ

  if (lik_init == "auto") {
    P0_try <- tryCatch(solve_lyapunov(TT1, QQ1), error = function(e) NULL)
    ok_stat <- !is.null(P0_try) && all(is.finite(P0_try)) &&
      min(Re(eigen((P0_try + t(P0_try)) / 2, symmetric = TRUE,
                   only.values = TRUE)$values)) > -1e-8
    lik_init <- if (ok_stat) "stationary" else "kappa"
    P0 <- if (ok_stat) P0_try else .build_P0(TT1, QQ1)
  } else if (lik_init == "stationary") {
    P0 <- tryCatch(solve_lyapunov(TT1, QQ1),
                   error = function(e) .build_P0(TT1, QQ1))
    if (anyNA(P0)) P0 <- .build_P0(TT1, QQ1)
  } else {
    P0 <- .build_P0(TT1, QQ1)
  }

  ## Regime probability: Pr[s_{t-1} = i | y_{1:t-1}]
  regime_prob <- ms_dr$pi0

  ## Initialise: s_0 = 0, P_0 for all regimes
  s0   <- numeric(n_state)
  Beta <- replicate(h, s0, simplify = FALSE)
  Pvar <- replicate(h, P0, simplify = FALSE)

  ## ---- Kim-Nelson loop (structural version) --------------------------------
  ##
  ## Key difference from shock-variance-only filter:
  ##   - Innovation uses FROM-regime i measurement matrices (ZZ_i, d_i)
  ##   - State prediction uses TO-regime j transition matrix (TT_j)
  ##   - Covariances use FROM-regime i P_i but TO-regime j QQ_j, HH_j, SS_j
  ##
  ## Following the comment block in ms-filter.R (dynhr lagged-state convention):
  ##   v_ij   = y_t - ZZ_i * b_i - d_i   (FROM-regime i measurement)
  ##   F_ij   = ZZ_i * P_i * ZZ_i' + HH_j (P_i from prior, HH_j to-regime)
  ##   K_ij   = (TT_j * P_i * ZZ_i' + SS_j) * F_ij^{-1}
  ##   b_hat  = TT_j * b_i + K_ij * v_ij
  ##   P_hat  = Joseph-form with TT_j - K_ij * ZZ_i and RR_j - K_ij * DD_i

  ll_const <- -0.5 * n_obs * log(2 * pi)
  me_diag  <- me_variance * diag(n_obs)
  loglik   <- 0
  ll_floor <- -1e300

  if (return_regime_probs)
    reg_prob_out <- matrix(0, h, n_T)

  for (t in seq_len(n_T)) {
    y_t    <- Y[, t]
    obs_ok <- is.finite(y_t)
    all_na <- !any(obs_ok)

    log_lik_joint <- matrix(-Inf, h, h)   # [from_i, to_j]
    beta_hat  <- array(0, c(n_state, h, h))
    p_hat     <- array(0, c(n_state, n_state, h, h))

    for (i in seq_len(h)) {
      b_i   <- Beta[[i]]
      P_i   <- Pvar[[i]]
      TT_j_list <- TT_list   # reference to avoid repeated indexing

      ## FROM-regime i measurement matrices
      ZZ_i  <- ZZ_list[[i]]
      DD_i  <- DD_list[[i]]
      d_i   <- d_list[[i]]
      tZZ_i <- t(ZZ_i)

      for (j in seq_len(h)) {
        TT_j  <- TT_j_list[[j]]
        RR_j  <- RR_list[[j]]
        DD_j  <- DD_list[[j]]
        cov_j <- regime_covs[[j]]
        QQ_j  <- cov_j$QQ
        HH_j  <- cov_j$HH
        SS_j  <- cov_j$SS

        ## Prediction: TO-regime j transition
        beta_pred <- drop(TT_j %*% b_i)
        P_pred    <- tcrossprod(TT_j %*% P_i, TT_j) + QQ_j
        P_pred    <- (P_pred + t(P_pred)) * 0.5

        if (all_na) {
          beta_hat[, i, j]    <- beta_pred
          p_hat[, , i, j]     <- P_pred
          log_lik_joint[i, j] <- log(P[i, j] * regime_prob[i])
        } else {
          ## Innovation (FROM-regime i): v = y_t - ZZ_i * b_i - d_i
          v_full  <- y_t - as.numeric(ZZ_i %*% b_i) - d_i
          v_obs   <- v_full
          if (any(!obs_ok)) v_obs[!obs_ok] <- 0

          ## Innovation covariance (FROM-regime i P_i, TO-regime j HH_j)
          F_ij  <- ZZ_i %*% P_i %*% tZZ_i + HH_j + me_diag
          F_ij  <- (F_ij + t(F_ij)) * 0.5

          if (any(!obs_ok)) {
            F_ij_obs <- F_ij[obs_ok, obs_ok, drop = FALSE]
            v_use    <- v_obs[obs_ok]
            n_obs_t  <- sum(obs_ok)
          } else {
            F_ij_obs <- F_ij
            v_use    <- v_obs
            n_obs_t  <- n_obs
          }

          Fc_ij <- tryCatch(chol(F_ij_obs), error = function(e) NULL)
          if (is.null(Fc_ij)) {
            log_lik_joint[i, j] <- -Inf
            beta_hat[, i, j]    <- beta_pred
            p_hat[, , i, j]     <- P_pred
            next
          }

          Fi_ij     <- chol2inv(Fc_ij)
          log_det_F <- 2 * sum(log(diag(Fc_ij)))
          ll_t      <- -0.5 * n_obs_t * log(2 * pi)
          quad      <- drop(crossprod(v_use, Fi_ij %*% v_use))
          lp_ij     <- ll_t - 0.5 * (log_det_F + quad)

          log_lik_joint[i, j] <- lp_ij + log(P[i, j] * regime_prob[i])

          ## Kalman gain: K_ij = (TT_j * P_i * ZZ_i' + SS_j) * F_ij^{-1}
          if (any(!obs_ok)) {
            ZZ_i_obs  <- ZZ_i[obs_ok, , drop = FALSE]
            DD_i_obs  <- DD_i[obs_ok, , drop = FALSE]
            SS_j_obs  <- RR_j %*% cov_j$Sigma_e %*% t(DD_i_obs)
            K_ij      <- (TT_j %*% P_i %*% t(ZZ_i_obs) + SS_j_obs) %*% Fi_ij
            b_upd     <- beta_pred + drop(K_ij %*% v_use)
            IKZ       <- TT_j - K_ij %*% ZZ_i_obs
            RmKD      <- RR_j - K_ij %*% DD_i_obs
          } else {
            K_ij  <- (TT_j %*% P_i %*% tZZ_i + SS_j) %*% Fi_ij
            b_upd <- beta_pred + drop(K_ij %*% v_use)
            IKZ   <- TT_j - K_ij %*% ZZ_i
            RmKD  <- RR_j - K_ij %*% DD_i
          }

          ## Joseph-form updated covariance
          P_upd <- tcrossprod(IKZ %*% P_i, IKZ) +
                   tcrossprod(RmKD %*% cov_j$Sigma_e, RmKD)
          P_upd <- (P_upd + t(P_upd)) * 0.5

          beta_hat[, i, j]  <- b_upd
          p_hat[, , i, j]   <- P_upd
        }
      }
    }

    ## Log-sum-exp for period t
    mx <- max(log_lik_joint)
    if (!is.finite(mx)) {
      loglik <- -Inf
      break
    }
    log_f_y <- mx + log(sum(exp(log_lik_joint - mx)))
    loglik  <- loglik + log_f_y
    if (loglik < ll_floor) { loglik <- -Inf; break }

    ## Hamilton filter: posterior joint regime probs
    prob_joint      <- exp(log_lik_joint - log_f_y)
    regime_prob_new <- colSums(prob_joint)
    regime_prob_new <- pmax(regime_prob_new, 1e-300)
    regime_prob_new <- regime_prob_new / sum(regime_prob_new)

    ## Kim collapsing: h^2 -> h
    Beta_new <- vector("list", h)
    Pvar_new <- vector("list", h)

    for (j in seq_len(h)) {
      if (regime_prob_new[j] < 1e-300) {
        Beta_new[[j]] <- Beta[[j]]
        Pvar_new[[j]] <- Pvar[[j]]
        next
      }
      wts <- prob_joint[, j] / regime_prob_new[j]
      b_j <- numeric(n_state)
      for (i in seq_len(h)) b_j <- b_j + wts[i] * beta_hat[, i, j]
      Beta_new[[j]] <- b_j

      P_j <- matrix(0, n_state, n_state)
      for (i in seq_len(h)) {
        diff_ij <- beta_hat[, i, j] - b_j
        P_j <- P_j + wts[i] * (p_hat[, , i, j] + tcrossprod(diff_ij))
      }
      Pvar_new[[j]] <- (P_j + t(P_j)) * 0.5
    }

    Beta        <- Beta_new
    Pvar        <- Pvar_new
    regime_prob <- regime_prob_new

    if (return_regime_probs)
      reg_prob_out[, t] <- regime_prob
  }

  list(
    loglik       = loglik,
    regime_probs = if (return_regime_probs) reg_prob_out else NULL,
    n_obs        = n_obs,
    n_T          = n_T
  )
}
