## R/obc-inversion-filter.R
## --------------------------------------------------------------------------
## Cuba-Borda / Guerrieri / Iacoviello / Zhong (2019, JAE 34:1073-1085)
## deterministic inversion-filter likelihood for OBC models.
##
## Requires an EXACTLY-IDENTIFIED system: n_obs == n_exo, DD square and
## well-conditioned at every period.
##
## Provides:
##   kalman_filter_obc_inversion() -- inversion-filter log-likelihood
## --------------------------------------------------------------------------


#' Cuba-Borda et al. (2019) inversion-filter log-likelihood for OBC models
#'
#' Given a regime path (from the existing OccBin guess-and-verify machinery),
#' inverts the observation equation each period to recover the shock:
#'
#'   eps_t = DD_{r_t}^{-1} (y_t - ZZ_{r_t} s_{t-1} - d_{r_t})
#'
#' and advances the state deterministically:
#'
#'   s_t = TT_{r_t} s_{t-1} + RR_{r_t} eps_t + c_state_{r_t}
#'
#' The log-likelihood is the change-of-variables Jacobian result:
#'
#'   log L = -T/2 * n * log(2*pi)
#'           - sum_t log|det DD_{r_t}|
#'           - 0.5 * sum_t eps_t' Sigma_e^{-1} eps_t
#'
#' Hard requirements (enforced with stop()):
#'   - n_obs == n_exo  (exactly identified)
#'   - DD must be square and non-singular at every period
#'
#' @param Y             n_obs x T observation matrix (n_obs must equal n_exo)
#' @param dr_slack      Slack-regime DecisionRules (from solve_perturbation)
#' @param regime_cache  R environment of per-regime policies built by
#'                      obc_ensure_policy() / obc_guess_verify()
#' @param model         dynhr_mod
#' @param params        Named numeric parameter vector
#' @param obs_vars      Character vector of observed variable names
#'                      (length must equal n_exo)
#' @param regime_path   Integer vector (length T): bitfield regime per period
#' @return List with:
#'   $loglik      scalar log-likelihood
#'   $shocks      n_exo x T matrix of recovered shocks
#'   $n_obs       integer
#'   $n_T         integer
#' @noRd
kalman_filter_obc_inversion <- function(Y, dr_slack, regime_cache,
                                         model, params, obs_vars,
                                         regime_path) {

  endo    <- dr_slack$endo_names
  exo     <- dr_slack$exo_names
  n_state <- length(dr_slack$state_idx)
  n_exo   <- length(exo)
  n_obs   <- length(obs_vars)

  ## Hard requirement: exactly identified
  if (n_obs != n_exo)
    stop(sprintf(
      "kalman_filter_obc_inversion: requires n_obs == n_exo (exactly identified). Got n_obs=%d, n_exo=%d. Use a sub-set of observables equal to n_exo.",
      n_obs, n_exo
    ))

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("kalman_filter_obc_inversion: observed variables not in model: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  ## Data orientation
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  if (length(regime_path) != n_T)
    stop(sprintf(
      "kalman_filter_obc_inversion: regime_path length (%d) != n_T (%d).",
      length(regime_path), n_T
    ))

  ## Shock covariance and its inverse
  Sigma_e        <- .get_shock_cov(model, exo, params)
  Sigma_e_inv    <- solve(Sigma_e)
  log_det_Sigma  <- log(abs(det(Sigma_e)))

  ## Observable steady-state means
  d <- dr_slack$ys[obs_vars]

  ## Pre-fetch slack policy (regime 0)
  pol_s <- get("0", envir = regime_cache, inherits = FALSE)  # noqa: unused var kept for clarity

  ## ---- Main loop -----------------------------------------------------------
  s        <- numeric(n_state)   # initial state = zero (deviation form)
  loglik   <- 0
  ll_const <- -0.5 * n_obs * log(2 * pi)
  shocks   <- matrix(0, n_exo, n_T)

  for (t in seq_len(n_T)) {
    regime_idx <- regime_path[t]
    key        <- as.character(regime_idx)

    pol <- get(key, envir = regime_cache, inherits = FALSE)

    ZZ      <- pol$ZZ      # n_obs x n_state
    DD      <- pol$DD      # n_obs x n_exo  (must be square)
    TT      <- pol$TT      # n_state x n_state
    RR      <- pol$RR      # n_state x n_exo
    c_state <- pol$c_state # n_state constant
    d_eff   <- d + pol$c_obs

    ## Invert: eps_t = DD^{-1} (y_t - ZZ s_{t-1} - d_eff)
    rhs <- Y[, t] - drop(ZZ %*% s) - d_eff

    ## Solve for eps_t and compute log|det DD|
    DD_det <- tryCatch(det(DD), error = function(e) 0)
    if (!is.finite(DD_det) || DD_det == 0)
      return(list(loglik = -Inf, shocks = NULL, n_obs = n_obs, n_T = n_T))

    eps_t <- tryCatch(drop(solve(DD, rhs)), error = function(e) NULL)
    if (is.null(eps_t) || !all(is.finite(eps_t)))
      return(list(loglik = -Inf, shocks = NULL, n_obs = n_obs, n_T = n_T))

    log_det_DD <- log(abs(DD_det))

    ## Log-likelihood contribution (Jacobian change-of-variables from y to eps):
    ##   p(y_t) = p(eps_t) * |det(d eps_t / d y_t)|
    ##          = N(eps_t; 0, Sigma_e) * |det DD^{-1}|
    ## => log p(y_t) = -n/2*log(2pi) - 1/2*log|det Sigma_e|
    ##                 - log|det DD| - 1/2 eps_t' Sigma_e^{-1} eps_t
    ll_t <- ll_const - 0.5 * log_det_Sigma - log_det_DD -
            0.5 * drop(crossprod(eps_t, Sigma_e_inv %*% eps_t))

    if (!is.finite(ll_t))
      return(list(loglik = -Inf, shocks = NULL, n_obs = n_obs, n_T = n_T))

    loglik <- loglik + ll_t

    ## Advance state deterministically
    s <- drop(TT %*% s) + drop(RR %*% eps_t) + c_state
    shocks[, t] <- eps_t
  }

  rownames(shocks) <- exo
  list(loglik = loglik, shocks = shocks, n_obs = n_obs, n_T = n_T)
}
