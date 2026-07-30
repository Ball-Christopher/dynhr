## R/hank-kalman.R
## --------------------------------------------------------------------------
## Kalman/posterior bridge for the general sequence-space HANK stack
## (R/hank-model.R), generalizing the scalar-shock MA companion form already
## built for the Krusell-Smith path in R/hank-estimation.R
## (hank_ma_state_space / hank_loglik_ss) to an arbitrary hank_model() with
## MULTIPLE named exogenous shocks, each following its own AR(1) process.
##
## CONSTRUCTION (truncated-MA companion form; matches hank_ma_state_space):
## The sequence-space solve gives impulse responses (MA coefficients), not a
## recursive law. For each exogenous shock z with AR(1) persistence rho_z, a
## unit innovation at t=0 produces the deterministic driving path
## dZ_t = rho_z^t, and hank_model_irf() gives the resulting deviation path of
## each observable -- this IS the shock's MA coefficient sequence Theta^z_s
## (s = 0, ..., T_h-1). Truncating at q <= T_h terms and packing a
## length-q shift register per shock (state s^z_t = (eps^z_t, ..., eps^z_{t-q+1}))
## gives an EXACT MA(q-1) state-space companion form; the only approximation
## is the truncation of the MA lag polynomial at q terms (same convention and
## same error source as hank_ma_state_space's `q`). q defaults to T_h (no
## truncation beyond the GE horizon itself). This construction is validated
## against the existing ABRS autocovariance likelihood in
## test-hank-kalman.R (oracle a): the two must agree to Kalman numerical
## precision when q = T_h (no truncation) and the gap must shrink as
## min(T_h, n_lags) grows when q is deliberately truncated below T_h.
##
## Shocks are packed as independent (Sigma_e diagonal in the AR-std
## parameterization used here); block-diagonal T/R and column-concatenated
## Z/D across shocks.
##
## NEAR-UNIT-ROOT CAVEAT: the truncation error of this companion form scales
## as rho^(2q)/(1 - rho^2) PER SHOCK, which explodes as any rho approaches 1
## (documented: -337 log-points at rho = 0.9992, q = 200).  For posteriors
## that visit high persistences use the exact-AR(1) stacked-covariance
## likelihood on the same object instead: hank_loglik_ar() in
## R/hank-kalman-ar.R (and likelihood = "exact_ar" in
## make_log_posterior_hank()), which keeps each shock's persistence tail in
## closed form.
## --------------------------------------------------------------------------


#' Truncated-MA state-space form of a general sequence-space HANK model
#'
#' Converts a linearised \code{\link{hank_model}} solution into an explicit
#' \code{\link{new_dsge_ss}} state space for a chosen set of observed
#' aggregates, given an AR(1) persistence/std spec for each exogenous shock.
#' Each shock's impulse response (\code{\link{hank_model_irf}} to a unit AR(1)
#' innovation) supplies its moving-average coefficients; the state stacks a
#' length-\code{q} shift register per shock (shocks enter independently, so
#' \code{T_mat}/\code{R_mat} are block-diagonal across shocks and
#' \code{Sigma_e} is diagonal).
#'
#' @param model A \code{\link{hank_model}} (e.g. from \code{\link{hank_ks_model}}
#'   or the \code{model} field of \code{\link{hank_nk_hank}}).
#' @param shock_specs Named list, one entry per \code{model$exogenous} shock,
#'   each \code{list(rho = <AR(1) persistence>, sigma = <innovation std>)}.
#'   Names must match \code{model$exogenous} exactly (order-independent).
#' @param observables Character vector of variable names (must be produced by
#'   the model, i.e. appear in \code{names(model$G)}) to treat as the
#'   observation vector, in order.
#' @param q Integer state length per shock (number of MA terms retained);
#'   default \code{model$T_h} (no truncation beyond the GE horizon).
#'
#' @return A lagged-timing \code{\link{new_dsge_ss}} object. Its
#'   \code{shock_names} field records \code{model$exogenous} order and
#'   \code{obs_names} records \code{observables}; \code{Theta_list},
#'   \code{rho_vec} and \code{q} record the per-shock MA coefficients,
#'   AR(1) persistences and truncation (consumed by
#'   \code{\link{hank_loglik_ar}}).
#' @seealso \code{\link{hank_kalman_loglik}}, \code{\link{hank_loglik_ar}},
#'   \code{\link{make_log_posterior_hank}}
#' @export
hank_state_space <- function(model, shock_specs, observables, q = NULL) {
  if (!inherits(model, "hank_model"))
    stop("hank_state_space: `model` must be a hank_model object.")
  exo <- model$exogenous
  if (!setequal(names(shock_specs), exo))
    stop("hank_state_space: `shock_specs` names must exactly match ",
         "model$exogenous = {", paste(exo, collapse = ", "), "}.")
  missing_obs <- setdiff(observables, names(model$G))
  if (length(missing_obs))
    stop("hank_state_space: observable(s) not produced by model: ",
         paste(missing_obs, collapse = ", "))

  T_h <- model$T_h
  if (is.null(q)) q <- T_h
  q <- min(q, T_h)
  n_obs   <- length(observables)
  n_shock <- length(exo)

  Theta_list <- .hank_theta_list(model, shock_specs, observables)

  ## Block-diagonal shift registers, one length-q block per shock.
  n_state <- q * n_shock
  TT <- matrix(0, n_state, n_state)
  RR <- matrix(0, n_state, n_shock)
  Z_lag <- matrix(0, n_obs, n_state)
  D_lag <- matrix(0, n_obs, n_shock)
  sigma_vec <- numeric(n_shock)
  rho_vec   <- numeric(n_shock)
  state_names <- character(n_state)

  for (k in seq_along(exo)) {
    z <- exo[k]
    idx <- ((k - 1L) * q + 1L):(k * q)
    if (q >= 2L) TT[idx[-1], idx[-q]][cbind(seq_len(q - 1L), seq_len(q - 1L))] <- 1
    RR[idx[1], k] <- 1
    Theta <- Theta_list[[z]]
    if (q >= 2L) for (s in 1:(q - 1L)) Z_lag[, idx[s]] <- Theta[s + 1L, ]
    D_lag[, k] <- Theta[1L, ]
    sigma_vec[k] <- shock_specs[[z]]$sigma
    rho_vec[k]   <- shock_specs[[z]]$rho
    state_names[idx] <- paste0(z, "_lag", seq_len(q) - 1L)
  }

  new_dsge_ss(T_mat = TT, R_mat = RR, Z_mat = Z_lag, D_mat = D_lag,
              Sigma_e = diag(sigma_vec^2, n_shock, n_shock),
              state_names = state_names, obs_names = observables,
              shock_names = exo, timing = "lagged",
              Theta_list = Theta_list, q = q,
              rho_vec = setNames(rho_vec, exo))
}


## Per-shock MA coefficients Theta^z (T_h x n_obs, one matrix per exogenous
## shock): a unit AR(1) innovation gives the deterministic driving path
## dZ_t = rho_z^t, fed through the full GE solve (hank_model_irf already
## accumulates through the whole block DAG to every variable in model$G).
## Shared by hank_state_space() and the exact_ar branch of
## make_log_posterior_hank() (which needs only Theta, not the companion-form
## matrices).  No validation; callers validate model/shock_specs/observables.
.hank_theta_list <- function(model, shock_specs, observables) {
  T_h <- model$T_h
  exo <- model$exogenous
  n_obs <- length(observables)
  Theta_list <- setNames(vector("list", length(exo)), exo)
  for (z in exo) {
    rho_z <- shock_specs[[z]]$rho
    dZ    <- setNames(list(rho_z^(seq_len(T_h) - 1L)), z)
    irf   <- hank_model_irf(model, dZ)
    Theta_list[[z]] <- matrix(sapply(observables, function(o) irf[[o]]),
                              T_h, n_obs,
                              dimnames = list(NULL, observables))
  }
  Theta_list
}


#' Kalman log-likelihood of aggregate HANK data on a general state space
#'
#' Runs dynhr's Gaussian Kalman filter core (\code{.kf_univariate_dispatch},
#' the same engine used by \code{kalman_filter}/\code{make_log_posterior})
#' directly on a \code{\link{hank_state_space}} object, without requiring a
#' \code{dr}-solution adapter. The initial state covariance is
#' \code{Sigma_e (x)-block per shock}, i.e. the stationary covariance of the
#' nilpotent shift-register state (exact, not approximate, since each block's
#' transition is strictly nilpotent).
#'
#' @param Y \code{T_data x n_obs} matrix (or data frame) of demeaned
#'   observations, columns in the same order as \code{ss$obs_names}.
#' @param ss A \code{\link{hank_state_space}} object.
#' @param me_var Measurement-error variance added to every observable's
#'   diagonal (default 0; needed when \code{n_obs > n_shock}, i.e. stochastic
#'   singularity).
#'
#' @return The scalar Gaussian log-likelihood.
#' @seealso \code{\link{hank_state_space}}, \code{\link{make_log_posterior_hank}}
#' @export
hank_kalman_loglik <- function(Y, ss, me_var = 0) {
  if (!inherits(ss, "dsge_ss"))
    stop("hank_kalman_loglik: `ss` must be a hank_state_space()/dsge_ss object.")
  Y <- as.matrix(Y)
  n_state <- ss$n_state
  q       <- ss$q
  n_shock <- ss$n_shock

  ## Stationary covariance of the block-nilpotent shift-register state: each
  ## shock's q-length register has stationary covariance sigma_z^2 * I_q
  ## (exact for a shift register fed by iid innovations), so the joint
  ## initial covariance is block-diagonal with sigma_z^2 * I_q blocks.
  P0 <- matrix(0, n_state, n_state)
  for (k in seq_len(n_shock)) {
    idx <- ((k - 1L) * q + 1L):(k * q)
    P0[idx, idx] <- ss$Sigma_e[k, k] * diag(q)
  }

  out <- .kf_univariate_dispatch(
    Y_minus_d = t(Y),                      # core expects n_obs x T
    ZZ = ss$Z_mat, TT = ss$T_mat, RR = ss$R_mat, DD = ss$D_mat,
    Sigma_e = ss$Sigma_e,
    s0 = rep(0, n_state),
    P_state = P0,
    P_inf_state = NULL,
    me_variance = me_var)
  out$loglik
}


#' Log-posterior factory for shock-process parameters of a HANK model
#'
#' Builds a \code{function(theta) -> list(logpost, loglik, logprior)} closure
#' (the interface expected by \code{rwmh}/\code{run_posterior_estimation}
#' style samplers) over a HANK shock-process parameter vector: one
#' persistence \code{rho_<shock>} and one std \code{sigma_<shock>} per
#' exogenous shock (household/steady-state parameters are held fixed --
#' re-solving the het-block steady state per draw is out of scope here).
#' Independent priors: Normal for each \code{rho} (soft-truncated to
#' \code{(-1, 1)} by returning \code{-Inf} outside the interval) and a
#' Half-Normal (implemented as \code{Normal(0, prior_sigma_sd)} folded to the
#' positive axis, i.e. \code{-Inf} for \code{sigma <= 0}) for each
#' \code{sigma}.
#'
#' @param model A \code{\link{hank_model}}.
#' @param Y \code{T_data x n_obs} matrix of demeaned observations.
#' @param observables Character vector of observable names (see
#'   \code{\link{hank_state_space}}).
#' @param q Optional truncation horizon passed to \code{\link{hank_state_space}}
#'   (\code{likelihood = "kalman"}) or \code{\link{hank_loglik_ar}}
#'   (\code{likelihood = "exact_ar"}).
#' @param me_var Measurement-error variance (see \code{\link{hank_kalman_loglik}}).
#' @param likelihood \code{"kalman"} (default) runs the Kalman filter on the
#'   truncated-MA companion form (\code{\link{hank_kalman_loglik}});
#'   \code{"exact_ar"} uses the exact-AR(1) stacked-covariance likelihood
#'   (\code{\link{hank_loglik_ar}}), which keeps each shock's AR(1) variance
#'   tail in closed form instead of truncating it with the MA lag polynomial
#'   -- prefer it whenever the \code{rho} posterior may visit near-unit-root
#'   values (see the representability bound in \code{\link{hank_loglik_ar}}:
#'   \eqn{|\rho|^{T_h}} should stay small, so bound/cap the persistences).
#'   A draw whose likelihood evaluation fails (e.g. a non-positive-definite
#'   stacked covariance at extreme \code{rho}) is rejected with
#'   \code{logpost = -Inf} rather than erroring the chain.
#' @param me_sd Per-observable measurement-error standard deviation for
#'   \code{likelihood = "exact_ar"} (scalar or length-\code{n_obs}); default
#'   \code{sqrt(me_var)}, so the two likelihood options describe the same
#'   measurement error unless overridden.
#' @param prior_rho_mean,prior_rho_sd Named numeric vectors (by shock name) or
#'   scalars (recycled) giving the Normal prior mean/sd for each shock's
#'   \code{rho}. Defaults \code{mean = 0.5, sd = 0.3}.
#' @param prior_sigma_sd Named numeric vector or scalar: half-Normal scale for
#'   each shock's \code{sigma}. Default \code{0.05}.
#' @param boundary Representability guard for \code{likelihood = "exact_ar"},
#'   applied per DRAW (the likelihoods' own \code{check_boundary} warns per
#'   call, which is unusable inside a sampler). \code{"warn"} (default) warns
#'   ONCE per closure, naming the offending shocks and the implied \code{rho}
#'   cap, and changes no posterior; \code{"reject"} makes
#'   \eqn{|\rho|^{T_h} \le} \code{boundary_tol} part of the PRIOR SUPPORT and
#'   returns \code{-Inf} outside it, exactly as the existing \eqn{|\rho| \ge 1}
#'   guard does; \code{"ignore"} restores the previous silent behaviour.
#'   Out here the sequence-space \code{Theta} is itself contaminated by the
#'   solve's terminal boundary, so the likelihood is WRONG rather than
#'   imprecise and no \code{q} repairs it -- see
#'   \code{\link{hank_theta_boundary_check}}.
#' @param boundary_tol Bound on \eqn{|\rho|^{T_h}}; default \code{1e-3}.
#'
#' @return A function \code{log_post_fn(theta)} where \code{theta} is a named
#'   numeric vector with entries \code{rho_<shock>} and \code{sigma_<shock>}
#'   for every \code{shock} in \code{model$exogenous}; returns
#'   \code{list(logpost, loglik, logprior)}.
#' @seealso \code{\link{hank_state_space}}, \code{\link{hank_kalman_loglik}},
#'   \code{\link{hank_loglik_ar}}, \code{rwmh}
#' @export
make_log_posterior_hank <- function(model, Y, observables, q = NULL,
                                    me_var = 0,
                                    likelihood = c("kalman", "exact_ar"),
                                    me_sd = NULL,
                                    prior_rho_mean = 0.5, prior_rho_sd = 0.3,
                                    prior_sigma_sd = 0.05,
                                    boundary = c("warn", "reject", "ignore"),
                                    boundary_tol = 1e-3) {
  likelihood <- match.arg(likelihood)
  boundary <- match.arg(boundary)
  if (!(is.numeric(boundary_tol) && length(boundary_tol) == 1L &&
        is.finite(boundary_tol) && boundary_tol > 0))
    stop("make_log_posterior_hank: `boundary_tol` must be a finite positive ",
         "scalar.")
  boundary_state <- new.env(parent = emptyenv())
  exo <- model$exogenous
  missing_obs <- setdiff(observables, names(model$G))
  if (length(missing_obs))
    stop("make_log_posterior_hank: observable(s) not produced by model: ",
         paste(missing_obs, collapse = ", "))
  if (is.null(me_sd)) me_sd <- sqrt(me_var)
  rep_named <- function(x, nm) {
    if (is.null(names(x))) setNames(rep(x, length.out = length(nm)), nm)
    else x[nm]
  }
  rho_mean <- rep_named(prior_rho_mean, exo)
  rho_sd   <- rep_named(prior_rho_sd, exo)
  sig_sd   <- rep_named(prior_sigma_sd, exo)
  ## Per-closure exact-AR cache: per-shock autocovariance slabs + the stacked
  ## gather index survive across draws (exactness-preserving; see the `cache`
  ## argument of hank_loglik_ar). A rho move re-derives that shock's Theta
  ## (different driving path) so its slab recomputes; a sigma-only move
  ## reuses every slab.
  ar_cache <- new.env(parent = emptyenv())

  function(theta) {
    rho_nm <- paste0("rho_", exo)
    sig_nm <- paste0("sigma_", exo)
    if (!all(c(rho_nm, sig_nm) %in% names(theta)))
      stop("make_log_posterior_hank: theta must have entries ",
           paste(c(rho_nm, sig_nm), collapse = ", "))
    rho   <- theta[rho_nm]; names(rho) <- exo
    sigma <- theta[sig_nm]; names(sigma) <- exo

    if (any(rho <= -1 | rho >= 1) || any(sigma <= 0))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = -Inf))

    logprior <- sum(stats::dnorm(rho, rho_mean, rho_sd, log = TRUE)) +
      sum(stats::dnorm(sigma, 0, sig_sd, log = TRUE) + log(2))  # half-normal

    shock_specs <- setNames(
      lapply(exo, function(z) list(rho = rho[[z]], sigma = sigma[[z]])), exo)
    if (likelihood == "kalman") {
      ss <- hank_state_space(model, shock_specs, observables, q = q)
      loglik <- hank_kalman_loglik(Y, ss, me_var = me_var)
    } else {
      ## Representability: reject/flag BEFORE spending a likelihood on a Theta
      ## the sequence-space solve cannot represent at this persistence.
      if (.hank_ar_boundary_gate(rho, model$T_h, boundary, boundary_tol,
                                 "make_log_posterior_hank", boundary_state))
        return(list(logpost = -Inf, loglik = NA_real_, logprior = logprior))
      ## exact_ar needs only the MA coefficients, not the (potentially large)
      ## q*n_shock-dimensional companion-form matrices.
      Theta_list <- .hank_theta_list(model, shock_specs, observables)
      loglik <- tryCatch(
        hank_loglik_ar(Y, Theta_list, rho = rho, sigma = sigma,
                       me_sd = me_sd, q = q, check_boundary = FALSE,
                       cache = ar_cache),
        error = function(e) -Inf)
      if (!is.finite(loglik))
        return(list(logpost = -Inf, loglik = loglik, logprior = logprior))
    }

    list(logpost = loglik + logprior, loglik = loglik, logprior = logprior)
  }
}
