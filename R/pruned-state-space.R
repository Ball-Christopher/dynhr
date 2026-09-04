## R/pruned-state-space.R
## --------------------------------------------------------------------------
## Order-2 pruned state-space object + Gaussian Kalman likelihood
##
## Implements the AFVRR (2018) pruned-SS estimation-grade SSM:
##   - pruned_state_space(dr, model, params) -- augmented-state constructor
##   - pruned_ss_moments(pss)                -- mean / Var / autocov
##   - pruned_ss_loglik(pss, Y, obs_vars)    -- Gaussian KF on augmented state
##
## References:
##   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F.
##     (2018). The pruned state-space system for non-linear DSGE models.
##     Review of Economic Studies, 85(1), 1-49. [AFVRR]
##
## The augmented state is:
##   xi_t = [ x1_t ; x2_t ; x1_t (x) x1_t ]   (d = 2*n_s + n_s^2)
##
## Transition:
##   xi_{t+1} = Tlin * xi_t + c + c_u + G * r_t
## Observation:
##   y_t = Dxi * xi_t + ys + ghss/2 + c_v + Gv * r_t
##
## where r_t = [eps_t ; eps_t(x)x1_t ; x1_t(x)eps_t ; eps_t(x)eps_t]
##
## Innovation covariance (at stationarity, a=0, P=Sigma_x):
##   Cr0  = Cov(r_t) = .order2_cov_r(0, Sigma_x, Sigma_e)
##   QQ   = G  * Cr0 * G'   (state noise, d x d)
##   HH   = Gv * Cr0 * Gv'  (obs noise,   n_obs x n_obs)
##   SS   = G  * Cr0 * Gv'  (cross,        d x n_obs)
##
## Cov(r_t) blocks (a = E[x1], M = P + aa', P = Cov(x1)):
##   eps                    : Sigma_e
##   eps(x)x1               : Sigma_e (x) a'   (cross with eps block: kron(Sigma_e, a'))
##   x1(x)eps               : a' (x) Sigma_e   (cross with eps block: kron(a', Sigma_e))
##   eps(x)eps (centered)   : Isserlis 4th-moment - vec(Sigma_e)vec(Sigma_e)'
##     = kron(Sigma_e, Sigma_e) + kron_comm(Sigma_e, Sigma_e)
##       (sum over all Wick contractions of the 4-index zero-mean Gaussian)
##   eps(x)x1, eps(x)x1     : kron(Sigma_e, M)
##   x1(x)eps, x1(x)eps     : kron(M, Sigma_e)
##   eps(x)x1, x1(x)eps     : C_ex1_x1e[i-s, k-u] = Sigma_e[i,u] * M[s,k]
##
## (At stationarity a=0 so the eps/eps(x)x1 and eps/x1(x)eps cross blocks
##  vanish; the kron(Sigma_e, M) and kron(M, Sigma_e) blocks remain because
##  M = P + 0 = Sigma_x != 0.)
##
## SUBTLE SPOT (Sigma_eta construction -- for orchestrator verification):
## See the detailed derivation in the comments at the end of this file.
## --------------------------------------------------------------------------


# ============================================================================
# Constructor
# ============================================================================

#' Build an order-2 pruned state-space object
#'
#' Extracts the AFVRR (2018) augmented linear state-space for an order-2
#' decision rule and returns a reusable \code{pruned_ss} object.  The
#' augmented state is \eqn{\xi_t = [x^{(1)}_t;\, x^{(2)}_t;\, x^{(1)}_t
#' \otimes x^{(1)}_t]} of dimension \eqn{d = 2 n_s + n_s^2}, where
#' \eqn{n_s} is the number of state variables.
#'
#' The machinery is a thin wrapper around the already-validated internal
#' \code{.order2_aug_system} helper in \code{stochsimul-monolith.R}.
#' \code{compute_moments_order2} calls the same internal, so
#' \code{pruned_ss_moments(pss)} is byte-identical to
#' \code{compute_moments_order2(dr, model, params)}.
#'
#' @param dr     A \code{DecisionRules2} object (output of
#'   \code{solve_perturbation(order = 2)}).
#' @param model  A parsed model object (output of \code{parse_mod}).
#' @param params Named numeric parameter vector.  Defaults to
#'   \code{model$param_values}.
#' @return An object of class \code{"pruned_ss"} with fields:
#'   \describe{
#'     \item{sys}{The raw augmented system list from \code{.order2_aug_system}.
#'       Carries \code{Tlin} (d x d transition), \code{G} (d x Dr selection),
#'       \code{Dxi} (n_endo x d observation map), \code{Gv} (n_endo x Dr
#'       direct noise), \code{cc}, \code{c_u} (constant drift in xi),
#'       \code{c_v} (constant drift in y), \code{Sigma_e}, and index vectors
#'       \code{ix1}, \code{ix2}, \code{ik}.}
#'     \item{Sigma_e}{n_exo x n_exo shock covariance.}
#'     \item{ys}{Named steady-state vector (n_endo).}
#'     \item{endo_names}{Character vector of endogenous variable names.}
#'     \item{exo_names}{Character vector of exogenous variable names.}
#'     \item{state_idx}{Integer index of state variables in the endo ordering.}
#'     \item{n_s}{Number of state variables.}
#'     \item{n_u}{Number of exogenous shocks.}
#'     \item{n_endo}{Number of endogenous variables.}
#'     \item{d}{Dimension of the augmented state (2*n_s + n_s^2).}
#'   }
#' @seealso \code{\link{pruned_ss_moments}}, \code{\link{pruned_ss_loglik}}
#' @export
pruned_state_space <- function(dr, model, params = NULL) {
  if (!inherits(dr, "DecisionRules2"))
    stop("pruned_state_space: dr must be a DecisionRules2 object.")

  if (is.null(params)) params <- model$param_values

  Sigma_e <- .get_shock_cov(model, dr$exo_names, params)
  sys     <- .order2_aug_system(dr, Sigma_e)

  structure(
    list(
      sys        = sys,
      Sigma_e    = Sigma_e,
      ys         = dr$ys,
      endo_names = dr$endo_names,
      exo_names  = dr$exo_names,
      state_idx  = dr$state_idx,
      n_s        = sys$n_s,
      n_u        = sys$n_u,
      n_endo     = sys$n_endo,
      d          = sys$d,
      dr         = dr,
      model      = model
    ),
    class = "pruned_ss"
  )
}


# ============================================================================
# Moments
# ============================================================================

#' Compute unconditional moments from a pruned state-space object
#'
#' Returns the same output as \code{\link{compute_moments_order2}} but accepts
#' a pre-built \code{\link{pruned_state_space}} object.  The computation
#' delegates entirely to \code{.order2_stationary_moments}, so results are
#' byte-identical to \code{compute_moments_order2}.
#'
#' @param pss  A \code{pruned_ss} object from \code{\link{pruned_state_space}}.
#' @param n_ar Number of autocorrelation lags (default 5).
#' @return A list with the same structure as \code{compute_moments_order2}:
#'   \code{mean}, \code{var_cov}, \code{std_dev}, \code{correlation},
#'   \code{autocorr}, \code{Sigma_e}, \code{Sigma_x}, \code{Var_x2},
#'   \code{mean_x2}.
#' @export
pruned_ss_moments <- function(pss, n_ar = 5L) {
  stopifnot(inherits(pss, "pruned_ss"))
  sys <- pss$sys

  st <- .order2_stationary_moments(sys)

  Sigma_x <- st$Sigma_x
  Var_x2  <- st$Var_x2
  mean_x2 <- st$mean_x2
  mn      <- st$mean
  names(mn) <- pss$endo_names
  Sigma_y <- st$var_cov
  rownames(Sigma_y) <- pss$endo_names
  colnames(Sigma_y) <- pss$endo_names

  ## Standard deviations and correlations
  variances <- pmax(diag(Sigma_y), 0)
  std_dev   <- sqrt(variances)
  names(std_dev) <- pss$endo_names

  sd_outer <- outer(std_dev, std_dev)
  sd_outer[sd_outer == 0] <- Inf
  corr_mat <- Sigma_y / sd_outer
  diag(corr_mat) <- 1
  rownames(corr_mat) <- pss$endo_names
  colnames(corr_mat) <- pss$endo_names

  ## Autocovariances (lag τ ≥ 1)
  n_endo <- pss$n_endo
  ghx    <- sys$hx  # n_s x n_s -- but we need the full ghx
  ## Reconstruct the full ghx from Dxi: Dxi = cbind(ghx_full, ghx_full, 0.5*ghxx_full)
  ## The first n_s columns of Dxi are ghx_full (n_endo x n_s).
  ## Actually Dxi is built from the full ghx, ghx, ghxx so extract correctly:
  ## Dxi[, ix1] = ghx[endo, ], Dxi[, ix2] = ghx[endo, ] (see .order2_aug_system)
  ## We only need ghx for the autocovariance propagation (same as compute_moments_order2).
  ## In compute_moments_order2 the full ghx is available directly from dr$ghx.
  ## Here we recover it from sys$Dxi: columns 1:n_s are ghx, n_s+1:2n_s are ghx again.
  n_s   <- sys$n_s
  ix1   <- sys$ix1
  ix2   <- sys$ix2
  ghx_full <- sys$Dxi[, ix1, drop = FALSE]   # n_endo x n_s

  ## S_sel: maps n_s x 1 state to n_s x n_endo (same as in compute_moments_order2)
  n_endo_full <- pss$n_endo
  sidx  <- pss$state_idx
  S_sel <- matrix(0, nrow = n_s, ncol = n_endo_full)
  for (i in seq_along(sidx)) S_sel[i, sidx[i]] <- 1

  autocorr <- array(0, dim = c(n_endo_full, n_endo_full, n_ar))
  endo <- pss$endo_names
  dimnames(autocorr) <- list(endo, endo, paste0("lag", seq_len(n_ar)))

  Gamma_prev_x1 <- ghx_full %*% Sigma_x %*% t(ghx_full)
  Gamma_prev_x2 <- ghx_full %*% Var_x2  %*% t(ghx_full)

  for (lag in seq_len(n_ar)) {
    Gamma_lag_x1 <- ghx_full %*% S_sel %*% Gamma_prev_x1
    Gamma_lag_x2 <- ghx_full %*% S_sel %*% Gamma_prev_x2
    Gamma_lag <- Gamma_lag_x1 + Gamma_lag_x2
    autocorr[, , lag] <- Gamma_lag / sd_outer
    Gamma_prev_x1 <- Gamma_lag_x1
    Gamma_prev_x2 <- Gamma_lag_x2
  }

  ## Third cumulant (skewness) -- delegate to compute_third_cumulant()
  skewness <- tryCatch(
    compute_third_cumulant(pss$dr, pss$model)$skewness,
    error = function(e) NULL
  )

  list(
    mean      = mn,
    var_cov   = Sigma_y,
    std_dev   = std_dev,
    correlation = corr_mat,
    autocorr  = autocorr,
    Sigma_e   = pss$Sigma_e,
    Sigma_x   = Sigma_x,
    Var_x2    = Var_x2,
    mean_x2   = mean_x2,
    skewness  = skewness
  )
}


# ============================================================================
# Gaussian Kalman likelihood on the augmented state
# ============================================================================

## ---------------------------------------------------------------------------
## Measurement-error floor guard (the near-degenerate-F hazard).
##
## When the model implies a NEAR-DEGENERATE one-step innovation covariance
## -- some linear combination of observables is almost perfectly predictable
## (smallest eigenvalue of the steady-state F orders of magnitude below its
## diagonal; e.g. a static rate that is nearly a combination of the other
## observables conditional on the past) -- even a "negligible" me_variance
## floor is LARGE relative to that eigenvalue. The floor then reshapes
## F^{-1} exactly in the direction where the likelihood's parameter
## discrimination is concentrated, producing a theta-dependent bias on data
## that do not carry that measurement error (2026-07-03: this artifact
## fully accounted for the P2d "misspecification bias" on rbc2shock, where
## min-eig(F) ~ 6e-9 vs the 1e-8 floor -- 158% -- while the floor was only
## 0.16% of the smallest F DIAGONAL, so no per-observable diagnostic can
## catch it; see ORDER3_PRUNED_SS_FOLLOWUP.md, MAJOR CORRECTION).
##
## Detector: iterate the correlated-noise Riccati to (near) steady state
## with HH_model (no floor) and compare me_variance against the smallest
## eigenvalue of the resulting F. Returns list(ratio = eigmin(F_0 + me I)/
## eigmin(F_0) = 1 + me/eigmin, loadings = |eigenvector|), or NULL when the
## me=0 recursion is not evaluable.
## ---------------------------------------------------------------------------
## Session memo for the detector: the result is a pure function of its
## inputs, and callers that leave me_floor_check = TRUE on a per-call path
## (e.g. kalman_filter called directly in a loop with me_variance > 0) would
## otherwise re-pay the full Riccati iteration on every call — measured at
## ~13x the cost of the small-model KF sweep itself (2026-08-05 perf-gate
## regression, introduced 7b2f217). Keyed on the full numeric content, so a
## hit is exact; bounded so a rogue per-draw caller cannot grow it without
## limit. Requires the optional 'digest' package; without it the detector
## simply recomputes (pre-memo behavior).
.me_floor_memo <- new.env(parent = emptyenv())

.pruned_me_floor_ratio <- function(Tlin, ZZ, QQ, HH_model, SS, Sxi0,
                                   me_variance, n_iter = 150L) {
  key <- NULL
  if (requireNamespace("digest", quietly = TRUE)) {
    key <- digest::digest(list(Tlin, ZZ, QQ, HH_model, SS, Sxi0,
                               me_variance, n_iter))
    hit <- get0(key, envir = .me_floor_memo, inherits = FALSE)
    if (!is.null(hit)) return(if (identical(hit, list())) NULL else hit)
  }
  P <- Sxi0
  res <- NULL
  for (k in seq_len(n_iter)) {
    Fm <- ZZ %*% P %*% t(ZZ) + HH_model
    Fi <- tryCatch(solve(Fm), error = function(e) NULL)
    if (is.null(Fi)) { P <- NULL; break }
    M <- Tlin %*% P %*% t(ZZ) + SS
    P_new <- Tlin %*% P %*% t(Tlin) + QQ - M %*% Fi %*% t(M)
    P_new <- (P_new + t(P_new)) * 0.5
    ## Converged-to-steady-state early exit: the detector only needs the
    ## fixed-point F, and stable systems typically converge in far fewer
    ## than n_iter steps.
    done <- max(abs(P_new - P)) <= 1e-12 * max(1, max(abs(P_new)))
    P <- P_new
    if (done) break
  }
  if (!is.null(P)) {
    F0 <- ZZ %*% P %*% t(ZZ) + HH_model
    eg <- eigen((F0 + t(F0)) * 0.5, symmetric = TRUE)
    emin <- min(eg$values)
    if (is.finite(emin) && emin > 0)
      res <- list(ratio    = 1 + me_variance / emin,
                  loadings = abs(eg$vectors[, which.min(eg$values)]))
  }
  if (!is.null(key)) {
    if (length(ls(.me_floor_memo, all.names = TRUE)) >= 64L)
      rm(list = ls(.me_floor_memo, all.names = TRUE), envir = .me_floor_memo)
    ## list() is the memo sentinel for a NULL (not-evaluable) result.
    assign(key, if (is.null(res)) list() else res, envir = .me_floor_memo)
  }
  res
}

.warn_me_floor_lock <- function(chk, obs_vars, me_variance, threshold = 1.5) {
  if (is.null(chk) || !is.finite(chk$ratio) || chk$ratio <= threshold)
    return(invisible(FALSE))
  main <- obs_vars[chk$loadings > 0.3]
  if (!length(main)) main <- obs_vars[which.max(chk$loadings)]
  warning(sprintf(paste0(
    "me_variance = %g is large relative to the smallest eigenvalue of the ",
    "model-implied (noise-free) one-step innovation covariance (inflates it ",
    "by factor %.2f): a linear combination of the observables (loading ",
    "mainly on %s) is nearly perfectly predictable by the model, so the ",
    "measurement noise you ASSUMED -- not the model -- dominates the ",
    "likelihood in that direction. me_variance is a genuine iid observation-",
    "noise variance (it enters both the innovation covariance and the state ",
    "update), so if your DATA really carry measurement error of this size, ",
    "ignore this warning. If the data were simulated withOUT measurement ",
    "error, this is a mis-specified noise model and it distorts the ",
    "likelihood (theta-dependently) -- set me_variance = 0. Disable this ",
    "check with options(dynhr.me_floor_check = FALSE)."),
    me_variance, chk$ratio, paste(main, collapse = ", ")),
    call. = FALSE)
  invisible(TRUE)
}

#' Gaussian Kalman filter log-likelihood on the pruned augmented state
#'
#' Runs a standard Kalman filter on the order-2 AFVRR augmented state
#' \eqn{\xi_t = [x^{(1)}_t;\, x^{(2)}_t;\, x^{(1)}_t \otimes x^{(1)}_t]}
#' using STATIONARY (time-invariant) noise covariances.  This is the
#' "cheap deterministic linear-KF on the AFVRR augmented pruned state"
#' described in the scope brief: cheaper than TPF, richer than the
#' unconditional-cumulant likelihood.
#'
#' The filter is initialized at the augmented stationary mean and covariance
#' (from \code{.order2_stationary_moments}).  Noise covariances are computed
#' at the stationary \eqn{\bar{a} = 0}, \eqn{\bar{P} = \Sigma_x} and held
#' constant throughout (the Gaussian approximation for the pruned system).
#'
#' The SSM used internally is the correlated-noise form:
#' \preformatted{
#'   xi_{t+1} = Tlin * xi_t + c_drift + G  * r_t
#'   y_t      = Dxi  * xi_t + d_y     + Gv * r_t
#' }
#' with \eqn{E[r_t r_t'] = Cr_0} (stationary), giving:
#' \preformatted{
#'   QQ = G  * Cr0 * G'   (d x d)
#'   HH = Gv * Cr0 * Gv'  (n_obs x n_obs)
#'   SS = G  * Cr0 * Gv'  (d x n_obs)
#' }
#' The correlated-noise KF converts to the uncorrelated form via the standard
#' substitution \eqn{\tilde{x}_{t+1} = \xi_{t+1} - S F^{-1} v_t} before each
#' step.
#'
#' @param pss      A \code{pruned_ss} object from
#'   \code{\link{pruned_state_space}}.
#' @param Y        Observation matrix (\code{n_obs x T} or \code{T x n_obs};
#'   transposed automatically when \code{ncol(Y) == n_obs}).  \code{NA}
#'   entries are handled with an exact PARTIAL measurement update: for each
#'   period only the non-\code{NA} observables enter the Kalman update (their
#'   rows/cols of \code{ZZ}/\code{HH}/\code{SS} and the corresponding
#'   log-density normalizing constant); a period with ALL observables missing
#'   falls back to a predict-only step (no update). Fully-observed and
#'   fully-missing periods are unaffected by this and match the previous
#'   behavior exactly.
#'   \code{Y} must be in LEVELS: the observation intercept is
#'   \code{ys + ghss/2 + c_v} (the order-2 mean), not zero.  Passing
#'   deviations from the steady state -- e.g. raw \code{simulate_model()} /
#'   \code{simulate_model_order2()} output, which are deviations -- makes
#'   every innovation carry the whole steady state and silently inflates
#'   \eqn{v_t^2/F_t} by orders of magnitude (F4-C: it turned an
#'   \eqn{O(\sigma^2)} model difference into a \eqn{\sigma}-invariant 2.8-nat
#'   offset against the linear Kalman filter).  Add \code{dr$ys[obs_vars]}
#'   to simulated deviations before filtering.
#' @param obs_vars Character vector of observed variable names (must be a
#'   subset of \code{pss$endo_names}).
#' @param me_variance  Scalar measurement-error jitter added to the diagonal
#'   of the innovation covariance (default 0). CAUTION: when the model
#'   implies a near-degenerate innovation covariance (some combination of
#'   observables is almost perfectly predictable one step ahead), even a
#'   tiny floor (1e-8) can dominate the smallest eigenvalue of F and bias
#'   likelihood comparisons on measurement-error-free data, theta-dependently.
#'   Match \code{me_variance} to what the data actually contain; the
#'   \code{me_floor_check} guard warns when the hazard is live.
#' @param me_floor_check Logical: when \code{me_variance > 0}, compare
#'   \code{me_variance} against the smallest eigenvalue of the filter's
#'   steady-state innovation covariance at \code{me_variance = 0} and warn
#'   if the floor inflates that eigenvalue by more than 50\%. Default
#'   \code{getOption("dynhr.me_floor_check", TRUE)}.
#' @return Scalar log-likelihood (numeric).  Returns \code{-Inf} on
#'   non-finite innovation covariance or other filter failure.
#'
#' @section Mixed-frequency data:
#' The exact partial-NA update makes \code{Y} a natural home for
#' mixed-frequency panels: give each observable its own row and set the
#' periods where it is unobserved to \code{NA}.  For example, a "monthly"
#' series observed every period alongside a "quarterly" series observed
#' only every third period:
#' \preformatted{
#'   T_n <- 120L
#'   Y <- rbind(monthly_series, quarterly_series)   # 2 x T_n, levels
#'   Y[2, -seq(3L, T_n, by = 3L)] <- NA              # quarterly: keep every 3rd
#'   ll <- pruned_ss_loglik(pss, Y, c("monthly_var", "quarterly_var"))
#' }
#' Each period's log-density uses only the observables that are non-\code{NA}
#' that period; a period with every observable missing contributes no update
#' (predict-only). See \code{Y} above for the exact partial-update semantics.
#'
#' @examples
#' \donttest{
#' mod_path <- system.file("extdata/models/rbc2shock.mod", package = "dynhr")
#' m   <- dynhr:::parse_mod(mod_path)
#' cm  <- dynhr:::compile_model(m, verbose = FALSE, max_order = 2L)
#' ss  <- dynhr:::solve_steady(cm, m$param_values, verbose = FALSE)
#' dr2 <- solve_perturbation(m, cm, ss$values, m$param_values,
#'                            order = 2L, verbose = FALSE)
#' pss <- pruned_state_space(dr2, m, m$param_values)
#'
#' ## Simulate data, then observe "c" every period ("monthly") and "y" only
#' ## every 3rd period ("quarterly").
#' T_n <- 120L
#' sim <- simulate_model_order2(dr2, n_periods = T_n, burn_in = 1000L,
#'                               model = m, pruning = TRUE)
#' Y <- t(sim[, c("c", "y")]) + dr2$ys[c("c", "y")]
#' Y[2, -seq(3L, T_n, by = 3L)] <- NA
#'
#' ll <- pruned_ss_loglik(pss, Y, c("c", "y"), me_variance = 1e-4)
#' ll
#' }
#' @seealso \code{\link{pruned_state_space}}, \code{\link{pruned_ss_moments}}
#' @export
pruned_ss_loglik <- function(pss, Y, obs_vars, me_variance = 0,
                             me_floor_check = getOption("dynhr.me_floor_check",
                                                        TRUE)) {
  stopifnot(inherits(pss, "pruned_ss"))
  sys <- pss$sys

  ## -- Observation index ------------------------------------------------------
  obs_idx <- match(obs_vars, pss$endo_names)
  if (any(is.na(obs_idx)))
    stop("pruned_ss_loglik: obs_vars not found in pss: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  n_obs <- length(obs_vars)

  ## -- Reshape Y to n_obs x T -------------------------------------------------
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## -- Augmented-state system matrices ----------------------------------------
  Tlin  <- sys$Tlin    # d x d
  G     <- sys$G       # d x Dr
  Dxi_f <- sys$Dxi     # n_endo x d  (full obs map)
  Gv_f  <- sys$Gv      # n_endo x Dr (full direct noise)
  d_dim <- sys$d

  ## Subset to observed variables
  ZZ  <- Dxi_f[obs_idx, , drop = FALSE]    # n_obs x d
  Gv  <- Gv_f[obs_idx,  , drop = FALSE]    # n_obs x Dr

  ## Observation intercept: ys[obs_vars] + ghss[obs_vars]/2 + c_v[obs_vars]
  ## In compute_moments_order2 the obs mean is: ys + Dxi*mu_xi + 0.5*ghss + c_v
  ## For the KF the intercept is ys + 0.5*ghss + c_v (the Dxi*xi part is state)
  d_y <- pss$ys[obs_vars] +
         0.5 * sys$ghss[obs_idx] +
         sys$c_v[obs_idx]

  ## State drift (absorbed into initial condition at stationarity)
  c_drift <- sys$cc + sys$c_u   # d-vector

  ## -- Stationary noise covariances -------------------------------------------
  ## Compute Sigma_x (first-order layer covariance)
  Sigma_x <- solve_lyapunov(sys$hx, sys$hu %*% sys$Sigma_e %*% t(sys$hu))

  ## Cr0 = Cov(r_t) at stationarity (a=0, P=Sigma_x)
  Cr0 <- .order2_cov_r(numeric(sys$n_s), Sigma_x, sys$Sigma_e)

  ## Augmented state noise covariance
  QQ  <- G  %*% Cr0 %*% t(G)            # d x d
  QQ  <- (QQ + t(QQ)) * 0.5             # symmetrize
  ## Observation noise covariance
  HH  <- Gv %*% Cr0 %*% t(Gv)          # n_obs x n_obs
  HH  <- (HH + t(HH)) * 0.5
  HH_model <- HH                        # model-implied part, pre-floor
  if (me_variance > 0) HH <- HH + me_variance * diag(n_obs)
  ## Cross-covariance state/observation noise
  SS  <- G  %*% Cr0 %*% t(Gv)          # d x n_obs

  ## -- Stationary initial conditions ------------------------------------------
  st  <- .order2_stationary_moments(sys)
  ## Stationary augmented-state mean
  mu0 <- as.numeric(solve(diag(d_dim) - Tlin, c_drift))
  ## Stationary augmented-state covariance (Lyapunov fixed point)
  Sxi0 <- solve_lyapunov(Tlin, QQ)
  Sxi0 <- (Sxi0 + t(Sxi0)) * 0.5

  ## -- Measurement-error floor guard ------------------------------------------
  if (me_variance > 0 && isTRUE(me_floor_check))
    .warn_me_floor_lock(
      .pruned_me_floor_ratio(Tlin, ZZ, QQ, HH_model, SS, Sxi0, me_variance),
      obs_vars, me_variance)

  ## -- Kalman filter (correlated noise, constant-gain) ------------------------
  ##
  ## Correlated-noise form: Cov(transition noise, obs noise) = SS.
  ## Standard KF requires uncorrelated noise.  Convert via:
  ##   xi_{t+1}' = xi_{t+1} - SS * F_t^{-1} * v_t
  ##   <=> predicting xi_{t+1} after absorbing the observation correlation.
  ##
  ## The textbook correlated-noise KF (Anderson & Moore 1979, Ch. 8):
  ##   Predict:   xi_{t|t-1} = Tlin * xi_{t-1|t-1} + c_drift
  ##              P_{t|t-1}  = Tlin * P_{t-1|t-1} * Tlin' + QQ
  ##   Innovation:v_t = y_t - d_y - ZZ * xi_{t|t-1}
  ##              F_t = ZZ * P_{t|t-1} * ZZ' + HH
  ##   Gain:      K_t = (Tlin * P_{t|t-1} * ZZ' + SS) * F_t^{-1}
  ##   Update:    xi_{t|t} = xi_{t|t-1} + (K_t - SS*F_t^{-1}*ZZ)*P_{t|t-1}^{-1}*...
  ##                  (equivalently via the Joseph form)
  ##   Actually the standard correlated form:
  ##     K_t = (P_{t|t-1} * ZZ' * F_t^{-1}) from state side; but here the
  ##     cross-covariance comes from the TRANSITION side (G*Cr0*Gv'), so:
  ##     M_t = Tlin * P_{t|t-1} * ZZ' + SS  (d x n_obs)
  ##     K_t = M_t * F_t^{-1}
  ##     xi_{t|t} = xi_{t|t-1} + P_{t|t-1} * ZZ' * F_t^{-1} * v_t
  ##     P_{t|t}  = P_{t|t-1} - P_{t|t-1} * ZZ' * F_t^{-1} * ZZ * P_{t|t-1}
  ##     P_{t+1|t} = Tlin*P_{t|t}*Tlin' + QQ - K_t*F_t*K_t' + K_t*M_t'+ M_t*K_t'
  ##               = ... (reduces to standard when SS=0)
  ##
  ## Use the equivalent "innovation form" (Shumway & Stoffer; DK ch. 4):
  ##   Define  A_t = Tlin - K_t * ZZ  (d x d)
  ##   Then    P_{t+1|t} = A_t * P_{t|t-1} * A_t' + QQ - SS*F_t^{-1}*SS'  ...
  ##   This is numerically standard but complex.
  ##
  ## Simplest correct implementation for the constant-coefficient case:
  ## Use the standard KF on the AUGMENTED-AUGMENTED model that absorbs SS
  ## via the Kalman recursion for the modified transition.  Or directly
  ## implement the correlated-noise KF update as above.
  ##
  ## We implement it straightforwardly following Anderson & Moore (1979):
  ##   F_t   = ZZ * P_t * ZZ' + HH
  ##   M_t   = Tlin * P_t * ZZ' + SS    (d x n_obs)
  ##   K_t   = M_t * F_t^{-1}           (d x n_obs  Kalman gain)
  ##   xi_{t+1|t} = Tlin * xi_t + c_drift + K_t * v_t  - SS*F_t^{-1}*v_t
  ##              = (Tlin - SS*F_t^{-1}*ZZ)*xi_t + c_drift + K_t*v_t
  ## Wait -- let's be precise.  The standard correlated-noise KF equations are:
  ##   Predict: xi_{t|t-1} = Tlin * xi_{t-1|t-1} + c_drift
  ##            P_{t|t-1}  = Tlin * P_{t-1|t-1} * Tlin' + QQ
  ##   Innovation: v_t = y_t - d_y - ZZ * xi_{t|t-1}
  ##               F_t = ZZ * P_{t|t-1} * ZZ' + HH
  ##   Filter update:
  ##     xi_{t|t} = xi_{t|t-1} + L_t * v_t
  ##       where L_t = P_{t|t-1}*ZZ'*F_t^{-1}  (filter gain, state update side)
  ##     P_{t|t}  = (I - L_t*ZZ) * P_{t|t-1}
  ##   One-step-ahead prediction for t+1:
  ##     xi_{t+1|t} = Tlin * xi_{t|t} + c_drift
  ##       = Tlin*(xi_{t|t-1} + L_t*v_t) + c_drift
  ##     P_{t+1|t}  = Tlin * P_{t|t} * Tlin' + QQ
  ##                   + Tlin*L_t*SS' + SS*L_t'*Tlin'  <-- cross terms from correlation
  ## This is getting tangled.  The cleanest derivation uses the Woodbury identity.
  ##
  ## CLEAN APPROACH (following AFVRR appendix directly):
  ## The correlated-noise model can be converted to uncorrelated by:
  ##   Define xi_t* = xi_t - Tlin^{-1} * SS * F_t^{-1} * v_t  (not practical here)
  ##
  ## Instead, use the EQUIVALENT UNCORRELATED REPRESENTATION:
  ## Define z_t = xi_t - SS * HH^{-1} * (y_t - d_y - ZZ * xi_{t|t-1})  ... also messy.
  ##
  ## The most robust approach for a one-shot implementation:
  ## Apply the result from Chib, Nardari & Shephard (2002) / DK (2012) sec. 4.3:
  ##   Modified transition: Tlin* = Tlin - SS * HH^{-1} * ZZ  (when HH invertible)
  ##   Modified noise cov:  QQ*   = QQ  - SS * HH^{-1} * SS'
  ##   Then run standard KF on the transformed model with:
  ##     QQ_tilde = QQ*
  ##     observation noise = HH (unchanged)
  ##     no cross-covariance
  ##
  ## This is exact when HH is invertible (which it is for the pruned SS with
  ## me_variance >= 0, since HH = Gv*Cr0*Gv' is PSD and the obs eq has
  ## additive Gaussian noise from the quadratic innovations).
  ##
  ## BUT: HH = Gv * Cr0 * Gv' may be singular (stochastic singularity of obs).
  ## So we use the general correlated-noise KF directly:
  .pruned_kf_correlated(Y, Tlin, ZZ, d_y, c_drift, QQ, HH, SS, mu0, Sxi0)
}


## Internal: Kalman filter for correlated state/observation noise.
##
## Model:
##   xi_{t+1} = Tlin * xi_t + c_drift + w_t   Cov(w_t) = QQ
##   y_t      = ZZ   * xi_t + d_y     + v_t   Cov(v_t) = HH
##   Cov(w_t, v_{t+1}) = 0  (future obs uncorrelated with current state noise)
##   Cov(w_t, v_t)      = SS  (current state noise correlated with current obs noise)
##
## Note: in dynhr's augmented system, the innovation r_t drives BOTH xi_{t+1}
## (via G) and y_t (via Gv).  However the observation y_t is measured AFTER
## xi_t is observed (timing: y_t depends on xi_t, and the noise entering xi_{t+1}
## is the SAME r_t that also enters y_t).  In the pruned-SS timing convention
## (AFVRR eq 3-4): y_t = Dxi*xi_t + Gv*r_t, xi_{t+1} = Tlin*xi_t + G*r_t + c.
## So w_t and v_t for the SAME t are correlated: Cov(G*r_t, Gv*r_t) = G*Cr0*Gv' = SS.
##
## The correlated KF equations (DK 2012 sec 4.3 / AM 1979 ch 7):
##   Predict:   xi_pred = Tlin * xi_upd + c_drift
##              P_pred  = Tlin * P_upd  * Tlin' + QQ
##   Innovation: v = y - d_y - ZZ * xi_pred
##               F = ZZ * P_pred * ZZ' + HH
##   Modified gain: K = (Tlin * P_pred * ZZ' + SS) * inv(F)
##     = M * inv(F)  where M = Tlin*P_pred*ZZ' + SS
##   Filter update:
##     xi_upd_new = xi_pred + (P_pred * ZZ' * inv(F)) * v
##                = xi_pred + L * v   where L = P_pred*ZZ'*inv(F)
##     P_upd_new  = P_pred - L * ZZ * P_pred
##   Next-step predict (rolled into next iteration):
##     P_{t+2|t+1} = Tlin * P_upd_new * Tlin' + QQ
##   CORRECTION: in the correlated model the one-step-ahead covariance update is:
##     P_pred_new = Tlin*P_upd_new*Tlin' + QQ - cross terms
##   The standard result (e.g. Harvey 1990 eq 3.2.3a for the correlated case):
##     P_{t+1|t} = Tlin*P_{t|t-1}*Tlin' + QQ - M*F^{-1}*M'
##   (This absorbs the cross-terms; derivation: P_{t+1|t} is the prediction
##    error covariance for xi_{t+1} given Y_{1:t}.  The optimal predictor is
##    Tlin*xi_{t|t} + c_drift + SS*F_t^{-1}*v_t, and P_{t+1|t} follows.)
##
## Actually the recursion that avoids two-step propagation is:
##   Compute P_pred as usual; then:
##   F = ZZ*P_pred*ZZ' + HH
##   M = Tlin*P_pred*ZZ' + SS   (d x n_obs)
##   K = M * inv(F)              (d x n_obs, the modified gain)
##   xi_upd = xi_pred + L*v  where L = P_pred*ZZ'*inv(F)
##   P_upd  = P_pred - L*ZZ*P_pred   (or: P_pred*(I - ZZ'*inv(F)*ZZ*P_pred))
##   And for the NEXT prediction:
##   xi_pred_new = Tlin*xi_upd + c_drift + K*v - L*v*...
##   Actually note: the one-step ahead prediction ALSO uses the innovation:
##   xi_pred_new = Tlin*xi_upd + c_drift
##               = Tlin*(xi_pred + L*v) + c_drift
##   P_pred_new  = Tlin*P_upd*Tlin' + QQ
##   BUT this double-counts the correlation: the proper P_pred_new is
##     P_pred_new = Tlin*P_upd*Tlin' + QQ  +  (correction for SS term)
##   The correction IS zero when computed via P_upd properly, ONLY IF
##   we use the correct one-step covariance formula.
##
## DEFINITIVE FORMULA (Harvey 1990, Section 3.2, equation 3.2.3a-c):
##   Let K_t = M_t * F_t^{-1} where M_t = Tlin*P_{t|t-1}*ZZ' + SS
##   P_{t+1|t} = Tlin*P_{t|t-1}*Tlin' + QQ - K_t*M_t'   -- (*)
##   xi_{t+1|t} = Tlin*xi_{t|t-1} + c_drift + K_t*v_t
##   xi_{t|t} = xi_{t|t-1} + P_{t|t-1}*ZZ'*F_t^{-1}*v_t
##
## (*) is valid when the correlated noise enters on the SAME time index.
## This is Harvey's "general model" with S_t = SS (his notation).
##
.pruned_kf_correlated <- function(Y, Tlin, ZZ, d_y, c_drift, QQ, HH, SS,
                                  mu0, Sxi0) {
  n_obs <- nrow(ZZ)
  n_T   <- ncol(Y)
  d_dim <- nrow(Tlin)

  xi <- mu0          # d-vector: current filtered state
  P  <- Sxi0         # d x d: current filtered covariance

  loglik <- 0

  for (t in seq_len(n_T)) {
    y_t <- Y[, t]

    ## Partial-NA measurement update: only the OBSERVED components of y_t
    ## enter the update; missing components are simply dropped (not treated
    ## as an all-or-nothing period). `o` indexes the observed observables.
    o <- which(!is.na(y_t))

    ## Fully-missing period: no measurement update at all; predict forward.
    if (length(o) == 0L) {
      xi <- as.numeric(Tlin %*% xi + c_drift)
      P  <- Tlin %*% P %*% t(Tlin) + QQ
      P  <- (P + t(P)) * 0.5
      next
    }

    n_t <- length(o)
    ll_const_t <- -0.5 * n_t * log(2 * pi)

    ZZ_o <- ZZ[o, , drop = FALSE]
    d_y_o <- d_y[o]
    HH_o <- HH[o, o, drop = FALSE]
    SS_o <- SS[, o, drop = FALSE]

    ## xi, P are carried as the one-step PREDICTION xi_{t|t-1}, P_{t|t-1}.
    ## Do NOT re-apply the predict step here: mu0/Sxi0 already equal xi_{1|0}/
    ## P_{1|0}, and the bottom of the loop produces xi_{t+1|t}/P_{t+1|t}.
    ## Predicting again would push xi/P through Tlin twice per period.
    v  <- y_t[o] - d_y_o - as.numeric(ZZ_o %*% xi)
    F  <- ZZ_o %*% P %*% t(ZZ_o) + HH_o
    F  <- (F + t(F)) * 0.5

    ## Log-det of F
    F_chol <- tryCatch(chol(F), error = function(e) NULL)
    if (is.null(F_chol)) return(-Inf)
    log_det_F <- 2 * sum(log(diag(F_chol)))
    F_inv <- chol2inv(F_chol)

    ## Log-likelihood contribution
    ll_t <- ll_const_t - 0.5 * log_det_F - 0.5 * sum(v * (F_inv %*% v))
    if (!is.finite(ll_t)) return(-Inf)
    loglik <- loglik + ll_t

    ## One-step-ahead update (Harvey 1990 eq 3.2.3a, correlated noise):
    ##   M = Tlin * P_{t|t-1} * ZZ_o' + SS_o   (d x n_t)
    ##   K = M * F_inv                          (d x n_t)  modified gain
    ##   xi_{t+1|t} = Tlin*xi_{t|t-1} + c_drift + K * v
    ##   P_{t+1|t}  = Tlin*P_{t|t-1}*Tlin' + QQ - K * M'
    M <- Tlin %*% P %*% t(ZZ_o) + SS_o      # d x n_t
    K <- M %*% F_inv                        # d x n_t
    xi <- as.numeric(Tlin %*% xi + c_drift + K %*% v)
    P  <- Tlin %*% P %*% t(Tlin) + QQ - K %*% t(M)
    P  <- (P + t(P)) * 0.5
  }

  loglik
}


# ============================================================================
# make_log_posterior dispatch (pruned likelihood)
# ============================================================================

#' Build a log-posterior function using the pruned-SS Gaussian likelihood
#'
#' Internal function called by \code{make_log_posterior} and
#' \code{make_posterior} when \code{likelihood = "pruned"}.  The likelihood
#' is the Gaussian KF on the order-2 AFVRR augmented state (\code{pruned_ss_loglik}).
#' An analytic gradient is NOT available; gradient-based optimizers fall back
#' to numerical finite differences (the closure is differentiable in principle,
#' but no tangent/adjoint pass exists yet for the augmented-state KF).
#'
#' @param model     Parsed model (from \code{parse_mod}).
#' @param data      Observation matrix (T x n_obs or n_obs x T).
#' @param prior_spec Prior spec data.frame from \code{prior_spec}.
#' @param obs_vars  Character vector of observed variable names.
#' @param compiled  dynhr_compiled (from \code{compile_model}).
#' @param me_variance Measurement-error jitter (default 0).
#' @param system_priors Named list of system-prior functions (default NULL).
#' @param power Power-posterior (generalised-Bayes) tempering exponent
#'   \eqn{\zeta}. \code{NULL} (default) resolves the \code{power_posterior}
#'   option once, at factory time (see \code{.resolve_power_posterior}); the
#'   default of 1 is bit-identical to the untempered posterior.
#' @return A closure \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}.
#' @noRd
make_log_posterior_pruned <- function(model, data, prior_spec, obs_vars,
                                       compiled, me_variance = 0,
                                       system_priors = NULL,
                                       power = NULL) {
  ## Resolve zeta ONCE here, not per draw -- the exponent is a property of the
  ## closure (see .resolve_power_posterior).
  power <- .resolve_power_posterior(power, "make_log_posterior_pruned")

  ## Data: ensure n_obs x T
  n_obs <- length(obs_vars)
  Y <- if (is.matrix(data) && nrow(data) == n_obs) data else t(data)

  ## Adapter over the shared closure builder (R/posterior-closure.R). What is
  ## specific to this branch: the order-2 lift + pruned-SS assembly (solve
  ## hook), the augmented-state Gaussian KF (loglik hook), NO cold retry after
  ## a failed warm-started steady-state solve, NO stationarity guard beyond
  ## the BK check, and a system prior that is added to $logpost WITHOUT
  ## folding into $logprior (mode "extra") against dr2 and the pruned
  ## Sigma_e -- all pinned by test-posterior-closure-parity.R.
  .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    solve_fn = function(model, compiled, sys_cache, ss, params, theta) {
      s1 <- .posterior_solve1(model, compiled, sys_cache, ss, params, "none")
      if (is.null(s1)) return(NULL)
      dr2 <- tryCatch(
        solve_perturbation_order2(model, compiled, ss, params, s1$dr,
                                  verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(dr2)) return(NULL)
      pss <- tryCatch(pruned_state_space(dr2, model, params),
                      error = function(e) NULL)
      if (is.null(pss)) return(NULL)
      list(dr = dr2, pss = pss)
    },
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      loglik <- tryCatch(
        pruned_ss_loglik(sol$pss, Y, obs_vars, me_variance = me_variance,
                         me_floor_check = me_floor_check),
        error = function(e) -Inf
      )
      if (!is.finite(loglik)) return(NULL)
      list(loglik = loglik, Sigma_e = sol$pss$Sigma_e)
    },
    power             = power,
    warm_retry        = FALSE,
    system_prior      = system_priors,
    system_prior_mode = "extra")
}


## ============================================================================
## Sigma_eta DERIVATION NOTE (for orchestrator verification)
## ============================================================================
##
## The augmented state xi_t = [x1_t; x2_t; x1_t (x) x1_t] evolves as:
##
##   xi_{t+1} = Tlin * xi_t + c_drift + G * r_t         (i)
##   y_t      = Dxi  * xi_t + d_y     + Gv * r_t        (ii)
##
## where the "raw innovation" vector is
##   r_t = [ eps_t ; eps_t(x)x1_t ; x1_t(x)eps_t ; eps_t(x)eps_t ]
##
## Cr0 = E[r_t r_t'] at stationarity (a=0, P=Sigma_x):
##
## Block structure of Cr0 (Dr x Dr, Dr = n_u + n_u*n_s + n_s*n_u + n_u*n_u):
##
##   Row/col groups (using stationary a=0, M = Sigma_x):
##   j1 = 1:n_u          (eps)
##   j2 = n_u+1 : n_u+n_u*n_s      (eps(x)x1)
##   j3 = n_u+n_u*n_s+1 : n_u+n_u*n_s+n_s*n_u  (x1(x)eps)
##   j4 = rest           (eps(x)eps, centered)
##
##   Cr0[j1, j1] = Sigma_e
##   Cr0[j1, j2] = kron(Sigma_e, a') = 0  (a=0)
##   Cr0[j1, j3] = kron(a', Sigma_e) = 0  (a=0)
##   Cr0[j2, j2] = kron(Sigma_e, M)       (= kron(Sigma_e, Sigma_x))
##   Cr0[j2, j3] = C_ex1_x1e  where C_ex1_x1e[(i-1)*n_s+j, (k-1)*n_u+l]
##                             = Sigma_e[i,l] * M[j,k]
##                             (Cov of (eps_i x1_j), (x1_k eps_l))
##   Cr0[j3, j3] = kron(M, Sigma_e)       (= kron(Sigma_x, Sigma_e))
##   Cr0[j4, j4] = E[(eps(x)eps)(eps(x)eps)'] - vec(Sigma_e)vec(Sigma_e)'
##               = .fourth_moment_gaussian(Sigma_e) - outer(vecSe, vecSe)
##
##   Isserlis theorem (Wick contractions) for 4th central moment of Gaussian:
##   E[eps_i eps_j eps_k eps_l] = Sigma_e[i,j]*Sigma_e[k,l]
##                               + Sigma_e[i,k]*Sigma_e[j,l]
##                               + Sigma_e[i,l]*Sigma_e[j,k]
##   So the CENTERED 4th moment (the [j4,j4] block) is:
##   E[eps(x)eps * eps'(x)eps'] - vecSe*vecSe'
##   = kron(Sigma_e, Sigma_e) + comm(n_u, n_u) * kron(Sigma_e, Sigma_e)
##   where comm(n_u, n_u) is the vec-permutation (commutation) matrix
##   [equivalently: .fourth_moment_gaussian(Sigma_e) - outer(vecSe, vecSe)]
##
## Given Cr0, the augmented state noise covariance is:
##   QQ = G * Cr0 * G'   (d x d)
##
## with G (d x Dr):
##   G[ix1, j1] = hu            (n_s x n_u)
##   G[ix2, j2] = hxu           (n_s x n_u*n_s)
##   G[ix2, j4] = 0.5*huu       (n_s x n_u^2)
##   G[ik,  j2] = kron(hu, hx)  (n_s^2 x n_u*n_s)
##   G[ik,  j3] = kron(hx, hu)  (n_s^2 x n_s*n_u)
##   G[ik,  j4] = kron(hu, hu)  (n_s^2 x n_u^2)
##
## Blocks of QQ = G*Cr0*G':
##
## (1) [ix1, ix1] = hu * Sigma_e * hu'
##     = the standard first-order state noise (same as in the linear KF)
##
## (2) [ix2, ix2]: the x2 layer innovation covariance.
##     G[ix2, j2] * Cr0[j2,j2] * G[ix2,j2]' + G[ix2,j2]*Cr0[j2,j3]*G[ix2,j3]'
##     + G[ix2,j3]*Cr0[j3,j2]*G[ix2,j2]' [but G[ix2,j3]=0]
##     + G[ix2,j4]*Cr0[j4,j4]*G[ix2,j4]'
##   = hxu*kron(Se,Sx)*hxu' + 0.5*huu*Cr0[j4,j4]*0.5*huu'
##     (plus hxu*Cr0[j2,j4]*huu' = 0 since Cr0[j2,j4] = 0)
##
## (3) [ik, ik]: the x1(x)x1 layer innovation covariance.
##     G[ik,j2]*kron(Se,Sx)*G[ik,j2]' + G[ik,j2]*C_ex1_x1e*G[ik,j3]'
##     + G[ik,j3]*C_ex1_x1e'*G[ik,j2]' + G[ik,j3]*kron(Sx,Se)*G[ik,j3]'
##     + G[ik,j4]*Cr0[j4,j4]*G[ik,j4]'
##   = kron(hu,hx)*kron(Se,Sx)*kron(hu,hx)' + 2*kron(hu,hx)*C_ex1_x1e*kron(hx,hu)'
##     + kron(hx,hu)*kron(Sx,Se)*kron(hx,hu)' + kron(hu,hu)*Cr0[j4,j4]*kron(hu,hu)'
##
##   THE ISSERLIS (4TH-MOMENT) TERMS: Cr0[j4,j4] appears in BOTH the x2
##   layer ([ix2,ix2]) and the x1(x)x1 layer ([ik,ik]).  In the x1(x)x1 block,
##   kron(hu,hu)*Cr0[j4,j4]*kron(hu,hu)' captures the variance of the quadratic
##   term huu*(eps(x)eps) -- this is the Isserlis contribution:
##     Var(huu*(eps(x)eps)) = huu*(E[(eps(x)eps)(eps(x)eps)'] - vecSe*vecSe')*huu'
##   For a scalar shock (n_u=1), Sigma_e = sigma^2,
##     Cr0[j4,j4] = E[eps^4] - sigma^4 = 3*sigma^4 - sigma^4 = 2*sigma^4
##   which is the variance of eps^2 (since Var(eps^2) = 2*sigma^4 for N(0,sigma^2)).
##
## (4) Cross blocks [ix1,ix2], [ix1,ik], [ix2,ik]: also non-zero in general;
##   all are computed by the matrix product G*Cr0*G' which handles them
##   automatically.  The [ix1,ix2] cross (first-order/second-order covariance)
##   is hu*Cr0[j1,j2]*hxu' = hu*0*hxu' = 0  (because Cr0[j1,j2]=kron(Se,a')=0
##   at stationarity with a=0).  Similarly [ix1,ik] = 0 at stationarity.
##   The [ix2,ik] cross is non-zero: it involves kron(hu,hx)*C_ex1_x1e etc.
##
## SUMMARY: QQ = G*Cr0*G' with Cr0 computed above at a=0, P=Sigma_x.
## The x1(x)x1 diagonal block includes Isserlis 4th-moment terms from
## kron(hu,hu)*Cr0[j4,j4]*kron(hu,hu)' and G[ik,j2]*kron(Se,Sx)*G[ik,j2]'
## (the "eps(x)x1" second-moment contribution to the kronecker layer).
## All cross-blocks are computed automatically by the G*Cr0*G' product.
##
## ============================================================================
