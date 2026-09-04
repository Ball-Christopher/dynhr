## R/method-of-moments.R
## --------------------------------------------------------------------------
## GMM / SMM estimation by moment matching (E2-A).
##
## The cumulant machinery in R/cumulant-likelihood.R already builds a
## CONTEMPORANEOUS model-implied moment map -- [mean; vec(Sigma_y); vec(c3);
## vec(c4)] -- and a Newey-West / analytic long-run covariance for it. What it
## did NOT have was
##
##   * autocovariances Gamma(h) = Cov(y_t, y_{t-h}) for h >= 1, which are the
##     moments that identify the PERSISTENCE parameters (an AR(1) rho is
##     invisible to the contemporaneous variance alone once sigma is free),
##   * a simulated (SMM) moment map,
##   * a J-test, sandwich standard errors, or a user-facing entry point.
##
## This file adds those four. It deliberately reuses, rather than re-derives:
##   * `.build_moment_vector()` (R/cumulant-gradient.R) for the contemporaneous
##     block -- extended here with a `lags` argument so that ONE function
##     defines the moment ordering for the loglik, the weight matrix and the
##     estimator,
##   * `compute_moments()` for Sigma_y / Sigma_state,
##   * `solution_derivatives()` for the analytic moment Jacobian,
##   * `.run_mode_finding()` (R/mode-orchestrate.R) for the optimiser,
##   * `.with_local_seed()` (R/tpf-likelihood.R) for SMM common random numbers.
##
## MOMENT ORDERING (one definition, used everywhere):
##   [ mean(obs) ]                                  if 1 %in% orders
##   [ vec(Sigma_y[obs, obs]) ]                     if 2 %in% orders
##   [ vec(c3) ]                                    if 3 %in% orders
##   [ vec(c4) ]                                    if 4 %in% orders
##   [ vec(Gamma(h1)) ... vec(Gamma(hH)) ]          for the lags >= 1, ascending
##
## AUTOCOVARIANCE RECURSION. In the lagged-state convention the decision rule
## is y_t = ghx s_{t-1} + ghu e_t with s = y[state_idx], so with the endo-space
## transition operator G = ghx %*% S (S the n_state x n_endo selection matrix)
##
##   Gamma(h) = Cov(y_t, y_{t-h}) = G Gamma(h-1),   Gamma(0) = Sigma_y.
##
## This is the SAME recursion compute_moments() uses for `autocorr`, and it is
## the correct one: the alternative closed form `Z hx^h Sigma_state Z'` DROPS
## the Cov(s_{t-1}, e_{t-h}) channel and is only valid when the observables
## carry no contemporaneous shock loading (ghu[obs, ] == 0).
## `.analytic_gmm_weight_matrix()` (R/cumulant-likelihood.R) uses this SAME
## recursion via `.mom_endo_transition()`; it used to carry the dropped-channel
## closed form, which made its Omega wrong at every lag h >= 1 on any model
## with a contemporaneously-loaded observable.
## --------------------------------------------------------------------------


# ============================================================================
# Moment-spec helpers
# ============================================================================

#' Normalise a `lags` specification to the strictly-positive lags, ascending.
#'
#' Lag 0 is the CONTEMPORANEOUS block and is expressed through `orders`
#' (order 2 == vec(Sigma_y)), so `lags = 0:4` and `lags = 1:4` describe the
#' same autocovariance block; `lags = 0` means "no autocovariances".
#' @noRd
.mom_lags <- function(lags) {
  if (is.null(lags) || length(lags) == 0L) return(integer(0))
  lg <- as.integer(lags)
  if (any(is.na(lg))) stop("`lags` must be a finite integer vector.",
                           call. = FALSE)
  if (any(lg < 0L)) stop("`lags` must be non-negative.", call. = FALSE)
  sort(unique(lg[lg >= 1L]))
}

#' Resolve and validate a `moments = list(orders =, lags =)` specification.
#' @noRd
.mom_spec <- function(moments) {
  if (is.null(moments)) moments <- list()
  if (!is.list(moments))
    stop("`moments` must be a list with elements `orders` and `lags`.",
         call. = FALSE)
  orders <- as.integer(moments$orders %||% 1:2)
  lags   <- .mom_lags(moments$lags %||% 0L)
  if (!length(orders) || any(is.na(orders)))
    stop("`moments$orders` must be a non-empty integer vector.", call. = FALSE)
  orders <- sort(unique(orders))
  if (any(!orders %in% 1:4))
    stop("`moments$orders` entries must be in 1:4.", call. = FALSE)
  ## Order 4 used to be refused here: the model side was n_obs x n_endo^3
  ## against an n_obs x n_obs^3 sample k-statistic, so the two only lined up
  ## when every endogenous variable was observed.  `.build_moment_vector()`
  ## now projects the model c4 onto the observables (.project_c4_obs), so the
  ## blocks match for any n_obs <= n_endo and order 4 is allowed.
  if (!length(orders) && !length(lags))
    stop("`moments` selects no moments at all.", call. = FALSE)
  list(orders = orders, lags = lags)
}

#' Human-readable names for every entry of the moment vector.
#' @noRd
.mom_moment_names <- function(obs_vars, orders, lags) {
  n <- length(obs_vars)
  nm <- character(0)
  if (1L %in% orders) nm <- c(nm, paste0("mean[", obs_vars, "]"))
  if (2L %in% orders) {
    ## column-major vec(): (i, j) at (j - 1) * n + i
    nm <- c(nm, as.vector(outer(obs_vars, obs_vars,
                                function(a, b) paste0("var[", a, ",", b, "]"))))
  }
  if (3L %in% orders) {
    ## row = i, col = (j - 1) * n + k, then as.numeric() (column-major)
    lab <- matrix("", n, n * n)
    for (i in seq_len(n)) for (j in seq_len(n)) for (k in seq_len(n))
      lab[i, (j - 1L) * n + k] <-
        paste0("c3[", obs_vars[i], ",", obs_vars[j], ",", obs_vars[k], "]")
    nm <- c(nm, as.vector(lab))
  }
  if (4L %in% orders) {
    ## row = i, col = (j - 1) * n^2 + (k - 1) * n + l, then as.numeric()
    lab <- matrix("", n, n * n * n)
    for (i in seq_len(n)) for (j in seq_len(n)) for (k in seq_len(n))
      for (l in seq_len(n))
        lab[i, (j - 1L) * n * n + (k - 1L) * n + l] <-
          paste0("c4[", obs_vars[i], ",", obs_vars[j], ",", obs_vars[k], ",",
                 obs_vars[l], "]")
    nm <- c(nm, as.vector(lab))
  }
  for (h in lags)
    nm <- c(nm, as.vector(outer(obs_vars, obs_vars, function(a, b)
      paste0("acov", h, "[", a, ",", b, "]"))))
  nm
}


# ============================================================================
# Model-implied autocovariances
# ============================================================================

#' Endo-space transition operator G with Gamma(h) = G Gamma(h-1).
#' @noRd
.mom_endo_transition <- function(dr) {
  n_endo <- length(dr$endo_names)
  st     <- dr$state_idx
  S <- matrix(0, length(st), n_endo)
  if (length(st)) S[cbind(seq_along(st), st)] <- 1
  dr$ghx %*% S
}

#' Model-implied observable autocovariances Gamma(h), h >= 1
#'
#' @param dr        Decision rules (order >= 1); only the first-order blocks
#'   enter, so an order-2 rule gives the same Gamma(h) as its order-1 core.
#' @param model,params As in \code{compute_moments()}.
#' @param obs_vars  Observable names.
#' @param lags      Strictly-positive lags (already normalised).
#' @param moments   Optional cached \code{compute_moments()} result.
#' @return Named list of n_obs x n_obs matrices, one per lag, or NULL when the
#'   stationary covariance is not finite (unit root / explosive draw).
#' @noRd
.mom_model_autocov <- function(dr, model, params, obs_vars, lags,
                               moments = NULL) {
  lags <- .mom_lags(lags)
  if (!length(lags)) return(list())
  mom <- moments %||% compute_moments(dr, model, params = params)
  Sig <- mom$var_cov
  obs_idx <- match(obs_vars, dr$endo_names)
  if (any(is.na(obs_idx))) return(NULL)
  if (!all(is.finite(Sig))) return(NULL)
  G <- .mom_endo_transition(dr)
  out <- vector("list", length(lags))
  names(out) <- paste0("lag", lags)
  Gam <- Sig
  done <- 0L
  for (i in seq_along(lags)) {
    while (done < lags[i]) {
      Gam  <- G %*% Gam
      done <- done + 1L
    }
    out[[i]] <- Gam[obs_idx, obs_idx, drop = FALSE]
    dimnames(out[[i]]) <- list(obs_vars, obs_vars)
  }
  out
}

#' Sample autocovariances Gamma_hat(h) = (T-h)^{-1} sum_t yc_t yc_{t-h}'
#'
#' The mean subtracted is the FULL-sample mean (the same one
#' \code{sample_cumulants()} uses), so the lag-0 element of this family and
#' the order-2 cumulant moment are computed off one centring.
#' @noRd
.mom_sample_autocov <- function(data, obs_vars, lags) {
  lags <- .mom_lags(lags)
  if (!length(lags)) return(list())
  Y  <- data[, obs_vars, drop = FALSE]
  T_obs <- nrow(Y)
  if (max(lags) >= T_obs)
    stop("method_of_moments(): requested lag ", max(lags),
         " with only ", T_obs, " observations.", call. = FALSE)
  Yc <- sweep(Y, 2L, colMeans(Y), FUN = "-")
  out <- vector("list", length(lags))
  names(out) <- paste0("lag", lags)
  for (i in seq_along(lags)) {
    h <- lags[i]
    A <- Yc[(h + 1L):T_obs, , drop = FALSE]
    B <- Yc[1L:(T_obs - h), , drop = FALSE]
    out[[i]] <- crossprod(A, B) / (T_obs - h)
    dimnames(out[[i]]) <- list(obs_vars, obs_vars)
  }
  out
}

#' Sample moment vector in the canonical ordering.
#' @noRd
.mom_sample_moment_vector <- function(data, obs_vars, orders, lags) {
  Y <- data[, obs_vars, drop = FALSE]
  m <- numeric(0)
  if (any(orders >= 1L)) {
    sc <- sample_cumulants(Y, max_order = max(2L, max(orders)))
    if (1L %in% orders) m <- c(m, unname(sc$mean[obs_vars]))
    if (2L %in% orders)
      m <- c(m, as.numeric(sc$var_cov[obs_vars, obs_vars, drop = FALSE]))
    if (3L %in% orders) m <- c(m, as.numeric(sc$c3))
    if (4L %in% orders) m <- c(m, as.numeric(sc$c4))
  }
  for (A in .mom_sample_autocov(Y, obs_vars, lags)) m <- c(m, as.numeric(A))
  m
}


# ============================================================================
# Per-observation moment contributions and the long-run covariance
# ============================================================================

#' T_eff x p matrix of per-observation raw moment contributions.
#'
#' Rows are the ALIGNED sample t = (H+1)..T with H = max(lags), so every
#' block is evaluated on one common index set (an autocovariance contribution
#' is undefined for t <= h). With no lags, H = 0 and the full sample is used,
#' which is what keeps the lags-free path identical to
#' \code{estimate_gmm_weight_matrix()}.
#' @noRd
.mom_contributions <- function(data, obs_vars, orders, lags) {
  Y  <- data[, obs_vars, drop = FALSE]
  n  <- length(obs_vars)
  T_obs <- nrow(Y)
  lags  <- .mom_lags(lags)
  H  <- if (length(lags)) max(lags) else 0L
  Yc <- sweep(Y, 2L, colMeans(Y), FUN = "-")
  idx <- (H + 1L):T_obs
  Te  <- length(idx)

  blocks <- list()
  if (1L %in% orders) blocks[[length(blocks) + 1L]] <- Y[idx, , drop = FALSE]
  if (2L %in% orders) {
    B <- matrix(0, Te, n * n)
    for (r in seq_len(Te)) {
      yc <- Yc[idx[r], ]
      B[r, ] <- as.numeric(outer(yc, yc))
    }
    blocks[[length(blocks) + 1L]] <- B
  }
  if (3L %in% orders) {
    B <- matrix(0, Te, n * n * n)
    for (r in seq_len(Te)) {
      yc <- Yc[idx[r], ]
      m3 <- matrix(0, n, n * n)
      for (a in seq_len(n)) for (b in seq_len(n)) for (cc in seq_len(n))
        m3[a, (b - 1L) * n + cc] <- yc[a] * yc[b] * yc[cc]
      B[r, ] <- as.numeric(m3)
    }
    blocks[[length(blocks) + 1L]] <- B
  }
  if (4L %in% orders) {
    ## The order-4 MOMENT is the fourth CUMULANT k4, not the raw fourth
    ## central moment, so the per-observation contribution must carry the
    ## Gaussian-pairing subtraction too -- otherwise .mom_omega() centres a
    ## raw-moment contribution on a cumulant and the long-run covariance is
    ## off by the 3 sigma^4 pairing.  Deterministic given the sample, so it is
    ## computed once, outside the row loop.
    Vh   <- crossprod(Yc) / max(1L, nrow(Yc) - 1L)     # n x n sample covariance
    pair <- numeric(n * n * n * n)
    for (i in seq_len(n)) for (j in seq_len(n))
      for (k in seq_len(n)) for (l in seq_len(n)) {
        ## flattened position of (row = i, col = (j-1)n^2 + (k-1)n + l), the
        ## sample_cumulants()$c4 layout under as.numeric().
        pos <- ((j - 1L) * n * n + (k - 1L) * n + l - 1L) * n + i
        pair[pos] <- Vh[i, j] * Vh[k, l] + Vh[i, k] * Vh[j, l] +
                     Vh[i, l] * Vh[j, k]
      }
    B <- matrix(0, Te, n * n * n * n)
    for (r in seq_len(Te)) {
      yc <- Yc[idx[r], ]
      ## y_i y_j y_k y_l is fully symmetric, so the 4-fold outer product
      ## flattened column-major IS the (i; j,k,l) layout above.
      B[r, ] <- as.numeric(outer(yc, outer(yc, outer(yc, yc)))) - pair
    }
    blocks[[length(blocks) + 1L]] <- B
  }
  for (h in lags) {
    B <- matrix(0, Te, n * n)
    for (r in seq_len(Te)) {
      B[r, ] <- as.numeric(outer(Yc[idx[r], ], Yc[idx[r] - h, ]))
    }
    blocks[[length(blocks) + 1L]] <- B
  }
  do.call(cbind, blocks)
}

#' Newey-West long-run covariance of the moment conditions (with lags).
#'
#' @param m_center p-vector the raw contributions are centred on (the
#'   model-implied moment vector at the current estimate; column means when
#'   NULL).
#' @return list(Omega, Omega_reg, bandwidth, T_eff, p)
#' @noRd
.mom_omega <- function(data, obs_vars, orders, lags, m_center = NULL,
                       bandwidth = NULL, ridge = 1e-6) {
  G_raw <- .mom_contributions(data, obs_vars, orders, lags)
  Te <- nrow(G_raw)
  p  <- ncol(G_raw)
  ctr <- if (is.null(m_center)) colMeans(G_raw) else m_center
  if (length(ctr) != p)
    stop(".mom_omega: centre has length ", length(ctr), " but there are ",
         p, " moments.", call. = FALSE)
  Gc <- sweep(G_raw, 2L, ctr, FUN = "-")

  if (is.null(bandwidth)) bandwidth <- floor(4 * (Te / 100)^(2 / 9))
  bandwidth <- as.integer(max(0L, bandwidth))

  Om <- crossprod(Gc) / Te
  if (bandwidth > 0L) {
    for (j in seq_len(bandwidth)) {
      if (j >= Te) break
      w  <- 1 - j / (bandwidth + 1)
      Gj <- crossprod(Gc[(j + 1L):Te, , drop = FALSE],
                      Gc[1L:(Te - j), , drop = FALSE]) / Te
      Om <- Om + w * (Gj + t(Gj))
    }
  }
  Om <- (Om + t(Om)) * 0.5

  dmax <- max(diag(Om))
  if (!is.finite(dmax) || dmax <= 0) dmax <- 1
  ## Omega_reg is only for the DIAGONAL weighting (whose 1/diag would divide
  ## by zero on a degenerate moment); the full inverse goes through
  ## .mom_omega_inverse(), which does its own scaled eigendecomposition, so no
  ## condition number is computed here -- that would be a second p x p eigen()
  ## per weight-matrix build for a number nothing reads.
  Om_reg <- Om + ridge * dmax * diag(p)
  list(Omega = Om, Omega_reg = Om_reg, bandwidth = bandwidth,
       T_eff = Te, p = p)
}

#' GMM criterion T * d' W d.
#'
#' The identity branch is written out rather than routed through the quadratic
#' form on purpose. `.cumulant_loglik()`'s identity path evaluates
#' `-0.5 * sum(delta^2) / length(delta) * T_obs`; the general quadratic form
#' associates the same products differently and lands one ulp away. Because
#' `-2 *` and `-0.5 *` are exact in binary floating point, computing the
#' identity criterion as `-2 *` that expression makes `-0.5 * criterion`
#' reproduce the cumulant log-likelihood BIT-FOR-BIT, which is the documented
#' contract (and the E2-A oracle (i) gate).
#' @noRd
.mom_criterion_value <- function(d, W, T_obs) {
  if (identical(attr(W, "method"), "identity"))
    return(-2 * (-0.5 * sum(d^2) / length(d) * T_obs))
  T_obs * as.numeric(crossprod(d, W %*% d))
}

#' Rank-truncated inverse of a long-run covariance, plus its numerical rank.
#'
#' The moment vector uses \code{vec(Sigma_y)}, not \code{vech}, so a symmetric
#' block contributes EXACTLY duplicated moments (\code{var[a,b]} and
#' \code{var[b,a]}). That is deliberate -- it is what makes the lags-free
#' moment vector identical to \code{.cumulant_loglik}'s -- but it makes
#' \eqn{\Omega} exactly singular, and a ridge-only inverse then puts a
#' \eqn{1/\mathrm{ridge}}-sized weight on a direction carrying no information
#' and, worse, reports \eqn{p - k} degrees of freedom for the J-test when the
#' true number of independent moments is smaller.  So: eigendecompose the
#' CORRELATION-scaled \eqn{\Omega} (scale-free, so the tolerance means the same
#' thing for a mean moment and a fourth-moment), drop the numerically-null
#' directions, and report the retained rank for the J-test's df.
#'
#' \eqn{\Omega = D C D \Rightarrow \Omega^{+} = D^{-1} C^{+} D^{-1}}.
#' @noRd
.mom_omega_inverse <- function(Omega, ridge = 1e-6, tol = 1e-9) {
  p <- nrow(Omega)
  d <- sqrt(pmax(diag(Omega), 0))
  d[!is.finite(d) | d <= 0] <- 1
  C <- Omega / outer(d, d)
  eg <- eigen((C + t(C)) * 0.5, symmetric = TRUE)
  lmax <- max(eg$values)
  if (!is.finite(lmax) || lmax <= 0)
    stop(".mom_omega_inverse: long-run covariance is not positive.",
         call. = FALSE)
  keep <- eg$values > tol * lmax
  V <- eg$vectors[, keep, drop = FALSE]
  lam <- eg$values[keep] + ridge * lmax
  Cinv <- V %*% diag(1 / lam, nrow = length(lam)) %*% t(V)
  W <- Cinv / outer(d, d)
  W <- (W + t(W)) * 0.5
  list(W = W, rank = sum(keep), p = p,
       condition_number = max(eg$values) / min(eg$values[keep]))
}

#' Symmetric inverse via Cholesky, with a pseudo-inverse fallback.
#' @noRd
.mom_sym_inverse <- function(M, what = "matrix") {
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(ch)) {
    Mi <- chol2inv(ch)
    return((Mi + t(Mi)) * 0.5)
  }
  eg <- eigen(M, symmetric = TRUE)
  tol <- max(dim(M)) * max(abs(eg$values)) * .Machine$double.eps
  keep <- abs(eg$values) > tol
  if (!any(keep))
    stop(".mom_sym_inverse: ", what, " is numerically zero.", call. = FALSE)
  V <- eg$vectors[, keep, drop = FALSE]
  Mi <- V %*% diag(1 / eg$values[keep], nrow = sum(keep)) %*% t(V)
  (Mi + t(Mi)) * 0.5
}


# ============================================================================
# Model solve at a candidate theta
# ============================================================================

#' Solve steady state + perturbation at `params`, returning NULL on failure.
#' @noRd
.mom_solve_at <- function(model, compiled, params, order = 1L) {
  ss <- tryCatch(solve_steady_state(model, compiled, params, verbose = FALSE),
                 error = function(e) NULL)
  if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)
  params <- ss$params %||% params
  ## suppressWarnings: an optimiser step INTO the non-stationary region is a
  ## normal event for a criterion that is +Inf there, and solve_perturbation()
  ## warns loudly on a Blanchard-Kahn violation. The rejection is handled by
  ## the `bk_satisfied` / spectral-radius tests below, so the warning would be
  ## pure noise -- hundreds of them per fit.
  dr <- suppressWarnings(tryCatch(
    solve_perturbation(model, compiled, ss$values, params,
                       order = as.integer(order), verbose = FALSE),
    error = function(e) NULL))
  if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)
  ## Stationarity: the whole moment map is a stationary-distribution object.
  sr <- tryCatch(
    max(Mod(eigen(dr$ghx[dr$state_idx, , drop = FALSE],
                  only.values = TRUE)$values)),
    error = function(e) Inf)
  if (!is.finite(sr) || sr >= 1) return(NULL)
  list(dr = dr, params = params, ss = ss$values)
}


# ============================================================================
# SMM: simulated moment map with common random numbers
# ============================================================================

#' Draw the fixed standardized innovations used by every SMM evaluation.
#'
#' Drawn ONCE, under \code{.with_local_seed()} so the caller's RNG stream is
#' restored, and reused at every theta -- the common-random-numbers condition
#' that makes the simulated criterion a smooth function of theta rather than a
#' fresh Monte-Carlo surface at each step.
#'
#' `byrow = TRUE` is load-bearing, not cosmetic: it makes period \eqn{t} draw
#' stream positions \eqn{(t-1)n_u + 1, \ldots, t n_u}, so at a fixed seed a
#' LONGER simulation is a strict EXTENSION of a shorter one rather than a
#' completely different path.  That is what lets `n_sim` be raised and the
#' remaining gap to the analytic GMM estimate read as simulation error
#' shrinking, instead of as one Monte-Carlo draw replaced by another.
#' @noRd
.mom_crn_draws <- function(n_rows, n_exo, seed) {
  .with_local_seed(seed, matrix(stats::rnorm(n_rows * n_exo), n_rows, n_exo,
                                byrow = TRUE))
}

#' Simulated moment vector at one theta (common random numbers).
#' @noRd
.mom_simulated_moments <- function(dr, model, params, obs_vars, orders, lags,
                                   U, n_periods, burn_in) {
  exo <- dr$exo_names
  Sig <- .get_shock_cov(model, exo, params)
  L <- tryCatch(t(chol(Sig)), error = function(e) NULL)
  if (is.null(L)) L <- .tpf_psd_sqrt(Sig)
  shocks <- U %*% t(L)
  sim <- tryCatch(
    .simulate_dr_any_order(dr, n_periods = n_periods, model = model,
                           burn_in = burn_in, shocks = shocks),
    error = function(e) NULL)
  if (is.null(sim)) return(NULL)
  lev <- attr(sim, "levels")
  if (is.null(lev)) return(NULL)
  colnames(lev) <- dr$endo_names
  if (!all(is.finite(lev[, obs_vars]))) return(NULL)
  .mom_sample_moment_vector(lev, obs_vars, orders, lags)
}


# ============================================================================
# Moment Jacobian dm(theta)/dtheta'
# ============================================================================

#' Analytic moment Jacobian for orders 1-2 plus autocovariances.
#'
#' Chains \code{solution_derivatives()} (d ghx, d ghu, d ys) through
#'   Sigma_s = hx Sigma_s hx' + hu Sigma_e hu'   (Lyapunov; differentiated by
#'                                                ONE further Lyapunov solve),
#'   Sigma_y = ghx Sigma_s ghx' + ghu Sigma_e ghu',
#'   Gamma(h) = G Gamma(h-1),   dGamma(h) = dG Gamma(h-1) + G dGamma(h-1).
#'
#' Parameters that are estimated shock standard deviations (a name in
#' \code{model$varexo_names} rather than in \code{model$param_values}) enter
#' ONLY through Sigma_e, so their solution derivatives are exactly zero and
#' only the d(Sigma_e) channel is formed for them.
#'
#' @return list(J, ok) with J the p x k Jacobian, or \code{ok = FALSE} plus a
#'   \code{reason} when the analytic route does not cover the request.
#' @noRd
.mom_jacobian_analytic <- function(model, compiled, dr, params, param_names,
                                   obs_vars, orders, lags, me_variance,
                                   h_rel = 1e-6) {
  fail <- function(reason) list(ok = FALSE, reason = reason)

  if (any(orders >= 3L))
    return(fail("cumulant orders >= 3 have no closed-form moment Jacobian"))
  if (.dr_perturbation_order(dr) > 1L)
    return(fail(paste0("the order-2 mean correction 0.5 * ghss has no ",
                       "closed-form parameter derivative here")))

  endo <- dr$endo_names
  exo  <- dr$exo_names
  n_endo <- length(endo)
  n_exo  <- length(exo)
  st     <- dr$state_idx
  n_s    <- length(st)
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx))) return(fail("obs_vars not found in the model"))

  ghx <- dr$ghx
  ghu <- dr$ghu
  hx  <- ghx[st, , drop = FALSE]
  hu  <- ghu[st, , drop = FALSE]
  Sig_e <- .get_shock_cov(model, exo, params)
  Sig_s <- .state_covariance(hx, hu, Sig_e)
  if (!all(is.finite(Sig_s))) return(fail("stationary state covariance is not finite"))
  Sig_y <- ghx %*% Sig_s %*% t(ghx) + ghu %*% Sig_e %*% t(ghu)
  Gop   <- .mom_endo_transition(dr)

  ## Split structural parameters (which move the solution) from estimated
  ## shock stds (which move only Sigma_e).
  struct <- param_names[param_names %in% names(model$param_values)]
  sd_par <- setdiff(param_names, struct)
  if (!all(sd_par %in% (model$varexo_names %||% character(0))))
    return(fail(paste0("parameter(s) ",
                       paste(setdiff(sd_par, model$varexo_names %||% character(0)),
                             collapse = ", "),
                       " are neither model parameters nor declared shocks")))

  sdv <- NULL
  if (length(struct)) {
    sdv <- tryCatch(
      suppressWarnings(
        solution_derivatives(model, compiled, dr, params, struct,
                             obs_vars = obs_vars, h_rel = h_rel)),
      error = function(e) NULL)
    if (is.null(sdv)) return(fail("solution_derivatives() failed"))
    bad <- struct[!vapply(struct, function(p) isTRUE(sdv$derivs[[p]]$ok),
                          logical(1))]
    if (length(bad))
      return(fail(paste0("solution_derivatives() failed for ",
                         paste(bad, collapse = ", "))))
  }

  lags <- .mom_lags(lags)
  p_len <- length(.mom_moment_names(obs_vars, orders, lags))
  J <- matrix(0, p_len, length(param_names),
              dimnames = list(.mom_moment_names(obs_vars, orders, lags),
                              param_names))

  for (j in seq_along(param_names)) {
    pnm <- param_names[j]
    if (pnm %in% struct) {
      d   <- sdv$derivs[[pnm]]
      dG  <- d$dG
      dH  <- d$dH
      dys <- d$dys
    } else {
      dG  <- matrix(0, n_endo, n_s)
      dH  <- matrix(0, n_endo, n_exo)
      dys <- setNames(numeric(n_endo), endo)
    }
    hstep <- max(h_rel * abs(params[[pnm]]), 1e-7)
    dSe <- .o2sd_dSigma_e(model, params, pnm, hstep, n_exo)

    dhx <- dG[st, , drop = FALSE]
    dhu <- dH[st, , drop = FALSE]

    ## dSigma_s solves the SAME Lyapunov operator with the differentiated RHS.
    dQ  <- dhu %*% Sig_e %*% t(hu) + hu %*% dSe %*% t(hu) +
           hu %*% Sig_e %*% t(dhu)
    rhs <- dhx %*% Sig_s %*% t(hx) + hx %*% Sig_s %*% t(dhx) + dQ
    dSs <- solve_lyapunov(hx, (rhs + t(rhs)) * 0.5)
    if (!all(is.finite(dSs)))
      return(fail("the differentiated Lyapunov solve is not finite"))

    dSy <- dG %*% Sig_s %*% t(ghx) + ghx %*% dSs %*% t(ghx) +
           ghx %*% Sig_s %*% t(dG) +
           dH %*% Sig_e %*% t(ghu) + ghu %*% dSe %*% t(ghu) +
           ghu %*% Sig_e %*% t(dH)

    col <- numeric(0)
    if (1L %in% orders) col <- c(col, unname(dys[obs_vars]))
    if (2L %in% orders)
      col <- c(col, as.numeric(dSy[obs_idx, obs_idx, drop = FALSE]))

    if (length(lags)) {
      dGop <- dG %*% {
        S <- matrix(0, n_s, n_endo)
        if (n_s) S[cbind(seq_len(n_s), st)] <- 1
        S
      }
      Gam  <- Sig_y
      dGam <- dSy
      done <- 0L
      for (h in lags) {
        while (done < h) {
          dGam <- dGop %*% Gam + Gop %*% dGam
          Gam  <- Gop %*% Gam
          done <- done + 1L
        }
        col <- c(col, as.numeric(dGam[obs_idx, obs_idx, drop = FALSE]))
      }
    }
    J[, j] <- col
  }
  list(ok = TRUE, J = J)
}

#' Central finite-difference moment Jacobian (full re-solve at each step).
#'
#' `moment_fn` is the SAME closure the criterion uses, so under SMM the common
#' random numbers are shared between the objective and its Jacobian.
#' @noRd
.mom_jacobian_fd <- function(theta, obs_vars, orders, lags, moment_fn,
                             h_rel = 1e-5) {
  k <- length(theta)
  m0 <- moment_fn(theta)
  if (is.null(m0)) return(NULL)
  J <- matrix(NA_real_, length(m0), k,
              dimnames = list(.mom_moment_names(obs_vars, orders, lags),
                              names(theta)))
  for (j in seq_len(k)) {
    h <- max(h_rel * abs(theta[[j]]), 1e-7)
    tp <- theta; tp[[j]] <- theta[[j]] + h
    tm <- theta; tm[[j]] <- theta[[j]] - h
    mp <- moment_fn(tp)
    mm <- moment_fn(tm)
    if (is.null(mp) || is.null(mm)) next
    J[, j] <- (mp - mm) / (2 * h)
  }
  J
}

#' Moment Jacobian with the analytic route first and a REPORTED FD fallback.
#'
#' The fallback is reported THREE ways -- a one-time warning, the returned
#' `$fell_back` / `$reason`, and the `print()` line -- because a silently
#' finite-differenced Jacobian degrades the standard errors without changing
#' anything a caller would otherwise notice.
#' @noRd
.mom_jacobian <- function(model, compiled, dr, params, theta, obs_vars,
                          orders, lags, me_variance, moment_fn,
                          jacobian = c("auto", "analytic", "fd"),
                          h_rel = 1e-6) {
  jacobian <- match.arg(jacobian)
  if (jacobian != "fd") {
    an <- .mom_jacobian_analytic(model, compiled, dr, params, names(theta),
                                 obs_vars, orders, lags, me_variance,
                                 h_rel = h_rel)
    if (isTRUE(an$ok))
      return(list(J = an$J, method = "analytic", fell_back = FALSE,
                  reason = NA_character_))
    if (jacobian == "analytic")
      stop("method_of_moments(jacobian = \"analytic\"): ", an$reason,
           call. = FALSE)
    msg <- paste0("method_of_moments(): the analytic moment Jacobian is not ",
                  "available here (", an$reason, "); falling back to central ",
                  "finite differences of the full re-solve.")
    .cumulant_warn_once("mom_jacobian_fd_fallback", msg)
    J <- .mom_jacobian_fd(theta, obs_vars, orders, lags, moment_fn,
                          h_rel = max(h_rel, 1e-5))
    return(list(J = J, method = "fd", fell_back = TRUE, reason = an$reason))
  }
  J <- .mom_jacobian_fd(theta, obs_vars, orders, lags, moment_fn,
                        h_rel = max(h_rel, 1e-5))
  list(J = J, method = "fd", fell_back = FALSE,
       reason = "requested by the caller")
}


# ============================================================================
# Main entry point
# ============================================================================

#' GMM / SMM estimation by moment matching
#'
#' Estimates structural parameters by minimising the GMM criterion
#' \deqn{Q(\theta) = T \, g(\theta)' W g(\theta), \qquad
#'       g(\theta) = \hat m - m(\theta)}
#' over a moment vector that combines the contemporaneous cumulants of
#' \code{\link{sample_cumulants}} (mean, variance, third cumulant) with the
#' autocovariances \eqn{\Gamma(h) = \mathrm{Cov}(y_t, y_{t-h})} for the
#' requested lags.  \code{method = "gmm"} uses the closed-form model-implied
#' moments; \code{method = "smm"} replaces them with simulated moments under
#' common random numbers.
#'
#' @section Relation to the cumulant likelihood:
#' With \code{moments = list(orders = 1:2, lags = 0)} and
#' \code{weighting = "identity"} the criterion is exactly
#' \eqn{-2\times} the \code{likelihood = "cumulant"} objective
#' (\code{.cumulant_loglik}), i.e. \code{-0.5 * result$criterion} reproduces
#' that log-likelihood bit-for-bit.  Non-zero \code{lags} extend the same
#' moment map; nothing else changes.
#'
#' @section Weighting:
#' \code{"identity"} (the default) uses \eqn{W = I_p / p}.  It is the default
#' deliberately: the efficient weight is asymptotically optimal but is a known
#' finite-sample pathology once \eqn{p} is large relative to \eqn{T} (the
#' too-many-moments problem, Donald & Newey 2001) -- see the
#' finite-sample caution in \code{\link{estimate_gmm_weight_matrix}}.
#' \code{"newey_west"} / \code{"optimal"} use \eqn{W = \hat\Omega^{-1}} with a
#' Newey-West HAC \eqn{\hat\Omega}; \code{"diagonal"} keeps only its diagonal.
#' \code{two_step = TRUE} always finishes with one efficient step at the
#' first-step estimate (skipped, with a warning, when \eqn{p \ge T}).
#'
#' @section Inference:
#' Standard errors come from the sandwich
#' \deqn{\widehat{\mathrm{Var}}(\hat\theta) = \tfrac{c}{T}
#'   (J'WJ)^{-1} J'W \hat\Omega W J (J'WJ)^{-1}}
#' with \eqn{J = \partial m(\theta)/\partial\theta'} and \eqn{c = 1} for GMM,
#' \eqn{c = 1 + 1/\tau} for SMM with \eqn{\tau = } \code{n_sim}.  The
#' over-identification statistic is \eqn{J = T g'\hat\Omega^{+}g} with
#' \eqn{r - k} degrees of freedom, where \eqn{r} is the NUMERICAL RANK of
#' \eqn{\hat\Omega} (\code{$omega_rank}).  The rank, not \eqn{p}, is the right
#' count: the moment vector uses \code{vec(Sigma_y)} rather than \code{vech},
#' so the symmetric contemporaneous block carries exactly duplicated entries
#' that are not independent over-identifying restrictions.
#'
#' @param model    A \code{dynhr_mod} (from \code{\link{parse_mod}}) or a
#'   \code{dynhr_solved} (from \code{\link{solve_model}}).
#' @param data     \eqn{T \times n_{obs}} matrix or data.frame whose columns
#'   include \code{obs_vars}.
#' @param obs_vars Character vector of observed variable names.  Defaults to
#'   the model's \code{varobs} declaration.
#' @param start    Named numeric vector of starting values.  Its NAMES select
#'   which quantities are estimated: a model parameter, or a declared shock
#'   (an estimated \code{stderr}).  Defaults to the prior means of the
#'   \code{estimated_params} block when the model has one.
#' @param moments  \code{list(orders = , lags = )}.  \code{orders} selects the
#'   contemporaneous cumulant orders (1 = mean, 2 = variance, 3 = third
#'   cumulant, 4 = fourth cumulant).  Orders 3-4 need an order-2 perturbation
#'   and add \eqn{n_{obs}^3} / \eqn{n_{obs}^4} moments each, so they are cheap
#'   only for small \eqn{n_{obs}}.  \code{lags} selects the autocovariance
#'   lags; lag 0 means "contemporaneous only" and adds nothing beyond
#'   \code{orders}.
#' @param method   \code{"gmm"} (closed-form model moments) or \code{"smm"}
#'   (simulated moments).
#' @param weighting First-step weight matrix; see the Weighting section.
#' @param two_step Run the efficient second step (default \code{TRUE}).
#' @param n_sim    SMM simulation length multiplier \eqn{\tau}: the simulated
#'   path has \code{n_sim * T} periods (default 10).
#' @param seed     Seed for the SMM common random numbers (default 20260902).
#'   The global RNG stream is restored afterwards.
#' @param compiled Optional \code{dynhr_compiled}; compiled from \code{model}
#'   when absent.
#' @param order    Perturbation order used for the moment map (1 or 2).
#'   Defaults to 2 when \code{orders} includes 3 or 4, else 1.
#' @param me_variance Measurement-error variance added to the diagonal of the
#'   model-implied contemporaneous covariance (default 0).  Being i.i.d. it
#'   does not enter \eqn{\Gamma(h)} for \eqn{h \ge 1}.
#' @param lower,upper Optional named numeric bounds for the optimiser.
#' @param optimizer Mode-finder passed to the package optimiser stack
#'   (default \code{"newrat"}; see \code{\link{find_mode}}).
#' @param n_iter   Optimiser iteration budget (default 2000).
#' @param bandwidth Newey-West lag truncation; \code{NULL} uses
#'   \eqn{\lfloor 4 (T/100)^{2/9}\rfloor}.
#' @param ridge    Ridge fraction for the \eqn{\hat\Omega} inversion.
#' @param jacobian \code{"auto"} (analytic where it reaches, else finite
#'   differences -- and the fallback is REPORTED, both as a one-time warning
#'   and in \code{$jacobian_method}), \code{"analytic"} (error if unavailable)
#'   or \code{"fd"}.
#' @param burn_in  SMM burn-in periods discarded before the simulated sample.
#' @param verbose  Print optimiser progress.
#'
#' @return An object of class \code{"dynhr_mom"}: a list with
#'   \code{estimate}, \code{se}, \code{vcov}, \code{criterion},
#'   \code{j_stat}, \code{j_df}, \code{j_pvalue}, \code{omega_rank},
#'   \code{weight_matrix},
#'   \code{moment_fit} (a data.frame with the sample and model moment, their
#'   difference and a standardized difference), \code{jacobian},
#'   \code{jacobian_method}, and the settings used.
#'
#' @references
#'   Hansen, L. P. (1982). Large sample properties of generalized method of
#'     moments estimators. \emph{Econometrica}, 50(4), 1029-1054.
#'   Duffie, D., & Singleton, K. J. (1993). Simulated moments estimation of
#'     Markov models of asset prices. \emph{Econometrica}, 61(4), 929-952.
#'   Donald, S. G., & Newey, W. K. (2001). Choosing the number of instruments.
#'     \emph{Econometrica}, 69(5), 1161-1191.
#'
#' @seealso \code{\link{estimate_gmm_weight_matrix}},
#'   \code{\link{make_log_posterior_cumulant}}, \code{\link{compute_moments}}
#' @examples
#' \donttest{
#' solved <- solve_model(system.file("extdata/models/rbc2shock.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' set.seed(1)
#' sim <- simulate_model(solved$dr, n_periods = 400L, model = solved$model)
#' Y   <- attr(sim, "levels")[, c("a", "b"), drop = FALSE]
#'
#' fit <- method_of_moments(solved, data = Y, obs_vars = c("a", "b"),
#'                          start = c(rho_a = 0.9, rho_b = 0.85),
#'                          moments = list(orders = 2L, lags = 0:2),
#'                          two_step = FALSE, verbose = FALSE)
#' fit
#' }
#' @export
method_of_moments <- function(model, data, obs_vars = NULL, start = NULL,
                              moments = list(orders = 1:2, lags = 0:4),
                              method = c("gmm", "smm"),
                              weighting = c("identity", "optimal",
                                            "newey_west", "diagonal"),
                              two_step = TRUE,
                              n_sim = 10L, seed = 20260902L,
                              compiled = NULL, order = NULL,
                              me_variance = 0,
                              lower = NULL, upper = NULL,
                              optimizer = "newrat", n_iter = 2000L,
                              bandwidth = NULL, ridge = 1e-6,
                              jacobian = c("auto", "analytic", "fd"),
                              burn_in = 100L, verbose = TRUE) {

  method    <- match.arg(method)
  weighting <- match.arg(weighting)
  jacobian  <- match.arg(jacobian)
  spec      <- .mom_spec(moments)
  orders    <- spec$orders
  lags      <- spec$lags

  ## ---- Unpack model / compiled -------------------------------------------
  if (inherits(model, "dynhr_solved")) {
    compiled <- compiled %||% model$compiled
    model    <- model$model
  }
  if (!inherits(model, "dynhr_mod"))
    stop("method_of_moments(): `model` must be a dynhr_mod or dynhr_solved.",
         call. = FALSE)
  if (is.null(compiled))
    compiled <- compile_model(model, verbose = FALSE,
                              max_order = if (any(orders >= 3L)) 2L else 1L)

  ## ---- Data / obs_vars ----------------------------------------------------
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    stop("method_of_moments(): `data` must be a numeric matrix.", call. = FALSE)
  if (is.null(obs_vars) || !length(obs_vars))
    obs_vars <- model$obs_vars %||% model$varobs_names
  if (is.null(obs_vars) || !length(obs_vars))
    stop("method_of_moments(): `obs_vars` not supplied and the model ",
         "declares no `varobs`.", call. = FALSE)
  if (!all(obs_vars %in% colnames(data)))
    stop("method_of_moments(): all `obs_vars` must be column names of `data`: ",
         paste(setdiff(obs_vars, colnames(data)), collapse = ", "),
         call. = FALSE)
  T_obs <- nrow(data)

  ## ---- Starting values ----------------------------------------------------
  if (is.null(start)) {
    ps <- tryCatch(extract_prior_spec(model), error = function(e) NULL)
    if (is.null(ps) || !nrow(ps))
      stop("method_of_moments(): supply `start` (the model has no ",
           "estimated_params block to take starting values from).",
           call. = FALSE)
    start <- setNames(ps$mean, ps$name)
  }
  if (is.null(names(start)) || any(!nzchar(names(start))))
    stop("method_of_moments(): `start` must be a NAMED numeric vector.",
         call. = FALSE)
  theta0 <- as.numeric(start)
  names(theta0) <- names(start)
  k <- length(theta0)

  known <- c(names(model$param_values), model$varexo_names %||% character(0))
  bad <- setdiff(names(theta0), known)
  if (length(bad))
    stop("method_of_moments(): `start` names not connected to the model ",
         "(neither a parameter nor a declared shock): ",
         paste(bad, collapse = ", "), call. = FALSE)

  solve_order <- as.integer(order %||% (if (any(orders >= 3L)) 2L else 1L))
  if (!solve_order %in% c(1L, 2L))
    stop("method_of_moments(): `order` must be 1 or 2.", call. = FALSE)
  if (any(orders >= 3L) && solve_order < 2L)
    stop("method_of_moments(): cumulant orders 3-4 need order = 2.",
         call. = FALSE)

  ## ---- Sample moments -----------------------------------------------------
  m_emp <- .mom_sample_moment_vector(data, obs_vars, orders, lags)
  p <- length(m_emp)
  mom_names <- .mom_moment_names(obs_vars, orders, lags)
  if (p != length(mom_names))
    stop(".mom internal: moment vector length ", p, " != ", length(mom_names),
         " names.", call. = FALSE)
  if (p < k)
    stop("method_of_moments(): ", p, " moments cannot identify ", k,
         " parameters. Add lags or orders.", call. = FALSE)

  ## ---- SMM common random numbers -----------------------------------------
  n_sim <- as.integer(n_sim)
  U <- NULL
  n_sim_periods <- NA_integer_
  if (method == "smm") {
    if (n_sim < 1L) stop("method_of_moments(): `n_sim` must be >= 1.",
                         call. = FALSE)
    n_sim_periods <- n_sim * T_obs
    U <- .mom_crn_draws(n_sim_periods + as.integer(burn_in),
                        length(model$varexo_names), seed)
  }

  ## ---- Model moment map ---------------------------------------------------
  cache <- new.env(parent = emptyenv())
  cache$dr <- NULL
  cache$params <- NULL

  model_moments <- function(theta, store = FALSE) {
    params <- .apply_theta_to_params(model, setNames(as.numeric(theta),
                                                     names(theta0)))
    sol <- .mom_solve_at(model, compiled, params, solve_order)
    if (is.null(sol)) return(NULL)
    m <- if (method == "gmm") {
      suppressWarnings(.build_moment_vector(sol$dr, model, sol$params,
                                            obs_vars, orders,
                                            me_variance = me_variance,
                                            lags = lags))
    } else {
      .mom_simulated_moments(sol$dr, model, sol$params, obs_vars, orders,
                             lags, U, n_sim_periods, as.integer(burn_in))
    }
    if (is.null(m) || length(m) != p || !all(is.finite(m))) return(NULL)
    if (store) { cache$dr <- sol$dr; cache$params <- sol$params }
    m
  }

  criterion <- function(theta, W) {
    m <- model_moments(theta)
    if (is.null(m)) return(NA_real_)
    d <- m_emp - m
    .mom_criterion_value(d, W, T_obs)
  }

  ## ---- Weight matrices ----------------------------------------------------
  W_identity <- diag(p) / p

  build_W <- function(kind, m_center) {
    if (kind == "identity") {
      return(structure(W_identity, method = "identity"))
    }
    om <- .mom_omega(data, obs_vars, orders, lags, m_center = m_center,
                     bandwidth = bandwidth, ridge = ridge)
    if (kind == "diagonal") {
      W <- diag(1 / diag(om$Omega_reg), nrow = p)
      rk <- p
    } else {
      inv <- .mom_omega_inverse(om$Omega, ridge = ridge)
      W  <- inv$W
      rk <- inv$rank
      if (!is.finite(inv$condition_number) || inv$condition_number > 1e12)
        warning("method_of_moments(): the retained long-run covariance is ",
                "ill-conditioned (condition number ",
                if (is.finite(inv$condition_number))
                  format(inv$condition_number, scientific = TRUE) else "Inf",
                "). Increase `ridge` or drop moments.", call. = FALSE)
    }
    attr(W, "method")    <- kind
    attr(W, "bandwidth") <- om$bandwidth
    attr(W, "rank")      <- rk
    attr(W, "Omega")     <- om$Omega
    W
  }

  ## The efficient weight needs a centre; before any estimate exists the model
  ## moments at the starting value are the natural one.
  m_start <- model_moments(theta0)
  if (is.null(m_start))
    stop("method_of_moments(): the model could not be solved at `start` ",
         "(steady state, Blanchard-Kahn or stationarity failure).",
         call. = FALSE)

  W1 <- build_W(weighting, if (weighting == "identity") NULL else m_start)

  ## ---- Optimisation -------------------------------------------------------
  lo <- setNames(rep(-Inf, k), names(theta0))
  up <- setNames(rep(Inf,  k), names(theta0))
  if (!is.null(lower)) lo[names(lower)] <- as.numeric(lower)
  if (!is.null(upper)) up[names(upper)] <- as.numeric(upper)
  bounds <- data.frame(name = names(theta0), lower = as.numeric(lo),
                       upper = as.numeric(up), stringsAsFactors = FALSE)

  optimise_step <- function(theta_start, W) {
    target <- function(theta) {
      q <- criterion(theta, W)
      list(logpost = if (is.na(q)) -Inf else -0.5 * q)
    }
    res <- .run_mode_finding(target, theta_start, bounds,
                             nm_maxit = as.integer(n_iter),
                             method = optimizer, verbose = verbose)
    res
  }

  fit1 <- optimise_step(theta0, W1)
  theta_hat <- fit1$theta_mode
  W_used <- W1
  steps <- 1L

  ## ---- Efficient second step ---------------------------------------------
  two_step_done <- FALSE
  if (isTRUE(two_step)) {
    if (p >= T_obs) {
      warning("method_of_moments(): skipping the efficient second step -- ",
              "p = ", p, " moments with T = ", T_obs, " observations makes ",
              "the long-run covariance rank-deficient (the too-many-moments ",
              "problem). Reduce `lags`/`orders` or keep identity weighting.",
              call. = FALSE)
    } else {
      m1 <- model_moments(theta_hat)
      if (is.null(m1)) {
        warning("method_of_moments(): first-step estimate is infeasible; ",
                "skipping the second step.", call. = FALSE)
      } else {
        W2 <- build_W("newey_west", m1)
        fit2 <- optimise_step(theta_hat, W2)
        theta_hat <- fit2$theta_mode
        W_used <- W2
        fit1 <- fit2
        steps <- 2L
        two_step_done <- TRUE
      }
    }
  }

  ## ---- Final moments, Jacobian, inference --------------------------------
  m_hat <- model_moments(theta_hat, store = TRUE)
  if (is.null(m_hat))
    stop("method_of_moments(): the model is infeasible at the final ",
         "estimate.", call. = FALSE)
  g_hat <- m_emp - m_hat
  Q_hat <- .mom_criterion_value(g_hat, W_used, T_obs)

  jac <- .mom_jacobian(model, compiled, cache$dr, cache$params, theta_hat,
                       obs_vars, orders, lags, me_variance,
                       moment_fn = function(th) model_moments(th),
                       jacobian = if (method == "smm" && jacobian == "auto")
                         "fd" else jacobian)
  J <- jac$J

  om_fin <- .mom_omega(data, obs_vars, orders, lags, m_center = m_hat,
                       bandwidth = bandwidth, ridge = ridge)
  Omega <- om_fin$Omega
  inv_fin <- .mom_omega_inverse(Omega, ridge = ridge)
  W_opt   <- inv_fin$W
  om_rank <- inv_fin$rank

  ## Sandwich. g(theta) = m_emp - m(theta), so dg/dtheta' = -J and the signs
  ## cancel throughout; c = 1 + 1/tau is the SMM simulation-noise inflation.
  scale_c <- if (method == "smm") 1 + 1 / n_sim else 1
  V <- matrix(NA_real_, k, k, dimnames = list(names(theta0), names(theta0)))
  se <- setNames(rep(NA_real_, k), names(theta0))
  if (!is.null(J) && all(is.finite(J))) {
    A  <- t(J) %*% W_used %*% J
    Ai <- tryCatch(.mom_sym_inverse(A, "J' W J"), error = function(e) NULL)
    if (!is.null(Ai)) {
      B <- t(J) %*% W_used %*% Omega %*% W_used %*% J
      V <- scale_c * (Ai %*% B %*% Ai) / T_obs
      V <- (V + t(V)) * 0.5
      dimnames(V) <- list(names(theta0), names(theta0))
      dv <- diag(V)
      se <- setNames(ifelse(dv > 0, sqrt(dv), NA_real_), names(theta0))
    }
  }

  ## ---- J-test -------------------------------------------------------------
  ## df = (rank of Omega) - k, NOT p - k: with vec() rather than vech() the
  ## symmetric contemporaneous block carries exactly duplicated moments, which
  ## are not independent over-identifying restrictions.
  j_stat <- T_obs * as.numeric(crossprod(g_hat, W_opt %*% g_hat)) / scale_c
  j_df   <- om_rank - k
  j_p    <- if (j_df > 0L) stats::pchisq(j_stat, df = j_df, lower.tail = FALSE)
            else NA_real_
  j_valid <- two_step_done || weighting %in% c("optimal", "newey_west")

  ## ---- Moment fit table ---------------------------------------------------
  om_sd <- sqrt(pmax(diag(Omega), 0) / om_fin$T_eff)
  fit_tbl <- data.frame(
    moment = mom_names,
    sample = m_emp,
    model  = m_hat,
    diff   = g_hat,
    t_stat = ifelse(om_sd > 0, g_hat / om_sd, NA_real_),
    stringsAsFactors = FALSE
  )
  rownames(fit_tbl) <- NULL

  out <- list(
    estimate       = theta_hat,
    se             = se,
    vcov           = V,
    criterion      = Q_hat,
    loglik_equiv   = -0.5 * Q_hat,
    j_stat         = j_stat,
    j_df           = j_df,
    j_pvalue       = j_p,
    j_valid        = j_valid,
    weight_matrix  = W_used,
    Omega          = Omega,
    omega_rank     = om_rank,
    moment_fit     = fit_tbl,
    moment_names   = mom_names,
    jacobian       = J,
    jacobian_method = jac$method,
    jacobian_fell_back = jac$fell_back,
    jacobian_reason = jac$reason,
    method         = method,
    weighting      = weighting,
    two_step       = two_step_done,
    steps          = steps,
    moments        = list(orders = orders, lags = lags),
    obs_vars       = obs_vars,
    n_obs          = T_obs,
    n_moments      = p,
    n_params       = k,
    n_sim          = if (method == "smm") n_sim else NA_integer_,
    seed           = if (method == "smm") seed else NA_integer_,
    bandwidth      = om_fin$bandwidth,
    optimizer      = fit1$method,
    convergence    = fit1$convergence,
    iterations     = fit1$iterations,
    start          = theta0
  )
  class(out) <- "dynhr_mom"
  out
}


# ============================================================================
# S3 methods
# ============================================================================

#' Print a method-of-moments fit
#'
#' @param x A \code{"dynhr_mom"} object from \code{\link{method_of_moments}}.
#' @param digits Significant digits for the estimate table.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.dynhr_mom <- function(x, digits = 4L, ...) {
  cat("=== dynhr method_of_moments ===\n")
  cat(sprintf("  Estimator      : %s (%s weighting%s)\n",
              toupper(x$method), x$weighting,
              if (isTRUE(x$two_step)) ", two-step efficient" else ""))
  cat(sprintf("  Moments        : orders %s | lags %s  ->  p = %d (rank %d)\n",
              paste(x$moments$orders, collapse = ","),
              if (length(x$moments$lags))
                paste(x$moments$lags, collapse = ",") else "none",
              x$n_moments, x$omega_rank))
  cat(sprintf("  Observables    : %s   (T = %d)\n",
              paste(x$obs_vars, collapse = ", "), x$n_obs))
  if (identical(x$method, "smm"))
    cat(sprintf("  Simulation     : n_sim = %d (tau), seed = %s\n",
                x$n_sim, format(x$seed)))
  cat(sprintf("  Criterion      : %.6g\n", x$criterion))
  cat(sprintf("  Moment Jacobian: %s%s\n", x$jacobian_method,
              if (isTRUE(x$jacobian_fell_back))
                paste0("  [fell back: ", x$jacobian_reason, "]") else ""))
  cat("\n  Estimates (sandwich SE):\n")
  tab <- data.frame(estimate = as.numeric(x$estimate),
                    se = as.numeric(x$se),
                    row.names = names(x$estimate))
  tab$t <- tab$estimate / tab$se
  print(round(tab, digits))
  cat(sprintf("\n  J-statistic    : %.4f  (df = %d, p = %.4f)%s\n",
              x$j_stat, x$j_df, x$j_pvalue,
              if (isTRUE(x$j_valid)) "" else
                "  [identity weighting: reference distribution approximate]"))
  invisible(x)
}

#' Summarise a method-of-moments fit
#'
#' Adds the moment-fit table (sample vs model moment, difference, and the
#' difference standardized by the HAC standard error of the moment) to what
#' \code{\link{print.dynhr_mom}} shows.
#'
#' @param object A \code{"dynhr_mom"} object.
#' @param n_moments Maximum number of moment rows to print (default 20).
#' @param digits Significant digits.
#' @param ... Ignored.
#' @return \code{object}, invisibly.
#' @export
summary.dynhr_mom <- function(object, n_moments = 20L, digits = 4L, ...) {
  print(object, digits = digits)
  cat("\n  Moment fit:\n")
  tb <- object$moment_fit
  n <- min(nrow(tb), as.integer(n_moments))
  show <- tb[seq_len(n), , drop = FALSE]
  show[, 2:5] <- lapply(show[, 2:5], function(z) signif(z, digits))
  print(show, row.names = FALSE)
  if (n < nrow(tb))
    cat(sprintf("  ... %d further moment(s) not shown.\n", nrow(tb) - n))
  invisible(object)
}
