## R/hank-kalman-ar.R
## --------------------------------------------------------------------------
## Exact-AR(1) variant of the truncated-MA stacked-covariance HANK likelihood
## (the ABRS autocovariance form of R/hank-estimation.R and the
## hank_state_space() bridge of R/hank-kalman.R).
##
## WHY.  The truncated-MA likelihoods (hank_loglik_aggregate, hank_loglik_ss,
## hank_state_space -> hank_kalman_loglik) keep the first q MA coefficients of
## each shock's impulse response and drop the rest.  For a shock with AR(1)
## persistence rho the dropped tail carries variance of order
## rho^(2q)/(1 - rho^2): negligible at tame persistences, catastrophic near
## the unit root.  Documented failure case (downstream application,
## 2026-07-12): at a posterior mode with
## rho = 0.9992 the plain truncated-MA likelihood (q = 200) was -337
## log-points off the exact Kalman filter of the equivalent .mod state space,
## at parameter values where a tame-rho oracle had previously matched to
## 3.4e-4.
##
## FIX (quasi-differencing).  The sequence-space Theta^z_s is the response of
## the observables to a unit innovation in shock z, i.e. to the anticipated-
## from-0 AR(1) driving path rho_z^t.  Quasi-differencing the MA coefficients,
##
##     Psi_s = Theta_s - rho_z * Theta_{s-1}        (Theta_{-1} = 0),
##
## cancels the shock's own geometric factor EXACTLY -- including the
## anticipation terms of the sequence-space solution: in the time-invariant
## region Psi_s equals the response to an UNANTICIPATED unit LEVEL pulse of
## the shock, which decays at the model's internal rate rather than at rho_z.
## The observable contribution of shock z is then the ARMA(1, q-1) process
##
##     y^z_t = rho_z y^z_{t-1} + sum_{s < q} Psi_s eps^z_{t-s},
##
## whose autocovariance has the closed form
##
##     G_z(k) = sig_z^2 * sum_d m_z(d) rho_z^{|k-d|} / (1 - rho_z^2),
##     m_z(d) = sum_j Psi_j Psi_{j-d}'          (finite cross-lag sums),
##
## i.e. EXACT in the shock persistence; the only truncation left is the
## internal-mode tail of Psi (roots of the model itself, not of the shock
## process).  At rho_z = 0 the formula reduces bit-exact to the plain
## truncated-MA autocovariance (0^0 = 1 keeps the d = k term).
##
## VALIDATION HISTORY (against an independent reference implementation
## of the same likelihood):
##   - reduces bit-exact to the plain truncated-MA formula at rho = 0;
##   - matches a brute-force MA(5000) analytic oracle to 6e-4 log-points at a
##     near-unit-root mode where the plain formula is off by +9.0;
##   - matches the exact .mod Kalman filter on rank_baseline_v3 at
##     0.97-capped persistences to -0.092 log-points (q = 300, T_h = 400,
##     10 observables), where the plain formula loses hundreds.
## The oracle tests in tests/testthat/test-hank-kalman-ar.R mirror all three.
##
## RESIDUAL REPRESENTABILITY BOUND (terminal-boundary contamination).  The
## exact-AR treatment is only as good as Theta itself.  The sequence-space
## solve imposes a terminal steady-state boundary at T_h; when a shock is so
## persistent that its level is still material at the horizon (|rho|^T_h not
## small), the retained Theta rows are contaminated WELL INSIDE the window
## and no q fixes it.  Measured on rank_baseline_v3 (T_h = 400, exact .mod
## MA coefficients as oracle -- the diag5_theta_boundary.R pattern):
## rho = 0.9992 gives 7e-2..1e-1 relative Theta error across rows 25..390
## (shock level still 0.73 at the boundary); rho = 0.985 gives ~3e-4;
## rho <= 0.96 is at machine precision.  Likelihood error scales as
## rho^(2q) + rho^T_h: at T_h = 400 a 0.985 cap still cost -11.3 log-points
## at the oracle, a 0.98 cap ~ -1.  Practical rule: keep |rho|^T_h below
## ~1e-3 (cap estimated persistences, e.g. at 0.98 for T_h = 400, or
## re-solve with a larger T_h).  hank_theta_boundary_check() computes the
## per-shock proxies; the likelihoods warn once per call when it is
## violated (check_boundary = FALSE to silence, e.g. inside samplers whose
## priors/bounds already enforce a cap).
## --------------------------------------------------------------------------


## Unit-innovation-variance exact-AR(1) autocovariance kernel of ONE shock:
## G(k), k = 0..n_lags, from the first q rows of Theta and persistence rho.
## Vectorized over lags: with M the q x n_obs^2 matrix of vec'd cross-lag
## sums m(d), the +/-d folds are two small matrix products,
##   G(k) = sum_{d>=0} [ m(d) w(k,d) + m(d)' w(k,-d) ] - m(0) w(k,0),
##   w(k,d) = rho^{|k-d|} / (1 - rho^2).
## No validation here; callers validate.
.hank_ar_autocov_kernel <- function(Theta, rho, n_lags, q) {
  n_obs <- ncol(Theta)
  Th  <- Theta[seq_len(q), , drop = FALSE]
  Psi <- Th - rho * rbind(rep(0, n_obs), Th[-q, , drop = FALSE])  # quasi-diff

  ## cross-lag sums m(d), d = 0..q-1 (m(-d) = t(m(d)))
  m <- lapply(0:(q - 1L), function(d) {
    idx <- seq_len(q - d)
    crossprod(Psi[idx + d, , drop = FALSE], Psi[idx, , drop = FALSE])
  })

  ks <- 0:n_lags; ds <- 0:(q - 1L)
  Wp <- outer(ds, ks, function(d, k) rho^abs(k - d)) / (1 - rho^2)
  Wm <- outer(ds, ks, function(d, k) rho^(k + d))    / (1 - rho^2)
  M  <- matrix(unlist(m), nrow = length(ds), byrow = TRUE)  # row d = vec(m(d))
  Gp <- crossprod(Wp, M)              # (n_lags+1) x n^2: sum_d w(k, d) vec(m(d))
  Gm <- crossprod(Wm, M)              # (n_lags+1) x n^2: sum_d w(k,-d) vec(m(d))

  G <- array(0, c(n_obs, n_obs, n_lags + 1L))
  m0 <- matrix(m[[1L]], n_obs, n_obs)
  for (k in seq_len(n_lags + 1L))
    G[, , k] <- matrix(Gp[k, ], n_obs, n_obs) +
      t(matrix(Gm[k, ], n_obs, n_obs)) - m0 * Wp[1L, k]
  G
}


## Gaussian loglik of the stacked demeaned sample under autocovariance array G
## (must already include any measurement error in G[,,1]), with exact
## missing-data handling by row deletion of the stacked system.
##
## Assembly is a single indexed gather from as.numeric(G): element (r, c) of
## the kept stacked covariance is G[i_r, j_c, |t_r - t_c| + 1], transposed
## when the column block is later -- exactly .hank_stacked_cov's block-fill
## semantics (G(t-s) for t >= s, t(G(s-t)) otherwise, then symmetrized), so
## the result is BIT-IDENTICAL to the old O(T_data^2) R block-fill loop
## (symmetrize-then-subset == subset-then-symmetrize elementwise). The gather
## index depends only on (T_data, n_obs, missing-data pattern); passing a
## `cache` environment reuses it across calls -- ported from the paper-side
## hank_abrs_loglik_ar(cache=), where this replaced the block-fill as the
## warm-eval bottleneck (0.33 s -> 0.017 s per likelihood eval).
## The gather index alone: it depends ONLY on (T_data, n_obs, missing pattern),
## so it survives every parameter move including a STRUCTURAL one (which
## invalidates every autocovariance slab but not this). Split out of
## .hank_stacked_loglik so hank_loglik_ar_grad() can share the one index --
## and the one factorization built from it -- instead of building a second.
.hank_stacked_index <- function(Y, cache = NULL) {
  Td <- nrow(Y); n_obs <- ncol(Y)
  yv <- as.numeric(t(Y))                       # time-major stacking
  keep <- is.finite(yv)
  if (!is.null(cache) && !is.null(cache$Sidx) &&
      identical(cache$keep, keep) && identical(cache$dims, c(Td, n_obs)))
    return(list(idx = cache$Sidx, keep = keep, yv = yv))
  N  <- Td * n_obs
  bT <- rep(seq_len(Td), each = n_obs)         # time block of each row/col
  bI <- rep.int(seq_len(n_obs), Td)            # within-block observable
  nn <- n_obs * n_obs
  Gidx <- matrix(0L, N, N)
  for (cc in seq_len(N)) {
    cB <- bT[cc]; cj <- bI[cc]
    k    <- abs(bT - cB)
    swap <- bT < cB                            # column block later -> transpose
    row_e <- ifelse(swap, cj, bI)
    col_e <- ifelse(swap, bI, cj)
    Gidx[, cc] <- k * nn + (col_e - 1L) * n_obs + row_e
  }
  idx <- Gidx[keep, keep, drop = FALSE]
  if (!is.null(cache)) {
    cache$Sidx <- idx; cache$keep <- keep; cache$dims <- c(Td, n_obs)
  }
  list(idx = idx, keep = keep, yv = yv)
}

## The stacked Gaussian loglik from an already-built Cholesky factor. Kept as
## one expression so every caller (hank_loglik_ar and hank_loglik_ar_grad)
## returns BIT-identical floats from the same (ch, yk).
.hank_stacked_loglik_from_chol <- function(ch, yk) {
  z <- backsolve(ch, yk, transpose = TRUE)
  -0.5 * (length(yk) * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(z^2))
}

.hank_stacked_loglik <- function(Y, G, cache = NULL) {
  Td <- nrow(Y); n_obs <- ncol(Y)
  if (dim(G)[1L] != n_obs)
    stop("stacked loglik: dim(G)[1] != ncol(Y)")
  if (dim(G)[3L] < Td)
    stop("autocovariance max lag < T_data - 1; increase T_h/n_lags")
  ix <- .hank_stacked_index(Y, cache)
  yk <- ix$yv[ix$keep]
  S  <- matrix(as.numeric(G)[ix$idx], length(yk), length(yk))
  S  <- (S + t(S)) / 2                         # symmetrize against round-off
  .hank_stacked_loglik_from_chol(chol(S), yk)
}


## Per-shock unscaled autocovariance slab, optionally cached. The slab
## sigma^-2 * Gamma_z depends only on (Theta_z, rho_z, n_lags, q_z) -- NOT on
## sigma_z or the measurement error -- so with a cache environment a
## one-coordinate optimizer/gradient perturbation recomputes at most ONE
## shock's slab (a sigma or me move recomputes none). Exactness-preserving:
## a cache hit returns the identical floats the kernel would produce.
.hank_ar_autocov_slab <- function(z, Theta, rho_z, n_lags, q_z, cache) {
  if (is.null(cache))
    return(.hank_ar_autocov_kernel(Theta, rho_z, n_lags, q_z))
  ent <- cache$A[[z]]
  if (is.null(ent) || ent$rho != rho_z || ent$q != q_z ||
      ent$n_lags != n_lags || !identical(ent$Theta, Theta)) {
    A <- .hank_ar_autocov_kernel(Theta, rho_z, n_lags, q_z)
    if (is.null(cache$A)) cache$A <- list()
    cache$A[[z]] <- list(rho = rho_z, q = q_z, n_lags = n_lags,
                         Theta = Theta, A = A)
  }
  cache$A[[z]]$A
}


## Warn (once per call) when |rho|^T_h says the sequence-space Theta is
## terminal-boundary contaminated; see the header and
## hank_theta_boundary_check().
.hank_theta_boundary_warn <- function(rho, T_h, tol, caller,
                                      silence = "check_boundary = FALSE") {
  lev <- abs(rho)^T_h
  bad <- which(lev > tol)
  if (!length(bad)) return(invisible(NULL))
  nm <- names(lev)[bad]
  if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("shock ", bad)
  warning(caller, ": |rho|^T_h = ",
          paste(sprintf("%.2e (%s)", lev[bad], nm), collapse = ", "),
          " exceeds ", format(tol),
          " -- the sequence-space Theta is terminal-boundary contaminated at",
          " this persistence (the shock level is still material at the solve",
          " horizon T_h = ", T_h, "), and the exact-AR quasi-differencing",
          " cannot repair Theta itself. Cap |rho| (<= ",
          sprintf("%.4f", tol^(1 / T_h)),
          " here) or re-solve with a larger T_h; see",
          " hank_theta_boundary_check(). Set ", silence, " to",
          " silence.", call. = FALSE)
  invisible(NULL)
}


## REPRESENTABILITY GATE for a posterior CLOSURE (as opposed to a single
## likelihood call). `check_boundary` on the likelihoods warns per call, which
## is unusable inside a sampler: it either floods the log or, as the shipped
## factories did, gets turned off and the bound stops existing. A closure needs
## a per-DRAW decision instead, so this returns "should this draw be rejected?"
## and warns at most once per closure:
##   "ignore" -- no bound (the old, silent behaviour)
##   "warn"   -- warn ONCE, naming the offending shocks and the implied cap,
##               then evaluate anyway (default: does not change any posterior)
##   "reject" -- treat |rho|^T_h <= tol as part of the PRIOR SUPPORT and return
##               -Inf, exactly as the existing |rho| >= 1 guard does
## Note what the bound is NOT: it is not a statement about stationarity. Theta
## itself is contaminated by the sequence-space solve's terminal boundary out
## here (see the file header), so the likelihood is wrong rather than merely
## imprecise, and no q or cache repairs it.
.hank_ar_boundary_gate <- function(rho, T_h, boundary, tol, caller, state) {
  if (identical(boundary, "ignore")) return(FALSE)
  if (!length(rho) || !is.finite(T_h)) return(FALSE)
  if (!any(abs(rho)^T_h > tol)) return(FALSE)
  if (identical(boundary, "reject")) return(TRUE)
  if (!isTRUE(state$warned)) {
    state$warned <- TRUE
    .hank_theta_boundary_warn(rho, T_h, tol, caller,
                              silence = "boundary = \"ignore\"")
  }
  FALSE
}


#' Exact-AR(1) autocovariance function from MA coefficients
#'
#' Exact-in-persistence counterpart of \code{\link{hank_autocov}} for a scalar
#' AR(1) shock: \code{Theta} are the unit-innovation MA coefficients of the
#' observables (the impulse responses to the anticipated-from-0 AR(1) driving
#' path \eqn{\rho^t}, e.g. from \code{\link{hank_ma_coefficients}}).  The
#' shock's geometric factor is removed exactly by quasi-differencing
#' (\eqn{\Psi_s = \Theta_s - \rho\,\Theta_{s-1}}) and restored in closed form,
#' \deqn{\Gamma(k) = \sigma^2 \sum_d m(d)\, \rho^{|k-d|} / (1-\rho^2), \qquad
#'       m(d) = \sum_j \Psi_j \Psi_{j-d}',}
#' so only the internal-mode tail of \eqn{\Psi} is truncated at \code{q}
#' terms -- not the shock persistence tail \eqn{\rho^{2q}} that the plain
#' truncated-MA \code{\link{hank_autocov}} loses.  Reduces bit-exact to
#' \code{\link{hank_autocov}} at \code{rho = 0}.
#'
#' @param Theta \code{T_h x n_obs} MA coefficients (see
#'   \code{\link{hank_ma_coefficients}}).
#' @param rho AR(1) persistence of the shock, \code{abs(rho) < 1}.
#' @param sigma_eps Standard deviation of the shock innovation.
#' @param n_lags Integer: maximum lag \eqn{k} to return.
#' @param q Integer: number of quasi-differenced MA terms retained (default all
#'   rows of \code{Theta}).
#'
#' @return An \code{n_obs x n_obs x (n_lags+1)} array with slice \code{k+1}
#'   giving \eqn{\Gamma(k)}.
#' @seealso \code{\link{hank_autocov}}, \code{\link{hank_loglik_aggregate_ar}},
#'   \code{\link{hank_loglik_ar}}
#' @export
hank_autocov_ar <- function(Theta, rho, sigma_eps, n_lags, q = NULL) {
  Theta <- as.matrix(Theta)
  if (abs(rho) >= 1)
    stop("hank_autocov_ar: `rho` must satisfy abs(rho) < 1.")
  T_h <- nrow(Theta)
  if (is.null(q)) q <- T_h
  q <- min(q, T_h)
  sigma_eps^2 * .hank_ar_autocov_kernel(Theta, rho, n_lags, q)
}


#' Exact-AR(1) Gaussian log-likelihood of aggregate HANK data (scalar shock)
#'
#' Drop-in exact-in-persistence variant of \code{\link{hank_loglik_aggregate}}:
#' the same stacked block-Toeplitz Gaussian likelihood, but with the
#' autocovariance built by \code{\link{hank_autocov_ar}} so the shock's AR(1)
#' variance tail is kept in closed form instead of being truncated with the MA
#' lag polynomial.  Use this whenever \code{rho} may approach the unit root:
#' the plain truncated-MA likelihood loses variance of order
#' \eqn{\rho^{2q}/(1-\rho^2)} (documented failure: \eqn{-337} log-points at
#' \eqn{\rho = 0.9992}, \eqn{q = 200}, vs. the exact Kalman filter), while
#' here the only truncation left is the internal-mode tail of the
#' quasi-differenced coefficients.  At \code{rho = 0} the two likelihoods
#' agree bit-exactly.
#'
#' Unlike \code{\link{hank_loglik_aggregate}}, missing observations
#' (\code{NA}/\code{NaN} entries of \code{Y}) are handled exactly, by row
#' deletion of the stacked system.
#'
#' @section Representability bound: the exact-AR treatment repairs the
#' likelihood formula, not \code{Theta} itself.  A sequence-space \code{Theta}
#' solved on horizon \code{T_h} with a terminal steady-state boundary is
#' contaminated well inside the window once \eqn{|\rho|^{T_h}} is no longer
#' small (measured reference: \eqn{\rho = 0.9992} at \code{T_h = 400} leaves
#' 7e-2..1e-1 relative error across rows 25..390; likelihood error scales as
#' \eqn{\rho^{2q} + \rho^{T_h}}).  Keep \eqn{|\rho|^{T_h}} below ~1e-3 (e.g.
#' cap \eqn{\rho} at 0.98 for \code{T_h = 400}); by default the function warns
#' when the bound is violated.  See \code{\link{hank_theta_boundary_check}}.
#'
#' @param Y \code{T_data x n_obs} matrix of demeaned aggregate observations
#'   (deviations from steady state), columns in the same order as
#'   \code{Theta}; \code{NA}/\code{NaN} entries allowed.
#' @param Theta MA coefficients (see \code{\link{hank_ma_coefficients}}).
#' @param rho AR(1) persistence of the shock, \code{abs(rho) < 1}.
#' @param sigma_eps Shock innovation standard deviation.
#' @param me_var Measurement-error variance added to every observable
#'   (default 0).
#' @param q Number of quasi-differenced MA terms retained (default all rows of
#'   \code{Theta}).
#' @param check_boundary Logical: warn when \eqn{|\rho|^{T_h}} exceeds
#'   \code{boundary_tol} (see the representability-bound section).  Set
#'   \code{FALSE} inside samplers whose priors/bounds already cap \code{rho}.
#' @param boundary_tol Threshold for the boundary warning (default
#'   \code{1e-3}).
#' @param cache Optional environment (\code{new.env()}) reused across calls
#'   with the same \code{Y} shape and missing-data pattern.  Two
#'   exactness-preserving layers (a cache hit returns identical floats to
#'   the uncached path): per-shock unscaled autocovariance slabs keyed on
#'   \code{(rho, q, Theta)} -- so a one-coordinate optimizer/gradient
#'   perturbation recomputes at most one shock's slab, and a
#'   \code{sigma}/measurement-error move recomputes none -- and the stacked
#'   covariance gather index keyed on \code{(T_data, n_obs, NA pattern)}.
#'   Intended for iterative callers (samplers, optimizers, per-coordinate
#'   finite differences); \code{\link{make_log_posterior_hank}} wires one in
#'   automatically for \code{likelihood = "exact_ar"}.  Default \code{NULL}
#'   (no caching).
#'
#' @return The scalar Gaussian log-likelihood.
#' @seealso \code{\link{hank_loglik_aggregate}}, \code{\link{hank_loglik_ar}},
#'   \code{\link{hank_theta_boundary_check}}
#' @export
hank_loglik_aggregate_ar <- function(Y, Theta, rho, sigma_eps, me_var = 0,
                                     q = NULL, check_boundary = TRUE,
                                     boundary_tol = 1e-3, cache = NULL) {
  Y <- as.matrix(Y)
  Theta <- as.matrix(Theta)
  n_obs <- ncol(Y)
  if (ncol(Theta) != n_obs)
    stop("hank_loglik_aggregate_ar: ncol(Theta) must equal ncol(Y).")
  if (abs(rho) >= 1)
    stop("hank_loglik_aggregate_ar: `rho` must satisfy abs(rho) < 1.")
  if (!is.null(cache) && !is.environment(cache))
    stop("hank_loglik_aggregate_ar: `cache` must be an environment ",
         "(e.g. new.env()) or NULL.")
  T_h <- nrow(Theta)
  if (is.null(q)) q <- T_h
  q <- min(q, T_h)
  if (isTRUE(check_boundary))
    .hank_theta_boundary_warn(rho, T_h, boundary_tol,
                              "hank_loglik_aggregate_ar")

  G <- sigma_eps^2 *
    .hank_ar_autocov_slab("shock", Theta, rho, nrow(Y) - 1L, q, cache)
  G[, , 1L] <- G[, , 1L] + diag(me_var, n_obs)
  .hank_stacked_loglik(Y, G, cache = cache)
}


#' Exact-AR(1) stacked-covariance log-likelihood of a multi-shock HANK model
#'
#' Exact-in-persistence alternative to running the Kalman filter on the
#' truncated-MA companion form (\code{\link{hank_state_space}} \eqn{\to}
#' \code{\link{hank_kalman_loglik}}): the Gaussian likelihood of the stacked
#' demeaned sample under the ABRS block-Toeplitz autocovariance, with each
#' shock's AR(1) persistence handled in closed form by quasi-differencing its
#' MA coefficients (\eqn{\Psi_s = \Theta_s - \rho_z\,\Theta_{s-1}}, which
#' cancels the shock's geometric factor exactly, anticipation terms included):
#' \deqn{\Gamma_z(k) = \sigma_z^2 \sum_d m_z(d)\, \rho_z^{|k-d|} /
#'       (1-\rho_z^2), \qquad m_z(d) = \sum_j \Psi_j \Psi_{j-d}'.}
#' The truncated-MA paths lose each shock's persistence variance tail
#' \eqn{\rho_z^{2q}/(1-\rho_z^2)} -- a documented \eqn{-337} log-point error
#' at a near-unit-root posterior mode (\eqn{\rho = 0.9992}, \code{q = 200})
#' vs. the exact Kalman filter -- while here the only truncation left is the
#' internal-mode tail of \eqn{\Psi_z}.  At all-zero persistences the two
#' formulations coincide (bit-exact against the plain ABRS autocovariance).
#'
#' Cost: one \code{(n_kept x n_kept)} Cholesky of the stacked covariance,
#' where \code{n_kept <= T_data * n_obs}; missing observations
#' (\code{NA}/\code{NaN}) are handled exactly by row deletion of the stacked
#' system.  Per-observable measurement error is supported (and required
#' whenever the observation set is stochastically singular, e.g. more
#' observables than shocks or observables tied by a static identity).
#'
#' @inheritSection hank_loglik_aggregate_ar Representability bound
#'
#' @param Y \code{T_data x n_obs} matrix (or data frame) of demeaned
#'   observations, columns in the order of \code{ss$obs_names} (for a
#'   \code{\link{hank_state_space}} object) or of the \code{Theta} columns;
#'   \code{NA}/\code{NaN} entries allowed.
#' @param ss Either a \code{\link{hank_state_space}} object (its stored
#'   \code{Theta_list}, per-shock \code{rho}/\code{sigma} and \code{q} are
#'   used; pass \code{rho}/\code{sigma}/\code{q} to override), or a named list
#'   of \code{T_h x n_obs} unit-innovation MA-coefficient matrices, one per
#'   shock (the response to the anticipated-from-0 AR(1) driving path
#'   \eqn{\rho_z^t} -- the same \code{Theta_list} that
#'   \code{\link{hank_state_space}} builds), in which case \code{rho} and
#'   \code{sigma} are required.
#' @param rho Named numeric vector of per-shock AR(1) persistences
#'   (\code{abs(rho) < 1}).  Defaults to the values stored in \code{ss}.
#' @param sigma Named numeric vector of per-shock innovation standard
#'   deviations.  Defaults to the values stored in \code{ss}.
#' @param me_sd Measurement-error standard deviation: scalar (recycled) or
#'   length-\code{n_obs} vector.  Default 0 (allowed only if the observation
#'   set is non-singular).
#' @param q Number of quasi-differenced MA terms retained per shock (default:
#'   \code{ss$q} for a state-space object, else all rows of \code{Theta}).
#' @inheritParams hank_loglik_aggregate_ar
#'
#' @return The scalar Gaussian log-likelihood.
#' @seealso \code{\link{hank_state_space}}, \code{\link{hank_kalman_loglik}},
#'   \code{\link{hank_loglik_aggregate_ar}},
#'   \code{\link{hank_theta_boundary_check}},
#'   \code{\link{make_log_posterior_hank}} (option
#'   \code{likelihood = "exact_ar"})
#' @export
hank_loglik_ar <- function(Y, ss, rho = NULL, sigma = NULL, me_sd = 0,
                           q = NULL, check_boundary = TRUE,
                           boundary_tol = 1e-3, cache = NULL) {
  Y <- as.matrix(Y)
  p <- .hank_ar_prepare(Y, ss, rho = rho, sigma = sigma, me_sd = me_sd, q = q,
                        check_boundary = check_boundary,
                        boundary_tol = boundary_tol, cache = cache)
  .hank_stacked_loglik(Y, p$G, cache = cache)
}


## Validation, (rho, sigma, q, Theta_list) normalization and the summed
## autocovariance array G -- everything hank_loglik_ar() does before the
## stacked Gaussian evaluation. Factored out so hank_loglik_ar_grad() shares
## it verbatim instead of calling hank_loglik_ar() and then re-deriving the
## same normalization -- the duplication that also cost it a second gather and
## Cholesky per call. Every message keeps the "hank_loglik_ar:" prefix
## deliberately: these describe THAT function's documented contract, and the
## grad path reached them through it until this refactor.
.hank_ar_prepare <- function(Y, ss, rho = NULL, sigma = NULL, me_sd = 0,
                             q = NULL, check_boundary = TRUE,
                             boundary_tol = 1e-3, cache = NULL) {
  Td <- nrow(Y); n_obs <- ncol(Y)
  if (!is.null(cache) && !is.environment(cache))
    stop("hank_loglik_ar: `cache` must be an environment ",
         "(e.g. new.env()) or NULL.")

  if (inherits(ss, "dsge_ss")) {
    Theta_list <- ss$Theta_list
    if (is.null(Theta_list))
      stop("hank_loglik_ar: `ss` carries no Theta_list; ",
           "build it with hank_state_space().")
    if (is.null(rho)) rho <- ss$rho_vec
    if (is.null(rho))
      stop("hank_loglik_ar: `ss` stores no per-shock rho ",
           "(object predates rho_vec storage?); pass `rho` explicitly, ",
           "named by shock.")
    if (is.null(sigma))
      sigma <- stats::setNames(sqrt(diag(ss$Sigma_e)), ss$shock_names)
    if (is.null(q)) q <- ss$q
  } else if (is.list(ss)) {
    Theta_list <- ss
    if (is.null(rho) || is.null(sigma))
      stop("hank_loglik_ar: with a raw Theta list, `rho` and `sigma` are ",
           "required (named numeric vectors, one entry per shock).")
  } else {
    stop("hank_loglik_ar: `ss` must be a hank_state_space()/dsge_ss object ",
         "or a named list of MA-coefficient matrices.")
  }

  shocks <- names(Theta_list)
  if (is.null(shocks) || any(!nzchar(shocks)))
    stop("hank_loglik_ar: the Theta list must be named by shock.")
  if (!all(shocks %in% names(rho)) || !all(shocks %in% names(sigma)))
    stop("hank_loglik_ar: `rho` and `sigma` must carry an entry for every ",
         "shock in {", paste(shocks, collapse = ", "), "}.")
  if (any(abs(rho[shocks]) >= 1))
    stop("hank_loglik_ar: every `rho` must satisfy abs(rho) < 1.")
  me_sd <- rep_len(me_sd, n_obs)

  T_h <- nrow(as.matrix(Theta_list[[1L]]))
  if (is.null(q)) q <- T_h
  if (isTRUE(check_boundary))
    .hank_theta_boundary_warn(rho[shocks], T_h, boundary_tol,
                              "hank_loglik_ar")

  G <- array(0, c(n_obs, n_obs, Td))
  for (z in shocks) {
    Theta <- as.matrix(Theta_list[[z]])
    if (ncol(Theta) != n_obs)
      stop("hank_loglik_ar: Theta for shock '", z, "' has ", ncol(Theta),
           " columns; Y has ", n_obs, ".")
    q_z <- min(q, nrow(Theta))
    G <- G + sigma[[z]]^2 *
      .hank_ar_autocov_slab(z, Theta, rho[[z]], Td - 1L, q_z, cache)
  }
  G[, , 1L] <- G[, , 1L] + diag(me_sd^2, n_obs)
  list(G = G, Theta_list = Theta_list, rho = rho, sigma = sigma, q = q,
       shocks = shocks, me_sd = me_sd, T_h = T_h)
}


#' Terminal-boundary representability check for sequence-space MA coefficients
#'
#' The sequence-space solve behind \code{\link{hank_state_space}} /
#' \code{\link{hank_ma_coefficients}} imposes a terminal steady-state boundary
#' at the horizon \code{T_h}.  When a shock's AR(1) persistence is high enough
#' that its level is still material at the horizon, the resulting
#' MA coefficients \eqn{\Theta} are contaminated \emph{well inside} the
#' retained window -- a defect no likelihood formula (including the exact-AR
#' variants \code{\link{hank_loglik_ar}} /
#' \code{\link{hank_loglik_aggregate_ar}}) can repair, because \eqn{\Theta}
#' itself is wrong.  Reference measurement (exact state-space MA coefficients
#' as oracle, \code{T_h = 400}): \eqn{\rho = 0.9992} leaves 7e-2..1e-1
#' relative \eqn{\Theta} error across rows 25..390 (shock level 0.73 at the
#' boundary); \eqn{\rho = 0.985} ~3e-4; \eqn{\rho \le 0.96} machine
#' precision.  Likelihood error scales as \eqn{\rho^{2q} + \rho^{T_h}}; at
#' \code{T_h = 400} a 0.985 persistence cap still cost \eqn{-11.3} log-points
#' against an exact-filter oracle, a 0.98 cap about \eqn{-1}.
#'
#' Two per-shock proxies are reported:
#' \itemize{
#'   \item \code{boundary_level} \eqn{= |\rho|^{T_h}}, the shock's driving
#'     level at the terminal boundary -- the contamination scale.  Flagged
#'     (\code{flag_boundary}) when it exceeds \code{tol_boundary}; the
#'     remedy is a persistence cap (\eqn{|\rho| \le}
#'     \code{tol_boundary^(1/T_h)}) or a larger solve horizon.  The proxy is
#'     conservative: models whose observables load only statically on the
#'     shock have no boundary error at any \eqn{\rho}.
#'   \item \code{psi_tail}, the max absolute quasi-differenced coefficient
#'     \eqn{\Psi} over the last tenth of the retained window relative to the
#'     overall max -- the \emph{internal-mode} truncation left in the
#'     exact-AR formula (the analog of \eqn{\rho^{2q}} for the plain
#'     truncated-MA one).  Flagged (\code{flag_tail}) above \code{tol_tail};
#'     the remedy is a larger \code{q}/\code{T_h}.  \code{NA} when only
#'     \code{rho}/\code{T_h} are supplied.
#' }
#'
#' @param ss A \code{\link{hank_state_space}} object, a named list of
#'   \code{T_h x n_obs} MA-coefficient matrices (as in
#'   \code{\link{hank_loglik_ar}}; requires \code{rho}), or \code{NULL} to
#'   check persistences alone via \code{rho} and \code{T_h}.  A bare numeric
#'   vector is accepted and treated as \code{rho}.
#' @param rho Named numeric vector of per-shock AR(1) persistences (required
#'   unless \code{ss} is a \code{hank_state_space} object that stores them).
#' @param T_h Solve horizon (required when no \eqn{\Theta} matrices are
#'   supplied; otherwise defaults to \code{nrow(Theta)}).
#' @param q Retained quasi-differenced terms for the \code{psi_tail} proxy
#'   (default: \code{ss$q} for a state-space object, else \code{T_h}).
#' @param tol_boundary Flag threshold for \code{boundary_level} (default
#'   \code{1e-3}; see above for its calibration).
#' @param tol_tail Flag threshold for \code{psi_tail} (default \code{1e-2}).
#'
#' @return A data frame with one row per shock: \code{shock}, \code{rho},
#'   \code{T_h}, \code{q}, \code{boundary_level}, \code{psi_tail},
#'   \code{flag_boundary}, \code{flag_tail}.
#' @seealso \code{\link{hank_loglik_ar}},
#'   \code{\link{hank_loglik_aggregate_ar}}, \code{\link{hank_state_space}}
#' @export
hank_theta_boundary_check <- function(ss = NULL, rho = NULL, T_h = NULL,
                                      q = NULL, tol_boundary = 1e-3,
                                      tol_tail = 1e-2) {
  Theta_list <- NULL
  if (inherits(ss, "dsge_ss")) {
    Theta_list <- ss$Theta_list
    if (is.null(rho)) rho <- ss$rho_vec
    if (is.null(q))   q   <- ss$q
  } else if (is.numeric(ss)) {
    if (is.null(rho)) rho <- ss              # bare rho vector in first slot
  } else if (is.list(ss)) {
    Theta_list <- ss
  }
  if (is.null(rho))
    stop("hank_theta_boundary_check: no per-shock `rho` available.")
  if (!is.null(Theta_list)) {
    if (is.null(T_h)) T_h <- nrow(as.matrix(Theta_list[[1L]]))
    if (is.null(names(rho)) && length(rho) == length(Theta_list))
      names(rho) <- names(Theta_list)
  }
  if (is.null(T_h))
    stop("hank_theta_boundary_check: `T_h` is required when no Theta ",
         "matrices are supplied.")
  if (is.null(q)) q <- T_h
  q <- min(q, T_h)

  shocks <- names(rho)
  if (is.null(shocks)) shocks <- paste0("shock", seq_along(rho))
  psi_tail <- rep(NA_real_, length(rho))
  if (!is.null(Theta_list)) {
    psi_tail <- vapply(seq_along(rho), function(i) {
      Theta <- Theta_list[[shocks[i]]]
      if (is.null(Theta)) return(NA_real_)
      Theta <- as.matrix(Theta)
      q_z <- min(q, nrow(Theta))
      Th  <- Theta[seq_len(q_z), , drop = FALSE]
      Psi <- Th - rho[[i]] * rbind(rep(0, ncol(Th)), Th[-q_z, , drop = FALSE])
      tail_rows <- max(1L, ceiling(0.9 * q_z)):q_z
      denom <- max(abs(Psi))
      if (denom == 0) 0 else max(abs(Psi[tail_rows, , drop = FALSE])) / denom
    }, numeric(1))
  }

  lev <- abs(rho)^T_h
  data.frame(shock = shocks,
             rho = as.numeric(rho),
             T_h = T_h,
             q = q,
             boundary_level = as.numeric(lev),
             psi_tail = psi_tail,
             flag_boundary = as.numeric(lev) > tol_boundary,
             flag_tail = !is.na(psi_tail) & psi_tail > tol_tail,
             row.names = NULL)
}
