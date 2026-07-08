## R/hank-panel.R
## --------------------------------------------------------------------------
## Per-household PANEL log-likelihood for the discount-factor-mixture HANK
## (P1 hand-off brief, references/level-vs-response-paper/PANEL_LIKELIHOOD_BRIEF.md,
## deliverables 1-3; deliverable 3 -- integrating out the latent AGGREGATE
## price path -- is hank_mixture_panel_loglik_marginal at the bottom of this
## file: plain Monte Carlo over prior path particles, with per-period SISR
## variance control explicitly left to P-panel research).
##
## THE KEY SIMPLIFICATION (the brief's framing): dynhr discretizes household
## state onto a FINITE (e, a) grid (n_e income states x n_a asset gridpoints,
## n_cell = n_e*n_a cells; see R/hank-distribution.R's cell-order convention).
## A household's latent trajectory therefore lives on a finite state space and
## its likelihood conditional on a known aggregate price path is an EXACT
## discrete-HMM forward filter -- a sparse matrix-vector recursion over
## n_cell states -- not a particle filter. The transition kernel for that HMM
## is EXACTLY hank_forward_operator()'s Lambda (income transition Pi composed
## with the type's savings-policy lottery weights onto the a-grid), so this
## file is pure PLUMBING atop existing machinery: no lottery/interpolation
## logic is re-derived here.
##
## STATE ORDERING: reuses R/hank-distribution.R's convention throughout --
## cell index (e,a) = (e-1)*n_a + a (income OUTER/slow, asset INNER/fast); a
## distribution/filter-mass vector has length n_e*n_a; Lambda is row-
## stochastic with d_next = t(Lambda) %*% d.
##
## FORWARD FILTER (scaled, log-space): standard scaled forward recursion,
##   alpha_1 = a0 .* b_1(y_1);  c_1 = sum(alpha_1);  alpha_1 <- alpha_1 / c_1
##   alpha_t = (t(Lambda_t) %*% alpha_{t-1}) .* b_t(y_t);  c_t = sum(alpha_t); ...
##   loglik = sum(log(c_t))
## which is algebraically identical to the unnormalized recursion but keeps
## every alpha_t a proper (unit-mass) distribution, avoiding underflow over
## long T (see the underflow certification test at T = 200).
##
## OBSERVATION MODEL: obs_spec is a character vector naming which household
## variables are observed, drawn from c("a", "c") (mirrors hank_state_space()'s
## `observables` character-vector convention, R/hank-kalman.R). me_var is the
## per-observable measurement-error VARIANCE: either a single scalar (recycled
## across obs_spec, mirroring hank_kalman_loglik()'s scalar me_var,
## R/hank-kalman.R) or a named numeric vector keyed by the same names as
## obs_spec (needed here because "a" and "c" live on very different scales and
## a single shared variance is rarely appropriate). Observation weights are
## Gaussian, independent across observed variables conditional on the (e,a)
## cell: b_t(y_t)[cell] = prod_{v in obs_spec} dnorm(y_t[v]; x_v(cell), sqrt(me_var[v])),
## where x_v(cell) is the type's period-t policy value (a_pol/c_pol, or the
## grid level a itself for "a") at that cell. A NA entry in y_t for a given
## variable/period contributes a uniform (weight-1) factor for that variable
## (standard missing-data handling for exact discrete HMMs; also used to
## build the forward-operator-consistency oracle test, which sets ALL
## observations missing so the filter reduces to pure kernel propagation).
## --------------------------------------------------------------------------


#' Gaussian observation weights for one period of a household HMM
#'
#' @param y_t Named numeric vector (subset of \code{names(obs_spec)}, possibly
#'   with \code{NA} entries for missing variables this period).
#' @param obs_spec Character vector of observed variable names, subset of
#'   \code{c("a", "c")}.
#' @param me_sd Named numeric vector of measurement-error STANDARD DEVIATIONS,
#'   keyed by the entries of \code{obs_spec}.
#' @param a_level Length-\code{n_cell} vector: the asset GRID level at each
#'   cell (only needed if \code{"a" \%in\% obs_spec}).
#' @param c_pol_vec Length-\code{n_cell} vector: the period's consumption
#'   policy at each cell, distribution-order (only needed if
#'   \code{"c" \%in\% obs_spec}).
#' @return Length-\code{n_cell} numeric vector of observation weights
#'   (product of per-variable Gaussian densities; uniform 1 where all
#'   observed variables are \code{NA} this period).
#' @keywords internal
.hank_hh_obs_weight <- function(y_t, obs_spec, me_sd, a_level, c_pol_vec) {
  n_cell <- length(a_level)
  w <- rep(1, n_cell)
  for (v in obs_spec) {
    yv <- y_t[[v]]
    if (is.null(yv) || is.na(yv)) next
    x <- if (v == "a") a_level else c_pol_vec
    w <- w * stats::dnorm(yv, mean = x, sd = me_sd[[v]])
  }
  w
}


#' Normalize \code{me_var} to a named per-\code{obs_spec} standard-deviation vector
#' @keywords internal
.hank_hh_me_sd <- function(me_var, obs_spec) {
  if (is.null(names(me_var))) {
    if (length(me_var) == 1L) {
      me_var <- setNames(rep(me_var, length(obs_spec)), obs_spec)
    } else if (length(me_var) == length(obs_spec)) {
      me_var <- setNames(me_var, obs_spec)
    } else {
      stop("hank_hh_hmm_loglik: unnamed 'me_var' must have length 1 or ",
           "length(obs_spec).")
    }
  }
  missing_v <- setdiff(obs_spec, names(me_var))
  if (length(missing_v))
    stop("hank_hh_hmm_loglik: 'me_var' is missing entries for: ",
         paste(missing_v, collapse = ", "))
  if (any(me_var[obs_spec] <= 0) || any(!is.finite(me_var[obs_spec])))
    stop("hank_hh_hmm_loglik: 'me_var' entries must be finite and > 0.")
  sqrt(me_var[obs_spec])
}


#' Resolve a household type's per-period policies over a horizon from an
#' \code{aggregate_path} spec
#'
#' Constant/stationary path (\code{aggregate_path = NULL} or a list with no
#' \code{r_path}/\code{w_path}): every period reuses the type's OWN steady-
#' state policy/forward-operator objects (\code{blocks_k$a}, \code{blocks_k$c},
#' \code{blocks_k$Lambda}) -- no re-solve. Time-varying path (\code{aggregate_path}
#' has \code{r_path}/\code{w_path}): the per-period policies come from
#' \code{\link{hank_value_transition}(blocks_k, r_path, w_path, T_h,
#' keep_policies = TRUE)} (or, if the caller only needs policies and not
#' values, this reuses the identical backward-EGM recursion
#' \code{\link{hank_td_nonlinear}} performs) -- reused verbatim rather than
#' re-deriving the backward EGM pass.
#'
#' @param blocks_k A \code{\link{hank_het_block}} (single type).
#' @param aggregate_path \code{NULL} (steady state) or a list with numeric
#'   \code{r_path}, \code{w_path} (length \code{T}, level paths; see
#'   \code{\link{hank_td_nonlinear}}).
#' @param T Integer horizon (number of periods observed for this household).
#' @return A list with \code{a_pol}, \code{c_pol} (each length-\code{T} lists
#'   of \code{n_e x n_a} matrices) and \code{Lambda} (length-\code{T} list of
#'   sparse forward operators, no-transpose convention).
#' @keywords internal
.hank_hh_type_path <- function(blocks_k, aggregate_path, T) {
  n_e <- blocks_k$n_e; n_a <- blocks_k$n_a
  is_stationary <- is.null(aggregate_path) ||
    (is.null(aggregate_path$r_path) && is.null(aggregate_path$w_path))
  if (is_stationary) {
    a_pol  <- rep(list(blocks_k$a), T)
    c_pol  <- rep(list(blocks_k$c), T)
    Lambda <- rep(list(blocks_k$Lambda), T)
    return(list(a_pol = a_pol, c_pol = c_pol, Lambda = Lambda))
  }
  r_path <- aggregate_path$r_path
  w_path <- aggregate_path$w_path
  if (is.null(r_path)) r_path <- rep(blocks_k$r, T)
  if (is.null(w_path)) w_path <- rep(blocks_k$w, T)
  stopifnot(length(r_path) == T, length(w_path) == T)
  vt <- hank_td_nonlinear(blocks_k, r_path = r_path, w_path = w_path,
                          T_h = T, keep_policies = TRUE)
  list(a_pol = vt$a_pol, c_pol = vt$c_pol, Lambda = vt$Lambda)
}


#' Per-household discrete-HMM trajectory log-likelihood (one candidate type)
#'
#' The EXACT log-likelihood of one household's observed trajectory
#' \code{y_i = y_{i,1:T}}, conditional on a known/candidate aggregate price
#' path and a fixed discount-factor type \code{beta_k}, computed by a scaled
#' forward filter over the household's discretized \code{(e, a)} state
#' space. Because dynhr discretizes household state onto a finite
#' \code{n_e * n_a}-cell grid, this is exact matrix-vector recursion (no
#' particle filter, no deconvolution): the transition kernel is exactly
#' \code{\link{hank_forward_operator}}'s \code{Lambda} (income transition
#' \code{Pi} composed with the type's savings-policy lottery weights), reused
#' verbatim rather than re-derived.
#'
#' The forward recursion (scaled to avoid underflow over long \code{T}):
#' \deqn{\alpha_1 \propto a_0 \odot b_1(y_1), \qquad
#'       \alpha_t \propto (\Lambda_{t-1}^\top \alpha_{t-1}) \odot b_t(y_t),}
#' with each \eqn{\alpha_t} renormalized to sum to 1 after multiplying by its
#' observation weight, and the log-likelihood accumulated as the sum of the
#' log-normalizing-constants \code{log(c_t)}; \code{sum(alpha_T) == 1}
#' identically, so \code{loglik} is recovered purely from the scaling
#' constants (see Details in the source for the exact algebra).
#'
#' \code{aggregate_path}: \code{NULL} (or a list with no \code{r_path}/
#' \code{w_path}) means a CONSTANT/steady-state path -- every period reuses
#' the type's own steady-state policy and forward operator
#' (\code{blocks_k$a}, \code{blocks_k$c}, \code{blocks_k$Lambda}), no re-solve.
#' A time-varying path is given as \code{list(r_path = ..., w_path = ...)}
#' (level paths, length \code{T}); the per-period policies then come from
#' \code{\link{hank_td_nonlinear}}'s \code{keep_policies = TRUE} backward-EGM
#' pass (already exposes \code{a_pol}/\code{c_pol}/\code{Lambda} per period --
#' reused, not re-derived).
#'
#' This function conditions on a KNOWN aggregate path; integrating out a
#' LATENT aggregate path via the existing particle-filter machinery
#' (\code{tpf_propagate_particles} / \code{ppf_likelihood}) is a distinct,
#' deferred deliverable -- \code{aggregate_path} is designed so a particle
#' (one aggregate-path draw) can be passed here as the observation density at
#' that particle without any change to this function's contract.
#'
#' @param y_i A \code{T x length(obs_spec)} numeric matrix (or data frame)
#'   with column names matching \code{obs_spec}; entries may be \code{NA} for
#'   missing observations (contribute a uniform observation weight for that
#'   variable/period).
#' @param obs_spec Character vector of observed household variable names,
#'   subset of \code{c("a", "c")} (asset holdings / consumption).
#' @param blocks_k A \code{\link{hank_het_block}}: the single candidate type
#'   this trajectory's likelihood is evaluated under.
#' @param Pi Numeric \code{n_e x n_e} income transition matrix (must equal
#'   \code{blocks_k$Pi}; passed explicitly per the brief's signature, and
#'   checked for consistency).
#' @param aggregate_path \code{NULL} for the steady-state path, or a list
#'   \code{list(r_path, w_path)} of length-\code{T} level paths for a
#'   time-varying path (see Details).
#' @param me_var Measurement-error variance: a scalar (recycled across
#'   \code{obs_spec}) or a named numeric vector keyed by \code{obs_spec}.
#' @param a0_dist Optional length-\code{n_e*n_a} initial distribution over
#'   \code{(e, a)} cells for period 1 (defaults to the type-conditional
#'   ergodic distribution \code{blocks_k$D}, i.e. \code{hank_mixture_dist}'s
#'   per-type ergodic distribution -- the correct initialization per the
#'   brief's "selection / initial conditions" identification point).
#'
#' @return A list with \code{loglik} (scalar) and \code{alpha_T} (the final
#'   filtered, normalized distribution over \code{(e,a)} cells, length
#'   \code{n_e*n_a}).
#' @seealso \code{\link{hank_forward_operator}}, \code{\link{hank_mixture_dist}},
#'   \code{\link{hank_td_nonlinear}}, \code{\link{hank_mixture_panel_loglik}}
#' @export
hank_hh_hmm_loglik <- function(y_i, obs_spec, blocks_k, Pi, aggregate_path,
                               me_var, a0_dist = NULL) {
  if (!inherits(blocks_k, "hank_het_block"))
    stop("hank_hh_hmm_loglik: 'blocks_k' must be a hank_het_block.")
  obs_spec <- match.arg(obs_spec, c("a", "c"), several.ok = TRUE)
  if (!isTRUE(all.equal(unname(Pi), unname(blocks_k$Pi))))
    stop("hank_hh_hmm_loglik: 'Pi' must equal blocks_k$Pi.")

  y_i <- as.matrix(y_i)
  if (!all(obs_spec %in% colnames(y_i)))
    stop("hank_hh_hmm_loglik: 'y_i' column names must include all of obs_spec.")
  T_h <- nrow(y_i)
  if (T_h < 1L) stop("hank_hh_hmm_loglik: 'y_i' must have at least one row.")

  n_e <- blocks_k$n_e; n_a <- blocks_k$n_a
  n_cell <- n_e * n_a
  me_sd <- .hank_hh_me_sd(me_var, obs_spec)

  a0 <- if (is.null(a0_dist)) blocks_k$D else a0_dist
  if (length(a0) != n_cell)
    stop("hank_hh_hmm_loglik: 'a0_dist' must have length n_e*n_a.")
  a0 <- a0 / sum(a0)

  a_level_cell <- rep(blocks_k$a_grid, each = n_e)  # distribution order: (e,a)

  paths <- .hank_hh_type_path(blocks_k, aggregate_path, T_h)

  loglik <- 0
  alpha <- a0
  for (t in seq_len(T_h)) {
    if (t > 1L) alpha <- as.numeric(Matrix::t(paths$Lambda[[t - 1L]]) %*% alpha)

    c_pol_vec <- .hank_mat_to_vec(paths$c_pol[[t]])
    y_t <- as.list(y_i[t, obs_spec, drop = TRUE])
    names(y_t) <- obs_spec
    w_t <- .hank_hh_obs_weight(y_t, obs_spec, me_sd, a_level_cell, c_pol_vec)

    alpha <- alpha * w_t
    ct <- sum(alpha)
    if (!is.finite(ct) || ct <= 0)
      stop("hank_hh_hmm_loglik: filter mass collapsed to zero at t = ", t,
           " (check me_var / grid coverage of observed values).")
    alpha <- alpha / ct
    loglik <- loglik + log(ct)
  }

  list(loglik = loglik, alpha_T = alpha)
}


#' Per-household PANEL log-likelihood under a discount-factor mixture
#'
#' Mixes \code{\link{hank_hh_hmm_loglik}} over the type quadrature
#' \code{\{beta_k, omega_k\}} of a discount-factor mixture \code{G_eta} (per
#' household: \code{log sum_k omega_k * p(y_i | beta_k, Theta)}, via
#' log-sum-exp for numerical stability) and sums the per-household
#' log-likelihoods across the panel (households are conditionally
#' independent given the aggregate path \code{Theta}, per the brief).
#'
#' @param panel A list of length \code{M}, one entry per household, each a
#'   \code{T_i x length(obs_spec)} numeric matrix/data frame as expected by
#'   \code{\link{hank_hh_hmm_loglik}}'s \code{y_i} (households may have
#'   different \code{T_i}).
#' @param obs_spec Character vector of observed household variable names,
#'   subset of \code{c("a", "c")}.
#' @param blocks_by_type A list of length \code{K} \code{\link{hank_het_block}}
#'   objects (one per mixture type; as returned by
#'   \code{\link{hank_mixture_blocks}}), sharing a grid.
#' @param omega Numeric length-\code{K} mixture weights, non-negative,
#'   summing to 1 (see \code{\link{hank_mixture_dist}}'s identical
#'   convention).
#' @param Pi Numeric \code{n_e x n_e} income transition matrix shared by every
#'   type (must equal each \code{blocks_by_type[[k]]$Pi}).
#' @param aggregate_path \code{NULL} for the steady-state path, or a list
#'   \code{list(r_path, w_path)} shared across households and types (see
#'   \code{\link{hank_hh_hmm_loglik}}). Conditions on a KNOWN aggregate path;
#'   integrating out a LATENT path via the particle filter is deliverable 3
#'   (deferred -- not implemented here). A future
#'   \code{hank_mixture_panel_loglik_marginal()} can call this function once
#'   per particle (one \code{aggregate_path} draw) and combine the per-
#'   particle household log-likelihoods with particle weights, without any
#'   change to this function's contract.
#' @param me_var Measurement-error variance: scalar or named vector keyed by
#'   \code{obs_spec} (see \code{\link{hank_hh_hmm_loglik}}).
#'
#' @return A list with \code{loglik} (scalar, summed over households),
#'   \code{loglik_i} (length-\code{M} per-household log-likelihoods), and
#'   \code{loglik_ik} (\code{M x K} matrix of per-household-per-type
#'   log-likelihoods, before mixing -- useful for posterior type
#'   probabilities / diagnostics).
#' @seealso \code{\link{hank_hh_hmm_loglik}}, \code{\link{hank_mixture_blocks}},
#'   \code{\link{hank_mixture_dist}}
#' @export
hank_mixture_panel_loglik <- function(panel, obs_spec, blocks_by_type, omega,
                                      Pi, aggregate_path, me_var) {
  .hank_mixture_check_omega(blocks_by_type, omega)
  M <- length(panel)
  K <- length(blocks_by_type)
  if (M < 1L) stop("hank_mixture_panel_loglik: 'panel' must have >= 1 household.")

  loglik_ik <- matrix(NA_real_, M, K)
  for (i in seq_len(M)) {
    for (k in seq_len(K)) {
      loglik_ik[i, k] <- hank_hh_hmm_loglik(
        y_i = panel[[i]], obs_spec = obs_spec, blocks_k = blocks_by_type[[k]],
        Pi = Pi, aggregate_path = aggregate_path, me_var = me_var,
        a0_dist = NULL)$loglik
    }
  }

  log_omega <- log(omega)
  loglik_i <- apply(loglik_ik, 1L, function(row) {
    m <- max(row + log_omega)
    m + log(sum(exp(row + log_omega - m)))
  })

  list(loglik = sum(loglik_i), loglik_i = loglik_i, loglik_ik = loglik_ik)
}


#' Draw aggregate price-path particles from a stationary AR(1) TFP prior
#' through the linear KS GE map
#'
#' @param ks A \code{\link{hank_ks_steady}} object.
#' @param rho,sigma AR(1) persistence and innovation std of the (single) TFP
#'   shock, deviations from \code{ks$Z}.
#' @param n_particles,T_h Number of path draws / horizon.
#' @param ge Optional precomputed \code{\link{hank_ks_ge_jacobian}(ks, T_h)}.
#' @return Length-\code{n_particles} list of \code{list(r_path, w_path)}
#'   level paths (each length \code{T_h}).
#' @keywords internal
.hank_panel_ks_path_draws <- function(ks, rho, sigma, n_particles, T_h,
                                      ge = NULL) {
  if (!inherits(ks, "hank_ks"))
    stop("hank_mixture_panel_loglik_marginal: 'ks' must be a hank_ks steady ",
         "state (see hank_ks_steady) when 'shock_specs' is an AR(1) spec.")
  if (is.null(ge)) ge <- hank_ks_ge_jacobian(ks, T_h)
  lapply(seq_len(n_particles), function(p) {
    ## Stationary AR(1) in deviations (the economy sits in its stochastic
    ## steady state before the sample starts -- consistent with the
    ## household filter's ergodic-D0 initialization).
    dZ <- numeric(T_h)
    dZ[1L] <- stats::rnorm(1L, 0, sigma / sqrt(1 - rho^2))
    if (T_h > 1L)
      for (t in 2L:T_h) dZ[t] <- rho * dZ[t - 1L] + stats::rnorm(1L, 0, sigma)
    irf <- hank_ks_linear_irf(ks, dZ, ge = ge)
    list(r_path = ks$r + irf$dr, w_path = ks$w + irf$dw)
  })
}


#' Panel log-likelihood with the latent aggregate path integrated out
#' (particle / Monte Carlo marginalization)
#'
#' The P1 panel brief's deliverable 3: the full-panel likelihood
#' \deqn{L(\eta) = E_\Theta\left[\prod_i L_i(\eta \mid \Theta)\right]}
#' with the common latent aggregate price path \eqn{\Theta = \{r_t, w_t\}}
#' integrated out by simple Monte Carlo over path particles: each particle is
#' one draw of the aggregate path from its prior; the per-particle panel
#' log-likelihood is \code{\link{hank_mixture_panel_loglik}} evaluated at
#' that path (households are conditionally independent given the path); and
#' the marginal is the equal-weight log-mean-exp across particles.
#'
#' The estimator of \eqn{L} is unbiased for any \code{n_particles}; the
#' returned log-likelihood (like every particle-filter log-likelihood) has a
#' downward finite-sample Jensen bias that vanishes as
#' \code{n_particles} grows. Monitor the returned \code{ess}: with many
#' households/periods the product likelihood concentrates on few paths and
#' the plain prior sampler degenerates -- per-period sequential importance
#' resampling (a Rao-Blackwellized particle filter re-using the household
#' HMM's per-period normalizing constants as the incremental weights) is the
#' known variance-control upgrade, deliberately NOT implemented here
#' (P-panel research scope).
#'
#' @param panel,obs_spec,blocks_by_type,omega,Pi,me_var As in
#'   \code{\link{hank_mixture_panel_loglik}}.
#' @param shock_specs The aggregate-path prior. Either
#'   \itemize{
#'     \item a \strong{function} \code{function(n_particles, T_h)} returning a
#'       length-\code{n_particles} list of \code{aggregate_path} objects
#'       (each \code{list(r_path, w_path)} level paths, or \code{NULL} for
#'       the steady-state path) -- fully model-agnostic; or
#'     \item a one-shock AR(1) spec in the package's \code{shock_specs}
#'       convention (\code{list(Z = list(rho = , sigma = ))}, see
#'       \code{\link{hank_state_space}}), in which case \code{ks} is required
#'       and TFP deviation paths are drawn from the stationary AR(1) and
#'       mapped to \code{(r, w)} level paths through the linear KS GE map
#'       \code{\link{hank_ks_linear_irf}} (Jacobian built once, or passed
#'       via \code{ge}).
#'   }
#' @param T_h Horizon of the drawn paths; defaults to the longest household
#'   trajectory in \code{panel}.
#' @param n_particles Number of aggregate-path particles.
#' @param ks,ge Only used with the AR(1) spec form of \code{shock_specs}:
#'   the \code{\link{hank_ks_steady}} object and (optionally) its
#'   precomputed \code{\link{hank_ks_ge_jacobian}}.
#' @param seed Optional integer; \code{set.seed(seed)} before drawing paths.
#'
#' @return A list with \code{loglik} (the marginal log-likelihood estimate),
#'   \code{loglik_p} (length-\code{n_particles} per-particle panel
#'   log-likelihoods), \code{ess} (effective sample size of the normalized
#'   particle weights, in \code{[1, n_particles]}), and \code{n_particles}.
#' @seealso \code{\link{hank_mixture_panel_loglik}},
#'   \code{\link{hank_ks_linear_irf}}
#' @export
hank_mixture_panel_loglik_marginal <- function(panel, obs_spec, blocks_by_type,
                                               omega, Pi, shock_specs, me_var,
                                               T_h = NULL, n_particles = 100L,
                                               ks = NULL, ge = NULL,
                                               seed = NULL) {
  if (is.null(T_h)) T_h <- max(vapply(panel, NROW, integer(1)))
  n_particles <- as.integer(n_particles)
  if (n_particles < 1L)
    stop("hank_mixture_panel_loglik_marginal: 'n_particles' must be >= 1.")
  if (!is.null(seed)) set.seed(seed)

  if (is.function(shock_specs)) {
    paths <- shock_specs(n_particles, T_h)
    if (!is.list(paths) || length(paths) != n_particles)
      stop("hank_mixture_panel_loglik_marginal: a function 'shock_specs' ",
           "must return a length-n_particles list of aggregate_path objects.")
  } else if (is.list(shock_specs) && length(shock_specs) == 1L &&
             all(c("rho", "sigma") %in% names(shock_specs[[1L]]))) {
    if (is.null(ks))
      stop("hank_mixture_panel_loglik_marginal: the AR(1) spec form of ",
           "'shock_specs' requires 'ks' (a hank_ks_steady object).")
    spec <- shock_specs[[1L]]
    paths <- .hank_panel_ks_path_draws(ks, rho = spec$rho, sigma = spec$sigma,
                                       n_particles = n_particles, T_h = T_h,
                                       ge = ge)
  } else {
    stop("hank_mixture_panel_loglik_marginal: 'shock_specs' must be a ",
         "function(n_particles, T_h) returning path draws, or a one-shock ",
         "list(<name> = list(rho = , sigma = )) AR(1) spec.")
  }

  loglik_p <- vapply(paths, function(pth) {
    hank_mixture_panel_loglik(panel = panel, obs_spec = obs_spec,
                              blocks_by_type = blocks_by_type, omega = omega,
                              Pi = Pi, aggregate_path = pth,
                              me_var = me_var)$loglik
  }, numeric(1))

  ## Equal-weight prior particles: log L-hat = logmeanexp(loglik_p).
  m <- max(loglik_p)
  loglik <- m + log(mean(exp(loglik_p - m)))
  w <- exp(loglik_p - m); w <- w / sum(w)
  ess <- 1 / sum(w^2)

  list(loglik = loglik, loglik_p = loglik_p, ess = ess,
       n_particles = n_particles)
}
