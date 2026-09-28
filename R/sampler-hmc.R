## R/sampler-hmc.R
## --------------------------------------------------------------------------
## Phase-2 split from nuts-monolith.R.
##
## Internals: .hmc_gradient, .hmc_leapfrog, .hmc_kinetic,
##             .hmc_sample_momentum, .hmc_find_stepsize
## Sampler:   dynhr_hmc() -- vanilla HMC with fixed step size
## Diagnostics: hmc_summary() -- summary table for HMC and NUTS results
## --------------------------------------------------------------------------

.hmc_gradient <- function(f, theta, method = "forward", eps = 1e-5) {
  nms <- names(theta)
  d   <- length(theta)

  # All finite-difference methods here use a RELATIVE step h = eps*max(|theta|,
  # 1e-3). A FIXED absolute step (numDeriv's "simple" default, h ~ 1e-4) is
  # wildly inaccurate for parameters whose magnitude is << 1 (e.g. a shock std
  # dev ~ 1e-3): on a Kalman log-likelihood the forward difference can come out
  # with the WRONG SIGN, which destroys leapfrog energy conservation -> NUTS
  # dual-averaging shrinks eps without bound, trees hit max depth, the sampler
  # hangs. The relative step matches numDeriv's Richardson gradient to ~1e-4.
  #
  #   "forward"    -- d+1 evals (reuses f(theta)); cheapest, and its O(h) bias is
  #                   negligible vs the gradient magnitude, so leapfrog energy
  #                   conservation is essentially identical to central. Default.
  #   "central"    -- 2d evals; O(h^2) accuracy. Use if forward looks too noisy.
  #   "Richardson" -- numDeriv extrapolation, ~4-6 evals/param; most accurate.
  if (identical(method, "Richardson") &&
      requireNamespace("numDeriv", quietly = TRUE)) {
    g <- numDeriv::grad(func = f, x = theta, method = "Richardson")
  } else if (identical(method, "central")) {
    g <- numeric(d)
    for (i in seq_len(d)) {
      h <- eps * max(abs(theta[i]), 1e-3)
      theta_p <- theta_m <- theta
      theta_p[i] <- theta[i] + h
      theta_m[i] <- theta[i] - h
      fp <- f(theta_p)
      fm <- f(theta_m)
      g[i] <- if (is.finite(fp) && is.finite(fm)) (fp - fm) / (2 * h) else 0
    }
  } else {
    f0 <- f(theta)
    if (!is.finite(f0)) f0 <- -1e300
    g <- numeric(d)
    for (i in seq_len(d)) {
      h <- eps * max(abs(theta[i]), 1e-3)
      theta_p <- theta
      theta_p[i] <- theta[i] + h
      fp <- f(theta_p)
      g[i] <- if (is.finite(fp)) (fp - f0) / h else 0
    }
  }

  g[!is.finite(g)] <- 0
  names(g) <- nms
  g
}


## ---------------------------------------------------------------------------
## W92: fused log-posterior + gradient.
##
## make_posterior_grad()'s Gaussian closures carry attr(grad_fn, "logpost_grad"),
## function(theta) -> list(logpost, grad), which returns the log-posterior the
## gradient pass already computed. With it a gradient sampler needs ONE
## evaluation per new position instead of a gradient call plus a separate
## log-posterior call, and -- since the value arrives with the gradient -- the
## gradient at a trajectory's current end point is CARRIED into the next
## leapfrog step instead of being recomputed there (.hmc_leapfrog() used to
## evaluate the gradient at both ends of every step: two gradients and one
## log-posterior per new point, now one fused evaluation).
##
## .hmc_fused_target() decides whether a sampler fuses. It refuses (NULL: the
## sampler runs its separate-call path exactly as before) unless grad_fn
## carries the companion AND the companion's log-posterior equals the
## sampler's own log_post_fn at the start point to 1e-10 relative, both
## finite. The check pairs the two closures: a companion built for another
## target (another power, lik_init, data, a system prior the gradient does not
## carry) differs there. When it fuses, the chain's target is the companion
## everywhere -- the start point's lp_curr included -- so the chain targets
## ONE deterministic function.
##
## In sampler space (eta under `transform`) the value and gradient are
## assembled exactly as make_transformed_logpost() and make_transformed_grad()
## assemble them (Jacobian added once to a finite value; ONE chain-rule
## application, W86), and a non-finite value is floored to -1e300 as the
## samplers' .lp_scalar() floors it.
##
## Returns NULL or list(vg = function(state) -> list(lp, grad), lp0, g0), with
## lp0 / g0 the value and gradient at `state_init`.
## ---------------------------------------------------------------------------
.hmc_fused_target <- function(log_post_fn, grad_fn, state_init, par_names,
                              transform = NULL, verbose = FALSE,
                              sampler = "sampler") {
  if (is.null(grad_fn)) return(NULL)
  vg_theta <- attr(grad_fn, "logpost_grad", exact = TRUE)
  if (!is.function(vg_theta)) return(NULL)
  ## theta-space list(logpost, grad) at sampler state `state` -> sampler space
  to_state <- if (is.null(transform)) {
    function(state, r) {
      lp <- r$logpost
      list(lp = if (is.finite(lp)) lp else -1e300, grad = r$grad)
    }
  } else {
    function(state, r) {
      lp <- r$logpost
      if (is.finite(lp)) lp <- lp + transform$log_jacobian(state)
      list(lp   = if (is.finite(lp)) lp else -1e300,
           grad = transform$dtheta_deta(state) * r$grad +
                  transform$dlog_jacobian(state))
    }
  }
  to_theta <- if (is.null(transform)) function(state) state
              else function(state) transform$to_constrained(state)
  vg <- function(state) {
    names(state) <- par_names
    to_state(state, vg_theta(to_theta(state)))
  }
  names(state_init) <- par_names
  theta0 <- to_theta(state_init)
  ref    <- log_post_fn(theta0)
  lp_ref <- if (is.list(ref)) ref$logpost else ref
  f0     <- vg_theta(theta0)
  ok <- length(lp_ref) == 1L && is.finite(lp_ref) &&
    length(f0$logpost) == 1L && is.finite(f0$logpost) &&
    abs(f0$logpost - lp_ref) <= 1e-10 * max(1, abs(lp_ref))
  if (!ok) {
    if (isTRUE(verbose))
      .dynhr_inform(sprintf(paste0(
        "%s: the gradient's fused log-posterior (%s) is not the sampler's ",
        "log-posterior (%s) at the start point; sampling with separate ",
        "log-posterior and gradient calls."),
        sampler, format(f0$logpost, digits = 15), format(lp_ref, digits = 15)),
        class = "dynhr_message_fuse_refused")
    return(NULL)
  }
  e0 <- to_state(state_init, f0)   # == vg(state_init), without re-evaluating
  list(vg = vg, lp0 = e0$lp, g0 = e0$grad)
}


#' Single leapfrog step
#'
#' Supports two mass-matrix modes:
#'   - Diagonal (legacy): pass \code{M_inv_diag} (numeric vector), leave
#'     \code{M_inv} NULL. Results are bit-identical to the pre-dense code.
#'   - Dense: pass \code{M_inv} (d×d matrix) and \code{chol_M} (upper
#'     Cholesky of M = solve(M_inv)), leave \code{M_inv_diag} NULL.
#'
#' @param theta Position (named numeric)
#' @param r Momentum (numeric, same length)
#' @param eps Step size (scalar)
#' @param grad_fn function(theta) -> gradient vector
#' @param M_inv_diag Inverse mass matrix diagonal (numeric vector; diagonal path)
#' @param M_inv Dense inverse mass matrix (d×d matrix; dense path)
#' @param g0 The gradient at \code{theta}, carried from the evaluation that
#'   produced it (W92 fused path; W94 the separate-call path too -- the start
#'   of a step is the previous step's end); NULL: computed with \code{grad_fn}.
#' @param vg_fn Fused path (W92): \code{function(theta) -> list(lp, grad)};
#'   when supplied the end point is evaluated with it (one call for value and
#'   gradient) and the result carries \code{lp}.
#' @return list(theta, r, g, lp) after one leapfrog step (\code{g} the
#'   gradient at the new position; \code{lp} its log density on the fused
#'   path, NULL otherwise), or NULL when the step left the well-defined
#'   region (non-finite gradient or position) -- the callers treat NULL as a
#'   divergence.
#' @noRd
.hmc_leapfrog <- function(theta, r, eps, grad_fn, M_inv_diag, M_inv = NULL,
                          g0 = NULL, vg_fn = NULL) {
  nms <- names(theta)
  g <- if (is.null(g0)) grad_fn(theta) else g0
  ## Bail out (NULL -> divergence upstream) BEFORE evaluating the gradient
  ## at a non-finite position: gradient closures call the log-posterior
  ## internally, and an NA/Inf theta crashes prior bounds checks rather
  ## than registering as a divergence. Non-finite gradients arise
  ## legitimately at the edge of the solvable region (BK violations, failed
  ## perturbed steady states in the implicit gradient's FD fallback).
  if (any(!is.finite(g))) return(NULL)
  r <- r + 0.5 * eps * g
  if (is.null(M_inv)) {
    ## Diagonal path (bit-identical to pre-dense behaviour)
    theta <- theta + eps * M_inv_diag * r
  } else {
    ## Dense path: θ += ε M⁻¹ r  (M_inv: d×d matrix, or a
    ## "dynhr_lowrank_metric" object applied in O(d k); sampler-metric.R)
    theta <- theta + eps * .metric_apply_inv(M_inv, r)
  }
  names(theta) <- nms
  if (any(!is.finite(theta))) return(NULL)
  lp <- NULL
  if (is.null(vg_fn)) {
    g <- grad_fn(theta)
  } else {
    e  <- vg_fn(theta)
    g  <- e$grad
    lp <- e$lp
  }
  if (any(!is.finite(g))) return(NULL)
  r <- r + 0.5 * eps * g
  list(theta = theta, r = r, g = g, lp = lp)
}


#' Kinetic energy
#'
#' Diagonal path: K(r) = 0.5 * sum(r^2 * M_inv_diag)
#' Dense path:    K(r) = 0.5 * rᵀ M⁻¹ r
#'
#' @param r Momentum vector
#' @param M_inv_diag Inverse mass diagonal (diagonal path; ignored if M_inv supplied)
#' @param M_inv Dense inverse mass matrix (dense path; NULL = diagonal)
#' @noRd
.hmc_kinetic <- function(r, M_inv_diag, M_inv = NULL) {
  if (is.null(M_inv)) {
    ## Diagonal path (bit-identical to pre-dense behaviour)
    0.5 * sum(r^2 * M_inv_diag)
  } else {
    ## Dense / low-rank path: 0.5 * rᵀ M⁻¹ r
    0.5 * sum(r * .metric_apply_inv(M_inv, r))
  }
}


#' Sample momentum from N(0, M)
#'
#' Diagonal path: r ~ N(0, diag(M_diag))  →  r = rnorm(d) * sqrt(M_diag)
#' Dense path:    r ~ N(0, M)             →  r = t(chol_M) %*% rnorm(d)
#'
#' Dense covariance derivation:
#'   Let U = chol_M (upper triangular, M = UᵀU).
#'   r = t(U) z  where z ~ N(0, I).
#'   Cov(r) = t(U) Cov(z) U = t(U) I U = t(U) U = M.  ✓
#'
#' @param d Dimension
#' @param M_diag Mass diagonal (diagonal path)
#' @param chol_M Upper Cholesky of M (dense path; NULL = diagonal), or a
#'   \code{"dynhr_lowrank_metric"} object (low-rank-plus-diagonal path).
#' @noRd
.hmc_sample_momentum <- function(d, M_diag, chol_M = NULL) {
  if (is.null(chol_M)) {
    ## Diagonal path (bit-identical to pre-dense behaviour)
    rnorm(d) * sqrt(M_diag)
  } else if (inherits(chol_M, "dynhr_lowrank_metric")) {
    ## Low-rank-plus-diagonal path (sampler-metric.R): p = diag(1/sd) A z
    .lowrank_sample_momentum(chol_M)
  } else {
    ## Dense path: r = t(U) z, U = chol_M (upper triangular)
    ## t(chol_M) is lower triangular; %*% rnorm(d) gives a draw from N(0,M)
    as.numeric(t(chol_M) %*% rnorm(d))
  }
}


#' Find a reasonable initial step size (Algorithm 4, Hoffman & Gelman 2014)
#'
#' Searches for eps such that acceptance probability ~ 0.5.
#' Supports diagonal (M_inv_diag / M_diag) and dense (M_inv / chol_M) paths.
#' Fused path (W92): with \code{vg_fn} (see \code{.hmc_fused_target()}) every
#' trial point is ONE value-and-gradient evaluation, and \code{lp0}/\code{g0}
#' (the value and gradient at \code{theta}, already known to the caller) are
#' not recomputed. On either path the gradient at \code{theta} (\code{g0},
#' or one \code{grad_fn} call) serves every trial (W94).
#' @noRd
.hmc_find_stepsize <- function(theta, lp_fn, grad_fn, M_inv_diag, M_diag,
                                M_inv = NULL, chol_M = NULL,
                                vg_fn = NULL, lp0 = NULL, g0 = NULL) {
  d <- length(theta)
  eps <- 1.0

  r <- .hmc_sample_momentum(d, M_diag, chol_M = chol_M)
  if (!is.null(vg_fn) && (is.null(lp0) || is.null(g0))) {
    e0 <- vg_fn(theta)
    lp0 <- e0$lp
    g0  <- e0$grad
  }
  if (is.null(vg_fn)) lp0 <- lp_fn(theta)
  ## W94: every trial below starts at theta -- the gradient there is taken
  ## ONCE (on the separate-call path too), not once per trial.
  if (is.null(g0)) g0 <- grad_fn(theta)
  H0 <- -lp0 + .hmc_kinetic(r, M_inv_diag, M_inv = M_inv)

  step <- .hmc_leapfrog(theta, r, eps, grad_fn, M_inv_diag, M_inv = M_inv,
                        g0 = g0, vg_fn = vg_fn)
  if (is.null(step) || any(!is.finite(step$theta))) {
    return(0.001)
  }

  lp1 <- if (is.null(vg_fn)) lp_fn(step$theta) else step$lp
  if (!is.finite(lp1)) lp1 <- -1e300
  H1 <- -lp1 + .hmc_kinetic(step$r, M_inv_diag, M_inv = M_inv)

  log_ratio <- -H1 + H0
  if (!is.finite(log_ratio)) return(0.001)

  a <- if (log_ratio > log(0.5)) 1 else -1

  for (k in 1:100) {
    eps_try <- eps * (2^a)
    step <- .hmc_leapfrog(theta, r, eps_try, grad_fn, M_inv_diag, M_inv = M_inv,
                          g0 = g0, vg_fn = vg_fn)
    if (is.null(step) || any(!is.finite(step$theta))) break

    lp1 <- if (is.null(vg_fn)) lp_fn(step$theta) else step$lp
    if (!is.finite(lp1)) break
    H1 <- -lp1 + .hmc_kinetic(step$r, M_inv_diag, M_inv = M_inv)
    log_ratio <- -H1 + H0
    if (!is.finite(log_ratio)) break

    if (a * log_ratio <= -a * log(2)) break
    eps <- eps_try
  }

  max(eps, 1e-10)
}


# ============================================================================
# dynhr_hmc() -- Basic HMC with fixed trajectory length
# ============================================================================

#' @param log_post_fn function(theta) -> list(logpost, loglik, logprior)
#' @param theta_init Named numeric vector of starting parameter values
#' @param n_draws Total draws including warmup
#' @param n_warmup Warmup draws (discarded from output$chain)
#' @param L Number of leapfrog steps per iteration
#' @param step_size Step size (NULL = auto-tune)
#' @param step_jitter Relative step-size jitter (default 0.2): each iteration
#'   integrates with \code{step_size * U(1 - step_jitter, 1 + step_jitter)}
#'   (Neal 2011, Sec. 5.4.2.2; Stan's \code{stepsize_jitter}). A fixed
#'   \code{L * step_size} can land near a period of the (preconditioned)
#'   target's Hamiltonian flow, and the chain then revisits near-identical
#'   orbits: on the two 2-d Gaussians of test-sampler-hmc.R (L = 25 / 20,
#'   3000 / 4000 draws), 10 of 50 seeds with no jitter missed a variance or
#'   sd by more than 15\% (worst 40\%), while jitter 0.2 held all 50 within
#'   12.5\%. The jitter is drawn independently
#'   of the state, so each iteration is still an exact (reversible) kernel.
#'   \code{0} disables it and reproduces the unjittered chain bit-for-bit.
#' @param grad_fn Optional analytical gradient function(theta) -> numeric
#' @param grad_method "simple" (forward diff) or "Richardson"
#' @param adapt_mass Whether to adapt diagonal mass matrix during warmup
#'   (once, at 70\% of warmup, from the draws after the first 20\%). Stan's
#'   convention: the INVERSE mass is set to the warmup posterior variance
#'   (before 0.9.3.50 the mass itself was, which squared the conditioning).
#'   \code{metric = "warmup_dense"} instead sets a Ledoit-Wolf dense inverse
#'   mass from the same window.
#' @param verbose Print progress messages
#' @param progressor progressr callback or NULL
#' @param chain_id Label for progress messages
#' @param transform Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), HMC runs
#'   in UNCONSTRAINED eta-space: `theta_init` is mapped to
#'   `eta0 = transform$to_unconstrained(theta_init)`, the target is
#'   `make_transformed_logpost(log_post_fn, transform, include_jacobian =
#'   TRUE)` (theta-space logpost + change-of-variables Jacobian, so eta draws
#'   are correctly distributed), and `grad_fn` (if supplied) is wrapped via
#'   `make_transformed_grad()`. If `grad_fn` is NULL, the existing numerical
#'   gradient differentiates the WRAPPED (Jacobian-included) eta-space
#'   target directly -- see `lp_scalar` below. `M_inv_diag`/`M_diag` (mass
#'   matrix) are interpreted in eta-space; the orchestrator is responsible
#'   for converting a theta-space inverse-Hessian diagonal via the delta
#'   method (see \code{run_posterior_estimation}'s `transform_params` path).
#'   The returned `chain`/`full_chain` are mapped back to theta-space via
#'   `to_constrained()`, and `logpost_trace`/`post_logpost` store the
#'   THETA-SPACE log-posterior (Jacobian subtracted back out). When NULL
#'   (default), behaviour is bit-identical to before.
#'
#' @return List compatible with rwmh() output + HMC-specific fields
#' @noRd
dynhr_hmc <- function(
    log_post_fn,
    theta_init,
    n_draws     = 5000L,
    n_warmup    = 2500L,
    L           = 25L,
    step_size   = NULL,
    grad_fn     = NULL,
    grad_method = "forward",
    adapt_mass  = TRUE,
    metric      = c("diagonal", "warmup_dense"),
    verbose     = TRUE,
    progressor  = NULL,
    chain_id    = NULL,
    transform   = NULL,
    M_inv       = NULL,
    chol_M      = NULL,
    step_jitter = 0.2
) {
  stopifnot(is.function(log_post_fn), is.numeric(theta_init))
  metric <- match.arg(metric)
  if (!is.numeric(step_jitter) || length(step_jitter) != 1L ||
      !is.finite(step_jitter) || step_jitter < 0 || step_jitter >= 1) {
    .dynhr_abort("`step_jitter` must be a single number in [0, 1).",
                 class = "dynhr_error_invalid_argument")
  }
  d <- length(theta_init)
  par_names <- names(theta_init)
  n_total <- n_draws + n_warmup

  # ---- Transformed (opt-in) target: operate on eta = to_unconstrained(theta)
  if (!is.null(transform)) {
    target_fn  <- make_transformed_logpost(log_post_fn, transform, include_jacobian = TRUE)
    state_init <- transform$to_unconstrained(theta_init)
    names(state_init) <- par_names
  } else {
    target_fn  <- log_post_fn
    state_init <- theta_init
  }

  # --- Scalar log-posterior wrapper (operates on eta when transformed) ---
  lp_scalar <- function(theta) {
    names(theta) <- par_names
    res <- target_fn(theta)
    val <- if (is.list(res)) res$logpost else res
    if (!is.finite(val)) -1e300 else val
  }

  # --- Gradient function ---
  # kernel_stats (F5): capture make_posterior_grad()'s fallback-usage
  # counters (attr(grad_fn, "kernel_stats")) from the caller-supplied
  # grad_fn BEFORE any transform-wrapping (which would drop the attribute).
  grad_fn_kernel_stats <- attr(grad_fn, "kernel_stats")
  if (is.null(grad_fn)) {
    grad <- function(theta) .hmc_gradient(lp_scalar, theta, method = grad_method)
  } else if (!is.null(transform)) {
    grad <- make_transformed_grad(grad_fn, transform)
  } else {
    grad <- grad_fn
  }
  # --- Fused value + gradient (W92): one evaluation per new position, the
  # gradient carried along the trajectory; NULL = separate calls as before.
  fz <- .hmc_fused_target(log_post_fn, grad_fn, state_init, par_names,
                          transform = transform, verbose = verbose,
                          sampler = "HMC")
  vg_fn <- if (is.null(fz)) NULL else fz$vg
  ## gradient at state_init when known (fused; W94 also the separate-call
  ## path once the step-size search has taken it)
  g_init <- if (is.null(fz)) NULL else fz$g0

  # --- Mass matrix ---
  # Dense path: M_inv (d×d) and chol_M (upper Cholesky of M = solve(M_inv))
  # are pre-supplied; we set M_diag/M_inv_diag to NULL-safe sentinels used
  # only by the diagonal path (they are ignored when M_inv is non-NULL).
  # Diagonal path: M_diag and M_inv_diag as before (identity to start).
  # When a dense metric is supplied we disable diagonal mass adaptation
  # (adapt_mass is overridden below) to preserve the caller-supplied metric.
  use_dense <- !is.null(M_inv) && is.matrix(M_inv)
  if (use_dense) {
    ## Dense: sentinel diagonals (never used in the hot path, but kept for
    ## the adapt_mass block that will be skipped)
    M_diag     <- rep(1, d)
    M_inv_diag <- NULL   ## signals: ignore diagonal path
  } else {
    M_diag     <- rep(1, d)
    M_inv_diag <- rep(1, d)
  }

  # --- Find initial step size ---
  if (is.null(step_size)) {
    if (is.null(fz)) g_init <- grad(state_init)   # W94: taken once
    step_size <- .hmc_find_stepsize(state_init, lp_scalar, grad,
                                    M_inv_diag, M_diag,
                                    M_inv = M_inv, chol_M = chol_M,
                                    vg_fn = vg_fn, lp0 = fz$lp0, g0 = g_init)
    if (verbose) .dynhr_inform(sprintf("HMC: auto step_size = %.4e", step_size))
  }

  # --- Dual-averaging step-size adaptation (Nesterov; Hoffman & Gelman 2014) ---
  # HMC runs L leapfrog steps per iteration, so a step size tuned for a SINGLE
  # step (.hmc_find_stepsize) is far too large for the full L-step trajectory --
  # the integrator goes unstable, every trajectory diverges, and the chain
  # freezes. We therefore adapt step_size during warmup to hit target_accept
  # over the actual L-step trajectory (identical machinery to dynhr_nuts()); the
  # found/supplied step_size is only the starting point. Frozen after warmup.
  target_accept <- 0.80
  mu_da    <- log(10 * step_size)   # target log step size (dual averaging)
  eps_bar  <- 1                     # averaged step size (log-scale tracking)
  H_bar    <- 0                     # averaged acceptance-stat shortfall
  gamma_da <- 0.05
  t0_da    <- 10
  kappa_da <- 0.75
  da_m     <- 0L                    # DA iteration count (reset at mass update)

  # --- Storage ---
  # `chain`/`full_chain` are always THETA-SPACE (mapped back via
  # to_constrained() when transformed). `state_chain` tracks the sampler
  # state (eta-space when transformed, theta-space otherwise) and is used
  # only for mass-matrix adaptation (which must operate in the space the
  # sampler actually moves in).
  chain         <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
  state_chain   <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
  logpost_trace <- numeric(n_total)
  accepted      <- logical(n_total)
  n_grad_evals  <- 0L

  state   <- state_init
  ## Fused (W92): the start point's value comes from the same function as
  ## every later one; g_curr is the gradient there, carried into the first
  ## leapfrog step of each trajectory.
  ## W94: on the separate-call path too the gradient at the current state is
  ## carried (NULL until first needed): the start point of a trajectory's
  ## first leapfrog step is the previous accepted end point, whose gradient
  ## the step that reached it already computed. Same values, half the calls.
  lp_curr <- if (is.null(fz)) lp_scalar(state) else fz$lp0
  g_curr  <- g_init
  trace_lp_curr <- if (!is.null(transform)) {
    lp_curr - transform$log_jacobian(state)
  } else {
    lp_curr
  }
  chain[1, ]        <- if (!is.null(transform)) transform$to_constrained(state) else state
  state_chain[1, ]  <- state
  logpost_trace[1]  <- trace_lp_curr
  accepted[1]       <- TRUE
  n_accept <- 0L
  t_start  <- Sys.time()

  for (i in 2:n_total) {
    r0 <- .hmc_sample_momentum(d, M_diag, chol_M = chol_M)
    H0 <- -lp_curr + .hmc_kinetic(r0, M_inv_diag, M_inv = M_inv)

    if (is.null(g_curr)) g_curr <- grad(state)
    state_prop <- state
    r_prop     <- r0
    g_prop     <- g_curr
    lp_last    <- NULL
    divergent  <- FALSE

    # Step-size jitter (state-independent, so the kernel stays exact); no RNG
    # draw when disabled, so step_jitter = 0 is the unjittered chain exactly.
    eps_i <- if (step_jitter > 0) step_size * (1 + step_jitter * (2 * runif(1) - 1)) else step_size
    for (l in seq_len(L)) {
      step <- .hmc_leapfrog(state_prop, r_prop, eps_i, grad,
                            M_inv_diag, M_inv = M_inv,
                            g0 = g_prop, vg_fn = vg_fn)
      if (is.null(step) || any(!is.finite(step$theta)) || any(!is.finite(step$r))) {
        divergent <- TRUE
        break
      }
      state_prop <- step$theta
      r_prop     <- step$r
      g_prop     <- step$g   # carried into the next step (W94: both paths)
      lp_last    <- step$lp  # NULL on the separate-call path
    }
    n_grad_evals <- n_grad_evals + L  # one new-point gradient per leapfrog step

    names(state_prop) <- par_names

    if (!divergent) {
      lp_prop   <- if (is.null(vg_fn)) lp_scalar(state_prop) else lp_last
      H1        <- -lp_prop + .hmc_kinetic(-r_prop, M_inv_diag, M_inv = M_inv)
      log_alpha <- -H1 + H0
      ## Metropolis acceptance probability, clamped to [0,1] -- the statistic
      ## the dual averager drives toward target_accept.
      accept_stat <- if (is.finite(log_alpha)) min(1, exp(log_alpha)) else 0

      if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        state   <- state_prop
        lp_curr <- lp_prop
        g_curr  <- g_prop
        n_accept <- n_accept + 1L
        accepted[i] <- TRUE
      } else {
        accepted[i] <- FALSE
      }
    } else {
      accept_stat <- 0      # divergent trajectory: zero acceptance
      accepted[i] <- FALSE
    }

    trace_lp_curr <- if (!is.null(transform)) {
      lp_curr - transform$log_jacobian(state)
    } else {
      lp_curr
    }
    chain[i, ]       <- if (!is.null(transform)) transform$to_constrained(state) else state
    state_chain[i, ] <- state
    logpost_trace[i] <- trace_lp_curr

    # --- Dual-averaging step-size update (warmup only) ---
    # Drive step_size toward target_accept over the full L-step trajectory.
    if (i <= n_warmup) {
      da_m      <- da_m + 1L
      w_da      <- 1 / (da_m + t0_da)
      H_bar     <- (1 - w_da) * H_bar + w_da * (target_accept - accept_stat)
      log_eps_m <- mu_da - (sqrt(da_m) / gamma_da) * H_bar
      step_size <- exp(log_eps_m)
      m_kappa   <- da_m^(-kappa_da)
      eps_bar   <- exp(m_kappa * log_eps_m + (1 - m_kappa) * log(eps_bar))
    }

    # --- Mass matrix adaptation at 70% of warmup (in sampler-state space) ---
    # Disabled when a dense metric is pre-supplied (use_dense): the caller's
    # metric is fixed and must not be overwritten by any adaptation.
    # On update we re-find an initial step size for the NEW metric and RESET the
    # dual averager (mu/eps_bar/H_bar/da_m) so it does not mix acceptance stats
    # across two mass geometries; it then re-converges over the rest of warmup.
    #
    # metric = "warmup_dense": at this same 70% checkpoint, compute the full
    # sample covariance of the window draws, apply Ledoit-Wolf analytic
    # shrinkage for stability (handles n_window < d), set M_inv/chol_M, and
    # re-find step size under the new dense metric.  The dense path is then
    # FROZEN for the remainder of warmup and the sampling phase (same ergodicity
    # guarantee as the diagonal adaptation: freeze-at-warmup-end is valid).
    # Falls back to the diagonal path when the window is too small.
    if (!use_dense && adapt_mass && i == floor(n_warmup * 0.7) && i > 50) {
      idx <- max(1, floor(n_warmup * 0.2)):i
      if (identical(metric, "warmup_dense")) {
        ## Dense path: Ledoit-Wolf shrunken covariance -> M_inv
        win_draws <- state_chain[idx, , drop = FALSE]
        dense_res <- .warmup_dense_metric(win_draws, verbose = verbose)
        if (!is.null(dense_res)) {
          M_inv      <- dense_res$M_inv
          chol_M     <- dense_res$chol_M
          M_diag     <- rep(1, d)     ## sentinel (hot path uses M_inv, not M_diag)
          M_inv_diag <- NULL          ## signals: ignore diagonal path
          use_dense  <- TRUE          ## freeze: skip further diagonal adaptation
          eps0_new   <- .hmc_find_stepsize(state, lp_scalar, grad,
                                           M_inv_diag, M_diag,
                                           M_inv = M_inv, chol_M = chol_M,
                                           vg_fn = vg_fn, lp0 = lp_curr,
                                           g0 = g_curr)
          mu_da     <- log(10 * eps0_new)
          eps_bar   <- 1
          H_bar     <- 0
          da_m      <- 0L
          step_size <- eps0_new
          if (verbose) .dynhr_inform(sprintf(
            "HMC: warmup_dense mass set (lambda=%.3f), reset step_size = %.4e",
            dense_res$lambda, eps0_new))
        } else {
          ## Fallback to diagonal when window is too small
          vars <- apply(state_chain[idx, , drop = FALSE], 2, var)
          vars[vars < 1e-12 | !is.finite(vars)] <- 1
          ## Stan's rule: the INVERSE mass is the posterior variance. This used
          ## to set the MASS to the variance (inverted), squaring the problem's
          ## conditioning instead of removing it (brief 23 W14 finding).
          M_inv_diag <- vars
          M_diag     <- 1 / vars
          eps0_new   <- .hmc_find_stepsize(state, lp_scalar, grad,
                                           M_inv_diag, M_diag,
                                           M_inv = NULL, chol_M = NULL,
                                           vg_fn = vg_fn, lp0 = lp_curr,
                                           g0 = g_curr)
          mu_da     <- log(10 * eps0_new)
          eps_bar   <- 1
          H_bar     <- 0
          da_m      <- 0L
          step_size <- eps0_new
          if (verbose) .dynhr_inform(sprintf(
            "HMC: warmup_dense fallback to diagonal, reset step_size = %.4e", eps0_new))
        }
      } else {
        ## Default diagonal adaptation
        vars <- apply(state_chain[idx, , drop = FALSE], 2, var)
        vars[vars < 1e-12 | !is.finite(vars)] <- 1
        ## Stan's rule: the INVERSE mass is the posterior variance. This used
        ## to set the MASS to the variance (inverted), squaring the problem's
        ## conditioning instead of removing it (brief 23 W14 finding).
        M_inv_diag <- vars
        M_diag     <- 1 / vars
        eps0_new   <- .hmc_find_stepsize(state, lp_scalar, grad,
                                          M_inv_diag, M_diag,
                                          M_inv = NULL, chol_M = NULL,
                                          vg_fn = vg_fn, lp0 = lp_curr,
                                          g0 = g_curr)
        mu_da     <- log(10 * eps0_new)
        eps_bar   <- 1
        H_bar     <- 0
        da_m      <- 0L
        step_size <- eps0_new
        if (verbose) .dynhr_inform(sprintf("HMC: adapted mass matrix, reset step_size = %.4e", eps0_new))
      }
    }

    # --- Freeze step size at end of warmup (use the dual-averaged value) ---
    if (i == n_warmup) {
      step_size <- eps_bar
      if (verbose) .dynhr_inform(sprintf("HMC: warmup complete, final step_size = %.4e", step_size))
    }

    # --- Progress ---
    if (i %% 500 == 0 || i == n_total) {
      rate    <- n_accept / (i - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / i * (n_total - i)
      ch_lab  <- if (is.null(chain_id)) "?" else as.character(chain_id)
      msg <- sprintf("HMC Ch%s %d/%d accept=%.0f%% lp=%.1f eps=%.3e ETA=%.0fs",
                      ch_lab, i, n_total, rate * 100, trace_lp_curr, step_size, eta)
      if (!is.null(progressor)) {
        progressor(message = msg, amount = 1)
      } else if (verbose) {
        .dynhr_inform(msg)
      }
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  post_chain   <- chain[(n_warmup + 1):n_total, , drop = FALSE]
  post_logpost <- logpost_trace[(n_warmup + 1):n_total]

  list(
    chain           = post_chain,
    full_chain      = chain,
    logpost_trace   = logpost_trace,
    post_logpost    = post_logpost,
    acceptance_rate = n_accept / (n_total - 1),
    step_size       = step_size,
    L               = L,
    mass_matrix     = if (use_dense && !is.null(M_inv)) M_inv else M_diag,
    n_draws         = as.integer(n_draws),
    n_burn          = as.integer(n_warmup),
    n_grad_evals    = n_grad_evals,
    elapsed_secs    = elapsed,
    sampler         = "hmc",
    kernel_stats    = if (is.null(grad_fn_kernel_stats)) NULL
                       else as.list(grad_fn_kernel_stats)
  )
}




hmc_summary <- function(result, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  cat(sprintf("\n=== %s Summary ===\n", toupper(result$sampler)))
  cat(sprintf("  Draws: %d post-warmup (%d warmup)\n", result$n_draws, result$n_burn))
  cat(sprintf("  Acceptance rate: %.1f%%\n", result$acceptance_rate * 100))
  cat(sprintf("  Step size: %.4e\n", result$step_size))
  cat(sprintf("  Wall time: %.1f sec (%.1f ms/draw)\n",
              result$elapsed_secs,
              result$elapsed_secs / result$n_draws * 1000))
  cat(sprintf("  Gradient evaluations: %d\n", result$n_grad_evals))

  if (result$sampler == "nuts") {
    cat(sprintf("  Mean tree depth: %.1f\n", result$mean_treedepth))
    cat(sprintf("  Divergent transitions: %d (%.1f%%)\n",
                result$n_divergent,
                result$n_divergent / result$n_draws * 100))
    if (result$n_divergent > 0) {
      cat("  WARNING: Divergences detected. Consider reparameterisation or smaller step_size.\n")
    }
  }

  # Parameter summary table
  ch <- result$chain
  qmat <- t(apply(ch, 2, quantile, probs = probs))
  means <- colMeans(ch)
  sds   <- apply(ch, 2, sd)

  out <- data.frame(
    mean = round(means, 4),
    sd   = round(sds, 4),
    round(qmat, 4),
    check.names = FALSE
  )

  cat("\nParameter estimates:\n")
  print(out)
  invisible(out)
}
