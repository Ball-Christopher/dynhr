## R/sampler-metric.R
## --------------------------------------------------------------------------
## Metric utilities for dense-metric HMC/NUTS (Stage 0, manifold-MCMC roadmap).
##
## softabs_metric(H, alpha): SoftAbs regularisation of a symmetric matrix H.
##   λ ↦ λ·coth(αλ)  (→|λ| as α→∞; floors near-zero λ at ≈1/α).
##   Returns list(G, G_inv, L, logdet) — all symmetric PD.
## --------------------------------------------------------------------------


#' Stable scalar SoftAbs map: λ · coth(α λ)
#'
#' Handles three regimes without overflow:
#'   * |α λ| > 20  →  |λ|   (asymptotic; coth → sign)
#'   * |α λ| < 1e-8 →  1/α  (near-zero; Taylor: λ coth(αλ) ≈ 1/α)
#'   * otherwise    →  λ / tanh(α λ)  (exact formula)
#'
#' @param lam eigenvalue (scalar, possibly negative or near-zero)
#' @param alpha softness parameter (large → |λ|)
#' @return softabs(λ) ≥ 0 (strictly positive: floored at 1/alpha for λ≈0)
#' @noRd
.softabs_scalar <- function(lam, alpha) {
  al <- alpha * lam
  if (abs(al) > 20) {
    ## Asymptotic: coth(x) → sign(x), so λ·coth(αλ) → |λ|
    abs(lam)
  } else if (abs(al) < 1e-8) {
    ## Near-zero Taylor: λ coth(αλ) = λ/(αλ) = 1/α + O((αλ)²/3)
    1 / alpha
  } else {
    ## Standard formula; tanh is never zero here (|αλ| ≥ 1e-8)
    lam / tanh(al)
  }
}


#' SoftAbs metric from a symmetric matrix H
#'
#' Eigendecomposes H = Q Λ Qᵀ (via eigen(H, symmetric = TRUE)) and maps each
#' eigenvalue λ to the SoftAbs value λ·coth(αλ), which:
#'   - equals |λ| for large |αλ| (absolute-value metric, PD by construction);
#'   - floors at 1/α for λ ≈ 0 (regularises flat directions);
#'   - is smooth and differentiable everywhere.
#'
#' Guarantees symmetric PD output regardless of H's definiteness.
#'
#' @param H  Symmetric matrix (d × d); typically the negative-log-posterior
#'   Hessian (positive entries mean curvature).  May be indefinite or
#'   near-singular — SoftAbs handles both.
#' @param alpha  Softness parameter (default 1e6).  As alpha → ∞ the map
#'   approaches |λ| (pure absolute value).  Smaller alpha gives a smoother
#'   but less tight lower bound.
#'
#' @return A named list:
#'   \item{G}{SoftAbs-regularised PD matrix (d × d)}
#'   \item{G_inv}{Inverse of G (d × d)}
#'   \item{L}{Upper Cholesky factor of G: \code{chol(G)}, so G = t(L) \%*\% L}
#'   \item{logdet}{log|G| = sum(log(softabs eigenvalues))}
#' @noRd
softabs_metric <- function(H, alpha = 1e6) {
  stopifnot(is.matrix(H), isSymmetric(H, tol = 1e-10))
  d <- nrow(H)

  ## Eigendecompose (symmetric = TRUE for stability / correct ordering)
  ev <- eigen(H, symmetric = TRUE)
  Q  <- ev$vectors           ## d × d orthogonal
  lam <- ev$values           ## length d (descending by default)

  ## SoftAbs map (vectorised over eigenvalues via sapply for clarity)
  lam_soft <- vapply(lam, .softabs_scalar, numeric(1L), alpha = alpha)
  ## Paranoia: ensure all values are strictly positive (should hold by math,
  ## but floating-point can produce tiny negatives near zero)
  lam_soft <- pmax(lam_soft, .Machine$double.eps)

  ## Reassemble G = Q diag(lam_soft) Qᵀ
  ## Using sweep is cleaner than Q %*% diag(lam_soft) %*% t(Q) and avoids
  ## allocating a full d×d diagonal matrix.
  G_raw <- Q %*% (sweep(t(Q), 1L, lam_soft, `*`))
  ## Symmetrize to kill numerical skew (||G - Gᵀ|| / ||G|| typically ~1e-16)
  G <- (G_raw + t(G_raw)) / 2

  ## Inverse: G⁻¹ = Q diag(1/lam_soft) Qᵀ
  G_inv_raw <- Q %*% (sweep(t(Q), 1L, 1 / lam_soft, `*`))
  G_inv <- (G_inv_raw + t(G_inv_raw)) / 2

  ## Cholesky of G (should succeed because G is PD)
  L <- tryCatch(chol(G), error = function(e) {
    ## Fallback: add a tiny nugget in case of floating-point non-PD
    chol(G + diag(1e-12 * max(diag(G)), d))
  })

  list(G = G, G_inv = G_inv, L = L, logdet = sum(log(lam_soft)))
}


# ============================================================================
# Ledoit-Wolf analytic shrinkage for warmup-dense mass matrix
# ============================================================================

#' Ledoit-Wolf analytic shrinkage toward scaled identity
#'
#' Computes the Ledoit-Wolf (2004) analytic shrinkage intensity lambda that
#' minimises the expected Frobenius loss between the shrunken sample
#' covariance and the true covariance, using the closed-form Ledoit & Wolf
#' (2004) analytic formula:
#'
#'   S_shrunk = (1 - lambda) * S + lambda * mu * I
#'
#' where mu = trace(S) / d (the scaled-identity target) and
#'
#'   lambda = min(1, max(0,
#'     ( (n - 2)/n * trace(S^2) + trace(S)^2 )
#'     / ( (n + 2) * (trace(S^2) - trace(S)^2 / d) )
#'   ))
#'
#' This is the Ledoit-Wolf 2004 closed-form shrinkage intensity for a Gaussian
#' population. It is invariant to the overall scale of S and is well-defined
#' for n >= 2 and d >= 2. For n > d the estimator is consistent; for n <= d
#' (more parameters than warmup draws) it regularises aggressively (lambda
#' can reach 1, collapsing to the scaled identity).
#'
#' @param S  Sample covariance matrix (d x d, symmetric PSD)
#' @param n  Number of observations used to compute S
#' @return List with:
#'   \item{S_shrunk}{Shrunken covariance matrix ((1-lambda)*S + lambda*mu*I)}
#'   \item{lambda}{Shrinkage intensity in [0,1]}
#'   \item{mu}{Target scale (trace(S)/d)}
#' @noRd
.ledoit_wolf_shrink <- function(S, n) {
  d  <- nrow(S)
  stopifnot(d >= 2L, n >= 2L)

  trS  <- sum(diag(S))
  trS2 <- sum(S * S)      # trace(S %*% S) = sum(S^2) element-wise for symmetric S

  mu   <- trS / d

  ## Ledoit-Wolf (2004) analytic formula
  ## numerator   = ((n-2)/n) * tr(S^2) + tr(S)^2
  ## denominator = (n+2)   * (tr(S^2) - tr(S)^2/d)
  numer <- ((n - 2) / n) * trS2 + trS^2
  denom <- (n + 2) * (trS2 - trS^2 / d)

  lambda <- if (denom <= 0 || !is.finite(denom)) {
    ## Degenerate: all eigenvalues equal (S = c*I) or rank-1 -> full shrinkage
    1
  } else {
    min(1, max(0, numer / denom))
  }

  S_shrunk <- (1 - lambda) * S + lambda * mu * diag(d)
  list(S_shrunk = S_shrunk, lambda = lambda, mu = mu)
}


#' Compute warm-start dense mass matrix from a window of sampler draws
#'
#' Given a matrix of within-window draws (rows = iterations, cols = params),
#' computes the sample covariance, applies Ledoit-Wolf analytic shrinkage for
#' numerical stability, and returns the inverse mass matrix M_inv = S_shrunk
#' (the estimated covariance, following the same convention as metric="hessian"
#' where M_inv = Sigma_prop) together with chol_M = chol(solve(M_inv)) for
#' use in the NUTS/HMC dense path.
#'
#' Falls back to NULL (diagonal path) when:
#'   - fewer than max(d + 1, 20) rows (window too short for reliable estimates)
#'   - S_shrunk is not positive-definite after shrinkage (should be rare)
#'   - any dimension is zero
#'
#' @param win_draws  Matrix (n_win x d) of draws in the window
#' @param verbose    Print a message on success/fallback
#' @return List(M_inv, chol_M, S_shrunk, lambda) or NULL on fallback
#' @noRd
.warmup_dense_metric <- function(win_draws, verbose = FALSE) {
  n_win <- nrow(win_draws)
  d     <- ncol(win_draws)
  if (d < 2L || n_win < max(d + 1L, 20L)) {
    if (verbose) .dynhr_inform(sprintf(
      "warmup_dense: window too small (n=%d, d=%d); falling back to diagonal.",
      n_win, d))
    return(NULL)
  }

  ## Sample covariance (n_win rows)
  S_raw <- cov(win_draws)          # (n_win-1) divisor; fine for Ledoit-Wolf
  ## Guard: replace any NA/Inf rows/cols with 1 on diagonal
  bad <- !is.finite(diag(S_raw))
  if (any(bad)) {
    S_raw[bad, ] <- 0
    S_raw[, bad] <- 0
    S_raw[cbind(which(bad), which(bad))] <- 1
  }

  ## Ledoit-Wolf shrinkage (use n_win - 1 = cov() denominator's effective n)
  lw    <- .ledoit_wolf_shrink(S_raw, n = n_win - 1L)
  S_shr <- lw$S_shrunk

  ## M_inv = S_shrunk (the estimated target covariance = inverse mass matrix).
  ## This follows the same convention as metric="hessian" where M_inv = Sigma_prop
  ## (the posterior covariance from mode-finding): M_inv encodes the TARGET
  ## COVARIANCE, so the mass M = solve(M_inv) = precision matrix.
  ## Momentum is drawn from N(0, M) = N(0, precision), which correctly
  ## whitens the geometry: Var(position) ~ Var(theta) ~ M_inv.
  ##
  ## chol_M = chol(M) = chol(solve(S_shrunk)) (upper triangular) is needed
  ## by .hmc_sample_momentum: r = t(chol_M) %*% z, Cov(r) = M.
  M_inv <- S_shr   ## M_inv = covariance estimate (convention: M_inv = Sigma)

  ## Symmetrise M_inv to kill floating-point skew
  M_inv <- (M_inv + t(M_inv)) / 2

  ## Verify M_inv is PD (should hold after Ledoit-Wolf; guard for edge cases)
  ev_inv <- tryCatch(eigen(M_inv, only.values = TRUE)$values, error = function(e) NULL)
  if (is.null(ev_inv) || any(ev_inv <= 0)) {
    if (verbose) .dynhr_inform("warmup_dense: S_shrunk not PD after shrinkage; falling back to diagonal.")
    return(NULL)
  }

  ## M = solve(M_inv) = precision matrix; chol_M = chol(M)
  M_for_chol <- tryCatch(solve(M_inv), error = function(e) NULL)
  if (is.null(M_for_chol) || !all(is.finite(M_for_chol))) {
    if (verbose) .dynhr_inform("warmup_dense: S_shrunk not invertible; falling back to diagonal.")
    return(NULL)
  }
  M_for_chol <- (M_for_chol + t(M_for_chol)) / 2   ## symmetrise before chol
  chol_M <- tryCatch(chol(M_for_chol), error = function(e) NULL)

  if (is.null(chol_M)) {
    if (verbose) .dynhr_inform("warmup_dense: M=solve(S_shrunk) not PD; falling back to diagonal.")
    return(NULL)
  }

  if (verbose) .dynhr_inform(sprintf(
    "warmup_dense: lambda=%.3f mu=%.4g n=%d d=%d -> dense mass matrix set.",
    lw$lambda, lw$mu, n_win, d))

  list(M_inv   = M_inv,
       chol_M  = chol_M,
       S_shrunk = S_shr,
       lambda  = lw$lambda,
       mu      = lw$mu)
}


# ============================================================================
# Low-rank-plus-diagonal inverse mass matrix (NUTS metric = "lowrank")
# ============================================================================
#
# Reference: J. Lao (2026), "The Universal Warmup Path: Automatic
# Preconditioner Selection for HMC", arXiv:2607.23788.  Its controller starts
# diagonal and, at dimension-derived window endpoints, promotes to a
# low-rank-plus-diagonal inverse mass matrix when the window evidence supports
# it (rank capacity k_cap = min(50, max(floor(d/2), 1)), pooled support floor
# N_min = 8 (k_cap + 1), first window ceil(N_min / M) per chain, 1.5x window
# growth, final 15% step size only).  The paper does not restate the
# low-rank-plus-diagonal estimator; it uses the Fisher-divergence (draws +
# scores) form of A. Seyboldt, E. L. Carlson & B. Carpenter (2026),
# "Preconditioning Hamiltonian Monte Carlo by minimizing Fisher Divergence",
# arXiv:2603.18845 -- .lowrank_estimate() below is their Algorithm 1
# (diagonal scale, joint span of whitened draws + scores, SPD geometric mean
# `spdm`, eigen-filter), with Lao's bulk-relative BBP-edge filter in place of
# their fixed `lambda <= 1/c or >= c` cut.  Their diagonal case (Theorem 2.2,
# eq. 6) is .fisher_diag_estimate() / NUTS metric = "fisher_diag".
#
# Convention (Lao eq. 1, and the dense path of this package):
#   G = M^{-1} = diag(sd) (I + U (Lambda - I) U') diag(sd),
#   p ~ N(0, M),  K(p) = p' G p / 2,  velocity = G p,
# with U (d x k) orthonormal columns and Lambda = diag(lam), lam > 0.  Since
# U' U = I the inverse is available in closed form,
#   M = diag(1/sd) (I + U (Lambda^{-1} - I) U') diag(1/sd),
# and A = I + U (Lambda^{-1/2} - I) U' satisfies A A = I + U (Lambda^{-1} - I) U',
# so p = diag(1/sd) A z, z ~ N(0, I), has covariance M.  Every operation is
# O(d k): the d x d matrices are never formed in the sampler hot path.
# ============================================================================


#' Construct a low-rank-plus-diagonal inverse-mass metric object
#'
#' @param sd  Positive numeric vector (length d): diagonal scale.
#' @param U   d x k matrix with orthonormal columns (k may be 0).
#' @param lam Positive numeric vector (length k): eigenvalues of the low-rank
#'   correction in the \code{sd}-whitened frame.
#' @return Object of class \code{"dynhr_lowrank_metric"}.
#' @noRd
.lowrank_metric <- function(sd, U = NULL, lam = numeric(0)) {
  sd <- as.numeric(sd)
  d  <- length(sd)
  if (is.null(U)) U <- matrix(0, d, 0L)
  U   <- as.matrix(U)
  lam <- as.numeric(lam)
  stopifnot(all(is.finite(sd)), all(sd > 0), nrow(U) == d,
            ncol(U) == length(lam), all(is.finite(lam)), all(lam > 0))
  structure(list(sd = sd, U = U, lam = lam, rank = length(lam)),
            class = "dynhr_lowrank_metric")
}


#' Apply the low-rank-plus-diagonal inverse mass: G r (velocity)
#' @noRd
.lowrank_apply_inv <- function(lr, r) {
  q <- lr$sd * as.numeric(r)
  if (lr$rank > 0L) {
    q <- q + as.numeric(lr$U %*% ((lr$lam - 1) * as.numeric(crossprod(lr$U, q))))
  }
  lr$sd * q
}


#' Draw momentum p ~ N(0, M) for a low-rank-plus-diagonal metric
#' @noRd
.lowrank_sample_momentum <- function(lr) {
  z <- rnorm(length(lr$sd))
  if (lr$rank > 0L) {
    z <- z + as.numeric(lr$U %*% ((lr$lam^(-0.5) - 1) * as.numeric(crossprod(lr$U, z))))
  }
  z / lr$sd
}


#' Explicit dense inverse mass G = M^{-1} (for reporting / tests only)
#' @noRd
.lowrank_dense_inv <- function(lr) {
  d <- length(lr$sd)
  core <- diag(d)
  if (lr$rank > 0L) core <- core + lr$U %*% ((lr$lam - 1) * t(lr$U))
  G <- core * outer(lr$sd, lr$sd)
  (G + t(G)) / 2
}


#' Explicit dense mass M (for reporting / tests only)
#' @noRd
.lowrank_dense_mass <- function(lr) {
  d <- length(lr$sd)
  core <- diag(d)
  if (lr$rank > 0L) core <- core + lr$U %*% ((1 / lr$lam - 1) * t(lr$U))
  M <- core / outer(lr$sd, lr$sd)
  (M + t(M)) / 2
}


#' Inverse-mass times momentum for any non-diagonal metric
#'
#' Dense matrix: \code{as.numeric(M_inv \%*\% r)} (unchanged expression, so the
#' dense path stays bit-identical).  Low-rank object: O(d k) product.
#' @noRd
.metric_apply_inv <- function(M_inv, r) {
  if (inherits(M_inv, "dynhr_lowrank_metric")) {
    .lowrank_apply_inv(M_inv, r)
  } else {
    as.numeric(M_inv %*% r)
  }
}


#' Rank capacity of the low-rank metric (Lao 2026, Sec. 3)
#' @noRd
.lowrank_rank_cap <- function(d) {
  as.integer(min(50L, max(floor(d / 2), 1L)))
}


#' Warmup schedule for metric = "lowrank" (Lao 2026, Sec. 3 / Fig. 1)
#'
#' One-step initialisation window; first mass-matrix window
#' \code{n1 = ceiling(N_min / n_chains)} with \code{N_min = 8 (k_cap + 1)};
#' nominal 1.5x growth; when the next grown window would not fit before the
#' step-size-only phase the last slow window absorbs the remainder; the final
#' 15\% of warmup adapts the step size only.  Same \code{data.frame(start,
#' end, type)} format as \code{.nuts_warmup_windows()}.
#'
#' @param n_warmup Warmup iterations (per chain).
#' @param d Dimension.
#' @param n_chains Chains pooled at each endpoint (1 for a single NUTS chain).
#' @noRd
.lowrank_warmup_windows <- function(n_warmup, d, n_chains = 1L) {
  n_warmup <- as.integer(n_warmup)
  if (n_warmup < 10L) {
    return(data.frame(start = 1L, end = n_warmup, type = "fast_init",
                      stringsAsFactors = FALSE))
  }
  k_cap <- .lowrank_rank_cap(d)
  n1    <- as.integer(ceiling(8 * (k_cap + 1) / n_chains))
  term_buffer <- as.integer(max(1L, floor(0.15 * n_warmup)))
  slow_start  <- 2L
  slow_end    <- n_warmup - term_buffer

  rows <- list(data.frame(start = 1L, end = 1L, type = "fast_init",
                          stringsAsFactors = FALSE))
  cur <- slow_start
  w   <- n1
  while (cur <= slow_end) {
    win_end <- cur + w - 1L
    next_w  <- as.integer(ceiling(1.5 * w))
    if (win_end >= slow_end || win_end + next_w > slow_end) win_end <- slow_end
    rows[[length(rows) + 1L]] <- data.frame(start = cur, end = win_end,
                                            type = "slow",
                                            stringsAsFactors = FALSE)
    cur <- win_end + 1L
    w   <- next_w
  }
  rows[[length(rows) + 1L]] <- data.frame(start = slow_end + 1L, end = n_warmup,
                                          type = "fast_final",
                                          stringsAsFactors = FALSE)
  do.call(rbind, rows)
}


#' SPD geometric mean A # B = A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}
#' @noRd
.spd_geometric_mean <- function(A, B) {
  ea <- eigen(A, symmetric = TRUE)
  va <- pmax(ea$values, .Machine$double.eps)
  A_h  <- ea$vectors %*% (sqrt(va) * t(ea$vectors))
  A_mh <- ea$vectors %*% ((1 / sqrt(va)) * t(ea$vectors))
  inner <- A_mh %*% B %*% A_mh
  inner <- (inner + t(inner)) / 2
  ei <- eigen(inner, symmetric = TRUE)
  inner_h <- ei$vectors %*% (sqrt(pmax(ei$values, 0)) * t(ei$vectors))
  G <- A_h %*% inner_h %*% A_h
  (G + t(G)) / 2
}


#' Fisher-divergence diagonal inverse mass from one warmup window
#'
#' Seyboldt, Carlson & Carpenter (2026, arXiv:2603.18845), Theorem 2.2 /
#' eq. (6): the diagonal affine map minimising the sample Fisher divergence to
#' N(0, I) has inverse mass \code{M^{-1} = diag(sigma^2)} with
#' \code{sigma_i^2 = sqrt(var(x_i) / var(g_i))} -- the geometric mean of the
#' draw variance and the inverse score variance.  For a Gaussian coordinate
#' with independent components both factors equal the variance, so it agrees
#' with the variance rule; on correlated targets it sits between the marginal
#' variance and the conditional variance \code{1 / var(g_i)}, penalising large
#' and small eigenvalues of the preconditioned covariance symmetrically.
#' Coordinates whose draw or score variance is degenerate (< 1e-12 or
#' non-finite) fall back to the draw variance, or 1 if that too is degenerate
#' (same guard as \code{.lowrank_estimate()}).
#'
#' @param draws N x d matrix of warmup states (sampler space), N >= 2.
#' @param grads N x d matrix of log-density gradients at those states.
#' @return Numeric length-d vector: the inverse-mass diagonal \code{sigma^2}.
#' @noRd
.fisher_diag_estimate <- function(draws, grads) {
  draws <- as.matrix(draws)
  grads <- as.matrix(grads)
  stopifnot(nrow(draws) >= 2L, identical(dim(draws), dim(grads)))
  var_x <- apply(draws, 2L, stats::var)
  var_g <- apply(grads, 2L, stats::var)
  ok_x  <- is.finite(var_x) & var_x >= 1e-12
  ok_g  <- is.finite(var_g) & var_g >= 1e-12
  sd2   <- sqrt(var_x / var_g)
  bad   <- !(ok_x & ok_g) | !is.finite(sd2)
  sd2[bad] <- ifelse(ok_x[bad], var_x[bad], 1)
  sd2
}


#' Estimate a low-rank-plus-diagonal inverse mass from one warmup window
#'
#' Fisher-divergence estimator (draws + scores; Seyboldt, Carlson & Carpenter
#' 2026, arXiv:2603.18845, Algorithm 1), within-chain version:
#' \enumerate{
#'   \item per-chain-centre draws x and scores g (Lao 2026 eq. 5, W);
#'   \item diagonal scale \code{sd_i^2 = sqrt(var(x_i) / var(g_i))};
#'   \item whiten \code{x~ = x / sd}, \code{g~ = g * sd}; restrict to the span
#'     Q of the whitened draws and scores (the full space when N - M >= d);
#'   \item \code{S = Cov(x~) # Cov(g~)^{-1}} (SPD geometric mean; equals the
#'     whitened covariance for a Gaussian target);
#'   \item eigen-decompose S; relative to the bulk level \code{c0} (geometric
#'     median of the eigenvalues) retain the directions with
#'     \code{|log(lambda / c0)| > log(tau)},
#'     \code{tau = max(cutoff, (1 + sqrt(d / N))^2)} (the BBP noise edge as a
#'     scale calibration, Lao Sec. 5.3), largest \code{|log|} first, at most
#'     \code{min(k_cap, floor(N / 8) - 1)} of them (support floor
#'     \code{N >= 8 (k + 1)});
#'   \item the non-retained directions share the geometric-mean bulk level
#'     \code{c}: \code{G = c diag(sd) (I + U (Lambda / c - I) U') diag(sd)}.
#' }
#' Rank 0 (promotion not supported) is the Fisher diagonal metric.
#'
#' @param draws N x d matrix of warmup states (sampler space).
#' @param grads N x d matrix of log-density gradients at those states.
#' @param chain Optional length-N chain labels (per-chain centring); NULL =
#'   one chain.
#' @param cutoff Minimum eigenvalue ratio for promotion (default 2).
#' @param max_rank Optional user rank cap (NULL = k_cap(d)).
#' @param gamma Ridge added to the projected covariances (default 1e-5).
#' @return list(metric, rank, eigenvalues, threshold, bulk, k_max, n)
#' @noRd
.lowrank_estimate <- function(draws, grads, chain = NULL, cutoff = 2,
                              max_rank = NULL, gamma = 1e-5) {
  draws <- as.matrix(draws)
  grads <- as.matrix(grads)
  N <- nrow(draws)
  d <- ncol(draws)
  stopifnot(nrow(grads) == N, ncol(grads) == d, N >= 3L, cutoff > 1)
  if (is.null(chain)) chain <- rep(1L, N)
  chain <- as.integer(factor(chain))
  n_ch  <- max(chain)

  xc <- draws
  gc <- grads
  for (cc in seq_len(n_ch)) {
    rows <- which(chain == cc)
    xc[rows, ] <- sweep(draws[rows, , drop = FALSE], 2L,
                        colMeans(draws[rows, , drop = FALSE]))
    gc[rows, ] <- sweep(grads[rows, , drop = FALSE], 2L,
                        colMeans(grads[rows, , drop = FALSE]))
  }
  dof <- N - n_ch

  var_x <- colSums(xc^2) / dof
  var_g <- colSums(gc^2) / dof
  ok_x  <- is.finite(var_x) & var_x >= 1e-12
  ok_g  <- is.finite(var_g) & var_g >= 1e-12
  sd2   <- sqrt(var_x / var_g)
  bad   <- !(ok_x & ok_g) | !is.finite(sd2)
  sd2[bad] <- ifelse(ok_x[bad], var_x[bad], 1)
  sd <- sqrt(sd2)

  xt <- sweep(xc, 2L, sd, "/")
  gt <- sweep(gc, 2L, sd, "*")
  xt[!is.finite(xt)] <- 0
  gt[!is.finite(gt)] <- 0

  ## Span of the whitened draws and scores (full space when well supported)
  if (dof >= d) {
    Q <- diag(d)
  } else {
    qq <- qr(cbind(t(xt), t(gt)))
    Q  <- qr.Q(qq)[, seq_len(qq$rank), drop = FALSE]
  }
  r_sub <- ncol(Q)
  xq <- xt %*% Q
  gq <- gt %*% Q
  Cx <- crossprod(xq) / dof + gamma * diag(r_sub)
  Cg <- crossprod(gq) / dof + gamma * diag(r_sub)
  S  <- .spd_geometric_mean(Cx, solve(Cg))
  es <- eigen(S, symmetric = TRUE)
  lam_all <- pmax(es$values, .Machine$double.eps)

  k_cap <- if (is.null(max_rank)) .lowrank_rank_cap(d) else
    as.integer(min(max_rank, .lowrank_rank_cap(d)))
  k_max <- as.integer(max(0L, min(k_cap, floor(N / 8) - 1L, r_sub - 1L)))
  tau   <- max(cutoff, (1 + sqrt(d / N))^2)
  c0    <- exp(stats::median(log(lam_all)))
  z     <- abs(log(lam_all / c0))
  cand  <- which(z > log(tau))
  keep  <- cand[order(-z[cand])]
  keep  <- sort(keep[seq_len(min(length(keep), k_max))])
  bulk  <- if (length(keep) > 0L) exp(mean(log(lam_all[-keep]))) else
    exp(mean(log(lam_all)))

  U   <- Q %*% es$vectors[, keep, drop = FALSE]
  lam <- lam_all[keep] / bulk
  list(metric      = .lowrank_metric(sd * sqrt(bulk), U, lam),
       rank        = length(keep),
       eigenvalues = lam_all / bulk,
       threshold   = tau,
       bulk        = bulk,
       k_max       = k_max,
       n           = N)
}
