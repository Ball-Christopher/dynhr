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
    if (verbose) message(sprintf(
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
    if (verbose) message("warmup_dense: S_shrunk not PD after shrinkage; falling back to diagonal.")
    return(NULL)
  }

  ## M = solve(M_inv) = precision matrix; chol_M = chol(M)
  M_for_chol <- tryCatch(solve(M_inv), error = function(e) NULL)
  if (is.null(M_for_chol) || !all(is.finite(M_for_chol))) {
    if (verbose) message("warmup_dense: S_shrunk not invertible; falling back to diagonal.")
    return(NULL)
  }
  M_for_chol <- (M_for_chol + t(M_for_chol)) / 2   ## symmetrise before chol
  chol_M <- tryCatch(chol(M_for_chol), error = function(e) NULL)

  if (is.null(chol_M)) {
    if (verbose) message("warmup_dense: M=solve(S_shrunk) not PD; falling back to diagonal.")
    return(NULL)
  }

  if (verbose) message(sprintf(
    "warmup_dense: lambda=%.3f mu=%.4g n=%d d=%d -> dense mass matrix set.",
    lw$lambda, lw$mu, n_win, d))

  list(M_inv   = M_inv,
       chol_M  = chol_M,
       S_shrunk = S_shr,
       lambda  = lw$lambda,
       mu      = lw$mu)
}
