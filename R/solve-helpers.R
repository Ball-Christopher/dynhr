## R/solve-helpers.R
## --------------------------------------------------------------------------
## Small numerical utilities used by the perturbation solver: safe SVD-based
## pseudoinverse, etc.
##
## Phase-1 split from perturbation-monolith.R (no logic changes).
## --------------------------------------------------------------------------

# ---- Utility: safe matrix inverse via SVD ----

#' SVD-based pseudoinverse with a RELATIVE singular-value cutoff
#'
#' @param M     Square or rectangular matrix.
#' @param rtol  RELATIVE singular-value cutoff: singular values satisfying
#'   \code{d <= rtol * max(d)} are treated as exact zeros.  Default
#'   \code{1e-12}.  A relative cutoff is scale-invariant --
#'   \code{.safe_inv(c * M) == .safe_inv(M) / c} for any finite \code{c != 0}
#'   -- which the previous ABSOLUTE cutoff (\code{d <= 1e-10}) was not: a
#'   perfectly invertible but badly scaled matrix (e.g. equations divided
#'   through by \code{1e-7}, so every singular value sits below \code{1e-10})
#'   had ALL of its singular values zeroed, silently handing the caller a zero
#'   "inverse".  On the QZ decision-rule path that produced a wrong decision
#'   rule still flagged \code{bk_satisfied = TRUE}.
#' @param tol   Optional ABSOLUTE cutoff escape hatch.  When non-\code{NULL} it
#'   overrides \code{rtol} and reproduces the historical behaviour
#'   (\code{tol = 1e-10} is the pre-2026-09 default).  No in-package caller
#'   uses it; it exists so a caller that genuinely knows the absolute noise
#'   floor of \code{M} can say so.
#' @param warn_label Optional character label.  When non-\code{NULL} and at
#'   least one singular value is truncated, ONE \code{warning()} per call is
#'   emitted naming the label and the number of truncated singular values.
#'   Used on the decision-rule call sites, where a truncation means the
#'   returned matrix is a PSEUDO-inverse and the resulting rule is not the
#'   unique solution of the linear system.
#' @return Pseudoinverse of M.
#' @noRd
.safe_inv <- function(M, rtol = 1e-12, tol = NULL, warn_label = NULL) {
  s <- svd(M)
  d <- s$d
  cutoff <- if (!is.null(tol)) {
    tol
  } else {
    d_max <- if (length(d)) max(d) else 0
    if (is.finite(d_max) && d_max > 0) rtol * d_max else 0
  }
  ## Numerically identical to ifelse(s$d > cutoff, 1/s$d, 0) but without the
  ## ifelse overhead; and scale rows of t(u) by d_inv directly instead of
  ## forming a diagonal matrix (diag(d) %*% X == d * X), saving an allocation
  ## and a matrix multiply.
  dropped <- d <= cutoff
  d_inv <- 1 / d
  d_inv[dropped] <- 0
  if (!is.null(warn_label) && any(dropped)) {
    n_drop <- sum(dropped)
    warning(sprintf(
      paste0("%s: .safe_inv() truncated %d of %d singular value%s at the ",
             "relative cutoff %.3g (max singular value %.3g); the result is a ",
             "PSEUDO-inverse, so the decision rule is NOT the unique solution ",
             "of the linear system."),
      warn_label, n_drop, length(d), if (n_drop == 1L) "" else "s",
      cutoff, if (length(d)) max(d) else NA_real_), call. = FALSE)
  }
  s$v %*% (d_inv * t(s$u))
}

# ---- The single discrete-Lyapunov solver ----------------------------------
#
# Consolidation note (D2, 2026-09-02). dynhr used to carry FOUR discrete
# Lyapunov solvers:
#   * solve_lyapunov()               R/stochsimul-monolith.R  (doubling, the
#                                    correct one -- relative tol, stability
#                                    gate, NaN contract, kron fallback)
#   * .solve_lyapunov()              R/backend-monolith.R     (direct kron only)
#   * a pskf copy and a tpf copy     (removed in wave 1)
# They disagreed on near-unit-root systems: the direct-kron version's rcond
# gate fires on highly non-normal but perfectly stable transition matrices
# (rcond(I - A (x) A) underflows machine eps while the Lyapunov equation
# itself is well posed), returning NaN where the doubling iteration converges
# fine; and its n == 1 shortcut Q / (1 - a^2) returned a NEGATIVE "variance"
# for |a| > 1 instead of NaN.
#
# There is now ONE implementation (below) and `.solve_lyapunov()` is a thin
# alias of it that preserves the NaN-on-nonstationary contract every caller
# already checks with `all(is.finite(.))`.
#' Solve the discrete Lyapunov equation X = A X A' + B
#'
#' Uses the doubling algorithm for efficiency.
#'
#' @param A Square matrix
#' @param B Symmetric positive semi-definite matrix
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Solution matrix X
#' @export
solve_lyapunov <- function(A, B, max_iter = 500L, tol = 1e-14) {
  n <- nrow(A)

  ## Fast path: if B is all zeros, solution is zero
  if (all(B == 0)) return(matrix(0, n, n))

  ## Doubling algorithm
  X <- B
  A_pow <- A
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    X_new <- X + A_pow %*% X %*% t(A_pow)
    if (any(!is.finite(X_new))) break
    diff <- max(abs(X_new - X))
    if (!is.finite(diff)) break
    ## RELATIVE convergence: a near-unit root gives a huge stationary covariance
    ## (entries ~ 1/(1-rho^2)), so the per-step increment can never fall below an
    ## ABSOLUTE 1e-14 (it plateaus at ~max|X| * machine-eps). An absolute test
    ## therefore never converges for near-unit-root systems -> the loop runs all
    ## max_iter steps and falls through to the O(n^6) kronecker solve (~0.5 s for
    ## n = 37). Scaling by max|X| makes it converge in the proper ~log2(mixing)
    ## steps for ANY stable A. (Harmless for well-damped A where max|X| ~ O(1).)
    if (diff < tol * max(1, max(abs(X_new)))) { converged <- TRUE; break }
    A_pow <- A_pow %*% A_pow
    if (any(!is.finite(A_pow))) break
    X <- X_new
  }
  if (converged) return(X)

  ## Stability gate before the O(n^6) vec/kronecker fallback.
  ##
  ## The discrete Lyapunov equation X = A X A' + B has a finite (PSD) solution
  ## only when A is stable (spectral radius < 1). The doubling loop above only
  ## fails to converge when A has a unit/explosive root: A^k does not decay, so
  ## the stationary covariance diverges and NO valid X exists. In that case the
  ## kronecker fallback is doubly bad -- it spends O(n^6) solving an (n^2 x n^2)
  ## system (e.g. ~0.5 s for n = 37) AND, because M = I - A (x) A is only
  ## *near*-singular for an explosive root (rcond just above machine-eps, so the
  ## guard below misses it), it returns a garbage non-PSD X instead of NaN.
  ##
  ## A cheap eigenvalue check (~0.2 ms for n = 37) short-circuits this: signal
  ## non-stationarity with NaN (the same contract as the singular-M guard
  ## below), letting the caller fall back to the exact-diffuse Kalman init /
  ## simulation-based moments. This is the dominant hot path when an estimation
  ## sampler probes BK-boundary draws whose decision rule is near-explosive
  ## (verified on NZSIM: ~473 ms/eval of wasted kronecker solves).
  if (max(Mod(eigen(A, only.values = TRUE)$values)) >= 1)
    return(matrix(NaN, n, n))

  ## Fallback: vec method  vec(X) = (I - A (x) A)^{-1} vec(B)
  I_n2 <- diag(n^2)
  AkA <- kronecker(A, A)
  M <- I_n2 - AkA
  # Check for singular system (unit roots from NN1 placeholders, etc.)
  if (rcond(M) < .Machine$double.eps) {
    # Return a matrix with NaN to signal non-stationarity; caller can
    # fall back to simulation-based welfare.
    # Unit root detected (e.g. NN1 placeholder equations). Not an error;
    # the caller can fall back to simulation-based computations.
    if (getOption("dynhr.warn_lyapunov", FALSE)) {
      warning("solve_lyapunov: system is singular (unit root detected). Returning NaN.")
    }
    return(matrix(NaN, n, n))
  }
  x_vec <- solve(M, as.vector(B))
  matrix(x_vec, nrow = n, ncol = n)
}

#' Internal alias of \code{\link{solve_lyapunov}}
#'
#' Historical name used by the gradient / Hessian / backend code paths.  It
#' used to be a separate direct-kronecker implementation in
#' \code{R/backend-monolith.R}; it is now a one-line alias so that there is a
#' SINGLE Lyapunov solver in the package.
#'
#' Contract (unchanged, and now strictly stronger): returns an \code{n x n}
#' matrix of \code{NaN} when no stationary solution exists (spectral radius of
#' \code{A} >= 1, or a singular \code{I - A (x) A}).  Every caller must test
#' the result with \code{all(is.finite(.))} -- a \code{tryCatch(..., error =)}
#' does NOT catch this, because no condition is signalled.
#'
#' @param A Square matrix.
#' @param Q Right-hand side (need not be positive semi-definite: the doubling
#'   iteration converges for any symmetric \code{Q} when \code{A} is stable,
#'   which the adjoint code paths rely on).
#' @param tol Accepted for backward compatibility and IGNORED (the previous
#'   direct-kronecker implementation ignored it too); the doubling iteration
#'   uses \code{solve_lyapunov()}'s own relative tolerance.
#' @return Solution matrix, or an all-\code{NaN} matrix (see Contract).
#' @noRd
.solve_lyapunov <- function(A, Q, tol = .LYAP_TOL) solve_lyapunov(A, Q)


# ---- Utility: damped-Newton (Armijo) line search ----

#' Backtracking Armijo line search for a stacked-path Newton step
#'
#' Shared by the perfect-foresight solver (\code{R/perfect-foresight-solve.R})
#' and the MCP path solver (\code{R/mcp-solve.R}); the two used to carry
#' byte-identical private copies (\code{.pf_line_search} /
#' \code{.mcp_line_search}).
#'
#' @param Y           T x n_endo path matrix (incumbent)
#' @param delta_vec   Numeric (T*n_endo): Newton direction
#' @param R_stack     Numeric (T*n_endo): current residual
#' @param J_stack     dgCMatrix: current Jacobian
#' @param theta_cur   Scalar: current merit value = 1/2 ||R||^2
#' @param fn_merit    Function to evaluate theta at a trial Y
#' @param sigma       Armijo parameter (default 1e-4)
#' @param max_ls      Maximum line-search iterations (default 20)
#' @param n_endo      Integer: number of endogenous variables
#' @param T           Integer: horizon
#' @return List with:
#'   $alpha      -- step length
#'   $theta_new  -- merit value at Y + alpha*delta
#'   $Y_new      -- updated path matrix (or NULL if no tried step reduced merit)
#'   $ls_iter    -- iterations used
#'   $accepted   -- logical: TRUE if the Armijo condition was met
#'
#' @details
#' On Armijo failure the search returns the smallest-merit trial point actually
#' evaluated (when it strictly improves on the incumbent), rather than
#' \code{NULL} with a discarded step. The previous behaviour -- return NULL,
#' caller applies a blind half Newton step -- could send the path to ~1e9 in a
#' single iteration and singularise the next Jacobian on cold starts (issue
#' M15). Returning the best vetted point guarantees monotone non-increase of
#' the merit function.
#' @noRd
.line_search <- function(Y, delta_vec, R_stack, J_stack,
                         theta_cur, fn_merit,
                         sigma = 1e-4, max_ls = 20L,
                         n_endo, T) {
  # Compute gradient: grad theta = J' . R
  grad <- as.numeric(Matrix::crossprod(J_stack, R_stack))
  directional_deriv <- sum(grad * delta_vec)

  # If the directional derivative is positive the Newton direction is not a
  # descent direction.  Fall back to steepest descent: Delta = -grad
  if (directional_deriv >= 0) {
    delta_vec <- -grad
    directional_deriv <- -sum(grad * grad)
  }

  alpha <- 1.0

  # Track the best (lowest-merit) trial point that strictly improves on the
  # incumbent, so an Armijo failure still yields a usable (merit-decreasing)
  # step rather than discarding all work (M15).
  best_theta <- theta_cur
  best_Y     <- NULL
  best_alpha <- 0

  for (ls_iter in seq_len(max_ls)) {
    # Trial point
    Y_trial <- Y
    for (t in seq_len(T)) {
      idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
      Y_trial[t, ] <- Y[t, ] + alpha * delta_vec[idx_t]
    }

    theta_new <- fn_merit(Y_trial)

    if (is.finite(theta_new) && theta_new < best_theta) {
      best_theta <- theta_new
      best_Y     <- Y_trial
      best_alpha <- alpha
    }

    # Armijo condition (isTRUE guards NaN merit values from invalid steps).
    if (isTRUE(theta_new <= theta_cur + sigma * alpha * directional_deriv)) {
      return(list(alpha = alpha, theta_new = theta_new,
                  Y_new = Y_trial, ls_iter = ls_iter, accepted = TRUE))
    }

    alpha <- alpha * 0.5
  }

  # Armijo never satisfied: return the best vetted (merit-decreasing) point
  # (may be NULL if none improved -- the caller decides).
  list(alpha = best_alpha, theta_new = best_theta,
       Y_new = best_Y, ls_iter = max_ls, accepted = FALSE)
}


# ---- Utility: overflow-safe log-sum-exp ----

#' Stable \code{log(sum(exp(x)))}
#'
#' Single implementation for the five open-coded
#' \code{max + log(sum(exp(x - max)))} sites that used to live in
#' mdd-calibration.R, sampler-smc.R, ms-filter.R, mdd-thames.R and
#' hank-panel.R.
#'
#' Contract: an all-\code{-Inf} input returns \code{-Inf}; an \code{NA} or
#' \code{NaN} anywhere in \code{x} returns \code{-Inf} (the degenerate-weight
#' convention the SMC/particle code has always used, so a broken likelihood
#' draw is rejected rather than poisoning the chain); \code{+Inf} propagates.
#'
#' @param x numeric vector of log-values.
#' @return scalar.
#' @noRd
.logsumexp <- function(x) {
  mx <- max(x)
  if (is.na(mx)) return(-Inf)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}
