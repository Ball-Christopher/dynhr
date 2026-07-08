## R/mode-trust-region.R
## --------------------------------------------------------------------------
## Trust-region Newton optimizer (Nocedal & Wright, Algorithms 4.1 + dogleg
## 4.3) with a positive-definite Hessian modification.  A deterministic
## workhorse for when gradients and Hessians are available (dynhr has analytic
## adjoint gradients/Hessians), complementing the stochastic/global mode
## finders (CMA-ES, JADE) and the quasi-Newton newrat/csminwel default.
##
## MINIMIZES fn.  To find a posterior MODE, pass the NEGATIVE log posterior
## (and its negative gradient/Hessian); see the example.
##
## This is a self-contained leaf optimizer: it is not threaded into the parallel
## run_mode_finding() dispatch (that is a possible follow-up), so it cannot
## introduce a silent fidelity fallback elsewhere.
## --------------------------------------------------------------------------


#' Positive-definite modification of a symmetric matrix (eigenvalue clamp)
#'
#' Returns \code{B} with its eigenvalues clamped to at least \code{tau} times
#' the largest eigenvalue (and at least \code{tau} in absolute terms), so the
#' trust-region model is strictly convex and the Newton step well-defined even
#' when the true Hessian is indefinite.
#'
#' @param B Symmetric numeric matrix.
#' @param tau Relative floor for eigenvalues (default 1e-8).
#' @return A symmetric positive-definite matrix.
#' @keywords internal
.tr_modify_pd <- function(B, tau = 1e-8) {
  B <- (B + t(B)) / 2
   e <- eigen(B, symmetric = TRUE)
  floor_val <- max(tau, tau * max(abs(e$values)))
  d <- pmax(e$values, floor_val)
  e$vectors %*% (d * t(e$vectors))
}


#' Dogleg trust-region subproblem step
#'
#' Approximately minimizes the quadratic model \eqn{g'p + 0.5 p'B p} subject to
#' \eqn{\|p\| \le \Delta} using the dogleg path between the Cauchy point and the
#' (PD) Newton step.
#'
#' @param g Gradient vector.
#' @param B Positive-definite model Hessian.
#' @param delta Trust-region radius.
#' @return The step vector \code{p}.
#' @keywords internal
.tr_dogleg <- function(g, B, delta) {
  ## Full Newton step.
  pN <- tryCatch(-solve(B, g), error = function(e) rep(NA_real_, length(g)))
  if (all(is.finite(pN)) && sqrt(sum(pN^2)) <= delta) return(pN)

  ## Cauchy (steepest-descent) step.
  gBg <- as.numeric(t(g) %*% B %*% g)
  gg  <- sum(g^2)
  if (gBg <= 0) {                       # shouldn't happen for PD B, but guard
    return(-delta * g / sqrt(gg))
  }
  pU <- -(gg / gBg) * g                 # unconstrained Cauchy minimizer
  if (sqrt(sum(pU^2)) >= delta)         # Cauchy point outside -> scale to bound
    return(-delta * g / sqrt(gg))

  if (!all(is.finite(pN)))              # no usable Newton step -> Cauchy
    return(pU)

  ## Dogleg: find tau in [0,1] s.t. ||pU + tau (pN - pU)|| = delta.
  d <- pN - pU
  a <- sum(d^2); b <- 2 * sum(pU * d); c <- sum(pU^2) - delta^2
  tau <- (-b + sqrt(max(b^2 - 4 * a * c, 0))) / (2 * a)
  pU + tau * d
}


#' Trust-region Newton minimizer
#'
#' Minimizes \code{fn} from \code{x0} using a dogleg trust-region method with a
#' positive-definite Hessian modification.
#'
#' @param fn Objective function to MINIMIZE: \code{fn(x, ...)} returns a scalar.
#' @param x0 Numeric starting vector.
#' @param gr Optional gradient \code{gr(x, ...)}; if \code{NULL}, a central
#'   finite-difference gradient is used.
#' @param he Optional Hessian \code{he(x, ...)}; if \code{NULL}, a
#'   finite-difference Hessian (via \code{numDeriv} if available) is used.
#' @param control List of tuning parameters: \code{delta0} (initial radius, 1),
#'   \code{delta_max} (max radius, 1e3), \code{eta} (acceptance threshold, 0.1),
#'   \code{tol_g} (gradient-norm tolerance, 1e-8), \code{tol_step} (step
#'   tolerance, 1e-12), \code{maxit} (500), \code{verbose} (FALSE).
#' @param ... Passed to \code{fn}/\code{gr}/\code{he}.
#'
#' @return A list with \code{par}, \code{value}, \code{gradient},
#'   \code{hessian}, \code{iterations}, \code{converged}, \code{delta}, and
#'   \code{convergence_message}.
#'
#' @examples
#' ## Rosenbrock (non-convex): converges to (1, 1).
#' f <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
#' fit <- mode_trust_region(f, x0 = c(-1.2, 1))
#' fit$par
#'
#' \dontrun{
#' ## Posterior mode: minimize the NEGATIVE log posterior.
#' lp  <- make_log_posterior(...)
#' fit <- mode_trust_region(function(t) -lp(t), x0 = theta_init)
#' }
#' @export
mode_trust_region <- function(fn, x0, gr = NULL, he = NULL,
                              control = list(), ...) {
  ctrl <- utils::modifyList(list(
    delta0 = 1, delta_max = 1e3, eta = 0.1,
    tol_g = 1e-8, tol_step = 1e-12, maxit = 500L, verbose = FALSE), control)

  grad_fn <- if (!is.null(gr)) function(x) gr(x, ...) else function(x)
    .tr_num_grad(fn, x, ...)
  hess_fn <- if (!is.null(he)) function(x) he(x, ...) else function(x) {
    if (requireNamespace("numDeriv", quietly = TRUE))
      numDeriv::hessian(function(z) fn(z, ...), x)
    else .tr_num_hess(fn, x, ...)
  }

  x  <- x0
  fx <- fn(x, ...)
  g  <- grad_fn(x)
  delta <- ctrl$delta0
  converged <- FALSE
  msg <- "maxit reached"
  it <- 0L

  for (it in seq_len(ctrl$maxit)) {
    if (sqrt(sum(g^2)) < ctrl$tol_g) {
      converged <- TRUE; msg <- "gradient tolerance"; break
    }
    B  <- .tr_modify_pd(hess_fn(x))
    p  <- .tr_dogleg(g, B, delta)
    np <- sqrt(sum(p^2))

    pred_red <- -(sum(g * p) + 0.5 * as.numeric(t(p) %*% B %*% p))
    fx_new   <- fn(x + p, ...)
    act_red  <- fx - fx_new
    rho <- if (pred_red > 0) act_red / pred_red else -Inf

    ## Radius update.
    if (rho < 0.25) {
      delta <- 0.25 * delta
    } else if (rho > 0.75 && np >= 0.99 * delta) {
      delta <- min(2 * delta, ctrl$delta_max)
    }
    ## Accept / reject.
    if (rho > ctrl$eta) {
      x <- x + p; fx <- fx_new; g <- grad_fn(x)
    }
    if (ctrl$verbose)
      cat(sprintf("  it %d: f=%.8g |g|=%.3e delta=%.3e rho=%.3f\n",
                  it, fx, sqrt(sum(g^2)), delta, rho))
    if (np < ctrl$tol_step && rho > ctrl$eta) {
      converged <- TRUE; msg <- "step tolerance"; break
    }
  }

  list(par = x, value = fx, gradient = g, hessian = hess_fn(x),
       iterations = it, converged = converged, delta = delta,
       convergence_message = msg)
}


#' Central finite-difference gradient (fallback)
#' @keywords internal
.tr_num_grad <- function(fn, x, ..., h = 1e-6) {
  n <- length(x); g <- numeric(n)
  for (i in seq_len(n)) {
    step <- h * max(1, abs(x[i]))
    xp <- x; xm <- x; xp[i] <- xp[i] + step; xm[i] <- xm[i] - step
    g[i] <- (fn(xp, ...) - fn(xm, ...)) / (2 * step)
  }
  g
}


#' Finite-difference Hessian (fallback when numDeriv is unavailable)
#' @keywords internal
.tr_num_hess <- function(fn, x, ..., h = 1e-4) {
  n <- length(x); H <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in i:n) {
    hi <- h * max(1, abs(x[i])); hj <- h * max(1, abs(x[j]))
    xpp <- x; xpp[i] <- xpp[i] + hi; xpp[j] <- xpp[j] + hj
    xpm <- x; xpm[i] <- xpm[i] + hi; xpm[j] <- xpm[j] - hj
    xmp <- x; xmp[i] <- xmp[i] - hi; xmp[j] <- xmp[j] + hj
    xmm <- x; xmm[i] <- xmm[i] - hi; xmm[j] <- xmm[j] - hj
    H[i, j] <- (fn(xpp, ...) - fn(xpm, ...) - fn(xmp, ...) + fn(xmm, ...)) /
      (4 * hi * hj)
    H[j, i] <- H[i, j]
  }
  H
}
