## R/param-transform.R
## --------------------------------------------------------------------------
## Unconstrained-parameter-transform layer (step 1 of 2).
##
## build_param_transform(prior_spec, par_names) returns an S3 object
## ("dynhr_param_transform") that maps each constrained parameter theta_j
## (living on its prior support) to an unconstrained eta_j (living on the
## whole real line), and back, together with the log-Jacobian of the
## transform and its derivatives -- the pieces needed by HMC/NUTS to sample
## in unconstrained space while still targeting the correct theta-space
## posterior.
##
## Per-parameter transform types (chosen from the resolved support (a, b)):
##   (-Inf,  Inf)  identity        eta = theta
##   (a,     Inf)  shifted log      eta = log(theta - a)
##   (-Inf,  b)    reflected log    eta = log(b - theta)
##   (a,     b)    scaled logit     eta = qlogis((theta - a) / (b - a))
##
## make_transformed_logpost() / make_transformed_grad() (step 2 helpers) wrap
## a theta-space log-posterior / gradient so samplers can operate on eta.
## --------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Support resolution
## ---------------------------------------------------------------------------

#' Resolve the (a, b) support for a single parameter from prior_spec
#'
#' Resolution order:
#'   1. Explicit finite lower/upper in prior_spec (non-NA), used as-is
#'      (one-sided is allowed: e.g. lower finite, upper = Inf).
#'   2. Distribution-implied support:
#'        beta                          -> (0, 1)
#'        gamma/inv_gamma/inv_gamma1/
#'          inv_gamma2                  -> (0, Inf)
#'        uniform                       -> (p1, p2)
#'        normal (or unknown)           -> (-Inf, Inf)
#'
#' In practice extract_prior_spec() always populates lower/upper (via
#' .default_lower / .default_upper, or directly from p1/p2 for "uniform"),
#' so step 1 fires for every distribution produced by the parser.  Step 2
#' is retained as a fallback for hand-built prior_spec data.frames (e.g.
#' test fixtures) that omit lower/upper or leave them NA.
#'
#' @param row Single-row data.frame / list with at least `distribution`,
#'   and optionally `p1`, `p2`, `lower`, `upper`.
#' @return list(a = lower bound, b = upper bound)
#' @noRd
.resolve_param_support <- function(row) {
  lo <- if (!is.null(row$lower)) row$lower else NA_real_
  hi <- if (!is.null(row$upper)) row$upper else NA_real_

  dist <- tolower(as.character(row$distribution %||% "normal"))
  p1   <- if (!is.null(row$p1)) row$p1 else NA_real_
  p2   <- if (!is.null(row$p2)) row$p2 else NA_real_

  ## Distribution-implied defaults (used to fill in NA lower/upper).
  default <- switch(dist,
    "beta" = c(0, 1),
    "gamma" =, "inv_gamma" =, "inv_gamma1" =, "inv_gamma2" = c(0, Inf),
    "uniform" = {
      lo_u <- if (!is.na(p1)) p1 else -Inf
      hi_u <- if (!is.na(p2)) p2 else Inf
      c(lo_u, hi_u)
    },
    c(-Inf, Inf)  # "normal" or unknown: real line
  )

  a <- if (!is.na(lo)) lo else default[1]
  b <- if (!is.na(hi)) hi else default[2]

  list(a = a, b = b)
}

#' Classify a (a, b) support into one of the four transform types
#' @noRd
.classify_support <- function(a, b) {
  a_finite <- is.finite(a)
  b_finite <- is.finite(b)
  if (!a_finite && !b_finite) {
    "identity"
  } else if (a_finite && !b_finite) {
    "log"
  } else if (!a_finite && b_finite) {
    "reflected_log"
  } else {
    "logit"
  }
}

## ---------------------------------------------------------------------------
## build_param_transform()
## ---------------------------------------------------------------------------

#' Build an unconstrained-parameter transform layer
#'
#' @param prior_spec data.frame (from extract_prior_spec()) with columns
#'   `name`, `distribution`, `p1`, `p2`, `lower`, `upper` (the latter four
#'   may be NA / absent for a given row).
#' @param par_names  Character vector of parameter names, in the order they
#'   appear in the parameter vector theta the samplers operate on.
#' @return An S3 list of class "dynhr_param_transform" with elements:
#'   \itemize{
#'     \item `types`: named character vector, one of "identity", "log",
#'       "reflected_log", "logit" per parameter
#'     \item `a`, `b`: named numeric vectors of resolved support bounds
#'     \item `to_unconstrained(theta)`: theta -> eta (named, vectorized)
#'     \item `to_constrained(eta)`: eta -> theta (named, vectorized)
#'     \item `log_jacobian(eta)`: scalar sum of log|d theta/d eta|
#'     \item `dlog_jacobian(eta)`: named vector d/deta_j log|J|
#'     \item `dtheta_deta(eta)`: named vector d theta_j / d eta_j
#'   }
#' @noRd
build_param_transform <- function(prior_spec, par_names) {
  n <- length(par_names)
  types <- setNames(character(n), par_names)
  a_vec <- setNames(numeric(n), par_names)
  b_vec <- setNames(numeric(n), par_names)

  for (j in seq_len(n)) {
    nm  <- par_names[j]
    idx <- match(nm, prior_spec$name)

    if (is.na(idx)) {
      warning(sprintf(
        "build_param_transform: parameter '%s' not found in prior_spec; using identity transform.",
        nm))
      types[nm] <- "identity"
      a_vec[nm] <- -Inf
      b_vec[nm] <-  Inf
      next
    }

    row <- prior_spec[idx, , drop = FALSE]
    sup <- .resolve_param_support(row)
    types[nm] <- .classify_support(sup$a, sup$b)
    a_vec[nm] <- sup$a
    b_vec[nm] <- sup$b
  }

  ## ---- elementwise scalar transforms (operate on a single (theta, type, a, b)) ----

  ## Relative clamp factor for boundary safety: 1e-10 of the local scale.
  .clamp_eps <- 1e-10

  .theta_to_eta1 <- function(theta, type, a, b) {
    switch(type,
      "identity" = theta,
      "log" = {
        d <- theta - a
        if (d <= 0) {
          scale <- max(abs(a), abs(theta), 1)
          d <- .clamp_eps * scale
        }
        log(d)
      },
      "reflected_log" = {
        d <- b - theta
        if (d <= 0) {
          scale <- max(abs(b), abs(theta), 1)
          d <- .clamp_eps * scale
        }
        log(d)
      },
      "logit" = {
        u <- (theta - a) / (b - a)
        if (u <= 0) u <- .clamp_eps
        if (u >= 1) u <- 1 - .clamp_eps
        qlogis(u)
      },
      stop("unknown transform type: ", type)
    )
  }

  .eta_to_theta1 <- function(eta, type, a, b) {
    switch(type,
      "identity" = eta,
      "log" = a + exp(eta),
      "reflected_log" = b - exp(eta),
      "logit" = a + (b - a) * plogis(eta),
      stop("unknown transform type: ", type)
    )
  }

  .log_jacobian1 <- function(eta, type, a, b) {
    switch(type,
      "identity" = 0,
      "log" = eta,
      "reflected_log" = eta,
      "logit" = log(b - a) + plogis(eta, log.p = TRUE) + plogis(-eta, log.p = TRUE),
      stop("unknown transform type: ", type)
    )
  }

  .dlog_jacobian1 <- function(eta, type) {
    switch(type,
      "identity" = 0,
      "log" = 1,
      "reflected_log" = 1,
      "logit" = 1 - 2 * plogis(eta),
      stop("unknown transform type: ", type)
    )
  }

  .dtheta_deta1 <- function(eta, type, a, b) {
    switch(type,
      "identity" = 1,
      "log" = exp(eta),
      "reflected_log" = -exp(eta),
      "logit" = (b - a) * plogis(eta) * plogis(-eta),
      stop("unknown transform type: ", type)
    )
  }

  ## ---- vectorized public API ----

  to_unconstrained <- function(theta) {
    if (is.null(names(theta))) names(theta) <- par_names
    eta <- setNames(numeric(n), par_names)
    for (j in seq_len(n)) {
      nm <- par_names[j]
      eta[nm] <- .theta_to_eta1(theta[[nm]], types[[nm]], a_vec[[nm]], b_vec[[nm]])
    }
    eta
  }

  to_constrained <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    theta <- setNames(numeric(n), par_names)
    for (j in seq_len(n)) {
      nm <- par_names[j]
      theta[nm] <- .eta_to_theta1(eta[[nm]], types[[nm]], a_vec[[nm]], b_vec[[nm]])
    }
    theta
  }

  log_jacobian <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    total <- 0
    for (j in seq_len(n)) {
      nm <- par_names[j]
      total <- total + .log_jacobian1(eta[[nm]], types[[nm]], a_vec[[nm]], b_vec[[nm]])
    }
    total
  }

  dlog_jacobian <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    out <- setNames(numeric(n), par_names)
    for (j in seq_len(n)) {
      nm <- par_names[j]
      out[nm] <- .dlog_jacobian1(eta[[nm]], types[[nm]])
    }
    out
  }

  dtheta_deta <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    out <- setNames(numeric(n), par_names)
    for (j in seq_len(n)) {
      nm <- par_names[j]
      out[nm] <- .dtheta_deta1(eta[[nm]], types[[nm]], a_vec[[nm]], b_vec[[nm]])
    }
    out
  }

  structure(
    list(
      par_names        = par_names,
      types            = types,
      a                = a_vec,
      b                = b_vec,
      to_unconstrained = to_unconstrained,
      to_constrained   = to_constrained,
      log_jacobian     = log_jacobian,
      dlog_jacobian    = dlog_jacobian,
      dtheta_deta      = dtheta_deta
    ),
    class = "dynhr_param_transform"
  )
}

## ---------------------------------------------------------------------------
## Step-2 helpers: wrap a theta-space log-posterior / gradient for eta-space
## samplers.  Defined here (alongside the transform) so step 2 can use them
## without further plumbing.
## ---------------------------------------------------------------------------

#' Wrap a theta-space log-posterior as an eta-space (unconstrained) function
#'
#' @param log_post_fn Function(theta) -> scalar logpost, or
#'   list(logpost = ..., loglik = ..., logprior = ..., ...)
#' @param tr          A "dynhr_param_transform" (from build_param_transform())
#' @param include_jacobian Logical; if TRUE (default) add `tr$log_jacobian(eta)`
#'   to the returned logpost (the change-of-variables correction needed when
#'   sampling eta with a target proportional to p(theta) * |d theta/d eta|).
#' @return function(eta) -> same shape as log_post_fn's return value (scalar
#'   or list), with `logpost` (and the bare scalar) adjusted by the Jacobian.
#'   If the underlying logpost is -Inf, -Inf is returned unchanged (no
#'   Jacobian added).
#' @noRd
make_transformed_logpost <- function(log_post_fn, tr, include_jacobian = TRUE) {
  function(eta) {
    if (is.null(names(eta))) names(eta) <- tr$par_names
    theta <- tr$to_constrained(eta)
    res   <- log_post_fn(theta)

    if (is.list(res)) {
      lp <- res$logpost
      if (!is.finite(lp)) return(res)
      if (include_jacobian) res$logpost <- lp + tr$log_jacobian(eta)
      res
    } else {
      lp <- res
      if (!is.finite(lp)) return(lp)
      if (include_jacobian) lp <- lp + tr$log_jacobian(eta)
      lp
    }
  }
}

#' Wrap a theta-space gradient as an eta-space gradient (chain rule)
#'
#' d/deta_j log p(theta(eta)) + log|J(eta)|
#'   = dtheta_j/deta_j * grad_fn(theta)_j + dlog_jacobian(eta)_j
#'
#' @param grad_fn Function(theta) -> named numeric vector, theta-space
#'   gradient of the (unadjusted) log-posterior.
#' @param tr      A "dynhr_param_transform" (from build_param_transform())
#' @return function(eta) -> named numeric vector, eta-space gradient of the
#'   Jacobian-adjusted log-posterior.
#' @noRd
make_transformed_grad <- function(grad_fn, tr) {
  function(eta) {
    if (is.null(names(eta))) names(eta) <- tr$par_names
    theta <- tr$to_constrained(eta)
    tr$dtheta_deta(eta) * grad_fn(theta) + tr$dlog_jacobian(eta)
  }
}

#' Delta-method conversion of a theta-space covariance to eta-space
#'
#' eta = g(theta)  =>  Var(eta) ~ J Var(theta) J'  with J = diag(deta/dtheta)
#'                            = diag(1/dtheta_deta) Var(theta) diag(1/dtheta_deta)
#'
#' i.e. Sigma_eta = D^{-1} Sigma_theta D^{-1}, D = diag(dtheta_deta(eta_at)).
#' D entries are floored at 1e-12 in absolute value before inverting to avoid
#' blow-up when a parameter sits very near a transform's asymptote (e.g.
#' eta -> -Inf for a "log" transform, where dtheta/deta -> 0).
#'
#' @param Sigma_theta Theta-space covariance matrix, or NULL.
#' @param tr          A "dynhr_param_transform" (from build_param_transform()).
#' @param theta_at    Theta-space point at which to evaluate the Jacobian
#'   (typically the mode or current chain position).
#' @return Eta-space covariance matrix with the same dimnames as
#'   \code{Sigma_theta}, or NULL if \code{Sigma_theta} is NULL / not a matrix.
#' @noRd
.cov_theta_to_eta <- function(Sigma_theta, tr, theta_at) {
  if (is.null(Sigma_theta) || !is.matrix(Sigma_theta)) return(NULL)
  eta_at <- tr$to_unconstrained(theta_at)
  d_vec  <- tr$dtheta_deta(eta_at)
  d_vec[abs(d_vec) < 1e-12] <- sign(d_vec[abs(d_vec) < 1e-12]) * 1e-12
  d_vec[d_vec == 0] <- 1e-12
  Dinv <- diag(1 / d_vec, nrow = length(d_vec))
  Sigma_eta <- Dinv %*% Sigma_theta %*% Dinv
  dimnames(Sigma_eta) <- dimnames(Sigma_theta)
  Sigma_eta
}

#' Delta-method conversion of a theta-space PRECISION/metric to eta-space
#'
#' Unlike a covariance, a precision-like object (Fisher information, Hessian
#' of neg-log-posterior) transforms with the JACOBIAN, not its inverse:
#'
#'   eta = g(theta)  =>  G_eta = J^{-T} G_theta J^{-1},  J = diag(deta/dtheta)
#'                             = diag(dtheta/deta) G_theta diag(dtheta/deta)
#'                             = D G_theta D,  D = diag(dtheta_deta(eta_at))
#'
#' i.e. G_eta = D^T G_theta D (D is diagonal so D^T = D). This is the inverse
#' rule of \code{.cov_theta_to_eta()}: consistent, since for a positive-definite
#' G_theta, \code{solve(.fim_theta_to_eta(G_theta, ...))} recovers
#' \code{.cov_theta_to_eta(solve(G_theta), ...)} exactly (see
#' test-nuts-whittle-metric.R). Same floor-guard on \code{dtheta_deta} as
#' \code{.cov_theta_to_eta()}, to avoid blow-up near a transform's asymptote.
#'
#' @param G_theta  Theta-space precision/metric matrix, or NULL.
#' @param tr       A "dynhr_param_transform" (from build_param_transform()).
#' @param theta_at Theta-space point at which to evaluate the Jacobian
#'   (typically the mode).
#' @return Eta-space precision/metric matrix with the same dimnames as
#'   \code{G_theta}, or NULL if \code{G_theta} is NULL / not a matrix.
#' @noRd
.fim_theta_to_eta <- function(G_theta, tr, theta_at) {
  if (is.null(G_theta) || !is.matrix(G_theta)) return(NULL)
  eta_at <- tr$to_unconstrained(theta_at)
  d_vec  <- tr$dtheta_deta(eta_at)
  d_vec[abs(d_vec) < 1e-12] <- sign(d_vec[abs(d_vec) < 1e-12]) * 1e-12
  d_vec[d_vec == 0] <- 1e-12
  D <- diag(d_vec, nrow = length(d_vec))
  G_eta <- D %*% G_theta %*% D
  dimnames(G_eta) <- dimnames(G_theta)
  G_eta
}
