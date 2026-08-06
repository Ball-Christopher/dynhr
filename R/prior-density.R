## R/prior-density.R
## --------------------------------------------------------------------------
## Phase-2 split from estimation-monolith.R.
##
## Per-distribution log-density functions for Bayesian estimation.
## These evaluate the log prior density for a single parameter value.
##
## Key convention: inv_gamma_pdf in Dynare.jl uses the custom InverseGamma1
## distribution (X = sqrt(Y), Y ~ InverseGamma).  IG1 != standard IG2.
## --------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## InverseGamma1 (IG1) helpers -- Dynare.jl's custom distribution.
## X ~ IG1(alpha, theta) means X = sqrt(Y) where Y ~ InverseGamma(alpha, theta).
##
## Log density:
##   log(2) + alpha*log(theta) - lgamma(alpha) - (2*alpha+1)*log(x) - theta/x^2
##
## Hyperparameters from (mean mu, std sig):
##   alpha: numerical root of
##     log(alpha-1) + log(sig^2+mu^2) - log(mu^2)
##       - 2*(lgamma(alpha) - lgamma(alpha-0.5)) = 0
##   theta = (alpha - 1) * (sig^2 + mu^2)
## ---------------------------------------------------------------------------

## alpha depends ONLY on the (fixed) prior hyperparameters (mu, sig), so the
## uniroot below is constant across every density evaluation -- memoize it.
## In an MCMC run .lp_ig1() is called once per IG parameter per draw; without
## this cache each call re-solves the root (profiled at ~6% of an RWMH draw).
.ig1_alpha_cache <- new.env(parent = emptyenv())
.ig1_alpha <- function(mu, sig) {
  key <- paste0(mu, "_", sig)
  hit <- .ig1_alpha_cache[[key]]
  if (!is.null(hit)) return(hit)
  sig2 <- sig^2; mu2 <- mu^2
  rhs_const <- log(sig2 + mu2) - log(mu2)
  obj <- function(a) log(a - 1) + rhs_const - 2 * (lgamma(a) - lgamma(a - 0.5))
  val <- uniroot(obj, interval = c(1 + 1e-9, 1e4), tol = 1e-12, extendInt = "upX")$root
  assign(key, val, envir = .ig1_alpha_cache)
  val
}

.lp_ig1 <- function(x, mu, sig) {
  if (x <= 0) return(-Inf)
  alpha <- .ig1_alpha(mu, sig)
  if (is.na(alpha)) return(-Inf)
  theta <- (alpha - 1) * (sig^2 + mu^2)
  log(2) + alpha * log(theta) - lgamma(alpha) - (2 * alpha + 1) * log(x) - theta / x^2
}

#' Evaluate log-prior density for a single parameter
#'
#' Supports: "beta", "gamma", "normal", "inv_gamma" (IG1 = Dynare.jl convention),
#' "inv_gamma2" (standard IG), "uniform".
#'
#' @param x   Parameter value
#' @param dist Distribution name (from prior_spec$distribution)
#' @param p1  First hyperparameter (mean, or lower for uniform)
#' @param p2  Second hyperparameter (std, or upper for uniform)
#' @param p3  Optional lower bound (-Inf by default)
#' @param p4  Optional upper bound (Inf by default)
#' @return Scalar log-density (finite or -Inf)
#' @noRd
log_prior_density <- function(x, dist, p1, p2, p3 = -Inf, p4 = Inf) {
  if (x < p3 || x > p4) return(-Inf)

  switch(.normalize_dist(dist),
    "beta" = {
      mu <- p1; sig <- p2
      a  <- mu * (mu*(1-mu)/sig^2 - 1)
      b  <- (1-mu) * (mu*(1-mu)/sig^2 - 1)
      if (a <= 0 || b <= 0) return(-Inf)
      dbeta(x, a, b, log = TRUE)
    },
    "gamma" = {
      mu <- p1; sig <- p2
      dgamma(x, shape = (mu/sig)^2, rate = mu/sig^2, log = TRUE)
    },
    "normal" = dnorm(x, mean = p1, sd = p2, log = TRUE),
    "inv_gamma" =, "inv_gamma1" = .lp_ig1(x, p1, p2),
    "inv_gamma2" = {
      mu <- p1; sig <- p2
      if (is.infinite(sig)) {
        nu <- 2.0001; s <- mu
      } else {
        nu <- 4 + 2 * (mu / sig)^2
        s  <- mu * (nu - 2) / nu
      }
      if (x <= 0) return(-Inf)
      (nu / 2) * log(nu * s / 2) - lgamma(nu / 2) -
        ((nu + 2) / 2) * log(x) - nu * s / (2 * x)
    },
    "uniform" = dunif(x, min = p1, max = p2, log = TRUE),
    ## Fail loud on an unrecognised distribution. Returning 0 (a flat prior)
    ## here silently turns a misspelled distribution name (e.g. "normial",
    ## "invgamma") into an IMPROPER flat prior on that parameter -- the prior is
    ## quietly dropped and the posterior is wrong with no error.
    stop("log_prior_density: unsupported prior distribution \"", dist, "\". ",
         "Supported: beta, gamma, normal, inv_gamma (= inv_gamma1), inv_gamma2, ",
         "uniform. A misspelled name would otherwise be given a silent flat prior.",
         call. = FALSE)
  )
}

#' Compute total log-prior density for a named parameter vector
#'
#' @param theta      Named numeric vector (parameter values to evaluate)
#' @param prior_spec Data.frame with columns: name, distribution, p1, p2,
#'                   lower, upper
#' @return Scalar log-prior (finite, or -Inf if any parameter is out of bounds)
#' @export
log_prior <- function(theta, prior_spec) {
  ## Hoist the data.frame column extractions and the distribution-name
  ## normalization out of the per-parameter loop: this function runs once
  ## per posterior evaluation, so per-row `$` dispatch and (especially)
  ## per-row .normalize_dist() are hot -- the latter's old trimws() cost
  ## ~200us/call on some Windows builds and dominated dynhr_benchmark()
  ## there (2026-08-05 diagnosis; see .normalize_dist's hot-path note).
  spec_name  <- prior_spec$name
  spec_dist  <- .normalize_dist(prior_spec$distribution)
  spec_p1    <- prior_spec$p1
  spec_p2    <- prior_spec$p2
  spec_lower <- prior_spec$lower
  spec_upper <- prior_spec$upper
  theta_nms  <- names(theta)

  lp <- 0
  for (i in seq_along(spec_name)) {
    nm <- spec_name[i]
    if (!(nm %in% theta_nms)) next
    x  <- theta[nm]

    ## Non-finite parameter values (NA/NaN/Inf, e.g. from a diverged HMC
    ## trajectory or a failed gradient) are outside every prior's support:
    ## return -Inf rather than letting `x < lo` evaluate to NA and crash.
    if (!is.finite(x)) return(-Inf)

    lo <- spec_lower[i]
    hi <- spec_upper[i]
    if (!is.na(lo) && x < lo) return(-Inf)
    if (!is.na(hi) && x > hi) return(-Inf)

    dist <- spec_dist[i]
    p1   <- spec_p1[i]
    p2   <- spec_p2[i]

    ll <- switch(dist,
      "inv_gamma" =, "inv_gamma1" = .lp_ig1(x, p1, p2),
      "inv_gamma2" = {
        if (x <= 0) { -Inf } else {
          shape <- (p1 / p2)^2 + 2
          scale <- p1 * (shape - 1)
          if (shape <= 2) { -log(x) }
          else { dgamma(1/x, shape = shape, rate = scale, log = TRUE) - 2*log(x) }
        }
      },
      "beta" = {
        if (x <= 0 || x >= 1) { -Inf } else {
          v <- p2^2
          a <- p1 * (p1 * (1 - p1) / v - 1)
          b <- (1 - p1) * (p1 * (1 - p1) / v - 1)
          ## Degenerate implied shape (a or b <= 0) is an INVALID beta, not a
          ## flat prior: return -Inf (mirrors log_prior_density()), don't drop it.
          if (a <= 0 || b <= 0) { -Inf } else { dbeta(x, a, b, log = TRUE) }
        }
      },
      "gamma" = {
        if (x <= 0) { -Inf }
        else { dgamma(x, shape = (p1/p2)^2, rate = p1/p2^2, log = TRUE) }
      },
      "normal"  = dnorm(x, mean = p1, sd = p2, log = TRUE),
      "uniform" = { if (x < p1 || x > p2) -Inf else -log(p2 - p1) },
      ## Fail loud on an unrecognised distribution. The posterior/likelihood
      ## closures call THIS function (not log_prior_density), so returning 0 (a
      ## flat prior) here silently turns a misspelled name into an improper flat
      ## prior on the hot path -- the bug the scoping report flagged, present in
      ## this copy after log_prior_density()/the gradient mirror were hardened.
      stop("log_prior: unsupported prior distribution \"", dist, "\". Supported: ",
           "beta, gamma, normal, inv_gamma (= inv_gamma1), inv_gamma2, uniform. ",
           "A misspelled name would otherwise be given a silent flat prior.",
           call. = FALSE)
    )

    if (!is.finite(ll)) return(-Inf)
    lp <- lp + ll
  }
  lp
}
