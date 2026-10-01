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

#' \code{alpha - 1} for the IG1 prior with mean \code{mu} and std \code{sig}
#'
#' Two fixes over the old \code{.ig1_alpha()}:
#'
#'  * The root is solved in \code{u = log(alpha - 1)} rather than in \code{alpha}. The
#'    old solve bracketed \code{alpha} on \code{[1 + 1e-9, 1e4]}, but a large \code{sig}
#'    drives \code{alpha - 1} far BELOW 1e-9 (sig = 1e5 with mu = 1.2 needs
#'    ~4.6e-11), and \code{extendInt = "upX"} only extends upward -- so the solve
#'    ERRORED for any diffuse-but-finite prior sd. In \code{u} the objective
#'    \code{u + log((sig^2 + mu^2)/mu^2) - 2*(lgamma(alpha) - lgamma(alpha - 1/2))}
#'    is increasing and unbounded below, so a root always exists.
#'  * \code{alpha - 1} is returned, NOT \code{alpha}. For sig >= ~1e9, \code{1 + (alpha - 1)}
#'    rounds to exactly 1 in double precision, and the old
#'    \code{theta = (alpha - 1) * (sig^2 + mu^2)} then collapsed to 0.
#' @noRd
.ig1_alpha_minus_1 <- function(mu, sig) {
  ## sd = Inf is the alpha -> 1 LIMIT, not an error: E[X^2] = theta/(alpha - 1)
  ## for X ~ IG1, so a finite sd requires alpha > 1 and sd = Inf pins alpha = 1.
  ## The old code fed uniroot an all-+Inf objective and died, which is why D6
  ## reported "no prior" for such parameters.
  if (is.infinite(sig)) return(0)
  key <- paste0(mu, "_", sig)
  hit <- .ig1_alpha_cache[[key]]
  if (!is.null(hit)) return(hit)
  rhs_const <- log(sig^2 + mu^2) - log(mu^2)
  obj <- function(u) {
    a <- 1 + exp(u)
    u + rhs_const - 2 * (lgamma(a) - lgamma(a - 0.5))
  }
  u <- uniroot(obj, interval = c(-700, log(1e4)), tol = 1e-14,
               extendInt = "upX")$root
  val <- exp(u)
  assign(key, val, envir = .ig1_alpha_cache)
  val
}

#' Shape alpha of the IG1 prior (see \code{.ig1_alpha_minus_1()} for the caveat:
#' for very diffuse priors alpha is numerically 1 and only \code{alpha - 1} carries
#' information, so \code{theta} must be built from that).
#' @noRd
.ig1_alpha <- function(mu, sig) 1 + .ig1_alpha_minus_1(mu, sig)

#' (alpha, theta) of the IG1 prior with mean \code{mu} and std \code{sig}
#'
#' The ONE mapping. \code{theta = (alpha - 1) * (sig^2 + mu^2)} is 0 * Inf at
#' sig = Inf, so the limit is taken in closed form. With alpha -> 1 the mean is
#' \code{sqrt(theta) * Gamma(alpha - 1/2)/Gamma(alpha) -> sqrt(theta * pi)}, so
#' matching the mean gives \code{theta = mu^2 / pi} -- the continuous limit of the
#' finite-sd parameterisation, which keeps the MEAN and sends the variance to
#' infinity. (It coincides with Dynare's \code{nu = 2, s = 2*mu^2/pi} convention for
#' an infinite-variance inverse gamma 1.)
#' @noRd
.ig1_params <- function(mu, sig) {
  am1 <- .ig1_alpha_minus_1(mu, sig)
  theta <- if (is.infinite(sig)) mu^2 / pi else am1 * (sig^2 + mu^2)
  list(alpha = 1 + am1, theta = theta)
}

.lp_ig1 <- function(x, mu, sig) {
  if (x <= 0) return(-Inf)
  ps <- .ig1_params(mu, sig)
  if (is.na(ps$alpha)) return(-Inf)
  log(2) + ps$alpha * log(ps$theta) - lgamma(ps$alpha) -
    (2 * ps$alpha + 1) * log(x) - ps$theta / x^2
}

#' (shape, scale) of the IG2 prior with mean \code{mu} and std \code{sig}
#'
#' The ONE mapping. \code{shape = (mu/sig)^2 + 2}, \code{scale = mu*(shape - 1)},
#' which for sig = Inf is simply \code{shape = 2}, \code{scale = mu} -- a PROPER inverse
#' gamma with mean \code{mu} and infinite variance (sd finite requires shape > 2, so
#' sd = Inf pins shape = 2). No special case is needed, and the two old ones
#' were both wrong and mutually inconsistent: \code{log_prior_density()} used
#' \code{nu = 2.0001, s = mu} (shape ~ 1, i.e. an INFINITE-MEAN prior) while
#' \code{log_prior()} -- the one the samplers call -- fell into a \code{shape <= 2} branch
#' and returned the IMPROPER Jeffreys density \code{-log(x)}.
#' @noRd
.ig2_params <- function(mu, sig) {
  shape <- (mu / sig)^2 + 2
  list(shape = shape, scale = mu * (shape - 1))
}

.lp_ig2 <- function(x, mu, sig) {
  if (x <= 0) return(-Inf)
  ps <- .ig2_params(mu, sig)
  ## a*log(b) - lgamma(a) - (a+1)*log(x) - b/x
  ps$shape * log(ps$scale) - lgamma(ps$shape) -
    (ps$shape + 1) * log(x) - ps$scale / x
}

## ---------------------------------------------------------------------------
## Dynare's generalised / shifted priors (PRIOR_3RD_PARAMETER / _4TH_)
## ---------------------------------------------------------------------------
## Dynare 7.1 `estimation/set_prior.m` + `estimation/priordens.m` +
## `estimation/prior_bounds.m` (and the `estimated_params` section of the
## reference manual) give p3/p4 a meaning that depends on the prior SHAPE:
##
##   BETA            generalised beta on [p3, p4] (defaults 0, 1); p1/p2 are the
##                   mean/sd of the parameter ON THAT INTERVAL
##                   (beta_specification(mu, sigma2, p3, p4); lpdfgbeta).
##   GAMMA           p3 is a SHIFT (lower end of the support, default 0); p1/p2
##                   are the mean/sd of the parameter, so the gamma law is fitted
##                   to mean p1 - p3 (gamma_specification: `mu = mu - lb`;
##                   priordens: lpdfgam(x - p3, ...)). p4 is not used.
##   INV_GAMMA(1|2)  p3 is a SHIFT, exactly as for GAMMA (priordens:
##                   lpdfig1(x - p3, ...) / lpdfig2(x - p3, ...)). p4 is not
##                   used. NOTE: Dynare 7.1's inverse_gamma_specification() for
##                   TYPE 2 builds `s = 2*mu*(1 + mu2/sigma2)` with the
##                   UNSHIFTED mu (type 1 uses `mu - lb` throughout), so a
##                   shifted IG2 in Dynare does not have mean p1. dynhr applies
##                   the shift consistently (the law of X - p3 is fitted to mean
##                   p1 - p3), which coincides with Dynare whenever p3 = 0.
##   NORMAL          p3/p4 are TRUNCATION bounds only (prior_bounds case 3);
##                   the density is the untruncated normal kernel.
##   UNIFORM         Dynare uniform_specification(): p1/p2 are the MEAN and
##                   STANDARD DEVIATION, so the support is
##                   [p1 - sqrt(3) p2, p1 + sqrt(3) p2]; when p3 AND p4 are
##                   both given they ARE the support (and p1/p2 are ignored /
##                   empty). extract_prior_spec() emits all four, as Dynare's
##                   bayestopt_ does: p1/p2 = mean/sd, p3/p4 = the bounds.
##
## In a prior spec the columns `p3`/`p4` carry exactly the generalisation
## parameters above (NA = the default); truncation lives in `lower`/`upper`.
## NA p3/p4 reproduce the ungeneralised laws BIT-FOR-BIT (x - 0, (x - 0)/1).

#' Generalisation parameters (a, b) of a prior, with Dynare's defaults
#'
#' \code{a} is the lower end of the natural support (beta: p3 or 0; gamma and both
#' inverse gammas: the shift p3 or 0) and \code{b} the upper end (beta: p4 or 1).
#' @noRd
.prior_gen_ab <- function(p3, p4) {
  a <- if (length(p3) == 0L || is.na(p3)) 0 else p3
  b <- if (length(p4) == 0L || is.na(p4)) 1 else p4
  c(a, b)
}

#' Support [a, b] of a UNIFORM prior -- Dynare's uniform_specification()
#'
#' \code{p3}/\code{p4} when BOTH are given (the bound form, \code{uniform_pdf, , , lb, ub}),
#' otherwise the mean/sd form \code{[p1 - sqrt(3) p2, p1 + sqrt(3) p2]}. A missing
#' or degenerate row yields a non-finite or empty interval; callers reject it.
#' @noRd
.uniform_ab <- function(p1, p2, p3 = NA_real_, p4 = NA_real_) {
  if (length(p3) == 1L && length(p4) == 1L && !is.na(p3) && !is.na(p4))
    return(c(p3, p4))
  if (length(p1) != 1L || length(p2) != 1L) return(c(NA_real_, NA_real_))
  h <- sqrt(3) * p2
  c(p1 - h, p1 + h)
}

#' Beta shapes (alpha, beta) of a (generalised) beta prior on [a, b] with mean
#' \code{p1} and sd \code{p2} -- Dynare's beta_specification().
#' @noRd
.beta_shapes <- function(p1, p2, a = 0, b = 1) {
  len <- b - a
  m <- (p1 - a) / len
  v <- (p2 / len)^2
  k <- m * (1 - m) / v - 1
  c(m * k, (1 - m) * k)
}

#' Log-density of ONE prior at ONE point, given an ALREADY-NORMALISED name
#'
#' \code{log_prior_density()} and \code{log_prior()} used to carry two independent
#' switch statements that disagreed (inv_gamma2 with sd = Inf most visibly).
#' Both now call this, so there is exactly one density per (dist, p1..p4).
#' Truncation bounds (\code{lower}/\code{upper}) are applied by the callers; \code{p3}/\code{p4}
#' are Dynare's generalisation parameters (see the block comment above):
#' beta support [p3, p4], gamma / inverse-gamma shift p3. NA = default.
#' @noRd
.lp_dist1 <- function(x, dist, p1, p2, p3 = NA_real_, p4 = NA_real_) {
  switch(dist,
    "beta" = {
      ## Hot path (log_prior runs this per parameter per draw): inline NA
      ## defaults rather than calling .prior_gen_ab().
      a <- if (is.na(p3)) 0 else p3
      b <- if (is.na(p4)) 1 else p4
      if (x <= a || x >= b) return(-Inf)
      sh <- .beta_shapes(p1, p2, a, b)
      ## A degenerate implied shape is an INVALID beta, not a flat prior.
      if (!is.finite(sh[1]) || !is.finite(sh[2]) || sh[1] <= 0 || sh[2] <= 0)
        return(-Inf)
      len <- b - a
      ## Dynare lpdfgbeta: log Beta((x - a)/(b - a)) - log(b - a).
      dbeta((x - a) / len, sh[1], sh[2], log = TRUE) - log(len)
    },
    "gamma" = {
      s <- if (is.na(p3)) 0 else p3
      y <- x - s
      if (y <= 0) return(-Inf)
      m <- p1 - s
      dgamma(y, shape = (m / p2)^2, rate = m / p2^2, log = TRUE)
    },
    "normal" = dnorm(x, mean = p1, sd = p2, log = TRUE),
    "inv_gamma" =, "inv_gamma1" = {
      s <- if (is.na(p3)) 0 else p3
      .lp_ig1(x - s, p1 - s, p2)
    },
    "inv_gamma2" = {
      s <- if (is.na(p3)) 0 else p3
      .lp_ig2(x - s, p1 - s, p2)
    },
    "uniform" = {
      ab <- .uniform_ab(p1, p2, p3, p4)
      if (!isTRUE(ab[1] < ab[2]) || !is.finite(ab[2] - ab[1]) ||
          x < ab[1] || x > ab[2]) -Inf else -log(ab[2] - ab[1])
    },
    ## Fail loud on an unrecognised distribution. Returning 0 (a flat prior)
    ## here silently turns a misspelled distribution name (e.g. "normial",
    ## "invgamma") into an IMPROPER flat prior on that parameter -- the prior is
    ## quietly dropped and the posterior is wrong with no error.
    .dynhr_abort(
      "unsupported prior distribution \"", dist, "\". Supported: beta, gamma, ",
      "normal, inv_gamma (= inv_gamma1), inv_gamma2, uniform. A misspelled ",
      "name would otherwise be given a silent flat prior.",
      class = "dynhr_error_unsupported_prior")
  )
}

#' The law \code{.lp_dist1()} scores, as (natural support, r, p, q) closures
#'
#' ONE place that turns (dist, p1..p4) into a sampler, a CDF and a quantile
#' function of the SAME distribution \code{.lp_dist1()} evaluates -- generalised
#' beta on [p3, p4], shift p3 for gamma / inverse gammas, IG1 = Dynare.jl's
#' \code{X = sqrt(Y)}, IG2 standard -- through the hyperparameter mappings the
#' density uses (\code{.beta_shapes}, \code{.ig1_params}, \code{.ig2_params}).
#' @param d Canonical (already normalised) distribution name.
#' @return list(lo, hi, r(n), p(x), q(u)); \code{p} is only called for x strictly
#'   inside (lo, hi) and \code{q} for u in [0, 1].
#' @noRd
.prior_law <- function(d, p1, p2, p3 = NA_real_, p4 = NA_real_, label = d) {
  bad_shape <- function(what)
    .dynhr_abort(
      "prior for \"", label, "\" (", d, ", p1 = ", p1, ", p2 = ", p2,
      "): ", what, " -- no proper distribution to draw from.",
      class = "dynhr_error_prior_spec")
  switch(d,
    "normal" = list(
      lo = -Inf, hi = Inf,
      r = function(n) rnorm(n, p1, p2),
      p = function(x) pnorm(x, p1, p2),
      q = function(u) qnorm(u, p1, p2)),
    "beta" = {
      ab <- .prior_gen_ab(p3, p4); a <- ab[1]; b <- ab[2]; len <- b - a
      sh <- .beta_shapes(p1, p2, a, b)
      if (!is.finite(sh[1]) || !is.finite(sh[2]) || sh[1] <= 0 || sh[2] <= 0)
        bad_shape(paste0("the implied beta shapes are not positive (sd too ",
                         "large for the mean on [", a, ", ", b, "])"))
      al <- sh[1]; be <- sh[2]
      list(lo = a, hi = b,
           r = function(n) a + len * rbeta(n, al, be),
           p = function(x) stats::pbeta((x - a) / len, al, be),
           q = function(u) a + len * stats::qbeta(u, al, be))
    },
    "gamma" = {
      s <- .prior_gen_ab(p3, NA_real_)[1]; m <- p1 - s
      shape <- (m / p2)^2; rate <- m / p2^2
      if (!(is.finite(shape) && shape > 0 && is.finite(rate) && rate > 0))
        bad_shape("the implied gamma shape/rate are not positive")
      list(lo = s, hi = Inf,
           r = function(n) s + rgamma(n, shape = shape, rate = rate),
           p = function(x) stats::pgamma(x - s, shape = shape, rate = rate),
           q = function(u) s + stats::qgamma(u, shape = shape, rate = rate))
    },
    "inv_gamma" =, "inv_gamma1" = {
      ## X - s = 1/sqrt(G), G ~ Gamma(alpha, rate = theta) (the square of an
      ## IG1 variable is IG2(alpha, theta)); X is DECREASING in G, so the
      ## tails swap in p/q.
      s <- .prior_gen_ab(p3, NA_real_)[1]
      ps <- .ig1_params(p1 - s, p2); al <- ps$alpha; th <- ps$theta
      if (!(is.finite(al) && al > 0 && is.finite(th) && th > 0))
        bad_shape("the implied inverse-gamma-1 parameters are not positive")
      list(lo = s, hi = Inf,
           r = function(n) s + 1 / sqrt(rgamma(n, shape = al, rate = th)),
           p = function(x) 1 - stats::pgamma(1 / (x - s)^2, shape = al, rate = th),
           q = function(u) s + 1 / sqrt(stats::qgamma(1 - u, shape = al, rate = th)))
    },
    "inv_gamma2" = {
      ## X - s = 1/G, G ~ Gamma(shape, rate = scale). The rate is the IG
      ## scale, NOT 1/scale (inverting it once made draws ~100x too large).
      s <- .prior_gen_ab(p3, NA_real_)[1]
      ps <- .ig2_params(p1 - s, p2); sh <- ps$shape; sc <- ps$scale
      if (!(is.finite(sh) && sh > 0 && is.finite(sc) && sc > 0))
        bad_shape("the implied inverse-gamma-2 parameters are not positive")
      list(lo = s, hi = Inf,
           r = function(n) s + 1 / rgamma(n, shape = sh, rate = sc),
           p = function(x) 1 - stats::pgamma(1 / (x - s), shape = sh, rate = sc),
           q = function(u) s + 1 / stats::qgamma(1 - u, shape = sh, rate = sc))
    },
    "uniform" = {
      ## Dynare convention: p1/p2 = mean/sd, or the support [p3, p4].
      ab <- .uniform_ab(p1, p2, p3, p4); a <- ab[1]; b <- ab[2]
      if (!(is.finite(a) && is.finite(b) && a < b))
        bad_shape(paste0("a uniform prior needs a finite support a < b (mean ",
                         "p1 and sd p2 > 0, or bounds p3 < p4); got [", a,
                         ", ", b, "]"))
      list(lo = a, hi = b,
           r = function(n) runif(n, a, b),
           p = function(x) (x - a) / (b - a),
           q = function(u) a + u * (b - a))
    },
    .dynhr_abort(
      "prior for \"", label, "\": unsupported prior distribution \"", d,
      "\". Supported: beta, gamma, normal, inv_gamma (= inv_gamma1), ",
      "inv_gamma2, uniform. (The pre-0.9.4 prior sampler drew an unrecognised ",
      "name -- including Dynare spellings such as \"beta_pdf\" -- from ",
      "N(p1, p2).)",
      class = "dynhr_error_unsupported_prior")
  )
}

#' Draw \code{n} values from EXACTLY the prior \code{.lp_dist1()} scores
#'
#' The ONE prior sampler. SMC, SBC, the prior predictive,
#' prior sensitivity, DIME, SMC^2 and profile-CI seeds all reach it through
#' \code{.smc_make_prior_sampler()}. It draws from the law of \code{.lp_dist1(x, dist,
#' p1, p2, p3, p4)} restricted to \code{[lower, upper]} -- the distribution
#' \code{log_prior()} scores (it rejects outside \code{[lower, upper]} and does not
#' renormalise, which is a constant):
#'
#'  * no effective truncation (\code{[lower, upper]} contains the natural support):
#'    a direct draw;
#'  * otherwise inverse-CDF on the truncated interval, \code{Q(U(F(lo), F(hi)))},
#'    exact for every continuous F -- never clamping (a point mass on the
#'    bound) and never a rescaled "generalised" beta standing in for a
#'    truncated one (the pre-fix SMC sampler did that, with a silent
#'    \code{max(v, 2)} shape clamp: sd 0.173 drawn vs 0.158 scored).
#'
#' The name is canonicalised with \code{.normalize_dist()} ("BETA_PDF", "beta_pdf"
#' and "beta" agree) and an unrecognised name is a classed error, never a
#' silent normal.
#'
#' @param n Number of draws.
#' @param dist Distribution name (any spelling \code{.normalize_dist()} accepts).
#' @param p1,p2,p3,p4 Hyperparameters as in the prior spec.
#' @param lower,upper Truncation bounds (NA = none).
#' @param label Parameter name for error messages.
#' @return Numeric vector of length \code{n}.
#' @noRd
.rprior_dist1 <- function(n, dist, p1, p2, p3 = NA_real_, p4 = NA_real_,
                          lower = -Inf, upper = Inf, label = dist) {
  d <- .normalize_dist(dist)
  law <- .prior_law(d, p1, p2, p3, p4, label = label)
  lo <- if (length(lower) == 0L || is.na(lower)) -Inf else lower
  hi <- if (length(upper) == 0L || is.na(upper))  Inf else upper
  cut_lo <- lo > law$lo
  cut_hi <- hi < law$hi
  if (!cut_lo && !cut_hi) return(law$r(n))
  p_lo <- if (cut_lo) law$p(lo) else 0
  p_hi <- if (cut_hi) law$p(hi) else 1
  if (!(p_hi > p_lo))
    .dynhr_abort(
      "prior for \"", label, "\": the bounds [", lo, ", ", hi, "] carry no ",
      "probability under the declared ", d, " distribution, so no draw is ",
      "possible.",
      class = "dynhr_error_empty_prior_support")
  law$q(runif(n, p_lo, p_hi))
}

#' Evaluate log-prior density for a single parameter
#'
#' Supports: "beta", "gamma", "normal", "inv_gamma" (IG1 = Dynare.jl convention),
#' "inv_gamma2" (standard IG), "uniform".
#'
#' @param x   Parameter value
#' @param dist Distribution name (from prior_spec$distribution)
#' @param p1  First hyperparameter (the mean; for uniform the Dynare MEAN)
#' @param p2  Second hyperparameter (the std; for uniform the Dynare STANDARD
#'   DEVIATION, support \code{[p1 - sqrt(3) p2, p1 + sqrt(3) p2]})
#' @param p3  Optional lower bound. \code{NA} means "no lower bound" (-Inf).
#' @param p4  Optional upper bound. \code{NA} means "no upper bound" (Inf).
#' @param gen_p3,gen_p4 Dynare generalisation parameters -- the prior spec's
#'   \code{p3}/\code{p4} COLUMNS (beta support, gamma / inverse-gamma shift; NA =
#'   default). Not to be confused with the truncation bounds \code{p3}/\code{p4} above.
#' @return Scalar log-density (finite or -Inf)
#' @noRd
log_prior_density <- function(x, dist, p1, p2, p3 = -Inf, p4 = Inf,
                              gen_p3 = NA_real_, gen_p4 = NA_real_) {
  ## NA bounds are what `extract_prior_spec()` writes for an unbounded
  ## side. `x < NA` is NA, and `if (NA)` ERRORS -- this function used to abort
  ## on any spec row with an NA bound. Treat NA as -Inf / Inf, exactly as
  ## `log_prior()` already did.
  if (!is.finite(x)) return(-Inf)
  if (!is.na(p3) && x < p3) return(-Inf)
  if (!is.na(p4) && x > p4) return(-Inf)

  .lp_dist1(x, .normalize_dist(dist), p1, p2, gen_p3, gen_p4)
}

#' Compute total log-prior density for a named parameter vector
#'
#' Hyperparameters follow Dynare: \code{p1}/\code{p2} are the prior MEAN and STANDARD
#' DEVIATION for every shape, and the optional columns \code{p3}/\code{p4} are Dynare's
#' PRIOR_3RD/4TH_PARAMETER (generalised-beta support, gamma / inverse-gamma
#' shift). For a UNIFORM prior this means (Dynare \code{uniform_specification()}):
#' the support is \code{[p3, p4]} when both are given, otherwise
#' \code{[p1 - sqrt(3) * p2, p1 + sqrt(3) * p2]} -- \code{p1}/\code{p2} are NOT the bounds.
#' A hand-built uniform on \code{[a, b]} is therefore written with \code{p3 = a,
#' p4 = b} (and \code{p1 = (a + b)/2}, \code{p2 = (b - a)/sqrt(12)}, which is what
#' [extract_prior_spec()] emits). \code{lower}/\code{upper} truncate every prior.
#'
#' @param theta      Named numeric vector (parameter values to evaluate)
#' @param prior_spec Data.frame with columns: name, distribution, p1, p2,
#'                   lower, upper (optionally p3, p4)
#' @return Scalar log-prior (finite, or -Inf if any parameter is out of bounds)
#' @export
log_prior <- function(theta, prior_spec) {
  ## The hyper-parameter conversions (beta shapes, inverse-gamma parameters,
  ## the uniform support, ...) depend only on the prior spec, so they are done
  ## once per spec (.prior_plan) and the density is evaluated per distribution
  ## family, vectorised over the parameters of that family.
  plan <- .prior_plan(prior_spec, names(theta))
  x <- as.numeric(theta)[plan$pos]
  ## Non-finite parameter values (NA/NaN/Inf, e.g. from a diverged HMC
  ## trajectory or a failed gradient) are outside every prior's support:
  ## return -Inf rather than letting `x < lo` evaluate to NA and crash.
  if (!all(is.finite(x))) return(-Inf)
  if (any(x < plan$lower, na.rm = TRUE) || any(x > plan$upper, na.rm = TRUE))
    return(-Inf)
  ll <- numeric(length(x))
  for (g in plan$groups) ll[g$k] <- .lp_group(g, x[g$k])
  if (!all(is.finite(ll))) return(-Inf)
  ## Accumulate in prior-spec order, exactly as the one-parameter-at-a-time
  ## evaluation did (sum() would add in extended precision and differ in the
  ## last bit).
  lp <- 0
  for (v in ll) lp <- lp + v
  lp
}

## --------------------------------------------------------------------------
## Prior plan: per-spec precomputation shared by log_prior() and
## .dlog_prior_grouped().
## --------------------------------------------------------------------------

## Small FIFO cache of plans: a posterior closure calls log_prior() with the
## SAME prior_spec and parameter names on every evaluation, so the plan is
## built once. A hit costs one identical() per slot (pointer-equal when the
## closure holds the spec object).
.prior_plan_cache <- new.env(parent = emptyenv())
.prior_plan_cache$slots <- list()

.prior_plan <- function(prior_spec, theta_names) {
  slots <- .prior_plan_cache$slots
  for (s in slots)
    if (identical(s$names, theta_names) && identical(s$spec, prior_spec))
      return(s$plan)
  plan <- .prior_plan_build(prior_spec, theta_names)
  slots <- c(list(list(spec = prior_spec, names = theta_names, plan = plan)),
             slots)
  .prior_plan_cache$slots <- slots[seq_len(min(length(slots), 8L))]
  plan
}

## Rows of the spec that name an entry of theta (a row that does not is
## skipped, as before), grouped by canonical distribution. `pos` is the
## position of each kept row in theta; `lower`/`upper` the truncation bounds
## (NA = none).
.prior_plan_build <- function(prior_spec, theta_names) {
  n_spec <- length(prior_spec$name)
  pos_all <- match(prior_spec$name, theta_names)
  rows <- which(!is.na(pos_all))
  dist <- .normalize_dist(prior_spec$distribution)[rows]
  p1 <- prior_spec$p1[rows]
  p2 <- prior_spec$p2[rows]
  ## Dynare generalisation parameters (beta support / gamma + IG shift); a
  ## hand-built spec without the columns gets the ungeneralised laws.
  p3 <- (prior_spec$p3 %||% rep(NA_real_, n_spec))[rows]
  p4 <- (prior_spec$p4 %||% rep(NA_real_, n_spec))[rows]
  groups <- list()
  for (d in unique(dist)) {
    k <- which(dist == d)
    groups[[d]] <- .prior_group(d, k, p1[k], p2[k], p3[k], p4[k])
  }
  list(pos = pos_all[rows], lower = prior_spec$lower[rows],
       upper = prior_spec$upper[rows], groups = groups)
}

## Precomputed constants of one distribution family over its rows `k`.
.prior_group <- function(d, k, p1, p2, p3, p4) {
  s <- ifelse(is.na(p3), 0, p3)
  switch(d,
    "normal" = list(d = d, k = k, p1 = p1, p2 = p2),
    "beta" = {
      a <- ifelse(is.na(p3), 0, p3)
      b <- ifelse(is.na(p4), 1, p4)
      len <- b - a
      sh1 <- sh2 <- numeric(length(k))
      for (i in seq_along(k)) {
        sh <- .beta_shapes(p1[i], p2[i], a[i], b[i])
        sh1[i] <- sh[1]; sh2[i] <- sh[2]
      }
      ## A degenerate implied shape is an INVALID beta, not a flat prior.
      valid <- is.finite(sh1) & is.finite(sh2) & sh1 > 0 & sh2 > 0
      list(d = d, k = k, a = a, b = b, len = len, loglen = log(len),
           sh1 = sh1, sh2 = sh2, valid = valid)
    },
    "gamma" = {
      m <- p1 - s
      list(d = d, k = k, s = s, shape = (m / p2)^2, rate = m / p2^2)
    },
    "inv_gamma" =, "inv_gamma1" = {
      alpha <- theta <- numeric(length(k))
      for (i in seq_along(k)) {
        ps <- .ig1_params(p1[i] - s[i], p2[i])
        alpha[i] <- ps$alpha; theta[i] <- ps$theta
      }
      ## c0 is the leading part of the IG1 log density, summed in the same
      ## left-to-right order as .lp_ig1().
      list(d = "inv_gamma", k = k, s = s, alpha = alpha, theta = theta,
           c0 = log(2) + alpha * log(theta) - lgamma(alpha),
           a2 = 2 * alpha + 1, th2 = 2 * theta, ok = !is.na(alpha))
    },
    "inv_gamma2" = {
      shape <- scale <- numeric(length(k))
      for (i in seq_along(k)) {
        ps <- .ig2_params(p1[i] - s[i], p2[i])
        shape[i] <- ps$shape; scale[i] <- ps$scale
      }
      list(d = d, k = k, s = s, shape = shape, scale = scale,
           c0 = shape * log(scale) - lgamma(shape), sh1 = shape + 1)
    },
    "uniform" = {
      lo <- hi <- numeric(length(k))
      for (i in seq_along(k)) {
        ab <- .uniform_ab(p1[i], p2[i], p3[i], p4[i])
        lo[i] <- ab[1]; hi[i] <- ab[2]
      }
      valid <- !is.na(lo < hi) & lo < hi & is.finite(hi - lo)
      list(d = d, k = k, lo = lo, hi = hi, valid = valid,
           nlen = -log(hi - lo))
    },
    ## Fail loud on an unrecognised distribution (a misspelled name would
    ## otherwise be given a silent flat prior): .lp_dist1() raises the error.
    .lp_dist1(0, d, 1, 1)
  )
}

## Log density of one family's parameters at x (same length as g$k).
.lp_group <- function(g, x) {
  switch(g$d,
    "normal" = dnorm(x, mean = g$p1, sd = g$p2, log = TRUE),
    "beta" = {
      ll <- rep(-Inf, length(x))
      ok <- g$valid & x > g$a & x < g$b
      if (any(ok))
        ll[ok] <- dbeta(((x - g$a) / g$len)[ok], g$sh1[ok], g$sh2[ok],
                        log = TRUE) - g$loglen[ok]
      ll
    },
    "gamma" = {
      y <- x - g$s
      ll <- rep(-Inf, length(x))
      ok <- y > 0
      if (any(ok))
        ll[ok] <- dgamma(y[ok], shape = g$shape[ok], rate = g$rate[ok],
                         log = TRUE)
      ll
    },
    "inv_gamma" = {
      y <- x - g$s
      ll <- rep(-Inf, length(x))
      ok <- y > 0 & g$ok
      if (any(ok))
        ll[ok] <- g$c0[ok] - g$a2[ok] * log(y[ok]) - g$theta[ok] / y[ok]^2
      ll
    },
    "inv_gamma2" = {
      y <- x - g$s
      ll <- rep(-Inf, length(x))
      ok <- y > 0
      if (any(ok))
        ll[ok] <- g$c0[ok] - g$sh1[ok] * log(y[ok]) - g$scale[ok] / y[ok]
      ll
    },
    "uniform" = ifelse(g$valid & x >= g$lo & x <= g$hi, g$nlen, -Inf)
  )
}

#' Analytic gradient of the total log-prior, grouped by distribution
#'
#' Same value as looping \code{.dlog_prior_density1()} over the spec rows
#' (zero for a parameter not in the spec, zero outside a prior's support),
#' with the hyper-parameter conversions taken from the cached prior plan and
#' the arithmetic vectorised over each distribution family.
#' @noRd
.dlog_prior_grouped <- function(theta, prior_spec) {
  g <- numeric(length(theta))
  names(g) <- names(theta)
  plan <- .prior_plan(prior_spec, names(theta))
  x <- as.numeric(theta)[plan$pos]
  d <- numeric(length(x))
  for (gr in plan$groups) d[gr$k] <- .dlp_group(gr, x[gr$k])
  g[plan$pos] <- d
  g
}

.dlp_group <- function(g, x) {
  out <- numeric(length(x))
  switch(g$d,
    "normal" = -(x - g$p1) / g$p2^2,
    "uniform" = out,
    "beta" = {
      ok <- g$valid & !(x <= g$a | x >= g$b)
      if (any(ok)) {
        y <- ((x - g$a) / g$len)[ok]
        out[ok] <- (((g$sh1[ok] - 1) / y - (g$sh2[ok] - 1) / (1 - y)) /
                      g$len[ok])
      }
      out
    },
    "gamma" = {
      y <- x - g$s
      ok <- !(y <= 0)
      if (any(ok))
        out[ok] <- (g$shape[ok] - 1) / y[ok] - g$rate[ok]
      out
    },
    "inv_gamma" = {
      y <- x - g$s
      ok <- !(y <= 0)
      if (any(ok))
        out[ok] <- -g$a2[ok] / y[ok] + g$th2[ok] / y[ok]^3
      out
    },
    "inv_gamma2" = {
      y <- x - g$s
      ok <- !(y <= 0)
      if (any(ok))
        out[ok] <- -g$sh1[ok] / y[ok] + g$scale[ok] / y[ok]^2
      out
    }
  )
}
