## R/system-priors.R
## --------------------------------------------------------------------------
## System priors (IRIS terminology; also "prior predictive restrictions").
##
## A system prior places a prior not on a *parameter* but on a *model feature*
## such as an impulse response, a variance share, or a spectral peak. Each
## feature is a scalar (or short vector) function of the solved model evaluated
## at a particular parameter draw; the prior is expressed as a log-density
## (possibly improper, possibly a hard indicator restriction).
##
## At each posterior evaluation the feature is computed once from the already-
## solved dr; there is NO re-solve. The accumulated log-density is added to the
## log-prior, so the posterior is:
##
##   log p(theta | data) = log p(data | theta) + log p(theta) + sum_k log f_k(feature_k(theta))
##
## where f_k is the k-th system prior density and feature_k operates on the
## state list(theta, model, dr, Sigma_e, params).
##
## References:
##   Andrle & Benes (2013). System priors: formulating priors about DSGE models'
##     properties. IMF Working Paper WP/13/257.
##   Irsova & Havranek (2015). Capital-Labor Substitution. Journal of Economic
##     Surveys (for practical variance-share restrictions).
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1.  Constructor: system_prior_spec()
## ---------------------------------------------------------------------------

#' Specify system priors for a DSGE model
#'
#' Constructs a list of system-prior entries. Each entry defines a \emph{feature
#' function} that maps the solved model state to a numeric scalar (or short
#' vector) and a density (or predicate) that acts as a soft or hard restriction
#' on that feature.
#'
#' A \emph{soft restriction} adds a log-density penalty: the feature is
#' evaluated, the log-density at that value is accumulated, and draws for which
#' the penalty is \code{-Inf} (e.g. outside a uniform prior's support) are
#' rejected.
#'
#' A \emph{hard restriction} is expressed as a predicate: if the predicate
#' returns \code{FALSE}, the entire draw gets log-posterior \code{-Inf}
#' (rejected) with no error. Hard restrictions are equivalent to an improper
#' uniform that is 0 inside the valid region and \code{-Inf} outside.
#'
#' The feature function receives a \code{state} list with fields:
#' \itemize{
#'   \item \code{theta}   -- named parameter vector (structural parameters)
#'   \item \code{model}   -- parsed \code{dynhr_mod}
#'   \item \code{dr}      -- solved decision rules (\code{ghx}, \code{ghu}, etc.)
#'   \item \code{Sigma_e} -- shock covariance matrix (n_exo x n_exo)
#'   \item \code{params}  -- full expanded named parameter vector
#' }
#'
#' Feature functions must be self-contained closures that do not reference
#' large external objects (multi-GB data frames, etc.) when the model is
#' estimated in parallel via \code{mirai}: they will be serialised to worker
#' daemons.
#'
#' @param ...  One or more system-prior entries created by
#'   \code{sp_entry()}, \code{sp_irf()}, \code{sp_variance_share()}, or
#'   \code{sp_custom()}.
#'
#' @return An object of class \code{"system_prior_spec"}: a named list of
#'   entries, each a list with \code{$feature} (function) and \code{$density}
#'   (list or function) and an optional \code{$label} (character).
#'
#' @examples
#' \dontrun{
#' ## Sign restriction on the IRF of output to a technology shock at horizon 1
#' sp <- system_prior_spec(
#'   sp_irf("y", "eps_a", horizon = 1L,
#'           density = list(type = "hard_sign", sign = "positive"))
#' )
#'
#' ## Soft Gaussian prior on the variance share of demand shocks in output
#' sp <- system_prior_spec(
#'   sp_variance_share("y", "eps_d",
#'     density = list(dist = "normal", p1 = 0.30, p2 = 0.10))
#' )
#' }
#'
#' @seealso \code{\link{sp_irf}}, \code{\link{sp_variance_share}},
#'   \code{\link{sp_custom}}, \code{make_log_posterior}
#' @export
system_prior_spec <- function(...) {
  entries <- list(...)
  ## Flatten if the user passed a list-of-entries rather than bare entries.
  ## Only flatten when the single argument looks like a list-of-entries (its
  ## first element is itself a list with a $feature key), not when it is a
  ## single entry passed directly.
  if (length(entries) == 1L && is.list(entries[[1L]]) &&
      !inherits(entries[[1L]], "sp_entry") &&
      !is.function(entries[[1L]]$feature))
    entries <- entries[[1L]]

  for (i in seq_along(entries)) {
    e <- entries[[i]]
    if (!is.list(e) || !is.function(e$feature))
      stop("system_prior_spec: entry ", i,
           " must have a $feature function. ",
           "Use sp_irf(), sp_variance_share(), or sp_custom().",
           call. = FALSE)
    if (is.null(e$density))
      stop("system_prior_spec: entry ", i, " has no $density.", call. = FALSE)
  }

  structure(entries, class = "system_prior_spec")
}


## ---------------------------------------------------------------------------
## 2.  Low-level entry builder
## ---------------------------------------------------------------------------

#' Build a single system-prior entry
#'
#' Internal constructor used by the convenience wrappers.  Not usually called
#' directly by users.
#'
#' @param feature  A function \code{function(state)} returning a numeric scalar
#'   (non-finite values yield \code{-Inf} for the whole draw).
#' @param density  Either:
#'   \itemize{
#'     \item A \emph{named list} with \code{$dist} (distribution name accepted
#'       by \code{\link{log_prior_density}}: "normal", "beta", "gamma",
#'       "inv_gamma", "uniform") and parameters \code{$p1}, \code{$p2}
#'       (optionally \code{$p3}, \code{$p4} for bounds). Or:
#'     \item \code{list(type = "hard_sign", sign = "positive")} /
#'       \code{"negative"} for a hard sign restriction (indicator
#'       \eqn{1\{x > 0\}}). Or:
#'     \item \code{list(type = "hard_predicate", predicate = function(x)
#'       logical(1))} for an arbitrary hard restriction. Or:
#'     \item A function \code{function(x) scalar log-density}.
#'   }
#' @param label  Optional character label for diagnostics.
#'
#' @return A list of class \code{"sp_entry"}.
#' @noRd
.sp_entry <- function(feature, density, label = NULL) {
  structure(
    list(feature = feature, density = density, label = label),
    class = "sp_entry"
  )
}


## ---------------------------------------------------------------------------
## 3.  Convenience constructors
## ---------------------------------------------------------------------------

#' System prior on an impulse response
#'
#' Creates a system-prior entry based on the impulse response of variable
#' \code{var} to shock \code{shock} at horizon \code{horizon}.  The feature
#' value is the IRF coefficient at that horizon (in model units, i.e. deviation
#' from steady state per one standard-deviation shock).
#'
#' The computation reuses \code{compute_irfs(dr, model, n_periods, params)}
#' from \code{R/stochsimul-monolith.R} — no re-solve.
#'
#' @param var      Endogenous variable name (character).
#' @param shock    Exogenous shock name (character).
#' @param horizon  Response horizon (integer >= 1).
#' @param density  Density specification (see \code{.sp_entry}).
#'   Common choices:
#'   \itemize{
#'     \item \code{list(type = "hard_sign", sign = "positive")} -- IRF > 0.
#'     \item \code{list(dist = "normal", p1 = 0.5, p2 = 0.2)} -- soft N(0.5,
#'       0.2²) prior on the IRF coefficient.
#'   }
#' @param n_periods Total IRF periods to compute (must be >= \code{horizon};
#'   default \code{horizon}).
#'
#' @return A list of class \code{"sp_entry"}.
#' @seealso \code{\link{system_prior_spec}}, \code{\link{sp_variance_share}}
#' @export
sp_irf <- function(var, shock, horizon, density, n_periods = NULL) {
  stopifnot(is.character(var),   length(var) == 1L)
  stopifnot(is.character(shock), length(shock) == 1L)
  horizon <- as.integer(horizon)
  stopifnot(horizon >= 1L)
  if (is.null(n_periods)) n_periods <- horizon
  n_periods <- max(as.integer(n_periods), horizon)

  ## Capture by value (not by reference to caller's frame)
  force(var); force(shock); force(horizon); force(n_periods)

  feature <- function(state) {
    dr     <- state$dr
    model  <- state$model
    params <- state$params
    irfs   <- compute_irfs(dr, model, n_periods = n_periods, params = params)
    mat    <- irfs[[shock]]
    if (is.null(mat))
      return(NA_real_)
    if (!(var %in% colnames(mat)))
      return(NA_real_)
    mat[horizon, var]
  }

  .sp_entry(feature, density,
            label = sprintf("irf(%s, %s, h=%d)", var, shock, horizon))
}


#' System prior on a variance share
#'
#' Creates a system-prior entry based on the share of the unconditional variance
#' of variable \code{var} that is attributable to shock \code{shock}.  The
#' feature value is in \eqn{[0, 1]} (fraction, NOT percentage).
#'
#' The computation reuses \code{compute_moments(dr, model, n_ar = 0L, params)}
#' from \code{R/stochsimul-monolith.R} — no re-solve, one Lyapunov solve per
#' shock.
#'
#' @param var    Endogenous variable name (character).
#' @param shock  Exogenous shock name (character).
#' @param density  Density specification (see \code{.sp_entry}).
#'   Common choices:
#'   \itemize{
#'     \item \code{list(dist = "beta", p1 = 0.50, p2 = 0.20)} -- soft prior
#'       that the shock explains about half of the variance.
#'     \item \code{list(dist = "uniform", p1 = 0, p2 = 1)} -- flat prior
#'       (effectively no restriction); mainly useful as a placeholder.
#'   }
#'
#' @return A list of class \code{"sp_entry"}.
#' @seealso \code{\link{system_prior_spec}}, \code{\link{sp_irf}}
#' @export
sp_variance_share <- function(var, shock, density) {
  stopifnot(is.character(var),   length(var) == 1L)
  stopifnot(is.character(shock), length(shock) == 1L)

  force(var); force(shock)

  feature <- function(state) {
    dr     <- state$dr
    model  <- state$model
    params <- state$params
    ## n_ar = 1L: compute_moments requires at least 1 lag (n_ar=0 triggers a
    ## 1:0 dimension bug in the autocorr array). We don't need autocorrelations;
    ## only var_decomp_pct is used.
    mom    <- compute_moments(dr, model, n_ar = 1L, params = params)
    ## var_decomp_pct is percentage; convert to fraction [0,1]
    pct_mat <- mom$var_decomp_pct
    if (is.null(pct_mat) || !(var %in% rownames(pct_mat)) ||
        !(shock %in% colnames(pct_mat)))
      return(NA_real_)
    pct_mat[var, shock] / 100
  }

  .sp_entry(feature, density,
            label = sprintf("vardec(%s, %s)", var, shock))
}


#' System prior on the frequency of a spectral peak
#'
#' Creates a system-prior entry whose feature value is the angular frequency
#' \eqn{\omega^* \in (0, \pi)} (radians per period) at which the marginal
#' spectral density of variable \code{var} is maximised.  A flat spectrum
#' (typical of white-noise processes) peaks at \eqn{\omega = 0}; a business-
#' cycle component peaks somewhere around \eqn{2\pi/32 \approx 0.196}
#' radians.
#'
#' \strong{Convention.}  The spectral density is computed via
#' \code{.whittle_spectral_density(omega, TT, RR, ZZ, DD, Sigma_e)} from
#' \code{R/whittle-likelihood.R}, using the lagged-state observation
#' convention \eqn{H = ZZ z (I - TT z)^{-1} RR + DD} where
#' \eqn{z = e^{-i\omega}}.  \code{TT}, \code{RR} are the state-to-state
#' and shock-to-state rows; \code{ZZ}, \code{DD} are the single row of
#' \code{dr$ghx} and \code{dr$ghu} corresponding to \code{var}.
#' Measurement error is ignored (\code{me_variance = 0}) since the feature
#' targets the model's theoretical spectral peak, not the observed data.
#'
#' \strong{Grid.}  The search uses \code{n_grid} equally-spaced frequencies in
#' \eqn{(0, \pi)}.  The default (200 points) resolves peaks to within
#' \eqn{\pi/200 \approx 0.016} radians, sufficient for identifying business-
#' cycle vs. lower-frequency peaks. Increase \code{n_grid} for finer
#' resolution if needed.
#'
#' @param var      Endogenous variable name (character).
#' @param density  Density specification (see \code{.sp_entry}).
#'   Typically a \code{list(dist = "normal", p1 = omega_target, p2 = sigma)}
#'   where \code{omega_target} is the desired peak frequency in radians, e.g.
#'   \code{2*pi/32} for a 32-period (business-cycle) peak.  Hard restrictions
#'   are also supported.
#' @param n_grid   Number of grid points in \eqn{(0, \pi)} for the peak search.
#'   Default 200L.
#'
#' @return A list of class \code{"sp_entry"}.
#' @seealso \code{\link{system_prior_spec}}, \code{\link{sp_irf}},
#'   \code{\link{sp_variance_share}}
#' @export
sp_spectral_peak <- function(var, density, n_grid = 200L) {
  stopifnot(is.character(var), length(var) == 1L)
  n_grid <- as.integer(n_grid)
  stopifnot(n_grid >= 10L)

  force(var); force(n_grid)

  feature <- function(state) {
    dr      <- state$dr
    Sigma_e <- state$Sigma_e

    endo      <- dr$endo_names
    state_idx <- dr$state_idx
    var_idx   <- match(var, endo)
    if (is.na(var_idx)) return(NA_real_)

    ghx <- dr$ghx
    ghu <- dr$ghu

    TT <- ghx[state_idx,            , drop = FALSE]   # n_state x n_state
    RR <- ghu[state_idx,            , drop = FALSE]   # n_state x n_exo
    ZZ <- ghx[var_idx,              , drop = FALSE]   # 1 x n_state
    DD <- ghu[var_idx,              , drop = FALSE]   # 1 x n_exo

    ## Grid over (0, pi), endpoints excluded to avoid DC and Nyquist artefacts.
    omegas <- seq(pi / (n_grid + 1L), pi * n_grid / (n_grid + 1L),
                  length.out = n_grid)

    ## Marginal spectral density of var: real part of the (1,1) element of S(omega)
    psd <- vapply(omegas, function(omega) {
      S <- dynhr:::.whittle_spectral_density(omega, TT, RR, ZZ, DD, Sigma_e,
                                              me_variance = 0)
      Re(S[1L, 1L])
    }, numeric(1L))

    ## Return the frequency of the maximum
    omegas[which.max(psd)]
  }

  .sp_entry(feature, density,
            label = sprintf("spectral_peak(%s)", var))
}


#' Custom system prior
#'
#' Most general entry point: supply any feature function and any density.
#'
#' @param feature  A function \code{function(state)} returning a numeric scalar.
#' @param density  A density specification (list or function; see
#'   \code{.sp_entry}).
#' @param label    Optional character label.
#'
#' @return A list of class \code{"sp_entry"}.
#' @seealso \code{\link{system_prior_spec}}
#' @export
sp_custom <- function(feature, density, label = NULL) {
  if (!is.function(feature))
    stop("sp_custom: `feature` must be a function.", call. = FALSE)
  .sp_entry(feature, density, label = label)
}


## ---------------------------------------------------------------------------
## 4.  Density evaluator
## ---------------------------------------------------------------------------

#' Evaluate a single density specification at a feature value
#'
#' @param density  A density specification as described in \code{.sp_entry}.
#' @param x        The feature value (numeric scalar).
#' @return Scalar log-density (finite or -Inf).
#' @noRd
.eval_sp_density <- function(density, x) {
  ## Hard sign restriction
  if (is.list(density) && identical(density$type, "hard_sign")) {
    sign_req <- density$sign %||% "positive"
    if (identical(sign_req, "positive") || identical(sign_req, "+"))
      return(if (x > 0) 0 else -Inf)
    if (identical(sign_req, "negative") || identical(sign_req, "-"))
      return(if (x < 0) 0 else -Inf)
    ## zero: exact zero is unlikely, so we treat it as feasible
    return(if (x == 0) 0 else -Inf)
  }

  ## Hard predicate restriction
  if (is.list(density) && identical(density$type, "hard_predicate")) {
    pred <- density$predicate
    if (!is.function(pred))
      stop(".eval_sp_density: hard_predicate$predicate must be a function.",
           call. = FALSE)
    result <- tryCatch(pred(x), error = function(e) FALSE)
    return(if (isTRUE(result)) 0 else -Inf)
  }

  ## Arbitrary log-density function
  if (is.function(density)) {
    val <- tryCatch(density(x), error = function(e) -Inf)
    ## +Inf log-density (e.g. degenerate Dirac-like) is treated as 0 penalty
    ## to avoid contaminating the posterior sum.
    return(if (is.nan(val) || is.na(val)) -Inf
           else if (val == Inf)  0
           else val)
  }

  ## Named-distribution list: {dist, p1, p2 [, p3, p4]}
  if (is.list(density) && !is.null(density$dist)) {
    p3 <- density$p3 %||% -Inf
    p4 <- density$p4 %||%  Inf
    val <- log_prior_density(x, density$dist, density$p1, density$p2, p3, p4)
    ## +Inf log-density (boundary of support for beta/uniform) is clamped to 0
    ## so the accumulated sum stays finite.
    if (isTRUE(val == Inf)) return(0)
    return(val)
  }

  warning(".eval_sp_density: unrecognised density specification; returning 0.",
          call. = FALSE)
  0
}


## ---------------------------------------------------------------------------
## 5.  System-prior evaluator (called inside log-posterior closures)
## ---------------------------------------------------------------------------

#' Evaluate all system priors and return their sum
#'
#' Called once per posterior evaluation, after the model has been solved.
#'
#' @param spec   A \code{system_prior_spec} (or \code{NULL}).
#' @param state  A named list with fields \code{theta}, \code{model},
#'   \code{dr}, \code{Sigma_e}, \code{params}.
#'
#' @return Scalar: sum of all system-prior log-densities. Returns \code{0} if
#'   \code{spec} is \code{NULL} or empty. Returns \code{-Inf} on the first
#'   non-finite feature value or density value without raising an error.
#' @noRd
.eval_system_priors <- function(spec, state) {
  if (is.null(spec) || length(spec) == 0L) return(0)

  lp_sum <- 0
  for (entry in spec) {
    ## Compute the feature, guarding against errors in user code
    fval <- tryCatch(
      entry$feature(state),
      error = function(e) {
        warning(".eval_system_priors: feature '",
                entry$label %||% "?", "' threw: ", conditionMessage(e),
                call. = FALSE)
        NA_real_
      }
    )

    ## Non-finite feature value: reject the draw cleanly
    if (!is.finite(fval)) return(-Inf)

    ## Evaluate density at feature value. +Inf is clamped to 0 inside
    ## .eval_sp_density; NA or -Inf short-circuits the whole draw.
    lp_k <- .eval_sp_density(entry$density, fval)
    if (is.na(lp_k) || lp_k == -Inf) return(-Inf)
    lp_sum <- lp_sum + lp_k
  }
  lp_sum
}


## ---------------------------------------------------------------------------
## 6.  Print method for diagnostics
## ---------------------------------------------------------------------------

#' @export
print.system_prior_spec <- function(x, ...) {
  n <- length(x)
  cat(sprintf("<system_prior_spec>  (%d %s)\n", n,
              if (n == 1L) "entry" else "entries"))
  for (i in seq_along(x)) {
    lbl <- x[[i]]$label %||% sprintf("entry %d", i)
    dns <- x[[i]]$density
    dns_str <- if (is.function(dns)) {
      "custom log-density fn"
    } else if (is.list(dns) && !is.null(dns$type)) {
      paste0("hard (", dns$type, ")")
    } else if (is.list(dns) && !is.null(dns$dist)) {
      paste0(dns$dist, "(", dns$p1, ", ", dns$p2, ")")
    } else {
      "?"
    }
    cat(sprintf("  [%d] %s  ->  %s\n", i, lbl, dns_str))
  }
  invisible(x)
}
