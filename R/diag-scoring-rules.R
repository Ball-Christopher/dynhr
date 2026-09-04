## R/diag-scoring-rules.R
## --------------------------------------------------------------------------
## Proper scoring rules for posterior-predictive forecast evaluation.
##
## Implements ensemble (sample-based) estimators for:
##   - CRPS  (Continuous Ranked Probability Score, univariate)
##   - Energy Score  (multivariate generalisation of CRPS)
##   - Variogram Score  (captures dependence structure)
##   - Logarithmic score (Gaussian plug-in from the draw mean/cov, or KDE)
##   - PIT values (probability integral transform, calibration diagnostic)
##   - Central-interval coverage at user-given levels
##
## plus the Diebold-Mariano test for two score series (dm_test()).
##
## ORIENTATION.  Every score here is NEGATIVELY oriented: LOWER IS BETTER.
## That includes the logarithmic score, which is reported as -log f(y) (the
## "ignorance score"), so that it can be averaged and differenced alongside
## CRPS/energy/variogram without a sign flip.
##
## References:
##   Gneiting & Raftery (2007) "Strictly proper scoring rules..."
##     JASA 102(477):359-378.
##   Scheuerer & Hamill (2015) "Variogram-based proper scoring rules"
##     Monthly Weather Review 143(4):1321-1334.
##   Diebold & Mariano (1995) "Comparing predictive accuracy"
##     JBES 13(3):253-263.
##   Harvey, Leybourne & Newbold (1997) "Testing the equality of prediction
##     mean squared errors" IJF 13(2):281-291.
##   Dawid (1984) / Diebold, Gunther & Tay (1998) -- PIT calibration.
##
## Public API: score_forecast(), dm_test()
## Internal helpers: .crps_ensemble(), .energy_score(), .variogram_score(),
##   .crps_gaussian(), .logs_gaussian_uni(), .logs_gaussian_joint(),
##   .logs_kde(), .pit_ensemble(), .pit_gaussian(), .coverage_ensemble(),
##   .coverage_gaussian()
## --------------------------------------------------------------------------


## ---- Internal helpers ----------------------------------------------------

#' Ensemble CRPS for a single variable
#'
#' Gneiting-Raftery (2007) eq. 21 (kernel representation):
#'   CRPS(F,y) = E|X - y| - 0.5 * E|X - X'|
#' where X, X' are iid draws from F.
#'
#' @param x Numeric vector of S ensemble draws.
#' @param y Scalar observation.
#' @return Scalar CRPS value (lower = better).
#' @noRd
.crps_single <- function(x, y) {
  S <- length(x)
  if (S < 2L) stop("CRPS requires >= 2 ensemble members")
  ## E|X - y|
  term1 <- mean(abs(x - y))
  ## 0.5 * E|X - X'| = sum of all pairwise |x_i - x_j| / S^2
  ## Efficient O(S log S) via sorted differences:
  xs  <- sort(x)
  idx <- seq_len(S)
  ## sum_{i<j} |x_i - x_j| = sum_i x_i * (2i - S - 1)  (sorted)
  ## Use double arithmetic to avoid integer overflow for large S.
  pw  <- sum(xs * (2.0 * idx - S - 1.0)) / (as.numeric(S) * S)
  term1 - pw
}

#' Ensemble CRPS for one or more variables
#'
#' @param X  Numeric matrix (S x n_obs) of ensemble draws.
#' @param y  Numeric vector of length n_obs of observations.
#' @return Named numeric vector of length n_obs.
#' @noRd
.crps_ensemble <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n_obs <- ncol(X)
  if (length(y) != n_obs)
    stop(sprintf(".crps_ensemble: ncol(X)=%d but length(y)=%d", n_obs, length(y)))
  vapply(seq_len(n_obs), function(j) .crps_single(X[, j], y[j]), numeric(1L))
}

#' Ensemble Energy Score (multivariate proper scoring rule)
#'
#' ES(F,y) = E||X - y|| - 0.5 * E||X - X'||
#' where ||.|| is the Euclidean norm over n_obs dimensions.
#' Reduces to CRPS when n_obs = 1.
#'
#' @param X  Numeric matrix (S x n_obs) of ensemble draws.
#' @param y  Numeric vector of length n_obs.
#' @return Scalar energy score.
#' @noRd
.energy_score <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  S     <- nrow(X)
  n_obs <- ncol(X)
  if (length(y) != n_obs)
    stop(sprintf(".energy_score: ncol(X)=%d but length(y)=%d", n_obs, length(y)))
  if (S < 2L) stop("Energy score requires >= 2 ensemble members")

  ## E||X_i - y||
  diffs1 <- sweep(X, 2L, y, `-`)
  term1  <- mean(sqrt(rowSums(diffs1^2)))

  ## 0.5 * E||X_i - X_j|| = (sum of pairwise Euclidean distances) / S^2.
  ## stats::dist() does the O(S^2) work in C (direct differences, so no
  ## Gram-matrix cancellation) and stores only the lower triangle -- far faster
  ## than the former R-level i-loop and lighter than a full S x S matrix.
  pw_sum <- sum(stats::dist(X))
  term2  <- pw_sum / (as.numeric(S) * S)  ## normalise by S^2; double avoids integer overflow

  term1 - term2
}

#' Ensemble Variogram Score
#'
#' VS_p(F,y) = sum_{k<l} w_kl * (|y_k - y_l|^p - mean_i |X_{i,k} - X_{i,l}|^p)^2
#'
#' Scheuerer & Hamill (2015).  Default weights w_kl = 1.
#' Default p = 0.5 (recommended for forecast evaluation).
#'
#' @param X  Numeric matrix (S x n_obs) of ensemble draws.
#' @param y  Numeric vector of length n_obs.
#' @param p  Order; default 0.5.
#' @param w  Weight matrix (n_obs x n_obs) or NULL for uniform weights.
#' @return Scalar variogram score.
#' @noRd
.variogram_score <- function(X, y, p = 0.5, w = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n_obs <- ncol(X)
  S     <- nrow(X)
  if (length(y) != n_obs)
    stop(sprintf(".variogram_score: ncol(X)=%d but length(y)=%d", n_obs, length(y)))
  if (n_obs < 2L) stop("Variogram score requires >= 2 observables")
  if (S < 2L) stop("Variogram score requires >= 2 ensemble members")
  if (is.null(w)) w <- matrix(1.0, n_obs, n_obs)

  total <- 0
  for (k in seq_len(n_obs - 1L)) {
    for (l in seq(k + 1L, n_obs)) {
      obs_diff  <- abs(y[k] - y[l])^p
      ens_diffs <- abs(X[, k] - X[, l])^p
      ens_mean  <- mean(ens_diffs)
      total     <- total + w[k, l] * (obs_diff - ens_mean)^2
    }
  }
  total
}


## ---- Closed-form Gaussian scores ----------------------------------------

#' Closed-form CRPS of a Gaussian predictive
#'
#' Gneiting & Raftery (2007) eq. 22:
#'   CRPS(N(mu, sigma^2), y)
#'     = sigma * ( w (2 Phi(w) - 1) + 2 phi(w) - 1/sqrt(pi) ),  w = (y-mu)/sigma
#'
#' Vectorised over `mu`, `sigma`, `y` (recycled by the usual R rules).
#'
#' @param mu    Numeric vector of predictive means.
#' @param sigma Numeric vector of predictive standard deviations (> 0).
#' @param y     Numeric vector of realisations.
#' @return Numeric vector of CRPS values.
#' @noRd
.crps_gaussian <- function(mu, sigma, y) {
  if (any(!is.finite(sigma)) || any(sigma <= 0))
    stop(".crps_gaussian: sigma must be finite and strictly positive",
         call. = FALSE)
  w <- (y - mu) / sigma
  sigma * (w * (2 * stats::pnorm(w) - 1) + 2 * stats::dnorm(w) -
             1 / sqrt(pi))
}

#' Univariate Gaussian logarithmic score (ignorance score, lower = better)
#'
#' @param mu,sigma,y As in \code{.crps_gaussian()}.
#' @return \code{-log dnorm(y, mu, sigma)}, vectorised.
#' @noRd
.logs_gaussian_uni <- function(mu, sigma, y) {
  if (any(!is.finite(sigma)) || any(sigma <= 0))
    stop(".logs_gaussian_uni: sigma must be finite and strictly positive",
         call. = FALSE)
  -stats::dnorm(y, mean = mu, sd = sigma, log = TRUE)
}

#' Joint (multivariate) Gaussian logarithmic score
#'
#' \eqn{-\log N(y; mu, Sigma)} computed from a Cholesky factor, with a
#' ridge-jitter fallback for a numerically non-PD \code{Sigma} (which a small
#' ensemble covariance routinely is).
#'
#' @param mu    Numeric vector (length k) of predictive means.
#' @param Sigma k x k predictive covariance.
#' @param y     Numeric vector (length k) of realisations.
#' @return Scalar joint log score, or \code{NA_real_} if \code{Sigma} cannot
#'   be factorised even after jittering.
#' @noRd
.logs_gaussian_joint <- function(mu, Sigma, y) {
  k     <- length(mu)
  Sigma <- as.matrix(Sigma)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  ch <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(ch)) {
    sc <- mean(diag(Sigma))
    if (!is.finite(sc) || sc <= 0) return(NA_real_)
    for (jit in c(1e-12, 1e-10, 1e-8, 1e-6)) {
      ch <- tryCatch(chol(Sigma + jit * sc * diag(k)), error = function(e) NULL)
      if (!is.null(ch)) break
    }
    if (is.null(ch)) return(NA_real_)
  }
  z       <- backsolve(ch, y - mu, transpose = TRUE)
  log_det <- 2 * sum(log(diag(ch)))
  0.5 * (k * log(2 * pi) + log_det + sum(z^2))
}

#' Kernel-density logarithmic score for one variable
#'
#' Gaussian-kernel density estimate of the predictive built from the ensemble,
#' evaluated at the realisation:
#'   f(y) = (1/S) sum_i dnorm(y; x_i, h)
#' with `h` the Silverman rule-of-thumb bandwidth (\code{stats::bw.nrd0}) by
#' default.  Unlike the Gaussian plug-in it does not assume normality, at the
#' cost of a bandwidth choice.
#'
#' @param x  Numeric vector of S ensemble draws.
#' @param y  Scalar realisation.
#' @param bw Bandwidth, or \code{NULL} for \code{bw.nrd0(x)}.
#' @return Scalar \code{-log f(y)}; \code{Inf} if the KDE underflows.
#' @noRd
.logs_kde <- function(x, y, bw = NULL) {
  x <- as.numeric(x)
  if (is.null(bw)) {
    bw <- tryCatch(stats::bw.nrd0(x), error = function(e) NA_real_)
    ## bw.nrd0 returns 0 for a degenerate ensemble; fall back to a tiny
    ## positive bandwidth so the score is Inf-or-finite, never NaN.
    if (!is.finite(bw) || bw <= 0)
      bw <- max(1e-8, 1e-8 * max(1, abs(mean(x))))
  }
  if (!is.finite(bw) || bw <= 0)
    stop(".logs_kde: bandwidth must be finite and strictly positive",
         call. = FALSE)
  ## log-sum-exp for numerical safety when y is deep in the tail.
  lk <- stats::dnorm(y, mean = x, sd = bw, log = TRUE)
  m  <- max(lk)
  if (!is.finite(m)) return(Inf)
  -(m + log(mean(exp(lk - m))))
}


## ---- Calibration: PIT and interval coverage ------------------------------

#' Randomised PIT of an ensemble forecast
#'
#' The probability integral transform of the realisation under the predictive.
#' With a finite ensemble the empirical CDF is a step function, so the
#' NON-randomised version \code{#\{x < y\}/S} can only take S+1 values and is
#' uniform on a lattice, not on (0,1).  The randomised version
#'
#'   u = ( #\{x_i < y\} + U (1 + #\{x_i = y\}) ) / (S + 1),   U ~ Unif(0,1)
#'
#' is EXACTLY Unif(0,1) whenever y is exchangeable with the S draws -- the
#' same rank-randomisation argument that makes SBC ranks exactly uniform
#' (see \code{.sbc_ranks()} in R/validate-sbc.R).  That exactness is what
#' makes a KS test on
#' the PITs an honest calibration check at moderate S.
#'
#' @param X         S x n_obs matrix of ensemble draws.
#' @param y         Numeric vector of length n_obs.
#' @param randomize Logical; use the randomised rank (default TRUE).
#' @return Numeric vector of length n_obs with values in [0, 1].
#' @noRd
.pit_ensemble <- function(X, y, randomize = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  S <- nrow(X)
  vapply(seq_len(ncol(X)), function(j) {
    xj      <- X[, j]
    n_below <- sum(xj < y[j])
    n_equal <- sum(xj == y[j])
    if (randomize) {
      (n_below + stats::runif(1L) * (1 + n_equal)) / (S + 1)
    } else {
      (n_below + 0.5 * n_equal) / S
    }
  }, numeric(1L))
}

#' PIT of a Gaussian predictive (exact, no ensemble noise)
#' @noRd
.pit_gaussian <- function(mu, sigma, y) stats::pnorm(y, mean = mu, sd = sigma)

#' Central prediction intervals and coverage indicators from an ensemble
#'
#' Equal-tailed intervals: the (1-level)/2 and (1+level)/2 sample quantiles.
#'
#' @param X      S x n_obs matrix of draws.
#' @param y      Numeric vector of length n_obs.
#' @param levels Numeric vector of nominal coverage levels in (0, 1).
#' @return List with \code{lower}, \code{upper}, \code{covered} (all
#'   n_levels x n_obs matrices) and \code{levels}.
#' @noRd
.coverage_ensemble <- function(X, y, levels) {
  X  <- as.matrix(X)
  y  <- as.numeric(y)
  nl <- length(levels)
  nj <- ncol(X)
  lo <- hi <- matrix(NA_real_, nl, nj)
  for (l in seq_len(nl)) {
    p <- c((1 - levels[l]) / 2, (1 + levels[l]) / 2)
    for (j in seq_len(nj)) {
      q <- stats::quantile(X[, j], probs = p, names = FALSE)
      lo[l, j] <- q[1L]
      hi[l, j] <- q[2L]
    }
  }
  .coverage_assemble(lo, hi, y, levels)
}

#' Central prediction intervals and coverage from Gaussian moments
#' @noRd
.coverage_gaussian <- function(mu, sigma, y, levels) {
  z  <- stats::qnorm((1 + levels) / 2)
  lo <- outer(z, sigma, function(a, b) -a * b) +
    matrix(mu, length(levels), length(mu), byrow = TRUE)
  hi <- outer(z, sigma, function(a, b) a * b) +
    matrix(mu, length(levels), length(mu), byrow = TRUE)
  .coverage_assemble(lo, hi, y, levels)
}

#' Shared assembly of the coverage return value
#' @noRd
.coverage_assemble <- function(lo, hi, y, levels) {
  lo <- matrix(lo, nrow = length(levels))
  hi <- matrix(hi, nrow = length(levels))
  yy <- matrix(y, nrow = length(levels), ncol = length(y), byrow = TRUE)
  cov_ind <- (yy >= lo) & (yy <= hi)
  storage.mode(cov_ind) <- "double"
  rn <- sprintf("%g%%", 100 * levels)
  dimnames(lo) <- dimnames(hi) <- dimnames(cov_ind) <- list(rn, NULL)
  list(lower = lo, upper = hi, covered = cov_ind, levels = levels)
}


## ---- Public API ----------------------------------------------------------

#' Proper scoring rules for posterior-predictive forecast evaluation
#'
#' Computes ensemble-based proper scoring rules (CRPS, energy score, variogram
#' score) comparing a set of predictive draws to realised observations.  Lower
#' scores indicate better-calibrated and sharper forecasts.
#'
#' @section Scoring rules:
#' \describe{
#'   \item{crps}{Continuous Ranked Probability Score (Gneiting & Raftery 2007,
#'     eq. 21).  Computed per observable; mean across observables also returned.
#'     Ensemble estimator: \eqn{E|X - y| - 0.5 E|X - X'|}.}
#'   \item{energy}{Energy Score (multivariate CRPS; Gneiting & Raftery 2007).
#'     Single scalar: \eqn{E\|X - y\| - 0.5 E\|X - X'\|} where \eqn{\|\cdot\|}
#'     is the Euclidean norm over all observables.  Equals CRPS when
#'     \code{n_obs = 1}.}
#'   \item{variogram}{Variogram Score of order \code{vs_p} (Scheuerer & Hamill
#'     2015).  Rewards correct dependence structure across observables.
#'     Requires \code{n_obs >= 2}.}
#'   \item{logs}{Logarithmic score, reported as \eqn{-\log f(y)} so that,
#'     like every other rule here, LOWER IS BETTER.  \code{logs_method =
#'     "gaussian"} uses the Gaussian plug-in built from the ensemble mean and
#'     covariance (per-observable marginals plus a joint multivariate value);
#'     \code{"kde"} replaces the per-observable marginals with a Gaussian
#'     kernel-density estimate, which does not assume normality.  With
#'     \code{predictive_moments} the Gaussian score is exact rather than a
#'     plug-in.}
#'   \item{pit}{Probability integral transform \eqn{u = F(y)} per observable.
#'     Calibrated forecasts have \eqn{u \sim U(0,1)}.  From an ensemble the
#'     randomised rank is used (exactly uniform under exchangeability); from
#'     \code{predictive_moments} it is \code{pnorm()} and needs no
#'     randomisation.}
#'   \item{coverage}{Equal-tailed central prediction intervals at each level
#'     in \code{coverage_levels}, with 0/1 coverage indicators.  Averaging the
#'     indicators over many forecast origins gives the empirical coverage
#'     rate; see \code{\link{forecast_backtest}}.}
#' }
#'
#' @section Gaussian predictive (\code{predictive_moments}):
#' For a linear Gaussian state space the predictive distribution is known in
#' closed form, and passing it through \code{predictive_moments =
#' list(mean = , cov = )} makes CRPS, the log score, the PIT and the coverage
#' intervals EXACT instead of Monte-Carlo estimates -- CRPS via
#' Gneiting & Raftery (2007) eq. 22.  The energy and variogram scores have no
#' closed form here and still require draws.
#'
#' @param predictive_draws Numeric matrix (S \eqn{\times} n_obs) of S draws
#'   from the posterior predictive distribution, where each column corresponds
#'   to one observable.  For a single observable, a numeric vector of length S
#'   is also accepted.  May be \code{NULL} when \code{predictive_moments} is
#'   supplied and only closed-form rules are requested.
#' @param observed Numeric vector of length n_obs giving the realised
#'   observations.
#' @param rules Character vector naming the rules to compute.  Any subset of
#'   \code{c("crps", "energy", "variogram", "logs", "pit", "coverage")}.
#'   Default: the first three (unchanged from earlier releases).
#' @param vs_p Variogram score order (default 0.5; recommended by
#'   Scheuerer & Hamill 2015).
#' @param vs_weights Optional weight matrix (n_obs \eqn{\times} n_obs) for the
#'   variogram score.  Default \code{NULL} uses uniform weights.
#' @param obs_vars Optional character vector of length n_obs for labelling
#'   per-observable CRPS values.
#' @param predictive_moments Optional \code{list(mean = , cov = )} giving the
#'   Gaussian predictive mean (length n_obs) and covariance
#'   (n_obs \eqn{\times} n_obs; a length-n_obs vector of variances is also
#'   accepted).  When supplied, the closed-form Gaussian scores are used --
#'   see the section above.
#' @param logs_method Either \code{"gaussian"} (default) or \code{"kde"}; see
#'   the \code{logs} entry above.  Ignored unless \code{"logs"} is requested.
#' @param coverage_levels Numeric vector of nominal central-interval levels in
#'   (0, 1).  Default \code{c(0.5, 0.9)}.  Ignored unless \code{"coverage"} is
#'   requested.
#' @param kde_bw Bandwidth for \code{logs_method = "kde"} (scalar, or one per
#'   observable).  Default \code{NULL} uses \code{stats::bw.nrd0()}.
#' @param pit_randomize Logical; use the randomised ensemble rank for the PIT
#'   (default \code{TRUE}).  Set \code{FALSE} for the plain empirical CDF.
#'   Has no effect on the \code{predictive_moments} path.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{crps}}{Numeric vector of per-observable CRPS (if requested).}
#'     \item{\code{crps_mean}}{Mean CRPS across observables.}
#'     \item{\code{energy}}{Scalar energy score (if requested).}
#'     \item{\code{variogram}}{Scalar variogram score (if requested).}
#'     \item{\code{logs}}{Per-observable log score (if requested).}
#'     \item{\code{logs_mean}}{Mean per-observable log score.}
#'     \item{\code{logs_joint}}{Joint multivariate Gaussian log score
#'       (\code{logs_method = "gaussian"} only; \code{NA} if the ensemble
#'       covariance is not factorisable).}
#'     \item{\code{pit}}{Per-observable PIT value in [0, 1].}
#'     \item{\code{coverage}}{List with \code{lower}, \code{upper},
#'       \code{covered} (n_levels \eqn{\times} n_obs) and \code{levels}.}
#'     \item{\code{rules}}{Character vector of rules computed.}
#'     \item{\code{n_draws}}{Number of ensemble members used (0 on the
#'       \code{predictive_moments}-only path).}
#'     \item{\code{n_obs}}{Number of observables.}
#'   }
#' @export
#'
#' @references
#' Gneiting, T., & Raftery, A. E. (2007). Strictly proper scoring rules,
#' prediction, and estimation. \emph{Journal of the American Statistical
#' Association}, 102(477), 359--378.
#'
#' Scheuerer, M., & Hamill, T. M. (2015). Variogram-based proper scoring
#' rules for probabilistic forecasts of multivariate quantities.
#' \emph{Monthly Weather Review}, 143(4), 1321--1334.
#'
#' Diebold, F. X., Gunther, T. A., & Tay, A. S. (1998). Evaluating density
#' forecasts with applications to financial risk management.
#' \emph{International Economic Review}, 39(4), 863--883.
#'
#' @seealso \code{\link{dm_test}} to compare two score series,
#'   \code{\link{forecast_backtest}} for the recursive expanding-window driver.
#'
#' @examples
#' set.seed(42)
#' S <- 1000L
#' mu <- c(0.5, -0.3)
#' sigma <- c(1.2, 0.8)
#' X <- cbind(rnorm(S, mu[1], sigma[1]),
#'            rnorm(S, mu[2], sigma[2]))
#' y_obs <- c(0.6, -0.5)
#' score_forecast(X, y_obs)
#'
#' ## Calibration diagnostics and the exact Gaussian predictive
#' score_forecast(X, y_obs, rules = c("logs", "pit", "coverage"))
#' score_forecast(NULL, y_obs, rules = c("crps", "logs", "pit"),
#'                predictive_moments = list(mean = mu, cov = sigma^2))
score_forecast <- function(predictive_draws = NULL,
                           observed,
                           rules      = c("crps", "energy", "variogram"),
                           vs_p       = 0.5,
                           vs_weights = NULL,
                           obs_vars  = NULL,
                           predictive_moments = NULL,
                           logs_method     = c("gaussian", "kde"),
                           coverage_levels = c(0.5, 0.9),
                           kde_bw          = NULL,
                           pit_randomize   = TRUE) {

  logs_method <- match.arg(logs_method)
  rules <- match.arg(rules,
                     c("crps", "energy", "variogram", "logs", "pit", "coverage"),
                     several.ok = TRUE)

  ## ---- Gaussian predictive moments (optional) ---------------------------
  pm <- NULL
  if (!is.null(predictive_moments)) {
    pm <- .validate_predictive_moments(predictive_moments)
    if (any(c("energy", "variogram") %in% rules) && is.null(predictive_draws))
      stop("score_forecast: the energy and variogram scores have no ",
           "closed form here -- supply 'predictive_draws' as well as ",
           "'predictive_moments', or drop those rules.", call. = FALSE)
  }

  ## Coerce and validate inputs
  have_draws <- !is.null(predictive_draws)
  if (!have_draws && is.null(pm))
    stop("score_forecast: supply 'predictive_draws' or 'predictive_moments'.",
         call. = FALSE)

  if (have_draws) {
    X <- if (is.null(dim(predictive_draws))) {
      matrix(as.numeric(predictive_draws), ncol = 1L)
    } else {
      as.matrix(predictive_draws)
    }
    S     <- nrow(X)
    n_obs <- ncol(X)
  } else {
    X     <- NULL
    S     <- 0L
    n_obs <- length(pm$mean)
  }
  y <- as.numeric(observed)

  if (length(y) != n_obs)
    stop(sprintf(
      "score_forecast: nrow/ncol mismatch -- predictive_draws has %d columns but observed has length %d",
      n_obs, length(y)))
  if (have_draws && S < 2L)
    stop("score_forecast: need at least 2 ensemble draws")
  if (!is.null(pm) && length(pm$mean) != n_obs)
    stop(sprintf(
      "score_forecast: predictive_moments$mean has length %d but there are %d observables",
      length(pm$mean), n_obs), call. = FALSE)

  if (!is.null(obs_vars) && length(obs_vars) != n_obs)
    stop("obs_names must have length n_obs")
  if (is.null(obs_vars))
    obs_vars <- if (have_draws && !is.null(colnames(X))) colnames(X)
                else paste0("obs_", seq_len(n_obs))

  out <- list(rules  = rules,
              n_draws = S,
              n_obs   = n_obs)

  if ("crps" %in% rules) {
    crps_vec <- if (is.null(pm)) .crps_ensemble(X, y)
                else .crps_gaussian(pm$mean, pm$sd, y)
    names(crps_vec) <- obs_vars
    out$crps        <- crps_vec
    out$crps_mean   <- mean(crps_vec)
  }

  if ("energy" %in% rules) {
    out$energy <- .energy_score(X, y)
  }

  if ("variogram" %in% rules) {
    if (n_obs < 2L) {
      warning("score_forecast: variogram score requires n_obs >= 2; skipping")
      out$variogram <- NA_real_
    } else {
      out$variogram <- .variogram_score(X, y, p = vs_p, w = vs_weights)
    }
  }

  ## ---- Logarithmic score ------------------------------------------------
  if ("logs" %in% rules) {
    if (!is.null(pm)) {
      logs_vec   <- .logs_gaussian_uni(pm$mean, pm$sd, y)
      logs_joint <- .logs_gaussian_joint(pm$mean, pm$cov, y)
    } else if (logs_method == "kde") {
      bw <- if (is.null(kde_bw)) rep(list(NULL), n_obs)
            else as.list(rep_len(as.numeric(kde_bw), n_obs))
      logs_vec <- vapply(seq_len(n_obs),
                         function(j) .logs_kde(X[, j], y[j], bw = bw[[j]]),
                         numeric(1L))
      logs_joint <- NA_real_
    } else {
      mu_hat <- colMeans(X)
      sd_hat <- apply(X, 2L, stats::sd)
      if (any(!is.finite(sd_hat)) || any(sd_hat <= 0))
        stop("score_forecast: degenerate ensemble -- a column of ",
             "predictive_draws has zero variance, so the Gaussian log score ",
             "is undefined. Use logs_method = 'kde' or supply ",
             "predictive_moments.", call. = FALSE)
      logs_vec   <- .logs_gaussian_uni(mu_hat, sd_hat, y)
      logs_joint <- if (n_obs > 1L) .logs_gaussian_joint(mu_hat, stats::cov(X), y)
                    else logs_vec
    }
    names(logs_vec) <- obs_vars
    out$logs        <- logs_vec
    out$logs_mean   <- mean(logs_vec)
    out$logs_joint  <- logs_joint
    out$logs_method <- if (is.null(pm)) logs_method else "gaussian_exact"
  }

  ## ---- PIT --------------------------------------------------------------
  if ("pit" %in% rules) {
    pit_vec <- if (is.null(pm)) .pit_ensemble(X, y, randomize = pit_randomize)
               else .pit_gaussian(pm$mean, pm$sd, y)
    names(pit_vec) <- obs_vars
    out$pit        <- pit_vec
  }

  ## ---- Interval coverage ------------------------------------------------
  if ("coverage" %in% rules) {
    lv <- as.numeric(coverage_levels)
    if (length(lv) < 1L || any(!is.finite(lv)) || any(lv <= 0) || any(lv >= 1))
      stop("score_forecast: coverage_levels must lie strictly inside (0, 1).",
           call. = FALSE)
    cvg <- if (is.null(pm)) .coverage_ensemble(X, y, lv)
           else .coverage_gaussian(pm$mean, pm$sd, y, lv)
    colnames(cvg$lower) <- colnames(cvg$upper) <-
      colnames(cvg$covered) <- obs_vars
    out$coverage <- cvg
  }

  class(out) <- c("dynhr_forecast_scores", "list")
  out
}

#' Validate and normalise a Gaussian `predictive_moments` argument
#'
#' Accepts \code{cov} as a full covariance matrix or as a vector of
#' variances, and returns both the covariance and the marginal sds.
#' @noRd
.validate_predictive_moments <- function(pm) {
  if (!is.list(pm) || is.null(pm$mean) || is.null(pm$cov))
    stop("score_forecast: predictive_moments must be a list with 'mean' and ",
         "'cov' elements.", call. = FALSE)
  mu <- as.numeric(pm$mean)
  k  <- length(mu)
  V  <- pm$cov
  V  <- if (is.null(dim(V))) {
    if (length(V) != k)
      stop(sprintf(paste0("score_forecast: predictive_moments$cov has length ",
                          "%d but 'mean' has length %d."), length(V), k),
           call. = FALSE)
    diag(as.numeric(V), nrow = k)
  } else {
    V <- as.matrix(V)
    if (nrow(V) != k || ncol(V) != k)
      stop(sprintf(paste0("score_forecast: predictive_moments$cov is %d x %d ",
                          "but 'mean' has length %d."), nrow(V), ncol(V), k),
           call. = FALSE)
    V
  }
  sdv <- sqrt(diag(V))
  if (any(!is.finite(sdv)) || any(sdv <= 0))
    stop("score_forecast: predictive_moments$cov must have strictly positive ",
         "diagonal entries.", call. = FALSE)
  list(mean = mu, cov = V, sd = sdv)
}

#' @export
print.dynhr_forecast_scores <- function(x, ...) {
  cat(sprintf("dynhr forecast scores  (S=%d draws, n_obs=%d)\n",
              x$n_draws, x$n_obs))
  if (!is.null(x$crps)) {
    cat(sprintf("  CRPS (mean): %.6g\n", x$crps_mean))
    if (x$n_obs > 1L) {
      cat(sprintf("  CRPS per obs: %s\n",
                  paste(sprintf("%s=%.4g", names(x$crps), x$crps), collapse = ", ")))
    }
  }
  if (!is.null(x$energy))
    cat(sprintf("  Energy score: %.6g\n", x$energy))
  if (!is.null(x$variogram))
    cat(sprintf("  Variogram score (p=default): %.6g\n", x$variogram))
  if (!is.null(x$logs)) {
    cat(sprintf("  Log score [%s] (mean): %.6g\n",
                x$logs_method %||% "gaussian", x$logs_mean))
    if (x$n_obs > 1L) {
      cat(sprintf("  Log score per obs: %s\n",
                  paste(sprintf("%s=%.4g", names(x$logs), x$logs),
                        collapse = ", ")))
      if (!is.null(x$logs_joint) && is.finite(x$logs_joint))
        cat(sprintf("  Log score (joint): %.6g\n", x$logs_joint))
    }
  }
  if (!is.null(x$pit))
    cat(sprintf("  PIT: %s\n",
                paste(sprintf("%s=%.4f", names(x$pit), x$pit),
                      collapse = ", ")))
  if (!is.null(x$coverage)) {
    cv <- x$coverage
    for (l in seq_along(cv$levels))
      cat(sprintf("  %g%% interval covered: %s\n", 100 * cv$levels[l],
                  paste(sprintf("%s=%s", colnames(cv$covered),
                                ifelse(cv$covered[l, ] > 0, "yes", "no")),
                        collapse = ", ")))
  }
  invisible(x)
}


## ---- Diebold-Mariano test ------------------------------------------------

#' Diebold-Mariano test of equal predictive accuracy
#'
#' Tests \eqn{H_0: E[d_t] = 0} for the loss differential
#' \eqn{d_t = s^{(1)}_t - s^{(2)}_t} formed from two series of scores over
#' the SAME forecast origins.  Both inputs must be negatively oriented
#' (lower = better), which every rule in \code{\link{score_forecast}} is, so
#' a NEGATIVE mean differential means model 1 forecasts better.
#'
#' @section Variance and small-sample correction:
#' \eqn{h}-step-ahead forecast errors are MA(\eqn{h-1}) dependent, so the
#' long-run variance of \eqn{d_t} is estimated by the HAC (Newey-West type)
#' estimator truncated at lag \eqn{h-1}:
#' \deqn{\hat V = \gamma_0 + 2 \sum_{k=1}^{L} w_k \gamma_k}
#' with \eqn{L = }\code{lag} (default \eqn{h-1}) and \eqn{w_k = 1}
#' (\code{kernel = "rectangular"}, the original Diebold-Mariano choice) or
#' \eqn{w_k = 1 - k/(L+1)} (\code{kernel = "bartlett"}, which guarantees a
#' non-negative \eqn{\hat V}).
#'
#' The raw statistic \eqn{DM = \bar d / \sqrt{\hat V / n}} is badly oversized
#' in small samples and at long horizons, so by default the
#' Harvey-Leybourne-Newbold (1997) correction is applied,
#' \deqn{DM^* = DM \sqrt{(n + 1 - 2h + h(h-1)/n)/n}}
#' and referred to a \eqn{t_{n-1}} distribution rather than the standard
#' normal.  Set \code{small_sample = FALSE} for the uncorrected
#' \eqn{N(0,1)} version.
#'
#' @param score1,score2 Numeric vectors of the same length: per-origin scores
#'   for model 1 and model 2 (lower = better).  Pairs with an \code{NA} in
#'   either series are dropped.
#' @param h Forecast horizon in periods (default 1).  Drives both the default
#'   HAC truncation lag and the small-sample correction.
#' @param lag HAC truncation lag.  Default \code{NULL} uses \code{h - 1}.
#' @param kernel \code{"rectangular"} (default) or \code{"bartlett"}.
#' @param alternative \code{"two.sided"} (default), \code{"less"} (model 1 is
#'   BETTER, i.e. \eqn{E[d] < 0}), or \code{"greater"}.
#' @param small_sample Apply the Harvey-Leybourne-Newbold correction and use
#'   the \eqn{t_{n-1}} reference distribution (default \code{TRUE}).
#' @return An object of class \code{"dynhr_dm_test"}: a list with
#'   \code{statistic}, \code{p_value}, \code{mean_diff}, \code{var_hac},
#'   \code{n}, \code{h}, \code{lag}, \code{df}, \code{kernel},
#'   \code{alternative} and \code{correction}.
#' @export
#'
#' @references
#' Diebold, F. X., & Mariano, R. S. (1995). Comparing predictive accuracy.
#' \emph{Journal of Business & Economic Statistics}, 13(3), 253--263.
#'
#' Harvey, D., Leybourne, S., & Newbold, P. (1997). Testing the equality of
#' prediction mean squared errors. \emph{International Journal of
#' Forecasting}, 13(2), 281--291.
#'
#' @seealso \code{\link{score_forecast}}, \code{\link{forecast_backtest}}
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' s_good <- rexp(n, rate = 2)          # lower scores = better model
#' s_bad  <- s_good + rexp(n, rate = 4) # strictly worse at every origin
#' dm_test(s_good, s_bad, h = 1, alternative = "less")
dm_test <- function(score1, score2,
                    h = 1L,
                    lag = NULL,
                    kernel = c("rectangular", "bartlett"),
                    alternative = c("two.sided", "less", "greater"),
                    small_sample = TRUE) {
  kernel      <- match.arg(kernel)
  alternative <- match.arg(alternative)

  s1 <- as.numeric(score1)
  s2 <- as.numeric(score2)
  if (length(s1) != length(s2))
    stop(sprintf(paste0("dm_test: score1 has length %d but score2 has ",
                        "length %d -- the two score series must cover the ",
                        "same forecast origins."), length(s1), length(s2)),
         call. = FALSE)

  d  <- s1 - s2
  ok <- is.finite(d)
  d  <- d[ok]
  n  <- length(d)
  if (n < 3L)
    stop(sprintf(paste0("dm_test: only %d usable (finite) paired ",
                        "observations -- need at least 3."), n), call. = FALSE)

  h <- as.integer(h)
  if (is.na(h) || h < 1L) stop("dm_test: h must be a positive integer.",
                               call. = FALSE)
  L <- if (is.null(lag)) h - 1L else as.integer(lag)
  if (is.na(L) || L < 0L) stop("dm_test: lag must be a non-negative integer.",
                               call. = FALSE)
  if (L > n - 1L) L <- n - 1L

  dbar <- mean(d)
  dc   <- d - dbar
  ## Autocovariances gamma_k, k = 0, ..., L (divisor n, as in DM 1995).
  gam <- function(k) if (k == 0L) sum(dc * dc) / n
                     else sum(dc[(k + 1L):n] * dc[1L:(n - k)]) / n
  V <- gam(0L)
  if (L >= 1L) {
    for (k in seq_len(L)) {
      w <- if (kernel == "bartlett") 1 - k / (L + 1) else 1
      V <- V + 2 * w * gam(k)
    }
  }
  if (!is.finite(V) || V <= 0)
    stop(sprintf(paste0("dm_test: the HAC long-run variance estimate is ",
                        "%.6g (not positive). Retry with ",
                        "kernel = \"bartlett\", which is guaranteed ",
                        "non-negative."), V), call. = FALSE)

  stat <- dbar / sqrt(V / n)
  if (isTRUE(small_sample)) {
    ## Harvey-Leybourne-Newbold (1997) eq. 9.
    adj <- (n + 1 - 2 * h + h * (h - 1) / n) / n
    if (adj <= 0)
      stop(sprintf(paste0("dm_test: the Harvey-Leybourne-Newbold correction ",
                          "factor is non-positive (n = %d, h = %d): the ",
                          "sample is too short for this horizon. Use ",
                          "small_sample = FALSE or shorten h."), n, h),
           call. = FALSE)
    stat <- stat * sqrt(adj)
    df   <- n - 1L
    pfun <- function(q) stats::pt(q, df = df)
  } else {
    df   <- Inf
    pfun <- stats::pnorm
  }

  p <- switch(alternative,
              two.sided = 2 * pfun(-abs(stat)),
              less      = pfun(stat),
              greater   = 1 - pfun(stat))

  structure(list(statistic   = stat,
                 p_value     = p,
                 mean_diff   = dbar,
                 var_hac     = V,
                 n           = n,
                 h           = h,
                 lag         = L,
                 df          = df,
                 kernel      = kernel,
                 alternative = alternative,
                 correction  = if (isTRUE(small_sample)) "HLN" else "none"),
            class = c("dynhr_dm_test", "list"))
}

#' @export
print.dynhr_dm_test <- function(x, ...) {
  cat("Diebold-Mariano test of equal predictive accuracy\n")
  cat(sprintf("  n = %d origins, h = %d, HAC lag = %d (%s), correction = %s\n",
              x$n, x$h, x$lag, x$kernel, x$correction))
  cat(sprintf("  mean loss differential (model1 - model2): %+.6g\n",
              x$mean_diff))
  cat(sprintf("  DM statistic: %.4f   df: %s\n", x$statistic,
              if (is.finite(x$df)) format(x$df) else "Inf (normal)"))
  cat(sprintf("  alternative: %s   p-value: %.4g\n",
              x$alternative, x$p_value))
  cat(sprintf("  -> %s\n",
              if (x$mean_diff < 0) "model 1 scores lower (better) on average"
              else "model 2 scores lower (better) on average"))
  invisible(x)
}
