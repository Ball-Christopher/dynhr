## R/diag-scoring-rules.R
## --------------------------------------------------------------------------
## Proper scoring rules for posterior-predictive forecast evaluation.
##
## Implements ensemble (sample-based) estimators for:
##   - CRPS  (Continuous Ranked Probability Score, univariate)
##   - Energy Score  (multivariate generalisation of CRPS)
##   - Variogram Score  (captures dependence structure)
##
## References:
##   Gneiting & Raftery (2007) "Strictly proper scoring rules..."
##     JASA 102(477):359-378.
##   Scheuerer & Hamill (2015) "Variogram-based proper scoring rules"
##     Monthly Weather Review 143(4):1321-1334.
##
## Public API: score_forecast()
## Internal helpers: .crps_ensemble(), .energy_score(), .variogram_score()
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
#' }
#'
#' @param predictive_draws Numeric matrix (S \eqn{\times} n_obs) of S draws
#'   from the posterior predictive distribution, where each column corresponds
#'   to one observable.  For a single observable, a numeric vector of length S
#'   is also accepted.
#' @param observed Numeric vector of length n_obs giving the realised
#'   observations.
#' @param rules Character vector naming the rules to compute.  Any subset of
#'   \code{c("crps", "energy", "variogram")}.  Default: all three.
#' @param vs_p Variogram score order (default 0.5; recommended by
#'   Scheuerer & Hamill 2015).
#' @param vs_weights Optional weight matrix (n_obs \eqn{\times} n_obs) for the
#'   variogram score.  Default \code{NULL} uses uniform weights.
#' @param obs_names Optional character vector of length n_obs for labelling
#'   per-observable CRPS values.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{crps}}{Numeric vector of per-observable CRPS (if requested).}
#'     \item{\code{crps_mean}}{Mean CRPS across observables.}
#'     \item{\code{energy}}{Scalar energy score (if requested).}
#'     \item{\code{variogram}}{Scalar variogram score (if requested).}
#'     \item{\code{rules}}{Character vector of rules computed.}
#'     \item{\code{n_draws}}{Number of ensemble members used.}
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
#' @examples
#' set.seed(42)
#' S <- 1000L
#' mu <- c(0.5, -0.3)
#' sigma <- c(1.2, 0.8)
#' X <- cbind(rnorm(S, mu[1], sigma[1]),
#'            rnorm(S, mu[2], sigma[2]))
#' y_obs <- c(0.6, -0.5)
#' score_forecast(X, y_obs)
score_forecast <- function(predictive_draws,
                           observed,
                           rules      = c("crps", "energy", "variogram"),
                           vs_p       = 0.5,
                           vs_weights = NULL,
                           obs_names  = NULL) {

  ## Coerce and validate inputs
  X <- if (is.null(dim(predictive_draws))) {
    matrix(as.numeric(predictive_draws), ncol = 1L)
  } else {
    as.matrix(predictive_draws)
  }
  y     <- as.numeric(observed)
  S     <- nrow(X)
  n_obs <- ncol(X)

  if (length(y) != n_obs)
    stop(sprintf(
      "score_forecast: nrow/ncol mismatch -- predictive_draws has %d columns but observed has length %d",
      n_obs, length(y)))
  if (S < 2L)
    stop("score_forecast: need at least 2 ensemble draws")

  rules <- match.arg(rules, c("crps", "energy", "variogram"), several.ok = TRUE)

  if (!is.null(obs_names) && length(obs_names) != n_obs)
    stop("obs_names must have length n_obs")
  if (is.null(obs_names))
    obs_names <- if (!is.null(colnames(X))) colnames(X) else paste0("obs_", seq_len(n_obs))

  out <- list(rules  = rules,
              n_draws = S,
              n_obs   = n_obs)

  if ("crps" %in% rules) {
    crps_vec        <- .crps_ensemble(X, y)
    names(crps_vec) <- obs_names
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

  class(out) <- c("dynhr_forecast_scores", "list")
  out
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
  invisible(x)
}
