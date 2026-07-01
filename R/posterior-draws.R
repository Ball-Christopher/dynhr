## R/posterior-draws.R
## --------------------------------------------------------------------------
## as_posterior_draws() -- resample-or-weight adapter for SMC output.
##
## SMC returns a weighted particle cloud; all downstream posterior consumers
## (diagnostics, Bayesian IRF, smoother) expect an equally-weighted draw
## matrix.  This adapter converts between the two representations by
## systematic resampling, reusing the same helper already used inside the
## SMC tempering loop (.smc_systematic_resample, sampler-smc.R).
##
## Usage:
##   draws <- as_posterior_draws(smc_result)
##   draws <- as_posterior_draws(smc_result, n_draws = 2000L, seed = 1L)
##
## The result is a plain n_draws x n_params matrix, compatible with every
## consumer that currently receives chains$chain.
## --------------------------------------------------------------------------


#' Convert an SMC result to a posterior draw matrix
#'
#' Converts the weighted particle cloud from \code{dynhr_smc()} or
#' \code{run_smc_mirai()} into an equally-weighted draw matrix by systematic
#' resampling.  The resulting matrix has the same column layout as an MCMC
#' \code{chain} matrix and is suitable for all posterior consumers
#' (diagnostics, Bayesian IRF, smoother, posterior-mean computation).
#'
#' When \code{smc_result$smc_weights} is absent (e.g. an older result that
#' pre-dates weight storage), the function assumes the particles are already
#' equally weighted (no resampling required) and returns the chain as-is or
#' subsampled uniformly to \code{n_draws}.
#'
#' @param smc_result  A list returned by \code{dynhr_smc()} or
#'   \code{run_smc_mirai()}.  Must contain at least \code{$chain} (N x d
#'   matrix) and optionally \code{$smc_weights} (length-N normalised weight
#'   vector).
#' @param n_draws     Number of draws to produce (default: number of
#'   particles, i.e. \code{nrow(smc_result$chain)}).
#' @param strategy    \code{"resample"} (default) performs systematic
#'   resampling from the weighted cloud.  \code{"direct"} returns the
#'   particle matrix as-is (or uniformly sub/supersampled), ignoring weights.
#'   Use \code{"direct"} only when \code{ess_target} was high enough that
#'   final weights are known to be approximately uniform.
#' @param seed        Optional integer seed for the resampling RNG.  When
#'   \code{NULL} (default), the current RNG state is used.
#' @param ess_warn    ESS fraction below which a warning is emitted (default
#'   \code{0.5}).  Set \code{NULL} to suppress.
#' @return Numeric matrix, \code{n_draws x d}, with column names matching the
#'   estimated parameter names.  Equally weighted (no weight attribute).
#' @noRd
as_posterior_draws <- function(smc_result,
                                n_draws   = NULL,
                                strategy  = c("resample", "direct"),
                                seed      = NULL,
                                ess_warn  = 0.5) {

  strategy <- match.arg(strategy)

  chain <- smc_result$chain
  if (!is.matrix(chain))
    chain <- as.matrix(chain)
  N <- nrow(chain)
  d <- ncol(chain)

  if (is.null(n_draws))
    n_draws <- N

  ## ---- Check / retrieve stored weights ------------------------------------
  w_norm <- smc_result$smc_weights
  has_weights <- !is.null(w_norm) && length(w_norm) == N && !all(w_norm == w_norm[1])

  ## ---- ESS warning --------------------------------------------------------
  if (has_weights && !is.null(ess_warn) && is.numeric(ess_warn)) {
    ess <- 1 / sum(w_norm^2)
    if (ess < ess_warn * n_draws) {
      warning(sprintf(
        paste0("as_posterior_draws: SMC effective sample size (ESS = %.0f) is ",
               "below %.0f%% of n_draws = %d (threshold = %.0f). The ",
               "resampled draws may not adequately represent the posterior. ",
               "Consider increasing n_particles or ess_target."),
        ess, ess_warn * 100, n_draws, ess_warn * n_draws
      ), call. = FALSE)
    }
  }

  ## ---- Apply strategy -----------------------------------------------------
  if (strategy == "direct" || !has_weights) {
    ## No-weight path: return as-is or uniformly subsample / re-index.
    if (n_draws == N) return(chain)
    if (!is.null(seed)) set.seed(seed)
    idx <- if (n_draws <= N) {
      sample.int(N, n_draws, replace = FALSE)
    } else {
      sample.int(N, n_draws, replace = TRUE)
    }
    return(chain[idx, , drop = FALSE])
  }

  ## Systematic resampling (low-variance, O(N))
  if (!is.null(seed)) set.seed(seed)
  idx <- .smc_systematic_resample(w_norm, n_draws)
  chain[idx, , drop = FALSE]
}


#' Check whether an SMC result needs the posterior-draws adapter
#'
#' Returns TRUE when \code{result} is an SMC result whose final particle
#' weights may be non-uniform (i.e. \code{smc_weights} was stored and is not
#' a uniform vector).  Used by the orchestration layer to decide whether to
#' call \code{as_posterior_draws()} before passing draws to consumers.
#'
#' @param result  Any sampler result list (rwmh, nuts, smc).
#' @return Logical scalar.
#' @noRd
.is_smc_weighted <- function(result) {
  if (is.null(result$sampler) || result$sampler != "smc")
    return(FALSE)
  w <- result$smc_weights
  if (is.null(w) || length(w) < 2L)
    return(FALSE)
  ## Non-trivial weights: coefficient of variation > 1e-6
  (sd(w) / max(mean(w), 1e-300)) > 1e-6
}
