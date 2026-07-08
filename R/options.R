## R/options.R
## Package-level option store for dynhr.
##
## dynhr_set_options() stores named values that override function argument
## defaults. Useful when one option (e.g. me_variance) must flow through many
## nested calls without threading it through every signature.
##
## Priority: explicit argument > global option > hard-coded default.

.dynhr_opts <- new.env(parent = emptyenv())

#' Set dynhr package-level options
#'
#' Named values stored here act as defaults for any function that calls
#' `.dynhr_opt(name)`. Explicit argument values always take precedence.
#'
#' Recognised options (non-exhaustive):
#' \describe{
#'   \item{`me_variance`}{Measurement-error variance vector or scalar, passed
#'     to `make_log_posterior()`, `run_mode_parallel()`, `run_mcmc_parallel()`.}
#'   \item{`perturb_scale`}{Starting-point perturbation scale for mode-finding.}
#'   \item{`nm_maxit`}{Iteration budget for each mode-finding chain.}
#'   \item{`seed_base`}{Base RNG seed for reproducible parallel runs.}
#'   \item{`proposal_cov_method`}{\code{"diagonal"} (default) or
#'     \code{"full"}; selects the Step 6 proposal-covariance strategy in
#'     \code{run_mode_finding()} -- see \code{proposal_cov}.}
#'   \item{`transform_params`}{\code{TRUE} (default); mode finding
#'     (\code{run_mode_finding()}) and the samplers
#'     (\code{run_posterior_estimation()}: RWMH, NUTS, MALA, ChEES, HMC) operate
#'     on an unconstrained reparameterisation built by
#'     \code{build_param_transform} instead of the raw, box-constrained
#'     parameter vector. Default \code{TRUE} because every gradient sampler
#'     step-collapses against the bounds of parameters near a unit root; the
#'     transform leaves the target invariant. Set \code{FALSE} to sample in raw
#'     theta-space -- see \code{\link{run_mode_finding}}.}
#'   \item{`rwmh_adapt_cov`}{\code{FALSE} (default); when \code{TRUE}, RWMH
#'     chains in \code{run_posterior_estimation()} use Haario et al. (2001)
#'     adaptive proposal covariance -- see \code{rwmh}'s `adapt_cov`.}
#'   \item{`rwmh_n_blocks`}{\code{1L} (default); number of randomized
#'     parameter blocks for RWMH chains in
#'     \code{run_posterior_estimation()} -- see \code{rwmh}'s
#'     `n_blocks`.}
#'   \item{`grad_method`}{\code{"hybrid"} (default) or \code{"implicit"};
#'     selects the \code{grad_method} passed to
#'     \code{\link{make_posterior_grad}} when \code{analytic_grad = TRUE} on
#'     the serial NUTS path of \code{run_posterior_estimation()} /
#'     \code{run_full_estimation()}. \code{"implicit"} uses the full
#'     implicit-differentiation gradient (Childers et al. 2022); has no effect
#'     on the parallel multi-chain (mirai) NUTS path, which always uses the
#'     numerical gradient.}
#'   \item{`power_posterior`}{Scalar in \code{(0, 1]} (default \code{1});
#'     power-posterior (generalised-Bayes) tempering exponent \eqn{\zeta}.
#'     The log-posterior returned by \code{make_log_posterior()} becomes
#'     \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)}: values below 1
#'     down-weight a potentially misspecified likelihood and widen the
#'     posterior (Bissiri, Holmes & Walker 2016, JRSS-B; Grünwald SafeBayes).
#'     The \code{$loglik} field always carries the \emph{raw} (untempered)
#'     log-likelihood; only \code{$logpost} is tempered. Default \code{1} is
#'     bit-identical to the standard Bayesian posterior and is compatible with
#'     every existing sampler without any changes.}
#'   \item{`allow_monge_metric`}{\code{FALSE} (default); the experimental
#'     \code{metric = "monge"} MALA option in \code{run_posterior_estimation()}
#'     errors unless this is \code{TRUE}, because the Monge metric collapses the
#'     proposal step on sharp / near-unit-root posteriors (it under-explored a
#'     tight direction to ~1\% of its variance in testing). It still warns when
#'     enabled; prefer \code{metric = "hessian"} (constant Laplace).}
#'   \item{`allow_whittle_fim_metric`}{\code{FALSE} (default); the experimental
#'     \code{metric = "whittle_fim"} NUTS option in
#'     \code{run_posterior_estimation()} errors unless this is \code{TRUE}.
#'     \code{"whittle_fim"} assembles a one-shot dense mass matrix from the
#'     Whittle (frequency-domain) Fisher information at the mode (no periodic
#'     recompute); it is unvalidated on the full estimation pipeline and costs
#'     one extra evaluation that can be expensive on large models. Falls back
#'     to the \code{"hessian"} metric on failure.}
#'   \item{`checkpoint_flush_every`}{\code{1000L} (default); the chunk size (in
#'     draws) at which checkpointed RWMH flushes its streaming buffer to disk and
#'     rewrites the restart state, when \code{run_posterior_estimation()} is
#'     called with \code{checkpoint_dir}. Smaller values bound RAM more tightly
#'     and lose less on an interruption, at the cost of more frequent I/O.}
#'   \item{`debug_kf_errors`}{\code{FALSE} (default); when \code{TRUE},
#'     a Kalman-filter error inside \code{make_log_posterior()} is RE-RAISED
#'     instead of being caught and treated as an infeasible draw
#'     (\code{loglik = -Inf}). The default keeps an infeasible parameter draw
#'     (singular covariance, unit-root with stationary init, ...) from crashing a
#'     long chain; set this \code{TRUE} when debugging to surface a genuine code
#'     bug that would otherwise be silently masked as a rejected draw.}
#' }
#'
#' @param ... Named option values.
#' @return Invisibly, the previous values of all modified options.
#' @export
dynhr_set_options <- function(...) {
  args <- list(...)
  if (length(args) == 0L) return(invisible(list()))
  nms <- names(args)
  if (is.null(nms) || any(nms == ""))
    stop("All arguments to dynhr_set_options() must be named.")
  prev <- lapply(nms, function(nm) {
    p <- if (exists(nm, envir = .dynhr_opts, inherits = FALSE))
      get(nm, envir = .dynhr_opts) else NULL
    assign(nm, args[[nm]], envir = .dynhr_opts)
    p
  })
  names(prev) <- nms
  invisible(prev)
}

#' Get all current dynhr package-level options
#' @return Named list of all currently set options.
#' @export
dynhr_get_options <- function() as.list(.dynhr_opts)

#' Reset dynhr options to package defaults
#'
#' @param ... Character names of options to reset. If empty, resets all.
#' @export
dynhr_reset_options <- function(...) {
  nms <- c(...)
  if (length(nms) == 0L) {
    rm(list = ls(.dynhr_opts), envir = .dynhr_opts)
  } else {
    for (nm in nms)
      if (exists(nm, envir = .dynhr_opts, inherits = FALSE))
        rm(list = nm, envir = .dynhr_opts)
  }
  invisible(NULL)
}

# Internal: resolve option with arg-value > global option > default priority.
# Usage: val <- .dynhr_opt("me_variance", arg_value, default = 0)
# - If arg_value is non-NULL, return it as-is.
# - Else if the option is set globally, return that.
# - Else return `default`.
.dynhr_opt <- function(name, arg_value = NULL, default = NULL) {
  if (!is.null(arg_value)) return(arg_value)
  if (exists(name, envir = .dynhr_opts, inherits = FALSE))
    return(get(name, envir = .dynhr_opts))
  default
}
