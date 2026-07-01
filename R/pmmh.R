## R/pmmh.R
## --------------------------------------------------------------------------
## pmmh() -- a thin, documented entry point for Particle Marginal
## Metropolis-Hastings. PMMH is not a separate sampler: it is a Random-Walk
## Metropolis-Hastings chain driven by an UNBIASED, NOISY particle-filter
## estimate of the marginal likelihood. By the pseudo-marginal argument
## (Andrieu & Roberts 2009; Andrieu, Doucet & Holenstein 2010) such a chain
## targets exactly the same posterior as one using the exact likelihood.
##
## In dynhr the unbiased particle likelihoods are:
##   - "tpf"  : Tempered Particle Filter (Herbst & Schorfheide 2019), for
##              nonlinear/order-2 state spaces.
##   - "ppf"  : OBC bootstrap particle filter (occasionally-binding
##              constraints).
##   - "copf" : OBC conditionally-optimal particle filter.
##
## So PMMH is exactly:
##
##   run_posterior_estimation(mode_result, methods = "RWMH",
##                            likelihood = "tpf" | "ppf" | "copf")
##
## with the likelihood already baked into mode_result$ctx (set when the mode
## was found). pmmh() simply forwards to run_posterior_estimation with
## methods = "RWMH" fixed, and is the discoverable, self-documenting handle on
## that capability. The correctness of the construction (noisy-loglik chain ->
## exact posterior) is proven end-to-end by tests/testthat/test-pmmh-oracle.R,
## which shows PMMH recovers the exact-Kalman posterior on a linear-Gaussian
## AR(1) within Monte-Carlo error.
## --------------------------------------------------------------------------

#' Particle Marginal Metropolis-Hastings (PMMH)
#'
#' Runs Particle Marginal Metropolis-Hastings: a Random-Walk Metropolis
#' chain whose acceptance ratio uses an \emph{unbiased particle-filter
#' estimate} of the marginal likelihood in place of the (intractable) exact
#' likelihood. Because the estimate is unbiased, the chain targets exactly the
#' same posterior as the exact-likelihood chain (the pseudo-marginal
#' guarantee of Andrieu, Doucet & Holenstein 2010).
#'
#' \code{pmmh()} is a thin convenience wrapper. It is \emph{identical} to
#' \preformatted{run_posterior_estimation(mode_result, methods = "RWMH", ...)}
#' when \code{mode_result$ctx$likelihood} is one of the unbiased particle
#' filters (\code{"tpf"}, \code{"ppf"}, \code{"copf"}). The likelihood is
#' chosen when the posterior mode is found (e.g. \code{find_mode(...,
#' likelihood = "tpf")}); \code{pmmh()} reads it from \code{mode_result} and
#' fixes the sampler to RWMH.
#'
#' @section Particle count and the variance preflight:
#' PMMH performance hinges on the Monte-Carlo variance of the log-likelihood
#' estimate at the mode. The rule of thumb (Doucet, Pitt, Deligiannidis &
#' Kohn 2015; the Dynare / Herbst-Schorfheide threshold) is to choose the
#' particle count \eqn{N} so that \eqn{\mathrm{sd}(\log \hat L) \approx 1}: too
#' few particles makes acceptance dominated by loglik noise; too many wastes
#' computation. \code{run_posterior_estimation} runs a preflight
#' (\code{\link{tpf_loglik_sd_preflight}} / the internal
#' \code{.tpf_pmcmc_preflight}) before sampling and \strong{warns with an
#' N-recommendation} when the SD exceeds 1. Set the particle count via
#' \code{nparticles} (and/or the filter's own option list:
#' \code{tpf_options$n_particles} for TPF, \code{obc_options$N} for OBC). Use
#' \code{seed = NULL} in the likelihood closure (the default): each proposal
#' must draw a \emph{fresh, independent} particle-filter estimate; a fixed seed
#' would reuse one pseudo-likelihood draw and break the PMMH invariance.
#'
#' @section Correlated pseudo-marginal (CPM):
#' Setting \code{tpf_options$cpm_rho_u} (a value in \eqn{(0,1)}) switches the
#' TPF path to the correlated pseudo-marginal sampler (Deligiannidis, Doucet &
#' Pitt 2018), which correlates consecutive loglik evaluations to cut variance
#' at fixed \eqn{N}. This is handled automatically by
#' \code{run_posterior_estimation}.
#'
#' @param mode_result A mode-finding result whose \code{ctx$likelihood} is
#'   \code{"tpf"}, \code{"ppf"}, or \code{"copf"}. (A Gaussian-likelihood
#'   \code{mode_result} would run ordinary exact-likelihood RWMH, not PMMH;
#'   \code{pmmh()} warns in that case.)
#' @param ... Passed through to \code{\link{run_posterior_estimation}} (e.g.
#'   \code{ndraws}, \code{nburn}, \code{nchains}, \code{nparticles},
#'   \code{Sigma_prop}, \code{parallel}, \code{checkpoint_dir}). The
#'   \code{methods} argument is fixed to \code{"RWMH"} and must not be passed.
#'
#' @return The list returned by \code{\link{run_posterior_estimation}}.
#'
#' @seealso \code{\link{run_posterior_estimation}},
#'   \code{\link{tpf_loglik_sd_preflight}},
#'   \code{\link{make_log_posterior_obc_ppf}}.
#'
#' @references
#' Andrieu, C., Doucet, A. & Holenstein, R. (2010). Particle Markov chain
#' Monte Carlo methods. \emph{JRSS-B} 72(3): 269-342.
#'
#' Doucet, A., Pitt, M.K., Deligiannidis, G. & Kohn, R. (2015). Efficient
#' implementation of Markov chain Monte Carlo when using an unbiased
#' likelihood estimator. \emph{Biometrika} 102(2): 295-313.
#'
#' Herbst, E. & Schorfheide, F. (2019). Tempered particle filtering.
#' \emph{J. Econometrics} 210(1): 26-44.
#'
#' @examples
#' \dontrun{
#' ## 1. Find the mode under a particle likelihood:
#' mode_res <- find_mode(model, data, prior_spec, obs_vars,
#'                       likelihood = "tpf", me_variance = 0.01,
#'                       tpf_options = list(n_particles = 2000L))
#'
#' ## 2. Run PMMH (RWMH over the unbiased TPF loglik). The variance preflight
#' ##    runs automatically and warns if more particles are needed.
#' post <- pmmh(mode_res, ndraws = 20000L, nburn = 5000L, nchains = 4L,
#'              nparticles = 2000L)
#'
#' ## Equivalent long form:
#' post <- run_posterior_estimation(mode_res, methods = "RWMH",
#'                                  ndraws = 20000L, nparticles = 2000L)
#' }
#' @export
pmmh <- function(mode_result, ...) {
  dots <- list(...)
  if ("methods" %in% names(dots)) {
    stop("pmmh(): do not pass `methods` -- PMMH fixes the sampler to RWMH. ",
         "Use run_posterior_estimation() directly for other samplers.")
  }

  lik <- tryCatch(mode_result$ctx$likelihood, error = function(e) NULL)
  if (is.null(lik) || !lik %in% c("tpf", "ppf", "copf")) {
    warning("pmmh(): mode_result$ctx$likelihood = ",
            if (is.null(lik)) "NULL" else sQuote(lik),
            " is not an unbiased particle likelihood (\"tpf\"/\"ppf\"/\"copf\"). ",
            "This will run ordinary exact-likelihood RWMH, NOT PMMH. Build the ",
            "mode_result with likelihood = \"tpf\" (or \"ppf\"/\"copf\") for PMMH.",
            call. = FALSE)
  }

  do.call(run_posterior_estimation,
          c(list(mode_result = mode_result, methods = "RWMH"), dots))
}
