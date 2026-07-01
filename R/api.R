## R/api.R
## --------------------------------------------------------------------------
## Public API for dynhr -- thin wrappers over the internal estimators.
##
## Exported verbs:
##   prior_spec()     -- extract prior spec data.frame from a parsed model
##   make_posterior() -- build a cached log-posterior closure
##   find_mode()      -- locate the posterior mode
##   mcmc()           -- RWMH sampler  -> dynhr_chains
##   smc()            -- SMC sampler   -> dynhr_chains
##   nuts()           -- NUTS sampler  -> dynhr_chains
##
## S3 class dynhr_chains:
##   new_dynhr_chains()     -- constructor (internal)
##   print.dynhr_chains()   -- compact header + parameter list
##   summary.dynhr_chains() -- posterior table (mean / sd / quantiles)
## --------------------------------------------------------------------------


# ============================================================================
# dynhr_chains S3 class
# ============================================================================

#' Wrap a raw sampler result as a dynhr_chains object
#'
#' Normalises the n_draws field to post-warmup count and stamps the class.
#' Works for rwmh(), dynhr_smc(), and dynhr_nuts() result lists.
#'
#' @param result List from rwmh(), dynhr_smc(), or dynhr_nuts()
#' @param sampler Character: "rwmh", "smc", or "nuts"
#' @return Object of class "dynhr_chains"
#' @noRd
new_dynhr_chains <- function(result, sampler) {
  if (is.null(result$sampler)) result$sampler <- sampler
  result$n_draws <- nrow(result$chain)
  class(result)  <- c("dynhr_chains", class(result))
  result
}


#' Print a dynhr_chains object
#'
#' Shows sampler, draw count, dimension, and key diagnostics.
#'
#' @param x   A \code{dynhr_chains} object
#' @param ... Ignored
#' @return \code{x}, invisibly
#' @export
print.dynhr_chains <- function(x, ...) {
  n_draws <- nrow(x$chain)
  n_par   <- ncol(x$chain)
  sampler <- toupper(x$sampler %||% "UNKNOWN")

  cat(sprintf("\n<dynhr_chains [%s]: %d draws x %d parameters>\n",
              sampler, n_draws, n_par))

  if (!is.null(x$acceptance_rate))
    cat(sprintf("  Acceptance rate : %.1f%%\n", x$acceptance_rate * 100))

  if (!is.null(x$n_divergent) && x$n_divergent > 0)
    cat(sprintf("  Divergences     : %d\n", x$n_divergent))

  if (!is.null(x$log_marginal_lik))
    cat(sprintf("  log p(Y|M)      : %.2f\n", x$log_marginal_lik))

  if (!is.null(x$n_stages))
    cat(sprintf("  Tempering stages: %d\n", x$n_stages))

  if (!is.null(x$elapsed_secs))
    cat(sprintf("  Elapsed         : %.1f sec\n", x$elapsed_secs))

  cat("  Parameters      :", paste(colnames(x$chain), collapse = ", "), "\n")
  invisible(x)
}


#' Posterior summary table for a dynhr_chains object
#'
#' @param object A \code{dynhr_chains} object
#' @param probs  Quantile probabilities (default 5 / 50 / 95 percent)
#' @param ...    Ignored
#' @return data.frame with columns mean, sd, and the requested quantiles
#'   (returned invisibly)
#' @export
summary.dynhr_chains <- function(object, probs = c(0.05, 0.50, 0.95), ...) {
  ch     <- object$chain
  means  <- colMeans(ch)
  sds    <- apply(ch, 2, sd)
  qmat   <- t(apply(ch, 2, quantile, probs = probs))

  out <- data.frame(
    mean = round(means, 4),
    sd   = round(sds,   4),
    round(qmat, 4),
    check.names = FALSE
  )

  sampler <- toupper(object$sampler %||% "?")
  cat(sprintf("\n<dynhr_chains [%s]: %d draws>\n\n", sampler, nrow(ch)))
  print(out)
  invisible(out)
}


# ============================================================================
# Public API functions
# ============================================================================

#' Extract prior specification from a parsed DSGE model
#'
#' Reads the \code{estimated_params} block from \code{model} and returns a
#' data.frame suitable for \code{\link{make_posterior}} and
#' \code{\link{find_mode}}.
#'
#' @param model dynhr_mod object from \code{\link{parse_mod}}
#' @return data.frame with columns: name, distribution, p1, p2, lower, upper,
#'   mean, std
#' @seealso \code{\link{make_posterior}}, \code{\link{find_mode}}
#' @export
prior_spec <- function(model) extract_prior_spec(model)


#' Build a cached log-posterior evaluator
#'
#' Constructs a closure that evaluates \eqn{\log p(\theta | Y)} for any
#' parameter vector \eqn{\theta}, caching the model structure to avoid
#' recompilation on each call.  Call this once before \code{\link{mcmc}},
#' \code{\link{smc}}, or \code{\link{nuts}}.
#'
#' @param model       dynhr_mod from \code{\link{parse_mod}}
#' @param data        Observation matrix (\eqn{T \times n_{\text{obs}}}),
#'   columns matching \code{obs_vars}
#' @param prior_spec  Prior spec data.frame from \code{\link{prior_spec}}
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled from \code{\link{compile_model}}
#' @param me_variance Measurement error variance added to the diagonal of the
#'   observation noise (default 0; use a small positive value for
#'   stochastically-singular models)
#' @param likelihood  Likelihood type: \code{"gaussian"} (Kalman filter,
#'   default) or \code{"cumulant"} (cumulant-matching, Mutschler 2015).
#'   The cumulant likelihood works with any perturbation order >= 1;
#'   orders >= 2 provide skewness/kurtosis content.
#' @param ...         Additional arguments forwarded to the internal
#'   likelihood constructor. For \code{likelihood = "cumulant"}, supports
#'   \code{order} (perturbation order, default 2), \code{cumulant_orders}
#'   (which orders to match, default 1:4), and \code{cumulant_weight}
#'   ("identity" or "precision").
#' @return A function \code{function(theta)} returning a named list
#'   \code{list(logpost, loglik, logprior)}
#' @seealso \code{\link{prior_spec}}, \code{\link{find_mode}}, \code{\link{mcmc}}
#' @export
make_posterior <- function(model, data, prior_spec, obs_vars, compiled,
                           me_variance = 0,
                           likelihood = c("gaussian", "cumulant", "pruned"),
                           ...) {
  likelihood <- match.arg(likelihood)
  make_log_posterior(model, data, prior_spec, obs_vars, compiled,
                     me_variance, likelihood = likelihood, ...)
}


#' Find the posterior mode
#'
#' Runs a multi-stage optimizer to locate \eqn{\arg\max_\theta \log p(\theta|Y)}.
#' The returned mode is the natural starting point for \code{\link{mcmc}} or
#' \code{\link{nuts}}, and the inverse Hessian at the mode provides the
#' initial proposal covariance for RWMH.
#'
#' @param log_post_fn Log-posterior closure from \code{\link{make_posterior}}
#' @param theta_init  Named numeric starting vector (e.g. prior means)
#' @param prior_spec  Prior spec from \code{\link{prior_spec}} (used for
#'   parameter bounds)
#' @param n_iter      Total optimizer iteration budget (default 10000)
#' @param method      Optimizer sequence. Default \code{"newrat"} (csminwel
#'   quasi-Newton, equivalent to Dynare \code{mode_compute = 4}): gradient-based
#'   with an infeasibility-backtracking line search, the most robust and
#'   cheapest choice on near-unit-root models. Other options:
#'   \code{"cmaes_newrat"} (CMA-ES global search then newrat polish -- use for
#'   multimodal posteriors or poor starts), \code{"combined"} (CMA-ES then
#'   L-BFGS-B; the former default), \code{"cmaes"}, \code{"nmkb"},
#'   \code{"jade"}, \code{"nelder"}, \code{"cmaes_nmkb"}, \code{"cmaes_jade"}
#' @param verbose     Print progress (default \code{TRUE})
#' @return Named list: \code{theta_mode}, \code{logpost}, \code{convergence},
#'   \code{iterations}, \code{method}
#' @seealso \code{\link{mcmc}}, \code{\link{nuts}}, \code{\link{make_posterior}}
#'
#' @references
#'   Hansen, N. (2016). The CMA evolution strategy: A tutorial.
#'     \emph{arXiv:1604.00772}.
#'   Nelder, J. A., & Mead, R. (1965). A simplex method for function
#'     minimization. \emph{The Computer Journal}, 7(4), 308-313.
#' @export
find_mode <- function(log_post_fn, theta_init, prior_spec,
                      n_iter = 10000L, method = "newrat", verbose = TRUE) {
  .run_mode_finding(log_post_fn, theta_init, prior_spec,
                    nm_maxit = n_iter, method = method, verbose = verbose)
}


#' Random Walk Metropolis-Hastings sampler
#'
#' Runs RWMH-MCMC with adaptive proposal scaling and returns a
#' \code{\link{dynhr_chains}} object.  A good starting point is the mode
#' from \code{\link{find_mode}} with the scaled inverse-Hessian as
#' \code{Sigma_prop}.
#'
#' @param log_post_fn Log-posterior closure from \code{\link{make_posterior}}
#' @param theta0      Named numeric starting vector
#' @param Sigma_prop  Proposal covariance matrix (\eqn{n_p \times n_p});
#'   scale by \eqn{2.38^2 / n_p} when using the mode Hessian
#' @param n_draws     Post-warmup draws to retain (default 2000)
#' @param n_warmup    Warmup draws for adaptive tuning, then discarded
#'   (default 1000)
#' @param ...         Additional arguments forwarded to the internal sampler:
#'   \code{scale}, \code{target_rate}, \code{adapt_every}, \code{verbose}
#' @return A \code{\link{dynhr_chains}} object whose \code{$chain} slot is the
#'   \eqn{n_{\text{draws}} \times n_p} post-warmup sample matrix
#' @seealso \code{\link{smc}}, \code{\link{nuts}}, \code{\link{find_mode}}
#'
#' @references
#'   Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., &
#'     Teller, E. (1953). Equation of state calculations by fast computing
#'     machines. \emph{Journal of Chemical Physics}, 21(6), 1087-1092.
#'   Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains
#'     and their applications. \emph{Biometrika}, 57(1), 97-109.
#'   Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2013).
#'     \emph{Bayesian Data Analysis} (3rd ed.). Chapman & Hall/CRC.
#' @export
mcmc <- function(log_post_fn, theta0, Sigma_prop,
                 n_draws = 2000L, n_warmup = 1000L, ...) {
  res <- rwmh(log_post_fn, theta0, Sigma_prop,
              n_draws = n_draws + n_warmup, n_burn = n_warmup, ...)
  new_dynhr_chains(res, "rwmh")
}


#' Sequential Monte Carlo sampler
#'
#' Runs SMC with adaptive likelihood tempering.  No mode estimate is required;
#' particles are initialised from the prior.  Returns a
#' \code{\link{dynhr_chains}} object.
#'
#' @param log_post_fn Log-posterior closure from \code{\link{make_posterior}}
#' @param prior_spec  Prior spec from \code{\link{prior_spec}}
#' @param n_particles Number of particles (default 2000; use 5000+ for
#'   production runs)
#' @param ...         Additional arguments forwarded to the internal sampler:
#'   \code{ess_target}, \code{n_mh_steps}, \code{mh_scale_factor},
#'   \code{parallel}, \code{verbose}
#' @return A \code{\link{dynhr_chains}} object.  The \code{$log_marginal_lik}
#'   slot contains the log marginal likelihood estimate
#'   \eqn{\log p(Y | \mathcal{M})}
#' @seealso \code{\link{mcmc}}, \code{\link{nuts}}
#'
#' @references
#'   Herbst, E. P., & Schorfheide, F. (2015). \emph{Bayesian Estimation of
#'     DSGE Models}. Princeton University Press. Chapter 10 (SMC).
#'   Chopin, N. (2002). A sequential particle filter method for static models.
#'     \emph{Biometrika}, 89(3), 539-552.
#' @export
smc <- function(log_post_fn, prior_spec, n_particles = 2000L, ...) {
  res <- dynhr_smc(log_post_fn, prior_spec = prior_spec,
                   n_particles = n_particles, ...)
  new_dynhr_chains(res, "smc")
}


#' No-U-Turn Sampler (NUTS)
#'
#' Runs NUTS with dual-averaging step-size adaptation and optional diagonal
#' mass-matrix adaptation.  Returns a \code{\link{dynhr_chains}} object.
#'
#' @param log_post_fn Log-posterior closure from \code{\link{make_posterior}}
#' @param theta0      Named numeric starting vector
#' @param n_draws     Post-warmup draws to retain (default 2000)
#' @param n_warmup    Warmup iterations for step-size / mass adaptation,
#'   then discarded (default 1000)
#' @param ...         Additional arguments forwarded to the internal sampler:
#'   \code{step_size}, \code{max_treedepth}, \code{target_accept},
#'   \code{adapt_mass}, \code{verbose}
#' @return A \code{\link{dynhr_chains}} object.  The \code{$n_divergent}
#'   slot reports divergent transitions (should be 0 in a well-tuned run)
#' @seealso \code{\link{mcmc}}, \code{\link{smc}}, \code{\link{find_mode}}
#'
#' @references
#'   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: Adaptively
#'     setting path lengths in Hamiltonian Monte Carlo.
#'     \emph{Journal of Machine Learning Research}, 15(1), 1593-1623.
#'   Neal, R. M. (2011). MCMC using Hamiltonian dynamics.
#'     In \emph{Handbook of Markov Chain Monte Carlo}. Chapman & Hall/CRC.
#' @export
nuts <- function(log_post_fn, theta0, n_draws = 2000L, n_warmup = 1000L, ...) {
  res <- dynhr_nuts(log_post_fn, theta0,
                    n_draws = n_draws, n_warmup = n_warmup, ...)
  new_dynhr_chains(res, "nuts")
}


#' DIME ensemble MCMC sampler
#'
#' Runs the Differential-Independence Mixture Ensemble (DIME) sampler
#' (Boehl 2022/2024).  Gradient-free, swarm-based, and robust to
#' multimodal posteriors.  Particles are initialised from the prior; no
#' mode estimate is required.  Returns a \code{\link{dynhr_chains}} object.
#'
#' @param log_post_fn Log-posterior closure from \code{\link{make_posterior}}
#' @param prior_spec  Prior spec from \code{\link{prior_spec}}
#' @param n_chain     Number of ensemble walkers (default NULL =
#'   max(5 * n_par, 20))
#' @param n_iter      Post-warmup iterations per walker (default 1000)
#' @param n_burn      Burn-in iterations (default 500)
#' @param ...         Additional arguments forwarded to \code{run_dime()}:
#'   \code{aimh_prob}, \code{sigma}, \code{rho}, \code{df}, \code{verbose}
#' @return A \code{\link{dynhr_chains}} object.  The \code{$chain} matrix has
#'   \eqn{n_{\text{iter}} \times n_{\text{chain}}} rows (all walkers, all
#'   post-warmup iterations).
#' @seealso \code{\link{mcmc}}, \code{\link{smc}}, \code{\link{nuts}}
#'
#' @references
#'   Boehl, G. (2022). DIME MCMC: A Swiss Army Knife for Bayesian Inference.
#'     SSRN Working Paper No. 4250395.
#' @export
dime <- function(log_post_fn, prior_spec, n_chain = NULL,
                 n_iter = 1000L, n_burn = 500L, ...) {
  res <- run_dime(log_post_fn, prior_spec = prior_spec,
                  n_chain = n_chain, n_iter = n_iter, n_burn = n_burn, ...)
  new_dynhr_chains(res, "dime")
}


# ============================================================================
# Diagnostic API aliases
# ============================================================================

#' Run the dynhr diagnostic battery
#'
#' Runs all available diagnostics given a set of estimation outputs and
#' returns a named list of \code{dynhr_diagnostic} objects.
#'
#' When the first argument is a \code{dynhr_posterior_result} object, all
#' required inputs are derived automatically.  In this mode, the \code{report}
#' argument controls output format.
#'
#' @param ...  Arguments forwarded to \code{run_all_diagnostics()}.  When the
#'   first argument is a \code{dynhr_posterior_result}, the \code{report},
#'   \code{report_file}, and \code{output_dir} arguments are available.
#' @return Named list of \code{dynhr_diagnostic} objects.
#' @seealso \code{\link{run_full_estimation}}, \code{\link{write_report}},
#'   \code{\link{solve_model}}, \code{\link{run_posterior_estimation}}
#' @export
run_diagnostics <- function(...) run_all_diagnostics(...)


# write_report() is defined in diag-report-llm.R (the full implementation)
# with format = c("llm", "html") dispatch.


# ============================================================================
# New simplified pipeline API (v0.3+)
# ============================================================================

#' Solve a DSGE model (parse -> compile -> steady -> perturbation)
#'
#' Convenience alias for \code{\link{solve_model}}.
#'
#' @param ...  Arguments forwarded to \code{\link{solve_model}}.
#' @return A \code{dynhr_solved} object.
#' @seealso \code{\link{solve_model}}, \code{\link{run_mode_finding}}
#' @export
solve_dsge <- function(...) solve_model(...)


#' Run mode-finding for a solved DSGE model
#'
#' Convenience alias for \code{\link{run_mode_finding}}.
#'
#' @param ...  Arguments forwarded to \code{\link{run_mode_finding}}.
#' @return A \code{dynhr_mode_result} object.
#' @seealso \code{\link{run_mode_finding}}, \code{\link{run_posterior_estimation}}
#' @export
estimate_mode <- function(...) run_mode_finding(...)


#' Run posterior estimation with multi-method dispatch
#'
#' Convenience alias for \code{\link{run_posterior_estimation}}.
#'
#' @param ...  Arguments forwarded to \code{\link{run_posterior_estimation}}.
#' @return A \code{dynhr_posterior_result} object.
#' @seealso \code{\link{run_posterior_estimation}}, \code{run_all_diagnostics}
#' @export
estimate_posterior <- function(...) run_posterior_estimation(...)


#' Two-stage SMC model-tempering (Mlikota & Schorfheide 2024)
#'
#' Automates the unbiased Bayes-factor workflow described in Mlikota &
#' Schorfheide (2024) "Sequential Monte Carlo with Model Tempering":
#'
#' \strong{Stage 1} runs standard likelihood-tempering SMC on the approximating
#' model M0 to obtain its posterior cloud and \eqn{\log Z_{M0}}.
#'
#' \strong{Stage 2} runs a model-tempering bridge SMC that starts from the M0
#' posterior cloud (\code{init_particles}) and bridges from M0 to M1.  The
#' result's \code{log_marginal_lik} equals \eqn{\log Z_{M1}} (the full MDD of
#' M1, not merely the ratio), because Stage 2 accumulates
#' \eqn{\log(Z_{M1}/Z_{M0})} and the wrapper adds \eqn{\log Z_{M0}}.
#'
#' ## M0 specification
#'
#' \describe{
#'   \item{\code{approx_loglik_fn} only (no \code{log_post_fn_M0}):}{
#'     \code{approx_loglik_fn(theta)} returns a scalar log-likelihood for M0.
#'     The wrapper constructs the full M0 log-posterior by combining it with
#'     the prior (which is shared between M0 and M1).}
#'   \item{Full \code{log_post_fn_M0} provided:}{
#'     Must return \code{list(logpost, loglik, logprior)} -- use this when M0
#'     has a different prior, or when a more efficient M0 evaluator exists.}
#' }
#'
#' ## Returned fields
#'
#' All fields from the Stage-2 (M1) SMC result are present for compatibility
#' with \code{\link[=smc]{smc()}} / \code{new_dynhr_chains()}, plus:
#' \itemize{
#'   \item \code{log_marginal_lik}: \eqn{\log Z_{M1}} (unbiased).
#'   \item \code{log_Z_M0}: \eqn{\log Z_{M0}} from Stage 1.
#'   \item \code{log_ratio}: \eqn{\log(Z_{M1}/Z_{M0})} from Stage 2.
#'   \item \code{stage1}: full Stage-1 SMC result for inspection.
#'   \item \code{marginal_valid}: always \code{TRUE} (M0 cloud supplied).
#' }
#'
#' @param log_post_fn_M1 function(theta) -> list(logpost, loglik, logprior).
#'   The M1 log-posterior.  Same interface as \code{\link{smc}}'s
#'   \code{log_post_fn}.
#' @param approx_loglik_fn function(theta) -> scalar.  M0 log-likelihood.
#'   Used in Stage 1 to build the M0 posterior and in Stage 2 as the bridge
#'   denominator.
#' @param prior_sampler function() -> named numeric vector drawn from the
#'   prior.  Required when \code{log_post_fn_M0 = NULL}.
#' @param prior_spec Prior specification (data.frame or named list).
#'   Alternative to \code{prior_sampler}; forwarded to both stages.
#' @param log_post_fn_M0 Optional.  Full M0 log-posterior function (same
#'   interface as \code{log_post_fn_M1}).  When \code{NULL} (default), built
#'   automatically from \code{approx_loglik_fn} + the shared prior.
#' @param n_particles Number of particles (same for both stages).
#' @param ess_target ESS ratio target for adaptive tempering.
#' @param n_mh_steps RWMH mutation steps per tempering stage.
#' @param seed_M0 RNG seed for Stage 1 (M0 run).  Default \code{1L}.
#' @param seed_M1 RNG seed for Stage 2 (M1 bridge run).  Default \code{2L}.
#' @param verbose Print progress messages.
#' @param ... Additional arguments forwarded to both internal
#'   \code{dynhr_smc()} calls: e.g. \code{mh_scale_factor},
#'   \code{mixture_weights}, \code{parallel}, \code{backend}.
#'
#' @return A \code{\link{dynhr_chains}} object (M1 posterior) with extra
#'   slots \code{log_Z_M0}, \code{log_ratio}, and \code{stage1}.
#'
#' @seealso \code{\link{smc}}, \code{\link{mcmc}}, \code{\link{nuts}}
#'
#' @references
#'   Mlikota, M., & Schorfheide, F. (2024). Sequential Monte Carlo with
#'     Model Tempering. \emph{Journal of Econometrics}.
#' @export
smc_model_tempered <- function(
    log_post_fn_M1,
    approx_loglik_fn,
    prior_sampler  = NULL,
    prior_spec     = NULL,
    log_post_fn_M0 = NULL,
    n_particles    = 2000L,
    ess_target     = 0.5,
    n_mh_steps     = 1L,
    seed_M0        = 1L,
    seed_M1        = 2L,
    verbose        = TRUE,
    ...
) {
  res <- dynhr_smc_model_tempered(
    log_post_fn_M1   = log_post_fn_M1,
    approx_loglik_fn = approx_loglik_fn,
    prior_sampler    = prior_sampler,
    prior_spec       = prior_spec,
    log_post_fn_M0   = log_post_fn_M0,
    n_particles      = n_particles,
    ess_target       = ess_target,
    n_mh_steps       = n_mh_steps,
    seed_M0          = seed_M0,
    seed_M1          = seed_M1,
    verbose          = verbose,
    ...
  )
  ch <- new_dynhr_chains(res, "smc")
  ## Attach wrapper-specific bookkeeping fields as attributes so they survive
  ## new_dynhr_chains() and remain inspectable.
  ch$log_Z_M0  <- res$log_Z_M0
  ch$log_ratio <- res$log_ratio
  ch$stage1    <- res$stage1
  ch
}

