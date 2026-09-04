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

  ## Delayed-acceptance chains (mcmc(screen_fn = )) carry the stage-1 kill
  ## rate and the count of expensive evaluations actually paid for.
  if (!is.null(x$screen_rate))
    cat(sprintf("  Screen rate     : %.1f%% (stage-1 rejects; %d expensive evals)\n",
                x$screen_rate * 100, x$n_expensive %||% NA_integer_))

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
  cat(sprintf("\n<dynhr_chains [%s]: %d draws>\n", sampler, nrow(ch)))
  if (!is.null(object$screen_rate))
    cat(sprintf("Screen rate: %.1f%% of proposals rejected at stage 1; %d expensive evaluations\n",
                object$screen_rate * 100, object$n_expensive %||% NA_integer_))
  cat("\n")
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
#' @seealso \code{\link{make_posterior}}, \code{\link{find_mode}},
#'   \code{\link{dynhr_model}} (the pipeline-object equivalent)
#' @examples
#' ## nk_demo.mod carries a nine-parameter estimated_params block
#' model  <- parse_mod(system.file("extdata/models/nk_demo.mod",
#'                                 package = "dynhr"), verbose = FALSE)
#' priors <- prior_spec(model)
#' priors[, c("name", "distribution", "p1", "p2")]
#'
#' ## prior means are the conventional starting point for find_mode()
#' theta0 <- setNames(priors$mean, priors$name)
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
#' @seealso \code{\link{prior_spec}}, \code{\link{find_mode}}, \code{\link{mcmc}},
#'   \code{\link{dm_posterior}} (the pipeline-object equivalent)
#' @examples
#' ## Model and data both ship with the package
#' model    <- parse_mod(system.file("extdata/models/nk_demo.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, verbose = FALSE)
#' priors   <- prior_spec(model)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#'
#' log_post <- make_posterior(model, data = Y, prior_spec = priors,
#'                            obs_vars = obs_vars, compiled = compiled)
#'
#' theta0 <- setNames(priors$mean, priors$name)
#' str(log_post(theta0))
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
#' @seealso \code{\link{mcmc}}, \code{\link{nuts}}, \code{\link{make_posterior}},
#'   \code{\link{dm_mode}} (the pipeline-object equivalent)
#'
#' @references
#'   Hansen, N. (2016). The CMA evolution strategy: A tutorial.
#'     \emph{arXiv:1604.00772}.
#'   Nelder, J. A., & Mead, R. (1965). A simplex method for function
#'     minimization. \emph{The Computer Journal}, 7(4), 308-313.
#' @examples
#' model    <- parse_mod(system.file("extdata/models/nk_demo.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, verbose = FALSE)
#' priors   <- prior_spec(model)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' log_post <- make_posterior(model, data = Y, prior_spec = priors,
#'                            obs_vars = obs_vars, compiled = compiled)
#'
#' theta0 <- setNames(priors$mean, priors$name)
#' mode <- find_mode(log_post, theta0, priors, n_iter = 200L, verbose = FALSE)
#' mode$logpost
#' round(mode$theta_mode, 4)
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
#' @param checkpoint_dir Optional directory for streaming checkpoints. When
#'   set, draws are streamed to per-chain files in \code{flush_every}-row
#'   chunks and a restart state (position, proposal scale and covariance,
#'   acceptance count, and -- crucially -- \code{.Random.seed}) is saved
#'   after every flush, so an interrupted run loses at most one flush
#'   window.
#' @param resume      Logical (default \code{FALSE}). With
#'   \code{checkpoint_dir} set, continue the saved run in that directory
#'   instead of starting fresh: pass a larger \code{n_draws} to extend a
#'   finished chain. The continuation is \emph{bit-identical} to having run
#'   the larger \code{n_draws} in one process (exact RNG and adaptation
#'   restore). Keep \code{n_warmup} at its original value: the saved warmup
#'   length governs the retained-draw split, while \code{n_warmup} still
#'   enters the total draw target \code{n_draws + n_warmup}. A checkpoint
#'   written under a different parameter configuration is refused. For bespoke samplers outside
#'   \code{mcmc()}, the same contract is available via
#'   \code{\link{mcmc_chain_state}} / \code{\link{mcmc_chain_restore}}.
#' @param flush_every Rows per checkpoint flush (default 1000; only used
#'   with \code{checkpoint_dir}).
#' @param screen_fn   Optional CHEAP log-posterior closure over the same
#'   parameters (e.g. an order-1 Kalman posterior when \code{log_post_fn} is
#'   an order-2 pruned or particle-filter posterior). When supplied, the
#'   sampler switches to \strong{delayed acceptance} (Christen & Fox 2005):
#'   each proposal is first accept/rejected against \code{screen_fn}, and only
#'   a survivor pays for an evaluation of \code{log_post_fn}, which is then
#'   accepted with the screen-corrected ratio
#'   \eqn{[\pi(\theta')/\pi(\theta)] \cdot [c(\theta)/c(\theta')]}. The chain
#'   targets \code{log_post_fn} \emph{exactly} however biased the screen is;
#'   a bad screen costs efficiency only. \code{screen_fn} must be
#'   deterministic in \eqn{\theta} (verified at \code{theta0}). The returned
#'   object gains \code{$screen_rate} (share of proposals killed by the
#'   screen) and \code{$n_expensive} (evaluations of \code{log_post_fn}).
#'   Pass \code{pseudo_marginal = TRUE} through \code{...} when
#'   \code{log_post_fn} is a particle filter: the sampler then refuses a
#'   fixed-seed closure and runs the PMCMC variance preflight.
#' @param ...         Additional arguments forwarded to the internal sampler:
#'   \code{scale}, \code{target_rate}, \code{adapt_every}, \code{verbose}
#'   (and, with \code{screen_fn}, \code{pseudo_marginal},
#'   \code{preflight_K}, \code{adapt_cov})
#' @return A \code{\link{dynhr_chains}} object whose \code{$chain} slot is the
#'   \eqn{n_{\text{draws}} \times n_p} post-warmup sample matrix
#' @seealso \code{\link{smc}}, \code{\link{nuts}}, \code{\link{find_mode}},
#'   \code{\link{mcmc_chain_state}}, \code{\link{dm_sample}} (the
#'   pipeline-object equivalent)
#'
#' @references
#'   Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., &
#'     Teller, E. (1953). Equation of state calculations by fast computing
#'     machines. \emph{Journal of Chemical Physics}, 21(6), 1087-1092.
#'   Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains
#'     and their applications. \emph{Biometrika}, 57(1), 97-109.
#'   Gelman, A., Carlin, J. B., Stern, H. S., & Rubin, D. B. (2013).
#'     \emph{Bayesian Data Analysis} (3rd ed.). Chapman & Hall/CRC.
#' @examples
#' ## Any function(theta) -> list(logpost, loglik, logprior) is a valid target;
#' ## in practice it comes from make_posterior().
#' log_post <- function(theta) {
#'   ll <- -0.5 * sum(((theta - c(0.5, -0.2)) / c(1, 0.5))^2)
#'   list(logpost = ll, loglik = ll, logprior = 0)
#' }
#'
#' set.seed(1)
#' chains <- mcmc(log_post, theta0 = c(a = 0, b = 0),
#'                Sigma_prop = diag(c(1, 0.25)),
#'                n_draws = 500L, n_warmup = 200L, verbose = FALSE)
#' chains
#' summary(chains)
#' @export
mcmc <- function(log_post_fn, theta0, Sigma_prop,
                 n_draws = 2000L, n_warmup = 1000L,
                 checkpoint_dir = NULL, resume = FALSE,
                 flush_every = 1000L, screen_fn = NULL, ...) {
  if (isTRUE(resume) && is.null(checkpoint_dir))
    stop("mcmc: resume = TRUE requires 'checkpoint_dir' (the directory of ",
         "the run to continue).", call. = FALSE)
  if (!is.null(screen_fn)) {
    ## Delayed acceptance: `screen_fn` is a CHEAP stand-in for `log_post_fn`
    ## used only to kill hopeless proposals before the expensive target is
    ## touched. The chain still targets `log_post_fn` exactly (see
    ## R/sampler-da.R); a biased screen costs efficiency, never correctness.
    ckpt <- NULL
    if (!is.null(checkpoint_dir)) {
      if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
      ckpt <- list(dir = checkpoint_dir, flush_every = flush_every,
                   resume = isTRUE(resume))
    }
    res <- rwmh_da(log_post_fn, screen_fn, theta0, Sigma_prop,
                   n_draws = n_draws + n_warmup, n_burn = n_warmup,
                   checkpoint = ckpt, ...)
    return(new_dynhr_chains(res, "rwmh_da"))
  }
  if (is.null(checkpoint_dir)) {
    res <- rwmh(log_post_fn, theta0, Sigma_prop,
                n_draws = n_draws + n_warmup, n_burn = n_warmup, ...)
  } else {
    if (!dir.exists(checkpoint_dir))
      dir.create(checkpoint_dir, recursive = TRUE)
    ckpt <- list(dir = checkpoint_dir, flush_every = flush_every,
                 resume = isTRUE(resume),
                 fingerprint = .ckpt_fingerprint(names(theta0)))
    res <- rwmh(log_post_fn, theta0, Sigma_prop,
                n_draws = n_draws + n_warmup, n_burn = n_warmup,
                checkpoint = ckpt, ...)
  }
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
#' @seealso \code{\link{mcmc}}, \code{\link{nuts}}, \code{\link{dm_sample}}
#'
#' @references
#'   Herbst, E. P., & Schorfheide, F. (2015). \emph{Bayesian Estimation of
#'     DSGE Models}. Princeton University Press. Chapter 10 (SMC).
#'   Chopin, N. (2002). A sequential particle filter method for static models.
#'     \emph{Biometrika}, 89(3), 539-552.
#' @examples
#' model    <- parse_mod(system.file("extdata/models/nk_demo.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, verbose = FALSE)
#' priors   <- prior_spec(model)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' log_post <- make_posterior(model, data = Y, prior_spec = priors,
#'                            obs_vars = obs_vars, compiled = compiled)
#'
#' set.seed(1)
#' ## n_particles = 100 keeps the example quick; use 2000+ in practice
#' chains <- smc(log_post, priors, n_particles = 100L, verbose = FALSE)
#' chains$log_marginal_lik
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
#' @seealso \code{\link{mcmc}}, \code{\link{smc}}, \code{\link{find_mode}},
#'   \code{\link{dm_sample}}
#'
#' @references
#'   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: Adaptively
#'     setting path lengths in Hamiltonian Monte Carlo.
#'     \emph{Journal of Machine Learning Research}, 15(1), 1593-1623.
#'   Neal, R. M. (2011). MCMC using Hamiltonian dynamics.
#'     In \emph{Handbook of Markov Chain Monte Carlo}. Chapman & Hall/CRC.
#' @examples
#' ## NUTS needs a gradient; without an analytic one it finite-differences the
#' ## target, so on a real DSGE posterior prefer make_posterior_grad(). Here a
#' ## cheap analytic target keeps the example fast.
#' log_post <- function(theta) {
#'   ll <- -0.5 * sum(((theta - c(0.5, -0.2)) / c(1, 0.5))^2)
#'   list(logpost = ll, loglik = ll, logprior = 0)
#' }
#'
#' set.seed(1)
#' chains <- nuts(log_post, theta0 = c(a = 0, b = 0),
#'                n_draws = 200L, n_warmup = 200L, verbose = FALSE)
#' chains$n_divergent          # 0 in a well-tuned run
#' summary(chains)
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
#' @seealso \code{\link{mcmc}}, \code{\link{smc}}, \code{\link{nuts}},
#'   \code{\link{dm_sample}}
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
#'
#' @section Coverage:
#' The battery spans D0-D41.  \code{run_diagnostics()} wires D0-D37 and D41 and
#' silently skips any stage whose inputs are absent, so it is safe to call with
#' whatever pieces you have.  Three members are exported standalone instead,
#' because they need only a model plus a parameter vector: D38
#' (\code{\link{diag_sloppiness}}), D39 (\code{\link{diag_stability_map}}) and
#' D40 (\code{\link{diag_near_unit_root}}).  See
#' \code{vignette("diagnostics")} for the stage table.
#'
#' @seealso \code{\link{run_full_estimation}}, \code{\link{write_report}},
#'   \code{\link{solve_model}}, \code{\link{run_posterior_estimation}},
#'   \code{\link{diag_sloppiness}}, \code{\link{diag_stability_map}},
#'   \code{\link{diag_near_unit_root}}, \code{\link{dm_diagnostics}}
#' @examples
#' model    <- parse_mod(system.file("extdata/models/rbc.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, verbose = FALSE)
#' steady   <- solve_steady(compiled, model$param_values,
#'                          endo_names = model$var_names,
#'                          exo_names  = model$varexo_names, verbose = FALSE)
#' dr <- solve_perturbation(model, compiled, steady$values,
#'                          model$param_values, verbose = FALSE)
#' set.seed(1)
#' Y <- as.matrix(simulate_model(dr, n_periods = 100L, model = model,
#'                               burn_in = 20L)[, "y", drop = FALSE])
#'
#' ## Pre-estimation stages only: no draws, so D5-D7 and the fit block skip.
#' diags <- suppressWarnings(
#'   run_diagnostics(model = model, compiled = compiled, dr = dr,
#'                   params = model$param_values, data = Y,
#'                   obs_names = "y", verbose = FALSE))
#' names(diags)
#' @export
run_diagnostics <- function(...) run_all_diagnostics(...)


# write_report() is defined in diag-report-llm.R (the full implementation)
# with format = c("llm", "html") dispatch.


# ============================================================================
# New simplified pipeline API (v0.3+)
# ============================================================================
#
# The three `function(...)` pass-through aliases that used to live here --
# solve_dsge() -> solve_model(), estimate_mode() -> run_mode_finding(),
# estimate_posterior() -> run_posterior_estimation() -- were REMOVED in the
# API-consistency wave. They added a second name for each verb with no
# signature of their own (so `args()`, tab-completion and the help pages were
# all empty `...`), and nothing in the package, the tests or the vignettes
# called them. Use the canonical names directly.


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

