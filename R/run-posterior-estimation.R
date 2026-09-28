## R/run-posterior-estimation.R
## --------------------------------------------------------------------------
## run_posterior_estimation() -- Multi-method posterior sampling dispatcher.
##
## Given a mode-finding result (dynhr_mode_result), runs one or more MCMC
## samplers in sequence, with automatic multi-chain convergence diagnostics
## and stoch_simul at the posterior mean.
##
## Dispatch rule:
##   - Vector arguments (methods, n_warmup, n_draws, n_chains) are paired with
##     each method. Recycled to length(methods).
##   - Scalar arguments apply to ALL methods.
## --------------------------------------------------------------------------


#' Run posterior estimation with multi-method dispatch
#'
#' Runs one or more MCMC samplers sequentially, starting each from the mode
#' (or the previous chain's final state).  Automatically computes convergence
#' diagnostics and stoch_simul (IRFs + moments) at the posterior mean of
#' each batch.
#'
#' @param mode_result  A \code{dynhr_mode_result} from
#'   \code{\link{run_mode_finding}}.
#' @param methods      Character vector of sampler names:
#'   \code{"RWMH"}, \code{"HMC"}, \code{"NUTS"}, \code{"SMC"},
#'   \code{"DSMH"} (\code{\link{dynhr_dsmh}}; parallel over chain groups on the
#'   mirai backend when \code{parallel = TRUE}; first tempering rung
#'   from the number of data rows).
#'   Default \code{c("RWMH")}.
#' @param n_warmup     Burn-in / warmup draws per method. Recycled to
#'   \code{length(methods)}.  Default \code{5000}.
#' @param n_draws      Post-warmup draws to retain per method. Recycled.
#'   Default \code{20000}.
#' @param n_chains     Number of MCMC chains per method (for RWMH). Recycled.
#'   Default \code{4}.
#' @param n_particles  Number of SMC particles (only for \code{"SMC"}), or
#'   the DSMH cloud size \eqn{N G} (only for \code{"DSMH"}).
#'   Recycled.  Default \code{2000}.
#' @param n_walkers    Number of ensemble walkers for the DIME sampler
#'   (\code{NULL} = \code{max(5 * n_par, 20)}).  Ignored for all other
#'   samplers.
#' @param parallel     Run multi-chain \code{"RWMH"}/\code{"NUTS"} batches and
#'   \code{"SMC"} in parallel (default \code{FALSE}). When \code{TRUE}, chains
#'   are dispatched to a persistent \code{mirai} daemon pool: standard
#'   Gaussian-likelihood mode results recompile the posterior once per daemon
#'   (\code{run_mcmc_mirai()}/\code{run_nuts_mirai()}/\code{run_smc_mirai()}
#'   via \code{.mirai_pool_init()}); OBC/PKF and cumulant-likelihood mode
#'   results instead ship the already-built log-posterior closure once via
#'   \code{.mirai_pool_closure()}. Otherwise chains run sequentially.
#' @param parallel_backend Parallel backend for RWMH chains: \code{"mirai"}
#'   (default).
#' @param n_cores      Worker (daemon) count when \code{parallel = TRUE}
#'   (\code{NULL} = auto-detect, capped at \code{n_chains}).
#' @param analytic_grad For the gradient samplers (\code{"NUTS"}, \code{"HMC"},
#'   \code{"MALA"}, \code{"CHEES"}), use the analytic gradient from
#'   \code{\link{make_posterior_grad}} (method: the \code{grad_method} option,
#'   default \code{"auto"}, exact \code{"adjoint_solution"} for the Gaussian
#'   likelihood) instead of a numerical one. Default \code{TRUE}: the exact
#'   gradient wherever the likelihood has one; a likelihood without one (the
#'   OBC filters, \code{"tpf"}, \code{"pskf"}, ...) keeps the numerical
#'   gradient without a warning (the verbose output names the gradient each
#'   stage used). \code{FALSE} always finite-differences the log-posterior.
#'   Applies to the serial samplers and to the parallel multi-chain NUTS path
#'   (each daemon builds the theta-space gradient; the sampler applies the
#'   eta chain rule once). The low-level \code{\link{nuts}} called without
#'   \code{grad_fn} still finite-differences; the runners are what supply
#'   the exact gradient.
#' @param Sigma_prop   Proposal covariance.  Defaults to the one stored in
#'   \code{mode_result$Sigma_prop}.  If a list, one per method.
#' @param run_stoch_simul  If \code{TRUE} (default), compute IRFs and
#'   theoretical moments at the posterior mean of each batch.
#' @param skip_mode_finding_check  If \code{FALSE} (default), abort before
#'   sampling when \code{mode_result} looks unusable (non-finite mode
#'   log-posterior, missing/all-NA \code{theta_mode}, or a \code{Sigma_prop}
#'   that looks like a Hessian-inversion fallback, e.g. a scaled identity or
#'   near-singular matrix). Set to \code{TRUE} to force MCMC despite these
#'   issues (e.g. with a hand-specified \code{Sigma_prop}, or for debugging).
#' @param nuts_timeout_seconds  Wall-clock timeout (seconds, default 300) for
#'   the \code{"NUTS"} sampler; if exceeded, NUTS sampling for that batch is
#'   abandoned and a diagnostic message is emitted.
#' @param transform_params  Default \code{TRUE}: run the random-walk
#'   and gradient-based samplers -- \code{"RWMH"}, \code{"NUTS"}, \code{"MALA"},
#'   \code{"CHEES"} and \code{"HMC"} (both the serial and parallel/mirai paths) -- in
#'   unconstrained eta-space via \code{build_param_transform} (built from
#'   \code{prior_spec}). \code{Sigma_prop} (theta-space) is converted to an
#'   eta-space covariance via the delta method (also for the correlated
#'   pseudo-marginal RWMH path, which then samples eta too), and the NUTS/ChEES
#'   \code{mass_diag} and the NUTS/MALA dense \code{metric = "hessian"} are
#'   similarly converted -- except at a mode on a prior
#'   bound, where \code{mode_result$Sigma_prop_eta} (see
#'   \code{\link{run_mode_finding}}) is used instead when no
#'   \code{Sigma_prop} is supplied (the delta method divides a bound
#'   parameter's variance by the squared distance to the bound there); draws are mapped back to
#'   theta-space before being returned, so downstream consumers (convergence
#'   diagnostics, stoch_simul, etc.) are unaffected. \strong{Strongly
#'   recommended for models with bounded parameters near a unit root}: the
#'   leapfrog/Langevin steps of the gradient samplers otherwise run into the
#'   bounds. Measured on a 68-parameter NZ DSGE (20 observables, T = 132;
#'   NUTS, adapted diagonal metric, 300 warmup + 300 draws, 2026-09-25 on the
#'   corrected mass adaptation): without the transform 246 of 300
#'   transitions were divergent and the minimum bulk ESS was 1.3; with it
#'   there were 0 divergences and a minimum bulk ESS of 57. May also be set
#'   globally via \code{dynhr_set_options(transform_params = TRUE)}.
#' @param rwmh_adapt_cov  Opt-in (default \code{FALSE}): forwarded to
#'   \code{rwmh}'s \code{adapt_cov} argument on the RWMH chains (serial and
#'   parallel/mirai paths) -- Haario et al. (2001) adaptive proposal
#'   covariance, frozen at the end of burn-in. May also be set globally via
#'   \code{dynhr_set_options(rwmh_adapt_cov = TRUE)}.
#' @param rwmh_n_blocks   Opt-in (default \code{1L}): forwarded to
#'   \code{rwmh}'s \code{n_blocks} argument on the RWMH chains (serial and
#'   parallel/mirai paths) -- randomized parameter blocking
#'   (Chib & Ramamurthy 2010 / Herbst & Schorfheide 2015 ch. 4). May also be
#'   set globally via \code{dynhr_set_options(rwmh_n_blocks = 2L)}.
#' @param metric       Metric for MALA, NUTS and HMC:
#'   \code{"diagonal"} (default; NUTS starts from the inverse-Hessian diagonal,
#'   HMC from the identity, and both ADAPT a diagonal inverse mass = the
#'   warmup posterior variance, Stan's convention; MALA uses the identity),
#'   \code{"warmup_dense"} (NUTS/HMC: diagonal windows, then one Ledoit-Wolf
#'   dense inverse mass from the final window's draws),
#'   \code{"fisher_diag"} (NUTS only, opt-in: diagonal inverse mass
#'   \eqn{\sqrt{\mathrm{var}(x)/\mathrm{var}(\nabla \log p)}} from each
#'   warmup window's draws and gradients; Seyboldt, Carlson & Carpenter 2026,
#'   arXiv:2603.18845),
#'   \code{"lowrank"} (NUTS only, opt-in: low-rank-plus-diagonal inverse mass
#'   from draws and gradients; Seyboldt et al. 2026 Algorithm 1 with the Lao
#'   2026 arXiv:2607.23788 schedule),
#'   \code{"hessian"} (dense metric from \code{Sigma_prop}; in eta-space,
#'   \eqn{D^{-1} \Sigma D^{-1}} with \eqn{D = diag(d\theta/d\eta)}, under
#'   \code{transform_params}),
#'   \code{"monge"} (position-dependent Monge metric \eqn{G = I + \alpha^2 g g^T},
#'   Stage 3a; MALA only), or \code{"whittle_fim"} (one-shot dense metric from
#'   the Whittle (frequency-domain) Fisher information at the mode; NUTS only).
#'   \strong{The \code{"whittle_fim"} metric is experimental and disabled by
#'   default}: it costs one extra ~seconds-to-minutes-scale evaluation at the
#'   mode (no periodic recompute -- frozen through warmup and sampling, same
#'   as \code{"hessian"}) and is unvalidated on the full estimation pipeline.
#'   It therefore errors unless
#'   \code{dynhr_set_options(allow_whittle_fim_metric = TRUE)} is set. On
#'   failure (degenerate FIM, non-invertible result) it falls back to the
#'   \code{"hessian"} metric (Sigma_prop, converted to eta-space exactly as
#'   for \code{metric = "hessian"}) with a message, or to the diagonal metric
#'   when Sigma_prop gives no positive-definite dense metric, mirroring the
#'   existing \code{"hessian"} branch's own not-PD fallback to diagonal.
#'   Measured result (2026-07-03, NZSIM 68 params): the NUTS warmup step
#'   size collapsed to ~1e-13 with divergences under this metric (the
#'   SoftAbs-floored FIM is badly scaled in weakly-identified directions),
#'   while the \code{"hessian"} metric warms up stably on the same
#'   posterior -- i.e. no measured win; prefer \code{"hessian"}.
#'   \strong{The \code{"monge"} metric is experimental and
#'   disabled by default}: \eqn{G} inflates along the gradient, so the proposal
#'   step collapses on sharp / near-unit-root posteriors -- in testing it
#'   under-explored a tight direction to ~1\% of its variance (ESS/draw ~0.002)
#'   where the constant \code{"hessian"} metric recovered it fully. It therefore
#'   errors unless \code{dynhr_set_options(allow_monge_metric = TRUE)} is set, and
#'   warns when used. Prefer \code{"hessian"} (constant Laplace) or \code{"diagonal"}.
#'   The parallel multi-chain NUTS path (\code{parallel = TRUE},
#'   \code{n_chains > 1}, mirai) takes the adapted metrics (\code{"diagonal"},
#'   \code{"warmup_dense"}, \code{"fisher_diag"}, \code{"lowrank"}); a fixed
#'   metric (\code{"hessian"}, \code{"whittle_fim"}) is not sent to the
#'   daemons, which then run \code{"diagonal"}.
#' @param monge_alpha  Softness parameter \eqn{\alpha \geq 0} for the Monge
#'   metric (default 1). Only used when \code{metric = "monge"}, which also
#'   requires \code{allow_monge_metric = TRUE} (see \code{metric}).
#'   May also be set globally via \code{dynhr_set_options(monge_alpha = 2)}.
#' @param checkpoint_dir  Optional directory enabling memory-streamed,
#'   restartable sampling. Supported for \code{"RWMH"}, \code{"NUTS"},
#'   \code{"MALA"}, \code{"CHEES"} and \code{"DIME"} (serial; plus the parallel
#'   mirai path for RWMH and NUTS). When set, each chain streams its draws to
#'   per-chain files (\code{chain_<id>.draws} / \code{.lp} / \code{.state.rds}) in
#'   flush-sized chunks, so RAM during sampling is bounded by the flush window
#'   rather than \code{n_draws * n_par}, and a restart state is saved after each
#'   flush. The files are per-chain, so parallel chains never collide. Flush size
#'   is set via \code{dynhr_set_options(checkpoint_flush_every = ...)} (default 1000).
#' @param resume  When \code{TRUE} and \code{checkpoint_dir} points at a prior
#'   run, continue each chain from its saved state -- RNG, position, scale and
#'   proposal covariance are restored exactly, so the continuation is identical
#'   to a single longer run -- adding \code{n_draws} more retained draws. The
#'   saved model / prior / parameter configuration must match (it is enforced).
#'   The checkpoint also records the spec's content hashes and the build: a
#'   resume whose target (model, data, likelihood, mode, sampler settings)
#'   differs, or that crosses a registered result change touching this run,
#'   is refused (see \code{\link{run_estimation}}, section Checkpoints).
#' @param on_mismatch  What a \code{resume} does when the checkpoint's target
#'   differs from this call's, or a registered result change between the two
#'   dynhr versions touches this run: \code{"refuse"} (default; an error) or
#'   \code{"warn"} (resume anyway; the result is then marked in
#'   \code{$provenance$integrity}, its run record and \code{print()}).
#'   A code or environment difference without such a change only warns.
#' @param verbose      Print progress messages.
#' @param seed         Optional integer. When non-\code{NULL},
#'   \code{set.seed(seed)} is called once at the very top of this function,
#'   \strong{before any mode-finding multistart or sampler runs} -- it seeds
#'   the entire call (chain dispersal draws, RWMH/NUTS/MALA/HMC/ChEES/SMC/DIME
#'   RNG, and any multistart perturbations triggered downstream), not just one
#'   sampler. This makes the whole \code{run_posterior_estimation()} call
#'   reproducible without the caller needing to wrap it in their own
#'   \code{set.seed()}. Default \code{NULL}: the ambient RNG stream is left
#'   untouched, i.e. identical to wrapping the call in \code{set.seed()}
#'   yourself (old behaviour). \strong{Note}: previously a \code{seed =}
#'   argument passed by a caller was silently absorbed by \code{...} and
#'   forwarded to samplers that do not have a \code{seed} parameter, so it was
#'   dropped without warning or error -- this explicit argument closes that
#'   gap. The parallel (mirai) RWMH, NUTS, SMC and DIME paths use it as their
#'   daemons' base seed (before 0.9.3.105 they ignored it and always used
#'   their fixed default, so different seeds gave identical parallel draws);
#'   with \code{seed = NULL} they keep that default.
#' @param ...          Additional arguments forwarded to each sampler
#'   (e.g. \code{target_accept = 0.25}, \code{adapt_every = 100}).
#'
#' @return An object of class \code{"dynhr_posterior_result"} containing:
#'   \describe{
#'     \item{\code{mode_result}}{The input mode result}
#'     \item{\code{chains}}{Named list of per-method chain results}
#'     \item{\code{pooled_draws}}{Combined draws from all chains/methods}
#'     \item{\code{convergence}}{Convergence diagnostics (R-hat, ESS)}
#'     \item{\code{posterior_mean}}{Posterior mean parameter vector}
#'     \item{\code{posterior_irfs}}{IRFs at posterior mean (if computed)}
#'     \item{\code{posterior_moments}}{Moments at posterior mean (if computed)}
#'     \item{\code{meta}}{Run metadata}
#'     \item{\code{resolved}}{What the sampler stages resolved at run time:
#'       \code{grad_method}, the analytic-gradient method each gradient-based
#'       stage used (named by method key, e.g. \code{"adjoint_solution"} for
#'       a requested \code{"auto"}), and \code{grad_method_requested}. The
#'       run record copies it.}
#'     \item{\code{run_record}}{A \code{dynhr_run_record}: resolved arguments,
#'       option snapshot, RNG state and provenance; replay it with
#'       \code{\link{dynhr_rerun}}.}
#'   }
#'
#' @examples
#' \dontrun{
#' mod  <- solve_model("my_model.mod")
#' mode <- run_mode_finding(mod, data, obs_vars = c("y", "pi", "r"))
#'
#' # Single method
#' post <- run_posterior_estimation(mode, methods = "RWMH", n_draws = 50000)
#'
#' # Multi-method dispatch
#' post <- run_posterior_estimation(mode,
#'   methods = c("RWMH", "HMC", "NUTS", "RWMH"),
#'   n_warmup = c(10000,  1000,  500,    500),
#'   n_draws  = 50000)
#' }
#'
#' @seealso \code{\link{solve_model}}, \code{\link{run_mode_finding}},
#'   \code{run_all_diagnostics}, \code{\link{dynhr_mcmc}}, \code{\link{nuts}},
#'   \code{\link{smc}}
#' @export
run_posterior_estimation <- function(mode_result,
                                     methods               = c("RWMH"),
                                     n_warmup              = 5000L,
                                     n_draws               = 20000L,
                                     n_chains              = 4L,
                                     n_particles           = 2000L,
                                     n_walkers             = NULL,
                                     parallel              = FALSE,
                                     parallel_backend      = "mirai",
                                     n_cores               = NULL,
                                     analytic_grad         = TRUE,
                                     Sigma_prop            = NULL,
                                     run_stoch_simul       = TRUE,
                                     skip_mode_finding_check = FALSE,
                                     nuts_timeout_seconds  = 300L,
                                     transform_params      = NULL,
                                     rwmh_adapt_cov        = NULL,
                                     rwmh_n_blocks         = NULL,
                                     metric                = c("diagonal", "hessian", "warmup_dense", "monge", "whittle_fim",
                                                               "lowrank", "fisher_diag"),
                                     monge_alpha           = NULL,
                                     checkpoint_dir        = NULL,
                                     resume                = FALSE,
                                     on_mismatch           = "refuse",
                                     verbose               = TRUE,
                                     seed                  = NULL,
                                     ...) {
  ## Retired spellings (0.9.3 renames) would otherwise travel through `...`
  ## to a sampler and die there as an unclassed "unused argument".
  .dynhr_reject_retired_args("run_posterior_estimation", ...names(),
                             names(sys.call()), .dynhr_retired_count_args)
  ## Own the message epoch for this run: repeat-suppressed warnings
  ## (`.dynhr_warn(once = TRUE)`) are keyed within it and re-arm for the
  ## next run, and the close reports what it suppressed. A nested call
  ## inherits this epoch rather than opening a second one.
  .dynhr_run_epoch <- .dynhr_epoch("run_posterior_estimation")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)
  ## Run record (R/run-record.R): resolved args, option snapshot and RNG
  ## state at ENTRY -- before the body touches any argument or the RNG.
  .rr <- .dynhr_rr_begin("run_posterior_estimation", environment(), list(...))
  ## E5 C2: a thin wrapper. The arguments become an estimation spec whose mode
  ## stage is `mode_result` (validate_spec() holds the cross-field checks) and
  ## the one spec runner samples it (.est_sampler_stage()).
  spec <- as_estimation_spec(.rr$args, entry = "run_posterior_estimation")
  ## The runner records the gradient method(s) the sampler stages resolved
  ## (.dynhr_rr_finish(result = )).
  invisible(.run_estimation_impl(spec, rr = .rr))
}

## Retired argument spellings of the estimation runners -> the current name
## (brief 32 P3). The 0.9.3 renames made run_posterior_estimation()'s counts
## follow the sampler convention, and unified data / obs_vars / me_variance
## across the likelihood entry points.
.dynhr_retired_count_args <- c(
  nburn = "n_warmup", ndraws = "n_draws", nchains = "n_chains",
  nparticles = "n_particles", nwalkers = "n_walkers")
.dynhr_retired_data_args <- c(
  Y = "data", obs_names = "obs_vars", observables = "obs_vars",
  me_var = "me_variance", me_sd = "me_variance")

## Abort (class dynhr_error_retired_argument) when a retired spelling in
## `retired` was supplied. `dot_names`: ...names() of the entry point (a
## retired name forwarded by a wrapper's `...`); `call_names`:
## names(sys.call()) (a retired name R would otherwise PARTIAL-match to a
## formal, e.g. me_var -> me_variance, so it never reaches `...`).
.dynhr_reject_retired_args <- function(fn, dot_names, call_names, retired) {
  nm <- unique(c(dot_names, call_names))
  hit <- intersect(nm[!is.na(nm) & nzchar(nm)], names(retired))
  if (!length(hit)) return(invisible(NULL))
  note <- ifelse(hit == "me_sd", ", which is a VARIANCE: pass me_sd^2", "")
  .dynhr_abort(fn, ": retired argument name",
               if (length(hit) > 1L) "s" else "", " ",
               paste0("`", hit, "` (use `", retired[hit], "`", note, ")",
                      collapse = ", "),
               ".", class = "dynhr_error_retired_argument")
}

## Upper-case run_posterior_estimation() method key of a sampler spec
## ("pmmh" is random-walk MH over a particle likelihood: key "RWMH").
.est_method_key <- function(s)
  if (identical(s$method, "pmmh")) "RWMH" else toupper(s$method)

## Particle-MCMC (PMMH) loglik-variance preflight at the mode, for the noisy
## unbiased-loglik filters: the TPF ("tpf") and the OBC particle filters
## ("ppf"/"copf"). .tpf_pmcmc_preflight is filter-agnostic (it evaluates the
## loglik K times with varied external seeds), so the same SD < 1
## (Herbst-Schorfheide) rule and N-recommendation apply. NULL when the
## likelihood is not a particle filter or the preflight is switched off.
.est_pmcmc_preflight <- function(log_post_fn, theta_mode, par_ctx, n_particles,
                                 verbose) {
  if (!par_ctx$likelihood %in% c("tpf", "ppf", "copf")) return(NULL)
  .vcat <- function(...) if (verbose) .dynhr_cat(...)
  ## TPF reads tpf_options; OBC PFs read obc_options (with the same
  ## pmcmc_preflight_* keys). Particle count lives in n_particles (TPF) or
  ## N (OBC); fall back to the sampler's `n_particles`.
  tpf_opts  <- if (identical(par_ctx$likelihood, "tpf"))
                 par_ctx$tpf_options %||% list()
               else par_ctx$obc_options %||% list()
  pf_K      <- tpf_opts$pmcmc_preflight_K    %||% 30L
  pf_skip   <- isTRUE(tpf_opts$pmcmc_preflight_skip)
  pf_n_part <- tpf_opts$n_particles %||% tpf_opts$N %||% n_particles
  if (pf_skip || pf_K <= 0L || is.null(log_post_fn) || is.null(theta_mode))
    return(NULL)
  .vcat(sprintf("-- TPF PMCMC preflight (K = %d) --\n", pf_K))
  ## .tpf_pmcmc_preflight calls set.seed(k) externally before each
  ## replicate. That only VARIES the filter noise when the closure was
  ## built with seed = NULL: the preflight detects the degenerate fixed-
  ## seed case and stops with an explicit "fixed seed" error.
  res <- .tpf_pmcmc_preflight(log_post_fn, theta_mode, K = pf_K,
                              verbose = verbose)
  if (!is.na(res$sd) && res$sd > 1) {
    n_needed <- ceiling(pf_n_part * res$n_needed_factor)
    .dynhr_warn(sprintf(
      "TPF loglik SD at mode = %.2f > 1 (Dynare threshold).", res$sd),
      "\nPMCMC acceptance dominated by loglik noise.",
      sprintf("\nCurrent n_particles = %.0f; to achieve SD < 1, raise to ~%.0f.",
              as.numeric(pf_n_part), as.numeric(n_needed)),
      call. = FALSE)
    res$n_needed <- n_needed
  } else if (!is.na(res$sd)) {
    .vcat(sprintf("  TPF loglik SD = %.3f (< 1 threshold, OK)\n", res$sd))
  }
  res
}

## The sampler stage of the spec runner (formerly the body of
## run_posterior_estimation()): runs the spec's sampler(s) in sequence from
## `mode_result`, pools the draws, computes convergence diagnostics and
## (outputs$stoch_simul) IRFs and moments at the posterior mean. Returns the
## dynhr_posterior_result without its run record.
.est_sampler_stage <- function(spec, mode_result) {
  samplers <- .spec_sampler_list(spec$sampler)
  cmp      <- spec$compute
  verbose  <- cmp$verbose
  parallel <- isTRUE(cmp$parallel)
  parallel_backend <- cmp$backend
  n_cores  <- cmp$n_cores
  seed     <- cmp$seed
  checkpoint_dir <- cmp$checkpoint_dir

  # Seed the WHOLE stage (chain-dispersal draws and every sampler's RNG)
  # before anything else touches the RNG stream. The parallel (mirai) paths
  # take the seed as their seed_base (E5 C2: they used to ignore it).
  if (!is.null(seed)) set.seed(seed)

  .vcat <- function(...) if (verbose) .dynhr_cat(...)

  solved      <- mode_result$solved
  model       <- solved$model
  compiled    <- solved$compiled
  log_post_fn <- mode_result$log_post_fn
  theta_mode  <- mode_result$theta_mode
  prior_spec  <- mode_result$prior_spec
  ## The proposal covariance: a sampler's own Sigma_prop, else the mode's.
  .sigma_of <- function(s) s$Sigma_prop %||% mode_result$Sigma_prop
  Sigma_prop  <- .sigma_of(Filter(function(s) "Sigma_prop" %in% names(s),
                                  samplers)[1L][[1L]])

  # ---- Checkpoint / restart (opt-in). When `checkpoint_dir` is set, each
  # chain streams its draws to chain_<id>.* files (RAM bounded by the flush
  # window, parallel-safe) and saves a restart state; `resume = TRUE` continues a
  # prior run with more draws. The config fingerprint forces the same parameters.
  sampler_checkpoint <- if (!is.null(checkpoint_dir)) {
    if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
    list(dir         = checkpoint_dir,
         flush_every = as.integer(.dynhr_opt("checkpoint_flush_every", default = 1000L)),
         resume      = isTRUE(cmp$resume),
         fingerprint = .ckpt_fingerprint(prior_spec$name, prior_spec))
  } else NULL

  # Opt-in unconstrained-parameter transform (eta-space sampling), per
  # sampler (sampler$transform_params). RWMH / NUTS / MALA / ChEES / HMC then
  # operate on eta = to_unconstrained(theta) (with the change-of-variables
  # Jacobian in the target), and results are mapped back to theta-space.
  # Sigma_prop (theta-space, from mode-finding) is converted to an ETA-SPACE
  # covariance via the delta method (.cov_theta_to_eta(), R/param-transform.R):
  #   Sigma_eta = D^{-1} Sigma_theta D^{-1}, D = diag(dtheta_deta(eta_mode)),
  # D floored at 1e-12 in absolute value. Shared by the serial and parallel
  # (mirai) paths. Default ON: in the constrained geometry of bounded
  # parameters near a unit root the gradient samplers fail -- NZSIM
  # re-measure 2026-09-25 (0.9.3.87, corrected mass adaptation, NUTS diag,
  # 300+300): 246/300 divergent, min bulk ESS 1.3 without the transform vs
  # 0 divergent, min ESS 57 with it. (The older "NUTS 5 draws, all divergent,
  # step 5.9e-11" evidence ran on the pre-0.9.3.50 inverted mass adaptation.)
  #
  # Mode on a prior bound (a KKT-active constrained mode): the delta method
  # divides a bound parameter's variance by (dtheta/deta)^2 ~ (distance to
  # the bound)^2, i.e. by ~1e-20 -> an eta proposal ~1e6+ too wide. There
  # run_mode_finding()'s Step 6 returns the correct eta-space proposal as
  # mode_result$Sigma_prop_eta (R/mode-hessian.R, "Step-6 proposal
  # covariance at a bound-constrained mode"); .bound_eta_of() hands it to
  # every theta->eta proposal / mass conversion below, provided the sampler
  # has no Sigma_prop of its own (the stored matrix belongs to the mode's
  # proposal) and the conversion point is the mode itself (a later stage
  # starts from the previous stage's last draw: the delta method there).
  # Interior modes carry no Sigma_prop_eta: every path is unchanged.
  .bound_eta_of <- function(s, theta_at) {
    Se <- mode_result$Sigma_prop_eta
    if (is.null(Se) || !is.null(s$Sigma_prop) ||
        !identical(theta_at, theta_mode)) return(NULL)
    Se
  }
  .transform_of <- function(s) {
    if (!isTRUE(s$transform_params)) return(list(pt = NULL, Sig_eta = NULL))
    pt <- build_param_transform(prior_spec, names(theta_mode))
    list(pt = pt,
         Sig_eta = .bound_eta_of(s, theta_mode) %||%
           .cov_theta_to_eta(.sigma_of(s), pt, theta_mode))
  }

  # Read the estimation context from mode_result; fall back to reconstructing
  # from legacy flat fields when mode_result pre-dates the ctx refactor.
  par_ctx <- if (!is.null(mode_result$ctx) &&
                 inherits(mode_result$ctx, "dynhr_estimation_context")) {
    mode_result$ctx
  } else {
    ctx_from_mode_result(mode_result)
  }

  # Inputs needed to provision a mirai daemon pool. The pool recompiles the
  # posterior per daemon via the standard Gaussian make_log_posterior, so it is
  # only valid when the mode result used neither OBC/PKF nor a non-standard
  # likelihood; otherwise the already-built closure is shipped.
  par_obs_vars    <- mode_result$obs_vars
  par_data        <- mode_result$data
  par_me_var      <- par_ctx$me_variance
  par_me_extra    <- par_ctx$me_extra
  par_shock_scale <- par_ctx$shock_scale
  par_standard    <- .ctx_is_standard_gaussian(par_ctx,
                       use_obc = !is.null(mode_result$obc_specs),
                       model   = model)
  ## Analytic gradients (NUTS / HMC / MALA / ChEES): analytic_grad (default
  ## TRUE) takes the exact gradient where the likelihood has one; where it has
  ## none the stage quietly falls back to numerical (the verbose output says
  ## so, below).
  grad_ok <- is.null(mode_result$obc_specs) &&
    .ctx_allows_analytic_gradient(par_ctx)

  # ------------------------------------------------------------------
  # Mode-finding quality check
  # ------------------------------------------------------------------
  # Run before any sampler starts.  Catches:
  #   (a) Non-finite logpost at mode     -> optimiser never converged
  #   (b) NULL / all-NA theta_mode       -> no starting point for MCMC
  #   (c) Sigma_prop is identity-scaled  -> Hessian inversion failed and
  #       a fallback diagonal was used; NUTS / HMC will stall or diverge
  #
  # Set skip_mode_finding_check = TRUE (mode$skip_check) to continue despite
  # these issues (e.g. a hand-specified Sigma_prop, or for debugging).
  if (!isTRUE(spec$mode$skip_check)) {

    # mode$logpost is the canonical field; mode$value is a legacy alias
    mode_logpost <- mode_result$mode$logpost %||%
                    mode_result$mode$value   %||%
                    mode_result$logpost      %||% NA_real_
    theta_mode_ok <- !is.null(theta_mode) && length(theta_mode) > 0L &&
                     any(is.finite(theta_mode))

    issues <- character(0)

    if (!is.finite(mode_logpost)) {
      issues <- c(issues, sprintf(
        "mode$value is %s (optimiser did not converge or Kalman filter returned -Inf/NaN).",
        format(mode_logpost)
      ))
    }
    if (!theta_mode_ok) {
      issues <- c(issues, "theta_mode is NULL, empty, or all-NA (no valid starting point).")
    }

    # Sigma_prop quality: check if it looks like a scaled identity (Hessian fallback)
    if (!is.null(Sigma_prop) && is.matrix(Sigma_prop) && nrow(Sigma_prop) > 1L) {
      diag_vals <- diag(Sigma_prop)
      off_diag  <- Sigma_prop[upper.tri(Sigma_prop)]
      # Identity-scaled fallback: all off-diagonal = 0, all diagonal equal
      is_scaled_identity <- (max(abs(off_diag)) < .Machine$double.eps * 1e4) &&
                            (diff(range(diag_vals)) / max(diag_vals, 1e-300) < 1e-10)
      if (is_scaled_identity) {
        issues <- c(issues,
          "Sigma_prop appears to be a scaled identity matrix: Hessian inversion failed and a fallback proposal was used. NUTS/HMC will likely stall.")
      }
      # Also warn if condition number is huge
      sv <- tryCatch(svd(Sigma_prop, nu = 0L, nv = 0L)$d, error = function(e) NULL)
      if (!is.null(sv) && min(sv) > 0) {
        kappa <- max(sv) / min(sv)
        if (kappa > 1e8) {
          issues <- c(issues, sprintf(
            "Sigma_prop condition number = %.2e (> 1e8): proposal covariance is nearly singular. NUTS may stall.",
            kappa
          ))
        }
      }
    }

    if (length(issues) > 0L) {
      msg <- paste0(
        "run_posterior_estimation: mode-finding check failed.\n",
        paste0("  - ", issues, collapse = "\n"),
        "\nRe-run run_mode_finding() with a different method or data, then retry.\n",
        "Pass skip_mode_finding_check = TRUE to force MCMC despite these issues."
      )
      stop(msg, call. = FALSE)
    }
  }

  .vcat("\n================================================================\n")
  .vcat("  run_posterior_estimation\n")
  .vcat("================================================================\n\n")

  # -------------------------------------------------------------------
  # Step 0.5: particle-MCMC (PMMH) loglik-variance preflight
  # Run before sampling so the SD warning appears early. Applies to any
  # noisy-unbiased-loglik filter used for PMMH: the TPF ("tpf") and the OBC
  # particle filters ("ppf"/"copf"). .tpf_pmcmc_preflight is filter-agnostic
  # (it just evaluates the loglik K times with varied external seeds), so the
  # same SD < 1 (Herbst-Schorfheide) rule and N-recommendation apply.
  # -------------------------------------------------------------------
  n_particles_1 <- Filter(Negate(is.null),
                          lapply(samplers, `[[`, "n_particles"))[1L][[1L]] %||% 2000L
  tpf_preflight_result <- .est_pmcmc_preflight(log_post_fn, theta_mode, par_ctx,
                                               n_particles_1, verbose)

  # -------------------------------------------------------------------
  # Step 1: The sampler sequence
  # -------------------------------------------------------------------
  n_methods <- length(samplers)
  unsupported <- vapply(samplers, function(s) identical(s$method, "smc2"),
                        logical(1))
  if (any(unsupported))
    .dynhr_abort("run_estimation: sampler \"smc2\" is not run by the spec runner ",
                 "yet; call smc2() directly.",
                 class = "dynhr_error_spec_unsupported")
  methods    <- vapply(samplers, .est_method_key, character(1))
  int_field  <- function(f) vapply(samplers, function(s)
    if (is.null(s[[f]])) NA_integer_ else as.integer(s[[f]]), integer(1))
  nburn_vec   <- int_field("n_warmup")
  ndraws_vec  <- int_field("n_draws")
  nchains_vec <- int_field("n_chains")

  .vcat(sprintf("  Methods: %s\n", paste(methods, collapse = " -> ")))
  .vcat(sprintf("  Burn-in: %s\n", paste(nburn_vec, collapse = ", ")))
  .vcat(sprintf("  Draws  : %s\n", paste(ndraws_vec, collapse = ", ")))
  .vcat(sprintf("  Chains : %s\n\n", paste(nchains_vec, collapse = ", ")))

  # -------------------------------------------------------------------
  # Step 2: Run samplers sequentially
  # -------------------------------------------------------------------
  chains_list <- list()
  current_theta <- theta_mode

  ## The gradient method each gradient-based stage actually used. "auto" (the
  ## grad_method option's default) is resolved by make_posterior_grad() once,
  ## at build time, and stamped on the closure as attr(, "grad_method"); it is
  ## printed here and kept on the stage result (chains[[m]]$resolved) and in
  ## the result's $resolved, which the run record copies.
  ##
  ## W94: the gradient is of the SAMPLED target, log_post_fn = the mode
  ## stage's objective (run_mode_finding() builds it from the likelihood
  ## spec's system prior and power_posterior, and its Step-5 gradient takes
  ## both, brief 28 S1). Before W94 neither reached this gradient: with a
  ## system prior NUTS / ChEES / MALA / HMC sampled log_post_fn with the
  ## gradient of another posterior (valid -- the MH step uses log_post_fn --
  ## but inefficient), and W92 refused the fused log-posterior there. Now the
  ## fused value IS log_post_fn's (the samplers still check it at the start
  ## point, .hmc_fused_target()). The infeasibility penalty is the one piece
  ## make_posterior_grad() cannot take: with it the objective is a finite
  ## penalty where the gradient's fused value is -Inf, so the companion is
  ## dropped and the sampler keeps separate calls.
  .stage_grad <- function(tag) {
    .vcat(sprintf("  [%s] building analytic gradient (%s)...\n", tag, grad_method))
    g <- make_posterior_grad(model, par_data, prior_spec, par_obs_vars,
                             compiled, me_variance = par_me_var,
                             me_extra    = par_me_extra,
                             shock_scale = par_shock_scale,
                             grad_method = grad_method,
                             likelihood  = par_ctx$likelihood,
                             freq_band   = par_ctx$freq_band,
                             lik_init    = par_ctx$lik_init %||% "auto",
                             power         = spec$likelihood$power_posterior,
                             system_priors = par_ctx$system_priors)
    .vcat(sprintf("  [%s] gradient method: %s\n", tag,
                  .grad_method_label(attr(g, "grad_method"), grad_method)))
    if (!is.null(par_ctx$infeasible_penalty) ||
        !is.null(spec$likelihood$infeasible_penalty))
      attr(g, "logpost_grad") <- NULL
    g
  }

  for (i in seq_len(n_methods)) {
    s  <- samplers[[i]]
    m  <- methods[i]
    nb <- nburn_vec[i]
    nd <- ndraws_vec[i]
    nc <- nchains_vec[i]
    np <- s$n_particles
    transform_params <- isTRUE(s$transform_params)
    tr <- .transform_of(s)
    param_transform <- tr$pt
    Sigma_prop_eta  <- tr$Sig_eta
    Sigma_prop      <- .sigma_of(s)
    metric          <- s$metric
    analytic_grad   <- isTRUE(s$analytic_grad) && grad_ok
    grad_method     <- s$grad_method
    ## set by the MALA / NUTS / ChEES / HMC branches when they build a gradient
    stage_grad_method <- NULL
    ## run_posterior_estimation()'s `...`, forwarded to the sampler
    sampler_args    <- s$extra

    .vcat(sprintf("--- Method %d/%d: %s (%d burn + %d draws) ---\n",
                  i, n_methods, m, nb, nd))
    if (m %in% c("NUTS", "HMC", "MALA", "CHEES") && !analytic_grad)
      .vcat(sprintf("  [%s] gradient: numerical (%s)\n", m,
                    if (!isTRUE(s$analytic_grad)) "analytic_grad = FALSE"
                    else "this likelihood has no analytic gradient"))

    chain_res <- switch(m,
      "RWMH" = {
        ## CPM dispatch: when cpm_rho_u is set and likelihood = "tpf", use
        ## rwmh_cpm (serial path only). Otherwise fall through to .run_rwmh_batch.
        ## log_post_fn here is the closure built at mode-finding time; its
        ## burn_in_init must already be 0L (the mode stage pins it whenever
        ## tpf_options$cpm_rho_u is set) since rwmh_cpm calls it with a non-NULL
        ## U_list on every step after the priming call.
        cpm_rho_u <- (par_ctx$tpf_options %||% list())$cpm_rho_u
        use_cpm   <- !is.null(cpm_rho_u) && is.numeric(cpm_rho_u) &&
                     is.finite(cpm_rho_u) && cpm_rho_u > 0 && cpm_rho_u < 1 &&
                     identical(par_ctx$likelihood, "tpf")
        if (use_cpm) {
          ## CPM preflight message
          if (!is.null(tpf_preflight_result) && !is.na(tpf_preflight_result$sd)) {
            sd_ll     <- tpf_preflight_result$sd
            sd_cpm_eff <- sd_ll * sqrt(2 * (1 - cpm_rho_u))
            .dynhr_inform(sprintf(
              "  CPM enabled (rho_u = %.2f). Effective loglik SD ~ %.3f (from %.3f).",
              cpm_rho_u, sd_cpm_eff, sd_ll))
            if (sd_ll < 0.5)
              .dynhr_warn("CPM requested but loglik SD = ", round(sd_ll, 3),
                      " < 0.5. Standard PMCMC already efficient; CPM overhead not needed.",
                      call. = FALSE)
          }
          if (nc > 1L)
            .dynhr_inform("  CPM: parallel chains not yet supported; running serial (nc=1).")
          ## transform_params: rwmh_cpm samples eta (with the Jacobian) using
          ## the eta-space Sigma_prop_eta, as .run_rwmh_batch does. W86: the
          ## transform used not to be passed, so the eta covariance drove a
          ## THETA-space random walk (sig near 0.02: proposals ~1/d = 50x too
          ## wide).
          cpm_Sigma <- if (transform_params && !is.null(Sigma_prop_eta)) Sigma_prop_eta else Sigma_prop
          cpm_res <- rwmh_cpm(
            log_post_fn  = log_post_fn,
            theta0       = current_theta,
            Sigma_prop   = cpm_Sigma,
            transform    = if (transform_params) param_transform,
            n_draws      = nd,
            n_burn       = nb,
            rho_u        = cpm_rho_u,
            verbose      = verbose,
            chain_id     = 1L
          )
          list(chains = list(cpm_res), chain_stats = NULL)
        } else {
        # transform_params: sample in eta-space using the eta-space
        # Sigma_prop_eta (delta-method conversion above); .run_rwmh_batch
        # maps the starting point on both the serial and parallel branches.
        rwmh_Sigma <- if (transform_params && !is.null(Sigma_prop_eta)) Sigma_prop_eta else Sigma_prop
        do.call(.run_rwmh_batch, c(list(
                                log_post_fn, current_theta, rwmh_Sigma,
                                prior_spec, n_chains = nc,
                                n_draws = nd, n_burn = nb,
                                verbose = verbose,
                                parallel = parallel,
                                parallel_backend = parallel_backend,
                                n_cores = n_cores,
                                parsed_model = if (par_standard) model else NULL,
                                Y = if (par_standard) par_data else NULL,
                                obs_names = par_obs_vars,
                                ctx = par_ctx,
                                me_variance = par_me_var,
                                me_extra = par_me_extra,
                                shock_scale = par_shock_scale,
                                use_closure = !par_standard,
                                transform = param_transform,
                                adapt_cov = isTRUE(s$adapt_cov),
                                n_blocks = s$n_blocks,
                                checkpoint = sampler_checkpoint),
                                if (!is.null(seed)) list(seed_base = seed),
                                sampler_args))
        }  # end else (non-CPM RWMH path)
      },
      "HMC"  = {
        ## Analytic gradient (W94: the same gate and builder as NUTS / MALA /
        ## ChEES; HMC used to get the numerical gradient whatever was asked).
        hmc_grad <- if (analytic_grad) .stage_grad("HMC") else NULL
        stage_grad_method <- attr(hmc_grad, "grad_method")
        ## metric: "warmup_dense" and "diagonal" are the HMC metrics.
        do.call(.run_hmc_batch, c(list(log_post_fn, current_theta,
                       n_draws = nd, n_warmup = nb,
                       metric = metric,
                       # parity with NUTS/MALA: run in eta-space so the
                       # leapfrog does not step-collapse against bounds
                       # (adapt_mass then tunes the diagonal mass in
                       # eta-space). dynhr_hmc maps draws back to theta.
                       transform = if (transform_params) param_transform else NULL,
                       verbose = verbose),
                       ## only when built: a grad_fn in `...` still reaches it
                       if (!is.null(hmc_grad)) list(grad_fn = hmc_grad),
                       sampler_args))
      },
      "MALA" = {
        # MALA: preconditioned Langevin with constant or position-dependent metric.
        # metric = "hessian": use G_inv = Sigma_prop (same as NUTS dense path).
        # metric = "diagonal": G = G_inv = identity (plain MALA).
        # metric = "monge":    position-dependent Monge metric (Stage 3a):
        #   G(theta) = I + alpha^2 * g*g', passed as metric_fn to dynhr_mala.
        mala_G_inv     <- NULL
        mala_G         <- NULL
        mala_metric_fn <- NULL

        if (identical(metric, "hessian") && !is.null(Sigma_prop) &&
            is.matrix(Sigma_prop) && all(is.finite(Sigma_prop))) {
          # G_inv is the covariance of the SAMPLER's coordinates: in eta-space
          # under transform_params (dynhr_mala runs on eta there) -- see
          # .dense_metric_basis(). Before W86 the theta-space Sigma_prop was
          # handed over as the eta metric: off by d_i d_j per entry.
          mala_db <- .dense_metric_basis(
            Sigma_prop, if (transform_params) param_transform, current_theta,
            if (transform_params && !is.null(param_transform))
              .bound_eta_of(s, current_theta))
          mala_G_inv_try <- mala_db$B  # Sigma_prop is the inverse-mass (M_inv)
          mala_G_try     <- tryCatch(solve(mala_G_inv_try), error = function(e) NULL)
          if (!is.null(mala_G_try)) {
            mala_G_inv <- mala_G_inv_try / tcrossprod(mala_db$dv)  # D^-1 B D^-1
            mala_G     <- mala_G_try * tcrossprod(mala_db$dv)      # D B^-1 D
            .vcat("  [MALA] using dense metric from Sigma_prop.\n")
          } else {
            .vcat("  [MALA] metric='hessian': Sigma_prop not invertible -- using identity.\n")
          }
        }

        # Analytic gradient (same gate as NUTS); also needed for the Monge metric_fn
        mala_grad <- if (analytic_grad) .stage_grad("MALA") else NULL
        stage_grad_method <- attr(mala_grad, "grad_method")

        # Monge metric: position-dependent G(theta) = I + alpha^2 g g'.
        # Requires a gradient function (analytic if available, else FD); the
        # metric_fn needs its own copy because it is called from .metric_at(),
        # which does NOT pass through grad_fn. validate_spec() gates it on the
        # allow_monge_metric option and warns about its step collapse.
        if (identical(metric, "monge")) {
          alpha_monge <- s$monge_alpha
          # If analytic grad is available reuse it; otherwise build FD closure.
          monge_grad_fn <- if (!is.null(mala_grad)) {
            mala_grad
          } else {
            par_names_local <- names(theta_mode)
            lp_for_grad <- function(th) {
              names(th) <- par_names_local
              res <- log_post_fn(th)
              val <- if (is.list(res)) res$logpost else res
              if (!is.finite(val)) return(-1e300)
              val
            }
            function(th) .hmc_gradient(lp_for_grad, th, method = "forward")
          }
          mala_metric_fn <- monge_metric_fn(monge_grad_fn, alpha = alpha_monge)
          .vcat(sprintf("  [MALA] Monge metric (alpha=%.4g).\n", alpha_monge))
        }

        do.call(.run_mala_batch, c(list(log_post_fn, current_theta,
                        n_draws    = nd,   n_warmup = nb,
                        G          = mala_G,
                        G_inv      = mala_G_inv,
                        metric_fn  = mala_metric_fn,
                        grad_fn    = mala_grad,
                        transform  = if (transform_params) param_transform else NULL,
                        chain_id   = 1L, checkpoint = sampler_checkpoint,
                        verbose    = verbose), sampler_args))
      },
      "NUTS" = {
        # Parallel multi-chain NUTS on a mirai pool. Standard Gaussian models
        # recompile the posterior per daemon; OBC/PKF and cumulant models ship
        # the already-built log_post_fn closure once instead.
        use_par_nuts <- parallel && nc > 1L &&
                        identical(parallel_backend, "mirai") &&
                        requireNamespace("mirai", quietly = TRUE)
        if (use_par_nuts) {
          # run_nuts_mirai() forwards the ADAPTED metrics ("diagonal",
          # "warmup_dense", "fisher_diag", "lowrank"); a fixed dense metric
          # ("hessian", "whittle_fim") is not shipped to daemons --
          # validate_spec() says so (dynhr_warning_metric_ignored).
          # Analytic/implicit gradient on the parallel path requires a
          # standard Gaussian model (par_standard): each daemon needs
          # .worker_model/.worker_cm/.worker_Y from .mirai_pool_init, which
          # only runs when parsed_model/Y are supplied (par_standard).
          par_analytic_grad <- analytic_grad && par_standard
          if (isTRUE(s$analytic_grad) && grad_ok && !par_standard)
            .vcat("  [NUTS] analytic_grad requires a standard Gaussian model (compiled per daemon); ",
                  "this OBC/PKF or cumulant run uses the numerical gradient.\n")
          else if (par_analytic_grad) {
            ## Each daemon builds its own closure with the same inputs (and the
            ## same default theta_ref), so it resolves "auto" the same way;
            ## resolve it here too, once, to report and record it.
            stage_grad_method <- .grad_method_resolved(
              model, par_data, prior_spec, par_obs_vars, compiled,
              grad_method = grad_method, likelihood = par_ctx$likelihood,
              me_extra = par_me_extra, shock_scale = par_shock_scale,
              lik_init = par_ctx$lik_init %||% "auto",
              me_variance = par_me_var)
            .vcat(sprintf("  [NUTS] parallel chains will build analytic gradients (%s) per daemon...\n",
                          .grad_method_label(stage_grad_method, grad_method)))
          }
          # The covariance behind the chains' initial mass (1/diag) and the
          # chain-2..N dispersion: ONE precedence, the serial NUTS branch's
          # (W87) --
          #   1. the sampler's own Sigma_prop,
          #   2. at a bound mode (no own Sigma_prop, chain start = the mode):
          #      mode_result$Sigma_prop_eta, already in eta-space,
          #   3. the mode's Sigma_prop (theta-space; delta-method to eta
          #      under transform_params).
          # All three carry the RWMH scale 2.38^2 / n_par. This used to take
          # mode_result$V_mode (UNSCALED) ahead of the sampler's own
          # Sigma_prop, and at a bound mode mixed it with the SCALED
          # Sigma_prop_eta. The mass scale is not neutral for NUTS: it moves
          # the initial step size and hence the whole adapted run. V_mode,
          # rescaled to the same 2.38^2 / n_par, is only the last resort when
          # the mode stored no Sigma_prop (the parallel path needs a matrix;
          # the serial branch then runs the identity mass).
          nuts_Sig <- Sigma_prop %||%
            (if (is.matrix(mode_result$V_mode))
               mode_result$V_mode * (2.38^2 / length(current_theta)))
          if (is.null(rownames(nuts_Sig)))
            rownames(nuts_Sig) <- colnames(nuts_Sig) <- names(current_theta)
          # transform_params: run_nuts_mirai's `Sigma_prop` (and hence its
          # mass_diag = 1/diag(Sigma_prop)) must be in ETA-SPACE -- convert
          # nuts_Sig via the same delta-method helper used for Sigma_prop_eta
          # above, evaluated at the current chain start `current_theta`
          # (mirrors the serial NUTS branch's `eta_mode_nuts`).
          # At a bound mode: mode_result$Sigma_prop_eta (see .bound_eta_of).
          if (transform_params && !is.null(param_transform)) {
            nuts_Sig_eta <- .bound_eta_of(s, current_theta) %||%
              .cov_theta_to_eta(nuts_Sig, param_transform, current_theta)
            if (!is.null(nuts_Sig_eta)) nuts_Sig <- nuts_Sig_eta
          }
          par_res <- do.call(run_nuts_mirai, c(list(
            parsed_model = if (par_standard) model else NULL,
            Y = if (par_standard) par_data else NULL,
            prior_spec = prior_spec,
            obs_names = par_obs_vars, theta_mode = current_theta,
            Sigma_prop = nuts_Sig, n_chains = nc, n_draws = nd, n_warmup = nb,
            n_cores = n_cores, me_variance = par_me_var,
            me_extra = par_me_extra,
            shock_scale = par_shock_scale,
            ctx = par_ctx,
            log_post_fn = if (par_standard) NULL else log_post_fn,
            transform = if (transform_params) param_transform else NULL,
            analytic_grad = par_analytic_grad,
            grad_method = grad_method,
            checkpoint = sampler_checkpoint,
            metric = if (metric %in% c("warmup_dense", "fisher_diag", "lowrank"))
                       metric else "diagonal",
            adapt = s$adapt,
            progress = verbose),
            if (!is.null(seed)) list(seed_base = seed)))
          list(chains = par_res$chains, chain_stats = par_res$chain_stats)
        } else {
        # NUTS warmup can stall indefinitely on steep / poorly-conditioned
        # posteriors (step-size adaptation keeps halving until it hits eps=0).
        # Wrap with setTimeLimit so a runaway warmup is caught cleanly.
        # Pre-condition with Sigma_prop: mass_diag = 1/diag(Sigma_prop) rescales
        # each parameter to unit posterior variance, fixing tiny initial step_size.
        nuts_mass <- if (!is.null(Sigma_prop) && is.matrix(Sigma_prop) &&
                         all(is.finite(diag(Sigma_prop))) &&
                         all(diag(Sigma_prop) > 0)) {
          1 / pmax(diag(Sigma_prop), 1e-12)
        } else NULL
        # transform_params: mass_diag must be supplied in ETA-SPACE. nuts_mass
        # above is the MASS 1/Var(theta); by the delta method
        # Var(eta) ~ Var(theta) / (dtheta/deta)^2, so the eta mass is
        #   m_eta = 1/Var(eta) = m_theta * dtheta_deta(eta_mode)^2
        # (0.9.3.131: this used to DIVIDE by d^2 -- the mass off by d^4 per
        # parameter, e.g. ~1e8 for a log-transformed parameter near 0.01).
        # dtheta_deta is floor-guarded (>= 1e-12 in absolute value) before
        # squaring to avoid blow-up near a transform's asymptote.
        # At a bound mode the eta-space diagonal is mode_result$Sigma_prop_eta's
        # (see .bound_eta_of): the delta method is ~1e6+ too wide there.
        Se_nuts <- if (transform_params && !is.null(param_transform))
          .bound_eta_of(s, current_theta)
        if (!is.null(Se_nuts) && !is.null(nuts_mass)) {
          nuts_mass <- 1 / pmax(diag(Se_nuts), 1e-12)
        } else if (transform_params && !is.null(param_transform) && !is.null(nuts_mass)) {
          eta_mode_nuts <- param_transform$to_unconstrained(current_theta)
          d_vec_nuts <- param_transform$dtheta_deta(eta_mode_nuts)
          d_vec_nuts[abs(d_vec_nuts) < 1e-12] <- 1e-12
          nuts_mass <- nuts_mass * d_vec_nuts^2
        }
        # --- Dense metric (metric = "hessian", and "whittle_fim"'s fallback) ---
        # Use Sigma_prop as the dense inverse-mass matrix M_inv, IN THE
        # SAMPLER'S SPACE. Sigma_prop is theta-space (mode-finding); under
        # transform_params dynhr_nuts runs on eta, so M_inv must be the eta
        # covariance D^-1 Sigma_prop D^-1 (D = diag(dtheta/deta) at the chain
        # start; Sigma_prop_eta at a bound mode) -- .dense_metric_basis().
        # W86: this path used to hand the THETA-space Sigma_prop over as the
        # eta M_inv (off by d_i d_j per entry, e.g. ~1e3 for a std near 0.02).
        # The inverse / Cholesky are taken on the theta-space matrix and scaled
        # analytically (M = D Sigma^-1 D, chol(M) = chol(Sigma^-1) D), as the
        # "whittle_fim" branch below does: the D^2 spread would otherwise be
        # folded into the condition number that solve()/chol() see.
        # The diagonal path uses 1/diag(Sigma_prop); this is the strict
        # generalisation to the full matrix. We do NOT use hessian_exact
        # directly.
        # metric = "whittle_fim" builds the same matrix first, silently: it is
        # that metric's fallback, already in the sampler's space (W87: the
        # fallback used to be unreachable -- the message said 'hessian' while
        # the run silently got the diagonal metric).
        nuts_M_inv  <- NULL
        nuts_chol_M <- NULL
        nuts_hess_msg <- identical(metric, "hessian")
        if (metric %in% c("hessian", "whittle_fim") && !is.null(Sigma_prop) &&
            is.matrix(Sigma_prop) && all(is.finite(Sigma_prop))) {
          # B is the inverse-mass up to the diagonal rescaling by dv.
          # M = solve(M_inv); chol_M = chol(M).
          nuts_db <- .dense_metric_basis(
            Sigma_prop, if (transform_params) param_transform, current_theta,
            if (transform_params && !is.null(param_transform))
              .bound_eta_of(s, current_theta))
          nuts_M_inv_try <- nuts_db$B
          nuts_M_try     <- tryCatch(solve(nuts_M_inv_try), error = function(e) NULL)
          if (!is.null(nuts_M_try)) {
            nuts_chol_M_try <- tryCatch(chol(nuts_M_try), error = function(e) NULL)
            if (!is.null(nuts_chol_M_try)) {
              nuts_M_inv  <- nuts_M_inv_try / tcrossprod(nuts_db$dv)       # D^-1 B D^-1
              nuts_chol_M <- nuts_chol_M_try *
                rep(nuts_db$dv, each = nrow(nuts_chol_M_try))              # chol(B^-1) D
              nuts_mass   <- NULL  # dense path overrides diagonal
              if (nuts_hess_msg)
                .vcat("  [NUTS] using dense mass matrix from Sigma_prop.\n")
            } else if (nuts_hess_msg) {
              .vcat("  [NUTS] metric='hessian': Sigma_prop not PD -- falling back to diagonal.\n")
            }
          } else if (nuts_hess_msg) {
            .vcat("  [NUTS] metric='hessian': Sigma_prop not invertible -- falling back to diagonal.\n")
          }
        }

        # --- One-shot dense metric from the Whittle FIM (metric = "whittle_fim") ---
        # Opt-in: validate_spec() errors unless allow_whittle_fim_metric = TRUE.
        # Assembled ONCE at the mode (current_theta), then frozen through warmup
        # and sampling -- same contract as "hessian". On any failure
        # (non-standard model, dss_list assembly, degenerate FIM, non-positive
        # transform Jacobian) the run keeps the "hessian" dense metric built
        # above (D^-1 Sigma_prop D^-1 under transform_params; Sigma_prop_eta at
        # a bound mode), or the diagonal metric when Sigma_prop gave none, and
        # says which. The fallback is NEVER routed through whittle_fim(): a
        # returned fallback_metric would go on to the theta->eta D-scaling below
        # and an eta-space matrix would be converted twice. whittle_fim() gets a
        # sentinel instead, which it returns as-is when the FIM is degenerate
        # (.whittle_fim_fallback); without it a degenerate FIM silently became
        # the identity metric.
        if (identical(metric, "whittle_fim")) {
          .wf_fallback <- function(why)
            .vcat(sprintf("  [NUTS] metric='whittle_fim': %s -- falling back to %s.\n",
                          why, if (!is.null(nuts_M_inv))
                                 "the 'hessian' metric (Sigma_prop)"
                               else "the diagonal metric"))
          if (!par_standard) {
            .wf_fallback("requires a standard Gaussian model")
          } else {
            T_obs_wf <- ncol(par_data)
            omega_grid_wf <- 2 * pi * seq_len(floor(T_obs_wf / 2)) / T_obs_wf
            asm <- tryCatch(
              .assemble_dss_list_at_mode(model, compiled, current_theta, par_obs_vars),
              error = function(e) .dynhr_reraise_bug(e, NULL)
            )
            if (is.null(asm) || !isTRUE(asm$ok)) {
              .wf_fallback("dss_list assembly failed")
            } else {
              fim <- tryCatch(
                whittle_fim(TT = asm$TT, RR = asm$RR, ZZ = asm$ZZ, DD = asm$DD,
                           Sigma_e = asm$Sigma_e, dss_list = asm$dss_list,
                           omega_grid = omega_grid_wf, T_obs = T_obs_wf,
                           prior_hess = NULL,
                           fallback_metric = list(degenerate = TRUE),
                           me_variance = par_me_var %||% 0),
                error = function(e) .dynhr_reraise_bug(e, NULL)
              )
              if (isTRUE(fim$degenerate)) {
                .wf_fallback("degenerate Whittle FIM")
              } else if (is.null(fim) || !is.matrix(fim$G_inv)) {
                .wf_fallback("whittle_fim() failed")
              } else {
                ## fim$G is the theta-space precision (metric); the dense NUTS
                ## mass matrix M IS the metric G (M_inv = G^{-1}), matching how
                ## the "hessian" branch feeds nuts_M_inv/nuts_chol_M (there
                ## Sigma_prop plays the role of M_inv directly).
                ##
                ## Transform to eta-space ANALYTICALLY from the factors
                ## whittle_fim() already returns (G_inv, L = chol(G), both
                ## SoftAbs-stabilised): with D = diag(dtheta_deta) > 0,
                ##   M      = G_eta        = D G D
                ##   M^{-1} = G_eta^{-1}   = D^{-1} G^{-1} D^{-1}
                ##   chol(G_eta)           = L D   (upper-tri x diagonal stays
                ##                                  upper-tri; diag > 0 for d > 0)
                ## Numerically re-solving/chol-ing G_eta instead is NOT viable
                ## on stiff posteriors: cond(G) ~ 1e11 on NZSIM-class models and
                ## the D^2 spread pushes cond(G_eta) past double precision --
                ## exactly the failure the 2026-07-03 NZSIM shootout hit.
                d_wf <- rep(1, nrow(fim$G))
                if (transform_params && !is.null(param_transform)) {
                  eta_wf <- param_transform$to_unconstrained(current_theta)
                  d_wf   <- param_transform$dtheta_deta(eta_wf)
                }
                if (all(is.finite(d_wf)) && all(d_wf > 0) &&
                    is.matrix(fim$G_inv) && is.matrix(fim$L)) {
                  nuts_M_inv  <- fim$G_inv * tcrossprod(1 / d_wf)   # D^-1 G^-1 D^-1
                  nuts_chol_M <- fim$L %*% diag(d_wf, nrow = length(d_wf))
                  dimnames(nuts_M_inv)  <- dimnames(fim$G)
                  nuts_mass   <- NULL
                  .vcat("  [NUTS] using dense mass matrix from the Whittle FIM (one-shot at the mode).\n")
                } else {
                  .wf_fallback("non-positive/non-finite transform Jacobian or missing FIM factors")
                }
              }
            }
          }
        }

        # Exact analytic gradient (compiled Kalman score for shock-std params +
        # numerical for the rest, or full implicit-differentiation gradient
        # when grad_method = "implicit") when requested on a standard Gaussian
        # model.
        nuts_grad <- if (analytic_grad) .stage_grad("NUTS") else NULL
        stage_grad_method <- attr(nuts_grad, "grad_method")
        t0_nuts <- proc.time()
        nuts_sec <- max(30L, as.integer(s$timeout))
        res <- tryCatch({
          setTimeLimit(elapsed = nuts_sec, transient = TRUE)
          do.call(dynhr_nuts, c(list(log_post_fn, current_theta,
                     n_draws = nd, n_warmup = nb,
                     mass_diag = nuts_mass, grad_fn = nuts_grad,
                     M_inv = nuts_M_inv, chol_M = nuts_chol_M,
                     ## Adapted metrics reach dynhr_nuts by name; the fixed
                     ## dense ones ("hessian"/"whittle_fim") arrive as
                     ## M_inv/chol_M above.
                     metric = if (metric %in% c("diagonal", "warmup_dense", "lowrank", "fisher_diag"))
                                metric else "diagonal",
                     transform = if (transform_params) param_transform else NULL,
                     chain_id = 1L, checkpoint = sampler_checkpoint), sampler_args))
        }, error = function(e) {
          elapsed <- round((proc.time() - t0_nuts)[["elapsed"]], 1)
          if (grepl("reached elapsed time limit|time limit|timed out",
                    conditionMessage(e), ignore.case = TRUE)) {
            .vcat(sprintf(
              "\n  [NUTS] Timed out after %.1fs (nuts_timeout_seconds = %d).\n",
              elapsed, nuts_sec))
            .vcat("  This usually means NUTS warmup is diverging (step-size -> 0).\n")
            .vcat("  Cause: ill-conditioned Sigma_prop (Hessian fallback) or very\n")
            .vcat("  steep / non-smooth posterior.  Try RWMH or SMC instead.\n")
            NULL
          } else {
            stop(e)   # re-raise unexpected errors
          }
        }, finally = {
          setTimeLimit(elapsed = Inf, transient = FALSE)
        })
        if (is.null(res)) {
          # Timeout path: return an empty result so the loop continues
          list(chains    = list(),
               chain_stats = data.frame(chain = integer(0), accept_rate = numeric(0),
                                        final_logpost = numeric(0),
                                        stringsAsFactors = FALSE))
        } else {
          list(chains = list(res), chain_stats = data.frame(
            chain = 1L, accept_rate = res$acceptance_rate %||% NA,
            final_logpost = tail(res$post_logpost %||% rep(NA, nd), 1),
            stringsAsFactors = FALSE
          ))
        }
        }  # end else (serial single-chain NUTS)
      },
      "SMC" = {
        # Standard Gaussian models recompile the posterior per daemon; OBC/PKF
        # and cumulant models ship the already-built log_post_fn closure once.
        use_par_smc <- parallel &&
                       identical(parallel_backend, "mirai") &&
                       requireNamespace("mirai", quietly = TRUE)
        res <- if (use_par_smc)
          do.call(run_smc_mirai, c(list(
                        parsed_model = if (par_standard) model else NULL,
                        Y = if (par_standard) par_data else NULL,
                        prior_spec = prior_spec, obs_names = par_obs_vars,
                        n_particles = np, n_cores = n_cores,
                        me_variance = par_me_var,
                        me_extra = par_me_extra,
                        shock_scale = par_shock_scale,
                        ctx = par_ctx,
                        log_post_fn = if (par_standard) NULL else log_post_fn,
                        verbose = verbose),
                        if (!is.null(seed)) list(seed_base = seed),
                        sampler_args))
        else
          do.call(dynhr_smc, c(list(log_post_fn, prior_spec = prior_spec,
                    n_particles = np), sampler_args))
        ## Resample the weighted particle cloud to an equally-weighted draw
        ## matrix so all downstream consumers (diagnostics, Bayesian IRF,
        ## smoother) receive a standard draw matrix.  The original particles
        ## and weights remain in $particles and $smc_weights.
        ## B5: re-index $post_logpost (and every other per-particle vector)
        ## with the SAME resampling indices as $chain, so THAMES and any other
        ## (draw, log-posterior) consumer pairs each draw with its own value.
        if (.is_smc_weighted(res)) {
          res <- .smc_equal_weight_result(res)
          if (verbose)
            .vcat(sprintf("  SMC: resampled %d particles to %d equally-weighted draws\n",
                          res$n_particles, res$n_draws))
        }
        list(chains = list(res), chain_stats = data.frame(
          chain = 1L, accept_rate = NA,
          final_logpost = NA,
          stringsAsFactors = FALSE
        ))
      },

      "DSMH" = {
        # Dynamic Striated MH (Waggoner-Wu-Zha 2016; Dynare 'dsmh'). With
        # parallel = TRUE on the mirai backend each stage's chain groups run on
        # a pool (dynhr_dsmh(n_cores =); seeded output is identical to serial).
        # Its cloud is already equally weighted (MH output), so no
        # resampling step. Dynare's first rung lambda1 = 1/(10 * n_par * T)
        # takes T from the data unless the caller supplied lambda1 /
        # lambda_schedule / n_obs through the sampler's extra arguments.
        dsmh_args <- sampler_args
        if (parallel && is.null(dsmh_args[["n_cores"]])) {
          if (identical(parallel_backend, "mirai") &&
              requireNamespace("mirai", quietly = TRUE))
            dsmh_args[["n_cores"]] <- .mirai_n_cores(n_cores, dsmh_args[["n_groups"]] %||% 10L)
          else
            .dynhr_inform("  DSMH: parallel runs need the mirai backend; running serially.")
        }
        if (is.null(dsmh_args[["n_obs"]]) && is.null(dsmh_args[["lambda1"]]) &&
            is.null(dsmh_args[["lambda_schedule"]]) && !is.null(par_data))
          dsmh_args[["n_obs"]] <- NROW(par_data)
        if (is.null(dsmh_args[["verbose"]])) dsmh_args[["verbose"]] <- verbose
        res <- do.call(dynhr_dsmh,
                       c(list(log_post_fn, prior_spec = prior_spec,
                              n_particles = np), dsmh_args))
        list(chains = list(res), chain_stats = data.frame(
          chain = 1L, accept_rate = res$accept_schedule[res$n_stages],
          final_logpost = NA,
          stringsAsFactors = FALSE
        ))
      },

      "CHEES" = {
        # ChEES-HMC (Hoffman, Radul & Sountsov 2021): fixed-length leapfrog
        # with dual-averaged step size AND trajectory-time adaptation via the
        # ChEES criterion.  Mirrors the NUTS serial branch.
        chees_mass <- if (!is.null(Sigma_prop) && is.matrix(Sigma_prop) &&
                          all(is.finite(diag(Sigma_prop))) &&
                          all(diag(Sigma_prop) > 0)) {
          1 / pmax(diag(Sigma_prop), 1e-12)
        } else NULL
        # transform_params: mass_diag in eta-space (same delta-method as NUTS;
        # at a bound mode mode_result$Sigma_prop_eta, see .bound_eta_of)
        Se_chees <- if (transform_params && !is.null(param_transform))
          .bound_eta_of(s, current_theta)
        if (!is.null(Se_chees) && !is.null(chees_mass)) {
          chees_mass <- 1 / pmax(diag(Se_chees), 1e-12)
        } else if (transform_params && !is.null(param_transform) && !is.null(chees_mass)) {
          eta_mode_chees <- param_transform$to_unconstrained(current_theta)
          d_vec_chees <- param_transform$dtheta_deta(eta_mode_chees)
          d_vec_chees[abs(d_vec_chees) < 1e-12] <- 1e-12
          chees_mass <- chees_mass * d_vec_chees^2   # the eta MASS (see the NUTS branch)
        }
        # Analytic gradient (same gate as NUTS)
        chees_grad <- if (analytic_grad) .stage_grad("ChEES") else NULL
        stage_grad_method <- attr(chees_grad, "grad_method")
        res <- do.call(dynhr_chees, c(list(log_post_fn, current_theta,
                           n_draws   = nd, n_warmup = nb,
                           mass_diag = chees_mass, grad_fn = chees_grad,
                           transform = if (transform_params) param_transform else NULL,
                           chain_id  = 1L, checkpoint = sampler_checkpoint,
                           verbose   = verbose), sampler_args))
        list(chains = list(res), chain_stats = data.frame(
          chain = 1L, accept_rate = res$acceptance_rate %||% NA,
          final_logpost = tail(res$post_logpost %||% rep(NA, nd), 1),
          stringsAsFactors = FALSE
        ))
      },

      "DIME" = {
        use_par_dime <- parallel &&
                        identical(parallel_backend, "mirai") &&
                        requireNamespace("mirai", quietly = TRUE)
        res <- if (use_par_dime)
          do.call(run_dime_mirai, c(list(
            parsed_model = if (par_standard) model else NULL,
            Y            = if (par_standard) par_data else NULL,
            prior_spec   = prior_spec, obs_names = par_obs_vars,
            n_chain      = s$n_walkers,
            n_iter       = nd, n_burn = nb,
            n_cores      = n_cores,
            me_variance  = par_me_var,
            me_extra     = par_me_extra,
            shock_scale  = par_shock_scale,
            ctx          = par_ctx,
            log_post_fn  = if (par_standard) NULL else log_post_fn,
            verbose      = verbose),
            if (!is.null(seed)) list(seed_base = seed)))
        else
          do.call(run_dime, c(list(log_post_fn, prior_spec = prior_spec,
                   n_chain = s$n_walkers,
                   n_iter  = nd, n_burn = nb,
                   checkpoint = sampler_checkpoint,
                   verbose = verbose), sampler_args))
        list(chains = list(res), chain_stats = data.frame(
          chain = 1L, accept_rate = res$acceptance_rate,
          final_logpost = tail(res$post_logpost[is.finite(res$post_logpost)], 1L),
          stringsAsFactors = FALSE
        ))
      }
    )

    if (!is.null(stage_grad_method))
      chain_res$resolved <- list(grad_method = stage_grad_method,
                                 grad_method_requested = grad_method)
    chains_list[[m]] <- chain_res

    # Update starting point for next method
    if (!is.null(chain_res$chains) && length(chain_res$chains) > 0) {
      last_chain <- chain_res$chains[[length(chain_res$chains)]]
      if (!is.null(last_chain$chain) && nrow(last_chain$chain) > 0)
        current_theta <- setNames(as.numeric(last_chain$chain[nrow(last_chain$chain), ]),
                                  colnames(last_chain$chain))
    }

    # Print chain stats
    if (verbose && nrow(chain_res$chain_stats) > 0) {
      print(chain_res$chain_stats, row.names = FALSE, digits = 3)
    }
    .vcat("\n")
  }

  # -------------------------------------------------------------------
  # Step 3: Pool draws and compute convergence
  # -------------------------------------------------------------------
  .vcat("-- Aggregating chains --\n")

  ## Pool each STORED run once. chains_list is keyed by method name, so a
  ## repeated method (methods = c("RWMH", "NUTS", "RWMH")) keeps only its
  ## last run; iterating `methods` itself visited that run twice and
  ## double-counted its draws.
  all_chains <- list()
  for (m in names(chains_list)) {
    if (!is.null(chains_list[[m]]$chains)) {
      for (ch in chains_list[[m]]$chains) {
        if (!is.null(ch$chain)) all_chains <- c(all_chains, list(ch$chain))
      }
    }
  }

  if (length(all_chains) > 0) {
    pooled_draws <- do.call(rbind, all_chains)
    .vcat(sprintf("  Pooled draws: %d x %d\n", nrow(pooled_draws), ncol(pooled_draws)))
  } else {
    pooled_draws <- NULL
    .vcat("  No valid chains produced.\n")
  }

  # Convergence diagnostics
  conv <- NULL
  if (length(all_chains) >= 2) {
    .vcat("  Computing convergence diagnostics...\n")
    conv <- .compute_convergence(all_chains)
    if (!is.null(conv$rhat)) {
      .vcat(sprintf("  R-hat: max=%.3f  (>1.05: %d/%d)\n",
                    max(conv$rhat, na.rm = TRUE),
                    sum(conv$rhat > 1.05, na.rm = TRUE),
                    length(conv$rhat)))
      .vcat(sprintf("  ESS:   median=%.0f  min=%.0f\n",
                    median(conv$ess, na.rm = TRUE),
                    min(conv$ess, na.rm = TRUE)))
    }
  }

  # Posterior mean
  posterior_mean <- if (!is.null(pooled_draws)) {
    setNames(colMeans(pooled_draws), colnames(pooled_draws))
  } else {
    theta_mode
  }

  # -------------------------------------------------------------------
  # Step 4: Stoch_simul at posterior mean
  # -------------------------------------------------------------------
  irfs    <- NULL
  moments <- NULL
  if (isTRUE(spec$outputs$stoch_simul)) {
    .vcat("-- Stoch_simul at posterior mean --\n")

    ## 0.9.4 (ledger A1): RE-SOLVE at the posterior mean through the
    ## likelihood's own pipeline (apply_theta_to_params() handles the
    ## shock-named entries), never the calibration-time `solved$dr`. If the
    ## posterior mean does not solve, both stay NULL.
    solve_state <- new.env(parent = emptyenv())
    solve_state$ss_warm <- NULL
    ## lik_init = "diffuse": moments/IRFs are well defined at a unit-root
    ## posterior mean, so do not apply the stationary-likelihood rejection.
    post_sol <- .solve_dr_for_theta(model, compiled,
                                    cache_system_structure(compiled),
                                    posterior_mean, solve_state,
                                    lik_init = "diffuse")
    if (is.null(post_sol)) {
      .dynhr_warn("run_posterior_estimation: the model does not solve at the ",
                  "posterior mean (steady state or Blanchard-Kahn failure); ",
                  "posterior_irfs and posterior_moments are NULL.")
    } else {
      dr <- post_sol$dr
      ## Clear the solve-time Sigma_e snapshot so `params` stays authoritative
      ## (the 0.9.3.7 asymmetry: compute_irfs() prefers dr$Sigma_e).
      dr$Sigma_e <- NULL
      params_post <- post_sol$params
      irfs <- compute_irfs(dr, model, params = params_post)
      moments <- compute_moments(dr, model, params = params_post)
      if (!is.null(irfs))
        .vcat(sprintf("  IRFs computed: %d shocks\n", length(irfs)))
      if (!is.null(moments))
        .vcat("  Theoretical moments computed.\n")
    }
  }

  # -------------------------------------------------------------------
  # Wrap up
  # -------------------------------------------------------------------
  .vcat("\n================================================================\n")
  .vcat("  Posterior estimation complete\n")
  .vcat("================================================================\n\n")

  ## Per-stage resolved gradient methods (stages that built an analytic
  ## gradient only), keyed like `chains`; the run record copies it.
  gm_res <- Filter(Negate(is.null), lapply(chains_list, `[[`, "resolved"))
  resolved <- list(
    grad_method = vapply(gm_res, `[[`, character(1), "grad_method"),
    grad_method_requested = vapply(gm_res, `[[`, character(1),
                                   "grad_method_requested"))

  result <- list(
    mode_result       = mode_result,
    chains            = chains_list,
    resolved          = resolved,
    pooled_draws      = pooled_draws,
    convergence       = conv,
    posterior_mean    = posterior_mean,
    posterior_irfs    = irfs,
    posterior_moments = moments,
    tpf_preflight     = tpf_preflight_result,
    meta = list(
      methods       = methods,
      n_warmup      = nburn_vec,
      n_draws       = ndraws_vec,
      n_chains      = nchains_vec,
      n_methods     = n_methods,
      seed          = seed   # provenance for the report (NULL when unseeded)
    )
  )
  class(result) <- c("dynhr_posterior_result", "list")
  result
}


#' Print method for dynhr_posterior_result
#' @noRd
#' @export
print.dynhr_posterior_result <- function(x, ...) {
  cat("\n<dynhr_posterior_result>\n")
  cat(sprintf("  Methods   : %s\n", paste(x$meta$methods, collapse = " -> ")))
  cat(sprintf("  Draws     : %d pooled\n",
              if (!is.null(x$pooled_draws)) nrow(x$pooled_draws) else 0))
  cat(sprintf("  Params    : %d estimated\n", nrow(x$mode_result$prior_spec)))
  if (!is.null(x$convergence$rhat)) {
    rh <- x$convergence$rhat
    cat(sprintf("  R-hat     : max = %.3f  (>1.05: %d / %d)\n",
                max(rh, na.rm = TRUE), sum(rh > 1.05, na.rm = TRUE), length(rh)))
  }
  if (!is.null(x$convergence$ess)) {
    es <- x$convergence$ess
    cat(sprintf("  ESS       : median = %.0f  min = %.0f\n",
                median(es, na.rm = TRUE), min(es, na.rm = TRUE)))
  }
  gl <- .grad_method_record_line(x$resolved)
  if (!is.null(gl)) cat(sprintf("  Gradient  : %s\n", gl))
  if (!is.null(x$posterior_irfs))
    cat(sprintf("  IRFs      : %d shocks at posterior mean\n", length(x$posterior_irfs)))
  if (!is.null(x$posterior_moments))
    cat("  Moments   : computed at posterior mean\n")
  il <- .est_integrity_lines(x$provenance$integrity)
  if (length(il)) cat(il, sep = "\n")
  invisible(x)
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Run a batch of RWMH chains
#'
#' When \code{parallel = TRUE} (and the necessary model inputs are supplied),
#' the chains run on a persistent mirai daemon pool via
#' \code{run_mcmc_mirai()}; otherwise they run sequentially.
#'
#' Two pool-provisioning paths: when \code{parsed_model}/\code{Y}/
#' \code{obs_names} are supplied, daemons recompile the standard Gaussian
#' posterior once via \code{.mirai_pool_init()}. When \code{use_closure = TRUE}
#' (OBC/PKF or cumulant models), the already-built \code{log_post_fn} closure
#' is shipped once via \code{.mirai_pool_closure()} instead.
#'
#' @param transform Optional "dynhr_param_transform" (opt-in). When non-NULL,
#'   \code{Sigma_prop} is interpreted as an eta-space covariance (the
#'   delta-method conversion from the theta-space proposal covariance is the
#'   caller's responsibility -- see \code{run_posterior_estimation}'s
#'   \code{Sigma_prop_eta}) and forwarded to \code{rwmh()}'s (serial path) or
#'   \code{run_mcmc_mirai()}'s (PARALLEL/mirai path) \code{transform} argument;
#'   extra chain starting points are dispersed in eta-space and mapped back
#'   via \code{to_constrained()} on both paths.
#' @param adapt_cov Opt-in (default \code{FALSE}); forwarded to
#'   \code{rwmh()}'s (serial path) or \code{run_mcmc_mirai()}'s
#'   (PARALLEL/mirai path) \code{adapt_cov} argument -- Haario et al. (2001)
#'   adaptive proposal covariance.
#' @param n_blocks Opt-in (default \code{1L}); forwarded to \code{rwmh()}'s
#'   (serial path) or \code{run_mcmc_mirai()}'s (PARALLEL/mirai path)
#'   \code{n_blocks} argument -- randomized parameter blocking.
#' @noRd
.run_rwmh_batch <- function(log_post_fn, theta_mode, Sigma_prop,
                             prior_spec, n_chains = 4L,
                             n_draws = 20000L, n_burn = 5000L,
                             verbose = TRUE,
                             parallel = FALSE, parallel_backend = "mirai",
                             n_cores = NULL, parsed_model = NULL, Y = NULL,
                             obs_names = NULL, me_variance = 0,
                             me_extra = NULL,
                             shock_scale = NULL,
                             use_closure = FALSE,
                             transform = NULL,
                             adapt_cov = FALSE, n_blocks = 1L,
                             seed_base = NULL, ctx = NULL, checkpoint = NULL, ...) {

  n_par     <- length(theta_mode)
  L <- .robust_chol(Sigma_prop, n_par)

  # Parallel path: a single mirai pool runs all chains. Either the parsed
  # model + data (standard Gaussian, recompiled per daemon) or a pre-built
  # log-posterior closure (OBC/PKF, cumulant; shipped once) is required.
  use_par <- isTRUE(parallel) && n_chains > 1L &&
             identical(parallel_backend, "mirai") &&
             requireNamespace("mirai", quietly = TRUE) &&
             (isTRUE(use_closure) ||
              (!is.null(parsed_model) && !is.null(Y) && !is.null(obs_names)))
  if (use_par) {
    if (verbose) .dynhr_cat(sprintf("    Parallel RWMH (mirai): %d chains\n", n_chains))
    par_res <- run_mcmc_mirai(
      parsed_model = parsed_model, Y = Y,
      prior_spec   = prior_spec, obs_names = obs_names,
      theta_mode   = theta_mode, Sigma_prop = Sigma_prop,
      n_chains     = n_chains, n_draws = n_draws, n_burn = n_burn,
      seed_base    = seed_base, n_cores = n_cores,
      me_variance  = me_variance,
      me_extra     = me_extra,
      shock_scale  = shock_scale,
      ctx          = ctx,
      log_post_fn  = if (isTRUE(use_closure)) log_post_fn else NULL,
      transform    = transform,
      adapt_cov    = adapt_cov,
      n_blocks     = n_blocks,
      checkpoint   = checkpoint,
      progress     = verbose)
    chain_list <- par_res$chains
    conv <- .compute_convergence(chain_list)
    combined <- if (!is.null(conv)) {
      new_dynhr_chains(list(
        chain           = conv$combined,
        acceptance_rate = mean(par_res$chain_stats$accept_rate, na.rm = TRUE),
        sampler         = "rwmh", n_chains = n_chains,
        chain_list      = chain_list, chain_stats = par_res$chain_stats), "rwmh")
    } else NULL
    return(list(chains = chain_list, chain_stats = par_res$chain_stats,
                combined = combined, convergence = conv))
  }

  chain_list  <- vector("list", n_chains)
  chain_stats <- data.frame(
    chain = integer(), accept_rate = numeric(),
    final_logpost = numeric(), stringsAsFactors = FALSE
  )

  for (ch in seq_len(n_chains)) {
    if (ch == 1L) {
      th0 <- theta_mode
    } else if (!is.null(transform)) {
      # Disperse in eta-space (Sigma_prop here is Sigma_prop_eta) and map
      # back to theta-space; eta has no boundary, so no clamping is needed.
      eta_mode <- transform$to_unconstrained(theta_mode)
      z   <- rnorm(n_par)
      eta0 <- eta_mode + 0.4 * as.numeric(L %*% z)
      names(eta0) <- names(theta_mode)
      th0 <- transform$to_constrained(eta0)
      if (!is.finite(log_post_fn(th0)$logpost)) th0 <- theta_mode
    } else {
      z   <- rnorm(n_par)
      th0 <- theta_mode + 0.4 * as.numeric(L %*% z)
      names(th0) <- names(theta_mode)
      for (i in seq_along(th0)) {
        th0[i] <- max(th0[i], prior_spec$lower[i] + 1e-8)
        th0[i] <- min(th0[i], prior_spec$upper[i] - 1e-8)
      }
      if (!is.finite(log_post_fn(th0)$logpost)) th0 <- theta_mode
    }

    if (verbose) .dynhr_cat(sprintf("    Chain %d/%d...\n", ch, n_chains))
    chain_list[[ch]] <- rwmh(log_post_fn, th0, Sigma_prop,
                               n_draws = n_draws + n_burn, n_burn = n_burn,
                               transform = transform, chain_id = ch,
                               adapt_cov = adapt_cov, n_blocks = n_blocks,
                               checkpoint = checkpoint, verbose = verbose, ...)

    if (!is.null(chain_list[[ch]]))
      chain_stats <- rbind(chain_stats, data.frame(
        chain         = ch,
        accept_rate   = chain_list[[ch]]$acceptance_rate,
        final_logpost = tail(chain_list[[ch]]$post_logpost, 1),
        stringsAsFactors = FALSE
      ))
  }

  # Compute convergence
  conv <- .compute_convergence(chain_list)

  # Build combined result
  combined <- if (!is.null(conv)) {
    list(
      chain           = conv$combined,
      acceptance_rate = mean(chain_stats$accept_rate, na.rm = TRUE),
      sampler         = "RWMH",
      n_chains        = n_chains,
      chain_list      = chain_list,
      chain_stats     = chain_stats
    )
  } else {
    list(
      chain           = NULL,
      acceptance_rate = NA,
      sampler         = "RWMH",
      n_chains        = 0L,
      chain_list      = chain_list,
      chain_stats     = chain_stats
    )
  }
  combined <- new_dynhr_chains(combined, "rwmh")

  list(chains = chain_list, chain_stats = chain_stats,
       combined = combined, convergence = conv)
}


#' Dense sampler metric in the sampler's coordinates (theta or eta)
#'
#' The dense \code{"hessian"} metric of NUTS / MALA is a covariance (inverse
#' mass). Under \code{transform_params} the sampler runs on
#' \eqn{\eta}, so the metric must be the eta covariance
#' \eqn{D^{-1} \Sigma_\theta D^{-1}}, \eqn{D = diag(d\theta/d\eta)} at
#' \code{theta_at} (floored at 1e-12 in absolute value, as
#' \code{.cov_theta_to_eta()}); at a mode on a prior bound it is the Step-6
#' eta proposal \code{Sigma_eta} itself. Returns the pieces instead of the
#' product so callers can invert / factor the theta-space matrix and rescale
#' analytically: \eqn{M^{-1} = B / (d d^T)}, \eqn{M = B^{-1} (d d^T)},
#' \eqn{chol(M) = chol(B^{-1}) D}.
#' @param Sigma_theta Theta-space covariance (Sigma_prop).
#' @param pt Parameter transform, or NULL (theta-space sampling).
#' @param theta_at Point at which D is evaluated (the chain start).
#' @param Sigma_eta Eta-space covariance at a bound mode, or NULL.
#' @return list(B = matrix, dv = numeric scaling vector; all 1 when no
#'   rescaling applies).
#' @noRd
.dense_metric_basis <- function(Sigma_theta, pt = NULL, theta_at = NULL,
                                Sigma_eta = NULL) {
  n <- nrow(Sigma_theta)
  if (!is.null(Sigma_eta)) return(list(B = Sigma_eta, dv = rep(1, n)))
  if (is.null(pt)) return(list(B = Sigma_theta, dv = rep(1, n)))
  dv <- as.numeric(pt$dtheta_deta(pt$to_unconstrained(theta_at)))
  small <- abs(dv) < 1e-12
  dv[small] <- sign(dv[small]) * 1e-12
  dv[dv == 0] <- 1e-12
  list(B = Sigma_theta, dv = dv)
}


#' Run a batch of HMC chains
#' @noRd
.run_hmc_batch <- function(log_post_fn, theta_mode,
                            n_draws = 2000L, n_warmup = 1000L,
                            verbose = TRUE, ...) {
  if (!exists("dynhr_hmc", mode = "function")) {
    if (verbose) .dynhr_cat("    HMC not available, falling back to NUTS...\n")
    return(.run_rwmh_batch(log_post_fn, theta_mode,
                            diag(length(theta_mode)), NULL,
                            n_chains = 1L, n_draws = n_draws,
                            n_burn = n_warmup, verbose = verbose, ...))
  }

  if (verbose) .dynhr_cat("    Running HMC...\n")
  res <- dynhr_hmc(log_post_fn, theta_mode,
                   n_draws = n_draws, n_warmup = n_warmup, ...)
  list(chains = list(res), chain_stats = data.frame(
    chain = 1L,
    accept_rate = res$acceptance_rate %||% NA,
    final_logpost = tail(res$post_logpost %||% rep(NA, n_draws), 1),
    stringsAsFactors = FALSE
  ))
}


#' Compute convergence diagnostics from a list of chain matrices
#' @noRd
.compute_convergence <- function(chains) {
  if (is.null(chains) || length(chains) == 0)
    return(list(rhat = NULL, ess = NULL, combined = NULL))

  # Handle chains that may be rwmh result lists or raw matrices
  chain_list <- lapply(chains, function(ch) {
    if (is.list(ch) && !is.null(ch$chain)) ch$chain else ch
  })

  valid <- which(!sapply(chain_list, is.null) &
                   sapply(chain_list, function(x) is.matrix(x) && nrow(x) > 0))
  if (length(valid) < 2) {
    combined <- if (length(valid) == 1) chain_list[[valid[1]]] else NULL
    return(list(rhat = NULL, ess = NULL, combined = combined))
  }

  chain_list <- chain_list[valid]
  combined <- do.call(rbind, chain_list)

  ## Chains of different lengths (a sampler sequence whose methods keep
  ## different numbers of draws, e.g. methods = c("RWMH", "NUTS") with
  ## n_draws = c(40, 30)): the split R-hat / ESS estimators need a common
  ## length, so each chain contributes its LAST n_min draws (its most
  ## converged part). `combined` still pools every draw. Before E5 C3 this
  ## crashed ("values must be length ..."). `n_used` records the length used.
  lens  <- vapply(chain_list, nrow, integer(1))
  n_min <- min(lens)
  if (any(lens != n_min))
    chain_list <- lapply(chain_list, function(m)
      m[seq.int(nrow(m) - n_min + 1L, nrow(m)), , drop = FALSE])

  ## One estimator for the whole package: the split-chain, rank-normalised,
  ## folded R-hat and the Geyer bulk ESS in diag-helpers.R (Vehtari et al.
  ## 2021), shared with D5 and sampler_diagnostics(). This used to be a third
  ## private implementation (plain unsplit R-hat, ESS on concatenated chains)
  ## that disagreed with what the diagnostics reported.
  p_names <- colnames(combined)
  if (is.null(p_names)) p_names <- paste0("theta_", seq_len(ncol(combined)))
  cs   <- .d5_convergence(chain_list, p_names)
  rhat <- setNames(cs$rhat,     cs$param)
  ess  <- setNames(cs$ess_bulk, cs$param)

  out <- list(rhat = rhat, ess = ess, combined = combined)
  if (any(lens != n_min)) out$n_used <- n_min
  out
}
