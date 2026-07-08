## R/run-posterior-estimation.R
## --------------------------------------------------------------------------
## run_posterior_estimation() -- Multi-method posterior sampling dispatcher.
##
## Given a mode-finding result (dynhr_mode_result), runs one or more MCMC
## samplers in sequence, with automatic multi-chain convergence diagnostics
## and stoch_simul at the posterior mean.
##
## Dispatch rule:
##   - Vector arguments (methods, nburn, ndraws, nchains) are paired with
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
#'   \code{"RWMH"}, \code{"HMC"}, \code{"NUTS"}, \code{"SMC"}.
#'   Default \code{c("RWMH")}.
#' @param nburn        Burn-in / warmup draws per method. Recycled to
#'   \code{length(methods)}.  Default \code{5000}.
#' @param ndraws       Post-warmup draws to retain per method. Recycled.
#'   Default \code{20000}.
#' @param nchains      Number of MCMC chains per method (for RWMH). Recycled.
#'   Default \code{4}.
#' @param nparticles   Number of SMC particles (only for \code{"SMC"}).
#'   Recycled.  Default \code{2000}.
#' @param nwalkers     Number of ensemble walkers for the DIME sampler
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
#'   (\code{NULL} = auto-detect, capped at \code{nchains}).
#' @param analytic_grad For \code{"NUTS"} on a standard Gaussian model, use the
#'   exact analytic gradient (\code{\link{make_posterior_grad}}) instead of a
#'   numerical one. Default \code{FALSE}. Applies to the serial NUTS path; the
#'   parallel multi-chain path uses the numerical gradient.
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
#'   eta-space covariance via the delta method, and the NUTS/ChEES
#'   \code{mass_diag} is similarly converted; draws are mapped back to
#'   theta-space before being returned, so downstream consumers (convergence
#'   diagnostics, stoch_simul, etc.) are unaffected. \strong{Strongly
#'   recommended for models with bounded parameters near a unit root}: in the
#'   constrained geometry every gradient sampler (MALA/ChEES/NUTS) otherwise
#'   suffers step-size collapse, which the eta-space transform resolves
#'   universally. May also be set globally via
#'   \code{dynhr_set_options(transform_params = TRUE)}.
#' @param rwmh_adapt_cov  Opt-in (default \code{FALSE}): forwarded to
#'   \code{rwmh}'s `adapt_cov` argument on the RWMH chains (serial and
#'   parallel/mirai paths) -- Haario et al. (2001) adaptive proposal
#'   covariance, frozen at the end of burn-in. May also be set globally via
#'   \code{dynhr_set_options(rwmh_adapt_cov = TRUE)}.
#' @param rwmh_n_blocks   Opt-in (default \code{1L}): forwarded to
#'   \code{rwmh}'s `n_blocks` argument on the RWMH chains (serial and
#'   parallel/mirai paths) -- randomized parameter blocking
#'   (Chib & Ramamurthy 2010 / Herbst & Schorfheide 2015 ch. 4). May also be
#'   set globally via \code{dynhr_set_options(rwmh_n_blocks = 2L)}.
#' @param metric       Metric for MALA and NUTS:
#'   \code{"diagonal"} (default, identity mass matrix),
#'   \code{"hessian"} (dense metric from \code{Sigma_prop}),
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
#'   \code{"hessian"} metric (Sigma_prop) with a message, mirroring the
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
#'   rather than \code{ndraws * n_par}, and a restart state is saved after each
#'   flush. The files are per-chain, so parallel chains never collide. Flush size
#'   is set via \code{dynhr_set_options(checkpoint_flush_every = ...)} (default 1000).
#' @param resume  When \code{TRUE} and \code{checkpoint_dir} points at a prior
#'   run, continue each chain from its saved state -- RNG, position, scale and
#'   proposal covariance are restored exactly, so the continuation is identical
#'   to a single longer run -- adding \code{ndraws} more retained draws. The
#'   saved model / prior / parameter configuration must match (it is enforced).
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
#'   gap.
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
#'   }
#'
#' @examples
#' \dontrun{
#' mod  <- solve_model("my_model.mod")
#' mode <- run_mode_finding(mod, data, obs_vars = c("y", "pi", "r"))
#'
#' # Single method
#' post <- run_posterior_estimation(mode, methods = "RWMH", ndraws = 50000)
#'
#' # Multi-method dispatch
#' post <- run_posterior_estimation(mode,
#'   methods = c("RWMH", "HMC", "NUTS", "RWMH"),
#'   nburn   = c(10000,  1000,  500,    500),
#'   ndraws  = 50000)
#' }
#'
#' @seealso \code{\link{solve_model}}, \code{\link{run_mode_finding}},
#'   \code{run_all_diagnostics}, \code{\link{mcmc}}, \code{\link{nuts}},
#'   \code{\link{smc}}
#' @export
run_posterior_estimation <- function(mode_result,
                                     methods               = c("RWMH"),
                                     nburn                 = 5000L,
                                     ndraws                = 20000L,
                                     nchains               = 4L,
                                     nparticles            = 2000L,
                                     nwalkers              = NULL,
                                     parallel              = FALSE,
                                     parallel_backend      = "mirai",
                                     n_cores               = NULL,
                                     analytic_grad         = FALSE,
                                     Sigma_prop            = NULL,
                                     run_stoch_simul       = TRUE,
                                     skip_mode_finding_check = FALSE,
                                     nuts_timeout_seconds  = 300L,
                                     transform_params      = NULL,
                                     rwmh_adapt_cov        = NULL,
                                     rwmh_n_blocks         = NULL,
                                     metric                = c("diagonal", "hessian", "warmup_dense", "monge", "whittle_fim"),
                                     monge_alpha           = NULL,
                                     checkpoint_dir        = NULL,
                                     resume                = FALSE,
                                     verbose               = TRUE,
                                     seed                  = NULL,
                                     ...) {

  # Seed the WHOLE run (mode-finding multistarts already happened upstream in
  # run_mode_finding(); this seeds chain-dispersal draws and every sampler's
  # RNG below) before anything else touches the RNG stream. Plain set.seed()
  # is used deliberately: withr is Suggests-only and not used anywhere else in
  # R/ (test-only), so we don't add a new runtime dependency for this.
  if (!is.null(seed)) set.seed(seed)

  .vcat <- function(...) if (verbose) cat(...)

  metric <- match.arg(metric)

  solved      <- mode_result$solved
  model       <- solved$model
  compiled    <- solved$compiled
  log_post_fn <- mode_result$log_post_fn
  theta_mode  <- mode_result$theta_mode
  prior_spec  <- mode_result$prior_spec
  Sigma_prop  <- Sigma_prop %||% mode_result$Sigma_prop

  # ------------------------------------------------------------------
  # Opt-in unconstrained-parameter transform (eta-space sampling)
  # ------------------------------------------------------------------
  # When enabled, RWMH / NUTS operate on eta = to_unconstrained(theta) (with
  # the change-of-variables Jacobian included in the target -- see
  # make_transformed_logpost()/dynhr_nuts()/rwmh()'s `transform` arg), and
  # results are mapped back to theta-space before being returned. Sigma_prop
  # (theta-space, from mode-finding) is converted to an ETA-SPACE covariance
  # via the delta method (.cov_theta_to_eta(), in R/param-transform.R):
  #
  #   eta = g(theta)  =>  Var(eta) ~ J Var(theta) J'  with J = diag(deta/dtheta)
  #                                 = diag(1/dtheta_deta) Var(theta) diag(1/dtheta_deta)
  #
  # i.e. Sigma_eta = D^{-1} Sigma_theta D^{-1}, D = diag(dtheta_deta(eta_mode)).
  # D entries are floored at 1e-12 in absolute value before inverting to
  # avoid blow-up when a parameter sits very near a transform's asymptote
  # (e.g. eta -> -Inf for a "log" transform, where dtheta/deta -> 0). This
  # conversion is shared by the serial and parallel (mirai) RWMH/NUTS paths.
  # Default ON: every gradient sampler (MALA/ChEES/NUTS) step-collapses in the
  # constrained geometry of bounded parameters near a unit root; the eta-space
  # transform fixes it universally and leaves the target invariant. Set
  # transform_params = FALSE (or dynhr_set_options) to sample in raw theta-space.
  transform_params <- isTRUE(.dynhr_opt("transform_params", transform_params,
                                         default = TRUE))

  # ---- RWMH opt-in upgrades (Haario adaptive covariance, randomized
  # blocking) -- resolved here and threaded through .run_rwmh_batch to both
  # the serial and parallel (mirai) RWMH paths. Defaults (FALSE / 1L) are
  # bit-identical to the pre-existing sampler.
  rwmh_adapt_cov <- isTRUE(.dynhr_opt("rwmh_adapt_cov", rwmh_adapt_cov,
                                       default = FALSE))
  rwmh_n_blocks  <- .dynhr_opt("rwmh_n_blocks", rwmh_n_blocks, default = 1L)

  # ---- Checkpoint / restart (opt-in, RWMH). When `checkpoint_dir` is set, each
  # chain streams its draws to chain_<id>.* files (RAM bounded by the flush
  # window, parallel-safe) and saves a restart state; `resume = TRUE` continues a
  # prior run with more draws. The config fingerprint forces the same parameters.
  sampler_checkpoint <- if (!is.null(checkpoint_dir)) {
    if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
    list(dir         = checkpoint_dir,
         flush_every = as.integer(.dynhr_opt("checkpoint_flush_every", default = 1000L)),
         resume      = isTRUE(resume),
         fingerprint = .ckpt_fingerprint(prior_spec$name, prior_spec))
  } else NULL

  param_transform <- NULL
  Sigma_prop_eta  <- NULL
  if (transform_params) {
    param_transform <- build_param_transform(prior_spec, names(theta_mode))
    Sigma_prop_eta  <- .cov_theta_to_eta(Sigma_prop, param_transform, theta_mode)
  }

  # Read the estimation context from mode_result; fall back to reconstructing
  # from legacy flat fields when mode_result pre-dates the ctx refactor.
  # Fixes: system_priors and freq_band silently dropped on parallel path.
  par_ctx <- if (!is.null(mode_result$ctx) &&
                 inherits(mode_result$ctx, "dynhr_estimation_context")) {
    mode_result$ctx
  } else {
    ctx_from_mode_result(mode_result)
  }

  # Inputs needed to provision a mirai daemon pool for parallel RWMH. The pool
  # recompiles the posterior per daemon via the standard Gaussian make_log_posterior,
  # so it is only valid when the mode result used neither OBC/PKF nor the
  # cumulant likelihood; otherwise .run_rwmh_batch falls back to the serial path.
  par_obs_vars    <- mode_result$obs_vars
  par_data        <- mode_result$data
  par_me_var      <- par_ctx$me_variance
  par_me_extra    <- par_ctx$me_extra
  par_shock_scale <- par_ctx$shock_scale
  par_standard    <- .ctx_is_standard_gaussian(par_ctx,
                       use_obc = !is.null(mode_result$obc_specs))

  # ------------------------------------------------------------------
  # Mode-finding quality check
  # ------------------------------------------------------------------
  # Run before any sampler starts.  Catches:
  #   (a) Non-finite logpost at mode     -> optimiser never converged
  #   (b) NULL / all-NA theta_mode       -> no starting point for MCMC
  #   (c) Sigma_prop is identity-scaled  -> Hessian inversion failed and
  #       a fallback diagonal was used; NUTS / HMC will stall or diverge
  #
  # Set skip_mode_finding_check = TRUE to continue despite these issues
  # (e.g. when using a hand-specified Sigma_prop, or for debugging).
  if (!isTRUE(skip_mode_finding_check)) {

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
  tpf_preflight_result <- NULL
  if (par_ctx$likelihood %in% c("tpf", "ppf", "copf")) {
    ## TPF reads tpf_options; OBC PFs read obc_options (with the same
    ## pmcmc_preflight_* keys). Particle count lives in n_particles (TPF) or
    ## N (OBC); fall back to the `nparticles` arg, then 1000.
    tpf_opts  <- if (identical(par_ctx$likelihood, "tpf"))
                   par_ctx$tpf_options %||% list()
                 else par_ctx$obc_options %||% list()
    pf_K      <- tpf_opts$pmcmc_preflight_K    %||% 30L
    pf_skip   <- isTRUE(tpf_opts$pmcmc_preflight_skip)
    pf_n_part <- tpf_opts$n_particles %||% tpf_opts$N %||% nparticles %||% 1000L

    if (!pf_skip && pf_K > 0L && !is.null(log_post_fn) &&
        !is.null(theta_mode)) {
      .vcat(sprintf("-- TPF PMCMC preflight (K = %d) --\n", pf_K))
      ## Use .tpf_pmcmc_preflight which calls set.seed(k) externally for each
      ## replicate, ensuring varied RNG even if the closure has a fixed seed.
      tpf_preflight_result <- .tpf_pmcmc_preflight(
        log_post_fn, theta_mode, K = pf_K, verbose = verbose)

      if (!is.na(tpf_preflight_result$sd) && tpf_preflight_result$sd > 1) {
        n_needed <- ceiling(pf_n_part * tpf_preflight_result$n_needed_factor)
        warning(sprintf(
          "TPF loglik SD at mode = %.2f > 1 (Dynare threshold).",
          tpf_preflight_result$sd),
          "\nPMCMC acceptance dominated by loglik noise.",
          sprintf("\nCurrent n_particles = %d; to achieve SD < 1, raise to ~%d.",
                  pf_n_part, n_needed),
          call. = FALSE)
        tpf_preflight_result$n_needed <- n_needed
      } else if (!is.na(tpf_preflight_result$sd)) {
        .vcat(sprintf("  TPF loglik SD = %.3f (< 1 threshold, OK)\n",
                      tpf_preflight_result$sd))
      }
    }
  }

  # -------------------------------------------------------------------
  # Step 1: Normalise vector arguments
  # -------------------------------------------------------------------
  n_methods <- length(methods)

  recycle <- function(x) {
    if (length(x) == 1L) rep(x, n_methods) else rep_len(x, n_methods)
  }

  methods    <- toupper(as.character(methods))
  nburn_vec  <- as.integer(recycle(nburn))
  ndraws_vec <- as.integer(recycle(ndraws))
  nchains_vec <- as.integer(recycle(nchains))
  npart_vec  <- as.integer(recycle(nparticles))
  ## nwalkers is scalar (DIME only) -- not recycled per method
  nwalkers_val <- nwalkers  # may be NULL (auto-sized in run_dime)

  valid_methods <- c("RWMH", "HMC", "NUTS", "MALA", "SMC", "DIME", "CHEES")
  bad <- setdiff(methods, valid_methods)
  if (length(bad) > 0)
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         ". Valid: ", paste(valid_methods, collapse = ", "))

  .vcat(sprintf("  Methods: %s\n", paste(methods, collapse = " -> ")))
  .vcat(sprintf("  Burn-in: %s\n", paste(nburn_vec, collapse = ", ")))
  .vcat(sprintf("  Draws  : %s\n", paste(ndraws_vec, collapse = ", ")))
  .vcat(sprintf("  Chains : %s\n\n", paste(nchains_vec, collapse = ", ")))

  # -------------------------------------------------------------------
  # Step 2: Run samplers sequentially
  # -------------------------------------------------------------------
  chains_list <- list()
  current_theta <- theta_mode

  for (i in seq_len(n_methods)) {
    m  <- methods[i]
    nb <- nburn_vec[i]
    nd <- ndraws_vec[i]
    nc <- nchains_vec[i]
    np <- npart_vec[i]

    .vcat(sprintf("--- Method %d/%d: %s (%d burn + %d draws) ---\n",
                  i, n_methods, m, nb, nd))

    chain_res <- switch(m,
      "RWMH" = {
        ## CPM dispatch: when cpm_rho_u is set and likelihood = "tpf", use
        ## rwmh_cpm (serial path only). Otherwise fall through to .run_rwmh_batch.
        cpm_rho_u <- (par_ctx$tpf_options %||% list())$cpm_rho_u
        use_cpm   <- !is.null(cpm_rho_u) && is.numeric(cpm_rho_u) &&
                     is.finite(cpm_rho_u) && cpm_rho_u > 0 && cpm_rho_u < 1 &&
                     identical(par_ctx$likelihood, "tpf")
        if (use_cpm) {
          ## CPM preflight message
          if (!is.null(tpf_preflight_result) && !is.na(tpf_preflight_result$sd)) {
            sd_ll     <- tpf_preflight_result$sd
            sd_cpm_eff <- sd_ll * sqrt(2 * (1 - cpm_rho_u))
            message(sprintf(
              "  CPM enabled (rho_u = %.2f). Effective loglik SD ~ %.3f (from %.3f).",
              cpm_rho_u, sd_cpm_eff, sd_ll))
            if (sd_ll < 0.5)
              warning("CPM requested but loglik SD = ", round(sd_ll, 3),
                      " < 0.5. Standard PMCMC already efficient; CPM overhead not needed.",
                      call. = FALSE)
          }
          if (nc > 1L)
            message("  CPM: parallel chains not yet supported; running serial (nc=1).")
          cpm_Sigma <- if (transform_params && !is.null(Sigma_prop_eta)) Sigma_prop_eta else Sigma_prop
          cpm_res <- rwmh_cpm(
            log_post_fn  = log_post_fn,
            theta0       = current_theta,
            Sigma_prop   = cpm_Sigma,
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
        # maps the starting point and falls back to the untransformed path
        # on the parallel mirai branch (transform = NULL there).
        rwmh_Sigma <- if (transform_params && !is.null(Sigma_prop_eta)) Sigma_prop_eta else Sigma_prop
        .run_rwmh_batch(log_post_fn, current_theta, rwmh_Sigma,
                                prior_spec, nchains = nc,
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
                                adapt_cov = rwmh_adapt_cov,
                                n_blocks = rwmh_n_blocks,
                                checkpoint = sampler_checkpoint, ...)
        }  # end else (non-CPM RWMH path)
      },
      "HMC"  = {
        ## metric: "warmup_dense" and "diagonal" are valid for HMC.
        ## "hessian" and "monge" are not passed to dynhr_hmc.
        hmc_metric_arg <- if (metric %in% c("diagonal", "warmup_dense")) metric else "diagonal"
        .run_hmc_batch(log_post_fn, current_theta,
                       n_draws = nd, n_warmup = nb,
                       metric = hmc_metric_arg,
                       # parity with NUTS/MALA: run in eta-space so the
                       # leapfrog does not step-collapse against bounds
                       # (adapt_mass then tunes the diagonal mass in
                       # eta-space). dynhr_hmc maps draws back to theta.
                       transform = if (transform_params) param_transform else NULL,
                       verbose = verbose, ...)
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
          mala_G_inv_try <- Sigma_prop  # Sigma_prop is the inverse-mass (M_inv)
          mala_G_try     <- tryCatch(solve(mala_G_inv_try), error = function(e) NULL)
          if (!is.null(mala_G_try)) {
            mala_G_inv <- mala_G_inv_try
            mala_G     <- mala_G_try
            .vcat("  [MALA] using dense metric from Sigma_prop.\n")
          } else {
            .vcat("  [MALA] metric='hessian': Sigma_prop not invertible -- using identity.\n")
          }
        }

        # Analytic gradient (same gate as NUTS); also needed for the Monge metric_fn
        mala_grad <- if (isTRUE(analytic_grad) &&
                         is.null(mode_result$obc_specs) &&
                         .ctx_allows_analytic_gradient(par_ctx)) {
          grad_method <- .dynhr_opt("grad_method", default = "hybrid")
          .vcat(sprintf("  [MALA] building analytic gradient (%s)...\n", grad_method))
          make_posterior_grad(model, par_data, prior_spec, par_obs_vars,
                              compiled, me_variance = par_me_var,
                              me_extra    = par_me_extra,
                              shock_scale = par_shock_scale,
                              grad_method = grad_method,
                              likelihood  = par_ctx$likelihood,
                              freq_band   = par_ctx$freq_band)
        } else NULL

        # Monge metric: position-dependent G(theta) = I + alpha^2 g g'.
        # Requires a gradient function (analytic if available, else FD).
        # Build the gradient closure for the metric's internal use (may differ
        # from mala_grad if analytic_grad=FALSE; dynhr_mala builds its own FD
        # grad when grad_fn=NULL, but the metric_fn needs its own copy because
        # it is called from .metric_at(), which does NOT pass through grad_fn).
        if (identical(metric, "monge")) {
          ## Monge metric is OFF by default. G = I + alpha^2 g g' inflates the
          ## metric along the gradient, so the proposal step shrinks by
          ## 1/(1 + alpha^2 |g|^2) wherever the gradient is large -- on sharp /
          ## near-unit-root posteriors the chain stalls. Verified independently:
          ## on an anisotropic target the tight direction recovered ~1% of its
          ## variance (ESS/draw ~0.002) vs full recovery under the constant
          ## Laplace metric. Keep it available, but behind an explicit opt-in.
          if (!isTRUE(.dynhr_opt("allow_monge_metric", default = FALSE))) {
            stop("metric = \"monge\" is experimental and disabled by default: the ",
                 "position-dependent Monge metric collapses the proposal step on ",
                 "sharp / near-unit-root posteriors and typically under-explores. ",
                 "Prefer metric = \"hessian\" (constant Laplace) or \"diagonal\". To ",
                 "use it anyway, call dynhr_set_options(allow_monge_metric = TRUE).",
                 call. = FALSE)
          }
          warning("metric = \"monge\": the Monge metric can collapse the proposal ",
                  "step on sharp / near-unit-root posteriors (under-exploration). ",
                  "Check mixing/ESS and prefer metric = \"hessian\" if the chain stalls.",
                  call. = FALSE)
          alpha_monge <- .dynhr_opt("monge_alpha", monge_alpha, default = 1)
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

        .run_mala_batch(log_post_fn, current_theta,
                        n_draws    = nd,   n_warmup = nb,
                        G          = mala_G,
                        G_inv      = mala_G_inv,
                        metric_fn  = mala_metric_fn,
                        grad_fn    = mala_grad,
                        transform  = if (transform_params) param_transform else NULL,
                        chain_id   = 1L, checkpoint = sampler_checkpoint,
                        verbose    = verbose, ...)
      },
      "NUTS" = {
        # The analytic gradient path (tangent/adjoint KF) has no me_extra or
        # shock_scale support: fall back to the numerical gradient.
        # Guard centralised via .ctx_allows_analytic_gradient(par_ctx).
        if (isTRUE(analytic_grad) &&
            !.ctx_allows_analytic_gradient(
              par_ctx %||% estimation_context())) {
          warning("analytic_grad ignored: the analytic gradient path does ",
                  "not support per-period me_extra (filter_tunes) or ",
                  "shock_scale (heteroskedastic_shocks). Using the numerical gradient.",
                  call. = FALSE)
          analytic_grad <- FALSE
        }

        # Parallel multi-chain NUTS on a mirai pool. Standard Gaussian models
        # recompile the posterior per daemon; OBC/PKF and cumulant models ship
        # the already-built log_post_fn closure once instead.
        use_par_nuts <- parallel && nc > 1L &&
                        identical(parallel_backend, "mirai") &&
                        requireNamespace("mirai", quietly = TRUE)
        if (use_par_nuts) {
          # Analytic/implicit gradient on the parallel path requires a
          # standard Gaussian model (par_standard): each daemon needs
          # .worker_model/.worker_cm/.worker_Y from .mirai_pool_init, which
          # only runs when parsed_model/Y are supplied (par_standard). OBC/PKF
          # and cumulant models ship a pre-built log_post_fn closure instead
          # (.mirai_pool_closure) and have no analytic gradient -- mirrors the
          # serial branch's `analytic_grad && par_standard` gate below.
          par_analytic_grad <- isTRUE(analytic_grad) && par_standard
          par_grad_method <- .dynhr_opt("grad_method", default = "hybrid")
          if (isTRUE(analytic_grad) && !par_standard)
            .vcat("  [NUTS] analytic_grad requires a standard Gaussian model (compiled per daemon); ",
                  "this OBC/PKF or cumulant run uses the numerical gradient.\n")
          else if (par_analytic_grad)
            .vcat(sprintf("  [NUTS] parallel chains will build analytic gradients (%s) per daemon...\n",
                          par_grad_method))
          nuts_Sig <- mode_result$V_mode %||% Sigma_prop
          if (is.null(rownames(nuts_Sig)))
            rownames(nuts_Sig) <- colnames(nuts_Sig) <- names(current_theta)
          # transform_params: run_nuts_mirai's `Sigma_prop` (and hence its
          # mass_diag = 1/diag(Sigma_prop)) must be in ETA-SPACE -- convert
          # nuts_Sig via the same delta-method helper used for Sigma_prop_eta
          # above, evaluated at the current chain start `current_theta`
          # (mirrors the serial NUTS branch's `eta_mode_nuts`).
          if (transform_params && !is.null(param_transform)) {
            nuts_Sig_eta <- .cov_theta_to_eta(nuts_Sig, param_transform, current_theta)
            if (!is.null(nuts_Sig_eta)) nuts_Sig <- nuts_Sig_eta
          }
          par_res <- run_nuts_mirai(
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
            grad_method = par_grad_method,
            checkpoint = sampler_checkpoint,
            progress = verbose)
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
        # transform_params: mass_diag must be supplied in ETA-SPACE. If
        # M_inv_diag (= nuts_mass above) approximates Var(theta), then by the
        # delta method Var(eta) ~ Var(theta) / (dtheta/deta)^2, i.e.
        #   m_inv_eta = m_inv_theta / dtheta_deta(eta_mode)^2
        # dtheta_deta is floor-guarded (>= 1e-12 in absolute value) before
        # squaring to avoid blow-up near a transform's asymptote.
        if (transform_params && !is.null(param_transform) && !is.null(nuts_mass)) {
          eta_mode_nuts <- param_transform$to_unconstrained(current_theta)
          d_vec_nuts <- param_transform$dtheta_deta(eta_mode_nuts)
          d_vec_nuts[abs(d_vec_nuts) < 1e-12] <- 1e-12
          nuts_mass <- nuts_mass / d_vec_nuts^2
        }
        # --- Dense metric (metric = "hessian") ---
        # Use Sigma_prop as the dense inverse-mass matrix M_inv.
        # RATIONALE: Sigma_prop is already (a) in the sampler's coordinate
        # space (theta or eta, depending on transform_params) and (b)
        # regularised by the mode-finding step.  The current diagonal path
        # uses 1/diag(Sigma_prop); this is the strict generalisation to the
        # full matrix.  We do NOT use hessian_exact directly to avoid
        # re-introducing the theta-vs-eta coordinate bug.
        nuts_M_inv  <- NULL
        nuts_chol_M <- NULL
        if (identical(metric, "hessian") && !is.null(Sigma_prop) &&
            is.matrix(Sigma_prop) && all(is.finite(Sigma_prop))) {
          # Sigma_prop is the inverse-mass (M_inv = Sigma_prop).
          # M = solve(M_inv); chol_M = chol(M).
          nuts_M_inv_try <- Sigma_prop
          nuts_M_try     <- tryCatch(solve(nuts_M_inv_try), error = function(e) NULL)
          if (!is.null(nuts_M_try)) {
            nuts_chol_M_try <- tryCatch(chol(nuts_M_try), error = function(e) NULL)
            if (!is.null(nuts_chol_M_try)) {
              nuts_M_inv  <- nuts_M_inv_try
              nuts_chol_M <- nuts_chol_M_try
              nuts_mass   <- NULL  # dense path overrides diagonal
              .vcat("  [NUTS] using dense mass matrix from Sigma_prop.\n")
            } else {
              .vcat("  [NUTS] metric='hessian': Sigma_prop not PD -- falling back to diagonal.\n")
            }
          } else {
            .vcat("  [NUTS] metric='hessian': Sigma_prop not invertible -- falling back to diagonal.\n")
          }
        }

        # --- One-shot dense metric from the Whittle FIM (metric = "whittle_fim") ---
        # Opt-in (mirrors the Monge gate exactly): errors unless
        # allow_whittle_fim_metric = TRUE. Assembled ONCE at the mode (current_theta),
        # then frozen through warmup and sampling -- same contract as "hessian".
        # On any failure (non-standard model, degenerate FIM, non-invertible
        # result) falls back to the "hessian" dense metric with a message,
        # rather than erroring the whole run.
        if (identical(metric, "whittle_fim")) {
          if (!isTRUE(.dynhr_opt("allow_whittle_fim_metric", default = FALSE))) {
            stop("metric = \"whittle_fim\" is experimental and disabled by default: ",
                 "it assembles a one-shot dense mass matrix from the Whittle ",
                 "(frequency-domain) Fisher information at the mode, which costs one ",
                 "extra evaluation that can be expensive on large models and is ",
                 "unvalidated on the full estimation pipeline. To use it anyway, call ",
                 "dynhr_set_options(allow_whittle_fim_metric = TRUE).",
                 call. = FALSE)
          }
          if (!par_standard) {
            .vcat("  [NUTS] metric='whittle_fim' requires a standard Gaussian model -- ",
                  "falling back to 'hessian'.\n")
          } else {
            T_obs_wf <- ncol(par_data)
            omega_grid_wf <- 2 * pi * seq_len(floor(T_obs_wf / 2)) / T_obs_wf
            asm <- tryCatch(
              .assemble_dss_list_at_mode(model, compiled, current_theta, par_obs_vars),
              error = function(e) NULL
            )
            if (is.null(asm) || !isTRUE(asm$ok)) {
              .vcat("  [NUTS] metric='whittle_fim': dss_list assembly failed -- falling back to 'hessian'.\n")
            } else {
              fallback_hess <- if (!is.null(nuts_M_inv)) {
                list(G = tryCatch(solve(nuts_M_inv), error = function(e) NULL),
                     G_inv = nuts_M_inv, L = nuts_chol_M, logdet = NA_real_)
              } else NULL
              fim <- tryCatch(
                whittle_fim(TT = asm$TT, RR = asm$RR, ZZ = asm$ZZ, DD = asm$DD,
                           Sigma_e = asm$Sigma_e, dss_list = asm$dss_list,
                           omega_grid = omega_grid_wf, T_obs = T_obs_wf,
                           prior_hess = NULL, fallback_metric = fallback_hess,
                           me_variance = par_me_var %||% 0),
                error = function(e) NULL
              )
              if (is.null(fim) || !is.matrix(fim$G_inv)) {
                .vcat("  [NUTS] metric='whittle_fim': whittle_fim() failed -- falling back to 'hessian'.\n")
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
                  .vcat("  [NUTS] metric='whittle_fim': non-positive/non-finite transform Jacobian or missing FIM factors -- falling back to 'hessian'.\n")
                }
              }
            }
          }
        }

        # Exact analytic gradient (compiled Kalman score for shock-std params +
        # numerical for the rest, or full implicit-differentiation gradient
        # when grad_method = "implicit") when requested on a standard Gaussian
        # model.
        nuts_grad <- if (isTRUE(analytic_grad) &&
                         is.null(mode_result$obc_specs) &&
                         .ctx_allows_analytic_gradient(par_ctx)) {
          grad_method <- .dynhr_opt("grad_method", default = "hybrid")
          .vcat(sprintf("  [NUTS] building analytic gradient (%s)...\n", grad_method))
          make_posterior_grad(model, par_data, prior_spec, par_obs_vars,
                              compiled, me_variance = par_me_var,
                              me_extra    = par_me_extra,
                              shock_scale = par_shock_scale,
                              grad_method = grad_method,
                              likelihood  = par_ctx$likelihood,
                              freq_band   = par_ctx$freq_band)
        } else NULL
        t0_nuts <- proc.time()
        nuts_sec <- max(30L, as.integer(nuts_timeout_seconds))
        res <- tryCatch({
          setTimeLimit(elapsed = nuts_sec, transient = TRUE)
          dynhr_nuts(log_post_fn, current_theta,
                     n_draws = nd, n_warmup = nb,
                     mass_diag = nuts_mass, grad_fn = nuts_grad,
                     M_inv = nuts_M_inv, chol_M = nuts_chol_M,
                     metric = if (metric %in% c("diagonal", "warmup_dense")) metric else "diagonal",
                     transform = if (transform_params) param_transform else NULL,
                     chain_id = 1L, checkpoint = sampler_checkpoint, ...)
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
          run_smc_mirai(parsed_model = if (par_standard) model else NULL,
                        Y = if (par_standard) par_data else NULL,
                        prior_spec = prior_spec, obs_names = par_obs_vars,
                        n_particles = np, n_cores = n_cores,
                        me_variance = par_me_var,
                        me_extra = par_me_extra,
                        shock_scale = par_shock_scale,
                        ctx = par_ctx,
                        log_post_fn = if (par_standard) NULL else log_post_fn,
                        verbose = verbose, ...)
        else
          dynhr_smc(log_post_fn, prior_spec = prior_spec,
                    n_particles = np, ...)
        ## Resample the weighted particle cloud to an equally-weighted draw
        ## matrix so all downstream consumers (diagnostics, Bayesian IRF,
        ## smoother) receive a standard draw matrix.  The original particles
        ## and weights remain in $particles and $smc_weights.
        if (.is_smc_weighted(res)) {
          res$chain <- as_posterior_draws(res)
          res$n_draws <- nrow(res$chain)
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

      "CHEES" = {
        # ChEES-HMC (Hoffman, Radul & Sountsov 2021): fixed-length leapfrog
        # with dual-averaged step size AND trajectory-time adaptation via the
        # ChEES criterion.  Mirrors the NUTS serial branch.
        chees_mass <- if (!is.null(Sigma_prop) && is.matrix(Sigma_prop) &&
                          all(is.finite(diag(Sigma_prop))) &&
                          all(diag(Sigma_prop) > 0)) {
          1 / pmax(diag(Sigma_prop), 1e-12)
        } else NULL
        # transform_params: mass_diag in eta-space (same delta-method as NUTS)
        if (transform_params && !is.null(param_transform) && !is.null(chees_mass)) {
          eta_mode_chees <- param_transform$to_unconstrained(current_theta)
          d_vec_chees <- param_transform$dtheta_deta(eta_mode_chees)
          d_vec_chees[abs(d_vec_chees) < 1e-12] <- 1e-12
          chees_mass <- chees_mass / d_vec_chees^2
        }
        # Analytic gradient (same gate as NUTS)
        chees_grad <- if (isTRUE(analytic_grad) &&
                          is.null(mode_result$obc_specs) &&
                          .ctx_allows_analytic_gradient(par_ctx)) {
          grad_method <- .dynhr_opt("grad_method", default = "hybrid")
          .vcat(sprintf("  [ChEES] building analytic gradient (%s)...\n", grad_method))
          make_posterior_grad(model, par_data, prior_spec, par_obs_vars,
                              compiled, me_variance = par_me_var,
                              me_extra    = par_me_extra,
                              shock_scale = par_shock_scale,
                              grad_method = grad_method,
                              likelihood  = par_ctx$likelihood,
                              freq_band   = par_ctx$freq_band)
        } else NULL
        res <- dynhr_chees(log_post_fn, current_theta,
                           n_draws   = nd, n_warmup = nb,
                           mass_diag = chees_mass, grad_fn = chees_grad,
                           transform = if (transform_params) param_transform else NULL,
                           chain_id  = 1L, checkpoint = sampler_checkpoint,
                           verbose   = verbose, ...)
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
          run_dime_mirai(
            parsed_model = if (par_standard) model else NULL,
            Y            = if (par_standard) par_data else NULL,
            prior_spec   = prior_spec, obs_names = par_obs_vars,
            n_chain      = nwalkers_val,
            n_iter       = nd, n_burn = nb,
            n_cores      = n_cores,
            me_variance  = par_me_var,
            me_extra     = par_me_extra,
            shock_scale  = par_shock_scale,
            ctx          = par_ctx,
            log_post_fn  = if (par_standard) NULL else log_post_fn,
            verbose      = verbose)
        else
          run_dime(log_post_fn, prior_spec = prior_spec,
                   n_chain = nwalkers_val,
                   n_iter  = nd, n_burn = nb,
                   checkpoint = sampler_checkpoint,
                   verbose = verbose, ...)
        list(chains = list(res), chain_stats = data.frame(
          chain = 1L, accept_rate = res$acceptance_rate,
          final_logpost = tail(res$post_logpost[is.finite(res$post_logpost)], 1L),
          stringsAsFactors = FALSE
        ))
      }
    )

    chains_list[[m]] <- chain_res

    # Update starting point for next method
    if (!is.null(chain_res$chains) && length(chain_res$chains) > 0) {
      last_chain <- chain_res$chains[[length(chain_res$chains)]]
      if (!is.null(last_chain$chain) && nrow(last_chain$chain) > 0)
        current_theta <- setNames(as.numeric(last_chain$chain[nrow(last_chain$chain), ]),
                                  colnames(last_chain$chain))
    }

    # Print chain stats
    if (nrow(chain_res$chain_stats) > 0) {
      print(chain_res$chain_stats, row.names = FALSE, digits = 3)
    }
    cat("\n")
  }

  # -------------------------------------------------------------------
  # Step 3: Pool draws and compute convergence
  # -------------------------------------------------------------------
  .vcat("-- Aggregating chains --\n")

  all_chains <- list()
  for (m in methods) {
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
  if (isTRUE(run_stoch_simul)) {
    .vcat("-- Stoch_simul at posterior mean --\n")
    params_post <- model$param_values
    for (nm in names(posterior_mean)) {
      if (nm %in% names(params_post)) params_post[[nm]] <- posterior_mean[[nm]]
    }

    dr <- solved$dr
    irfs <- compute_irfs(dr, model, params = params_post)
    moments <- compute_moments(dr, model, params = params_post)
    if (!is.null(irfs))
      .vcat(sprintf("  IRFs computed: %d shocks\n", length(irfs)))
    if (!is.null(moments))
      .vcat("  Theoretical moments computed.\n")
  }

  # -------------------------------------------------------------------
  # Wrap up
  # -------------------------------------------------------------------
  .vcat("\n================================================================\n")
  .vcat("  Posterior estimation complete\n")
  .vcat("================================================================\n\n")

  result <- list(
    mode_result       = mode_result,
    chains            = chains_list,
    pooled_draws      = pooled_draws,
    convergence       = conv,
    posterior_mean    = posterior_mean,
    posterior_irfs    = irfs,
    posterior_moments = moments,
    tpf_preflight     = tpf_preflight_result,
    meta = list(
      methods       = methods,
      nburn         = nburn_vec,
      ndraws        = ndraws_vec,
      nchains       = nchains_vec,
      n_methods     = n_methods
    )
  )
  class(result) <- c("dynhr_posterior_result", "list")
  invisible(result)
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
  if (!is.null(x$posterior_irfs))
    cat(sprintf("  IRFs      : %d shocks at posterior mean\n", length(x$posterior_irfs)))
  if (!is.null(x$posterior_moments))
    cat("  Moments   : computed at posterior mean\n")
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
#'   \code{run_mcmc_mirai()}'s (PARALLEL/mirai path) `transform` argument;
#'   extra chain starting points are dispersed in eta-space and mapped back
#'   via `to_constrained()` on both paths.
#' @param adapt_cov Opt-in (default \code{FALSE}); forwarded to
#'   \code{rwmh()}'s (serial path) or \code{run_mcmc_mirai()}'s
#'   (PARALLEL/mirai path) `adapt_cov` argument -- Haario et al. (2001)
#'   adaptive proposal covariance.
#' @param n_blocks Opt-in (default \code{1L}); forwarded to \code{rwmh()}'s
#'   (serial path) or \code{run_mcmc_mirai()}'s (PARALLEL/mirai path)
#'   `n_blocks` argument -- randomized parameter blocking.
#' @noRd
.run_rwmh_batch <- function(log_post_fn, theta_mode, Sigma_prop,
                             prior_spec, nchains = 4L,
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
                             seed_base = 42L, ctx = NULL, checkpoint = NULL, ...) {

  n_par     <- length(theta_mode)
  L <- .robust_chol(Sigma_prop, n_par)

  # Parallel path: a single mirai pool runs all chains. Either the parsed
  # model + data (standard Gaussian, recompiled per daemon) or a pre-built
  # log-posterior closure (OBC/PKF, cumulant; shipped once) is required.
  use_par <- isTRUE(parallel) && nchains > 1L &&
             identical(parallel_backend, "mirai") &&
             requireNamespace("mirai", quietly = TRUE) &&
             (isTRUE(use_closure) ||
              (!is.null(parsed_model) && !is.null(Y) && !is.null(obs_names)))
  if (use_par) {
    if (verbose) cat(sprintf("    Parallel RWMH (mirai): %d chains\n", nchains))
    par_res <- run_mcmc_mirai(
      parsed_model = parsed_model, Y = Y,
      prior_spec   = prior_spec, obs_names = obs_names,
      theta_mode   = theta_mode, Sigma_prop = Sigma_prop,
      n_chains     = nchains, n_draws = n_draws, n_burn = n_burn,
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
        sampler         = "rwmh", n_chains = nchains,
        chain_list      = chain_list, chain_stats = par_res$chain_stats), "rwmh")
    } else NULL
    return(list(chains = chain_list, chain_stats = par_res$chain_stats,
                combined = combined, convergence = conv))
  }

  chain_list  <- vector("list", nchains)
  chain_stats <- data.frame(
    chain = integer(), accept_rate = numeric(),
    final_logpost = numeric(), stringsAsFactors = FALSE
  )

  for (ch in seq_len(nchains)) {
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

    if (verbose) cat(sprintf("    Chain %d/%d...\n", ch, nchains))
    chain_list[[ch]] <- rwmh(log_post_fn, th0, Sigma_prop,
                               n_draws = n_draws + n_burn, n_burn = n_burn,
                               transform = transform, chain_id = ch,
                               adapt_cov = adapt_cov, n_blocks = n_blocks,
                               checkpoint = checkpoint, ...)

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
      n_chains        = nchains,
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


#' Run a batch of HMC chains
#' @noRd
.run_hmc_batch <- function(log_post_fn, theta_mode,
                            n_draws = 2000L, n_warmup = 1000L,
                            verbose = TRUE, ...) {
  if (!exists("dynhr_hmc", mode = "function")) {
    if (verbose) cat("    HMC not available, falling back to NUTS...\n")
    return(.run_rwmh_batch(log_post_fn, theta_mode,
                            diag(length(theta_mode)), NULL,
                            nchains = 1L, n_draws = n_draws,
                            n_burn = n_warmup, verbose = verbose, ...))
  }

  if (verbose) cat("    Running HMC...\n")
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
  n_post <- nrow(chain_list[[1]])
  n_par  <- ncol(chain_list[[1]])
  m      <- length(chain_list)
  combined <- do.call(rbind, chain_list)

  chain_means <- sapply(chain_list, colMeans, simplify = "matrix")
  if (is.list(chain_means)) chain_means <- do.call(cbind, chain_means)
  grand_mean  <- rowMeans(chain_means)

  rhat <- ess <- setNames(numeric(n_par), colnames(combined))

  for (j in seq_len(n_par)) {
    W     <- mean(sapply(chain_list, function(ch) var(ch[, j])))
    B     <- n_post / (m - 1) * sum((chain_means[j, ] - grand_mean[j])^2)
    V_hat <- (1 - 1/n_post) * W + (1/n_post) * B
    rhat[j] <- sqrt(V_hat / max(W, 1e-20))

    x       <- combined[, j]
    n_total <- length(x)
    if (var(x) < 1e-20) { ess[j] <- n_total; next }
    max_lag <- min(500, n_total / 3)
    ac      <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
    pairs   <- seq(1, length(ac) - 1, by = 2)
    tau     <- 0
    for (k in pairs) {
      if (ac[k] + ac[k+1] < 0) break
      tau <- tau + ac[k] + ac[k+1]
    }
    ess[j] <- n_total / (1 + 2 * tau)
  }

  list(rhat = rhat, ess = ess, combined = combined)
}
