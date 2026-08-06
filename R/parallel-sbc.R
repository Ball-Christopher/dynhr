## R/parallel-sbc.R
## ---------------------------------------------------------------------------
## Mirai-parallel driver for SBC replications.
##
## `run_sbc_mirai()` mirrors the serial for-loop in dynhr_sbc() but dispatches
## replications over a mirai daemon pool.  Called by dynhr_sbc() when n_cores
## is non-NULL; the serial loop is the NULL fallback.
##
## Key design decisions:
##   - dynhr_compiled (list of R closures) cannot travel through mori
##     shared memory.  We ship the lightweight parsed model and recompile on
##     each daemon inside everywhere(), mirroring .mirai_pool_init()
##     (R/parallel-mirai.R:251-262).  Startup cost: ~0.5-2 s per daemon.
##   - prior_sampler is a closure that closes only over prior_spec (a
##     data.frame); it is safe to ship via .args.
##   - Per-replication seeds: seed_base + i, matching the serial loop in
##     dynhr_sbc(), so serial and parallel give bit-identical results for the
##     same seed.
##   - RNGkind is forced to Mersenne-Twister on each daemon to match the main
##     session (daemons default to L'Ecuyer-CMRG).
## ---------------------------------------------------------------------------


#' Parallel SBC replications via mirai
#'
#' Dispatches \code{n_replications} calls to
#' \code{.sbc_one_replication} across a mirai daemon pool, collects
#' results, and returns a list of per-replication outputs identical in
#' structure to the serial loop in \code{dynhr_sbc}.
#'
#' @param model        Parsed dynhr_mod (lightweight; recompiled on daemons).
#' @param prior_spec   Prior specification data.frame from
#'   \code{extract_prior_spec}.
#' @param prior_sampler Closure (\code{function() -> named numeric}); returned
#'   by \code{.smc_make_prior_sampler}.
#' @param obs_vars     Character vector of observed variable names.
#' @param T_obs        Number of observation periods to simulate.
#' @param presample    Pre-sample burn-in periods (integer).
#' @param n_replications Number of SBC replications.
#' @param n_draws,n_burn,thin Sampler settings.
#' @param sampler      Currently only \code{"rwmh"}.
#' @param me_variance  Measurement error variance (scalar or named vector).
#' @param lik_init     Kalman filter initial conditions method.
#' @param transform_params Logical; apply unconstrained reparameterisation.
#' @param adapt_cov    Logical; adapt RWMH proposal covariance.
#' @param seed_base    Base seed; replication \code{i} uses
#'   \code{seed_base + i}.
#' @param n_cores      Number of daemon workers.
#' @param verbose      Logical; print per-replication progress.
#' @return A list of length \code{n_replications}; each element is either
#'   \code{list(ok=TRUE, ranks=<named integer>, L_effective=<int>)} or
#'   \code{list(ok=FALSE, reason=<character>)}.
#' @noRd
run_sbc_mirai <- function(model, prior_spec, prior_sampler,
                          obs_vars, T_obs, presample,
                          n_replications, n_draws, n_burn, thin,
                          sampler, me_variance, lik_init,
                          transform_params, adapt_cov,
                          seed_base, n_cores, verbose = TRUE,
                          likelihood = "gaussian", order = 1L) {

  ## ---- daemon pool -------------------------------------------------------
  ## BLAS/OpenMP pinned to 1 per daemon (env inherited at spawn, restored on
  ## exit) -- prevents OpenBLAS-build oversubscription; see
  ## .mirai_pin_blas_threads and the note in .smc_pool_setup.
  .restore_blas <- .mirai_pin_blas_threads()
  on.exit(.restore_blas(), add = TRUE)
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(NULL), add = TRUE)

  ## Load dynhr on every daemon and recompile the model once per daemon.
  ## `<<-` is required so the binding reaches the daemon's global environment
  ## (plain `<-` stays in the everywhere() expression's local frame and is not
  ## visible to subsequent mirai_map tasks via bare name lookup).
  mirai::everywhere(
    {
      suppressMessages(library(dynhr))
      .cmpl        <- utils::getFromNamespace("compile_model", "dynhr")
      .worker_cm  <<- .cmpl(.worker_model, verbose = FALSE)
      .worker_model <<- .worker_model
    },
    .args = list(.worker_model = model)
  )

  ## ---- dispatch ----------------------------------------------------------
  ## Build a task wrapper; all per-replication args are shipped via .args.
  ## seed_base + .i matches the serial loop: seed + i.
  task_args <- list(
    .prior_spec       = prior_spec,
    .obs_vars         = obs_vars,
    .T_obs            = T_obs,
    .presample        = presample,
    .n_draws          = n_draws,
    .n_burn           = n_burn,
    .thin             = thin,
    .sampler          = sampler,
    .me_variance      = me_variance,
    .lik_init         = lik_init,
    .transform_params = transform_params,
    .adapt_cov        = adapt_cov,
    .seed_base        = seed_base,
    .verbose          = verbose,
    .likelihood       = likelihood,
    .order            = order
  )

  sbc_task <- function(.i,
             .prior_spec, .obs_vars, .T_obs, .presample,
             .n_draws, .n_burn, .thin, .sampler, .me_variance, .lik_init,
             .transform_params, .adapt_cov, .seed_base, .verbose,
             .likelihood, .order) {
      ## Force Mersenne-Twister to match the main session's RNG.
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      set.seed(.seed_base + .i)

      ## Retrieve the per-daemon globals set by everywhere().
      .m  <- get(".worker_model", envir = globalenv())
      .cm <- get(".worker_cm",   envir = globalenv())

      ## Rebuild prior_sampler on the daemon from prior_spec so the closure
      ## environment is entirely local (no cross-session frame references).
      .ps <- dynhr:::.smc_make_prior_sampler(.prior_spec)

      dynhr:::.sbc_one_replication(
        i                = .i,
        model            = .m,
        compiled         = .cm,
        prior_spec       = .prior_spec,
        prior_sampler    = .ps,
        obs_vars         = .obs_vars,
        T_obs            = .T_obs,
        presample        = .presample,
        n_draws          = .n_draws,
        n_burn           = .n_burn,
        thin             = .thin,
        sampler          = .sampler,
        me_variance      = .me_variance,
        lik_init         = .lik_init,
        transform_params = .transform_params,
        adapt_cov        = .adapt_cov,
        likelihood       = .likelihood,
        order            = .order,
        seed             = .seed_base + .i,
        verbose          = .verbose
      )
  }
  ## Sever the task env (it takes all data via .args and uses dynhr::: / daemon
  ## globals) so mirai_map does not serialise this frame's heavy locals (the
  ## parsed model etc.) with every one of the n_replications tasks.
  environment(sbc_task) <- asNamespace("dynhr")
  results <- mirai::mirai_map(
    seq_len(n_replications), sbc_task, .args = task_args
  )[]

  ## ---- collect -----------------------------------------------------------
  ## results is a list of mirai value objects (already resolved by []).
  lapply(results, function(r) {
    if (inherits(r, "error")) {
      list(ok = FALSE, reason = conditionMessage(r))
    } else {
      r
    }
  })
}
