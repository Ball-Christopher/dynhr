## R/smc2.R
## --------------------------------------------------------------------------
## dynhr_smc2() -- SMC^2: Sequential Monte Carlo over parameters (theta),
## wrapping a NOISY but UNBIASED particle-filter likelihood estimate (the
## order-3 Tempered Particle Filter, make_log_posterior_tpf(), or the
## measurement-side SV Rao-Blackwellized particle filter,
## make_log_posterior_sv_rbpf()).
##
## This is pseudo-marginal SMC on the EXTENDED space (theta, u) where u is
## the particle filter's internal randomness (Chopin, Jacob & Papaspiliopoulos
## 2013 "SMC^2"; Andrieu, Doucet & Holenstein 2010 for the pseudo-marginal
## argument PMMH borrows from). Because the PF loglik estimate is unbiased,
## the final-stage (lambda = 1) marginal of the extended-space target IS the
## exact posterior, and the telescoped SMC normalising constant is a valid
## (noisy but unbiased) estimate of the evidence -- even though intermediate
## tempered targets pi_lambda ∝ prior * Lhat^lambda are NOT equal to
## prior * L^lambda for any single u (they only equal it in expectation).
##
## Key observation exploited here: dynhr_smc()'s existing tempering +
## mutation loop (R/sampler-smc.R) is ALREADY pseudo-marginal-correct by
## construction, for any log_post_fn that returns a fresh, unbiased loglik
## estimate on every call:
##   - stage-0 particles get ONE fresh evaluation each (their stored loglik);
##   - incremental tempering weights use ONLY the stored log_liks vector
##     (never re-evaluated between stages, even across a resample);
##   - the RWMH mutation step computes the incumbent's tempered log-target
##     from its STORED loglik (ll_i, carried in from the previous stage /
##     stage 0), evaluates the PROPOSAL fresh via log_post_fn(theta_prop) --
##     the ONE PF evaluation per proposal -- and on accept promotes the
##     proposal's fresh loglik to be the new stored value. The incumbent is
##     NEVER re-run ("Monte Carlo within Metropolis", which would be invalid,
##     never happens).
## So dynhr_smc2() does not reimplement tempering/resampling/mutation; it
## composes the noisy-likelihood closure (this file's only new logic) and
## forwards to dynhr_smc() unmodified, exactly per the brief's "reuse the
## existing machinery in R/sampler-smc.R" instruction.
## --------------------------------------------------------------------------

## Allow-lists of factory-specific pass-through arguments accepted via
## `likelihood_args`. Deliberately narrow and deliberately SEPARATE from
## dynhr_smc2's own top-level arguments: both make_log_posterior_tpf() and
## make_log_posterior_sv_rbpf() have their OWN `n_particles` argument (the
## number of PARTICLE-FILTER particles), which is a different quantity from
## dynhr_smc2's `n_particles` (the number of THETA particles in the outer
## SMC population) -- reusing a single flat `n_particles=` argument for both
## would silently conflate the two counts. Routing PF-level tuning through a
## named list mirrors the codebase's existing ctx$tpf_options convention
## (see the "tpf" branch of make_log_posterior() in R/posterior.R).
## `mh_scale` was REMOVED from this list (and from the whole TPF stack) on
## 2026-09-02: it had been inert since the mutation step was corrected to hold
## the ancestor state fixed, so it tuned nothing while looking like a knob.
## Passing it now trips the allow-list check below, which is the point --
## a silently ignored tuning parameter is worse than an error.
.smc2_tpf_allow <- c("n_particles", "ess_target", "n_mh",
                      "max_stages_u", "order", "burn_in_init")
.smc2_sv_allow  <- c("n_particles", "stochastic_volatility", "power")

#' Reject `likelihood_args` keys the chosen factory does not accept.
#'
#' The allow-list used to be applied with \code{intersect()} alone, which
#' SILENTLY DROPPED anything unrecognised -- so a typo, or a knob that had
#' been retired (\code{mh_scale}), looked like it was tuning the filter while
#' doing nothing at all. That is the exact failure mode this list exists to
#' prevent, so an unknown key is now an error naming the accepted set.
#' \code{seed} is checked separately and earlier (it has its own message).
#' @noRd
.smc2_check_allowed <- function(args, allow, which) {
  if (!length(args)) return(invisible(NULL))
  nm  <- names(args)
  bad <- setdiff(nm[nzchar(nm)], c(allow, "seed"))
  if (length(bad))
    stop("dynhr_smc2: likelihood_args ",
         paste(sQuote(bad), collapse = ", "),
         if (length(bad) > 1L) " are not accepted" else " is not accepted",
         " for likelihood = ", sQuote(which), ". Accepted: ",
         paste(sQuote(allow), collapse = ", "), ".", call. = FALSE)
  if (any(!nzchar(nm)) || is.null(nm))
    stop("dynhr_smc2: every element of likelihood_args must be named.",
         call. = FALSE)
  invisible(NULL)
}


#' Compose the noisy, unbiased particle-filter log-posterior used by
#' dynhr_smc2()
#'
#' Thin dispatcher over make_log_posterior_tpf() / make_log_posterior_sv_rbpf()
#' (both @export'd factories owned elsewhere -- this file never redefines
#' their internals). Always forces `seed = NULL` in the factory call: SMC^2's
#' pseudo-marginal validity requires a FRESH, independent particle-filter draw
#' at every evaluation (mirrors pmmh()'s and both factories' own seed = NULL
#' convention for valid PMMH/SMC^2; see their roxygen "seed" sections). A
#' caller-supplied `likelihood_args$seed` is rejected up front rather than
#' silently overridden, so a misunderstanding is loud, not silent.
#'
#' @return function(theta) -> list(logpost, loglik, logprior), suitable as
#'   dynhr_smc()'s `log_post_fn`.
#' @noRd
.smc2_make_loglik <- function(likelihood, model, data, prior_spec, obs_vars,
                               compiled, me_variance, system_priors,
                               likelihood_args) {
  if (!is.null(likelihood_args$seed)) {
    stop("dynhr_smc2: likelihood_args$seed must not be set. SMC^2's pseudo-",
         "marginal validity requires a FRESH, independent particle-filter ",
         "draw at every evaluation -- a fixed seed makes the loglik ",
         "deterministic and breaks the extended-space SMC argument (the same ",
         "reason pmmh() and both factories require seed = NULL for MCMC/SMC ",
         "use).", call. = FALSE)
  }

  if (identical(likelihood, "tpf")) {
    if (!is.numeric(me_variance) || length(me_variance) != 1L ||
        !is.finite(me_variance) || me_variance <= 0) {
      stop("dynhr_smc2: likelihood = \"tpf\" requires a positive scalar ",
           "`me_variance` (the TPF's tempering instrument); got ",
           if (is.null(me_variance)) "NULL" else me_variance, ".",
           call. = FALSE)
    }
    .smc2_check_allowed(likelihood_args, .smc2_tpf_allow, "tpf")
    extra <- likelihood_args[intersect(names(likelihood_args), .smc2_tpf_allow)]
    do.call(make_log_posterior_tpf,
            c(list(model = model, data = data, prior_spec = prior_spec,
                   obs_vars = obs_vars, compiled = compiled,
                   me_variance = me_variance, system_priors = system_priors,
                   seed = NULL),
              extra))
  } else if (identical(likelihood, "sv_rbpf")) {
    .smc2_check_allowed(likelihood_args, .smc2_sv_allow, "sv_rbpf")
    extra <- likelihood_args[intersect(names(likelihood_args), .smc2_sv_allow)]
    do.call(make_log_posterior_sv_rbpf,
            c(list(model = model, data = data, prior_spec = prior_spec,
                   obs_vars = obs_vars, compiled = compiled,
                   me_variance = if (is.null(me_variance)) 0 else me_variance,
                   seed = NULL),
              extra))
  } else {
    stop("dynhr_smc2: unknown `likelihood` ", sQuote(as.character(likelihood)[1]),
         " -- must be \"tpf\" (order-3 Tempered Particle Filter, via ",
         "make_log_posterior_tpf()) or \"sv_rbpf\" (stochastic-volatility ",
         "Rao-Blackwellized particle filter, via make_log_posterior_sv_rbpf()).",
         call. = FALSE)
  }
}


#' SMC^2 -- Sequential Monte Carlo over parameters wrapping a noisy unbiased
#' particle-filter likelihood
#'
#' Likelihood-tempered SMC over the parameter vector theta, where the
#' likelihood at each theta is a NOISY but UNBIASED estimate from a particle
#' filter (order-3 Tempered Particle Filter or the SV-on-shocks
#' Rao-Blackwellized particle filter), with pseudo-marginal (PMMH-style)
#' Metropolis-Hastings mutation moves. This is the parameter-SMC analogue of
#' \code{\link{pmmh}()}: exactly as PMMH is RWMH driven by an unbiased
#' particle loglik, \code{dynhr_smc2()} is \code{dynhr_smc()} driven by
#' one.
#'
#' \strong{Correctness (why this is valid).} Each theta-particle carries its
#' own stored log-likelihood estimate, drawn fresh when the particle was
#' created (stage 0) or last accepted a mutation proposal. The incumbent's
#' stored estimate is NEVER re-evaluated -- re-running the particle filter on
#' the current particle during a mutation sweep ("Monte Carlo within
#' Metropolis") would bias the sampler. Tempering increments and the
#' evidence-telescoping product use the stored estimates; RWMH mutation moves
#' evaluate a fresh particle-filter estimate ONLY at the proposal, and promote
#' it to be the new stored value on acceptance. This IS the extended-space
#' (theta, u) SMC construction (Chopin, Jacob & Papaspiliopoulos 2013; the
#' pseudo-marginal argument of Andrieu, Doucet & Holenstein 2010): because
#' the PF estimate is unbiased, the lambda = 1 marginal is the exact
#' posterior and the accumulated log-evidence is unbiased, even though
#' intermediate tempered targets are not literally prior * L^lambda for any
#' fixed realisation of the filter's randomness.
#'
#' All tempering, adaptive-ESS scheduling, resampling and mutation machinery
#' is reused UNCHANGED from \code{dynhr_smc()} (R/sampler-smc.R); this
#' function's only job is composing the noisy-likelihood closure via
#' \code{make_log_posterior_tpf()} / \code{make_log_posterior_sv_rbpf()} and
#' forwarding to it. The returned list has the same shape as
#' \code{dynhr_smc()}'s (same field names), so the internal helpers
#' \code{dynhr:::smc_summary()}/\code{dynhr:::smc_plot_diagnostics()}
#' (unexported) work on it directly.
#'
#' @param model dynhr_mod from \code{\link{parse_mod}}.
#' @param data Observation matrix in the factories' native n_obs x T
#'   orientation (columns = time periods; rows = \code{obs_vars}, in order).
#' @param prior_spec Prior specification (data.frame or named list); also
#'   used to build the theta-particle prior sampler via the same
#'   \code{.smc_make_prior_sampler()} dynhr_smc() itself uses.
#' @param obs_vars Character vector of observed variable names.
#' @param compiled dynhr_compiled from \code{\link{compile_model}}.
#' @param likelihood \code{"tpf"} (order-3 Tempered Particle Filter) or
#'   \code{"sv_rbpf"} (SV-on-shocks Rao-Blackwellized particle filter).
#' @param me_variance Measurement-error variance. Required (positive scalar)
#'   for \code{likelihood = "tpf"} (the TPF's tempering instrument); optional
#'   for \code{"sv_rbpf"} (default 0, matching
#'   \code{make_log_posterior_sv_rbpf}'s own default).
#' @param system_priors Optional system priors list forwarded to
#'   \code{make_log_posterior_tpf} (ignored for \code{"sv_rbpf"}, which has no
#'   such argument).
#' @param likelihood_args Named list of extra arguments forwarded to the
#'   chosen factory. For \code{"tpf"}: any of \code{n_particles} (PF particle
#'   count), \code{ess_target} (the FILTER's internal resampling threshold --
#'   distinct from this function's own \code{ess_target}, which governs the
#'   OUTER theta-tempering schedule), \code{n_mh},
#'   \code{max_stages_u}, \code{order}, \code{burn_in_init} (e.g. the
#'   SBC-certified order-3 config \code{list(n_particles = 150, n_mh = 1,
#'   burn_in_init = 50)}). For \code{"sv_rbpf"}: \code{n_particles},
#'   \code{stochastic_volatility}, \code{power}. \code{seed} must NOT be set
#'   here (errors loudly) -- see Correctness above.
#' @param n_particles Number of THETA particles in the outer SMC population
#'   (not to be confused with the PF particle count inside
#'   \code{likelihood_args}).
#' @param ess_target Target ESS ratio for the OUTER (theta) adaptive
#'   tempering schedule (0.5-0.99).
#' @param n_mh_steps RWMH mutation steps per outer tempering stage.
#' @param mh_scale_factor,mut_target,mixture_weights Forwarded to
#'   \code{dynhr_smc} unchanged (mutation proposal scale / target
#'   acceptance / optional Herbst-Schorfheide 3-component mixture proposal).
#' @param parallel Evaluate the per-theta-particle-filter loglik calls (both
#'   stage-0 initials and RWMH mutation proposals) across a persistent mirai
#'   daemon pool, via \code{dynhr_smc}'s own \code{parallel = TRUE}
#'   path (R/sampler-smc.R: \code{.smc_pool_setup}/\code{.smc_pmap}) --
#'   forwarded unchanged. This is where SMC^2's cost lives (one particle-filter
#'   run per theta-particle per stage), so parallelising it is the payoff.
#'   Pseudo-marginal validity is preserved: each worker task evaluates ONLY a
#'   stage-0 initial draw or a fresh RWMH proposal, never a stored incumbent
#'   (see \code{dynhr_smc}'s \code{.eval_particle}/\code{.mutate_one}), and
#'   \code{.smc_pmap} seeds every task deterministically from
#'   \code{seed_base}/the mutation-stage counter/particle index -- a fixed-but-
#'   distinct seed per (stage, particle) task, which still gives an
#'   independent particle-filter draw \code{u} for every proposal (the
#'   pseudo-marginal requirement), not a seed reused across evaluations of the
#'   same theta. If \pkg{mirai} (or the requested \code{backend}) is
#'   unavailable, this mirrors \code{dynhr_smc}'s own behaviour exactly: no
#'   error is raised, and the run silently falls back to
#'   \pkg{future.apply} (if installed) or plain serial evaluation.
#' @param backend Parallel backend, forwarded to \code{dynhr_smc}
#'   (\code{"mirai"}, the default and only backend combined with SMC^2's
#'   pseudo-marginal seeding guarantees above; \code{"future"} is accepted for
#'   parity with \code{dynhr_smc} but not specifically verified here).
#' @param seed_base Integer seed base forwarded to \code{dynhr_smc}. NOTE:
#'   it governs only the per-task seeding of the parallel (\code{mirai})
#'   evaluation path; in the default serial mode, and for the host-side
#'   resampling uniforms in either mode, it has no effect -- SMC^2 runs are
#'   NOT reproducible end-to-end from \code{seed_base} (and the
#'   particle-filter randomness inside each loglik evaluation must stay
#'   fresh regardless; see Correctness above).
#' @param verbose Print progress.
#' @param progressor Optional progressr callback, forwarded to
#'   \code{dynhr_smc}.
#'
#' @return A \code{\link{dynhr_chains}} object (so \code{print()},
#'   \code{summary()} and \code{plot()} work exactly as for
#'   \code{\link{mcmc}} / \code{\link{smc}} / \code{\link{nuts}} /
#'   \code{\link{dime}}).  Every field of \code{dynhr_smc()}'s return value is
#'   retained unchanged (particles/chain, smc_weights, log_liks,
#'   log_marginal_lik, lambda_schedule, ess_schedule, accept_schedule,
#'   n_eval, ...), plus \code{$sampler = "smc2"} and \code{$likelihood}
#'   recording which particle filter was used.
#'
#' @references
#' Chopin, N., Jacob, P.E. & Papaspiliopoulos, O. (2013). SMC^2: an efficient
#' algorithm for sequential analysis of state space models. \emph{JRSS-B}
#' 75(3): 397-426.
#'
#' Andrieu, C., Doucet, A. & Holenstein, R. (2010). Particle Markov chain
#' Monte Carlo methods. \emph{JRSS-B} 72(3): 269-342.
#'
#' Herbst, E. & Schorfheide, F. (2014). Sequential Monte Carlo sampling for
#' DSGE models. \emph{J. Applied Econometrics} 29(7): 1073-1098.
#'
#' @seealso \code{dynhr_smc}, \code{\link{pmmh}},
#'   \code{\link{make_log_posterior_tpf}}, \code{\link{make_log_posterior_sv_rbpf}}
#'
#' @examples
#' \donttest{
#' ## Deliberately tiny: 20 theta-particles around a 100-particle inner TPF.
#' ## A real run needs hundreds of each, and the TPF needs me_variance > 0.
#' model    <- parse_mod(system.file("extdata/models/nk_demo.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, max_order = 2L, verbose = FALSE)
#' priors   <- prior_spec(model)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#'
#' fit <- dynhr_smc2(model, Y, priors, obs_vars, compiled,
#'                   likelihood      = "tpf",
#'                   me_variance     = 1e-4,
#'                   n_particles     = 20L,
#'                   likelihood_args = list(n_particles = 100L),
#'                   verbose = FALSE)
#' fit
#' }
#' @export
dynhr_smc2 <- function(
    model, data, prior_spec, obs_vars, compiled,
    likelihood         = c("tpf", "sv_rbpf"),
    me_variance        = NULL,
    system_priors      = NULL,
    likelihood_args    = list(),
    n_particles        = 200L,
    ess_target         = 0.5,
    n_mh_steps         = 1L,
    mh_scale_factor    = 0.5,
    mut_target         = 0.25,
    mixture_weights    = NULL,
    parallel           = FALSE,
    backend            = "mirai",
    seed_base          = 1L,
    verbose            = TRUE,
    progressor         = NULL
) {
  if (!is.character(likelihood) || length(likelihood) < 1L ||
      !likelihood[1] %in% c("tpf", "sv_rbpf")) {
    stop("dynhr_smc2: `likelihood` must be \"tpf\" or \"sv_rbpf\", got ",
         if (is.character(likelihood) && length(likelihood) >= 1L)
           sQuote(likelihood[1]) else class(likelihood)[1], ".", call. = FALSE)
  }
  likelihood <- likelihood[1]

  loglik_fn <- .smc2_make_loglik(likelihood, model, data, prior_spec,
                                  obs_vars, compiled, me_variance,
                                  system_priors, likelihood_args)

  prior_sampler <- .smc_make_prior_sampler(prior_spec)

  out <- dynhr_smc(
    log_post_fn      = loglik_fn,
    prior_sampler    = prior_sampler,
    n_particles      = n_particles,
    ess_target       = ess_target,
    n_mh_steps       = n_mh_steps,
    mh_scale_factor  = mh_scale_factor,
    mut_target       = mut_target,
    mixture_weights  = mixture_weights,
    ## Forwarded unchanged to dynhr_smc()'s own mirai daemon-pool path
    ## (.smc_pool_setup()/.smc_pmap(), R/sampler-smc.R). Verified pseudo-
    ## marginal-safe: .eval_particle()/.mutate_one() (the two task closures
    ## shipped to workers) evaluate ONLY stage-0 initial draws / fresh RWMH
    ## proposals -- the stored incumbent loglik (ll_i) is closed over from the
    ## driver process and never re-run on a worker -- and .smc_pmap() seeds
    ## every task deterministically from seed_base + stage*n_particles + i,
    ## giving each (stage, particle) task its own independent PF draw. See
    ## dynhr_smc2()'s `parallel` roxygen for the full argument.
    parallel         = parallel,
    backend          = backend,
    seed_base        = seed_base,
    verbose          = verbose,
    progressor       = progressor
  )

  out$sampler    <- "smc2"
  out$likelihood <- likelihood
  ## C6 (API consistency): every sampler entry point returns a dynhr_chains
  ## object. new_dynhr_chains() only stamps the class and normalises
  ## $n_draws to nrow($chain) -- no field of dynhr_smc()'s return value is
  ## dropped or renamed, so the fields test-smc2.R pins are all still there.
  new_dynhr_chains(out, "smc2")
}
