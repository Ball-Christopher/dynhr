## R/estimation-context.R
## --------------------------------------------------------------------------
## estimation_context() -- bundle all per-estimation options into one object.
##
## Motivation: options that currently thread through ~6 entry points x
## serial+parallel paths (me_variance, likelihood, lik_init, me_extra,
## shock_scale, freq_band, system_priors, tpf_options, gradient_policy).
## A new option costs one new field here, not a package-wide sweep.
##
## ROADMAP Tier 7 item 2.
## ROADMAP Tier 8 item 10 (Phase 1): plan compiles into ctx.
## --------------------------------------------------------------------------


#' Create an estimation context object
#'
#' Bundles all per-estimation options that thread through the estimation
#' pipeline (mode-finding, posterior estimation, parallel runners) into a
#' single \code{estimation_context} object.  Pass it to
#' \code{make_log_posterior}, \code{\link{run_mode_finding}},
#' \code{\link{run_posterior_estimation}}, \code{\link{run_full_estimation}},
#' or \code{\link{dynhr_sbc}} to override the default individual arguments.
#'
#' When \code{plan} is supplied the plan is compiled into the context:
#' \code{plan_to_filter_tunes} resolves the filter-tunes spec and
#' \code{plan_to_shock_scale_spec} resolves the shock-scale spec.  Both specs
#' are stored on the plan object (as attributes) and the plan itself is stored
#' in \code{$plan} for provenance.  The model-aware entry points
#' (\code{\link{run_full_estimation}}, \code{\link{run_mode_finding}}) consume
#' these specs via the existing helpers (\code{.resolve_filter_tunes},
#' \code{.expand_observables_for_tunes}, \code{.build_shock_scale_matrix}).
#'
#' @param me_variance  Measurement-error variance (scalar, default 0).
#'   Required to be \code{> 0} when \code{likelihood = "tpf"}.
#' @param likelihood   Likelihood type: \code{"gaussian"} (Kalman filter,
#'   default), \code{"cumulant"}, \code{"whittle"}, \code{"tpf"},
#'   \code{"pskf"} (Pruned Skewed Kalman Filter), \code{"student_t"}
#'   (Gaussian KF recursions with multivariate Student-t per-period
#'   log-likelihood; requires \code{student_df}), or \code{"pkf"}
#'   (OBC Penalty Kalman Filter used by \code{run_mode_finding} /
#'   \code{run_full_estimation} when the model has occasionally-binding
#'   constraints).
#' @param lik_init     Kalman filter \code{P0} initialisation (default
#'   \code{"auto"}).  Forwarded to \code{\link{kalman_filter}} on the
#'   Gaussian path.  Ignored for cumulant/whittle/tpf likelihoods.
#' @param me_extra     \code{n_obs x T} matrix of per-period extra measurement
#'   error variances (resolved from a \code{filter_tunes} block); \code{NULL}
#'   = no filter tunes.  Must be \code{NULL} when \code{plan} is supplied.
#' @param shock_scale  \code{n_exo x T} matrix of per-period shock standard
#'   deviation scale factors (resolved from a \code{heteroskedastic_shocks}
#'   block); \code{NULL} = no heteroskedastic shocks.  Must be \code{NULL}
#'   when \code{plan} is supplied.
#' @param freq_band    Numeric(2) \code{c(lo, hi)} in radians.  Only
#'   meaningful when \code{likelihood = "whittle"}; a warning is issued if
#'   a non-default band is supplied with another likelihood.  Default
#'   \code{c(0, pi)}.
#' @param system_priors  Named list of system-prior functions, or \code{NULL}.
#'   Closures in this list must only close over scalar constants (IRF
#'   horizons, target values, distribution parameters), not large
#'   environments -- see \code{\link{validate_context}}.
#' @param tpf_options  Named list of options forwarded to
#'   \code{\link{make_log_posterior_tpf}} when \code{likelihood = "tpf"}.
#'   Typical keys: \code{n_particles}, \code{ess_target}, \code{n_mh},
#'   \code{seed}, \code{cpm_rho_u} (CPM AR(1) correlation, default \code{NULL}
#'   = disabled; set to a value in (0, 1) to enable the correlated
#'   pseudo-marginal RWMH sampler -- see \code{rwmh_cpm}).
#' @param obc_options  Named list of options for the OBC PKF path
#'   (\code{likelihood = "pkf"}).  Reserved keys: \code{max_inner},
#'   \code{proposal} (for PPF; \code{"bootstrap"} or \code{"copf"}),
#'   \code{N_particles} (for PPF).  Default \code{list()}.
#' @param obc_specs   Pre-parsed OBC specs list from \code{obc_parse_tags()}.
#'   When \code{NULL} (default), \code{conditional_forecast()} will call
#'   \code{obc_parse_tags()} internally.  Set by \code{run_mode_finding()} /
#'   \code{run_full_estimation()} to avoid re-parsing on each forecast call.
#' @param gradient_policy  \code{"auto"} (default), \code{"analytic"}, or
#'   \code{"numerical"}.  Controls the gradient path used by
#'   \code{make_posterior_grad}; \code{me_extra} and \code{shock_scale} are
#'   now supported by the analytic gradient path (Tier 7, ROADMAP item 3).
#' @param plan         Optional \code{dynhr_plan} object.  When supplied the
#'   plan is compiled into the context: the filter-tunes and shock-scale specs
#'   are resolved and attached to the stored plan.  The plan is also stored in
#'   \code{$plan} for provenance.  Supplying explicit \code{me_extra} or
#'   \code{shock_scale} alongside \code{plan} is an error.
#' @param sample_start Optional character date string (e.g. \code{"1990q1"})
#'   required when the plan contains date-literal tune periods.  Passed to
#'   \code{plan_to_filter_tunes}.
#' @param ms_spec  Optional \code{\link{ms_dsge_spec}} object for
#'   Markov-switching DSGE estimation (shock-variance switching only).  When
#'   non-\code{NULL} the likelihood is evaluated via \code{\link{ms_kim_filter}}
#'   (Kim-Nelson GPB(2) filter) instead of the standard Kalman filter.
#'   Incompatible with \code{likelihood != "gaussian"}, \code{me_extra},
#'   \code{shock_scale}, and analytic gradient.  Default \code{NULL}
#'   (non-MS, standard path).  Supply exactly one of \code{ms_spec} or
#'   \code{ms_struct_spec}, or neither.
#' @param ms_struct_spec  Optional \code{\link{ms_struct_spec}} object for
#'   structural Markov-switching DSGE estimation (full structural-parameter
#'   switching across regimes).  When non-\code{NULL} the posterior evaluator
#'   calls \code{\link{solve_ms_perturbation}} at each draw to compute
#'   regime-specific decision rules, then evaluates the likelihood via
#'   \code{\link{ms_kim_filter_struct}}.  Incompatible with \code{ms_spec},
#'   \code{me_extra}, \code{shock_scale}, and analytic gradient.  Default
#'   \code{NULL} (non-structural-MS path).
#' @param student_df  Positive finite scalar; degrees of freedom for the
#'   multivariate Student-t per-period log-likelihood.  Required when
#'   \code{likelihood = "student_t"}, ignored otherwise.
#' @return An object of class \code{"dynhr_estimation_context"}.
#' @seealso \code{\link{validate_context}}, \code{ctx_from_mode_result}
#' @export
estimation_context <- function(
    me_variance     = 0,
    likelihood      = c("gaussian", "cumulant", "whittle", "tpf", "pskf",
                        "student_t", "pkf", "ppf", "copf", "pruned"),
    lik_init        = "auto",
    me_extra        = NULL,
    shock_scale     = NULL,
    freq_band       = c(0, pi),
    system_priors   = NULL,
    tpf_options     = list(),
    obc_options     = list(),
    obc_specs       = NULL,
    gradient_policy = c("auto", "analytic", "numerical"),
    plan            = NULL,
    sample_start    = NULL,
    ms_spec         = NULL,
    ms_struct_spec  = NULL,
    student_df      = NULL
) {
  likelihood      <- match.arg(likelihood)
  gradient_policy <- match.arg(gradient_policy)

  ## ---- Plan compilation (Tier 8 item 10 Phase 1) ----------------------
  ## When a plan is supplied: run the two estimation adapters
  ## (plan_to_filter_tunes + plan_to_shock_scale_spec) and attach the
  ## resolved specs as attributes on the plan object.  The model-aware entry
  ## points consume these specs via .resolve_filter_tunes /
  ## .expand_observables_for_tunes / .build_shock_scale_matrix, exactly as
  ## the legacy inline guard blocks did.  We do NOT resolve to matrices here
  ## because that requires the parsed model (varexo_names, T, etc.).
  ## The plan is stored in $plan for provenance.
  stored_plan <- NULL
  if (!is.null(plan)) {
    if (!inherits(plan, "dynhr_plan"))
      stop("estimation_context: 'plan' must be a dynhr_plan object.", call. = FALSE)
    if (!is.null(me_extra))
      stop("estimation_context: supply either 'plan' or 'me_extra', not both.",
           call. = FALSE)
    if (!is.null(shock_scale))
      stop("estimation_context: supply either 'plan' or 'shock_scale', not both.",
           call. = FALSE)
    .check_plan_tpf_compat(plan, likelihood = likelihood)
    ## Resolve specs and attach as attributes so entry points can pick them up.
    attr(plan, ".filter_tunes_spec") <- plan_to_filter_tunes(plan,
                                          sample_start = sample_start)
    attr(plan, ".shock_scale_spec")  <- plan_to_shock_scale_spec(plan)
    stored_plan <- plan
  }

  ## Soft warning: pkf with me_variance = 0 is technically allowed but near-
  ## singular for n_obs == n_state models.
  if (likelihood == "pkf" && isTRUE(me_variance == 0)) {
    warning(
      "estimation_context: likelihood = \"pkf\" with me_variance = 0 may ",
      "cause near-singular F matrices when n_obs >= n_state. ",
      "Consider me_variance = 1e-8 or larger.",
      call. = FALSE
    )
  }

  ## Validation: student_t requires student_df > 0
  if (likelihood == "student_t") {
    if (is.null(student_df))
      stop("estimation_context: likelihood = \"student_t\" requires student_df ",
           "(degrees of freedom, a positive scalar).", call. = FALSE)
    if (!is.numeric(student_df) || length(student_df) != 1L ||
        !is.finite(student_df) || student_df <= 0)
      stop("estimation_context: student_df must be a positive finite scalar.",
           call. = FALSE)
  }

  ## Validation: tpf requires me_variance > 0
  if (likelihood == "tpf") {
    if (!is.numeric(me_variance) || length(me_variance) != 1L ||
        !is.finite(me_variance) || me_variance <= 0) {
      stop(
        "estimation_context: likelihood = \"tpf\" requires me_variance > 0.\n",
        "The Tempered Particle Filter uses me_variance as its tempering ",
        "instrument; me_variance = 0 is not allowed.\n",
        "Provide me_variance > 0 (e.g. 1e-4 * var(data))."
      )
    }
  }

  ## Validation: freq_band outside c(0,pi) with non-whittle likelihood
  if (!identical(likelihood, "whittle") &&
      (!isTRUE(all.equal(freq_band, c(0, pi))))) {
    warning(
      "estimation_context: freq_band != c(0, pi) has no effect when ",
      "likelihood != \"whittle\".",
      call. = FALSE
    )
  }

  ## Validation: exactly one of ms_spec / ms_struct_spec (or neither)
  if (!is.null(ms_spec) && !is.null(ms_struct_spec))
    stop("estimation_context: supply at most one of 'ms_spec' (shock-variance ",
         "switching) or 'ms_struct_spec' (structural switching), not both.",
         call. = FALSE)

  ## Validation: ms_spec compatibility checks (shock-variance switching)
  if (!is.null(ms_spec)) {
    if (!inherits(ms_spec, "ms_dsge_spec"))
      stop("estimation_context: ms_spec must be an ms_dsge_spec object.",
           call. = FALSE)
    if (!identical(likelihood, "gaussian"))
      stop("estimation_context: ms_spec is only compatible with ",
           "likelihood = \"gaussian\" (Kim-Nelson filter).", call. = FALSE)
    if (!is.null(me_extra))
      stop("estimation_context: ms_spec is incompatible with me_extra.",
           call. = FALSE)
    if (!is.null(shock_scale))
      stop("estimation_context: ms_spec is incompatible with shock_scale; ",
           "use shock_scales in the ms_dsge_spec instead.", call. = FALSE)
    ## Analytic gradient is not yet implemented for the MS likelihood
    if (identical(gradient_policy, "analytic"))
      stop("estimation_context: gradient_policy = \"analytic\" is not yet ",
           "supported for MS-DSGE (ms_spec != NULL). Use \"numerical\" or \"auto\".",
           call. = FALSE)
  }

  ## Validation: ms_struct_spec compatibility checks (structural switching)
  if (!is.null(ms_struct_spec)) {
    if (!inherits(ms_struct_spec, "ms_struct_spec"))
      stop("estimation_context: ms_struct_spec must be an ms_struct_spec object.",
           call. = FALSE)
    if (!identical(likelihood, "gaussian"))
      stop("estimation_context: ms_struct_spec is only compatible with ",
           "likelihood = \"gaussian\" (Kim-Nelson filter).", call. = FALSE)
    if (!is.null(me_extra))
      stop("estimation_context: ms_struct_spec is incompatible with me_extra.",
           call. = FALSE)
    if (!is.null(shock_scale))
      stop("estimation_context: ms_struct_spec is incompatible with shock_scale.",
           call. = FALSE)
    if (identical(gradient_policy, "analytic"))
      stop("estimation_context: gradient_policy = \"analytic\" is not yet ",
           "supported for structural MS-DSGE (ms_struct_spec != NULL). ",
           "Use \"numerical\" or \"auto\".", call. = FALSE)
  }

  structure(
    list(
      me_variance     = me_variance,
      likelihood      = likelihood,
      lik_init        = lik_init,
      me_extra        = me_extra,
      shock_scale     = shock_scale,
      freq_band       = freq_band,
      system_priors   = system_priors,
      tpf_options     = tpf_options,
      obc_options     = obc_options,
      obc_specs       = obc_specs,
      gradient_policy = gradient_policy,
      plan            = stored_plan,
      ms_spec         = ms_spec,
      ms_struct_spec  = ms_struct_spec,
      student_df      = student_df
    ),
    class = c("dynhr_estimation_context", "list")
  )
}


#' @export
#' @noRd
print.dynhr_estimation_context <- function(x, ...) {
  cat("<dynhr_estimation_context>\n")
  cat(sprintf("  likelihood    : %s\n", x$likelihood))
  cat(sprintf("  me_variance   : %g\n", x$me_variance))
  cat(sprintf("  lik_init      : %s\n", x$lik_init))
  cat(sprintf("  freq_band     : [%.4g, %.4g]\n", x$freq_band[1], x$freq_band[2]))
  cat(sprintf("  me_extra      : %s\n",
              if (is.null(x$me_extra)) "NULL"
              else sprintf("matrix [%d x %d]", nrow(x$me_extra), ncol(x$me_extra))))
  cat(sprintf("  shock_scale   : %s\n",
              if (is.null(x$shock_scale)) "NULL"
              else sprintf("matrix [%d x %d]", nrow(x$shock_scale), ncol(x$shock_scale))))
  cat(sprintf("  system_priors : %s\n",
              if (is.null(x$system_priors)) "NULL"
              else sprintf("%d prior(s)", length(x$system_priors))))
  cat(sprintf("  tpf_options   : %s\n",
              if (length(x$tpf_options) == 0) "default"
              else paste(names(x$tpf_options), collapse = ", ")))
  cat(sprintf("  obc_options   : %s\n",
              if (length(x$obc_options) == 0) "default"
              else paste(names(x$obc_options), collapse = ", ")))
  cat(sprintf("  obc_specs     : %s\n",
              if (is.null(x$obc_specs)) "NULL (parsed on demand)"
              else sprintf("%d constraint(s)", length(x$obc_specs))))
  cat(sprintf("  gradient_policy: %s\n", x$gradient_policy))
  cat(sprintf("  plan          : %s\n",
              if (is.null(x$plan)) "NULL" else "<dynhr_plan>"))
  cat(sprintf("  ms_spec       : %s\n",
              if (is.null(x$ms_spec)) "NULL"
              else sprintf("<ms_dsge_spec> %d regimes", x$ms_spec$n_regimes)))
  cat(sprintf("  ms_struct_spec: %s\n",
              if (is.null(x$ms_struct_spec)) "NULL"
              else sprintf("<ms_struct_spec> %d regimes (structural)", x$ms_struct_spec$n_regimes)))
  if (!is.null(x$student_df))
    cat(sprintf("  student_df    : %g\n", x$student_df))
  invisible(x)
}


#' Check whether a context allows the analytic gradient path
#'
#' Compatibility test: returns \code{TRUE} when nothing in the context rules
#' the analytic (tangent/adjoint KF) gradient path out -- standard Gaussian
#' likelihood and gradient_policy not explicitly \code{"numerical"}.  The
#' analytic gradient now supports per-period \code{me_extra} and
#' \code{shock_scale} inputs (ROADMAP Tier 7 item 3), so those are no longer
#' blocking conditions.  Whether to USE the analytic path is the caller's
#' \code{analytic_grad} flag; this predicate only answers whether it is allowed.
#' Replaces three duplicated inline guards in the runner files.
#'
#' @param ctx An \code{estimation_context} object.
#' @return Logical scalar.
#' @noRd
.ctx_allows_analytic_gradient <- function(ctx) {
  ## Unit-root draws already return logpost = -Inf before the gradient is
  ## evaluated (the spectral_radius >= 1 guard in make_log_posterior_whittle),
  ## so the near-unit-root regime is self-limiting for whittle.
  ## cumulant: make_posterior_grad(likelihood = "cumulant") mirrors the forward
  ## order-2 solve and dispatches to cumulant_loglik_grad (Tier 14 B2).
  ## pskf/tpf are excluded: those likelihoods have no analytic gradient path.
  ## Note: "pruned" is NOT listed here -- analytic gradient is not yet
  ## implemented for the pruned-SS KF (numerical FD fallback applies).
  !identical(ctx$gradient_policy, "numerical") &&
    ctx$likelihood %in% c("gaussian", "whittle", "cumulant")
}


#' Check whether a context is the "standard Gaussian" fast path
#'
#' The standard Gaussian fast path allows the daemon to recompile the posterior
#' once per worker (via \code{.mirai_pool_init}) instead of shipping a pre-built
#' closure.  This replaces three slightly-different \code{par_standard}
#' computations in the runner files.
#'
#' @param ctx     An \code{estimation_context} object.
#' @param use_obc Logical; \code{TRUE} if the model has OBC/PKF constraints.
#' @return Logical scalar.
#' @noRd
.ctx_is_standard_gaussian <- function(ctx, use_obc = FALSE) {
  !isTRUE(use_obc) &&
    identical(ctx$likelihood, "gaussian") &&
    is.null(ctx$me_extra) &&
    is.null(ctx$shock_scale) &&
    isTRUE(all.equal(ctx$freq_band, c(0, pi)))
}


#' Validate a context for mirai-safe serialization
#'
#' Checks that the \code{system_priors} closures in \code{ctx} do not close
#' over large non-base environments.  Closures that capture large objects
#' (e.g. a full estimation frame or an \code{R6} object) will cause mirai to
#' serialize the entire environment, inflating the task payload and potentially
#' breaking serialization.  Safe closures close over scalar constants and base
#' environments only.
#'
#' @param ctx An \code{estimation_context} object.
#' @return \code{ctx} invisibly, after issuing warnings for any unsafe closure.
#' @export
validate_context <- function(ctx) {
  stopifnot(inherits(ctx, "dynhr_estimation_context"))

  sp <- ctx$system_priors
  if (is.null(sp)) return(invisible(ctx))

  safe_envs <- c("R_GlobalEnv", "base", "package:dynhr", "",
                 grep("^package:", search(), value = TRUE))

  for (nm in names(sp)) {
    fn <- sp[[nm]]
    if (!is.function(fn)) next
    env_nm <- tryCatch(environmentName(environment(fn)),
                       error = function(e) NA_character_)
    if (!is.na(env_nm) && !(env_nm %in% safe_envs)) {
      warning(
        "validate_context: system_priors[[\"", nm, "\"]] closes over a ",
        "non-base environment ('", env_nm, "'). This may cause mirai ",
        "serialization overhead or failures. Ensure the closure only captures ",
        "scalar constants.",
        call. = FALSE
      )
    }
  }

  invisible(ctx)
}


#' Reconstruct an estimation context from a legacy mode result
#'
#' When \code{run_posterior_estimation()} receives a \code{dynhr_mode_result}
#' that pre-dates the context refactor (i.e. has no \code{$ctx} field), this
#' helper reconstructs a context from the legacy flat fields.
#'
#' @param mr A \code{dynhr_mode_result} object.
#' @return An \code{estimation_context} object.
#' @noRd
ctx_from_mode_result <- function(mr) {
  estimation_context(
    me_variance   = mr$me_variance   %||% 0,
    likelihood    = mr$likelihood    %||% "gaussian",
    lik_init      = mr$lik_init      %||% "auto",
    me_extra      = mr$me_extra      %||% NULL,
    shock_scale   = mr$shock_scale   %||% NULL,
    freq_band     = mr$freq_band     %||% c(0, pi),
    system_priors = mr$system_priors %||% NULL,
    tpf_options   = mr$tpf_options   %||% list()
  )
}
