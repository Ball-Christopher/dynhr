## R/posterior.R
## --------------------------------------------------------------------------
## Phase-2 split from estimation-monolith.R.
##
## make_log_posterior(): closure factory called once before MCMC.
## Returns a function(theta) -> list(logpost, loglik, logprior).
## --------------------------------------------------------------------------

## Helper: build shock covariance matrix Sigma_e from model parameters.
## Handles both diagonal variances and cross-shock correlations from the
## `shocks;` block (parsed by parse_shocks_block() in R/parse-blocks.R).
##
## corr entries: if the row has a corr_expr column (non-NA), re-evaluate the
## expression against the current params so that an estimated correlation
## parameter rho_ab tracks theta during MCMC. The cov entry is then
##   cov_ij = rho_ij * sd_i * sd_j
## where sd_i / sd_j come from the ALREADY re-evaluated stderr vector, so the
## off-diagonal also picks up any moving stderr parameters. If |corr| > 1 at
## some theta the resulting matrix will not be PSD; chol() inside the Kalman
## filter will fail and the posterior will return -Inf naturally -- consistent
## with how negative variances from variance_expr are handled (no stop()).
.get_shock_cov <- function(model, exo_names, params) {
  n_exo   <- length(exo_names)
  stderr  <- .get_shock_stderr(model, exo_names, params)
  Sigma_e <- diag(stderr^2, nrow = n_exo)
  rownames(Sigma_e) <- colnames(Sigma_e) <- exo_names

  ## Incorporate cross-shock correlations from the parsed shocks block.
  corr_df <- model$shocks$correlations
  if (!is.null(corr_df) && nrow(corr_df) > 0) {
    ## Build a parameter eval env once (lazy, only if corr_expr present).
    penv <- NULL
    get_penv <- function() {
      if (is.null(penv))
        penv <<- list2env(as.list(params), parent = baseenv())
      penv
    }
    for (k in seq_len(nrow(corr_df))) {
      r <- corr_df[k, ]
      i <- match(r$var1, exo_names)
      j <- match(r$var2, exo_names)
      if (is.na(i) || is.na(j)) next

      ## Distinguish cov-type rows (absolute covariance; from `var a,b=expr`)
      ## from corr-type rows (correlation ratio; from `corr a,b=expr`).
      has_cov_expr <- "cov_expr" %in% names(r) && !is.na(r$cov_expr) &&
                      nzchar(as.character(r$cov_expr))
      has_cov_val  <- "cov" %in% names(r) && !is.na(r$cov)

      if (has_cov_expr || has_cov_val) {
        ## cov-type: re-evaluate the expression; fall back to parse-time value.
        cov_val <- NA_real_
        if (has_cov_expr) {
          cov_val <- tryCatch(
            .eval_cached_expr(as.character(r$cov_expr), get_penv()),
            error = function(e) NA_real_
          )
        }
        if (is.na(cov_val) || !is.finite(cov_val)) cov_val <- r$cov
        if (is.finite(cov_val))
          Sigma_e[i, j] <- Sigma_e[j, i] <- cov_val
      } else {
        ## corr-type: re-evaluate corr_expr; fall back to parse-time snapshot.
        rho <- NA_real_
        if ("corr_expr" %in% names(r) && !is.na(r$corr_expr) &&
            nzchar(r$corr_expr)) {
          rho <- tryCatch(
            .eval_cached_expr(as.character(r$corr_expr), get_penv()),
            error = function(e) NA_real_
          )
        }
        if (is.na(rho) || !is.finite(rho)) rho <- r$corr

        if (is.finite(rho)) {
          Sigma_e[i, j] <- Sigma_e[j, i] <- rho * stderr[i] * stderr[j]
        }
      }
    }
  }

  Sigma_e
}

#' Shock covariance matrix implied by a parsed model
#'
#' Public accessor that builds the \code{n_exo x n_exo} structural-shock
#' covariance matrix \eqn{\Sigma_e} from a model's \code{shocks} block,
#' honouring \code{stderr}, \code{var a,b = cov}, and \code{corr a,b}
#' (including expression-valued entries re-evaluated against \code{params}).
#' This is the same matrix the Kalman filter, Whittle likelihood, and the
#' spectral/identification diagnostics (D13, D15, D23, D24) use internally, so
#' it is the correct value to pass as their \code{sigma_e=} argument.
#'
#' @param model A parsed model object (from \code{\link{parse_mod}}).
#' @param params Named numeric parameter vector to evaluate expression-valued
#'   covariance/correlation entries against. Defaults to
#'   \code{model$param_values}.
#' @param exo_names Character vector of exogenous shock names, in the order the
#'   rows/columns of \eqn{\Sigma_e} should take. Defaults to
#'   \code{model$varexo_names}.
#' @return A named \code{n_exo x n_exo} covariance matrix.
#' @examples
#' \dontrun{
#' m  <- parse_mod("model.mod")
#' Se <- shock_cov(m)
#' diag_cer(dr = dr, data = Y, model = m, sigma_e = Se)
#' }
#' @export
shock_cov <- function(model, params = NULL, exo_names = NULL) {
  if (is.null(params))    params    <- model$param_values
  if (is.null(exo_names)) exo_names <- model$varexo_names
  .get_shock_cov(model, exo_names, params)
}

## Apply an estimation draw `theta` onto the model parameter vector.
##
## Model parameters are overwritten in place. Crucially, a `theta` entry whose
## name is NOT a model parameter but IS an exogenous shock name is also injected
## into `params` under the shock name: this is the Dynare `stderr <shock>` /
## `corr <a>,<b>` estimated_params convention, where the estimated quantity is a
## shock standard deviation (or correlation), not a structural parameter.
## `.get_shock_stderr` / `.get_shock_cov` then pick it up (Priority 0 / corr
## lookup). Without this the estimated shock std is silently dropped and the
## likelihood is computed with the frozen parse-time stderr (P0 bug: the shock
## stds have zero effect on the loglik). For models that estimate stds the
## Dynare-rebuilt `sig_*` way, those ARE model params and take the first branch.
##
## @param model  dynhr_mod.
## @param theta  Named numeric draw (from prior_spec ordering).
## @param params Optional starting vector (default model$param_values).
## @return params with theta applied (structural params + injected shock stds).
#' Apply an estimated parameter vector to a model's parameter values
#'
#' Merges a named vector of ESTIMATED quantities (a posterior draw, a
#' posterior mode, a swept point) into a model's parameter vector, returning
#' the result. This is the conversion that sits between an estimation result
#' and the solver: the output is what you pass as
#' \code{solve_steady_state(model, params = ...)}.
#'
#' It exists as a public route because
#' \code{\link{set_param_values}} cannot do this job: that function ERRORS on
#' any name outside \code{model$param_names}, whereas an estimated vector
#' routinely also carries SHOCK STANDARD DEVIATIONS under Dynare's
#' \code{stderr <shock>} convention. Those arrive named for the shock, not for
#' a parameter, and are injected under the shock's own name so the downstream
#' shock-covariance builder picks them up as that shock's standard deviation.
#'
#' @param model A parsed model (see \code{\link{parse_mod}}).
#' @param theta Named numeric vector of estimated values. Names may be model
#'   parameters, shock names (an estimated \code{stderr <shock>}), or
#'   correlation entries such as \code{"corr e_a,e_b"}.
#' @param params Optional named parameter vector to merge INTO. \code{NULL}
#'   (default) starts from \code{model$param_values}. Pass an existing vector
#'   to layer several updates without going back to the model each time.
#'
#' @return The merged named numeric parameter vector. \code{model} itself is
#'   NOT modified.
#'
#' @section Names that are deliberately passed over:
#' Entries of \code{theta} that are neither a model parameter nor a declared
#' shock -- in practice \code{corr <a>,<b>} entries -- are left OUT of the
#' returned vector by design: a correlation is not a parameter value and is
#' consumed separately when the shock covariance is assembled. An entry that is
#' none of the three is also passed over here rather than raised, because the
#' builder-time guard in \code{\link{make_log_posterior}} is the layer that
#' rejects an unusable prior target, and duplicating that check here would make
#' the same mistake fail in two places with different messages. **If you are
#' calling this directly, compare \code{names(theta)} against the returned
#' names when you need to be sure nothing was dropped silently.**
#'
#' @seealso \code{\link{set_param_values}} (parameters only, strict),
#'   \code{\link{solve_steady_state}}, \code{\link{make_log_posterior}}
#' @examples
#' \donttest{
#' mod <- system.file("extdata", "models", "rbc", "rbc.mod", package = "dynhr")
#' if (nzchar(mod)) {
#'   m  <- parse_mod(mod)
#'   p  <- apply_theta_to_params(m, c(alpha = 0.33))
#'   ss <- solve_steady_state(m, params = p)
#'   ss$converged
#' }
#' }
#' @export
apply_theta_to_params <- function(model, theta, params = NULL) {
  if (is.null(params)) params <- model$param_values
  exo <- model$varexo_names %||% character(0)
  for (nm in names(theta)) {
    if (nm %in% names(params)) {
      params[nm] <- theta[nm]
    } else if (nm %in% exo) {
      ## Estimated shock std (Dynare `stderr <shock>`): inject under the shock
      ## name so .get_shock_stderr Priority 0 uses it as that shock's stderr.
      params[nm] <- theta[nm]
    }
    ## else: a corr <a>,<b> entry or an unknown name — left for .get_shock_cov's
    ## correlation handling / the builder-time validation guard to deal with.
  }
  params
}

## Internal alias, retained so the many existing internal call sites keep
## working. A direct binding, NOT a wrapper: this runs once per posterior
## evaluation, so an extra frame would be pure overhead, and two
## implementations could drift.
.apply_theta_to_params <- apply_theta_to_params

## Validate that every estimated prior maps to something the likelihood uses.
## Fail loudly (rather than silently dropping) when a prior name is neither a
## model parameter, nor a shock name (estimated std), nor a `corr a,b` entry.
## This closes the P0 silent-drop class: an estimated quantity that the
## likelihood cannot apply must error at builder time, not be ignored.
.validate_prior_targets <- function(model, prior_spec, where = "make_log_posterior") {
  if (is.null(prior_spec) || is.null(prior_spec$name)) return(invisible(NULL))
  pn  <- names(model$param_values) %||% character(0)
  exo <- model$varexo_names %||% character(0)
  is_corr <- grepl("[, ]", prior_spec$name) |
    grepl("^corr[_.]", prior_spec$name, ignore.case = TRUE)
  unknown <- prior_spec$name[!(prior_spec$name %in% pn) &
                               !(prior_spec$name %in% exo) & !is_corr]
  if (length(unknown) > 0L) {
    stop(sprintf(
      paste0("%s: %d estimated prior(s) are not connected to the likelihood ",
             "(not a model parameter, an exogenous shock std, or a correlation): ",
             "%s. A `stderr <shock>` prior must name a declared shock; a ",
             "structural prior must name a declared parameter. Fix the ",
             "estimated_params block or rename the prior."),
      where, length(unknown), paste(unknown, collapse = ", ")),
      call. = FALSE)
  }
  invisible(NULL)
}

## Shared solve pipeline: theta -> (dr, params) or NULL on infeasibility.
##
## Factored out of make_log_posterior's gaussian-path closure so that
## make_loglik_contrib() (R/posterior.R) can reuse the EXACT same
## steady-state-solve -> BK-check -> stationarity-guard pipeline rather than
## re-deriving a subtly different one (the source of the make_posterior vs
## hand-rolled-pipeline loglik discrepancy on nk_small this closes).
##
## `state` is a mutable environment holding the per-closure warm-start cache
## (`state$ss_warm`), so repeated calls across theta draws keep warm-starting
## the steady-state solve exactly like make_log_posterior does.
##
## @param model, compiled, sys_cache As in make_log_posterior.
## @param theta      Named numeric draw.
## @param state      environment with a `ss_warm` field (mutable cache).
## @param shock_scale Passed through only to decide whether a near-unit-root
##   draw must be rejected outright (heteroskedastic shocks are incompatible
##   with the diffuse phase) -- mirrors make_log_posterior's guard exactly.
## @param lik_init   As in kalman_filter(); used only to resolve which init
##   would be "in force" for the stationarity guard (mirrors
##   make_log_posterior's resolution so the SAME draws are rejected).
## @return list(dr = <decision rules>, params = <params>) on success, or
##   NULL on infeasibility (BK violation or unit root under a stationary
##   init).
.solve_dr_for_theta <- function(model, compiled, sys_cache, theta, state,
                                lik_init = "auto", shock_scale = NULL) {
  params <- .apply_theta_to_params(model, theta)

  ss_result <- solve_steady_state(model, compiled, params,
                                  y0 = state$ss_warm, verbose = FALSE)
  if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
    ## Retry once from the cold initval-based guess before declaring
    ## infeasible (mirrors make_log_posterior's warm-start retry).
    if (!is.null(state$ss_warm))
      ss_result <- solve_steady_state(model, compiled, params,
                                      verbose = FALSE)
    if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
      state$ss_warm <- NULL
      return(NULL)
    }
  }
  state$ss_warm <- ss_result$ss

  ## Re-derive any steady_state_model-computed parameter (Tier 13 #1 fix).
  params <- ss_result$params %||% params
  sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
  dr  <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
  if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)

  ## Stationarity guard identical to make_log_posterior's gaussian path.
  ns <- length(dr$state_idx)
  ev <- dr$eigenvalues
  spectral_radius <- if (!is.null(ev) && length(ev) >= ns)
    max(Mod(ev[seq_len(ns)]))
  else
    max(Mod(eigen(dr$ghx[dr$state_idx, , drop = FALSE],
                  only.values = TRUE)$values))
  init_in_force <- if (identical(lik_init, "auto"))
    (if (spectral_radius > 1 - 1e-6) "diffuse" else "stationary")
  else lik_init
  if (spectral_radius >= 1 &&
      (identical(init_in_force, "stationary") || !is.null(shock_scale))) {
    return(NULL)
  }

  list(dr = dr, params = params)
}

#' Create a cached log-posterior evaluator for MCMC
#'
#' Call this ONCE before MCMC. Returns a closure that evaluates the
#' log-posterior for any parameter vector theta, using pre-cached
#' system structure to avoid repeated compilation.
#'
#' @param model       dynhr_mod (from parse_mod())
#' @param data        Observation matrix (T x n_obs), columns = obs_vars
#' @param prior_spec  Prior specification data.frame (from extract_prior_spec())
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled (from compile_model())
#' @param me_variance Measurement error variance (default 0). Use 0 for
#'   exactly-identified models (n_obs == n_shocks); HH = DD*Sigma_e*DD'
#'   already regularises F_t. Set > 0 for stochastically-singular models.
#' @param likelihood  Likelihood type: \code{"gaussian"} (Kalman filter,
#'   default), \code{"cumulant"} (cumulant-matching, Mutschler 2015),
#'   \code{"whittle"} (frequency-domain Whittle likelihood), or
#'   \code{"pskf"} (Pruned Skewed Kalman Filter for skew-normal shocks).
#'   The cumulant likelihood works with any perturbation order >= 1;
#'   orders >= 2 provide skewness/kurtosis content.
#'   The Whittle likelihood requires a stationary, complete (no NA) panel;
#'   it supports band-restricted estimation via \code{freq_band}.
#'   The PSKF uses the CSN state-space recursion; set \code{skew} in the
#'   shocks block (or \code{estimated_params}) to specify shock shapes.
#' @param lik_init    Kalman filter \code{P0} initialization, forwarded to
#'   \code{\link{kalman_filter}} on the \code{likelihood = "gaussian"} path
#'   (default \code{"auto"}: stationary models use the Lyapunov \code{P0};
#'   unit-root models use the exact diffuse initialization). See
#'   \code{\link{kalman_filter}} for \code{"stationary"}, \code{"diffuse"},
#'   and \code{"kappa"}. Ignored when \code{likelihood = "cumulant"} or
#'   \code{likelihood = "whittle"}.
#' @param freq_band   Numeric(2) \code{c(lo, hi)} in radians (only used when
#'   \code{likelihood = "whittle"}). Restricts the Whittle sum to Fourier
#'   frequencies in \code{(lo, hi]}. Default \code{c(0, pi)} uses all
#'   frequencies. See \code{\link{whittle_business_cycle_band}}.
#' @param power        Power-posterior (generalised-Bayes) tempering exponent
#'   \eqn{\zeta \in (0, 1]}. The returned log-posterior is
#'   \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)}: setting
#'   \eqn{\zeta < 1} down-weights a (potentially misspecified) likelihood and
#'   widens the posterior in a calibration-consistent way (Bissiri, Holmes &
#'   Walker 2016, JRSS-B; Grünwald SafeBayes). Default \code{1} is
#'   bit-identical to the untempered posterior. The \code{$loglik} field in
#'   the returned list always carries the \emph{raw} (untempered) log-likelihood
#'   so that SMC tempering, marginal-likelihood estimates, and other callers that
#'   need the true likelihood are unaffected. Only \code{$logpost} is tempered.
#'   Can also be set globally via
#'   \code{dynhr_set_options(power_posterior = 0.5)}.
#' @param pruned_order For \code{likelihood = "pruned"}: the perturbation order
#'   of the AFVRR pruned state space, \code{2L} (default) or \code{3L}. Order 3
#'   uses \code{pruned_ss_loglik3} (skewness/kurtosis content via the cubic
#'   augmented state). Ignored for other likelihoods.
#' @param ...         Additional arguments passed to the cumulant likelihood
#'   constructor when \code{likelihood = "cumulant"} (e.g. \code{order},
#'   \code{cumulant_orders}, \code{cumulant_weight}), or to
#'   \code{make_log_posterior_tpf} when \code{likelihood = "tpf"} (e.g.
#'   \code{n_particles}, \code{ess_target}, \code{n_mh}, \code{seed}).
#' @param me_extra Optional \code{n_obs x T} matrix of ADDITIONAL per-period
#'   measurement-error variance, added on top of \code{me_variance}. Rows are
#'   observables in \code{obs_vars} order and columns are periods, so a row
#'   swap is a different model -- see the multi-value testing note in
#'   \code{CLAUDE.md}.
#' @param shock_scale Optional \code{n_exo x T} matrix of KNOWN per-period
#'   shock standard-deviation scale factors (the deterministic-volatility
#'   path). Rows are shocks in declaration order.
#' @param system_priors Optional \code{system_prior_spec} object placing priors
#'   on model-implied quantities (moments, IRF features) rather than on
#'   parameters directly; its log-density is added to the parameter prior.
#' @param infeasible_penalty Optional soft penalty replacing the hard
#'   \code{-Inf} returned where the model cannot be solved, which lets
#'   gradient-based samplers see a slope out of an infeasible region instead of
#'   a cliff. \code{NULL} (default) is OFF and byte-identical to the
#'   unpenalised behaviour. A single number is taken as
#'   \code{list(scale = x, floor = -1e6)}; a list must supply \code{scale} and
#'   may supply \code{floor}.
#' @param ctx Optional \code{\link{estimation_context}} bundling
#'   \code{me_variance}, \code{likelihood}, \code{lik_init}, \code{me_extra},
#'   \code{shock_scale}, \code{freq_band}, \code{system_priors},
#'   \code{infeasible_penalty} and the Markov-switching specs. **When supplied,
#'   its fields OVERRIDE the individual arguments silently**, so pass a context
#'   or the individual arguments, not both.
#' @param student_df Positive finite scalar: degrees of freedom for
#'   \code{likelihood = "student_t"}. Ignored by the other likelihoods.
#'
#' @return A closure \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}, where \code{theta} is the named
#'   vector of ESTIMATED parameters described by \code{prior_spec}.
#'
#'   The closure captures a compiled model, so it holds external pointers and
#'   **does not survive \code{saveRDS}/\code{readRDS} or transport to a worker
#'   process**. Rebuild it in the target session by calling this function again
#'   with the same arguments rather than serialising the closure.
#'
#' @seealso \code{\link{estimation_context}}, \code{\link{extract_prior_spec}},
#'   \code{\link{log_prior}}
#' @export
make_log_posterior <- function(model, data, prior_spec, obs_vars = NULL,
                               compiled, me_variance = 0,
                               likelihood = c("gaussian", "cumulant",
                                              "whittle", "tpf", "pskf",
                                              "student_t", "pruned",
                                              "ppf", "copf", "sv_rbpf"),
                               lik_init = "auto",
                               me_extra = NULL,
                               shock_scale = NULL,
                               freq_band = c(0, pi),
                               system_priors = NULL,
                               infeasible_penalty = NULL,
                               ctx = NULL,
                               power = NULL,
                               student_df = NULL,
                               pruned_order = 2L,
                               ...) {
  ## When a ctx is supplied, unpack its fields over the individual args.
  ## Individual args supplied alongside ctx are silently overridden by ctx.
  ms_spec_ctx        <- NULL
  ms_struct_spec_ctx <- NULL
  if (!is.null(ctx) && inherits(ctx, "dynhr_estimation_context")) {
    me_variance        <- ctx$me_variance
    likelihood         <- ctx$likelihood
    lik_init           <- ctx$lik_init
    me_extra           <- ctx$me_extra
    shock_scale        <- ctx$shock_scale
    freq_band          <- ctx$freq_band
    system_priors      <- ctx$system_priors
    infeasible_penalty <- ctx$infeasible_penalty %||% infeasible_penalty
    ms_spec_ctx        <- ctx$ms_spec        # shock-variance MS path
    ms_struct_spec_ctx <- ctx$ms_struct_spec # structural MS path (new)
    student_df         <- ctx$student_df %||% student_df
    pruned_order       <- ctx$pruned_order %||% pruned_order
    ## tpf_options are merged into ... via do.call below when likelihood="tpf"
  }
  likelihood <- match.arg(likelihood)
  if (!pruned_order %in% c(2L, 3L))
    stop("make_log_posterior: `pruned_order` must be 2 or 3.")

  ## Default obs_vars from the model's varobs declaration (parse_mod also
  ## exposes it as model$obs_vars) — pathological-DSGE paper gap #4.
  if (is.null(obs_vars) || length(obs_vars) == 0L) {
    obs_vars <- model$obs_vars %||% model$varobs_names
    if (is.null(obs_vars) || length(obs_vars) == 0L)
      stop("make_log_posterior: `obs_vars` not supplied and the model ",
           "declares no `varobs`; pass obs_vars= or add a varobs line to ",
           "the .mod.", call. = FALSE)
  }

  ## Power-posterior (generalised-Bayes) tempering exponent zeta in (0, 1].
  ## Resolves: explicit arg > global option `power_posterior` > default 1.
  ## Default 1 is bit-identical to the untempered posterior.
  power <- .dynhr_opt("power_posterior", power, default = 1)
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0) {
    stop("make_log_posterior: `power` must be a finite scalar in (0, 1].")
  }
  if (power > 1) {
    warning("make_log_posterior: `power` > 1 produces a 'cold' (over-confident) ",
            "posterior. This is valid but unusual; set power <= 1 for standard ",
            "generalised-Bayes tempering.")
  }

  ## Continuous Blanchard-Kahn violation measure for the infeasibility penalty.
  ## The discrete count of excess unstable roots has NO gradient (it is constant
  ## between integer changes), so it gives the optimizer no direction toward
  ## feasibility. Instead measure how far the BOUNDARY-crossing generalized
  ## eigenvalues are from the unit circle: BK requires exactly `model$n_forward`
  ## eigenvalues with modulus > 1 (infinite roots count as unstable); a violation
  ## means some eigenvalue is on the wrong side. The penalty is the sum of squared
  ## distances of the misclassified eigenvalues from |lambda| = 1 -- continuous in
  ## the parameters, so as the optimizer moves the offending eigenvalue across the
  ## unit circle the penalty falls smoothly to zero (then the real likelihood
  ## takes over). Falls back to the discrete excess count when eigenvalues are
  ## unavailable. (User-requested: penalise eigenvalue crossings, not the count.)
  .bk_eigen_violation <- function(dr, model) {
    ev     <- if (!is.null(dr)) dr$eigenvalues else NULL
    target <- model$n_forward                     # required # of unstable roots
    if (is.null(ev) || length(ev) == 0L || is.null(target) || !is.finite(target)) {
      nu <- if (!is.null(dr)) dr$n_unstable else NULL
      return(if (!is.null(nu) && !is.null(target))
               as.numeric(abs(nu - target)) + 1e-6 else 1)
    }
    mm  <- Mod(ev)
    n_inf <- sum(!is.finite(mm))                  # infinite roots are unstable
    mm  <- sort(mm[is.finite(mm)], decreasing = TRUE)
    need <- max(as.integer(target) - n_inf, 0L)   # finite unstable roots needed
    n   <- length(mm)
    if (n == 0L) return(as.numeric(abs(target - n_inf)) + 1e-6)
    k   <- min(need, n)
    top <- if (k >= 1L) mm[seq_len(k)] else numeric(0)   # these should be > 1
    bot <- if (n > k)  mm[(k + 1L):n] else numeric(0)    # these should be <= 1
    viol <- sum(pmax(1 - top, 0)^2) + sum(pmax(bot - 1, 0)^2)
    ## Guarantee a strictly-positive penalty whenever BK actually failed, even in
    ## the rare case the count matches but the ordering is degenerate.
    if (viol <= 0) viol <- 1e-6
    viol
  }

  ## Validate and normalise infeasible_penalty.
  ## Accepted forms:
  ##   NULL          -> default OFF (byte-identical to existing behaviour)
  ##   numeric(1)    -> treat as list(scale = x, floor = -1e6)
  ##   list(scale=x) -> use scale; floor defaults to -1e6
  .penalty_cfg <- NULL
  if (!is.null(infeasible_penalty)) {
    if (is.numeric(infeasible_penalty) && length(infeasible_penalty) == 1L) {
      infeasible_penalty <- list(scale = infeasible_penalty)
    }
    if (!is.list(infeasible_penalty) || is.null(infeasible_penalty$scale)) {
      stop("make_log_posterior: infeasible_penalty must be NULL, a numeric scale, ",
           "or list(scale = <number>, floor = <number>).")
    }
    .penalty_cfg <- list(
      scale = infeasible_penalty$scale,
      floor = infeasible_penalty$floor %||% -1e6
    )
  }

  ## Fail loudly if any estimated prior cannot be applied to the likelihood
  ## (P0 guard: a `stderr <shock>` prior must name a declared shock; a
  ## structural prior must name a declared parameter). Prevents silently
  ## optimizing a likelihood in which some estimated parameters are inert.
  .validate_prior_targets(model, prior_spec, where = "make_log_posterior")
  ## NB: full prior_spec structure validation lives in extract_prior_spec() (the
  ## canonical producer), NOT here -- callers legitimately pass minimal / alternate-
  ## schema specs straight to make_log_posterior (e.g. the tpf early-stop tests),
  ## and log_prior() still fails loud at evaluation on an unknown distribution.

  if (likelihood == "tpf") {
    ## Hard stop on me_variance <= 0: the TPF requires positive measurement
    ## error variance as its tempering instrument (Landmine 1).
    if (!is.numeric(me_variance) || length(me_variance) != 1L ||
        !is.finite(me_variance) || me_variance <= 0) {
      stop(
        "make_log_posterior: likelihood = \"tpf\" requires me_variance > 0.\n",
        "The Tempered Particle Filter uses measurement error variance as its\n",
        "tempering instrument; me_variance = 0 is not allowed.\n",
        "Provide me_variance > 0 (e.g. 1e-4 * var(data))."
      )
    }
    ## Hard stop if model has filter_tunes: per-period extra ME variances are
    ## a Kalman-filter feature and are incompatible with TPF.
    ## model$filter_tunes is either a bare data.frame of tunes or a spec list
    ## whose rows live in $tunes. (NB a data.frame IS a list — test
    ## is.data.frame first.)
    ft <- model$filter_tunes
    ft_nrow <- if (is.data.frame(ft)) nrow(ft) else nrow(ft$tunes)
    if (!is.null(ft_nrow) && ft_nrow > 0) {
      stop(
        "make_log_posterior: likelihood = \"tpf\" is incompatible with ",
        "filter_tunes. Per-period measurement error adjustments (filter_tunes) ",
        "are a time-domain Kalman filter feature; TPF uses a fixed me_variance."
      )
    }
    ## Hard stop if shock_scale is supplied: heteroskedastic shocks require
    ## per-period Kalman filter updates and are incompatible with TPF.
    if (!is.null(shock_scale)) {
      stop(
        "make_log_posterior: likelihood = \"tpf\" is incompatible with ",
        "shock_scale (heteroskedastic_shocks). Use the Gaussian likelihood."
      )
    }
    ## data must be n_obs x T (TPF convention) — transpose if needed.
    ## make_log_posterior callers pass data T x n_obs; tpf expects n_obs x T.
    data_tpf <- if (ncol(data) == length(obs_vars)) t(data) else data
    ## Merge ctx$tpf_options (if any) with ... ; explicit ... wins over ctx.
    tpf_extra_args <- if (!is.null(ctx) && length(ctx$tpf_options) > 0L) {
      ## only use ctx$tpf_options keys not already in ...
      dots <- list(...)
      modifyList(ctx$tpf_options, dots)
    } else list(...)
    return(do.call(make_log_posterior_tpf,
                   c(list(model, data_tpf, prior_spec, obs_vars,
                          compiled, me_variance = me_variance,
                          system_priors = system_priors),
                     tpf_extra_args)))
  }

  if (likelihood == "sv_rbpf") {
    ## Measurement-side stochastic volatility on the shocks: Rao-Blackwellized
    ## particle filter over the latent AR(1) log-variance states, with the DSGE
    ## states integrated analytically via kf_step (see R/sv-rbpf.R). The SV spec
    ## already lives on `model$stochastic_volatility` (resolved by the runner);
    ## the factory's default stochastic_volatility = NULL uses it.
    if (!is.null(shock_scale))
      stop("make_log_posterior: likelihood = \"sv_rbpf\" is incompatible with a ",
           "deterministic shock_scale (heteroskedastic_shocks). The SV filter ",
           "supplies the shock scale from the latent volatility path itself.",
           call. = FALSE)
    ## data must be n_obs x T (particle-filter convention) — transpose if needed.
    data_sv <- if (ncol(data) == length(obs_vars)) t(data) else data
    ## Forward ONLY the SV factory's own tuning args from ... (run_full_estimation
    ## forwards its whole ... here, so a blind splat would pass unrelated args
    ## like n_iter into the strict factory signature and error).
    dots <- list(...)
    sv_extra_args <- dots[intersect(names(dots),
                                    c("n_particles", "seed", "stochastic_volatility"))]
    return(do.call(make_log_posterior_sv_rbpf,
                   c(list(model = model, data = data_sv, prior_spec = prior_spec,
                          obs_vars = obs_vars, compiled = compiled,
                          me_variance = me_variance, power = power),
                     sv_extra_args)))
  }

  if (likelihood == "ppf" || likelihood == "copf") {
    ## OBC particle filter (occasionally-binding-constraints). Bootstrap PF
    ## ("ppf") or conditionally-optimal PF ("copf"); both yield an UNBIASED
    ## log-marginal-likelihood estimate, so RWMH over this closure is PMMH for
    ## OBC models. This branch mirrors the "tpf" path: the same me_variance > 0
    ## guard (PF weights degenerate at me_variance = 0), the same
    ## incompatibility with filter_tunes / shock_scale (per-period Kalman
    ## features), and the same data transpose to n_obs x T.
    if (!is.numeric(me_variance) || length(me_variance) != 1L ||
        !is.finite(me_variance) || me_variance <= 0) {
      stop(
        "make_log_posterior: likelihood = \"", likelihood, "\" requires ",
        "me_variance > 0.\nThe OBC particle filter weights p(y_t | ...) ",
        "degenerate at me_variance = 0; provide me_variance > 0 (e.g. 1e-4)."
      )
    }
    ft      <- model$filter_tunes
    ft_nrow <- if (is.data.frame(ft)) nrow(ft) else nrow(ft$tunes)
    if (!is.null(ft_nrow) && ft_nrow > 0) {
      stop(
        "make_log_posterior: likelihood = \"", likelihood, "\" is incompatible ",
        "with filter_tunes. Per-period measurement error adjustments are a ",
        "Kalman filter feature; the OBC PF uses a fixed me_variance."
      )
    }
    if (!is.null(shock_scale)) {
      stop(
        "make_log_posterior: likelihood = \"", likelihood, "\" is incompatible ",
        "with shock_scale (heteroskedastic_shocks). Use the Gaussian likelihood."
      )
    }
    ## make_log_posterior_obc_ppf accepts data either orientation, but pass the
    ## n_obs x T convention explicitly for parity with the tpf branch.
    data_ppf <- if (ncol(data) == length(obs_vars)) t(data) else data
    ## The obc_ppf factory has a NARROWER interface than the Gaussian/TPF
    ## paths: it takes only (model, data, prior_spec, obs_vars, compiled,
    ## specs, me_variance, N, proposal, regime_guess, seed). It does NOT accept
    ## system_priors, power, lik_init, ctx, me_extra -- pass only its args. Any
    ## of {N, specs, proposal, regime_guess, seed} may be supplied via ... .
    dots      <- list(...)
    ppf_allow <- c("specs", "N", "proposal", "regime_guess", "seed")
    ppf_extra <- dots[intersect(names(dots), ppf_allow)]
    ## likelihood == "copf" selects the conditionally-optimal proposal unless
    ## the caller already named one explicitly in ... .
    if (is.null(ppf_extra$proposal) && likelihood == "copf")
      ppf_extra$proposal <- "copf"
    return(do.call(make_log_posterior_obc_ppf,
                   c(list(model, data_ppf, prior_spec, obs_vars, compiled,
                          me_variance = me_variance),
                     ppf_extra)))
  }

  if (likelihood == "pskf") {
    ## Pruned Skewed Kalman Filter. Deterministic likelihood (no particles);
    ## analytic gradient not supported (brief: .ctx_allows_analytic_gradient
    ## returns FALSE for pskf since likelihood != "gaussian").
    ## Merge any cut_tol from ... for exposed pruning parameter.
    pskf_dots <- list(...)
    cut_tol   <- if (!is.null(pskf_dots$cut_tol)) pskf_dots$cut_tol else 0.01
    max_q     <- if (!is.null(pskf_dots$max_q)) pskf_dots$max_q else 5L
    return(make_log_posterior_pskf(
      model, data, prior_spec, obs_vars,
      compiled,
      me_variance   = me_variance,
      system_priors = system_priors,
      cut_tol       = cut_tol,
      max_q         = max_q
    ))
  }

  if (likelihood == "student_t") {
    ## Student-t innovation likelihood: Gaussian KF recursions + multivariate-t
    ## per-period log-density. Analytic gradient not supported (no adjoint path).
    if (is.null(student_df))
      stop("make_log_posterior: likelihood = \"student_t\" requires student_df ",
           "(degrees of freedom, a positive scalar).", call. = FALSE)
    nu <- student_df
    return(local({
      nu_ <- nu; mv_ <- me_variance; li_ <- lik_init
      sys_cache_ <- cache_system_structure(compiled)
      ss_warm_ <- NULL
      function(theta) {
        lp <- log_prior(theta, prior_spec)
        if (!is.finite(lp))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
        params <- .apply_theta_to_params(model, theta)
        ss_result <- solve_steady_state(model, compiled, params,
                                        y0 = ss_warm_, verbose = FALSE)
        if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
          ss_warm_ <<- NULL
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
        }
        ss_warm_ <<- ss_result$ss
        params <- ss_result$params %||% params
        sys <- extract_system_matrices_fast(sys_cache_, ss_result$ss, params)
        dr  <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
        if (is.null(dr) || !isTRUE(dr$bk_satisfied))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
        kf <- tryCatch(
          kalman_filter_student_t(data, dr, model, params, obs_vars,
                                  student_df  = nu_,
                                  me_variance = mv_,
                                  lik_init    = li_),
          error = function(e) {
            ## Same masking issue as the gaussian path: a genuine bug in the
            ## Student-t filter is indistinguishable from an infeasible draw here.
            ## dynhr_set_options(debug_kf_errors = TRUE) RE-RAISES for debugging.
            if (isTRUE(.dynhr_opt("debug_kf_errors", default = FALSE))) stop(e)
            NULL
          }
        )
        if (is.null(kf) || !is.finite(kf$loglik))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
        loglik <- kf$loglik
        sp_lp  <- if (!is.null(system_priors)) {
          .eval_system_priors(
            system_priors,
            list(theta   = theta,
                 model   = model,
                 dr      = dr,
                 Sigma_e = .get_shock_cov(model, model$varexo_names, params),
                 params  = params))
        } else 0
        list(logpost = lp + power * loglik + sp_lp,
             loglik  = loglik,
             logprior = lp)
      }
    }))
  }

  if (likelihood == "cumulant") {
    return(make_log_posterior_cumulant(model, data, prior_spec, obs_vars,
                                        compiled, me_variance = me_variance,
                                        system_priors = system_priors,
                                        ...))
  }

  if (likelihood == "pruned") {
    ## Pruned-SS Gaussian KF on the AFVRR augmented state (order 2 by default,
    ## order 3 when pruned_order = 3L).  Analytic gradient NOT available;
    ## gradient-based optimizers use numerical FD.
    if (pruned_order == 3L)
      return(make_log_posterior_pruned3(model, data, prior_spec, obs_vars,
                                        compiled, me_variance = me_variance,
                                        system_priors = system_priors))
    return(make_log_posterior_pruned(model, data, prior_spec, obs_vars,
                                      compiled, me_variance = me_variance,
                                      system_priors = system_priors))
  }

  if (likelihood == "whittle") {
    extra <- list(...)
    debias_arg <- if (!is.null(extra$debias)) extra$debias else TRUE
    return(make_log_posterior_whittle(model, data, prior_spec, obs_vars,
                                       compiled, me_variance = me_variance,
                                       freq_band = freq_band,
                                       system_priors = system_priors,
                                       debias = debias_arg))
  }

  ## ---- Structural MS-DSGE branch ------------------------------------------
  ## When ms_struct_spec_ctx is present, solve per-regime decision rules at
  ## each draw (solve_ms_perturbation) then evaluate the structural Kim-Nelson
  ## filter (ms_kim_filter_struct). The single-regime steady-state solve and
  ## solve_perturbation are bypassed; the solver is the main per-draw cost.
  ##
  ## Per-regime parameter convention:
  ##   The spec carries params_by_regime (baseline param vectors). At each draw,
  ##   theta is applied to EACH regime's params via .apply_theta_to_params, so
  ##   estimated structural parameters affect all regimes identically (global
  ##   params) or selectively (regime params override at spec construction time).
  ##   If the spec was built with regime-specific overrides, the estimated theta
  ##   still applies on top of the baseline -- this is the correct convention for
  ##   "estimate a subset of params while holding regime differences fixed."
  ##
  ## Warm-start cache: a per-closure `ms_dr_warm` stores the last converged
  ##   MsDecisionRules. The solver uses it as a starting point when params
  ##   change only slightly across MCMC draws.  Per-closure (not global), so
  ##   each parallel chain keeps its own cache.
  if (!is.null(ms_struct_spec_ctx)) {
    if (is.null(compiled$lead_lag_incidence) &&
        !is.null(compiled$model$lead_lag_incidence))
      compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

    sys_cache_ms <- cache_system_structure(compiled)
    ms_ss_warm   <- NULL            # per-draw ss warm-start (regime 1 only; others follow)

    return(local({
      ## Capture spec fields once (avoid repeated $ lookups in the hot path)
      P_ms     <- ms_struct_spec_ctx$transition
      n_reg    <- ms_struct_spec_ctx$n_regimes
      pbr_base <- ms_struct_spec_ctx$params_by_regime  # baseline param lists
      sbr_pre  <- ms_struct_spec_ctx$ss_by_regime      # NULL = solve per draw
      mv_ms    <- me_variance
      li_ms    <- lik_init
      sp_ms    <- system_priors
      pw_ms    <- power

      ## Per-closure warm starts (per-chain state)
      ms_ss_warm_  <- NULL   # length-n_reg list of per-regime ss warm starts
      ms_dr_warm_  <- NULL   # last converged MsDecisionRules (unused currently;
                             ## solver always runs from QZ init -- warm-starting
                             ## the iteration is future work)

      function(theta) {
        lp <- log_prior(theta, prior_spec)
        if (!is.finite(lp))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

        ## Apply theta to each regime's params (structural params update all
        ## regimes; regime-distinguishing params were fixed at spec build time).
        params_list <- lapply(pbr_base, function(pr_s)
          .apply_theta_to_params(model, theta, pr_s))

        ## Per-regime steady states: use pre-supplied or solve per draw.
        if (!is.null(sbr_pre)) {
          ss_list <- sbr_pre
        } else {
          ss_list <- vector("list", n_reg)
          for (s in seq_len(n_reg)) {
            y0_s <- if (!is.null(ms_ss_warm_)) ms_ss_warm_[[s]] else NULL
            ss_s <- solve_steady_state(model, compiled, params_list[[s]],
                                       y0 = y0_s, verbose = FALSE)
            if (is.null(ss_s) || !isTRUE(ss_s$converged)) {
              ## Retry cold
              if (!is.null(y0_s))
                ss_s <- solve_steady_state(model, compiled, params_list[[s]],
                                           verbose = FALSE)
              if (is.null(ss_s) || !isTRUE(ss_s$converged)) {
                ms_ss_warm_ <<- NULL
                return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
              }
            }
            ## Sync any SSM-computed params
            params_list[[s]] <- ss_s$params %||% params_list[[s]]
            ss_list[[s]] <- ss_s$ss
          }
          ms_ss_warm_ <<- ss_list
        }

        ## Solve coupled MS decision rules
        ms_dr <- tryCatch(
          solve_ms_perturbation(
            model            = model,
            compiled         = compiled,
            ss_by_regime     = ss_list,
            params_by_regime = params_list,
            P                = P_ms,
            tol              = 1e-10,
            max_iter         = 500L,
            verbose          = FALSE
          ),
          error = function(e) NULL
        )
        if (is.null(ms_dr))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

        ## BK check: all regimes must satisfy BK
        if (!all(vapply(ms_dr$dr, function(d) isTRUE(d$bk_satisfied), logical(1))))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

        ## MSS check
        if (!isTRUE(ms_dr$mss_ok))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

        ## Kim-Nelson structural filter
        ## obs_vars and data use the same T x n_obs convention as the standard path.
        ## ms_kim_filter_struct expects data as n_obs x T.
        kf <- tryCatch(
          ms_kim_filter_struct(
            t(data), ms_dr, model,
            params     = params_list[[1L]],   # common params for obs equation
            obs_vars   = obs_vars,
            me_variance = mv_ms,
            lik_init    = li_ms
          ),
          error = function(e) {
            if (isTRUE(.dynhr_opt("debug_kf_errors", default = FALSE))) stop(e)
            NULL
          }
        )
        if (is.null(kf) || !is.finite(kf$loglik))
          return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

        ## System priors (use regime-1 dr as the representative)
        if (!is.null(sp_ms)) {
          sp_lp <- .eval_system_priors(
            sp_ms,
            list(theta   = theta,
                 model   = model,
                 dr      = ms_dr$dr[[1L]],
                 Sigma_e = .get_shock_cov(model, model$varexo_names, params_list[[1L]]),
                 params  = params_list[[1L]]))
          if (!is.finite(sp_lp))
            return(list(logpost = -Inf, loglik = kf$loglik, logprior = lp))
          lp <- lp + sp_lp
        }

        list(logpost = pw_ms * kf$loglik + lp, loglik = kf$loglik, logprior = lp)
      }
    }))
  }

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  sys_cache <- cache_system_structure(compiled)

  ## Warm-start cache for the steady-state solve. Consecutive MCMC proposals
  ## are close in parameter space, so the previously converged steady state is
  ## an excellent initial guess: Newton converges in a handful of iterations
  ## and the expensive optim fallback is never triggered. For models with a
  ## closed-form steady_state_model block this is a no-op (analytical path is
  ## hit first); for numerically-solved models it is the dominant per-draw
  ## saving (see inst/benchmarks/posterior.R). Per-closure state, so each
  ## parallel chain keeps its own warm start.
  ss_warm <- NULL
  ## me-floor hazard guard (see R/pruned-state-space.R and kalman_filter()'s
  ## me_floor_check arg): warn at most once per closure, not once per draw.
  .me_floor_checked <- FALSE

  function(theta) {
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    params <- .apply_theta_to_params(model, theta)

    ss_result <- solve_steady_state(model, compiled, params,
                                    y0 = ss_warm, verbose = FALSE)
    if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
      ## The warm guess may have been misleading (e.g. a large jump in the
      ## proposal). Retry once from the cold initval-based guess before
      ## declaring the draw infeasible.
      if (!is.null(ss_warm))
        ss_result <- solve_steady_state(model, compiled, params,
                                        verbose = FALSE)
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
        ss_warm <<- NULL
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
      }
    }
    ss_warm <<- ss_result$ss

    ## Re-derive any steady_state_model-computed parameter so the linearization
    ## point uses the consistent (not stale) p_c (no-op for non-SSM-parameter
    ## models). Without this the dynamic system is built at an invalid steady
    ## state (F != 0) with the wrong p_c, silently biasing the loglik (Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
    dr  <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied)) {
      if (is.null(.penalty_cfg)) {
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
      }
      ## Penalty mode: measure BK violation as the number of excess unstable
      ## roots beyond the Blanchard-Kahn rank condition.  The penalty is
      ## monotone in the violation count so the optimizer is pushed back toward
      ## feasibility.  loglik stays -Inf (correctness); only logpost carries the
      ## gradient-restoring penalty.
      bk_viol <- .bk_eigen_violation(dr, model)
      logpost_pen <- .penalty_cfg$floor -
        .penalty_cfg$scale * (1 + bk_viol)
      return(list(logpost = logpost_pen, loglik = -Inf, logprior = lp,
                  infeasible = TRUE, violation = bk_viol))
    }

    ## Stationarity guard. The stable generalized eigenvalues from the QZ
    ## solve ARE the eigenvalues of the state-transition block
    ## ghx[state_idx, ] (verified bit-equal to eigen(ghx_state) across the
    ## model corpus to ~1e-13), so reuse them instead of a fresh O(n^3)
    ## eigen() call (~100 us/draw on the feasible path). They are ordered
    ## stable-first and, given bk_satisfied, n_stable == n_state. This is NOT
    ## redundant with the BK check: the QZ select tolerance (1+1e-6) can admit
    ## a near-unit root that makes the stationary (Lyapunov-initialised) Kalman
    ## filter invalid, so we still reject |lambda| >= 1 here.
    ns <- length(dr$state_idx)
    ev <- dr$eigenvalues
    spectral_radius <- if (!is.null(ev) && length(ev) >= ns)
      max(Mod(ev[seq_len(ns)]))
    else
      max(Mod(eigen(dr$ghx[dr$state_idx, , drop = FALSE],
                    only.values = TRUE)$values))
    ## A unit root makes the stationary Lyapunov P0 invalid (NaN), but the
    ## exact-diffuse initialization handles it (Koopman-Durbin / DK 2012 ch.5).
    ## Mirror kalman_filter's lik_init resolution: "auto" switches to the
    ## diffuse init for a (near-)unit root, and "diffuse"/"kappa" are explicit
    ## diffuse requests -- in all three the draw is filterable, so reject ONLY
    ## when the init in force is the stationary Lyapunov one. The diffuse
    ## forward loglik is validated against the closed-form local-level
    ## likelihood to ~1e-12 (test-kalman-diffuse.R), so un-gating it produces a
    ## correct (not silently-wrong) posterior. Exception: a diffuse phase is
    ## unsupported with heteroskedastic shock_scale (kalman_filter stop()s), so
    ## unit-root draws are still rejected in that case.
    init_in_force <- if (identical(lik_init, "auto"))
      (if (spectral_radius > 1 - 1e-6) "diffuse" else "stationary")
    else lik_init
    if (spectral_radius >= 1 &&
        (identical(init_in_force, "stationary") || !is.null(shock_scale))) {
      if (is.null(.penalty_cfg)) {
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
      }
      ## Penalty mode: violation = sum of how far each stable eigenvalue exceeds
      ## 1 in modulus.  Continuous across the feasibility boundary (zero exactly
      ## at |lambda| = 1); monotone — larger eigenvalues give a larger penalty.
      viol <- if (is.finite(spectral_radius)) max(spectral_radius - 1, 0) else 1
      logpost_pen <- .penalty_cfg$floor -
        .penalty_cfg$scale * (1 + viol)
      return(list(logpost = logpost_pen, loglik = -Inf, logprior = lp,
                  infeasible = TRUE, violation = viol))
    }

    ## MS-DSGE path: use Kim-Nelson filter when ms_spec is present.
    ## The non-MS (default) path is bit-identical when ms_spec is NULL.
    if (!is.null(ms_spec_ctx)) {
      kf <- tryCatch(
        ms_kim_filter(t(data), dr, model, params, obs_vars,
                      ms_spec     = ms_spec_ctx,
                      me_variance = me_variance,
                      lik_init    = lik_init),
        error = function(e) {
          if (isTRUE(.dynhr_opt("debug_kf_errors", default = FALSE))) stop(e)
          NULL  # infeasible draw -> -Inf (see the non-MS branch's note)
        }
      )
    } else {
      kf <- tryCatch(
        kalman_filter(data, dr, model, params, obs_vars,
                      return_filtered = FALSE, me_variance = me_variance,
                      lik_init = lik_init, me_extra = me_extra,
                      shock_scale = shock_scale,
                      me_floor_check = !.me_floor_checked &&
                        isTRUE(getOption("dynhr.me_floor_check", TRUE))),
        error = function(e) {
          ## Lyapunov / inv_sympd / other KF failures must not propagate: they
          ## indicate an infeasible parameter draw (singular covariance, unit-root
          ## with stationary init, etc.). Return NULL -> -Inf rather than crashing
          ## the optimizer chain (RC1b fix). A genuine code bug is indistinguishable
          ## from an infeasible draw here and would be silently masked as a rejected
          ## draw; dynhr_set_options(debug_kf_errors = TRUE) RE-RAISES instead, so a
          ## masked bug surfaces during debugging.
          if (isTRUE(.dynhr_opt("debug_kf_errors", default = FALSE))) stop(e)
          NULL
        }
      )
      .me_floor_checked <<- TRUE   # guard once per closure, not per MCMC draw
    }
    if (is.null(kf) || !is.finite(kf$loglik))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## System priors: penalty terms on model features (IRF signs, variance
    ## shares, etc.). Evaluated once from the already-solved dr; no re-solve.
    if (!is.null(system_priors)) {
      sp_lp <- .eval_system_priors(
        system_priors,
        list(theta   = theta,
             model   = model,
             dr      = dr,
             Sigma_e = .get_shock_cov(model, model$varexo_names, params),
             params  = params))
      if (!is.finite(sp_lp))
        return(list(logpost = -Inf, loglik = kf$loglik, logprior = lp))
      lp <- lp + sp_lp
    }

    ## `power` tempers only the logpost; $loglik carries the raw likelihood so
    ## SMC tempering, marginal-likelihood estimators, and diagnostics see the
    ## true likelihood unchanged.
    list(logpost = power * kf$loglik + lp, loglik = kf$loglik, logprior = lp)
  }
}

#' Build a per-period log-likelihood-contribution closure
#'
#' \code{robust_confidence_set()} (Andrews-Mikusheva LM2,
#' \code{\link{robust_confidence_set}}) and \code{run_diagnostics(loglik_contrib_fn
#' = )} need a \code{theta -> } length-\eqn{T} per-period log-likelihood
#' closure. Hand-reconstructing \code{solve_perturbation -> kalman_filter}
#' outside \code{make_log_posterior} risks a likelihood that differs
#' from the posterior's (a config mismatch in \code{me_variance} handling,
#' observable steady-state constants, or \code{lik_init} routing). This
#' builder reuses \code{make_log_posterior}'s EXACT steady-state-solve,
#' Blanchard-Kahn / stationarity feasibility guards, and
#' \code{kalman_filter()} call (via the shared internal
#' \code{.solve_dr_for_theta()} helper) so that
#' \code{sum(make_loglik_contrib(...)(theta))} equals the \code{$loglik}
#' that \code{make_log_posterior}/\code{\link{make_posterior}} would
#' return for the same \code{theta}, to floating-point precision.
#'
#' Only the Gaussian (Kalman filter) likelihood exposes per-period
#' contributions; other likelihoods (cumulant, whittle, pskf, tpf, pruned,
#' ...) do not exist as a single per-period decomposition in the same sense
#' and this builder errors clearly rather than approximate one.
#'
#' @param model       dynhr_mod (from \code{\link{parse_mod}})
#' @param data        Observation matrix (\eqn{T \times n_{\text{obs}}}),
#'   columns matching \code{obs_vars}
#' @param prior_spec  Prior specification data.frame (from
#'   \code{extract_prior_spec}). Unused for the returned
#'   likelihood-only contributions (no prior term is added), but
#'   \code{theta} must still supply every name the likelihood needs (the
#'   builder validates unknown prior targets the same way
#'   \code{make_log_posterior} does, in case the caller shares a
#'   \code{prior_spec} with \code{theta}'s naming). Pass \code{NULL} to skip
#'   this cross-check entirely.
#' @param obs_vars    Character vector of observed variable names. Defaults
#'   to \code{model$obs_vars} (populated from a \code{varobs} declaration)
#'   when \code{NULL}, mirroring \code{make_log_posterior}'s default.
#' @param compiled    dynhr_compiled (from \code{\link{compile_model}})
#' @param me_variance Measurement error variance (default 0), identical
#'   convention to \code{make_log_posterior} / \code{\link{kalman_filter}}.
#' @param lik_init    Kalman filter \code{P0} initialization, forwarded to
#'   \code{\link{kalman_filter}} (default \code{"auto"}); see
#'   \code{make_log_posterior}.
#' @param me_extra    \code{n_obs x T} matrix of additional per-observable,
#'   per-period measurement-error variance (filter_tunes soft tunes); see
#'   \code{\link{kalman_filter}}.
#' @param shock_scale \code{n_exo x T} heteroskedastic shock-scale matrix;
#'   see \code{\link{kalman_filter}}.
#' @param likelihood  Likelihood type. Only \code{"gaussian"} (the default)
#'   is supported; any other value errors ("gaussian only for now") rather
#'   than silently approximating a per-period decomposition it does not have.
#' @param ...         Accepted for signature parity with
#'   \code{make_log_posterior}; unused (errors are raised instead of
#'   silently ignoring arguments that would change the likelihood, except
#'   arguments consumed by the gaussian path already listed above).
#' @return A function \code{function(theta)} returning a numeric vector of
#'   length \code{nrow(data)}: the per-period Gaussian log-likelihood
#'   contribution (prediction-error decomposition), \strong{likelihood
#'   only} (no prior term). Returns a length-\code{nrow(data)} vector of
#'   \code{-Inf} for an out-of-domain \code{theta} (steady-state solve
#'   failure, Blanchard-Kahn violation, or a unit root under a stationary
#'   \code{lik_init}) rather than erroring, matching
#'   \code{\link{robust_confidence_set}}'s contract.
#' @seealso \code{make_log_posterior}, \code{\link{make_posterior}},
#'   \code{\link{robust_confidence_set}}
#' @export
make_loglik_contrib <- function(model, data, prior_spec = NULL, obs_vars = NULL,
                                compiled, me_variance = 0,
                                lik_init = "auto",
                                me_extra = NULL,
                                shock_scale = NULL,
                                likelihood = "gaussian",
                                ...) {
  if (!identical(likelihood, "gaussian"))
    stop("make_loglik_contrib: likelihood = \"", likelihood, "\" is not ",
         "supported -- gaussian only for now. Per-period log-likelihood ",
         "contributions are only exposed by the Kalman-filter (gaussian) ",
         "path; the other likelihoods (cumulant, whittle, pskf, tpf, ",
         "pruned, student_t, ppf/copf, MS-DSGE) do not have a per-period ",
         "decomposition wired up and this builder will not silently ",
         "approximate one.", call. = FALSE)

  ## Default obs_vars from the model's varobs declaration, identical to
  ## make_log_posterior (pathological-DSGE paper gap #4 / v9024).
  if (is.null(obs_vars) || length(obs_vars) == 0L) {
    obs_vars <- model$obs_vars %||% model$varobs_names
    if (is.null(obs_vars) || length(obs_vars) == 0L)
      stop("make_loglik_contrib: `obs_vars` not supplied and the model ",
           "declares no `varobs`; pass obs_vars= or add a varobs line to ",
           "the .mod.", call. = FALSE)
  }

  ## Same fail-loud guard as make_log_posterior (only when a prior_spec is
  ## actually supplied -- the contributions are likelihood-only, so a caller
  ## may reasonably pass prior_spec = NULL and skip this cross-check).
  if (!is.null(prior_spec))
    .validate_prior_targets(model, prior_spec, where = "make_loglik_contrib")

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  sys_cache <- cache_system_structure(compiled)
  n_T       <- nrow(data)

  ## Per-closure mutable state: warm-started steady state (mirrors
  ## make_log_posterior's ss_warm cache) and the me-floor-check once-guard.
  state             <- new.env(parent = emptyenv())
  state$ss_warm     <- NULL
  .me_floor_checked <- FALSE

  function(theta) {
    solved <- .solve_dr_for_theta(model, compiled, sys_cache, theta, state,
                                  lik_init = lik_init, shock_scale = shock_scale)
    if (is.null(solved)) return(rep(-Inf, n_T))

    kf <- tryCatch(
      kalman_filter(data, solved$dr, model, solved$params, obs_vars,
                    return_filtered = FALSE, me_variance = me_variance,
                    return_ll_contrib = TRUE,
                    lik_init = lik_init, me_extra = me_extra,
                    shock_scale = shock_scale,
                    me_floor_check = !.me_floor_checked &&
                      isTRUE(getOption("dynhr.me_floor_check", TRUE))),
      error = function(e) {
        if (isTRUE(.dynhr_opt("debug_kf_errors", default = FALSE))) stop(e)
        NULL
      }
    )
    .me_floor_checked <<- TRUE

    if (is.null(kf) || is.null(kf$loglik_contrib) ||
        length(kf$loglik_contrib) != n_T || !is.finite(kf$loglik))
      return(rep(-Inf, n_T))

    kf$loglik_contrib
  }
}

#' Print posterior summary table
#' @noRd
print_posterior <- function(mcmc, prior_spec) {
  chain <- mcmc$chain
  cat("\nPOSTERIOR SUMMARY\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("%-15s %10s %10s %10s %10s %10s\n",
              "Parameter", "Prior Mean", "Post Mean", "Post Median", "5%", "95%"))
  cat(strrep("-", 80), "\n")
  for (i in seq_len(ncol(chain))) {
    nm      <- colnames(chain)[i]
    pr_mean <- prior_spec$mean[prior_spec$name == nm]
    if (length(pr_mean) == 0) pr_mean <- NA
    x <- chain[, i]
    cat(sprintf("%-15s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                nm, pr_mean, mean(x), median(x),
                quantile(x, 0.05), quantile(x, 0.95)))
  }
  cat(strrep("-", 80), "\n")
}

#' @noRd
print.EstimationResult <- function(x, ...) {
  cat("=== DSGE Estimation Result ===\n")
  cat("Observed vars:   ", paste(x$obs_vars, collapse = ", "), "\n")
  cat("Estimated params:", ncol(x$mcmc$chain), "\n")
  cat("MCMC draws:      ", x$mcmc$n_draws, "(burn-in:", x$mcmc$n_burn, ")\n")
  cat("Acceptance rate: ", sprintf("%.1f%%", x$mcmc$acceptance_rate * 100), "\n")
  cat("==============================\n")
  invisible(x)
}
