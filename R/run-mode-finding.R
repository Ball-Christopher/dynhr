## R/run-mode-finding.R
## --------------------------------------------------------------------------
## run_mode_finding() -- Unified entry point for posterior mode-finding.
##
## Given a solved model (dynhr_solved), data, and observable variable names,
## this function:
##   1. Extracts the prior specification (auto eps_* -> sig_* rename)
##   2. Builds the log-posterior closure
##   3. Runs multi-stage mode-finding
##   4. Computes the Hessian-based proposal covariance
##   5. Returns a comprehensive result object ready for estimation
## --------------------------------------------------------------------------


## The `...` names of run_mode_finding() that go to the OPTIMISER (the formals
## of .run_mode_finding() after its three positional arguments); every other
## `...` name goes to the log-posterior constructor. One definition for the
## runner's split and for the spec's check of mode$extra.
.rmf_optimiser_arg_names <- function()
  setdiff(names(formals(.run_mode_finding)),
          c("log_post_fn", "theta_init", "prior_spec"))

## The keys the mode stage reads from `mode_options` (mode$options of a spec).
## `verbose` is accepted and ignored (the spec's compute$verbose governs).
.mode_option_keys <- c("use_analytic_grad", "use_analytic_hess", "h0_method",
                       "record_curvature", "verbose")

#' Validate and clip a caller-supplied \code{theta_init} against a prior spec
#'
#' Matching is BY NAME against \code{priors$name} (never by position: the prior
#' spec's row order is the \code{estimated_params} order after the \code{eps_* -> sig_*}
#' rename, which a caller has no reason to reproduce). The returned vector is
#' reordered to the spec's own order and carries the spec's names, so it is a
#' drop-in replacement for \code{setNames(priors$mean, priors$name)}.
#'
#' Values are clipped strictly INSIDE \code{[lower, upper]}: a start sitting exactly
#' on a bound makes the eta-transform infinite and makes the bounded optimiser
#' stages start on the feasibility boundary. The nudge is \code{1e-8 * (upper -
#' lower)} for a finite support, and \code{1e-8 * max(1, |bound|)} when the other
#' side is infinite.
#'
#' @param theta_init Named numeric vector from the caller.
#' @param priors     Prior spec data.frame from \code{extract_prior_spec()}.
#' @return Named numeric vector, length/order/names matching \code{priors$name}.
#' @noRd
.rmf_check_theta_init <- function(theta_init, priors) {
  want <- as.character(priors$name)

  if (!is.numeric(theta_init) || !is.null(dim(theta_init)))
    .dynhr_abort("run_mode_finding: `theta_init` must be a named numeric ",
                 "vector (got ", class(theta_init)[1L], ").")
  nm <- names(theta_init)
  if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm)))
    .dynhr_abort("run_mode_finding: `theta_init` must be NAMED; its names are ",
                 "matched against the prior spec's parameter names (",
                 paste(want, collapse = ", "), ").")
  if (anyDuplicated(nm))
    .dynhr_abort("run_mode_finding: `theta_init` has duplicated name(s): ",
                 paste(unique(nm[duplicated(nm)]), collapse = ", "), ".")

  missing_nm <- setdiff(want, nm)
  extra_nm   <- setdiff(nm, want)
  if (length(missing_nm) || length(extra_nm))
    .dynhr_abort(
      "run_mode_finding: `theta_init` names do not match the estimated ",
      "parameters. ",
      if (length(missing_nm))
        paste0("Missing: ", paste(missing_nm, collapse = ", "), ". ") else "",
      if (length(extra_nm))
        paste0("Unknown: ", paste(extra_nm, collapse = ", "), ". ") else "",
      "Expected exactly: ", paste(want, collapse = ", "), ".")

  th <- as.numeric(theta_init[want])
  bad <- !is.finite(th)
  if (any(bad))
    .dynhr_abort("run_mode_finding: `theta_init` has non-finite value(s) for: ",
                 paste(want[bad], collapse = ", "),
                 ". NA/NaN/Inf starts are not allowed.")

  lo <- suppressWarnings(as.numeric(priors$lower))
  hi <- suppressWarnings(as.numeric(priors$upper))
  if (length(lo) != length(want)) lo <- rep(-Inf, length(want))
  if (length(hi) != length(want)) hi <- rep( Inf, length(want))
  lo[is.na(lo)] <- -Inf
  hi[is.na(hi)] <-  Inf

  for (i in seq_along(th)) {
    l <- lo[i]; u <- hi[i]
    if (!(is.finite(l) || is.finite(u))) next
    eps <- if (is.finite(l) && is.finite(u)) 1e-8 * (u - l) else
      1e-8 * max(1, abs(if (is.finite(l)) l else u))
    if (is.finite(l) && th[i] <= l + eps) th[i] <- l + eps
    if (is.finite(u) && th[i] >= u - eps) th[i] <- u - eps
  }

  setNames(th, want)
}

#' Default mode-finding start: INITVAL / estimated_params_init, else prior mean
#'
#' Dynare starts the optimiser at the \code{estimated_params}
#' INITVAL (overridden by an \code{estimated_params_init} block; \code{set_prior.m} falls
#' back to the prior mean only where INITVAL is NaN). \code{extract_prior_spec()}
#' records that value in \code{priors$init}; before this fix it was read nowhere and
#' every run started at the prior means.
#'
#' Per parameter: \code{init} when finite and inside \code{[lower, upper]} (nudged off a
#' bound exactly as \code{.rmf_check_theta_init()} does), otherwise the prior mean.
#' A finite \code{init} OUTSIDE the bounds (Dynare's check_prior_bounds() would
#' error) falls back to the prior mean with a classed warning.
#'
#' @param priors Prior spec from \code{extract_prior_spec()}.
#' @return list(theta = named numeric in spec order, source = label).
#' @noRd
.rmf_default_theta_init <- function(priors) {
  nm    <- as.character(priors$name)
  start <- setNames(as.numeric(priors$mean), nm)
  init  <- suppressWarnings(as.numeric(priors$init %||% rep(NA_real_, length(nm))))
  if (length(init) != length(nm) || !any(is.finite(init)))
    return(list(theta = start, source = "prior means"))

  lo <- suppressWarnings(as.numeric(priors$lower)); lo[is.na(lo)] <- -Inf
  hi <- suppressWarnings(as.numeric(priors$upper)); hi[is.na(hi)] <-  Inf
  fin    <- is.finite(init)
  inside <- fin & init >= lo & init <= hi
  if (any(fin & !inside))
    .dynhr_warn(
      "run_mode_finding: INITVAL outside the prior bounds for ",
      paste0(nm[fin & !inside], " (", init[fin & !inside], " not in [",
             lo[fin & !inside], ", ", hi[fin & !inside], "])", collapse = "; "),
      "; starting those parameter(s) at the prior mean instead.",
      class = "dynhr_warning_init_outside_bounds")
  if (!any(inside))
    return(list(theta = start, source = "prior means"))

  clipped <- .rmf_check_theta_init(setNames(init[inside], nm[inside]),
                                   priors[inside, , drop = FALSE])
  start[inside] <- clipped
  src <- if (all(inside)) "INITVAL (estimated_params / estimated_params_init)"
         else sprintf("INITVAL for %d of %d parameters, prior means otherwise",
                      sum(inside), length(nm))
  list(theta = start, source = src)
}


## The Kalman P0 in force at `theta` under `lik_init` (.grad_init_in_force:
## "stationary", "diffuse", "kappa" or "reject"), for choosing the newrat H0
## seed: posterior_hessian() -- the analytic seed -- refuses every init but
## "stationary". NA when the model does not solve at theta (the seed builder
## then fails and falls back on its own, as before).
.rmf_h0_init_in_force <- function(model, compiled, theta, lik_init) {
  sol <- .grad_solve_dr(model, compiled, cache_system_structure(compiled), theta)
  if (is.null(sol)) return(NA_character_)
  dr <- sol$dr
  .grad_init_in_force(dr$ghx[dr$state_idx, , drop = FALSE],
                      dr$ghu[dr$state_idx, , drop = FALSE],
                      .get_shock_cov(model, model$varexo_names, sol$params),
                      lik_init, dr = dr)
}


#' Run posterior mode-finding for a solved DSGE model
#'
#' Takes a solved model (from \code{\link{solve_model}}), observation data,
#' and observable variable names, and runs the full mode-finding pipeline:
#' prior extraction -> log-posterior construction -> multi-stage optimisation
#' -> Hessian-based proposal covariance.
#'
#' @param solved      A \code{dynhr_solved} object from \code{\link{solve_model}}.
#' @param data        Observation matrix (\eqn{T \times n_{obs}}), column names
#'   matching \code{obs_vars}.  May also be a path to a CSV file.
#' @param obs_vars    Character vector of observable variable names.
#' @param n_iter      Total optimizer iteration budget (default 10000).
#' @param method      Optimizer sequence (default \code{"newrat"}: the
#'   csminwel quasi-Newton optimizer, equivalent to Dynare's
#'   \code{mode_compute = 4}). Newrat uses the analytic gradient (built
#'   automatically for Gaussian/cumulant/Whittle likelihoods) and an
#'   infeasibility-backtracking line search, which is markedly more robust and
#'   far cheaper on near-unit-root models than the gradient-free CMA-ES global
#'   search: it needs hundreds of evaluations where CMA-ES needs thousands and
#'   tends to dwell on the feasibility boundary. For a multimodal posterior or
#'   a poor starting point, use \code{"cmaes_newrat"} (CMA-ES global search to
#'   escape poor starts, then newrat polish) or \code{"combined"} (the former
#'   default: CMA-ES then L-BFGS-B). See \code{\link{find_mode}} for the full
#'   list of options.
#' @param me_variance Measurement error variance added to the observation
#'   noise diagonal. \code{NULL} (default) reads the \code{me_variance}
#'   option (0 unless set), as the estimation spec's
#'   \code{likelihood$me_variance} does.  Use a small positive value for
#'   stochastically singular models.
#' @param likelihood  Likelihood type: \code{"gaussian"} (Kalman filter,
#'   default), \code{"cumulant"} (cumulant-matching, Mutschler 2015),
#'   \code{"whittle"} (frequency-domain Whittle likelihood),
#'   \code{"pruned"} (Gaussian KF on the AFVRR pruned state space; see
#'   \code{pruned_order}), \code{"pskf"} (Pruned Skewed Kalman Filter for
#'   skew-normal shocks; pass \code{cut_tol} via \code{...}), or
#'   \code{"student_t"} (Gaussian-KF recursions with a multivariate Student-t
#'   per-period density; requires \code{student_df} passed via \code{...}).
#'   The Whittle path requires a stationary, complete panel; use
#'   \code{freq_band} (via \code{posterior_options}) for band-restricted
#'   estimation. The pruned, pskf, and student_t paths are all deterministic
#'   (no stochastic-particle noise) but have no analytic gradient (FD
#'   fallback); pruned additionally runs the multi-start step via the
#'   closure-shipped parallel path. \code{"tpf"}/\code{"ppf"}/\code{"copf"}
#'   (particle-filter likelihoods with a noisy unbiased loglik estimate) are
#'   deliberately NOT accepted by this function -- their estimation noise
#'   breaks deterministic optimizers (Nelder-Mead/CMA-ES/newrat all assume a
#'   fixed objective at repeated evaluations of the same point); use
#'   \code{\link{run_full_estimation}} or particle MCMC (PMMH) via
#'   \code{\link{run_posterior_estimation}} instead. (An estimation spec run
#'   by \code{\link{run_estimation}} does accept them and, when its mode
#'   stage runs, warns that the optimiser works on a noisy objective.)
#' @param pruned_order Integer, \code{2L} (default) or \code{3L}: AFVRR
#'   pruned state-space order used when \code{likelihood = "pruned"}
#'   (ignored otherwise; \code{3L} with another likelihood is an error).
#' @param data_col_map Optional named character vector mapping model observable
#'   names to CSV column names when they differ.
#' @param mode_options  List of additional options for mode-finding,
#'   e.g. \code{list(verbose = TRUE)}.
#' @param posterior_options  List of additional options forwarded to the
#'   internal log-posterior constructor (e.g. \code{list(order = 2)} for
#'   cumulant likelihood).
#' @param parallel  Run mode-finding (Step 5) as a multi-start search and the
#'   Hessian-based proposal covariance (Step 6) on a mirai daemon pool
#'   (default \code{FALSE}). Step 5 dispatches \code{n_starts} dispersed
#'   starts via \code{run_mode_mirai()} and keeps the best; Step 6 ships the
#'   log-posterior closure to the daemons as-is. Both apply to ANY likelihood
#'   (Gaussian, OBC/PKF, cumulant).
#'
#'   \strong{Performance:} this is the single biggest lever for the Step-6
#'   covariance Hessian. The default \code{FALSE} path computes the
#'   finite-difference Hessian serially (\code{num_hessian}: \eqn{O(n_{par}^2)}
#'   likelihood evaluations); \code{parallel = TRUE} distributes those across
#'   the daemon pool (\code{num_hessian_mirai}, roughly an \eqn{n}-fold
#'   speed-up on \eqn{n} cores). It is left opt-in by default rather than
#'   auto-enabled because spinning up and tearing down a mirai daemon pool has
#'   a fixed setup cost and a lifecycle to manage, and is undesirable in CI or
#'   embedded contexts. For a large parameter vector (where the Hessian
#'   dominates wall time) prefer \code{parallel = TRUE}.
#' @param n_cores   Daemon count when \code{parallel = TRUE} (\code{NULL} =
#'   auto).
#' @param n_starts  Number of dispersed starting points for the parallel
#'   multi-start mode search in Step 5 (\code{NULL} = one per daemon).
#'   Ignored when \code{parallel = FALSE}.
#' @param transform_params  Run Step 5 mode-finding in unconstrained eta-space
#'   (default \code{TRUE}): via
#'   \code{build_param_transform} (built from \code{priors}) and
#'   \code{.run_mode_finding}'s \code{transform} argument -- see there for
#'   the exact semantics (Jacobian-free objective, eta-space bounds for
#'   bounded optimizer stages, \code{theta_mode} returned in theta-space as
#'   always). Applies to both the serial path and the PARALLEL multi-start
#'   path (\code{run_mode_mirai}'s \code{transform} argument): each
#'   daemon's \code{.run_mode_finding} call receives the transform, and
#'   dispersed multi-start points are jittered in eta-space (no bound
#'   clamping) before being mapped back to theta-space starts. May also be
#'   set globally via \code{dynhr_set_options(transform_params = TRUE)}. Has
#'   no effect on Step 6 (proposal covariance), which is always computed in
#'   theta-space here; \code{run_posterior_estimation()} performs its own
#'   delta-method conversion when sampling in eta-space. Exception: at a
#'   mode on a prior bound (a constrained, KKT mode, where that conversion
#'   divides by \eqn{(d\theta/d\eta)^2 \approx} (distance to the bound)^2),
#'   Step 6 also returns \code{Sigma_prop_eta}, see Value.
#' @param proposal_cov_method  Step 6 proposal-covariance strategy:
#'   \code{"diagonal"} (default) keeps the existing
#'   \code{build_sigma_prop} behaviour (full-Hessian-with-cap,
#'   empirical, or prior-std fallback). \code{"full"} instead uses
#'   \code{proposal_cov(method = "full")}: the inverse of the full
#'   negative Hessian at the mode with eigenvalue repair (Dynare /
#'   Herbst & Schorfheide 2015 style), capturing posterior correlations
#'   for RWMH. Falls back to the diagonal univariate-curvature proposal
#'   with a warning if the full Hessian is unusable. May also be set via
#'   \code{dynhr_set_options(proposal_cov_method = "full")}.
#' @param filter_tunes Call-level tune override.  \code{NULL} (default) uses
#'   the \code{filter_tunes} block from the \code{.mod} file.  Pass a
#'   \code{filter_tunes_spec} object (from \code{\link{filter_tunes}()})
#'   to override the mod-file block, or \code{FALSE} to disable tunes.
#' @param heteroskedastic_shocks  Call-level heteroskedastic-shocks override.
#'   \code{NULL} (default) uses the \code{heteroskedastic_shocks} block from
#'   the \code{.mod} file if present.  Mutually exclusive with \code{plan}.
#' @param stochastic_volatility  Call-level stochastic-volatility override.
#'   \code{NULL} (default) uses the \code{stochastic_volatility} block from
#'   the \code{.mod} file if present; supply one to attach or replace the
#'   SV specification for this call only, leaving the model object untouched.
#' @param plan  A \code{\link{dynhr_plan}} object bundling
#'   \code{filter_tunes} and \code{heteroskedastic_shocks} overrides.
#'   Mutually exclusive with \code{filter_tunes} and
#'   \code{heteroskedastic_shocks}.
#' @param use_exact_hessian  Opt-in (default \code{NULL}: the
#'   \code{use_exact_hessian} option, else \code{FALSE}): build the RWMH
#'   proposal covariance from the analytic exact posterior Hessian
#'   (\code{posterior_hessian}) computed at the mode, instead of a
#'   finite-difference Hessian. Only applied for the standard Gaussian KF
#'   likelihood with no OBC, \code{me_extra}, heteroskedastic shocks, diffuse
#'   initialisation, or missing data; falls back to the numerical Hessian
#'   otherwise. The Hessian is cached in \code{result$hessian_exact}. Can also
#'   be set globally via \code{dynhr_set_options(use_exact_hessian = TRUE)}.
#'
#'   \strong{Why the default is \code{FALSE}.} The analytic Hessian is not
#'   uniformly faster than finite differences -- there is a crossover with
#'   model size. \code{posterior_hessian} forms second-order solution
#'   derivatives at a cost of \eqn{O(n_{par}^2 \cdot n_{state}^3)}, whereas
#'   the finite-difference Hessian is \eqn{O(n_{par}^2)} \emph{cheap} forward
#'   filter evaluations. For a small parameter vector and/or a long sample the
#'   analytic route wins; for a large state and many parameters (e.g.
#'   \eqn{n_{par} \gtrsim 60}, \eqn{n_{state} \gtrsim 70}) the finite-
#'   difference Hessian is the faster choice, so it is the default. To speed up
#'   the (default) numerical path use \code{parallel = TRUE}; do not reach for
#'   \code{use_exact_hessian = TRUE} expecting a speed-up on a large model.
#'   (Separately, the option \code{use_analytic_hess}, default \code{TRUE},
#'   uses the same analytic Hessian only to \emph{seed} the newrat optimiser's
#'   initial curvature \code{H0}; it does not set the reported covariance. On a
#'   large model that seed is itself expensive and of uncertain benefit -- set
#'   \code{mode_options = list(use_analytic_hess = FALSE)} (one run) or
#'   \code{dynhr_set_options(use_analytic_hess = FALSE)} (the session; the
#'   spec field \code{mode$analytic_h0}) to fall back to the cheap default
#'   \code{H0} if mode-finding setup is the bottleneck.)
#' @param verbose  Print progress messages (default TRUE).
#' @param theta_init  Optional named numeric vector of starting values for
#'   Step 4, in place of the default start. Names are matched against the prior
#'   spec's \code{name} column (the canonical keys \code{extract_prior_spec()}
#'   produces -- \code{sig_*} after the \code{eps_*} rename, \code{corr a,b},
#'   \code{skew x}); a vector with missing or unknown names, or carrying any
#'   \code{NA}/non-finite entry, is an error. Each value is clipped strictly
#'   inside the prior spec's \code{[lower, upper]} (nudged by
#'   \eqn{10^{-8}} times the support width when it sits on a bound), so the
#'   bounded optimiser stages and the eta-transform stay finite. The supplied
#'   start is used by EVERY path that reads the start: the serial optimiser,
#'   the parallel multi-start portfolio (where chain 1 is always the start
#'   itself), and the Hessian evaluated at the start. Default \code{NULL}
#'   starts each parameter at its \code{estimated_params} INITVAL (as
#'   overridden by an \code{estimated_params_init} block; the \code{init}
#'   column of \code{extract_prior_spec()}) when that is finite and inside
#'   \code{[lower, upper]}, and at its prior mean otherwise -- Dynare's rule.
#'
#'   The motivating use is \code{\link{diag_prior_sensitivity}}, which asks a
#'   LOCAL question -- "does the mode move when the prior is flattened?" --
#'   and so starts its flat-prior search at the informative mode rather than
#'   at the midpoint of each flat uniform support.
#' @param ...      Additional arguments forwarded to the mode-finder.
#'
#' @return An object of class \code{"dynhr_mode_result"} containing:
#'   \describe{
#'     \item{\code{solved}}{The input \code{dynhr_solved} object}
#'     \item{\code{data}}{Observation matrix}
#'     \item{\code{obs_vars}}{Observable variable names}
#'     \item{\code{prior_spec}}{Prior specification data.frame}
#'     \item{\code{log_post_fn}}{Log-posterior closure}
#'     \item{\code{theta_init}}{Starting parameter vector (INITVAL / prior
#'       means, or the clipped \code{theta_init=} argument when one was
#'       supplied)}
#'     \item{\code{theta_mode}}{Posterior mode parameter vector}
#'     \item{\code{mode}}{Full mode-finding result list}
#'     \item{\code{Sigma_prop}}{Proposal covariance for MCMC (theta-space).
#'       When the mode sits within the finite-difference stencil of a prior
#'       bound (a constrained mode: the posterior still rises toward the
#'       bound), those parameters are differentiated one-sided (inward) and
#'       decoupled from the rest: the other parameters get the inverse
#'       curvature of the posterior with the bound parameters held at the
#'       mode, and each bound parameter the inverse of its one-sided inward
#'       curvature -- a local step scale, not a posterior variance (its
#'       marginal is truncated and gradient-dominated; proposals across the
#'       bound are rejected by the prior). A non-concave inward direction
#'       falls back to the prior variance.}
#'     \item{\code{Sigma_prop_eta}}{Present only at such a bound mode: the
#'       eta-space (unconstrained, \code{transform_params}) proposal
#'       covariance. Interior parameters: the usual delta-method conversion
#'       of \code{Sigma_prop}; each bound parameter: the inverse negative
#'       curvature of the Jacobian-adjusted eta log-posterior at its own
#'       conditional mode (an O(1) variance), decoupled.}
#'     \item{\code{V_mode}}{Estimated variance-covariance at mode}
#'     \item{\code{me_variance},\code{likelihood}}{Settings used}
#'     \item{\code{run_record}}{A \code{dynhr_run_record}: resolved arguments,
#'       option snapshot, RNG state and provenance; replay it with
#'       \code{\link{dynhr_rerun}}.}
#'   }
#'
#' @examples
#' ## Model and data both ship with the package
#' solved   <- solve_model(system.file("extdata/models/nk_demo.mod",
#'                                     package = "dynhr"), verbose = FALSE)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#'
#' ## n_iter is capped here to keep the example quick
#' mode <- run_mode_finding(solved, Y, obs_vars = obs_vars,
#'                          n_iter = 200L, verbose = FALSE)
#' mode
#' round(mode$theta_mode, 4)
#'
#' ## The proposal covariance feeds straight into dynhr_mcmc()
#' dim(mode$Sigma_prop)
#'
#' @seealso \code{\link{solve_model}}, \code{\link{run_posterior_estimation}},
#'   \code{\link{find_mode}}, \code{\link{make_posterior}}
#' @export
run_mode_finding <- function(solved,
                             data,
                             obs_vars = NULL,
                             n_iter           = 10000L,
                             method           = "newrat",
                             me_variance      = NULL,
                             likelihood       = c("gaussian", "cumulant", "whittle",
                                                  "pruned", "pskf", "student_t",
                                                  "sv_rbpf"),
                             pruned_order     = 2L,
                             data_col_map     = NULL,
                             mode_options     = list(),
                             posterior_options = list(),
                             parallel         = FALSE,
                             n_cores          = NULL,
                             n_starts         = NULL,
                             proposal_cov_method = NULL,
                             transform_params = NULL,
                             filter_tunes     = NULL,
                             heteroskedastic_shocks = NULL,
                             stochastic_volatility = NULL,
                             plan             = NULL,
                             use_exact_hessian = NULL,
                             verbose          = TRUE,
                             theta_init       = NULL,
                             ...) {
  ## Own the message epoch for this run: repeat-suppressed warnings
  ## (`.dynhr_warn(once = TRUE)`) are keyed within it and re-arm for the
  ## next run, and the close reports what it suppressed. A nested call
  ## inherits this epoch rather than opening a second one.
  .dynhr_run_epoch <- .dynhr_epoch("run_mode_finding")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)
  ## Run record (R/run-record.R): resolved args, option snapshot and RNG
  ## state at ENTRY -- before the body touches any argument or the RNG.
  .rr <- .dynhr_rr_begin("run_mode_finding", environment(), list(...))
  ## A thin wrapper. The arguments become an estimation spec (the
  ## single home of the cross-field checks, validate_spec()) and the one spec
  ## runner does the work (mode stage: .est_mode_stage()).
  spec <- as_estimation_spec(.rr$args, entry = "run_mode_finding")
  .run_estimation_impl(spec, rr = .rr)
}

## Log-posterior-constructor extras (likelihood$extra / mode$extra) that
## make_posterior_grad() takes too: the mode stage forwards them to both.
.mode_grad_extras <- c("cumulant_orders", "cumulant_weight", "weight_matrix",
                       "debias")
## Extras that change the objective of `likelihood` but that the analytic
## gradient cannot take: with one of them the mode stage keeps the numerical
## gradient rather than differentiate another posterior. (The per-daemon pool
## gradient, .mirai_bind_worker_grad, is built only for the standard Gaussian
## posterior -- see grad_on_pool -- so it never sees these extras.)
.mode_grad_blocking_extras <- function(likelihood)
  switch(likelihood,
         cumulant = c("order", "h"),
         character(0))

## The mode stage of the spec runner (formerly the body of run_mode_finding()).
## `spec`: a validated dynhr_estimation_spec; `inp`: the prepared inputs
## (.est_inputs(): model, compiled, data, solved). Returns list(result = the
## dynhr_mode_result without its run record, plus the resolved pieces the
## later stages of a full run need: model, compiled, priors, data, obs_vars,
## me_extra, shock_scale, use_obc, obc_specs). `proposal = FALSE` skips Step 6
## (the proposal covariance: Sigma_prop / V_mode / hessian_exact stay NULL)
## for a full run that samples nothing. `optimise = FALSE` (prior-initialised
## samplers, see .spec_mode_stage_needed()) builds the posterior and its
## context but skips Steps 5-6: the result has no mode (theta_mode, mode,
## Sigma_prop NULL; meta$mode_skipped TRUE).
.est_mode_stage <- function(spec, inp, proposal = TRUE, optimise = TRUE) {
  lik <- spec$likelihood
  md  <- spec$mode
  cmp <- spec$compute
  n_iter              <- md$n_iter
  method              <- md$method
  me_variance         <- lik$me_variance
  likelihood          <- lik$type
  pruned_order        <- lik$pruned_order
  mode_options        <- md$options
  parallel            <- cmp$parallel
  n_cores             <- cmp$n_cores
  n_starts            <- md$n_starts
  proposal_cov_method <- md$proposal_cov
  transform_params    <- isTRUE(md$transform_params)
  use_exact_hessian   <- isTRUE(md$exact_hessian)
  verbose             <- cmp$verbose
  theta_init          <- md$theta_init
  plan                <- lik$plan
  filter_tunes        <- lik$filter_tunes
  heteroskedastic_shocks <- lik$heteroskedastic_shocks
  freq_band           <- lik$freq_band
  lik_init            <- lik$lik_init
  filter_method       <- lik$filter_method %||% "auto"
  system_priors       <- lik$system_priors
  ## run_mode_finding()'s `...`: log-posterior-constructor extras (student_df
  ## is a likelihood field of the spec) and optimiser arguments.
  ## Names the optimiser takes (its formals after the three positional
  ## arguments) go to the optimiser only; the rest to the constructor.
  .opt_names <- .rmf_optimiser_arg_names()
  lp_dots <- c(lik$extra, md$extra[setdiff(names(md$extra), .opt_names)],
               if (!is.null(lik$student_df)) list(student_df = lik$student_df))

  ## ONE daemon pool spans the whole standard parallel run: the seeded portfolio
  ## (Step 5) provisions it; the Hessian (Step 6) re-binds .worker_lp on it
  ## instead of recompiling the model on every daemon again. shared_pool_ok flags
  ## that a live caller-owned pool exists; it is torn down exactly once on exit
  ## (covers the error path too). build_sigma_prop reads it via pool_ready.
  shared_pool_ok <- FALSE
  on.exit(if (isTRUE(shared_pool_ok)) {
    try(mirai::daemons(NULL), silent = TRUE)
  }, add = TRUE)

  .vcat <- function(...) if (verbose) .dynhr_cat(...)

  model    <- inp$model
  compiled <- inp$compiled
  solved   <- inp$solved
  obs_vars <- spec$obs_vars

  ## Apply unified plan= if supplied (the plan / filter_tunes /
  ## heteroskedastic_shocks conflicts are checked by likelihood_spec()).
  ## Plan adaptation routes through estimation_context()
  ## so that the compiled specs live in ctx$plan for provenance.
  if (!is.null(plan)) {
    .plan_ctx <- estimation_context(plan = plan,
                                    sample_start = model$sample_start,
                                    likelihood   = likelihood)
    filter_tunes           <- attr(.plan_ctx$plan, ".filter_tunes_spec")
    heteroskedastic_shocks <- attr(.plan_ctx$plan, ".shock_scale_spec")
    rm(.plan_ctx)
  }

  ## Apply call-level filter_tunes and heteroskedastic_shocks overrides.
  model <- .resolve_filter_tunes(model, filter_tunes)
  model <- .resolve_heteroskedastic_shocks(model, heteroskedastic_shocks)
  model <- .resolve_stochastic_volatility(model, lik$stochastic_volatility)

  .vcat("\n================================================================\n")
  .vcat("  run_mode_finding\n")
  .vcat("================================================================\n\n")

  # -------------------------------------------------------------------
  # Step 1: Data (loaded from the spec by .est_inputs())
  # -------------------------------------------------------------------
  .vcat("-- Step 1: Load data --\n")
  data <- inp$data
  .vcat(sprintf("  Data: %d x %d  (obs: %s)\n",
                nrow(data), ncol(data), paste(colnames(data), collapse = ", ")))

  ## Expand observables for filter_tunes (no-op when no tunes present).
  tunes_exp <- .expand_observables_for_tunes(model, obs_vars, data)
  obs_vars  <- tunes_exp$obs_vars
  data      <- tunes_exp$Y
  me_extra  <- tunes_exp$me_extra

  ## Build shock_scale matrix from heteroskedastic_shocks block.
  shock_scale_mat <- .build_shock_scale_matrix(
    model, model$varexo_names, nrow(data))

  # -------------------------------------------------------------------
  # Step 2: Extract prior specification (auto eps_* -> sig_* rename)
  # -------------------------------------------------------------------
  .vcat("-- Step 2: Extract prior specification --\n")
  priors <- spec$model$prior_spec %||% extract_prior_spec(model, verbose = verbose)
  .vcat(sprintf("  %d parameters to estimate\n", nrow(priors)))

  # -------------------------------------------------------------------
  # Step 3: Build log-posterior closure
  # -------------------------------------------------------------------
  .vcat("-- Step 3: Build log-posterior --\n")

  ## OBC: likelihood$obc, else auto-detected mcp tags. validate_spec() has
  ## already refused an OBC model with a likelihood only the Kalman filter
  ## implements (the single home of that check).
  use_obc <- .spec_obc_active(spec)

  obc_specs <- NULL
  if (use_obc) {
    ## Announced, not silent: an OBC model is estimated through the OBC/PKF
    ## log-posterior (per the fail-loud-on-fidelity-downgrade convention).
    .vcat("  OBC model detected -- using PKF log-posterior\n")
    .dynhr_inform("run_mode_finding: model has OBC (mcp=) tagged equations; ",
            "auto-switching likelihood from the default \"", likelihood,
            "\" to the OBC/PKF log-posterior (standard Kalman-filter ",
            "likelihoods cannot represent occasionally-binding constraints). ",
            "Pass posterior_options = list(obc_filter = \"ppf\"|\"copf\") to ",
            "use a particle filter instead.")
    obc_specs   <- obc_parse_tags(model)
    log_post_fn <- make_log_posterior_obc_pkf(
      model, data, priors, obs_vars, compiled,
      specs       = obc_specs,
      me_variance = me_variance,
      max_inner   = lik$obc_max_inner
    )
  } else {
    ## CPM routing (v1 limitation, documented on make_log_posterior_tpf's
    ## `burn_in_init` @param): when likelihood$tpf_options$cpm_rho_u is set,
    ## the sampler stage will later dispatch this closure to rwmh_cpm(),
    ## which calls it with a non-NULL U_list on every step after the priming
    ## call. burn_in_init > 0 hard-errors against a supplied U_list (the CPM
    ## slot layout has no burn-in slots), so this is the build site that must
    ## pin burn_in_init = 0L for that case -- unless the caller already passed
    ## burn_in_init explicitly.
    .mf_dots <- lp_dots
    .cpm_rho_u <- lik$tpf_options$cpm_rho_u %||% NULL
    if (identical(likelihood, "tpf") && !is.null(.cpm_rho_u) &&
        is.null(.mf_dots$burn_in_init)) {
      .mf_dots$burn_in_init <- 0L
    }
    ## The target, ONE list: the objective gets it all and the Step-5
    ## gradient gets the part make_posterior_grad() takes (the
    ## gradient used to miss the extras and the tempering exponent passed
    ## through them, and so differentiated another posterior). The exponent
    ## is the typed field (validate_spec() refuses it as an `extra`).
    lp_target <- list(
      me_variance        = me_variance,
      likelihood         = likelihood,
      pruned_order       = pruned_order,
      lik_init           = lik_init,
      filter_method      = filter_method,
      me_extra           = me_extra,
      shock_scale        = shock_scale_mat,
      freq_band          = freq_band,
      system_priors      = system_priors,
      infeasible_penalty = lik$infeasible_penalty,
      power              = lik$power_posterior)
    ## the PSKF CDF setting is a typed field (it changes results); only the
    ## pskf factory takes it
    if (identical(likelihood, "pskf")) lp_target$pskf_cdf <- lik$pskf_cdf
    log_post_fn <- do.call(make_log_posterior,
      c(list(model, data, priors, obs_vars, compiled), lp_target, .mf_dots))
  }

  # -------------------------------------------------------------------
  # Step 4: Starting values -- INITVAL / prior means, or a caller-supplied vector
  # -------------------------------------------------------------------
  ## `theta_init` is assigned ONCE here and every downstream path (serial
  ## optimiser, parallel multi-start portfolio, mirai Hessian at the start,
  ## the recorded result field) reads that one object -- nothing re-derives a
  ## start from `priors$mean`. Overriding it here is therefore sufficient.
  if (is.null(theta_init)) {
    .def_init    <- .rmf_default_theta_init(priors)
    theta_init   <- .def_init$theta
    .init_source <- .def_init$source
    .vcat(sprintf("-- Step 4: Initialise from %s --\n", .init_source))
  } else {
    .vcat("-- Step 4: Initialise from supplied theta_init --\n")
    theta_init <- .rmf_check_theta_init(theta_init, priors)
    .init_source <- "supplied theta_init"
  }

  # Evaluate at the start
  if (isTRUE(optimise)) {
    lp_cal <- log_post_fn(theta_init)
    .vcat(sprintf("  Log-posterior at %s: %.4f\n", .init_source,
                  lp_cal$logpost %||% -Inf))
  } else {
    .vcat("-- Steps 5-6 skipped: the sampler(s) start from the prior and use no mode --\n")
  }

  # -------------------------------------------------------------------
  # Step 5: Mode-finding
  # -------------------------------------------------------------------
  .vcat(sprintf("-- Step 5: Mode-finding (method = '%s', max_iter = %d) --\n",
                method, n_iter))

  mo <- modifyList(list(), mode_options)

  ## Build a temporary ctx to use .ctx_is_standard_gaussian (avoids duplication).
  .tmp_mode_ctx <- estimation_context(
    me_variance  = me_variance,
    likelihood   = likelihood,
    me_extra     = me_extra,
    shock_scale  = shock_scale_mat,
    freq_band    = freq_band,
    pruned_order = pruned_order,
    student_df   = lik$student_df
  )
  par_standard <- .ctx_is_standard_gaussian(.tmp_mode_ctx, use_obc = use_obc,
                                            model = model)
  rm(.tmp_mode_ctx)
  ## Decouple the parallel MULTI-START (Step 5) from the parallel HESSIAN
  ## (Step 6, use_par_hess below -- independent).
  ##
  ## For the H0-seed-dependent methods (newrat, cmaes_newrat) + parallel we
  ## run a SEEDED PORTFOLIO (run_mode_mirai): the H0 seed is computed ONCE on
  ## the host (a parallel FD Hessian at the start) and shipped to the daemons;
  ## chain 1 = newrat @ theta_init with the seed = bit-identical to the serial
  ## run; chains 2+ cycle dispersed starts and the global searchers. The
  ## gradient-free / global methods use the plain parallel multi-start.
  seed_dependent_method <- method %in% c("newrat", "cmaes_newrat")
  use_portfolio <- isTRUE(parallel) && seed_dependent_method && par_standard &&
    requireNamespace("mirai", quietly = TRUE)
  use_par_mode <- isTRUE(parallel) && !seed_dependent_method &&
    requireNamespace("mirai", quietly = TRUE)
  if (isTRUE(parallel) && seed_dependent_method && !par_standard)
    .vcat(sprintf(paste0("  [parallel] method '%s' with a non-standard likelihood: ",
                         "running serially + parallel Hessian only (Step 6).\n"),
                  method))

  ## Multi-start dispersal of the parallel paths: the spec's mode$perturb_scale
  ## and compute$seed (else the seed_base option).
  mode_seed_base <- cmp$seed %||% .dynhr_opt("seed_base", default = 42L)

  # Opt-in unconstrained-parameter transform (eta-space mode finding). Built
  # from `priors` (the prior_spec already in scope) -- see
  # build_param_transform() and .run_mode_finding()'s `transform` argument.
  # Applies to both the serial and parallel (mirai multi-start) paths.
  mode_transform <- NULL
  if (transform_params) {
    mode_transform <- build_param_transform(priors, names(theta_init))
  }

  ## Step 6 differences the EXACT gradient when the mode stage has one
  ## (2n gradient calls instead of the 2n(n+1)-evaluation lp stencil). The
  ## serial Step-5 path builds it (grad_fn below) and Step 6 reuses it; the
  ## per-daemon pool path (no host gradient) builds it on each daemon under
  ## the same eligibility (grad_on_pool).
  grad_fn <- NULL
  grad_on_pool <- isTRUE(mo$use_analytic_grad %||% TRUE) && !use_obc &&
    identical(likelihood, "gaussian") && par_standard

  if (!isTRUE(optimise)) {
    mode_res <- NULL
  } else if (use_portfolio) {
    ## Seeded portfolio: compute the shared H0 seed once on the host (parallel
    ## FD Hessian at the prior mean -- option B), then run the multi-start
    ## portfolio with the seed shipped to every newrat chain.
    ##
    ## POOL SHARING: the seed and the portfolio share ONE daemon pool. The seed
    ## (interior prior mean, many FD evals) wants the fast "stationary" filter;
    ## the mode-search wants the robust "auto". Rather than build + tear down two
    ## pools (recompiling the model on every daemon each time), we compile once
    ## here and let run_mode_mirai RE-BIND .worker_lp to "auto"
    ## (pool_ready = TRUE). The pool is torn down once, in the finally block.
    .vcat("  Parallel seeded portfolio (mirai): shared H0 seed + multi-start (chain 1 = serial replica)\n")
    np_ <- length(theta_init)
    ncs <- .mirai_n_cores(n_cores, np_ * (np_ + 1L) %/% 2L)
    pool_sh <- tryCatch(
      .mirai_pool_init(ncs, model, data, priors, obs_vars, me_variance,
                       me_extra = me_extra, shock_scale = shock_scale_mat,
                       ## the H0 seed is the curvature of the SAME posterior
                       ## (system prior included) the portfolio optimises
                       system_priors = system_priors,
                       lik_init = "stationary",
                       filter_method = filter_method),
      error = function(e) {
        ## A worker/host version skew is not "pool unavailable": every stage
        ## would hit it again. Re-raise it.
        if (inherits(e, "dynhr_error_worker_version_skew")) stop(e)
        ## Nor is a programming error: tear the half-started pool down, then
        ## re-raise it.
        bug <- .dynhr_is_programming_error(e)
        if (!bug)
          .vcat(sprintf("  (shared daemon pool unavailable: %s; each stage self-provisions)\n",
                        conditionMessage(e)))
        try(mirai::daemons(NULL), silent = TRUE)
        if (bug) stop(e)
        NULL
      })
    pool_ok <- !is.null(pool_sh)
    ## Keep this pool ALIVE for Step 6 (build_sigma_prop re-binds .worker_lp on
    ## it). Torn down once by this stage's on.exit. We deliberately do NOT
    ## tear it down here.
    shared_pool_ok <- pool_ok
    par_mode <- tryCatch({
      H0_seed <- if (pool_ok) tryCatch({
          H <- -num_hessian_mirai(NULL, theta_init, h = 1e-4, n_cores = ncs,
                                  verbose = FALSE, pool_ready = TRUE)
          if (any(!is.finite(H))) NULL else H
        }, error = function(e) {
          if (.dynhr_is_programming_error(e)) stop(e)
          .vcat(sprintf("  (shared H0 seed unavailable: %s; chains run unseeded)\n",
                        conditionMessage(e))); NULL
        }) else NULL
      if (!is.null(H0_seed)) .vcat("  Shared H0 seed built; launching portfolio.\n")
      ## The portfolio runs in the SAME space as the serial path (eta-space when
      ## transform_params = TRUE); mode_task passes the THETA-space gradient and
      ## lets .run_mode_finding apply the chain rule, exactly like the host path.
      run_mode_mirai(
        parsed_model   = model, Y = data, prior_spec = priors, obs_names = obs_vars,
        theta_init     = theta_init, n_chains = n_starts, nm_maxit = n_iter,
        method         = method, n_cores = if (pool_ok) ncs else n_cores,
        perturb_scale  = md$perturb_scale, seed_base = mode_seed_base,
        me_variance    = me_variance,
        me_extra       = me_extra, shock_scale = shock_scale_mat,
        ## the daemons' .worker_lp is re-bound with it (pool_ready): without
        ## it the chains optimised the posterior WITHOUT the system prior
        ##
        system_priors  = system_priors,
        transform      = mode_transform,
        analytic_grad  = isTRUE(mo$use_analytic_grad %||% TRUE),
        likelihood     = likelihood, freq_band = freq_band,
        ## the stage's own lik_init: run_mode_mirai re-binds .worker_lp to it
        ## (pool_ready) and the chains build their gradients with it.
        newrat_H0_seed = H0_seed, lik_init = lik_init,
        filter_method = filter_method, pool_ready = pool_ok,
        progress = verbose)
    })
    ## Pool stays alive for Step 6; torn down once by the on.exit above.
    mode_res <- par_mode$best
    mode_res$multistart <- par_mode
  } else if (use_par_mode) {
    .vcat(sprintf("  Parallel multi-start mode-finding (mirai)%s\n",
                  if (!par_standard) " [closure-shipped]" else ""))
    par_mode <- run_mode_mirai(
      parsed_model = if (par_standard) model else NULL,
      Y            = if (par_standard) data  else NULL,
      prior_spec   = priors,
      obs_names    = if (par_standard) obs_vars else NULL,
      theta_init   = theta_init,
      n_chains     = n_starts,
      nm_maxit     = n_iter,
      method       = method,
      n_cores      = n_cores,
      perturb_scale = md$perturb_scale,
      seed_base    = mode_seed_base,
      me_variance  = me_variance,
      me_extra     = me_extra,
      shock_scale  = shock_scale_mat,
      ## the per-daemon .worker_lp (standard path) is built with it; the
      ## shipped closure (non-standard path) already carries it
      system_priors = system_priors,
      log_post_fn  = if (par_standard) NULL else log_post_fn,
      transform    = mode_transform,
      ## Build the analytic gradient on each daemon so the parallel newrat /
      ## cmaes_newrat chains converge (FD-only newrat is far too slow). Only the
      ## standard recompile-per-daemon path can build it. Opt out via
      ## mode_options$use_analytic_grad = FALSE.
      analytic_grad = par_standard && isTRUE(mo$use_analytic_grad %||% TRUE),
      likelihood   = likelihood,
      ## The pool's .worker_lp and the chains' gradients use this init; it
      ## defaulted to "auto" whatever the stage's lik_init was.
      lik_init     = lik_init,
      filter_method = filter_method,
      progress     = verbose
    )
    mode_res <- par_mode$best
    mode_res$multistart <- par_mode
  } else {
    ## P2: build the analytic posterior gradient for the L-BFGS-B polish stage,
    ## so it needs ONE gradient eval per step instead of FD's (n+1) objective
    ## evals. When transform_params = TRUE, the chain rule in .run_mode_finding
    ## converts the theta-space grad_fn into an eta-space gradient
    ## (RC4 fix: gate no longer requires is.null(mode_transform)).
    ## combined_optimize keeps the best point across stages, so even a poor
    ## gradient-polished step never regresses the CMA-ES result. Opt out with
    ## mode_options$use_analytic_grad = FALSE.
    grad_fn <- NULL
    use_analytic_grad <- isTRUE(mo$use_analytic_grad %||% TRUE)
    ## Methods that use a gradient: combined, nelder, cmaes_jade (L-BFGS-B
    ## polish), and the newrat/cmaes_newrat csminwel stages.
    methods_using_grad <- c("combined", "nelder", "cmaes_jade",
                            "newrat", "cmaes_newrat")
    if (use_analytic_grad && !use_obc &&
        method %in% methods_using_grad &&
        likelihood %in% c("gaussian", "cumulant", "whittle")) {
      ## The SAME target as log_post_fn: lp_target's gradient
      ## arguments (power_posterior, system prior, ...) plus the extras
      ## make_posterior_grad() takes. An extra that changes the objective but
      ## that the gradient cannot take leaves the numerical gradient.
      grad_blockers <- intersect(names(.mf_dots), .mode_grad_blocking_extras(likelihood))
      grad_fn <- if (length(grad_blockers)) {
        .vcat(sprintf("  (analytic gradient unavailable: extra argument(s) %s; using numerical FD)\n",
                      paste(grad_blockers, collapse = ", ")))
        NULL
      } else tryCatch(
        do.call(make_posterior_grad, c(
          list(model, data, priors, obs_vars, compiled, verbose = FALSE),
          lp_target[intersect(names(lp_target), names(formals(make_posterior_grad)))],
          .mf_dots[intersect(names(.mf_dots), .mode_grad_extras)])),
        error = function(e) {
          if (.dynhr_is_programming_error(e)) stop(e)
          .vcat(sprintf("  (analytic gradient unavailable: %s; using numerical FD)\n",
                        conditionMessage(e)))
          NULL
        })
      if (!is.null(grad_fn)) {
        if (method %in% c("newrat", "cmaes_newrat")) {
          .vcat("  Using analytic gradient for the newrat/csminwel stage\n")
        } else {
          .vcat("  Using analytic gradient for the L-BFGS-B polish stage\n")
        }
      }
    }
    ## PSKF has no analytic gradient, and optim()'s internal L-BFGS-B finite
    ## differences cannot be kept on one pruning selection: a step across a
    ## selection switch returns jump / (2 * 1e-3). The polish stage gets the
    ## frozen-selection difference instead (.pskf_mode_grad_fn). It is the
    ## optimiser's gradient only: Step 6 must not difference it again (each
    ## call freezes the selection at its own point), so `grad_fn` stays NULL
    ## and the Step-6 Hessian is the frozen num_hessian / numDeriv stencil.
    ## The newrat / cmaes_newrat stages difference through csminwel's own
    ## numerical gradient, which is frozen there.
    mode_grad_fn <- grad_fn
    if (is.null(grad_fn)) {
      mode_grad_fn <- .pskf_mode_polish_grad(likelihood, method, log_post_fn,
                                             priors)
      if (!is.null(mode_grad_fn))
        .vcat("  Using frozen-selection finite differences for the L-BFGS-B polish stage\n")
    }

    ## Build a theta-space analytic posterior Hessian supplier for the newrat
    ## initial H0. Only available for standard Gaussian KF (same gate as
    ## use_exact_hessian). Falls back to csminwel default (1e-4 * I) if
    ## unavailable/fails -- .run_mode_finding handles the NULL case.
    hessian_fn <- NULL
    ## mode_options$use_analytic_hess (this run) over the spec field
    ## mode$analytic_h0 (default: the use_analytic_hess option at build time)
    use_analytic_hess <- isTRUE(mo$use_analytic_hess %||% md$analytic_h0 %||% TRUE)
    h0_method <- mo$h0_method %||% "auto"   # "auto" | "numderiv" | "analytic"
    ## Shared eligibility for an analytic-quality H0 seed (standard Gaussian KF).
    h0_eligible <- method %in% c("newrat", "cmaes_newrat") && !use_obc &&
      identical(likelihood, "gaussian") &&
      is.null(me_extra) && is.null(shock_scale_mat) && !anyNA(data)

    ## Option B: when a parallel pool is available, seed newrat's H0 from a
    ## PARALLEL finite-difference Hessian at theta_init (an interior point:
    ## a supplied start is clipped strictly inside the bounds) instead of the
    ## analytic posterior_hessian. ~20x faster (NZSIM: 205s -> 9s on 16 cores)
    ## and an equivalent seed: csminwel BFGS-refines H0. "auto" uses it
    ## whenever a pool exists; opt out with mode_options$h0_method = "analytic".
    want_parallel_h0 <- h0_eligible && h0_method %in% c("auto", "numderiv") &&
      isTRUE(parallel) && par_standard && requireNamespace("mirai", quietly = TRUE)
    if (want_parallel_h0) {
      hessian_fn <- local({
        .lp <- log_post_fn; .model <- model; .data <- data; .priors <- priors
        .obs <- obs_vars;  .mev <- me_variance; .nc <- n_cores
        .fm <- filter_method
        function(theta) {
          ncs <- .mirai_n_cores(.nc, length(theta) * (length(theta) + 1L) %/% 2L)
          sh <- .mirai_pool_init(ncs, .model, .data, .priors, .obs, .mev,
                                 lik_init = "stationary", filter_method = .fm)
          on.exit({ mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }, add = TRUE)
          ## .lp is unused on the pool_ready path (daemons evaluate .worker_lp);
          ## pass NULL so num_hessian_mirai's task closure never serialises this
          ## heavy posterior closure (see the note there).
          -num_hessian_mirai(NULL, theta, h = 1e-4, n_cores = ncs,
                             verbose = FALSE, pool_ready = TRUE)
        }
      })
      .vcat("  Will build newrat initial H0 from a parallel finite-difference Hessian (option B)\n")
    } else if (use_analytic_hess && h0_eligible &&
               !is.na(h0_init <- .rmf_h0_init_in_force(
                 model, compiled, theta_init, lik_init)) &&
               h0_init != "stationary") {
      ## posterior_hessian()'s curvature kernel starts from the Lyapunov P0,
      ## so it refuses an init that puts another P0 in force at the start
      ## point (lik_init = "kappa"; "diffuse" or "auto" at a unit root). The
      ## seed is then csminwel's default H0 -- chosen here, not by catching
      ## the refusal.
      .vcat(sprintf(paste0("  (analytic newrat H0 not used: lik_init = \"%s\" puts the ",
                           "\"%s\" P0 in force at the start; csminwel default H0)\n"),
                    lik_init, h0_init))
    } else if (use_analytic_hess && h0_eligible) {
      hessian_fn <- tryCatch({
        local({
          .model    <- model
          .compiled <- compiled
          .priors   <- priors
          .obs_vars <- obs_vars
          .data     <- data
          .me_var   <- me_variance
          .sys_pr   <- system_priors
          .lik_init <- lik_init
          function(theta) {
            pm <- .model$param_values
            pm[names(theta)] <- theta
            ss_h <- tryCatch(
              solve_steady_state(.model, .compiled, pm, verbose = FALSE),
              error = function(e) .dynhr_reraise_bug(e, NULL))
            if (is.null(ss_h) || !isTRUE(ss_h$converged))
              stop("steady state did not converge in hessian_fn")
            pm2 <- ss_h$params %||% pm
            dr_h <- solve_perturbation(.model, .compiled, ss_h$ss, pm2,
                                       verbose = FALSE)
            ## posterior_hessian returns the Hessian of log-posterior (negative
            ## curvature matrix); negate it to get neg-logpost Hessian (pos-def
            ## at mode) for .make_pd and csminwel's H0_inv.
            ## check_mode = FALSE: this seeds newrat's INITIAL H0 at whatever
            ## theta the optimizer is at (typically prior means) -- by
            ## construction not a critical point, so the mode-criticality
            ## guard would always warn spuriously here. The at-mode Hessian
            ## call below (use_exact_hessian) keeps the auto-guard.
            H_logpost <- posterior_hessian(
              .model, .compiled, dr_h, pm2,
              names(theta), .obs_vars, t(.data),
              me_variance  = .me_var,
              include_prior = TRUE, prior_spec = .priors,
              system_priors = .sys_pr,
              lik_init = .lik_init,
              ## the seed only needs a good curvature, and the "loop" path
              ## (full second-order solution derivatives) is several times
              ## slower than the d2X-free methods at a Smets-Wouters-sized
              ## model
              t2_method = "auto",
              check_mode = FALSE)
            ## Convert logpost Hessian -> neg-logpost Hessian (flip sign).
            -H_logpost
          }
        })
      }, error = function(e) {
        if (.dynhr_is_programming_error(e)) stop(e)
        .vcat(sprintf("  (analytic Hessian for newrat H0 unavailable: %s)\n",
                      conditionMessage(e)))
        NULL
      })
      if (!is.null(hessian_fn))
        .vcat("  Will build newrat initial H0 from analytic posterior Hessian\n")
    }

    ## `...` is dual-purposed: log-posterior-constructor extras (e.g.
    ## student_df for likelihood="student_t", cut_tol for "pskf") were already
    ## consumed above by make_log_posterior(); strip them here so they don't
    ## also leak into .run_mode_finding()'s optimizer-args `...` (which has no
    ## such formals and errors on an unused argument).
    ## `...` goes to both the log-posterior constructor and the optimiser;
    ## the optimiser takes only its own named options (a likelihood extra
    ## such as cumulant_orders is not one of them).
    dots_optim <- md$extra[intersect(names(md$extra),
                                     names(formals(.run_mode_finding)))]
    ## Opt-in (default FALSE, off by default): stash the newrat/csminwel final
    ## BFGS inverse-Hessian (H0 seed updated by every curvature pair collected
    ## along the optimiser trajectory) in mode_res$H_bfgs, THETA-space. Purely
    ## a measurement/diagnostic hook -- see .run_mode_finding's
    ## `record_curvature` param; no effect on mode-finding itself. Set via
    ## mode_options$record_curvature = TRUE.
    mode_res <- do.call(.run_mode_finding, c(
      list(
        log_post_fn, theta_init, priors,
        nm_maxit    = n_iter,
        method      = method,
        transform   = mode_transform,
        grad_fn     = mode_grad_fn,
        hessian_fn  = hessian_fn,
        verbose     = verbose,
        record_curvature = isTRUE(mo$record_curvature %||% FALSE)
      ),
      dots_optim
    ))
  }

  if (isTRUE(optimise) && (is.null(mode_res) || !is.finite(mode_res$logpost)))
    .dynhr_abort("Mode-finding failed or returned non-finite log-posterior.",
                 class = "dynhr_error_mode_failed")

  theta_mode <- mode_res$theta_mode
  if (isTRUE(optimise))
    .vcat(sprintf("  Log-posterior at mode: %.4f\n", mode_res$logpost))

  # -------------------------------------------------------------------
  # Step 6: Proposal covariance (via Hessian or prior std)
  # -------------------------------------------------------------------
  .vcat("-- Step 6: Build proposal covariance --\n")

  ## Curvature-aware estimation: when use_exact_hessian = TRUE and
  ## the model is the standard Gaussian KF (no OBC/PKF, cumulant/whittle,
  ## me_extra, heteroskedastic shocks, diffuse init, or missing data), compute
  ## the analytic posterior Hessian at the mode (re-solving the decision rule
  ## there) and feed it to the proposal-covariance builder in place of the
  ## finite-difference Hessian. Opt-in; silent fallback to numerical on any
  ## failure or inapplicable model. Cached in the result for NUTS/Laplace reuse.
  hessian_exact <- NULL
  .allows_exact_hess <- identical(likelihood, "gaussian") &&
    !isTRUE(use_obc) && is.null(me_extra) && is.null(shock_scale_mat) &&
    !anyNA(data) &&
    ## posterior_hessian() has no curvature kernel for these inits (it
    ## refuses them; see .rmf_h0_init_in_force)
    !(lik_init %in% c("diffuse", "kappa"))
  if (!isTRUE(proposal) || !isTRUE(optimise)) {
    .vcat("  (skipped: no sampler stage uses it)\n")
  } else if (isTRUE(use_exact_hessian) && .allows_exact_hess) {
    hessian_exact <- tryCatch({
      pm <- model$param_values
      pm[names(theta_mode)] <- theta_mode
      ss_m <- solve_steady_state(model, compiled, pm, verbose = FALSE)
      if (!isTRUE(ss_m$converged))
        stop("steady state did not converge at the mode")
      ## Re-derive SSM-computed params so the exact Hessian uses the consistent
      ## (not stale) p_c (no-op for non-SSM-parameter models).
      pm <- ss_m$params %||% pm
      dr_m <- solve_perturbation(model, compiled, ss_m$ss, pm, verbose = FALSE)
      ## At a bound-constrained (KKT) mode the gradient is non-zero along the
      ## active bounds BY CONSTRUCTION, so posterior_hessian()'s
      ## mode-criticality guard would warn about a legitimate mode; skip it
      ## there (build_sigma_prop decouples those coordinates).
      .bnd6 <- .step6_bounds(priors, names(theta_mode))
      .h6   <- max(1e-4, 1e-4 * abs(theta_mode)) * pmax(1, abs(theta_mode))
      .kkt6 <- any(.step6_bound_active(theta_mode, .bnd6$lower, .bnd6$upper,
                                       2 * .h6))
      H <- posterior_hessian(model, compiled, dr_m, pm, names(theta_mode),
                             obs_vars, t(data), me_variance = me_variance,
                             include_prior = TRUE, prior_spec = priors,
                             system_priors = system_priors,
                             lik_init = lik_init,
                             t2_method = "auto",
                             check_mode = if (.kkt6) FALSE else NULL)
      if (any(!is.finite(H))) stop("non-finite entries in exact Hessian")
      .vcat("  Computed exact posterior Hessian at mode (analytic adjoint).\n")
      H
    }, error = function(e) {
      if (.dynhr_is_programming_error(e)) stop(e)
      .vcat(sprintf("  Exact Hessian unavailable (%s); using numerical.\n",
                    conditionMessage(e)))
      NULL
    })
  } else if (isTRUE(use_exact_hessian)) {
    .vcat("  use_exact_hessian = TRUE but the model is not standard-Gaussian; using numerical.\n")
  }

  V_mode_hess <- NULL
  ## Eta-space proposal covariance, set only at a bound-constrained mode
  ## (R/mode-hessian.R, "Step-6 proposal covariance at a bound-constrained
  ## mode"); NULL otherwise, and then not stored.
  Sigma_prop_eta <- NULL
  if (!isTRUE(proposal) || !isTRUE(optimise)) {
    Sigma_prop <- NULL
  } else if (identical(proposal_cov_method, "full")) {
    .vcat("  Using full-Hessian proposal covariance (method = 'full').\n")
    Sigma_prop <- proposal_cov(
      lp_fn      = log_post_fn,
      theta_mode = theta_mode,
      prior_spec = priors,
      verbose    = verbose,
      method     = "full",
      grad_fn    = grad_fn
    )
    Sigma_prop_eta <- attr(Sigma_prop, "Sigma_eta")
    attr(Sigma_prop, "Sigma_eta") <- NULL
  } else {
    .bsp <- build_sigma_prop(
      return_V     = TRUE,
      lp_fn        = log_post_fn,
      theta_mode   = theta_mode,
      prior_spec   = priors,
      verbose      = verbose,
      parallel     = parallel,
      n_cores      = n_cores,
      hess_exact   = hessian_exact,
      par_standard = par_standard,
      parsed_model = model,
      data         = data,
      obs_vars     = obs_vars,
      me_variance  = me_variance,
      me_extra     = me_extra,
      shock_scale  = shock_scale_mat,
      lik_init     = lik_init,
      filter_method = filter_method,
      ## Re-bind .worker_lp on the live portfolio pool instead of recompiling the
      ## model on every daemon again (Step 5 -> Step 6 pool sharing).
      pool_ready   = shared_pool_ok,
      ## The Hessian from the exact gradient (host / per daemon)
      grad_fn      = grad_fn,
      grad_on_pool = grad_on_pool,
      system_priors = system_priors
    )
    Sigma_prop <- .bsp$Sigma
    V_mode_hess <- .bsp$V_mode
    Sigma_prop_eta <- .bsp$Sigma_eta
  }

  # -------------------------------------------------------------------
  # Wrap up
  # -------------------------------------------------------------------
  .vcat("\n================================================================\n")
  .vcat("  Mode-finding complete\n")
  .vcat("================================================================\n\n")

  ## Build estimation context to store in the result: the ctx carries all
  ## options so they can be read back without loss (system_priors, freq_band).
  ## Store plan provenance in mode_ctx$plan.
  ## When use_obc = TRUE, stamp the OBC filter type (likelihood$obc_filter) so
  ## that conditional_forecast() dispatches to the correct terminal-state path
  ## (ppf/copf stamping). The mode stage always
  ## builds make_log_posterior_obc_pkf.
  mode_ctx_likelihood <- if (use_obc) lik$obc_filter else likelihood
  mode_ctx <- estimation_context(
    me_variance   = me_variance,
    likelihood    = mode_ctx_likelihood,
    lik_init      = lik_init,
    me_extra      = me_extra,
    shock_scale   = shock_scale_mat,
    freq_band     = freq_band,
    system_priors = system_priors,
    tpf_options   = lik$tpf_options,
    obc_specs     = obc_specs,
    student_df    = lik$student_df,
    ## pruned_order=3 is only valid with likelihood="pruned" or "tpf"; the
    ## OBC branch rewrites the ctx likelihood, so fall back to the default there.
    pruned_order  = if (mode_ctx_likelihood %in% c("pruned", "tpf")) pruned_order else 2L
  )
  if (!is.null(plan)) mode_ctx$plan <- plan

  result <- list(
    solved         = solved,
    data           = data,
    obs_vars       = obs_vars,
    prior_spec     = priors,
    log_post_fn    = log_post_fn,
    theta_init     = theta_init,
    theta_mode     = theta_mode,
    mode           = mode_res,
    Sigma_prop     = Sigma_prop,
    V_mode         = mode_res$V_mode %||% V_mode_hess,
    hessian_exact  = hessian_exact,
    log_marglik_laplace = NA_real_,
    me_variance    = me_variance,
    me_extra       = me_extra,
    shock_scale    = shock_scale_mat,
    likelihood     = likelihood,
    freq_band      = freq_band,
    system_priors  = system_priors,
    obc_specs      = obc_specs,
    ctx            = mode_ctx,
    meta           = list(
      method     = method,
      n_iter     = n_iter,
      use_obc    = use_obc
    )
  )
  ## Laplace marginal likelihood: free once the exact posterior
  ## Hessian is available. NA when the Hessian is absent or -H is not PD.
  if (!is.null(hessian_exact))
    result$log_marglik_laplace <- laplace_log_marglik(result)
  if (!is.null(Sigma_prop_eta)) result$Sigma_prop_eta <- Sigma_prop_eta
  if (!isTRUE(optimise)) result$meta$mode_skipped <- TRUE

  class(result) <- c("dynhr_mode_result", "list")
  list(result = result, model = model, compiled = compiled, priors = priors,
       data = data, obs_vars = obs_vars, me_extra = me_extra,
       shock_scale = shock_scale_mat, use_obc = use_obc,
       obc_specs = obc_specs, log_post_fn = log_post_fn)
}


#' Print method for dynhr_mode_result
#' @noRd
#' @export
print.dynhr_mode_result <- function(x, ...) {
  cat("\n<dynhr_mode_result>\n")
  cat(sprintf("  Estimated params: %d\n", nrow(x$prior_spec)))
  cat(sprintf("  Observed vars  : %s\n", paste(x$obs_vars, collapse = ", ")))
  if (isTRUE(x$meta$mode_skipped)) {
    cat("  Mode           : not run (prior-initialised sampler)\n")
  } else {
    cat(sprintf("  Log-posterior  : %.4f\n", x$mode$logpost))
    cat(sprintf("  Method         : %s\n", x$meta$method))
  }
  if (!is.null(x$obc_specs))
    cat(sprintf("  OBC constraints: %d\n", length(x$obc_specs)))
  il <- .est_integrity_lines(x$provenance$integrity)
  if (length(il)) cat(il, sep = "\n")
  invisible(x)
}


#' Build a proposal covariance matrix for MCMC
#'
#' Constructs \code{Sigma_prop} from the negative inverse Hessian at the
#' mode, scaled by the optimal RWMH factor \eqn{2.38^2 / n_{par}}.  Falls
#' back to a diagonal matrix from prior standard deviations when the
#' Hessian is unavailable or non-positive-definite.
#'
#' @param lp_fn       Log-posterior closure.
#' @param theta_mode  Mode parameter vector.
#' @param prior_spec  Prior spec data.frame.
#' @param pooled_draws Optional pooled draw matrix for empirical covariance.
#' @param verbose     Print progress.
#' @param parallel    Evaluate the numerical Hessian on a mirai daemon pool
#'   (default \code{FALSE}). The log-posterior closure \code{lp_fn} is
#'   shipped to the daemons as-is, so this works for ANY likelihood
#'   (Gaussian, OBC/PKF, cumulant) -- unlike the MCMC/NUTS/SMC mirai paths,
#'   which require a parsed model to recompile per daemon.
#' @param n_cores     Daemon count for \code{parallel = TRUE} (\code{NULL} =
#'   auto, capped at the number of Hessian evaluation tasks).
#' @return Proposal covariance matrix.
#' @noRd
build_sigma_prop <- function(lp_fn, theta_mode, prior_spec,
                              pooled_draws = NULL, verbose = TRUE,
                              parallel = FALSE, n_cores = NULL,
                              hess_exact = NULL,
                              ## For the parallel Hessian's per-daemon model
                              ## re-init (avoids shipping the dead-pointer
                              ## closure); only used when par_standard = TRUE.
                              par_standard = FALSE, parsed_model = NULL,
                              data = NULL, obs_vars = NULL, me_variance = 0,
                              me_extra = NULL, shock_scale = NULL,
                              ## pool_ready = TRUE: the caller has already stood up
                              ## a daemon pool (model compiled per daemon) and owns
                              ## teardown; we only RE-BIND .worker_lp to the filter
                              ## we need (stationary/auto) instead of a fresh pool.
                              pool_ready = FALSE,
                              ## The log-posterior's lik_init: the daemon
                              ## filter when "stationary" disagrees with lp_fn
                              ## at the mode.
                              lik_init = "auto",
                              ## The log-posterior's Kalman recursion (value
                              ## path), for the same daemon filter.
                              filter_method = "auto",
                              ## return_V = TRUE: return list(Sigma, V_mode) so
                              ## the caller can keep the regularised inverse
                              ## Hessian (documented as result$V_mode) instead
                              ## of discarding it. V_mode is NULL on the
                              ## empirical-covariance path.
                              return_V = FALSE,
                              ## The mode stage's theta-space gradient
                              ## (make_posterior_grad()). When usable the
                              ## Hessian is the central difference of it
                              ## (grad_hessian(): 2n calls) instead of the
                              ## 2n(n+1)-evaluation lp stencil.
                              grad_fn = NULL,
                              ## The per-daemon (par_standard) pool path
                              ## may build the gradient on each daemon from
                              ## the daemon's model (the Step-5 eligibility:
                              ## analytic gradient not opted out, no OBC,
                              ## no blocking extras) and difference it there.
                              grad_on_pool = FALSE,
                              ## The system prior the daemon posterior and
                              ## gradient carry on the pool path.
                              system_priors = NULL) {
  n_par <- length(theta_mode)

  .vcat <- function(...) if (verbose) .dynhr_cat(...)
  .out <- function(Sigma, V_mode = NULL, Sigma_eta = NULL) {
    if (!isTRUE(return_V)) return(Sigma)
    res <- list(Sigma = Sigma, V_mode = V_mode)
    if (!is.null(Sigma_eta)) res$Sigma_eta <- Sigma_eta
    res
  }

  # Priority 1: empirical covariance from pooled draws
  if (!is.null(pooled_draws) && nrow(pooled_draws) > n_par) {
    .vcat("  Using empirical covariance from pooled draws.\n")
    Sigma <- cov(pooled_draws)
    opt_scale <- 2.38^2 / n_par
    rownames(Sigma) <- colnames(Sigma) <- names(theta_mode)
    return(.out(Sigma * opt_scale))
  }

  # Priority 2: Hessian at mode. When an exact posterior Hessian (analytic
  # adjoint, log-posterior i.e. include_prior = TRUE) is supplied it is used
  # directly in place of the finite-difference Hessian; all downstream repair
  # (ill-conditioning fallback, prior-scale cap, PD regularisation) is identical.
  #
  # Bound-constrained (KKT) mode: coordinates whose central stencil (num_hessian
  # reaches 2 h_i) leaves the prior support are differentiated one-sided and
  # decoupled (.step6_bound_hessian, R/mode-hessian.R); the stencil for the
  # rest never moves them. No such coordinate -> the path below is unchanged.
  #
  # With an exact gradient (grad_fn on the host, or grad_on_pool on the
  # per-daemon pool) the Hessian is the symmetrised central difference of the
  # gradient (grad_hessian(): 2n calls, R/mode-hessian.R) -- the same
  # bound-active selection, the bound coordinates' one-sided curvature from
  # their own gradient component. A gradient stencil with a non-finite entry
  # falls back to the log-posterior stencil below.
  h <- max(1e-4, 1e-4 * abs(theta_mode))
  h_vec  <- h * pmax(1, abs(theta_mode))
  bnds   <- .step6_bounds(prior_spec, names(theta_mode))
  at_bnd <- .step6_bound_active(theta_mode, bnds$lower, bnds$upper, 2 * h_vec)
  coords <- which(!at_bnd)
  os_grad <- NULL
  .grad_ok <- function(gh) {
    ok <- all(is.finite(gh$hess))
    if (!ok)
      .vcat("  (gradient stencil has non-finite entries; using the log-posterior stencil)\n")
    ok
  }
  if (!is.null(hess_exact)) {
    .vcat("  Using exact posterior Hessian at mode (analytic adjoint).\n")
    hess <- hess_exact
  } else {
    hess <- NULL
    if (.step6_grad_usable(grad_fn)) {
      .vcat(sprintf("  Estimating Hessian at mode (central differences of the analytic gradient, %d calls)...\n",
                    2L * length(coords) + if (any(at_bnd)) 1L + 2L * sum(at_bnd) else 0L))
      gh <- grad_hessian(grad_fn, theta_mode, coords = coords, at_bnd = at_bnd,
                         bnds = bnds)
      if (.grad_ok(gh)) { hess <- gh$hess; os_grad <- gh$os }
    }
    use_par_hess <- isTRUE(parallel) && length(coords) > 1L &&
                    requireNamespace("mirai", quietly = TRUE)
    if (is.null(hess) && use_par_hess && par_standard) {
      ## Re-initialise the model natively on each daemon (.mirai_pool_init) and
      ## reuse that pool, instead of shipping the log-posterior closure (whose
      ## compiled-model Rcpp external pointers serialise to DEAD pointers on the
      ## daemon -> a ~45x-slower per-eval fallback: 1908s vs 166s on NZSIM).
      ##
      ## Prefer lik_init = "stationary" for the daemon evals: at a stationary
      ## mode (including the near-unit-root NZSIM mode) it is exact -- identical
      ## likelihood to "auto"/"diffuse" -- but a fresh daemon built with "auto"
      ## re-runs the eigen decision every eval and resolves a near-unit-root
      ## mode to the slow EXACT-DIFFUSE filter (~500 ms/eval vs stationary's
      ## ~5-20 ms). Validate it ONCE against the host log-posterior at the mode;
      ## only a genuine unit-root model (where the filters disagree) needs the
      ## diffuse path, in which case re-init with "auto". Non-finite entries
      ## (e.g. boundary -Inf perturbations) flow to the shared make_pd handler
      ## below, exactly as the serial numDeriv path does -- no separate re-run.
      .vcat("  Estimating Hessian at mode (mirai, parallel; per-daemon re-init)...\n")
      n_cores_h <- .mirai_n_cores(n_cores,
                                  length(coords) * (length(coords) + 1L) %/% 2L)
      host_lp_mode <- lp_fn(theta_mode)$logpost %||% -Inf
      ## Set the daemon filter: on a SHARED pool (pool_ready, caller-owned) just
      ## re-bind .worker_lp -- no model recompile; otherwise stand up a fresh
      ## pool here. Returns the share handle (NULL on the rebind path). The
      ## system prior is part of the target (the host lp_fn carries it).
      .set_filter <- function(lik) {
        if (isTRUE(pool_ready)) {
          .mirai_rebind_worker_lp(prior_spec, obs_vars, me_variance = me_variance,
                                  me_extra = me_extra, shock_scale = shock_scale,
                                  system_priors = system_priors,
                                  lik_init = lik, filter_method = filter_method)
          NULL
        } else {
          .mirai_pool_init(n_cores_h, parsed_model, data, prior_spec, obs_vars,
                           me_variance, me_extra = me_extra,
                           shock_scale = shock_scale,
                           system_priors = system_priors, lik_init = lik,
                           filter_method = filter_method)
        }
      }
      fit_filter <- "stationary"
      sh <- .set_filter("stationary")
      stat_lp_mode <- tryCatch({
        .tm <- theta_mode
        mirai::mirai({
          lpf <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
          r <- lpf(.tm); if (is.list(r)) r$logpost else r
        }, .tm = .tm)[]
      }, error = function(e) .dynhr_reraise_bug(e, NA_real_))
      if (!is.finite(stat_lp_mode) ||
            abs(stat_lp_mode - host_lp_mode) > 1e-4 * max(1, abs(host_lp_mode))) {
        .vcat(sprintf("  Stationary filter disagrees at the mode; re-initialising with lik_init = %s.\n",
                      lik_init))
        ## Only tear down a pool we OWN; a shared pool is just re-bound.
        if (!isTRUE(pool_ready)) { mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }
        sh <- .set_filter(lik_init)
        fit_filter <- lik_init
      }
      .vcat(sprintf("    (filter = %s, %d daemons%s)\n", fit_filter, n_cores_h,
                    if (isTRUE(pool_ready)) ", shared pool" else ""))
      ## The exact gradient of the daemon posterior, built per daemon from
      ## its model (.mirai_bind_worker_grad), differenced over the pool. Not
      ## used when the daemon gradient resolves to "hybrid" or is unavailable.
      if (isTRUE(grad_on_pool)) {
        gm <- .mirai_bind_worker_grad(prior_spec, obs_vars,
                                      me_variance = me_variance,
                                      me_extra = me_extra,
                                      shock_scale = shock_scale,
                                      system_priors = system_priors)
        if (!is.na(gm) && !identical(gm, "hybrid")) {
          .vcat("    (central differences of the analytic gradient on the pool)\n")
          gh <- grad_hessian_mirai(theta_mode, coords = coords, at_bnd = at_bnd,
                                   bnds = bnds, verbose = verbose)
          if (.grad_ok(gh)) { hess <- gh$hess; os_grad <- gh$os }
        }
      }
      if (is.null(hess))
        hess <- tryCatch(
          num_hessian_mirai(lp_fn, theta_mode, h = h, n_cores = n_cores_h,
                            verbose = verbose, pool_ready = TRUE,
                            coords = coords),
          error = function(e) .dynhr_reraise_bug(e,
            matrix(NA_real_, n_par, n_par)))
      ## Tear down only a pool we created here; the caller owns a shared pool.
      if (!isTRUE(pool_ready)) { mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }
    } else if (is.null(hess) && use_par_hess) {
      ## Non-standard likelihood (me_extra / shock_scale / non-Gaussian): the
      ## per-daemon re-init path does not cover it, so ship the closure.
      .vcat("  Estimating Hessian at mode (mirai, parallel)...\n")
      hess <- num_hessian_mirai(lp_fn, theta_mode, h = h, n_cores = n_cores,
                                verbose = verbose, coords = coords)
    } else if (is.null(hess)) {
      .vcat("  Estimating Hessian at mode (numDeriv)...\n")
      hess <- num_hessian(lp_fn, theta_mode, h = h, coords = coords)
    }
  }
  bh <- NULL
  if (any(at_bnd)) {
    bh <- .step6_bound_hessian(hess, lp_fn, theta_mode, at_bnd, h_vec, bnds,
                               prior_spec, keep_diag = !is.null(hess_exact),
                               verbose = verbose, os = os_grad)
    hess <- bh$hess
  }
  .pc   <- .proposal_cov_from_hessian(hess, prior_spec, theta_mode, verbose)
  ## Bound-active mode: also the eta-space proposal (see R/mode-hessian.R).
  Sigma_eta <- if (!is.null(bh))
    .step6_eta_cov(.pc$Sigma, bh, theta_mode, prior_spec, scale = 2.38^2 / n_par)
  .out(.pc$Sigma, .pc$V_mode, Sigma_eta)
}



#' Posterior-curvature proposal covariance from a Hessian at the mode
#'
#' Shared by \code{run_mode_finding()} and the one-call estimation entry
#' points (\code{run_full_estimation()}, \code{estimate-runner.R}). Extracted
#' 2026-09-04: the two samplers used to test \code{mode_res$V_mode}, which
#' \code{.run_mode_finding()} never sets, so the Hessian branch was DEAD and
#' every RWMH proposal silently came from prior variances. On a well-identified
#' posterior that froze the chain outright (0\% acceptance, zero posterior
#' variance) while the inverse-Hessian proposal sampled at 23.5\%.
#'
#' @param hess Hessian of the LOG-posterior at the mode (negative definite
#'   there), n_par x n_par.
#' @param prior_spec Prior specification; \code{$std} supplies the per-
#'   parameter scale used to cap genuinely flat eigen-directions.
#' @param theta_mode Named mode vector (supplies the dimnames).
#' @param verbose Passed through to the \code{.vcat()} progress notes.
#' @return \code{list(V_mode, Sigma)} — the regularised inverse-Hessian
#'   posterior covariance and the scaled RWMH proposal covariance.
#' @noRd
.proposal_cov_from_hessian <- function(hess, prior_spec, theta_mode,
                                       verbose = TRUE) {
  n_par <- length(theta_mode)
  ## `.vcat` is a LOCAL closure in each caller (run-mode-finding.R:251, :945),
  ## not a package-level function, so it must be re-made here or the
  ## regularisation branches below would fail to find it.
  .vcat <- function(...) if (verbose) .dynhr_cat(...)

  if (any(!is.finite(hess)) ||
      !(is.finite(rcond(-hess)) && rcond(-hess) > .Machine$double.eps)) {
    ## Non-finite entries or ill-conditioned: use .make_pd to preserve finite
    ## local curvature and assign prior-var to bad coordinates, instead of
    ## discarding all curvature information and falling back to a pure
    ## prior-std diagonal.
    if (any(!is.finite(hess))) {
      .vcat("  Hessian has non-finite values; applying .make_pd regularisation.\n")
    } else {
      .vcat("  Hessian is ill-conditioned; applying .make_pd regularisation.\n")
    }
    pv <- prior_spec$std^2
    pv[!is.finite(pv) | pv <= 0] <- 1
    ## .make_pd expects the negative-Hessian (positive at mode); hess here IS
    ## the negative log-posterior Hessian (returned by num_hessian as H of logpost,
    ## so it should be negative definite at the mode; we negate for .make_pd).
    V_mode <- .make_pd(-hess, cond_target = 100, prior_var = pv)
  } else {
    V_mode <- solve(-hess)
  }

  # Cap genuinely flat eigen-directions at the prior scale.
  #
  # RATIONALE FOR EIGEN-BASIS CAPPING (vs. old per-coordinate marginal cap):
  #
  # The old approach computed shrink = sqrt(pmin(1, prior_var / diag(V_mode)))
  # and applied Ds %*% V_mode %*% Ds. This caps the MARGINAL variance of each
  # parameter. A parameter that is weakly identified *marginally* but tightly
  # constrained *conditionally* (via a strong posterior correlation with another
  # parameter) has a large marginal variance but a tiny conditional variance.
  # Capping its marginal variance to the prior scale collapses the full-Hessian
  # conditional structure, widening the proposal far beyond the posterior ridge
  # (e.g. NZSIM's phipi: marginal ±1.9 vs. Dynare-posterior ±0.01).
  #
  # The correct criterion for "flat direction" is weak curvature in the
  # EIGEN-BASIS of the Hessian (i.e., a large eigenvalue of V_mode = H^{-1}).
  # The prior variance in that eigen-direction is q'*diag(prior_var)*q, where q
  # is the eigenvector. An eigen-direction is only "flat" if its posterior
  # variance (eigenvalue) exceeds the prior variance *in that same direction*
  # — this cannot happen for a tightly conditionally-identified parameter
  # because its tight conditional constraints live in the small-eigenvalue
  # (well-identified) directions of V_mode.
  #
  # Result: the cap now preserves off-diagonal correlations and conditional
  # precision, while still bounding genuinely uninformative eigen-directions
  # (e.g. uniform-prior flat directions where the Hessian curvature is ~0).
  prior_var <- prior_spec$std^2
  prior_var[!is.finite(prior_var) | prior_var <= 0] <- Inf  # no cap if undefined
  eig_v <- eigen(V_mode, symmetric = TRUE)
  Q   <- eig_v$vectors           # columns are eigenvectors
  lam <- eig_v$values            # eigenvalues of V_mode (>= 0 after .make_pd)
  # Prior variance in each eigen-direction: q_j' diag(prior_var) q_j
  prior_var_dir <- colSums(Q^2 * prior_var)   # length n_par
  prior_var_dir[!is.finite(prior_var_dir)] <- Inf  # no cap where prior is undefined
  lam_capped <- pmin(lam, prior_var_dir)
  n_capped <- sum(lam_capped < lam - 1e-12 * pmax(abs(lam), 1))
  if (n_capped > 0L) {
    V_mode <- Q %*% diag(lam_capped, nrow = n_par) %*% t(Q)
    V_mode <- (V_mode + t(V_mode)) / 2          # re-symmetrise after reconstruction
    .vcat(sprintf("  Capped %d eigen-direction(s) at the prior scale (flat-likelihood directions).\n",
                  n_capped))
  }

  opt_scale <- 2.38^2 / n_par
  Sigma <- V_mode * opt_scale
  rownames(Sigma) <- colnames(Sigma) <- names(theta_mode)

  # Ensure positive definiteness
  eig <- eigen(Sigma, symmetric = TRUE)
  if (min(eig$values) <= 0) {
    .vcat("  Regularising non-positive-definite covariance.\n")
    Sigma <- eig$vectors %*% diag(pmax(eig$values, 1e-6), nrow = n_par) %*% t(eig$vectors)
  }

  list(V_mode = V_mode, Sigma = Sigma)
}

#' Numerical Hessian via central differences
#'
#' @param coords Coordinates to differentiate (default all). Entries outside
#'   \code{coords} x \code{coords} are left 0 -- used by Step 6 to keep the
#'   stencil off a bound-active coordinate (see .step6_bound_hessian).
#'
#' A PSKF log-posterior is differenced on the pruning selection made at
#' \code{theta} (evaluated first, so it records; see .pskf_freeze_open): a
#' second difference across a selection switch would otherwise return
#' jump / h^2. Functions that run no PSKF filter are unaffected.
#' @noRd
num_hessian <- function(fn, theta, h = 1e-4, coords = seq_along(theta)) {
  n <- length(theta)
  H <- matrix(0, nrow = n, ncol = n)
  fr <- .pskf_freeze_open()
  on.exit(.pskf_freeze_close(fr), add = TRUE)
  fn <- .pskf_freeze_wrap(fn)
  f0 <- fn(theta)
  f0_val <- if (is.list(f0)) f0$logpost else f0

  for (a in seq_along(coords)) {
    i <- coords[a]
    for (j in coords[a:length(coords)]) {
      th_pp <- th_pm <- th_mp <- th_mm <- theta
      hi <- h * max(1, abs(theta[i]))
      hj <- h * max(1, abs(theta[j]))
      th_pp[i] <- th_pp[i] + hi; th_pp[j] <- th_pp[j] + hj
      th_pm[i] <- th_pm[i] + hi; th_pm[j] <- th_pm[j] - hj
      th_mp[i] <- th_mp[i] - hi; th_mp[j] <- th_mp[j] + hj
      th_mm[i] <- th_mm[i] - hi; th_mm[j] <- th_mm[j] - hj

      f_pp <- fn(th_pp); f_pp <- if (is.list(f_pp)) f_pp$logpost else f_pp
      f_pm <- fn(th_pm); f_pm <- if (is.list(f_pm)) f_pm$logpost else f_pm
      f_mp <- fn(th_mp); f_mp <- if (is.list(f_mp)) f_mp$logpost else f_mp
      f_mm <- fn(th_mm); f_mm <- if (is.list(f_mm)) f_mm$logpost else f_mm

      H[i, j] <- H[j, i] <- (f_pp - f_pm - f_mp + f_mm) / (4 * hi * hj)
    }
  }
  H
}


#' Parallel numerical Hessian via central differences on a mirai daemon pool
#'
#' Drop-in replacement for \code{\link{num_hessian}} that ships the
#' log-posterior closure \code{fn} to a mirai daemon pool once (via
#' \code{\link{.mirai_pool_closure}}) and evaluates each parameter pair's
#' four-point central-difference stencil in parallel. Unlike the MCMC/NUTS/SMC
#' mirai paths, this needs no parsed model to recompile per daemon -- \code{fn}
#' is used as-is, so it works for Gaussian, OBC/PKF, and cumulant likelihoods
#' alike.
#'
#' @param fn       log-posterior function: theta -> list(logpost, ...) or scalar.
#' @param theta    Named parameter vector (the mode).
#' @param h        Base step-size fraction (as in \code{\link{num_hessian}}).
#' @param n_cores  Daemon count (\code{NULL} = auto, capped at the number of
#'   parameter pairs).
#' @param verbose  Print timing.
#' @param coords   Coordinates to differentiate (default all), as in
#'   \code{\link{num_hessian}}; other entries are left 0.
#' @return n x n Hessian matrix (unnamed dimnames, as \code{num_hessian}).
#' @noRd
num_hessian_mirai <- function(fn, theta, h = 1e-4, n_cores = NULL,
                              verbose = TRUE, pool_ready = FALSE,
                              coords = seq_along(theta)) {
  n <- length(theta)
  par_names <- names(theta)

  nc <- length(coords)
  pairs <- vector("list", nc * (nc + 1L) %/% 2L)
  k <- 0L
  for (a in seq_len(nc)) for (j in coords[a:nc]) {
    k <- k + 1L
    pairs[[k]] <- c(coords[a], j)
  }

  n_cores <- .mirai_n_cores(n_cores, length(pairs))
  if (verbose) .dynhr_cat(sprintf("    %d Hessian evaluations on %d daemons\n",
                           length(pairs) * 4L, n_cores))

  t0 <- proc.time()
  ## pool_ready = TRUE: the caller has already stood up the daemon pool with
  ## `.worker_lp` set (e.g. via .mirai_pool_init, which re-COMPILES the model
  ## natively per daemon) and manages teardown. Preferred for compiled-model
  ## posteriors: shipping the closure here (.mirai_pool_closure) serialises the
  ## model's Rcpp external pointers to dead pointers on the daemon, forcing a
  ## slow per-eval fallback (~45x slower -- 1908s vs 166s on NZSIM).
  if (!pool_ready) {
    .mirai_pool_closure(n_cores, fn)
    on.exit(mirai::daemons(NULL), add = TRUE)
  }

  ## PSKF: the pruning selection the log-posterior makes at theta, recorded
  ## here in the main process (fn is at hand) and replayed by every daemon
  ## evaluation (.num_hessian_stencil seeds its frozen-selection scope with
  ## it). Seeding needs no change to the posterior closure's interface, which
  ## an explicit keep_override argument would. An empty list -- not a PSKF
  ## posterior, or pool_ready, where the daemons hold a posterior recompiled
  ## per daemon (standard Gaussian models only) -- leaves the stencil as
  ## before.
  keep_paths <- if (!pool_ready && is.function(fn))
    .pskf_record_centre(fn, theta) else list()

  ## .num_hessian_stencil is a namespace function, so mirai_map serialises
  ## its environment as a reference to the daemon's dynhr namespace, not as
  ## this frame. A closure over this frame would ship `fn` -- the passed
  ## log-posterior -- with EVERY one of the n(n+1)/2 tasks; for a real
  ## make_posterior() closure (heavy adjoint/model state) that was ~50x
  ## slower per eval inside a full run_mode_finding than standalone (NZSIM:
  ## 270s vs 5s).
  raw <- mirai::mirai_map(
    pairs, .num_hessian_stencil,
    .args = list(.theta = theta, .h = h, .par_names = par_names,
                 .keep_paths = keep_paths)
  )[]

  H <- matrix(0, nrow = n, ncol = n)
  for (k in seq_along(pairs)) {
    p <- pairs[[k]]
    v <- raw[[k]]
    if (inherits(v, "miraiError") || inherits(v, "errorValue")) v <- NA_real_
    H[p[1L], p[2L]] <- H[p[2L], p[1L]] <- v
  }

  if (verbose) .dynhr_cat(sprintf("    Hessian: %.1f sec\n",
                           (proc.time() - t0)[["elapsed"]]))
  H
}

#' One four-point central-difference Hessian entry on a mirai daemon
#'
#' Evaluates the daemon's \code{.worker_lp} on the stencil of the parameter
#' pair \code{pair} (see \code{num_hessian_mirai}). \code{.keep_paths}
#' (from \code{.pskf_record_centre} in the main process) seeds a
#' frozen-selection scope, so a PSKF posterior is differenced on the pruning
#' selection made at the centre; an empty list evaluates plainly.
#' @noRd
.num_hessian_stencil <- function(pair, .theta, .h, .par_names,
                                 .keep_paths = list()) {
  lp <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
  .eval <- function(th) {
    names(th) <- .par_names
    r <- lp(th)
    if (is.list(r)) r$logpost else r
  }
  if (length(.keep_paths)) {
    fr <- .pskf_freeze_open(.keep_paths)
    on.exit(.pskf_freeze_close(fr), add = TRUE)
    .eval <- .pskf_freeze_wrap(.eval)
  }
  i <- pair[1L]; j <- pair[2L]
  hi <- .h * max(1, abs(.theta[i]))
  hj <- .h * max(1, abs(.theta[j]))
  th_pp <- th_pm <- th_mp <- th_mm <- .theta
  th_pp[i] <- th_pp[i] + hi; th_pp[j] <- th_pp[j] + hj
  th_pm[i] <- th_pm[i] + hi; th_pm[j] <- th_pm[j] - hj
  th_mp[i] <- th_mp[i] - hi; th_mp[j] <- th_mp[j] + hj
  th_mm[i] <- th_mm[i] - hi; th_mm[j] <- th_mm[j] - hj
  f_pp <- .eval(th_pp); f_pm <- .eval(th_pm)
  f_mp <- .eval(th_mp); f_mm <- .eval(th_mm)
  (f_pp - f_pm - f_mp + f_mm) / (4 * hi * hj)
}


#' Parallel Step-6 Hessian from the exact gradient on a mirai daemon pool
#'
#' The pool counterpart of \code{grad_hessian()} (R/mode-hessian.R): the
#' same evaluation points (\code{.step6_grad_plan}) and the same assembly
#' (\code{.step6_grad_finish}), the gradient calls spread over the daemons.
#' The caller has stood up the pool (\code{.mirai_pool_init}) and bound each
#' daemon's gradient with \code{.mirai_bind_worker_grad()}; a daemon without
#' one (or a failed task) yields an NA gradient, i.e. a non-finite Hessian
#' the caller treats as "gradient unavailable".
#'
#' @param theta Named parameter vector (the mode).
#' @param h Base relative step, as in \code{grad_hessian()}.
#' @param coords,at_bnd,bnds As in \code{grad_hessian()}.
#' @param verbose Print timing.
#' @return list(hess, os) as \code{grad_hessian()}.
#' @noRd
grad_hessian_mirai <- function(theta, h = .step6_grad_h,
                               coords = seq_along(theta),
                               at_bnd = rep(FALSE, length(theta)),
                               bnds = NULL, verbose = TRUE) {
  n    <- length(theta)
  h_g  <- h * pmax(1, abs(theta))
  plan <- .step6_grad_plan(theta, coords, at_bnd, h_g, bnds)
  if (verbose) .dynhr_cat(sprintf("    %d gradient evaluations on the pool\n",
                                  length(plan$points)))
  t0 <- proc.time()
  task <- function(th, .ev) {
    gfn <- get0(".worker_grad", envir = globalenv(), inherits = FALSE)
    if (!is.function(gfn)) return(NULL)
    .ev(gfn, th)
  }
  ## Severed from this frame for the reason given in num_hessian_mirai(): the
  ## task needs nothing from it. `.ev` is the host's .step6_grad_eval (its
  ## namespace environment resolves to the daemon's dynhr).
  environment(task) <- globalenv()
  raw <- mirai::mirai_map(plan$points, task,
                          .args = list(.ev = .step6_grad_eval))[]
  G <- lapply(raw, function(v)
    if (is.numeric(v) && length(v) == n) as.numeric(v) else rep(NA_real_, n))
  if (verbose) .dynhr_cat(sprintf("    Hessian: %.1f sec\n",
                                  (proc.time() - t0)[["elapsed"]]))
  .step6_grad_finish(plan, G)
}
