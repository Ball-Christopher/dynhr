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
#'   noise diagonal (default 0).  Use a small positive value for
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
#'   deliberately NOT accepted here -- their estimation noise breaks
#'   deterministic optimizers (Nelder-Mead/CMA-ES/newrat all assume a fixed
#'   objective at repeated evaluations of the same point); use
#'   \code{\link{run_full_estimation}} or particle MCMC (PMMH) via
#'   \code{\link{run_posterior_estimation}} instead.
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
#'   \code{.run_mode_finding}'s `transform` argument -- see there for
#'   the exact semantics (Jacobian-free objective, eta-space bounds for
#'   bounded optimizer stages, `theta_mode` returned in theta-space as
#'   always). Applies to both the serial path and the PARALLEL multi-start
#'   path (\code{run_mode_mirai}'s `transform` argument): each
#'   daemon's \code{.run_mode_finding} call receives the transform, and
#'   dispersed multi-start points are jittered in eta-space (no bound
#'   clamping) before being mapped back to theta-space starts. May also be
#'   set globally via \code{dynhr_set_options(transform_params = TRUE)}. Has
#'   no effect on Step 6 (proposal covariance), which is always computed in
#'   theta-space here; \code{run_posterior_estimation()} performs its own
#'   delta-method conversion when sampling in eta-space.
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
#' @param use_exact_hessian  Opt-in (default \code{FALSE}): build the RWMH
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
#'   \code{dynhr_set_options(use_analytic_hess = FALSE)} to fall back to the
#'   cheap default \code{H0} if mode-finding setup is the bottleneck.)
#' @param verbose  Print progress messages (default TRUE).
#' @param ...      Additional arguments forwarded to the mode-finder.
#'
#' @return An object of class \code{"dynhr_mode_result"} containing:
#'   \describe{
#'     \item{\code{solved}}{The input \code{dynhr_solved} object}
#'     \item{\code{data}}{Observation matrix}
#'     \item{\code{obs_vars}}{Observable variable names}
#'     \item{\code{prior_spec}}{Prior specification data.frame}
#'     \item{\code{log_post_fn}}{Log-posterior closure}
#'     \item{\code{theta_init}}{Starting parameter vector (prior means)}
#'     \item{\code{theta_mode}}{Posterior mode parameter vector}
#'     \item{\code{mode}}{Full mode-finding result list}
#'     \item{\code{Sigma_prop}}{Proposal covariance for MCMC}
#'     \item{\code{V_mode}}{Estimated variance-covariance at mode}
#'     \item{\code{me_variance},\code{likelihood}}{Settings used}
#'   }
#'
#' @examples
#' \dontrun{
#' mod  <- solve_model("my_model.mod")
#' data <- read.csv("data.csv")
#' mode <- run_mode_finding(mod, data, obs_vars = c("y", "pi", "r"))
#' }
#'
#' @seealso \code{\link{solve_model}}, \code{\link{run_posterior_estimation}},
#'   \code{\link{find_mode}}, \code{\link{make_posterior}}
#' @export
run_mode_finding <- function(solved,
                             data,
                             obs_vars = NULL,
                             n_iter           = 10000L,
                             method           = "newrat",
                             me_variance      = 0,
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
                             use_exact_hessian = FALSE,
                             verbose          = TRUE,
                             ...) {

  ## Capture BEFORE match.arg() collapses the default vector to its first
  ## element: this is the only reliable "did the caller pass likelihood=
  ## explicitly?" signal (mirrors the missing()/sentinel convention used
  ## elsewhere in this codebase for optional args with non-NULL defaults).
  likelihood_explicit <- !missing(likelihood)
  likelihood <- match.arg(likelihood)
  ## For the Whittle path, freq_band can be passed via posterior_options.
  ## Extract it here so it is threaded into make_log_posterior().
  freq_band <- posterior_options$freq_band %||% c(0, pi)

  ## ONE daemon pool spans the whole standard parallel run: the seeded portfolio
  ## (Step 5) provisions it; the Hessian (Step 6) re-binds .worker_lp on it
  ## instead of recompiling the model on every daemon again. shared_pool_ok flags
  ## that a live caller-owned pool exists; it is torn down exactly once on exit
  ## (covers the error path too). build_sigma_prop reads it via pool_ready.
  shared_pool_ok <- FALSE
  on.exit(if (isTRUE(shared_pool_ok)) {
    try(mirai::daemons(NULL), silent = TRUE)
  }, add = TRUE)

  proposal_cov_method <- .dynhr_opt("proposal_cov_method", proposal_cov_method,
                                     default = "diagonal")
  proposal_cov_method <- match.arg(proposal_cov_method, c("diagonal", "full"))
  transform_params <- isTRUE(.dynhr_opt("transform_params", transform_params,
                                         default = TRUE))
  use_exact_hessian <- isTRUE(.dynhr_opt("use_exact_hessian", use_exact_hessian,
                                          default = FALSE))
  .vcat <- function(...) if (verbose) cat(...)

  model    <- solved$model
  compiled <- solved$compiled

  ## Default obs_vars from the model's varobs declaration (same rule as
  ## make_log_posterior) — pathological-DSGE paper gap #4.
  if (is.null(obs_vars) || length(obs_vars) == 0L) {
    obs_vars <- model$obs_vars %||% model$varobs_names
    if (is.null(obs_vars) || length(obs_vars) == 0L)
      stop("run_mode_finding: `obs_vars` not supplied and the model ",
           "declares no `varobs`; pass obs_vars= or add a varobs line to ",
           "the .mod.", call. = FALSE)
  }

  ## Apply unified plan= if supplied.  Error if both plan and individual args.
  ## Tier 8 item 10: plan adaptation now routes through estimation_context()
  ## so that the compiled specs live in ctx$plan for provenance.
  if (!is.null(plan)) {
    if (!inherits(plan, "dynhr_plan"))
      stop("run_mode_finding: 'plan' must be a dynhr_plan object.", call. = FALSE)
    if (!is.null(filter_tunes))
      stop("run_mode_finding: supply either 'plan' or 'filter_tunes', not both.",
           call. = FALSE)
    if (!is.null(heteroskedastic_shocks))
      stop("run_mode_finding: supply either 'plan' or 'heteroskedastic_shocks', not both.",
           call. = FALSE)
    ## Compile plan inside estimation_context (TPF check + spec resolution).
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
  model <- .resolve_stochastic_volatility(model, stochastic_volatility)

  .vcat("\n================================================================\n")
  .vcat("  run_mode_finding\n")
  .vcat("================================================================\n\n")

  # -------------------------------------------------------------------
  # Step 1: Load data (if path provided)
  # -------------------------------------------------------------------
  .vcat("-- Step 1: Load data --\n")
  if (is.character(data)) {
    data_raw <- read.csv(data)
    if (!is.null(data_col_map)) {
      for (mod_nm in names(data_col_map)) {
        dat_nm <- data_col_map[[mod_nm]]
        if (dat_nm %in% names(data_raw) && !(mod_nm %in% names(data_raw)))
          names(data_raw)[names(data_raw) == dat_nm] <- mod_nm
      }
    }
    if (!all(obs_vars %in% names(data_raw)))
      stop(sprintf("obs_vars not found in data: %s",
                   paste(setdiff(obs_vars, names(data_raw)), collapse = ", ")))
    data <- as.matrix(data_raw[, obs_vars])
  }
  if (is.null(colnames(data))) colnames(data) <- obs_vars
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
  priors <- extract_prior_spec(model, verbose = verbose)
  .vcat(sprintf("  %d parameters to estimate\n", nrow(priors)))

  # -------------------------------------------------------------------
  # Step 3: Build log-posterior closure
  # -------------------------------------------------------------------
  .vcat("-- Step 3: Build log-posterior --\n")

  # Auto-detect OBC: any equation carries an MCP tag?
  use_obc <- any(vapply(model$equations,
                         function(e) {
                           tg <- e$tag
                           !is.null(tg) && !is.na(tg) &&
                             grepl("^\\s*mcp\\s*=", tg, perl = TRUE)
                         }, logical(1L)))

  obc_specs <- NULL
  if (use_obc) {
    ## OBC models can only be evaluated through the OBC/PKF (or PPF/COPF,
    ## chosen later via posterior_options$obc_filter) log-posterior -- the
    ## standard Gaussian/cumulant/whittle/pruned closures do not model the
    ## occasionally-binding constraint. If the caller left `likelihood` at
    ## its default, silently switching to PKF is a reasonable convenience,
    ## but it must be ANNOUNCED (not silent) per this package's fail-loud-on
    ## -fidelity-downgrade convention. If the caller EXPLICITLY asked for an
    ## incompatible likelihood, that is very likely a mistake (they think
    ## they are estimating e.g. "cumulant" but PKF is silently substituted
    ## instead) -- stop() with an actionable message rather than overriding
    ## their choice without telling them.
    if (likelihood_explicit) {
      stop(sprintf(paste0(
        "run_mode_finding: model has OBC (mcp=) tagged equations, which are ",
        "only supported by the OBC/PKF log-posterior -- but likelihood = ",
        "\"%s\" was explicitly requested. These are incompatible: OBC models ",
        "must be estimated via the occasionally-binding-constraint filter ",
        "(PKF by default; PPF/COPF via posterior_options$obc_filter). Either ",
        "drop the `likelihood` argument (default) to use PKF automatically, ",
        "or remove the mcp= tags from the model if you intend to estimate it ",
        "as a standard linear model."),
        likelihood), call. = FALSE)
    }
    .vcat("  OBC model detected -- using PKF log-posterior\n")
    message("run_mode_finding: model has OBC (mcp=) tagged equations; ",
            "auto-switching likelihood from the default \"", likelihood,
            "\" to the OBC/PKF log-posterior (standard Kalman-filter ",
            "likelihoods cannot represent occasionally-binding constraints). ",
            "Pass posterior_options = list(obc_filter = \"ppf\"|\"copf\") to ",
            "use a particle filter instead.")
    obc_specs   <- obc_parse_tags(model)
    max_inner   <- posterior_options$max_inner %||% 10L
    log_post_fn <- make_log_posterior_obc_pkf(
      model, data, priors, obs_vars, compiled,
      specs       = obc_specs,
      me_variance = me_variance,
      max_inner   = max_inner
    )
  } else {
    ## CPM routing (v1 limitation, documented on make_log_posterior_tpf's
    ## `burn_in_init` @param): when posterior_options$tpf_options$cpm_rho_u is
    ## set, run_posterior_estimation() will later dispatch this closure to
    ## rwmh_cpm(), which calls it with a non-NULL U_list on every step after
    ## the priming call. burn_in_init > 0 hard-errors against a supplied
    ## U_list (the CPM slot layout has no burn-in slots), so this is the
    ## build site that must pin burn_in_init = 0L for that case -- unless the
    ## caller already passed burn_in_init explicitly via `...`.
    .mf_dots <- list(...)
    .cpm_rho_u <- posterior_options$tpf_options$cpm_rho_u %||% NULL
    if (identical(likelihood, "tpf") && !is.null(.cpm_rho_u) &&
        is.null(.mf_dots$burn_in_init)) {
      .mf_dots$burn_in_init <- 0L
    }
    log_post_fn <- do.call(make_log_posterior,
      c(list(model, data, priors, obs_vars, compiled,
             me_variance        = me_variance,
             likelihood         = likelihood,
             pruned_order       = pruned_order,
             lik_init           = posterior_options$lik_init %||% "auto",
             me_extra           = me_extra,
             shock_scale        = shock_scale_mat,
             freq_band          = freq_band,
             system_priors      = posterior_options$system_priors %||% NULL,
             infeasible_penalty = posterior_options$infeasible_penalty %||% NULL),
        .mf_dots))
  }

  # -------------------------------------------------------------------
  # Step 4: Starting values from prior means
  # -------------------------------------------------------------------
  .vcat("-- Step 4: Initialise from prior means --\n")
  theta_init <- setNames(priors$mean, priors$name)

  # Evaluate at calibration
  lp_cal <- log_post_fn(theta_init)
  .vcat(sprintf("  Log-posterior at prior means: %.4f\n", lp_cal$logpost %||% -Inf))

  # -------------------------------------------------------------------
  # Step 5: Mode-finding
  # -------------------------------------------------------------------
  .vcat(sprintf("-- Step 5: Mode-finding (method = '%s', max_iter = %d) --\n",
                method, n_iter))

  mo <- modifyList(list(), mode_options)

  ## Build a temporary ctx to use .ctx_is_standard_gaussian (avoids duplication).
  ## student_t needs student_df threaded through here too (estimation_context()
  ## validates it at construction) -- pulled from `...` the same way as the
  ## final mode_ctx below.
  .tmp_dots <- list(...)
  .tmp_mode_ctx <- estimation_context(
    me_variance  = me_variance,
    likelihood   = likelihood,
    me_extra     = me_extra,
    shock_scale  = shock_scale_mat,
    freq_band    = freq_band,
    pruned_order = pruned_order,
    student_df   = .tmp_dots$student_df %||% NULL
  )
  rm(.tmp_dots)
  par_standard <- .ctx_is_standard_gaussian(.tmp_mode_ctx, use_obc = use_obc)
  rm(.tmp_mode_ctx)
  ## Decouple the parallel MULTI-START (Step 5) from the parallel HESSIAN
  ## (Step 6, use_par_hess below -- independent).
  ##
  ## The parallel multi-start re-initialises the model per daemon and CANNOT
  ## carry newrat's analytic H0 seed (posterior_hessian, ~216 s -- absurd to
  ## rebuild on every daemon). A seed-less parallel newrat chain is ~10x slower
  ## than the serial newrat (measured: 5687 s vs 570 s on NZSIM) and, on a
  ## razor-thin near-unit-root surface, its eta-space/FD steps stall at the start
  ## -- so parallel = TRUE used to return the UNOPTIMISED prior mean for the
  ## (default) newrat method. For these H0-seed-dependent methods, run Step 5
  ## SERIALLY (with the seed + analytic gradient) and parallelise ONLY the Step 6
  ## Hessian. The multi-start parallelism still applies to the gradient-free /
  ## global methods (cmaes, jade, nmkb, combined, ...), which need no seed.
  ## UPDATE: the H0 seed is now computed ONCE on the host (a parallel FD Hessian
  ## at the prior mean, ~9 s) and shipped to the daemons as a cheap matrix, so the
  ## "seed not replicable per daemon -> serial" constraint no longer holds. For
  ## seeded methods + parallel we now run a SEEDED PORTFOLIO (run_mode_mirai):
  ## chain 1 = newrat @ theta_init with the seed = bit-identical to the serial
  ## run; chains 2+ cycle dispersed starts and the global searchers, all newrat
  ## chains sharing the seed. Same wall time as the serial newrat (chains run
  ## concurrently) but far more robust (multi-start + multi-algo, keep-best).
  seed_dependent_method <- method %in% c("newrat", "cmaes_newrat")
  use_portfolio <- isTRUE(parallel) && seed_dependent_method && par_standard &&
    requireNamespace("mirai", quietly = TRUE)
  use_par_mode <- isTRUE(parallel) && !seed_dependent_method &&
    requireNamespace("mirai", quietly = TRUE)
  if (isTRUE(parallel) && seed_dependent_method && !par_standard)
    .vcat(sprintf(paste0("  [parallel] method '%s' with a non-standard likelihood: ",
                         "running serially + parallel Hessian only (Step 6).\n"),
                  method))

  # Opt-in unconstrained-parameter transform (eta-space mode finding). Built
  # from `priors` (the prior_spec already in scope) -- see
  # build_param_transform() and .run_mode_finding()'s `transform` argument.
  # Applies to both the serial and parallel (mirai multi-start) paths.
  mode_transform <- NULL
  if (transform_params) {
    mode_transform <- build_param_transform(priors, names(theta_init))
  }

  if (use_portfolio) {
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
                       lik_init = "stationary"),
      error = function(e) {
        .vcat(sprintf("  (shared daemon pool unavailable: %s; each stage self-provisions)\n",
                      conditionMessage(e)))
        try(mirai::daemons(NULL), silent = TRUE); NULL
      })
    pool_ok <- !is.null(pool_sh)
    ## Keep this pool ALIVE for Step 6 (build_sigma_prop re-binds .worker_lp on
    ## it). Torn down once by run_mode_finding's on.exit. We deliberately do NOT
    ## tear it down here.
    shared_pool_ok <- pool_ok
    par_mode <- tryCatch({
      H0_seed <- if (pool_ok) tryCatch({
          H <- -num_hessian_mirai(NULL, theta_init, h = 1e-4, n_cores = ncs,
                                  verbose = FALSE, pool_ready = TRUE)
          if (any(!is.finite(H))) NULL else H
        }, error = function(e) {
          .vcat(sprintf("  (shared H0 seed unavailable: %s; chains run unseeded)\n",
                        conditionMessage(e))); NULL
        }) else NULL
      if (!is.null(H0_seed)) .vcat("  Shared H0 seed built; launching portfolio.\n")
      ## The portfolio runs in the SAME space as the serial path (eta-space when
      ## transform_params = TRUE). The earlier "eta breaks on the daemon" was a
      ## double chain-rule bug in mode_task (it pre-transformed grad_fn AND passed
      ## transform, so .run_mode_finding chain-ruled twice -> newrat stuck at the
      ## prior mean); now fixed (mode_task passes the THETA-space gradient and lets
      ## .run_mode_finding apply the chain rule, exactly like the host path). The
      ## H0 seed is theta-space; .run_mode_finding chain-rules it to eta too. Running
      ## in eta-space recovers ~248 nats of mode quality on NZSIM (9060 vs 8812).
      run_mode_mirai(
        parsed_model   = model, Y = data, prior_spec = priors, obs_names = obs_vars,
        theta_init     = theta_init, n_chains = n_starts, nm_maxit = n_iter,
        method         = method, n_cores = if (pool_ok) ncs else n_cores,
        me_variance    = me_variance,
        me_extra       = me_extra, shock_scale = shock_scale_mat,
        transform      = mode_transform,
        analytic_grad  = isTRUE(mo$use_analytic_grad %||% TRUE),
        likelihood     = likelihood, freq_band = freq_band,
        newrat_H0_seed = H0_seed, lik_init = "auto", pool_ready = pool_ok,
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
      me_variance  = me_variance,
      me_extra     = me_extra,
      shock_scale  = shock_scale_mat,
      log_post_fn  = if (par_standard) NULL else log_post_fn,
      transform    = mode_transform,
      ## Build the analytic gradient on each daemon so the parallel newrat /
      ## cmaes_newrat chains converge (FD-only newrat is far too slow). Only the
      ## standard recompile-per-daemon path can build it. Opt out via
      ## mode_options$use_analytic_grad = FALSE.
      analytic_grad = par_standard && isTRUE(mo$use_analytic_grad %||% TRUE),
      likelihood   = likelihood,
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
      grad_fn <- tryCatch(
        make_posterior_grad(
          model, data, priors, obs_vars, compiled,
          me_variance = me_variance, me_extra = me_extra,
          shock_scale = shock_scale_mat, likelihood = likelihood,
          freq_band = freq_band, verbose = FALSE),
        error = function(e) {
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

    ## Build a theta-space analytic posterior Hessian supplier for the newrat
    ## initial H0. Only available for standard Gaussian KF (same gate as
    ## use_exact_hessian). Falls back to csminwel default (1e-4 * I) if
    ## unavailable/fails -- .run_mode_finding handles the NULL case.
    hessian_fn <- NULL
    use_analytic_hess <- isTRUE(mo$use_analytic_hess %||% TRUE)
    h0_method <- mo$h0_method %||% "auto"   # "auto" | "numderiv" | "analytic"
    ## Shared eligibility for an analytic-quality H0 seed (standard Gaussian KF).
    h0_eligible <- method %in% c("newrat", "cmaes_newrat") && !use_obc &&
      identical(likelihood, "gaussian") &&
      is.null(me_extra) && is.null(shock_scale_mat) && !anyNA(data)

    ## Option B: when a parallel pool is available, seed newrat's H0 from a
    ## PARALLEL finite-difference Hessian at theta_init (the prior mean -- an
    ## interior point where the FD Hessian is finite) instead of the analytic
    ## posterior_hessian. ~20x faster (NZSIM: 205s -> 9s on 16 cores) and an
    ## equivalent seed: csminwel BFGS-refines H0, and the two H0 matrices agree
    ## to <1% Frobenius. mode-orchestrate .make_pd-regularises + inverts the
    ## result and falls back to csminwel's default 1e-4*I on any failure, so an
    ## indefinite/approximate seed is safe. "auto" uses it whenever a pool
    ## exists; opt out with mode_options$h0_method = "analytic".
    want_parallel_h0 <- h0_eligible && h0_method %in% c("auto", "numderiv") &&
      isTRUE(parallel) && par_standard && requireNamespace("mirai", quietly = TRUE)
    if (want_parallel_h0) {
      hessian_fn <- local({
        .lp <- log_post_fn; .model <- model; .data <- data; .priors <- priors
        .obs <- obs_vars;  .mev <- me_variance; .nc <- n_cores
        function(theta) {
          ncs <- .mirai_n_cores(.nc, length(theta) * (length(theta) + 1L) %/% 2L)
          sh <- .mirai_pool_init(ncs, .model, .data, .priors, .obs, .mev,
                                 lik_init = "stationary")
          on.exit({ mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }, add = TRUE)
          ## .lp is unused on the pool_ready path (daemons evaluate .worker_lp);
          ## pass NULL so num_hessian_mirai's task closure never serialises this
          ## heavy posterior closure (see the note there).
          -num_hessian_mirai(NULL, theta, h = 1e-4, n_cores = ncs,
                             verbose = FALSE, pool_ready = TRUE)
        }
      })
      .vcat("  Will build newrat initial H0 from a parallel finite-difference Hessian (option B)\n")
    } else if (use_analytic_hess && h0_eligible) {
      hessian_fn <- tryCatch({
        local({
          .model    <- model
          .compiled <- compiled
          .priors   <- priors
          .obs_vars <- obs_vars
          .data     <- data
          .me_var   <- me_variance
          function(theta) {
            pm <- .model$param_values
            pm[names(theta)] <- theta
            ss_h <- tryCatch(
              solve_steady_state(.model, .compiled, pm, verbose = FALSE),
              error = function(e) NULL)
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
              check_mode = FALSE)
            ## Convert logpost Hessian -> neg-logpost Hessian (flip sign).
            -H_logpost
          }
        })
      }, error = function(e) {
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
    dots_optim <- list(...)
    dots_optim[c("student_df", "cut_tol")] <- NULL
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
        grad_fn     = grad_fn,
        hessian_fn  = hessian_fn,
        verbose     = verbose,
        record_curvature = isTRUE(mo$record_curvature %||% FALSE)
      ),
      dots_optim
    ))
  }

  if (is.null(mode_res) || !is.finite(mode_res$logpost))
    stop("Mode-finding failed or returned non-finite log-posterior.")

  theta_mode <- mode_res$theta_mode
  .vcat(sprintf("  Log-posterior at mode: %.4f\n", mode_res$logpost))

  # -------------------------------------------------------------------
  # Step 6: Proposal covariance (via Hessian or prior std)
  # -------------------------------------------------------------------
  .vcat("-- Step 6: Build proposal covariance --\n")

  ## Tier 11 #2 (curvature-aware estimation): when use_exact_hessian = TRUE and
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
    !identical(posterior_options$lik_init %||% "auto", "diffuse")
  if (isTRUE(use_exact_hessian) && .allows_exact_hess) {
    hessian_exact <- tryCatch({
      pm <- model$param_values
      pm[names(theta_mode)] <- theta_mode
      ss_m <- solve_steady_state(model, compiled, pm, verbose = FALSE)
      if (!isTRUE(ss_m$converged))
        stop("steady state did not converge at the mode")
      ## Re-derive SSM-computed params so the exact Hessian uses the consistent
      ## (not stale) p_c (no-op for non-SSM-parameter models; Tier 13 #1).
      pm <- ss_m$params %||% pm
      dr_m <- solve_perturbation(model, compiled, ss_m$ss, pm, verbose = FALSE)
      H <- posterior_hessian(model, compiled, dr_m, pm, names(theta_mode),
                             obs_vars, t(data), me_variance = me_variance,
                             include_prior = TRUE, prior_spec = priors)
      if (any(!is.finite(H))) stop("non-finite entries in exact Hessian")
      .vcat("  Computed exact posterior Hessian at mode (analytic adjoint).\n")
      H
    }, error = function(e) {
      .vcat(sprintf("  Exact Hessian unavailable (%s); using numerical.\n",
                    conditionMessage(e)))
      NULL
    })
  } else if (isTRUE(use_exact_hessian)) {
    .vcat("  use_exact_hessian = TRUE but the model is not standard-Gaussian; using numerical.\n")
  }

  if (identical(proposal_cov_method, "full")) {
    .vcat("  Using full-Hessian proposal covariance (method = 'full').\n")
    Sigma_prop <- proposal_cov(
      lp_fn      = log_post_fn,
      theta_mode = theta_mode,
      prior_spec = priors,
      verbose    = verbose,
      method     = "full"
    )
  } else {
    Sigma_prop <- build_sigma_prop(
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
      ## Re-bind .worker_lp on the live portfolio pool instead of recompiling the
      ## model on every daemon again (Step 5 -> Step 6 pool sharing).
      pool_ready   = shared_pool_ok
    )
  }

  # -------------------------------------------------------------------
  # Wrap up
  # -------------------------------------------------------------------
  .vcat("\n================================================================\n")
  .vcat("  Mode-finding complete\n")
  .vcat("================================================================\n\n")

  ## Build estimation context to store in the result.
  ## This fixes the system_priors and freq_band bugs in run_posterior_estimation:
  ## the ctx carries all options so they can be read back without loss.
  ## Tier 8 item 10: store plan provenance in mode_ctx$plan.
  ## When use_obc = TRUE, stamp the OBC filter type that was actually used so
  ## that conditional_forecast() dispatches to the correct terminal-state path
  ## (Tier 10 item 5; Tier 15 §C: ppf/copf stamping).
  ## run_mode_finding always calls make_log_posterior_obc_pkf; mode_ctx$likelihood
  ## is "pkf".  Callers using PPF/COPF can override via posterior_options$obc_filter.
  mode_ctx_likelihood <- if (use_obc) {
    ft <- posterior_options$obc_filter %||% "pkf"
    match.arg(as.character(ft), c("pkf", "ppf", "copf"))
  } else likelihood
  ## student_t needs student_df threaded through to estimation_context() too
  ## (it validates student_df at construction); it is passed to
  ## make_log_posterior() above via `...`, picked back up here the same way.
  dots_ <- list(...)
  mode_ctx <- estimation_context(
    me_variance   = me_variance,
    likelihood    = mode_ctx_likelihood,
    lik_init      = posterior_options$lik_init %||% "auto",
    me_extra      = me_extra,
    shock_scale   = shock_scale_mat,
    freq_band     = freq_band,
    system_priors = posterior_options$system_priors %||% NULL,
    tpf_options   = posterior_options$tpf_options %||% list(),
    obc_specs     = obc_specs,
    student_df    = dots_$student_df %||% NULL,
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
    V_mode         = mode_res$V_mode %||% NULL,
    hessian_exact  = hessian_exact,
    log_marglik_laplace = NA_real_,
    me_variance    = me_variance,
    me_extra       = me_extra,
    shock_scale    = shock_scale_mat,
    likelihood     = likelihood,
    freq_band      = freq_band,
    system_priors  = posterior_options$system_priors %||% NULL,
    obc_specs      = obc_specs,
    ctx            = mode_ctx,
    meta           = list(
      method     = method,
      n_iter     = n_iter,
      use_obc    = use_obc
    )
  )
  ## Laplace marginal likelihood (Tier 11 #2): free once the exact posterior
  ## Hessian is available. NA when the Hessian is absent or -H is not PD.
  if (!is.null(hessian_exact))
    result$log_marglik_laplace <- laplace_log_marglik(result)

  class(result) <- c("dynhr_mode_result", "list")
  invisible(result)
}


#' Print method for dynhr_mode_result
#' @noRd
#' @export
print.dynhr_mode_result <- function(x, ...) {
  cat("\n<dynhr_mode_result>\n")
  cat(sprintf("  Estimated params: %d\n", nrow(x$prior_spec)))
  cat(sprintf("  Observed vars  : %s\n", paste(x$obs_vars, collapse = ", ")))
  cat(sprintf("  Log-posterior  : %.4f\n", x$mode$logpost))
  cat(sprintf("  Method         : %s\n", x$meta$method))
  if (!is.null(x$obc_specs))
    cat(sprintf("  OBC constraints: %d\n", length(x$obc_specs)))
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
                              pool_ready = FALSE) {
  n_par <- length(theta_mode)

  .vcat <- function(...) if (verbose) cat(...)

  # Priority 1: empirical covariance from pooled draws
  if (!is.null(pooled_draws) && nrow(pooled_draws) > n_par) {
    .vcat("  Using empirical covariance from pooled draws.\n")
    Sigma <- cov(pooled_draws)
    opt_scale <- 2.38^2 / n_par
    rownames(Sigma) <- colnames(Sigma) <- names(theta_mode)
    return(Sigma * opt_scale)
  }

  # Priority 2: Hessian at mode. When an exact posterior Hessian (analytic
  # adjoint, log-posterior i.e. include_prior = TRUE) is supplied it is used
  # directly in place of the finite-difference Hessian; all downstream repair
  # (ill-conditioning fallback, prior-scale cap, PD regularisation) is identical.
  if (!is.null(hess_exact)) {
    .vcat("  Using exact posterior Hessian at mode (analytic adjoint).\n")
    hess <- hess_exact
  } else {
    h <- max(1e-4, 1e-4 * abs(theta_mode))
    use_par_hess <- isTRUE(parallel) && n_par > 1L &&
                    requireNamespace("mirai", quietly = TRUE)
    if (use_par_hess && par_standard) {
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
      n_cores_h <- .mirai_n_cores(n_cores, n_par * (n_par + 1L) %/% 2L)
      host_lp_mode <- lp_fn(theta_mode)$logpost %||% -Inf
      ## Set the daemon filter: on a SHARED pool (pool_ready, caller-owned) just
      ## re-bind .worker_lp -- no model recompile; otherwise stand up a fresh
      ## pool here. Returns the share handle (NULL on the rebind path).
      .set_filter <- function(lik) {
        if (isTRUE(pool_ready)) {
          .mirai_rebind_worker_lp(prior_spec, obs_vars, me_variance = me_variance,
                                  me_extra = me_extra, shock_scale = shock_scale,
                                  lik_init = lik)
          NULL
        } else {
          .mirai_pool_init(n_cores_h, parsed_model, data, prior_spec, obs_vars,
                           me_variance, me_extra = me_extra,
                           shock_scale = shock_scale, lik_init = lik)
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
      }, error = function(e) NA_real_)
      if (!is.finite(stat_lp_mode) ||
            abs(stat_lp_mode - host_lp_mode) > 1e-4 * max(1, abs(host_lp_mode))) {
        .vcat("  Stationary filter disagrees at the mode; re-initialising with auto (diffuse).\n")
        ## Only tear down a pool we OWN; a shared pool is just re-bound.
        if (!isTRUE(pool_ready)) { mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }
        sh <- .set_filter("auto")
        fit_filter <- "auto"
      }
      .vcat(sprintf("    (filter = %s, %d daemons%s)\n", fit_filter, n_cores_h,
                    if (isTRUE(pool_ready)) ", shared pool" else ""))
      hess <- tryCatch(
        num_hessian_mirai(lp_fn, theta_mode, h = h, n_cores = n_cores_h,
                          verbose = verbose, pool_ready = TRUE),
        error = function(e) matrix(NA_real_, n_par, n_par))
      ## Tear down only a pool we created here; the caller owns a shared pool.
      if (!isTRUE(pool_ready)) { mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }
    } else if (use_par_hess) {
      ## Non-standard likelihood (me_extra / shock_scale / non-Gaussian): the
      ## per-daemon re-init path does not cover it, so ship the closure.
      .vcat("  Estimating Hessian at mode (mirai, parallel)...\n")
      hess <- num_hessian_mirai(lp_fn, theta_mode, h = h, n_cores = n_cores,
                                verbose = verbose)
    } else {
      .vcat("  Estimating Hessian at mode (numDeriv)...\n")
      hess <- num_hessian(lp_fn, theta_mode, h = h)
    }
  }
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

  Sigma
}


#' Numerical Hessian via central differences
#' @noRd
num_hessian <- function(fn, theta, h = 1e-4) {
  n <- length(theta)
  H <- matrix(0, nrow = n, ncol = n)
  f0 <- fn(theta)
  f0_val <- if (is.list(f0)) f0$logpost else f0

  for (i in seq_len(n)) {
    for (j in i:n) {
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
#' @return n x n Hessian matrix (unnamed dimnames, as \code{num_hessian}).
#' @noRd
num_hessian_mirai <- function(fn, theta, h = 1e-4, n_cores = NULL,
                              verbose = TRUE, pool_ready = FALSE) {
  n <- length(theta)
  par_names <- names(theta)

  pairs <- vector("list", n * (n + 1L) %/% 2L)
  k <- 0L
  for (i in seq_len(n)) for (j in i:n) {
    k <- k + 1L
    pairs[[k]] <- c(i, j)
  }

  n_cores <- .mirai_n_cores(n_cores, length(pairs))
  if (verbose) cat(sprintf("    %d Hessian evaluations on %d daemons\n",
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

  stencil <- function(pair, .theta, .h, .par_names) {
    lp <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    .eval <- function(th) {
      names(th) <- .par_names
      r <- lp(th)
      if (is.list(r)) r$logpost else r
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
  ## CRITICAL: sever the stencil's environment. It needs nothing from this
  ## frame (it takes .theta/.h/.par_names via .args and fetches the worker
  ## log-posterior via get0(".worker_lp", globalenv()) on the daemon). Left
  ## attached, mirai_map serialises this frame -- including `fn`, the passed
  ## log-posterior -- with EVERY one of the n(n+1)/2 tasks. For a real
  ## make_posterior() closure (heavy adjoint/model state) that is ~50x slower
  ## per eval inside a full run_mode_finding than standalone (NZSIM: 270s vs
  ## 5s). globalenv() ships nothing.
  environment(stencil) <- globalenv()

  raw <- mirai::mirai_map(
    pairs, stencil,
    .args = list(.theta = theta, .h = h, .par_names = par_names)
  )[]

  H <- matrix(0, nrow = n, ncol = n)
  for (k in seq_along(pairs)) {
    p <- pairs[[k]]
    v <- raw[[k]]
    if (inherits(v, "miraiError") || inherits(v, "errorValue")) v <- NA_real_
    H[p[1L], p[2L]] <- H[p[2L], p[1L]] <- v
  }

  if (verbose) cat(sprintf("    Hessian: %.1f sec\n",
                           (proc.time() - t0)[["elapsed"]]))
  H
}
