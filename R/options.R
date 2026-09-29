## R/options.R
## Package-level option store for dynhr.
##
## dynhr_set_options() stores named values that override function argument
## defaults. Useful when one option (e.g. me_variance) must flow through many
## nested calls without threading it through every signature.
##
## Priority: explicit argument > global option > registry default.

.dynhr_opts <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------------
# The option registry: the ONE place an option's default lives.
# ---------------------------------------------------------------------------
# Every name read through `.dynhr_opt("<name>")` anywhere in R/ must have an
# entry here (guarded by tests/testthat/test-fix-0926-run-record.R, which
# parses R/). Each entry carries
#   default          the value used when neither an argument nor
#                    dynhr_set_options() supplies one;
#   changes_results  TRUE when a non-default value can change a numerical
#                    result (a posterior, a mode, a draw) rather than only
#                    diagnostics, guards or I/O;
#   description      one line, for dynhr_get_options() and the run record.
# `.dynhr_opt(name)` without a `default =` falls back to the registry default.
# Call sites that still pass `default =` must pass exactly this value (also
# guarded). The run record (R/run-record.R) snapshots EVERY registered option
# via dynhr_get_options(effective = TRUE), so a registered option is what makes
# a global setting replayable -- register new options here first.
# Keep in sync with the \describe{} list in dynhr_set_options()'s roxygen.
.dynhr_option_registry <- list(
  me_variance = list(
    default = 0, changes_results = TRUE,
    description = "Measurement-error variance for the PSOCK/future parallel mode and MCMC helpers"),
  perturb_scale = list(
    default = 0.5, changes_results = TRUE,
    description = "Multi-start perturbation scale (fraction of prior sd) for run_mode_parallel"),
  nm_maxit = list(
    default = 5000L, changes_results = TRUE,
    description = "Per-chain iteration budget for run_mode_parallel"),
  mode_method = list(
    default = "cmaes_nmkb", changes_results = TRUE,
    description = "Optimiser for run_mode_parallel (NOT run_mode_finding's method)"),
  seed_base = list(
    default = 42L, changes_results = TRUE,
    description = "Base RNG seed for run_mode_parallel and the parallel (mirai) mode starts and RWMH / NUTS chains when no seed is passed"),
  proposal_cov_method = list(
    default = "diagonal", changes_results = TRUE,
    description = "run_mode_finding Step 6 proposal covariance: diagonal or full"),
  transform_params = list(
    default = TRUE, changes_results = TRUE,
    description = "Mode finding and samplers work in unconstrained eta-space"),
  use_exact_hessian = list(
    default = FALSE, changes_results = TRUE,
    description = "run_mode_finding: RWMH proposal from the analytic posterior Hessian"),
  use_analytic_hess = list(
    default = TRUE, changes_results = TRUE,
    description = "run_mode_finding newrat: seed the initial curvature H0 from the analytic posterior Hessian (standard Gaussian Kalman filter)"),
  rwmh_adapt_cov = list(
    default = FALSE, changes_results = TRUE,
    description = "run_posterior_estimation RWMH: Haario adaptive covariance"),
  rwmh_n_blocks = list(
    default = 1L, changes_results = TRUE,
    description = "run_posterior_estimation RWMH: number of randomized parameter blocks"),
  grad_method = list(
    default = "auto", changes_results = TRUE,
    description = "Analytic posterior gradient method for NUTS/MALA/ChEES: auto, hybrid, implicit, adjoint or adjoint_solution"),
  nuts_adapt = list(
    default = "independent", changes_results = TRUE,
    description = "Parallel multi-chain NUTS warmup: independent or pooled"),
  monge_alpha = list(
    default = 1, changes_results = TRUE,
    description = "Monge-metric MALA alpha when the monge_alpha argument is NULL"),
  power_posterior = list(
    default = 1, changes_results = TRUE,
    description = "Power-posterior tempering exponent zeta in (0, 1]"),
  pskf_cdf = list(
    default = "accurate", changes_results = TRUE,
    description = "PSKF multivariate-normal CDF evaluation: accurate (deterministic C++ tilted separation-of-variables lattice, log-scale error below max(1e-5, 1e-7|log p|); compensation up to dim 5) or fast (pre-0.9.4 Mendell-Elston compensation, plain Miwa)"),
  allow_monge_metric = list(
    default = FALSE, changes_results = FALSE,
    description = "Permit the experimental metric = \"monge\" MALA option"),
  allow_whittle_fim_metric = list(
    default = FALSE, changes_results = FALSE,
    description = "Permit the experimental metric = \"whittle_fim\" NUTS option"),
  checkpoint_flush_every = list(
    default = 1000L, changes_results = FALSE,
    description = "Checkpointed RWMH flush chunk size (draws)"),
  debug_kf_errors = list(
    default = FALSE, changes_results = FALSE,
    description = "Re-raise Kalman-filter errors inside the log-posterior instead of loglik = -Inf"),
  verbosity = list(
    default = "info", changes_results = FALSE,
    description = "Message verbosity level (see dynhr_set_verbosity)")
)

#' Set dynhr package-level options
#'
#' Named values stored here act as defaults for any function that calls
#' \code{.dynhr_opt(name)}. Explicit argument values always take precedence:
#' the priority is explicit argument, then the value set here, then the
#' registered default listed below.
#'
#' Every option dynhr reads is registered in one internal table, which holds
#' its default and whether it can change a numerical result.
#' \code{dynhr_get_options(effective = TRUE)} returns all of them with their
#' effective values, and the run record attached to every estimation result
#' (see \code{\link{dynhr_rerun}}) stores that snapshot so a run made under
#' global options can be replayed exactly. Setting a name that is not in the
#' list below is allowed but has no effect.
#'
#' Registered options (default in parentheses):
#' \describe{
#'   \item{\code{me_variance}}{(\code{0}) Measurement-error variance used by the
#'     internal PSOCK/future parallel helpers \code{run_mode_parallel()} and
#'     \code{run_mcmc_parallel()} when their \code{me_variance} argument is
#'     \code{NULL}. It does \emph{not} reach \code{make_log_posterior()},
#'     \code{run_mode_finding()} or \code{run_full_estimation()}, which take
#'     \code{me_variance} as an argument.}
#'   \item{\code{perturb_scale}}{(\code{0.5}) Starting-point perturbation scale
#'     (fraction of the prior standard deviation) for the multi-start
#'     \code{run_mode_parallel()}.}
#'   \item{\code{nm_maxit}}{(\code{5000L}) Iteration budget for each
#'     \code{run_mode_parallel()} chain.}
#'   \item{\code{mode_method}}{(\code{"cmaes_nmkb"}) Optimiser used by
#'     \code{run_mode_parallel()} when its \code{method} argument is
#'     \code{NULL}. Not read by \code{run_mode_finding()} (whose \code{method}
#'     defaults to \code{"newrat"}) nor by \code{run_full_estimation()}
#'     (\code{mode_method}, default \code{"newrat"}).}
#'   \item{\code{seed_base}}{(\code{42L}) Base RNG seed for the multi-start
#'     \code{run_mode_parallel()}, and for the parallel (mirai) mode
#'     starts and RWMH / NUTS chains when no seed is passed (chain k seeds
#'     with \code{seed_base + k}); a spec's \code{compute$seed} sets it for
#'     the run.}
#'   \item{\code{use_exact_hessian}}{(\code{FALSE}) \code{run_mode_finding()}
#'     builds the RWMH proposal covariance from the analytic exact posterior
#'     Hessian -- see that function's \code{use_exact_hessian} argument.}
#'   \item{\code{use_analytic_hess}}{(\code{TRUE}) \code{run_mode_finding()}
#'     with \code{method = "newrat"} / \code{"cmaes_newrat"} on the standard
#'     Gaussian Kalman likelihood seeds the optimiser's initial curvature
#'     \code{H0} from the analytic posterior Hessian (it does not set the
#'     reported covariance). \code{FALSE} uses csminwel's default \code{H0}.
#'     The spec field \code{mode$analytic_h0}; a
#'     \code{mode_options$use_analytic_hess} entry overrides it for one run.}
#'   \item{\code{monge_alpha}}{(\code{1}) Monge-metric scale used by the
#'     \code{metric = "monge"} MALA option of
#'     \code{run_posterior_estimation()} when its \code{monge_alpha} argument
#'     is \code{NULL}.}
#'   \item{\code{verbosity}}{(\code{"info"}) dynhr message level; set it with
#'     \code{\link{dynhr_set_verbosity}}.}
#'   \item{\code{proposal_cov_method}}{\code{"diagonal"} (default) or
#'     \code{"full"}; selects the Step 6 proposal-covariance strategy in
#'     \code{run_mode_finding()} -- see \code{proposal_cov}.}
#'   \item{\code{transform_params}}{\code{TRUE} (default); mode finding
#'     (\code{run_mode_finding()}) and the samplers
#'     (\code{run_posterior_estimation()}: RWMH, NUTS, MALA, ChEES, HMC) operate
#'     on an unconstrained reparameterisation built by
#'     \code{build_param_transform} instead of the raw, box-constrained
#'     parameter vector. Default \code{TRUE} because the gradient samplers
#'     run into the bounds of parameters near a unit root: on a 68-parameter
#'     NZ DSGE, NUTS without the transform had 246 of 300 transitions
#'     divergent (minimum bulk ESS 1.3) versus 0 divergent (minimum bulk ESS
#'     57) with it -- see \code{\link{run_posterior_estimation}}. The
#'     transform leaves the target invariant. Set \code{FALSE} to sample in raw
#'     theta-space -- see \code{\link{run_mode_finding}}.}
#'   \item{\code{rwmh_adapt_cov}}{\code{FALSE} (default); when \code{TRUE}, RWMH
#'     chains in \code{run_posterior_estimation()} use Haario et al. (2001)
#'     adaptive proposal covariance -- see \code{rwmh}'s \code{adapt_cov}.}
#'   \item{\code{rwmh_n_blocks}}{\code{1L} (default); number of randomized
#'     parameter blocks for RWMH chains in
#'     \code{run_posterior_estimation()} -- see \code{rwmh}'s
#'     \code{n_blocks}.}
#'   \item{\code{grad_method}}{\code{"auto"} (default), \code{"hybrid"},
#'     \code{"implicit"}, \code{"adjoint"} or \code{"adjoint_solution"};
#'     selects the \code{grad_method} passed to
#'     \code{\link{make_posterior_grad}} when \code{analytic_grad = TRUE} (the
#'     default) for the NUTS, HMC, MALA and ChEES samplers of
#'     \code{run_posterior_estimation()} / \code{run_full_estimation()} /
#'     \code{run_estimation()}.
#'     \code{"auto"} picks the fastest exact method the likelihood supports
#'     (\code{"adjoint_solution"} for the Gaussian likelihood,
#'     \code{"implicit"} for Whittle / pruned and for cumulant (except
#'     \code{"adjoint_solution"} for cumulant orders that include 3 but not 4),
#'     \code{"hybrid"} only
#'     where no analytic Kalman kernel covers the data); the method it chose
#'     is \code{attr(grad_fn, "grad_method")}. \code{"hybrid"} (the previous
#'     default) finite-differences the decision-rule-moving parameters; the
#'     other three are exact -- see \code{\link{make_posterior_grad}}.}
#'   \item{\code{nuts_adapt}}{\code{"independent"} (default) or
#'     \code{"pooled"}; warmup adaptation on the parallel multi-chain (mirai)
#'     NUTS path of \code{run_posterior_estimation()} /
#'     \code{run_full_estimation()}. \code{"pooled"} warms all chains up
#'     together (Lao 2026): one shared step size and one metric
#'     (\code{"diagonal"}, \code{"fisher_diag"} or \code{"lowrank"}) fitted to
#'     the pooled, per-chain-centred warmup draws, then frozen identically for
#'     every chain. Not available with \code{metric = "warmup_dense"} or a
#'     sampler checkpoint.}
#'   \item{\code{power_posterior}}{Scalar in \code{(0, 1]} (default \code{1});
#'     power-posterior (generalised-Bayes) tempering exponent \eqn{\zeta}.
#'     The log-posterior returned by \code{make_log_posterior()} becomes
#'     \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)}: values below 1
#'     down-weight a potentially misspecified likelihood and widen the
#'     posterior (Bissiri, Holmes & Walker 2016, JRSS-B; Grünwald SafeBayes).
#'     The \code{$loglik} field always carries the \emph{raw} (untempered)
#'     log-likelihood; only \code{$logpost} is tempered. Default \code{1} is
#'     bit-identical to the standard Bayesian posterior and is compatible with
#'     every existing sampler without any changes.}
#'   \item{\code{pskf_cdf}}{\code{"accurate"} (default) or \code{"fast"};
#'     how the skew-normal (PSKF) likelihood, \code{likelihood = "pskf"} and
#'     \code{make_log_posterior_pskf_order2()}, evaluates its multivariate-normal
#'     CDFs. \code{"accurate"} evaluates the pruning compensation with the
#'     deterministic CDF up to dimension 5, and every 3- to 7-dimensional CDF
#'     with a deterministic C++ separation-of-variables lattice rule (Genz
#'     reordering, minimax tilting, log scale) whose error is below
#'     max(1e-5, 1e-7 |log p|); 2-dimensional CDFs are computed on the log
#'     scale. \code{"fast"} is the evaluation before dynhr 0.9.4
#'     (Mendell-Elston for the compensation, plain Miwa): about 3 to 4 times cheaper under
#'     multi-shock skew and approximate: measured 3e-4 nat off an exact
#'     likelihood with one skewed shock (T = 20), 0.002 nat with two skewed
#'     shocks (T = 100) and 0.3 nat with three (T = 60). (Before dynhr 0.9.4.8
#'     a Mendell-Elston sign error put it 0.14 to 18 nats off.) Resolved once
#'     when the log-posterior is built.}
#'   \item{\code{allow_monge_metric}}{\code{FALSE} (default); the experimental
#'     \code{metric = "monge"} MALA option in \code{run_posterior_estimation()}
#'     errors unless this is \code{TRUE}, because the Monge metric collapses the
#'     proposal step on sharp / near-unit-root posteriors (it under-explored a
#'     tight direction to ~1\% of its variance in testing). It still warns when
#'     enabled; prefer \code{metric = "hessian"} (constant Laplace).}
#'   \item{\code{allow_whittle_fim_metric}}{\code{FALSE} (default); the experimental
#'     \code{metric = "whittle_fim"} NUTS option in
#'     \code{run_posterior_estimation()} errors unless this is \code{TRUE}.
#'     \code{"whittle_fim"} assembles a one-shot dense mass matrix from the
#'     Whittle (frequency-domain) Fisher information at the mode (no periodic
#'     recompute); it is unvalidated on the full estimation pipeline and costs
#'     one extra evaluation that can be expensive on large models. Falls back
#'     to the \code{"hessian"} metric on failure.}
#'   \item{\code{checkpoint_flush_every}}{\code{1000L} (default); the chunk size (in
#'     draws) at which checkpointed RWMH flushes its streaming buffer to disk and
#'     rewrites the restart state, when \code{run_posterior_estimation()} is
#'     called with \code{checkpoint_dir}. Smaller values bound RAM more tightly
#'     and lose less on an interruption, at the cost of more frequent I/O.}
#'   \item{\code{debug_kf_errors}}{\code{FALSE} (default); when \code{TRUE},
#'     a Kalman-filter error inside \code{make_log_posterior()} is RE-RAISED
#'     instead of being caught and treated as an infeasible draw
#'     (\code{loglik = -Inf}). The default keeps an infeasible parameter draw
#'     (singular covariance, unit-root with stationary init, ...) from crashing a
#'     long chain; set this \code{TRUE} when debugging to surface a genuine code
#'     bug that would otherwise be silently masked as a rejected draw.}
#' }
#'
#' The options live in dynhr's own store, NOT in R's \code{options()}:
#' \code{getOption("dynhr.allow_monge_metric")} is always \code{NULL}. Read
#' them with \code{dynhr_get_options()} (the values set) or
#' \code{dynhr_get_options(effective = TRUE)} (every registered option with
#' the value in force).
#'
#' @param ... Named option values. A \code{NULL} value un-sets the option, so
#'   its registered default applies again (as \code{dynhr_reset_options()}).
#' @return Invisibly, a named list of the previous values of the modified
#'   options, \code{NULL} for an option that was not set (its registered
#'   default was in force). \code{do.call(dynhr_set_options, old)} restores
#'   them.
#' @examples
#' old <- dynhr_set_options(allow_monge_metric = TRUE)
#' dynhr_get_options(effective = TRUE)$allow_monge_metric
#' do.call(dynhr_set_options, old)   # back to the registered default
#' dynhr_get_options(effective = TRUE)$allow_monge_metric
#' @export
dynhr_set_options <- function(...) {
  args <- list(...)
  if (length(args) == 0L) return(invisible(list()))
  nms <- names(args)
  if (is.null(nms) || any(nms == ""))
    stop("All arguments to dynhr_set_options() must be named.")
  prev <- lapply(nms, function(nm) {
    p <- if (exists(nm, envir = .dynhr_opts, inherits = FALSE))
      get(nm, envir = .dynhr_opts) else NULL
    ## NULL un-sets (brief 32 P3c): the previous value of an unset option is
    ## returned as NULL, so the restore idiom
    ## old <- dynhr_set_options(x = v); do.call(dynhr_set_options, old)
    ## used to STORE NULL, and .dynhr_opt(x) then returned NULL instead of the
    ## registered default (a NULL grad_method, seed_base, ...).
    if (is.null(args[[nm]])) {
      if (exists(nm, envir = .dynhr_opts, inherits = FALSE))
        rm(list = nm, envir = .dynhr_opts)
    } else {
      assign(nm, args[[nm]], envir = .dynhr_opts)
    }
    p
  })
  names(prev) <- nms
  invisible(prev)
}

#' Get dynhr package-level options
#'
#' @param effective \code{FALSE} (default) returns only the options set with
#'   \code{\link{dynhr_set_options}} -- the form to save and later restore with
#'   \code{do.call(dynhr_set_options, old)}. \code{TRUE} returns EVERY
#'   registered option (see \code{\link{dynhr_set_options}} for the list) with
#'   its effective value: the value set, or else its registered default. Any
#'   set option that is not registered is appended after the registered ones.
#'   This is the snapshot stored in an estimation result's run record.
#' @return Named list. With \code{effective = TRUE} it carries an attribute
#'   \code{"set"}: the names of the options whose value came from
#'   \code{dynhr_set_options()} rather than from the registered default.
#' @examples
#' eff <- dynhr_get_options(effective = TRUE)
#' eff$power_posterior
#' attr(eff, "set")
#' @export
dynhr_get_options <- function(effective = FALSE) {
  set_opts <- as.list(.dynhr_opts)
  if (!isTRUE(effective)) return(set_opts)
  reg <- names(.dynhr_option_registry)
  out <- lapply(.dynhr_option_registry, `[[`, "default")
  for (nm in names(set_opts)) out[nm] <- list(set_opts[[nm]])
  extra <- sort(as.character(setdiff(names(set_opts), reg)), method = "radix")
  out <- out[c(reg, extra)]
  attr(out, "set") <- sort(as.character(names(set_opts)), method = "radix")
  out
}

#' Reset dynhr options to package defaults
#'
#' @param ... Character names of options to reset. If empty, resets all.
#' @export
dynhr_reset_options <- function(...) {
  nms <- c(...)
  if (length(nms) == 0L) {
    rm(list = ls(.dynhr_opts), envir = .dynhr_opts)
  } else {
    for (nm in nms)
      if (exists(nm, envir = .dynhr_opts, inherits = FALSE))
        rm(list = nm, envir = .dynhr_opts)
  }
  invisible(NULL)
}

# Internal: resolve option with arg-value > global option > default priority.
# Usage: val <- .dynhr_opt("me_variance", arg_value)
# - If arg_value is non-NULL, return it as-is.
# - Else if the option is set globally, return that.
# - Else return `default` when one is passed, otherwise the registry default
#   (`.dynhr_option_registry`; NULL for an unregistered name). A passed
#   `default` must equal the registry's -- guarded by
#   test-fix-0926-run-record.R.
.dynhr_opt <- function(name, arg_value = NULL, default) {
  if (!is.null(arg_value)) return(arg_value)
  if (exists(name, envir = .dynhr_opts, inherits = FALSE))
    return(get(name, envir = .dynhr_opts))
  if (missing(default)) return(.dynhr_option_registry[[name]]$default)
  default
}

# Internal: resolve the power-posterior (generalised-Bayes) tempering exponent
# ONCE, at log-posterior FACTORY time.
#
# Every `make_log_posterior_<x>()` factory takes `power = NULL` meaning "resolve
# `power_posterior` from the option store now". Resolving at factory time (not
# per evaluation) is what makes the exponent a property of the closure: an MCMC
# run cannot silently change target distribution half-way through because some
# other code called `dynhr_set_options(power_posterior = )`, the resolution cost
# leaves the per-draw hot path, and -- since a mirai daemon has its own
# `.dynhr_opts` (see `.dynhr_daemon_state`) -- a closure built on the host keeps
# the host's exponent wherever it is evaluated.
#
# `where` names the calling factory in the error message.
.resolve_power_posterior <- function(power, where) {
  power <- .dynhr_opt("power_posterior", power, default = 1)
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0)
    stop(where, ": `power` must be a finite scalar in (0, 1].", call. = FALSE)
  if (power > 1)
    .dynhr_warn(where, ": `power` > 1 produces a 'cold' (over-confident) ",
            "posterior. This is valid but unusual; set power <= 1 for ",
            "standard generalised-Bayes tempering.", call. = FALSE)
  power
}

# ---------------------------------------------------------------------------
# Shipping package options to parallel workers
# ---------------------------------------------------------------------------
# `.dynhr_opts` is a namespace-private ENVIRONMENT, not base `options()`, so a
# mirai daemon (a fresh R process that merely `library(dynhr)`s) starts with an
# EMPTY option store: every `dynhr_set_options(power_posterior = 0.5,
# me_variance = 1e-4, debug_kf_errors = TRUE, ...)` set in the host session was
# silently invisible to the workers, and the parallel path then evaluated a
# DIFFERENT posterior from the serial one. Base `options(dynhr.* = )` (e.g.
# `dynhr.use_rcpp`, `dynhr.me_floor_check`) has the same problem: R does not
# inherit the host's options into a spawned process.
#
# `.dynhr_daemon_state()` snapshots both stores on the host; the snapshot rides
# in the `.args` of every `everywhere()` pool-init block and is replayed there
# by `.dynhr_daemon_apply()`. Snapshot at pool-init time and apply ONCE PER
# DAEMON -- never per task -- so the cost is O(n_daemons), not O(n_draws).

# ---------------------------------------------------------------------------
# Worker / host version skew
# ---------------------------------------------------------------------------
# Every worker (mirai daemon, PSOCK node, future worker) is a fresh R process
# that `library(dynhr)`s the INSTALLED package. Under `devtools::load_all()`
# with a stale install -- or with R_LIBS pointing somewhere else -- the workers
# silently ran DIFFERENT code from the host session (this landed in the suite
# register three times). The snapshot therefore also carries the host's build
# identity, and `.dynhr_daemon_apply()` aborts with class
# `dynhr_error_worker_version_skew` when the worker's build differs:
#   * package version    -- always compared;
#   * R-code fingerprint  -- compared when the host is a pkgload dev load
#                            (a GIT_COMMIT stamp then says nothing about the
#                            code actually loaded) or has no GIT_COMMIT stamp;
#   * GIT_COMMIT stamp    -- compared otherwise (installed host).
# The check function and the fingerprint function travel INSIDE the snapshot
# as function objects: a namespace closure serialises as a reference to that
# namespace, so the BODY is the host's while names resolve in the worker's own
# dynhr -- a worker whose install predates this check is therefore still
# checked. Both bodies must stay self-contained (base R + tools only, no dynhr
# helpers).

# Hash of every function in namespace `ns` (name + deparsed closure). Deparsed
# WITHOUT srcrefs, so a load_all() namespace (srcrefs kept) and an installed,
# byte-compiled one built from the same sources hash identically.
.dynhr_code_fingerprint <- function(ns) {
  ## radix sort: locale-INDEPENDENT order (host and worker may differ in
  ## LC_COLLATE -- testthat alone switches it to C).
  nms <- sort(ls(ns, all.names = TRUE, sorted = FALSE), method = "radix")
  ctl <- c("keepNA", "keepInteger", "niceNames", "showAttributes")
  txt <- vector("list", length(nms))
  for (k in seq_along(nms)) {
    f <- get(nms[[k]], envir = ns, inherits = FALSE)
    if (is.function(f)) txt[[k]] <- c(nms[[k]], deparse(f, control = ctl))
  }
  tf <- tempfile("dynhr-fp-")
  on.exit(unlink(tf), add = TRUE)
  writeLines(unlist(txt, use.names = FALSE), tf, useBytes = TRUE)
  unname(tools::md5sum(tf))
}

# First line of the GIT_COMMIT stamp under a package root (installed layout,
# or <source>/inst under pkgload); NA when absent.
.dynhr_git_stamp_at <- function(path) {
  for (f in file.path(path, c("GIT_COMMIT", file.path("inst", "GIT_COMMIT")))) {
    if (file.exists(f)) {
      txt <- trimws(readLines(f, n = 1L, warn = FALSE))
      if (length(txt) == 1L && nzchar(txt)) return(txt)
    }
  }
  NA_character_
}

# Host fingerprint cache for an INSTALLED namespace (its code does not change
# under us); a dev load is always re-hashed.
.dynhr_build_cache <- new.env(parent = emptyenv())

# The host's build identity, shipped as `.dynhr_daemon_state()$build`.
.dynhr_host_build <- function() {
  ns     <- asNamespace("dynhr")
  path   <- getNamespaceInfo(ns, "path")
  dev    <- exists(".__DEVTOOLS__", envir = ns, inherits = FALSE)
  commit <- .dynhr_git_stamp_at(path)
  fp <- NA_character_
  if (dev || is.na(commit)) {
    if (!dev && !is.null(.dynhr_build_cache$fingerprint)) {
      fp <- .dynhr_build_cache$fingerprint
    } else {
      fp <- .dynhr_code_fingerprint(ns)
      if (!dev) .dynhr_build_cache$fingerprint <- fp
    }
  }
  list(version = unname(getNamespaceVersion(ns)), git_commit = commit,
       dev_load = dev, path = path, fingerprint = fp,
       fingerprint_fn = .dynhr_code_fingerprint,
       check = .dynhr_daemon_skew)
}

# Compare a shipped host build identity with THIS process's dynhr. Returns NULL
# when they match, otherwise the error message. Self-contained (see the header
# note): it runs with the host's body inside the worker's namespace.
.dynhr_daemon_skew <- function(build) {
  ns     <- asNamespace("dynhr")
  w_ver  <- unname(getNamespaceVersion(ns))
  w_path <- getNamespaceInfo(ns, "path")
  w_commit <- NA_character_
  for (f in file.path(w_path, c("GIT_COMMIT", file.path("inst", "GIT_COMMIT")))) {
    if (is.na(w_commit) && file.exists(f)) {
      txt <- trimws(readLines(f, n = 1L, warn = FALSE))
      if (length(txt) == 1L && nzchar(txt)) w_commit <- txt
    }
  }
  why <- NULL
  if (!identical(as.character(build$version), w_ver)) {
    why <- "the package versions differ"
  } else if (!is.na(build$fingerprint)) {
    if (!identical(build$fingerprint_fn(ns), build$fingerprint))
      why <- paste0("the versions match but the R code differs",
                    if (isTRUE(build$dev_load))
                      " (the host is a development load via pkgload / devtools::load_all())")
  } else if (!identical(build$git_commit, w_commit)) {
    why <- "the versions match but the builds come from different commits (GIT_COMMIT)"
  }
  if (is.null(why)) return(NULL)
  desc <- function(ver, commit, dev, path)
    paste0("dynhr ", ver,
           if (!is.na(commit)) paste0(" (GIT_COMMIT ", substr(commit, 1L, 12L), ")"),
           if (isTRUE(dev)) " [dev load]", " at ", path)
  paste0("dynhr worker version skew: ", why, ".\n",
         "  host  : ", desc(build$version, build$git_commit, build$dev_load,
                            build$path), "\n",
         "  worker: ", desc(w_ver, w_commit, FALSE, w_path), "\n",
         "Parallel workers load the INSTALLED package, so they would run ",
         "different code from this session. Fix: install the current source ",
         "(`R CMD INSTALL --preclean <package dir>`), or start R with R_LIBS ",
         "pointing at a library that holds a current build.")
}

# Snapshot the host's dynhr option state for shipping to a worker process.
# Returns list(opts  = <as.list(.dynhr_opts)>,
#              base  = <the dynhr.* base options>,
#              build = <.dynhr_host_build(): version / commit / fingerprint>).
.dynhr_daemon_state <- function() {
  base_all <- options()
  list(opts  = as.list(.dynhr_opts),
       base  = base_all[grep("^dynhr\\.", names(base_all))],
       build = .dynhr_host_build())
}

# Replay a `.dynhr_daemon_state()` snapshot in the current process.
# Idempotent, and AUTHORITATIVE for the package store: `.dynhr_opts` is cleared
# before the snapshot is written, so re-initialising a LIVE pool (e.g.
# `.mirai_rebind_worker_lp`) cannot leave a stale value from a previous run
# behind. Base options are merged, not cleared -- the daemon's own R defaults
# for unrelated `dynhr.*` keys are none of our business.
#
# FIRST it checks this process's build against the host's (`state$build`, see
# `.dynhr_daemon_skew`) and aborts with class `dynhr_error_worker_version_skew`
# on a mismatch. The condition is built inline rather than via `.dynhr_abort()`
# because a host-shipped copy of this body may run in a worker namespace that
# predates R/conditions.R.
.dynhr_daemon_apply <- function(state) {
  if (is.null(state) || !is.list(state)) return(invisible(NULL))
  build <- state$build
  if (is.list(build) && is.function(build$check)) {
    msg <- build$check(build)
    if (!is.null(msg))
      stop(structure(
        class = c("dynhr_error_worker_version_skew", "dynhr_error", "error",
                  "condition"),
        list(message = msg, call = NULL)))
  }
  if (!is.null(state$opts)) {
    old <- ls(.dynhr_opts, all.names = TRUE)
    if (length(old)) rm(list = old, envir = .dynhr_opts)
    if (length(state$opts)) list2env(state$opts, envir = .dynhr_opts)
  }
  if (length(state$base)) do.call(options, state$base)
  invisible(NULL)
}


#' Situation report: which dynhr is loaded, and what parallel workers will load
#'
#' Prints the facts needed to tell whether a parallel run executes the code you
#' think it does. Parallel workers (mirai daemons, PSOCK nodes, future workers)
#' are fresh R processes that \code{library(dynhr)} the INSTALLED package from
#' \code{.libPaths()}; under \code{devtools::load_all()} that can be an older build than
#' the one in your session. dynhr refuses to start such a pool (error class
#' \code{dynhr_error_worker_version_skew}); \code{dynhr_sitrep()} shows the same
#' comparison before you start one.
#'
#' Reported: the loaded dynhr version, its \code{GIT_COMMIT} stamp, whether it was
#' loaded via pkgload (\code{devtools::load_all()}), the version and stamp of the
#' installed copy that workers would load (flagged when it differs), the BLAS
#' and LAPACK R is linked against, the core count, and the mirai daemon status.
#'
#' @return Invisibly, a list of class \code{"dynhr_sitrep"} holding the reported
#'   fields; printing it reproduces the report.
#' @examples
#' dynhr_sitrep()
#' @export
dynhr_sitrep <- function() {
  ns   <- asNamespace("dynhr")
  path <- getNamespaceInfo(ns, "path")
  inst_path <- NA_character_
  for (lib in .libPaths()) {
    if (is.na(inst_path) && file.exists(file.path(lib, "dynhr", "DESCRIPTION")))
      inst_path <- file.path(lib, "dynhr")
  }
  inst_version <- NA_character_
  inst_commit  <- NA_character_
  if (!is.na(inst_path)) {
    inst_version <- unname(read.dcf(file.path(inst_path, "DESCRIPTION"),
                                    fields = "Version")[1L, 1L])
    inst_commit  <- .dynhr_git_stamp_at(inst_path)
  }
  si <- utils::sessionInfo()
  st <- mirai::status()
  out <- structure(list(
    version           = unname(getNamespaceVersion(ns)),
    git_commit        = .dynhr_git_stamp_at(path),
    dev_load          = exists(".__DEVTOOLS__", envir = ns, inherits = FALSE),
    path              = path,
    installed_path    = inst_path,
    installed_version = inst_version,
    installed_commit  = inst_commit,
    blas              = if (is.null(si$BLAS)) NA_character_ else si$BLAS,
    lapack            = if (is.null(si$LAPACK)) NA_character_ else si$LAPACK,
    cores             = parallel::detectCores(),
    mirai_connections = as.integer(st$connections),
    mirai_daemons     = as.character(st$daemons)
  ), class = "dynhr_sitrep")
  print(out)
  invisible(out)
}

#' @export
print.dynhr_sitrep <- function(x, ...) {
  na_or <- function(v, alt = "none")
    if (length(v) != 1L || is.na(v)) alt else as.character(v)
  skew <- !is.na(x$installed_version) &&
    (!identical(x$installed_version, x$version) ||
       (!x$dev_load && !identical(x$installed_commit, x$git_commit)))
  worker_line <- if (is.na(x$installed_version)) {
    "NOT INSTALLED -- parallel workers cannot load dynhr"
  } else {
    paste0("dynhr ", x$installed_version,
           if (!is.na(x$installed_commit))
             paste0(" (GIT_COMMIT ", substr(x$installed_commit, 1L, 12L), ")"),
           " at ", x$installed_path,
           if (skew) "\n                <-- DIFFERS from the loaded build: parallel runs will abort"
           else if (x$dev_load) "\n                (same version; R-code identity is checked at pool start)")
  }
  cat("dynhr situation report\n",
      "  loaded      : dynhr ", x$version,
      if (x$dev_load) " [dev load via pkgload]", " at ", x$path, "\n",
      "  GIT_COMMIT  : ", na_or(x$git_commit), "\n",
      "  workers load: ", worker_line, "\n",
      "  BLAS        : ", na_or(x$blas, "unknown"), "\n",
      "  LAPACK      : ", na_or(x$lapack, "unknown"), "\n",
      "  cores       : ", na_or(x$cores, "unknown"), "\n",
      "  mirai       : ",
      if (isTRUE(x$mirai_connections > 0L))
        paste0(x$mirai_connections, " daemon connection(s) at ",
               paste(x$mirai_daemons, collapse = ", "))
      else "no daemons running", "\n",
      sep = "")
  invisible(x)
}
