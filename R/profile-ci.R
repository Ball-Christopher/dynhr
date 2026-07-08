## R/profile-ci.R
## --------------------------------------------------------------------------
## profile_ci() -- profile-likelihood/posterior confidence intervals with
## FEASIBILITY-SEEDED nuisance re-optimization (pathological-DSGE paper,
## DYNHR_GAPS.md gap #5, referee(2) A2).
##
## THE MODE-SEEDING TRAP: a naive profile re-optimizes the nuisance
## parameters at each grid value of the profiled parameter, starting FROM
## THE MODE's nuisance values. That start can be Blanchard-Kahn INFEASIBLE
## (indeterminate) at the fixed grid value even though it was feasible at
## the mode itself -- the model's determinacy region is generally curved in
## parameter space (e.g. the "generalized Taylor principle" is not simply
## psi1 >= 1; see replication/R/10d_weakid_profile2.R in the
## pathological-dsge-paper repo). A local optimizer seeded there either
## fails outright or converges to the nearest feasible point, silently
## reporting the WRONG bound of the profile confidence set -- often
## overstating identification strength.
##
## The fix: at each grid point, re-optimize from (a) the mode start AND
## (b) `n_seeds` additional starts drawn from the prior and SCREENED for
## BK-feasibility (and stationarity) at that fixed grid value BEFORE
## optimizing; keep the best (highest logpost) result across all starts
## that converge. When no start is feasible, the grid point is recorded as
## `infeasible = TRUE` (not silently -Inf-best): omit it from the profile
## set exactly like a real infeasible region, without contaminating the
## comparison of which starts DID work.
## --------------------------------------------------------------------------


## Internal: Blanchard-Kahn feasibility screen at a given full parameter
## vector `params` (named numeric, model$param_values-shaped). Mirrors
## 10d_weakid_profile2.R's `bk_at()`: solve steady state, then perturbation,
## and check `dr$bk_satisfied`. Returns FALSE (not an error) on any failure
## so it is safe to call inside a screening loop.
.profile_ci_bk_feasible <- function(model, compiled, params) {
  ss <- tryCatch(solve_steady_state(model, compiled, params, verbose = FALSE),
                error = function(e) NULL)
  if (is.null(ss) || !isTRUE(ss$converged)) return(FALSE)
  ## An indeterminate/explosive draw is an EXPECTED outcome of the
  ## feasibility screen (that is the whole point of screening), not a
  ## warning-worthy event for the caller -- suppress the
  ## .solve_from_system() "BK violation" warning here (mirrors 10d's
  ## bk_at(), which treats it as a plain FALSE return).
  dr <- tryCatch(
    suppressWarnings(solve_perturbation(model, compiled, ss$ss,
                                        params = ss$params %||% params,
                                        order = 1L, verbose = FALSE)),
    error = function(e) NULL
  )
  if (is.null(dr)) return(FALSE)
  isTRUE(dr$bk_satisfied)
}


## Internal: build the model-based per-grid-point objective and feasibility
## screen from (model, data, prior_spec, obs_vars, compiled). Returns a list
## with $logpost_fn(theta_named_full) -> scalar (-Inf on any failure, never
## errors) and $feasible_fn(theta_named_full) -> logical.
.profile_ci_model_objective <- function(model, data, prior_spec, obs_vars,
                                        compiled, me_variance, likelihood,
                                        ...) {
  lp_fn <- make_log_posterior(model, data, prior_spec, obs_vars = obs_vars,
                              compiled = compiled, me_variance = me_variance,
                              likelihood = likelihood, ...)
  logpost_fn <- function(theta) {
    ## Same rationale as .profile_ci_bk_feasible(): an infeasible draw during
    ## re-optimization is expected (the optimizer is actively probing
    ## infeasible neighbourhoods), not a warning-worthy event.
    r <- tryCatch(suppressWarnings(lp_fn(theta)), error = function(e) NULL)
    if (is.null(r) || !is.list(r) || !is.finite(r$logpost)) return(-Inf)
    r$logpost
  }
  feasible_fn <- function(theta) {
    params <- .apply_theta_to_params(model, theta)
    .profile_ci_bk_feasible(model, compiled, params)
  }
  list(logpost_fn = logpost_fn, feasible_fn = feasible_fn)
}


## Internal: draw `n` feasibility-screened nuisance starting vectors at a
## fixed profiled-parameter value. `free_names` are the nuisance parameters
## being re-optimized (all estimated params except `param`). Draws come from
## `prior_sampler()` (a `.smc_make_prior_sampler(prior_spec)` closure); each
## candidate is screened with `feasible_fn` before being kept. Gives up after
## `max_tries` attempts (mirrors 10d's `tries < 600` budget) -- returns
## however many feasible starts were found (possibly zero; the caller must
## handle that as "no feasible seed").
.profile_ci_seed_starts <- function(param, v, free_names, prior_sampler,
                                    theta_mode, feasible_fn, n_seeds,
                                    max_tries = 600L) {
  starts <- list()
  tries  <- 0L
  while (length(starts) < n_seeds && tries < max_tries) {
    tries <- tries + 1L
    draw  <- prior_sampler()
    th    <- theta_mode
    th[names(draw)] <- draw
    th[param] <- v
    if (isTRUE(feasible_fn(th))) {
      starts[[length(starts) + 1L]] <- th[free_names]
    }
  }
  starts
}


## Internal: re-optimize the nuisance parameters at one grid value `v`,
## trying the mode-start plus feasibility-screened seed starts, and return
## the best result. `objective_fn(theta_full)` is a SCALAR objective to be
## MAXIMIZED (log-posterior or log-likelihood); csminwel minimizes, so the
## wrapper negates.
##
## Returns list(value = best logpost (or -Inf), theta = best full theta (or
## NULL), convergence, feasible = TRUE/FALSE (was ANY start feasible),
## won_by = "mode" | "seed" | NA, n_seeds_tried, n_seeds_feasible).
.profile_ci_optimize_at <- function(param, v, free_names, theta_mode,
                                    objective_fn, feasible_fn, prior_sampler,
                                    n_seeds, n_iter, optimizer_crit = 1e-6) {
  th_mode_full        <- theta_mode
  th_mode_full[param]  <- v
  mode_start_feasible <- isTRUE(feasible_fn(th_mode_full))

  seed_starts <- if (n_seeds > 0L && !is.null(prior_sampler)) {
    .profile_ci_seed_starts(param, v, free_names, prior_sampler, theta_mode,
                            feasible_fn, n_seeds)
  } else list()

  starts <- list()
  start_labels <- character(0)
  if (mode_start_feasible) {
    starts[[length(starts) + 1L]] <- th_mode_full[free_names]
    start_labels <- c(start_labels, "mode")
  }
  if (length(seed_starts) > 0L) {
    starts <- c(starts, seed_starts)
    start_labels <- c(start_labels, rep("seed", length(seed_starts)))
  }

  if (length(starts) == 0L) {
    return(list(value = -Inf, theta = NULL, convergence = NA_integer_,
               feasible = FALSE, won_by = NA_character_,
               n_seeds_tried = length(seed_starts),
               n_seeds_feasible = length(seed_starts)))
  }

  neg_obj <- function(xf) {
    th <- th_mode_full
    th[free_names] <- xf
    val <- objective_fn(th)
    if (!is.finite(val)) 1e10 else -val
  }

  best_val   <- -Inf
  best_theta <- NULL
  best_label <- NA_character_
  for (i in seq_along(starts)) {
    x0 <- as.numeric(starts[[i]])
    r  <- tryCatch(
      csminwel(neg_obj, x0 = x0, H0 = diag(length(free_names)) * 1e-2,
              crit = optimizer_crit, nit = n_iter, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(r) || !is.finite(r$value)) next
    val <- -r$value
    if (val > best_val) {
      best_val   <- val
      th <- th_mode_full
      th[free_names] <- r$x
      best_theta <- th
      best_label <- start_labels[i]
    }
  }

  list(value = best_val, theta = best_theta,
      convergence = if (is.null(best_theta)) NA_integer_ else 0L,
      feasible = TRUE, won_by = best_label,
      n_seeds_tried = length(seed_starts),
      n_seeds_feasible = length(seed_starts))
}


## Internal: turn a logical "in-set" vector over a sorted grid into a list
## of [lo, hi] intervals (runs of consecutive TRUE grid points), handling a
## DISCONNECTED profile set. Singleton points become a zero-width interval
## [v, v].
.profile_ci_runs <- function(grid, in_set) {
  if (!any(in_set)) return(list())
  idx <- which(in_set)
  ## Split idx into runs of consecutive grid POSITIONS (not values), so a
  ## non-uniform grid is still handled correctly.
  brk <- c(0, which(diff(idx) != 1L), length(idx))
  runs <- list()
  for (i in seq_len(length(brk) - 1L)) {
    seg <- idx[(brk[i] + 1L):brk[i + 1L]]
    runs[[length(runs) + 1L]] <- c(lo = grid[min(seg)], hi = grid[max(seg)])
  }
  runs
}


#' Profile-likelihood/posterior confidence interval with feasibility-seeded
#' nuisance re-optimization
#'
#' Profiles the log-posterior (or a user-supplied objective) over a grid of
#' values of one parameter, re-optimizing the remaining ("nuisance")
#' parameters at each grid point. Returns the profile curve, per-point
#' feasibility/convergence diagnostics, and the level-set confidence interval
#' (or disconnected set of intervals) from the \eqn{\chi^2(1)}
#' likelihood-ratio cutoff.
#'
#' \strong{The mode-seeding trap this closes:} a naive profile re-optimizes
#' nuisance parameters starting FROM THE POSTERIOR MODE's nuisance values at
#' every grid point. That start can be Blanchard-Kahn \emph{infeasible}
#' (indeterminate) at a grid value far from the mode even when it was
#' feasible at the mode itself, because a DSGE model's determinacy region is
#' generally a curved manifold in parameter space (e.g. the "generalized
#' Taylor principle" boundary is not the naive \eqn{\psi_1 = 1} line). A
#' local optimizer seeded there either fails or converges to the nearest
#' feasible sliver, silently reporting the WRONG profile value and
#' overstating identification strength. \code{profile_ci()} additionally
#' seeds the nuisance re-optimization from \code{n_seeds} prior draws that
#' are SCREENED for feasibility at the fixed grid value \emph{before}
#' optimizing, and reports (\code{$won_by}) whether the winning start at each
#' grid point was the mode start or a feasibility-screened seed -- the
#' paper's headline diagnostic for when the naive profile would have been
#' wrong.
#'
#' Two calling conventions are supported:
#' \itemize{
#'   \item \strong{Model-based}: supply \code{model}, \code{data},
#'     \code{prior_spec}, \code{obs_vars}, \code{compiled} (as for
#'     \code{\link{make_posterior}}); feasibility is the model's
#'     Blanchard-Kahn / stationarity check.
#'   \item \strong{Objective-based} (for testing or non-DSGE use):
#'     supply \code{objective_fn(theta)} (a named-vector -> scalar
#'     log-posterior/log-likelihood to maximize) and, optionally,
#'     \code{feasible_fn(theta)} (a named-vector -> logical feasibility
#'     screen; default always-feasible) and \code{prior_sampler()} (a
#'     no-argument closure returning a named numeric draw used for
#'     feasibility-seeded starts).
#' }
#'
#' @param param      Character scalar: the name of the parameter to profile.
#' @param model,data,prior_spec,obs_vars,compiled Model-based inputs; see
#'   \code{\link{make_posterior}}. Required unless \code{objective_fn} is
#'   supplied.
#' @param mode       Named numeric vector: the (unconstrained) posterior
#'   mode, e.g. \code{run_mode_finding(...)$theta_mode}. Required (supplies
#'   both the nuisance start and the free-parameter names via its own
#'   names, together with \code{param}).
#' @param objective_fn Optional function(theta_named_vector) -> scalar to
#'   maximize, bypassing the model/data/prior_spec/compiled path entirely.
#'   When supplied, \code{model} etc. are not required.
#' @param feasible_fn Optional function(theta_named_vector) -> logical,
#'   used only with \code{objective_fn}. Default: always \code{TRUE} (no
#'   feasibility screening -- appropriate for objectives with no
#'   feasibility constraint).
#' @param prior_sampler Optional no-argument closure returning a named
#'   numeric draw, used to generate feasibility-screened seed starts. For
#'   the model-based path this defaults to
#'   \code{.smc_make_prior_sampler(prior_spec)}; for the \code{objective_fn}
#'   path it must be supplied explicitly if \code{n_seeds > 0}.
#' @param free       Character vector of nuisance parameter names to
#'   re-optimize at each grid point. Defaults to
#'   \code{setdiff(names(mode), param)}.
#' @param grid       Numeric vector of explicit grid values for \code{param},
#'   OR \code{NULL} to auto-generate \code{n_grid} points spanning
#'   \code{mode[param] +/- grid_width_sd * prior_spec$std} (model-based
#'   path) -- supply \code{grid} explicitly when using \code{objective_fn}
#'   without a \code{prior_spec}.
#' @param n_grid     Number of auto-generated grid points (default 15L).
#'   Ignored when \code{grid} is supplied.
#' @param grid_width_sd Half-width of the auto-generated grid in prior
#'   standard deviations (default 4).
#' @param level      Confidence level for the \eqn{\chi^2(1)} likelihood-ratio
#'   cutoff (default 0.95).
#' @param n_seeds    Number of feasibility-screened nuisance seed starts to
#'   try at each grid point, IN ADDITION to the mode-start (default 8L).
#'   Set to 0 to reproduce the naive (mode-seeded-only) profile.
#' @param n_iter     csminwel iteration budget per re-optimization
#'   (default 200L).
#' @param seed       Optional integer RNG seed (sets \code{set.seed()} once
#'   at the top for reproducible seed draws).
#' @param cores      Number of cores for \code{parallel::mclapply} over grid
#'   points (default 1L; not supported on Windows, where it silently
#'   falls back to serial).
#' @param me_variance,likelihood,... Forwarded to
#'   \code{make_log_posterior} on the model-based path (ignored with
#'   \code{objective_fn}).
#' @param verbose    Print a progress dot per grid point (default TRUE).
#'
#' @return An object of class \code{"dynhr_profile_ci"}, a list with:
#'   \describe{
#'     \item{\code{grid}}{The grid of \code{param} values, sorted.}
#'     \item{\code{profile}}{Profile log-posterior at each grid point
#'       (\code{-Inf} where infeasible).}
#'     \item{\code{drop}}{\code{profile - max(profile, na.rm=TRUE)}, the
#'       usual "nats below max" profile display.}
#'     \item{\code{feasible}}{Logical: was ANY start (mode or seed) feasible
#'       at this grid point.}
#'     \item{\code{won_by}}{Character: \code{"mode"} or \code{"seed"} --
#'       which start produced the best value at this grid point (\code{NA}
#'       if infeasible). \strong{The paper's headline diagnostic}: any
#'       \code{"seed"} entries are grid points where the naive mode-seeded
#'       profile would have gotten the wrong answer.}
#'     \item{\code{theta_at}}{List of full parameter vectors at the
#'       optimum for each grid point (\code{NULL} where infeasible).}
#'     \item{\code{intervals}}{Data frame with columns \code{lo}, \code{hi}:
#'       the (possibly multiple, i.e. disconnected) grid intervals inside
#'       the level-set confidence region.}
#'     \item{\code{connected}}{Logical: \code{TRUE} iff \code{intervals} has
#'       exactly one row.}
#'     \item{\code{ci}}{Numeric \code{c(lo, hi)}: the outer envelope of
#'       \code{intervals} (\code{min(lo), max(hi)}), for convenient use when
#'       the caller does not care about disconnectedness.}
#'     \item{\code{level}, \code{crit}}{The requested level and the
#'       \eqn{\chi^2(1)} cutoff used (\code{qchisq(level, df = 1) / 2} nats
#'       below the max).}
#'     \item{\code{mode_value}}{The profile value at \code{mode[param]}
#'       (approximately the unprofiled mode logpost, used as the
#'       normalisation constant).}
#'   }
#'
#' @seealso \code{\link{run_mode_finding}}, \code{\link{make_posterior}},
#'   \code{\link{robust_confidence_set}}
#' @references
#'   Referee(2) A2 correction to the pathological-DSGE weak-identification
#'   case study: replication \code{10d_weakid_profile2.R} (naive mode-seeded
#'   profile vs. feasibility-seeded true profile on \code{nk_small}).
#' @export
profile_ci <- function(param,
                       model = NULL, data = NULL, prior_spec = NULL,
                       obs_vars = NULL, compiled = NULL,
                       mode,
                       objective_fn = NULL, feasible_fn = NULL,
                       prior_sampler = NULL,
                       free = NULL,
                       grid = NULL, n_grid = 15L, grid_width_sd = 4,
                       level = 0.95, n_seeds = 8L, n_iter = 200L,
                       seed = NULL, cores = 1L,
                       me_variance = 0, likelihood = "gaussian",
                       verbose = TRUE, ...) {

  if (missing(mode) || is.null(mode) || is.null(names(mode)) ||
      any(!nzchar(names(mode))))
    stop("profile_ci: `mode` must be a named numeric vector (e.g. ",
         "run_mode_finding(...)$theta_mode).", call. = FALSE)
  if (!is.character(param) || length(param) != 1L || !nzchar(param))
    stop("profile_ci: `param` must be a single parameter name.", call. = FALSE)
  if (!(param %in% names(mode)))
    stop("profile_ci: `param` = \"", param, "\" is not a name in `mode`.",
         call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  use_objective <- !is.null(objective_fn)
  if (!use_objective) {
    if (is.null(model) || is.null(data) || is.null(prior_spec) ||
        is.null(compiled))
      stop("profile_ci: supply either `objective_fn` or all of `model`, ",
           "`data`, `prior_spec`, `compiled` (obs_vars defaults from the ",
           "model's varobs).", call. = FALSE)
    built <- .profile_ci_model_objective(model, data, prior_spec, obs_vars,
                                         compiled, me_variance, likelihood, ...)
    objective_fn <- built$logpost_fn
    if (is.null(feasible_fn)) feasible_fn <- built$feasible_fn
    if (is.null(prior_sampler) && n_seeds > 0L)
      prior_sampler <- .smc_make_prior_sampler(prior_spec)
  } else {
    if (is.null(feasible_fn)) feasible_fn <- function(theta) TRUE
  }

  free_names <- free %||% setdiff(names(mode), param)
  if (length(free_names) == 0L)
    stop("profile_ci: no nuisance parameters to re-optimize (`free` is ",
         "empty after excluding `param`). Nothing to profile if `mode` has ",
         "only one parameter -- the 'profile' is then just the raw ",
         "objective curve; call objective_fn directly instead.",
         call. = FALSE)

  ## ---- Grid construction -------------------------------------------------
  if (is.null(grid)) {
    if (use_objective)
      stop("profile_ci: `grid` must be supplied explicitly when using ",
           "`objective_fn` without a model-based `prior_spec` (no prior sd ",
           "to auto-scale the grid width).", call. = FALSE)
    sd_p <- prior_spec$std[match(param, prior_spec$name)]
    if (is.na(sd_p) || !is.finite(sd_p) || sd_p <= 0)
      stop("profile_ci: could not read a finite prior std for `param` = \"",
           param, "\" to auto-build the grid; supply `grid` explicitly.",
           call. = FALSE)
    lo <- mode[[param]] - grid_width_sd * sd_p
    hi <- mode[[param]] + grid_width_sd * sd_p
    lo_bound <- prior_spec$lower[match(param, prior_spec$name)]
    hi_bound <- prior_spec$upper[match(param, prior_spec$name)]
    if (!is.na(lo_bound)) lo <- max(lo, lo_bound)
    if (!is.na(hi_bound)) hi <- min(hi, hi_bound)
    grid <- seq(lo, hi, length.out = n_grid)
  }
  grid <- sort(unique(as.numeric(grid)))
  n_g  <- length(grid)

  ## ---- Per-grid-point re-optimization ------------------------------------
  one_point <- function(v) {
    if (isTRUE(verbose)) cat(".")
    .profile_ci_optimize_at(param, v, free_names, mode, objective_fn,
                            feasible_fn, prior_sampler, n_seeds, n_iter)
  }

  use_parallel <- is.numeric(cores) && cores > 1L &&
    .Platform$OS.type != "windows"
  results <- if (use_parallel) {
    parallel::mclapply(grid, one_point, mc.cores = cores)
  } else {
    lapply(grid, one_point)
  }
  if (isTRUE(verbose)) cat("\n")

  profile_val <- vapply(results, `[[`, numeric(1), "value")
  feasible_v  <- vapply(results, `[[`, logical(1), "feasible")
  won_by_v    <- vapply(results, function(r) r$won_by %||% NA_character_,
                        character(1))
  theta_at    <- lapply(results, `[[`, "theta")
  n_seeds_tried    <- vapply(results, `[[`, integer(1), "n_seeds_tried")
  n_seeds_feasible <- vapply(results, `[[`, integer(1), "n_seeds_feasible")

  if (!any(feasible_v))
    stop("profile_ci: NO grid point had a feasible start (mode or seed). ",
         "Widen `grid`, increase `n_seeds`, or check that `mode` is itself ",
         "feasible.", call. = FALSE)

  finite_val <- profile_val[feasible_v]
  max_val    <- max(finite_val, na.rm = TRUE)
  drop_v     <- ifelse(feasible_v, profile_val - max_val, NA_real_)

  ## chi^2(1) likelihood-ratio cutoff, expressed as nats below the max
  ## (LR = 2*(max - v) <= qchisq(level, 1)  <=>  v - max >= -crit_nats).
  crit_nats <- stats::qchisq(level, df = 1L) / 2

  in_set <- feasible_v & is.finite(drop_v) & (drop_v >= -crit_nats)
  intervals_raw <- .profile_ci_runs(grid, in_set)
  intervals <- if (length(intervals_raw) == 0L) {
    data.frame(lo = numeric(0), hi = numeric(0))
  } else {
    do.call(rbind.data.frame, lapply(intervals_raw, function(r)
      data.frame(lo = unname(r["lo"]), hi = unname(r["hi"]))))
  }

  ci <- if (nrow(intervals) == 0L) c(NA_real_, NA_real_)
       else c(min(intervals$lo), max(intervals$hi))

  mode_idx   <- which.min(abs(grid - mode[[param]]))
  mode_value <- profile_val[mode_idx]

  structure(list(
    param       = param,
    grid        = grid,
    profile     = profile_val,
    drop        = drop_v,
    feasible    = feasible_v,
    won_by      = won_by_v,
    theta_at    = theta_at,
    n_seeds_tried    = n_seeds_tried,
    n_seeds_feasible = n_seeds_feasible,
    intervals   = intervals,
    connected   = nrow(intervals) <= 1L,
    ci          = ci,
    level       = level,
    crit        = crit_nats,
    mode_value  = mode_value,
    mode_param  = mode[[param]]
  ), class = "dynhr_profile_ci")
}


#' @export
print.dynhr_profile_ci <- function(x, ...) {
  cat(sprintf("Profile CI for '%s' at level %.3g\n", x$param, x$level))
  cat(sprintf("  mode value at %s = %.6g: %.4f\n", x$param, x$mode_param,
             x$mode_value))
  if (nrow(x$intervals) == 0L) {
    cat("  level set: EMPTY (no grid point within the LR cutoff)\n")
  } else {
    cat(sprintf("  connected: %s\n", x$connected))
    for (i in seq_len(nrow(x$intervals))) {
      cat(sprintf("  interval %d: [%.6g, %.6g]\n", i, x$intervals$lo[i],
                 x$intervals$hi[i]))
    }
  }
  n_seed_wins <- sum(x$won_by == "seed", na.rm = TRUE)
  if (n_seed_wins > 0L) {
    cat(sprintf(paste0("  NOTE: %d/%d grid point(s) were won by a ",
                       "feasibility-seeded start (not the mode start) -- ",
                       "the naive mode-seeded profile would have been ",
                       "wrong there.\n"), n_seed_wins, length(x$won_by)))
  }
  invisible(x)
}
