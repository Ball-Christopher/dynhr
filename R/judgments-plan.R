## R/judgments-plan.R
## --------------------------------------------------------------------------
## Unified judgments / plan layer — ROADMAP Tier 7 item 7.
##
## A dynhr_plan bundles three kinds of model interventions in one object:
##
##   plan_tune(var, periods, values, stderr)    -- in-sample filter tunes
##   plan_condition(var, horizons, values, ...) -- out-of-sample conditions
##   plan_shock_scale(var, periods, scales)     -- in-sample shock scaling
##
## REFERENCE-FRAME DISCIPLINE (critical):
##   - In-sample entries (plan_tune, plan_shock_scale) use SAMPLE PERIOD
##     indexing: integers 1..T relative to sample_start, or Dynare-style
##     date literals (e.g. "1990q1") that are resolved at adapter time.
##   - Out-of-sample entries (plan_condition) use FORECAST HORIZON offsets:
##     1-based integers from the forecast origin (terminal filter state).
##   These two reference frames are kept separate and are NEVER mixed.
##
## Design: WRAP, not replace.  dynhr_plan() is a user-facing front end;
## adapters convert it to the three existing internal representations:
##   plan_to_filter_tunes(plan, sample_start)  -> filter_tunes_spec
##   plan_to_conditions(plan)                  -> data.frame + attributes
##   plan_to_shock_scale_spec(plan)            -> het_shocks_spec
##
## ROADMAP Tier 7 item 7 · implementation 2026-06-12
## --------------------------------------------------------------------------


## ===========================================================================
## Entry constructors
## ===========================================================================

#' Specify a single in-sample filter-tune entry for a dynhr_plan
#'
#' Creates a tune entry for use inside \code{\link{dynhr_plan}()}.  The entry
#' targets one endogenous variable over a set of sample periods.
#'
#' \strong{Reference frame}: \code{periods} are SAMPLE period indices (1-based
#' integers relative to \code{sample_start}), OR Dynare-style date literals
#' (e.g. \code{"1990q1"}) which are resolved when the plan is used in
#' estimation (a \code{sample_start} must then be supplied to
#' \code{plan_to_filter_tunes()}).  This reference frame is distinct from the
#' forecast-horizon offsets used by \code{\link{plan_condition}()}.
#'
#' @param var      Name of the endogenous variable to tune (character scalar).
#'   Must be an unobserved variable (not already in \code{obs_names}).
#' @param periods  Integer vector of sample periods (1-based), or character
#'   vector of Dynare-style date literals.
#' @param values   Numeric vector of tuned values, one per period.  A scalar
#'   is recycled across all periods.
#' @param stderr   Numeric scalar giving measurement-noise standard deviation
#'   at the tuned periods (soft tune).  \code{NULL} (default) = hard tune:
#'   variable is pinned exactly.
#'
#' @return A \code{"plan_tune_entry"} object for use in
#'   \code{\link{dynhr_plan}()}.
#'
#' @seealso \code{\link{dynhr_plan}}, \code{\link{plan_condition}},
#'   \code{\link{plan_shock_scale}}
#' @export
plan_tune <- function(var, periods, values, stderr = NULL) {
  if (!is.character(var) || length(var) != 1L || nchar(var) == 0L)
    stop("plan_tune: 'var' must be a non-empty character scalar.", call. = FALSE)
  if (length(periods) == 0L)
    stop(sprintf("plan_tune: var '%s' has zero periods.", var), call. = FALSE)

  ## Accept date literals (character) or integers
  date_literals <- is.character(periods) && !all(grepl("^-?[0-9]+$", periods))
  if (!date_literals) {
    periods <- as.integer(periods)
    if (any(!is.finite(periods)))
      stop(sprintf("plan_tune: var '%s' has non-finite periods.", var),
           call. = FALSE)
  }

  values <- as.numeric(values)
  if (length(values) == 1L && length(periods) > 1L)
    values <- rep(values, length(periods))
  if (length(values) != length(periods))
    stop(sprintf(
      "plan_tune: var '%s' has %d periods but %d values.",
      var, length(periods), length(values)), call. = FALSE)

  if (!is.null(stderr)) {
    stderr <- as.numeric(stderr)
    if (length(stderr) == 1L && length(periods) > 1L)
      stderr <- rep(stderr, length(periods))
    if (length(stderr) != length(periods))
      stop(sprintf(
        "plan_tune: var '%s' has %d periods but %d stderr values.",
        var, length(periods), length(stderr)), call. = FALSE)
    if (any(!is.finite(stderr) | stderr < 0))
      stop(sprintf("plan_tune: var '%s' stderr must be finite and non-negative.",
                   var), call. = FALSE)
  }

  entry <- list(var = var, periods = periods, values = values,
                stderr = stderr, date_literals = date_literals)
  class(entry) <- "plan_tune_entry"
  entry
}


#' Specify a single out-of-sample condition entry for a dynhr_plan
#'
#' Creates a condition entry for use inside \code{\link{dynhr_plan}()}.  The
#' entry pins (or softly nudges) an observable variable at one or more forecast
#' horizons.
#'
#' \strong{Reference frame}: \code{horizons} are FORECAST HORIZON offsets
#' (1-based integers from the forecast origin, i.e. the terminal filter state
#' \eqn{s_{T|T}}).  This reference frame is distinct from the sample period
#' indices used by \code{\link{plan_tune}()} and
#' \code{\link{plan_shock_scale}()}.
#'
#' \strong{Anticipated vs unanticipated}: all condition entries in one plan
#' must share the same \code{type}.  Mixing is not supported because the
#' Waggoner-Zha anticipated path (stacked QR) and the unanticipated path
#' (period-by-period) produce materially different shock paths even for
#' identical conditions.  Supply separate plans if you need both.
#'
#' \strong{Soft conditions}: when \code{method = "soft"}, the \code{type}
#' argument is ignored by the underlying \code{.soft_forecast()} engine
#' (forward-KF unanticipated path).  A warning is issued if
#' \code{type = "anticipated"} is combined with \code{method = "soft"}.
#'
#' @param var      Name of the observable variable to condition on (character
#'   scalar).  Must be in \code{obs_names} at forecast time.
#' @param horizons Integer vector of forecast horizons (1-based offsets from
#'   the forecast origin).
#' @param values   Numeric vector of conditioned values, one per horizon.  A
#'   scalar is recycled.
#' @param type     Character; one of \code{"unanticipated"} (default) or
#'   \code{"anticipated"}.  Ignored when \code{method = "soft"}.
#' @param method   Character; one of \code{"hard"} (default) or
#'   \code{"soft"}.
#' @param stderr   Numeric scalar giving the pseudo-observation noise standard
#'   deviation for soft conditions.  Ignored when \code{method = "hard"}.
#'   \code{NULL} = use the default in \code{.soft_forecast()}.
#'
#' @return A \code{"plan_condition_entry"} object for use in
#'   \code{\link{dynhr_plan}()}.
#'
#' @seealso \code{\link{dynhr_plan}}, \code{\link{plan_tune}},
#'   \code{\link{plan_shock_scale}}
#' @export
plan_condition <- function(var, horizons, values,
                           type   = c("unanticipated", "anticipated"),
                           method = c("hard", "soft"),
                           stderr = NULL) {
  if (!is.character(var) || length(var) != 1L || nchar(var) == 0L)
    stop("plan_condition: 'var' must be a non-empty character scalar.", call. = FALSE)

  type   <- match.arg(type)
  method <- match.arg(method)

  ## Warn when type = "anticipated" is combined with method = "soft"
  ## (type is ignored by .soft_forecast; L5 landmine)
  if (method == "soft" && type == "anticipated")
    warning(
      "plan_condition: type = \"anticipated\" is ignored when method = \"soft\". ",
      "The soft-conditioning engine (.soft_forecast) always uses a forward-KF ",
      "unanticipated path regardless of type.",
      call. = FALSE)

  horizons <- as.integer(horizons)
  if (length(horizons) == 0L)
    stop(sprintf("plan_condition: var '%s' has zero horizons.", var), call. = FALSE)
  if (any(!is.finite(horizons)) || any(horizons < 1L))
    stop(sprintf("plan_condition: var '%s' horizons must be finite positive integers.",
                 var), call. = FALSE)

  values <- as.numeric(values)
  if (length(values) == 1L && length(horizons) > 1L)
    values <- rep(values, length(horizons))
  if (length(values) != length(horizons))
    stop(sprintf(
      "plan_condition: var '%s' has %d horizons but %d values.",
      var, length(horizons), length(values)), call. = FALSE)

  if (!is.null(stderr)) {
    stderr <- as.numeric(stderr)
    if (length(stderr) == 1L && length(horizons) > 1L)
      stderr <- rep(stderr, length(horizons))
    if (length(stderr) != length(horizons))
      stop(sprintf(
        "plan_condition: var '%s' has %d horizons but %d stderr values.",
        var, length(horizons), length(stderr)), call. = FALSE)
  }

  entry <- list(var = var, horizons = horizons, values = values,
                type = type, method = method, stderr = stderr)
  class(entry) <- "plan_condition_entry"
  entry
}


#' Specify a single shock-scaling entry for a dynhr_plan
#'
#' Creates a shock-scale entry for use inside \code{\link{dynhr_plan}()}.  The
#' entry applies multiplicative scale factors to a shock's standard deviation
#' over a set of sample periods.
#'
#' \strong{Reference frame}: \code{periods} are SAMPLE period indices (1-based
#' integers relative to \code{sample_start}).  This reference frame is
#' distinct from the forecast-horizon offsets used by
#' \code{\link{plan_condition}()}.
#'
#' \strong{Out-of-sample shock scaling}: \code{conditional_forecast()} takes
#' a shock covariance matrix \code{Q} but does not support per-horizon scale
#' vectors in the current version.  If a plan containing shock-scale entries
#' is passed to \code{conditional_forecast()}, a warning is issued and the
#' shock-scale entries are ignored on the forecast side.
#'
#' @param var     Name of the exogenous shock to scale (character scalar).
#'   Must match a name in \code{dr$exo_names}.
#' @param periods Integer vector of sample periods (1-based).
#' @param scales  Numeric vector of scale factors (multiplicative on shock
#'   std), one per period.  A scalar is recycled across all periods.  Must be
#'   positive and finite.
#'
#' @return A \code{"plan_shock_scale_entry"} object for use in
#'   \code{\link{dynhr_plan}()}.
#'
#' @seealso \code{\link{dynhr_plan}}, \code{\link{plan_tune}},
#'   \code{\link{plan_condition}}
#' @export
plan_shock_scale <- function(var, periods, scales) {
  if (!is.character(var) || length(var) != 1L || nchar(var) == 0L)
    stop("plan_shock_scale: 'var' must be a non-empty character scalar.",
         call. = FALSE)

  periods <- as.integer(periods)
  if (length(periods) == 0L)
    stop(sprintf("plan_shock_scale: var '%s' has zero periods.", var), call. = FALSE)
  if (any(!is.finite(periods)))
    stop(sprintf("plan_shock_scale: var '%s' has non-finite periods.", var),
         call. = FALSE)

  scales <- as.numeric(scales)
  if (length(scales) == 1L && length(periods) > 1L)
    scales <- rep(scales, length(periods))
  if (length(scales) != length(periods))
    stop(sprintf(
      "plan_shock_scale: var '%s' has %d periods but %d scales.",
      var, length(periods), length(scales)), call. = FALSE)
  if (any(!is.finite(scales) | scales <= 0))
    stop(sprintf(
      "plan_shock_scale: var '%s' scales must be finite and positive.", var),
      call. = FALSE)

  entry <- list(var = var, periods = periods, scales = scales)
  class(entry) <- "plan_shock_scale_entry"
  entry
}


## ===========================================================================
## Plan constructor
## ===========================================================================

#' Build a unified judgments plan for a DSGE model
#'
#' Collects in-sample filter tunes, out-of-sample forecast conditions, and
#' shock-scale entries into a single \code{dynhr_plan} object.  The plan can
#' then be passed to \code{\link{run_full_estimation}()},
#' \code{\link{run_mode_finding}()}, or \code{\link{conditional_forecast}()}
#' in place of the individual \code{filter_tunes}, \code{heteroskedastic_shocks},
#' and \code{conditions} arguments.
#'
#' \strong{Reference-frame discipline}:
#' \itemize{
#'   \item \emph{In-sample entries} (\code{\link{plan_tune}},
#'     \code{\link{plan_shock_scale}}) use \strong{sample period indexing}:
#'     1-based integers relative to \code{sample_start}, or Dynare-style date
#'     literals resolved at adapter time.
#'   \item \emph{Out-of-sample entries} (\code{\link{plan_condition}}) use
#'     \strong{forecast horizon offsets}: 1-based from the terminal filter
#'     state \eqn{s_{T|T}}.
#'   \item These two reference frames are never mixed.
#' }
#'
#' \strong{Composition rules}:
#' \itemize{
#'   \item Duplicate tune entries (same var + any overlapping period) are an
#'     error.
#'   \item Duplicate shock-scale entries (same shock + any overlapping period)
#'     are an error.
#'   \item Mixed anticipated/unanticipated hard conditions in one plan are an
#'     error (the Waggoner-Zha anticipated path and the unanticipated path
#'     produce materially different shock paths).
#'   \item Shock-scale entries do not conflict with condition entries: scales
#'     operate during filtering (in-sample); conditions operate during
#'     forecasting (out-of-sample).
#' }
#'
#' \strong{TPF incompatibility}: a plan with tune or shock-scale entries
#' cannot be used with \code{likelihood = "tpf"}.  The incompatibility is
#' checked at integration time in \code{run_full_estimation()} /
#' \code{run_mode_finding()}.
#'
#' @param ...  One or more entries created by \code{\link{plan_tune}()},
#'   \code{\link{plan_condition}()}, or \code{\link{plan_shock_scale}()}.
#'
#' @return An object of class \code{"dynhr_plan"} with three sub-lists:
#'   \describe{
#'     \item{\code{$tunes}}{List of \code{plan_tune_entry} objects.}
#'     \item{\code{$conditions}}{List of \code{plan_condition_entry} objects.}
#'     \item{\code{$shock_scales}}{List of \code{plan_shock_scale_entry} objects.}
#'   }
#'
#' @examples
#' \dontrun{
#' ## Pin unobserved variable 'y' at 0.01 for periods 20-22; condition on
#' ## obs_y at horizons 1-4; scale shock eps during periods 15-18.
#' p <- dynhr_plan(
#'   plan_tune("y", 20:22, 0.01),
#'   plan_condition("obs_y", horizons = 1:4, values = 0.05),
#'   plan_shock_scale("eps", 15:18, 0.5)
#' )
#'
#' ## Use in estimation (replaces filter_tunes + heteroskedastic_shocks args):
#' result <- run_full_estimation(model = ..., data = ..., obs_vars = ...,
#'                               plan = p)
#'
#' ## Use in conditional forecast (conditions section only):
#' cf <- conditional_forecast(model, dr, Y, plan = p, horizon = 8L)
#' }
#'
#' @seealso \code{\link{plan_tune}}, \code{\link{plan_condition}},
#'   \code{\link{plan_shock_scale}},
#'   \code{\link{plan_to_filter_tunes}}, \code{\link{plan_to_conditions}},
#'   \code{\link{plan_to_shock_scale_spec}}
#' @export
dynhr_plan <- function(...) {
  entries <- list(...)
  ## Flatten a single list-of-entries argument
  if (length(entries) == 1L && is.list(entries[[1L]]) &&
      !inherits(entries[[1L]], c("plan_tune_entry", "plan_condition_entry",
                                 "plan_shock_scale_entry")))
    entries <- entries[[1L]]

  tunes        <- list()
  conditions   <- list()
  shock_scales <- list()

  for (i in seq_along(entries)) {
    e <- entries[[i]]
    if (inherits(e, "plan_tune_entry")) {
      tunes <- c(tunes, list(e))
    } else if (inherits(e, "plan_condition_entry")) {
      conditions <- c(conditions, list(e))
    } else if (inherits(e, "plan_shock_scale_entry")) {
      shock_scales <- c(shock_scales, list(e))
    } else {
      stop(sprintf(
        "dynhr_plan: argument %d is not a plan_tune(), plan_condition(), or plan_shock_scale() entry (got class '%s').",
        i, class(e)[1L]), call. = FALSE)
    }
  }

  ## ---- Composition validation -------------------------------------------

  ## Duplicate tune: same var + overlapping integer periods
  if (length(tunes) >= 2L) {
    for (i in seq_along(tunes)) {
      for (j in seq_along(tunes)) {
        if (j <= i) next
        if (tunes[[i]]$var == tunes[[j]]$var &&
            !isTRUE(tunes[[i]]$date_literals) &&
            !isTRUE(tunes[[j]]$date_literals)) {
          overlap <- intersect(tunes[[i]]$periods, tunes[[j]]$periods)
          if (length(overlap) > 0L)
            stop(sprintf(
              "dynhr_plan: duplicate tune entries for var '%s' at period(s) %s.",
              tunes[[i]]$var, paste(overlap, collapse = ", ")), call. = FALSE)
        }
      }
    }
  }

  ## Duplicate shock_scale: same var + overlapping periods
  if (length(shock_scales) >= 2L) {
    for (i in seq_along(shock_scales)) {
      for (j in seq_along(shock_scales)) {
        if (j <= i) next
        if (shock_scales[[i]]$var == shock_scales[[j]]$var) {
          overlap <- intersect(shock_scales[[i]]$periods, shock_scales[[j]]$periods)
          if (length(overlap) > 0L)
            stop(sprintf(
              "dynhr_plan: duplicate shock_scale entries for shock '%s' at period(s) %s.",
              shock_scales[[i]]$var, paste(overlap, collapse = ", ")), call. = FALSE)
        }
      }
    }
  }

  ## Mixed anticipated/unanticipated hard conditions (L3 landmine)
  if (length(conditions) >= 2L) {
    hard_conds <- Filter(function(e) e$method == "hard", conditions)
    if (length(hard_conds) >= 2L) {
      types <- vapply(hard_conds, function(e) e$type, character(1L))
      if (length(unique(types)) > 1L)
        stop(
          "dynhr_plan: cannot mix anticipated and unanticipated hard conditions in one plan. ",
          "The Waggoner-Zha anticipated path (stacked QR, R/conditional-forecast.R:675-698) ",
          "and the unanticipated path (period-by-period, R/conditional-forecast.R:699-703) ",
          "produce materially different shock paths even for identical conditions. ",
          "Supply separate plans.",
          call. = FALSE)
    }
  }

  plan <- list(tunes = tunes, conditions = conditions, shock_scales = shock_scales)
  class(plan) <- "dynhr_plan"
  plan
}


## ===========================================================================
## print / summary methods
## ===========================================================================

#' @export
print.dynhr_plan <- function(x, ...) {
  n_t <- length(x$tunes)
  n_c <- length(x$conditions)
  n_s <- length(x$shock_scales)

  cat(sprintf(
    "<dynhr_plan>  [%d tune(s) | %d condition(s) | %d shock_scale(s)]\n",
    n_t, n_c, n_s))

  if (n_t > 0L) {
    cat("\nIn-sample tunes (sample period reference frame):\n")
    for (e in x$tunes) {
      kind <- if (is.null(e$stderr)) "hard" else "soft"
      pstr <- if (length(e$periods) <= 6L) paste(e$periods, collapse = ", ")
              else paste(c(as.character(e$periods[1:5]), "..."), collapse = ", ")
      cat(sprintf("  var %-12s  periods [%s]  %s\n", e$var, pstr, kind))
    }
  }

  if (n_c > 0L) {
    cat("\nOut-of-sample conditions (forecast horizon reference frame):\n")
    for (e in x$conditions) {
      hstr <- if (length(e$horizons) <= 6L) paste(e$horizons, collapse = ", ")
              else paste(c(e$horizons[1:5], "..."), collapse = ", ")
      cat(sprintf("  var %-12s  horizons [%s]  type=%s  method=%s\n",
                  e$var, hstr, e$type, e$method))
    }
  }

  if (n_s > 0L) {
    cat("\nShock scales (sample period reference frame):\n")
    for (e in x$shock_scales) {
      pstr <- if (length(e$periods) <= 6L) paste(e$periods, collapse = ", ")
              else paste(c(e$periods[1:5], "..."), collapse = ", ")
      cat(sprintf("  shock %-10s  periods [%s]\n", e$var, pstr))
    }
  }

  invisible(x)
}


#' @export
summary.dynhr_plan <- function(object, ...) {
  n_t <- length(object$tunes)
  n_c <- length(object$conditions)
  n_s <- length(object$shock_scales)

  cat(sprintf(
    "dynhr_plan summary: %d tune(s), %d condition(s), %d shock_scale(s)\n",
    n_t, n_c, n_s))

  if (n_t > 0L) {
    vars      <- vapply(object$tunes, function(e) e$var, character(1L))
    kinds     <- vapply(object$tunes,
                        function(e) if (is.null(e$stderr)) "hard" else "soft",
                        character(1L))
    n_periods <- vapply(object$tunes, function(e) length(e$periods), integer(1L))
    cat("\nIn-sample tunes:\n")
    df <- data.frame(var = vars, type = kinds, n_periods = n_periods,
                     stringsAsFactors = FALSE)
    print(df, row.names = FALSE)
  }

  if (n_c > 0L) {
    vars    <- vapply(object$conditions, function(e) e$var, character(1L))
    types   <- vapply(object$conditions, function(e) e$type, character(1L))
    methods <- vapply(object$conditions, function(e) e$method, character(1L))
    n_hor   <- vapply(object$conditions, function(e) length(e$horizons), integer(1L))
    cat("\nOut-of-sample conditions:\n")
    df <- data.frame(var = vars, type = types, method = methods,
                     n_horizons = n_hor, stringsAsFactors = FALSE)
    print(df, row.names = FALSE)
  }

  if (n_s > 0L) {
    shocks    <- vapply(object$shock_scales, function(e) e$var, character(1L))
    n_periods <- vapply(object$shock_scales, function(e) length(e$periods), integer(1L))
    cat("\nShock scales:\n")
    df <- data.frame(shock = shocks, n_periods = n_periods,
                     stringsAsFactors = FALSE)
    print(df, row.names = FALSE)
  }

  invisible(object)
}


## ===========================================================================
## Adapters
## ===========================================================================

#' Convert a dynhr_plan to a filter_tunes_spec
#'
#' Extracts the in-sample tune entries from a \code{dynhr_plan} and assembles
#' them into a \code{filter_tunes_spec} suitable for passing to
#' \code{\link{run_full_estimation}()} or \code{\link{run_mode_finding}()} via
#' the \code{filter_tunes} argument.
#'
#' Date-literal periods are resolved via \code{sample_start} using
#' \code{.date_to_period()}.  If any tune entry uses date literals and
#' \code{sample_start} is \code{NULL}, an error is raised.
#'
#' @param plan         A \code{dynhr_plan} object.
#' @param sample_start Optional character scalar; date literal for sample
#'   period 1 (e.g. \code{"1990q1"}).  Required when any tune entry's
#'   \code{periods} contain date literals.
#'
#' @return A \code{filter_tunes_spec} object (same class as returned by
#'   \code{\link{filter_tunes}()}).  Returns an empty spec if the plan has no
#'   tune entries.
#'
#' @seealso \code{\link{dynhr_plan}}, \code{\link{filter_tunes}}
#' @export
plan_to_filter_tunes <- function(plan, sample_start = NULL) {
  stopifnot(inherits(plan, "dynhr_plan"))

  if (length(plan$tunes) == 0L)
    return(filter_tunes())

  entries <- lapply(plan$tunes, function(e) {
    periods <- e$periods
    ## Resolve date literals if needed
    if (isTRUE(e$date_literals)) {
      if (is.null(sample_start))
        stop(sprintf(
          "plan_to_filter_tunes: tune for var '%s' uses date literals but sample_start is NULL.",
          e$var), call. = FALSE)
      periods <- vapply(as.character(periods),
                        function(s) .date_to_period(s, sample_start),
                        integer(1L))
    }
    tune(e$var, periods, e$values, stderr = e$stderr)
  })

  do.call(filter_tunes, entries)
}


#' Convert a dynhr_plan to a conditions data.frame for conditional_forecast
#'
#' Extracts the out-of-sample condition entries from a \code{dynhr_plan} and
#' assembles them into a data.frame with columns \code{var}, \code{horizon},
#' \code{value}, and optionally \code{stderr}.  The consensus \code{type} and
#' \code{method} are returned as attributes.
#'
#' @param plan  A \code{dynhr_plan} object.
#'
#' @return A data.frame with columns \code{var} (character), \code{horizon}
#'   (integer), \code{value} (numeric), \code{stderr} (numeric or \code{NA}).
#'   Two attributes are set:
#'   \describe{
#'     \item{\code{attr(., "type")}}{Consensus type (\code{"unanticipated"} or
#'       \code{"anticipated"}).}
#'     \item{\code{attr(., "method")}}{Consensus method (\code{"hard"} or
#'       \code{"soft"}).}
#'   }
#'   Returns an empty data.frame (zero rows) with default type/method if the
#'   plan has no condition entries.
#'
#' @seealso \code{\link{dynhr_plan}}, \code{\link{conditional_forecast}}
#' @export
plan_to_conditions <- function(plan) {
  stopifnot(inherits(plan, "dynhr_plan"))

  if (length(plan$conditions) == 0L) {
    df <- data.frame(var = character(0), horizon = integer(0),
                     value = numeric(0), stderr = numeric(0),
                     stringsAsFactors = FALSE)
    attr(df, "type")   <- "unanticipated"
    attr(df, "method") <- "hard"
    return(df)
  }

  rows <- lapply(plan$conditions, function(e) {
    n      <- length(e$horizons)
    se_vec <- if (!is.null(e$stderr)) e$stderr else rep(NA_real_, n)
    data.frame(var     = rep(e$var, n),
               horizon = e$horizons,
               value   = e$values,
               stderr  = se_vec,
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  ## Determine consensus type / method
  ## (mixed-anticipated guard already fires in dynhr_plan())
  all_types   <- vapply(plan$conditions, function(e) e$type,   character(1L))
  all_methods <- vapply(plan$conditions, function(e) e$method, character(1L))
  consensus_type   <- if (length(unique(all_types))   == 1L) all_types[1L]   else "unanticipated"
  consensus_method <- if (length(unique(all_methods)) == 1L) all_methods[1L] else "hard"

  attr(df, "type")   <- consensus_type
  attr(df, "method") <- consensus_method
  df
}


#' Convert a dynhr_plan to a het_shocks_spec for heteroskedastic shocks
#'
#' Extracts the shock-scale entries from a \code{dynhr_plan} and assembles
#' them into a \code{het_shocks_spec} suitable for passing to
#' \code{\link{run_full_estimation}()} or \code{\link{run_mode_finding}()} via
#' the \code{heteroskedastic_shocks} argument.
#'
#' @param plan  A \code{dynhr_plan} object.
#'
#' @return A \code{het_shocks_spec} object (same class as returned by
#'   \code{\link{heteroskedastic_shocks}()}).  Returns an empty spec if the
#'   plan has no shock-scale entries.
#'
#' @seealso \code{\link{dynhr_plan}}, \code{\link{heteroskedastic_shocks}}
#' @export
plan_to_shock_scale_spec <- function(plan) {
  stopifnot(inherits(plan, "dynhr_plan"))

  if (length(plan$shock_scales) == 0L)
    return(heteroskedastic_shocks())

  entries <- lapply(plan$shock_scales, function(e)
    shock_scale_entry(e$var, e$periods, e$scales))

  do.call(heteroskedastic_shocks, entries)
}


## ===========================================================================
## TPF incompatibility check (internal helper)
## ===========================================================================

#' Check a plan for TPF incompatibility (internal)
#'
#' Called at integration entry points when a plan is supplied.  Stops with a
#' clear error if the plan contains tunes or shock scales AND the likelihood
#' is "tpf".
#'
#' @param plan       A \code{dynhr_plan} object.
#' @param likelihood Character scalar likelihood type.
#' @noRd
.check_plan_tpf_compat <- function(plan, likelihood = NULL) {
  has_tunes  <- length(plan$tunes) > 0L
  has_scales <- length(plan$shock_scales) > 0L
  if (!has_tunes && !has_scales) return(invisible(NULL))
  if (!identical(likelihood, "tpf")) return(invisible(NULL))

  what <- c(if (has_tunes) "filter tunes" else NULL,
            if (has_scales) "shock scales" else NULL)
  stop(sprintf(
    "dynhr_plan with %s is not compatible with likelihood = \"tpf\". ",
    paste(what, collapse = " and ")),
    "The Tempered Particle Filter (R/posterior.R) hard-stops on filter_tunes and ",
    "shock_scale. Remove the plan's tunes/scales or use likelihood = \"gaussian\".",
    call. = FALSE)
}
