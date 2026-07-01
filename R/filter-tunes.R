## R/filter-tunes.R
## --------------------------------------------------------------------------
## Implicit-observable-expansion adapter for the filter_tunes block
## (Dynare #2020). Tuned (otherwise-unobserved) endogenous variables join the
## observables; the data matrix gains a column that is NA except at the tune
## periods, where it holds the tuned value; the measurement-error covariance
## is extended so that:
##   - hard tunes (no stderr)  -> ME variance 0 at tune periods (exact)
##   - soft tunes (stderr > 0) -> ME variance stderr^2 at tune periods
##
## Downstream: kalman_filter()'s standard ("missing data") path already drops
## NA observables per period; soft, per-period heteroskedastic ME variances
## are threaded through via the `me_extra` argument (see R/kalman-filter.R).
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1.  Constructor: filter_tunes() and tune()
## ---------------------------------------------------------------------------

#' Specify call-level filter tunes for a DSGE model
#'
#' Constructs a \code{filter_tunes_spec} object -- the same structure produced
#' by the \code{filter_tunes ... end;} block in a \code{.mod} file -- so that
#' tunes can be supplied or overridden at estimation time without editing the
#' model file.
#'
#' Each tune targets a single endogenous variable over a set of periods.
#' Hard tunes (\code{stderr = NULL}) pin the variable exactly at
#' \code{values}; soft tunes (\code{stderr > 0}) add Gaussian measurement
#' noise with standard deviation \code{stderr} at the tuned periods.
#'
#' Scalar recycling: a scalar \code{values} or \code{stderr} is broadcast
#' across all periods when \code{periods} has length > 1, matching the
#' \code{filter_tunes} block behaviour.
#'
#' @param ...  One or more tune entries created by \code{\link{tune}()}.
#'
#' @return An object of class \code{"filter_tunes_spec"}: a list with one
#'   element \code{$tunes}, a data.frame with columns \code{var} (character),
#'   \code{periods} (list of integer vectors), \code{values} (list of numeric
#'   vectors), \code{stderr} (list of numeric vectors or \code{NULL} for hard
#'   tunes). This is identical in structure to the \code{$filter_tunes} slot
#'   produced by \code{\link{parse_mod}()}.
#'
#' Pass this object to the \code{filter_tunes} argument of
#' \code{\link{run_full_estimation}()} or \code{\link{run_mode_finding}()}
#' to override the \code{.mod}-file block. Pass \code{FALSE} to ignore the
#' mod-file block entirely.
#'
#' @examples
#' \dontrun{
#' ## Hard tune: pin z at 0.1 for periods 3-5
#' tunes <- filter_tunes(
#'   tune("z", 3:5, 0.1)
#' )
#'
#' ## Mixed: hard tune on z, soft tune on a (per-period value, global stderr)
#' tunes <- filter_tunes(
#'   tune("z", 3:5, 0.1),
#'   tune("a", 7L, -0.2, stderr = 0.05)
#' )
#'
#' ## Override the mod-file block at estimation time:
#' result <- run_full_estimation(
#'   model        = my_model,
#'   data         = Y,
#'   obs_vars     = c("y", "pi", "r"),
#'   filter_tunes = tunes
#' )
#'
#' ## Disable the mod-file block entirely:
#' result <- run_full_estimation(
#'   model        = my_model,
#'   data         = Y,
#'   obs_vars     = c("y", "pi", "r"),
#'   filter_tunes = FALSE
#' )
#' }
#'
#' @seealso \code{\link{tune}}, \code{\link{run_full_estimation}},
#'   \code{\link{run_mode_finding}}
#' @export
filter_tunes <- function(...) {
  entries <- list(...)
  ## Flatten a single list-of-entries argument (same idiom as system_prior_spec).
  if (length(entries) == 1L && is.list(entries[[1L]]) &&
      !inherits(entries[[1L]], "tune_entry"))
    entries <- entries[[1L]]

  for (i in seq_along(entries)) {
    e <- entries[[i]]
    if (!inherits(e, "tune_entry"))
      stop(sprintf(
        "filter_tunes: argument %d is not a tune() entry (got %s).",
        i, class(e)[1L]), call. = FALSE)
  }

  if (length(entries) == 0L) {
    empty <- data.frame(var = character(0L), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$values  <- list()
    empty$stderr  <- list()
    out <- list(tunes = empty)
    class(out) <- "filter_tunes_spec"
    return(out)
  }

  df <- data.frame(
    var = vapply(entries, function(e) e$var, character(1L)),
    stringsAsFactors = FALSE
  )
  df$periods <- lapply(entries, function(e) e$periods)
  df$values  <- lapply(entries, function(e) e$values)
  df$stderr  <- lapply(entries, function(e) e$stderr)

  out <- list(tunes = df)
  class(out) <- "filter_tunes_spec"
  out
}


#' Specify a single tune entry
#'
#' Helper for \code{\link{filter_tunes}()} defining one tune.
#'
#' @param var     Name of the endogenous variable to tune (character scalar).
#' @param periods Integer vector of sample periods (1-based).
#' @param values  Numeric vector of tuned values, one per period.  A scalar
#'   is recycled across all periods.
#' @param stderr  Numeric scalar (or vector, one per period) giving the
#'   standard deviation of measurement noise at the tuned periods (soft tune).
#'   \code{NULL} (default) = hard tune: variable is pinned exactly.
#'
#' @return A \code{"tune_entry"} list for use inside \code{\link{filter_tunes}()}.
#'
#' @examples
#' \dontrun{
#' filter_tunes(
#'   tune("z", 3:5, 0.1),
#'   tune("a", 7L, -0.2, stderr = 0.05)
#' )
#' }
#'
#' @seealso \code{\link{filter_tunes}}
#' @export
tune <- function(var, periods, values, stderr = NULL) {
  if (!is.character(var) || length(var) != 1L || nchar(var) == 0L)
    stop("tune: 'var' must be a non-empty character scalar.", call. = FALSE)

  periods <- as.integer(periods)
  if (length(periods) == 0L)
    stop(sprintf("tune: var '%s' has zero periods.", var), call. = FALSE)
  if (any(!is.finite(periods)))
    stop(sprintf("tune: var '%s' has non-finite periods.", var), call. = FALSE)

  values <- as.numeric(values)
  if (length(values) == 1L && length(periods) > 1L)
    values <- rep(values, length(periods))
  if (length(values) != length(periods))
    stop(sprintf(
      "tune: var '%s' has %d periods but %d values.",
      var, length(periods), length(values)), call. = FALSE)

  if (!is.null(stderr)) {
    stderr <- as.numeric(stderr)
    if (length(stderr) == 1L && length(periods) > 1L)
      stderr <- rep(stderr, length(periods))
    if (length(stderr) != length(periods))
      stop(sprintf(
        "tune: var '%s' has %d periods but %d stderr values.",
        var, length(periods), length(stderr)), call. = FALSE)
    if (any(!is.finite(stderr) | stderr < 0))
      stop(sprintf("tune: var '%s' stderr must be finite and non-negative.", var),
           call. = FALSE)
  }

  entry <- list(var = var, periods = periods, values = values, stderr = stderr)
  class(entry) <- "tune_entry"
  entry
}


## ---------------------------------------------------------------------------
## 2.  Internal: resolve call-level filter_tunes arg for estimation entry points
## ---------------------------------------------------------------------------

#' Resolve call-level filter_tunes argument for estimation
#'
#' Given the \code{filter_tunes} user argument and the parsed model, return
#' the model with \code{$filter_tunes} set appropriately:
#'   \code{NULL}  -> use mod-file block unchanged
#'   \code{FALSE} -> replace with empty tunes (disable)
#'   \code{filter_tunes_spec} -> override block with call-level spec
#'
#' @param model         Parsed \code{dynhr_mod}.
#' @param filter_tunes  User-supplied argument.
#' @return Possibly-modified \code{model} list.
#' @noRd
.resolve_filter_tunes <- function(model, filter_tunes) {
  if (is.null(filter_tunes))
    return(model)
  if (isFALSE(filter_tunes)) {
    empty <- data.frame(var = character(0L), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$values  <- list()
    empty$stderr  <- list()
    model$filter_tunes <- list(tunes = empty)
    return(model)
  }
  if (inherits(filter_tunes, "filter_tunes_spec")) {
    model$filter_tunes <- filter_tunes
    return(model)
  }
  stop(paste0(
    "filter_tunes must be NULL (use mod-file block), FALSE (disable), ",
    "or a filter_tunes_spec object from filter_tunes()."),
    call. = FALSE)
}


## ---------------------------------------------------------------------------
## 3.  .expand_observables_for_tunes (internal adapter)
## ---------------------------------------------------------------------------

#' Expand observables for filter_tunes
#'
#' @param model     A `dynhr_mod` (from `parse_mod()`); uses
#'   `model$filter_tunes$tunes` and `model$var_names`.
#' @param obs_names Character vector of observed variable names already in
#'   `Y` (in column order).
#' @param Y         `T x n_obs` data matrix (or data.frame), columns matching
#'   `obs_names`.
#' @return If `model$filter_tunes$tunes` has zero rows, returns
#'   `list(obs_vars = obs_names, Y = as.matrix(Y), me_extra = NULL)`
#'   unchanged (cheap no-op). Otherwise a list with:
#'   \describe{
#'     \item{obs_vars}{`obs_names` with the tuned variables appended (in
#'       first-appearance order; a variable with multiple tune rows appears
#'       once).}
#'     \item{Y}{`T x length(obs_vars)` matrix; the appended columns are `NA`
#'       except at the tune periods, which hold the tuned values.}
#'     \item{me_extra}{`length(obs_vars) x T` matrix of additional
#'       measurement-error variances to add to the diagonal of `H` at each
#'       period (0 for hard tunes and for all original observables; `stderr^2`
#'       for soft tunes at their tune periods). Pass to
#'       `kalman_filter(..., me_extra = me_extra)` and
#'       `kalman_smoother(..., me_extra = me_extra)`.}
#'   }
#' @noRd
.expand_observables_for_tunes <- function(model, obs_names, Y) {
  Y <- as.matrix(Y)
  tunes <- model$filter_tunes$tunes
  if (is.null(tunes) || nrow(tunes) == 0L)
    return(list(obs_vars = obs_names, Y = Y, me_extra = NULL))

  n_T <- nrow(Y)

  ## ---- Guards -----------------------------------------------------------
  for (i in seq_len(nrow(tunes))) {
    v  <- tunes$var[i]
    pr <- tunes$periods[[i]]
    if (!(v %in% model$var_names))
      stop(sprintf(
        "filter_tunes: '%s' is not an endogenous variable in this model.",
        v), call. = FALSE)
    if (v %in% obs_names)
      stop(sprintf(
        "filter_tunes: '%s' is already an observed variable; tunes apply only to unobserved endogenous variables.",
        v), call. = FALSE)
    if (any(pr < 1L) || any(pr > n_T))
      stop(sprintf(
        "filter_tunes: var %s has tune periods outside the sample (1..%d): %s",
        v, n_T, paste(pr[pr < 1L | pr > n_T], collapse = ", ")), call. = FALSE)
  }

  ## ---- New observables: one column per distinct tuned var ----------------
  tuned_vars <- unique(tunes$var)
  obs_vars   <- c(obs_names, tuned_vars)
  n_obs_orig <- length(obs_names)
  n_obs_new  <- length(obs_vars)

  Y_new <- matrix(NA_real_, nrow = n_T, ncol = n_obs_new)
  Y_new[, seq_len(n_obs_orig)] <- Y
  colnames(Y_new) <- obs_vars

  me_extra <- matrix(0, nrow = n_obs_new, ncol = n_T)

  for (i in seq_len(nrow(tunes))) {
    v       <- tunes$var[i]
    col_idx <- n_obs_orig + match(v, tuned_vars)
    pr      <- tunes$periods[[i]]
    vals    <- tunes$values[[i]]
    se      <- tunes$stderr[[i]]

    Y_new[pr, col_idx] <- vals
    if (!is.null(se))
      me_extra[col_idx, pr] <- se^2
  }

  list(obs_vars = obs_vars, Y = Y_new, me_extra = me_extra)
}
