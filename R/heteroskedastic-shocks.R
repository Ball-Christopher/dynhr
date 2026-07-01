## R/heteroskedastic-shocks.R
## --------------------------------------------------------------------------
## Call-level constructors and internal helpers for the heteroskedastic_shocks
## block (time-varying shock standard deviations).
##
## Syntax in .mod files:
##   heteroskedastic_shocks;
##     var e_a; periods 5, 6, 7; scales 0.1, 0.1, 0.5;
##     var e_r; periods 2020q1:2020q4; scales 3.0;
##   end;
##
## The `scales` values are MULTIPLICATIVE factors on shock standard deviation,
## so the effective shock covariance at period t is:
##   Sigma_e_t = diag(s_t) Sigma_e diag(s_t)
## where s_t is the column of the shock_scale matrix for period t (1 = baseline).
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1.  Constructor: heteroskedastic_shocks() and shock_scale_entry()
## ---------------------------------------------------------------------------

#' Specify call-level heteroskedastic shocks for a DSGE model
#'
#' Constructs a \code{het_shocks_spec} object -- the same structure produced by
#' the \code{heteroskedastic_shocks ... end;} block in a \code{.mod} file -- so
#' that time-varying shock standard deviations can be supplied or overridden at
#' estimation time without editing the model file.
#'
#' Each entry targets a single exogenous shock over a set of periods.  The
#' \code{scales} values are multiplicative factors on the shock standard
#' deviation (1 = baseline).  A scalar \code{scales} is broadcast across all
#' periods when \code{periods} has length > 1.
#'
#' @param ...  One or more entries created by \code{\link{shock_scale_entry}()}.
#'
#' @return An object of class \code{"het_shocks_spec"}: a list with one element
#'   \code{$scales}, a data.frame with columns \code{var} (character),
#'   \code{periods} (list of integer vectors), \code{scales} (list of numeric
#'   vectors).  This is identical in structure to the \code{$heteroskedastic_shocks}
#'   slot produced by \code{\link{parse_mod}()}.
#'
#' Pass this object to the \code{heteroskedastic_shocks} argument of
#' \code{\link{run_full_estimation}()} or \code{\link{run_mode_finding}()}
#' to override the \code{.mod}-file block.  Pass \code{FALSE} to ignore the
#' mod-file block entirely.
#'
#' @examples
#' \dontrun{
#' ## Scale e_a by 0.1 during periods 5-7
#' hs <- heteroskedastic_shocks(
#'   shock_scale_entry("e_a", 5:7, 0.1)
#' )
#'
#' ## Multiple shocks, scalar recycling
#' hs <- heteroskedastic_shocks(
#'   shock_scale_entry("e_a", 5:7, c(0.1, 0.1, 0.5)),
#'   shock_scale_entry("e_r", 10:13, 3.0)
#' )
#'
#' ## Override the mod-file block at estimation time:
#' result <- run_full_estimation(
#'   model               = my_model,
#'   data                = Y,
#'   obs_vars            = c("y", "pi", "r"),
#'   heteroskedastic_shocks = hs
#' )
#'
#' ## Disable the mod-file block entirely:
#' result <- run_full_estimation(
#'   model               = my_model,
#'   data                = Y,
#'   obs_vars            = c("y", "pi", "r"),
#'   heteroskedastic_shocks = FALSE
#' )
#' }
#'
#' @seealso \code{\link{shock_scale_entry}}, \code{\link{run_full_estimation}},
#'   \code{\link{run_mode_finding}}
#' @export
heteroskedastic_shocks <- function(...) {
  entries <- list(...)
  ## Flatten a single list-of-entries argument (same idiom as filter_tunes).
  if (length(entries) == 1L && is.list(entries[[1L]]) &&
      !inherits(entries[[1L]], "shock_scale_entry"))
    entries <- entries[[1L]]

  for (i in seq_along(entries)) {
    e <- entries[[i]]
    if (!inherits(e, "shock_scale_entry"))
      stop(sprintf(
        "heteroskedastic_shocks: argument %d is not a shock_scale_entry() (got %s).",
        i, class(e)[1L]), call. = FALSE)
  }

  if (length(entries) == 0L) {
    empty <- data.frame(var = character(0L), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$scales  <- list()
    out <- list(scales = empty)
    class(out) <- "het_shocks_spec"
    return(out)
  }

  df <- data.frame(
    var = vapply(entries, function(e) e$var, character(1L)),
    stringsAsFactors = FALSE
  )
  df$periods <- lapply(entries, function(e) e$periods)
  df$scales  <- lapply(entries, function(e) e$scales)

  out <- list(scales = df)
  class(out) <- "het_shocks_spec"
  out
}


#' Specify a single shock-scale entry
#'
#' Helper for \code{\link{heteroskedastic_shocks}()} defining the scale factors
#' for one exogenous shock over a set of periods.
#'
#' @param var     Name of the exogenous shock (character scalar).
#' @param periods Integer vector of sample periods (1-based).
#' @param scales  Numeric vector of scale factors (multiplicative on shock std),
#'   one per period.  A scalar is recycled across all periods.  Must be
#'   positive and finite.
#'
#' @return A \code{"shock_scale_entry"} list for use inside
#'   \code{\link{heteroskedastic_shocks}()}.
#'
#' @examples
#' \dontrun{
#' heteroskedastic_shocks(
#'   shock_scale_entry("e_a", 5:7, 0.1),
#'   shock_scale_entry("e_r", 10:13, 3.0)
#' )
#' }
#'
#' @seealso \code{\link{heteroskedastic_shocks}}
#' @export
shock_scale_entry <- function(var, periods, scales) {
  if (!is.character(var) || length(var) != 1L || nchar(var) == 0L)
    stop("shock_scale_entry: 'var' must be a non-empty character scalar.",
         call. = FALSE)

  periods <- as.integer(periods)
  if (length(periods) == 0L)
    stop(sprintf("shock_scale_entry: var '%s' has zero periods.", var),
         call. = FALSE)
  if (any(!is.finite(periods)))
    stop(sprintf("shock_scale_entry: var '%s' has non-finite periods.", var),
         call. = FALSE)

  scales <- as.numeric(scales)
  if (length(scales) == 1L && length(periods) > 1L)
    scales <- rep(scales, length(periods))
  if (length(scales) != length(periods))
    stop(sprintf(
      "shock_scale_entry: var '%s' has %d periods but %d scales.",
      var, length(periods), length(scales)), call. = FALSE)
  if (any(!is.finite(scales) | scales <= 0))
    stop(sprintf("shock_scale_entry: var '%s' scales must be finite and positive.",
                 var), call. = FALSE)

  entry <- list(var = var, periods = periods, scales = scales)
  class(entry) <- "shock_scale_entry"
  entry
}


## ---------------------------------------------------------------------------
## 2.  Internal: resolve call-level heteroskedastic_shocks arg
## ---------------------------------------------------------------------------

#' Resolve call-level heteroskedastic_shocks argument for estimation
#'
#' @param model                 Parsed \code{dynhr_mod}.
#' @param heteroskedastic_shocks User-supplied argument.
#' @return Possibly-modified \code{model} list.
#' @noRd
.resolve_heteroskedastic_shocks <- function(model, heteroskedastic_shocks) {
  if (is.null(heteroskedastic_shocks))
    return(model)
  if (isFALSE(heteroskedastic_shocks)) {
    empty <- data.frame(var = character(0L), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$scales  <- list()
    model$heteroskedastic_shocks <- list(scales = empty)
    return(model)
  }
  if (inherits(heteroskedastic_shocks, "het_shocks_spec")) {
    model$heteroskedastic_shocks <- heteroskedastic_shocks
    return(model)
  }
  stop(paste0(
    "heteroskedastic_shocks must be NULL (use mod-file block), FALSE (disable), ",
    "or a het_shocks_spec object from heteroskedastic_shocks()."),
    call. = FALSE)
}


## ---------------------------------------------------------------------------
## 3.  .build_shock_scale_matrix (internal adapter)
## ---------------------------------------------------------------------------

#' Build the shock_scale matrix from model$heteroskedastic_shocks
#'
#' @param model     A \code{dynhr_mod}; uses
#'   \code{model$heteroskedastic_shocks$scales} and checks against
#'   \code{dr$exo_names} (Landmine 9: row order = dr$exo_names, NOT
#'   model$varexo_names).
#' @param exo_names Character vector of shock names in dr$exo_names order.
#' @param n_T       Number of sample periods.
#' @return If \code{model$heteroskedastic_shocks} is NULL or has zero rows,
#'   returns \code{NULL} (identity -- no heteroskedasticity). Otherwise an
#'   \code{n_exo x T} matrix of scale factors (1 = baseline).
#' @noRd
.build_shock_scale_matrix <- function(model, exo_names, n_T) {
  spec <- model$heteroskedastic_shocks
  if (is.null(spec)) return(NULL)
  df <- spec$scales
  if (is.null(df) || nrow(df) == 0L) return(NULL)

  n_exo <- length(exo_names)
  mat   <- matrix(1, nrow = n_exo, ncol = n_T)
  rownames(mat) <- exo_names

  for (i in seq_len(nrow(df))) {
    v  <- df$var[i]
    pr <- df$periods[[i]]
    sc <- df$scales[[i]]
    ri <- match(v, exo_names)
    if (is.na(ri))
      stop(sprintf(
        ".build_shock_scale_matrix: shock '%s' in heteroskedastic_shocks is not ",
        "in dr$exo_names (%s).",
        v, paste(exo_names, collapse = ", ")), call. = FALSE)
    ## Clamp periods to [1, n_T] with a warning on out-of-range.
    bad <- pr < 1L | pr > n_T
    if (any(bad))
      warning(sprintf(
        ".build_shock_scale_matrix: shock '%s' has periods outside [1, %d]: %s (ignored).",
        v, n_T, paste(pr[bad], collapse = ", ")), call. = FALSE)
    ok_pr <- pr[!bad]; ok_sc <- sc[!bad]
    if (length(ok_pr) > 0L) mat[ri, ok_pr] <- ok_sc
  }

  mat
}
