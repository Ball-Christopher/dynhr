## R/time-orientation.R
## --------------------------------------------------------------------------
## ONE convention for every time-indexed matrix the package takes or returns:
## ROWS ARE PERIODS (T x n). Covariance stacks are the one exception (n x n x T,
## a stack of matrices); path / draw / quantile arrays put periods first.
##
## Inputs used to be documented n x T in some places and T x n in others, and
## several sites GUESSED the orientation from the shape (`if (nrow(Y) != n_obs)
## Y <- t(Y)`): a sample whose length equalled the number of observables was
## then read transposed -- kalman_filter() scored a correctly oriented 2 x 2
## sample at -6639.3 instead of -1926.7. .orient_time_matrix() is the single
## rule every entry point uses instead:
##   * names decide: columns named by the series -> T x n as given; rows named
##     by the series -> the old n x T orientation, transposed (with a one-time
##     message of class dynhr_message_orientation_transposed);
##   * unnamed: an unambiguous T x n shape is accepted; the old n x T shape is
##     an error with a t() hint (never silently transposed); a square matrix
##     is an error asking for column names.
## --------------------------------------------------------------------------

#' Resolve a time-indexed input to T x n (rows = periods)
#'
#' @param x Numeric matrix, data frame or (when `n_series == 1`) vector.
#' @param series Character vector of the series (column) names the caller
#'   expects, in the caller's order; or an integer count when the series are
#'   unnamed.
#' @param arg Argument name, for messages.
#' @param caller Function name, for messages.
#' @param subset Logical: may the matrix carry only SOME of the series (by
#'   name)? The result then has just those columns; the caller expands.
#'   Unnamed matrices must always carry every series.
#' @param n_periods Optional expected number of periods (rows); checked when
#'   given.
#' @return Numeric T x n matrix; when `series` is character its columns are
#'   named (and ordered as `series`, or as the subset found).
#' @noRd
.orient_time_matrix <- function(x, series, arg, caller, subset = FALSE,
                                n_periods = NULL) {
  nm <- if (is.character(series)) series else NULL
  n  <- if (is.character(series)) length(series) else as.integer(series)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    if (n != 1L)
      .dynhr_abort(caller, ": `", arg, "` must be a matrix with one column ",
                   "per series (T x ", n, ", rows = periods).",
                   class = "dynhr_error_orientation")
    x <- matrix(x, ncol = 1L, dimnames = list(NULL, nm))
  }
  if (!is.numeric(x) && !is.logical(x))
    .dynhr_abort(caller, ": `", arg, "` must be numeric.",
                 class = "dynhr_error_orientation")
  storage.mode(x) <- "double"
  rn <- rownames(x); cn <- colnames(x)
  named_by <- function(labels)
    !is.null(nm) && !is.null(labels) && length(labels) > 0L &&
      !anyNA(labels) && all(labels %in% nm) && !anyDuplicated(labels) &&
      (subset || length(labels) == n)
  col_hit <- named_by(cn)
  row_hit <- named_by(rn)
  if (col_hit && row_hit && !identical(nrow(x), ncol(x)))
    row_hit <- FALSE
  if (col_hit && row_hit)
    .dynhr_abort(caller, ": `", arg, "` has the series names on BOTH its rows ",
                 "and its columns, so its orientation is ambiguous. Use column ",
                 "names only (T x n, rows = periods).",
                 class = "dynhr_error_orientation")
  if (!col_hit && row_hit) {
    .dynhr_inform(caller, ": `", arg, "` is named by series on its ROWS (the ",
                  "old n x T orientation); it was transposed. The convention ",
                  "is T x n -- rows are periods, columns the named series.",
                  once = TRUE, key = paste0("orientation:", caller, ":", arg),
                  class = "dynhr_message_orientation_transposed")
    x <- t(x)
    col_hit <- TRUE
  }
  if (col_hit) {
    keep <- if (is.null(nm)) colnames(x) else nm[nm %in% colnames(x)]
    x <- x[, keep, drop = FALSE]
  } else {
    if (!is.null(cn) && !is.null(nm) && length(cn) && !all(cn %in% nm)) {
      bad <- setdiff(cn, nm)
      .dynhr_abort(caller, ": `", arg, "` has column names that are not ",
                   "series of this call (", paste(utils::head(bad, 5),
                   collapse = ", "), if (length(bad) > 5) ", ..." else "",
                   "); expected ", paste(nm, collapse = ", "), ".",
                   class = "dynhr_error_orientation")
    }
    if (nrow(x) == ncol(x) && ncol(x) == n && n > 1L)
      .dynhr_abort(caller, ": `", arg, "` is square (", n, " x ", n, "), so ",
                   "its orientation cannot be read from its shape. Give it ",
                   "column names", if (!is.null(nm)) paste0(" (",
                   paste(nm, collapse = ", "), ")") else "",
                   ": the convention is T x n, rows = periods.",
                   class = "dynhr_error_orientation")
    if (ncol(x) != n) {
      hint <- if (nrow(x) == n)
        paste0(" It looks like the old n x T orientation: transpose it ",
               "with t() (the convention is T x n, rows = periods).") else ""
      .dynhr_abort(caller, ": `", arg, "` has ", ncol(x), " column(s); ",
                   "expected ", n, " (one per series, T x n, rows = periods).",
                   hint, class = "dynhr_error_orientation")
    }
    if (!is.null(nm)) colnames(x) <- nm
  }
  if (!is.null(n_periods) && nrow(x) != n_periods)
    .dynhr_abort(caller, ": `", arg, "` has ", nrow(x), " row(s); expected ",
                 n_periods, " (one per period).",
                 class = "dynhr_error_orientation")
  .check_time_values(x, arg, caller)
  x
}

#' Refuse values a time-indexed input cannot mean
#'
#' Observation data (`data` / `Y`): NA -- and NaN, MATLAB / Dynare / IRIS's
#' missing-value marker, which data exported from them carry (the IRIS fixture
#' in test-trend-balanced-growth.R) -- is a missing observation; Inf and -Inf
#' are not. The filters' `!is.finite()` test used to treat them as missing too
#' (same loglik as NA) while the smoother returned NaN.
#' `me_extra` (variances) must be finite and non-negative, `shock_scale`
#' finite: an NA / NaN there used to fail deep inside with "missing value
#' where TRUE/FALSE needed", and a negative variance gave loglik -Inf.
#' @noRd
.check_time_values <- function(x, arg, caller) {
  if (arg %in% c("data", "Y")) {
    if (!is.double(x) || !anyNA(x) && all(is.finite(x))) return(invisible(x))
    bad <- is.infinite(x)
    if (any(bad)) {
      k <- which(bad, arr.ind = TRUE)[1L, ]
      .dynhr_abort(caller, ": `", arg, "` has ", sum(bad), " non-finite ",
                   "value(s) (first: ", format(x[k[1L], k[2L]]), " at period ",
                   k[1L], ", column ", colnames(x)[k[2L]] %||% k[2L], "). Only ",
                   "NA (or NaN) marks a missing observation; Inf / -Inf are refused.",
                   class = "dynhr_error_nonfinite_data")
    }
  } else if (arg %in% c("me_extra", "shock_scale")) {
    if (!all(is.finite(x)))
      .dynhr_abort(caller, ": `", arg, "` must be finite everywhere (it has ",
                   sum(!is.finite(x)), " NA / NaN / Inf value(s)).",
                   class = "dynhr_error_bad_argument")
    if (identical(arg, "me_extra") && any(x < 0))
      .dynhr_abort(caller, ": `me_extra` holds measurement-error VARIANCES and ",
                   "must be >= 0 (minimum ", min(x), ").",
                   class = "dynhr_error_bad_argument")
  }
  invisible(x)
}
