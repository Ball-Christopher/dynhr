## R/parse-heteroskedastic-shocks.R
## --------------------------------------------------------------------------
## Parser for the `heteroskedastic_shocks ... end;` block:
##
##   heteroskedastic_shocks;
##     var e_a; periods 5, 6, 7; scales 0.1, 0.1, 0.5;
##     var e_r; periods 2020q1:2020q4; scales 3.0;
##   end;
##
## Structure mirrors parse-filter-tunes.R:34-155.
## Each `var NAME; periods P; scales S;` group becomes one row of a
## data.frame with list-columns `periods`/`scales`.
## `scales` values are MULTIPLICATIVE factors on shock std (1 = baseline).
## Scalar recycling: a single value is broadcast across all periods.
## --------------------------------------------------------------------------

#' Parse the heteroskedastic_shocks block into a structured list
#'
#' @param body         Body text of the heteroskedastic_shocks block (between
#'   \code{heteroskedastic_shocks;} and \code{end;}).
#' @param sample_start Optional date literal (e.g. "1990q1") giving the date
#'   of sample period 1, used to resolve date-literal periods. May be
#'   \code{NULL} if the block uses only integer periods.
#' @return A list with one element: \code{scales}, a data.frame with columns
#'   \code{var} (character), \code{periods} (list of integer vectors),
#'   \code{scales} (list of numeric vectors). One row per \code{var} statement.
#' @noRd
parse_heteroskedastic_shocks_block <- function(body, sample_start = NULL) {
  entries <- list()

  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  cur_var     <- NULL
  cur_periods <- NULL
  cur_scales  <- NULL

  flush <- function() {
    if (is.null(cur_var)) return(invisible())
    if (is.null(cur_periods))
      stop(sprintf("heteroskedastic_shocks: var %s has no 'periods' statement.",
                   cur_var), call. = FALSE)
    if (is.null(cur_scales))
      stop(sprintf("heteroskedastic_shocks: var %s has no 'scales' statement.",
                   cur_var), call. = FALSE)
    ## Scalar-value recycling: a single scale is broadcast across all periods.
    if (length(cur_scales) == 1L && length(cur_periods) > 1L)
      cur_scales <- rep(cur_scales, length(cur_periods))
    if (length(cur_scales) != length(cur_periods))
      stop(sprintf(
        "heteroskedastic_shocks: var %s has %d periods but %d scales.",
        cur_var, length(cur_periods), length(cur_scales)), call. = FALSE)
    if (any(!is.finite(cur_scales) | cur_scales <= 0))
      stop(sprintf(
        "heteroskedastic_shocks: var %s scales must be finite and positive.",
        cur_var), call. = FALSE)

    entries[[length(entries) + 1L]] <<- list(
      var     = cur_var,
      periods = cur_periods,
      scales  = cur_scales
    )
    cur_var     <<- NULL
    cur_periods <<- NULL
    cur_scales  <<- NULL
  }

  safe_eval_list <- function(spec) {
    spec <- trimws(spec)
    spec <- sub(",\\s*$", "", spec)        # trailing comma tolerance
    toks <- strsplit(spec, "[,\\s]+", perl = TRUE)[[1]]
    toks <- toks[nchar(toks) > 0]
    vapply(toks, function(tok) {
      val <- tryCatch(eval(parse(text = tok)), error = function(e) NA_real_)
      as.numeric(val)
    }, numeric(1), USE.NAMES = FALSE)
  }

  for (s in stmts) {
    if (grepl("^\\s*end\\s*$", s, ignore.case = TRUE)) next

    m_var <- regmatches(s, regexec(
      "^\\s*var\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*$", s, perl = TRUE))[[1]]
    if (length(m_var) > 0 && nchar(m_var[1]) > 0) {
      flush()
      cur_var <- m_var[2]
      next
    }

    m_periods <- regmatches(s, regexec(
      "^\\s*periods?\\s+(.+)$", s, perl = TRUE))[[1]]
    if (length(m_periods) > 0 && nchar(m_periods[1]) > 0) {
      if (is.null(cur_var))
        stop("heteroskedastic_shocks: 'periods' statement without a preceding 'var'.",
             call. = FALSE)
      cur_periods <- .expand_period_list(m_periods[2], sample_start)
      next
    }

    m_scales <- regmatches(s, regexec(
      "^\\s*scales?\\s+(.+)$", s, perl = TRUE))[[1]]
    if (length(m_scales) > 0 && nchar(m_scales[1]) > 0) {
      if (is.null(cur_var))
        stop("heteroskedastic_shocks: 'scales' statement without a preceding 'var'.",
             call. = FALSE)
      cur_scales <- safe_eval_list(m_scales[2])
      next
    }

    stop(sprintf("heteroskedastic_shocks: unrecognised statement '%s'.", s),
         call. = FALSE)
  }
  flush()

  if (length(entries) == 0L) {
    empty <- data.frame(var = character(0), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$scales  <- list()
    return(list(scales = empty))
  }

  df <- data.frame(var = vapply(entries, function(x) x$var, character(1)),
                   stringsAsFactors = FALSE)
  df$periods <- lapply(entries, function(x) x$periods)
  df$scales  <- lapply(entries, function(x) x$scales)
  list(scales = df)
}
