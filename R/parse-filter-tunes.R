## R/parse-filter-tunes.R
## --------------------------------------------------------------------------
## Parser for the `filter_tunes ... end;` block (Dynare #2020 / #2030 syntax):
##
##   filter_tunes;
##     var c;
##     periods 1:4;
##     values 1.02, 1.03, 1.01, 1.00;
##
##     var k;
##     periods 5;
##     values 10.5;
##     stderr 0.01;
##   end;
##
## Each `var NAME; periods P; values V; [stderr S;]` group becomes one row of
## a tibble-like data.frame with list-columns `periods`/`values`/`stderr`.
## Hard tunes have `stderr = NULL` (NA list entry); soft tunes have a stderr
## scalar (constant across periods) or vector (one value per tuned period).
## --------------------------------------------------------------------------

#' Parse the filter_tunes block into a structured list
#'
#' @param body         Body text of the filter_tunes block (between
#'   `filter_tunes;` and `end;`).
#' @param sample_start  Optional date literal (e.g. "1990q1") giving the date
#'   of sample period 1, used to resolve date-literal periods. May be `NULL`
#'   if the block uses only integer periods.
#' @return A list with one element: `tunes`, a data.frame with columns
#'   `var` (character), `periods` (list of integer vectors), `values` (list
#'   of numeric vectors), `stderr` (list of numeric vectors or `NULL` for
#'   hard tunes). One row per `var` statement.
#' @noRd
parse_filter_tunes_block <- function(body, sample_start = NULL) {
  tunes <- list()

  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  cur_var     <- NULL
  cur_periods <- NULL
  cur_values  <- NULL
  cur_stderr  <- NULL

  flush <- function() {
    if (is.null(cur_var)) return(invisible())
    if (is.null(cur_periods))
      stop(sprintf("filter_tunes: var %s has no 'periods' statement.",
                    cur_var), call. = FALSE)
    if (is.null(cur_values))
      stop(sprintf("filter_tunes: var %s has no 'values' statement.",
                    cur_var), call. = FALSE)
    ## Scalar-value recycling: a single value is broadcast across all periods
    ## (matches Dynare behaviour for constant tunes over a range).
    if (length(cur_values) == 1L && length(cur_periods) > 1L)
      cur_values <- rep(cur_values, length(cur_periods))
    if (length(cur_values) != length(cur_periods))
      stop(sprintf(
        "filter_tunes: var %s has %d periods but %d values.",
        cur_var, length(cur_periods), length(cur_values)), call. = FALSE)
    if (!is.null(cur_stderr) && length(cur_stderr) != 1L &&
        length(cur_stderr) != length(cur_periods))
      stop(sprintf(
        "filter_tunes: var %s has %d periods but %d stderr values ",
        cur_var, length(cur_periods), length(cur_stderr)),
        call. = FALSE)
    if (!is.null(cur_stderr) && length(cur_stderr) == 1L &&
        length(cur_periods) > 1L)
      cur_stderr <- rep(cur_stderr, length(cur_periods))

    tunes[[length(tunes) + 1L]] <<- list(
      var     = cur_var,
      periods = cur_periods,
      values  = cur_values,
      stderr  = cur_stderr
    )
    cur_var     <<- NULL
    cur_periods <<- NULL
    cur_values  <<- NULL
    cur_stderr  <<- NULL
  }

  safe_eval_list <- function(spec) {
    spec <- trimws(spec)
    spec <- sub(",\\s*$", "", spec)        # trailing comma tolerance (#2030)
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
        stop("filter_tunes: 'periods' statement without a preceding 'var'.",
             call. = FALSE)
      cur_periods <- .expand_period_list(m_periods[2], sample_start)
      next
    }

    m_values <- regmatches(s, regexec(
      "^\\s*values?\\s+(.+)$", s, perl = TRUE))[[1]]
    if (length(m_values) > 0 && nchar(m_values[1]) > 0) {
      if (is.null(cur_var))
        stop("filter_tunes: 'values' statement without a preceding 'var'.",
             call. = FALSE)
      cur_values <- safe_eval_list(m_values[2])
      next
    }

    m_stderr <- regmatches(s, regexec(
      "^\\s*stderr\\s+(.+)$", s, perl = TRUE))[[1]]
    if (length(m_stderr) > 0 && nchar(m_stderr[1]) > 0) {
      if (is.null(cur_var))
        stop("filter_tunes: 'stderr' statement without a preceding 'var'.",
             call. = FALSE)
      cur_stderr <- safe_eval_list(m_stderr[2])
      next
    }

    stop(sprintf("filter_tunes: unrecognised statement '%s'.", s),
         call. = FALSE)
  }
  flush()

  if (length(tunes) == 0L) {
    empty <- data.frame(var = character(0), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$values  <- list()
    empty$stderr  <- list()
    return(list(tunes = empty))
  }

  df <- data.frame(var = vapply(tunes, function(x) x$var, character(1)),
                   stringsAsFactors = FALSE)
  df$periods <- lapply(tunes, function(x) x$periods)
  df$values  <- lapply(tunes, function(x) x$values)
  df$stderr  <- lapply(tunes, function(x) x$stderr)
  list(tunes = df)
}
