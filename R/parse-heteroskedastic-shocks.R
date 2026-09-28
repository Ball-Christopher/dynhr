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
#' @param first_obs Dataset row of sample period 1 (Dynare's
#'   \code{estimation(first_obs=)}, default 1), or \code{NA} when
#'   \code{first_obs} was given as a date, so that its row is unknown.
#' @return A list with one element: \code{scales}, a data.frame with columns
#'   \code{var} (character), \code{periods} (list of integer vectors),
#'   \code{scales} (list of numeric vectors). One row per \code{var} statement.
#'   Periods are SAMPLE periods (1 = the first row of the data passed to the
#'   estimation).
#'
#' @details Period semantics follow Dynare 7 (reference manual,
#' \code{heteroskedastic_shocks}; \code{dynare_estimation_init.m}):
#' \itemize{
#'   \item an integer period indexes the ORIGINAL dataset, so sample period =
#'     \code{p - first_obs + 1};
#'   \item a date (\code{2020Q1}, \code{2020Q1:2020Q4}; Dynare 7) indexes the
#'     dataset's dates, so sample period = \code{d - date(first_obs) + 1}.
#'     The data's date index is known here only when \code{first_obs} is
#'     itself a date; otherwise a date aborts with class
#'     \code{dynhr_error_mod_date_unresolved} (as \code{observation_trends}
#'     does for a date-valued \code{first_obs}).  Likewise an integer period
#'     cannot be placed when \code{first_obs} is a date.
#' }
#' @noRd
parse_heteroskedastic_shocks_block <- function(body, sample_start = NULL,
                                               first_obs = 1L) {
  entries <- list()

  ## One `periods` statement -> sample periods (see @details).
  expand_periods <- function(spec) {
    spec <- sub(",\\s*$", "", trimws(spec))
    toks <- strsplit(spec, "[,\\s]+", perl = TRUE)[[1]]
    toks <- toks[nchar(toks) > 0L]
    out  <- integer(0)
    for (tok in toks) {
      parts  <- trimws(strsplit(tok, ":", fixed = TRUE)[[1]])
      is_int <- grepl("^-?[0-9]+$", parts)
      if (all(is_int)) {
        if (is.na(first_obs))
          .dynhr_abort(
            "heteroskedastic_shocks: integer period '", tok, "' indexes the ",
            "original dataset, but estimation(first_obs=) is a date, so the ",
            "dataset row of the first observation is unknown without the data ",
            "file. Give the periods as dates, or first_obs as an integer.",
            class = c("dynhr_error_mod_date_unresolved",
                      "dynhr_error_mod_syntax"))
        out <- c(out, .expand_period_token(tok, context = "heteroskedastic_shocks") -
                   as.integer(first_obs) + 1L)
      } else if (any(is_int)) {
        .dynhr_abort("heteroskedastic_shocks: period range '", tok, "' mixes ",
                     "an integer and a date.", class = "dynhr_error_mod_syntax")
      } else {
        out <- c(out, .expand_period_token(
          tok, sample_start, context = "heteroskedastic_shocks",
          no_start_why = paste0(
            "A date is resolved against the data's date index, which dynhr ",
            "knows only from a date-valued estimation(first_obs=...) (the ",
            "date of the first observation).")))
      }
    }
    out
  }

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

  sandbox <- .dynhr_sandbox_env()
  safe_eval_list <- function(spec) {
    spec <- trimws(spec)
    spec <- sub(",\\s*$", "", spec)        # trailing comma tolerance
    toks <- strsplit(spec, "[,\\s]+", perl = TRUE)[[1]]
    toks <- toks[nchar(toks) > 0]
    ## A-SEC (0.9.4): each token is evaluated in the .mod allowlist sandbox
    ## (numeric literals + elementary arithmetic only), never in the caller
    ## frame.  A disallowed call aborts with dynhr_error_unsafe_mod_expression;
    ## text that is not R, or names an unbound symbol, stays NA as before.
    vapply(toks, function(tok) {
      val <- .dynhr_sandbox_eval(tok, sandbox, .dynhr_safe_fn_names,
                                 context = "the heteroskedastic_shocks value")
      if (is.null(val) || !is.numeric(val) || length(val) != 1L) NA_real_
      else as.numeric(val)
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
      cur_periods <- expand_periods(m_periods[2])
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
