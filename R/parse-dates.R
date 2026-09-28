## R/parse-dates.R
## --------------------------------------------------------------------------
## Minimal, dependency-free date-literal support for the period lists of the
## filter_tunes (Dynare #2020), heteroskedastic_shocks (Dynare 7 dates
## syntax) and shock_paths blocks: parse "1990q1" / "1990m3" / "1990" style
## literals into a (year, subperiod, freq) triple, and convert a date literal
## to an integer sample period given a sample-start date.
##
## Supported literals (case-insensitive, Dynare's `dates` syntax):
##   YYYYqQ   quarterly,  Q in 1..4
##   YYYYmM   monthly,    M in 1..12
##   YYYYsS / YYYYhS   semi-annual (bi-annual), S in 1..2  (Dynare 7)
##   YYYYy / YYYYa     annual
##   YYYY     annual (only as a sample start: in a period list a bare integer
##            is an integer period, as in Dynare)
## --------------------------------------------------------------------------

#' Parse a Dynare-style date literal
#'
#' @param str Character scalar, e.g. "1990q1", "1990m3", "1990".
#' @return A list with `year` (integer), `period` (integer subperiod, 1-based;
#'   always 1 for annual), and `freq` (integer: 4 = quarterly, 12 = monthly,
#'   2 = semi-annual, 1 = annual), or `NULL` if `str` is not a recognised
#'   date literal.
#' @noRd
.parse_dynare_date <- function(str) {
  str <- trimws(str)
  m <- regmatches(str, regexec(
    "^([0-9]{4})([qQmMsShH])([0-9]{1,2})$", str, perl = TRUE))[[1]]
  if (length(m) > 0 && nchar(m[1]) > 0) {
    year   <- as.integer(m[2])
    letter <- tolower(m[3])
    sub    <- as.integer(m[4])
    freq   <- switch(letter, q = 4L, m = 12L, s = 2L, h = 2L)
    if (sub < 1L || sub > freq) return(NULL)
    return(list(year = year, period = sub, freq = freq))
  }
  m2 <- regmatches(str, regexec("^([0-9]{4})[yYaA]?$", str, perl = TRUE))[[1]]
  if (length(m2) > 0 && nchar(m2[1]) > 0) {
    return(list(year = as.integer(m2[2]), period = 1L, freq = 1L))
  }
  NULL
}

#' Convert a date literal to an integer sample period
#'
#' Sample period 1 corresponds to `sample_start`. Both `date_str` and
#' `sample_start` must use the same frequency (year/half/quarter/month).
#'
#' @param date_str     Date literal, e.g. "1991q2".
#' @param sample_start Date literal for period 1 of the sample, e.g. "1990q1".
#' @param context      Block name used in error messages.
#' @return Integer sample period (1-based), or `NULL` if either literal fails
#'   to parse.  Mismatched frequencies abort (`dynhr_error_mod_syntax`).
#' @noRd
.date_to_period <- function(date_str, sample_start, context = "filter_tunes") {
  d0 <- .parse_dynare_date(sample_start)
  d1 <- .parse_dynare_date(date_str)
  if (is.null(d0) || is.null(d1)) return(NULL)
  if (d0$freq != d1$freq)
    .dynhr_abort(sprintf(
      "%s: date '%s' frequency does not match sample_start '%s'.",
      context, date_str, sample_start), class = "dynhr_error_mod_syntax")
  freq <- d0$freq
  offset <- (d1$year - d0$year) * freq + (d1$period - d0$period)
  as.integer(offset + 1L)
}

#' Inverse of `.date_to_period()`: sample period(s) -> Dynare date literal(s)
#'
#' @param period       Integer vector of sample periods (1 = `sample_start`).
#' @param sample_start Date literal for sample period 1.
#' @return Character vector of date literals in Dynare's upper-case form
#'   (`1990Q1`, `1990M3`, `1990S2`, `1990Y`).
#' @noRd
.period_to_date <- function(period, sample_start) {
  d0 <- .parse_dynare_date(sample_start)
  if (is.null(d0))
    .dynhr_abort("could not parse date literal '", sample_start, "'.",
                 class = "dynhr_error_mod_syntax")
  f   <- d0$freq
  idx <- (d0$year * f + (d0$period - 1L)) + (as.integer(period) - 1L)
  yr  <- idx %/% f
  sub <- idx %% f + 1L
  if (f == 1L) return(paste0(yr, "Y"))
  letter <- switch(as.character(f), "2" = "S", "4" = "Q", "12" = "M")
  paste0(yr, letter, sub)
}

#' Expand a "lo:hi" or single-value period token into an integer sequence
#'
#' Accepts plain integers ("3", "3:7") and Dynare date literals
#' ("1990q1", "1990q1:1991q4"). Date literals require `sample_start`.
#'
#' @param token        Character scalar, a single period token (no commas).
#' @param sample_start Optional date literal for sample period 1. Required
#'   only if `token` contains date literals.
#' @param context      Block name used in error messages.
#' @param no_start_why Optional sentence appended to the error raised when a
#'   date literal meets a `NULL` `sample_start` (says where the start date
#'   would come from).
#' @return Integer vector of sample periods.
#' @noRd
.expand_period_token <- function(token, sample_start = NULL,
                                 context = "filter_tunes",
                                 no_start_why = NULL) {
  token <- trimws(token)
  parts <- strsplit(token, ":", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  if (length(parts) == 1L) {
    p <- parts[1]
  } else if (length(parts) == 2L) {
    lo_str <- parts[1]; hi_str <- parts[2]
  } else {
    .dynhr_abort(sprintf("%s: malformed period range '%s'.", context, token),
                 class = "dynhr_error_mod_syntax")
  }

  to_int <- function(s) {
    if (grepl("^-?[0-9]+$", s)) return(as.integer(s))
    if (is.null(.parse_dynare_date(s)))
      .dynhr_abort(sprintf("%s: could not parse date literal '%s'.",
                           context, s), class = "dynhr_error_mod_syntax")
    if (is.null(sample_start))
      .dynhr_abort(sprintf(paste0(
        "%s: date literal '%s' used in periods without an ",
        "estimation sample start date; use integer periods instead.%s"),
        context, s, if (is.null(no_start_why)) "" else paste0(" ", no_start_why)),
        class = c("dynhr_error_mod_date_unresolved", "dynhr_error_mod_syntax"))
    .date_to_period(s, sample_start, context = context)
  }

  if (length(parts) == 1L) {
    return(to_int(p))
  }
  lo <- to_int(lo_str); hi <- to_int(hi_str)
  if (hi < lo)
    .dynhr_abort(sprintf("%s: empty period range '%s' (hi < lo).", context,
                         token), class = "dynhr_error_mod_syntax")
  seq.int(lo, hi)
}

#' Expand a comma- or space-separated list of period tokens
#'
#' Tolerates trailing commas (Dynare #2030).
#'
#' @param spec         Character scalar, e.g. "3, 5:7, 10" or "1990q1:1991q4".
#' @param sample_start Optional date literal for sample period 1.
#' @param context      Block name used in error messages.
#' @param no_start_why See `.expand_period_token()`.
#' @return Integer vector of sample periods (concatenation of all tokens, in
#'   the order given; not deduplicated or sorted).
#' @noRd
.expand_period_list <- function(spec, sample_start = NULL,
                                context = "filter_tunes",
                                no_start_why = NULL) {
  spec <- trimws(spec)
  spec <- sub(",\\s*$", "", spec)          # trailing comma tolerance (#2030)
  if (!nchar(spec)) return(integer(0))
  ## Split on commas and/or whitespace.
  toks <- strsplit(spec, "[,\\s]+", perl = TRUE)[[1]]
  toks <- toks[nchar(toks) > 0]
  out <- integer(0)
  for (tok in toks)
    out <- c(out, .expand_period_token(tok, sample_start, context = context,
                                       no_start_why = no_start_why))
  out
}
