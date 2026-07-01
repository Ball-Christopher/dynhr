## R/parse-dates.R
## --------------------------------------------------------------------------
## Minimal, dependency-free date-literal support for the filter_tunes block
## (Dynare #2020). Only what filter_tunes needs: parse "1990q1" / "1990m3" /
## "1990" style literals into a (year, subperiod, freq) triple, and convert a
## date literal to an integer sample period given a sample-start date.
##
## Supported literals:
##   YYYYqQ   quarterly,  Q in 1..4
##   YYYYmM   monthly,    M in 1..12
##   YYYY     annual
## --------------------------------------------------------------------------

#' Parse a Dynare-style date literal
#'
#' @param str Character scalar, e.g. "1990q1", "1990m3", "1990".
#' @return A list with `year` (integer), `period` (integer subperiod, 1-based;
#'   always 1 for annual), and `freq` (integer: 4 = quarterly, 12 = monthly,
#'   1 = annual), or `NULL` if `str` is not a recognised date literal.
#' @noRd
.parse_dynare_date <- function(str) {
  str <- trimws(str)
  m <- regmatches(str, regexec(
    "^([0-9]{4})([qQmM])([0-9]{1,2})$", str, perl = TRUE))[[1]]
  if (length(m) > 0 && nchar(m[1]) > 0) {
    year   <- as.integer(m[2])
    letter <- tolower(m[3])
    sub    <- as.integer(m[4])
    freq   <- if (letter == "q") 4L else 12L
    if (sub < 1L || sub > freq) return(NULL)
    return(list(year = year, period = sub, freq = freq))
  }
  m2 <- regmatches(str, regexec("^([0-9]{4})$", str, perl = TRUE))[[1]]
  if (length(m2) > 0 && nchar(m2[1]) > 0) {
    return(list(year = as.integer(m2[2]), period = 1L, freq = 1L))
  }
  NULL
}

#' Convert a date literal to an integer sample period
#'
#' Sample period 1 corresponds to `sample_start`. Both `date_str` and
#' `sample_start` must use the same frequency (year/quarter/month).
#'
#' @param date_str     Date literal, e.g. "1991q2".
#' @param sample_start Date literal for period 1 of the sample, e.g. "1990q1".
#' @return Integer sample period (1-based), or `NULL` if either literal fails
#'   to parse or the frequencies disagree.
#' @noRd
.date_to_period <- function(date_str, sample_start) {
  d0 <- .parse_dynare_date(sample_start)
  d1 <- .parse_dynare_date(date_str)
  if (is.null(d0) || is.null(d1)) return(NULL)
  if (d0$freq != d1$freq)
    stop(sprintf(
      "filter_tunes: date '%s' frequency does not match sample_start '%s'.",
      date_str, sample_start), call. = FALSE)
  freq <- d0$freq
  offset <- (d1$year - d0$year) * freq + (d1$period - d0$period)
  as.integer(offset + 1L)
}

#' Expand a "lo:hi" or single-value period token into an integer sequence
#'
#' Accepts plain integers ("3", "3:7") and Dynare date literals
#' ("1990q1", "1990q1:1991q4"). Date literals require `sample_start`.
#'
#' @param token        Character scalar, a single period token (no commas).
#' @param sample_start Optional date literal for sample period 1. Required
#'   only if `token` contains date literals.
#' @return Integer vector of sample periods.
#' @noRd
.expand_period_token <- function(token, sample_start = NULL) {
  token <- trimws(token)
  parts <- strsplit(token, ":", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  if (length(parts) == 1L) {
    p <- parts[1]
  } else if (length(parts) == 2L) {
    lo_str <- parts[1]; hi_str <- parts[2]
  } else {
    stop(sprintf("filter_tunes: malformed period range '%s'.", token),
         call. = FALSE)
  }

  to_int <- function(s) {
    if (grepl("^-?[0-9]+$", s)) return(as.integer(s))
    if (is.null(sample_start))
      stop(sprintf(paste0(
        "filter_tunes: date literal '%s' used in periods without an ",
        "estimation sample start date; use integer periods instead."),
        s), call. = FALSE)
    p <- .date_to_period(s, sample_start)
    if (is.null(p))
      stop(sprintf("filter_tunes: could not parse date literal '%s'.", s),
           call. = FALSE)
    p
  }

  if (length(parts) == 1L) {
    return(to_int(p))
  }
  lo <- to_int(lo_str); hi <- to_int(hi_str)
  if (hi < lo)
    stop(sprintf("filter_tunes: empty period range '%s' (hi < lo).", token),
         call. = FALSE)
  seq.int(lo, hi)
}

#' Expand a comma- or space-separated list of period tokens
#'
#' Tolerates trailing commas (Dynare #2030).
#'
#' @param spec         Character scalar, e.g. "3, 5:7, 10" or "1990q1:1991q4".
#' @param sample_start Optional date literal for sample period 1.
#' @return Integer vector of sample periods (concatenation of all tokens, in
#'   the order given; not deduplicated or sorted).
#' @noRd
.expand_period_list <- function(spec, sample_start = NULL) {
  spec <- trimws(spec)
  spec <- sub(",\\s*$", "", spec)          # trailing comma tolerance (#2030)
  if (!nchar(spec)) return(integer(0))
  ## Split on commas and/or whitespace.
  toks <- strsplit(spec, "[,\\s]+", perl = TRUE)[[1]]
  toks <- toks[nchar(toks) > 0]
  out <- integer(0)
  for (tok in toks)
    out <- c(out, .expand_period_token(tok, sample_start))
  out
}
