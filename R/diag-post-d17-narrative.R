## R/diag-post-d17-narrative.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D17 narrative identification
## --------------------------------------------------------------------------

## ---- Period parsing ---------------------------------------------------------
## A period label is mapped to the CALENDAR INTERVAL it covers, so that a
## quarterly range selects every observation dated inside it whatever the
## data's own dating convention (quarter-start, quarter-end, monthly, ...).
## Until 0.9.4 a range end "2008-Q4" meant the single day 2008-10-01: monthly
## data lost Nov/Dec, and end-of-quarter dates (2008-12-31) matched nothing.
##
## Accepted labels: "2008-Q4", "2008Q4", "2008q4", "2008:Q4", "2008 Q4";
## "2008-M03", "2008M3", "2008m03", "2008-03"; any ISO date "2008-10-15"
## (an exact day). Returns list(start, end) of Date vectors, NA where the
## label is not recognised.
.d17_period_bounds <- function(x) {
  x <- trimws(as.character(x))
  n <- length(x)
  start <- rep(as.Date(NA), n)
  end   <- rep(as.Date(NA), n)
  first_of <- function(yr, mo) as.Date(sprintf("%04d-%02d-01", yr, mo))
  last_of  <- function(yr, mo) {
    nx_yr <- ifelse(mo == 12L, yr + 1L, yr)
    nx_mo <- ifelse(mo == 12L, 1L, mo + 1L)
    first_of(nx_yr, nx_mo) - 1
  }

  q_re <- "^([0-9]{4})\\s*[-:/ ]?\\s*[Qq]([1-4])$"
  m_re <- "^([0-9]{4})\\s*(?:[-:/ ]?\\s*[Mm]|-)([0-9]{1,2})$"
  is_q <- grepl(q_re, x)
  is_m <- !is_q & grepl(m_re, x, perl = TRUE)
  if (any(is_q)) {
    yr <- as.integer(sub(q_re, "\\1", x[is_q]))
    qq <- as.integer(sub(q_re, "\\2", x[is_q]))
    start[is_q] <- first_of(yr, 3L * qq - 2L)
    end[is_q]   <- last_of(yr, 3L * qq)
  }
  if (any(is_m)) {
    yr <- as.integer(sub(m_re, "\\1", x[is_m], perl = TRUE))
    mo <- as.integer(sub(m_re, "\\2", x[is_m], perl = TRUE))
    ok <- mo >= 1L & mo <= 12L
    idx <- which(is_m)[ok]
    start[idx] <- first_of(yr[ok], mo[ok])
    end[idx]   <- last_of(yr[ok], mo[ok])
  }
  rest <- !is_q & !is_m
  if (any(rest)) {
    d <- as.Date(x[rest], optional = TRUE)
    start[rest] <- d
    end[rest]   <- d
  }
  list(start = start, end = end)
}

## Observation dates -> one Date per row (the start of the labelled period for
## period labels; decimal-year numerics as from time(ts) map to the start of
## the month they fall in).
.d17_obs_dates <- function(dates) {
  if (inherits(dates, "Date")) return(dates)
  if (inherits(dates, "POSIXt")) return(as.Date(dates))
  if (is.numeric(dates)) {
    yr <- floor(dates + 1e-9)
    mo <- floor((dates - yr) * 12 + 1e-6) + 1L
    return(as.Date(sprintf("%04d-%02d-01", as.integer(yr), as.integer(mo))))
  }
  .d17_period_bounds(dates)$start
}

.d17_sign <- function(s, ep_name) {
  if (is.null(s) || length(s) != 1L || is.na(s))
    .dynhr_abort("d17_narrative_identification: episode '", ep_name,
                 "' has no `sign` (expected \"positive\" or \"negative\").")
  s <- tolower(trimws(as.character(s)))
  if (s %in% c("positive", "pos", "+", "up", "increase")) return("positive")
  if (s %in% c("negative", "neg", "-", "down", "decrease")) return("negative")
  .dynhr_abort("d17_narrative_identification: episode '", ep_name,
               "' has sign '", s, "'; expected \"positive\" or \"negative\".")
}

#' D17. Narrative identification (with category and historical decomp support)
#'
#' Tests whether the model's structural shock contributions at known historical
#' episodes have the expected sign and magnitude. Supports both individual
#' shocks and shock categories (from @dynhr:shock_categories), using
#' historical decomposition to compute the combined contribution of a
#' category to a target variable.
#'
#' \strong{Scale.} The quantity tested is the contribution of the shock (or
#' the SUM over a category's shocks) to the target variable, divided by the
#' full-sample standard deviation of that same contribution series. Each
#' shock's own size therefore enters through its smoothed values; the score is
#' invariant to rescaling a shock and is not contaminated by other shocks.
#'
#' \strong{Dating.} Row \code{t} of the historical decomposition is the
#' variable DATED \code{t} and is matched to \code{dates[t]}. A range endpoint
#' given as a period label ("2008-Q4", "2008M10") covers the whole period, so a
#' quarterly range selects every monthly, quarter-start or quarter-end
#' observation inside it.
#'
#' \strong{Tested value:} the peak, i.e. the largest-|value| standardised
#' contribution in the episode window. (Until 0.9.4 the rule was "peak if
#' |peak| >= |window mean|, else the mean", which always selects the peak.)
#' The window mean is reported alongside; review episodes where peak and mean
#' have opposite signs.
#'
#' @param hist_decomp      Output of historical_decomposition() (a list with
#'                         \code{$contributions}, T x n_endo matrices).
#' @param smoothed_shocks  Optional T x n_shock matrix; only used to check that
#'                         its row count matches the decomposition.
#' @param shock_names      unused (shocks are resolved against the
#'                         decomposition's own component names)
#' @param endo_names       character vector of endogenous variable names
#'                         (default: column names of the contributions)
#' @param dates            Date, POSIXct, decimal-year numeric, or period-label
#'                         character vector of length T
#' @param episodes         list of narrative episodes (from .extract_narratives
#'                         or manual): fields name, date_range (length 2),
#'                         shock, variable, sign, optional min_sd, description
#' @param shock_categories named list mapping category -> shock names
#'                         (from meta$shock_categories)
#' @param sd_threshold     Default minimum |contribution/sd| (default 1.0)
#' @return dynhr_diagnostic list
#'
#' @references Antolin-Diaz, J., & Rubio-Ramirez, J. F. (2018). Narrative sign
#'   restrictions for SVARs. \emph{American Economic Review}, 108(10), 2802-2829.
#' @noRd
d17_narrative_identification <- function(hist_decomp,
                                         smoothed_shocks = NULL,
                                         shock_names = NULL,
                                         endo_names = NULL,
                                         dates,
                                         episodes,
                                         shock_categories = list(),
                                         sd_threshold = 1.0,
                                         meta = NULL) {

  ## ---- Validate inputs up front ----
  contribs <- hist_decomp$contributions
  if (!is.list(contribs) || length(contribs) == 0L || is.null(names(contribs)))
    .dynhr_abort("d17_narrative_identification: `hist_decomp$contributions` ",
                 "must be a named list of T x n_endo matrices ",
                 "(historical_decomposition() output).")
  contribs <- lapply(contribs, as.matrix)
  TT <- nrow(contribs[[1L]])
  if (!all(vapply(contribs, nrow, integer(1)) == TT))
    .dynhr_abort("d17_narrative_identification: contribution matrices have ",
                 "differing row counts.")
  if (length(dates) != TT)
    .dynhr_abort("d17_narrative_identification: `dates` has length ",
                 length(dates), " but the decomposition has ", TT,
                 " periods; row t must be dated dates[t].")
  if (!is.null(smoothed_shocks) && NROW(smoothed_shocks) != TT)
    .dynhr_abort("d17_narrative_identification: `smoothed_shocks` has ",
                 NROW(smoothed_shocks), " rows but the decomposition has ",
                 TT, ".")
  col_names <- colnames(contribs[[1L]])
  if (is.null(col_names)) col_names <- endo_names
  if (is.null(col_names) || length(col_names) != ncol(contribs[[1L]]))
    .dynhr_abort("d17_narrative_identification: cannot name the decomposition's ",
                 "columns; supply `endo_names` (length ", ncol(contribs[[1L]]),
                 ").")
  if (!is.numeric(sd_threshold) || length(sd_threshold) != 1L ||
      !is.finite(sd_threshold) || sd_threshold < 0)
    .dynhr_abort("d17_narrative_identification: `sd_threshold` must be a ",
                 "single non-negative number.")
  if (is.null(shock_categories)) shock_categories <- list()

  dates_parsed <- .d17_obs_dates(dates)
  if (anyNA(dates_parsed))
    .dynhr_abort("d17_narrative_identification: could not parse `dates` ",
                 "(first unparseable: '", as.character(dates[which(is.na(dates_parsed))[1L]]),
                 "'). Use Date, decimal years, or labels like \"2008-Q4\" / \"2008-M10\".")

  ## Validate and normalise every episode before evaluating any.
  episodes <- lapply(seq_along(episodes), function(i) {
    ep <- episodes[[i]]
    nm <- ep$name %||% sprintf("episode_%d", i)
    if (is.null(ep$date_range) || length(ep$date_range) != 2L)
      .dynhr_abort("d17_narrative_identification: episode '", nm,
                   "' needs `date_range` = c(start, end).")
    b <- .d17_period_bounds(ep$date_range)
    if (anyNA(b$start))
      .dynhr_abort("d17_narrative_identification: episode '", nm,
                   "' has an unparseable date_range (",
                   paste(ep$date_range, collapse = " : "), ").")
    if (b$start[1L] > b$end[2L])
      .dynhr_abort("d17_narrative_identification: episode '", nm,
                   "' date_range starts after it ends.")
    if (is.null(ep$shock) || !nzchar(ep$shock))
      .dynhr_abort("d17_narrative_identification: episode '", nm,
                   "' has no `shock`.")
    if (is.null(ep$variable) || !nzchar(ep$variable))
      .dynhr_abort("d17_narrative_identification: episode '", nm,
                   "' has no `variable`.")
    min_sd <- ep$min_sd %||% sd_threshold
    if (!is.numeric(min_sd) || length(min_sd) != 1L || !is.finite(min_sd))
      .dynhr_abort("d17_narrative_identification: episode '", nm,
                   "' has a non-numeric `min_sd`.")
    list(name = nm, description = ep$description %||% "",
         date_range = ep$date_range, from = b$start[1L], to = b$end[2L],
         shock = ep$shock, variable = ep$variable,
         sign = .d17_sign(ep$sign, nm), min_sd = min_sd)
  })

  ## ---- Resolve shock reference to decomposition components ----
  resolve_shocks <- function(shock_ref) {
    ## A component of the decomposition itself (a shock, or a group built by
    ## historical_decomposition(shock_groups = ...)).
    if (shock_ref %in% setdiff(names(contribs), "initial"))
      return(list(type = "individual", shocks = shock_ref, missing = character(0)))
    if (shock_ref %in% names(shock_categories)) {
      cat_shocks <- trimws(unlist(strsplit(as.character(shock_categories[[shock_ref]]), ",")))
      cat_shocks <- cat_shocks[nzchar(cat_shocks)]
      return(list(type = "category",
                  shocks = intersect(cat_shocks, names(contribs)),
                  missing = setdiff(cat_shocks, names(contribs))))
    }
    list(type = "unknown", shocks = character(0), missing = shock_ref)
  }

  skip <- function(ep, reason) {
    list(name = ep$name, status = "SKIP", reason = reason,
         description = ep$description, shock_ref = ep$shock,
         target_var = ep$variable, date_range = ep$date_range)
  }

  ## ---- Evaluate each episode ----
  checks <- lapply(episodes, function(ep) {

    rows <- which(dates_parsed >= ep$from & dates_parsed <= ep$to)
    if (length(rows) == 0)
      return(skip(ep, "date range outside sample"))

    shock_info <- resolve_shocks(ep$shock)
    if (length(shock_info$shocks) == 0) {
      why <- if (shock_info$type == "category")
        sprintf("no member of category '%s' in the decomposition (%s)",
                ep$shock, paste(shock_info$missing, collapse = ","))
      else sprintf("shock/category '%s' not found", ep$shock)
      return(skip(ep, why))
    }

    var_idx <- match(ep$variable, col_names)
    if (is.na(var_idx))
      return(skip(ep, sprintf("variable '%s' not found", ep$variable)))

    ## Combined contribution of the shock(s) to the target variable.
    member_series <- lapply(shock_info$shocks,
                            function(s) contribs[[s]][, var_idx])
    names(member_series) <- shock_info$shocks
    combined <- Reduce(`+`, member_series)
    full_sd  <- stats::sd(combined, na.rm = TRUE)
    if (!is.finite(full_sd) || full_sd < 1e-15)
      return(skip(ep, "zero variance in contribution"))

    series_sd   <- combined / full_sd
    window_norm <- series_sd[rows]
    if (!all(is.finite(window_norm)))
      return(skip(ep, "non-finite contribution in window"))

    peak_idx <- which.max(abs(window_norm))
    peak_sd  <- window_norm[peak_idx]
    mean_sd  <- mean(window_norm)

    check_val   <- peak_sd
    actual_sign <- if (check_val > 0) "positive" else if (check_val < 0) "negative" else "zero"
    sign_ok     <- actual_sign == ep$sign
    mag_ok      <- abs(check_val) >= ep$min_sd

    status <- if (sign_ok && mag_ok) "PASS" else if (sign_ok) "WEAK" else "FAIL"

    list(
      name            = ep$name,
      description     = ep$description,
      shock_ref       = ep$shock,
      shock_type      = shock_info$type,
      shock_members   = shock_info$shocks,
      missing_members = shock_info$missing,
      target_var      = ep$variable,
      date_range      = ep$date_range,
      rows            = rows,
      expected_sign   = ep$sign,
      actual_sign     = actual_sign,
      min_sd          = ep$min_sd,
      check_sd        = check_val,
      peak_sd         = peak_sd,
      mean_sd         = mean_sd,
      peak_raw        = combined[rows[peak_idx]],
      peak_date       = as.character(dates[rows[peak_idx]]),
      contrib_sd      = full_sd,
      sign_ok         = sign_ok,
      mag_ok          = mag_ok,
      status          = status,
      series_sd       = series_sd,
      member_contribs = lapply(member_series, function(v) v[rows])
    )
  })

  ## ---- Summary ----
  n_total  <- length(checks)
  statuses <- vapply(checks, `[[`, character(1), "status")
  n_pass   <- sum(statuses == "PASS")
  n_weak   <- sum(statuses == "WEAK")
  n_fail   <- sum(statuses == "FAIL")
  n_skip   <- sum(statuses == "SKIP")

  ## Nothing evaluated -> not a verdict (was FALSE, i.e. a FAIL badge).
  pass <- if (n_total == n_skip) NA else n_fail == 0

  detail_of <- function(ch) {
    if (ch$status == "SKIP") return(ch$reason)
    sprintf("peak %.2fsd at %s (window mean %.2fsd); expected %s, |peak| >= %.2f%s",
            ch$peak_sd, ch$peak_date, ch$mean_sd, ch$expected_sign, ch$min_sd,
            if (length(ch$missing_members))
              sprintf("; members not in decomposition: %s",
                      paste(ch$missing_members, collapse = ","))
            else "")
  }

  lines <- sprintf(
    "D17 Narrative identification: %d episodes, %d PASS, %d WEAK, %d FAIL, %d SKIP.",
    n_total, n_pass, n_weak, n_fail, n_skip)
  for (ch in checks) {
    if (ch$status == "SKIP") {
      lines <- c(lines, sprintf("  [SKIP] %s: %s", ch$name, ch$reason))
    } else {
      type_tag <- if (ch$shock_type == "category") {
        sprintf("%s={%s}", ch$shock_ref, paste(ch$shock_members, collapse = ","))
      } else {
        ch$shock_ref
      }
      lines <- c(lines, sprintf("  [%s] %s: %s -> %s | %s", ch$status, ch$name,
                                type_tag, ch$target_var, detail_of(ch)))
    }
  }

  ## ---- Plot: standardised contribution with episode windows ----
  plots <- list()
  evaluated <- Filter(function(ch) ch$status != "SKIP", checks)
  if (requireNamespace("ggplot2", quietly = TRUE) && length(evaluated) > 0L) {
    panel_of <- function(ch) sprintf("%s -> %s", ch$shock_ref, ch$target_var)
    panels   <- vapply(evaluated, panel_of, character(1))
    keep     <- !duplicated(panels)
    plot_df <- do.call(rbind, lapply(evaluated[keep], function(ch)
      data.frame(date = dates_parsed, value_sd = ch$series_sd,
                 panel = panel_of(ch), stringsAsFactors = FALSE)))

    ## Shade exactly the observations the check used, padded half a period
    ## either side so a one-period episode is visible.
    half <- if (TT > 1L) as.numeric(stats::median(diff(dates_parsed))) / 2 else 15
    ep_df <- do.call(rbind, lapply(evaluated, function(ch) {
      data.frame(panel  = panel_of(ch),
                 xmin   = dates_parsed[min(ch$rows)] - half,
                 xmax   = dates_parsed[max(ch$rows)] + half,
                 thr    = if (ch$expected_sign == "positive") ch$min_sd else -ch$min_sd,
                 check  = ch$check_sd,
                 xmid   = dates_parsed[ch$rows[which.max(abs(ch$series_sd[ch$rows]))]],
                 verdict = ch$status,
                 stringsAsFactors = FALSE)
    }))
    ep_df$verdict <- factor(ep_df$verdict, levels = c("PASS", "WEAK", "FAIL"))
    verdict_cols <- c(PASS = dynhr_colours$green, WEAK = dynhr_colours$orange,
                      FAIL = dynhr_colours$red)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = date, y = value_sd)) +
      ggplot2::geom_rect(
        data = ep_df,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                     fill = verdict),
        inherit.aes = FALSE, alpha = 0.2) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.4) +
      ggplot2::geom_segment(
        data = ep_df,
        ggplot2::aes(x = xmin, xend = xmax, y = thr, yend = thr),
        inherit.aes = FALSE, linetype = "dashed", colour = "grey20") +
      ggplot2::geom_point(
        data = ep_df,
        ggplot2::aes(x = xmid, y = check, colour = verdict),
        inherit.aes = FALSE, size = 2) +
      ggplot2::scale_fill_manual(values = verdict_cols, drop = TRUE,
                                 name = "Episode") +
      ggplot2::scale_colour_manual(values = verdict_cols, drop = TRUE,
                                   guide = "none") +
      ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 2) +
      ggplot2::labs(
        title    = "D17: Narrative episodes in the historical decomposition",
        subtitle = paste("Shock/category contribution to the target variable,",
                         "in SDs of that contribution;\nshading = episode",
                         "window, dashed = required sign x min_sd, dot = tested value"),
        x = NULL, y = "Contribution / its sample SD"
      ) +
      theme_dynhr()

    plots$narrative_decomp <- .apply_meta(p, meta)
  }

  badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
  ep_lines <- vapply(checks, function(r)
    sprintf("    %s: %s (%s)", r$name, r$status, detail_of(r)), character(1))
  action <- if (is.na(pass))
    "No episode could be evaluated; check date ranges, shock/category and variable names."
  else if (n_fail > 0)
    sprintf("%d episode(s) FAILED (wrong sign). %d episode(s) WEAK (insufficient magnitude). Check shock category assignments and sign conventions.",
            n_fail, n_weak)
  else if (n_weak > 0)
    sprintf("0 failed, %d episode(s) WEAK (contribution below min_sd but correct sign). Weak != failed: model broadly consistent but response may be muted.",
            n_weak)
  else
    "All narrative restrictions satisfied."

  .make_result(
    result  = list(checks = checks, n_pass = n_pass, n_weak = n_weak,
                   n_fail = n_fail, n_skip = n_skip),
    pass    = pass,
    plots   = plots,
    summary = paste(lines, collapse = "\n"),
    llm_summary = paste(c(
      sprintf("D17 | Narrative Identification | %s", badge),
      sprintf("  episodes=%d pass=%d weak=%d fail=%d skip=%d",
              n_total, n_pass, n_weak, n_fail, n_skip),
      ep_lines,
      sprintf("  action: %s", action)
    ), collapse = "\n")
  )
}
