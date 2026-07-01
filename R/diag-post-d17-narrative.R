## R/diag-post-d17-narrative.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D17 narrative identification
## --------------------------------------------------------------------------

#' D17. Narrative identification (with category and historical decomp support)
#'
#' Tests whether the model's structural shock contributions at known historical
#' episodes have the expected sign and magnitude. Supports both individual
#' shocks and shock categories (from @dynhr:shock_categories), using
#' historical decomposition to compute the combined contribution of a
#' category to a target variable.
#'
#' \strong{Sign selection rule:} sign classification uses the peak contribution
#' if \code{|peak| >= |mean contribution over the episode window|}, and the
#' mean contribution otherwise.  Episodes where the peak and mean have
#' opposite signs may therefore be classified differently depending on which
#' dominates in absolute value.  For spike episodes (single-period peak), the
#' peak dominates; for sustained episodes, the mean dominates.  Users should
#' review episodes with opposite-sign peak and mean to confirm that the
#' classification matches economic intuition.
#'
#' @param hist_decomp      Output of historical_decomposition()
#' @param smoothed_shocks  T x n_shock matrix (for individual shock checks)
#' @param shock_names      character vector of shock names
#' @param endo_names       character vector of endogenous variable names
#' @param dates            Date or character vector of length T
#' @param episodes         list of narrative episodes (from .extract_narratives or manual)
#' @param shock_categories named list mapping category -> shock names
#'                         (from meta$shock_categories)
#' @param sd_threshold     Default minimum |contribution/sd| (default 1.0)
#' @return dynhr_diagnostic list
#'
#' @references Antolin-Diaz, J., & Rubio-Ramirez, J. F. (2018). Narrative sign
#'   restrictions for SVARs. \emph{American Economic Review}, 108(10), 2802-2829.
#' @noRd
d17_narrative_identification <- function(hist_decomp,
                                         smoothed_shocks,
                                         shock_names,
                                         endo_names,
                                         dates,
                                         episodes,
                                         shock_categories = list(),
                                         sd_threshold = 1.0,
                                         meta = NULL) {

  TT    <- nrow(smoothed_shocks)
  n_shk <- ncol(smoothed_shocks)

  if (is.null(colnames(smoothed_shocks))) {
    colnames(smoothed_shocks) <- shock_names
  }

  ## ---- Parse dates ----
  parse_qdate <- function(qstr) {
    parts <- strsplit(qstr, "-Q|q")[[1]]
    yr <- as.integer(parts[1])
    qq <- as.integer(parts[2])
    as.Date(sprintf("%d-%02d-01", yr, (qq - 1) * 3 + 1))
  }

  if (is.character(dates) && length(dates) > 0 && grepl("Q|q", dates[1])) {
    dates_parsed <- as.Date(sapply(dates, parse_qdate))
  } else {
    dates_parsed <- as.Date(dates)
  }

  ## ---- Helper: get rows for a date range ----
  get_rows <- function(date_range) {
    if (is.null(date_range) || length(date_range) < 2) return(integer(0))
    if (grepl("Q|q", date_range[1])) {
      d_start <- parse_qdate(date_range[1])
      d_end   <- parse_qdate(date_range[2])
    } else {
      d_start <- as.Date(date_range[1])
      d_end   <- as.Date(date_range[2])
    }
    which(dates_parsed >= d_start & dates_parsed <= d_end)
  }

  ## ---- Resolve shock reference to individual shock names ----
  resolve_shocks <- function(shock_ref) {
    # Check if it's a category name
    if (shock_ref %in% names(shock_categories)) {
      cat_shocks <- shock_categories[[shock_ref]]
      # Clean up: trim whitespace, split if comma-separated string
      if (length(cat_shocks) == 1 && grepl(",", cat_shocks)) {
        cat_shocks <- trimws(strsplit(cat_shocks, ",")[[1]])
      } else {
        cat_shocks <- trimws(cat_shocks)
      }
      return(list(type = "category", name = shock_ref, shocks = cat_shocks))
    }
    # Otherwise treat as individual shock
    if (shock_ref %in% shock_names) {
      return(list(type = "individual", name = shock_ref, shocks = shock_ref))
    }
    return(list(type = "unknown", name = shock_ref, shocks = character(0)))
  }

  ## ---- Evaluate each episode ----
  checks <- lapply(episodes, function(ep) {

    rows <- get_rows(ep$date_range)
    if (length(rows) == 0) {
      return(list(name = ep$name, status = "SKIP",
                  reason = "date range outside sample",
                  description = ep$description %||% ""))
    }

    shock_info <- resolve_shocks(ep$shock)
    if (shock_info$type == "unknown" || length(shock_info$shocks) == 0) {
      return(list(name = ep$name, status = "SKIP",
                  reason = sprintf("shock/category '%s' not found", ep$shock),
                  description = ep$description %||% ""))
    }

    target_var <- ep$variable %||% "y"
    var_idx    <- match(target_var, endo_names)
    if (is.na(var_idx)) {
      return(list(name = ep$name, status = "SKIP",
                  reason = sprintf("variable '%s' not found", target_var),
                  description = ep$description %||% ""))
    }

    ## Compute combined contribution of shock(s) to target variable
    ## using historical decomposition
    combined_contrib <- rep(0, TT)
    member_contribs  <- list()

    for (shk in shock_info$shocks) {
      if (shk %in% names(hist_decomp$contributions)) {
        contrib_mat <- hist_decomp$contributions[[shk]]
        shk_contrib <- contrib_mat[, var_idx]
        combined_contrib <- combined_contrib + shk_contrib
        member_contribs[[shk]] <- shk_contrib[rows]
      }
    }

    ## Extract window values and normalise
    window_vals <- combined_contrib[rows]
    full_sd     <- sd(combined_contrib, na.rm = TRUE)

    if (full_sd < 1e-15) {
      return(list(name = ep$name, status = "SKIP",
                  reason = "zero variance in contribution",
                  description = ep$description %||% ""))
    }

    window_norm <- window_vals / full_sd

    ## Find peak (largest absolute value in window)
    peak_idx  <- which.max(abs(window_norm))
    peak_sd   <- window_norm[peak_idx]
    peak_raw  <- window_vals[peak_idx]
    peak_date <- as.character(dates[rows[peak_idx]])

    ## Also compute mean contribution in window (for sustained episodes)
    mean_sd <- mean(window_norm, na.rm = TRUE)

    ## Check sign (use peak for spike episodes, mean for sustained)
    expected_sign <- tolower(ep$sign %||% "negative")
    # Use the value with larger absolute magnitude
    check_val   <- if (abs(peak_sd) >= abs(mean_sd)) peak_sd else mean_sd
    actual_sign <- if (check_val > 0) "positive" else "negative"
    sign_ok     <- actual_sign == expected_sign

    ## Check magnitude
    min_sd <- ep$min_sd %||% sd_threshold
    mag_ok <- abs(check_val) >= min_sd

    ## Overall status
    status <- if (sign_ok && mag_ok) "PASS"
    else if (sign_ok && !mag_ok) "WEAK"
    else "FAIL"

    list(
      name            = ep$name,
      description     = ep$description %||% "",
      shock_ref       = ep$shock,
      shock_type      = shock_info$type,
      shock_members   = shock_info$shocks,
      target_var      = target_var,
      date_range      = ep$date_range,
      rows            = rows,
      expected_sign   = expected_sign,
      actual_sign     = actual_sign,
      peak_sd         = round(peak_sd, 2),
      mean_sd         = round(mean_sd, 2),
      peak_raw        = peak_raw,
      peak_date       = peak_date,
      sign_ok         = sign_ok,
      mag_ok          = mag_ok,
      status          = status,
      member_contribs = member_contribs
    )
  })

  ## ---- Summary ----
  n_total  <- length(checks)
  statuses <- sapply(checks, `[[`, "status")
  n_pass   <- sum(statuses == "PASS")
  n_weak   <- sum(statuses == "WEAK")
  n_fail   <- sum(statuses == "FAIL")
  n_skip   <- sum(statuses == "SKIP")

  pass <- n_fail == 0 && n_skip < n_total

  lines <- c(
    sprintf("D17 Narrative identification: %d episodes, %d PASS, %d WEAK, %d FAIL, %d SKIP.",
            n_total, n_pass, n_weak, n_fail, n_skip)
  )
  for (ch in checks) {
    if (ch$status == "SKIP") {
      lines <- c(lines, sprintf("  [SKIP] %s: %s", ch$name, ch$reason))
    } else {
      type_tag <- if (ch$shock_type == "category") {
        sprintf("%s={%s}", ch$shock_ref, paste(ch$shock_members, collapse = ","))
      } else {
        ch$shock_ref
      }
      lines <- c(lines, sprintf(
        "  [%s] %s: %s -> %s | peak=%.2fsd mean=%.2fsd at %s (expected %s)",
        ch$status, ch$name, type_tag, ch$target_var,
        ch$peak_sd, ch$mean_sd, ch$peak_date, ch$expected_sign
      ))
    }
  }

  ## ---- Plot: contribution time series with episode markers ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE) && length(dates_parsed) == TT) {

    # Collect all unique (shock_ref, target_var) pairs from episodes
    plot_panels <- unique(data.frame(
      shock_ref = sapply(checks, function(ch) ch$shock_ref %||% NA),
      target    = sapply(checks, function(ch) ch$target_var %||% NA),
      stringsAsFactors = FALSE
    ))
    plot_panels <- plot_panels[complete.cases(plot_panels), ]

    if (nrow(plot_panels) > 0) {
      plot_df_list <- list()

      for (i in seq_len(nrow(plot_panels))) {
        shock_ref <- plot_panels$shock_ref[i]
        target    <- plot_panels$target[i]
        var_idx   <- match(target, endo_names)
        if (is.na(var_idx)) next

        shock_info <- resolve_shocks(shock_ref)
        if (length(shock_info$shocks) == 0) next

        combined <- rep(0, TT)
        for (shk in shock_info$shocks) {
          if (shk %in% names(hist_decomp$contributions)) {
            combined <- combined + hist_decomp$contributions[[shk]][, var_idx]
          }
        }

        full_sd <- sd(combined, na.rm = TRUE)
        if (full_sd < 1e-15) next

        plot_df_list[[length(plot_df_list) + 1]] <- data.frame(
          date     = dates_parsed,
          value_sd = combined / full_sd,
          panel    = sprintf("%s -> %s", shock_ref, target),
          stringsAsFactors = FALSE
        )
      }

      if (length(plot_df_list) > 0) {
        plot_df <- do.call(rbind, plot_df_list)

        # Build episode shading data frame with panel column for facet-scoped rects
        ep_rect_rows <- list()
        for (ch in checks) {
          if (ch$status == "SKIP" || is.null(ch$date_range)) next
          panel_name <- sprintf("%s -> %s", ch$shock_ref, ch$target_var)
          if (grepl("Q|q", ch$date_range[1])) {
            d_start <- parse_qdate(ch$date_range[1])
            d_end   <- parse_qdate(ch$date_range[2])
          } else {
            d_start <- as.Date(ch$date_range[1])
            d_end   <- as.Date(ch$date_range[2])
          }
          rect_col <- switch(ch$status,
                             PASS = dynhr_colours$green,
                             WEAK = dynhr_colours$orange,
                             FAIL = dynhr_colours$red,
                             dynhr_colours$grey)
          ep_rect_rows[[length(ep_rect_rows) + 1L]] <- data.frame(
            panel  = panel_name,
            xmin   = d_start,
            xmax   = d_end,
            fill   = rect_col,
            stringsAsFactors = FALSE
          )
        }

        p <- ggplot2::ggplot(plot_df,
                             ggplot2::aes(x = date, y = value_sd)) +
          ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.4) +
          ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
          ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 2) +
          ggplot2::labs(
            title    = "D17: Historical Decomposition -- Narrative Episodes",
            subtitle = "Category/shock contributions to target variable (SD units)",
            x = NULL, y = "Contribution (SD units)"
          ) +
          theme_dynhr()

        # Add episode shading (facet-scoped: rects only render in matching panel)
        if (length(ep_rect_rows) > 0L) {
          ep_df <- do.call(rbind, ep_rect_rows)
          for (fc in unique(ep_df$fill)) {
            sub_df <- ep_df[ep_df$fill == fc, , drop = FALSE]
            p <- p + ggplot2::geom_rect(
              data = sub_df,
              ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE,
              fill = fc, alpha = 0.15
            )
          }
        }

        plots$narrative_decomp <- .apply_meta(p, meta)
      }
    }
  }

  .make_result(
    result  = list(checks = checks, n_pass = n_pass, n_weak = n_weak,
                   n_fail = n_fail, n_skip = n_skip),
    pass    = pass,
    plots   = plots,
    summary = paste(lines, collapse = "\n"),
    llm_summary = {
      episode_results <- checks
      badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
      ep_lines <- vapply(episode_results, function(r) {
        detail <- r$detail %||% r$reason %||% NULL
        sprintf("    %s: %s%s", r$name, r$status,
                if (!is.null(detail)) sprintf(" (%s)", detail) else "")
      }, character(1))
      n_ep    <- length(episode_results)
      n_pass_ <- sum(vapply(episode_results, function(r) r$status == "PASS", logical(1)))
      n_weak_ <- sum(vapply(episode_results, function(r) r$status == "WEAK", logical(1)))
      n_fail_ <- sum(vapply(episode_results, function(r) r$status == "FAIL", logical(1)))
      paste(c(
        sprintf("D17 | Narrative Identification | %s", badge),
        sprintf("  episodes=%d pass=%d weak=%d fail=%d", n_ep, n_pass_, n_weak_, n_fail_),
        paste(ep_lines, collapse = "\n"),
        sprintf("  action: %s",
                if (n_fail_ > 0)
                  sprintf("%d episode(s) FAILED (wrong sign/magnitude). %d episode(s) WEAK (insufficient magnitude). Check shock category assignments and sign conventions.",
                          n_fail_, n_weak_)
                else if (n_weak_ > 0)
                  sprintf("0 failed, %d episode(s) WEAK (contribution below sd_threshold but correct sign). Weak != failed: model broadly consistent but response may be muted.",
                          n_weak_)
                else
                  "All narrative restrictions satisfied.")
      ), collapse = "\n")
    }
  )
}
