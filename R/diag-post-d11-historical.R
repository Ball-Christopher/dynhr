## R/diag-post-d11-historical.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D11 historical decomposition diagnostic.
## NOTE: The exported historical_decomposition() is defined in smoother-monolith.R.
## The matrix-level helper was removed 2026-09 (unused).
## --------------------------------------------------------------------------

#' D11. Historical decomposition
#'
#' Produces stacked area plots of the historical decomposition, with
#' NZ macroeconomic event markers overlaid for narrative context.
#'
#' @param hd_data      Data frame with columns: date, variable, shock, value.
#'                     Or a list of matrices (variable -> shocks x time).
#'                     Or the native dynhr smoother format: a list with
#'                     \code{$contributions} (named list of T x n_endo matrices,
#'                     one per shock), \code{$endo_names}, and \code{$exo_names}.
#' @param dates        Date vector (required if hd_data is a list of matrices)
#' @param var_names    Character vector (optional)
#' @param shock_names  Character vector (optional)
#' @param show_events  Logical -- annotate NZ macroeconomic events (default TRUE)
#' @return dynhr_diagnostic list
#'
#' @references Banbura, M., Giannone, D., & Reichlin, L. (2010). Large Bayesian vector
#'   auto regressions. \emph{Journal of Applied Econometrics}, 25(1), 71-92.
#' @noRd
d11_historical_decomposition <- function(hd_data,
                                         dates       = NULL,
                                         var_names   = NULL,
                                         shock_names = NULL,
                                         show_events = TRUE,
                                         max_vars    = 9L,
                                         meta        = NULL) {
    # Native dynhr smoother format: a list with $contributions (one [endo x time]
    # matrix per shock), $total, $endo_names, $exo_names. Convert it into the
    # per-variable [shock x time] list that the generic path below consumes.
    if (is.list(hd_data) && !is.null(hd_data$contributions) &&
        !is.null(hd_data$endo_names)) {
      contribs <- hd_data$contributions
      en <- hd_data$endo_names
      sn <- hd_data$exo_names %||% names(contribs) %||%
            paste0("shock_", seq_along(contribs))
      total <- hd_data$total
      # Select the most dynamic variables (largest total variance) when there
      # are more than max_vars, so the output stays readable.
      keep <- seq_along(en)
      if (!is.null(var_names)) {
        keep <- which(en %in% var_names)
      } else if (length(en) > max_vars && !is.null(total)) {
        keep <- order(apply(total, 1, stats::var), decreasing = TRUE)[seq_len(max_vars)]
      }
      hd_data <- stats::setNames(lapply(keep, function(vi) {
        mat <- do.call(rbind, lapply(contribs, function(C) C[vi, ]))
        rownames(mat) <- sn
        mat
      }), en[keep])
      var_names <- en[keep]
    }

    # Convert list-of-matrices to long data frame
    if (is.list(hd_data) && !is.data.frame(hd_data)) {
      # Explicit early return instead of stop() so orchestrator receives a clean result
      if (is.null(dates)) {
        return(.make_result(
          pass    = NA,
          plots   = list(),
          summary = "D11 skipped: dates required for list-of-matrices hd input."
        ))
      }

      dfs <- list()
      for (v in names(hd_data)) {
        mat <- hd_data[[v]]

        # Guard: ensure matrix column count matches dates length
        n_dates  <- length(dates)
        n_tcols  <- if (is.matrix(mat) || is.data.frame(mat)) ncol(mat) else length(mat)
        if (is.na(n_dates) || is.na(n_tcols) || n_dates < 1L || n_tcols < 1L) {
          return(.make_result(
            pass    = NA,
            plots   = list(),
            summary = sprintf(
              "D11 skipped: variable '%s' has %s columns but dates has length %s. Ensure hd matrices and dates are consistent.",
              v, as.character(n_tcols), as.character(n_dates))
          ))
        }
        if (n_tcols != n_dates) {
          return(.make_result(
            pass    = NA,
            plots   = list(),
            summary = sprintf(
              "D11 skipped: variable '%s' matrix has %d columns but dates has length %d. Cannot align time periods.",
              v, n_tcols, n_dates)
          ))
        }

        s_names <- rownames(mat) %||% paste0("shock_", seq_len(nrow(mat)))
        for (s in seq_len(nrow(mat))) {
          dfs[[length(dfs) + 1]] <- data.frame(
            date     = dates,
            variable = v,
            shock    = s_names[s],
            value    = mat[s, ]
          )
        }
      }
      hd_data <- do.call(rbind, dfs)
    }

    if (!inherits(hd_data$date, "Date")) {
      hd_data$date <- as.Date(hd_data$date)
    }

    if (is.null(var_names))   var_names   <- unique(hd_data$variable)
    if (is.null(shock_names)) shock_names <- unique(hd_data$shock)

    n_shocks <- length(shock_names)

    # ---- Fill palette --------------------------------------------------
    # "initial" (the shock-free trajectory from the smoothed initial state)
    # and "constraint" (accumulated OBC binding intercepts) are components of
    # the decomposition, NOT structural shocks: they take neutral greys, and
    # the shock palette is sized by the shock rows ALONE.  Sizing it by every
    # row would recolour every shock the moment an "initial" row appears --
    # the same decomposition would come out in different colours depending on
    # whether the initial condition was supplied.
    nonshock_fill <- c(initial = dynhr_na_fill, constraint = dynhr_na_colour)
    shock_only    <- setdiff(shock_names, names(nonshock_fill))
    fill_values   <- c(
      stats::setNames(.dynhr_palette_fun(dynhr_palette_light)(length(shock_only)),
                      shock_only),
      nonshock_fill[intersect(names(nonshock_fill), shock_names)]
    )

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # Build per-variable plots stored under hd_<v> keys.
    # Also accumulate all data for a single consolidated 3x3 faceted overview.
    for (v in var_names) {
      hd_v <- hd_data[hd_data$variable == v, ]

      # Use geom_col (stat="identity") so positive and negative contributions
      # stack correctly on their respective sides of zero.
      p_base <- ggplot2::ggplot(hd_v,
                                ggplot2::aes(x = date, y = value, fill = shock)) +
        # width in days on a Date axis; ~85 keeps a small gap at quarterly spacing
        ggplot2::geom_col(position = "stack", alpha = 0.85, width = 85) +
        ggplot2::geom_hline(yintercept = 0,
                            colour = dynhr_colours$grey, linewidth = 0.4) +
        ggplot2::scale_fill_manual(name = "Shock", values = fill_values,
                                   na.value = dynhr_na_fill) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(title = sprintf("D11: Historical decomposition -- %s", v),
             x = NULL, y = "Contribution (% dev from SS)")

      # Add NZ event markers if requested
      if (show_events) {
        events <- .nz_events()
        date_range <- range(hd_v$date, na.rm = TRUE)
        events <- events[events$date >= date_range[1] & events$date <= date_range[2], ]
        if (nrow(events) > 0) {
          p_events <- .add_nz_event_markers(p_base, events)
        } else {
          p_events <- p_base
        }
      } else {
        p_events <- p_base
      }

      # Store only under the hd_<v> key (with event markers / meta)
      plots[[paste0("hd_", v)]] <- .apply_meta(p_events, meta)
    }

    # Consolidated overview: one 3x3 faceted plot for quick cross-variable
    # comparison (replaces the need to flip between 9 separate PNGs).
    # Event markers are not added here since they would clutter the small panels.
    n_cols_overview <- min(3L, length(var_names))
    p_overview <- ggplot2::ggplot(
      hd_data[hd_data$variable %in% var_names, ],
      ggplot2::aes(x = date, y = value, fill = shock)
    ) +
      ggplot2::geom_col(position = "stack", alpha = 0.80, width = 85) +
      ggplot2::geom_hline(yintercept = 0,
                          colour = dynhr_colours$grey, linewidth = 0.3) +
      ggplot2::scale_fill_manual(name = "Shock", values = fill_values,
                                 na.value = dynhr_na_fill) +
      ggplot2::facet_wrap(~ variable, scales = "free_y",
                          ncol = n_cols_overview) +
      # Quarterly x-axis: "2010 Q1" format on a Date axis
      ggplot2::scale_x_date(
        date_breaks = "5 years",
        date_labels = "%Y"
      ) +
      theme_dynhr_compact() +
      ggplot2::theme(
        axis.text.x   = ggplot2::element_text(size = ggplot2::rel(0.75)),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(0.4, "lines")
      ) +
      ggplot2::labs(
        title    = "D11: Historical decomposition overview",
        subtitle = sprintf("%d variables, %d shocks", length(var_names), length(shock_names)),
        x = NULL, y = "Contribution (% dev from SS)")
    attr(p_overview, "dynhr_fig_height") <- 9.0
    plots$hd_overview <- .apply_meta(p_overview, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(hd_data = hd_data, var_names = var_names, shock_names = shock_names),
      pass    = NA,  # Informational
      plots   = plots,
      summary = sprintf(
        "D11 Historical decomposition: %d variables, %d shocks, %d periods (%s to %s).%s",
        length(var_names), n_shocks, length(unique(hd_data$date)),
        {d <- min(hd_data$date); sprintf("%d-Q%d", as.integer(format(d, "%Y")),
                                          (as.integer(format(d, "%m")) - 1L) %/% 3L + 1L)},
        {d <- max(hd_data$date); sprintf("%d-Q%d", as.integer(format(d, "%Y")),
                                          (as.integer(format(d, "%m")) - 1L) %/% 3L + 1L)},
        ifelse(show_events, " NZ event markers included.", "")
      ),
      llm_summary = {
        n_per   <- length(unique(hd_data$date))
        dr      <- range(hd_data$date)
        paste(c(
          "D11 | Historical Decomposition | INFO",
          sprintf("  variables=%d shocks=%d periods=%d",
                  length(var_names), length(shock_names), n_per),
          sprintf("  date_range=%s to %s",
                  format(dr[1], "%Y-%m"), format(dr[2], "%Y-%m")),
          "  action: Informational only. Inspect plots for dominant shock contributors by episode."
        ), collapse = "\n")
      }
    )
}
