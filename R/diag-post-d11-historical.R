## R/diag-post-d11-historical.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D11 historical decomposition diagnostic.
## NOTE: The exported historical_decomposition() is defined in smoother-monolith.R.
## The internal helper here is .compute_hd_from_matrices() to avoid shadowing.
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
    pal <- if (n_shocks <= length(dynhr_palette)) dynhr_palette[seq_len(n_shocks)]
    else colorRampPalette(dynhr_palette)(n_shocks)

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
        scale_fill_dynhr_light(name = "Shock") +
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
      scale_fill_dynhr_light(name = "Shock") +
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


#' Historical decomposition via Kalman smoother
#'
#' Decomposes the path of each endogenous variable into contributions
#' from each structural shock, using smoothed shocks and the state-space
#' transition matrices.
#'
#' For x_{t+1} = F x_t + G eps_t, y_t = H x_t:
#'   y_t = H * sum_{s=1}^{t} F^{t-s} G[:,j] eps_{j,s}  for each shock j
#'
#' @param smoothed_shocks  T x n_shock matrix (from kalman_smoother)
#' @param F_mat            n_endo x n_endo state transition matrix
#' @param G_mat            n_endo x n_shock shock impact matrix
#' @param endo_names       character vector of endogenous variable names
#' @param shock_names      character vector of shock names
#' @return List with:
#'   $contributions  named list of T x n_endo matrices, one per shock
#'   $total          T x n_endo matrix (sum of all contributions; should equal the smoothed states)
#' @noRd
.compute_hd_from_matrices <- function(smoothed_shocks, F_mat, G_mat,
                                      endo_names, shock_names) {

  TT    <- nrow(smoothed_shocks)
  n_end <- nrow(F_mat)
  n_shk <- ncol(G_mat)

  ## Pre-allocate: one T x n_endo matrix per shock
  contributions <- setNames(
    lapply(seq_len(n_shk), function(j) matrix(0, TT, n_end)),
    shock_names
  )

  ## Run forward simulation for each shock independently
  for (j in seq_len(n_shk)) {
    x_j <- rep(0, n_end)   # state attributable to shock j
    g_j <- G_mat[, j]      # impact column for shock j

    for (t in seq_len(TT)) {
      eps_jt <- smoothed_shocks[t, j]
      if (is.na(eps_jt)) eps_jt <- 0
      x_j <- F_mat %*% x_j + g_j * eps_jt
      contributions[[j]][t, ] <- as.numeric(x_j)
    }
    colnames(contributions[[j]]) <- endo_names
  }

  ## Total (should reconstruct smoothed states)
  total <- Reduce(`+`, contributions)
  colnames(total) <- endo_names

  list(contributions = contributions, total = total)
}
