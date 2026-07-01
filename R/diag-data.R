## R/diag-data.R
## --------------------------------------------------------------------------
## Phase-3+ addition.
##
## diag_data() -- plots observable time series and summary statistics for the
## estimation dataset.  Useful as a first-pass sanity check before estimation.
## --------------------------------------------------------------------------

#' Data diagnostic: plot observable time series and summary statistics
#'
#' Produces time-series plots of every observable variable, a correlation
#' heatmap, and a summary statistics table.  When \code{dates} are provided
#' the x-axis is in calendar time; otherwise in observation index.
#'
#' @param data       Numeric matrix or data frame (T x n_obs).  Columns are
#'   observables; rows are time periods.
#' @param obs_names  Character vector of observable names.  Defaults to
#'   \code{colnames(data)}.
#' @param dates      Optional Date (or coercible) vector of length T for the
#'   x-axis.  NULL uses integer index.
#' @param meta       A \code{\link{diag_meta}} object for plot provenance.
#' @return dynhr_diagnostic list with plots: \code{time_series},
#'   \code{distribution}, \code{correlation}.
#' @noRd
diag_data <- function(data,
                      obs_names = NULL,
                      dates     = NULL,
                      meta      = NULL) {
    if (!is.matrix(data)) data <- as.matrix(data)
    T_obs  <- nrow(data)
    n_obs  <- ncol(data)

    if (is.null(obs_names)) {
      obs_names <- if (!is.null(colnames(data))) colnames(data)
      else paste0("obs_", seq_len(n_obs))
    }
    colnames(data) <- obs_names

    # Build x-axis: dates or integer index
    x_vals <- if (!is.null(dates)) {
      as.Date(dates)
    } else {
      seq_len(T_obs)
    }
    use_dates <- inherits(x_vals, "Date")

    # ---- Summary statistics ------------------------------------------------
    sds     <- apply(data, 2, sd,     na.rm = TRUE)
    means   <- apply(data, 2, mean,   na.rm = TRUE)
    medians <- apply(data, 2, median, na.rm = TRUE)
    mins    <- apply(data, 2, min,    na.rm = TRUE)
    maxs    <- apply(data, 2, max,    na.rm = TRUE)
    n_miss  <- apply(data, 2, function(x) sum(is.na(x)))

    summary_tbl <- data.frame(
      observable = obs_names,
      n          = T_obs - n_miss,
      missing    = n_miss,
      mean       = round(means,   4),
      sd         = round(sds,     4),
      median     = round(medians, 4),
      min        = round(mins,    4),
      max        = round(maxs,    4),
      stringsAsFactors = FALSE
    )

    # ---- Correlation matrix ------------------------------------------------
    cor_mat <- cor(data, use = "pairwise.complete.obs")

    plots <- list()

    if (requireNamespace("ggplot2", quietly = TRUE)) {

      # -- Time-series panel -------------------------------------------------
      long_df <- do.call(rbind, lapply(seq_len(n_obs), function(j) {
        data.frame(
          x   = x_vals,
          y   = data[, j],
          obs = obs_names[j],
          stringsAsFactors = FALSE
        )
      }))

      p_ts <- ggplot2::ggplot(long_df,
                              ggplot2::aes(x = x, y = y)) +
        geom_dynhr_zero() +
        ggplot2::geom_line(colour = dynhr_colours$mid_blue,
                           linewidth = 0.4, na.rm = TRUE) +
        ggplot2::facet_wrap(~ obs, scales = "free_y",
                            ncol = min(3L, ceiling(sqrt(n_obs)))) +
        # Show x-axis date labels on ALL panels (not just the bottom row) so
        # the time context is readable at a glance in any row of the grid.
        { if (use_dates)
            ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y")
          else
            ggplot2::scale_x_continuous(n.breaks = 4L) } +
        theme_dynhr_compact() +
        ggplot2::theme(
          # Rotate x labels slightly so they fit without collision in narrow panels
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                              size = ggplot2::rel(0.75))
        ) +
        ggplot2::labs(
          title    = "Data: Observable time series",
          subtitle = sprintf("%d observations  x  %d variables  x  %d missing",
                             T_obs, n_obs, sum(n_miss)),
          x = if (use_dates) NULL else "Observation"
        )

      p_ts <- .apply_meta(p_ts, meta)
      plots$time_series <- p_ts

      # -- Distribution panel (density per observable) -----------------------
      dens_df <- do.call(rbind, lapply(seq_len(n_obs), function(j) {
        v <- data[, j]
        v <- v[is.finite(v)]
        if (length(v) < 4L) return(NULL)
        kde <- density(v, n = 256)
        data.frame(x = kde$x, density = kde$y,
                   obs = obs_names[j], stringsAsFactors = FALSE)
      }))
      dens_df$obs <- factor(dens_df$obs, levels = rev(obs_names))

      p_dist <- ggplot2::ggplot(dens_df,
                                ggplot2::aes(x = x, fill = obs, colour = obs)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = density),
                              alpha = 0.55, colour = NA) +
        ggplot2::geom_line(ggplot2::aes(y = density), linewidth = 0.4) +
        ggplot2::facet_grid(rows   = ggplot2::vars(obs),
                            scales = "free",
                            switch = "y") +
        # Expand x range so dense distributions (e.g. q_obs) are not clipped.
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.05)) +
        scale_fill_dynhr() +
        scale_colour_dynhr() +
        theme_dynhr_diagnostic() +
        ggplot2::theme(
          legend.position  = "none",
          axis.title.y     = ggplot2::element_blank(),
          axis.text.y      = ggplot2::element_blank(),
          axis.ticks.y     = ggplot2::element_blank(),
          strip.text.y.left = ggplot2::element_text(
            angle = 0, hjust = 1, size = ggplot2::rel(0.80), face = "bold"),
          panel.spacing    = ggplot2::unit(0.15, "lines")
        ) +
        ggplot2::labs(title    = "Data: Observable distributions (ridge)",
                      subtitle = "Empirical kernel density per variable",
                      x        = "Value")

      p_dist <- .apply_meta(p_dist, meta)
      plots$distribution <- p_dist

      # -- Correlation heatmap -----------------------------------------------
      if (n_obs > 1L) {
        cor_long <- do.call(rbind, lapply(seq_len(n_obs), function(i)
          do.call(rbind, lapply(seq_len(n_obs), function(j)
            data.frame(var1 = obs_names[i], var2 = obs_names[j],
                       r    = cor_mat[i, j], stringsAsFactors = FALSE)
          ))
        ))

        p_cor <- ggplot2::ggplot(cor_long,
                                 ggplot2::aes(x = var2, y = var1, fill = r)) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
          ggplot2::geom_text(ggplot2::aes(
            label  = sprintf("%.2f", r),
            colour = ifelse(abs(r) > 0.6, "white", "black")),
            size = 3.0, show.legend = FALSE) +
          scale_fill_dynhr_sunset(limits = c(-1, 1), name = "r") +
          ggplot2::scale_colour_identity() +
          theme_dynhr_diagnostic() +
          ggplot2::theme(
            axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid   = ggplot2::element_blank()
          ) +
          ggplot2::labs(title    = "Data: Pairwise correlations",
                        x = NULL, y = NULL)

        p_cor <- .apply_meta(p_cor, meta)
        plots$correlation <- p_cor
      }
    }

    # Pass: no missing data and all variables have variance
    all_finite  <- all(n_miss == 0)
    all_varying <- all(sds > 0)
    pass <- all_finite && all_varying

    .make_result(
      result  = list(summary = summary_tbl, cor = cor_mat,
                     T_obs = T_obs, n_obs = n_obs),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "Data: %d obs x %d vars. Missing: %d. %s%s",
        T_obs, n_obs, sum(n_miss),
        ifelse(all_finite && all_varying, "PASS.",
               paste(c(if (!all_finite)  "Missing data detected.",
                        if (!all_varying) "Zero-variance variable(s) detected."),
                     collapse = " ")),
        if (!is.null(meta$model_name))
          sprintf(" [%s]", meta$model_name) else ""
      )
    )
}
