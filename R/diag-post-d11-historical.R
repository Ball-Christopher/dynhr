## R/diag-post-d11-historical.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D11 historical decomposition diagnostic.
## NOTE: The exported historical_decomposition() is defined in smoother-monolith.R.
## The matrix-level helper was removed 2026-09 (unused).
## --------------------------------------------------------------------------

## Orientation of a native decomposition's matrices: "time_endo" (T x n_endo,
## what historical_decomposition() returns) or "endo_time" (n_endo x T, what
## historical_decomposition_obc() returns). Stated by `$orientation` when
## present; otherwise read from the dimnames, and refused when ambiguous.
.d11_orientation <- function(hd) {
  if (!is.null(hd$orientation)) {
    if (!hd$orientation %in% c("time_endo", "endo_time"))
      .dynhr_abort(sprintf("D11: unknown decomposition orientation '%s'.",
                           hd$orientation))
    return(hd$orientation)
  }
  tot <- hd$total
  en  <- hd$endo_names
  if (!is.null(colnames(tot)) && is.null(rownames(tot))) return("time_endo")
  if (!is.null(rownames(tot)) && is.null(colnames(tot))) return("endo_time")
  if (!is.null(en)) {
    if (ncol(tot) == length(en) && nrow(tot) != length(en)) return("time_endo")
    if (nrow(tot) == length(en) && ncol(tot) != length(en)) return("endo_time")
  }
  .dynhr_abort(paste0(
    "D11: cannot tell whether the decomposition matrices are T x n_endo or ",
    "n_endo x T. Set `$orientation` (\"time_endo\" or \"endo_time\")."))
}

## Native dynhr decomposition -> per-variable [component x T] list plus the
## series the components must add up to. `actual` is the smoother's own path
## (`$smoothed`, rebuilt from its states) when available -- an independent
## series -- and otherwise `$total`, which is the component sum by
## construction and so checks nothing.
.d11_from_native <- function(hd, var_names, max_vars) {
  orient <- .d11_orientation(hd)
  as_te  <- function(M) if (orient == "time_endo") M else t(M)
  contribs <- lapply(hd$contributions, function(C) as_te(as.matrix(C)))
  total    <- as_te(as.matrix(hd$total))
  en <- colnames(total) %||% hd$endo_names
  if (is.null(en) || length(en) != ncol(total))
    .dynhr_abort(sprintf("D11: decomposition has %d variables but %d names.",
                         ncol(total), length(en)))
  for (nm in names(contribs)) {
    if (!identical(dim(contribs[[nm]]), dim(total)))
      .dynhr_abort(sprintf("D11: component '%s' is %s, total is %s.", nm,
                           paste(dim(contribs[[nm]]), collapse = " x "),
                           paste(dim(total), collapse = " x ")))
  }
  sn <- names(contribs)
  if (is.null(sn))
    sn <- hd$exo_names %||% paste0("shock_", seq_along(contribs))

  independent <- is.matrix(hd$smoothed) && identical(dim(hd$smoothed), dim(total))
  actual <- if (independent) hd$smoothed else total
  colnames(actual) <- en

  keep <- seq_along(en)
  if (!is.null(var_names)) {
    keep <- which(en %in% var_names)
    if (length(keep) == 0L) keep <- seq_along(en)
  }
  if (length(keep) > max_vars) {
    v <- apply(actual[, keep, drop = FALSE], 2L, stats::var)
    keep <- keep[order(v, decreasing = TRUE)[seq_len(max_vars)]]
  }

  mats <- stats::setNames(lapply(keep, function(vi) {
    mat <- do.call(rbind, lapply(contribs, function(C) C[, vi]))
    rownames(mat) <- sn
    mat
  }), en[keep])

  ## Adding-up, per period, over the plotted variables.
  resid_by_period <- if (independent) {
    comp_sum <- Reduce(`+`, contribs)
    apply(abs(comp_sum[, keep, drop = FALSE] - actual[, keep, drop = FALSE]),
          1L, max)
  } else NULL

  list(mats = mats, actual = actual[, keep, drop = FALSE],
       var_names = en[keep], shock_names = sn,
       resid_by_period = resid_by_period,
       tol = hd$adding_up_tol %||% 1e-8,
       transition_residual = hd$transition_residual,
       transition_ok = hd$transition_ok)
}

## "2003-Q2" label for a Date, or the period number otherwise.
.d11_period_label <- function(x) {
  if (inherits(x, "Date"))
    sprintf("%d-Q%d", as.integer(format(x, "%Y")),
            (as.integer(format(x, "%m")) - 1L) %/% 3L + 1L)
  else sprintf("period %d", as.integer(x))
}

#' D11. Historical decomposition
#'
#' Stacked-bar plots of the historical decomposition, one bar segment per
#' shock plus the initial-condition (and OBC constraint) components, with the
#' series they add up to drawn as a line on top. Positive and negative
#' contributions stack on their own side of zero, so in every period the line
#' sits at (sum of the bars above zero) + (sum of the bars below zero).
#'
#' For a native decomposition carrying the smoother's own path
#' (\code{$smoothed}, i.e. \code{historical_decomposition()} called with the
#' whole smoother result), the adding-up of the components to that path is
#' checked in every period and reported; a failure, or a failed
#' \code{$transition_ok}, sets \code{pass = FALSE}. Otherwise the diagnostic
#' is informational (\code{pass = NA}).
#'
#' @param hd_data      Data frame with columns: date, variable, shock, value.
#'                     Or a list of matrices (variable -> shocks x time).
#'                     Or a native dynhr decomposition
#'                     (\code{historical_decomposition()} /
#'                     \code{historical_decomposition_obc()}): a list with
#'                     \code{$contributions}, \code{$total} and optionally
#'                     \code{$smoothed}, \code{$orientation},
#'                     \code{$endo_names}, \code{$exo_names}.
#' @param dates        Date vector, one per period. Optional for list input:
#'                     without it the x axis is the period index.
#' @param var_names    Variables to plot (optional; unknown names ignored).
#' @param shock_names  Character vector (optional)
#' @param show_events  Logical -- annotate NZ macroeconomic events on a Date
#'                     axis (default TRUE)
#' @param max_vars     Most-variable variables kept when there are more.
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
    native <- NULL
    actual_long <- NULL

    if (is.list(hd_data) && !is.data.frame(hd_data) &&
        is.list(hd_data$contributions) && is.matrix(hd_data$total)) {
      native    <- .d11_from_native(hd_data, var_names, max_vars)
      hd_data   <- native$mats
      var_names <- native$var_names
    }

    if (is.list(hd_data) && !is.data.frame(hd_data)) {
      n_T <- unique(vapply(hd_data, function(m) ncol(as.matrix(m)), integer(1)))
      if (length(n_T) != 1L || n_T < 1L) {
        return(.make_result(
          pass    = NA,
          plots   = list(),
          summary = "D11 skipped: hd matrices do not share one period count."
        ))
      }
      if (!is.null(dates) && length(dates) != n_T) {
        return(.make_result(
          pass    = NA,
          plots   = list(),
          summary = sprintf(
            "D11 skipped: hd has %d periods but dates has length %d. Cannot align time periods.",
            n_T, length(dates))
        ))
      }
      x_vals <- if (is.null(dates)) seq_len(n_T) else as.Date(dates)

      dfs <- list()
      for (v in names(hd_data)) {
        mat <- as.matrix(hd_data[[v]])
        s_names <- rownames(mat) %||% paste0("shock_", seq_len(nrow(mat)))
        for (s in seq_len(nrow(mat))) {
          dfs[[length(dfs) + 1L]] <- data.frame(
            date = x_vals, variable = v, shock = s_names[s],
            value = mat[s, ], stringsAsFactors = FALSE)
        }
        act_v <- if (!is.null(native)) native$actual[, v] else colSums(mat)
        actual_long[[v]] <- data.frame(date = x_vals, variable = v,
                                       value = as.numeric(act_v),
                                       stringsAsFactors = FALSE)
      }
      hd_data     <- do.call(rbind, dfs)
      actual_long <- do.call(rbind, actual_long)
    } else {
      if (!inherits(hd_data$date, "Date")) hd_data$date <- as.Date(hd_data$date)
      actual_long <- stats::aggregate(value ~ date + variable, data = hd_data,
                                      FUN = sum)
    }

    if (is.null(var_names))   var_names   <- unique(hd_data$variable)
    if (is.null(shock_names)) shock_names <- unique(hd_data$shock)
    is_date <- inherits(hd_data$date, "Date")

    # ---- Adding-up check -------------------------------------------------
    resid <- if (!is.null(native$resid_by_period)) native$resid_by_period else NULL
    add_up <- NULL
    pass   <- NA
    if (!is.null(resid)) {
      scale <- max(1, max(abs(native$actual)))
      rel   <- max(resid) / scale
      ok    <- rel <= native$tol
      add_up <- list(residual_by_period = resid, residual = max(resid),
                     relative = rel, ok = ok,
                     worst_period = which.max(resid),
                     transition_residual = native$transition_residual,
                     transition_ok = native$transition_ok)
      if (!ok || isFALSE(native$transition_ok)) pass <- FALSE
    }

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
    hd_data$shock <- factor(hd_data$shock, levels = names(fill_values))

    actual_label <- if (!is.null(add_up)) "Smoothed series" else "Sum of components"
    x_lab   <- if (is_date) NULL else "Period"
    y_lab   <- "Contribution (deviation from steady state)"
    bar_w   <- if (is_date) 0.9 * min(diff(sort(unique(as.numeric(hd_data$date)))), 91) else 0.9
    if (!is.finite(bar_w)) bar_w <- if (is_date) 80 else 0.9

    .layers <- function(p, lw, act) {
      p +
        # Stacked bars: ggplot2 stacks positive and negative values on their
        # own side of zero, so the line below equals top-of-positives plus
        # bottom-of-negatives in every period.
        ggplot2::geom_col(position = "stack", alpha = 0.9, width = bar_w) +
        ggplot2::geom_hline(yintercept = 0,
                            colour = dynhr_colours$grey, linewidth = 0.4) +
        ggplot2::geom_line(data = act,
                           ggplot2::aes(x = date, y = value,
                                        linetype = actual_label),
                           inherit.aes = FALSE, colour = "black",
                           linewidth = lw) +
        ggplot2::scale_fill_manual(name = "Component", values = fill_values,
                                   drop = FALSE, na.value = dynhr_na_fill) +
        ggplot2::scale_linetype_manual(name = NULL,
                                       values = stats::setNames("solid", actual_label))
    }

    subtitle <- if (!is.null(add_up)) {
      if (add_up$ok)
        sprintf("Bars add up to the smoothed series in every period (max residual %.2g)",
                add_up$residual)
      else
        sprintf("ADDING-UP FAILED: max residual %.3g at %s",
                add_up$residual, .d11_period_label(
                  sort(unique(hd_data$date))[add_up$worst_period]))
    } else {
      "Line = sum of the plotted components"
    }

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    for (v in var_names) {
      hd_v <- hd_data[hd_data$variable == v, ]
      act_v <- actual_long[actual_long$variable == v, ]
      p_base <- .layers(
        ggplot2::ggplot(hd_v, ggplot2::aes(x = date, y = value, fill = shock)),
        0.6, act_v)
      p_base <- p_base +
        theme_dynhr_diagnostic() +
        ggplot2::labs(title = sprintf("D11: Historical decomposition -- %s", v),
                      subtitle = subtitle, x = x_lab, y = y_lab)

      p_events <- p_base
      if (show_events && is_date) {
        events <- .nz_events()
        date_range <- range(hd_v$date, na.rm = TRUE)
        events <- events[events$date >= date_range[1] & events$date <= date_range[2], ]
        if (nrow(events) > 0) p_events <- .add_nz_event_markers(p_base, events)
      }
      plots[[paste0("hd_", v)]] <- .apply_meta(p_events, meta)
    }

    # Consolidated overview (no event markers: they clutter small panels).
    n_cols_overview <- min(3L, length(var_names))
    p_overview <- .layers(
      ggplot2::ggplot(hd_data[hd_data$variable %in% var_names, ],
                      ggplot2::aes(x = date, y = value, fill = shock)),
      0.4, actual_long[actual_long$variable %in% var_names, ])
    p_overview <- p_overview +
      ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = n_cols_overview) +
      theme_dynhr_compact() +
      ggplot2::theme(
        axis.text.x   = ggplot2::element_text(size = ggplot2::rel(0.75)),
        # free_y panels: without tick labels their scales cannot be read.
        axis.text.y   = ggplot2::element_text(size = ggplot2::rel(0.7)),
        axis.ticks.y  = ggplot2::element_line(colour = dynhr_colours$grey),
        axis.title.y  = ggplot2::element_text(size = ggplot2::rel(0.8), angle = 90),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(0.4, "lines")
      ) +
      ggplot2::labs(
        title    = "D11: Historical decomposition overview",
        subtitle = sprintf("%d variables, %d components. %s",
                           length(var_names), length(shock_names), subtitle),
        x = x_lab, y = y_lab)
    attr(p_overview, "dynhr_fig_height") <- min(9.0, 2.5 + 2.2 * ceiling(length(var_names) / 3))
    plots$hd_overview <- .apply_meta(p_overview, meta)

    }  # end requireNamespace guard

    periods <- sort(unique(hd_data$date))
    n_per   <- length(periods)
    rng     <- sprintf("%s to %s", .d11_period_label(periods[1]),
                       .d11_period_label(periods[n_per]))
    add_txt <- if (is.null(add_up)) {
      " Adding-up not checked (no independent smoothed series supplied)."
    } else {
      sprintf(" Adding-up to the smoothed series: max residual %.3g (%s)%s.",
              add_up$residual, if (add_up$ok) "OK" else "FAILED",
              if (is.null(add_up$transition_residual)) "" else
                sprintf("; transition residual %.3g", add_up$transition_residual))
    }
    events_txt <- if (show_events && is_date) " NZ event markers included." else ""
    badge <- if (isFALSE(pass)) "FAIL" else "INFO"

    .make_result(
      result  = list(hd_data = hd_data, actual = actual_long,
                     var_names = var_names, shock_names = shock_names,
                     adding_up = add_up),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "D11 Historical decomposition: %d variables, %d components, %d periods (%s).%s%s",
        length(var_names), length(shock_names), n_per, rng, add_txt, events_txt),
      llm_summary = paste(c(
        sprintf("D11 | Historical Decomposition | %s", badge),
        sprintf("  variables=%d components=%d periods=%d",
                length(var_names), length(shock_names), n_per),
        sprintf("  range=%s", rng),
        if (!is.null(add_up))
          sprintf("  adding_up_residual=%.3g ok=%s transition_residual=%s",
                  add_up$residual, add_up$ok,
                  format(add_up$transition_residual %||% NA_real_, digits = 3)),
        if (isFALSE(pass))
          "  action: components do not add up to the smoothed series; check the smoother inputs ($transition_residual_by_period)."
        else
          "  action: Informational only. Inspect plots for dominant shock contributors by episode."
      ), collapse = "\n")
    )
}
