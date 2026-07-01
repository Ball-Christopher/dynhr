## R/diag-shock-dominant.R
## --------------------------------------------------------------------------
## Phase-3+ addition.
##
## diag_shock_dominant() -- per-period dominant shock identification and
## shock decomposition summary from historical_decomposition() output.
## --------------------------------------------------------------------------

#' Dominant shock identification and shock decomposition diagnostic
#'
#' For each time period and selected variable, identifies which structural
#' shock (or shock category) made the largest absolute contribution to the
#' model-implied path.  Produces:
#' \itemize{
#'   \item A stacked contribution bar chart with the dominant shock
#'     highlighted per period.
#'   \item A "dominant-shock timeline" showing which shock dominated at each
#'     date.
#'   \item A summary table: fraction of periods dominated by each shock.
#' }
#'
#' @param hist_decomp    Output of \code{historical_decomposition()}: a list
#'   with \code{$contributions} (named list of T?--n_endo matrices, one per
#'   shock) and \code{$total} (T?--n_endo matrix).
#' @param var_sel        Character vector of variable names to analyse.
#'   Default: all variables in \code{hist_decomp$total}.
#' @param dates          Optional Date vector of length T for the x-axis.
#' @param shock_categories Named list mapping category label to shock names.
#'   When provided, contributions of grouped shocks are summed before
#'   computing dominance; category names appear in plots and tables.
#' @param top_n_shocks   Number of shocks/categories to show in the stacked
#'   chart; remaining are collapsed into "Other" (default 6).
#' @param meta           A \code{\link{diag_meta}} object for plot provenance.
#' @return dynhr_diagnostic list.  \code{$result$dominance} is a named list
#'   (one entry per variable) of data frames with columns
#'   \code{period}, \code{dominant_shock}, \code{share}.
#' @noRd
diag_shock_dominant <- function(hist_decomp,
                                 var_sel          = NULL,
                                 dates            = NULL,
                                 shock_categories = list(),
                                 top_n_shocks     = 6L,
                                 meta             = NULL) {
  
    total_mat   <- hist_decomp$total
    contrib_raw <- hist_decomp$contributions   # named list of T?--n_endo matrices
    T_obs       <- nrow(total_mat)
    endo_names  <- colnames(total_mat)
    shock_names <- names(contrib_raw)

    if (is.null(var_sel)) var_sel <- endo_names
    var_sel <- intersect(var_sel, endo_names)

    x_vals    <- if (!is.null(dates)) as.Date(dates)
                 else seq_len(T_obs)
    use_dates <- inherits(x_vals, "Date")

    # ---- Resolve shock categories ------------------------------------------
    # If shock_categories supplied, build aggregated contributions
    if (length(shock_categories) > 0) {
      agg_contrib <- lapply(names(shock_categories), function(cat_nm) {
        members <- shock_categories[[cat_nm]]
        # Clean up: split comma-separated strings
        if (length(members) == 1L && grepl(",", members))
          members <- trimws(strsplit(members, ",")[[1]])
        members <- intersect(trimws(members), shock_names)
        if (length(members) == 0L) return(NULL)
        # Sum matrices
        Reduce("+", contrib_raw[members])
      })
      names(agg_contrib) <- names(shock_categories)
      # Keep uncategorised shocks individually
      all_cat_members <- unlist(lapply(shock_categories, function(m) {
        if (length(m) == 1L && grepl(",", m)) trimws(strsplit(m, ",")[[1]]) else trimws(m)
      }))
      uncategorised <- setdiff(shock_names, all_cat_members)
      contrib_use <- c(Filter(Negate(is.null), agg_contrib),
                       contrib_raw[uncategorised])
    } else {
      contrib_use <- contrib_raw
    }
    label_names <- names(contrib_use)

    # ---- Per-variable analysis ---------------------------------------------
    dominance_list <- list()
    plots          <- list()

    for (v in var_sel) {
      v_idx <- match(v, endo_names)
      if (is.na(v_idx)) next

      # T ?-- n_labels matrix of contributions for variable v
      contrib_v <- do.call(cbind, lapply(contrib_use, function(mat) mat[, v_idx]))
      colnames(contrib_v) <- label_names

      # Dominant shock per period: largest absolute contribution
      dom_idx   <- apply(abs(contrib_v), 1, which.max)
      dom_shock <- label_names[dom_idx]
      dom_share <- abs(contrib_v[cbind(seq_len(T_obs), dom_idx)]) /
                   (rowSums(abs(contrib_v)) + 1e-300)

      dominance_list[[v]] <- data.frame(
        period         = seq_len(T_obs),
        x              = x_vals,
        dominant_shock = dom_shock,
        share          = dom_share,
        stringsAsFactors = FALSE
      )

      # ---- Plots -----------------------------------------------------------
      if (requireNamespace("ggplot2", quietly = TRUE)) {

        # Identify top-N shocks by total absolute contribution
        total_abs  <- colSums(abs(contrib_v))
        top_labels <- names(sort(total_abs, decreasing = TRUE))[seq_len(min(top_n_shocks, length(label_names)))]
        other_idx  <- !(label_names %in% top_labels)

        contrib_plot <- contrib_v
        if (any(other_idx)) {
          contrib_plot <- cbind(contrib_v[, top_labels, drop = FALSE],
                                Other = rowSums(contrib_v[, other_idx, drop = FALSE]))
          plot_labels <- c(top_labels, "Other")
        } else {
          contrib_plot <- contrib_v[, top_labels, drop = FALSE]
          plot_labels  <- top_labels
        }

        n_plot_labels <- length(plot_labels)
        pal <- if (n_plot_labels <= length(dynhr_palette))
          dynhr_palette[seq_len(n_plot_labels)]
        else colorRampPalette(dynhr_palette)(n_plot_labels)
        names(pal) <- plot_labels

        # Long format for ggplot
        long_df <- do.call(rbind, lapply(seq_along(plot_labels), function(si) {
          sh <- plot_labels[si]
          data.frame(x = x_vals, value = contrib_plot[, sh],
                     shock = sh, stringsAsFactors = FALSE)
        }))
        long_df$shock <- factor(long_df$shock, levels = rev(plot_labels))

        # Panel 1: stacked contribution chart
        p_stack <- ggplot2::ggplot(long_df, ggplot2::aes(x = x, y = value, fill = shock)) +
          ggplot2::geom_col(position = "stack", width = if (use_dates) NULL else 1,
                            alpha = 0.85) +
          ggplot2::geom_hline(yintercept = 0, colour = "white", linewidth = 0.3) +
          ggplot2::scale_fill_manual(values = rev(pal), name = "Shock") +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title    = sprintf("Shock decomposition: %s", v),
            subtitle = sprintf("Top %d shocks/categories * stacked contributions",
                               min(top_n_shocks, length(label_names))),
            x = if (use_dates) NULL else "Period",
            y = sprintf("Contribution to %s", v)
          )
        p_stack <- .apply_meta(p_stack, meta)
        plots[[paste0("stacked_", v)]] <- p_stack

        # Panel 2: dominant-shock timeline (colour = shock, y = share)
        dom_df <- dominance_list[[v]]
        dom_df$dominant_shock <- factor(dom_df$dominant_shock,
                                        levels = c(top_labels, setdiff(label_names, top_labels)))

        # Colour map for dominant shocks (use same palette as above)
        dom_pal <- c(pal, setNames(rep(dynhr_colours$grey,
                                       length(setdiff(label_names, top_labels))),
                                   setdiff(label_names, top_labels)))

        p_dom <- ggplot2::ggplot(dom_df, ggplot2::aes(x = x, y = share,
                                                        fill = dominant_shock,
                                                        colour = dominant_shock)) +
          ggplot2::geom_col(width = if (use_dates) NULL else 1) +
          ggplot2::scale_fill_manual(values = dom_pal, name = "Dominant shock") +
          ggplot2::scale_colour_manual(values = dom_pal, guide = "none") +
          ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", x * 100),
                                      limits = c(0, 1)) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title    = sprintf("Dominant shock: %s", v),
            subtitle = "Fraction of total abs. contribution explained by dominant shock",
            x        = if (use_dates) NULL else "Period",
            y        = "Dominant share"
          )
        p_dom <- .apply_meta(p_dom, meta)
        plots[[paste0("dominant_", v)]] <- p_dom
      }
    }

    # ---- Dominance summary table -------------------------------------------
    dom_summary <- do.call(rbind, lapply(names(dominance_list), function(v) {
      d <- dominance_list[[v]]
      freq <- table(d$dominant_shock) / nrow(d)
      data.frame(
        variable = v,
        shock    = names(freq),
        pct_dominant = as.numeric(freq) * 100,
        stringsAsFactors = FALSE
      )
    }))
    dom_summary <- dom_summary[order(dom_summary$variable, -dom_summary$pct_dominant), ]

    # Top dominant per variable for the summary string
    top_per_var <- do.call(rbind, lapply(names(dominance_list), function(v) {
      d <- dominance_list[[v]]
      freq <- sort(table(d$dominant_shock), decreasing = TRUE)
      data.frame(variable = v,
                 top_shock = names(freq)[1],
                 pct = as.numeric(freq[1]) / nrow(d) * 100,
                 stringsAsFactors = FALSE)
    }))

    summary_txt <- paste(
      sprintf("Dominant shocks: %d variable(s). Top: %s",
              length(var_sel),
              paste(sprintf("%s->%s(%.0f%%)", top_per_var$variable,
                            top_per_var$top_shock, top_per_var$pct),
                    collapse = ", ")),
      collapse = "\n"
    )

    .make_result(
      result  = list(dominance     = dominance_list,
                     dom_summary   = dom_summary,
                     var_sel       = var_sel,
                     label_names   = label_names),
      pass    = NA,   # informational
      plots   = plots,
      summary = summary_txt
    )
}
