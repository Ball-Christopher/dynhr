## R/diag-post-d10-vardec.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D10 variance decomposition; .get_vd_long() helper
## --------------------------------------------------------------------------

.get_vd_long <- function(vd_pct) {
  if (is.list(vd_pct) && !is.data.frame(vd_pct)) {
    if (!is.null(vd_pct$moments) && is.list(vd_pct$moments))
      vd_pct <- vd_pct$moments
    if (is.list(vd_pct) && !is.null(vd_pct$var_decomp_pct))
      vd_pct <- vd_pct$var_decomp_pct
  }
  if (!is.matrix(vd_pct) && !is.data.frame(vd_pct))
    stop("vd_pct must be a matrix [n_var x n_shock]")
  vd_pct <- as.matrix(vd_pct)
  var_names   <- rownames(vd_pct) %||% paste0("V", seq_len(nrow(vd_pct)))
  shock_names <- colnames(vd_pct) %||% paste0("S", seq_len(ncol(vd_pct)))
  rows <- list()
  for (i in seq_along(var_names))
    for (j in seq_along(shock_names))
      rows[[length(rows) + 1L]] <- data.frame(
        horizon = Inf, variable = var_names[i],
        shock = shock_names[j], share = vd_pct[i, j] / 100,
        stringsAsFactors = FALSE
      )
  do.call(rbind, rows)
}


#' D10. Variance decomposition
#'
#' Displays the unconditional forecast error variance decomposition and, when
#' foreign shocks are present, checks whether their share of output variance
#' falls within the empirically plausible range for small open-economy DSGE
#' models.
#'
#' Pass gate (foreign-shock check): the share of output variance attributed to
#' foreign shocks must lie in [0.15, 0.60].  This range is based on the
#' empirical range documented by Kamber et al. (2016) and Justiniano & Preston
#' (2010) for SOE DSGE models (typically 20-50\%).  The earlier gate of
#' [0.05, 0.70] was an excessively wide 65pp window that accepted almost any
#' parameterisation; the tighter gate catches models with trivial foreign
#' influence (<15\%) or overwhelming foreign dominance (>60\%).
#' \code{pass = NA} when no foreign shocks are identified.
#'
#' @param vd_data       Variance decomposition data.  Accepted formats:
#'   \enumerate{
#'     \item A matrix [n_var x n_shock] of variance shares (rows sum to 100\%
#'       or 1.0 per variable).
#'     \item A data frame with columns \code{variable}, \code{shock},
#'       \code{share} (shares in [0, 1]).
#'     \item A model moment object with \code{$moments$var_decomp_pct} slot.
#'   }
#' @param var_names     Character vector of variable names (optional; inferred
#'   from \code{vd_data} when possible).
#' @param shock_names   Character vector of shock names (optional; inferred
#'   from \code{vd_data} when possible).
#' @param foreign_shocks Character vector of shock names to treat as foreign
#'   (for the [0.15, 0.60] output-variance gate).  If NULL, the diagnostic
#'   auto-detects common foreign shock name patterns
#'   (e.g. \code{eps_ys}, \code{eps_ps}, \code{eps_rs}).
#' @param meta          Optional metadata list for plot annotation.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing \code{vd_long} (long-format data frame),
#'     \code{foreign_check} (list with \code{foreign_share}, \code{benchmark},
#'     \code{plausible}, and \code{message}; NULL when no foreign shocks).}
#'   \item{pass}{Logical -- foreign share in [0.15, 0.60]; \code{NA} when
#'     no foreign shocks are present.}
#'   \item{plots}{ggplot2 stacked bar chart and heatmap of variance shares.}
#'   \item{summary}{Human-readable summary with foreign share statistics.}
#'
#' @references
#'   Sims, C. A. (1980). Macroeconomics and reality. \emph{Econometrica},
#'     48(1), 1-48.
#'   Kamber, G., Morley, J., & Wong, B. (2016). Assessing the periphery's role
#'     in the global business cycle in a data-rich environment.
#'     \emph{North American Journal of Economics and Finance}, 38, 171-188.
#'   Justiniano, A., & Preston, B. (2010). Can structural small open-economy
#'     models account for the influence of foreign disturbances?
#'     \emph{Journal of International Economics}, 81(1), 61-74.
#' @noRd
d10_variance_decomposition <- function(vd_data, var_names = NULL,
                                       shock_names = NULL,
                                       foreign_shocks = NULL,
                                       meta = NULL) {

    vd_long <- if (is.data.frame(vd_data) &&
                   all(c("variable","shock","share") %in% names(vd_data))) {
      vd_data
    } else {
      .get_vd_long(vd_data)
    }

    all_vars   <- unique(vd_long$variable)
    all_shocks <- unique(vd_long$shock)

    if (is.null(foreign_shocks))
      foreign_shocks <- intersect(all_shocks,
                                  c("eps_ys","eps_ps","eps_rs","eps_ys_","eps_ps_","eps_rs_"))

    foreign_check <- NULL
    if (length(foreign_shocks) > 0 && "y" %in% all_vars) {
      y_rows <- vd_long[vd_long$variable == "y", ]
      foreign_share <- sum(y_rows$share[y_rows$shock %in% foreign_shocks])
      # Gate: [0.15, 0.60] per Kamber et al. (2016) and Justiniano & Preston (2010).
      # Empirical SOE DSGE models assign 20–50% of output variance to foreign
      # shocks.  The earlier [0.05, 0.70] window was a 65 pp acceptance band that
      # passed almost every reasonable model.  The tighter [0.15, 0.60] window
      # catches models with trivial foreign influence (<15%) or overwhelming
      # foreign dominance (>60%).
      foreign_check <- list(
        foreign_share = foreign_share,
        benchmark     = "[0.15, 0.60] (Kamber et al. 2016; Justiniano & Preston 2010)",
        plausible     = foreign_share > 0.15 && foreign_share < 0.60,
        message       = sprintf("Foreign share of y variance: %.1f%% (plausible range 15%%--60%%)",
                                foreign_share * 100)
      )
    }

    pass <- if (!is.null(foreign_check)) foreign_check$plausible else NA
    summary_lines <- c(
      sprintf("D10 Variance decomposition: %d variables, %d shocks",
              length(all_vars), length(all_shocks)),
      if (!is.null(foreign_check)) paste(" ", foreign_check$message)
    )

    var_names   <- all_vars
    shock_names <- all_shocks
    vd_mat <- tapply(vd_long$share, list(vd_long$shock, vd_long$variable), sum)
    vd_mat[is.na(vd_mat)] <- 0
    foreign_share <- if (!is.null(foreign_check)) foreign_check$foreign_share else NULL

    # Row-sum sanity check: each variable's shares should sum to ~1
    # vd_mat rows=shocks, cols=variables; column sums = shares per variable
    col_sums <- colSums(vd_mat)
    bad_vars <- names(col_sums)[abs(col_sums - 1) > 0.005]
    if (length(bad_vars) > 0) {
      summary_lines <- c(summary_lines,
        sprintf("  WARNING: row shares do not sum to 1 for: %s",
                paste(bad_vars, collapse = ", ")))
    }

    # Build 100%-stacked bar chart + heatmap table
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # vd_long already has variable, shock, share columns
      vd_plot_df <- vd_long[, c("variable", "shock", "share")]
      colnames(vd_plot_df) <- c("Variable", "Shock", "Share")

      p_vd <- ggplot2::ggplot(vd_plot_df,
                              ggplot2::aes(x = Variable, y = Share, fill = Shock)) +
        ggplot2::geom_col(position = "fill") +
        ggplot2::geom_text(
          ggplot2::aes(
            label  = ifelse(Share > 0.02, sprintf("%.0f%%", Share * 100), "")
          ),
          colour = "grey15",   # dark text reads well on the light pastel fills
          position = ggplot2::position_fill(vjust = 0.5),
          size = 2.8,
          show.legend = FALSE
        ) +
        scale_fill_dynhr_light() +
        ggplot2::scale_y_continuous(
          expand = c(0, 0),
          labels = scales::percent
        ) +
        theme_dynhr() +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "D10: Variance decomposition",
          x = NULL, y = "Share"
        )
      plots$vardec <- .apply_meta(p_vd, meta)

      # Heatmap / coloured-table alternative
      p_tbl <- ggplot2::ggplot(vd_plot_df,
                               ggplot2::aes(x = Shock, y = Variable, fill = Share)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::geom_text(
          ggplot2::aes(
            label  = ifelse(Share > 0.02, sprintf("%.0f%%", Share * 100), ""),
            # cividis is DARK at low share, light/yellow at high share, so low
            # cells need white text and high cells need dark text.
            colour = Share < 0.5
          ),
          size = 2.8,
          show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(
          values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
        scale_fill_dynhr_cividis(labels = scales::percent) +
        theme_dynhr() +
        ggplot2::theme(
          panel.grid      = ggplot2::element_blank(),
          axis.line       = ggplot2::element_blank(),
          axis.ticks      = ggplot2::element_blank(),
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(
          title = "D10: Variance decomposition (heatmap)",
          x = "Shock", y = "Variable", fill = "Share"
        )
      plots$vardec_table <- .apply_meta(p_tbl, meta)
    }

    .make_result(
      result  = list(vd_long = vd_long, foreign_check = foreign_check),
      pass    = pass,
      plots   = plots,
      summary = paste(summary_lines, collapse = "\n"),
      llm_summary = {
        badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
        n_var <- length(var_names); n_shk <- length(shock_names)
        top_df <- which(vd_mat == max(vd_mat), arr.ind = TRUE)
        top_shock <- rownames(vd_mat)[top_df[1,1]]
        top_var   <- colnames(vd_mat)[top_df[1,2]]
        paste(c(
          sprintf("D10 | Variance Decomposition | %s", badge),
          sprintf("  variables=%d shocks=%d", n_var, n_shk),
          if (!is.null(foreign_share))
            sprintf("  foreign_share_y=%.2f (benchmark ~0.33)", foreign_share),
          sprintf("  largest_contributor: %s->%s=%.2f", top_shock, top_var, max(vd_mat)),
          sprintf("  action: %s",
                  if (isTRUE(pass)) "Foreign shock share within plausible range [15%, 60%] (Kamber et al. 2016)."
                  else sprintf("Foreign shock share %.1f%% outside [15%%, 60%%]. %s",
                               foreign_share * 100,
                               if (!is.null(foreign_share) && foreign_share < 0.15)
                                 "Too little foreign influence -- check open-economy transmission channels."
                               else
                                 "Foreign shocks dominate output -- consider whether domestic shocks are properly identified."))
        ), collapse = "\n")
      }
    )
}
