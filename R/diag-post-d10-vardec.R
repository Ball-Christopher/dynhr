## R/diag-post-d10-vardec.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D10 variance decomposition; .get_vd_long() / .d10_parse_input() helpers
## --------------------------------------------------------------------------

## Matrix [n_var x n_shock] -> long data frame with shares in [0, 1].
## `scale` is what one variable's row sums to: 100 for percentages, 1 for
## fractions, NULL to infer it from the row sums. (The matrix is documented as
## either; until 0.9.4 it was always divided by 100, so a fraction matrix came
## out 100x too small.)
.get_vd_long <- function(vd_pct, scale = NULL) {
  if (!is.matrix(vd_pct) && !is.data.frame(vd_pct))
    .dynhr_abort("vd_pct must be a matrix [n_var x n_shock]")
  vd_pct <- as.matrix(vd_pct)
  if (is.null(scale)) {
    rs <- rowSums(vd_pct)
    rs <- rs[is.finite(rs) & rs > 0]
    scale <- if (length(rs) > 0L && all(abs(rs - 1) < 5e-3)) 1
             else if (length(rs) > 0L && all(abs(rs - 100) < 0.5)) 100
             else if (any(vd_pct > 1 + 1e-8, na.rm = TRUE)) 100
             else 1
  }
  var_names   <- rownames(vd_pct) %||% paste0("V", seq_len(nrow(vd_pct)))
  shock_names <- colnames(vd_pct) %||% paste0("S", seq_len(ncol(vd_pct)))
  ## Row-major (variable, then shock): as.vector(t(.)) walks the rows.
  data.frame(
    horizon  = Inf,
    variable = rep(var_names, each = length(shock_names)),
    shock    = rep(shock_names, times = length(var_names)),
    share    = as.vector(t(vd_pct)) / scale,
    stringsAsFactors = FALSE
  )
}

## Normalise every accepted D10 input to a long data frame, plus whatever side
## information (Sigma_e, absolute contributions, var_cov) the input carries so
## the correlated-shock caveat can be quantified.
.d10_parse_input <- function(vd_data) {
  info <- list(vd_long = NULL, Sigma_e = NULL, var_decomp = NULL,
               var_cov = NULL, horizon_note = NULL)
  if (is.list(vd_data) && !is.data.frame(vd_data)) {
    if (!is.null(vd_data$moments) && is.list(vd_data$moments))
      vd_data <- vd_data$moments
    if (is.null(vd_data$var_decomp_pct))
      .dynhr_abort("D10: list input has no $var_decomp_pct ",
                   "(or $moments$var_decomp_pct).")
    info$Sigma_e    <- vd_data$Sigma_e
    info$var_decomp <- vd_data$var_decomp
    info$var_cov    <- vd_data$var_cov
    info$vd_long    <- .get_vd_long(vd_data$var_decomp_pct, scale = 100)
    return(info)
  }
  if (is.data.frame(vd_data) &&
      all(c("variable", "shock", "share") %in% names(vd_data))) {
    vd_long <- as.data.frame(vd_data, stringsAsFactors = FALSE)
    vd_long$variable <- as.character(vd_long$variable)
    vd_long$shock    <- as.character(vd_long$shock)
    ## A multi-horizon FEVD in long form: summing across horizons (as the old
    ## code did) gives shares > 1. Use the unconditional horizon when present,
    ## else the longest one.
    if ("horizon" %in% names(vd_long)) {
      hz <- unique(vd_long$horizon)
      if (length(hz) > 1L) {
        h_use <- if (any(is.infinite(hz))) Inf else max(hz)
        vd_long <- vd_long[vd_long$horizon == h_use, , drop = FALSE]
        info$horizon_note <- sprintf(
          "  Input has %d horizons; using horizon %s.", length(hz),
          format(h_use))
      }
    }
    info$vd_long <- vd_long
    return(info)
  }
  info$vd_long <- .get_vd_long(vd_data)
  info
}

## Fraction of each variable's variance that the per-shock contributions do
## NOT account for.  Since 0.9.4 compute_moments() orthogonalises correlated
## shocks with the Cholesky factor of Sigma_e, so this is zero (to rounding) by
## construction; it is kept as a cheap sanity check on whatever decomposition
## the caller actually handed in.  NULL when the inputs are unavailable.
.d10_cov_gap <- function(var_decomp, var_cov) {
  if (is.null(var_decomp) || is.null(var_cov)) return(NULL)
  vars <- intersect(rownames(var_decomp), rownames(var_cov))
  if (length(vars) == 0L) return(NULL)
  tot  <- diag(as.matrix(var_cov))[vars]
  part <- rowSums(as.matrix(var_decomp)[vars, , drop = FALSE])
  ok   <- is.finite(tot) & is.finite(part) & tot > 0
  gap  <- 1 - part[ok] / tot[ok]
  gap
}


#' D10. Variance decomposition
#'
#' Displays the unconditional forecast error variance decomposition and, when
#' foreign shocks are present, checks whether their share of output variance
#' falls within the empirically plausible range for small open-economy DSGE
#' models.
#'
#' Pass gate (foreign-shock check): the share of \code{output_var} variance
#' attributed to foreign shocks must lie in [0.15, 0.60].  This range is based
#' on the empirical range documented by Kamber et al. (2016) and Justiniano &
#' Preston (2010) for SOE DSGE models (typically 20-50\%).
#' \code{pass = NA} when no foreign shocks are identified.
#'
#' Correlated shocks: since 0.9.4 \code{compute_moments()} orthogonalises the
#' shocks with the lower Cholesky factor of \code{Sigma_e} in the DECLARED
#' shock order (Dynare's \code{stoch_simul} convention), so the shares sum to
#' 100\% with no unattributed covariance term.  When the input carries a
#' non-diagonal \code{Sigma_e}, D10 reports an INFO note that the split of the
#' shared variance between correlated shocks depends on the \code{varexo}
#' order, and (when \code{var_cov} is available) the unattributed fraction per
#' variable, which is now zero up to rounding.
#'
#' @param vd_data       Variance decomposition data.  Accepted formats:
#'   \enumerate{
#'     \item A matrix [n_var x n_shock] of variance shares (rows sum to 100\%
#'       or 1.0 per variable; the scale is inferred from the row sums).
#'     \item A data frame with columns \code{variable}, \code{shock},
#'       \code{share} (shares in [0, 1]) and optionally \code{horizon}
#'       (the unconditional / longest horizon is used).
#'     \item A \code{compute_moments()} list (\code{$var_decomp_pct}) or an
#'       object with a \code{$moments$var_decomp_pct} slot.
#'   }
#' @param var_names,shock_names Ignored (names are taken from
#'   \code{vd_data}); kept for call compatibility.
#' @param foreign_shocks Character vector of shock names to treat as foreign.
#'   Names not present in \code{vd_data} are dropped with a warning.  If NULL,
#'   common names are auto-detected (\code{eps_ys}, \code{eps_ps},
#'   \code{eps_rs}).
#' @param output_var    Name of the output variable for the foreign gate.
#' @param meta          Optional metadata list for plot annotation.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing \code{vd_long} (long-format data frame,
#'     shares in [0, 1]), \code{foreign_check} (NULL when no foreign shocks),
#'     \code{row_sums}, \code{no_decomp_vars} (variables with no finite
#'     decomposition, e.g. nonstationary), \code{correlated_shocks} and
#'     \code{cov_gap}.}
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
                                       output_var = "y",
                                       meta = NULL) {

    parsed  <- .d10_parse_input(vd_data)
    vd_long <- parsed$vd_long

    all_vars   <- unique(vd_long$variable)
    all_shocks <- unique(vd_long$shock)

    ## Per-variable row sums; a variable whose shares are all zero / non-finite
    ## has no decomposition (compute_moments() zeroes nonstationary rows).
    row_sums <- vapply(all_vars, function(v) {
      s <- vd_long$share[vd_long$variable == v]
      if (any(!is.finite(s))) NA_real_ else sum(s)
    }, numeric(1))
    no_decomp_vars <- all_vars[is.na(row_sums) | row_sums == 0]
    dec_vars <- setdiff(all_vars, no_decomp_vars)
    bad_vars <- dec_vars[abs(row_sums[dec_vars] - 1) > 0.005]

    gate <- c(0.15, 0.60)
    if (is.null(foreign_shocks)) {
      foreign_shocks <- intersect(all_shocks,
                                  c("eps_ys","eps_ps","eps_rs","eps_ys_","eps_ps_","eps_rs_"))
    } else {
      missing_fs <- setdiff(foreign_shocks, all_shocks)
      if (length(missing_fs) > 0L)
        .dynhr_warn(sprintf(
          "D10: foreign_shocks not in the decomposition, ignored: %s (shocks: %s)",
          paste(missing_fs, collapse = ", "), paste(all_shocks, collapse = ", ")))
      foreign_shocks <- intersect(foreign_shocks, all_shocks)
    }

    foreign_check <- NULL
    gate_note <- NULL
    if (length(foreign_shocks) > 0 && !(output_var %in% dec_vars)) {
      gate_note <- sprintf(
        "  Foreign-share gate not applied: output variable '%s' %s.", output_var,
        if (output_var %in% all_vars) "has no finite decomposition"
        else "is not in the decomposition")
    }
    if (length(foreign_shocks) > 0 && output_var %in% dec_vars) {
      y_rows <- vd_long[vd_long$variable == output_var, ]
      foreign_share <- sum(y_rows$share[y_rows$shock %in% foreign_shocks])
      foreign_check <- list(
        foreign_share  = foreign_share,
        foreign_shocks = foreign_shocks,
        benchmark      = "[0.15, 0.60] (Kamber et al. 2016; Justiniano & Preston 2010)",
        plausible      = foreign_share >= gate[1] && foreign_share <= gate[2],
        message        = sprintf(
          "Foreign share (%s) of %s variance: %.1f%% (plausible range 15%%--60%%)",
          paste(foreign_shocks, collapse = ", "), output_var, foreign_share * 100)
      )
    }

    ## Correlated shocks: shares ignore the covariance terms.
    Se <- parsed$Sigma_e
    correlated <- !is.null(Se) && length(Se) > 1L &&
      any(abs(Se[row(Se) != col(Se)]) > 1e-12 * max(abs(diag(Se)), 1e-300))
    cov_gap <- if (correlated) .d10_cov_gap(parsed$var_decomp, parsed$var_cov) else NULL

    pass <- if (!is.null(foreign_check)) foreign_check$plausible else NA
    summary_lines <- c(
      sprintf("D10 Variance decomposition: %d variables, %d shocks",
              length(all_vars), length(all_shocks)),
      parsed$horizon_note,
      if (!is.null(foreign_check)) paste(" ", foreign_check$message),
      gate_note,
      if (length(bad_vars) > 0)
        sprintf("  WARNING: shares do not sum to 1 for: %s",
                paste(sprintf("%s (%.3f)", bad_vars, row_sums[bad_vars]),
                      collapse = ", ")),
      if (length(no_decomp_vars) > 0)
        sprintf("  No finite decomposition (e.g. nonstationary): %s",
                paste(no_decomp_vars, collapse = ", ")),
      if (correlated)
        paste0("  INFO: Sigma_e has non-zero correlations; shocks are ",
               "Cholesky-orthogonalised in the DECLARED shock order (Dynare's ",
               "convention), so the shares sum to 100% but the split of the ",
               "shared variance depends on the varexo order",
               if (length(cov_gap) > 0)
                 sprintf(". Largest unattributed share: %s %.2f%% (0 by ",
                         names(cov_gap)[which.max(abs(cov_gap))],
                         100 * cov_gap[which.max(abs(cov_gap))])
               else "",
               if (length(cov_gap) > 0) "construction)." else ".")
    )

    ## shocks x variables matrix over the variables that have a decomposition
    vd_dec <- vd_long[vd_long$variable %in% dec_vars, , drop = FALSE]
    vd_mat <- if (nrow(vd_dec) > 0L)
      tapply(vd_dec$share, list(vd_dec$shock, vd_dec$variable), sum) else NULL
    if (!is.null(vd_mat)) vd_mat[is.na(vd_mat)] <- 0
    foreign_share <- if (!is.null(foreign_check)) foreign_check$foreign_share else NULL

    # Build 100%-stacked bar chart + heatmap table
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE) && nrow(vd_dec) > 0L) {
      vd_plot_df <- vd_dec[, c("variable", "shock", "share")]
      colnames(vd_plot_df) <- c("Variable", "Shock", "Share")
      ## Declaration order: shocks left-to-right / in the legend, variables
      ## top-to-bottom (coord_flip draws the first level at the bottom).
      vd_plot_df$Shock    <- factor(vd_plot_df$Shock, levels = all_shocks)
      vd_plot_df$Variable <- factor(vd_plot_df$Variable, levels = rev(dec_vars))

      subtitle <- paste(c(
        "Unconditional (h = Inf) share of each variable's variance, by shock",
        if (!is.null(foreign_check))
          sprintf("foreign share of %s = %.0f%% (gate 15%%-60%%: %s)", output_var,
                  100 * foreign_share,
                  if (isTRUE(pass)) "PASS" else "FAIL"),
        if (length(no_decomp_vars) > 0)
          sprintf("omitted (no finite decomposition): %s",
                  paste(no_decomp_vars, collapse = ", ")),
        if (correlated) "correlated shocks: covariance terms unassigned"
      ), collapse = "\n")

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
          subtitle = subtitle,
          x = NULL, y = "Share of variance"
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
        scale_fill_dynhr_cividis(
          labels = scales::percent, limits = c(0, 1),
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          guide = ggplot2::guide_colourbar(
            barwidth = grid::unit(12, "lines"))) +
        theme_dynhr() +
        ggplot2::theme(
          panel.grid      = ggplot2::element_blank(),
          axis.line       = ggplot2::element_blank(),
          axis.line.x     = ggplot2::element_blank(),
          axis.line.y     = ggplot2::element_blank(),
          axis.ticks      = ggplot2::element_blank(),
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(
          title = "D10: Variance decomposition (heatmap)",
          subtitle = subtitle,
          x = "Shock", y = "Variable", fill = "Share of variance"
        )
      plots$vardec_table <- .apply_meta(p_tbl, meta)
    }

    .make_result(
      result  = list(vd_long = vd_long, foreign_check = foreign_check,
                     row_sums = row_sums, no_decomp_vars = no_decomp_vars,
                     correlated_shocks = correlated, cov_gap = cov_gap),
      pass    = pass,
      plots   = plots,
      summary = paste(summary_lines, collapse = "\n"),
      llm_summary = {
        badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
        top_line <- NULL
        if (!is.null(vd_mat) && length(vd_mat) > 0L) {
          top_df <- which(vd_mat == max(vd_mat), arr.ind = TRUE)
          top_line <- sprintf("  largest_contributor: %s->%s=%.2f",
                              rownames(vd_mat)[top_df[1, 1]],
                              colnames(vd_mat)[top_df[1, 2]], max(vd_mat))
        }
        action <- if (is.null(foreign_share))
          "No foreign shocks / output variable identified; informational only."
        else if (isTRUE(pass))
          "Foreign shock share within plausible range [15%, 60%] (Kamber et al. 2016)."
        else sprintf("Foreign shock share %.1f%% outside [15%%, 60%%]. %s",
                     foreign_share * 100,
                     if (foreign_share < gate[1])
                       "Too little foreign influence -- check open-economy transmission channels."
                     else
                       "Foreign shocks dominate output -- consider whether domestic shocks are properly identified.")
        paste(c(
          sprintf("D10 | Variance Decomposition | %s", badge),
          sprintf("  variables=%d shocks=%d", length(all_vars), length(all_shocks)),
          if (!is.null(foreign_share))
            sprintf("  foreign_share_%s=%.2f (benchmark [0.15, 0.60], Kamber et al. 2016)",
                    output_var, foreign_share),
          top_line,
          if (length(bad_vars) > 0)
            sprintf("  shares_not_summing_to_1: %s", paste(bad_vars, collapse = ", ")),
          if (correlated) "  caveat: correlated shocks, covariance terms unassigned",
          sprintf("  action: %s", action)
        ), collapse = "\n")
      }
    )
}
