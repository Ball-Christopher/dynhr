## R/diag-mcmc-d5-convergence.R
## --------------------------------------------------------------------------
## D5 MCMC convergence diagnostics (enhanced Phase 4/5)
##
## Implements bayesplot-equivalent diagnostics:
##   - Rank-normalised R-hat (Vehtari et al. 2021)
##   - Bulk-ESS and Tail-ESS (separately reported)
##   - Trace plots, rank overlay, ACF bar, density overlay by chain
##   - Posterior intervals (forest plot)
##   - NUTS-specific: treedepth histogram, BFMI, divergence count
##   - LLM-friendly summary with actionable advice
## --------------------------------------------------------------------------

#' D5. MCMC convergence diagnostics
#'
#' Computes rank-normalised R-hat, Bulk-ESS, and Tail-ESS per parameter,
#' generates trace plots, rank overlay histograms, ACF bar charts, per-chain
#' density overlays, and a posterior intervals (forest) plot.  If NUTS
#' metadata is present (treedepths, divergences, energy), also produces
#' HMC-specific diagnostics (BFMI, treedepth histogram).
#'
#' @param draws          Matrix (n_draws x n_params) from one chain, or a
#'   \code{dynhr_chains} object.
#' @param param_names    Character vector (optional; inferred from colnames).
#' @param chains_list    List of draw matrices for multi-chain diagnostics.
#'   Ignored when \code{draws} is a \code{dynhr_chains} object.
#' @param nuts_meta      List with NUTS diagnostics: \code{treedepths},
#'   \code{divergences}, \code{n_divergent}, \code{step_size}, and optionally
#'   \code{energy_trace}.  Auto-detected from a \code{dynhr_chains} object.
#' @section Badge (0.9.4):
#' Three-tier, following Vehtari et al. (2021):
#' \describe{
#'   \item{FAIL}{any \eqn{\hat R \ge 1.05}, or \eqn{\hat R} not computable
#'     (non-finite or constant draws -- a stuck sampler must never be reported
#'     as converged). 1.05 is a practical convention, NOT a Vehtari et al.
#'     number: their paper defines the single cutoff 1.01. It is the value the
#'     R-hat plot has always drawn as its "not converged" line.}
#'   \item{WARN}{any \eqn{1.01 < \hat R < 1.05} ("keep sampling / investigate"),
#'     or any bulk- or tail-ESS below \code{ess_target} /
#'     \code{ess_tail_target}.}
#'   \item{PASS}{otherwise.}
#' }
#' There is no ESS FAIL tier: a low ESS means the Monte Carlo error is larger
#' than you want, not that the draws are from the wrong distribution.
#'
#' @param ess_target     Minimum Bulk-ESS threshold. \code{NULL} (default) uses
#'   Vehtari et al. (2021)'s own scaling rule, \code{100 * n_chains} (400 for
#'   the standard four chains). The pre-0.9.4 flat default of 1000 silently
#'   assumed ten chains and appears nowhere in the literature.
#' @param ess_tail_target Minimum Tail-ESS threshold. \code{NULL} (default) is
#'   again \code{100 * n_chains}; Vehtari et al. apply the same rule to bulk-
#'   and tail-ESS.
#' @param acf_max_params Max parameters in ACF plot, sorted by worst lag-1
#'   (default 12).
#' @return dynhr_diagnostic with \code{$plots}, \code{$summary}, and
#'   \code{$llm_summary}.
#' @references
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Burkner, P.-C.
#'   (2021). Rank-normalization, folding, and localization: an improved
#'   \eqn{\hat R} for assessing convergence of MCMC. \emph{Bayesian Analysis},
#'   16(2), 667-718.
#' @noRd
d5_mcmc_convergence <- function(draws,
                                param_names      = NULL,
                                chains_list      = NULL,
                                nuts_meta        = NULL,
                                ess_target       = NULL,
                                ess_tail_target  = NULL,
                                acf_max_params   = 12L,
                                meta             = NULL) {

    # --- Unpack dynhr_chains -----------------------------------------------
    if (inherits(draws, "dynhr_chains")) {
      ch_obj <- draws
      draws  <- ch_obj$chain
      if (is.null(nuts_meta) && !is.null(ch_obj$treedepths)) {
        nuts_meta <- list(
          treedepths   = ch_obj$treedepths,
          divergences  = ch_obj$divergences %||% rep(FALSE, nrow(draws)),
          n_divergent  = ch_obj$n_divergent %||% 0L,
          step_size    = ch_obj$step_size,
          energy_trace = ch_obj$energy_trace
        )
      }
    }

    draws   <- as.matrix(draws)
    n_draws <- nrow(draws)
    n_par   <- ncol(draws)

    if (is.null(param_names))
      param_names <- colnames(draws) %||% paste0("theta_", seq_len(n_par))
    colnames(draws) <- param_names

    # --- Chains list -------------------------------------------------------
    if (is.null(chains_list)) {
      chains_list_full <- list(draws)
    } else {
      chains_list_full <- lapply(chains_list, function(m) {
        m2 <- as.matrix(m)
        colnames(m2) <- param_names
        m2
      })
    }
    n_chains <- length(chains_list_full)
    is_multi <- n_chains >= 2L

    ## Vehtari et al. (2021): bulk- and tail-ESS should both exceed
    ## 100 x (number of chains). A single chain is split in half here, so the
    ## effective chain count for the rule is max(n_chains, 2) -- which also
    ## keeps the one-chain default at the paper's familiar 200-per-split floor
    ## rather than a degenerate 100.
    if (is.null(ess_target))      ess_target      <- 100 * max(n_chains, 2L)
    if (is.null(ess_tail_target)) ess_tail_target <- 100 * max(n_chains, 2L)

    dims <- vapply(chains_list_full, function(m) c(nrow(m), ncol(m)),
                   numeric(2))
    if (any(dims[2, ] != n_par))
      .dynhr_abort("D5: every chain in `chains_list` must have ", n_par,
                   " columns (one per parameter).")
    if (length(unique(dims[1, ])) > 1L)
      .dynhr_abort("D5: chains in `chains_list` have unequal lengths (",
                   paste(dims[1, ], collapse = ", "),
                   "); R-hat and ESS need equal-length chains.")

    # --- Convergence statistics (Vehtari et al. 2021) ----------------------
    # A single chain is split in half, so split-R-hat is available too.
    conv <- .d5_convergence(chains_list_full, param_names)

    # Lag-1 ACF per parameter (first chain)
    acf1 <- vapply(seq_len(n_par), function(j) {
      y <- chains_list_full[[1]][, j]
      if (length(y) < 2L || !all(is.finite(y)) || stats::var(y) == 0)
        return(NA_real_)
      acf(y, lag.max = 1L, plot = FALSE)$acf[2, 1, 1]
    }, numeric(1))
    names(acf1) <- param_names

    # --- Badge (PASS / WARN / FAIL) ----------------------------------------
    # NA R-hat (non-finite or constant draws) is NOT computable and stays a
    # hard failure: a stuck sampler must not be reported as converged. A low
    # ESS is a Monte-Carlo-error statement, so it only ever warns.
    ess_bulk_fail <- is.na(conv$ess_bulk) | conv$ess_bulk < ess_target
    ess_tail_fail <- is.na(conv$ess_tail) | conv$ess_tail < ess_tail_target
    rhat_fail     <- is.na(conv$rhat) | conv$rhat >= 1.05
    rhat_warn     <- !is.na(conv$rhat) & conv$rhat > 1.01 & conv$rhat < 1.05
    pass <- !any(rhat_fail)
    warn <- pass && any(rhat_warn | ess_bulk_fail | ess_tail_fail)
    status <- if (!pass) "FAIL" else if (warn) "WARN" else "PASS"

    # --- Plots -------------------------------------------------------------
    plots <- list()
    USE_GG <- requireNamespace("ggplot2", quietly = TRUE)

    if (USE_GG) {
      col_pass <- dynhr_colours$green   %||% "#107C10"
      col_warn <- dynhr_colours$orange  %||% "#D18A00"
      col_fail <- dynhr_colours$red     %||% "#A80000"
      col_blue <- dynhr_colours$mid_blue %||% "#1B7CB6"

      # Derive a parameter-group column from the name prefix for faceting.
      .param_group <- function(nms) {
        prefixes <- c("sig", "sigma", "rho", "phi", "b")
        grp <- sub("_.*", "", nms)
        # normalise to a known group where possible
        grp <- ifelse(grepl("^sig", grp, ignore.case = TRUE), "sigma",
               ifelse(grepl("^rho", grp, ignore.case = TRUE), "rho",
               ifelse(grepl("^phi", grp, ignore.case = TRUE), "phi",
               ifelse(grepl("^b$",  grp, ignore.case = TRUE), "b", grp))))
        grp
      }

      # (a) ESS: Bulk + Tail bars -------------------------------------------
      ess_df <- rbind(
        data.frame(param = param_names, ess = conv$ess_bulk,
                   type = "Bulk", stringsAsFactors = FALSE),
        data.frame(param = param_names, ess = conv$ess_tail,
                   type = "Tail", stringsAsFactors = FALSE)
      )
      ess_df$param  <- factor(ess_df$param, levels = param_names)
      ess_df$group  <- .param_group(as.character(ess_df$param))
      ess_df$status <- ifelse(
        ess_df$type == "Bulk",
        ifelse(ess_df$ess >= ess_target,      "Adequate", "Low"),
        ifelse(ess_df$ess >= ess_tail_target, "Adequate", "Low")
      )
      # Not computable (non-finite / constant draws): draw a zero-height
      # "Low" bar rather than silently dropping the parameter.
      ess_df$not_comp <- is.na(ess_df$ess)
      ess_df$status[ess_df$not_comp] <- "Low"
      ess_df$ess[ess_df$not_comp] <- 0
      tgt_df <- data.frame(type   = c("Bulk", "Tail"),
                           target = c(ess_target, ess_tail_target),
                           stringsAsFactors = FALSE)

      p_ess <- ggplot2::ggplot(
        ess_df,
        ggplot2::aes(x = param, y = ess,
                     fill = interaction(type, status))
      ) +
        ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
        ggplot2::geom_text(
          data = ess_df[ess_df$not_comp, , drop = FALSE],
          ggplot2::aes(x = param, y = 0), inherit.aes = FALSE,
          label = "not computable", vjust = -0.5, size = 3,
          colour = col_fail) +
        ggplot2::geom_hline(data = tgt_df,
                            ggplot2::aes(yintercept = target),
                            linetype = "dashed", colour = col_fail,
                            linewidth = 0.5) +
        ggplot2::scale_fill_manual(
          values = c("Bulk.Adequate" = col_blue, "Bulk.Low" = col_fail,
                     "Tail.Adequate" = col_pass, "Tail.Low" = col_warn),
          labels = c("Bulk.Adequate" = "Bulk OK",  "Bulk.Low" = "Bulk low",
                     "Tail.Adequate" = "Tail OK",  "Tail.Low" = "Tail low"),
          name = NULL
        ) +
        # facet_grid with space="free_x" prevents single-bar groups from getting
        # oversized panels; small groups share the column width proportionally.
        ggplot2::facet_grid(type ~ group, scales = "free_x", space = "free_x") +
        theme_dynhr_diagnostic() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          # Use a plain font for strip text to avoid "phi"->"phl" / "sig"->"slg"
          # ligature substitution that Source Sans 3 applies to these prefixes.
          strip.text  = ggplot2::element_text(family = "sans")
        ) +
        ggplot2::labs(
          title    = "D5: Effective sample size (Bulk and Tail)",
          subtitle = sprintf("Bulk target >= %g  |  Tail target >= %g",
                              ess_target, ess_tail_target),
          x = NULL, y = "ESS")
      attr(p_ess, "dynhr_fig_height") <- 7.0
      plots$ess <- .apply_meta(p_ess, meta)

      # (b) Rank-normalised R-hat -------------------------------------------
      if (!all(is.na(conv$rhat))) {
        rhat_status <- ifelse(is.na(conv$rhat), "Unknown",
                        ifelse(conv$rhat < 1.01, "Converged (< 1.01)",
                          ifelse(conv$rhat < 1.05, "Marginal (1.01-1.05)",
                                 "Not converged (>= 1.05)")))
        rhat_df <- data.frame(
          param  = factor(param_names, levels = param_names),
          rhat   = conv$rhat,
          status = rhat_status,
          group  = .param_group(param_names),
          stringsAsFactors = FALSE
        )
        rhat_df <- rhat_df[!is.na(rhat_df$rhat), , drop = FALSE]
        y_top <- max(1.06, max(rhat_df$rhat) + 0.005)
        p_rhat <- ggplot2::ggplot(
          rhat_df, ggplot2::aes(x = param, y = rhat, fill = status)
        ) +
          ggplot2::geom_segment(
            ggplot2::aes(xend = param, y = 1, yend = rhat),
            colour = "grey60", linewidth = 0.5) +
          ggplot2::geom_point(shape = 21, size = 3, colour = "grey20") +
          ggplot2::geom_hline(yintercept = 1, colour = "grey40",
                              linewidth = 0.3) +
          ggplot2::coord_cartesian(ylim = c(min(1, min(rhat_df$rhat)), y_top)) +
          ggplot2::geom_hline(yintercept = 1.05, linetype = "dashed",
                              colour = col_fail, linewidth = 0.5) +
          ggplot2::geom_hline(yintercept = 1.01, linetype = "dotted",
                              colour = col_warn, linewidth = 0.4) +
          ggplot2::scale_fill_manual(
            values = c("Converged (< 1.01)"      = col_pass,
                       "Marginal (1.01-1.05)"    = col_warn,
                       "Not converged (>= 1.05)" = col_fail,
                       "Unknown"                 = "grey70"),
            name = NULL
          ) +
          ggplot2::facet_wrap(~ group, scales = "free_x") +
          theme_dynhr_diagnostic() +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            strip.text  = ggplot2::element_text(family = "sans")
          ) +
          ggplot2::labs(
            title    = "D5: Rank-normalised R-hat (Vehtari et al. 2021)",
            subtitle = sprintf(
              "Max of bulk and folded split-R-hat, %d chain(s) split in half  |  Dashed = 1.05 (FAIL)  |  Dotted = 1.01 (Vehtari et al. 2021 cutoff; 1.01-1.05 = WARN)",
              n_chains),
            x = NULL, y = expression(hat(R)[rank]))
        plots$rhat <- .apply_meta(p_rhat, meta)
      }

      # (c) Trace plots -------------------------------------------------------
      n_show   <- min(n_par, 20L)
      thin_idx <- unique(round(seq(1L, n_draws,
                                    length.out = min(n_draws, 2000L))))
      chain_pal <- .chain_palette(n_chains)

      trace_df <- do.call(rbind, lapply(seq_len(n_chains), function(ci) {
        m      <- chains_list_full[[ci]]
        # Clamp thin_idx to this chain's row count — prevents subscript OOB
        # when chains are shorter than the pooled draw count.
        idx_m  <- thin_idx[thin_idx <= nrow(m)]
        if (length(idx_m) == 0L) idx_m <- seq_len(min(nrow(m), 10L))
        thin <- m[idx_m, seq_len(n_show), drop = FALSE]
        data.frame(
          iter  = rep(idx_m, n_show),
          value = as.vector(thin),
          param = rep(param_names[seq_len(n_show)], each = length(idx_m)),
          chain = sprintf("Chain %d", ci),
          stringsAsFactors = FALSE
        )
      }))
      trace_df$param <- factor(trace_df$param,
                                levels = param_names[seq_len(n_show)])

      # Tick formatter: 50,000 -> "50k", 1,500,000 -> "1.5M".
      .k_label <- function(x) {
        ifelse(abs(x) >= 1e6,
               sprintf("%.1fM", x / 1e6),
               ifelse(abs(x) >= 1e3,
                      sprintf("%gk", x / 1e3),
                      as.character(x)))
      }
      p_trace <- ggplot2::ggplot(
        trace_df,
        ggplot2::aes(x = iter, y = value, colour = chain, group = chain)
      ) +
        ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
        ggplot2::facet_wrap(~ param, scales = "free_y", ncol = 4L) +
        ggplot2::scale_colour_manual(values = chain_pal, name = NULL) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          override.aes = list(linewidth = 1.2, alpha = 1))) +
        ggplot2::scale_x_continuous(n.breaks = 3L, labels = .k_label) +
        ggplot2::scale_y_continuous(n.breaks = 3L) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(strip.text = ggplot2::element_text(
          size = ggplot2::rel(0.85))) +
        ggplot2::labs(
          title    = "D5: MCMC trace plots",
          subtitle = sprintf("%d chain(s)  |  %s draws",
                              n_chains, format(n_draws, big.mark = ",")),
          x = "Iteration", y = "Value")
      # Hint to the report templates: tall facet grid needs more height.
      p_trace <- .apply_meta(p_trace, meta)
      attr(p_trace, "dynhr_fig_height") <- 9.5
      plots$trace <- p_trace

      # (d) Rank overlay (multi-chain only) -----------------------------------
      if (is_multi) {
        ro_df <- .rank_overlay_df(chains_list_full, param_names,
                                   max_params = min(n_par, 16L))
        if (!is.null(ro_df) && nrow(ro_df) > 0) {
          ro_df$param <- factor(ro_df$param,
                                levels = intersect(param_names, ro_df$param))
          p_rank <- ggplot2::ggplot(
            ro_df,
            ggplot2::aes(x = bin_mid / (n_chains * nrow(chains_list_full[[1]])),
                         y = freq,
                         colour = chain, group = chain)
          ) +
            ggplot2::geom_line(linewidth = 0.5) +
            ggplot2::geom_point(size = 1) +
            ggplot2::geom_hline(
              yintercept = unique(ro_df$expected),
              linetype = "dashed", colour = "grey50", linewidth = 0.3
            ) +
            ggplot2::facet_wrap(~ param, scales = "free_y", ncol = 4L) +
            ggplot2::scale_x_continuous(limits = c(0, 1),
                                        breaks = c(0, 0.5, 1)) +
            ggplot2::scale_colour_manual(
              values = .chain_palette(n_chains), name = NULL) +
            theme_dynhr_diagnostic() +
            ggplot2::theme(strip.text = ggplot2::element_text(
              size = ggplot2::rel(0.75))) +
            ggplot2::labs(
              title = "D5: Rank overlay (Vehtari 2021)",
              subtitle = sprintf(paste(
                "Draws per rank bin (10 bins); each chain should scatter",
                "around %s (dashed) if chains mix"),
                format(unique(ro_df$expected), big.mark = ",")),
              x = "Rank among all pooled draws (fraction)", y = "Draws in bin")
          plots$rank_overlay <- .apply_meta(p_rank, meta)
        }
      }

      # (e) ACF bar chart (top N worst lag-1) ---------------------------------
      acf_order  <- names(sort(abs(acf1), decreasing = TRUE))
      acf_params <- acf_order[seq_len(min(length(acf_order), acf_max_params))]
      n_chain1   <- nrow(chains_list_full[[1]])
      max_lag    <- min(20L, floor(n_chain1 / 5L))
      if (length(acf_params) > 0L && max_lag >= 1L) {

      acf_df <- do.call(rbind, lapply(acf_params, function(pn) {
        a <- acf(chains_list_full[[1]][, pn], lag.max = max_lag,
                 plot = FALSE)$acf[-1, 1, 1]
        data.frame(param = pn, lag = seq_along(a), acf = a,
                   stringsAsFactors = FALSE)
      }))
      # Worst-mixing parameter at the TOP of the table (tiles plot bottom-up).
      acf_df$param <- factor(acf_df$param, levels = rev(acf_params))

      p_acf <- ggplot2::ggplot(
        acf_df, ggplot2::aes(x = factor(lag), y = param, fill = acf)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::geom_text(
          ggplot2::aes(
            label  = sub("^(-?)0\\.", "\\1.", sprintf("%.2f", acf)),
            colour = abs(acf) > 0.5
          ),
          size = 2.7
        ) +
        ggplot2::scale_colour_manual(
          values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
        scale_fill_dynhr_sunset(limits = c(-1, 1), name = "ACF") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.line  = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank()) +
        ggplot2::labs(
          title    = sprintf("D5: Autocorrelation by lag -- top %d by worst lag-1",
                              length(acf_params)),
          subtitle = sprintf(
            "Chain 1 (%s draws)  |  |ACF| > %.3f (= 2/sqrt(n)) at many lags indicates slow mixing",
            format(n_chain1, big.mark = ","), 2 / sqrt(n_chain1)),
          x = "Lag", y = NULL)
      p_acf <- .apply_meta(p_acf, meta)
      attr(p_acf, "dynhr_fig_height") <- max(4, 0.42 * length(acf_params) + 2)
      plots$acf <- p_acf
      }

      # (f) Density overlay by chain (multi-chain) ----------------------------
      if (is_multi) {
        n_dens  <- min(n_par, 16L)
        dens_df <- do.call(rbind, lapply(seq_len(n_chains), function(ci) {
          m <- chains_list_full[[ci]]
          do.call(rbind, lapply(seq_len(n_dens), function(j) {
            d <- density(m[, j], n = 256L)
            data.frame(x = d$x, y = d$y,
                       param = param_names[j],
                       chain = sprintf("Chain %d", ci),
                       stringsAsFactors = FALSE)
          }))
        }))
        dens_df$param <- factor(dens_df$param,
                                 levels = param_names[seq_len(n_dens)])

        p_dens <- ggplot2::ggplot(
          dens_df,
          ggplot2::aes(x = x, y = y,
                       colour = chain, group = chain)
        ) +
          ggplot2::geom_line(linewidth = 0.5, alpha = 0.8) +
          ggplot2::scale_colour_manual(
            values = .chain_palette(n_chains), name = NULL) +
          ggplot2::facet_wrap(~ param, scales = "free", ncol = 4L) +
          ggplot2::scale_x_continuous(n.breaks = 3L) +
          theme_dynhr_diagnostic() +
          ggplot2::theme(strip.text = ggplot2::element_text(
            size = ggplot2::rel(0.75))) +
          ggplot2::labs(
            title    = "D5: Per-chain posterior densities",
            subtitle = "Overlapping chains = convergence; separation = non-convergence",
            x = "Value", y = "Density")
        plots$density_overlay <- .apply_meta(p_dens, meta)
      }

      # (g) Posterior intervals (forest plot) -- faceted by param group --------
      pool   <- do.call(rbind, chains_list_full)
      int_df <- do.call(rbind, lapply(param_names, function(pn) {
        x <- pool[, pn]
        data.frame(param = pn,
                   med = median(x),
                   q05 = quantile(x, 0.05), q25 = quantile(x, 0.25),
                   q75 = quantile(x, 0.75), q95 = quantile(x, 0.95),
                   group = .param_group(pn),
                   stringsAsFactors = FALSE)
      }))
      int_df$param <- factor(int_df$param, levels = rev(param_names))

      p_int <- ggplot2::ggplot(int_df,
                                          ggplot2::aes(y = param)) +
        ggplot2::geom_segment(
          ggplot2::aes(x = q05, xend = q95, yend = param),
          colour = col_blue, linewidth = 1.0, alpha = 0.35) +
        ggplot2::geom_segment(
          ggplot2::aes(x = q25, xend = q75, yend = param),
          colour = col_blue, linewidth = 2.8, alpha = 0.55) +
        ggplot2::geom_point(ggplot2::aes(x = med),
                            colour = col_blue, size = 2.0) +
        # One panel per parameter group with its own x range, so a
        # tightly-estimated parameter is not squashed to a dot by a diffuse one.
        ## One column stacks every group; at 68 parameters that is a 74 in
        ## figure that no page can hold and a capped one collapses to nothing
        ## legible (NZSIM-scale check).  Widen instead: ~20 facets per column.
        ggplot2::facet_wrap(~ group, scales = "free",
                            ncol = max(1L, ceiling(length(unique(int_df$group)) / 20))) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(
          # Suppress "phi"->"phl" / "sig"->"slg" font ligature in strip labels
          strip.text = ggplot2::element_text(family = "sans")
        ) +
        ggplot2::labs(
          title    = "D5: Posterior intervals (all parameters)",
          subtitle = sprintf("Pooled %d chain(s)  |  point = median, thick = 50%%, thin = 90%% interval  |  x-axis free per group",
                             n_chains),
          x = "Parameter value", y = NULL)
      p_int <- .apply_meta(p_int, meta)
      ## Capped: at 68 parameters the uncapped formula asked for 74 in, which
      ## ggsave() refuses (> 50 in) and no page or screen can show anyway.
      n_grp_int <- length(unique(int_df$group))
      ncol_int  <- max(1L, ceiling(n_grp_int / 20))
      attr(p_int, "dynhr_fig_height") <-
        min(30, max(4, (0.35 * n_par + 0.9 * n_grp_int) / ncol_int + 1.5))
      plots$intervals <- p_int

      # (h) NUTS-specific plots -----------------------------------------------
      if (!is.null(nuts_meta))
        plots <- c(plots,
                   .d5_nuts_plots(nuts_meta, col_fail, col_warn, col_pass,
                                  diag_meta = meta))
    }

    # --- LLM summary -------------------------------------------------------
    llm_summary <- .d5_llm(conv, acf1, n_draws, n_chains, n_par,
                             ess_bulk_fail, ess_tail_fail, rhat_fail,
                             rhat_warn, ess_target, ess_tail_target, status,
                             nuts_meta, param_names)

    .make_result(
      result = list(convergence  = conv,
                    acf1         = acf1,
                    n_chains     = n_chains,
                    n_draws      = n_draws,
                    nuts_meta    = nuts_meta),
      pass        = pass,
      warn        = warn,
      plots       = plots,
      summary     = .d5_console(conv, acf1, n_draws, n_chains,
                                 n_par, status, nuts_meta),
      llm_summary = llm_summary
    )
}


#' D21. Bayesian KPS posterior precision updating indicator
#'
#' Koop, Pesaran & Smith (2013, JBES) indicator. For an identified parameter
#' the posterior precision grows linearly in the sample size T, so
#' \code{precision / T} converges to a nonzero constant and the elasticity
#' \eqn{d \log(\mathrm{precision}) / d \log T} tends to 1; for an unidentified
#' parameter precision converges to a constant, \code{precision / T -> 0} and
#' the elasticity tends to 0.
#'
#' Precision is the inverse of each parameter's MARGINAL posterior variance.
#' Because posterior precision is (approximately) prior precision + T x
#' information, the log-log slope over the whole grid is pulled below 1 by the
#' prior at small T even for identified parameters. The badge therefore uses
#' the TERMINAL elasticity, between the two largest sample sizes, which is the
#' closest available estimate of the asymptotic rate. The full-range OLS slope
#' and the KPS ratio \code{precision / T} are reported alongside.
#'
#' Draws are supplied either as
#' - \code{draws_by_T}: list of posterior draw matrices (n_draws x n_params),
#'   one per sample size, or
#' - \code{kps_runner_fn}: function(T) returning such a draw matrix.
#' Columns are matched by NAME when the matrices carry column names, so the
#' column order may differ between sample sizes.
#'
#' @param draws_by_T List of draw matrices (n_draws x n_params). If
#'   \code{sample_sizes} is NULL, the names must be the sample sizes.
#' @param sample_sizes Numeric vector of (>= 3 distinct, positive) sample sizes
#'   aligned to the list order.
#' @param kps_runner_fn Optional function(T) -> draw matrix.
#' @param param_names Optional parameter names. If the draw matrices have
#'   column names, these select (and order) the columns by name; otherwise
#'   they label the columns positionally.
#' @param slope_threshold Minimum terminal elasticity considered adequately
#'   identified (1 = rate-T updating, 0 = no updating). \strong{The default
#'   0.80 is a package choice with no literature source}: Koop, Pesaran &
#'   Smith (2013) propose the precision-growth concept and its two limiting
#'   cases (elasticity -> 1 identified, -> 0 unidentified) but name no numeric
#'   cutoff. Move it if your grid of sample sizes is short or widely spaced.
#' @param meta Optional plot metadata.
#' @return dynhr_diagnostic object; \code{result} holds \code{sample_sizes}
#'   (sorted), \code{precision} (n_params x n_T), \code{precision_ratio}
#'   (precision / T), \code{elasticity} (terminal, used for the badge),
#'   \code{slope} (full-range OLS log-log slope) and \code{weak_params}.
#' @section Why the terminal elasticity:
#' KPS's claim is about the RATE at which precision grows, so the statistic
#' must be scale-free in the units of precision -- which a log-log elasticity
#' is and a levels fit \eqn{\mathrm{precision} = p_0 + cT} is not. The
#' full-range OLS log-log slope is mechanically pulled below 1 by the prior's
#' fixed contribution at small T, so the terminal (two-largest-T) elasticity is
#' the least biased available estimate of the asymptotic rate; the full-range
#' slope is reported alongside for reference only.
#' @references
#'   Koop, G., Pesaran, M. H., & Smith, R. P. (2013). On identification of
#'   Bayesian DSGE models. \emph{Journal of Business & Economic Statistics},
#'   31(3), 300-314.
#' @noRd
d21_kps_precision_update <- function(draws_by_T      = NULL,
                                     sample_sizes    = NULL,
                                     kps_runner_fn   = NULL,
                                     param_names     = NULL,
                                     slope_threshold = 0.80,
                                     meta            = NULL) {

    if (is.null(draws_by_T) && is.null(kps_runner_fn)) {
      return(.make_result(
        pass = NA,
        summary = paste(
          "D21 KPS posterior precision update: SKIPPED -- provide draws_by_T",
          "or kps_runner_fn(T) to evaluate precision growth with sample size."
        )
      ))
    }

    if (is.null(draws_by_T)) {
      if (!is.function(kps_runner_fn)) {
        .dynhr_abort("D21: `kps_runner_fn` must be a function(T).")
      }
      if (is.null(sample_sizes) || length(sample_sizes) < 3L) {
        .dynhr_abort("D21: when using kps_runner_fn, provide at least 3 sample_sizes.")
      }
      draws_by_T <- lapply(sample_sizes, kps_runner_fn)
      names(draws_by_T) <- as.character(sample_sizes)
    }
    if (is.matrix(draws_by_T) || is.data.frame(draws_by_T) || !is.list(draws_by_T)) {
      .dynhr_abort("D21: `draws_by_T` must be a list of draw matrices, one per sample size.")
    }

    nm <- names(draws_by_T)
    nm_num <- if (!is.null(nm) && all(grepl("^[0-9]+(\\.[0-9]+)?$", nm)))
      as.numeric(nm) else NULL
    if (is.null(sample_sizes)) {
      if (is.null(nm_num)) {
        .dynhr_abort("D21: sample_sizes not supplied and could not be inferred ",
                     "from names(draws_by_T).")
      }
      sample_sizes <- nm_num
    }
    sample_sizes <- as.numeric(sample_sizes)
    if (length(draws_by_T) != length(sample_sizes)) {
      .dynhr_abort("D21: draws_by_T (", length(draws_by_T), ") and sample_sizes (",
                   length(sample_sizes), ") must have the same length.")
    }
    if (!is.null(nm_num) && any(nm_num != sample_sizes)) {
      .dynhr_abort("D21: names(draws_by_T) (", paste(nm, collapse = ", "),
                   ") disagree with sample_sizes (",
                   paste(sample_sizes, collapse = ", "), ").")
    }
    if (length(sample_sizes) < 3L || any(!is.finite(sample_sizes)) ||
        any(sample_sizes <= 0) || anyDuplicated(sample_sizes)) {
      .dynhr_abort("D21: need at least 3 distinct, positive, finite sample sizes; got ",
                   paste(sample_sizes, collapse = ", "), ".")
    }

    ord <- order(sample_sizes)
    sample_sizes <- sample_sizes[ord]
    mats <- lapply(draws_by_T[ord], as.matrix)

    # Resolve parameter names, then align every matrix BY NAME when possible.
    cn1 <- colnames(mats[[1]])
    if (is.null(param_names)) {
      param_names <- cn1 %||% paste0("theta_", seq_len(ncol(mats[[1]])))
    }
    param_names <- as.character(param_names)
    if (anyDuplicated(param_names)) {
      .dynhr_abort("D21: duplicated parameter names: ",
                   paste(unique(param_names[duplicated(param_names)]), collapse = ", "))
    }
    n_par <- length(param_names)
    for (k in seq_along(mats)) {
      m <- mats[[k]]
      if (!is.numeric(m) || nrow(m) < 2L) {
        .dynhr_abort("D21: draws at T = ", sample_sizes[k],
                     " must be a numeric matrix with at least 2 draws.")
      }
      cn <- colnames(m)
      if (!is.null(cn)) {
        miss <- setdiff(param_names, cn)
        if (length(miss) > 0L) {
          .dynhr_abort("D21: draws at T = ", sample_sizes[k],
                       " lack column(s): ", paste(miss, collapse = ", "))
        }
        m <- m[, param_names, drop = FALSE]
      } else if (ncol(m) != n_par) {
        .dynhr_abort("D21: draws at T = ", sample_sizes[k], " have ", ncol(m),
                     " unnamed columns; expected ", n_par, ".")
      } else {
        colnames(m) <- param_names
      }
      if (any(!is.finite(m))) {
        .dynhr_abort("D21: draws at T = ", sample_sizes[k],
                     " contain non-finite values.")
      }
      mats[[k]] <- m
    }

    post_var <- vapply(mats, function(m) apply(m, 2L, stats::var),
                       numeric(n_par))
    post_var <- matrix(post_var, nrow = n_par,
                       dimnames = list(param_names, as.character(sample_sizes)))
    if (any(post_var <= 0)) {
      bad <- which(post_var <= 0, arr.ind = TRUE)
      .dynhr_abort("D21: zero posterior variance for ",
                   paste(unique(param_names[bad[, 1]]), collapse = ", "),
                   " -- precision undefined (fixed parameter?).")
    }
    precision <- 1 / post_var
    precision_ratio <- sweep(precision, 2L, sample_sizes, "/")

    logT <- log(sample_sizes)
    K <- length(sample_sizes)
    slope <- vapply(seq_len(n_par), function(i) {
      unname(stats::coef(stats::lm(log(precision[i, ]) ~ logT))[2])
    }, numeric(1))
    elasticity <- (log(precision[, K]) - log(precision[, K - 1L])) /
      (logT[K] - logT[K - 1L])
    names(slope) <- names(elasticity) <- param_names

    weak <- param_names[elasticity < slope_threshold]
    pass <- length(weak) == 0L

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      lab <- sprintf("%s (elasticity %.2f)", param_names, elasticity)
      names(lab) <- param_names
      n_col <- if (n_par <= 3L) n_par else min(4L, ceiling(sqrt(n_par)))
      status <- ifelse(elasticity < slope_threshold, "weak updating", "rate-T updating")
      growth_df <- data.frame(
        T = rep(sample_sizes, each = n_par),
        parameter = rep(param_names, times = K),
        precision = as.vector(precision),
        status = rep(status, times = K),
        stringsAsFactors = FALSE
      )
      ref_df <- do.call(rbind, lapply(c(1, 0), function(b) data.frame(
        T = rep(sample_sizes, each = n_par),
        parameter = rep(param_names, times = K),
        precision = as.vector(outer(precision[, K], (sample_sizes / sample_sizes[K])^b)),
        reference = if (b == 1) "slope 1 (identified)" else "slope 0 (no updating)",
        stringsAsFactors = FALSE
      )))
      growth_df$parameter <- factor(lab[growth_df$parameter], levels = lab)
      ref_df$parameter <- factor(lab[ref_df$parameter], levels = lab)

      p_kps <- ggplot2::ggplot(growth_df, ggplot2::aes(x = T, y = precision)) +
        ggplot2::geom_line(data = ref_df,
                           ggplot2::aes(linetype = reference),
                           colour = dynhr_colours$grey, linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(colour = status), linewidth = 0.7) +
        ggplot2::geom_point(ggplot2::aes(colour = status), size = 1.8) +
        ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = n_col) +
        ggplot2::scale_x_log10(breaks = sample_sizes,
                               labels = format(sample_sizes, trim = TRUE,
                                               scientific = FALSE)) +
        ggplot2::scale_y_log10() +
        ggplot2::scale_colour_manual(
          values = c("rate-T updating" = dynhr_colours$mid_blue,
                     "weak updating"   = dynhr_colours$red),
          name = NULL, drop = FALSE) +
        ggplot2::scale_linetype_manual(
          values = c("slope 1 (identified)" = "dashed",
                     "slope 0 (no updating)" = "dotted"),
          name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(legend.position = "bottom",
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
          title = "D21: KPS posterior precision updating",
          subtitle = sprintf(paste0(
            "Precision 1/Var vs sample size; reference lines pass through the largest-T point.\n",
            "Badge: elasticity between the two largest T >= %.2f"), slope_threshold),
          x = "Sample size T (log scale)",
          y = "Posterior precision 1/Var (log scale)"
        )
      attr(p_kps, "dynhr_fig_height") <- max(4.5, 2.4 * ceiling(n_par / n_col) + 2)
      plots$kps_precision_growth <- .apply_meta(p_kps, meta)
    }

    fmt_list <- function(x) paste(sprintf("%s=%.2f", names(x), x), collapse = ", ")
    .make_result(
      result = list(
        sample_sizes    = sample_sizes,
        precision       = precision,
        precision_ratio = precision_ratio,
        elasticity      = elasticity,
        slope           = slope,
        slope_threshold = slope_threshold,
        weak_params     = weak
      ),
      pass = pass,
      plots = plots,
      summary = sprintf(paste0(
        "D21 KPS precision update: %d parameters over %d sample sizes (T = %s). ",
        "Terminal elasticity range [%.2f, %.2f] (threshold %.2f). %s"),
        n_par, K, paste(sample_sizes, collapse = ", "),
        min(elasticity), max(elasticity), slope_threshold,
        if (pass) "PASS -- all parameters update at close to rate T."
        else sprintf("FAIL -- weak precision updating for: %s", paste(weak, collapse = ", "))
      ),
      llm_summary = paste(c(
        sprintf("D21 | KPS Precision Updating | %s", if (pass) "PASS" else "FAIL"),
        sprintf("  params=%d sample_sizes=%s threshold(terminal elasticity)=%.2f",
                n_par, paste(sample_sizes, collapse = ","), slope_threshold),
        sprintf("  terminal_elasticity: %s", fmt_list(elasticity)),
        sprintf("  full_range_slope: %s", fmt_list(slope)),
        sprintf("  precision/T at T=%s: %s", sample_sizes[K],
                fmt_list(precision_ratio[, K])),
        if (length(weak) > 0) sprintf("  weak_updating: %s", paste(head(weak, 6), collapse = ", ")),
        sprintf("  action: %s",
                if (pass) "Posterior precision grows at close to rate T for every parameter."
                else "Precision of the listed parameters is not growing with T (precision/T -> 0 suggests non-identification); review local/spectral identification and prior dominance.")
      ), collapse = "\n")
    )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# The Vehtari et al. (2021) split-R-hat / ESS estimators D5 reports used to
# live here as `.d5_*`. 0.9.4 (ledger A5) MOVED them to R/diag-helpers.R and
# DELETED the wrong shared implementations they were written to replace, so
# there is now exactly one ESS/R-hat implementation in the package. The names
# are unchanged (`.d5_split`, `.d5_zscale`, `.d5_degenerate`, `.d5_rhat_basic`,
# `.d5_autocov`, `.d5_ess_basic`, `.d5_convergence`).

#' Rank overlay data.frame for ggplot
#' @noRd
.rank_overlay_df <- function(chains_list, param_names, max_params = 16L,
                              n_bins = 10L) {
  n_chains <- length(chains_list)
  n        <- nrow(chains_list[[1]])
  n_par    <- min(ncol(chains_list[[1]]), max_params)
  total_N  <- n * n_chains
  expected <- n / n_bins

  do.call(rbind, lapply(seq_len(n_par), function(j) {
    all_vals <- unlist(lapply(chains_list, function(m) m[, j]))
    pooled_r <- rank(all_vals, ties.method = "first")
    breaks   <- seq(0, total_N, length.out = n_bins + 1L)
    do.call(rbind, lapply(seq_len(n_chains), function(ci) {
      chain_r <- pooled_r[((ci - 1L) * n + 1L):(ci * n)]
      h       <- cut(chain_r, breaks = breaks, include.lowest = TRUE)
      freq    <- tabulate(h, nbins = n_bins)
      bin_mid <- (breaks[-1] + breaks[-(n_bins + 1L)]) / 2
      data.frame(param    = param_names[j],
                 bin_mid  = bin_mid,
                 freq     = freq,
                 chain    = sprintf("Chain %d", ci),
                 expected = expected,
                 stringsAsFactors = FALSE)
    }))
  }))
}

#' Chain colour palette (up to 8 chains)
#' @noRd
.chain_palette <- function(n) {
  base <- c(dynhr_palette_vibrant, unname(tol_vibrant["grey"]), "#000000")
  rep_len(base, max(1L, n))
}

#' NUTS-specific plots
#' @noRd
.d5_nuts_plots <- function(nuts_meta, col_fail, col_warn, col_pass,
                             diag_meta = NULL) {
  plots <- list()
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(plots)

  if (!is.null(nuts_meta$treedepths)) {
    td_df <- data.frame(depth = as.integer(nuts_meta$treedepths))
    max_td <- max(td_df$depth, na.rm = TRUE)
    p_td <- ggplot2::ggplot(
      td_df, ggplot2::aes(x = depth)
    ) +
      ggplot2::geom_bar(fill = dynhr_colours$mid_blue %||% "#1B7CB6",
                        width = 0.7) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title    = "D5 (NUTS): Tree depth per iteration",
        subtitle = sprintf("Max reached = %d  |  Mean = %.1f  |  High max-treedepth = too-small step size (NUTS hitting max depth)",
                           max_td, mean(td_df$depth, na.rm = TRUE)),
        x = "Tree depth", y = "Count")
    plots$nuts_treedepth <- .apply_meta(p_td, diag_meta)
  }

  if (!is.null(nuts_meta$energy_trace) && length(nuts_meta$energy_trace) > 2L) {
    bfmi_v <- .bfmi(nuts_meta$energy_trace)
    if (!is.na(bfmi_v)) {
      bfmi_df <- data.frame(
        chain = "Chain 1",
        bfmi  = bfmi_v,
        status = if (bfmi_v < 0.3) "Low (< 0.3)" else "OK (>= 0.3)"
      )
      p_bfmi <- ggplot2::ggplot(
        bfmi_df,
        ggplot2::aes(x = chain, y = bfmi, fill = status)
      ) +
        ggplot2::geom_col(width = 0.35) +
        ggplot2::geom_hline(yintercept = 0.3, linetype = "dashed",
                            colour = col_fail, linewidth = 0.5) +
        ggplot2::scale_fill_manual(
          values = c("Low (< 0.3)" = col_fail,
                     "OK (>= 0.3)" = col_pass),
          name = NULL
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title    = "D5 (NUTS): BFMI -- Bayesian Fraction of Missing Information",
          subtitle = "BFMI < 0.3 = posterior tails not well explored; consider reparameterisation",
          x = NULL, y = "BFMI")
      plots$nuts_bfmi <- .apply_meta(p_bfmi, diag_meta)
    }
  }

  plots
}

#' Console-friendly one-line summary
#' @noRd
.d5_console <- function(conv, acf1, n_draws, n_chains, n_par, status, nuts_meta) {
  rhat_str  <- if (!all(is.na(conv$rhat)))
    sprintf(" Rank-normalised split R-hat max=%.3f.", max(conv$rhat, na.rm = TRUE))
  else " R-hat: not computable."
  bad <- conv$param[is.na(conv$rhat) | is.na(conv$ess_bulk) | is.na(conv$ess_tail)]
  bad_str <- if (length(bad))
    sprintf(" Not computable (non-finite or constant draws): %s.",
            paste(bad, collapse = ", "))
  else ""
  acf_ok <- acf1[is.finite(acf1)]
  acf_str <- if (length(acf_ok))
    sprintf(" Lag-1 ACF worst: %.2f (%s).", max(abs(acf_ok)),
            names(acf_ok)[which.max(abs(acf_ok))])
  else ""
  .min_or_na <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
  nuts_str  <- if (!is.null(nuts_meta) && (nuts_meta$n_divergent %||% 0L) > 0)
    sprintf(" NUTS: %d divergences.", nuts_meta$n_divergent)
  else ""
  sprintf(
"D5 [%s]: %d chain(s) x %s draws, %d params. Bulk-ESS min=%.0f median=%.0f. Tail-ESS min=%.0f.%s%s%s%s",
    status,
    n_chains, format(n_draws, big.mark = ","), n_par,
    .min_or_na(conv$ess_bulk),
    stats::median(conv$ess_bulk, na.rm = TRUE),
    .min_or_na(conv$ess_tail),
    rhat_str, bad_str, acf_str,
    nuts_str
  )
}

#' LLM-friendly D5 summary with actionable advice
#' @noRd
.d5_llm <- function(conv, acf1, n_draws, n_chains, n_par,
                     ess_bulk_fail, ess_tail_fail, rhat_fail,
                     rhat_warn, ess_target, ess_tail_target, status,
                     nuts_meta, param_names) {

  acf1 <- acf1[is.finite(acf1)]
  if (!length(acf1)) acf1 <- c(`(none)` = 0)
  .mn <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
  bad <- param_names[is.na(conv$rhat) | is.na(conv$ess_bulk) |
                     is.na(conv$ess_tail)]
  worst_bulk <- param_names[which.min(conv$ess_bulk)]
  worst_tail <- param_names[which.min(conv$ess_tail)]
  worst_rhat <- if (!all(is.na(conv$rhat))) param_names[which.max(conv$rhat)]
                else NA_character_
  worst_acf1 <- names(which.max(abs(acf1)))
  top5_acf   <- head(paste0(
    names(sort(acf1, decreasing = TRUE)),
    "=", sprintf("%.2f", sort(acf1, decreasing = TRUE))
  ), 5)

  # NUTS numbers
  nuts_lines <- character(0)
  if (!is.null(nuts_meta)) {
    n_div <- nuts_meta$n_divergent %||% 0L
    nuts_lines <- sprintf(
      "  nuts: divergences=%d (%.2f%%) step_size=%.3e mean_treedepth=%.1f",
      n_div, 100 * n_div / max(1L, n_draws),
      nuts_meta$step_size %||% NA_real_,
      mean(nuts_meta$treedepths %||% NA_real_, na.rm = TRUE)
    )
    if (!is.null(nuts_meta$energy_trace)) {
      bv <- .bfmi(nuts_meta$energy_trace)
      nuts_lines <- c(nuts_lines,
        sprintf("  bfmi=%.3f %s", bv,
                if (!is.na(bv) && bv < 0.3) "[WARN: <0.3 -- consider reparameterisation]"
                else "[OK]"))
    }
  }

  # Action advice
  issues <- character(0)
  if (length(bad))
    issues <- c(issues, sprintf(
      "%s: R-hat/ESS not computable (non-finite or constant draws) -- sampler stuck or parameter fixed; drop it or inspect the chain",
      paste(bad, collapse = ", ")))
  if (any(rhat_fail & !is.na(conv$rhat)))
    issues <- c(issues, sprintf(
      "%s: rhat=%.3f>=1.05 -- chains not converged; increase draws or check for multi-modality",
      worst_rhat, max(conv$rhat, na.rm = TRUE)
    ))
  else if (any(rhat_warn))
    issues <- c(issues, sprintf(
      "%s: rhat=%.3f>1.01 -- above the Vehtari et al. (2021) cutoff but below 1.05; keep sampling before using the draws",
      worst_rhat, max(conv$rhat, na.rm = TRUE)
    ))
  if (any(ess_bulk_fail & !is.na(conv$ess_bulk)))
    issues <- c(issues, sprintf(
      "%s: ess_bulk=%.0f<target_%g -- high autocorrelation; tune proposal or increase draws",
      worst_bulk, .mn(conv$ess_bulk), ess_target
    ))
  if (any(ess_tail_fail & !is.na(conv$ess_tail)) &&
      !any(ess_bulk_fail & !is.na(conv$ess_bulk)))
    issues <- c(issues, sprintf(
      "%s: ess_tail=%.0f<target_%g -- tails poorly sampled; check parameter boundaries",
      worst_tail, .mn(conv$ess_tail), ess_tail_target
    ))
  if (max(abs(acf1)) > 0.5 && length(issues) == 0)
    issues <- c(issues, sprintf(
      "%s: acf1=%.2f>0.5 -- tune proposal scale or reparameterise",
      worst_acf1, max(abs(acf1))
    ))
  if (!is.null(nuts_meta)) {
    n_div <- nuts_meta$n_divergent %||% 0L
    if (n_div > 0)
      issues <- c(issues, sprintf(
        "NUTS: %d divergences (%.1f%%) -- reparameterise or reduce step_size",
        n_div, 100 * n_div / max(1L, n_draws)
      ))
    bv <- if (!is.null(nuts_meta$energy_trace)) .bfmi(nuts_meta$energy_trace)
          else NA_real_
    if (!is.na(bv) && bv < 0.3)
      issues <- c(issues, sprintf(
        "NUTS: bfmi=%.3f<0.3 -- heavy tails; try non-centred parameterisation",
        bv
      ))
  }
  action <- if (length(issues) == 0)
    "All convergence metrics within acceptable range."
  else paste(issues, collapse = "; ")

  paste(c(
    sprintf("D5 | MCMC Convergence | %s", status),
    sprintf("  chains=%d draws=%s params=%d",
            n_chains, format(n_draws, big.mark = ","), n_par),
    sprintf("  ess_bulk: min=%.0f (%s) median=%.0f below_target_%g=%d/%d",
            .mn(conv$ess_bulk), worst_bulk,
            median(conv$ess_bulk, na.rm = TRUE),
            ess_target, sum(ess_bulk_fail & !is.na(conv$ess_bulk)), n_par),
    sprintf("  ess_tail: min=%.0f (%s) median=%.0f below_target_%g=%d/%d",
            .mn(conv$ess_tail), worst_tail,
            median(conv$ess_tail, na.rm = TRUE),
            ess_tail_target, sum(ess_tail_fail & !is.na(conv$ess_tail)), n_par),
    if (!all(is.na(conv$rhat)))
      sprintf("  rhat (rank-normalised split, max of bulk/folded): max=%.3f (%s) above_1.05=%d/%d above_1.01=%d/%d",
              max(conv$rhat, na.rm = TRUE),
              if (!is.na(worst_rhat)) worst_rhat else "?",
              sum(rhat_fail & !is.na(conv$rhat)), sum(!is.na(conv$rhat)),
              sum(!is.na(conv$rhat) & conv$rhat >= 1.01, na.rm = TRUE),
              sum(!is.na(conv$rhat)))
    else
      "  rhat: NA (not computable)",
    if (length(bad))
      sprintf("  not_computable: %s", paste(bad, collapse = ", ")),
    sprintf("  acf1 (chain 1): worst=%.2f (%s) top5: %s",
            acf1[worst_acf1], worst_acf1,
            paste(top5_acf, collapse = ", ")),
    nuts_lines,
    sprintf("  action: %s", action)
  ), collapse = "\n")
}
