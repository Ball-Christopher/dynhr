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
#' @param ess_target     Minimum Bulk-ESS threshold (default 1000).
#' @param ess_tail_target Minimum Tail-ESS threshold (default 400).
#' @param acf_max_params Max parameters in ACF plot, sorted by worst lag-1
#'   (default 12).
#' @return dynhr_diagnostic with \code{$plots}, \code{$summary}, and
#'   \code{$llm_summary}.
#' @noRd
d5_mcmc_convergence <- function(draws,
                                param_names      = NULL,
                                chains_list      = NULL,
                                nuts_meta        = NULL,
                                ess_target       = 1000,
                                ess_tail_target  = 400,
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

    # --- Convergence statistics --------------------------------------------
    if (is_multi) {
      conv <- .convergence_summary(chains_list_full)
    } else {
      ess_b <- vapply(seq_len(n_par), function(j) {
        z <- .rank_normalise(list(draws[, j]))[[1]]
        .effective_sample_size(z)
      }, numeric(1))
      ess_t <- vapply(seq_len(n_par), function(j) {
        x   <- draws[, j]
        q05 <- quantile(x, 0.05); q95 <- quantile(x, 0.95)
        min(.effective_sample_size(as.numeric(x <= q05)),
            .effective_sample_size(as.numeric(x <= q95)))
      }, numeric(1))
      conv <- data.frame(param = param_names, rhat = NA_real_,
                         ess_bulk = ess_b, ess_tail = ess_t,
                         stringsAsFactors = FALSE)
    }

    # Lag-1 ACF per parameter (first chain)
    acf1 <- vapply(seq_len(n_par), function(j)
      acf(draws[, j], lag.max = 1L, plot = FALSE)$acf[2, 1, 1],
      numeric(1))
    names(acf1) <- param_names

    # --- Pass/fail ---------------------------------------------------------
    ess_bulk_fail <- conv$ess_bulk < ess_target
    ess_tail_fail <- conv$ess_tail < ess_tail_target
    rhat_fail     <- !is.na(conv$rhat) & conv$rhat >= 1.05
    pass <- !any(ess_bulk_fail | ess_tail_fail | rhat_fail, na.rm = TRUE)

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
      tgt_df <- data.frame(type   = c("Bulk", "Tail"),
                           target = c(ess_target, ess_tail_target),
                           stringsAsFactors = FALSE)

      p_ess <- ggplot2::ggplot(
        ess_df,
        ggplot2::aes(x = param, y = ess,
                     fill = interaction(type, status))
      ) +
        ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
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
          subtitle = sprintf("Bulk target >= %d  |  Tail target >= %d",
                              ess_target, ess_tail_target),
          x = NULL, y = "ESS")
      attr(p_ess, "dynhr_fig_height") <- 7.0
      plots$ess <- .apply_meta(p_ess, meta)

      # (b) Rank-normalised R-hat -------------------------------------------
      if (is_multi && !all(is.na(conv$rhat))) {
        rhat_status <- ifelse(is.na(conv$rhat), "Unknown",
                        ifelse(conv$rhat < 1.05, "Converged",
                          ifelse(conv$rhat < 1.1, "Marginal", "Not converged")))
        rhat_df <- data.frame(
          param  = factor(param_names, levels = param_names),
          rhat   = conv$rhat,
          status = rhat_status,
          group  = .param_group(param_names),
          stringsAsFactors = FALSE
        )
        p_rhat <- ggplot2::ggplot(
          rhat_df, ggplot2::aes(x = param, y = rhat, fill = status)
        ) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::geom_hline(yintercept = 1.05, linetype = "dashed",
                              colour = col_fail, linewidth = 0.5) +
          ggplot2::geom_hline(yintercept = 1.01, linetype = "dotted",
                              colour = col_warn, linewidth = 0.4) +
          ggplot2::scale_fill_manual(
            values = c("Converged"    = col_pass,
                       "Marginal"     = col_warn,
                       "Not converged"= col_fail,
                       "Unknown"      = "grey70"),
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
            subtitle = "Dashed = 1.05 warning  |  Dotted = 1.01 strict  |  Goal: < 1.01",
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
          p_rank <- ggplot2::ggplot(
            ro_df,
            ggplot2::aes(x = bin_mid, y = freq,
                         colour = chain, group = chain)
          ) +
            ggplot2::geom_line(linewidth = 0.5) +
            ggplot2::geom_point(size = 1) +
            ggplot2::geom_hline(
              yintercept = unique(ro_df$expected),
              linetype = "dashed", colour = "grey50", linewidth = 0.3
            ) +
            ggplot2::facet_wrap(~ param, scales = "free_y", ncol = 4L) +
            ggplot2::scale_colour_manual(
              values = .chain_palette(n_chains), name = NULL) +
            theme_dynhr_diagnostic() +
            ggplot2::theme(strip.text = ggplot2::element_text(
              size = ggplot2::rel(0.75))) +
            ggplot2::labs(
              title = "D5: Rank overlay (Vehtari 2021)",
              subtitle = paste("Each chain's pooled-rank frequency should be",
                               "uniform (dashed)"),
              x = "Pooled rank", y = "Frequency")
          plots$rank_overlay <- .apply_meta(p_rank, meta)
        }
      }

      # (e) ACF bar chart (top N worst lag-1) ---------------------------------
      acf_order  <- names(sort(abs(acf1), decreasing = TRUE))
      acf_params <- acf_order[seq_len(min(length(acf_order), acf_max_params))]
      max_lag    <- min(20L, floor(n_draws / 5L))
      ci_band    <- 2 / sqrt(n_draws)

      acf_df <- do.call(rbind, lapply(acf_params, function(pn) {
        a <- acf(draws[, pn], lag.max = max_lag,
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
          subtitle = sprintf("Cells coloured by ACF; values persistently above the +/-2/sqrt(%s) band indicate slow mixing",
                              format(n_draws, big.mark = ",")),
          x = "Lag", y = NULL)
      p_acf <- .apply_meta(p_acf, meta)
      attr(p_acf, "dynhr_fig_height") <- max(4, 0.42 * length(acf_params) + 2)
      plots$acf <- p_acf

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
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                            colour = "grey50", linewidth = 0.3) +
        ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
        theme_dynhr_compact() +
        ggplot2::theme(
          # Suppress "phi"->"phl" / "sig"->"slg" font ligature in strip labels
          strip.text = ggplot2::element_text(family = "sans")
        ) +
        ggplot2::labs(
          title    = "D5: Posterior intervals (all parameters)",
          subtitle = "Median  +  50% CrI (thick)  +  90% CrI (thin)  |  Faceted by parameter group",
          x = "Parameter value", y = NULL)
      plots$intervals <- .apply_meta(p_int, meta)

      # (h) NUTS-specific plots -----------------------------------------------
      if (!is.null(nuts_meta))
        plots <- c(plots,
                   .d5_nuts_plots(nuts_meta, col_fail, col_warn, col_pass,
                                  diag_meta = meta))
    }

    # --- LLM summary -------------------------------------------------------
    llm_summary <- .d5_llm(conv, acf1, n_draws, n_chains, n_par,
                             ess_bulk_fail, ess_tail_fail, rhat_fail,
                             ess_target, ess_tail_target, pass, nuts_meta,
                             param_names)

    .make_result(
      result = list(convergence  = conv,
                    acf1         = acf1,
                    n_chains     = n_chains,
                    n_draws      = n_draws,
                    nuts_meta    = nuts_meta),
      pass        = pass,
      plots       = plots,
      summary     = .d5_console(conv, acf1, n_draws, n_chains,
                                 n_par, pass, nuts_meta),
      llm_summary = llm_summary
    )
}


#' D21. Bayesian KPS posterior precision updating indicator
#'
#' Implements the Koop-Pesaran-Smith idea using either:
#' - \code{draws_by_T}: named list of posterior draw matrices for increasing T, or
#' - \code{kps_runner_fn}: function(T) returning a draw matrix.
#'
#' For each parameter, fits log(precision) on log(T). A slope near 1 indicates
#' regular identification; substantially below 1 indicates weak identification.
#'
#' @param draws_by_T Named list of draw matrices (n_draws x n_params).
#' @param sample_sizes Optional numeric vector of sample sizes aligned to list order.
#' @param kps_runner_fn Optional function(T) -> draw matrix.
#' @param param_names Optional parameter names.
#' @param slope_threshold Minimum slope considered adequately identified.
#' @return dynhr_diagnostic object.
#' @noRd
d21_kps_precision_update <- function(draws_by_T   = NULL,
                                     sample_sizes = NULL,
                                     kps_runner_fn = NULL,
                                     param_names   = NULL,
                                     slope_threshold = 0.80,
                                     meta          = NULL) {

    if (is.null(draws_by_T) && is.null(kps_runner_fn)) {
      return(.make_result(
        pass = NA,
        summary = paste(
          "D21 KPS posterior precision update: SKIPPED -- provide draws_by_T",
          "or kps_runner_fn(T) to evaluate precision growth with sample size."
        )
      ))
    }

    if (is.null(draws_by_T) && !is.null(kps_runner_fn)) {
      if (is.null(sample_sizes) || length(sample_sizes) < 3L) {
        stop("When using kps_runner_fn, provide at least 3 sample_sizes.")
      }
      draws_by_T <- lapply(sample_sizes, kps_runner_fn)
      names(draws_by_T) <- as.character(sample_sizes)
    }

    if (is.null(sample_sizes)) {
      nm <- names(draws_by_T)
      if (!is.null(nm) && all(grepl("^[0-9]+$", nm))) {
        sample_sizes <- as.numeric(nm)
      } else {
        stop("sample_sizes not supplied and could not be inferred from names(draws_by_T).")
      }
    }
    if (length(draws_by_T) != length(sample_sizes)) {
      stop("draws_by_T and sample_sizes must have the same length.")
    }

    mats <- lapply(draws_by_T, as.matrix)
    n_par <- ncol(mats[[1]])
    if (is.null(param_names)) {
      param_names <- colnames(mats[[1]]) %||% paste0("theta_", seq_len(n_par))
    }
    for (i in seq_along(mats)) {
      if (ncol(mats[[i]]) != n_par) stop("All draws matrices must have same number of columns.")
      colnames(mats[[i]]) <- param_names
    }

    precision <- sapply(mats, function(m) 1 / (apply(m, 2, stats::sd)^2 + 1e-16))
    if (!is.matrix(precision)) precision <- matrix(precision, nrow = n_par)
    rownames(precision) <- param_names

    logT <- log(as.numeric(sample_sizes))
    slope <- vapply(seq_len(n_par), function(i) {
      y <- log(pmax(precision[i, ], 1e-16))
      coef(stats::lm(y ~ logT))[2]
    }, numeric(1))
    names(slope) <- param_names

    weak <- names(slope)[slope < slope_threshold]
    pass <- length(weak) == 0

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      growth_df <- do.call(rbind, lapply(seq_along(sample_sizes), function(k) {
        data.frame(
          T = sample_sizes[k],
          parameter = param_names,
          precision = precision[, k],
          stringsAsFactors = FALSE
        )
      }))
      growth_df$parameter <- factor(growth_df$parameter, levels = param_names)

      p_kps <- ggplot2::ggplot(
        growth_df, ggplot2::aes(x = T, y = precision, colour = parameter)
      ) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(size = 1.2) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        scale_colour_dynhr_vibrant() +
        theme_dynhr_diagnostic() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(
          title = "D21: Posterior precision growth by sample size",
          subtitle = "log(precision)-log(T) slope near 1 indicates stronger identification",
          x = "Sample size T (log scale)",
          y = "Posterior precision 1/sd^2 (log scale)"
        )
      plots$kps_precision_growth <- .apply_meta(p_kps, meta)
    }

    .make_result(
      result = list(
        sample_sizes = sample_sizes,
        precision = precision,
        slope = slope,
        weak_params = weak
      ),
      pass = pass,
      plots = plots,
      summary = sprintf(
        "D21 KPS precision update: %d parameters over %d sample sizes. slope range [%.2f, %.2f]. %s",
        n_par, length(sample_sizes), min(slope), max(slope),
        if (pass) "PASS -- all slopes above threshold."
        else sprintf("FAIL -- weak precision updating for: %s", paste(weak, collapse = ", "))
      ),
      llm_summary = paste(c(
        sprintf("D21 | KPS Precision Updating | %s", if (pass) "PASS" else "FAIL"),
        sprintf("  params=%d sample_sizes=%s slope_threshold=%.2f",
                n_par, paste(sample_sizes, collapse = ","), slope_threshold),
        sprintf("  slope_min=%.2f (%s) slope_median=%.2f slope_max=%.2f",
                min(slope), names(which.min(slope)), median(slope), max(slope)),
        if (length(weak) > 0) sprintf("  weak_updating: %s", paste(head(weak, 6), collapse = ", ")),
        sprintf("  action: %s",
                if (pass) "Posterior precision scales appropriately with data size."
                else "Some parameters do not exhibit near-T precision scaling; review local/spectral identification and prior dominance.")
      ), collapse = "\n")
    )
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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
  base <- c("#1B7CB6", "#E8540A", "#0E8040", "#9B59B6",
            "#D4AC0D", "#2E86C1", "#C0392B", "#27AE60")
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
.d5_console <- function(conv, acf1, n_draws, n_chains, n_par, pass, nuts_meta) {
  rhat_str  <- if (!all(is.na(conv$rhat)))
    sprintf(" Rhat max=%.3f.", max(conv$rhat, na.rm = TRUE))
  else " Rhat: N/A (single chain)."
  nuts_str  <- if (!is.null(nuts_meta) && (nuts_meta$n_divergent %||% 0L) > 0)
    sprintf(" NUTS: %d divergences.", nuts_meta$n_divergent)
  else ""
  sprintf(
    "D5 [%s]: %d chain(s) x %s draws, %d params. Bulk-ESS min=%.0f median=%.0f. Tail-ESS min=%.0f.%s Lag-1 ACF worst: %.2f (%s).%s",
    if (pass) "PASS" else "FAIL",
    n_chains, format(n_draws, big.mark = ","), n_par,
    min(conv$ess_bulk, na.rm = TRUE),
    median(conv$ess_bulk, na.rm = TRUE),
    min(conv$ess_tail, na.rm = TRUE),
    rhat_str,
    max(abs(acf1)), names(which.max(abs(acf1))),
    nuts_str
  )
}

#' LLM-friendly D5 summary with actionable advice
#' @noRd
.d5_llm <- function(conv, acf1, n_draws, n_chains, n_par,
                     ess_bulk_fail, ess_tail_fail, rhat_fail,
                     ess_target, ess_tail_target, pass, nuts_meta,
                     param_names) {

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
  if (any(rhat_fail, na.rm = TRUE))
    issues <- c(issues, sprintf(
      "%s: rhat=%.3f>1.05 -- chains not converged; increase draws or check for multi-modality",
      worst_rhat, max(conv$rhat, na.rm = TRUE)
    ))
  if (any(ess_bulk_fail, na.rm = TRUE))
    issues <- c(issues, sprintf(
      "%s: ess_bulk=%.0f<target_%d -- high autocorrelation; tune proposal or increase draws",
      worst_bulk, min(conv$ess_bulk, na.rm = TRUE), ess_target
    ))
  if (any(ess_tail_fail, na.rm = TRUE) && !any(ess_bulk_fail, na.rm = TRUE))
    issues <- c(issues, sprintf(
      "%s: ess_tail=%.0f<target_%d -- tails poorly sampled; check parameter boundaries",
      worst_tail, min(conv$ess_tail, na.rm = TRUE), ess_tail_target
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
    sprintf("D5 | MCMC Convergence | %s", if (pass) "PASS" else "FAIL"),
    sprintf("  chains=%d draws=%s params=%d",
            n_chains, format(n_draws, big.mark = ","), n_par),
    sprintf("  ess_bulk: min=%.0f (%s) median=%.0f below_target_%d=%d/%d",
            min(conv$ess_bulk, na.rm = TRUE), worst_bulk,
            median(conv$ess_bulk, na.rm = TRUE),
            ess_target, sum(ess_bulk_fail, na.rm = TRUE), n_par),
    sprintf("  ess_tail: min=%.0f (%s) median=%.0f below_target_%d=%d/%d",
            min(conv$ess_tail, na.rm = TRUE), worst_tail,
            median(conv$ess_tail, na.rm = TRUE),
            ess_tail_target, sum(ess_tail_fail, na.rm = TRUE), n_par),
    if (!all(is.na(conv$rhat)))
      sprintf("  rhat: max=%.3f (%s) above_1.05=%d/%d above_1.01=%d/%d",
              max(conv$rhat, na.rm = TRUE),
              if (!is.na(worst_rhat)) worst_rhat else "?",
              sum(rhat_fail, na.rm = TRUE), sum(!is.na(conv$rhat)),
              sum(!is.na(conv$rhat) & conv$rhat >= 1.01, na.rm = TRUE),
              sum(!is.na(conv$rhat)))
    else
      "  rhat: NA (single chain -- provide chains_list for multi-chain R-hat)",
    sprintf("  acf1: worst=%.2f (%s) top5: %s",
            acf1[worst_acf1], worst_acf1,
            paste(top5_acf, collapse = ", ")),
    nuts_lines,
    sprintf("  action: %s", action)
  ), collapse = "\n")
}
