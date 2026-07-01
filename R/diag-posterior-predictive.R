## R/diag-posterior-predictive.R
## --------------------------------------------------------------------------
## Posterior predictive checks (PPC) for DSGE models.
##
## diag_posterior_predictive():
##   Draws theta from the posterior, solves the model at each draw,
##   collects model-implied moments, and compares the posterior-predictive
##   distribution to data sample moments.
##
## Analogous to D4 (prior predictive) but uses posterior draws.  This is
## the standard Bayesian model-checking tool (Gelman et al. 1996; BDA3).
##
## For DSGE specifically: moment-based PPC checks that the posterior places
## substantial probability mass on data-consistent model behaviour.  The
## "test statistic" T(y) is the empirical sample moment; the posterior
## predictive p-value (PPP) is Pr(T(yrep) >= T(y) | y).  PPP near 0 or 1
## flags systematic misfit.
##
## LLM summary: compact key=value block with PPP values per moment +
## action advice if any PPP is outside [0.05, 0.95].
## --------------------------------------------------------------------------

#' Posterior predictive checks for a DSGE model
#'
#' Draws \code{n_draws} parameter vectors uniformly from the posterior chains
#' (with replacement), solves the model at each draw, collects model-implied
#' moments, and computes posterior-predictive p-values (PPP) by comparing
#' the posterior-predictive distribution of each moment to the corresponding
#' empirical data value.
#'
#' A PPP near 0 means the model systematically under-produces that moment
#' relative to the data; near 1 means it over-produces.  PPP outside
#' \code{[ppp_warn_lo, ppp_warn_hi]} are flagged.
#'
#' @param chains        \code{dynhr_chains} object or draw matrix
#'   (n_draws x n_params) from the estimated posterior.
#' @param model_solve_fn Function \code{theta -> named numeric} returning
#'   model-implied moments (same interface as used by D4).
#' @param data_moments  Named numeric vector of empirical sample moments.
#' @param moment_names  Character vector (optional; inferred from
#'   \code{data_moments}).
#' @param n_ppc         Number of posterior draws to use for PPC
#'   (default 500; capped at number of available draws).
#' @param ppp_warn_lo   Lower PPP threshold for flagging (default 0.05).
#' @param ppp_warn_hi   Upper PPP threshold for flagging (default 0.95).
#' @param meta          Optional \code{\link{diag_meta}} provenance object.
#' @return A \code{dynhr_diagnostic} object with \code{$result$ppp_values},
#'   \code{$result$post_moments}, \code{$plots}, \code{$summary}, and
#'   \code{$llm_summary}.
#' @export
diag_posterior_predictive <- function(chains,
                                       model_solve_fn,
                                       data_moments,
                                       moment_names  = NULL,
                                       n_ppc         = 500L,
                                       ppp_warn_lo   = 0.05,
                                       ppp_warn_hi   = 0.95,
                                       meta          = NULL) {

    draws <- if (inherits(chains, "dynhr_chains")) chains$chain
             else as.matrix(chains)

    if (is.null(moment_names))
      moment_names <- names(data_moments) %||%
                      paste0("m_", seq_along(data_moments))
    n_mom <- length(data_moments)

    # Sample from posterior (with replacement if n_ppc > n_draws)
    n_avail  <- nrow(draws)
    n_use    <- min(n_ppc, n_avail)
    idx      <- sample.int(n_avail, size = n_use, replace = (n_use > n_avail))
    post_sub <- draws[idx, , drop = FALSE]

    # --- Compute model moments at each posterior draw ----------------------
    post_moments <- matrix(NA_real_, nrow = n_use, ncol = n_mom)
    colnames(post_moments) <- moment_names
    n_success <- 0L

    for (i in seq_len(n_use)) {
      theta <- post_sub[i, ]
      mom   <- model_solve_fn(theta)
      if (!is.null(mom) && all(is.finite(mom))) {
        n_success <- n_success + 1L
        post_moments[n_success, ] <- mom[seq_len(n_mom)]
      }
    }

    post_moments <- post_moments[seq_len(n_success), , drop = FALSE]

    if (n_success < 10L) {
      return(.make_result(
        pass    = NA,
        summary = sprintf(
          "diag_posterior_predictive: only %d/%d draws produced valid moments. Posterior may be near a singularity.",
          n_success, n_use)
      ))
    }

    # --- Posterior predictive p-values -------------------------------------
    ppp <- vapply(seq_len(n_mom), function(j) {
      mean(post_moments[, j] >= data_moments[j], na.rm = TRUE)
    }, numeric(1))
    names(ppp) <- moment_names

    # Flag: extreme PPP values
    flagged  <- (ppp < ppp_warn_lo) | (ppp > ppp_warn_hi)
    n_flag   <- sum(flagged)
    pass     <- n_flag == 0L

    # --- Coverage fraction -------------------------------------------------
    # Fraction of moments where data value lies within 90% PPI of model
    lo90 <- apply(post_moments, 2, quantile, 0.05, na.rm = TRUE)
    hi90 <- apply(post_moments, 2, quantile, 0.95, na.rm = TRUE)
    in90 <- (data_moments >= lo90) & (data_moments <= hi90)

    # --- Plots -------------------------------------------------------------
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE) && n_success >= 2L) {

      # Long-form posterior-predictive draws
      pp_long <- do.call(rbind, lapply(seq_len(n_mom), function(j) {
        data.frame(moment = moment_names[j],
                   value  = post_moments[, j],
                   stringsAsFactors = FALSE)
      }))
      pp_long$moment <- factor(pp_long$moment, levels = moment_names)

      data_df <- data.frame(
        moment = factor(moment_names, levels = moment_names),
        value  = as.numeric(data_moments)
      )

      # Colour data line: red if flagged, green if in 90% PPI
      data_df$colour <- ifelse(
        flagged[match(data_df$moment, moment_names)], "Flagged", "In 90% PPI"
      )
      col_in  <- dynhr_colours$green  %||% "#107C10"
      col_out <- dynhr_colours$red    %||% "#A80000"

      plots$ppc_density <- ggplot2::ggplot(pp_long,
                                            ggplot2::aes(x = value)) +
        ggplot2::geom_density(
          fill  = dynhr_colours$light_blue %||% "#D6E9F8",
          colour = dynhr_colours$mid_blue  %||% "#1B7CB6",
          alpha = 0.5, linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data_df,
          ggplot2::aes(xintercept = value, colour = colour),
          linewidth = 0.8
        ) +
        ggplot2::scale_colour_manual(
          values = c("Flagged" = col_out, "In 90% PPI" = col_in),
          name = NULL
        ) +
        ggplot2::facet_wrap(~ moment, scales = "free", ncol = 3L) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title    = "Posterior predictive checks",
          subtitle = sprintf(
            "%d valid draws  |  vertical line = data  |  %d/%d moments outside [%.0f%%, %.0f%%] PPI",
            n_success, n_flag, n_mom,
            ppp_warn_lo * 100, ppp_warn_hi * 100
          ),
          x = "Moment value", y = "Posterior-predictive density"
        )
      if (!is.null(meta)) plots$ppc_density <- .apply_meta(plots$ppc_density, meta)

      # PPP bar chart
      ppp_df <- data.frame(
        moment = factor(moment_names, levels = moment_names),
        ppp    = ppp,
        status = ifelse(flagged, "Flagged", "OK")
      )
      plots$ppc_ppp <- ggplot2::ggplot(
        ppp_df, ggplot2::aes(x = moment, y = ppp, fill = status)
      ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_hline(yintercept = ppp_warn_lo, linetype = "dashed",
                            colour = col_out, linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = ppp_warn_hi, linetype = "dashed",
                            colour = col_out, linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted",
                            colour = "grey50", linewidth = 0.3) +
        ggplot2::scale_fill_manual(
          values = c("Flagged" = col_out, "OK" = col_in), name = NULL
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
          angle = 45, hjust = 1)) +
        ggplot2::labs(
          title    = "Posterior predictive p-values (PPP)",
          subtitle = sprintf(
            "Dashed lines = warning thresholds [%.2f, %.2f] | PPP near 0.5 = good fit",
            ppp_warn_lo, ppp_warn_hi
          ),
          x = NULL, y = "PPP"
        )
      if (!is.null(meta)) plots$ppc_ppp <- .apply_meta(plots$ppc_ppp, meta)
    }

    # --- Summaries ---------------------------------------------------------
    flag_names <- moment_names[flagged]
    ppp_str    <- paste(sprintf("%s=%.2f", moment_names, ppp), collapse = ", ")
    cov_pct    <- round(100 * mean(in90), 1)

    console_summary <- sprintf(
      "PPC [%s]: %d/%d draws valid. %d/%d moments outside [%.2f, %.2f] PPP. Data in 90%% PPI: %d/%d (%.0f%%). %s",
      if (pass) "PASS" else "FAIL",
      n_success, n_use, n_flag, n_mom, ppp_warn_lo, ppp_warn_hi,
      sum(in90), n_mom, cov_pct,
      if (!pass) paste("Flagged:", paste(flag_names, collapse = ", ")) else ""
    )

    # Build PPP table for LLM
    ppp_table <- paste(vapply(seq_len(n_mom), function(j) {
      flag_str <- if (flagged[j]) " [FLAG]" else ""
      sprintf("    %s: ppp=%.2f in90=%-3s model_median=%.3f data=%.3f%s",
              moment_names[j], ppp[j],
              if (in90[j]) "YES" else "NO",
              median(post_moments[, j], na.rm = TRUE),
              data_moments[j],
              flag_str)
    }, character(1)), collapse = "\n")

    action <- if (pass) {
      sprintf("All %d moments within acceptable PPP range [%.2f, %.2f]. Posterior fit looks reasonable.",
              n_mom, ppp_warn_lo, ppp_warn_hi)
    } else {
      flagged_detail <- vapply(flag_names, function(fn) {
        p <- ppp[fn]
        if (p < ppp_warn_lo)
          sprintf("%s: ppp=%.2f -- model systematically under-produces this moment; check model structure or prior", fn, p)
        else
          sprintf("%s: ppp=%.2f -- model systematically over-produces this moment; check shock variances or steady state", fn, p)
      }, character(1))
      paste(flagged_detail, collapse = "; ")
    }

    llm_summary <- paste(c(
      sprintf("PPC | Posterior Predictive Checks | %s",
              if (pass) "PASS" else "FAIL"),
      sprintf("  draws_used=%d valid=%d moments=%d",
              n_use, n_success, n_mom),
      sprintf("  data_in_90pct_ppi=%d/%d (%.0f%%)",
              sum(in90), n_mom, cov_pct),
      "  ppp_values (ppp near 0.5 = good fit, <0.05 or >0.95 = flagged):",
      ppp_table,
      sprintf("  action: %s", action)
    ), collapse = "\n")

    .make_result(
      result = list(
        ppp_values    = ppp,
        post_moments  = post_moments,
        data_moments  = data_moments,
        flagged       = flagged,
        n_success     = n_success,
        n_ppc         = n_use,
        coverage_90   = in90
      ),
      pass        = pass,
      plots       = plots,
      summary     = console_summary,
      llm_summary = llm_summary
    )
}
