## R/diag-post-d16-subsample.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D16 subsample stability (forest plots)
## --------------------------------------------------------------------------

#' D16. Subsample stability
#'
#' Compares posterior parameter estimates across the full sample and
#' one or more subsample splits. Produces forest plots showing posterior
#' medians and credible intervals for each subsample.
#'
#' @param results_full  Matrix or data frame (n_draws x n_params) -- full sample
#' @param results_sub   Named list of matrices/data frames -- subsample draws.
#'                      e.g. list("Pre-GFC" = draws1, "Post-GFC" = draws2)
#' @param param_names   Character vector (optional)
#' @param split_labels  Character vector -- labels for the subsamples
#' @param ci_level      Credible interval level (default 0.90)
#' @param mean_sd_tol   Maximum allowable |subsample mean - full mean| in units
#'   of the full-sample posterior SD (default 2.0). Parameters whose subsample
#'   mean drifts more than this many posterior SDs from the full-sample mean
#'   are flagged even if their CIs overlap, because CI overlap is vacuously
#'   satisfied by wide (poorly-identified) posteriors.  Gate:
#'   pass = all(CI_overlap) AND all(|sub_mean - full_mean| <= mean_sd_tol * full_sd).
#'   Reference: Lubik & Schorfheide (2004) discuss parameter stability across
#'   sub-samples as a test for indeterminacy.
#' @return dynhr_diagnostic list with forest plot
#' @references Lubik, T. A., & Schorfheide, F. (2004). Testing for indeterminacy:
#'   An application to US monetary policy. \emph{American Economic Review}, 94(1), 190-217.
#' @noRd
d16_subsample_stability <- function(results_full,
                                    results_sub,
                                    param_names  = NULL,
                                    split_labels = NULL,
                                    ci_level     = 0.90,
                                    mean_sd_tol  = 2.0,
                                    meta         = NULL) {

    results_full <- as.matrix(results_full)
    n_par <- ncol(results_full)
    if (is.null(param_names)) {
      param_names <- if (!is.null(colnames(results_full))) colnames(results_full)
      else paste0("theta_", seq_len(n_par))
    }

    if (is.null(split_labels)) split_labels <- names(results_sub)
    if (is.null(split_labels)) split_labels <- paste0("Sub_", seq_along(results_sub))

    alpha <- (1 - ci_level) / 2

    # Compute summary stats for each sample
    .summarise_draws <- function(draws, label) {
      draws <- as.matrix(draws)
      data.frame(
        Parameter = param_names,
        Median    = apply(draws, 2, median),
        Lo        = apply(draws, 2, quantile, probs = alpha),
        Hi        = apply(draws, 2, quantile, probs = 1 - alpha),
        Sample    = label
      )
    }

    summary_list <- list(.summarise_draws(results_full, "Full sample"))
    for (i in seq_along(results_sub)) {
      summary_list[[i + 1]] <- .summarise_draws(results_sub[[i]], split_labels[i])
    }
    summary_df <- do.call(rbind, summary_list)

    # Gate requires BOTH:
    #   (a) CI overlap: subsample CI overlaps full-sample CI
    #   (b) Mean stability: |subsample mean - full mean| <= mean_sd_tol * full SD
    # Criterion (b) prevents wide (poorly-identified) posteriors from trivially
    # passing via CI overlap alone.  mean_sd_tol = 2.0 by default.
    full_summary <- summary_list[[1]]
    full_sd <- apply(as.matrix(results_full), 2, sd)
    overlap_check   <- rep(TRUE, n_par)
    mean_drift_check <- rep(TRUE, n_par)

    for (i in seq_along(results_sub)) {
      sub_summary <- summary_list[[i + 1]]
      sub_draws   <- as.matrix(results_sub[[i]])
      sub_means   <- colMeans(sub_draws)
      full_means  <- colMeans(as.matrix(results_full))
      for (j in seq_len(n_par)) {
        full_lo <- full_summary$Lo[j]
        full_hi <- full_summary$Hi[j]
        sub_lo  <- sub_summary$Lo[j]
        sub_hi  <- sub_summary$Hi[j]
        # (a) CI overlap
        if (sub_hi < full_lo || sub_lo > full_hi) {
          overlap_check[j] <- FALSE
        }
        # (b) Mean stability
        sd_j <- max(full_sd[j], 1e-16)
        if (abs(sub_means[j] - full_means[j]) > mean_sd_tol * sd_j) {
          mean_drift_check[j] <- FALSE
        }
      }
    }

    stable_check    <- overlap_check & mean_drift_check
    pass            <- all(stable_check)
    unstable_params <- param_names[!stable_check]
    # Distinguish the two failure modes for reporting
    ci_fail_params   <- param_names[!overlap_check]
    drift_fail_params <- param_names[!mean_drift_check]

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    n_samples <- length(unique(summary_df$Sample))
    pal <- if (n_samples <= 3) c(dynhr_colours$dark_blue, dynhr_colours$orange, dynhr_colours$teal)
    else dynhr_palette[seq_len(n_samples)]

    summary_df$Parameter <- factor(summary_df$Parameter, levels = rev(param_names))

    # Facet by parameter with free x-scales so small- and large-scale
    # parameters each fill their own panel -- avoids cramping on shared axis.
    p_for <- ggplot2::ggplot(summary_df,
                           ggplot2::aes(x = Median, y = Sample,
                               colour = Sample, shape = Sample)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lo, xmax = Hi),
                     height = 0.35, linewidth = 0.5) +
      ggplot2::facet_wrap(~ Parameter, scales = "free_x",
                          ncol = min(4L, ceiling(sqrt(n_par)))) +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_x_continuous(n.breaks = 3L) +
      theme_dynhr_compact() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.75)),
        axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.70)),
        legend.position = "none"
      ) +
      ggplot2::labs(title = "D16: Subsample stability -- parameter comparison",
           subtitle = sprintf("%d%% credible intervals | %d subsamples | free x-scale per parameter",
                              round(ci_level * 100), length(results_sub)),
           x = "Parameter value")
    plots$forest <- .apply_meta(p_for, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(summary_df       = summary_df,
                     unstable_params  = unstable_params,
                     overlap_check    = overlap_check,
                     mean_drift_check = mean_drift_check,
                     ci_fail_params   = ci_fail_params,
                     drift_fail_params = drift_fail_params),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "D16 Subsample stability: %d params, %d subsamples. %s",
        n_par, length(results_sub),
        if (pass) {
          "PASS -- all subsample CIs overlap and means within 2 posterior SD of full sample."
        } else {
          parts <- c()
          if (length(ci_fail_params) > 0)
            parts <- c(parts, sprintf("CI non-overlap: %s", paste(ci_fail_params, collapse = ", ")))
          if (length(drift_fail_params) > 0)
            parts <- c(parts, sprintf("mean drift >%.1fSD: %s", mean_sd_tol, paste(drift_fail_params, collapse = ", ")))
          sprintf("FAIL -- %d unstable parameter(s). %s",
                  length(unstable_params), paste(parts, collapse = "; "))
        }
      ),
      llm_summary = {
        n_params   <- n_par
        badge      <- if (isTRUE(pass)) "PASS" else "FAIL"
        n_unstable <- length(unstable_params)
        paste(c(
          sprintf("D16 | Subsample Stability | %s", badge),
          sprintf("  params=%d unstable=%d/%d mean_sd_tol=%.1f", n_params, n_unstable, n_params, mean_sd_tol),
          if (length(ci_fail_params) > 0)
            sprintf("  ci_nonoverlap: %s", paste(ci_fail_params, collapse = ", ")),
          if (length(drift_fail_params) > 0)
            sprintf("  mean_drift_>%.1fSD: %s", mean_sd_tol, paste(drift_fail_params, collapse = ", ")),
          sprintf("  action: %s",
                  if (isTRUE(pass))
                    "All parameter estimates stable across subsamples (CI overlap + mean stability)."
                  else if (length(drift_fail_params) > 0)
                    sprintf("%s have subsample means drifting >%.1f full-sample posterior SD. Wide CIs do not mask this instability. Check for structural breaks or consider time-varying parameter model.",
                            paste(head(drift_fail_params, 3), collapse = ", "), mean_sd_tol)
                  else
                    sprintf("%s CI non-overlap across subsamples. Check for structural breaks.",
                            paste(head(ci_fail_params, 3), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}
