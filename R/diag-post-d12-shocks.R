## R/diag-post-d12-shocks.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D12 smoothed shock properties (Ljung-Box, ACF)
## --------------------------------------------------------------------------

#' D12. Smoothed shock properties
#'
#' Tests whether the smoothed (filtered) structural shocks from the Kalman
#' smoother are approximately i.i.d. N(0,1). Three separate gate criteria
#' must all pass:
#' \enumerate{
#'   \item \strong{Serial correlation}: Ljung-Box p-value > 0.05 (no
#'     significant autocorrelation up to \code{max_lag} lags).
#'   \item \strong{Mean}: |mean| < 0.1 (shocks centred near zero).
#'   \item \strong{Variance}: |var - 1| < 0.3 (shocks have unit variance,
#'     within 30 percentage points).
#' }
#'
#' @param shocks       Matrix (T x n_shocks) -- smoothed shock series
#' @param shock_names  Character vector (optional)
#' @param max_lag      Maximum lag for Ljung-Box and ACF (default 8)
#' @return dynhr_diagnostic list
#'
#' @references Ljung, G. M., & Box, G. E. P. (1978). On a measure of lack of fit in
#'   time series models. \emph{Biometrika}, 65(2), 297-303.
#' @noRd
d12_smoothed_shocks <- function(shocks,
                                shock_names = NULL,
                                max_lag = 8,
                                meta = NULL) {

    shocks <- as.matrix(shocks)
    n_t <- nrow(shocks)
    n_s <- ncol(shocks)
    if (is.null(shock_names)) {
      shock_names <- if (!is.null(colnames(shocks))) colnames(shocks)
      else paste0("eps_", seq_len(n_s))
    }

    # Ljung-Box for each shock
    lb_results <- lapply(seq_len(n_s), function(j) {
      res <- .ljung_box(shocks[, j], max_lag = max_lag)
      res$shock <- shock_names[j]
      res
    })
    names(lb_results) <- shock_names

    lb_pass_vec  <- sapply(lb_results, function(x) isTRUE(x$pass))
    mean_pass_vec <- abs(colMeans(shocks)) < 0.1
    var_pass_vec  <- abs(apply(shocks, 2, var) - 1) < 0.3
    all_pass <- all(lb_pass_vec & mean_pass_vec & var_pass_vec)

    # Mean and variance of each shock (should be ~0 and ~1)
    shock_stats <- data.frame(
      Shock    = shock_names,
      Mean     = colMeans(shocks),
      Variance = apply(shocks, 2, var),
      LB_Stat  = sapply(lb_results, `[[`, "statistic"),
      LB_pval  = sapply(lb_results, `[[`, "p_value"),
      LB_pass  = sapply(lb_results, function(x) isTRUE(x$pass))
    )

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # (a) ACF plots for each shock
    acf_data_list <- list()
    for (j in seq_len(n_s)) {
      acf_obj <- acf(shocks[, j], lag.max = max_lag, plot = FALSE)
      acf_data_list[[j]] <- data.frame(
        Lag  = as.numeric(acf_obj$lag[-1]),
        ACF  = as.numeric(acf_obj$acf[-1]),
        Shock = shock_names[j]
      )
    }
    acf_df <- do.call(rbind, acf_data_list)
    ci <- qnorm(0.975) / sqrt(n_t)

    p_acf <- ggplot2::ggplot(acf_df, ggplot2::aes(x = Lag, y = ACF)) +
      geom_dynhr_zero() +
      ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
                 colour = dynhr_colours$red, linewidth = 0.4) +
      ggplot2::geom_col(fill = dynhr_colours$mid_blue, width = 0.5) +
      ggplot2::facet_wrap(~ Shock, ncol = 3) +
      theme_dynhr_compact() +
      ggplot2::labs(title = "D12: ACF of smoothed shocks",
           subtitle = sprintf("Dashed lines = 95%% CI (+/-%.3f)", ci),
           x = "Lag")
    plots$acf <- .apply_meta(p_acf, meta)

    # (b) Time series of shocks -- standardised to unit variance per shock so
    #     the +/-2 reference lines are always visible regardless of raw scale.
    shock_sds <- apply(shocks, 2, sd)
    shock_sds[shock_sds < 1e-12] <- 1  # guard against zero-variance shocks
    shocks_std <- sweep(shocks, 2, shock_sds, "/")

    shocks_df <- as.data.frame(shocks_std)
    colnames(shocks_df) <- shock_names
    shock_long <- reshape2::melt(shocks_df,
                                 measure.vars   = shock_names,
                                 variable.name  = "Shock",
                                 value.name     = "Value")
    shock_long$Time <- rep(seq_len(n_t), n_s)

    p_ss <- ggplot2::ggplot(shock_long, ggplot2::aes(x = Time, y = Value)) +
      geom_dynhr_zero() +
      ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted",
                 colour = dynhr_colours$red, linewidth = 0.4) +
      ggplot2::facet_wrap(~ Shock, scales = "free_y", ncol = 3) +
      theme_dynhr_compact() +
      ggplot2::labs(title = "D12: Smoothed structural shocks (standardised)",
           subtitle = "Each shock divided by its own sd (unit variance). Dotted lines at +/-2.",
           x = "Time")
    plots$shock_series <- .apply_meta(p_ss, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(shock_stats = shock_stats, lb_results = lb_results),
      pass    = all_pass,
      plots   = plots,
      summary = sprintf(
        "D12 Smoothed shocks: %d shocks, %d periods. Ljung-Box (lag=%d): %d/%d pass (p > 0.05). %s",
        n_s, n_t, max_lag, sum(shock_stats$LB_pass), n_s,
        ifelse(all_pass, "Shocks appear serially uncorrelated.",
               sprintf("Serial correlation detected in: %s",
                       paste(shock_names[!shock_stats$LB_pass], collapse = ", ")))
      ),
      llm_summary = {
        pass     <- all_pass
        n_shocks <- n_s
        badge    <- if (isTRUE(pass)) "PASS" else "FAIL"
        lb_pass  <- vapply(lb_results, function(r) isTRUE(r$pass), logical(1))
        lb_pvals <- vapply(lb_results, function(r) r$p_value %||% NA_real_, numeric(1))
        paste(c(
          sprintf("D12 | Smoothed Shocks | %s", badge),
          sprintf("  shocks=%d ljung_box_pass=%d/%d (p>0.05)",
                  n_shocks, sum(lb_pass, na.rm = TRUE), n_shocks),
          if (any(!lb_pass, na.rm = TRUE)) {
            bad <- shock_names[!lb_pass & !is.na(lb_pass)]
            sprintf("  serial_correlation: %s",
                    paste(sprintf("%s(p=%.3f)", bad, lb_pvals[!lb_pass & !is.na(lb_pass)]),
                          collapse = ", "))
          },
          sprintf("  action: %s",
                  if (isTRUE(pass))
                    "Smoothed shocks show no serial correlation. Model fit looks adequate."
                  else {
                    bad <- shock_names[!lb_pass & !is.na(lb_pass)]
                    sprintf("%s: serial correlation in smoothed shocks suggests missing persistent component. Consider adding trend or AR state.",
                            paste(bad, collapse = ", "))
                  })
        ), collapse = "\n")
      }
    )
}
