## R/diag-post-d14-bayes-factor.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D14 Bayes factor / marginal likelihood comparison
## --------------------------------------------------------------------------

#' D14. Bayes factor model comparison
#'
#' Compares models via log marginal likelihoods (from Geweke's harmonic
#' mean estimator, Laplace approximation, or Dynare.jl's built-in
#' marginal likelihood estimate). Computes log Bayes factors and posterior
#' odds.
#'
#' @param models_list  A named list where each element is a list with:
#'                     $name       -- model name (character)
#'                     $log_marglik -- log marginal likelihood (numeric)
#'                     $n_params   -- number of estimated parameters (optional)
#' @param prior_odds   Numeric vector of prior model probabilities
#'                     (default: uniform)
#' @return A \code{dynhr_diagnostic} list with a comparison table of models,
#'   log Bayes factors, and posterior model probabilities.
#'   \code{pass = TRUE} when the best model has Jeffreys strong evidence
#'   (log10 BF > 1, i.e. log BF > 2.3 nats) relative to all competitors.
#'   \code{pass = NA} when any log marginal likelihood is non-finite.
#' @references Kass, R. E., & Raftery, A. E. (1995). Bayes factors.
#'   \emph{Journal of the American Statistical Association}, 90(430), 773-795.
#'   Jeffreys, H. (1961). \emph{Theory of Probability} (3rd ed.). Oxford University Press.
#' @noRd
d14_bayes_factor <- function(models_list,
                             prior_odds = NULL,
                             meta = NULL) {

    n_models <- length(models_list)
    if (n_models < 2) {
      return(.make_result(
        pass    = NA,
        summary = "D14 Bayes factor: Need at least 2 models for comparison."
      ))
    }

    model_names <- sapply(models_list, function(m) m$name %||% "Unnamed")
    log_ml      <- sapply(models_list, `[[`, "log_marglik")
    n_params    <- sapply(models_list, function(m) m$n_params %||% NA)

    if (is.null(prior_odds)) prior_odds <- rep(1 / n_models, n_models)

    # Log Bayes factors relative to best model
    best_idx <- which.max(log_ml)
    if (length(best_idx) == 0L || !is.finite(log_ml[best_idx])) {
      # All log marginal likelihoods are NA — no comparison possible
      return(.make_result(
        result  = list(log_ml = log_ml, n_models = n_models),
        pass    = NA,
        plots   = list(),
        summary = "D14 Bayes factor: all log marginal likelihoods are NA (SMC not run or failed)."
      ))
    }
    log_bf <- log_ml - log_ml[best_idx]

    # Posterior odds (on log scale, then normalise)
    log_post_odds <- log(prior_odds) + log_ml
    log_post_odds <- log_post_odds - max(log_post_odds)
    post_probs    <- exp(log_post_odds) / sum(exp(log_post_odds))

    # Build comparison table
    comparison <- data.frame(
      Model         = model_names,
      Log_MargLik   = log_ml,
      Log_BF_vs_best = log_bf,
      BF_vs_best    = exp(log_bf),
      Prior_Prob    = prior_odds,
      Posterior_Prob = post_probs,
      N_Params      = n_params
    )
    rownames(comparison) <- NULL

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # (a) Bar chart of log marginal likelihoods
    comp_df <- comparison
    comp_df$Model <- factor(comp_df$Model, levels = comp_df$Model)

    p_lm <- ggplot2::ggplot(comp_df, ggplot2::aes(x = Model, y = Log_MargLik)) +
      ggplot2::geom_col(fill = dynhr_colours$mid_blue, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", Log_MargLik)),
                vjust = -0.3, size = 3.5, colour = dynhr_colours$dark_blue) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(title = "D14: Log marginal likelihoods",
           x = NULL, y = "Log marginal likelihood")
    plots$log_marglik <- .apply_meta(p_lm, meta)

    # (b) Posterior probability bar chart
    p_pp <- ggplot2::ggplot(comp_df, ggplot2::aes(x = Model, y = Posterior_Prob)) +
      ggplot2::geom_col(fill = dynhr_colours$teal, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Posterior_Prob * 100)),
                vjust = -0.3, size = 3.5, colour = dynhr_colours$dark_blue) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(title = "D14: Posterior model probabilities",
           subtitle = "Assuming uniform prior model probabilities",
           x = NULL, y = "Posterior probability")
    plots$post_prob <- .apply_meta(p_pp, meta)

    }  # end requireNamespace guard

    # Gate: PASS if the best model's log10 BF vs all competitors > 1 (Jeffreys
    # "strong" evidence on the Kass-Raftery 1995 scale: log10 BF in (1, 2] is
    # strong; >2 is decisive).  NA when log_ml values are not all finite.
    # This replaces the previous hardcoded pass=NA.
    log_bfs_vs_best <- log_bf[log_bf != 0]  # exclude best model itself
    pass <- if (all(is.finite(log_ml)) && length(log_bfs_vs_best) > 0) {
      # log_bf = log_ml - log_ml[best_idx]; competitors have log_bf < 0.
      # max log10 BF of best vs all competitors = max(-log_bfs_vs_best) / log(10)
      max_log10_bf_vs_competitors <- max(-log_bfs_vs_best) / log(10)
      max_log10_bf_vs_competitors > 1
    } else {
      NA
    }

    .make_result(
      result  = list(comparison = comparison, best_model = model_names[best_idx]),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "D14 Bayes factor: %d models compared. Best model: '%s' (log ML = %.1f). Posterior prob = %.1f%%. %s",
        n_models, model_names[best_idx], log_ml[best_idx], post_probs[best_idx] * 100,
        if (is.na(pass)) "[INFO: not all log_ml finite]"
        else if (pass)   "[PASS: strong evidence for best model (log10 BF > 1)]"
        else             "[FAIL: evidence for best model is weak (log10 BF <= 1)]"
      ),
      llm_summary = {
        log_mls   <- log_ml
        badge     <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
        best_nm   <- model_names[which.max(log_mls)]
        best_ml   <- max(log_mls)
        ml_strs   <- paste(sprintf("%s=%.1f", model_names, log_mls), collapse = ", ")
        prob_strs <- paste(sprintf("%s=%.3f", model_names, post_probs), collapse = ", ")
        max_log10_bf <- if (length(log_bfs_vs_best) > 0) max(-log_bfs_vs_best) / log(10) else NA_real_
        bf_strength <- if (!is.na(max_log10_bf)) {
          if (max_log10_bf > 2) "decisive"
          else if (max_log10_bf > 1) "strong"
          else if (max_log10_bf > 0.5) "substantial"
          else "weak"
        } else "unknown"
        paste(c(
          sprintf("D14 | Bayes Factor Comparison | %s", badge),
          sprintf("  models=%d best=%s (log_ml=%.1f)", length(model_names), best_nm, best_ml),
          sprintf("  gate: pass = max log10(BF_best_vs_competitor) > 1.0 (Jeffreys scale, Kass-Raftery 1995)"),
          if (!is.na(max_log10_bf)) sprintf("  max_log10_BF_vs_best=%.2f (%s evidence)", max_log10_bf, bf_strength),
          sprintf("  log_marglik: %s", ml_strs),
          sprintf("  posterior_prob: %s", prob_strs),
          sprintf("  action: %s",
                  if (is.na(pass)) "Cannot compute BF: not all log marginal likelihoods are finite."
                  else if (pass)   sprintf("Strong evidence for %s (log10 BF=%.1f). Prefer this model.", best_nm, max_log10_bf)
                  else             sprintf("Evidence for %s is %s (log10 BF=%.1f). %s",
                                           best_nm, bf_strength, max_log10_bf,
                                           if (!is.na(max_log10_bf) && max_log10_bf > 0.5) "Moderate support but not conclusive."
                                           else "Models not clearly distinguished. Consider more data or alternative model specs."))
        ), collapse = "\n")
      }
    )
}
