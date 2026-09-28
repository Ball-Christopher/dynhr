## R/diag-post-d14-bayes-factor.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D14 Bayes factor / marginal likelihood comparison
## --------------------------------------------------------------------------

## Jeffreys (1961) evidence category for log10 BF >= 0, as tabulated by
## Kass & Raftery (1995, Sec. 3.2): (0, 1/2] bare mention, (1/2, 1]
## substantial, (1, 2] strong, > 2 decisive.
.d14_jeffreys_category <- function(log10_bf) {
  ifelse(!is.finite(log10_bf), "unknown",
  ifelse(log10_bf > 2,   "decisive",
  ifelse(log10_bf > 1,   "strong",
  ifelse(log10_bf > 0.5, "substantial", "not worth more than a bare mention"))))
}

## Long table of every per-model marginal-likelihood estimate supplied as
## models_list[[i]]$mdd_table (a marginal_likelihoods() data frame): Model,
## Estimator, Log_MDD, SE. NULL when no model carries one. Reporting only --
## D14's comparison and badge use $log_marglik alone.
.d14_estimator_table <- function(models_list, model_names) {
  parts <- lapply(seq_along(models_list), function(i) {
    tb <- if (is.list(models_list[[i]])) models_list[[i]]$mdd_table else NULL
    if (!is.data.frame(tb) || !all(c("estimator", "log_mdd") %in% names(tb)) ||
        nrow(tb) == 0L)
      return(NULL)
    data.frame(Model     = model_names[i],
               Estimator = as.character(tb$estimator),
               Log_MDD   = as.numeric(tb$log_mdd),
               SE        = if ("se" %in% names(tb)) as.numeric(tb$se)
                           else NA_real_,
               stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) return(NULL)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' D14. Bayes factor model comparison
#'
#' Compares models via log marginal likelihoods supplied by the caller (e.g.
#' \code{\link{laplace_log_marglik}}, \code{\link{thames_mdd}}, the SMC
#' evidence, or Dynare's estimate; D14 itself does not estimate them). All
#' arithmetic is on the natural-log scale: log Bayes factors
#' \eqn{\ln BF_{m,\mathrm{best}} = \ln p(Y|m) - \ln p(Y|\mathrm{best})} and
#' posterior model probabilities
#' \eqn{p(m|Y) = \exp\{\ln\pi_m + \ln p(Y|m) - \mathrm{LSE}_k(\ln\pi_k + \ln p(Y|k))\}}.
#'
#' @param models_list  A list whose elements are lists with
#'                     $name        -- model name (character; falls back to
#'                                     the list name, then "Model i"),
#'                     $log_marglik -- natural-log marginal likelihood
#'                                     (numeric scalar, required),
#'                     $n_params    -- number of estimated parameters (optional),
#'                     $mdd_table   -- optional \code{\link{marginal_likelihoods}}
#'                                     table (estimator / log_mdd / se). Reported
#'                                     next to the comparison as
#'                                     \code{result$estimators}; it never
#'                                     changes \code{log_marglik} or the badge.
#' @param prior_odds   Prior model probabilities (non-negative, one per model;
#'                     normalised to sum to 1). Default: uniform.
#' @param log_marglik_se Optional Monte Carlo standard errors of the log
#'   marginal likelihoods, on the natural-log scale. A named numeric vector
#'   (matched to the model names; missing models get \code{NA}) or an unnamed
#'   vector of length \code{length(models_list)} taken positionally. SMC and
#'   bridge-sampling estimators supply one; a Laplace approximation or a
#'   single-chain modified harmonic mean generally does not, so the default
#'   \code{NULL} keeps the point-estimate gate.
#' @param se_multiple  Width of the simulation-noise band, in combined SEs
#'   (default 2). \strong{A package choice}: 2 SE is the standard
#'   normal-approximation "not distinguishable" convention used for bridge- and
#'   SMC-estimated log-evidences; Kass & Raftery specify no such multiplier.
#' @param meta         Plot caption metadata (see \code{.apply_meta}).
#' @return A \code{dynhr_diagnostic} list. \code{result$comparison} holds, per
#'   model, the log ML, the natural-log and log10 Bayes factor against the best
#'   model, the prior and posterior model probabilities. \code{result} also
#'   carries \code{best_model}, \code{runner_up}, \code{log10_bf_runner_up}
#'   (log10 BF of the best model against its CLOSEST competitor) and its
#'   Jeffreys \code{evidence} category.
#'   \code{pass = TRUE} when the best model has at least Jeffreys "strong"
#'   evidence (log10 BF > 1, i.e. ln BF > 2.303) against EVERY competitor,
#'   i.e. against the runner-up; a tie for best is \code{FALSE}.
#'   \code{pass = NA} when fewer than two models are given, any log marginal
#'   likelihood is non-finite, or (with \code{log_marglik_se}) the best-versus-
#'   runner-up gap is inside the simulation-noise band. The gate is on the
#'   Bayes factor, so it does not depend on \code{prior_odds}.
#' @section Simulation noise:
#' Marginal-likelihood estimates are themselves Monte Carlo quantities.
#' Comparing two models whose estimated log evidences differ by less than the
#' noise in those estimates is not a comparison at all -- the standard
#' cautionary practice is to re-run with different seeds and check stability.
#' When \code{log_marglik_se} is supplied and the best-versus-runner-up gap is
#' within \code{se_multiple} times the combined SE
#' \eqn{\sqrt{se_{best}^2 + se_{runner}^2}}, the badge is downgraded to INFO
#' with an explicit message. This is a noise check LAYERED ON TOP of, and
#' independent of, the Kass-Raftery evidential scale, whose 0.5 / 1 / 2 log10
#' bands are unchanged.
#' @references Kass, R. E., & Raftery, A. E. (1995). Bayes factors.
#'   \emph{Journal of the American Statistical Association}, 90(430), 773-795.
#'   Jeffreys, H. (1961). \emph{Theory of Probability} (3rd ed.). Oxford University Press.
#' @noRd
d14_bayes_factor <- function(models_list,
                             prior_odds = NULL,
                             log_marglik_se = NULL,
                             se_multiple = 2,
                             meta = NULL) {

    if (!is.list(models_list)) {
      .dynhr_abort("d14_bayes_factor: `models_list` must be a list of model lists.")
    }
    n_models <- length(models_list)
    if (n_models < 2) {
      return(.make_result(
        pass    = NA,
        summary = "D14 Bayes factor: Need at least 2 models for comparison."
      ))
    }

    list_names <- names(models_list) %||% rep("", n_models)
    model_names <- vapply(seq_len(n_models), function(i) {
      m  <- models_list[[i]]
      nm <- if (is.list(m)) m$name else NULL
      if (is.character(nm) && length(nm) == 1L && !is.na(nm) && nzchar(nm)) nm
      else if (!is.na(list_names[i]) && nzchar(list_names[i])) list_names[i]
      else sprintf("Model %d", i)
    }, character(1))
    model_names <- make.unique(model_names, sep = "_")

    log_ml <- vapply(seq_len(n_models), function(i) {
      m <- models_list[[i]]
      v <- if (is.list(m)) m$log_marglik else NULL
      if (length(v) == 1L && is.logical(v) && is.na(v)) return(NA_real_)
      if (!is.numeric(v) || length(v) != 1L) {
        .dynhr_abort("d14_bayes_factor: model '", model_names[i],
                     "' needs a numeric scalar `log_marglik`.")
      }
      as.numeric(v)
    }, numeric(1))
    n_params <- vapply(models_list, function(m) {
      k <- if (is.list(m)) m$n_params else NULL
      if (is.numeric(k) && length(k) == 1L) as.numeric(k) else NA_real_
    }, numeric(1))

    ## MC standard errors of the log-ML estimates: named -> matched by model
    ## name, unnamed -> positional. Everything else is NA (no noise gate).
    se_vec <- stats::setNames(rep(NA_real_, n_models), model_names)
    if (!is.null(log_marglik_se)) {
      if (!is.numeric(log_marglik_se))
        .dynhr_abort("d14_bayes_factor: `log_marglik_se` must be numeric.")
      if (!is.null(names(log_marglik_se))) {
        hit <- intersect(model_names, names(log_marglik_se))
        se_vec[hit] <- as.numeric(log_marglik_se[hit])
      } else if (length(log_marglik_se) == n_models) {
        se_vec[] <- as.numeric(log_marglik_se)
      } else {
        .dynhr_abort("d14_bayes_factor: unnamed `log_marglik_se` must have ",
                     n_models, " elements (one per model).")
      }
      if (any(!is.na(se_vec) & se_vec < 0))
        .dynhr_abort("d14_bayes_factor: `log_marglik_se` must be non-negative.")
    }
    if (length(se_multiple) != 1L || !is.finite(se_multiple) || se_multiple < 0)
      .dynhr_abort("d14_bayes_factor: `se_multiple` must be a single ",
                   "non-negative number.")

    uniform_prior <- is.null(prior_odds)
    if (uniform_prior) prior_odds <- rep(1, n_models)
    if (!is.numeric(prior_odds) || length(prior_odds) != n_models ||
        any(!is.finite(prior_odds)) || any(prior_odds < 0) || sum(prior_odds) <= 0) {
      .dynhr_abort("d14_bayes_factor: `prior_odds` must be ", n_models,
                   " finite non-negative numbers with a positive sum.")
    }
    prior_prob <- prior_odds / sum(prior_odds)
    estimators <- .d14_estimator_table(models_list, model_names)

    finite <- is.finite(log_ml)
    if (!all(finite)) {
      comparison <- data.frame(Model = model_names, Log_MargLik = log_ml,
                               Prior_Prob = prior_prob, N_Params = n_params)
      return(.make_result(
        result  = list(comparison = comparison, log_ml = log_ml,
                       n_models = n_models, estimators = estimators),
        pass    = NA,
        plots   = list(),
        summary = sprintf(
          "D14 Bayes factor: non-finite log marginal likelihood for %s (estimator not run or failed); no comparison made.",
          paste0("'", model_names[!finite], "'", collapse = ", "))
      ))
    }

    # Log Bayes factors (natural log) relative to the best model
    best_idx  <- which.max(log_ml)
    log_bf    <- log_ml - log_ml[best_idx]
    log10_bf  <- log_bf / log(10)
    runner_idx <- seq_len(n_models)[-best_idx][which.max(log_ml[-best_idx])]
    # log10 BF of best vs its CLOSEST competitor (0 on a tie)
    log10_bf_runner <- (log_ml[best_idx] - log_ml[runner_idx]) / log(10)
    evidence <- .d14_jeffreys_category(log10_bf_runner)

    # Posterior model probabilities via log-sum-exp (zero prior -> -Inf -> 0)
    lw   <- log(prior_prob) + log_ml
    lmax <- max(lw)
    post_probs <- exp(lw - (lmax + log(sum(exp(lw - lmax)))))

    comparison <- data.frame(
      Model            = model_names,
      Log_MargLik      = log_ml,
      Log_BF_vs_best   = log_bf,
      Log10_BF_vs_best = log10_bf,
      BF_vs_best       = exp(log_bf),
      Prior_Prob       = prior_prob,
      Posterior_Prob   = post_probs,
      N_Params         = n_params
    )
    rownames(comparison) <- NULL

    pass <- log10_bf_runner > 1

    ## Simulation-noise gate: a gap smaller than the Monte Carlo error of the
    ## two log-ML ESTIMATES is not evidence either way, whichever side of the
    ## Kass-Raftery bands it happens to land on.
    gap_nats     <- log_ml[best_idx] - log_ml[runner_idx]
    combined_se  <- sqrt(se_vec[[best_idx]]^2 + se_vec[[runner_idx]]^2)
    noise_band   <- se_multiple * combined_se
    within_noise <- isTRUE(is.finite(noise_band) && gap_nats < noise_band)
    if (within_noise) pass <- NA
    badge <- .badge_str(list(pass = pass, errored = FALSE))
    noise_note <- if (within_noise) sprintf(
      paste0(" [INFO: the %.2f-nat gap is inside the %.1f x combined-SE ",
             "simulation-noise band (%.2f nats; SE %.2f / %.2f): the two ",
             "models are not distinguishable by these log-ML ESTIMATES, ",
             "whatever the Bayes factor says. Re-run the estimator with ",
             "different seeds/tuning.]"),
      gap_nats, se_multiple, noise_band,
      se_vec[[best_idx]], se_vec[[runner_idx]]) else ""

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # (a) Evidence against each model relative to the best, on the Jeffreys
    #     log10 scale, with the category boundaries (gate at 1) drawn.
    ev_df <- data.frame(
      Model   = factor(model_names, levels = rev(model_names)),
      against = -log10_bf,
      is_best = seq_len(n_models) == best_idx,
      label   = ifelse(seq_len(n_models) == best_idx,
                       sprintf("best  (ln ML = %.2f)", log_ml),
                       sprintf("%.2f  (ln BF = %.2f)", -log10_bf, -log_bf))
    )
    # Jeffreys boundaries; staggered label heights so they never collide.
    thr <- data.frame(x = c(0.5, 1, 2),
                      lab = c("substantial", "strong (gate)", "decisive"),
                      lty = c("dotted", "dashed", "dotted"),
                      vj  = c(-0.4, -1.8, -0.4))
    # pseudo-log axis: linear through the Jeffreys range (0-2), compressed
    # beyond it so a decisive outlier does not squash the boundaries.
    x_trans <- scales::pseudo_log_trans(sigma = 0.5, base = 10)
    x_max   <- max(2.5, max(ev_df$against) * 3)
    x_brks  <- c(0, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    x_brks  <- x_brks[x_brks <= x_max]
    p_bf <- ggplot2::ggplot(ev_df, ggplot2::aes(x = against, y = Model)) +
      ggplot2::geom_vline(data = thr, ggplot2::aes(xintercept = x),
                          linetype = thr$lty, colour = dynhr_colours$grey) +
      ggplot2::geom_text(data = thr, ggplot2::aes(x = x, y = Inf, label = lab),
                         inherit.aes = FALSE, vjust = thr$vj, hjust = 0.5,
                         size = 3, colour = dynhr_colours$grey) +
      ggplot2::geom_col(ggplot2::aes(fill = is_best), width = 0.6,
                        show.legend = FALSE) +
      ggplot2::geom_point(data = ev_df[ev_df$is_best, , drop = FALSE],
                          colour = dynhr_colours$teal, size = 3) +
      ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1, size = 3.3,
                         colour = dynhr_colours$dark_blue) +
      ggplot2::scale_fill_manual(values = c(`FALSE` = dynhr_colours$mid_blue,
                                            `TRUE` = dynhr_colours$teal)) +
      ggplot2::scale_x_continuous(trans = x_trans, breaks = x_brks,
                                  limits = c(0, x_max),
                                  expand = ggplot2::expansion(mult = c(0, 0))) +
      ggplot2::coord_cartesian(clip = "off") +
      theme_dynhr_diagnostic() +
      ggplot2::theme(plot.title.position = "plot",
                     plot.subtitle = ggplot2::element_text(
                       margin = ggplot2::margin(b = 26))) +
      ggplot2::labs(
        title = "D14: Bayes factor of the best model against each model",
        subtitle = sprintf(
          "Best '%s' vs runner-up '%s': log10 BF = %.3f (%s) -> %s\nGate: PASS when log10 BF > 1 against the runner-up",
          model_names[best_idx], model_names[runner_idx], log10_bf_runner,
          evidence, badge),
        x = "log10 BF(best : model), pseudo-log axis   [ln BF = 2.303 x log10 BF]",
        y = NULL)
    plots$log_bf <- .apply_meta(p_bf, meta)

    # (b) Posterior probability bars (horizontal: long model names stay legible)
    pp_df <- data.frame(Model = factor(model_names, levels = rev(model_names)),
                        Posterior_Prob = post_probs,
                        label = ifelse(post_probs < 0.0005 & post_probs > 0, "<0.1%",
                                       sprintf("%.1f%%", post_probs * 100)))
    prior_txt <- if (uniform_prior) "Uniform prior model probabilities" else
      paste(strwrap(paste0("Prior model probabilities: ",
                           paste(sprintf("%s = %.3g", model_names, prior_prob),
                                 collapse = ", ")), width = 95),
            collapse = "\n")
    p_pp <- ggplot2::ggplot(pp_df, ggplot2::aes(x = Posterior_Prob, y = Model)) +
      ggplot2::geom_col(fill = dynhr_colours$teal, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.15, size = 3.3,
                         colour = dynhr_colours$dark_blue) +
      ggplot2::scale_x_continuous(labels = scales::percent_format(),
                                  limits = c(0, 1.1), breaks = seq(0, 1, 0.25),
                                  expand = ggplot2::expansion(mult = c(0, 0))) +
      theme_dynhr_diagnostic() +
      ggplot2::theme(plot.title.position = "plot") +
      ggplot2::labs(title = "D14: Posterior model probabilities",
           subtitle = prior_txt,
           x = "Posterior probability", y = NULL)
    plots$post_prob <- .apply_meta(p_pp, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(comparison = comparison,
                     best_model = model_names[best_idx],
                     runner_up  = model_names[runner_idx],
                     log10_bf_runner_up = log10_bf_runner,
                     evidence   = evidence,
                     log_marglik_se = se_vec,
                     gap_nats   = gap_nats,
                     combined_se = combined_se,
                     se_multiple = se_multiple,
                     within_noise_band = within_noise,
                     estimators = estimators),
      pass    = pass,
      plots   = plots,
      summary = paste0(sprintf(
        "D14 Bayes factor: %d models compared. Best model: '%s' (log ML = %.2f), posterior prob = %.1f%%. log10 BF vs runner-up '%s' = %.2f (%s). %s",
        n_models, model_names[best_idx], log_ml[best_idx], post_probs[best_idx] * 100,
        model_names[runner_idx], log10_bf_runner, evidence,
        if (isTRUE(pass)) "[PASS: strong evidence for best model (log10 BF > 1)]"
        else if (is.na(pass)) ""
        else "[FAIL: evidence for best model is not strong (log10 BF <= 1)]"
      ), noise_note),
      llm_summary = paste(c(
        sprintf("D14 | Bayes Factor Comparison | %s", badge),
        sprintf("  models=%d best=%s (log_ml=%.2f) runner_up=%s",
                n_models, model_names[best_idx], log_ml[best_idx], model_names[runner_idx]),
        "  gate: pass = log10(BF best vs runner-up) > 1.0 (Jeffreys 'strong', Kass-Raftery 1995); INFO when the gap is inside se_multiple x combined SE",
        sprintf("  log_marglik_se: %s (combined=%s, band=%s nats, within_noise=%s)",
                paste(sprintf("%s=%s", model_names,
                              ifelse(is.na(se_vec), "NA", sprintf("%.3f", se_vec))),
                      collapse = ", "),
                if (is.finite(combined_se)) sprintf("%.3f", combined_se) else "NA",
                if (is.finite(noise_band)) sprintf("%.3f", noise_band) else "NA",
                within_noise),
        sprintf("  log10_BF_best_vs_runner_up=%.2f (ln BF=%.2f; %s)",
                log10_bf_runner, log10_bf_runner * log(10), evidence),
        sprintf("  log_marglik: %s",
                paste(sprintf("%s=%.2f", model_names, log_ml), collapse = ", ")),
        sprintf("  posterior_prob: %s",
                paste(sprintf("%s=%.3f", model_names, post_probs), collapse = ", ")),
        if (!is.null(estimators))
          sprintf("  mdd_estimators (reported, not gated): %s",
                  paste(sprintf("%s/%s=%s", estimators$Model,
                                estimators$Estimator,
                                ifelse(is.finite(estimators$Log_MDD),
                                       sprintf("%.2f", estimators$Log_MDD),
                                       "NA")),
                        collapse = ", ")),
        sprintf("  action: %s",
                if (is.na(pass))
                  sprintf("The %.2f-nat gap between %s and %s is inside the %.1f x combined-SE band (%.2f nats): re-estimate the marginal likelihoods with different seeds before preferring either.",
                          gap_nats, model_names[best_idx], model_names[runner_idx],
                          se_multiple, noise_band)
                else if (pass) sprintf("Strong evidence for %s (log10 BF=%.2f vs %s). Prefer this model.",
                                  model_names[best_idx], log10_bf_runner, model_names[runner_idx])
                else sprintf("Evidence for %s over %s is %s (log10 BF=%.2f). %s",
                             model_names[best_idx], model_names[runner_idx], evidence,
                             log10_bf_runner,
                             if (log10_bf_runner > 0.5) "Moderate support but not conclusive."
                             else "Models not clearly distinguished. Consider more data or alternative model specs."))
      ), collapse = "\n")
    )
}
