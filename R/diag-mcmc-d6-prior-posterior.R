## R/diag-mcmc-d6-prior-posterior.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R; ridge-density rewrite Phase-3+.
##
## D6 posterior vs prior -- ridge density layout
## --------------------------------------------------------------------------

#' D6. Posterior vs prior updating (ridge density)
#'
#' Displays prior and posterior distributions as stacked ridge densities: one
#' row per parameter, each with its own x-axis (so parameters with very
#' different scales remain readable).  Parameters where the posterior closely
#' resembles the prior (overlap > 0.80) are flagged and sorted to the top.
#'
#' @param draws           Matrix (n_draws x n_params) -- posterior draws
#' @param prior_density_fn Function: (x, param_name) -> prior density value(s).
#'   Vectorised over \code{x}.
#' @param param_names     Character vector (optional; taken from \code{colnames(draws)})
#' @param meta            A \code{\link{diag_meta}} object for plot provenance
#'   (model name, data hash, date caption).  Pass \code{NULL} to suppress.
#' @return dynhr_diagnostic list
#' @references Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches
#'   to the calculation of posterior moments. In J. M. Bernardo et al. (eds),
#'   \emph{Bayesian Statistics 4}, Oxford University Press.
#'
#' @noRd
d6_posterior_vs_prior <- function(draws,
                                  prior_density_fn,
                                  param_names = NULL,
                                  meta        = NULL) {

    draws <- as.matrix(draws)
    n_par <- ncol(draws)
    if (is.null(param_names)) {
      param_names <- if (!is.null(colnames(draws))) colnames(draws)
      else paste0("theta_", seq_len(n_par))
    }

    # ---- Per-parameter KDE and overlap score --------------------------------
    overlap_scores <- setNames(numeric(n_par), param_names)
    density_rows   <- vector("list", n_par)

    for (j in seq_len(n_par)) {
      post_kde   <- density(draws[, j], n = 512)
      dx         <- diff(post_kde$x[1:2])

      # Modestly extend the grid (+/-15% of KDE range) so the prior is not
      # hard-clipped to the posterior support, but not so far that the panel
      # collapses for wide-support priors.
      kde_range  <- diff(range(post_kde$x))
      prior_x    <- seq(min(post_kde$x) - 0.15 * kde_range,
                        max(post_kde$x) + 0.15 * kde_range, length.out = 512L)
      prior_dx   <- diff(prior_x[1:2])
      # Evaluate elementwise: some prior_density_fn implementations are not
      # vectorised (they use scalar if() on bounds), which errors on a grid.
      .pf <- function(xs) vapply(xs, function(xx) {
        v <- tryCatch(prior_density_fn(xx, param_names[j]), error = function(e) NA_real_)
        if (length(v) != 1L || !is.finite(v)) NA_real_ else v
      }, numeric(1))
      prior_vals_wide <- .pf(prior_x)

      # Normalise prior over the grid
      prior_norm_wide <- prior_vals_wide / (sum(prior_vals_wide, na.rm = TRUE) * prior_dx + 1e-300)

      # Evaluate prior on the posterior KDE grid for overlap computation
      prior_vals <- .pf(post_kde$x)
      prior_norm_on_kde <- prior_vals / (sum(prior_vals_wide, na.rm = TRUE) * prior_dx + 1e-300)
      prior_norm_on_kde[!is.finite(prior_norm_on_kde)] <- 0

      # Overlap coefficient: integral of min(posterior, prior) on KDE grid
      overlap_scores[j] <- sum(pmin(post_kde$y, prior_norm_on_kde)) * dx

      density_rows[[j]] <- data.frame(
        x            = c(post_kde$x, prior_x),
        density      = c(post_kde$y, prior_norm_wide),
        Distribution = rep(c("Posterior", "Prior"),
                           c(length(post_kde$x), length(prior_x))),
        Parameter    = param_names[j],
        stringsAsFactors = FALSE
      )
    }

    all_df <- do.call(rbind, density_rows)

    uninformative <- param_names[overlap_scores > 0.80]

    # ---- Ridge layout -------------------------------------------------------
    # Sort: uninformative (flagged) at top so they are immediately visible.
    ordered_params <- c(
      param_names[param_names %in% uninformative],
      param_names[!param_names %in% uninformative]
    )
    # facet_grid with free scales gives each parameter its own x-axis --
    # essential because parameter values span very different magnitudes.
    all_df$Parameter <- factor(all_df$Parameter, levels = rev(ordered_params))

    flag_label <- if (length(uninformative) > 0)
      sprintf("Flagged (overlap > 0.80): %s", paste(uninformative, collapse = ", "))
    else "All parameters updated by data"

    plots <- list()

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Layer order: prior ribbon first, then posterior on top so both visible.
      prior_df     <- all_df[all_df$Distribution == "Prior",     ]
      posterior_df <- all_df[all_df$Distribution == "Posterior", ]

      p <- ggplot2::ggplot(mapping = ggplot2::aes(x = x)) +
        # --- Prior (bottom layer) ---
        ggplot2::geom_ribbon(data = prior_df,
                             ggplot2::aes(ymin = 0, ymax = density,
                                          fill = Distribution),
                             alpha = 0.45, colour = NA) +
        ggplot2::geom_line(data = prior_df,
                           ggplot2::aes(y = density, colour = Distribution),
                           linewidth = 0.45) +
        # --- Posterior (top layer, slightly lower alpha so both visible) ---
        ggplot2::geom_ribbon(data = posterior_df,
                             ggplot2::aes(ymin = 0, ymax = density,
                                          fill = Distribution),
                             alpha = 0.40, colour = NA) +
        ggplot2::geom_line(data = posterior_df,
                           ggplot2::aes(y = density, colour = Distribution),
                           linewidth = 0.45) +
        # Small-multiples grid: facet_wrap frees BOTH axes per panel (row-only
        # facet_grid would share one x-axis and crush small-scale parameters),
        # and a multi-column layout keeps each panel tall enough to read.
        ggplot2::facet_wrap(~ Parameter, ncol = 4L, scales = "free") +
        ggplot2::scale_x_continuous(
          n.breaks = 3L,
          labels   = function(x) formatC(x, format = "g", digits = 2)
        ) +
        ggplot2::scale_fill_manual(
          values = c("Posterior" = dynhr_colours$dark_blue,
                     "Prior"     = dynhr_colours$orange), name = NULL) +
        ggplot2::scale_colour_manual(
          values = c("Posterior" = dynhr_colours$dark_blue,
                     "Prior"     = dynhr_colours$orange), name = NULL) +
        theme_dynhr_compact() +
        ggplot2::labs(
          title    = "D6: Prior vs Posterior (ridge density)",
          subtitle = flag_label,
          x        = "Parameter value"
        )

      p <- .apply_meta(p, meta)
      attr(p, "dynhr_fig_height") <- max(5, ceiling(n_par / 4) * 1.7)
      plots$ridge_prior_posterior <- p
    }

    pass            <- length(uninformative) == 0
    n_uninformative <- length(uninformative)

    .make_result(
      result  = list(overlap_scores = overlap_scores, uninformative = uninformative),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "D6 Posterior vs prior: %d params. Overlap range [%.2f, %.2f]. %s",
        n_par, min(overlap_scores), max(overlap_scores),
        ifelse(pass, "PASS -- all posteriors updated by data.",
               sprintf("FAIL -- %d uninformative: %s",
                       length(uninformative), paste(uninformative, collapse = ", ")))
      ),
      llm_summary = {
        n_par   <- length(overlap_scores)
        badge   <- if (pass) "PASS" else "FAIL"
        worst5_uninf <- head(sort(overlap_scores, decreasing = TRUE), 5)
        most_inf     <- head(sort(overlap_scores), 5)
        paste(c(
          sprintf("D6 | Prior vs Posterior | %s", badge),
          sprintf("  params=%d uninformative_gt0.80=%d/%d",
                  n_par, n_uninformative, n_par),
          sprintf("  overlap_min=%.2f (%s) overlap_max=%.2f (%s) overlap_median=%.2f",
                  min(overlap_scores), names(which.min(overlap_scores)),
                  max(overlap_scores), names(which.max(overlap_scores)),
                  median(overlap_scores)),
          if (length(most_inf) > 0)
            sprintf("  most_informative: %s",
                    paste(sprintf("%s=%.2f", names(most_inf), most_inf),
                          collapse = ", ")),
          if (n_uninformative > 0)
            sprintf("  uninformative: %s",
                    paste(sprintf("%s=%.2f", names(worst5_uninf), worst5_uninf),
                          collapse = ", ")),
          sprintf("  action: %s",
                  if (n_uninformative == 0)
                    "All parameters updated by data. Prior specification looks reasonable."
                  else
                    sprintf("%s: data barely moves prior (overlap>0.80). Check identification or widen/narrow prior.",
                            paste(names(head(sort(overlap_scores, decreasing = TRUE),
                                              min(3, n_uninformative))), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}
