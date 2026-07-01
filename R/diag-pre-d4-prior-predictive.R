## R/diag-pre-d4-prior-predictive.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D4 prior predictive checks
## --------------------------------------------------------------------------

#' D4. Prior predictive checks (Fernandez-Villaverde & Guerron-Quintana 2020)
#'
#' Draws parameter vectors from the joint prior, solves the model at each
#' draw, collects summary statistics, and compares the prior-predictive
#' distribution of those statistics to data sample values.
#'
#' @param model_solve_fn  Function: theta -> named numeric vector of moments
#' @param prior_draw_fn   Function: () -> theta (single draw from joint prior)
#' @param data_moments    Named numeric vector of sample moments from the data
#' @param moment_names    Character vector (optional, inferred from data_moments)
#' @param n_draws         Number of prior draws (default 500)
#' @return dynhr_diagnostic list with density plots
#' @references Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive
#'   assessment of model fitness via realized discrepancies.
#'   \emph{Statistica Sinica}, 6(4), 733-807.
#'   Fernandez-Villaverde, J., & Guerron-Quintana, P. A. (2020). Estimating DSGE
#'   models: Recent advances and future challenges.
#'   \emph{Annual Review of Economics}, 13, 229-255.
#' @noRd
d4_prior_predictive <- function(model_solve_fn,
                                prior_draw_fn,
                                data_moments,
                                moment_names = NULL,
                                n_draws = 500,
                                meta = NULL) {

    if (is.null(moment_names)) moment_names <- names(data_moments)
    n_mom <- length(data_moments)
    if (is.null(moment_names)) moment_names <- paste0("m_", seq_len(n_mom))

    # Draw from prior and collect moments
    prior_moments <- matrix(NA_real_, nrow = n_draws, ncol = n_mom)
    # Guard: dimnames assignment deferred until after row/col counts are confirmed
    n_success <- 0
    n_skip    <- 0L

    # Detect if prior_draw_fn accepts an argument n (matrix sampler style).
    # Try calling with n=1; if that fails, fall back to no-argument call.
    .call_prior_fn <- function(fn) {
      tryCatch(fn(1L), error = function(e) tryCatch(fn(), error = function(e2) NULL))
    }

    for (i in seq_len(n_draws)) {
      # --- Normalise prior_draw_fn() return shape explicitly ---
      # Accept both scalar-theta samplers () -> named numeric and
      # matrix samplers (n) -> matrix (n x n_par); take first row.
      raw_draw <- .call_prior_fn(prior_draw_fn)
      # Preserve parameter names: model_solve_fn maps theta BY NAME, so an
      # unnamed draw silently reuses the baseline calibration for every draw
      # (degenerate prior predictive).
      nms <- if (is.matrix(raw_draw) || is.data.frame(raw_draw)) {
        colnames(raw_draw)
      } else {
        names(raw_draw)
      }
      if (is.matrix(raw_draw) || is.data.frame(raw_draw)) {
        raw_draw <- raw_draw[1L, ]
      }
      theta_draw <- stats::setNames(as.numeric(raw_draw), nms)
      if (!is.numeric(theta_draw) || !all(is.finite(theta_draw)) || length(theta_draw) == 0L) {
        n_skip <- n_skip + 1L
        next
      }

      # --- Validate model output shape explicitly ---
      mom <- tryCatch(model_solve_fn(theta_draw), error = function(e) NULL)
      if (is.null(mom) || !is.numeric(mom) || length(mom) != n_mom || !all(is.finite(mom))) {
        n_skip <- n_skip + 1L
        next
      }

      n_success <- n_success + 1L
      prior_moments[n_success, ] <- as.numeric(mom)
    }

    # Trim to successful draws; guard the dimnames assignment against mismatch
    prior_moments <- prior_moments[seq_len(n_success), , drop = FALSE]
    if (nrow(prior_moments) == n_success && ncol(prior_moments) == n_mom) {
      colnames(prior_moments) <- moment_names
    }

    if (n_success < 10L) {
      return(.make_result(
        pass    = NA,
        plots   = list(),
        summary = sprintf(
          "D4 Prior predictive: insufficient valid prior draws (%d / %d valid, %d skipped due to non-conforming draw shape or non-finite model output). Prior may be badly calibrated or prior_draw_fn returns unexpected shapes.",
          n_success, n_draws, n_skip)
      ))
    }

    # Degeneracy guard: if the prior-predictive moments barely vary across
    # draws, the sampler is effectively broken (e.g. unnamed theta ignored by
    # model_solve_fn) -- p-values and densities would be meaningless.
    mom_sd <- apply(prior_moments, 2, stats::sd)
    if (all(mom_sd < 1e-12)) {
      return(.make_result(
        pass    = NA,
        plots   = list(),
        summary = sprintf(
          "D4 Prior predictive: DEGENERATE -- all %d prior draws produce identical moments. prior_draw_fn/model_solve_fn wiring is broken (check that draws are NAMED parameter vectors).",
          n_success)
      ))
    }

    # Check: what fraction of prior-predictive distribution is <= data value?
    p_vals <- numeric(n_mom)
    for (j in seq_len(n_mom)) {
      p_vals[j] <- mean(prior_moments[, j] <= data_moments[j])
    }
    names(p_vals) <- moment_names

    # Pass if at least 50% of moments fall within [0.025, 0.975] of
    # prior-predictive distribution.
    # The old all-or-nothing criterion (pass = all(in_interval)) was a
    # constant-fail for miscalibrated priors: any single outlier moment caused
    # FAIL with no nuance.  The 50% majority rule distinguishes "a few moments
    # extreme" (calibration hint) from "all moments extreme" (prior failure).
    # p_extreme counts p in {0, 1} -- moments where ALL prior draws are above
    # (p=0) or below (p=1) the data; these are the most egregious failures.
    in_interval <- (p_vals >= 0.025) & (p_vals <= 0.975)
    pass        <- mean(in_interval) >= 0.5
    p_extreme   <- sum(p_vals == 0 | p_vals == 1)

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # Density plots for each moment
    pp_long <- reshape2::melt(as.data.frame(prior_moments), variable.name = "Moment",
                              value.name = "Value")
    data_df <- data.frame(Moment = factor(moment_names, levels = moment_names),
                          Value  = data_moments)

    p_pp <- ggplot2::ggplot(pp_long, ggplot2::aes(x = Value)) +
      ggplot2::geom_density(fill = dynhr_colours$light_blue, colour = dynhr_colours$mid_blue,
                   alpha = 0.4, linewidth = 0.5) +
      ggplot2::geom_vline(data = data_df, ggplot2::aes(xintercept = Value),
                 colour = dynhr_colours$red, linewidth = 0.7, linetype = "solid") +
      ggplot2::facet_wrap(~ Moment, scales = "free", ncol = 3) +
      ggplot2::scale_x_continuous(n.breaks = 3L) +
      ggplot2::scale_y_continuous(n.breaks = 3L) +
      theme_dynhr_compact() +
      ggplot2::theme(panel.spacing.y = ggplot2::unit(1.0, "lines"),
                     # plain family: display-font ligatures garble strip labels
                     strip.text = ggplot2::element_text(family = "sans")) +
      ggplot2::labs(title = "D4: Prior predictive checks",
           subtitle = sprintf("%d / %d valid draws (%d skipped) | Red line = sample data value",
                              n_success, n_draws, n_skip),
           x = "Moment value", y = NULL)
    # Taller canvas so the per-moment panels are not vertically compressed.
    attr(p_pp, "dynhr_fig_height") <- max(7, ceiling(n_mom / 3) * 2.0)
    plots$prior_predictive <- .apply_meta(p_pp, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(prior_moments = prior_moments, p_values = p_vals,
                     n_success = n_success, n_draws = n_draws,
                     p_extreme = p_extreme, in_interval = in_interval),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "D4 Prior predictive: %d/%d valid draws (%d skipped). Data in 95%% interval: %d/%d moments (pass=>=50%%). p_extreme=%d. %s",
        n_success, n_draws, n_skip, sum(in_interval), n_mom, p_extreme,
        if (pass) "PASS."
        else sprintf("FAIL -- %d/%d moments outside 95%% interval. Outliers: %s",
                     n_mom - sum(in_interval), n_mom,
                     paste(moment_names[!in_interval], collapse = ", "))
      ),
      llm_summary = {
        badge <- if (pass) "PASS" else "FAIL"
        # Annotate direction for extreme p-values and wrap to <=100 chars/line
        pval_parts <- sprintf("%s=%.2f%s",
                              moment_names, p_vals,
                              ifelse(in_interval, "",
                                     ifelse(p_vals == 1, "[above all draws]",
                                            ifelse(p_vals == 0, "[below all draws]",
                                                   "[FLAG]"))))
        # Break the p_values line at ~90 chars
        .wrap_pvals <- function(parts, width = 90) {
          lines <- character(0); cur <- ""
          for (p in parts) {
            candidate <- if (nchar(cur) == 0) p else paste0(cur, ", ", p)
            if (nchar(candidate) > width && nchar(cur) > 0) {
              lines <- c(lines, cur); cur <- p
            } else {
              cur <- candidate
            }
          }
          if (nchar(cur) > 0) lines <- c(lines, cur)
          lines
        }
        pval_lines <- .wrap_pvals(pval_parts)
        paste(c(
          sprintf("D4 | Prior Predictive Checks | %s", badge),
          sprintf("  draws=%d valid=%d moments=%d in_interval=%d/%d p_extreme=%d",
                  n_draws, n_success, n_mom, sum(in_interval), n_mom, p_extreme),
          sprintf("  gate: pass = mean(in_95pct) >= 0.50 (majority criterion)"),
          sprintf("  p_values:"),
          paste0("    ", pval_lines),
          sprintf("  action: %s",
                  if (pass)
                    "Prior predictive encompasses majority of observed moments. Check flagged moments."
                  else if (p_extreme == n_mom)
                    "All moments are extreme (p in {0,1}): prior is grossly miscalibrated. Reset prior location and scale."
                  else {
                    flagged <- moment_names[!in_interval]
                    sprintf("%d/%d moments outside 95%% interval (%d extreme). Moments: %s. Tighten prior or check model parameterisation.",
                            length(flagged), n_mom, p_extreme,
                            paste(head(flagged, 5), collapse = ", "))
                  })
        ), collapse = "\n")
      }
    )
}
