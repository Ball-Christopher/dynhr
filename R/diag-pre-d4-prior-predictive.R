## R/diag-pre-d4-prior-predictive.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D4 prior predictive checks
## --------------------------------------------------------------------------

#' D4. Prior predictive checks (Fernandez-Villaverde & Guerron-Quintana 2020)
#'
#' Draws parameter vectors from the joint prior, evaluates the model's summary
#' statistics at each draw, and locates each sample (data) statistic inside the
#' prior-predictive distribution of that statistic.
#'
#' For moment \eqn{j} with data value \eqn{d_j} and prior-predictive draws
#' \eqn{m_j^{(i)}}, the reported tail probability is the mid-rank
#' \deqn{p_j = \Pr(m_j < d_j) + \tfrac12 \Pr(m_j = d_j),}
#' so a statistic that equals every draw (a constant moment) gets
#' \eqn{p_j = 0.5}, not 1. A moment is "in the 95\% interval" when
#' \eqn{0.025 \le p_j \le 0.975}; it is "extreme" when the data lie strictly
#' above (\eqn{p_j = 1}) or strictly below (\eqn{p_j = 0}) every draw. The
#' badge PASSes when at least half of the moments are in the interval.
#'
#' Moments are matched BY NAME: when both \code{data_moments} and the output
#' of \code{model_solve_fn} are named, the model output is subset and
#' reordered to \code{names(data_moments)} (so \code{model_solve_fn} may
#' return more moments than the data carry). Unnamed output is matched by
#' position and must have the same length as \code{data_moments}.
#'
#' @param model_solve_fn  Function: theta -> (named) numeric vector of moments.
#' @param prior_draw_fn   Function returning one draw from the joint prior as a
#'   NAMED numeric vector. A function with no formals is called as
#'   \code{prior_draw_fn()}; otherwise as \code{prior_draw_fn(1L)} (matrix
#'   sampler style, first row used).
#' @param data_moments    Named numeric vector of sample moments from the data.
#'   Non-finite entries are dropped (and reported).
#' @param moment_names    Character vector used only when \code{data_moments}
#'   is unnamed and has the same length.
#' @param n_draws         Number of prior draws (default 500).
#' @param meta            Optional plot metadata (see \code{.apply_meta}).
#' @return dynhr_diagnostic list with a density plot.
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

    if (is.list(data_moments) || !is.numeric(data_moments) ||
        length(data_moments) == 0L)
      .dynhr_abort("D4: `data_moments` must be a non-empty numeric vector ",
                   "(got ", class(data_moments)[1L], ").")
    data_moments <- stats::setNames(as.numeric(data_moments), names(data_moments))
    if (is.null(names(data_moments))) {
      names(data_moments) <- if (!is.null(moment_names) &&
                                 length(moment_names) == length(data_moments))
        as.character(moment_names) else paste0("m_", seq_along(data_moments))
    }
    # Non-finite data moments carry no information and would turn the badge
    # into NA; drop them and say so.
    dropped_nf <- names(data_moments)[!is.finite(data_moments)]
    data_moments <- data_moments[is.finite(data_moments)]
    if (length(data_moments) == 0L)
      .dynhr_abort("D4: every entry of `data_moments` is non-finite.")
    moment_names <- names(data_moments)
    n_mom <- length(data_moments)

    prior_moments <- matrix(NA_real_, nrow = n_draws, ncol = n_mom,
                            dimnames = list(NULL, moment_names))
    n_success  <- 0L
    n_skip     <- 0L
    first_fail <- NULL

    # Call convention fixed once from the signature. The old code retried
    # fn() after any failed fn(1L), which hid genuine sampler errors. A prior
    # sampler that errors is broken, not unlucky: its error propagates.
    if (!is.function(prior_draw_fn))
      .dynhr_abort("D4: `prior_draw_fn` must be a function.")
    no_arg <- !is.primitive(prior_draw_fn) && length(formals(prior_draw_fn)) == 0L
    .draw <- function() if (no_arg) prior_draw_fn() else prior_draw_fn(1L)

    for (i in seq_len(n_draws)) {
      raw_draw <- .draw()
      if (is.null(raw_draw)) {
        if (is.null(first_fail)) first_fail <- "prior_draw_fn returned NULL"
        n_skip <- n_skip + 1L; next
      }
      # Preserve parameter names: model_solve_fn maps theta BY NAME, so an
      # unnamed draw silently reuses the baseline calibration for every draw.
      if (is.matrix(raw_draw) || is.data.frame(raw_draw)) {
        nms <- colnames(raw_draw)
        raw_draw <- unlist(raw_draw[1L, , drop = TRUE])
      } else {
        nms <- names(raw_draw)
      }
      theta_draw <- suppressWarnings(stats::setNames(as.numeric(raw_draw), nms))
      if (length(theta_draw) == 0L || !all(is.finite(theta_draw))) {
        if (is.null(first_fail)) first_fail <- "prior draw empty or non-finite"
        n_skip <- n_skip + 1L; next
      }

      # A solve that legitimately fails at a prior draw (BK violation, no
      # steady state) must RETURN non-finite moments; those draws are counted
      # in n_skip below.  An error propagates (no tryCatch by convention).
      mom <- model_solve_fn(theta_draw)
      if (!is.numeric(mom)) {
        if (is.null(first_fail)) first_fail <- "model_solve_fn returned a non-numeric value"
        n_skip <- n_skip + 1L; next
      }
      # Align BY NAME when possible: positional matching silently pairs the
      # wrong moments when the orders differ.
      mnm <- names(mom)
      if (!is.null(mnm) && !anyNA(mnm) && all(nzchar(mnm))) {
        miss <- setdiff(moment_names, mnm)
        if (length(miss)) {
          if (is.null(first_fail)) first_fail <- sprintf("model_solve_fn output lacks data moment(s): %s",
                             paste(utils::head(miss, 5L), collapse = ", "))
          n_skip <- n_skip + 1L; next
        }
        mom <- mom[moment_names]
      } else if (length(mom) != n_mom) {
        if (is.null(first_fail)) first_fail <- sprintf(
          "unnamed model_solve_fn output has length %d, data_moments has %d",
          length(mom), n_mom)
        n_skip <- n_skip + 1L; next
      }
      if (!all(is.finite(mom))) {
        if (is.null(first_fail)) first_fail <- "model_solve_fn returned non-finite moments"
        n_skip <- n_skip + 1L; next
      }

      n_success <- n_success + 1L
      prior_moments[n_success, ] <- as.numeric(mom)
    }
    prior_moments <- prior_moments[seq_len(n_success), , drop = FALSE]

    drop_note <- if (length(dropped_nf))
      sprintf(" Dropped non-finite data moment(s): %s.",
              paste(dropped_nf, collapse = ", "))
    else ""
    fail_note <- if (!is.null(first_fail))
      sprintf(" First failure: %s.", first_fail) else ""

    if (n_success < 10L) {
      return(.make_result(
        pass    = NA,
        plots   = list(),
        summary = sprintf(
          "D4 Prior predictive: insufficient valid prior draws (%d / %d valid, %d skipped).%s%s",
          n_success, n_draws, n_skip, fail_note, drop_note)
      ))
    }

    # Degeneracy guard: if no moment varies across draws, the draws never
    # reach the moments -- p-values and densities would be meaningless.
    mom_sd    <- apply(prior_moments, 2, stats::sd)
    mom_scale <- pmax(1, abs(colMeans(prior_moments)))
    const_mom <- mom_sd <= 1e-10 * mom_scale
    if (all(const_mom)) {
      return(.make_result(
        pass    = NA,
        plots   = list(),
        summary = sprintf(
          "D4 Prior predictive: DEGENERATE -- all %d prior draws produce identical moments. The drawn parameters never reach model_solve_fn's output (unnamed draws, or a solve function that reuses a fixed decision rule).%s",
          n_success, drop_note)
      ))
    }

    # Mid-rank tail probability (ties count half).
    p_below <- colMeans(sweep(prior_moments, 2, data_moments, `<`))
    p_equal <- colMeans(sweep(prior_moments, 2, data_moments, `==`))
    p_vals  <- p_below + 0.5 * p_equal
    names(p_vals) <- moment_names
    above_all <- p_below == 1                # data strictly above every draw
    below_all <- (p_below + p_equal) == 0    # data strictly below every draw

    # Pass if at least 50% of moments fall within [0.025, 0.975] of the
    # prior-predictive distribution (majority rule: "a few moments extreme" is
    # a calibration hint, "most moments extreme" is a prior failure).
    in_interval <- (p_vals >= 0.025) & (p_vals <= 0.975)
    pass        <- mean(in_interval) >= 0.5
    p_extreme   <- sum(above_all | below_all)
    status <- ifelse(in_interval, "inside",
                     ifelse(above_all, "above all draws",
                            ifelse(below_all, "below all draws", "outside 95%")))
    names(status) <- moment_names
    const_names <- moment_names[const_mom]
    const_note <- if (length(const_names))
      sprintf(" Constant across draws (the prior does not move them): %s.",
              paste(const_names, collapse = ", "))
    else ""

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      strip_lab <- sprintf("%s  (p = %.2f)", moment_names, p_vals)
      mfac <- function(x) factor(x, levels = moment_names, labels = strip_lab)
      vary <- moment_names[!const_mom]
      pp_long <- data.frame(
        Moment = rep(moment_names, each = n_success),
        Value  = as.vector(prior_moments), stringsAsFactors = FALSE)
      pp_dens <- pp_long[pp_long$Moment %in% vary, , drop = FALSE]
      pp_dens$Moment <- mfac(pp_dens$Moment)
      q <- apply(prior_moments, 2, stats::quantile, probs = c(0.025, 0.975),
                 names = FALSE)
      band_df <- data.frame(Moment = mfac(moment_names), lo = q[1L, ], hi = q[2L, ])
      band_df <- band_df[!const_mom, , drop = FALSE]
      lev_in  <- "Data: inside 95% interval"
      lev_out <- "Data: outside 95% interval"
      data_df <- data.frame(
        Moment = mfac(moment_names), Value = as.numeric(data_moments),
        Status = factor(ifelse(in_interval, lev_in, lev_out),
                        levels = c(lev_in, lev_out)))
      const_df <- data.frame(Moment = mfac(const_names),
                             Value  = unname(colMeans(prior_moments)[const_mom]))

      p_pp <- ggplot2::ggplot() +
        ggplot2::geom_density(data = pp_dens, ggplot2::aes(x = .data$Value),
                              fill = dynhr_colours$light_blue,
                              colour = dynhr_colours$mid_blue,
                              alpha = 0.4, linewidth = 0.5) +
        ggplot2::geom_vline(data = band_df, ggplot2::aes(xintercept = .data$lo),
                            colour = dynhr_colours$mid_blue, linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_vline(data = band_df, ggplot2::aes(xintercept = .data$hi),
                            colour = dynhr_colours$mid_blue, linetype = "dashed",
                            linewidth = 0.4)
      if (nrow(const_df)) {
        # A constant moment has no density; mark the point mass instead of
        # leaving an empty panel.
        p_pp <- p_pp +
          ggplot2::geom_vline(data = const_df,
                              ggplot2::aes(xintercept = .data$Value),
                              colour = dynhr_colours$mid_blue, linewidth = 2.2,
                              alpha = 0.35) +
          ggplot2::geom_text(data = const_df,
                             ggplot2::aes(x = .data$Value, y = 0),
                             label = "all draws identical", vjust = -0.6,
                             hjust = -0.05, size = 3, colour = dynhr_colours$grey)
      }
      p_pp <- p_pp +
        ggplot2::geom_vline(data = data_df,
                            ggplot2::aes(xintercept = .data$Value,
                                         colour = .data$Status),
                            linewidth = 0.8) +
        ggplot2::scale_colour_manual(
          values = stats::setNames(c(dynhr_colours$teal, dynhr_colours$red),
                                   c(lev_in, lev_out)),
          drop = FALSE, name = NULL) +
        ggplot2::facet_wrap(~ Moment, scales = "free", ncol = 3) +
        ggplot2::scale_x_continuous(
          n.breaks = 4L, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        theme_dynhr_compact() +
        ggplot2::theme(panel.spacing.y = ggplot2::unit(1.0, "lines"),
                       panel.spacing.x = ggplot2::unit(1.2, "lines"),
                       legend.position = "bottom",
                       # plain family: display-font ligatures garble strip labels
                       strip.text = ggplot2::element_text(family = "sans")) +
        ggplot2::labs(
          title = "D4: Prior predictive checks",
          subtitle = sprintf(paste0(
            "%d / %d valid prior draws. Shaded: prior-predictive density; ",
            "dashed: its central 95%% interval.\n%d / %d data moments inside ",
            "(PASS needs >= 50%%). p = P(draw < data), ties counted half."),
            n_success, n_draws, sum(in_interval), n_mom),
          x = "Moment value (each panel on its own scale)", y = NULL)
      attr(p_pp, "dynhr_fig_height") <- max(7, ceiling(n_mom / 3) * 2.0 + 1)
      plots$prior_predictive <- .apply_meta(p_pp, meta)
    }

    flagged <- moment_names[!in_interval]
    flag_txt <- paste(sprintf("%s (%s)", flagged, status[flagged]), collapse = ", ")
    .make_result(
      result  = list(prior_moments = prior_moments, p_values = p_vals,
                     data_moments = data_moments,
                     n_success = n_success, n_draws = n_draws, n_skip = n_skip,
                     p_extreme = p_extreme, in_interval = in_interval,
                     status = status, constant_moments = const_names,
                     dropped_moments = dropped_nf),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "D4 Prior predictive: %d/%d valid draws (%d skipped). Data in 95%% interval: %d/%d moments (pass needs >= 50%%). Extreme (outside every draw): %d. %s%s%s%s",
        n_success, n_draws, n_skip, sum(in_interval), n_mom, p_extreme,
        if (pass && !length(flagged)) "PASS."
        else if (pass) sprintf("PASS; flagged: %s.", flag_txt)
        else sprintf("FAIL -- outside: %s.", flag_txt),
        const_note, drop_note, fail_note),
      llm_summary = {
        badge <- if (pass) "PASS" else "FAIL"
        pval_parts <- sprintf("%s=%.3f%s", moment_names, p_vals,
                              ifelse(in_interval, "", sprintf("[%s]", status)))
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
        paste(c(
          sprintf("D4 | Prior Predictive Checks | %s", badge),
          sprintf("  draws=%d valid=%d skipped=%d moments=%d in_interval=%d/%d extreme=%d",
                  n_draws, n_success, n_skip, n_mom, sum(in_interval), n_mom, p_extreme),
          "  gate: pass = mean(0.025 <= p <= 0.975) >= 0.50; p = P(draw < data) + P(draw = data)/2",
          "  p_values:",
          paste0("    ", .wrap_pvals(pval_parts)),
          if (length(const_names))
            sprintf("  constant across draws: %s", paste(const_names, collapse = ", ")),
          if (length(dropped_nf))
            sprintf("  dropped non-finite data moments: %s",
                    paste(dropped_nf, collapse = ", ")),
          sprintf("  action: %s",
                  if (pass && !length(flagged))
                    "Prior predictive encompasses every observed moment."
                  else if (pass)
                    sprintf("Prior predictive encompasses the majority of observed moments; check %s.",
                            paste(utils::head(flagged, 5L), collapse = ", "))
                  else if (p_extreme == n_mom)
                    "Every moment lies outside all prior draws: prior is grossly miscalibrated. Reset prior location and scale."
                  else
                    sprintf("%d/%d moments outside the 95%% interval (%d extreme): %s. Revisit the prior or the model parameterisation.",
                            length(flagged), n_mom, p_extreme,
                            paste(utils::head(flagged, 5L), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}
