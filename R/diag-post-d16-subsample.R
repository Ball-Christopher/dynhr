## R/diag-post-d16-subsample.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D16 subsample stability (forest plot + standardised-shift plot)
## --------------------------------------------------------------------------

#' D16. Subsample stability
#'
#' Compares the posterior of each parameter on one or more subsamples with the
#' full-sample posterior.
#'
#' **Statistic.** The subsample data are a SUBSET of the full-sample data, so
#' the two posteriors are not independent: the full-sample posterior mean
#' already contains the subsample's information. Under a stable, correctly
#' specified model (and Bernstein-von Mises asymptotics) the full-sample mean
#' is the efficient estimator and the subsample mean an inefficient one, so --
#' exactly as in Hausman (1978) -- the sampling variance of their difference is
#' the DIFFERENCE of the posterior variances:
#' \deqn{z_{ij} = (\bar\theta^{sub_i}_j - \bar\theta^{full}_j) /
#'   \sqrt{\max(V^{sub_i}_j - V^{full}_j, 0) + V^{sub_i}_j/ESS^{sub_i}_j +
#'   V^{full}_j/ESS^{full}_j}}
#' which is approximately N(0, 1) under stability. The last two terms are the
#' Monte Carlo variance of the two posterior means (ESS from the draws). For a
#' subsample holding a fraction f of the data, the difference has standard
#' deviation \eqn{\sigma_{full}\sqrt{(1-f)/f}}, so a threshold stated in
#' full-sample posterior SDs (the pre-0.9.4 gate) false-alarms heavily for
#' short subsamples (f = 0.25: 25\% per parameter).
#'
#' **Gate.** pass = every \eqn{|z_{ij}| \le} \code{z_tol} AND every subsample
#' credible interval overlaps the full-sample one. By default \code{z_tol} is
#' the Bonferroni critical value for a family-wise false-alarm rate of
#' \code{family_alpha} over all (subsample, parameter) pairs. The CI-overlap
#' check alone is not a usable test (wide posteriors overlap trivially and,
#' for nested data, it essentially never fires); it is reported for reference.
#'
#' @section Precondition: the subsamples must be NESTED:
#' Hausman's (1978) variance-of-the-difference formula
#' \eqn{\mathrm{Var}(\hat\theta_{ineff} - \hat\theta_{eff}) =
#' \mathrm{Var}(\hat\theta_{ineff}) - \mathrm{Var}(\hat\theta_{eff})} holds
#' only when one estimator is a strict efficiency improvement on the other --
#' here, when the subsample's data are a SUBSET of the full sample's. This is a
#' precondition, not a technicality: for a non-nested split (two disjoint
#' periods, a rolling window compared with another window) the two posteriors
#' are independent, the variance of the difference is their SUM, and
#' subtracting instead of adding can drive the denominator to zero or negative
#' -- the statistic is then not merely conservative, it is undefined. D16 does
#' not police this, because draw matrices carry no sample ranges: \code{
#' results_sub} is a list of draws and nothing in it says which periods
#' produced them. What D16 CAN detect it does detect -- a subsample posterior
#' markedly TIGHTER than the full-sample one (variance < 0.9x) is impossible
#' under nesting and is flagged (\code{sub_tighter}); its z then carries only
#' the Monte Carlo variance and is unreliable.
#'
#' For a non-nested or rolling split, use a structural-break test built for
#' unknown break points instead -- the sup-Wald / sup-LM family of Andrews
#' (1993) -- not this statistic.
#'
#' @param results_full  Matrix or data frame (n_draws x n_params) -- full sample
#' @param results_sub   Named list of matrices/data frames -- subsample draws,
#'   e.g. \code{list("Pre-GFC" = draws1, "Post-GFC" = draws2)}. If both the
#'   full-sample and a subsample matrix carry column names, columns are matched
#'   BY NAME; otherwise by position.
#' @param param_names   Character vector (optional; defaults to the column
#'   names of \code{results_full})
#' @param split_labels  Character vector -- labels for the subsamples
#' @param ci_level      Credible interval level (default 0.90)
#' @param z_tol         Threshold on |z| (see Details). NULL (default) = the
#'   Bonferroni value \code{qnorm(1 - family_alpha / (2 * n_sub * n_par))}.
#' @param family_alpha  Family-wise false-alarm rate used when \code{z_tol} is
#'   NULL (default 0.05).
#' @param meta          Optional \code{dynhr_diag_meta} for plot captions.
#' @return dynhr_diagnostic list with forest and z-statistic plots. The
#'   \code{result$stability_df} table has one row per (subsample, parameter).
#' @references Hausman, J. A. (1978). Specification tests in econometrics.
#'   \emph{Econometrica}, 46(6), 1251-1271.
#'
#'   Andrews, D. W. K. (1993). Tests for parameter instability and structural
#'   change with unknown change point. \emph{Econometrica}, 61(4), 821-856.
#'   (The appropriate alternative for NON-nested / rolling subsample splits,
#'   where Hausman's variance-difference form does not apply.)
#'
#'   Lubik, T. A., & Schorfheide, F. (2004). Testing for indeterminacy:
#'   An application to US monetary policy. \emph{American Economic Review},
#'   94(1), 190-217.
#' @noRd
d16_subsample_stability <- function(results_full,
                                    results_sub,
                                    param_names  = NULL,
                                    split_labels = NULL,
                                    ci_level     = 0.90,
                                    z_tol        = NULL,
                                    family_alpha = 0.05,
                                    meta         = NULL) {

    results_full <- as.matrix(results_full)
    n_par <- ncol(results_full)
    if (!is.list(results_sub) || length(results_sub) == 0L)
      .dynhr_abort("D16: `results_sub` must be a non-empty list of draw matrices.")
    n_sub <- length(results_sub)
    if (is.null(param_names)) {
      param_names <- if (!is.null(colnames(results_full))) colnames(results_full)
      else paste0("theta_", seq_len(n_par))
    }
    if (length(param_names) != n_par)
      .dynhr_abort(sprintf("D16: %d param_names for %d full-sample columns.",
                           length(param_names), n_par))

    if (is.null(split_labels)) split_labels <- names(results_sub)
    if (is.null(split_labels)) split_labels <- paste0("Sub_", seq_len(n_sub))
    if (length(split_labels) != n_sub)
      .dynhr_abort(sprintf("D16: %d split_labels for %d subsamples.",
                           length(split_labels), n_sub))
    empty_lab <- is.na(split_labels) | !nzchar(split_labels)
    split_labels[empty_lab] <- paste0("Sub_", which(empty_lab))
    # Labels key the plot palette/factor levels: keep them unique and distinct
    # from the "Full sample" row.
    split_labels <- make.unique(c("Full sample", as.character(split_labels)),
                                sep = "_")[-1L]

    # Align every subsample's columns with the full sample (by name when both
    # have names, otherwise by position).
    full_cn <- colnames(results_full)
    results_sub <- lapply(seq_len(n_sub), function(i) {
      m <- as.matrix(results_sub[[i]])
      if (!is.null(full_cn) && !is.null(colnames(m))) {
        miss <- setdiff(full_cn, colnames(m))
        if (length(miss) > 0L)
          .dynhr_abort(sprintf("D16: subsample '%s' lacks column(s): %s",
                               split_labels[i], paste(miss, collapse = ", ")))
        m <- m[, full_cn, drop = FALSE]
      } else if (ncol(m) != n_par) {
        .dynhr_abort(sprintf("D16: subsample '%s' has %d columns, full sample %d.",
                             split_labels[i], ncol(m), n_par))
      }
      m
    })
    all_draws <- c(list(results_full), results_sub)
    for (m in all_draws) {
      if (nrow(m) < 2L || any(!is.finite(m)))
        .dynhr_abort("D16: every draw matrix needs >= 2 rows and finite values.")
    }

    alpha <- (1 - ci_level) / 2
    n_tests <- n_sub * n_par
    if (is.null(z_tol)) z_tol <- stats::qnorm(1 - family_alpha / (2 * n_tests))

    .summarise_draws <- function(draws, label) {
      data.frame(
        Parameter = param_names,
        Median    = apply(draws, 2, stats::median),
        Lo        = apply(draws, 2, stats::quantile, probs = alpha, names = FALSE),
        Hi        = apply(draws, 2, stats::quantile, probs = 1 - alpha, names = FALSE),
        Sample    = label,
        stringsAsFactors = FALSE
      )
    }
    # Per-column ESS (capped at the number of draws); falls back to n for a
    # degenerate (constant) column, whose variance is zero anyway.
    .col_ess <- function(draws) {
      vapply(seq_len(ncol(draws)), function(j) {
        e <- .d5_ess_basic(draws[, j, drop = FALSE])
        if (is.finite(e) && e > 0) min(e, nrow(draws)) else nrow(draws)
      }, numeric(1))
    }

    summary_list <- c(list(.summarise_draws(results_full, "Full sample")),
                      lapply(seq_len(n_sub), function(i)
                        .summarise_draws(results_sub[[i]], split_labels[i])))
    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL

    full_summary <- summary_list[[1]]
    full_mean <- colMeans(results_full)
    full_var  <- apply(results_full, 2, stats::var)
    full_ess  <- .col_ess(results_full)

    stab_rows <- vector("list", n_sub)
    for (i in seq_len(n_sub)) {
      sub_draws <- results_sub[[i]]
      sub_summary <- summary_list[[i + 1]]
      sub_mean <- colMeans(sub_draws)
      sub_var  <- apply(sub_draws, 2, stats::var)
      sub_ess  <- .col_ess(sub_draws)
      mc_var   <- sub_var / sub_ess + full_var / full_ess
      se_diff  <- sqrt(pmax(sub_var - full_var, 0) + mc_var)
      diff     <- sub_mean - full_mean
      z <- ifelse(se_diff > 0, diff / se_diff,
                  ifelse(diff == 0, 0, sign(diff) * Inf))
      overlap <- !(sub_summary$Hi < full_summary$Lo |
                   sub_summary$Lo > full_summary$Hi)
      stab_rows[[i]] <- data.frame(
        Sample      = split_labels[i],
        Parameter   = param_names,
        mean_sub    = sub_mean,
        mean_full   = full_mean,
        sd_sub      = sqrt(sub_var),
        sd_full     = sqrt(full_var),
        se_diff     = se_diff,
        z           = z,
        ci_overlap  = overlap,
        drift_ok    = abs(z) <= z_tol,
        sub_tighter = sub_var < 0.9 * full_var,
        stringsAsFactors = FALSE
      )
    }
    stability_df <- do.call(rbind, stab_rows)
    rownames(stability_df) <- NULL

    by_par <- function(v) vapply(param_names, function(p)
      all(v[stability_df$Parameter == p]), logical(1), USE.NAMES = FALSE)
    overlap_check    <- by_par(stability_df$ci_overlap)
    mean_drift_check <- by_par(stability_df$drift_ok)
    stable_check     <- overlap_check & mean_drift_check
    pass             <- all(stable_check)
    unstable_params  <- param_names[!stable_check]
    ci_fail_params   <- param_names[!overlap_check]
    drift_fail_params <- param_names[!mean_drift_check]
    tighter_rows <- stability_df[stability_df$sub_tighter, , drop = FALSE]
    tighter_txt <- if (nrow(tighter_rows) > 0L)
      sprintf("subsample posterior tighter than full sample (z unreliable, is the subsample nested?): %s",
              paste(sprintf("%s/%s", tighter_rows$Sample, tighter_rows$Parameter),
                    collapse = ", "))
    else NULL
    max_abs_z <- max(abs(stability_df$z))

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    sample_levels <- c("Full sample", split_labels)
    pal <- stats::setNames(
      c(dynhr_primary_colour,
        rep_len(setdiff(dynhr_palette, dynhr_primary_colour), n_sub)),
      sample_levels)
    shapes <- stats::setNames(rep_len(c(16, 17, 15, 18, 3, 4, 8), n_sub + 1L),
                              sample_levels)

    summary_df$Parameter <- factor(summary_df$Parameter, levels = param_names)
    summary_df$Sample    <- factor(summary_df$Sample, levels = rev(sample_levels))
    full_band <- data.frame(
      Parameter = factor(param_names, levels = param_names),
      Lo = full_summary$Lo, Hi = full_summary$Hi, Median = full_summary$Median)

    p_for <- ggplot2::ggplot(summary_df,
                           ggplot2::aes(x = Median, y = Sample,
                               colour = Sample, shape = Sample)) +
      ggplot2::geom_rect(data = full_band,
                         ggplot2::aes(xmin = Lo, xmax = Hi,
                                      ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE,
                         fill = dynhr_na_fill, alpha = 0.6) +
      ggplot2::geom_vline(data = full_band,
                          ggplot2::aes(xintercept = Median),
                          colour = dynhr_na_colour, linetype = "dashed",
                          linewidth = 0.4) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = Lo, xmax = Hi),
                     orientation = "y", width = 0.35, linewidth = 0.5) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~ Parameter, scales = "free_x",
                          ncol = min(4L, ceiling(sqrt(n_par)))) +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_shape_manual(values = shapes) +
      ggplot2::scale_x_continuous(n.breaks = 3L) +
      theme_dynhr_compact() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.75)),
        axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.70)),
        panel.spacing.x = ggplot2::unit(1.2, "lines"),
        legend.position = "none"
      ) +
      ggplot2::labs(title = "D16: Subsample stability -- posterior medians and credible intervals",
           subtitle = sprintf("%d%% CIs | shaded band + dashed line = full-sample CI and median | free x-scale",
                              round(ci_level * 100)),
           x = "Parameter value")
    plots$forest <- .apply_meta(p_for, meta)

    zdf <- stability_df
    zdf$Parameter <- factor(zdf$Parameter, levels = rev(param_names))
    zdf$Sample    <- factor(zdf$Sample, levels = split_labels)
    zdf$Status    <- factor(ifelse(zdf$drift_ok, "within threshold", "beyond threshold"),
                            levels = c("within threshold", "beyond threshold"))
    z_lim <- max(z_tol * 1.15, min(max(abs(zdf$z[is.finite(zdf$z)]), 0) * 1.05, 50))
    zdf$z_plot <- pmax(pmin(zdf$z, z_lim), -z_lim)
    p_z <- ggplot2::ggplot(zdf, ggplot2::aes(x = z_plot, y = Parameter,
                                             colour = Sample, shape = Status)) +
      ggplot2::annotate("rect", xmin = -z_tol, xmax = z_tol,
                        ymin = -Inf, ymax = Inf,
                        fill = dynhr_na_fill, alpha = 0.5) +
      ggplot2::geom_vline(xintercept = c(-z_tol, z_tol),
                          colour = dynhr_colours$red, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, colour = dynhr_na_colour) +
      ggplot2::geom_point(size = 2.4,
                          position = ggplot2::position_dodge(width = 0.5)) +
      ggplot2::scale_colour_manual(values = pal[split_labels], name = "Subsample") +
      ggplot2::guides(colour = ggplot2::guide_legend(order = 1),
                      shape  = ggplot2::guide_legend(order = 2)) +
      ggplot2::scale_shape_manual(values = c("within threshold" = 16,
                                             "beyond threshold" = 4),
                                  name = NULL, drop = TRUE) +
      ggplot2::coord_cartesian(xlim = c(-z_lim, z_lim)) +
      theme_dynhr() +
      ggplot2::labs(
        title = "D16: Standardised subsample shift (Hausman z)",
        subtitle = sprintf(
          "z = (sub mean - full mean) / sqrt(V_sub - V_full + MC var) ~ N(0,1) if stable | dashed = +/-%.2f%s",
          z_tol, if (any(abs(zdf$z) > z_lim)) " | off-scale values clipped" else ""),
        x = "z (standard errors of the difference)", y = NULL)
    plots$zstat <- .apply_meta(p_z, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(summary_df        = summary_df,
                     stability_df      = stability_df,
                     unstable_params   = unstable_params,
                     overlap_check     = overlap_check,
                     mean_drift_check  = mean_drift_check,
                     ci_fail_params    = ci_fail_params,
                     drift_fail_params = drift_fail_params,
                     z_tol             = z_tol,
                     max_abs_z         = max_abs_z),
      pass    = pass,
      plots   = plots,
      summary = paste(c(sprintf(
        "D16 Subsample stability: %d params, %d subsamples, max|z| = %.2f (threshold %.2f). %s",
        n_par, n_sub, max_abs_z, z_tol,
        if (pass) {
          "PASS -- all subsample CIs overlap the full-sample CI and no standardised shift exceeds the threshold."
        } else {
          parts <- character(0)
          if (length(ci_fail_params) > 0)
            parts <- c(parts, sprintf("CI non-overlap: %s", paste(ci_fail_params, collapse = ", ")))
          if (length(drift_fail_params) > 0)
            parts <- c(parts, sprintf("|z| > %.2f: %s", z_tol, paste(drift_fail_params, collapse = ", ")))
          sprintf("FAIL -- %d unstable parameter(s). %s",
                  length(unstable_params), paste(parts, collapse = "; "))
        }), tighter_txt), collapse = " NOTE: "),
      llm_summary = {
        badge      <- if (isTRUE(pass)) "PASS" else "FAIL"
        paste(c(
          sprintf("D16 | Subsample Stability | %s", badge),
          sprintf("  params=%d subsamples=%d unstable=%d/%d max|z|=%.2f z_tol=%.2f",
                  n_par, n_sub, length(unstable_params), n_par, max_abs_z, z_tol),
          if (length(ci_fail_params) > 0)
            sprintf("  ci_nonoverlap: %s", paste(ci_fail_params, collapse = ", ")),
          if (length(drift_fail_params) > 0)
            sprintf("  z_beyond_%.2f: %s", z_tol, paste(drift_fail_params, collapse = ", ")),
          if (!is.null(tighter_txt)) sprintf("  note: %s", tighter_txt),
          sprintf("  action: %s",
                  if (isTRUE(pass))
                    "Parameter estimates stable across subsamples (Hausman-standardised shift + CI overlap)."
                  else if (length(drift_fail_params) > 0)
                    sprintf("%s shift by more than %.2f standard errors of the subsample-vs-full difference. Check for structural breaks or consider a time-varying parameter model.",
                            paste(head(drift_fail_params, 3), collapse = ", "), z_tol)
                  else
                    sprintf("%s CI non-overlap across subsamples. Check for structural breaks.",
                            paste(head(ci_fail_params, 3), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}
