## R/diag-post-d18-welfare.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D18 welfare plausibility stub
## --------------------------------------------------------------------------

#' D18. Welfare plausibility
#'
#' Uses the Ramsey workflow output, when available, to report welfare metrics.
#' If Ramsey results are absent, returns an informational diagnostic.
#'
#' The welfare gap is reported as a percentage of the absolute steady-state
#' welfare value, making it model-unit-independent.  By default the diagnostic
#' is informational only (\code{pass = NA}).  An optional \code{gap_threshold}
#' (percent) can be supplied via \code{...}; if provided, the diagnostic FAILs
#' when \code{|gap / steady_welfare| * 100 > gap_threshold}.
#'
#' @param ramsey_result Optional \code{dynhr_ramsey_result} or
#'   \code{dynhr_ramsey_result2}.  Must contain \code{$welfare$gap_vs_steady}.
#' @param ... Optional named arguments passed through.  Recognised:
#'   \describe{
#'     \item{\code{gap_threshold}}{Numeric: maximum acceptable welfare gap as a
#'       percentage of steady-state welfare.  When provided, \code{pass = FALSE}
#'       if \code{|gap_pct| > gap_threshold}.  No default (informational only).}
#'     \item{\code{meta}}{Metadata list for plot annotation.}
#'   }
#' @return A \code{dynhr_diagnostic} list.  \code{pass = NA} by default
#'   (informational); \code{pass = TRUE/FALSE} only when \code{gap_threshold}
#'   is supplied via \code{...}.
#' @noRd
d18_welfare_plausibility <- function(ramsey_result = NULL, ...) {
  # Accept any Ramsey object (dynhr_ramsey_result OR dynhr_ramsey_result2) that
  # carries a welfare block with a steady-state gap.
  .has_welfare <- !is.null(ramsey_result) &&
    !is.null(ramsey_result$welfare) &&
    !is.null(ramsey_result$welfare$gap_vs_steady)
  if (!.has_welfare) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = paste(
        "D18 Welfare plausibility: Ramsey output not provided.",
        "Run run_full_estimation(..., run_ramsey = TRUE) to enable welfare checks."
      ),
      llm_summary = paste(
        "[INFO] D18 Welfare plausibility",
        "status=skipped reason=no_ramsey_result",
        "action: rerun with run_ramsey=TRUE"
      )
    ))
  }

  w <- ramsey_result$welfare
  # Guard against zero-length slots that arise when Ramsey solves only
  # partially (e.g. objective_mean is computed but steady_value is not yet
  # filled in).  These cause "argument is of length zero" errors downstream.
  .safe_num <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    as.numeric(x[[1L]])
  }
  gap              <- .safe_num(w$gap_vs_steady)
  unconditional_v  <- .safe_num(w$unconditional_value)
  steady_v         <- .safe_num(w$steady_value)
  discount_v       <- .safe_num(w$discount)
  objective_mean_v <- .safe_num(w$objective_mean)

  # D18 is demoted to ALWAYS-INFO (pass = NA).
  # The gap is reported relative to |steady-state welfare| so that it is
  # model-unit-independent.  A model-specific threshold can be supplied via
  # gap_threshold; if provided, FAIL when |gap_pct| > gap_threshold (percent).
  gap_threshold <- list(...)$gap_threshold
  pass <- if (!is.null(gap_threshold) && is.finite(gap) && is.finite(steady_v) &&
               abs(steady_v) > 1e-16) {
    abs(gap / steady_v) * 100 <= gap_threshold
  } else {
    NA  # default: informational only
  }

  # Report gap as % of |steady-state welfare|
  gap_pct <- if (is.finite(gap) && is.finite(steady_v) && abs(steady_v) > 1e-16)
    gap / abs(steady_v) * 100 else NA_real_

  # Plot: unconditional (stochastic) vs deterministic-steady-state welfare.
  meta <- list(...)$meta
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      is.finite(unconditional_v) && is.finite(steady_v)) {
    wdf <- data.frame(
      Measure = factor(c("Unconditional\n(stochastic)", "Deterministic\nsteady state"),
                       levels = c("Unconditional\n(stochastic)", "Deterministic\nsteady state")),
      Welfare = c(unconditional_v, steady_v))
    gap_label <- if (is.finite(gap_pct)) sprintf("Gap = %.2f%% of |SS welfare|", gap_pct)
                 else sprintf("Gap = %.4g", gap)
    p_w <- ggplot2::ggplot(wdf, ggplot2::aes(x = Measure, y = Welfare, fill = Measure)) +
      ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4g", Welfare)),
                         vjust = -0.4, size = 3, colour = "grey15") +
      ggplot2::scale_fill_manual(values = c(dynhr_colours$mid_blue, dynhr_colours$teal)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.15))) +
      theme_dynhr() +
      ggplot2::labs(
        title = "D18: Welfare -- unconditional vs steady state",
        subtitle = sprintf("%s (discount = %.3f)", gap_label, discount_v %||% NA),
        x = NULL, y = "Planner welfare")
    plots$welfare <- .apply_meta(p_w, meta)
  }

  badge <- if (is.na(pass)) "INFO" else if (isTRUE(pass)) "PASS" else "FAIL"

  .make_result(
    result = list(
      welfare_unconditional = unconditional_v,
      welfare_steady        = steady_v,
      welfare_gap           = gap,
      welfare_gap_pct       = gap_pct,
      objective_mean        = objective_mean_v,
      discount              = discount_v
    ),
    pass = pass,
    plots = plots,
    summary = sprintf(
      "D18 Welfare plausibility [%s]: gap vs steady = %.6f (%s) (discount=%.3f).",
      badge,
      if (is.finite(gap)) gap else NA_real_,
      if (is.finite(gap_pct)) sprintf("%.2f%% of |SS welfare|", gap_pct) else "SS welfare unavailable",
      if (is.finite(discount_v)) discount_v else NA_real_
    ),
    llm_summary = {
      gap_note <- if (is.finite(gap_pct))
        sprintf("gap_pct_of_ss=%.2f%% (guidance: gaps >10%% may indicate scaling issues)", gap_pct)
      else
        "gap_pct_of_ss=NA (steady_value unavailable)"
      paste(c(
        sprintf("D18 | Welfare Plausibility | %s", badge),
        sprintf("  welfare_unconditional=%s welfare_steady=%s gap=%s discount=%s",
                if (is.finite(unconditional_v)) sprintf("%.6f", unconditional_v) else "NA",
                if (is.finite(steady_v))         sprintf("%.6f", steady_v)        else "NA",
                if (is.finite(gap))              sprintf("%.6f", gap)             else "NA",
                if (is.finite(discount_v))       sprintf("%.3f", discount_v)      else "NA"),
        sprintf("  %s", gap_note),
        sprintf("  action: %s",
                if (!is.null(gap_threshold))
                  sprintf("gap_threshold=%.1f%%. %s.",
                          gap_threshold,
                          if (isTRUE(pass)) "Gap within threshold." else "Gap exceeds threshold -- inspect planner objective scaling.")
                else
                  "Informational: no threshold set. Use gap_threshold= argument to activate gate.")
      ), collapse = "\n")
    }
  )
}
