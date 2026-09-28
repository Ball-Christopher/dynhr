## R/diag-post-d18-welfare.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D18 welfare plausibility
## --------------------------------------------------------------------------

#' D18. Welfare plausibility
#'
#' Reads the welfare block of a Ramsey result and reports the stochastic
#' welfare gap
#' \deqn{\Delta W = W^{unc} - W^{ss},\qquad W = E\sum_t \beta^t u_t,}
#' where \eqn{W^{unc} = E[u]/(1-\beta)} is unconditional (stochastic) welfare
#' and \eqn{W^{ss} = u(\bar y)/(1-\beta)} is welfare at the deterministic
#' steady state.  \strong{Sign convention:} \eqn{\Delta W > 0} is a welfare
#' GAIN from uncertainty, \eqn{\Delta W < 0} a LOSS.
#'
#' Reported measures:
#' \describe{
#'   \item{\code{welfare_gap}}{\eqn{\Delta W}, in discounted-sum objective units.
#'     Recomputed as \eqn{W^{unc} - W^{ss}} whenever both components are
#'     available; a stored \code{gap_vs_steady} that disagrees is flagged.}
#'   \item{\code{welfare_gap_per_period}}{\eqn{(1-\beta)\Delta W = E[u] - u(\bar y)},
#'     in per-period objective units.}
#'   \item{\code{cev_pct}}{Consumption-equivalent variation in PERCENT, only when
#'     the planner objective is literally \code{log(x)} for one variable
#'     \code{x}: then \eqn{\lambda = \exp((1-\beta)\Delta W) - 1} is the
#'     permanent proportional change in \code{x} that is welfare-equivalent to
#'     the gap.  \code{NA} for any other objective (a CEV needs the utility form).}
#'   \item{\code{welfare_gap_pct}}{\eqn{100\,\Delta W / |W^{ss}|}.  Invariant to
#'     rescaling the objective but NOT to adding a constant to it (and
#'     undefined when \eqn{W^{ss}=0}, e.g. \code{log(c)} with \eqn{\bar c = 1});
#'     treat it as a scaling heuristic, not a welfare measure.}
#' }
#'
#' A first-order (\code{order = 1}) solution omits the second-order shift in
#' the means of the endogenous variables, so its welfare gap -- and welfare
#' rankings built on it -- can have the wrong sign; the summary says so.
#'
#' Accepts \code{dynhr_ramsey_result} (\code{ramsey_policy()}) and
#' \code{dynhr_ramsey_result2} (\code{ramsey_model()}, which stores
#' \code{welfare$unconditional} when it simulated).  The \code{nn1} path at
#' \code{order = 1} stores the steady value as a placeholder for the
#' unconditional one; that is treated as "not computed".
#'
#' @param ramsey_result Optional \code{dynhr_ramsey_result} or
#'   \code{dynhr_ramsey_result2}.
#' @param gap_threshold Optional non-negative number: maximum acceptable
#'   \code{|welfare_gap_pct|}.  When supplied the diagnostic PASSes/FAILs on it
#'   (\code{pass = NA} if the percentage is undefined).  Default \code{NULL}:
#'   informational only.
#' @param meta Optional \code{dynhr_diag_meta} for plot captions.
#' @return A \code{dynhr_diagnostic}.  \code{pass = FALSE} whenever the welfare
#'   block is internally inconsistent (stored gap != components, or discount
#'   outside (0, 1)), otherwise as described for \code{gap_threshold}.
#' @noRd
d18_welfare_plausibility <- function(ramsey_result = NULL, gap_threshold = NULL,
                                     meta = NULL) {
  if (!is.null(gap_threshold) &&
      !(is.numeric(gap_threshold) && length(gap_threshold) == 1L &&
        is.finite(gap_threshold) && gap_threshold >= 0)) {
    .dynhr_abort("D18: `gap_threshold` must be a single finite non-negative number (percent).")
  }

  .skip <- function(reason, action) {
    .make_result(
      result  = NULL, pass = NA, plots = list(),
      summary = sprintf("D18 Welfare plausibility: skipped (%s).", reason),
      llm_summary = paste(
        "[INFO] D18 Welfare plausibility",
        sprintf("status=skipped reason=%s", gsub(" ", "_", reason)),
        sprintf("action: %s", action), sep = "\n"))
  }
  if (is.null(ramsey_result) || is.null(ramsey_result$welfare)) {
    return(.skip("no Ramsey result",
                 "rerun with run_ramsey=TRUE or pass ramsey_result="))
  }

  w <- ramsey_result$welfare
  # Zero-length / NULL slots (partial Ramsey solves) become NA.
  .num <- function(x) {
    if (is.null(x) || length(x) == 0L || !is.numeric(x)) return(NA_real_)
    as.numeric(x[[1L]])
  }
  rmeta      <- ramsey_result$meta
  order      <- .num(rmeta$order)
  steady_v   <- .num(w$steady_value)
  uncond_v   <- .num(if (!is.null(w$unconditional_value)) w$unconditional_value
                     else w$unconditional)
  gap_stored <- .num(w$gap_vs_steady)
  discount_v <- .num(w$discount)
  obj_mean_v <- .num(w$objective_mean)

  # nn1 order-1: there IS no stochastic welfare (certainty equivalence), and
  # since 0.9.4 the solver reports NA plus a machine-readable reason instead of
  # copying the steady-state value in as a placeholder.  Keep the belt-and-
  # braces override for results produced by older versions.
  reason      <- w$unconditional_reason
  placeholder <- (identical(rmeta$method, "nn1") && isTRUE(order == 1)) ||
    (is.character(reason) && length(reason) == 1L && !is.na(reason) &&
     identical(reason, "order1_certainty_equivalent"))
  if (placeholder) uncond_v <- NA_real_

  have_parts <- is.finite(uncond_v) && is.finite(steady_v)
  gap <- if (have_parts) uncond_v - steady_v else gap_stored
  if (!is.finite(gap)) {
    return(.skip(
      if (placeholder)
        paste0("order-1 Ramsey solution is certainty-equivalent: the ",
               "unconditional mean IS the steady state, so there is no ",
               "stochastic welfare gap to report")
      else "no stochastic welfare in the Ramsey result",
      "solve the Ramsey problem with a stochastic welfare evaluation (e.g. ramsey_policy(order = 2))"))
  }

  inconsistent <- have_parts && is.finite(gap_stored) &&
    abs(gap_stored - gap) > 1e-8 * max(1, abs(uncond_v), abs(steady_v))
  discount_ok <- is.finite(discount_v) && discount_v > 0 && discount_v < 1
  discount_bad <- is.finite(discount_v) && !discount_ok

  gap_pp <- if (discount_ok) (1 - discount_v) * gap else NA_real_
  gap_pct <- if (is.finite(steady_v) && abs(steady_v) > 1e-16)
    100 * gap / abs(steady_v) else NA_real_

  obj_text <- ramsey_result$objective$text
  log_obj <- FALSE
  if (is.character(obj_text) && length(obj_text) == 1L && !is.na(obj_text)) {
    # The parser stores the objective wrapped in parentheses: "(log(c))".
    txt <- gsub("\\s+", "", obj_text)
    while (grepl("^\\(.*\\)$", txt)) txt <- substr(txt, 2L, nchar(txt) - 1L)
    log_obj <- grepl("^log\\([A-Za-z_.][A-Za-z0-9_.]*\\)$", txt)
  }
  cev_pct <- if (log_obj && is.finite(gap_pp)) 100 * expm1(gap_pp) else NA_real_

  pass <- if (inconsistent || discount_bad) {
    FALSE
  } else if (!is.null(gap_threshold) && is.finite(gap_pct)) {
    abs(gap_pct) <= gap_threshold
  } else {
    NA
  }
  badge <- if (is.na(pass)) "INFO" else if (isTRUE(pass)) "PASS" else "FAIL"

  direction <- if (gap > 0) "gain" else if (gap < 0) "loss" else "none"
  .f <- function(x, fmt = "%.6g") if (is.finite(x)) sprintf(fmt, x) else "NA"
  order_note <- if (isTRUE(order == 1))
    "order=1: first-order welfare omits the second-order mean shift; the gap's sign and welfare rankings may be wrong -- use order=2"
  else NULL
  problems <- c(
    if (inconsistent) sprintf("stored gap_vs_steady=%s != unconditional - steady=%s",
                              .f(gap_stored), .f(gap)),
    if (discount_bad) sprintf("discount=%s is outside (0, 1)", .f(discount_v)),
    if (!is.null(gap_threshold) && !is.finite(gap_pct) && !inconsistent && !discount_bad)
      "gap_threshold not applied: |steady-state welfare| is zero or unavailable"
  )

  # ---- Plot: the gap itself (levels are affine-arbitrary) ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    per_period <- is.finite(gap_pp)
    val <- if (per_period) gap_pp else gap
    lab <- if (is.finite(cev_pct)) sprintf("%.4g  (CEV %+.3g%%)", val, cev_pct)
           else sprintf("%.4g", val)
    gdf <- data.frame(
      what = "Stochastic - steady state",
      value = val,
      sign = factor(if (val >= 0) "Gain" else "Loss", levels = c("Gain", "Loss")))
    subtitle <- paste0(
      sprintf("%s: discounted gap = %s", toupper(direction), .f(gap, "%.4g")),
      if (is.finite(gap_pct)) sprintf(" (%.3g%% of abs SS welfare)", gap_pct) else "",
      sprintf("; discount = %s; order = %s", .f(discount_v, "%.4g"), .f(order, "%d")),
      if (!is.null(order_note)) "\nOrder-1 welfare: sign may be wrong (no second-order mean shift)" else "")
    p_w <- ggplot2::ggplot(gdf, ggplot2::aes(x = value, y = what, fill = sign)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5) +
      ggplot2::geom_col(width = 0.5, show.legend = FALSE) +
      ggplot2::geom_label(ggplot2::aes(x = value / 2, label = lab),
                          fill = "white", colour = "grey15", size = 3.5,
                          linewidth = 0) +
      ggplot2::scale_fill_manual(values = c(Gain = unname(tol_vibrant["teal"]),
                                            Loss = unname(tol_vibrant["orange"])),
                                 drop = FALSE) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
      ggplot2::expand_limits(x = 0) +
      theme_dynhr() +
      ggplot2::labs(
        title = "D18: Welfare gap from uncertainty (right of 0 = gain)",
        subtitle = subtitle,
        x = if (per_period)
          "(1 - discount) x [W(unconditional) - W(steady state)]  (per-period objective units)"
        else "W(unconditional) - W(steady state)  (discounted objective units)",
        y = NULL)
    plots$welfare <- .apply_meta(p_w, meta)
  }

  summary <- paste0(
    sprintf("D18 Welfare plausibility [%s]: welfare %s from uncertainty, gap = %s (per period %s",
            badge, direction, .f(gap), .f(gap_pp)),
    if (is.finite(cev_pct)) sprintf(", CEV %+.4f%%", cev_pct) else "",
    if (is.finite(gap_pct)) sprintf(", %.2f%% of |SS welfare|", gap_pct) else "",
    sprintf(") (discount=%s, order=%s).", .f(discount_v, "%.4g"), .f(order, "%d")),
    if (length(problems)) paste0(" ", paste(problems, collapse = "; "), ".") else "",
    if (!is.null(order_note)) paste0(" Caveat: ", order_note, ".") else "")

  action <- if (length(problems) && !is.na(pass) && !pass)
    "inspect the Ramsey welfare block (discount / gap bookkeeping)"
  else if (!is.null(order_note))
    "re-evaluate welfare at order=2 before ranking policies"
  else if (!is.null(gap_threshold))
    sprintf("gap_threshold=%.4g%%: %s", gap_threshold,
            if (isTRUE(pass)) "gap within threshold" else if (isFALSE(pass))
              "gap exceeds threshold -- inspect planner objective scaling"
            else "threshold not applicable")
  else "none (informational; pass gap_threshold= to gate on |gap_pct|)"

  llm_summary <- paste(c(
    sprintf("D18 | Welfare Plausibility | %s", badge),
    sprintf("  welfare_unconditional=%s welfare_steady=%s gap=%s direction=%s",
            .f(uncond_v), .f(steady_v), .f(gap), direction),
    sprintf("  gap_per_period=%s cev_pct=%s gap_pct_of_ss=%s discount=%s order=%s",
            .f(gap_pp), .f(cev_pct), .f(gap_pct), .f(discount_v, "%.4g"), .f(order, "%d")),
    if (length(problems)) sprintf("  problem: %s", problems),
    if (!is.null(order_note)) sprintf("  caveat: %s", order_note),
    sprintf("  action: %s", action)
  ), collapse = "\n")

  .make_result(
    result = list(
      welfare_unconditional  = uncond_v,
      welfare_steady         = steady_v,
      welfare_gap            = gap,
      welfare_gap_per_period = gap_pp,
      welfare_gap_pct        = gap_pct,
      cev_pct                = cev_pct,
      direction              = direction,
      objective_mean         = obj_mean_v,
      discount               = discount_v,
      order                  = order,
      gap_consistent         = !inconsistent,
      discount_valid         = discount_ok
    ),
    pass = pass,
    plots = plots,
    summary = summary,
    llm_summary = llm_summary
  )
}
