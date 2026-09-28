## R/diag-mcmc-d7-mode-robustness.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D7 mode-finding robustness across starts
## --------------------------------------------------------------------------

#' D7. Mode-finding robustness
#'
#' Evaluates whether the posterior mode is robust across multiple dispersed
#' starting values. Uses the output from run_mode_parallel() which already
#' runs N chains from perturbed initial conditions.
#'
#' A start whose log-posterior is non-finite, or at/below the infeasible
#' floor (\code{-1e19}; the log-posterior's \code{-Inf} substitute is
#' \code{-1e20}), is a FAILED start: it is counted and reported but takes no
#' part in the best mode, the basins, the spread or the plot's gap axis.
#'
#' Pass criterion: >= 2 converged starts reach the best mode (within
#' \code{tol} nats of it). Fewer than 2 starts supplied gives \code{NA}
#' (nothing to compare); >= 2 starts of which fewer than 2 converged is a
#' FAIL. The \code{gap_warn} threshold is reported in the summary and drawn
#' on the plot but does \strong{not} affect the pass/fail gate.
#' Basins are formed greedily on the log-posterior value alone (a start joins
#' the highest unassigned mode within \code{tol} below it); parameter
#' agreement inside the best basin is reported as \code{param_spread}
#' (max - min per parameter, in the parameters' own units) but is not gated.
#'
#' @section Why \code{param_spread} is reported, not gated:
#' D7 answers one question -- is there more than one mode? -- and the
#' log-posterior value is the criterion that answers it: two starts within a
#' few nats of the same optimum ARE at the same mode, by definition. Starts
#' that agree on the log-posterior but disagree on the parameters are sitting
#' on a FLAT RIDGE, which is a weak-identification symptom, not a multimodality
#' one. Gating D7 on the spread would double-count that problem and FAIL
#' well-behaved-but-weakly-identified models on the wrong diagnostic; the
#' rank-based identification diagnostics (D1, D20, D27) are where a ridge
#' belongs, so cross-check a large \code{param_spread} against those. This
#' matches how Dynare separates the two concerns operationally
#' (\code{mode_check} slices for curvature; re-running \code{mode_compute} from
#' other starts for multimodality) -- no inspected implementation gates on a
#' scaled parameter distance within a basin.
#'
#' @param mode_multistart  Output of run_mode_parallel() (uses \code{$chains}),
#'   a list with \code{$results} (per-start results), or a bare list of
#'   per-start results; each result carries \code{$logpost} (or
#'   \code{$value}) and \code{$theta_mode} (or \code{$par}).
#' @param tol       Tolerance in nats for "same basin" (default 3.0)
#' @param gap_warn  Flag any converged start more than this many nats below
#'   the best (default 50)
#' @param meta      Optional diag_meta() caption metadata.
#' @return dynhr_diagnostic list
#'
#' @noRd
d7_mode_robustness <- function(mode_multistart, tol = 3.0, gap_warn = 50,
                                meta = NULL) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    .dynhr_abort("`tol` must be a single finite number >= 0.")
  if (!is.numeric(gap_warn) || length(gap_warn) != 1L ||
      !is.finite(gap_warn) || gap_warn < 0)
    .dynhr_abort("`gap_warn` must be a single finite number >= 0.")
  if (!is.list(mode_multistart))
    .dynhr_abort("`mode_multistart` must be a list (run_mode_parallel() output).")

  ## ---- Extract per-start results ----
  if (!is.null(mode_multistart[["results"]])) {
    chain_results <- mode_multistart[["results"]]
  } else if (!is.null(mode_multistart[["chains"]])) {
    chain_results <- mode_multistart[["chains"]]
  } else {
    chain_results <- mode_multistart[vapply(mode_multistart, function(x)
      is.list(x) && (!is.null(x$logpost) || !is.null(x$value)), logical(1))]
  }

  n_chains <- length(chain_results)
  if (n_chains < 2) {
    return(.make_result(
      result  = list(n_chains = n_chains),
      pass    = NA,
      plots   = list(),
      summary = sprintf("D7 Mode robustness: only %d start(s) supplied -- need >= 2 for comparison.",
                        n_chains)
    ))
  }

  logposts <- vapply(chain_results, function(ch) {
    v <- if (is.list(ch) && !is.null(ch$logpost)) ch$logpost
         else if (is.list(ch) && !is.null(ch$value)) ch$value
         else NA_real_
    if (is.numeric(v) && length(v) == 1L) as.numeric(v) else NA_real_
  }, numeric(1))
  thetas <- lapply(chain_results, function(ch) {
    if (!is.list(ch)) NULL
    else if (!is.null(ch$theta_mode)) ch$theta_mode
    else ch$par
  })

  ## ---- Failed starts: non-finite or at the infeasible floor ----
  ok       <- is.finite(logposts) & logposts > -1e19
  failed   <- which(!ok)
  n_failed <- length(failed)
  n_valid  <- sum(ok)

  if (n_valid > 0L) {
    best_idx <- which(ok)[which.max(logposts[ok])]
    best_lp  <- logposts[best_idx]
  } else {
    best_idx <- NA_integer_
    best_lp  <- NA_real_
  }
  gap_from_best <- ifelse(ok, best_lp - logposts, NA_real_)  # >= 0; NA = failed
  in_best_basin <- ok & !is.na(gap_from_best) & gap_from_best <= tol
  n_at_best     <- sum(in_best_basin)

  ## ---- Parameter spread within the best basin ----
  param_spread     <- NULL
  max_spread_param <- NA_character_
  max_spread_val   <- NA_real_
  basin_th <- Filter(Negate(is.null), thetas[in_best_basin])
  if (length(basin_th) >= 2L &&
      length(unique(vapply(basin_th, length, integer(1)))) == 1L) {
    basin_thetas <- do.call(rbind, lapply(basin_th, as.numeric))
    colnames(basin_thetas) <- names(basin_th[[1L]])
    param_spread <- apply(basin_thetas, 2, function(col) max(col) - min(col))
    if (length(param_spread) > 0L && all(is.finite(param_spread))) {
      k <- which.max(param_spread)
      max_spread_param <- if (is.null(names(param_spread))) sprintf("theta[%d]", k)
                          else names(param_spread)[k]
      max_spread_val   <- param_spread[[k]]
    }
  }

  ## ---- Basins: greedy on logpost, converged starts only ----
  basins   <- list()
  assigned <- !ok
  for (i in order(logposts, decreasing = TRUE)) {
    if (assigned[i]) next
    members <- which(!assigned & logposts <= logposts[i] &
                     logposts[i] - logposts <= tol)
    basins <- c(basins, list(list(center = logposts[i], members = members,
                                  n = length(members))))
    assigned[members] <- TRUE
  }
  n_basins <- length(basins)

  chain_tbl <- data.frame(
    chain   = seq_len(n_chains),
    logpost = round(logposts, 2),
    gap     = round(gap_from_best, 2),
    status  = ifelse(!ok, "failed", ifelse(in_best_basin, "best", "other")),
    stringsAsFactors = FALSE
  )

  pass       <- n_at_best >= 2L
  spread     <- if (n_valid > 0L) best_lp - min(logposts[ok]) else NA_real_
  far_chains <- which(ok & gap_from_best > gap_warn)

  ## ---- Summary ----
  lines <- c(
    sprintf("D7 Mode robustness: %d starts (%d converged, %d failed), %d distinct basin(s) (tol=%.1f nats).",
            n_chains, n_valid, n_failed, n_basins, tol),
    if (n_valid > 0L) sprintf("  Best logpost: %.2f (start %d)", best_lp, best_idx)
    else "  No start converged to a finite log-posterior.",
    sprintf("  Starts at best mode: %d/%d -- %s",
            n_at_best, n_chains,
            if (pass) "PASS"
            else if (n_valid < 2L) "FAIL (fewer than 2 starts converged)"
            else "FAIL (best mode found by only 1 start)"),
    if (n_valid > 0L)
      sprintf("  Converged logpost range: [%.2f, %.2f] (spread=%.2f nats)",
              min(logposts[ok]), best_lp, spread),
    if (n_failed > 0L)
      sprintf("  Failed starts: %s", paste(sprintf("s%d", failed), collapse = ", ")),
    if (!is.na(max_spread_val))
      sprintf("  Max param spread at best mode: %s (%.4g, parameter units)",
              max_spread_param, max_spread_val),
    if (length(far_chains) > 0L)
      sprintf("  WARNING: %d start(s) >%.0f nats below best: %s",
              length(far_chains), gap_warn,
              paste(sprintf("s%d (gap=%.1f)", far_chains, gap_from_best[far_chains]),
                    collapse = ", "))
  )

  ## ---- Plot: gap below the best mode, per start ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$mode_comparison <- .apply_meta(
      .d7_plot_gaps(gap_from_best, ok, in_best_basin, tol, gap_warn,
                    best_lp, best_idx, n_at_best),
      meta)
  }

  badge <- if (pass) "PASS" else "FAIL"
  action <- if (pass && n_at_best == n_chains)
    "All starts converged to the same mode. Mode finding is robust."
  else if (pass)
    sprintf("%d/%d starts reached the best mode; %d found other local modes, %d failed. The best mode is reproducible, but inspect the other basins.",
            n_at_best, n_chains, n_valid - n_at_best, n_failed)
  else if (n_valid < 2L)
    sprintf("Only %d/%d starts converged. Check starting values, prior support and the optimiser budget.",
            n_valid, n_chains)
  else
    sprintf("Best mode found by only 1 start (%d basins). Posterior may be multimodal: try more starts (n_starts=16+) or an SMC sampler.",
            n_basins)
  llm <- paste(c(
    sprintf("D7 | Mode Robustness | %s", badge),
    sprintf("  n_starts=%d n_converged=%d n_failed=%d n_basins=%d best_logpost=%s",
            n_chains, n_valid, n_failed, n_basins,
            if (n_valid > 0L) sprintf("%.2f", best_lp) else "NA"),
    sprintf("  n_at_best=%d/%d (%.0f%%)", n_at_best, n_chains,
            100 * n_at_best / n_chains),
    if (!is.na(max_spread_val))
      sprintf("  max_param_spread=%s=%.4g", max_spread_param, max_spread_val),
    sprintf("  action: %s", action)
  ), collapse = "\n")

  .make_result(
    result = list(
      n_chains      = n_chains,
      n_converged   = n_valid,
      n_failed      = n_failed,
      failed        = failed,
      logposts      = logposts,
      best_idx      = best_idx,
      best_logpost  = best_lp,
      n_at_best     = n_at_best,
      n_basins      = n_basins,
      basins        = basins,
      gap_from_best = gap_from_best,
      spread        = spread,
      chain_table   = chain_tbl,
      param_spread  = param_spread
    ),
    pass        = pass,
    plots       = plots,
    summary     = paste(lines, collapse = "\n"),
    llm_summary = llm
  )
}

## Lollipop plot of each start's gap below the best mode on a pseudo-log
## axis (0 = best), with the tol (basin) and gap_warn reference lines.
## Failed starts are drawn as crosses above the largest gap.
#' @noRd
.d7_plot_gaps <- function(gap_from_best, ok, in_best_basin, tol, gap_warn,
                          best_lp, best_idx, n_at_best) {
  n_chains <- length(gap_from_best)
  n_failed <- sum(!ok)
  lv <- c("At best mode (within tol)", "Other local mode", "Failed start")
  top <- max(c(gap_from_best[ok], gap_warn, tol, 1)) * 2
  df <- data.frame(
    chain  = factor(seq_len(n_chains)),
    gap    = ifelse(ok, gap_from_best, top),
    status = factor(ifelse(!ok, lv[3], ifelse(in_best_basin, lv[1], lv[2])),
                    levels = lv)
  )
  ref_cols <- unname(tol_vibrant[c("teal", "red")])
  ref_df <- data.frame(y   = c(tol, gap_warn),
                       lab = c(sprintf("tol = %g nats (same basin)", tol),
                               sprintf("gap_warn = %g nats", gap_warn)))
  cols <- stats::setNames(unname(tol_vibrant[c("blue", "orange", "grey")]), lv)
  shapes <- stats::setNames(c(16, 16, 4), lv)
  subtitle <- if (any(ok))
    sprintf("Best log-posterior %.2f (start %d) | %d/%d starts within tol%s",
            best_lp, best_idx, n_at_best, n_chains,
            if (n_failed > 0L) sprintf(" | %d failed (crosses at top)", n_failed)
            else "")
  else sprintf("All %d starts failed (crosses at top)", n_chains)

  ggplot2::ggplot(df, ggplot2::aes(x = chain, y = gap)) +
    ggplot2::scale_x_discrete(limits = levels(df$chain)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::geom_hline(data = ref_df, ggplot2::aes(yintercept = y),
                        linetype = c("dashed", "dotted"), colour = ref_cols,
                        linewidth = 0.6) +
    ggplot2::annotate("text", x = 0.5, y = ref_df$y, label = ref_df$lab,
                      hjust = 0, vjust = -0.4, size = 3.2, colour = ref_cols) +
    ggplot2::geom_segment(data = df[ok, , drop = FALSE],
                          ggplot2::aes(xend = chain, y = 0, yend = gap,
                                       colour = status),
                          linewidth = 0.8, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(colour = status, shape = status), size = 3.5) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::scale_y_continuous(
      transform = scales::pseudo_log_trans(sigma = 1, base = 10),
      breaks = c(0, 1, 3, 10, 30, 100, 300, 1e3, 1e4, 1e5, 1e6, 1e7),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE,
                                  trim = TRUE, drop0trailing = TRUE),
      expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(
      title    = "D7: Mode-finding robustness across starts",
      subtitle = subtitle,
      x = "Start (1 = supplied initial value)",
      y = "Nats below best mode (pseudo-log scale)",
      colour = NULL, shape = NULL) +
    theme_dynhr_diagnostic()
}
