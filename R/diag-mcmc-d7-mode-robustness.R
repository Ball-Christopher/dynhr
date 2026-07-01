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
#' Pass criterion: >= 2 chains find the best mode (within \code{tol} nats).
#' The \code{gap_warn} threshold is reported in the diagnostic summary but
#' does \strong{not} affect the pass/fail gate; it is a user-readable flag
#' for unusually large mode spreads.
#' The spread of converged modes is also reported for diagnosis.
#'
#' @param mode_multistart  Output of run_mode_parallel() or readRDS("mode_multistart.rds")
#'                         Expected structure: list with $results (list of per-chain results,
#'                         each with $logpost and $theta_mode) and $best
#' @param tol              Tolerance in nats for "same basin" (default 3.0)
#' @param gap_warn         Warn if any chain is > this many nats below best (default 50)
#' @return dynhr_diagnostic list
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#'   using multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'
#' @noRd
d7_mode_robustness <- function(mode_multistart, tol = 3.0, gap_warn = 50,
                                meta = NULL) {

  ## ---- Extract per-chain logpost values ----
  # Handle different structures from run_mode_parallel
  if (!is.null(mode_multistart$results)) {
    chain_results <- mode_multistart$results
  } else if (!is.null(mode_multistart$chains)) {
    chain_results <- mode_multistart$chains
  } else {
    # Try to extract from top-level list elements
    chain_results <- mode_multistart[sapply(mode_multistart, function(x)
      is.list(x) && !is.null(x$logpost))]
  }

  n_chains <- length(chain_results)
  if (n_chains < 2) {
    return(.make_result(
      result  = list(n_chains = n_chains),
      pass    = NA,
      plots   = list(),
      summary = sprintf("D7 Mode robustness: only %d chain -- need >= 2 for comparison.", n_chains)
    ))
  }

  ## Extract logpost and theta from each chain
  logposts <- vapply(chain_results, function(ch) {
    if (!is.null(ch$logpost)) ch$logpost
    else if (!is.null(ch$value)) ch$value
    else NA_real_
  }, numeric(1))

  thetas <- lapply(chain_results, function(ch) {
    if (!is.null(ch$theta_mode)) ch$theta_mode
    else if (!is.null(ch$par)) ch$par
    else NULL
  })

  best_idx   <- which.max(logposts)
  best_lp    <- logposts[best_idx]
  best_theta <- thetas[[best_idx]]

  ## ---- Classify chains by basin ----
  gap_from_best <- best_lp - logposts   # positive = worse than best
  in_best_basin <- gap_from_best <= tol
  n_at_best     <- sum(in_best_basin)

  ## ---- Parameter spread at best basin ----
  if (n_at_best >= 2 && !is.null(best_theta)) {
    basin_thetas <- do.call(rbind, thetas[in_best_basin])
    if (!is.null(basin_thetas) && nrow(basin_thetas) >= 2) {
      param_spread <- apply(basin_thetas, 2, function(col) max(col) - min(col))
      max_spread_param <- names(which.max(param_spread))
      max_spread_val   <- max(param_spread)
    } else {
      param_spread     <- NULL
      max_spread_param <- "unknown"
      max_spread_val   <- NA
    }
  } else {
    param_spread     <- NULL
    max_spread_param <- "N/A"
    max_spread_val   <- NA
  }

  ## ---- Detect distinct local modes ----
  # Cluster chains by logpost proximity
  sorted_lp <- sort(logposts, decreasing = TRUE)
  basins     <- list()
  assigned   <- rep(FALSE, n_chains)

  for (i in order(logposts, decreasing = TRUE)) {
    if (assigned[i]) next
    basin_members <- which(abs(logposts - logposts[i]) <= tol & !assigned)
    basins <- c(basins, list(list(
      center  = logposts[i],
      members = basin_members,
      n       = length(basin_members)
    )))
    assigned[basin_members] <- TRUE
  }

  n_basins <- length(basins)

  ## ---- Build summary table ----
  chain_tbl <- data.frame(
    chain    = seq_len(n_chains),
    logpost  = round(logposts, 2),
    gap      = round(gap_from_best, 2),
    at_best  = ifelse(in_best_basin, "YES", ""),
    stringsAsFactors = FALSE
  )

  ## ---- Pass/fail logic ----
  # Pass if: best basin found by >= 2 chains
  # Warn if: large spread or many distinct basins
  pass <- n_at_best >= 2

  ## ---- Build summary string ----
  lines <- c(
    sprintf("D7 Mode robustness: %d chains, %d distinct basins (tol=%.1f nats).",
            n_chains, n_basins, tol),
    sprintf("  Best logpost: %.2f (chain %d)", best_lp, best_idx),
    sprintf("  Chains at best basin: %d/%d %s",
            n_at_best, n_chains,
            if (pass) "PASS" else "FAIL -- mode found by only 1 chain"),
    sprintf("  Logpost range: [%.2f, %.2f] (spread=%.2f nats)",
            min(logposts), max(logposts), max(logposts) - min(logposts))
  )

  if (!is.na(max_spread_val)) {
    lines <- c(lines, sprintf(
      "  Max param spread at best basin: %s (%.4f)",
      max_spread_param, max_spread_val
    ))
  }

  # Flag chains far from best
  far_chains <- which(gap_from_best > gap_warn)
  if (length(far_chains) > 0) {
    lines <- c(lines, sprintf(
      "  WARNING: %d chain(s) >%.0f nats from best: %s",
      length(far_chains), gap_warn,
      paste(sprintf("ch%d (gap=%.1f)", far_chains, gap_from_best[far_chains]),
            collapse = ", ")
    ))
  }

  ## ---- Plot: logpost by chain ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plt_df <- data.frame(
      chain   = factor(seq_len(n_chains)),
      logpost = logposts,
      basin   = ifelse(in_best_basin, "Best basin", "Other")
    )
    p_mc <- ggplot2::ggplot(plt_df,
                                             ggplot2::aes(x = chain, y = logpost, fill = basin)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_hline(yintercept = best_lp - tol, linetype = "dashed",
                          colour = dynhr_colours$red, linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = c("Best basin" = dynhr_colours$mid_blue,
                                            "Other" = dynhr_colours$grey)) +
      ggplot2::labs(title = "D7: Mode-Finding Robustness",
                    subtitle = sprintf("Best=%.2f | %d/%d chains in best basin (tol=%.1f nats)",
                                       best_lp, n_at_best, n_chains, tol),
                    x = "Chain", y = "Log-posterior at mode", fill = NULL) +
      theme_dynhr_diagnostic()
    plots$mode_comparison <- .apply_meta(p_mc, meta)
  }

  .make_result(
    result = list(
      n_chains        = n_chains,
      logposts        = logposts,
      best_idx        = best_idx,
      best_logpost    = best_lp,
      n_at_best       = n_at_best,
      n_basins        = n_basins,
      basins          = basins,
      gap_from_best   = gap_from_best,
      chain_table     = chain_tbl,
      param_spread    = param_spread
    ),
    pass    = pass,
    plots   = plots,
    summary = paste(lines, collapse = "\n"),
    llm_summary = {
      badge <- if (pass) "PASS" else "FAIL"
      paste(c(
        sprintf("D7 | Mode Robustness | %s", badge),
        sprintf("  n_starts=%d n_basins=%d best_logpost=%.2f",
                n_chains, n_basins, best_lp),
        sprintf("  n_at_best=%d/%d (%.0f%%)",
                n_at_best, n_chains, 100 * n_at_best / max(1, n_chains)),
        if (!is.null(max_spread_param) && max_spread_param != "N/A")
          sprintf("  max_param_spread=%s=%.4f", max_spread_param, max_spread_val),
        sprintf("  action: %s",
                if (pass)
                  "All starts converged to same mode. Mode finding is robust."
                else if (n_basins > 1)
                  sprintf("Multiple modes detected (%d basins). Posterior may be multimodal. Try SMC sampler or more starts.",
                          n_basins)
                else
                  sprintf("Only %d/%d starts reached best mode. Try more starts (n_starts=16+) or check starting values.",
                          n_at_best, n_chains))
      ), collapse = "\n")
    }
  )
}
