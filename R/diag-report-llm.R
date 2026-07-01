## R/diag-report-llm.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## write_llm_report() LLM-readable diagnostic report writer
## --------------------------------------------------------------------------

#' Write a consolidated LLM-friendly diagnostic report
#'
#' Produces a single markdown file containing all diagnostic results in a
#' structured, machine-parseable format optimised for copy-pasting into an
#' LLM prompt to facilitate iterative model development.
#'
#' Output structure:
#'   1. Model metadata (dimensions, observables, parent model)
#'   2. Pass/fail dashboard (one-line summary per diagnostic)
#'   3. Delta table vs parent model (logpost, D9 ratios, etc.)
#'   4. Posterior mode estimates with auto-flagging (boundary, >2sd from prior)
#'   5. D5 detail: ESS + R-hat per parameter
#'   6. D9 detail: model vs data SD ratios
#'   7. D10 detail: variance decomposition table
#'   8. D8 detail: IRF benchmark pass/fail
#'   9. Auto-detected issues (boundaries, low ESS, high R-hat, D9 failures)
#'  10. Suggested LLM prompt template
#'
#' @param model_name     Character -- model variant (e.g., "03e_habit")
#' @param output_dir     Character -- path to write the report
#' @param mode_result    List: $theta_mode (named vec), $logpost (numeric)
#' @param prior_spec     Data frame: name, distribution, p1, p2, lower, upper
#' @param d5,d6,d8,d9,d10  dynhr_diagnostic objects (NULL = skip section)
#' @param model_meta     List: $n_endo, $n_exo, $n_obs, $n_params, $obs_names, $n_eqs
#' @param chain_stats    Data frame from MCMC parallel output
#' @param parent_model   Character -- name of parent model for delta tracking
#' @param parent_results List: $logpost, $d9_ratios (named vec), etc.
#' @param notes          Character vector -- free-form notes
#' @return Invisible path to written markdown file
#' @noRd
write_llm_report <- function(model_name,
                             output_dir,
                             mode_result    = NULL,
                             prior_spec     = NULL,
                             d5 = NULL, d6 = NULL,
                             d8 = NULL, d9 = NULL,
                             d10 = NULL,
                             model_meta     = NULL,
                             data_summary   = NULL,
                             chain_stats    = NULL,
                             benchmarks     = NULL,
                             parent_model   = NULL,
                             parent_results = NULL,
                             notes          = NULL) {

  out_path <- file.path(output_dir, paste0(model_name, "_llm_report.md"))
  L <- character(0)
  add   <- function(...) L <<- c(L, paste0(...))
  blank <- function()    L <<- c(L, "")

  ## -- Header --------------------------------------------------------
  add("# Diagnostic Report: ", model_name)
  add("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
  add("Format: LLM-readable model diagnostic summary")
  blank()

  ## -- Model metadata ------------------------------------------------
  add("## Model Structure")
  if (!is.null(model_meta)) {
    add("- Endogenous variables: ", model_meta$n_endo)
    add("- Exogenous shocks: ", model_meta$n_exo)
    add("- Observables: ", model_meta$n_obs,
        " [", paste(model_meta$obs_names, collapse = ", "), "]")
    add("- Estimated parameters: ", model_meta$n_params)
    add("- Equations: ", model_meta$n_eqs)
  }
  if (!is.null(parent_model))
    add("- Parent model: ", parent_model)
  blank()

  ## -- Dashboard -----------------------------------------------------
  add("## Diagnostic Dashboard")
  blank()
  add("| Diagnostic | Status | Key metric |")
  add("|---|---|---|")

  .st <- function(dx) {
    if (is.null(dx))           return("SKIP")
    if (isTRUE(dx$errored))    return("ERROR")
    if (is.na(dx$pass))        return("INFO")
    if (isTRUE(dx$pass))       return("PASS")
    "FAIL"
  }

  ## D5
  d5m <- ""
  if (!is.null(d5) && !is.null(d5$result$ess)) {
    min_e <- min(d5$result$ess, na.rm = TRUE)
    max_r <- if (!is.null(d5$result$rhat))
      max(d5$result$rhat$rhat, na.rm = TRUE) else NA
    d5m <- sprintf("min ESS=%d, max Rhat=%.3f", round(min_e),
                   if (is.na(max_r)) NA_real_ else max_r)
  }
  add("| D5 MCMC Convergence | ", .st(d5), " | ", d5m, " |")

  ## D6
  d6m <- ""
  if (!is.null(d6) && !is.null(d6$result$overlap_scores)) {
    ov <- d6$result$overlap_scores
    d6m <- sprintf("overlap [%.2f, %.2f]", min(ov), max(ov))
  }
  add("| D6 Prior vs Posterior | ", .st(d6), " | ", d6m, " |")

  ## D8
  d8m <- ""
  if (!is.null(d8) && !is.null(d8$result$checks)) {
    np <- sum(sapply(d8$result$checks, function(x) x$status == "PASS"))
    d8m <- sprintf("%d/%d benchmarks", np, length(d8$result$checks))
  }
  add("| D8 IRF Plausibility | ", .st(d8), " | ", d8m, " |")

  ## D9
  d9m <- ""
  if (!is.null(d9) && !is.null(d9$result$sd_table)) {
    tbl <- d9$result$sd_table
    rc  <- intersect(c("ratio", "model_data_ratio"), names(tbl))[1]
    if (!is.na(rc)) {
      np <- sum(tbl[[rc]] >= 0.5 & tbl[[rc]] <= 2.0, na.rm = TRUE)
      d9m <- sprintf("%d/%d obs in [0.5, 2.0]", np, nrow(tbl))
    }
  }
  add("| D9 Moment Matching | ", .st(d9), " | ", d9m, " |")

  ## D10
  add("| D10 Variance Decomp | ", .st(d10), " | |")
  blank()

  ## -- Delta from parent ---------------------------------------------
  if (!is.null(parent_results)) {
    add("## Delta from Parent (", parent_model, ")")
    blank()
    add("| Metric | ", parent_model, " | ", model_name, " | Change |")
    add("|---|---|---|---|")
    ## A dynhr_mode_result carries the mode log-posterior at $mode$logpost, NOT
    ## at the top level; reading mode_result$logpost returned NULL, so this delta
    ## row silently vanished (NULL - x -> numeric(0)). Read the nested value,
    ## keeping a top-level fallback for inner-optimiser-list inputs.
    mlp <- mode_result$logpost %||% mode_result$mode$logpost
    if (!is.null(parent_results$logpost) && !is.null(mlp)) {
      d <- mlp - parent_results$logpost
      add(sprintf("| Log-posterior (mode) | %.1f | %.1f | %+.1f |",
                  parent_results$logpost, mlp, d))
    }
    blank()
  }

  ## -- Mode estimates ------------------------------------------------
  if (!is.null(mode_result) && !is.null(prior_spec)) {
    add("## Posterior Mode Estimates")
    blank()
    add("| Parameter | Prior | Prior mean | Mode | Delta | Flag |")
    add("|---|---|---|---|---|---|")

    ## Fallback must be the nested theta vector ($mode$theta_mode), NOT the whole
    ## inner $mode list (which is not a named parameter vector).
    theta <- mode_result$theta_mode %||% mode_result$mode$theta_mode

    for (i in seq_len(nrow(prior_spec))) {
      nm  <- prior_spec$name[i]
      lo  <- prior_spec$lower[i]
      hi  <- prior_spec$upper[i]
      pm  <- prior_spec$p1[i]
      psd <- prior_spec$p2[i]
      val <- if (nm %in% names(theta)) as.numeric(theta[[nm]]) else NA_real_

      flag <- ""
      if (!is.na(val) && !is.na(lo) && !is.na(hi)) {
        rng <- hi - lo
        if (val < lo + 0.02 * rng)       flag <- "AT LOWER BOUND"
        else if (val > hi - 0.02 * rng)  flag <- "AT UPPER BOUND"
        else if (abs(val - pm) > 2 * psd) flag <- ">2sd from prior"
      }

      add(sprintf("| %s | %s(%.4f, %.4f) | %.4f | %.6f | %+.4f | %s |",
                  nm, prior_spec$distribution[i], pm, psd,
                  pm, val, val - pm, flag))
    }
    blank()
  }

  ## -- D5 detail -----------------------------------------------------
  if (!is.null(d5) && !is.null(d5$result$ess)) {
    add("## D5 MCMC Convergence Detail")
    blank()
    add("| Parameter | ESS | R-hat | Flag |")
    add("|---|---|---|---|")
    ess <- d5$result$ess
    rv  <- if (!is.null(d5$result$rhat)) d5$result$rhat$rhat else NULL
    for (nm in names(ess)) {
      e <- round(ess[[nm]])
      r <- if (!is.null(rv) && nm %in% names(rv)) sprintf("%.3f", rv[[nm]]) else "--"
      fl <- ""
      if (e < 400) fl <- "LOW ESS"
      if (!is.null(rv) && nm %in% names(rv) && rv[[nm]] > 1.05)
        fl <- paste(fl, "HIGH RHAT")
      add(sprintf("| %s | %d | %s | %s |", nm, e, r, trimws(fl)))
    }
    blank()
  }

  ## -- D9 detail -----------------------------------------------------
  if (!is.null(d9) && !is.null(d9$result$sd_table)) {
    add("## D9 Moment Matching Detail")
    blank()
    tbl <- d9$result$sd_table
    add("| Observable | Data SD | Model SD | Ratio | Status |")
    add("|---|---|---|---|---|")
    rc <- intersect(c("ratio", "model_data_ratio"), names(tbl))[1]
    for (i in seq_len(nrow(tbl))) {
      r   <- tbl[i, ]
      rat <- if (!is.na(rc)) r[[rc]] else NA
      st  <- if (!is.na(rat) && rat >= 0.5 && rat <= 2.0) "PASS" else "FAIL"
      vn  <- r$observable %||% r$variable %||% rownames(tbl)[i]
      dsd <- r$data_sd %||% r$sd_data %||% NA
      msd <- r$model_sd %||% r$sd_model %||% NA
      add(sprintf("| %s | %.5f | %.5f | %.2fx | %s |", vn, dsd, msd, rat, st))
    }
    blank()
  }

  ## -- D10 detail ----------------------------------------------------
  if (!is.null(d10) && !is.null(d10$result$vd_long)) {
    add("## D10 Variance Decomposition (%, unconditional)")
    blank()
    vd <- d10$result$vd_long
    vd_w <- reshape(vd[, c("variable","shock","share")],
                    idvar = "variable", timevar = "shock", direction = "wide")
    names(vd_w) <- sub("^share\\.", "", names(vd_w))
    sc <- setdiff(names(vd_w), "variable")
    add(paste(c("| Variable", sc, "|"), collapse = " | "))
    add(paste(c("|---", rep("---", length(sc)), "|"), collapse = "|"))
    for (i in seq_len(nrow(vd_w))) {
      vals <- sprintf("%.1f", as.numeric(vd_w[i, sc]) * 100)
      add(paste(c("|", vd_w$variable[i], vals, "|"), collapse = " | "))
    }
    blank()
  }

  ## -- D8 detail -----------------------------------------------------
  if (!is.null(d8) && !is.null(d8$result$checks)) {
    add("## D8 IRF Benchmark Checks")
    blank()
    add("| Benchmark | Status | Shock | Variable | Peak | Horizon |")
    add("|---|---|---|---|---|---|")
    for (chk in d8$result$checks) {
      add(sprintf("| %s | %s | %s | %s | %.4f | %d |",
                  chk$benchmark, chk$status,
                  chk$shock %||% "--", chk$variable %||% "--",
                  chk$peak_value %||% NA_real_, chk$peak_horizon %||% NA_integer_))
    }
    blank()
  }

  ## -- Chain stats ---------------------------------------------------
  if (!is.null(chain_stats)) {
    add("## MCMC Chain Summary")
    blank()
    add("```")
    add(paste(capture.output(print(chain_stats, row.names = FALSE, digits = 3)),
              collapse = "\n"))
    add("```")
    blank()
  }

  ## -- Auto-detected issues ------------------------------------------
  add("## Auto-detected Issues")
  blank()
  issues <- character(0)

  # Boundaries
  if (!is.null(mode_result) && !is.null(prior_spec)) {
    ## Fallback must be the nested theta vector ($mode$theta_mode), NOT the whole
    ## inner $mode list (which is not a named parameter vector).
    theta <- mode_result$theta_mode %||% mode_result$mode$theta_mode
    for (i in seq_len(nrow(prior_spec))) {
      nm <- prior_spec$name[i]; lo <- prior_spec$lower[i]; hi <- prior_spec$upper[i]
      val <- if (nm %in% names(theta)) as.numeric(theta[[nm]]) else NA
      if (!is.na(val) && !is.na(lo) && !is.na(hi)) {
        rng <- hi - lo
        if (val < lo + 0.02 * rng)
          issues <- c(issues, sprintf(
            "- **%s = %.4f** at lower bound (%.4f). Consider widening or reparameterising.", nm, val, lo))
        if (val > hi - 0.02 * rng)
          issues <- c(issues, sprintf(
            "- **%s = %.4f** at upper bound (%.4f). Consider widening or reparameterising.", nm, val, hi))
      }
    }
  }

  # Low ESS
  if (!is.null(d5) && !is.null(d5$result$ess))
    for (nm in names(d5$result$ess))
      if (d5$result$ess[[nm]] < 400)
        issues <- c(issues, sprintf(
          "- **%s** ESS = %d (< 400). Poor mixing.", nm, round(d5$result$ess[[nm]])))

  # High R-hat
  if (!is.null(d5) && !is.null(d5$result$rhat))
    for (nm in names(d5$result$rhat$rhat))
      if (d5$result$rhat$rhat[[nm]] > 1.05)
        issues <- c(issues, sprintf(
          "- **%s** R-hat = %.3f (> 1.05). Chains not converged.", nm, d5$result$rhat$rhat[[nm]]))

  # D9 failures
  if (!is.null(d9) && !is.null(d9$result$sd_table)) {
    tbl <- d9$result$sd_table
    rc  <- intersect(c("ratio", "model_data_ratio"), names(tbl))[1]
    if (!is.na(rc))
      for (i in seq_len(nrow(tbl))) {
        r <- tbl[[rc]][i]
        v <- tbl$observable[i] %||% tbl$variable[i] %||% rownames(tbl)[i]
        if (!is.na(r) && (r < 0.5 || r > 2.0))
          issues <- c(issues, sprintf(
            "- **%s** model/data SD ratio = %.2fx (outside [0.5, 2.0]).", v, r))
      }
  }

  if (length(issues) == 0) add("No issues auto-detected.")
  else for (iss in issues) add(iss)
  blank()

  ## -- Notes ---------------------------------------------------------
  if (!is.null(notes) && length(notes) > 0) {
    add("## Session Notes")
    for (n in notes) add(n)
    blank()
  }

  ## -- Prompt template -----------------------------------------------
  add("## Suggested LLM Prompt")
  blank()
  add("```")
  add("I am iterating on a DSGE model (", model_name, ") estimated with dynhr.")
  add("The diagnostic report above shows the current state.")
  add("Key issues to address:")
  if (length(issues) > 0) for (iss in issues) add("  ", iss)
  else add("  (none auto-detected -- review D9/D10 for economic plausibility)")
  add("Please suggest targeted model changes to fix these issues,")
  add("keeping the model structure minimal and well-identified.")
  add("```")
  blank()

  ## -- Write ---------------------------------------------------------
  writeLines(L, out_path)
  message(sprintf("[dynhr] LLM report: %s (%d lines)", out_path, length(L)))
  invisible(out_path)
}


# ---------------------------------------------------------------------------
#' Render a dynhr diagnostic report
#'
#' Dispatches on \code{format}:
#' \describe{
#'   \item{\code{"llm"}}{Writes a compact structured plain-text Markdown file
#'     using the \code{$llm_summary} fields from every diagnostic.  No
#'     external dependencies.}
#'   \item{\code{"html"}}{Renders the package's Quarto template
#'     (\code{inst/templates/report.qmd}) to a self-contained HTML file with
#'     tabsets, callout blocks, and (optionally) interactive plotly IRFs.
#'     Requires the \code{quarto} R package and a Quarto CLI installation.}
#'   \item{\code{"pdf"}}{Renders the package's Typst template
#'     (\code{inst/templates/report-pdf.qmd}) to PDF via Quarto's Typst
#'     engine.  No LaTeX required.  Requires the \code{quarto} R package
#'     and a Quarto CLI (>= 1.4) installation.}
#' }
#'
#' @param diagnostics A \code{dynhr_diagnostic_suite} (named list of
#'   \code{dynhr_diagnostic} objects, as returned by
#'   \code{\link{run_diagnostics}}).
#' @param file        Output file path.  Default: \code{"dynhr_report.md"}
#'   for \code{format = "llm"}, \code{"dynhr_report.html"} for \code{"html"}.
#' @param format      One of \code{"llm"} (default), \code{"html"}, or \code{"pdf"}.
#' @param model_name  Optional character label embedded in the report header.
#' @param irfs        Optional named list of IRF matrices for interactive
#'   plots in the HTML report.
#' @param ...         Additional arguments passed to the underlying writer.
#' @return \code{file}, invisibly.
#' @export
write_report <- function(diagnostics,
                          file       = NULL,
                          format     = "llm",
                          model_name = NULL,
                          irfs       = NULL,
                          ...) {
  format <- match.arg(format, c("llm", "html", "pdf"))

  # Auto file extension
  if (is.null(file)) {
    file <- switch(format,
                   html = "dynhr_report.html",
                   pdf  = "dynhr_report.pdf",
                   "dynhr_report.md")
  }

  if (format == "llm") {
    # Use the new suite-based LLM formatter if we have a diagnostic suite;
    # fall back to the legacy component-based writer for old-style calls.
    if (inherits(diagnostics, "dynhr_diagnostic_suite") ||
        (is.list(diagnostics) &&
         any(vapply(diagnostics,
                    function(x) inherits(x, "dynhr_diagnostic"),
                    logical(1))))) {
      txt <- format_llm_report(diagnostics, model_name = model_name)
      writeLines(txt, file)
      message("[dynhr] LLM report written to ", file)
      return(invisible(file))
    }
    # Legacy path
    write_llm_report(diagnostics, file, ...)
    return(invisible(file))
  }

  # ---- HTML or PDF via Quarto ----
  if (!requireNamespace("quarto", quietly = TRUE))
    stop(paste(
      "The 'quarto' R package is required for HTML/PDF output.",
      "Install it with: install.packages('quarto').",
      "Also ensure the Quarto CLI is installed (https://quarto.org)."
    ), call. = FALSE)

  is_pdf      <- identical(format, "pdf")
  tpl_name    <- if (is_pdf) "report-pdf.qmd" else "report.qmd"
  quarto_fmt  <- if (is_pdf) "typst" else "html"
  out_ext     <- if (is_pdf) "pdf" else "html"
  fmt_label   <- toupper(format)

  tpl <- system.file("templates", tpl_name, package = "dynhr")
  if (!file.exists(tpl))
    stop(sprintf("Quarto template not found at inst/templates/%s", tpl_name),
         call. = FALSE)

  # Write report to a temp directory, then move to requested path.
  # The quarto R package YAML-serialises execute_params, so it can't
  # transport ggplots / NA / nested lists.  We spool the suite (and IRFs)
  # to RDS and pass file paths as string parameters instead.
  tmp_dir  <- tempfile("dynhr_report_")
  dir.create(tmp_dir)
  tmp_qmd  <- file.path(tmp_dir, tpl_name)
  file.copy(tpl, tmp_qmd)

  diags_rds <- file.path(tmp_dir, "diags.rds")
  saveRDS(diagnostics, diags_rds)
  irfs_rds <- ""
  if (!is.null(irfs) && length(irfs) > 0L) {
    # Pre-render an IRF ggplot per shock so the Quarto subprocess (which
    # may not have dynhr installed) doesn't need to call theme_dynhr_*.
    plots <- lapply(names(irfs), function(sh) {
      mat <- as.data.frame(irfs[[sh]])
      if (nrow(mat) == 0L || ncol(mat) == 0L) return(NULL)
      mat$horizon <- seq_len(nrow(mat))
      vars <- setdiff(names(mat), "horizon")
      long <- do.call(rbind, lapply(vars, function(v)
        data.frame(horizon = mat$horizon, variable = v,
                   value = mat[[v]], stringsAsFactors = FALSE)
      ))
      long$variable <- factor(long$variable, levels = vars)
      p <- ggplot2::ggplot(long,
              ggplot2::aes(x = horizon, y = value)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                            linewidth = 0.25) +
        ggplot2::geom_line(colour = dynhr_primary_colour,
                           linewidth = 0.5) +
        ggplot2::facet_wrap(~ variable, scales = "free_y") +
        theme_dynhr() +
        ggplot2::labs(title = sprintf("IRF: shock = %s", sh),
                      x = "Horizon (quarters)", y = "Response")
      attr(p, "dynhr_fig_height") <- 7.5
      p
    })
    names(plots) <- names(irfs)
    plots <- plots[!vapply(plots, is.null, logical(1))]
    irfs_rds <- file.path(tmp_dir, "irfs.rds")
    saveRDS(plots, irfs_rds)
  }

  # Pre-render the LLM summary AND the executive summary here so the template
  # doesn't need to call into the dynhr package (Quarto runs in a fresh R
  # subprocess that may only have a stale *installed* dynhr, not the
  # devtools::load_all() development copy).
  llm_txt <- file.path(tmp_dir, "llm.txt")
  writeLines(format_llm_report(diagnostics,
                                model_name = model_name %||% "model"),
             llm_txt)

  exec_txt <- file.path(tmp_dir, "exec.txt")
  writeLines(format_executive_summary(diagnostics,
                                       model_name = model_name %||% "model",
                                       unicode    = FALSE),
             exec_txt)

  render_args <- list(
    input          = tmp_qmd,
    output_format  = quarto_fmt,
    execute_params = list(
      diags_rds  = diags_rds,
      irfs_rds   = irfs_rds,
      llm_txt    = llm_txt,
      exec_txt   = exec_txt,
      model_name = model_name %||% "model"
    ),
    output_file    = basename(file),
    quiet          = !isTRUE(getOption("dynhr.report.verbose", FALSE))
    )
    do.call(quarto::quarto_render, render_args)

    out_path <- file.path(tmp_dir, basename(file))
    if (!file.exists(out_path)) {
      # quarto may emit using the qmd basename + format extension
      out_path <- file.path(tmp_dir,
                            sub("\\.qmd$", paste0(".", out_ext),
                                basename(tmp_qmd)))
    }
    if (!file.exists(out_path))
      stop("Quarto produced no output file in ", tmp_dir)

    file.copy(out_path, file, overwrite = TRUE)
    message("[dynhr] ", fmt_label, " report written to ", file)

  invisible(file)
}
