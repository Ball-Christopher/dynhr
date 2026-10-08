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

  ## Single source of truth for the badge (including the 0.9.4 WARN level);
  ## "SKIP" is this report's own extra state for a diagnostic that never ran.
  .st <- function(dx) {
    if (is.null(dx)) return("SKIP")
    .badge_str(dx)
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
      ## Read D9's own per-row verdict when present (ok / marginal / outside /
      ## non-finite) instead of recomputing it from the ratio.
      np <- if (!is.null(tbl$status)) sum(tbl$status %in% c("ok", "marginal"))
            else sum(tbl[[rc]] >= 0.5 & tbl[[rc]] <= 2.0, na.rm = TRUE)
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
      st  <- if (!is.null(r$status)) toupper(r$status)
             else if (!is.na(rat) && rat >= 0.5 && rat <= 2.0) "PASS" else "FAIL"
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
        outside <- if (!is.null(tbl$status)) identical(tbl$status[i], "outside")
                   else (!is.na(r) && (r < 0.5 || r > 2.0))
        if (outside)
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
  .dynhr_inform(sprintf("[dynhr] LLM report: %s (%d lines)", out_path, length(L)))
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
#'     tabsets, callout blocks and pre-rendered IRF facet plots.  Requires the
#'     \code{quarto} R package and a Quarto CLI installation.}
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
#' @param irfs        Optional IRF input for the report's IRF section: a named
#'   list of \code{horizon x variable} matrices, or a \code{solve_model()} /
#'   \code{stoch_simul()} result, whose \code{$irfs} element is unwrapped.
#'   Entries that are not matrices are dropped with a warning naming them.
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
      .dynhr_inform("[dynhr] LLM report written to ", file)
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

  ## Normalise the output extension UP FRONT. quarto_render() writes
  ## `output_file` with the format's own extension appended when the two
  ## disagree, so write_report(file = "x.html", format = "pdf") produced
  ## "x.html.pdf" in the temp dir and then aborted with "Quarto produced no
  ## output file", leaving nothing at `file`.
  cur_ext <- tolower(tools::file_ext(file))
  if (!identical(cur_ext, out_ext)) {
    old_file <- file
    file <- if (nzchar(cur_ext)) sub(paste0("\\.", cur_ext, "$"),
                                     paste0(".", out_ext), file)
            else paste0(file, ".", out_ext)
    .dynhr_warn("write_report(): format = \"", format,
                "\" writes a .", out_ext, " file; output path changed from \"",
                old_file, "\" to \"", file, "\".")
  }

  ## `irfs` may be a solve_model()/stoch_simul() result rather than the bare
  ## named list of matrices -- that is the most natural call, and D8 already
  ## unwraps it this way (.get_irfs_long()).
  if (!is.null(irfs) && !is.null(irfs[["irfs", exact = TRUE]]))
    irfs <- irfs[["irfs", exact = TRUE]]
  if (!is.null(irfs) && length(irfs) > 0L) {
    keep <- vapply(irfs, is.matrix, logical(1))
    if (any(!keep))
      .dynhr_warn("write_report(): dropping non-matrix irfs entr",
                  if (sum(!keep) == 1L) "y: " else "ies: ",
                  paste(names(irfs)[!keep] %||%
                          which(!keep), collapse = ", "))
    irfs <- irfs[keep]
  }

  # Write report to a temp directory, then move to requested path.
  # The quarto R package YAML-serialises execute_params, so it can't
  # transport ggplots / NA / nested lists.  We spool the suite (and IRFs)
  # to RDS and pass file paths as string parameters instead.
  tmp_dir  <- tempfile("dynhr_report_")
  dir.create(tmp_dir)
  ## Every render used to leak a multi-MB directory for the rest of the
  ## session (and a FAILED render leaked one too). Keep it only when the user
  ## explicitly asks to inspect the intermediates.
  if (!isTRUE(getOption("dynhr.report.keep_tmp", FALSE)))
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  tmp_qmd  <- file.path(tmp_dir, tpl_name)
  file.copy(tpl, tmp_qmd)

  diags_rds <- file.path(tmp_dir, "diags.rds")
  saveRDS(diagnostics, diags_rds)

  ## The writer/template contract: classification, ordering, colours,
  ## explanations, actions and provenance are computed HERE and read ONLY
  ## from this file by the templates.
  rmeta <- .report_meta(diagnostics, model_name = model_name, irfs = irfs)
  meta_rds <- file.path(tmp_dir, "report-meta.rds")
  saveRDS(rmeta, meta_rds)

  ## dynhr is not necessarily installed in the Quarto subprocess, so the one
  ## escaping helper travels as source.
  escape_r <- file.path(tmp_dir, "escape.R")
  writeLines(.report_escape_source(), escape_r)

  ## Same for the figure height / pagination rule: ONE copy, in the package.
  plotutil_r <- file.path(tmp_dir, "plotutil.R")
  writeLines(.report_plotutil_source(), plotutil_r)

  irfs_rds <- ""
  if (!is.null(irfs) && length(irfs) > 0L) {
    # Pre-render an IRF ggplot per shock so the Quarto subprocess (which
    # may not have dynhr installed) doesn't need to call theme_dynhr_*.
    plots <- .report_irf_plots(
      irfs, horizon_lab = .irf_horizon_label(rmeta$provenance$frequency))
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
      meta_rds   = meta_rds,
      escape_r   = escape_r,
      plotutil_r = plotutil_r,
      ## Escaped here, at the boundary: an unescaped model name reached the
      ## Pandoc title and "model_*a*_$x$_<v1>" rendered as "model_a_x_".
      model_name = rmeta$model_name
    ),
    output_file    = basename(file),
    quiet          = !isTRUE(getOption("dynhr.report.verbose", FALSE))
    )
    do.call(quarto::quarto_render, render_args)

    out_path <- file.path(tmp_dir, basename(file))
    if (!file.exists(out_path)) {
      ## Fall back to whatever quarto actually emitted with the right
      ## extension -- a basename-derived guess only works when the template
      ## name and the requested name happen to agree.
      cand <- list.files(tmp_dir, pattern = paste0("\\.", out_ext, "$"),
                         full.names = TRUE)
      if (length(cand) > 0L) out_path <- cand[[1L]]
    }
    if (!file.exists(out_path))
      stop("Quarto produced no output file in ", tmp_dir)

    file.copy(out_path, file, overwrite = TRUE)
    .dynhr_inform("[dynhr] ", fmt_label, " report written to ", file)

  invisible(file)
}


# ---------------------------------------------------------------------------
# report-meta.rds: the writer/template contract
# ---------------------------------------------------------------------------
# Everything classification-, ordering-, colour-, explanation- or provenance-
# related is computed HERE, in the process that holds the loaded dynhr, and
# serialised into the tmp dir next to diags.rds. The Quarto templates read only
# this list; they never touch `$pass` / `$warn` / `$errored` themselves.
#
# The reason is the 2026-09-17 report review's blocker A1: report-pdf.qmd had
# re-implemented the badge classification as `sum(isTRUE(r$pass))`, which
# counts a WARN (encoded as pass = TRUE + warn = TRUE) as a PASS. The same PDF
# printed "9 PASS | 2 FAIL" in its summary table and "8 PASS | 1 WARN | 2 FAIL"
# in the executive summary two pages earlier. A classification that exists in
# two places will disagree; this one exists in one.

#' Human-readable sample span of the estimation data
#' @noRd
.diag_sample_span <- function(data, dates = NULL) {
  lab <- function(v) {
    if (length(v) == 0L) return(NA_character_)
    paste(format(v[[1L]]), "to", format(v[[length(v)]]))
  }
  if (!is.null(dates) && length(dates) > 0L) return(lab(dates))
  if (is.null(data)) return(NA_character_)
  rn <- rownames(data)
  if (!is.null(rn) && length(rn) > 0L) return(lab(rn))
  NA_character_
}

#' Data frequency implied by a vector of dates ("quarterly"/"monthly"/...)
#' @noRd
.diag_frequency <- function(dates) {
  if (is.null(dates) || length(dates) < 3L) return(NA_character_)
  d <- suppressWarnings(as.numeric(diff(as.Date(dates))))
  d <- d[is.finite(d)]
  if (length(d) == 0L) return(NA_character_)
  med <- stats::median(d)
  if      (med >= 0.5   && med <= 1.5)   "daily"
  else if (med >= 6     && med <= 8)     "weekly"
  else if (med >= 27    && med <= 32)    "monthly"
  else if (med >= 88    && med <= 94)    "quarterly"
  else if (med >= 175   && med <= 190)   "semiannual"
  else if (med >= 355   && med <= 372)   "annual"
  else NA_character_
}

#' Axis label for the IRF horizon, given a frequency
#' @noRd
.irf_horizon_label <- function(freq) {
  if (length(freq) != 1L || is.na(freq)) return("Horizon")
  unit <- switch(as.character(freq),
                 daily      = "days",
                 weekly     = "weeks",
                 monthly    = "months",
                 quarterly  = "quarters",
                 semiannual = "half-years",
                 annual     = "years",
                 NULL)
  if (is.null(unit)) "Horizon" else sprintf("Horizon (%s)", unit)
}

.REPORT_IRF_UNITS <- paste0(
  "one-standard-deviation shock; responses in the model's own units ",
  "(deviations from steady state)")

.REPORT_LEGEND <- paste0(
  "PASS ok | WARN soft failure, pass = TRUE with a warning flag | ",
  "FAIL gate failed | ERROR did not run | ",
  "INFO reported, not gated")

#' Build the report metadata list consumed by the Quarto templates
#'
#' @param diagnostics A \code{dynhr_diagnostic_suite} (or any named list of
#'   \code{dynhr_diagnostic} objects; non-diagnostic entries are ignored).
#' @param model_name  Optional character label. Escaped for the Pandoc title.
#' @param irfs        Optional named list of IRF matrices (used only to record
#'   which shocks the report covers).
#' @param ...         Reserved.
#' @return A named list with exactly the fields in the RF1/RF2 contract:
#'   \code{badges}, \code{counts}, \code{colours}, \code{order},
#'   \code{disp_id}, \code{explanations}, \code{actions}, \code{provenance},
#'   \code{model_name}, \code{irf_units}, \code{legend}.
#' @noRd
.report_meta <- function(diagnostics, model_name = NULL, irfs = NULL, ...) {
  diag_items <- Filter(function(r) inherits(r, "dynhr_diagnostic"), diagnostics)
  nms <- names(diag_items) %||% character(0)

  badges <- if (length(nms) == 0L) stats::setNames(character(0), character(0))
            else vapply(diag_items, .badge_str, character(1))

  cnt <- .badge_counts(diag_items)
  ## The contract's `counts` sums to the number of diagnostics, so the
  ## `total` element .badge_counts() carries is dropped rather than left in
  ## to be double-counted by a template that sums the vector.
  counts <- c(PASS = cnt[["PASS"]], WARN = cnt[["WARN"]], FAIL = cnt[["FAIL"]],
              ERROR = cnt[["ERROR"]], INFO = cnt[["INFO"]])
  storage.mode(counts) <- "integer"

  lv <- c("PASS", "WARN", "FAIL", "ERROR", "INFO")
  colours <- stats::setNames(vapply(lv, .badge_colour, character(1)), lv)

  prov_in <- attr(diagnostics, "provenance")
  provenance <- list(
    dynhr_version = as.character(utils::packageVersion("dynhr")),
    git_commit    = .dynhr_git_stamp_at(getNamespaceInfo(asNamespace("dynhr"), "path")),
    model_name    = prov_in$model_name  %||% model_name %||% NA_character_,
    model_file    = prov_in$model_file  %||% NA_character_,
    n_obs         = prov_in$n_obs       %||% NA_integer_,
    obs_vars      = prov_in$obs_vars    %||% NA_character_,
    data_hash     = prov_in$data_hash   %||% NA_character_,
    sample_span   = prov_in$sample_span %||% NA_character_,
    frequency     = prov_in$frequency   %||% NA_character_,
    seed          = prov_in$seed        %||% NA,
    sampler       = prov_in$sampler     %||% NA_character_,
    n_draws       = prov_in$n_draws     %||% NA_integer_,
    n_warmup      = prov_in$n_warmup    %||% NA,
    n_chains      = prov_in$n_chains    %||% NA_integer_,
    n_shocks_irf  = if (is.null(irfs)) NA_integer_ else length(irfs),
    rendered_at   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )

  list(
    badges       = badges,
    counts       = counts,
    colours      = colours,
    order        = .diag_order(diag_items),
    disp_id      = stats::setNames(
                     vapply(nms, function(n) .diag_meta_for(n)$disp_id,
                            character(1)), nms),
    explanations = stats::setNames(lapply(nms, .diag_explanation_for), nms),
    actions      = stats::setNames(
                     vapply(nms, function(n)
                              .diag_action_for(n, badges[[n]]),
                            character(1)), nms),
    provenance   = provenance,
    model_name   = .escape_for(model_name %||% "model", target = "title"),
    irf_units    = .REPORT_IRF_UNITS,
    legend       = .REPORT_LEGEND
  )
}


#' Source text of the escaping helper, for the Quarto subprocess
#'
#' dynhr may not be installed in the environment Quarto spawns, so the
#' templates cannot call \code{dynhr:::.escape_for()}. The writer drops this
#' file beside the qmd and the templates \code{source()} it.
#' @noRd
.report_escape_source <- function() {
  c("## Written by dynhr:::write_report(); do not edit.",
    "## Source of .escape_for() -- the ONLY way text may reach a sink in the",
    "## report templates. It depends on nothing outside base R.",
    paste0(".escape_for <- ", paste(deparse(.escape_for), collapse = "\n")),
    "")
}


#' Source text of the figure-sizing / pagination helpers, for the Quarto
#' subprocess
#'
#' Same reason as \code{.report_escape_source()}: dynhr may not be installed
#' where Quarto runs, so the helpers travel as source rather than being
#' called as \code{dynhr:::}. Both templates \code{source()} this file and
#' neither keeps a second copy of the height rule.
#' @noRd
.report_plotutil_source <- function() {
  fns <- c(".fig_height_for", ".facet_var_names", ".aes_var_name",
           ".paginate_key", ".paginate_clone_layer", ".paginate_page",
           ".paginate_finish", ".paginate_plot", ".paginate_facets",
           ".paginate_discrete")
  c("## Written by dynhr:::write_report(); do not edit.",
    "## Figure height rule + many-panel pagination for the report templates.",
    "## Depends on base R and ggplot2 only.",
    unlist(lapply(fns, function(nm)
      paste0(nm, " <- ", paste(deparse(get(nm)), collapse = "\n")))),
    "")
}


#' Pre-render one IRF facet plot per shock
#'
#' @param irfs        Named list of \code{horizon x variable} matrices.
#' @param horizon_lab X-axis label (see \code{.irf_horizon_label}).
#' @return Named list of ggplots, one per shock with usable data.
#' @noRd
.report_irf_plots <- function(irfs, horizon_lab = "Horizon") {
  plots <- lapply(names(irfs), function(sh) {
    mat <- as.data.frame(irfs[[sh]])
    if (nrow(mat) == 0L || ncol(mat) == 0L) return(NULL)
    mat$horizon <- seq_len(nrow(mat))
    vars <- setdiff(names(mat), "horizon")
    long <- do.call(rbind, lapply(vars, function(v)
      data.frame(horizon = mat$horizon, variable = v,
                 value = mat[[v]], stringsAsFactors = FALSE)
    ))

    ## scales = "free_y" with no zero tolerance draws a response of size
    ## 1e-17 as a full-amplitude hump whose 7-significant-digit scientific
    ## axis labels ("8.237751e-18 ... -1.729928e-16") blow out the facet
    ## width and squeeze every neighbouring panel. Flag those facets and give
    ## them a symmetric, human-sized range instead.
    fin  <- long$value[is.finite(long$value)]
    gmax <- if (length(fin)) max(abs(fin)) else 0
    tol  <- 1e-12 * gmax
    vmax <- vapply(split(long$value, long$variable),
                   function(v) {
                     v <- v[is.finite(v)]
                     if (length(v)) max(abs(v)) else 0
                   }, numeric(1))
    zero_vars <- names(vmax)[gmax > 0 & vmax <= tol]

    lab_of <- function(v)
      ifelse(v %in% zero_vars, paste0(v, " (numerically zero)"), v)
    long$variable <- factor(lab_of(long$variable), levels = lab_of(vars))

    p <- ggplot2::ggplot(long, ggplot2::aes(x = horizon, y = value)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                          linewidth = 0.25) +
      ggplot2::geom_line(colour = dynhr_primary_colour, linewidth = 0.5)
    if (length(zero_vars) > 0L) {
      ## Fixed symmetric limits for the flagged facets ONLY: geom_blank
      ## widens their scale without touching the data or the other panels.
      span  <- 1e-3 * gmax
      blank <- data.frame(
        horizon  = rep(long$horizon[1L], 2L * length(zero_vars)),
        variable = factor(rep(lab_of(zero_vars), each = 2L),
                          levels = levels(long$variable)),
        value    = rep(c(-span, span), times = length(zero_vars)))
      p <- p + ggplot2::geom_blank(data = blank)
    }
    p <- p +
      ggplot2::facet_wrap(~ variable, scales = "free_y") +
      ggplot2::scale_y_continuous(
        labels = function(v) formatC(v, format = "fg", digits = 3,
                                     drop0trailing = TRUE)) +
      theme_dynhr() +
      ggplot2::labs(
        title    = sprintf("IRF: shock = %s", sh),
        subtitle = if (length(zero_vars) > 0L)
          paste0(.REPORT_IRF_UNITS,
                 ". Panels marked \"numerically zero\" respond by less than ",
                 "1e-12 of the largest response to this shock.")
        else .REPORT_IRF_UNITS,
        x = horizon_lab,
        y = "Response (deviation from steady state)")
    attr(p, "dynhr_fig_height") <- 7.5
    p
  })
  names(plots) <- names(irfs)
  plots[!vapply(plots, is.null, logical(1))]
}
