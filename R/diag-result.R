## R/diag-result.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## .make_result() constructor for dynhr_diagnostic S3 class;
## format_llm_report() suite -> LLM-ready text;
## print methods
##
## LLM summary design
## ------------------
## Every dynhr_diagnostic carries a $llm_summary field: a compact,
## structured text block (~5-15 lines) encoding:
##   - diagnostic name + status
##   - key numeric outcomes as key=value pairs on one line
##   - a single "action:" line with concrete next steps
##
## format_llm_report() concatenates all $llm_summary fields from a suite,
## adds a one-line dashboard, and produces a document that can be pasted
## into an LLM context at ~200-400 tokens total (vs thousands for
## decoding visual output).
## --------------------------------------------------------------------------

#' Standardised result constructor
#'
#' @param result      Raw diagnostic output (list, data.frame, matrix, etc.)
#' @param pass        Logical: TRUE = pass, FALSE = fail, NA = informational only
#' @param plots       Named list of ggplot objects
#' @param summary     Human-readable description string (shown in console)
#' @param llm_summary Compact structured text for LLM consumption (NULL = auto)
#' @param warn        Logical: TRUE marks a SOFT failure. \code{warn = TRUE}
#'   always travels with \code{pass = TRUE} -- see \code{\link{.badge_str}} for
#'   why the WARN level is encoded that way rather than as a third value of
#'   \code{pass}.
#' @return A list with class "dynhr_diagnostic"
#' @noRd
.make_result <- function(result      = NULL,
                         pass        = NA,
                         plots       = list(),
                         summary     = "",
                         llm_summary = NULL,
                         errored     = FALSE,
                         warn        = FALSE) {
  warn <- isTRUE(warn)
  ## A WARN is a pass that wants attention. Encoding it as pass = TRUE keeps
  ## every `identical(pass, FALSE)` gate in the package (notably
  ## run_all_diagnostics(halt_on_identification_fail = TRUE)) from halting on a
  ## soft failure, which is the whole point of the level.
  if (warn) pass <- TRUE
  if (is.null(llm_summary)) {
    badge <- .badge_str(list(pass = pass, errored = errored, warn = warn))
    llm_summary <- sprintf("[%s] %s", badge, summary)
  }
  out <- list(
    result      = result,
    pass        = pass,
    warn        = warn,
    plots       = plots,
    summary     = summary,
    llm_summary = llm_summary,
    errored     = errored
  )
  class(out) <- "dynhr_diagnostic"
  out
}


# ---------------------------------------------------------------------------
#' Format a diagnostic suite as LLM-ready plain text
#'
#' Concatenates the \code{$llm_summary} field from every
#' \code{dynhr_diagnostic} in \code{suite} into a single document optimised
#' for pasting into an LLM context.  The output is intentionally terse:
#' key=value pairs on single lines, no markdown tables, no HTML.
#'
#' @param suite       Named list of class \code{"dynhr_diagnostic_suite"},
#'   or any named list of \code{"dynhr_diagnostic"} objects.
#' @param model_name  Optional character label for the model.
#' @param include_plots Logical: append a "plots available" inventory line
#'   for each diagnostic (default FALSE).
#' @return Character scalar -- the full LLM-ready report.
#' @export
format_llm_report <- function(suite, model_name = NULL, include_plots = FALSE) {
  stopifnot(is.list(suite))

  lines <- character(0)
  add   <- function(...) lines <<- c(lines, paste0(...))
  sep   <- function()    lines <<- c(lines, paste(rep("-", 60), collapse = ""))

  # Header
  add("## dynhr diagnostic report")
  if (!is.null(model_name)) add("model: ", model_name)
  add("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"))

  # Count outcomes
  diag_items <- Filter(function(r) inherits(r, "dynhr_diagnostic"), suite)
  ## Same ordering rule as the executive summary and the rendered report
  ## (severity, then group, then importance) -- see .diag_order().
  diag_items <- diag_items[.diag_order(diag_items)]
  cnt <- .badge_counts(diag_items)
  add("status: ", cnt[["PASS"]], " PASS | ", cnt[["WARN"]], " WARN | ",
      cnt[["FAIL"]], " FAIL | ", cnt[["ERROR"]], " ERROR | ",
      cnt[["INFO"]], " INFO | total=", cnt[["total"]])

  # One-line dashboard
  sep()
  add("## dashboard")
  for (nm in names(diag_items)) {
    r     <- diag_items[[nm]]
    add(sprintf("  [%-5s] %s", .badge_str(r), nm))
  }

  # Full llm_summary per diagnostic
  sep()
  add("## diagnostic details")
  nms <- names(diag_items)
  for (i in seq_along(nms)) {
    nm <- nms[[i]]
    r  <- diag_items[[nm]]
    # Use $llm_summary if present, else fall back to $summary
    txt <- r$llm_summary %||% r$summary %||% "(no summary)"
    add(txt)
    if (include_plots && length(r$plots) > 0)
      add("  plots: ", paste(names(r$plots), collapse = ", "))
    sep()
  }

  paste(lines, collapse = "\n")
}


# ---------------------------------------------------------------------------
#' Write an LLM-ready diagnostic report to file
#'
#' Calls \code{\link{format_llm_report}} and writes the result.  Can also
#' accept a \code{dynhr_diagnostic_suite} as its first argument (new-style),
#' or the legacy positional arguments (old-style, preserved for back-compat).
#'
#' @param diagnostics \code{dynhr_diagnostic_suite} or character model name
#'   (old-style).
#' @param file        Output path (default \code{"dynhr_report.md"}).
#' @param model_name  Optional label override.
#' @param ...         Additional arguments (ignored for suite input; forwarded
#'   for old-style input).
#' @return \code{file}, invisibly.
#' @export
write_llm_suite_report <- function(diagnostics, file = "dynhr_report.md",
                                    model_name = NULL, ...) {
  if (inherits(diagnostics, "dynhr_diagnostic_suite") ||
      (is.list(diagnostics) &&
       any(vapply(diagnostics, function(x) inherits(x, "dynhr_diagnostic"),
                  logical(1))))) {
    txt <- format_llm_report(diagnostics, model_name = model_name)
    writeLines(txt, file)
    .dynhr_inform("[dynhr] LLM report written to ", file)
    return(invisible(file))
  }
  # Legacy path: diagnostics is model_name (character)
  write_llm_report(diagnostics, file, ...)
}


# ---------------------------------------------------------------------------
# Pass / fail badge (HTML) for RMarkdown/Quarto output
# ---------------------------------------------------------------------------

#' Plain-text status badge for a diagnostic result
#'
#' Single source of truth for the ERROR/INFO/PASS/WARN/FAIL classification used
#' by the console print methods, the orchestrator summary, the LLM report, and
#' the executive summary (previously this branch was inlined in 7 places).
#'
#' WARN is the SOFT-failure level added in 0.9.4: "this needs attention" as
#' opposed to FAIL's "do not trust these results".  It is carried by a separate
#' \code{$warn} flag on top of \code{pass = TRUE} rather than as a third value
#' of \code{$pass}, so that code testing \code{identical(pass, FALSE)} (e.g.
#' \code{run_all_diagnostics(halt_on_identification_fail = TRUE)}) never halts
#' on a WARN, and code testing \code{is.na(pass)} never mistakes it for INFO.
#' The split mirrors how Stan and ArviZ separate convergence warnings from
#' hard errors; it is an engineering convention, not a statistical one.
#'
#' @param r A \code{dynhr_diagnostic} (uses \code{$errored}, \code{$pass} and
#'   \code{$warn}).
#' @return Character scalar: "ERROR", "INFO", "PASS", "WARN", or "FAIL".
#' @noRd
.badge_str <- function(r) {
  if (isTRUE(r$errored)) "ERROR"
  else if (is.na(r$pass)) "INFO"
  else if (isTRUE(r$pass)) (if (isTRUE(r$warn)) "WARN" else "PASS")
  else "FAIL"
}

#' Badge counts for a list of diagnostics
#'
#' @param diag_items List of \code{dynhr_diagnostic} objects.
#' @return Named integer vector with elements PASS, WARN, FAIL, ERROR, INFO,
#'   total (in that order).
#' @noRd
.badge_counts <- function(diag_items) {
  b <- vapply(diag_items, .badge_str, character(1))
  c(PASS  = sum(b == "PASS"),  WARN  = sum(b == "WARN"),
    FAIL  = sum(b == "FAIL"),  ERROR = sum(b == "ERROR"),
    INFO  = sum(b == "INFO"),  total = length(b))
}

#' Colour for a badge level (Paul Tol vibrant; WARN = amber/orange)
#' @noRd
.badge_colour <- function(badge) {
  switch(badge,
         PASS  = unname(tol_vibrant["teal"]),
         WARN  = unname(tol_vibrant["orange"]),
         FAIL  = unname(tol_vibrant["red"]),
         ERROR = unname(tol_vibrant["magenta"]),
         unname(tol_vibrant["grey"]))
}


# ---------------------------------------------------------------------------
# S3 print methods
# ---------------------------------------------------------------------------

#' Print method for dynhr_diagnostic_suite
#' @param x   A \code{dynhr_diagnostic_suite} object (list of diagnostics).
#' @param ... Currently unused.
#' @export
#' @method print dynhr_diagnostic_suite
print.dynhr_diagnostic_suite <- function(x, ...) {
  cat("dynhr Diagnostic Suite\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")

  diag_items <- Filter(function(r) inherits(r, "dynhr_diagnostic"), x)
  for (nm in names(diag_items)) {
    cat(sprintf("  [%-5s] %s\n", .badge_str(diag_items[[nm]]), nm))
  }

  cnt <- .badge_counts(diag_items)
  cat(sprintf("\nTotal: %d PASS, %d WARN, %d FAIL, %d ERROR, %d INFO\n",
              cnt[["PASS"]], cnt[["WARN"]], cnt[["FAIL"]],
              cnt[["ERROR"]], cnt[["INFO"]]))
  invisible(x)
}


#' Print method for individual dynhr_diagnostic
#' @param x   A \code{dynhr_diagnostic} object.
#' @param ... Currently unused.
#' @export
#' @method print dynhr_diagnostic
print.dynhr_diagnostic <- function(x, ...) {
  cat(sprintf("[%s] %s\n", .badge_str(x), x$summary))
  if (length(x$plots) > 0)
    cat(sprintf("  Plots: %s\n", paste(names(x$plots), collapse = ", ")))
  invisible(x)
}


# ---------------------------------------------------------------------------
# .DIAG_META — display metadata for the executive summary
# ---------------------------------------------------------------------------
# Keyed by the name used in the results list (same key as the orchestrator).
# Each entry: disp_id, group (A/B/C/D), importance (integer, lower = more
# important), action (one-line recommended fix for FAIL/ERROR).
# Diagnostics absent from this table fall back to group="?", importance=99L.
# ---------------------------------------------------------------------------

.DIAG_META <- list(
  model_summary    = list(disp_id = "D-A1", group = "A", importance =  1L,
    action = "Fix BK conditions or model structure before proceeding."),
  data             = list(disp_id = "D-A2", group = "A", importance =  2L,
    action = "Check observable transformations and sample span."),
  d19              = list(disp_id = "D-A3", group = "A", importance =  3L,
    action = "Review second-order solution quality; consider pruning."),
  expectations     = list(disp_id = "D-A4", group = "A", importance =  4L,
    action = "Fix @dynhr:expectations block constraints."),
  d0               = list(disp_id = "D-A5", group = "A", importance =  5L,
    action = paste0("The static Jacobian is singular. WARN = the solved ",
                    "system has a unit root, which makes the singularity ",
                    "expected (Dynare says the same); FAIL = a genuinely ",
                    "redundant or missing equation. Check the listed ",
                    "collinear equations and unpinned variables.")),
  d40              = list(disp_id = "D-A6", group = "A", importance =  6L,
    action = paste0("A state root is unit/explosive (FAIL) or merely ",
                    "long-lived relative to the sample (WARN). For a unit ",
                    "root use lik_init = 'auto'/'diffuse' and do not rely on ",
                    "unconditional moments; for a near-unit root treat ",
                    "moment-based diagnostics and long-run responses with ",
                    "caution.")),

  d1               = list(disp_id = "D-B1",  group = "B", importance =  1L,
    action = "Check Jacobian rank; fix collinear/unidentified parameters."),
  d23              = list(disp_id = "D-B2",  group = "B", importance =  2L,
    action = "Spectral rank deficient; add observables, restrict priors, or supply a dr-solve function as model_solve_fn."),
  d24              = list(disp_id = "D-B3",  group = "B", importance =  3L,
    action = "Global identification weak; check KL profile per parameter."),
  d20              = list(disp_id = "D-B4",  group = "B", importance =  4L,
    action = "Fisher matrix rank-deficient; CRLBs unreliable."),
  d25              = list(disp_id = "D-B5",  group = "B", importance =  5L,
    action = "Higher-order gain unavailable; recheck D1 baseline rank."),
  d3               = list(disp_id = "D-B6",  group = "B", importance =  6L,
    action = "Parameters with highest mu* are most influential on moments; consider estimating them or restricting priors."),
  d22              = list(disp_id = "D-B7",  group = "B", importance =  7L,
    action = "Check which observables identify which parameters."),
  d37              = list(disp_id = "D-B14", group = "B", importance = 14L,
    action = paste0("The Komunjer-Ng (A,B,C,D) rank condition fails: theta ",
                    "is not identified from the transfer function. Read the ",
                    "null space for the unidentified directions; a ",
                    "non-minimal state space is the usual cause.")),
  d26              = list(disp_id = "D-B8",  group = "B", importance =  8L,
    action = paste0("WARN = estimates move more than 1-for-1 with a ",
                    "calibrated value (|elasticity| > elasticity_tol); ",
                    "FAIL = precision collapses or identification is lost ",
                    "across the perturbation grid. Estimate the listed ",
                    "calibrated parameters, or report results across ",
                    "calibrations.")),
  d4               = list(disp_id = "D-B9",  group = "B", importance =  9L,
    action = "Prior predictive does not cover data; tighten/widen priors."),
  d30              = list(disp_id = "D-B10", group = "B", importance = 10L,
    action = "Precision SVD available; inspect sv_comparison for borderline singular values."),
  d27              = list(disp_id = "D-B11", group = "B", importance = 11L,
    action = "OBC regime changes identification; check piecewise params."),
  d28              = list(disp_id = "D-B12", group = "B", importance = 12L,
    action = "Regime-switching identification; check regime-dependent params."),

  d5               = list(disp_id = "D-C1",  group = "C", importance =  1L,
    action = "Increase chain length, tune step size, or reparameterise."),
  d7               = list(disp_id = "D-C2",  group = "C", importance =  2L,
    action = "Multiple basins found; rerun from more starting points."),
  d6               = list(disp_id = "D-C3",  group = "C", importance =  3L,
    action = "Check posteriors that did not update (overlap near 1)."),
  prior_sensitivity = list(disp_id = "D-C4", group = "C", importance =  4L,
    action = "Informative priors driving results; try flat priors."),
  d21              = list(disp_id = "D-C5",  group = "C", importance =  5L,
    action = "Posterior precision not increasing with T; weak identification."),

  d8               = list(disp_id = "D-D1",  group = "D", importance =  1L,
    action = paste0("Default benchmarks check response SIGNS only; peak ",
                    "timing/magnitude windows are opt-in. Review the ",
                    "calibration against the flagged benchmark, or drop the ",
                    "window you supplied.")),
  d9               = list(disp_id = "D-D2",  group = "D", importance =  2L,
    action = "Adjust shock sizes or model structure to match data SD."),
  d10              = list(disp_id = "D-D3",  group = "D", importance =  3L,
    action = "Review variance decomposition for implausible dominance."),
  d11              = list(disp_id = "D-D4",  group = "D", importance =  4L,
    action = "Inspect historical decomposition plots for plausibility."),
  d12              = list(disp_id = "D-D5",  group = "D", importance =  5L,
    action = "Serial correlation in shocks; check for missing dynamics."),
  d41              = list(disp_id = "D-D23", group = "D", importance = 23L,
    action = paste0("The Kalman-filter innovations are not white: the ",
                    "likelihood's own residuals are autocorrelated, so the ",
                    "state space is missing dynamics. A standardised ",
                    "variance far from 1 instead points at the measurement-",
                    "error or shock scaling.")),
  bayesian_irf     = list(disp_id = "D-D6",  group = "D", importance =  6L,
    action = "Review posterior IRF credible bands for implausible ranges."),
  d13              = list(disp_id = "D-D7",  group = "D", importance =  7L,
    action = "Cross-eq restrictions rejected; review model structure."),
  d15              = list(disp_id = "D-D8",  group = "D", importance =  8L,
    action = "DSGE-VAR tightness suggests misspecification; check lambda."),
  d29              = list(disp_id = "D-D9",  group = "D", importance =  9L,
    action = paste0("The data's variance/autocovariance moments pin the ",
                    "listed parameters weakly (sampling-weighted ",
                    "minimum-distance t-ratios, |theta|/SE). Add informative ",
                    "observables or calibrate them; D20 with ",
                    "weighting = \"sampling\" is the same information in ",
                    "rank form. The Stock-Wright S statistic is reported, ",
                    "not gated.")),
  d16              = list(disp_id = "D-D10", group = "D", importance = 10L,
    action = "Subsample instability; check for structural breaks."),
  d14              = list(disp_id = "D-D11", group = "D", importance = 11L,
    action = "Bayes factor comparison; consider model averaging."),
  model_comparison = list(disp_id = "D-D12", group = "D", importance = 12L,
    action = "Multi-model comparison; review relative log-BF."),
  posterior_predictive = list(disp_id = "D-D13", group = "D", importance = 13L,
    action = "Posterior predictive p-values outside [0.05, 0.95]."),
  d17              = list(disp_id = "D-D14", group = "D", importance = 14L,
    action = "Review sign/magnitude of shock contribution at episode; consider adjusting shock calibration or priors."),
  shock_dominant   = list(disp_id = "D-D15", group = "D", importance = 15L,
    action = "Review dominant shock shares for implausible dominance."),
  d18              = list(disp_id = "D-D16", group = "D", importance = 16L,
    action = "Welfare gap exceeds threshold (|gap/steady_welfare| * 100 > gap_threshold); check Ramsey/optimal policy calibration."),
  d31              = list(disp_id = "D-D17", group = "D", importance = 17L,
    action = paste0("The constrained path violates a bound, or the ",
                    "complementarity solve did not converge. Bounds are ",
                    "compared in LEVELS (steady state + deviation): check ",
                    "the units of `bound`, and that a bound given as a ",
                    "parameter name carries the intended value.")),
  d32              = list(disp_id = "D-D18", group = "D", importance = 18L,
    action = "OBC binding episodes; review binding-period summary."),

  # Deep-parameter block (D33-D36 + passport synthesis)
  d33              = list(disp_id = "D-B13", group = "B", importance = 13L,
    action = "Declare reduced_form= maps or justify borrowed identification."),
  d34              = list(disp_id = "D-D19", group = "D", importance = 19L,
    action = "Private-block parameters unstable; check Lucas-critique exposure."),
  d35              = list(disp_id = "D-D20", group = "D", importance = 20L,
    action = "Soft parameters; compare sandwich vs Hessian standard errors."),
  d36              = list(disp_id = "D-D21", group = "D", importance = 21L,
    action = "Calibrated parameters in tension with data; consider estimating."),
  deep_passport    = list(disp_id = "D-D22", group = "D", importance = 22L,
    action = "Review deep-parameter passport grades.")
)

## An ERROR means the diagnostic never ran, so its .DIAG_META action -- which
## explains how to fix a FAILING result -- is not merely unhelpful but wrong:
## the review found "[ERROR] d12 / Error in solve(): system is exactly singular
## / ACTION: Serial correlation in shocks; check for missing dynamics."
.DIAG_ERROR_ACTION <- paste0(
  "This diagnostic did not complete -- fix the error above before reading ",
  "its action.")

#' Recommended-action line for a diagnostic, given its badge
#'
#' @param nm    Diagnostic name as used in the results list.
#' @param badge Badge from \code{\link{.badge_str}}.
#' @return Character scalar ("" when there is no action to give).
#' @noRd
.diag_action_for <- function(nm, badge) {
  if (identical(badge, "ERROR")) return(.DIAG_ERROR_ACTION)
  .diag_meta_for(nm)$action %||% ""
}

#' Look up display metadata for a diagnostic by its result-list name
#' @noRd
.diag_meta_for <- function(nm) {
  m <- .DIAG_META[[nm]]
  if (is.null(m)) list(disp_id = nm, group = "?", importance = 99L, action = "")
  else m
}


# ---------------------------------------------------------------------------
# One ordering rule for the whole report
# ---------------------------------------------------------------------------
# The 2026-09-17 report review (D1) found the front page and the detail
# sections disagreed about what mattered: the executive summary sorted by
# .DIAG_META group then importance, while both Quarto templates rendered in
# orchestrator list order, and the HTML tab strip put ERROR *after* PASS.
#
# Severity leads, because "did not run" is strictly worse than "ran and
# failed", and PASS is the section nobody needs to read. Within a severity
# level the executive summary's existing rule (group, then importance) is kept,
# with the name as a deterministic tiebreak so two diagnostics that share a
# group and importance never swap places between renders.

.BADGE_RANK <- c(ERROR = 1L, FAIL = 2L, WARN = 3L, INFO = 4L, PASS = 5L)
.GROUP_RANK <- c(A = 1L, B = 2L, C = 3L, D = 4L, "?" = 5L)

#' Rank of a group letter, 5L for anything unregistered
#' @noRd
.group_rank <- function(g) {
  r <- unname(.GROUP_RANK[match(as.character(g), names(.GROUP_RANK))])
  if (length(r) != 1L || is.na(r)) 5L else r
}

#' Order a set of diagnostics: severity, then group, then importance, then name
#'
#' @param diag_items Named list of \code{dynhr_diagnostic} objects.
#' @return Character vector of names in report order.
#' @noRd
.diag_order <- function(diag_items) {
  nms <- names(diag_items)
  if (is.null(nms) || length(nms) == 0L) return(character(0))
  badges <- vapply(diag_items, .badge_str, character(1))
  brank  <- unname(.BADGE_RANK[badges])
  brank[is.na(brank)] <- 6L
  metas  <- lapply(nms, .diag_meta_for)
  nms[order(brank,
            vapply(metas, function(m) .group_rank(m$group),      integer(1)),
            vapply(metas, function(m) as.integer(m$importance), integer(1)),
            nms)]
}
