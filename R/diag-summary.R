## R/diag-summary.R
## --------------------------------------------------------------------------
## Executive summary for the dynhr diagnostic suite.
##
## Public API:
##   summary.dynhr_diagnostic_suite()  -- S3 method; prints + returns invisibly
##   write_executive_summary()         -- writes markdown file
##
## Internal:
##   format_executive_summary()        -- returns character vector of lines
## --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# Collapse a (possibly multi-line) summary to one line and truncate with an
# ellipsis at a word boundary where possible.
.diag_truncate <- function(x, width = 120L) {
  x <- gsub("[[:space:]]+", " ", trimws(x))
  if (nchar(x) <= width) return(x)
  out <- substr(x, 1L, width - 3L)
  sp  <- regexpr(" [^ ]*$", out)
  if (sp > width %/% 2) out <- substr(out, 1L, sp - 1L)
  paste0(out, "...")
}


# ---------------------------------------------------------------------------
# format_executive_summary()
# ---------------------------------------------------------------------------

#' Format the executive summary as a character vector (one element per line)
#'
#' Used by \code{summary.dynhr_diagnostic_suite} (console) and
#' \code{write_executive_summary} (markdown file) and by the Quarto templates.
#'
#' @param suite      A \code{dynhr_diagnostic_suite} (named list).
#' @param model_name Optional character label for the model.
#' @param n_critical Integer: maximum Critical findings shown (default 6).
#' @param n_watch    Integer: maximum Watch list items shown (default 8).
#' @param unicode    Logical: not currently used; reserved for box-drawing
#'   variants. Pass \code{FALSE} for plain markdown output (default \code{TRUE}
#'   for both to keep output identical and portable).
#' @return Character vector, one element per output line.
#' @noRd
format_executive_summary <- function(suite,
                                      model_name = NULL,
                                      n_critical = 6L,
                                      n_watch    = 8L,
                                      unicode    = TRUE) {

  diag_items <- Filter(function(r) inherits(r, "dynhr_diagnostic"), suite)

  # ---- Badge and metadata per item ----------------------------------------
  items <- lapply(names(diag_items), function(nm) {
    r <- diag_items[[nm]]
    m <- .diag_meta_for(nm)
    badge <- .badge_str(r)
    list(nm           = nm,
         r            = r,
         badge        = badge,
         disp_id      = m$disp_id,
         group        = m$group,
         importance   = m$importance,
         action       = m$action,
         short_summary = .diag_truncate(r$summary %||% "", 120L))
  })

  # ---- Global counts -------------------------------------------------------
  n_pass <- sum(vapply(items, function(x) x$badge == "PASS",  logical(1)))
  n_fail <- sum(vapply(items, function(x) x$badge == "FAIL",  logical(1)))
  n_err  <- sum(vapply(items, function(x) x$badge == "ERROR", logical(1)))
  n_info <- sum(vapply(items, function(x) x$badge == "INFO",  logical(1)))

  # ---- Group scorecard: PASS/(PASS+FAIL) -----------------------------------
  group_labels <- c(A = "Pre-solve / structure",
                    B = "Identification & sensitivity",
                    C = "Estimation / MCMC",
                    D = "Post-estimation")
  sc <- lapply(names(group_labels), function(g) {
    gi <- Filter(function(x) x$group == g, items)
    gp <- sum(vapply(gi, function(x) x$badge == "PASS",              logical(1)))
    gf <- sum(vapply(gi, function(x) x$badge %in% c("FAIL", "ERROR"), logical(1)))
    list(g        = g,
         label    = group_labels[[g]],
         pass     = gp,
         testable = gp + gf,
         n_fail   = gf)
  })
  names(sc) <- names(group_labels)

  # ---- Verdict: FAIL if A/B fail, WARN if C/D fail, else PASS --------------
  ab_fails <- sum(vapply(
    Filter(function(x) x$group %in% c("A", "B"), items),
    function(x) x$badge %in% c("FAIL", "ERROR"), logical(1)))
  verdict <- if (ab_fails > 0L || n_err > 0L) "FAIL"
             else if (n_fail > 0L)             "WARN"
             else                              "PASS"

  # ---- Sort fails by group order then importance ---------------------------
  group_rank <- c(A = 1L, B = 2L, C = 3L, D = 4L, "?" = 5L)
  fails <- Filter(function(x) x$badge %in% c("FAIL", "ERROR"), items)
  if (length(fails) > 0L) {
    ord <- order(
      vapply(fails, function(x) group_rank[[x$group]] %||% 5L, integer(1)),
      vapply(fails, function(x) x$importance,                   integer(1))
    )
    fails <- fails[ord]
  }

  critical        <- head(fails, n_critical)
  remaining_fails <- if (length(fails) > n_critical) tail(fails, length(fails) - n_critical) else list()

  # INFO items at high importance (rank <= 4 within group)
  infos <- Filter(function(x) x$badge == "INFO" && x$importance <= 4L, items)
  # Sort infos by group then importance too
  if (length(infos) > 0L) {
    io <- order(
      vapply(infos, function(x) group_rank[[x$group]] %||% 5L, integer(1)),
      vapply(infos, function(x) x$importance,                   integer(1))
    )
    infos <- infos[io]
  }

  watch  <- head(c(remaining_fails, infos), n_watch)
  passes <- Filter(function(x) x$badge == "PASS", items)
  if (length(passes) > 0L) {
    po <- order(
      vapply(passes, function(x) group_rank[[x$group]] %||% 5L, integer(1)),
      vapply(passes, function(x) x$importance,                   integer(1))
    )
    passes <- passes[po]
  }

  # ---- Assemble lines ------------------------------------------------------
  lines <- character(0L)
  add   <- function(...) lines <<- c(lines, paste0(...))
  sep   <- function(ch = "-", n = 80L) add(strrep(ch, n))

  sep("=")
  add("           dynhr DIAGNOSTIC EXECUTIVE SUMMARY")
  sep("=")
  if (!is.null(model_name))
    add(sprintf("  Model:  %s", model_name))
  add(sprintf("  Run:    %s", format(Sys.Date(), "%Y-%m-%d")))
  sep("=")
  add("")

  verdict_line <- switch(verdict,
    FAIL = "!! FAIL !!",
    WARN = "!! WARN !!",
    PASS = "   PASS   ")
  add(sprintf("OVERALL VERDICT: %s", verdict_line))
  add(sprintf("  %d PASS  |  %d FAIL  |  %d ERROR  |  %d INFO  (%d diagnostics total)",
              n_pass, n_fail, n_err, n_info, length(items)))
  add("")

  add("GROUP SCORECARD   (PASS / testable; INFO excluded from denominator)")
  add("")
  for (s in sc) {
    status_label <- if (s$n_fail > 0L) "[FAIL]"
                    else if (s$pass < s$testable) "[WARNING]"
                    else "[PASS]"
    add(sprintf("  %s  %-36s  %d / %d  %s",
                s$g, s$label, s$pass, s$testable, status_label))
  }

  # One-line verdict rationale
  add("")
  if (verdict == "FAIL") {
    worst_g <- if (sc[["A"]]$n_fail > 0L) "A" else "B"
    add(sprintf("  Verdict rule: FAIL because Group %s has %d failure(s) (pipeline blocker).",
                worst_g, sc[[worst_g]]$n_fail))
  } else if (verdict == "WARN") {
    add("  Verdict rule: WARN because C/D diagnostics have failures (no group-A/B blockers).")
  } else {
    add("  Verdict rule: PASS -- no failures in any group.")
  }
  add("")

  sep()
  add("CRITICAL FINDINGS   (top by group then importance rank)")
  sep()
  add("")
  if (length(critical) == 0L) {
    add("  (none)")
    add("")
  } else {
    for (x in critical) {
      add(sprintf("[%s] %s  %s", x$badge, x$disp_id, x$nm))
      s <- x$short_summary
      if (nchar(s) > 0L) {
        # Indent continuation lines of multi-line summaries
        slines <- strsplit(s, "\n", fixed = TRUE)[[1]]
        for (sl in slines) add(sprintf("  %s", sl))
      }
      if (nchar(x$action) > 0L)
        add(sprintf("  ACTION: %s", x$action))
      add("")
    }
  }

  sep()
  add("WATCH LIST   (remaining FAILs + high-importance INFOs)")
  sep()
  add("")
  if (length(watch) == 0L) {
    add("  (none)")
  } else {
    for (x in watch) {
      s80 <- .diag_truncate(x$short_summary, 80L)
      add(sprintf("[%s] %-7s  %s -- %s",
                  x$badge, x$disp_id, x$nm, s80))
    }
  }
  add("")

  sep()
  if (length(passes) > 0L) {
    pass_ids <- paste(vapply(passes, function(x) x$disp_id, character(1)),
                      collapse = "  ")
    # Wrap at ~75 chars
    max_w <- 72L
    words  <- strsplit(pass_ids, "  ")[[1]]
    cur    <- "PASSING  "
    for (w in words) {
      if (nchar(cur) + nchar(w) + 2L > max_w && cur != "PASSING  ") {
        add(cur)
        cur <- paste0("         ", w)
      } else {
        cur <- if (cur == "PASSING  ") paste0(cur, w) else paste0(cur, "  ", w)
      }
    }
    add(cur)
  } else {
    add("PASSING   (none)")
  }
  sep("=")

  lines
}


# ---------------------------------------------------------------------------
# S3 summary method
# ---------------------------------------------------------------------------

#' Summary method for a dynhr diagnostic suite
#'
#' Prints a one-to-two-page executive summary to the console and returns the
#' character vector of lines invisibly (suitable for \code{capture.output}).
#'
#' @param object A \code{dynhr_diagnostic_suite}.
#' @param ...    Currently unused.
#' @param model_name Optional character label override for the model name.
#' @param n_critical Integer: maximum Critical findings shown (default 6).
#' @param n_watch    Integer: maximum Watch list items shown (default 8).
#' @return The character vector of lines, invisibly.
#' @export
summary.dynhr_diagnostic_suite <- function(object, ...,
                                            model_name = NULL,
                                            n_critical = 6L,
                                            n_watch    = 8L) {
  lines <- format_executive_summary(object,
                                     model_name = model_name,
                                     n_critical = n_critical,
                                     n_watch    = n_watch,
                                     unicode    = TRUE)
  cat(paste(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}


# ---------------------------------------------------------------------------
# write_executive_summary()
# ---------------------------------------------------------------------------

#' Write the executive summary as a Markdown file
#'
#' Produces a plain Markdown document (no HTML, no ggplot) that can be
#' opened in any text editor or committed alongside the diagnostic suite RDS.
#'
#' @param suite      A \code{dynhr_diagnostic_suite}.
#' @param file       Output path (default \code{"executive_summary.md"}).
#' @param model_name Optional character label for the model.
#' @param ...        Additional arguments forwarded to
#'   \code{format_executive_summary} (e.g. \code{n_critical}, \code{n_watch}).
#' @return \code{file}, invisibly.
#' @export
write_executive_summary <- function(suite,
                                     file       = "executive_summary.md",
                                     model_name = NULL,
                                     ...) {
  lines <- format_executive_summary(suite,
                                     model_name = model_name,
                                     unicode    = FALSE,
                                     ...)
  writeLines(lines, file)
  message("[dynhr] Executive summary written to ", file)
  invisible(file)
}
