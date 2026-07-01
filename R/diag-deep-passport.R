## R/diag-deep-passport.R
## --------------------------------------------------------------------------
## Deep-Parameter Passport — the synthesising deliverable of the DP module.
##
## Collapses the (now fourteen) identification / sensitivity diagnostics plus
## the genuinely-new deep-parameter diagnostics into ONE row per parameter: a
## verdict across the four questions a parameter must answer to be called
## "deep" — Identified? Primitive? Invariant? Robust? — and an A-F deepness
## grade. Auxiliary shock parameters are shown but graded "-".
##
## It is intentionally *defensive*: every axis degrades to NA when its source
## diagnostic is absent (D33/D34/D35 not built yet) or is known-unreliable
## (D24's vacuous KL — see DEEP_PARAMETER_DIAGNOSTICS_PLAN.md §2). It therefore
## runs today on a suite containing only D1/D6 and fills in automatically as the
## remaining axes land.
## --------------------------------------------------------------------------


# ---- per-axis extractors (best-effort, NA when unavailable) ---------------

# Identified (Q1): a parameter is flagged weak if a *reliable* local-ID
# diagnostic names it. We deliberately exclude D24 (global KL) until its
# theta-re-solve bug is fixed, since in its current form it flags everything.
.axis_identified <- function(suite, params) {
  out <- stats::setNames(rep(NA, length(params)), params)
  if (is.null(suite)) return(out)
  weak <- character(0); assessed <- FALSE
  for (nm in c("d1", "d20")) {
    d <- suite[[nm]]
    wp <- d$result$weak_params
    if (!is.null(wp)) { weak <- c(weak, as.character(wp)); assessed <- TRUE }
  }
  if (assessed) {
    out[] <- TRUE
    out[params %in% unique(weak)] <- FALSE
  }
  out
}

# Informed (Q1): posterior moves away from the prior. Uses D6's overlap
# coefficient (overlap > 0.80 == uninformative), else D21's KPS slope.
.axis_informed <- function(suite, params) {
  out <- stats::setNames(rep(NA, length(params)), params)
  if (is.null(suite)) return(out)
  ov <- suite$d6$result$overlap_scores
  if (!is.null(ov)) {
    hit <- intersect(params, names(ov))
    out[hit] <- ov[hit] < 0.80
  }
  # Fallback / supplement: KPS contraction slope (slope < threshold -> not informed)
  sl <- suite$d21$result$slopes %||% suite$d21$result$slope
  if (!is.null(sl) && !is.null(names(sl))) {
    hit <- intersect(params[is.na(out)], names(sl))
    if (length(hit)) out[hit] <- sl[hit] > 0.2
  }
  out
}

# Generic axis from a not-yet-built diagnostic that, by contract, exposes
# `$result$passport_axis`: a named logical keyed by parameter (TRUE = good).
.axis_from <- function(suite, diag_name, params) {
  out <- stats::setNames(rep(NA, length(params)), params)
  ax <- suite[[diag_name]]$result$passport_axis
  if (!is.null(ax) && !is.null(names(ax))) {
    hit <- intersect(params, names(ax))
    out[hit] <- as.logical(ax[hit])
  }
  out
}


# ---------------------------------------------------------------------------
#' Deep-Parameter Passport
#'
#' Synthesises a per-parameter "deepness" scorecard from the deep-parameter
#' taxonomy and the available identification / deep diagnostics.  Each deep
#' parameter is scored on five axes — identified, informed, structural,
#' invariant, robust — and given an A-F grade; auxiliary shock parameters are
#' listed but graded \code{"-"}.
#'
#' @param deep_spec  A \code{\link{build_deep_spec}} object, or a parsed model
#'   (from which a taxonomy is built).
#' @param suite      Optional \code{dynhr_diagnostic_suite} (the output of
#'   \code{run_all_diagnostics()}); axes are read from it where present.
#' @param draws      Optional posterior draw matrix (\eqn{n \times p}) with
#'   named columns; adds posterior median/sd context columns.
#' @param axes       Optional named list of per-parameter logical vectors that
#'   override the suite-derived axes, e.g.
#'   \code{list(structural = c(theta_H = FALSE, alpha = TRUE))}.
#' @param meta       Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic} whose \code{$result$passport} is the tidy
#'   scorecard \code{data.frame}.
#' @export
deep_parameter_passport <- function(deep_spec,
                                     suite = NULL,
                                     draws = NULL,
                                     axes  = NULL,
                                     meta  = NULL) {
  tryCatch({
    # Accept a model and build the taxonomy on the fly.
    if (!inherits(deep_spec, "dynhr_deep_spec")) {
      deep_spec <- build_deep_spec(model = deep_spec)
    }
    if (nrow(deep_spec) == 0) {
      return(.make_result(
        pass = NA,
        summary = "Deep-Parameter Passport: no parameters to classify."))
    }

    params <- deep_spec$param

    # --- assemble the five axes (suite-derived, then explicit overrides) ---
    A <- list(
      identified = .axis_identified(suite, params),
      informed   = .axis_informed(suite, params),
      structural = .axis_from(suite, "d33", params),
      invariant  = .axis_from(suite, "d34", params),
      robust     = .axis_from(suite, "d35", params),
      calibrated = .axis_from(suite, "d36", params)
    )
    if (!is.null(axes)) {
      for (nm in intersect(names(axes), names(A))) {
        ax <- axes[[nm]]
        hit <- intersect(params, names(ax))
        A[[nm]][hit] <- as.logical(ax[hit])
      }
    }

    # --- posterior context columns (optional) ---
    pmed <- psd <- stats::setNames(rep(NA_real_, length(params)), params)
    if (!is.null(draws)) {
      draws <- as.matrix(draws)
      hit <- intersect(params, colnames(draws))
      if (length(hit)) {
        pmed[hit] <- apply(draws[, hit, drop = FALSE], 2, stats::median)
        psd[hit]  <- apply(draws[, hit, drop = FALSE], 2, stats::sd)
      }
    }

    # --- grade each row ---
    grade <- vapply(seq_along(params), function(i) {
      .deepness_grade(
        identified = A$identified[i], informed = A$informed[i],
        structural = A$structural[i], invariant = A$invariant[i],
        robust     = A$robust[i],     calibrated = A$calibrated[i],
        is_deep    = deep_spec$is_deep[i])
    }, character(1))

    passport <- data.frame(
      param      = params,
      class      = deep_spec$class,
      partition  = deep_spec$partition,
      is_deep    = deep_spec$is_deep,
      identified = A$identified,
      informed   = A$informed,
      structural = A$structural,
      invariant  = A$invariant,
      robust     = A$robust,
      calibrated = A$calibrated,
      reduced_form = deep_spec$reduced_form,
      post_median  = unname(pmed),
      post_sd      = unname(psd),
      source     = deep_spec$source,
      grade      = grade,
      stringsAsFactors = FALSE
    )
    # order: deep first, then by grade (worst first, unassessed "?" after the
    # graded ones), auxiliary last
    ord <- order(!passport$is_deep,
                 match(passport$grade, c("F","E","D","C","B","A","?","-")),
                 passport$param)
    passport <- passport[ord, , drop = FALSE]
    rownames(passport) <- NULL

    # --- coverage of the assessment ---
    deep_idx   <- which(passport$is_deep)
    axis_names <- c("identified", "informed", "structural", "invariant",
                    "robust", "calibrated")
    n_assessed <- vapply(axis_names, function(a) sum(!is.na(passport[deep_idx, a])),
                         integer(1))
    not_yet <- axis_names[n_assessed == 0]

    # --- plot ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plots$passport <- .plot_passport(passport, meta = meta)
    }

    # --- summaries ---
    n_deep <- length(deep_idx)
    grade_tab <- table(factor(passport$grade[deep_idx],
                              levels = c("A","B","C","D","E","F","?")))
    n_unassessed <- sum(passport$grade[deep_idx] == "?")
    flagged <- function(axis) {
      f <- passport$param[passport$is_deep & identical_false(passport[[axis]])]
      f[!is.na(f)]
    }
    summary_txt <- sprintf(
      "Deep-Parameter Passport: %d deep params graded (%s). %s",
      n_deep,
      paste(sprintf("%s:%d", names(grade_tab), as.integer(grade_tab)),
            collapse = " "),
      if (length(not_yet))
        sprintf("Axes not yet assessed: %s.", paste(not_yet, collapse = ", "))
      else "All axes assessed.")

    n_auto <- sum(deep_spec$source == "auto" & deep_spec$is_deep)

    llm <- paste(c(
      "DEEP-PARAMETER PASSPORT | INFO",
      sprintf("  deep=%d auxiliary=%d", n_deep, nrow(passport) - n_deep),
      sprintf("  grades: %s",
              paste(sprintf("%s=%d", names(grade_tab), as.integer(grade_tab)),
                    collapse = " ")),
      if (n_unassessed > 0)
        sprintf("  unassessed(?): %d deep param(s) have no axis assessed (e.g. calibrated, data-silent)",
                n_unassessed),
      {
        worst <- passport$param[passport$is_deep &
                                passport$grade %in% c("D","E","F")]
        if (length(worst))
          sprintf("  flagged(D-F): %s", paste(utils::head(worst, 8), collapse = ", "))
        else NULL
      },
      if (length(flagged("identified")))
        sprintf("  weak_id: %s", paste(flagged("identified"), collapse = ", ")),
      if (length(flagged("structural")))
        sprintf("  borrowed: %s", paste(flagged("structural"), collapse = ", ")),
      if (length(flagged("invariant")))
        sprintf("  non_invariant: %s", paste(flagged("invariant"), collapse = ", ")),
      if (length(flagged("robust")))
        sprintf("  soft: %s", paste(flagged("robust"), collapse = ", ")),
      if (length(flagged("calibrated")))
        sprintf("  calibration_tension: %s", paste(flagged("calibrated"), collapse = ", ")),
      if (length(not_yet))
        sprintf("  not_yet_assessed: %s", paste(not_yet, collapse = ", ")),
      if (n_auto > 0)
        sprintf("  warning: %d deep param(s) auto-classified; declare @dynhr:deep to confirm",
                n_auto),
      "  action: grades synthesise D1/D6 (+D33/D34/D35 when present); a D/E/F grade means at least one deepness axis failed."
    ), collapse = "\n")

    .make_result(
      result  = list(passport = passport, axes = A, not_assessed = not_yet),
      pass    = NA,
      plots   = plots,
      summary = summary_txt,
      llm_summary = llm
    )
  }, error = function(e) {
    .make_result(pass = NA, errored = TRUE,
                 summary = paste("Deep-Parameter Passport: ERROR --",
                                 conditionMessage(e)))
  })
}

# small helper: vectorised "is identically FALSE" that tolerates NA
identical_false <- function(v) !is.na(v) & (v == FALSE)


# ---------------------------------------------------------------------------
# Render the passport as a coloured param x axis tile grid with a grade column.
# ---------------------------------------------------------------------------
.plot_passport <- function(passport, meta = NULL) {
  axes <- c("identified", "informed", "structural", "invariant", "robust",
            "calibrated")
  axis_lab <- c(identified = "ID", informed = "Info", structural = "Struct",
                invariant = "Invar", robust = "Robust", calibrated = "Calib")

  # Drop rows where ALL axes are NA (fully "not assessed" deep params) to keep
  # the heatmap compact.  Record the count for the subtitle.
  deep_rows <- passport[passport$is_deep, , drop = FALSE]
  assessed_mask <- if (nrow(deep_rows) > 0)
    vapply(seq_len(nrow(deep_rows)), function(i)
      any(!is.na(deep_rows[i, axes])), logical(1))
  else
    logical(0)
  passport_plot <- rbind(
    deep_rows[assessed_mask, , drop = FALSE],
    passport[!passport$is_deep, , drop = FALSE]
  )
  n_unassessed_dropped <- sum(!assessed_mask)

  long <- do.call(rbind, lapply(axes, function(a) {
    v <- passport_plot[[a]]
    status <- ifelse(!passport_plot$is_deep, "n/a (aux)",
              ifelse(is.na(v), "not assessed",
              ifelse(v, "pass", "fail")))
    data.frame(param = passport_plot$param, axis = unname(axis_lab[a]),
               status = status, stringsAsFactors = FALSE)
  }))

  lvl <- rev(passport_plot$param)  # top row = best (deep, A) after the ordering
  long$param <- factor(long$param, levels = lvl)
  long$axis  <- factor(long$axis, levels = unname(axis_lab))
  long$status <- factor(long$status,
                        levels = c("pass", "fail", "not assessed", "n/a (aux)"))

  fill_vals <- c(
    "pass"         = dynhr_colours$teal,
    "fail"         = dynhr_colours$red,
    "not assessed" = dynhr_colours$grey,
    "n/a (aux)"    = "#E8E8E8")

  grade_df <- data.frame(param = factor(passport_plot$param, levels = lvl),
                         grade = passport_plot$grade)

  # Row height: at least 0.32 inches per row so labels are legible.
  n_rows <- length(lvl)
  fig_h  <- max(4.0, n_rows * 0.32 + 2.0)

  subtitle_txt <- paste0(
    "Per-parameter deepness: ID / Informed / Structural / Invariant / Robust  ->  grade",
    if (n_unassessed_dropped > 0)
      sprintf("  |  %d fully-unassessed deep param(s) not shown", n_unassessed_dropped)
    else "")

  p <- ggplot2::ggplot(long, ggplot2::aes(x = axis, y = param, fill = status)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(data = grade_df,
                       ggplot2::aes(x = "Calib", y = param, label = grade),
                       inherit.aes = FALSE, hjust = -1.4, fontface = "bold",
                       size = 3.2, colour = dynhr_colours$dark_blue) +
    ggplot2::scale_fill_manual(values = fill_vals, drop = FALSE, name = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title    = "Deep-Parameter Passport",
      subtitle = subtitle_txt,
      x = NULL, y = NULL) +
    theme_dynhr() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 28, 5.5, 5.5),
      # Ensure row labels stay readable when there are many rows
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.85))
    )
  attr(p, "dynhr_fig_height") <- fig_h

  .apply_meta(p, meta)
}
