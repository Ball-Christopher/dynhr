## R/diag-pre-prior-sensitivity.R
## --------------------------------------------------------------------------
## Prior Sensitivity Diagnostic — compares mode estimates under informative
## priors vs flat (max-entropy / uniform) priors over the same support.
##
## Parameters whose posterior mode shifts substantially under flat priors
## are "prior-driven" — the data alone does not pin them down.
##
## Usage:
##   Added to run_all_diagnostics() when the model has an estimated_params
##   block and the data is available.
## --------------------------------------------------------------------------


#' Prior sensitivity: informative vs flat (uniform) priors
#'
#' Re-runs mode-finding with uniform priors over the same support as the
#' original informative priors, then compares the two mode vectors.
#' Parameters whose relative shift exceeds `threshold` are flagged as
#' "prior-driven."
#'
#' @param solved       A \code{dynhr_solved} object.
#' @param data         Observation matrix (T x n_obs).
#' @param obs_vars     Observable variable names.
#' @param mode_inf     The \code{dynhr_mode_result} from informative-prior run.
#' @param me_variance  Measurement error variance (same as used for informative).
#' @param threshold    Relative shift threshold for flagging (default 0.50,
#'   i.e. a 50% change in parameter value relative to its informative-mode
#'   magnitude is flagged).
#' @param n_iter       Optimizer iterations for flat-prior mode finding
#'   (default 500).  Fewer iterations are acceptable since we only need
#'   a rough comparison.
#' @param verbose      Print progress messages.
#' @param meta         Optional [diag_meta()] list applied to the plots
#'   (model name, sample size); `NULL` derives a default from the model.
#' @return A \code{dynhr_diagnostic} list.
#' @export
diag_prior_sensitivity <- function(solved,
                                    data,
                                    obs_vars,
                                    mode_inf,
                                    me_variance    = 0,
                                    threshold      = 0.50,
                                    n_iter         = 500L,
                                    verbose        = TRUE,
                                    meta           = NULL) {

  .vcat <- function(...) if (verbose) cat(...)

    .vcat("[prior_sensitivity] Building flat-prior model...\n")

    mod_file <- solved$model$source_file
    if (is.null(mod_file) || !file.exists(mod_file)) {
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: model source file not found.",
        errored = TRUE
      ))
    }

    # ---- 1. Create flat-prior .mod file ----
    mod_text <- readLines(mod_file)
    est_start <- grep("estimated_params;", mod_text, ignore.case = TRUE)[1]
    if (is.na(est_start)) {
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: no estimated_params block found.",
        errored = TRUE
      ))
    }

    # Find the end of estimated_params block
    est_end_candidates <- grep("^end;", mod_text)
    est_end <- est_end_candidates[est_end_candidates > est_start][1]
    if (is.na(est_end)) {
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: cannot find end of estimated_params block.",
        errored = TRUE
      ))
    }

    # Parse the existing estimated_params to extract distribution types and bounds
    est_block <- mod_text[(est_start + 1):(est_end - 1)]
    est_block <- est_block[!grepl("^//|^#|^$", est_block)]

    # Build flat-prior replacement block
    flat_lines <- c("estimated_params;", "")
    for (line in est_block) {
      parts <- strsplit(trimws(line), ",")[[1]]
      if (length(parts) < 2) next

      # stderr lines: stderr NAME, distribution, p1, p2, lower, upper
      # param lines:  NAME,    distribution, p1, p2, lower, upper
      has_stderr <- grepl("stderr", parts[1])
      
      if (has_stderr) {
        # stderr eps_a, inv_gamma_pdf, 0.010, 0.006, 0.0005, 0.50
        name <- trimws(gsub("stderr", "", parts[1]))
        if (length(parts) >= 6) {
          lower <- trimws(parts[length(parts) - 1])
          upper <- trimws(gsub(";", "", parts[length(parts)]))
          flat_lines <- c(flat_lines,
            sprintf("  stderr %s, uniform_pdf, , , %s, %s;", name, lower, upper))
        } else {
          flat_lines <- c(flat_lines,
            sprintf("  stderr %s, uniform_pdf, , ;", name))
        }
      } else {
        # rho_a, beta_pdf, 0.90, 0.04, 0.50, 0.999
        name <- trimws(parts[1])
        if (length(parts) >= 6) {
          lower <- trimws(parts[length(parts) - 1])
          upper <- trimws(gsub(";", "", parts[length(parts)]))
          flat_lines <- c(flat_lines,
            sprintf("  %s, uniform_pdf, , , %s, %s;", name, lower, upper))
        } else {
          flat_lines <- c(flat_lines,
            sprintf("  %s, uniform_pdf, , ;", name))
        }
      }
    }
    flat_lines <- c(flat_lines, "", "end;")

    # Write flat-prior model
    flat_mod_file <- tempfile(fileext = "_flat.mod")
    new_text <- c(
      mod_text[1:(est_start - 1)],
      flat_lines,
      mod_text[(est_end + 1):length(mod_text)]
    )
    writeLines(new_text, flat_mod_file)
    .vcat(sprintf("[prior_sensitivity] Flat-prior model: %s\n", basename(flat_mod_file)))

    # ---- 2. Solve and run mode-finding with flat priors ----
    .vcat("[prior_sensitivity] Solving flat-prior model...\n")
    solved_flat <- solve_model(flat_mod_file, verbose = FALSE)

    if (is.null(solved_flat) || !inherits(solved_flat, "dynhr_solved") ||
        is.null(solved_flat$model) || is.null(solved_flat$compiled)) {
      unlink(flat_mod_file)
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: flat-prior model solve returned invalid result.",
        errored = TRUE
      ))
    }
    .vcat("[prior_sensitivity] Flat model BK satisfied.\n")

    .vcat("[prior_sensitivity] Running mode-finding with flat priors...\n")
    mode_flat <- run_mode_finding(solved_flat, data, obs_vars = obs_vars,
                                   n_iter = n_iter, method = "nmkb",
                                   me_variance = me_variance,
                                   verbose = FALSE)
    if (is.null(mode_flat) || is.null(mode_flat$mode) ||
        !is.finite(mode_flat$mode$logpost %||% NA_real_)) {
      unlink(flat_mod_file)
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: flat-prior mode finding returned invalid result.",
        errored = TRUE
      ))
    }
    .vcat(sprintf("[prior_sensitivity] Flat mode logpost = %.2f\n", 
                  mode_flat$mode$logpost))

    # ---- 3. Compare modes ----
    theta_inf <- mode_inf$mode$theta_mode
    theta_flat <- mode_flat$mode$theta_mode
    all_params <- union(names(theta_inf), names(theta_flat))

    comparisons <- data.frame(
      param = character(),
      informative = numeric(),
      flat = numeric(),
      abs_shift = numeric(),
      rel_shift = numeric(),
      flagged = logical(),
      stringsAsFactors = FALSE
    )

    for (p in all_params) {
      v_inf <- if (p %in% names(theta_inf)) theta_inf[p] else NA_real_
      v_flat <- if (p %in% names(theta_flat)) theta_flat[p] else NA_real_
      abs_shift <- abs(v_inf - v_flat)
      ref_val <- max(abs(v_inf), 1e-10)
      rel_shift <- abs_shift / ref_val
      flagged <- rel_shift > threshold
      comparisons <- rbind(comparisons, data.frame(
        param = p,
        informative = round(v_inf, 6),
        flat = round(v_flat, 6),
        abs_shift = round(abs_shift, 6),
        rel_shift = round(rel_shift, 3),
        flagged = flagged,
        stringsAsFactors = FALSE
      ))
    }

    n_flagged <- sum(comparisons$flagged)
    all_pass <- n_flagged == 0

    # ---- 4. Build result ----
    result <- list(
      comparisons = comparisons,
      theta_inf = theta_inf,
      theta_flat = theta_flat,
      n_flagged = n_flagged,
      threshold = threshold,
      logpost_inf = mode_inf$mode$logpost,
      logpost_flat = mode_flat$mode$logpost
    )

    # Summary
    if (all_pass) {
      summary_str <- sprintf(
        "Prior sensitivity: %d/%d parameters within threshold (%.0f%%). All parameters are data-driven.",
        nrow(comparisons) - n_flagged, nrow(comparisons), threshold * 100)
    } else {
      flagged_params <- comparisons$param[comparisons$flagged]
      summary_str <- sprintf(
        "Prior sensitivity: %d/%d parameters flagged (threshold=%.0f%%). Prior-driven: %s.",
        n_flagged, nrow(comparisons), threshold * 100,
        paste(flagged_params, collapse = ", "))
    }

    # ---- 5. Plot ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      comp_plot <- comparisons
      comp_plot$param <- factor(comp_plot$param, levels = rev(comp_plot$param))
      comp_plot$status <- ifelse(comp_plot$flagged, "Prior-driven", "Data-driven")

      p <- ggplot2::ggplot(comp_plot) +
        ggplot2::geom_segment(
          ggplot2::aes(x = informative, y = param, 
                       xend = flat, yend = param,
                       colour = status),
          linewidth = 1.2, alpha = 0.8) +
        ggplot2::geom_point(
          ggplot2::aes(x = informative, y = param, colour = "Informative"),
          size = 2.5) +
        ggplot2::geom_point(
          ggplot2::aes(x = flat, y = param, colour = "Flat (uniform)"),
          size = 2.5, shape = 17) +
        ggplot2::scale_colour_manual(
          values = c("Informative" = "#0077BB",
                     "Flat (uniform)" = "#EE7733",
                     "Data-driven" = "#009988",
                     "Prior-driven" = "#CC3311"),
          name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "Prior Sensitivity: Informative vs Flat Priors",
          subtitle = sprintf("Flags parameters with |shift| > %.0f%% of informative-mode magnitude",
                             threshold * 100),
          x = "Parameter value (mode)", y = NULL)

      plots$mode_comparison <- .apply_meta(p,
        meta %||% diag_meta(model_name = solved$model$source_file %||% "model"))

      # Relative shift bar chart
      comp_plot2 <- comp_plot
      comp_plot2$param <- factor(comp_plot2$param, levels = comp_plot2$param[order(comp_plot2$rel_shift)])
      p2 <- ggplot2::ggplot(comp_plot2) +
        ggplot2::geom_col(ggplot2::aes(x = rel_shift, y = param, fill = status),
                          alpha = 0.8) +
        ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                           colour = "#CC3311", linewidth = 0.5) +
        ggplot2::scale_fill_manual(
          values = c("Data-driven" = "#009988", "Prior-driven" = "#CC3311"),
          name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "Prior Sensitivity: Relative Shift by Parameter",
          subtitle = sprintf("Threshold = %.0f%% (dashed line). Bars to the right are prior-driven.",
                             threshold * 100),
          x = sprintf("|shift| / |informative mode|"), y = NULL)

      plots$relative_shift <- .apply_meta(p2,
        meta %||% diag_meta(model_name = solved$model$source_file %||% "model"))
    }

    .make_result(
      result  = result,
      pass    = all_pass,
      plots   = plots,
      summary = summary_str,
      llm_summary = sprintf(
        "[%s] PriorSensitivity flagged=%d/%d threshold=%.2f",
        if (all_pass) "PASS" else "FAIL",
        n_flagged, nrow(comparisons), threshold)
    )
}
