## R/diag-pre-d26-calibration-sensitivity.R
## --------------------------------------------------------------------------
## Phase A implementation: D26 Calibration Sensitivity Diagnostic.
##
## Implements Iskrev (2019b): reports how the identification strength of
## estimated parameters changes when calibrated parameters are perturbed.
##
## Algorithm:
##   1. Accept list of calibrated θ_c (fixed) and estimated θ_e (free)
##   2. For each calibrated parameter, perturb by ±10% (or user-specified grid)
##   3. At each perturbation, re-run D20 (Fisher information strength) for
##      estimated parameters
##   4. Report sensitivity matrix; flag fragile parameters whose identification
##      strength is sensitive to calibration assumptions
## --------------------------------------------------------------------------

#' D26. Calibration sensitivity diagnostic (Iskrev 2019b)
#'
#' Assesses how the identification strength of estimated parameters changes
#' when calibrated (fixed) parameters are perturbed. For each calibrated
#' parameter, the diagnostic re-evaluates the Fisher information-based
#' strength index (as in D20) at perturbed values of that parameter while
#' holding all others fixed.
#'
#' A parameter is flagged as "fragile" if perturbing any calibrated parameter
#' by the specified relative amount causes the estimated-parameter identification
#' strength to vary substantially.  The criterion is RELATIVE:
#'
#'   fragile_i = (s_min_i < 0.1 * s_max_possible) OR (s_max_i / s_min_i > fragility_tol)
#'
#' where s_min_i and s_max_i are the minimum and maximum |s_i| across all
#' calibration perturbations for estimated parameter i.
#'
#' The earlier \code{strength_threshold = 1.0} CRLB-style gate has been removed
#' because it conflates calibration sensitivity (D26's purpose) with intrinsic
#' identification weakness already captured by D1 and D20.  Using |s_i| < 1 as
#' the criterion caused 100% of estimated parameters to be flagged on typical
#' DSGE models (e.g. all 20/20 on the 03j model), making the diagnostic a
#' constant-fail with no discriminating power.
#'
#' The new gate:
#'   - \code{s_min < 0.1}: order-of-magnitude weakness in the perturbed baseline;
#'     indicates a calibrated parameter is so important that removing it destroys
#'     identification.
#'   - \code{s_max / s_min > fragility_tol}: identification strength varies by
#'     more than a factor of \code{fragility_tol} (default 10) across perturbations;
#'     indicates the estimated parameter's identifiability depends strongly on the
#'     exact calibration.
#'
#' @param model_solve_fn Function: theta -> named numeric vector of moments.
#'   Accepts the full parameter vector (calibrated + estimated).
#' @param theta_c Numeric vector of calibrated parameter values.
#' @param theta_e Numeric vector of estimated parameter values.
#' @param param_names_c Character vector of calibrated parameter names.
#' @param param_names_e Character vector of estimated parameter names.
#' @param perturbation_grid Named list of perturbation vectors for each
#'   calibrated parameter. If NULL, defaults to c(-0.10, 0, +0.10) relative
#'   perturbations.
#' @param eps Step size for finite differences in D20-style Jacobian.
#' @param strength_threshold (Deprecated) Not used in the gate; retained for
#'   backward compatibility.  The old CRLB-style |s_i| < 1.0 threshold has
#'   been replaced by the relative \code{fragility_tol} criterion.
#' @param fragility_tol Ratio threshold: flag if s_max / s_min > fragility_tol
#'   across perturbations.  Default 10 (order-of-magnitude swing indicates
#'   strong calibration sensitivity).  A parameter is also flagged if its
#'   minimum strength across all perturbations satisfies s_min < 0.1 (absolute
#'   near-zero strength), regardless of the ratio.  Both criteria are applied
#'   independently; a parameter is "fragile" if either fires.
#' @return dynhr_diagnostic list with results including:
#'   - sensitivity_matrix: array of strength indices (n_est x n_cal x n_pert)
#'   - fragile_params: estimated parameters flagged as fragile
#'   - fragility_table: data.frame with per-parameter sensitivity statistics
#' @references
#'   Iskrev, N. (2019b). What to expect when you're calibrating...
#'     *Journal of Economic Dynamics and Control*, forthcoming.
#' @noRd
d26_calibration_sensitivity <- function(model_solve_fn,
                                        theta_c,
                                        theta_e,
                                        param_names_c = NULL,
                                        param_names_e = NULL,
                                        perturbation_grid = NULL,
                                        eps = 1e-5,
                                        strength_threshold = 1.0,
                                        fragility_tol = 10.0,
                                        meta = NULL) {
  n_cal <- length(theta_c)
  n_est <- length(theta_e)

  if (is.null(param_names_c)) param_names_c <- names(theta_c) %||%
    paste0("theta_c_", seq_len(n_cal))
  if (is.null(param_names_e)) param_names_e <- names(theta_e) %||%
    paste0("theta_e_", seq_len(n_est))

  names(theta_c) <- param_names_c
  names(theta_e) <- param_names_e

  # Default perturbation grid: -10%, baseline, +10%
  if (is.null(perturbation_grid)) {
    perturbation_grid <- c(-0.10, 0, 0.10)
  }
  n_pert <- length(perturbation_grid)

  ## Wrap the solve-dependent computation: a failing model_solve_fn (e.g. a
  ## model that does not solve at a perturbed calibration) must yield a graceful
  ## ERROR diagnostic, not abort the whole run_all_diagnostics() pass.
  tryCatch({

    # Full parameter vector builder
    .build_full_theta <- function(th_c, th_e) {
      full <- c(th_c, th_e)
      names(full) <- c(param_names_c, param_names_e)
      full
    }

    # ---------------------------------------------------------------
    # Baseline: compute strength at unperturbed calibration
    # ---------------------------------------------------------------
    theta_full <- .build_full_theta(theta_c, theta_e)
    J_baseline <- .numerical_jacobian(
      function(th_e_only) {
        full <- c(theta_c, th_e_only)
        names(full) <- c(param_names_c, param_names_e)
        model_solve_fn(full)
      },
      theta_e, eps = eps
    )
    colnames(J_baseline) <- param_names_e

    # Baseline Fisher strength
    I_raw <- crossprod(J_baseline)
    I_reg <- I_raw + diag(1e-10, n_est)
    I_inv <- solve(I_reg)
    se_baseline <- sqrt(pmax(diag(I_inv), 0))
    s_baseline <- as.numeric(theta_e) / pmax(se_baseline, 1e-16)
    names(s_baseline) <- param_names_e

    # ---------------------------------------------------------------
    # 2. Perturb each calibrated parameter and recompute strengths
    # ---------------------------------------------------------------
    # Strength array: n_est x n_cal x n_pert
    strength_array <- array(NA_real_,
                           dim = c(n_est, n_cal, n_pert),
                           dimnames = list(
                             estimated = param_names_e,
                             calibrated = param_names_c,
                             perturbation = paste0("p", perturbation_grid)
                           ))

    # Also store the perturbation values
    pert_values <- matrix(NA_real_, nrow = n_cal, ncol = n_pert,
                          dimnames = list(param_names_c, paste0("p", perturbation_grid)))

    for (j in seq_len(n_cal)) {
      for (k in seq_len(n_pert)) {
        # Perturb calibrated parameter j by relative amount
        rel_change <- perturbation_grid[k]
        theta_c_pert <- theta_c
        theta_c_pert[j] <- theta_c[j] * (1 + rel_change)
        pert_values[j, k] <- theta_c_pert[j]

        # Compute Jacobian of estimated parameters at this calibration
        J_pert <- .numerical_jacobian(
          function(th_e_only) {
            full <- c(theta_c_pert, th_e_only)
            names(full) <- c(param_names_c, param_names_e)
            model_solve_fn(full)
          },
          theta_e, eps = eps
        )
        if (is.null(J_pert) || any(!is.finite(J_pert))) J_pert <- NULL

        if (!is.null(J_pert)) {
          colnames(J_pert) <- param_names_e
          I_raw_k <- crossprod(J_pert)
          I_reg_k <- I_raw_k + diag(1e-10, n_est)
          I_inv_k <- solve(I_reg_k)
          se_k <- sqrt(pmax(diag(I_inv_k), 0))
          s_k <- as.numeric(theta_e) / pmax(se_k, 1e-16)
          strength_array[, j, k] <- s_k
        }
      }
    }

    # ---------------------------------------------------------------
    # 3. Compute fragility statistics
    # ---------------------------------------------------------------
    fragility_list <- list()

    for (i in seq_len(n_est)) {
      # Strength min, max, range across all perturbations — operate on |s|
      s_all <- abs(as.numeric(strength_array[i, , ]))
      s_all <- s_all[is.finite(s_all)]
      if (length(s_all) == 0) next

      s_min <- min(s_all)
      s_max <- max(s_all)
      s_range <- s_max - s_min

      # Which calibrated parameters cause the largest swing?
      max_swing <- 0
      worst_cal <- NA_character_
      for (j in seq_len(n_cal)) {
        s_j <- abs(as.numeric(strength_array[i, j, ]))
        s_j <- s_j[is.finite(s_j)]
        if (length(s_j) >= 2) {
          swing <- diff(range(s_j))
          if (swing > max_swing) {
            max_swing <- swing
            worst_cal <- param_names_c[j]
          }
        }
      }

      # Relative fragility criterion (replaces CRLB |s_i| < 1.0 gate):
      #   (a) s_min < 0.1: order-of-magnitude weakness in any perturbation
      #   (b) s_max / s_min > fragility_tol: strength varies by factor > tol
      # This separates calibration sensitivity from intrinsic identification
      # weakness (the latter is already measured by D1 and D20).
      ratio <- s_max / max(s_min, 1e-16)
      fragile <- (s_min < 0.1) || (ratio > fragility_tol)
      fragility_list[[param_names_e[i]]] <- list(
        s_baseline = s_baseline[i],
        s_min = s_min,
        s_max = s_max,
        s_range = s_range,
        ratio = ratio,
        max_swing = max_swing,
        worst_calibrated = worst_cal,
        fragile = fragile
      )
    }

    # Build fragility table
    fragility_table <- do.call(rbind, lapply(seq_len(n_est), function(i) {
      fl <- fragility_list[[param_names_e[i]]]
      if (is.null(fl)) return(NULL)
      data.frame(
        parameter = param_names_e[i],
        s_baseline = fl$s_baseline,
        s_min = fl$s_min,
        s_max = fl$s_max,
        s_range = fl$s_range,
        ratio = fl$ratio,
        worst_calibrated = fl$worst_calibrated,
        fragile = fl$fragile,
        stringsAsFactors = FALSE
      )
    }))

    fragile_params <- param_names_e[vapply(fragility_list, function(fl)
      isTRUE(fl$fragile), logical(1))]
    pass <- length(fragile_params) == 0

    # ---------------------------------------------------------------
    # 4. Build plots
    # ---------------------------------------------------------------
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        requireNamespace("reshape2", quietly = TRUE)) {

      # (a) Strength heatmap: estimated params x calibrated params
      # (Take mean absolute strength across perturbations)
      strength_mean <- apply(abs(strength_array), c(1, 2), mean, na.rm = TRUE)
      dimnames(strength_mean) <- list(param_names_e, param_names_c)

      s_long <- reshape2::melt(strength_mean)
      colnames(s_long) <- c("Estimated", "Calibrated", "MeanAbsStrength")

      # Cap the colour scale at the 90th percentile so that a single
      # dominant row (high |s_i|) does not collapse all other rows into
      # a uniform dark colour.  Values above the cap are squished to the
      # maximum colour rather than shown as NA.
      strength_values <- s_long$MeanAbsStrength
      scale_upper <- stats::quantile(strength_values, 0.90, na.rm = TRUE)
      # Guard: if all values are equal (degenerate case), keep limits = NULL
      if (!is.finite(scale_upper) || scale_upper <= 0) scale_upper <- NULL

      p_sh <- ggplot2::ggplot(
        s_long,
        ggplot2::aes(x = Calibrated, y = Estimated, fill = MeanAbsStrength)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        scale_fill_dynhr_cividis(
          limits  = if (!is.null(scale_upper)) c(0, scale_upper) else NULL,
          oob     = scales::squish,
          n.breaks = 4L,
          name    = "|s_i|",
          guide   = ggplot2::guide_colourbar(barwidth = 12, barheight = 0.5)
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                              size = ggplot2::rel(0.75)),
          axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.75))
        ) +
        ggplot2::labs(
          title = "D26: Identification strength sensitivity to calibration",
          subtitle = sprintf("Mean |s_i| across %d perturbations per calibrated param",
                             n_pert),
          x = "Calibrated parameter", y = "Estimated parameter"
        )
      plots$strength_heatmap <- .apply_meta(p_sh, meta)

      # (b) Fragility bar chart (if any fragile params)
      if (nrow(fragility_table) > 0) {
        ft_plot <- fragility_table
        ft_plot$fragile_label <- ifelse(ft_plot$fragile, "Fragile", "Robust")

        p_fb <- ggplot2::ggplot(
          ft_plot,
          ggplot2::aes(x = reorder(parameter, s_range),
                       y = s_range, fill = fragile_label)
        ) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::coord_flip() +
          ggplot2::scale_fill_manual(
            values = c("Fragile" = dynhr_colours$red,
                       "Robust" = dynhr_colours$mid_blue),
            name = NULL
          ) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D26: Identification strength range across calibrations",
            subtitle = sprintf("Threshold: |s_i| < %.1f or range ratio > %.1f",
                               strength_threshold, fragility_tol),
            x = NULL, y = "Range of |s_i| across perturbations"
          )
        plots$fragility_bars <- .apply_meta(p_fb, meta)
      }
    }

    # ---------------------------------------------------------------
    # 5. Build summary
    # ---------------------------------------------------------------
    summary_text <- sprintf(
      "D26 Calibration sensitivity: %d estimated params, %d calibrated. %s%s",
      n_est, n_cal,
      if (pass) "PASS -- no fragile parameters."
      else sprintf("FAIL -- %d fragile parameter(s): %s.",
                    length(fragile_params),
                    paste(fragile_params, collapse = ", ")),
      if (nrow(fragility_table) > 0) {
        sprintf(" Min |s_i| = %.2f, Max range = %.2f.",
                min(fragility_table$s_min, na.rm = TRUE),
                max(fragility_table$s_range, na.rm = TRUE))
      } else ""
    )

    .make_result(
      result = list(
        strength_array = strength_array,
        fragility_table = fragility_table,
        fragile_params = fragile_params,
        perturbation_values = pert_values,
        s_baseline = s_baseline,
        se_baseline = se_baseline,
        J_baseline = J_baseline
      ),
      pass = pass,
      plots = plots,
      summary = summary_text,
      llm_summary = {
        badge <- if (pass) "PASS" else "FAIL"
        paste(c(
          sprintf("D26 | Calibration Sensitivity (Iskrev 2019b) | %s", badge),
          sprintf("  n_est=%d n_cal=%d n_pert=%d fragility_tol=%.0f",
                  n_est, n_cal, n_pert, fragility_tol),
          sprintf("  gate: s_min<0.1 OR s_max/s_min>%.0f (NOT |s_i|<1 CRLB -- see D1/D20)",
                  fragility_tol),
          if (nrow(fragility_table) > 0)
            sprintf("  max_ratio=%.1f min_s_min=%.3f",
                    max(fragility_table$ratio, na.rm = TRUE),
                    min(fragility_table$s_min, na.rm = TRUE)),
          if (length(fragile_params) > 0)
            sprintf("  fragile: %s", paste(fragile_params, collapse = ", ")),
          sprintf("  action: %s",
                  if (pass)
                    "Estimated parameters are robust to calibration perturbations (relative criterion)."
                  else
                    sprintf("Fragile parameters detected (%s). Perturbing calibrated parameters causes >%.0fx swing in identification strength. Perform sensitivity analysis on the listed calibrated parameters.",
                            paste(head(fragile_params, 5), collapse = ", "),
                            fragility_tol))
        ), collapse = "\n")
      }
    )
  },
  error = function(e) {
    msg <- conditionMessage(e)
    .make_result(
      result  = list(error = msg),
      pass    = NA,
      plots   = list(),
      summary = sprintf("D26 ERROR: calibration sensitivity failed: %s", msg),
      llm_summary = sprintf(
        "D26 | Calibration Sensitivity (Iskrev 2019b) | ERROR\n  %s", msg)
    )
  })
}
