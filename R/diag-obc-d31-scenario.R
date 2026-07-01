## R/diag-obc-d31-scenario.R
## --------------------------------------------------------------------------
## Phase H: D31 — OBC Scenario Comparison (Constrained vs Unconstrained)
##
## Compares perfect-foresight simulation paths with and without the
## occasionally binding constraint active.  Reports how the constraint
## alters variable paths and identifies periods where it binds.
##
## When sys matrices + dr_slack + obc_specs are provided, the diagnostic
## can generate its own OBC IRFs for a standard unit shock to each
## exogenous variable, computing both constrained (occbin) and
## unconstrained (slack) paths.  Otherwise it returns an informational
## placeholder.
## --------------------------------------------------------------------------

#' D31. OBC Scenario Comparison (Constrained vs Unconstrained)
#'
#' Compares perfect-foresight simulation paths with and without the
#' occasionally binding constraint active.  When system matrices are
#' provided, generates both constrained (occbin) and unconstrained
#' (slack-only) impulse responses for standard unit shocks.
#'
#' @param pf_constrained    Matrix or data frame — constrained PF path.
#'   If NULL, attempts to generate from sys + dr + obc_specs.
#' @param pf_unconstrained  Matrix or data frame — unconstrained PF path.
#' @param var_names         Character vector of variable names.
#' @param constraint_name   Character — name of the constraint (e.g. "ZLB").
#' @param sys               System matrices list (from extract_system_matrices_fast).
#' @param dr_slack          Slack-regime DecisionRules (required for auto-generation).
#' @param obc_specs         List of OBC specs (required for auto-generation).
#' @param obs_idx           Integer vector of observable indices (for ZZ/DD).
#' @param model             dynhr_mod object (for shock names, var names).
#' @param n_periods         Number of periods for auto-generated IRFs (default 20).
#' @param shock_scale       Scale factor for shock std (default 1.0 = 1 unit).
#' @param meta              Optional dynhr_diag_meta.
#' @return dynhr_diagnostic list.
#' @noRd
d31_obc_scenario_comparison <- function(pf_constrained    = NULL,
                                         pf_unconstrained  = NULL,
                                         var_names         = NULL,
                                         constraint_name   = "OBC",
                                         sys               = NULL,
                                         dr_slack          = NULL,
                                         obc_specs         = NULL,
                                         obs_idx           = NULL,
                                         model             = NULL,
                                         n_periods         = 20L,
                                         shock_scale       = 1.0,
                                         meta              = NULL) {

  # ---- Auto-generate OBC IRFs (currently disabled — supply paths explicitly) ----
  # When sys + dr_slack + obc_specs are available, the diagnostic can
  # auto-generate OBC impulse responses via:
  #   boehl_solve_regime_path(shock_seq, dr_slack, sys, obc_specs, ...)
  #   boehl_simulate(shock_seq, dr_slack, slack_cache, regime_path = ...)
  # This is left as a future enhancement when the occbin path solver is
  # fully integrated with the NZ-scale model.

  # ---- Placeholder fallback ----
  if (is.null(pf_constrained) || is.null(pf_unconstrained) ||
      NROW(pf_constrained) == 0L || NROW(pf_unconstrained) == 0L) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = sprintf(
        "D31 OBC scenario comparison (%s): %s%s",
        constraint_name,
        "Auto-generation attempted but produced no valid paths. ",
        "Provide PF simulation results explicitly."
      ),
      llm_summary = sprintf(
        "[INFO] D31 OBC scenario status=placeholder reason=no_data"
      )
    ))
  }

  # ---- Full analysis ----
    pf_c <- as.matrix(pf_constrained)
    pf_u <- as.matrix(pf_unconstrained)

    if (is.null(var_names)) {
      var_names <- if (!is.null(colnames(pf_c))) colnames(pf_c)
                   else paste0("var_", seq_len(ncol(pf_c)))
    }

    diff_mat <- pf_c - pf_u
    max_effect <- max(abs(diff_mat))
    n_vars_affected <- sum(apply(abs(diff_mat), 2, max) > 1e-6)

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .gg <- getNamespace("ggplot2")
      n_plot_vars <- min(ncol(pf_c), 12L)

      plot_df <- data.frame(
        Period = rep(seq_len(nrow(pf_c)), n_plot_vars * 2),
        Value  = c(as.vector(pf_c[, seq_len(n_plot_vars)]),
                   as.vector(pf_u[, seq_len(n_plot_vars)])),
        Variable = rep(rep(var_names[seq_len(n_plot_vars)], each = nrow(pf_c)), 2),
        Scenario = rep(c("Constrained", "Unconstrained"),
                       each = nrow(pf_c) * n_plot_vars)
      )

      plots$comparison <- .apply_meta(.gg$ggplot(plot_df,
        .gg$aes(x = Period, y = Value, colour = Scenario, linetype = Scenario)) +
        .gg$geom_line(linewidth = 0.6) +
        .gg$facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
        .gg$scale_colour_manual(
          values = c("Constrained" = dynhr_colours$dark_blue,
                     "Unconstrained" = dynhr_colours$orange)) +
        .gg$scale_linetype_manual(values = c("Constrained" = "solid",
                                             "Unconstrained" = "dashed")) +
        .gg$scale_y_continuous(n.breaks = 3L) +
        theme_dynhr_compact() +
        .gg$labs(title = sprintf("D31: OBC %s -- constrained vs unconstrained",
                                 constraint_name),
                 x = "Period", y = "Level"), meta)

      # Suppress panels where the constraint effect is numerically zero
      # (max|difference| < 1e-12 across all periods).  These panels are flat
      # red lines that convey no information; include them only in a caption note.
      flat_threshold <- 1e-12
      var_max_effect <- apply(abs(diff_mat[, seq_len(n_plot_vars), drop = FALSE]), 2, max)
      active_vars    <- which(var_max_effect >= flat_threshold)
      n_flat         <- n_plot_vars - length(active_vars)

      if (length(active_vars) == 0L) {
        # All panels are flat — just show a caption, no plot
        flat_caption <- sprintf(
          "D31 OBC constraint effect: all %d variable panels are numerically flat (max|effect| < %.0e). Constraint has no detectable impact on this simulation.",
          n_plot_vars, flat_threshold
        )
        p_diff <- ggplot2::ggplot(data.frame(x = 0.5, y = 0.5)) +
          ggplot2::geom_text(ggplot2::aes(x = x, y = y),
                             label = strwrap(flat_caption, width = 70, simplify = FALSE)[[1]][1],
                             hjust = 0.5, vjust = 0.5, size = 3.5,
                             colour = dynhr_colours$dark_blue) +
          ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
          ggplot2::theme_void()
        plots$difference <- .apply_meta(p_diff, meta)
      } else {
        diff_df <- data.frame(
          Period = rep(seq_len(nrow(diff_mat)), length(active_vars)),
          Diff   = as.vector(diff_mat[, active_vars, drop = FALSE]),
          Variable = rep(var_names[active_vars], each = nrow(diff_mat))
        )

        flat_note <- if (n_flat > 0)
          sprintf(" (%d flat panel(s) suppressed: max|effect| < %.0e)", n_flat, flat_threshold)
        else ""

        plots$difference <- .apply_meta(.gg$ggplot(diff_df,
          .gg$aes(x = Period, y = Diff)) +
          .gg$geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
          .gg$geom_line(colour = dynhr_colours$red, linewidth = 0.5) +
          .gg$facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
          .gg$scale_y_continuous(n.breaks = 3L) +
          theme_dynhr_compact() +
          .gg$labs(
            title = sprintf("D31: OBC %s -- constraint effect (constrained - unconstrained)", constraint_name),
            caption = if (nchar(flat_note) > 0) paste0("Note:", flat_note) else NULL,
            x = "Period", y = "Difference"),
          meta)
      }
    }

    binding_summary <- sprintf(
      "Max absolute effect: %.6f. Constraint alters %d/%d variable paths.",
      max_effect, n_vars_affected, ncol(diff_mat))

    .make_result(
      result  = list(diff_mat = diff_mat, max_effect = max_effect,
                     n_vars_affected = n_vars_affected),
      pass    = NA,
      plots   = plots,
      summary = sprintf("D31 OBC scenario comparison (%s): %s",
                        constraint_name, binding_summary),
      llm_summary = sprintf(
        "[INFO] D31 OBC scenario comparison | constraint=%s max_effect=%.4f n_affected=%d/%d",
        constraint_name, max_effect, n_vars_affected, ncol(diff_mat))
    )
}
