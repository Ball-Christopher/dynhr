## R/diag-obc-d32-binding.R
## --------------------------------------------------------------------------
## Phase H: D32 -- OBC Binding Periods Summary
##
## Analyses when and for how long occasionally binding constraints were
## active in a simulation or historical episode.  Reports binding episodes,
## their durations, and which constraints bind in each period.
##
## When sys + dr_slack + obc_specs are provided, can auto-generate a
## regime path from a standard unit-shock IRF.
## --------------------------------------------------------------------------

#' D32. OBC Binding Periods Summary
#'
#' Analyses when and for how long occasionally binding constraints were
#' active.  Reports binding episodes, their durations, and per-constraint
#' breakdown.  When system matrices are available, auto-generates a regime
#' path from occbin IRF computation.
#'
#' @param regime_path     Integer vector (length T). 0 = slack, >0 = binding.
#' @param constraint_names Character vector mapping regime indices to names.
#'   An element named "0" gives the slack-regime label.  Default "Slack".
#' @param dates           Optional Date vector (length T) for x-axis.
#' @param var_paths       Optional matrix (n_endo x T) of endogenous paths.
#' @param var_names       Optional character vector of variable names.
#' @param constraint_var  Character: name of constrained variable for annotation.
#' @param sys             System matrices (for auto-generation).
#' @param dr_slack        Slack-regime DecisionRules (for auto-generation).
#' @param obc_specs       List of OBC specs (for auto-generation).
#' @param obs_idx         Observable indices (for auto-generation).
#' @param model           dynhr_mod object.
#' @param n_periods       Number of periods for auto-generated IRF (default 20).
#' @param meta            Optional dynhr_diag_meta.
#' @return dynhr_diagnostic list with binding episode table and plots.
#' @noRd
d32_obc_binding_summary <- function(regime_path      = NULL,
                                     constraint_names = NULL,
                                     dates            = NULL,
                                     var_paths        = NULL,
                                     var_names        = NULL,
                                     constraint_var   = NULL,
                                     sys              = NULL,
                                     dr_slack         = NULL,
                                     obc_specs        = NULL,
                                     obs_idx          = NULL,
                                     model            = NULL,
                                     n_periods        = 20L,
                                     meta             = NULL) {

  # ---- Auto-generate regime path if data provided ----
  if (is.null(regime_path) && !is.null(sys) && !is.null(dr_slack) &&
      !is.null(obc_specs) && length(obc_specs) > 0L) {
    if (is.null(obs_idx) && !is.null(model)) {
      n_endo <- nrow(dr_slack$ghx)
      obs_idx <- seq_len(n_endo)
    }
    n_shock <- ncol(dr_slack$ghu) %||% 0L
    if (n_shock > 0) {
      shock_seq <- matrix(0, nrow = n_shock, ncol = n_periods)
      shock_seq[1, 1] <- 1.0  # unit shock to first exo variable

      # Solve occbin path -- returns regime path from guess-and-verify
      occbin_res <- boehl_solve_regime_path(shock_seq, dr_slack, sys, obc_specs,
                                             obs_idx = obs_idx, max_iter = 30L)
      if (!is.null(occbin_res)) {
        regime_path <- occbin_res$regime_path
        if (!is.null(occbin_res$paths) && is.null(var_paths)) {
          var_paths <- occbin_res$paths
        }
      }
    }
      if (is.null(constraint_var) && length(obc_specs) > 0L) {
        constraint_var <- obc_specs[[1]]$var_name %||% constraint_var
      }
      if (is.null(var_names) && !is.null(model)) {
        var_names <- model$var_names
      }
  }

  # ---- Placeholder ----
  if (is.null(regime_path)) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = paste(
        "D32 OBC binding summary: NOT YET AVAILABLE.",
        "Provide regime_path (integer vector, 0=slack, non-zero=binding)",
        "or supply sys + dr_slack + obc_specs for auto-generation."
      ),
      llm_summary = "[INFO] D32 OBC binding status=placeholder reason=no_data"
    ))
  }

  # ---- Full analysis ----
  regime_path <- as.integer(regime_path)
  T_total <- length(regime_path)
  is_binding <- regime_path != 0L

  if (!any(is_binding)) {
    return(.make_result(
      result  = list(n_episodes = 0L, total_binding = 0L,
                     binding_frac = 0, episodes = data.frame()),
      pass    = TRUE,
      plots   = list(),
      summary = "D32 OBC binding summary: No binding episodes. (PASS)",
      llm_summary = "[PASS] D32 OBC binding status=no_binding"
      ))
    }

    # Default constraint names
    if (is.null(constraint_names)) {
      n_specs <- length(obc_specs %||% list())
      if (n_specs > 0) {
        cnames <- vapply(obc_specs, function(s) s$name %||% "OBC",
                         character(1))
        constraint_names <- setNames(as.list(seq_along(cnames)), cnames)
      } else {
        constraint_names <- c("0" = "Slack", "1" = "Binding")
      }
    }

    # Identify binding episodes via rle
    runs <- rle(is_binding)
    run_ends <- cumsum(runs$lengths)
    run_starts <- c(1L, run_ends[-length(run_ends)] + 1L)
    binding_run_idx <- which(runs$values)
    n_episodes <- length(binding_run_idx)

    episodes <- data.frame(
      episode     = seq_len(n_episodes),
      start_idx   = run_starts[binding_run_idx],
      end_idx     = run_ends[binding_run_idx],
      duration    = runs$lengths[binding_run_idx],
      start_date  = if (!is.null(dates))
                      dates[run_starts[binding_run_idx]] else as.Date(NA),
      end_date    = if (!is.null(dates))
                      dates[run_ends[binding_run_idx]] else as.Date(NA),
      regime_val  = NA_integer_,
      regime_name = NA_character_,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(n_episodes)) {
      ep_path <- regime_path[episodes$start_idx[i]:episodes$end_idx[i]]
      ep_bind <- ep_path[ep_path != 0L]
      tbl <- sort(table(ep_bind), decreasing = TRUE)
      dominant_regime <- as.integer(names(tbl)[1L])
      episodes$regime_val[i] <- dominant_regime
      nm <- constraint_names[[as.character(dominant_regime)]]
      episodes$regime_name[i] <- nm %||% sprintf("Constraint_%d", dominant_regime)
    }

    total_binding <- sum(is_binding)
    binding_frac <- total_binding / T_total

    # Degenerate case: always binding (100%)
    # A solid block of one colour conveys no information about regime variation.
    # Replace with an explanatory annotation plot.
    always_binding <- all(is_binding)

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .gg <- getNamespace("ggplot2")

      if (always_binding) {
        degenerate_msg <- sprintf(
          "D32: Constraint is binding in ALL %d periods (binding_frac = 1.0).\nThe timeline is a solid block -- no regime variation to display.\nThis may indicate the floor level is too loose or the shock is too large.",
          T_total
        )
        plots$regime_timeline <- .apply_meta(
          ggplot2::ggplot() +
            ggplot2::annotate(
              "text", x = 0.5, y = 0.5,
              label = degenerate_msg,
              hjust = 0.5, vjust = 0.5, size = 3.8,
              colour = dynhr_colours$red
            ) +
            ggplot2::theme_void(),
          meta
        )
      } else {
      regime_df <- data.frame(
        Period = seq_len(T_total),
        Regime = factor(regime_path,
                        levels = sort(unique(regime_path)))
      )
      plots$regime_timeline <- .apply_meta(.gg$ggplot(regime_df,
        .gg$aes(x = Period, y = 1, fill = Regime)) +
        .gg$geom_tile(height = 0.8) +
        .gg$scale_fill_manual(
          values = c("0" = dynhr_colours$dark_blue, "1" = dynhr_colours$red,
                     "2" = dynhr_colours$orange,    "3" = dynhr_colours$green)) +
        .gg$scale_y_continuous(breaks = NULL) +
        .gg$labs(title = "D32: OBC binding regime timeline",
                 x = "Period", y = NULL, fill = "Regime") +
        theme_dynhr_diagnostic() +
        .gg$theme(axis.text.y = ggplot2::element_blank(),
                  panel.grid = ggplot2::element_blank()), meta)
      }  # end else (not always_binding)

      if (!is.null(var_paths) && !is.null(constraint_var)) {
        var_idx <- match(constraint_var, var_names)
        if (!is.na(var_idx) && var_idx <= nrow(var_paths)) {
          cv_df <- data.frame(
            Period = seq_len(ncol(var_paths)),
            Value  = as.numeric(var_paths[var_idx, ]),
            Binding = is_binding
          )
          plots$constraint_var <- .apply_meta(.gg$ggplot(cv_df,
            .gg$aes(x = Period, y = Value, colour = Binding)) +
            .gg$geom_line(linewidth = 0.6) +
            .gg$geom_hline(yintercept = if (length(obc_specs) > 0)
                             obc_specs[[1]]$bound else 0,
                           linetype = "dashed", colour = "grey50") +
            .gg$scale_colour_manual(values = c("FALSE" = dynhr_colours$dark_blue,
                                               "TRUE" = dynhr_colours$red)) +
            theme_dynhr_diagnostic() +
            .gg$labs(title = sprintf("D32: Constrained variable (%s) with binding periods",
                                     constraint_var),
                     x = "Period", y = constraint_var), meta)
        }
      }
    }

    # Summary
    degenerate_note <- if (always_binding) {
      " [DEGENERATE: always binding -- check floor level or shock size.]"
    } else ""

    summary_str <- sprintf(
      "D32 OBC binding summary: %d episode(s), %d/%d periods binding (%.1f%%, binding_frac=%.3f).%s",
      n_episodes, total_binding, T_total, 100 * binding_frac, binding_frac, degenerate_note
    )

    .make_result(
      result  = list(episodes = episodes, n_episodes = n_episodes,
                     total_binding = total_binding,
                     binding_frac = binding_frac,
                     always_binding = always_binding),
      pass    = NA,
      plots   = plots,
      summary = summary_str,
      llm_summary = sprintf(
        "[INFO] D32 OBC binding summary | n_episodes=%d binding_frac=%.3f%s",
        n_episodes, binding_frac,
        if (always_binding) " degenerate=always_binding" else "")
    )
}
