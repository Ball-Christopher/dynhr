## R/diag-obc.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## OBC diagnostic stubs: d_obc_binding_summary()
## --------------------------------------------------------------------------


#' OBC Diagnostic: Binding periods summary
#'
#' Analyses when and for how long occasionally binding constraints were
#' active in a simulation or historical episode.  Reports binding episodes,
#' their durations, and which constraints bind in each period.
#'
#' @param regime_path Integer vector (length T). 0 = slack, non-zero indices
#'   indicate which constraint(s) are binding.
#' @param constraint_names Character vector mapping regime indices to human-
#'   readable names (e.g. \code{c("ZLB", "LTV")}). An element named "0" gives
#'   the slack-regime label.  Defaults to generic "Constraint k" names.
#' @param dates Optional Date vector (length T) for x-axis in plots.
#' @param var_paths Optional matrix (n_endo x T) of endogenous variable paths,
#'   used to overlay constraint-violating variables in plots.
#' @param var_names Optional character vector (length n_endo) of variable names
#'   corresponding to \code{var_paths} rows.
#' @param constraint_var Character: name of the variable being constrained
#'   (e.g. "r" for ZLB, "b" for LTV). Used for annotation in plots.
#' @param meta Optional \code{dynhr_diag_meta} from \code{diag_meta()}.
#' @return \code{dynhr_diagnostic} list with binding episode table and plots.
#' @noRd
d_obc_binding_summary <- function(regime_path      = NULL,
                                   constraint_names = NULL,
                                   dates            = NULL,
                                   var_paths        = NULL,
                                   var_names        = NULL,
                                   constraint_var   = NULL,
                                   meta             = NULL) {
  
    # ---- 1. Check prerequisites ----
    if (is.null(regime_path)) {
      return(.make_result(
        result = NULL, pass = NA, plots = list(),
        summary = "OBC binding summary: provide regime_path (integer vector, 0=slack, non-zero=binding)."
      ))
    }

    regime_path <- as.integer(regime_path)
    T_total <- length(regime_path)

    # ---- 2. Identify episodes ----
    is_binding <- regime_path != 0L
    if (!any(is_binding)) {
      return(.make_result(
        result = list(n_episodes = 0L, total_binding = 0L, binding_frac = 0, episodes = data.frame()),
        pass = TRUE, plots = list(),
        summary = "OBC binding summary: No binding episodes detected. Constraint never binds (PASS).",
        llm_summary = "[PASS] OBC Binding Summary status=no_binding action: The constraint never binds."
      ))
    }

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
      start_date  = if (!is.null(dates)) dates[run_starts[binding_run_idx]] else as.Date(NA),
      end_date    = if (!is.null(dates)) dates[run_ends[binding_run_idx]] else as.Date(NA),
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
      if (!is.null(constraint_names) && !is.null(names(constraint_names))) {
        nm <- constraint_names[as.character(dominant_regime)]
        episodes$regime_name[i] <- if (!is.na(nm)) nm else sprintf("Constraint_%d", dominant_regime)
      } else if (!is.null(constraint_names) && dominant_regime <= length(constraint_names)) {
        episodes$regime_name[i] <- constraint_names[dominant_regime]
      } else {
        episodes$regime_name[i] <- sprintf("Constraint_%d", dominant_regime)
      }
    }

    total_binding <- sum(is_binding)
    binding_frac <- total_binding / T_total
    avg_duration <- mean(episodes$duration)
    max_duration <- max(episodes$duration)

    # ---- 3. Per-regime breakdown ----
    regime_types <- unique(episodes$regime_name)
    n_eps <- vapply(regime_types, function(rt) sum(episodes$regime_name == rt), integer(1))
    tot_p <- vapply(regime_types, function(rt) sum(episodes$duration[episodes$regime_name == rt]), integer(1))
    avg_d <- vapply(regime_types, function(rt) mean(episodes$duration[episodes$regime_name == rt]), numeric(1))
    max_d <- vapply(regime_types, function(rt) max(episodes$duration[episodes$regime_name == rt]), integer(1))
    per_regime <- data.frame(
      constraint = regime_types, n_episodes = n_eps,
      total_periods = tot_p, avg_duration = avg_d, max_duration = max_d,
      stringsAsFactors = FALSE
    )

    # ---- 4. Pass/fail ----
    pass <- if (binding_frac < 0.25 && avg_duration < 10) TRUE
            else if (binding_frac > 0.75) FALSE
            else NA

    # ---- 5. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .gg <- getNamespace("ggplot2")

      regime_df <- data.frame(
        Period  = seq_len(T_total),
        Binding = factor(is_binding, levels = c(FALSE, TRUE), labels = c("Slack", "Binding"))
      )
      if (!is.null(dates) && length(dates) == T_total) regime_df$Period <- dates

      p_regime <- .gg$ggplot(regime_df, .gg$aes(x = .data[["Period"]], y = .data[["Binding"]])) +
        .gg$geom_step(colour = "#003B73", linewidth = 0.6) +
        .gg$theme_minimal(base_size = 10) +
        .gg$labs(title = "OBC: Binding episodes",
                 subtitle = sprintf("%d episode(s), %.1f%% binding, avg duration %.1f",
                                    n_episodes, binding_frac * 100, avg_duration),
                 x = if (!is.null(dates)) "Date" else "Period", y = "Constraint state")
      plots$regime_path <- p_regime

      if (n_episodes > 0L) {
        dur_df <- episodes
        dur_df$label <- if (!is.null(dates)) format(dur_df$start_date, "%Y-%m")
                        else sprintf("Ep_%d", dur_df$episode)
        p_dur <- .gg$ggplot(dur_df, .gg$aes(x = .data[["label"]], y = .data[["duration"]],
                                             fill = .data[["regime_name"]])) +
          .gg$geom_col(width = 0.6) + .gg$coord_flip() +
          .gg$theme_minimal(base_size = 10) +
          .gg$labs(title = "OBC: Episode durations", x = NULL, y = "Duration")
        plots$episode_durations <- p_dur
      }
    }

    # ---- 6. Summary ----
    result <- list(
      n_episodes = n_episodes, total_binding = total_binding,
      binding_frac = binding_frac, avg_duration = avg_duration,
      max_duration = max_duration, episodes = episodes, per_regime = per_regime
    )

    detail <- if (nrow(per_regime) > 0L) {
      paste(apply(per_regime, 1L, function(r) {
        sprintf("%s: %s ep, %s per, avg %.1f",
                r["constraint"], r["n_episodes"], r["total_periods"], as.numeric(r["avg_duration"]))
      }), collapse = " | ")
    } else ""

    summary_str <- sprintf(
      "OBC binding summary: %d episode(s), %.1f%% binding. %s",
      n_episodes, binding_frac * 100,
      if (total_binding > 0) sprintf("Avg duration %.1f, max %d. %s", avg_duration, max_duration, detail)
      else "Constraint never binds."
    )

    llm_summary <- paste(c(
      sprintf("OBC | Binding Summary | %s",
              if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "FAIL" else "INFO"),
      sprintf("  episodes=%d binding_frac=%.3f avg_dur=%.1f max_dur=%d",
              n_episodes, binding_frac, avg_duration, max_duration),
      sprintf("  action: %s",
              if (n_episodes == 0L) "No binding episodes."
              else if (isTRUE(pass)) "Binding episodes are short and infrequent."
              else if (isFALSE(pass)) "Constraint binds most of the time."
              else sprintf("%d binding episode(s) detected.", n_episodes))
    ), collapse = "\n")

    .make_result(result = result, pass = pass, plots = plots,
                 summary = summary_str, llm_summary = llm_summary)
}
