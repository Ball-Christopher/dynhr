## R/diag-model-compare.R
## --------------------------------------------------------------------------
## Phase-3+ addition.
##
## model_comparison() -- compare two or more estimated models.
## Covers: parameter posteriors (ridge), IRF comparison, moment comparison,
## model evidence (log marginal likelihoods / Bayes factors).
## --------------------------------------------------------------------------

#' Compare two or more estimated dynhr models
#'
#' Accepts a named list of model result objects and produces side-by-side
#' comparisons of posterior parameter estimates, impulse responses, second
#' moments, and log marginal likelihoods.
#'
#' Each element of \code{models} should be a named list with any subset of:
#' \describe{
#'   \item{\code{draws}}{Posterior draws matrix (n_draws ?-- n_params) with
#'     named columns.}
#'   \item{\code{mode}}{List with \code{$theta_mode} and \code{$logpost}.}
#'   \item{\code{irfs}}{Named list of IRF matrices (shock -> T?--n_vars), or a
#'     long data frame with columns \code{horizon, variable, shock, value}.}
#'   \item{\code{moments}}{List with \code{$std_dev} (named numeric).}
#'   \item{\code{log_marglik}}{Scalar log marginal likelihood.}
#' }
#'
#' @param models      Named list of model result objects (see above).  At
#'   least two entries required.
#' @param param_sel   Character vector of parameter names to compare.
#'   \code{NULL} (default) uses the union of parameter names across all models.
#' @param shock_sel   Shock names to compare in IRF panels (\code{NULL} = all).
#' @param var_sel     Variable names for IRF and moment panels (\code{NULL} = all,
#'   capped at 8 for readability).
#' @param obs_names   Observable names for moment comparison (\code{NULL} = all).
#' @param data_moments Named numeric vector of empirical SDs for moment
#'   comparison.  \code{NULL} suppresses the data reference line.
#' @param ci_level    Credible interval level for the parameter forest plot
#'   (default 0.90).
#' @param irf_periods Number of horizons in IRF comparison (default 20).
#' @param meta        A \code{\link{diag_meta}} object for plot provenance.
#' @return dynhr_diagnostic list with plots:
#'   \code{parameter_ridge}, \code{parameter_forest}, \code{irfs_<shock>},
#'   \code{moments}, \code{model_evidence}.
#' @noRd
model_comparison <- function(models,
                              param_sel     = NULL,
                              shock_sel     = NULL,
                              var_sel       = NULL,
                              obs_names     = NULL,
                              data_moments  = NULL,
                              ci_level      = 0.90,
                              irf_periods   = 20L,
                              meta          = NULL) {

    if (length(models) < 2L)
      stop("model_comparison() requires at least 2 models.")
    model_names <- names(models)
    if (is.null(model_names))
      model_names <- paste0("Model_", seq_along(models))

    alpha <- (1 - ci_level) / 2

    # ---- Colour palette per model ------------------------------------------
    n_models <- length(models)
    model_pal <- setNames(
      if (n_models <= length(dynhr_palette)) dynhr_palette[seq_len(n_models)]
      else colorRampPalette(dynhr_palette)(n_models),
      model_names
    )

    plots <- list()

    # ========================================================================
    # 1. Parameter comparison -- ridge density + forest plot
    # ========================================================================
    draws_list <- Filter(Negate(is.null),
                         lapply(models, function(m) m$draws))

    if (length(draws_list) >= 2L) {
      # Union of parameter names
      all_params <- unique(unlist(lapply(draws_list, colnames)))
      if (!is.null(param_sel)) all_params <- intersect(param_sel, all_params)

      # ---- Ridge density ---------------------------------------------------
      dens_rows <- do.call(rbind, lapply(model_names, function(mn) {
        d <- models[[mn]]$draws
        if (is.null(d)) return(NULL)
        d <- as.matrix(d)
        pars <- intersect(all_params, colnames(d))
        do.call(rbind, lapply(pars, function(p) {
          x <- d[, p]
          x <- x[is.finite(x)]
          if (length(x) < 4L) return(NULL)
          kde <- density(x, n = 256)
          data.frame(x = kde$x, density = kde$y,
                     model = mn, param = p,
                     stringsAsFactors = FALSE)
        }))
      }))

      if (!is.null(dens_rows) && nrow(dens_rows) > 0) {
        n_p <- length(all_params)
        dens_rows$param <- factor(dens_rows$param, levels = rev(all_params))

        p_ridge <- ggplot2::ggplot(dens_rows,
                                   ggplot2::aes(x     = x,
                                                fill  = model,
                                                colour = model)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = density),
                               alpha = 0.40, colour = NA) +
          ggplot2::geom_line(ggplot2::aes(y = density), linewidth = 0.40) +
          ggplot2::facet_grid(rows   = ggplot2::vars(param),
                              scales = "free",
                              switch = "y") +
          ggplot2::scale_fill_manual(values = model_pal, name = NULL) +
          ggplot2::scale_colour_manual(values = model_pal, name = NULL) +
          theme_dynhr_diagnostic() +
          ggplot2::theme(
            axis.title.y     = ggplot2::element_blank(),
            axis.text.y      = ggplot2::element_blank(),
            axis.ticks.y     = ggplot2::element_blank(),
            strip.text.y.left = ggplot2::element_text(
              angle = 0, hjust = 1,
              size  = ggplot2::rel(if (n_p > 15) 0.65 else 0.80)),
            panel.spacing    = ggplot2::unit(0.15, "lines")
          ) +
          ggplot2::labs(
            title    = "Model comparison: posterior parameter distributions",
            subtitle = sprintf("%d models * %d parameters", n_models, n_p),
            x = "Parameter value"
          )
        p_ridge <- .apply_meta(p_ridge, meta)
        plots$parameter_ridge <- p_ridge
      }

      # ---- Forest plot (median + CI per model) -----------------------------
      forest_rows <- do.call(rbind, lapply(model_names, function(mn) {
        d <- models[[mn]]$draws
        if (is.null(d)) return(NULL)
        d <- as.matrix(d)
        pars <- intersect(all_params, colnames(d))
        do.call(rbind, lapply(pars, function(p) {
          x <- d[, p]
          x <- x[is.finite(x)]
          if (length(x) < 4L) return(NULL)
          data.frame(
            model  = mn,
            param  = p,
            median = median(x),
            lo     = quantile(x, alpha),
            hi     = quantile(x, 1 - alpha),
            stringsAsFactors = FALSE
          )
        }))
      }))

      if (!is.null(forest_rows) && nrow(forest_rows) > 0) {
        forest_rows$param <- factor(forest_rows$param,
                                    levels = rev(all_params))
        p_forest <- ggplot2::ggplot(forest_rows,
                                    ggplot2::aes(x = median, y = param,
                                                 colour = model, shape = model)) +
          ggplot2::geom_point(
            position = ggplot2::position_dodge(width = 0.55), size = 2) +
          ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi),
            position = ggplot2::position_dodge(width = 0.55),
            orientation = "y", width = 0.2, linewidth = 0.5) +
          ggplot2::scale_colour_manual(values = model_pal, name = NULL) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title    = "Model comparison: parameter estimates",
            subtitle = sprintf("%d%% credible intervals", round(ci_level * 100)),
            x = "Parameter value", y = NULL
          )
        p_forest <- .apply_meta(p_forest, meta)
        plots$parameter_forest <- p_forest
      }
    }

    # ========================================================================
    # 2. IRF comparison
    # ========================================================================
    irf_list <- Filter(Negate(is.null), lapply(models, function(m) m$irfs))

    if (length(irf_list) >= 2L) {
      # Normalise all to long format
      long_irfs <- do.call(rbind, lapply(model_names, function(mn) {
        irfs <- models[[mn]]$irfs
        if (is.null(irfs)) return(NULL)
        df <- .get_irfs_long(irfs)
        if (is.null(df)) return(NULL)
        df <- df[df$horizon <= irf_periods, ]
        df$model <- mn
        df
      }))

      if (!is.null(long_irfs) && nrow(long_irfs) > 0) {
        shocks_avail <- unique(long_irfs$shock)
        vars_avail   <- unique(long_irfs$variable)
        if (!is.null(shock_sel)) shocks_avail <- intersect(shock_sel, shocks_avail)
        if (!is.null(var_sel))   vars_avail   <- intersect(var_sel,   vars_avail)
        vars_avail <- head(vars_avail, 8L)

        for (sh in shocks_avail) {
          df_sh <- long_irfs[long_irfs$shock == sh &
                               long_irfs$variable %in% vars_avail, ]
          if (nrow(df_sh) == 0L) next

          p_irf <- ggplot2::ggplot(df_sh,
                                   ggplot2::aes(x = horizon, y = value,
                                                colour = model,
                                                linetype = model)) +
            ggplot2::geom_line(linewidth = 0.65, na.rm = TRUE) +
            ggplot2::geom_hline(yintercept = 0,
                                colour = dynhr_colours$grey, linewidth = 0.25) +
            ggplot2::facet_wrap(~ variable, scales = "free_y",
                                ncol = min(3L, length(vars_avail))) +
            ggplot2::scale_colour_manual(values = model_pal, name = NULL) +
            ggplot2::scale_linetype_manual(
              values = setNames(c("solid","dashed","dotted","dotdash",
                                  "longdash","twodash")[seq_len(n_models)],
                                model_names), name = NULL) +
            theme_dynhr_diagnostic() +
            ggplot2::labs(
              title    = sprintf("Model comparison: IRF -- shock = %s", sh),
              x = "Horizon (quarters)", y = "Response"
            )
          p_irf <- .apply_meta(p_irf, meta)
          plots[[paste0("irfs_", sh)]] <- p_irf
        }
      }
    }

    # ========================================================================
    # 3. Moment comparison
    # ========================================================================
    moment_rows <- do.call(rbind, lapply(model_names, function(mn) {
      mom <- models[[mn]]$moments
      if (is.null(mom)) return(NULL)
      sds <- if (!is.null(mom$std_dev)) mom$std_dev
             else if (!is.null(mom$moments$std_dev)) mom$moments$std_dev
             else NULL
      if (is.null(sds)) return(NULL)
      if (!is.null(obs_names)) sds <- sds[intersect(obs_names, names(sds))]
      data.frame(model = mn, variable = names(sds), sd = as.numeric(sds),
                 stringsAsFactors = FALSE)
    }))

    if (!is.null(moment_rows) && nrow(moment_rows) > 0) {
      # Add data reference if provided
      ref_rows <- if (!is.null(data_moments)) {
        nm_sel <- intersect(names(data_moments), unique(moment_rows$variable))
        data.frame(model = "Data", variable = nm_sel,
                   sd = as.numeric(data_moments[nm_sel]),
                   stringsAsFactors = FALSE)
      } else NULL

      all_moment_df <- rbind(moment_rows, ref_rows)
      all_moment_df$model <- factor(all_moment_df$model,
                                    levels = c(model_names, if (!is.null(ref_rows)) "Data"))
      mom_pal_ext <- c(model_pal, if (!is.null(ref_rows)) c(Data = dynhr_colours$grey))

      p_mom <- ggplot2::ggplot(all_moment_df,
                               ggplot2::aes(x = sd, y = variable,
                                            colour = model, shape = model)) +
        ggplot2::geom_point(
          position = ggplot2::position_dodge(width = 0.55), size = 2.5) +
        ggplot2::scale_colour_manual(values = mom_pal_ext, name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title    = "Model comparison: theoretical standard deviations",
          subtitle = if (!is.null(ref_rows)) "Grey = data SD" else NULL,
          x = "Standard deviation", y = NULL
        )
      p_mom <- .apply_meta(p_mom, meta)
      plots$moments <- p_mom
    }

    # ========================================================================
    # 4. Model evidence table
    # ========================================================================
    log_mls <- vapply(model_names, function(mn) {
      v <- models[[mn]]$log_marglik
      if (is.null(v)) NA_real_ else as.numeric(v)
    }, numeric(1))

    evidence_tbl <- NULL
    if (any(is.finite(log_mls))) {
      best_lml    <- max(log_mls, na.rm = TRUE)
      log_bfs     <- log_mls - best_lml
      # Posterior model probabilities (uniform prior)
      n_finite    <- sum(is.finite(log_bfs))
      log_pp      <- ifelse(is.finite(log_bfs), log_bfs, -1e300)
      log_pp      <- log_pp - max(log_pp)
      pp          <- exp(log_pp) / sum(exp(log_pp))

      evidence_tbl <- data.frame(
        model        = model_names,
        log_marglik  = log_mls,
        log_BF       = log_bfs,
        post_prob    = pp,
        stringsAsFactors = FALSE
      )

      if (requireNamespace("ggplot2", quietly = TRUE)) {
        ev_df <- evidence_tbl[is.finite(evidence_tbl$log_marglik), ]
        ev_df$model <- factor(ev_df$model, levels = ev_df$model)
        p_ev <- ggplot2::ggplot(ev_df, ggplot2::aes(x = model, y = log_BF,
                                                      fill = model)) +
          ggplot2::geom_col(width = 0.6) +
          ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", log_BF),
                                           vjust = -0.3),
                             size = 3.5, colour = dynhr_colours$dark_blue) +
          ggplot2::scale_fill_manual(values = model_pal, guide = "none") +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title    = "Model comparison: log Bayes factors vs best model",
            subtitle = "Reference = 0 (best model)",
            x = NULL, y = "Log Bayes factor"
          )
        p_ev <- .apply_meta(p_ev, meta)
        plots$model_evidence <- p_ev
      }
    }

    # ========================================================================
    # Summary
    # ========================================================================
    n_draws_avail <- sum(sapply(models, function(m) !is.null(m$draws)))
    n_irfs_avail  <- sum(sapply(models, function(m) !is.null(m$irfs)))
    n_mom_avail   <- sum(sapply(models, function(m) !is.null(m$moments)))

    ev_str <- if (!is.null(evidence_tbl) && any(is.finite(evidence_tbl$log_BF))) {
      best_nm <- evidence_tbl$model[which.max(evidence_tbl$log_marglik)]
      sprintf(" Best evidence: %s.", best_nm)
    } else ""

    .make_result(
      result  = list(evidence_tbl  = evidence_tbl,
                     model_names   = model_names,
                     param_sel     = param_sel,
                     n_draws_avail = n_draws_avail,
                     n_irfs_avail  = n_irfs_avail),
      pass    = NA,
      plots   = plots,
      summary = sprintf(
        "Model comparison: %d models. Draws: %d/%d * IRFs: %d/%d * Moments: %d/%d.%s",
        n_models, n_draws_avail, n_models,
        n_irfs_avail,  n_models,
        n_mom_avail,   n_models,
        ev_str
      )
    )
}


# ---------------------------------------------------------------------------
#' Compare an extended model to its parent (genealogy comparison)
#'
#' A focused wrapper around \code{model_comparison()} designed for
#' model-genealogy pipelines where a \emph{child} model adds parameters or
#' equations to a \emph{parent} model.  The function:
#'
#' \enumerate{
#'   \item Overlays posteriors for shared parameters.
#'   \item Reports whether the child's posteriors for shared parameters differ
#'     substantially from the parent (shift > 1 posterior SD in either model).
#'   \item Reports the log Bayes factor \eqn{\log p(Y|\text{child}) -
#'     \log p(Y|\text{parent})} if marginal likelihoods are available.
#' }
#'
#' @param chains        \code{dynhr_chains} (or draws matrix) from the child
#'   model.
#' @param parent_chains \code{dynhr_chains} (or draws matrix) from the parent
#'   model.
#' @param child_name    Human-readable label for the child model (default
#'   \code{"child"}).
#' @param parent_name   Human-readable label for the parent model (default
#'   \code{"parent"}).
#' @param param_sel     Parameter names to compare; \code{NULL} (default) uses
#'   the intersection of parameter names between the two models.
#' @param shift_threshold Posterior-SD multiples above which a shift is
#'   flagged as substantial (default 1.0).
#' @param meta          Optional provenance descriptor from
#'   \code{\link{diag_meta}}.
#' @return A \code{dynhr_diagnostic} object with pass = TRUE if the child
#'   improves on the parent (positive log BF) or there is no evidence of
#'   regression in shared parameters, pass = FALSE if the child model's
#'   shared parameters shift substantially relative to the parent and
#'   marginal likelihood decreases, and pass = NA when marginal likelihoods
#'   are not available.
#' @seealso \code{\link{run_diagnostics}}, \code{model_comparison}
#' @export
diag_parent_comparison <- function(chains, parent_chains,
                                    child_name  = "child",
                                    parent_name = "parent",
                                    param_sel   = NULL,
                                    shift_threshold = 1.0,
                                    meta = NULL) {

    # Extract draw matrices
    child_draws  <- if (inherits(chains,        "dynhr_chains")) chains$chain
                    else as.matrix(chains)
    parent_draws <- if (inherits(parent_chains, "dynhr_chains")) parent_chains$chain
                    else as.matrix(parent_chains)

    # Shared parameters
    shared_params <- intersect(colnames(child_draws), colnames(parent_draws))
    if (!is.null(param_sel))
      shared_params <- intersect(param_sel, shared_params)

    if (length(shared_params) == 0L) {
      return(.make_result(
        pass    = NA,
        summary = sprintf(
          "diag_parent_comparison [%s vs %s]: no shared parameters to compare",
          child_name, parent_name)
      ))
    }

    # Compute posterior summaries for shared parameters
    .post_summary <- function(draws, pnames) {
      means <- colMeans(draws[, pnames, drop = FALSE])
      sds   <- apply(draws[, pnames, drop = FALSE], 2L, sd)
      list(mean = means, sd = sds)
    }

    child_post  <- .post_summary(child_draws,  shared_params)
    parent_post <- .post_summary(parent_draws, shared_params)

    # Flag parameters with shifts > shift_threshold * avg-SD
    avg_sd <- (child_post$sd + parent_post$sd) / 2
    shifts <- abs(child_post$mean - parent_post$mean) / avg_sd
    flagged <- names(shifts[shifts > shift_threshold])

    # Extract log marginal likelihoods if available
    log_ml_child  <- if (inherits(chains,        "dynhr_chains")) chains$log_marginal_lik  else NULL
    log_ml_parent <- if (inherits(parent_chains, "dynhr_chains")) parent_chains$log_marginal_lik else NULL

    log_bf <- if (!is.null(log_ml_child) && !is.null(log_ml_parent))
      log_ml_child - log_ml_parent else NA_real_

    # Build parameter comparison table
    param_table <- data.frame(
      param       = shared_params,
      parent_mean = round(parent_post$mean[shared_params], 4),
      parent_sd   = round(parent_post$sd[shared_params],   4),
      child_mean  = round(child_post$mean[shared_params],  4),
      child_sd    = round(child_post$sd[shared_params],    4),
      shift_sds   = round(shifts[shared_params],            2),
      flagged     = shared_params %in% flagged,
      stringsAsFactors = FALSE
    )

    # Verdict
    has_ml <- !is.na(log_bf)
    if (has_ml) {
      bf_10   <- log_bf / log(10)
      bf_str  <- sprintf("  log10 BF = %.2f (%s)\n",
                         bf_10,
                         if      (bf_10 > 2)   "decisive evidence for child"
                         else if (bf_10 > 1)   "strong evidence for child"
                         else if (bf_10 > 0.5) "substantial evidence for child"
                         else if (bf_10 > 0)   "barely worth mentioning"
                         else if (bf_10 > -0.5)"barely worth mentioning against"
                         else if (bf_10 > -1)  "substantial evidence for parent"
                         else if (bf_10 > -2)  "strong evidence for parent"
                         else                  "decisive evidence for parent")
      pass <- log_bf > 0
    } else {
      bf_str <- "  Log marginal likelihood not available (run SMC for exact BF)\n"
      pass   <- NA
    }

    summary_txt <- paste0(
      sprintf("diag_parent_comparison [%s vs %s]: %d shared parameters\n",
              child_name, parent_name, length(shared_params)),
      bf_str,
      if (length(flagged) > 0)
        sprintf("  Substantial shifts (> %.1f SD): %s\n",
                shift_threshold, paste(flagged, collapse = ", "))
      else
        "  No substantial parameter shifts detected\n"
    )

    # Summary plot (only if ggplot2 available)
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      df_plot <- rbind(
        data.frame(model = parent_name, param = shared_params,
                   mean  = parent_post$mean[shared_params],
                   sd    = parent_post$sd[shared_params],
                   stringsAsFactors = FALSE),
        data.frame(model = child_name,  param = shared_params,
                   mean  = child_post$mean[shared_params],
                   sd    = child_post$sd[shared_params],
                   stringsAsFactors = FALSE)
      )
      df_plot$model <- factor(df_plot$model,
                               levels = c(parent_name, child_name))
      model_pal <- c("#1B7CB6", "#E8540A")
      names(model_pal) <- c(parent_name, child_name)

      p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = mean, y = param, colour = model,
                     xmin = mean - 2 * sd, xmax = mean + 2 * sd)
      ) +
        ggplot2::geom_pointrange(
          ggplot2::aes(shape = model),
          position = ggplot2::position_dodge(width = 0.4),
          size = 0.6, linewidth = 0.65
        ) +
        ggplot2::scale_colour_manual(values = model_pal, name = NULL) +
        ggplot2::scale_shape_manual(
          values = setNames(c(16L, 17L), c(parent_name, child_name)),
          name = NULL
        ) +
        ggplot2::labs(
          title    = sprintf("Parent vs child: %s vs %s", parent_name, child_name),
          subtitle = sprintf("%d shared parameters -- mean +/- 2 SD",
                             length(shared_params)),
          x = "Posterior estimate",
          y = NULL
        )
      if (exists("theme_dynhr_diagnostic"))
        p <- p + theme_dynhr_diagnostic()

      p <- if (!is.null(meta)) .apply_meta(p, meta) else p
      plots$parameter_forest <- p
    }

    .make_result(
      pass    = pass,
      summary = summary_txt,
      result  = list(param_table   = param_table,
                     log_bf        = log_bf,
                     flagged       = flagged,
                     child_name    = child_name,
                     parent_name   = parent_name),
      plots   = plots
    )

}
