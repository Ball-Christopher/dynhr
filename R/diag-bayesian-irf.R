## R/diag-bayesian-irf.R
## --------------------------------------------------------------------------
## Phase-3+ addition.
##
## diag_bayesian_irf() -- posterior IRF credible bands by re-solving the
## model at a subsample of MCMC draws.
## --------------------------------------------------------------------------

#' Bayesian IRF credible bands from posterior draws
#'
#' Re-solves the model at a subsample of posterior draws and collects impulse
#' responses.  Plots the posterior median plus credible bands (default 10-90%
#' and 16-84%) against the point estimate at the mode.
#'
#' Computational note: each draw requires a full perturbation solve
#' (\code{stoch_simul}).  Use \code{n_subsample} to trade resolution for
#' speed.  500 draws at a 24-variable model takes ~2-3 minutes.
#'
#' @param model       dynhr_mod from \code{parse_mod()} (used for parameter
#'   names and stoch_simul).
#' @param compiled    dynhr_compiled from \code{compile_model()}.
#' @param draws       Posterior draws matrix (n_draws ?-- n_params); column
#'   names must match estimated parameter names.
#' @param irf_periods Number of impulse-response horizons (default 20).
#' @param n_subsample Number of draws to use (default 400; set lower for speed).
#' @param ci_bands    Two-element numeric vector of lower/upper quantile
#'   probabilities for the outer credible band (default \code{c(0.10, 0.90)}).
#' @param inner_bands Inner credible band (default \code{c(0.16, 0.84)}).
#' @param shock_sel   Character vector of shock names to plot.  \code{NULL}
#'   (default) plots all shocks.
#' @param var_sel     Character vector of variable names to plot.  \code{NULL}
#'   (default) plots all endogenous variables (capped at 12 for readability).
#' @param meta        A \code{\link{diag_meta}} object for plot provenance.
#' @return dynhr_diagnostic list with \code{$result$irf_array}
#'   (n_subsample ?-- n_periods ?-- n_vars ?-- n_shocks array) and per-shock plots.
#' @noRd
diag_bayesian_irf <- function(model,
                               compiled,
                               draws,
                               irf_periods  = 20L,
                               n_subsample  = 400L,
                               ci_bands     = c(0.10, 0.90),
                               inner_bands  = c(0.16, 0.84),
                               shock_sel    = NULL,
                               var_sel      = NULL,
                               meta         = NULL) {

    draws    <- as.matrix(draws)
    n_total  <- nrow(draws)
    n_use    <- min(n_subsample, n_total)
    draw_idx <- if (n_use < n_total) sample(n_total, n_use) else seq_len(n_total)

    par_names   <- colnames(draws)
    model_pars  <- names(model$param_values)
    endo_names  <- model$var_names
    exo_names   <- model$varexo_names

    if (is.null(shock_sel)) shock_sel <- exo_names
    if (is.null(var_sel))   var_sel   <- head(endo_names, 12L)

    shock_sel <- intersect(shock_sel, exo_names)
    var_sel   <- intersect(var_sel,   endo_names)

    n_shk  <- length(shock_sel)
    n_vars <- length(var_sel)
    n_T    <- as.integer(irf_periods)

    # Storage: [draw, period, variable, shock]
    irf_array <- array(
      NA_real_,
      dim      = c(n_use, n_T, n_vars, n_shk),
      dimnames = list(NULL,
                      paste0("h", seq_len(n_T)),
                      var_sel,
                      shock_sel)
    )

    n_ok    <- 0L
    n_fail  <- 0L
    m_work  <- model   # working copy

    for (ki in seq_len(n_use)) {
      theta_k <- draws[draw_idx[ki], ]

      # Update model parameters
      pv <- m_work$param_values
      for (nm in par_names) {
        if (nm %in% model_pars) pv[[nm]] <- theta_k[nm]
      }
      m_work$param_values <- pv

      sim_k <- stoch_simul(m_work, verbose = FALSE)

      if (is.null(sim_k) || is.null(sim_k$irfs)) { n_fail <- n_fail + 1L; next }

      irfs_k <- sim_k$irfs   # named list by shock name

      for (si in seq_along(shock_sel)) {
        sh <- shock_sel[si]
        if (!sh %in% names(irfs_k)) next
        mat_k <- irfs_k[[sh]]   # n_periods ?-- n_endo (or n_endo ?-- n_periods)
        if (ncol(mat_k) != length(endo_names)) mat_k <- t(mat_k)
        if (nrow(mat_k) < n_T || ncol(mat_k) != length(endo_names)) next
        colnames(mat_k) <- endo_names
        for (vi in seq_along(var_sel)) {
          v <- var_sel[vi]
          if (v %in% colnames(mat_k))
            irf_array[ki, , vi, si] <- mat_k[seq_len(n_T), v]
        }
      }
      n_ok <- n_ok + 1L
    }

    m_work$param_values <- model$param_values   # restore

    if (n_ok == 0L) {
      return(.make_result(
        pass    = NA,
        summary = sprintf(
          "Bayesian IRF: all %d draws failed to solve. Check model stability.", n_use)
      ))
    }

    # ---- Compute quantile bands -------------------------------------------
    probs_all <- sort(unique(c(ci_bands, inner_bands, 0.50)))

    quant_list <- lapply(seq_along(shock_sel), function(si) {
      lapply(seq_along(var_sel), function(vi) {
        mat_si_vi <- irf_array[, , vi, si]   # n_use ?-- n_T
        finite_rows <- rowSums(is.finite(mat_si_vi)) == n_T
        if (sum(finite_rows) < 3L) return(NULL)
        apply(mat_si_vi[finite_rows, , drop = FALSE], 2,
              quantile, probs = probs_all, na.rm = TRUE)
      })
    })

    # ---- Plots -------------------------------------------------------------
    plots <- list()

    if (requireNamespace("ggplot2", quietly = TRUE)) {

      for (si in seq_along(shock_sel)) {
        sh <- shock_sel[si]
        plot_rows <- list()

        for (vi in seq_along(var_sel)) {
          v  <- var_sel[vi]
          qm <- quant_list[[si]][[vi]]
          if (is.null(qm)) next

          horizons <- seq_len(n_T)
          median_v <- qm["50%", ]
          lo_out   <- qm[sprintf("%.0f%%", ci_bands[1]    * 100), ]
          hi_out   <- qm[sprintf("%.0f%%", ci_bands[2]    * 100), ]
          lo_in    <- qm[sprintf("%.0f%%", inner_bands[1] * 100), ]
          hi_in    <- qm[sprintf("%.0f%%", inner_bands[2] * 100), ]

          plot_rows[[v]] <- data.frame(
            horizon  = rep(horizons, 5),
            value    = c(median_v, lo_out, hi_out, lo_in, hi_in),
            band     = rep(c("median","outer_lo","outer_hi","inner_lo","inner_hi"),
                          each = n_T),
            variable = v,
            stringsAsFactors = FALSE
          )
        }

        if (length(plot_rows) == 0L) next
        df_sh <- do.call(rbind, plot_rows)

        # Pivot to wide for ribbon geoms
        df_wide <- merge(
          df_sh[df_sh$band == "median",   c("horizon","value","variable")],
          df_sh[df_sh$band == "outer_lo", c("horizon","value","variable")],
          by = c("horizon","variable"), suffixes = c("_med","_lo_out")
        )
        df_wide <- merge(df_wide,
          df_sh[df_sh$band == "outer_hi", c("horizon","value","variable")],
          by = c("horizon","variable"))
        names(df_wide)[names(df_wide) == "value"] <- "hi_out"
        df_wide <- merge(df_wide,
          df_sh[df_sh$band == "inner_lo", c("horizon","value","variable")],
          by = c("horizon","variable"))
        names(df_wide)[names(df_wide) == "value"] <- "lo_in"
        df_wide <- merge(df_wide,
          df_sh[df_sh$band == "inner_hi", c("horizon","value","variable")],
          by = c("horizon","variable"))
        names(df_wide)[names(df_wide) == "value"] <- "hi_in"

        p <- ggplot2::ggplot(df_wide, ggplot2::aes(x = horizon)) +
          # Outer band (10-90 or user-defined)
          ggplot2::geom_ribbon(ggplot2::aes(ymin = value_lo_out, ymax = hi_out),
                               fill  = dynhr_colours$light_blue,
                               alpha = 0.30) +
          # Inner band (16-84)
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lo_in, ymax = hi_in),
                               fill  = dynhr_colours$mid_blue,
                               alpha = 0.45) +
          # Median
          ggplot2::geom_line(ggplot2::aes(y = value_med),
                             colour    = dynhr_colours$dark_blue,
                             linewidth = 0.7) +
          ggplot2::geom_hline(yintercept = 0,
                              colour    = dynhr_colours$grey,
                              linewidth = 0.3) +
          ggplot2::facet_wrap(~ variable, scales = "free_y",
                              ncol = min(3L, n_vars)) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title    = sprintf("Bayesian IRF: shock = %s", sh),
            subtitle = sprintf(
              "%d draws * outer %d-%d%% * inner %d-%d%%",
              n_ok,
              round(ci_bands[1] * 100),    round(ci_bands[2] * 100),
              round(inner_bands[1] * 100), round(inner_bands[2] * 100)),
            x = "Horizon (quarters)", y = "Response"
          )
        p <- .apply_meta(p, meta)
        plots[[paste0("bayesian_irf_", sh)]] <- p
      }
    }

    pass <- n_fail / n_use < 0.20   # fewer than 20% of draws failed

    .make_result(
      result  = list(irf_array  = irf_array,
                     n_ok       = n_ok,
                     n_fail     = n_fail,
                     shock_sel  = shock_sel,
                     var_sel    = var_sel,
                     probs      = probs_all),
      pass    = pass,
      plots   = plots,
      summary = sprintf(
        "Bayesian IRF: %d/%d draws solved (%d shocks x %d vars x %d horizons). %s",
        n_ok, n_use, n_shk, n_vars, n_T,
        ifelse(pass, "PASS.", sprintf("FAIL: %d draws failed.", n_fail))
      )
    )
}
