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
#' and 16-84%).
#'
#' Each IRF is the response to a ONE-STANDARD-DEVIATION shock at that draw:
#' estimated shock stds (\code{stderr <shock>} columns, named after the
#' shock) scale the unit-shock \code{ghu} per shock, so shock-std uncertainty
#' is part of the bands.
#'
#' Computational note: each draw requires a steady-state + first-order solve
#' (the likelihood's own pipeline, reusing \code{compiled}).  Use
#' \code{n_subsample} to trade resolution for speed.
#'
#' @param model       dynhr_mod from \code{parse_mod()}.
#' @param compiled    dynhr_compiled from \code{compile_model()}.
#' @param draws       Posterior draws matrix (n_draws x n_params); column
#'   names must match estimated parameter names (structural parameters or
#'   shock names for estimated shock stds).
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
    if (is.null(par_names))
      .dynhr_abort("diag_bayesian_irf: `draws` needs column names (the estimated parameter names).")
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

    ## Each draw is re-solved with the likelihood's own pipeline
    ## (.d8_irfs_at_theta). Until 0.9.4 this called stoch_simul() on a model
    ## copy with only STRUCTURAL params overwritten: estimated `stderr <shock>`
    ## draws (shock-named columns) were dropped, so every band was scaled by
    ## the calibrated shock stds, and `compiled` was ignored (a recompile per
    ## draw).
    if (is.null(compiled$lead_lag_incidence) &&
        !is.null(compiled$model$lead_lag_incidence))
      compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
    sys_cache <- cache_system_structure(compiled)
    state <- new.env(parent = emptyenv())
    state$ss_warm <- NULL

    for (ki in seq_len(n_use)) {
      theta_k <- draws[draw_idx[ki], ]
      names(theta_k) <- par_names

      irfs_k <- .d8_irfs_at_theta(model, compiled, theta_k, n_periods = n_T,
                                  sys_cache = sys_cache, state = state)
      if (is.null(irfs_k)) { n_fail <- n_fail + 1L; next }

      for (si in seq_along(shock_sel)) {
        mat_k <- irfs_k[[shock_sel[si]]]    # n_T x n_endo, columns named
        if (is.null(mat_k)) next
        irf_array[ki, , , si] <- mat_k[, var_sel, drop = FALSE]
      }
      n_ok <- n_ok + 1L
    }

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

        ## Quantile rows are looked up by POSITION in probs_all: quantile()'s
        ## row names keep decimals (0.025 -> "2.5%") while the old
        ## sprintf("%.0f%%") lookup gave "2%", so any non-integer-percent band
        ## errored.
        q_row <- function(p) which.min(abs(probs_all - p))
        ## A variable that does not respond to this shock (|band| at rounding
        ## level everywhere) is omitted: its panel would be pure 1e-17 noise.
        resp_scale <- max(vapply(quant_list[[si]], function(qm)
          if (is.null(qm)) 0 else max(abs(qm)), numeric(1)))
        n_flat <- 0L
        for (vi in seq_along(var_sel)) {
          v  <- var_sel[vi]
          qm <- quant_list[[si]][[vi]]
          if (is.null(qm)) next
          if (max(abs(qm)) <= 1e-10 * resp_scale) { n_flat <- n_flat + 1L; next }

          plot_rows[[v]] <- data.frame(
            horizon  = seq_len(n_T),
            variable = v,
            med      = qm[q_row(0.50), ],
            lo_out   = qm[q_row(ci_bands[1]), ],
            hi_out   = qm[q_row(ci_bands[2]), ],
            lo_in    = qm[q_row(inner_bands[1]), ],
            hi_in    = qm[q_row(inner_bands[2]), ],
            stringsAsFactors = FALSE
          )
        }

        if (length(plot_rows) == 0L) next
        df_wide  <- do.call(rbind, plot_rows)
        n_panels <- length(plot_rows)

        p <- ggplot2::ggplot(df_wide, ggplot2::aes(x = horizon)) +
          # Outer band (10-90 or user-defined)
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lo_out, ymax = hi_out),
                               fill  = dynhr_colours$light_blue,
                               alpha = 0.35) +
          # Inner band (16-84)
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lo_in, ymax = hi_in),
                               fill  = dynhr_colours$mid_blue,
                               alpha = 0.45) +
          ggplot2::geom_hline(yintercept = 0,
                              colour    = dynhr_colours$grey,
                              linewidth = 0.3) +
          # Median
          ggplot2::geom_line(ggplot2::aes(y = med),
                             colour    = dynhr_colours$dark_blue,
                             linewidth = 0.7) +
          ggplot2::facet_wrap(~ variable, scales = "free_y",
                              ncol = min(3L, n_panels)) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title    = sprintf("Bayesian IRF: shock = %s", sh),
            subtitle = paste0(sprintf(
              "Posterior median; bands %g-%g%% (light), %g-%g%% (dark); %d draws",
              ci_bands[1] * 100,    ci_bands[2] * 100,
              inner_bands[1] * 100, inner_bands[2] * 100, n_ok),
              if (n_flat > 0L)
                sprintf("; %d non-responding variable(s) omitted", n_flat)
              else ""),
            x = "Horizon (1 = impact)",
            y = "Response to a 1 s.d. shock"
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
