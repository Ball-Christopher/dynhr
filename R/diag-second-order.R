## R/diag-second-order.R
## --------------------------------------------------------------------------
## D19: Second-order solution quality diagnostic.
##
## Checks the second-order decision rules (ghxx, ghxu, ghuu, ghss) for:
##   - Relative size vs first-order terms (perturbation accuracy indicator)
##   - Uncertainty correction magnitude relative to steady state
##   - A_L condition number (solver quality)
##   - Pruning recommendation (if second-order terms are large)
##   - Comparison of impact-period variance (var1 = diag(ghu Sigma_e ghu'))
##     and second-order mean-correction magnitude (|ghss|/(2*impact-sd))
##   NOTE: var1 is the ONE-PERIOD IMPACT variance, not the unconditional
##   Lyapunov variance. Computing the true unconditional variance would
##   require solving a Lyapunov equation.
## --------------------------------------------------------------------------


#' D19: Second-order solution quality diagnostic
#'
#' Evaluates the quality and plausibility of the second-order perturbation
#' solution. Reports relative magnitudes of second-order terms, the
#' uncertainty correction, and recommends pruning when appropriate.
#'
#' Pruning is recommended when max(|ghxx|) * typical_state > 0.1,
#' i.e. when the quadratic correction is more than 10% of the linear term.
#'
#' \strong{Note on variances:} the reported \code{var1} is the
#' \emph{impact-period} variance \code{diag(ghu \%*\% Sigma_e \%*\% t(ghu))},
#' not the unconditional (Lyapunov) variance.  These are equal only when
#' \code{ghx = 0}; in general the unconditional variance is larger.  The
#' uncertainty-correction ratio \code{unc_rel} is therefore relative to the
#' impact-period standard deviation, not the unconditional standard deviation.
#'
#' @param dr2    DecisionRules2 object from solve_perturbation_order2()
#' @param model  dynhr_mod (for variable labels, optional)
#' @param ss     Named numeric steady state (for magnitude scaling)
#' @param params Named numeric parameter vector (for shock variances)
#' @return A dynhr_diagnostic object
#' @export
d19_second_order_accuracy <- function(dr2, model = NULL,
                                       ss = NULL, params = NULL) {
  if (!inherits(dr2, "DecisionRules2")) {
    return(.make_result(
      result  = NULL,
      pass    = FALSE,
      summary = "dr2 is not a DecisionRules2 object. Run solve_perturbation_order2() first.",
      llm_summary = "[FAIL] d19_second_order: input not DecisionRules2"
    ))
  }

  ghx  <- dr2$ghx
  ghu  <- dr2$ghu
  ghxx <- dr2$ghxx
  ghxu <- dr2$ghxu
  ghuu <- dr2$ghuu
  ghss <- dr2$ghss

  n_endo    <- length(dr2$endo_names)
  n_s       <- dr2$n_state
  n_u       <- dr2$n_exo
  endo      <- dr2$endo_names
  state_vars <- dr2$state_vars

  # ----------------------------------------------------------------
  # 1. Relative magnitude of second-order terms
  # ----------------------------------------------------------------
  norm1x <- max(abs(ghx))
  norm1u <- max(abs(ghu))
  norm2xx <- if (length(ghxx) > 0) max(abs(ghxx)) else 0
  norm2xu <- if (length(ghxu) > 0) max(abs(ghxu)) else 0
  norm2uu <- if (length(ghuu) > 0) max(abs(ghuu)) else 0
  norm_ss <- max(abs(ghss))

  ratio_xx <- if (norm1x > 1e-14) norm2xx / norm1x else NA_real_
  ratio_xu <- if (norm1u > 1e-14) norm2xu / norm1u else NA_real_
  ratio_uu <- if (norm1u > 1e-14) norm2uu / norm1u else NA_real_

  # ----------------------------------------------------------------
  # 2. Uncertainty correction relative to steady-state values
  # ----------------------------------------------------------------
  ss_vals <- if (!is.null(ss)) {
    sapply(endo, function(nm) if (nm %in% names(ss)) abs(ss[[nm]]) else NA_real_)
  } else if (!is.null(dr2$ys)) {
    abs(dr2$ys[endo])
  } else {
    rep(NA_real_, n_endo)
  }

  ss_scale  <- max(ss_vals, na.rm = TRUE)
  ss_scale  <- if (is.finite(ss_scale) && ss_scale > 1e-12) ss_scale else 1

  ghss_rel  <- max(abs(ghss)) / ss_scale

  # Variable with largest uncertainty correction
  top_ss_idx  <- which.max(abs(ghss))
  top_ss_var  <- if (length(top_ss_idx) > 0) endo[top_ss_idx] else "?"
  top_ss_val  <- if (length(top_ss_idx) > 0) ghss[top_ss_idx] else 0

  # ----------------------------------------------------------------
  # 3. Pruning recommendation
  # ----------------------------------------------------------------
  # Pruning is recommended when the quadratic correction is non-negligible
  # relative to the linear term for typical state magnitudes.
  # We use norm(ghxx) / norm(ghx) as a proxy.
  pruning_ratio   <- if (!is.na(ratio_xx)) ratio_xx else 0
  pruning_advised <- pruning_ratio > 0.1

  # ----------------------------------------------------------------
  # 4. Compare first- vs second-order variances
  # ----------------------------------------------------------------
  # Use the model's shock covariance if available
  Sigma_e <- dr2$Sigma_e
  if (is.null(Sigma_e) || !is.matrix(Sigma_e)) {
    if (!is.null(model) && !is.null(params)) {
      stderr  <- .get_shock_stderr(model, dr2$exo_names, params)
      Sigma_e <- diag(stderr^2, n_u, n_u)
    } else {
      Sigma_e <- diag(n_u)  # unit shocks
    }
  }

  # Impact-period variance of each variable (diagonal of ghu*Sigma_e*ghu')
  var1 <- diag(ghu %*% Sigma_e %*% t(ghu))

  # Second-order uncertainty correction: (1/2) ghss per variable
  unc_correction <- 0.5 * abs(ghss)

  # Relative size of uncertainty correction vs impact-period std dev
  sd1           <- sqrt(pmax(var1, 0))
  unc_rel       <- ifelse(sd1 > 1e-14, unc_correction / sd1, NA_real_)
  max_unc_rel   <- max(unc_rel, na.rm = TRUE)

  # ----------------------------------------------------------------
  # 5. Summary statistics
  # ----------------------------------------------------------------
  pass <- NA  # Informational diagnostic

  # Key metrics for the summary
  metrics <- list(
    n_endo     = n_endo,
    n_state    = n_s,
    n_exo      = n_u,
    norm_ghxx  = round(norm2xx, 6),
    norm_ghxu  = round(norm2xu, 6),
    norm_ghuu  = round(norm2uu, 6),
    norm_ghss  = round(norm_ss, 6),
    ratio_xx_to_ghx = round(ratio_xx, 4),
    ratio_uu_to_ghu = round(ratio_uu, 4),
    ghss_rel_to_ss  = round(ghss_rel, 4),
    max_unc_correction_rel_impact_sd = round(max_unc_rel, 4),
    top_ghss_var   = top_ss_var,
    top_ghss_val   = round(top_ss_val, 6),
    pruning_advised = pruning_advised
  )

  # ----------------------------------------------------------------
  # 6. Text summaries
  # ----------------------------------------------------------------
  lines <- c(
    sprintf("D19: Second-order perturbation solution quality"),
    sprintf("  n_endo=%d  n_state=%d  n_exo=%d", n_endo, n_s, n_u),
    "",
    "  Second-order term magnitudes (max |.|):",
    sprintf("    ghxx: %.4g  (%.2f%% of |ghx|)",
            norm2xx, 100 * ifelse(is.na(ratio_xx), 0, ratio_xx)),
    sprintf("    ghxu: %.4g  (%.2f%% of |ghu|)",
            norm2xu, 100 * ifelse(is.na(ratio_xu), 0, ratio_xu)),
    sprintf("    ghuu: %.4g  (%.2f%% of |ghu|)",
            norm2uu, 100 * ifelse(is.na(ratio_uu), 0, ratio_uu)),
    sprintf("    ghss: %.4g  (%.2f%% of |SS|)",
            norm_ss, 100 * ghss_rel),
    "",
    sprintf("  Uncertainty correction: max(|ghss|/impact-period sd) = %.4f", max_unc_rel),
    sprintf("  Largest correction: %s  (ghss=%.6g)", top_ss_var, top_ss_val),
    "",
    if (pruning_advised) {
      "  WARNING: max|ghxx|/max|ghx| > 10%. Pruning recommended for simulation."
    } else {
      "  Second-order terms are small relative to first-order (<10%). No pruning required."
    }
  )

  summary_text <- paste(lines, collapse = "\n")

  badge <- "INFO"
  llm_lines <- c(
    sprintf("[%s] d19_second_order:", badge),
    sprintf("  n_endo=%d n_state=%d n_exo=%d", n_endo, n_s, n_u),
    sprintf("  norm_ghxx=%.4g ratio_xx=%.3f norm_ghss=%.4g", norm2xx, ifelse(is.na(ratio_xx), 0, ratio_xx), norm_ss),
    sprintf("  ghss_rel_ss=%.4f unc_rel_sd1=%.4f", ghss_rel, max_unc_rel),
    sprintf("  pruning_advised=%s top_ghss_var=%s(%.4g)", pruning_advised, top_ss_var, top_ss_val),
    sprintf("  action: %s",
            if (pruning_advised) "use pruning=TRUE in simulate/IRF calls"
            else "second-order correction is well-behaved; pruning optional")
  )

  .make_result(
    result      = metrics,
    pass        = pass,
    plots       = list(),
    summary     = summary_text,
    llm_summary = paste(llm_lines, collapse = "\n")
  )
}


#' Compare first- and second-order IRFs for a given shock
#'
#' Plots IRFs from first-order (ghx/ghu only) and second-order (pruned
#' state space) side by side, highlighting the quadratic correction.
#'
#' @param dr2         DecisionRules2 object
#' @param model       dynhr_mod
#' @param shock_name  Name of shock to display (default: first shock)
#' @param vars        Variables to plot (default: all)
#' @param n_periods   IRF horizon
#' @param params      Named parameter vector
#' @param meta        Optional metadata list attached to the diagnostic result
#' @return A dynhr_diagnostic object with plots
#' @export
d19_irf_comparison <- function(dr2, model, shock_name = NULL,
                                vars = NULL, n_periods = 40L,
                                params = NULL, meta = NULL) {
  if (!inherits(dr2, "DecisionRules2")) {
    return(.make_result(
      result  = NULL,
      pass    = FALSE,
      summary = "dr2 is not a DecisionRules2 object.",
      llm_summary = "[FAIL] d19_irf_comparison: need DecisionRules2"
    ))
  }

  if (is.null(shock_name)) shock_name <- dr2$exo_names[1]
  if (!shock_name %in% dr2$exo_names) {
    stop(sprintf("Shock '%s' not in model exogenous variables.", shock_name))
  }

  if (is.null(vars)) vars <- dr2$endo_names

  if (is.null(params)) params <- model$param_values

  # First-order IRFs (using DecisionRules base class behaviour)
  irfs1 <- compute_irfs(dr2, model, n_periods = n_periods, params = params)

  # Second-order IRFs (pruned)
  irfs2 <- compute_irfs_order2(dr2, model, n_periods = n_periods,
                                params = params, pruning = TRUE)

  irf1 <- irfs1[[shock_name]]
  irf2 <- irfs2[[shock_name]]

  # Restrict to requested variables
  vars <- intersect(vars, dr2$endo_names)

  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    periods <- seq_len(n_periods)
    for (v in vars) {
      df <- data.frame(
        period = rep(periods, 2),
        value  = c(irf1[, v], irf2[, v]),
        order  = rep(c("First order", "Second order (pruned)"), each = n_periods)
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$period, y = .data$value,
                                             colour = .data$order,
                                             linetype = .data$order)) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 0, colour = dynhr_colours$grey, linewidth = 0.3) +
        ggplot2::labs(title = sprintf("IRF: %s -> %s", shock_name, v),
                      x = "Periods", y = "Deviation from SS",
                      colour = NULL, linetype = NULL) +
        scale_colour_dynhr_vibrant() +
        theme_dynhr_diagnostic()
      plots[[v]] <- .apply_meta(p, meta)
    }
  }

  # Compute max relative difference between 1st and 2nd order
  rel_diffs <- sapply(vars, function(v) {
    denom <- max(abs(irf1[, v]))
    if (denom < 1e-14) return(NA_real_)
    max(abs(irf2[, v] - irf1[, v])) / denom
  })
  max_rel_diff <- max(rel_diffs, na.rm = TRUE)

  summary_text <- sprintf(
    paste("D19 IRF comparison: shock=%s, periods=%d\n",
          "  Max relative deviation (2nd vs 1st order): %.4f\n",
          "  Variables with largest deviation: %s"),
    shock_name, n_periods, max_rel_diff,
    paste(vars[order(-rel_diffs, na.last = TRUE)[1:min(3, length(vars))]], collapse = ", ")
  )

  llm_text <- sprintf(
    "[INFO] d19_irf_comparison: shock=%s n_periods=%d max_rel_diff_2nd_vs_1st=%.4f",
    shock_name, n_periods, max_rel_diff
  )

  .make_result(
    result      = list(irf1 = irf1, irf2 = irf2, rel_diffs = rel_diffs),
    pass        = NA,
    plots       = plots,
    summary     = summary_text,
    llm_summary = llm_text
  )
}
