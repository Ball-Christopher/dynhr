## R/diag-deep-d34-invariance.R
## --------------------------------------------------------------------------
## D34. Policy-partitioned invariance test (the operational Lucas critique).
##
## The Lucas critique says a model is only structural if its *private-sector*
## deep parameters (preferences, technology, frictions) stay put when the
## *policy rule* changes. D16 (subsample stability) tests every estimated
## parameter symmetrically; D34 uses the @dynhr:deep partition to ask the
## sharper, directional question:
##
##   - PRIVATE deep block: must be invariant across the policy regimes
##     (overlapping credible intervals). A non-invariant private parameter is a
##     literal Lucas-critique violation -> the "deep" parameter is absorbing the
##     policy change and is not deep.
##   - POLICY block: is *expected* to move when the regime changes. A policy
##     block that does NOT move is itself suspicious (mis-dated break, or the
##     policy did not actually change) and is reported as a caveat, not a pass.
##   - AUXILIARY (shock) block: reported for context only -- shock volatilities
##     moving across regimes is stochastic-volatility, not a Lucas violation.
##
## A complementary super-exogeneity test (Engle, Hendry & Richard 1983) asks
## whether policy innovations load on the private structural residuals (they
## should not, under invariance).
##
## EFFICIENCY: D34 consumes the SAME inputs as D16 -- the full-sample `draws`
## and the named list of per-regime draws `results_sub`. The estimation
## pipeline therefore re-estimates each regime once and both diagnostics read
## those draws; D34 adds the partition, the directional interpretation and the
## super-exogeneity test on top.
##
## Exposes `$result$passport_axis` (named logical over PRIVATE params,
## TRUE = invariant) for the Deep-Parameter Passport's "invariant" column.
##
## References:
##   Lucas, R. E. (1976). Econometric policy evaluation: a critique.
##   Engle, R. F., Hendry, D. F., & Richard, J.-F. (1983). Exogeneity.
##     Econometrica, 51(2), 277-304.
##   Fernandez-Villaverde, J., & Rubio-Ramirez, J. F. (2008). How structural
##     are structural parameters? NBER Macroeconomics Annual 2007, 22, 83-137.
##   Inoue, A., & Rossi, B. (2011). Identifying the sources of instabilities in
##     macroeconomic fluctuations. Review of Economics and Statistics, 93(4).
## --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#' D34. Policy-partitioned invariance test (operational Lucas critique)
#'
#' @param model       Parsed model (for the @dynhr:deep policy/private
#'   partition).  Optional if \code{deep_spec} is supplied.
#' @param deep_spec   Optional \code{\link{build_deep_spec}} (built from
#'   \code{model} and the draw column names otherwise).
#' @param draws       Full-sample posterior draws (\eqn{n \times p}) with named
#'   columns -- the same object D16 receives.
#' @param results_sub Named list of per-regime posterior draw matrices, e.g.
#'   \code{list("Pre-1990" = d1, "Post-1990" = d2)} -- the same object D16
#'   receives. Two or more regimes give an across-regime test; one regime falls
#'   back to a regime-vs-full comparison.
#' @param param_names Optional parameter names (defaults to \code{colnames(draws)}).
#' @param ci_level    Credible-interval level for the overlap test (default 0.90).
#' @param super_exog  Optional list for the Engle-Hendry-Richard test:
#'   \code{list(shocks = <T x k matrix with named columns, or named list of
#'   per-regime matrices>, policy = <character names of policy innovations>,
#'   regime = <length-T factor; required when shocks is a single matrix>)}.
#' @param meta        Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic}; \code{$result$passport_axis} is a named
#'   logical over the private block (TRUE = invariant) for the Passport.
#' @noRd
d34_policy_invariance <- function(model       = NULL,
                                  deep_spec   = NULL,
                                  draws       = NULL,
                                  results_sub = NULL,
                                  param_names = NULL,
                                  ci_level    = 0.90,
                                  super_exog  = NULL,
                                  meta        = NULL) {
  tryCatch({
    if (is.null(draws))
      return(.make_result(pass = NA,
        summary = "D34 policy invariance: no full-sample draws supplied."))
    draws <- as.matrix(draws)
    if (is.null(param_names))
      param_names <- colnames(draws) %||% paste0("theta_", seq_len(ncol(draws)))
    colnames(draws) <- param_names

    if (is.null(results_sub) || length(results_sub) == 0)
      return(.make_result(pass = NA,
        summary = paste("D34 policy invariance: no per-regime draws (results_sub);",
                        "supply >=1 regime split (shared with D16).")))
    if (is.null(names(results_sub)))
      names(results_sub) <- paste0("Regime_", seq_along(results_sub))

    if (is.null(deep_spec))
      deep_spec <- build_deep_spec(model = model, param_names = param_names)

    # Partition the *estimated* parameters.
    ds  <- deep_spec[match(param_names, deep_spec$param), , drop = FALSE]
    ds$param <- param_names
    ds$class[is.na(ds$class)] <- "unknown"
    block <- ifelse(ds$class == "policy", "policy",
             ifelse(ds$is_deep %in% TRUE, "private", "auxiliary"))
    block[is.na(block)] <- "auxiliary"
    names(block) <- param_names
    private_p <- param_names[block == "private"]
    policy_p  <- param_names[block == "policy"]

    # Samples for the across-regime test: the regimes themselves if >=2, else
    # the single regime against the full sample.
    if (length(results_sub) >= 2) {
      test_draws <- results_sub
    } else {
      test_draws <- c(list("Full sample" = draws), results_sub)
    }
    summ_test <- lapply(names(test_draws), function(nm)
      .summarise_draws_ci(test_draws[[nm]], nm, param_names, ci_level))
    names(summ_test) <- names(test_draws)

    # Per-parameter invariance: do all pairs of test-sample CIs overlap?
    n_s <- length(summ_test)
    invariant <- vapply(seq_along(param_names), function(j) {
      lo <- vapply(summ_test, function(s) s$lo[j], numeric(1))
      hi <- vapply(summ_test, function(s) s$hi[j], numeric(1))
      ok <- TRUE
      for (a in seq_len(n_s - 1)) for (b in (a + 1):n_s)
        if (.ci_disjoint(lo[a], hi[a], lo[b], hi[b])) ok <- FALSE
      ok
    }, logical(1))
    names(invariant) <- param_names

    private_unstable <- private_p[!invariant[private_p]]   # Lucas violations
    policy_moved     <- policy_p[!invariant[policy_p]]     # expected
    policy_static    <- policy_p[invariant[policy_p]]      # suspicious

    # --- optional super-exogeneity test ---
    se <- NULL
    if (!is.null(super_exog)) {
      se <- tryCatch(.d34_super_exogeneity(super_exog, alpha = 1 - ci_level),
                     error = function(e)
                       list(error = conditionMessage(e), rejected = NA))
    }
    se_rejected <- if (!is.null(se)) isTRUE(se$rejected) else FALSE

    # Verdict: the Lucas critique holds iff no private parameter is
    # non-invariant and super-exogeneity (if tested) is not rejected.
    n_private <- length(private_p)
    pass <- if (n_private == 0 && is.null(se)) NA
            else (length(private_unstable) == 0 && !se_rejected)

    passport_axis <- if (n_private > 0)
      stats::setNames(invariant[private_p], private_p) else logical(0)

    # --- plot ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_df <- do.call(rbind, c(
        list(.summarise_draws_ci(draws, "Full sample", param_names, ci_level)),
        lapply(names(results_sub), function(nm)
          .summarise_draws_ci(results_sub[[nm]], nm, param_names, ci_level))))
      plot_df$block  <- block[plot_df$param]
      plot_df$status <- invariant[plot_df$param]
      plots$forest <- tryCatch(.plot_d34_forest(plot_df, ci_level, meta),
                               error = function(e) NULL)
    }

    summary_txt <- sprintf(
      "D34 Policy invariance: %d private / %d policy / %d aux params, %d regime(s). %s%s",
      n_private, length(policy_p), sum(block == "auxiliary"), length(results_sub),
      if (is.na(pass)) "No private deep parameters estimated (nothing to test)."
      else if (isTRUE(pass)) "PASS -- private block invariant across regimes."
      else sprintf("FAIL -- Lucas violation: %s",
                   paste(c(private_unstable,
                           if (se_rejected) "super-exogeneity rejected"),
                         collapse = ", ")),
      if (length(policy_static) > 0)
        sprintf(" (caveat: policy params %s did not move -- check the break date)",
                paste(policy_static, collapse = ", ")) else "")

    badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
    llm <- paste(c(
      sprintf("D34 | Policy Invariance (Lucas) | %s", badge),
      sprintf("  private=%d policy=%d aux=%d regimes=%d ci=%.0f%%",
              n_private, length(policy_p), sum(block == "auxiliary"),
              length(results_sub), 100 * ci_level),
      if (length(private_unstable))
        sprintf("  lucas_violation(private non-invariant): %s",
                paste(private_unstable, collapse = ", ")),
      if (length(policy_moved))
        sprintf("  policy_responded(expected): %s",
                paste(policy_moved, collapse = ", ")),
      if (length(policy_static))
        sprintf("  policy_static(caveat): %s", paste(policy_static, collapse = ", ")),
      if (!is.null(se) && is.null(se$error))
        sprintf("  super_exogeneity: %s (%d/%d private eqns load on policy)",
                if (isTRUE(se$rejected)) "REJECTED" else "ok",
                sum(is.finite(se$table$p) & se$table$p < (1 - ci_level), na.rm = TRUE),
                nrow(se$table)),
      sprintf("  action: %s",
              if (isTRUE(pass))
                "Private deep parameters are invariant to the policy shift; the model survives the Lucas critique on this split."
              else if (is.na(pass))
                "Estimate at least one private deep parameter (not just policy/shock params) to run the Lucas test."
              else sprintf("%s move with the policy regime -- they are absorbing the policy change, not structural. Re-specify or treat as regime-specific.",
                           paste(utils::head(c(private_unstable,
                                 if (se_rejected) "(super-exog)"), 3), collapse = ", ")))
    ), collapse = "\n")

    .make_result(
      result = list(invariant = invariant, block = block,
                    private_unstable = private_unstable,
                    policy_moved = policy_moved, policy_static = policy_static,
                    super_exogeneity = se, passport_axis = passport_axis,
                    summaries = summ_test),
      pass    = pass,
      plots   = plots,
      summary = summary_txt,
      llm_summary = llm)
  }, error = function(e) {
    .make_result(pass = NA, errored = TRUE,
                 summary = paste("D34 policy invariance: ERROR --",
                                 conditionMessage(e)))
  })
}


# Engle-Hendry-Richard super-exogeneity: do policy innovations enter the
# private structural residuals? For each private shock e_p, F-test the full
# model  e_p ~ regime * (policy innovations)  against the reduced  e_p ~ regime.
# Joint significance of the policy terms rejects super-exogeneity / invariance.
.d34_super_exogeneity <- function(super_exog, alpha = 0.10) {
  shocks <- super_exog$shocks
  if (is.list(shocks) && !is.data.frame(shocks) && !is.matrix(shocks)) {
    reg_lab <- rep(names(shocks), vapply(shocks, NROW, integer(1)))
    shocks  <- do.call(rbind, lapply(shocks, as.matrix))
    regime  <- factor(reg_lab)
  } else {
    shocks <- as.matrix(shocks)
    regime <- super_exog$regime
    if (is.null(regime)) stop("super_exog$regime is required when shocks is a single matrix")
    regime <- as.factor(regime)
  }
  all_names <- colnames(shocks)
  if (is.null(all_names)) stop("super_exog$shocks needs column names")
  pol_names <- intersect(super_exog$policy, all_names)
  if (length(pol_names) == 0) stop("no policy-shock columns matched super_exog$policy")
  priv_names <- setdiff(all_names, pol_names)
  if (length(priv_names) == 0) stop("no private-shock columns remain")
  if (nlevels(regime) < 2) stop("super_exog$regime needs >= 2 levels")

  Pn <- make.names(pol_names)
  rows <- lapply(priv_names, function(pn) {
    dat <- data.frame(y = shocks[, pn], regime = regime)
    for (k in seq_along(pol_names)) dat[[Pn[k]]] <- shocks[, pol_names[k]]
    f_full <- stats::as.formula(
      paste0("y ~ regime * (", paste(Pn, collapse = " + "), ")"))
    m_full <- tryCatch(stats::lm(f_full, dat), error = function(e) NULL)
    m_red  <- tryCatch(stats::lm(y ~ regime, dat), error = function(e) NULL)
    if (is.null(m_full) || is.null(m_red))
      return(data.frame(shock = pn, F = NA_real_, p = NA_real_))
    an <- tryCatch(stats::anova(m_red, m_full), error = function(e) NULL)
    if (is.null(an) || nrow(an) < 2)
      return(data.frame(shock = pn, F = NA_real_, p = NA_real_))
    data.frame(shock = pn, F = an$F[2], p = an$`Pr(>F)`[2])
  })
  tab <- do.call(rbind, rows)
  # Bonferroni across private equations.
  rejected <- any(is.finite(tab$p) & tab$p < alpha / max(1L, nrow(tab)))
  list(table = tab, rejected = rejected, alpha = alpha,
       policy = pol_names, private = priv_names)
}


# Partitioned forest: facet by block (Private / Policy / Auxiliary) with the
# private panel front-and-centre; regimes dodged within each parameter row.
.plot_d34_forest <- function(plot_df, ci_level, meta) {
  blk_levels <- c("private", "policy", "auxiliary")
  blk_labels <- c(private = "Private (must be invariant)",
                  policy  = "Policy (expected to move)",
                  auxiliary = "Auxiliary / shocks")
  plot_df$block <- factor(blk_labels[plot_df$block],
                          levels = unname(blk_labels[blk_levels]))
  ord <- unique(plot_df$param[order(match(plot_df$block, levels(plot_df$block)))])
  plot_df$param <- factor(plot_df$param, levels = rev(ord))

  samples <- unique(plot_df$sample)
  pal <- if (length(samples) <= 3)
    stats::setNames(c(dynhr_colours$grey, dynhr_colours$dark_blue, dynhr_colours$orange)[seq_along(samples)], samples)
  else stats::setNames(dynhr_palette[seq_along(samples)], samples)

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = median, y = param,
                                    colour = sample)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 1.9) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi),
                            position = ggplot2::position_dodge(width = 0.6),
                            height = 0.3, linewidth = 0.45) +
    ggplot2::facet_grid(rows = ggplot2::vars(block),
                        scales = "free_y", space = "free_y",
                        # Wrap long strip labels so "must be invariant" is not truncated
                        labeller = ggplot2::label_wrap_gen(width = 20)) +
    # Limit tick density to avoid overlapping x-axis labels in narrow panels
    ggplot2::scale_x_continuous(n.breaks = 3L) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::labs(
      title    = "D34: Policy-partitioned invariance (operational Lucas critique)",
      subtitle = sprintf("%d%% CIs | private block should be invariant; policy block is expected to move",
                         round(100 * ci_level)),
      x = "Parameter value", y = NULL)
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}
