## R/diag-second-order.R
## --------------------------------------------------------------------------
## D19: Second-order solution quality diagnostic.
##
## For the second-order decision rule
##   y_t - y* = ghx x + ghu u + 0.5 ghxx (x (x) x) + ghxu (u (x) x)
##              + 0.5 ghuu (u (x) u) + 0.5 ghss,      x = x_{t-1}, u = u_t,
## reports, per endogenous variable:
##   - nonlinearity ratio  sd(quadratic term) / sd(linear term), with
##     z = (x_{t-1}, u_t) ~ N(0, diag(Sigma_x, Sigma_e)) and Sigma_x the
##     unconditional FIRST-ORDER state covariance (Lyapunov). Both moments
##     are exact for Gaussian z: Var(G z) = G S G',
##     Var(z' A z) = 2 tr(A S A S).  Dimensionless and scale-aware (the old
##     max|ghxx|/max|ghx| ratio ignored the shock size and mixed units).
##   - risk (uncertainty) correction 0.5 * ghss, relative to the unconditional
##     first-order sd and to |steady state|.
## Sigma_e is dr2$Sigma_e -- the covariance ghss was solved with -- so the
## risk correction and the sd it is compared against are on the same scale.
## There is NO identity-covariance default (ghu excludes Sigma_e).
## --------------------------------------------------------------------------


## Shock covariance for D19: dr2$Sigma_e (what ghss was solved with), else the
## model's covariance at `params`. Returns NULL when neither is available.
.d19_sigma_e <- function(dr2, model, params) {
  n_u <- length(dr2$exo_names)
  Se  <- dr2$Sigma_e
  if (is.matrix(Se) && all(dim(Se) == c(n_u, n_u)) && all(is.finite(Se)))
    return(Se)
  if (is.null(model)) return(NULL)
  if (is.null(params)) params <- model$param_values
  .get_shock_cov(model, dr2$exo_names, params)
}


#' D19: Second-order solution quality diagnostic
#'
#' Evaluates how much the second-order perturbation solution departs from the
#' first-order one, per endogenous variable, on the scale of the model's own
#' shocks.
#'
#' Write the decision rule as \eqn{y - y^* = G z + z' A z + \tfrac12 g_{\sigma\sigma}}
#' with \eqn{z = (x_{t-1}, u_t)}, \eqn{G = [g_x, g_u]} and \eqn{A} built from
#' \code{ghxx}, \code{ghxu}, \code{ghuu}. With \eqn{z \sim N(0, S)},
#' \eqn{S = \mathrm{diag}(\Sigma_x, \Sigma_e)} and \eqn{\Sigma_x} the
#' unconditional first-order state covariance:
#' \itemize{
#'   \item \code{nl_ratio} = \eqn{\sqrt{2\,\mathrm{tr}(ASAS)} / \sqrt{G S G'}},
#'     the standard deviation of the quadratic term relative to that of the
#'     linear term. Values above \code{threshold} flag a variable whose
#'     second-order dynamics are material (first-order IRFs/moments are then a
#'     poor guide; simulate with the pruned second-order system).
#'   \item \code{risk_corr} = \eqn{\tfrac12 g_{\sigma\sigma}}, the constant
#'     uncertainty correction, reported relative to the unconditional
#'     first-order sd (\code{risk_rel_sd}) and to \eqn{|y^*|}
#'     (\code{risk_rel_ss}, \code{NA} where the steady state is zero).
#' }
#' If the first-order state transition has a unit root, \eqn{\Sigma_x} does not
#' exist; the state block is then set to zero (impact-period moments only) and
#' the summary says so.
#'
#' The shock covariance is \code{dr2$Sigma_e} (the covariance \code{ghss} was
#' solved with); only when that is absent is it rebuilt from \code{model} and
#' \code{params}. There is no identity-covariance fallback: without either,
#' the diagnostic returns a FAIL result.
#'
#' @param dr2       DecisionRules2 object from solve_perturbation_order2()
#' @param model     dynhr_mod (only used when \code{dr2$Sigma_e} is absent)
#' @param ss        Named numeric steady state (default \code{dr2$ys})
#' @param params    Named numeric parameter vector (only used with \code{model})
#' @param threshold Nonlinearity-ratio threshold above which second-order
#'   effects are reported as material (default 0.1)
#' @param meta      Optional metadata list attached to the plot
#' @return A dynhr_diagnostic object (informational: \code{pass = NA}) whose
#'   \code{result$by_variable} holds the per-variable table and whose
#'   \code{plots$nonlinearity} shows it.
#' @keywords internal
#' @export
d19_second_order_accuracy <- function(dr2, model = NULL,
                                       ss = NULL, params = NULL,
                                       threshold = 0.1, meta = NULL) {
  if (!inherits(dr2, "DecisionRules2")) {
    return(.make_result(
      result  = NULL,
      pass    = FALSE,
      summary = "dr2 is not a DecisionRules2 object. Run solve_perturbation_order2() first.",
      llm_summary = "[FAIL] d19_second_order: input not DecisionRules2"
    ))
  }

  endo      <- dr2$endo_names
  exo       <- dr2$exo_names
  state_idx <- dr2$state_idx
  n_endo    <- length(endo)
  n_s       <- length(state_idx)
  n_u       <- length(exo)
  n_z       <- n_s + n_u

  Sigma_e <- .d19_sigma_e(dr2, model, params)
  if (is.null(Sigma_e)) {
    return(.make_result(
      result  = NULL,
      pass    = FALSE,
      summary = paste("D19: no shock covariance (dr2$Sigma_e absent and no model",
                      "supplied); refusing to assume unit shocks."),
      llm_summary = "[FAIL] d19_second_order: no Sigma_e available"
    ))
  }

  ghx  <- matrix(dr2$ghx, n_endo, n_s)
  ghu  <- matrix(dr2$ghu, n_endo, n_u)
  ghss <- as.numeric(dr2$ghss)
  ghxx <- if (length(dr2$ghxx)) matrix(dr2$ghxx, n_endo, n_s * n_s) else matrix(0, n_endo, n_s * n_s)
  ghxu <- if (length(dr2$ghxu)) matrix(dr2$ghxu, n_endo, n_s * n_u) else matrix(0, n_endo, n_s * n_u)
  ghuu <- if (length(dr2$ghuu)) matrix(dr2$ghuu, n_endo, n_u * n_u) else matrix(0, n_endo, n_u * n_u)

  # ----------------------------------------------------------------
  # 1. Unconditional first-order state covariance
  # ----------------------------------------------------------------
  stationary <- TRUE
  Sigma_x    <- matrix(0, n_s, n_s)
  if (n_s > 0L) {
    hx <- ghx[state_idx, , drop = FALSE]
    hu <- ghu[state_idx, , drop = FALSE]
    Sx <- .solve_lyapunov(hx, hu %*% Sigma_e %*% t(hu))
    if (all(is.finite(Sx))) Sigma_x <- (Sx + t(Sx)) / 2 else stationary <- FALSE
  }
  S <- matrix(0, n_z, n_z)
  if (n_s > 0L) S[seq_len(n_s), seq_len(n_s)] <- Sigma_x
  S[n_s + seq_len(n_u), n_s + seq_len(n_u)] <- Sigma_e

  # ----------------------------------------------------------------
  # 2. Per-variable linear / quadratic standard deviations
  # ----------------------------------------------------------------
  # ghxx/ghxu/ghuu columns are column-major with the FIRST index fastest:
  # ghxu column (s, u) multiplies u_u * x_s, so matrix(ghxu[i, ], n_s, n_u)
  # is Hxu[s, u] and the term is x' Hxu u (no 1/2).
  G    <- cbind(ghx, ghu)
  var1 <- pmax(rowSums((G %*% S) * G), 0)
  var2 <- numeric(n_endo)
  for (i in seq_len(n_endo)) {
    Hxx <- matrix(ghxx[i, ], n_s, n_s)
    Hxu <- matrix(ghxu[i, ], n_s, n_u)
    Huu <- matrix(ghuu[i, ], n_u, n_u)
    A   <- 0.5 * rbind(cbind((Hxx + t(Hxx)) / 2, Hxu),
                       cbind(t(Hxu), (Huu + t(Huu)) / 2))
    AS  <- A %*% S
    var2[i] <- max(2 * sum(AS * t(AS)), 0)
  }
  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  nl_ratio <- ifelse(sd1 > 1e-14, sd2 / sd1, ifelse(sd2 > 1e-14, Inf, NA_real_))

  # ----------------------------------------------------------------
  # 3. Risk correction 0.5 * ghss
  # ----------------------------------------------------------------
  if (is.list(ss) && !is.null(ss$values)) ss <- ss$values
  if (is.null(ss)) ss <- dr2$ys
  ss_vals <- if (is.null(ss)) rep(NA_real_, n_endo) else
    unname(ss[endo])
  ss_vals <- as.numeric(ss_vals)
  risk_corr   <- 0.5 * ghss
  risk_rel_sd <- ifelse(sd1 > 1e-14, abs(risk_corr) / sd1, NA_real_)
  risk_rel_ss <- ifelse(is.finite(ss_vals) & abs(ss_vals) > 1e-12,
                        abs(risk_corr) / abs(ss_vals), NA_real_)

  by_var <- data.frame(
    variable    = endo,
    sd_linear   = sd1,
    sd_quad     = sd2,
    nl_ratio    = nl_ratio,
    risk_corr   = risk_corr,
    risk_rel_sd = risk_rel_sd,
    risk_rel_ss = risk_rel_ss,
    stringsAsFactors = FALSE
  )

  .max_or_na <- function(v) if (any(is.finite(v))) max(v[is.finite(v)]) else NA_real_
  .which_max <- function(v) {
    if (!any(!is.na(v))) return(NA_character_)
    endo[which.max(replace(v, is.na(v), -Inf))]
  }
  max_nl       <- if (any(is.infinite(nl_ratio))) Inf else .max_or_na(nl_ratio)
  top_nl_var   <- .which_max(nl_ratio)
  max_risk_sd  <- .max_or_na(risk_rel_sd)
  max_risk_ss  <- .max_or_na(risk_rel_ss)
  top_risk_var <- if (n_endo > 0L) endo[which.max(abs(risk_corr))] else NA_character_
  top_risk_val <- if (n_endo > 0L) max(abs(risk_corr)) * sign(risk_corr[which.max(abs(risk_corr))]) else NA_real_
  material     <- isTRUE(max_nl > threshold)
  material_vars <- endo[!is.na(nl_ratio) & nl_ratio > threshold]

  metrics <- list(
    n_endo      = n_endo,
    n_state     = n_s,
    n_exo       = n_u,
    stationary  = stationary,
    threshold   = threshold,
    norm_ghxx   = if (length(ghxx)) max(abs(ghxx)) else 0,
    norm_ghxu   = if (length(ghxu)) max(abs(ghxu)) else 0,
    norm_ghuu   = if (length(ghuu)) max(abs(ghuu)) else 0,
    norm_ghss   = if (length(ghss)) max(abs(ghss)) else 0,
    max_nl_ratio       = max_nl,
    top_nl_var         = top_nl_var,
    max_risk_rel_sd    = max_risk_sd,
    max_risk_rel_ss    = max_risk_ss,
    top_risk_var       = top_risk_var,
    top_risk_corr      = top_risk_val,
    nonlinearity_material = material,
    material_vars      = material_vars,
    by_variable        = by_var,
    Sigma_e            = Sigma_e,
    Sigma_x            = Sigma_x
  )

  # ----------------------------------------------------------------
  # 4. Plot: per-variable nonlinearity ratio and relative risk correction
  # ----------------------------------------------------------------
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE) && n_endo > 0L) {
    lab_nl   <- "Nonlinearity: sd(quadratic) / sd(linear)"
    lab_risk <- "Risk correction: |0.5 ghss| / sd(linear)"
    ord <- endo[order(ifelse(is.finite(nl_ratio), nl_ratio, -1))]
    df <- data.frame(
      variable = factor(rep(endo, 2), levels = ord),
      measure  = factor(rep(c(lab_nl, lab_risk), each = n_endo),
                        levels = c(lab_nl, lab_risk)),
      value    = c(ifelse(is.finite(nl_ratio), nl_ratio, NA_real_), risk_rel_sd),
      flag     = c(ifelse(!is.na(nl_ratio) & nl_ratio > threshold,
                          "Above threshold", "Below threshold"),
                   rep("No threshold (risk)", n_endo)),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$value), , drop = FALSE]
    thr_df <- data.frame(measure = factor(lab_nl, levels = levels(df$measure)),
                         x = threshold)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$variable)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$value,
                                         yend = .data$variable),
                            colour = dynhr_na_colour, linewidth = 0.4) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$flag), size = 2.4) +
      ggplot2::geom_vline(data = thr_df, ggplot2::aes(xintercept = .data$x),
                          linetype = "dashed", colour = dynhr_colours$red) +
      ggplot2::facet_wrap(~ measure, scales = "free_x") +
      ggplot2::scale_colour_manual(
        values = c("Above threshold" = dynhr_colours$red,
                   "Below threshold" = dynhr_primary_colour,
                   "No threshold (risk)" = dynhr_colours$teal),
        drop = FALSE) +
      ggplot2::labs(
        title    = "D19: second-order solution quality",
        subtitle = sprintf(paste0("Unconditional first-order sd%s; dashed line = ",
                                  "threshold %.2g"),
                           if (stationary) "" else " (unit root: impact-period only)",
                           threshold),
        x = "Ratio (dimensionless)", y = NULL, colour = NULL) +
      theme_dynhr_diagnostic()
    plots$nonlinearity <- .apply_meta(p, meta)
  }

  # ----------------------------------------------------------------
  # 5. Text summaries
  # ----------------------------------------------------------------
  fmt <- function(v) if (is.na(v)) "NA" else sprintf("%.4g", v)
  lines <- c(
    "D19: Second-order perturbation solution quality",
    sprintf("  n_endo=%d  n_state=%d  n_exo=%d", n_endo, n_s, n_u),
    if (!stationary)
      "  NOTE: first-order state transition has a unit root; sds are impact-period only.",
    "",
    sprintf("  Nonlinearity sd(quadratic)/sd(linear): max = %s  (%s)",
            fmt(max_nl), top_nl_var),
    sprintf("  Risk correction 0.5*ghss: largest |.| = %s  (%s)",
            fmt(top_risk_val), top_risk_var),
    sprintf("    max |0.5 ghss| / sd(linear) = %s;  max |0.5 ghss| / |SS| = %s",
            fmt(max_risk_sd), fmt(max_risk_ss)),
    "",
    if (material) {
      sprintf(paste0("  Second-order effects are material (ratio > %.2g) for: %s.\n",
                     "  First-order IRFs/moments understate them; simulate with the ",
                     "pruned second-order system."),
              threshold, paste(material_vars, collapse = ", "))
    } else {
      sprintf("  Second-order effects are small (all ratios <= %.2g).", threshold)
    }
  )
  summary_text <- paste(lines, collapse = "\n")

  llm_lines <- c(
    "[INFO] d19_second_order:",
    sprintf("  n_endo=%d n_state=%d n_exo=%d stationary=%s",
            n_endo, n_s, n_u, stationary),
    sprintf("  max_nl_ratio=%s (%s) threshold=%.3g material=%s",
            fmt(max_nl), top_nl_var, threshold, material),
    sprintf("  risk_corr=%s (%s) rel_sd=%s rel_ss=%s",
            fmt(top_risk_val), top_risk_var, fmt(max_risk_sd), fmt(max_risk_ss)),
    sprintf("  action: %s",
            if (material) "use the pruned order-2 system for IRFs/simulation/moments"
            else "second-order effects small; first-order results are a fair guide")
  )

  .make_result(
    result      = metrics,
    pass        = NA,
    plots       = plots,
    summary     = summary_text,
    llm_summary = paste(llm_lines, collapse = "\n")
  )
}


#' Compare first- and second-order IRFs for a given shock
#'
#' Plots IRFs from the first-order rule (ghx/ghu only) and the pruned
#' second-order system for the SAME one-standard-deviation shock vector, as
#' deviations from the deterministic steady state net of the constant risk
#' correction (so the two coincide when the model is linear).
#'
#' The shock vector is column \code{k} of the lower Cholesky factor of
#' \code{dr2$Sigma_e} (the model's covariance at \code{params} when that is
#' absent), times \code{shock_size} -- the \code{compute_irfs()} convention.
#' A shock with zero variance gets a unit impulse of size \code{shock_size},
#' as in \code{compute_irfs()}. The second-order response is the difference of
#' two pruned simulations (\code{simulate_model_order2()}) with and without the
#' impulse, which removes the \code{ghss} drift exactly.
#'
#' @param dr2         DecisionRules2 object
#' @param model       dynhr_mod
#' @param shock_name  Name of shock to display (default: first shock)
#' @param vars        Variables to plot (default: all)
#' @param n_periods   IRF horizon
#' @param params      Named parameter vector (default \code{model$param_values})
#' @param shock_size  Shock size in standard deviations (default 1)
#' @param meta        Optional metadata list attached to the diagnostic result
#' @return A dynhr_diagnostic object with \code{plots$irf_comparison}
#' @keywords internal
#' @export
d19_irf_comparison <- function(dr2, model, shock_name = NULL,
                                vars = NULL, n_periods = 40L,
                                params = NULL, shock_size = 1, meta = NULL) {
  if (!inherits(dr2, "DecisionRules2")) {
    return(.make_result(
      result  = NULL,
      pass    = FALSE,
      summary = "dr2 is not a DecisionRules2 object.",
      llm_summary = "[FAIL] d19_irf_comparison: need DecisionRules2"
    ))
  }
  exo  <- dr2$exo_names
  endo <- dr2$endo_names
  if (is.null(shock_name)) shock_name <- exo[1]
  if (!shock_name %in% exo)
    .dynhr_abort(sprintf("Shock '%s' not in model exogenous variables.", shock_name))
  n_periods <- as.integer(n_periods)
  if (n_periods < 1L) .dynhr_abort("n_periods must be >= 1.")
  if (is.null(vars)) vars <- endo
  vars <- intersect(vars, endo)
  if (length(vars) == 0L) .dynhr_abort("None of `vars` are endogenous variables.")
  if (is.null(params)) params <- model$param_values

  Sigma_e <- .d19_sigma_e(dr2, model, params)
  n_u <- length(exo)
  k   <- match(shock_name, exo)
  ev  <- eigen((Sigma_e + t(Sigma_e)) / 2, symmetric = TRUE, only.values = TRUE)$values
  L   <- if (min(ev) > 1e-14 * max(1, max(ev))) t(chol(Sigma_e)) else
    diag(sqrt(pmax(diag(Sigma_e), 0)), nrow = n_u)
  eps <- L[, k] * shock_size
  if (all(eps == 0) && shock_size != 0) {
    eps <- numeric(n_u)
    eps[k] <- shock_size
  }

  # First order: the package's standard IRF, fed the same covariance.
  dr_first <- dr2
  dr_first$Sigma_e <- Sigma_e
  # A10 (0.9.4): compute_irfs() now gives caller `params` PRECEDENCE over
  # dr$Sigma_e (the 0.9.3.7 Kalman-path rule, applied to both IRF orders).
  # D19 resolves the covariance ITSELF in .d19_sigma_e() -- which prefers
  # dr2$Sigma_e -- and has already pinned it on dr_first, so passing `params`
  # on as well would make compute_irfs() rebuild a DIFFERENT covariance from
  # the .mod and break the shock-vector check just below.  Pass params = NULL
  # so dr_first$Sigma_e wins; `Sigma_e` already reflects `params` whenever
  # dr2 carries no covariance of its own.
  irf1 <- compute_irfs(dr_first, model, n_periods = n_periods,
                       shock_size = shock_size, params = NULL)[[shock_name]]
  if (!isTRUE(all.equal(as.numeric(irf1[1, ]),
                        as.numeric(dr2$ghu %*% eps), tolerance = 1e-10)))
    .dynhr_abort("d19_irf_comparison: first-order impact does not match ghu %*% eps ",
                 "(shock-vector convention drifted from compute_irfs()).")

  # Second order (pruned): impulse path minus no-impulse path, same eps.
  E1 <- matrix(0, n_periods, n_u)
  E1[1, ] <- eps
  E0 <- matrix(0, n_periods, n_u)
  old <- options(dynhr.warn_near_unit_root = FALSE)
  on.exit(options(old), add = TRUE)
  irf2 <- simulate_model_order2(dr2, n_periods = n_periods, shocks = E1,
                                model = model, burn_in = 0L, pruning = TRUE) -
          simulate_model_order2(dr2, n_periods = n_periods, shocks = E0,
                                model = model, burn_in = 0L, pruning = TRUE)
  attr(irf2, "levels") <- NULL
  irf1 <- as.matrix(irf1)
  irf2 <- as.matrix(irf2)
  colnames(irf2) <- endo
  rownames(irf1) <- rownames(irf2) <- paste0("t", seq_len(n_periods))

  rel_diffs <- vapply(vars, function(v) {
    denom <- max(abs(irf1[, v]))
    if (denom < 1e-14) return(NA_real_)
    max(abs(irf2[, v] - irf1[, v])) / denom
  }, numeric(1))
  max_rel_diff <- if (any(!is.na(rel_diffs))) max(rel_diffs, na.rm = TRUE) else NA_real_
  top_vars <- vars[order(-rel_diffs, na.last = NA)]
  top_vars <- intersect(top_vars, vars[!is.na(rel_diffs) & rel_diffs > 1e-8])
  top_vars <- utils::head(top_vars, 3L)

  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    lab1 <- "First order"
    lab2 <- "Second order (pruned)"
    df <- data.frame(
      period   = rep(rep(seq_len(n_periods), 2), length(vars)),
      value    = unlist(lapply(vars, function(v) c(irf1[, v], irf2[, v]))),
      order    = factor(rep(rep(c(lab1, lab2), each = n_periods), length(vars)),
                        levels = c(lab1, lab2)),
      variable = factor(rep(vars, each = 2L * n_periods), levels = vars)
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$period, y = .data$value,
                                          colour = .data$order,
                                          linetype = .data$order)) +
      geom_dynhr_zero() +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ variable, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = function(l) unique(round(pretty(l)))) +
      ggplot2::scale_y_continuous(n.breaks = 4,
                                  labels = function(b) formatC(b, format = "g", digits = 2)) +
      ggplot2::labs(
        title    = sprintf("D19: IRFs to a %s-sd %s shock, first vs second order",
                           format(shock_size), shock_name),
        subtitle = sprintf(paste0("Deviation from steady state, net of the constant ",
                                  "risk correction 0.5 ghss. Max gap / max |1st order| ",
                                  "= %s%s"),
                           if (is.na(max_rel_diff)) "NA" else sprintf("%.3g", max_rel_diff),
                           if (length(top_vars)) sprintf(" (%s)", top_vars[1]) else ""),
        x = "Periods after impact (1 = impact)", y = "Deviation from steady state",
        colour = NULL, linetype = NULL) +
      scale_colour_dynhr_vibrant() +
      theme_dynhr_diagnostic()
    plots$irf_comparison <- .apply_meta(p, meta)
  }

  summary_text <- paste0(
    sprintf("D19 IRF comparison: shock=%s (impulse sd-scaled, size %s), periods=%d\n",
            shock_name, format(shock_size), n_periods),
    sprintf("  Max relative deviation (2nd vs 1st order): %s\n",
            if (is.na(max_rel_diff)) "NA" else sprintf("%.4f", max_rel_diff)),
    sprintf("  Variables with largest deviation: %s",
            if (length(top_vars)) paste(top_vars, collapse = ", ") else "none"))

  llm_text <- sprintf(
    "[INFO] d19_irf_comparison: shock=%s n_periods=%d max_rel_diff_2nd_vs_1st=%s",
    shock_name, n_periods,
    if (is.na(max_rel_diff)) "NA" else sprintf("%.4f", max_rel_diff))

  .make_result(
    result      = list(irf1 = irf1, irf2 = irf2, rel_diffs = rel_diffs,
                       shock = setNames(eps, exo)),
    pass        = NA,
    plots       = plots,
    summary     = summary_text,
    llm_summary = llm_text
  )
}
