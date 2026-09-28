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
##     (homogeneity test of the regime posterior means not rejected). A
##     non-invariant private parameter is a literal Lucas-critique violation
##     -> the "deep" parameter is absorbing the policy change and is not deep.
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
##   Hausman, J. A. (1978). Specification tests in econometrics.
##     Econometrica, 46(6), 1251-1271.
##   Inoue, A., & Rossi, B. (2011). Identifying the sources of instabilities in
##     macroeconomic fluctuations. Review of Economics and Statistics, 93(4).
## --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#' D34. Policy-partitioned invariance test (operational Lucas critique)
#'
#' **Statistic.** For each parameter the regime posterior means are compared
#' with a homogeneity (Cochran-Q / Wald) statistic. With K >= 2 regimes the
#' regimes are taken to be NON-OVERLAPPING samples, so (Bernstein-von Mises)
#' their posterior means are independent with variances \eqn{V_i}:
#' \deqn{Q_j = \sum_i w_{ij} (\bar\theta_{ij} - \tilde\theta_j)^2,\quad
#'   w_{ij} = 1 / (V_{ij} + V_{ij}/ESS_{ij})}
#' with \eqn{\tilde\theta_j} the \eqn{w}-weighted mean, and
#' \eqn{Q_j \sim \chi^2_{K-1}} under invariance (for K = 2,
#' \eqn{Q = z^2}). The \eqn{V/ESS} term is the Monte Carlo variance of each
#' posterior mean. Overlapping regimes make the test conservative (the shared
#' data correlate the means positively); a shared informative prior does too.
#' With ONE regime the regime is a subset of the full sample, so -- as in D16
#' -- the Hausman variance of the difference is \eqn{V_{sub} - V_{full}}
#' (plus Monte Carlo variance), and \eqn{Q = z^2 \sim \chi^2_1}.
#' Credible-interval overlap is reported (\code{ci_overlap}) but is not a
#' test: its size depends on the relative widths and on K, and for the nested
#' one-regime case it essentially never fires.
#'
#' **D34 regimes vs D16 nested subsamples.** The two diagnostics answer
#' different questions and take their samples differently, so pick by the
#' question, not by which inputs you happen to have.
#' \itemize{
#'   \item \strong{D16} splits the sample into NESTED subsamples (a sub-period
#'     against the full sample) and tests EVERY estimated parameter
#'     symmetrically, with the Hausman variance \eqn{V_{sub} - V_{full}}. It
#'     asks "is the estimate stable?" and needs no economic story about the
#'     break.
#'   \item \strong{D34} splits the sample into NON-OVERLAPPING policy REGIMES
#'     and tests them DIRECTIONALLY against the deep-parameter partition: the
#'     private block is expected to be invariant (its Q statistics carry the
#'     verdict), while the policy block is expected to MOVE and is reported as
#'     context, never as a failure. It needs a substantive regime split --
#'     a genuine, dateable change in the policy rule -- and the
#'     \code{@dynhr:deep} partition that says which parameters are policy.
#'     Handing D34 arbitrary sub-periods makes its directional reading
#'     meaningless even though the arithmetic still runs.
#' }
#' D34's one-regime fallback is deliberately D16's nested Hausman comparison,
#' so a single regime gives the D16 statistic read through the policy/private
#' partition.
#'
#' **Why the Engle-Hendry-Richard super-exogeneity test is direct-call only.**
#' \code{super_exog} is never built by \code{run_all_diagnostics()}; the
#' caller must supply the regime vector and the shock matrices. Two reasons.
#' (i) The classical EHR / Hendry-Santos F-test assumes iid, correctly-sized
#' residuals under the null. Smoothed DSGE shocks are neither: they are
#' serially correlated and variance-shrunk even under a correctly specified
#' model (the same Durbin-Koopman law-of-total-variance property that keeps
#' Ljung-Box out of D12's gate), so the F-test's nominal size does not hold on
#' them. (ii) Only the caller knows which dates constitute a policy regime
#' change; Favero & Hendry (1992) show that only LOCATION shifts in the
#' marginal process are informative for super-exogeneity, so an automatically
#' chosen break would usually be testing the wrong thing. No DSGE toolchain
#' (Dynare, IRIS, DSGE.jl, RISE) automates this test either.
#'
#' **Multiple testing.** The verdict family is the private-block invariance
#' tests plus (when supplied) the super-exogeneity equations; each is tested
#' at \code{family_alpha / n_family} (Bonferroni). The policy block ("moved")
#' and the auxiliary block (context) are each Bonferroni-corrected within
#' their own block at \code{family_alpha}.
#'
#' @param model       Parsed model (for the @dynhr:deep policy/private
#'   partition).  Optional if \code{deep_spec} is supplied.
#' @param deep_spec   Optional \code{\link{build_deep_spec}} (built from
#'   \code{model} and the draw column names otherwise).
#' @param draws       Full-sample posterior draws (\eqn{n \times p}) with named
#'   columns -- the same object D16 receives.
#' @param results_sub Named list of per-regime posterior draw matrices, e.g.
#'   \code{list("Pre-1990" = d1, "Post-1990" = d2)} -- the same object D16
#'   receives. Columns are matched to \code{draws} BY NAME when both carry
#'   column names (otherwise by position). Two or more (non-overlapping)
#'   regimes give an across-regime test; one regime falls back to a nested
#'   regime-vs-full Hausman comparison.
#' @param param_names Optional parameter names (defaults to \code{colnames(draws)}).
#' @param ci_level    Credible-interval level for the plotted/reported
#'   intervals (default 0.90). Not used by the verdict.
#' @param family_alpha Family-wise false-alarm rate of the verdict (default
#'   0.05); see Details.
#' @param super_exog  Optional list for the Engle-Hendry-Richard test:
#'   \code{list(shocks = <T x k matrix with named columns, or named list of
#'   per-regime matrices>, policy = <character names of policy innovations>,
#'   regime = <length-T factor; required when shocks is a single matrix>)}.
#' @param meta        Optional \code{\link{diag_meta}} provenance descriptor.
#' @references
#'   Engle, R. F., Hendry, D. F. & Richard, J.-F. (1983). Exogeneity.
#'   \emph{Econometrica}, 51(2), 277-304.
#'   Favero, C. & Hendry, D. F. (1992). Testing the Lucas critique: a review.
#'   \emph{Econometric Reviews}, 11(3), 265-306.
#'   Hausman, J. A. (1978). Specification tests in econometrics.
#'   \emph{Econometrica}, 46(6), 1251-1271.
#' @return A \code{dynhr_diagnostic}; \code{$result$invariance_table} holds
#'   one row per parameter (statistic, df, p, threshold, verdict);
#'   \code{$result$passport_axis} is a named logical over the private block
#'   (TRUE = invariant) for the Passport.
#' @noRd
d34_policy_invariance <- function(model        = NULL,
                                  deep_spec    = NULL,
                                  draws        = NULL,
                                  results_sub  = NULL,
                                  param_names  = NULL,
                                  ci_level     = 0.90,
                                  family_alpha = 0.05,
                                  super_exog   = NULL,
                                  meta         = NULL) {
  if (is.null(draws))
    return(.make_result(pass = NA,
      summary = "D34 policy invariance: no full-sample draws supplied."))
  if (is.null(results_sub) || length(results_sub) == 0)
    return(.make_result(pass = NA,
      summary = paste("D34 policy invariance: no per-regime draws (results_sub);",
                      "supply >=1 regime split (shared with D16).")))
  if (!is.list(results_sub) || is.data.frame(results_sub))
    .dynhr_abort("D34: `results_sub` must be a list of per-regime draw matrices.")
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !(ci_level > 0 && ci_level < 1))
    .dynhr_abort("D34: `ci_level` must be a single number in (0, 1).")
  if (!is.numeric(family_alpha) || length(family_alpha) != 1L ||
      !(family_alpha > 0 && family_alpha < 1))
    .dynhr_abort("D34: `family_alpha` must be a single number in (0, 1).")

  draws <- as.matrix(draws)
  full_cn <- colnames(draws)
  if (is.null(param_names))
    param_names <- full_cn %||% paste0("theta_", seq_len(ncol(draws)))
  if (length(param_names) != ncol(draws))
    .dynhr_abort(sprintf("D34: %d param_names for %d full-sample columns.",
                         length(param_names), ncol(draws)))
  colnames(draws) <- param_names

  reg_names <- names(results_sub)
  if (is.null(reg_names)) reg_names <- rep("", length(results_sub))
  empty <- is.na(reg_names) | !nzchar(reg_names)
  reg_names[empty] <- paste0("Regime_", which(empty))
  reg_names <- make.unique(c("Full sample", reg_names), sep = "_")[-1L]

  # Align every regime's columns with the full sample: BY NAME when both carry
  # names (a positional match silently compares different parameters when the
  # column orders differ), otherwise by position.
  results_sub <- lapply(seq_along(results_sub), function(i) {
    m <- as.matrix(results_sub[[i]])
    if (!is.null(full_cn) && !is.null(colnames(m))) {
      miss <- setdiff(full_cn, colnames(m))
      if (length(miss) > 0L)
        .dynhr_abort(sprintf("D34: regime '%s' lacks column(s): %s",
                             reg_names[i], paste(miss, collapse = ", ")))
      m <- m[, full_cn, drop = FALSE]
    } else if (ncol(m) != length(param_names)) {
      .dynhr_abort(sprintf("D34: regime '%s' has %d columns, full sample %d.",
                           reg_names[i], ncol(m), length(param_names)))
    }
    colnames(m) <- param_names
    m
  })
  names(results_sub) <- reg_names
  for (m in c(list(draws), results_sub)) {
    if (nrow(m) < 4L || any(!is.finite(m)))
      .dynhr_abort("D34: every draw matrix needs >= 4 rows and finite values.")
  }

  if (is.null(deep_spec))
    deep_spec <- build_deep_spec(model = model, param_names = param_names)

  # Partition the *estimated* parameters.
  ds  <- deep_spec[match(param_names, deep_spec$param), , drop = FALSE]
  ds$class[is.na(ds$class)] <- "unknown"
  block <- ifelse(ds$class == "policy", "policy",
           ifelse(ds$is_deep %in% TRUE, "private", "auxiliary"))
  block[is.na(block)] <- "auxiliary"
  names(block) <- param_names
  private_p <- param_names[block == "private"]
  policy_p  <- param_names[block == "policy"]
  aux_p     <- param_names[block == "auxiliary"]
  n_private <- length(private_p)
  nested    <- length(results_sub) < 2L

  # --- per-parameter homogeneity statistic ---
  test_draws <- if (nested) c(list("Full sample" = draws), results_sub) else results_sub
  summ_test <- lapply(names(test_draws), function(nm)
    .summarise_draws_ci(test_draws[[nm]], nm, param_names, ci_level))
  names(summ_test) <- names(test_draws)
  mom <- lapply(test_draws, .d34_moments)
  means <- vapply(mom, function(x) x$mean, numeric(length(param_names)))
  vars  <- vapply(mom, function(x) x$var,  numeric(length(param_names)))
  mcvar <- vapply(mom, function(x) x$mcvar, numeric(length(param_names)))
  means <- matrix(means, nrow = length(param_names))
  vars  <- matrix(vars,  nrow = length(param_names))
  mcvar <- matrix(mcvar, nrow = length(param_names))
  stat <- vapply(seq_along(param_names), function(j) {
    if (nested) {
      # column 1 = full sample, column 2 = the nested regime (Hausman). The
      # variance DIFFERENCE is itself a Monte Carlo estimate whose error is
      # large relative to it at low ESS (Var(V_hat) ~ 2 V^2 / ESS); one MC
      # standard error of the difference is added so noisy draws do not
      # shrink the denominator and over-reject.
      vd_se <- sqrt(2 * (vars[j, 1] * mcvar[j, 1] + vars[j, 2] * mcvar[j, 2]))
      v <- max(vars[j, 2] - vars[j, 1], 0) + vd_se + mcvar[j, 1] + mcvar[j, 2]
      d <- means[j, 2] - means[j, 1]
      return(if (v > 0) d^2 / v else if (d == 0) 0 else Inf)
    }
    v <- vars[j, ] + mcvar[j, ]
    # a constant (zero-variance) regime: any difference in means is decisive
    if (any(v <= 0)) return(if (diff(range(means[j, ])) == 0) 0 else Inf)
    w <- 1 / v
    pooled <- sum(w * means[j, ]) / sum(w)
    sum(w * (means[j, ] - pooled)^2)
  }, numeric(1))
  df_q <- if (nested) 1L else length(results_sub) - 1L
  pval <- stats::pchisq(stat, df = df_q, lower.tail = FALSE)
  names(stat) <- names(pval) <- param_names

  ci_overlap <- vapply(seq_along(param_names), function(j) {
    lo <- vapply(summ_test, function(s) s$lo[j], numeric(1))
    hi <- vapply(summ_test, function(s) s$hi[j], numeric(1))
    max(lo) <= min(hi)
  }, logical(1))

  # --- optional super-exogeneity test (errors propagate: a malformed spec
  # must not be silently dropped into a PASS) ---
  se <- NULL
  if (!is.null(super_exog))
    se <- .d34_super_exogeneity(super_exog, alpha = family_alpha)
  n_se <- if (is.null(se)) 0L else nrow(se$table)

  # --- Bonferroni thresholds ---
  n_family  <- n_private + n_se
  alpha_fam <- if (n_family > 0) family_alpha / n_family else NA_real_
  alpha_blk <- stats::setNames(rep(NA_real_, length(param_names)), param_names)
  alpha_blk[private_p] <- alpha_fam
  if (length(policy_p)) alpha_blk[policy_p] <- family_alpha / length(policy_p)
  if (length(aux_p))    alpha_blk[aux_p]    <- family_alpha / length(aux_p)
  invariant <- pval > alpha_blk
  names(invariant) <- param_names

  if (!is.null(se)) {
    se$alpha_each <- alpha_fam
    se$table$reject <- se$table$p <= alpha_fam
    se$rejected <- any(se$table$reject)
  }
  se_rejected <- !is.null(se) && isTRUE(se$rejected)

  private_unstable <- private_p[!invariant[private_p]]   # Lucas violations
  policy_moved     <- policy_p[!invariant[policy_p]]     # expected
  policy_static    <- policy_p[invariant[policy_p]]      # suspicious

  # Verdict: the Lucas critique holds iff no private parameter is
  # non-invariant and super-exogeneity (if tested) is not rejected.
  pass <- if (n_family == 0) NA
          else (length(private_unstable) == 0 && !se_rejected)

  passport_axis <- if (n_private > 0)
    stats::setNames(invariant[private_p], private_p) else logical(0)

  z_eq <- stats::qnorm(pval / 2, lower.tail = FALSE)
  z_crit <- stats::qnorm(alpha_blk / 2, lower.tail = FALSE)
  inv_tab <- data.frame(
    param = param_names, block = unname(block),
    stat = unname(stat), df = df_q, p = unname(pval),
    z_equiv = unname(z_eq), alpha = unname(alpha_blk),
    z_crit = unname(z_crit), invariant = unname(invariant),
    ci_overlap = ci_overlap, stringsAsFactors = FALSE)

  # --- plots ---
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_df <- do.call(rbind, c(
      list(.summarise_draws_ci(draws, "Full sample", param_names, ci_level)),
      lapply(names(results_sub), function(nm)
        .summarise_draws_ci(results_sub[[nm]], nm, param_names, ci_level))))
    plot_df$block <- unname(block[plot_df$param])
    plots$forest <- .plot_d34_forest(plot_df, names(results_sub), ci_level, meta)
    plots$zstat  <- .plot_d34_zstat(inv_tab, nested, family_alpha, meta)
  }

  fmt_p <- function(p) paste(sprintf("%s (p=%.2g)", p, pval[p]), collapse = ", ")
  summary_txt <- sprintf(
    "D34 Policy invariance: %d private / %d policy / %d aux params, %d regime(s)%s; family alpha %.2g (Bonferroni over %d test(s)). %s%s",
    n_private, length(policy_p), length(aux_p), length(results_sub),
    if (nested) " (nested regime-vs-full Hausman test)" else "",
    family_alpha, n_family,
    if (is.na(pass)) "No private deep parameters estimated (nothing to test)."
    else if (isTRUE(pass)) "PASS -- private block invariant across regimes."
    else sprintf("FAIL -- Lucas violation: %s",
                 paste(c(if (length(private_unstable)) fmt_p(private_unstable),
                         if (se_rejected) "super-exogeneity rejected"),
                       collapse = ", ")),
    if (length(policy_static) > 0)
      sprintf(" (caveat: policy params %s did not move -- check the break date)",
              paste(policy_static, collapse = ", ")) else "")

  badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
  llm <- paste(c(
    sprintf("D34 | Policy Invariance (Lucas) | %s", badge),
    sprintf("  private=%d policy=%d aux=%d regimes=%d%s test=chi2(df=%d) family_alpha=%.2g n_family=%d",
            n_private, length(policy_p), length(aux_p), length(results_sub),
            if (nested) " (nested-vs-full)" else "", df_q, family_alpha, n_family),
    if (length(private_unstable))
      sprintf("  lucas_violation(private non-invariant): %s", fmt_p(private_unstable)),
    if (length(policy_moved))
      sprintf("  policy_responded(expected): %s", fmt_p(policy_moved)),
    if (length(policy_static))
      sprintf("  policy_static(caveat): %s", paste(policy_static, collapse = ", ")),
    if (!is.null(se))
      sprintf("  super_exogeneity: %s (%d/%d private eqns load on policy at p <= %.2g)",
              if (se_rejected) "REJECTED" else "ok",
              sum(se$table$reject), n_se, alpha_fam),
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
                  invariance_table = inv_tab,
                  private_unstable = private_unstable,
                  policy_moved = policy_moved, policy_static = policy_static,
                  super_exogeneity = se, passport_axis = passport_axis,
                  nested = nested, summaries = summ_test),
    pass    = pass,
    plots   = plots,
    summary = summary_txt,
    llm_summary = llm)
}


# Posterior mean, variance and Monte Carlo variance of the mean (variance /
# ESS, ESS capped at the number of draws) for each column of a draw matrix.
.d34_moments <- function(m) {
  n <- nrow(m)
  v <- apply(m, 2, stats::var)
  ess <- vapply(seq_len(ncol(m)), function(j) {
    e <- if (v[j] > 0) .d5_ess_basic(m[, j, drop = FALSE]) else NA_real_
    if (is.finite(e) && e > 0) min(e, n) else n
  }, numeric(1))
  list(mean = colMeans(m), var = unname(v), mcvar = unname(v) / ess)
}


# Engle-Hendry-Richard super-exogeneity: do policy innovations enter the
# private structural residuals? For each private shock e_p and each regime r,
# F-test  e_p ~ 1 + (policy innovations)  against  e_p ~ 1  on regime r's
# rows. The per-regime tests are independent (disjoint rows), so their
# p-values are combined with Fisher's method (-2 sum log p ~ chi2_{2R}), which
# stays exact when the residual variance differs across regimes (a pooled
# regime-interacted F test does not: it over-rejects under a volatility
# break). `rejected` uses a Bonferroni threshold alpha / (number of equations).
.d34_super_exogeneity <- function(super_exog, alpha = 0.10) {
  if (!is.list(super_exog) || is.null(super_exog$shocks))
    .dynhr_abort("D34: `super_exog` must be a list with `shocks` and `policy`.")
  shocks <- super_exog$shocks
  if (is.list(shocks) && !is.data.frame(shocks)) {
    if (is.null(names(shocks)) || any(!nzchar(names(shocks))) ||
        anyDuplicated(names(shocks)))
      .dynhr_abort("D34: super_exog$shocks as a list needs unique regime names.")
    mats <- lapply(shocks, as.matrix)
    cn <- colnames(mats[[1]])
    if (is.null(cn))
      .dynhr_abort("D34: super_exog$shocks needs column names.")
    mats <- lapply(names(mats), function(nm) {
      m <- mats[[nm]]
      if (is.null(colnames(m)) || !setequal(colnames(m), cn))
        .dynhr_abort(sprintf("D34: super_exog$shocks[['%s']] columns differ from the first regime's.", nm))
      m[, cn, drop = FALSE]                   # align BY NAME
    })
    regime <- factor(rep(names(shocks), vapply(mats, nrow, integer(1))),
                     levels = names(shocks))
    shocks <- do.call(rbind, mats)
  } else {
    shocks <- as.matrix(shocks)
    regime <- super_exog$regime
    if (is.null(regime))
      .dynhr_abort("D34: super_exog$regime is required when shocks is a single matrix.")
    if (length(regime) != nrow(shocks))
      .dynhr_abort(sprintf("D34: super_exog$regime has length %d, shocks have %d rows.",
                           length(regime), nrow(shocks)))
    regime <- droplevels(as.factor(regime))
  }
  if (!is.numeric(shocks) || any(!is.finite(shocks)) || anyNA(regime))
    .dynhr_abort("D34: super_exog$shocks and $regime must be finite / non-missing.")
  all_names <- colnames(shocks)
  if (is.null(all_names))
    .dynhr_abort("D34: super_exog$shocks needs column names.")
  pol_names <- intersect(super_exog$policy, all_names)
  if (length(pol_names) == 0)
    .dynhr_abort("D34: no policy-shock columns matched super_exog$policy.")
  priv_names <- setdiff(all_names, pol_names)
  if (length(priv_names) == 0)
    .dynhr_abort("D34: super_exog has no private-shock columns.")
  if (nlevels(regime) < 2)
    .dynhr_abort("D34: super_exog$regime needs >= 2 levels.")
  q <- length(pol_names)
  n_r <- table(regime)
  if (any(n_r < q + 3L))
    .dynhr_abort(sprintf("D34: every regime needs >= %d rows for the super-exogeneity F test.",
                         q + 3L))

  # Per-regime QR of [1, policy]; reused for every private equation.
  qrs <- lapply(levels(regime), function(lv) {
    X <- cbind(1, shocks[regime == lv, pol_names, drop = FALSE])
    qx <- qr(X)
    if (qx$rank < ncol(X))
      .dynhr_abort(sprintf("D34: policy innovations are collinear/constant in regime '%s'.", lv))
    qx
  })
  names(qrs) <- levels(regime)

  rows <- lapply(priv_names, function(pn) {
    p_r <- vapply(levels(regime), function(lv) {
      y <- shocks[regime == lv, pn]
      rss_red  <- sum((y - mean(y))^2)
      rss_full <- sum(qr.resid(qrs[[lv]], y)^2)
      df2 <- length(y) - q - 1L
      if (rss_full <= 0) return(if (rss_red > 0) 0 else 1)
      Fst <- ((rss_red - rss_full) / q) / (rss_full / df2)
      stats::pf(Fst, q, df2, lower.tail = FALSE)
    }, numeric(1))
    fisher <- -2 * sum(log(pmax(p_r, .Machine$double.xmin)))
    data.frame(shock = pn, stat = fisher, df = 2L * length(p_r),
               p = stats::pchisq(fisher, 2L * length(p_r), lower.tail = FALSE),
               min_regime_p = min(p_r), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  rejected <- any(tab$p <= alpha / nrow(tab))
  list(table = tab, rejected = rejected, alpha = alpha,
       policy = pol_names, private = priv_names)
}


.d34_block_labels <- c(private   = "Private (must be invariant)",
                       policy    = "Policy (expected to move)",
                       auxiliary = "Auxiliary / shocks (context)")

# Forest: one panel per parameter (free x -- parameters live on very different
# scales), panels ordered private -> policy -> auxiliary, block named in the
# strip; regimes on the y axis, full sample shown in grey for reference.
.plot_d34_forest <- function(plot_df, regimes, ci_level, meta) {
  blk_short <- c(private = "private", policy = "policy", auxiliary = "aux")
  blk_rank <- match(plot_df$block, names(blk_short))
  ord <- unique(plot_df$param[order(blk_rank)])
  panel_lab <- stats::setNames(
    sprintf("%s [%s]", ord, blk_short[plot_df$block[match(ord, plot_df$param)]]), ord)
  plot_df$panel <- factor(panel_lab[plot_df$param], levels = unname(panel_lab))
  lv <- c("Full sample", regimes)
  plot_df$sample <- factor(plot_df$sample, levels = rev(lv))
  reg_cols <- rep_len(setdiff(dynhr_palette, dynhr_na_colour), length(regimes))
  pal <- stats::setNames(c(dynhr_na_colour, reg_cols), lv)
  n_par <- length(ord)

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = median, y = sample, colour = sample)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi),
                           orientation = "y", width = 0.3, linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ panel, scales = "free_x",
                        ncol = min(4L, ceiling(sqrt(n_par)))) +
    ggplot2::scale_x_continuous(n.breaks = 3L) +
    ggplot2::scale_colour_manual(values = pal, breaks = lv, name = NULL) +
    theme_dynhr_compact() +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(1.2, "lines"),
                   legend.position = "bottom") +
    ggplot2::labs(
      title    = "D34: Policy invariance -- posterior by regime",
      subtitle = sprintf("median and %d%% CI | [private] must not move, [policy] is expected to | free x-scale",
                         round(100 * ci_level)),
      x = "Parameter value", y = NULL)
  .apply_meta(p, meta)
}

# Test statistic: equivalent |z| of each parameter's homogeneity test, one
# panel per block, with the block's Bonferroni critical value drawn.
.plot_d34_zstat <- function(tab, nested, family_alpha, meta) {
  tab$block_lab <- factor(.d34_block_labels[tab$block],
                          levels = unname(.d34_block_labels))
  tab$param <- factor(tab$param, levels = rev(unique(tab$param[order(
    match(tab$block, names(.d34_block_labels)))])))
  tab$Verdict <- factor(ifelse(tab$invariant, "invariant", "moved"),
                        levels = c("invariant", "moved"))
  finite_z <- tab$z_equiv[is.finite(tab$z_equiv)]
  crit <- unique(tab[, c("block_lab", "z_crit")])
  z_lim <- max(c(crit$z_crit, finite_z, 1), na.rm = TRUE) * 1.1
  z_lim <- min(z_lim, max(c(crit$z_crit * 3, 10), na.rm = TRUE))
  clipped <- any(tab$z_equiv > z_lim)
  tab$z_plot <- pmin(tab$z_equiv, z_lim)

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = z_plot, y = param)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = z_plot, yend = param),
                          colour = dynhr_na_colour, linewidth = 0.4) +
    ggplot2::geom_vline(data = crit, ggplot2::aes(xintercept = z_crit),
                        colour = dynhr_colours$red, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = Verdict, shape = Verdict), size = 2.6) +
    ggplot2::facet_grid(rows = ggplot2::vars(block_lab), scales = "free_y",
                        space = "free_y", drop = TRUE,
                        labeller = ggplot2::label_wrap_gen(width = 14)) +
    ggplot2::scale_colour_manual(values = c(invariant = dynhr_colours$mid_blue,
                                            moved = dynhr_colours$orange),
                                 name = NULL, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(invariant = 16, moved = 17),
                                name = NULL, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, z_lim)) +
    theme_dynhr() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0)) +
    ggplot2::labs(
      title = "D34: Regime homogeneity test per parameter",
      subtitle = sprintf(
        "%s | dashed = Bonferroni cut-off, family alpha %.2g%s",
        if (nested) "nested regime-vs-full Hausman |z|"
        else "|z| equivalent of the chi-square test across regimes",
        family_alpha, if (clipped) " | clipped" else ""),
      x = "|z| (equivalent standard errors)", y = NULL)
  .apply_meta(p, meta)
}
