## dynhr-package.R
## --------------------------------------------------------------------------
## Package-level documentation and roxygen tags.
## --------------------------------------------------------------------------

#' dynhr: DSGE Modelling Toolkit
#'
#' Parses Dynare `.mod` files, solves the resulting models by first-order
#' perturbation, runs Bayesian estimation (random-walk Metropolis-Hastings,
#' sequential Monte Carlo, NUTS), and reports a battery of identification,
#' convergence, fit, and narrative diagnostics.
#'
#' The public API is exported and documented. Over 90 functions are
#' exported covering parsing, solving, estimation, filtering, diagnostics,
#' optimal policy (Ramsey / OSR / discretionary), OBC, welfare analysis,
#' and reporting. Functions still marked internal (no `@export`) remain
#' accessible via `dynhr:::fn`.
#'
#' @keywords internal
#' @useDynLib dynhr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import parallel
#' @importFrom stats acf approx coef complete.cases cor cov dbeta density dgamma dnorm dunif integrate lm.fit optim pchisq plogis pnorm predict qlogis qnorm quantile rbeta rcauchy rchisq reorder reshape rgamma rnorm runif sd setNames sigma simulate spline uniroot var
#' @importFrom utils capture.output head modifyList read.csv tail
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline hist legend lines mtext par points polygon text
"_PACKAGE"

## Package-level constants for numerical thresholds.
## Single source of truth — update here to change everywhere.

## Kalman filter: log-likelihood values below this threshold trigger a
## failure (explosive or degenerate observation variance).
.KF_LL_MIN  <- -1e8

## Linear algebra: default SVD tolerance for pseudo-inverse.
.PINV_TOL   <- 1e-12

## Lyapunov / Riccati: convergence tolerance.
.LYAP_TOL   <- 1e-12
.DARE_TOL   <- 1e-14

## Jacobian / Hessian regularisation: small diagonal ridge added when
## matrices are near-singular.
.JAC_REG    <- 1e-8

## Diffuse Kalman initialisation: multiplier on identity for unit-root
## states when the Lyapunov solution is unavailable.
.DIFFUSE_SCALE <- 1e6

## suppress R CMD check NOTE about "no visible binding for global variable"
## introduced by ggplot2 / data.table NSE inside the monolith files. These
## are temporary; once each file is split they get proper
## `@importFrom rlang .data` annotations.
utils::globalVariables(c(
  ".", ".data", ".N", ".SD",
  "value", "variable", "shock", "horizon", "period", "share",
  "parameter", "draw", "chain", "iteration",
  "x", "y", "ymin", "ymax", "lower", "upper", "mean", "median",
  "ess", "rhat", "logpost", "loglik", "logprior",
  "dynhr_colours", "dynhr_palette",
  ## data.table NSE helpers / functions used unqualified in the transform and
  ## gradient-primitive monoliths (data.table is in Suggests; these calls are
  ## reached only when it is installed):
  ":=", "set", "copy", "data.table", "fread", "rbindlist", "as.IDate",
  "shift", "i.y_trend",
  ## ggplot2 / diagnostic-plot and summary-table NSE column names:
  "ACF", "Coef", "Hi", "Lag", "Lo", "Log_MargLik", "Model", "Moment",
  "Regime", "SD", "Source", "Time", "Type", "Value", "Variable", "Verdict",
  "axis", "basin", "bfmi", "bin", "bin_mid", "block", "calibrated", "corr",
  "count", "depth", "flat", "freq", "grade", "h", "hi", "implied", "index",
  "iter", "lab", "label", "lag", "lambda", "ll", "lo", "loading", "log_BF",
  "med", "meta", "model", "obs", "param", "q05", "q25", "q75", "q95", "r",
  "set", "sloppy", "softness", "status", "target", "theta", "type",
  "value_med", "var1", "var2", "verdict", "abs_E", "abs_loading",
  "ln_norm", "ln_trend",
  ## NZSIM neutral-rate / data-prep series names (transform-monolith):
  ".worker_model", ".x", "b_", "b_ratio", "b_ratio_trend", "c_", "c_sh",
  "c_sh_trend", "c_trend", "dp_", "dp_actual", "dp_trend", "g_", "g_sh",
  "g_sh_trend", "g_trend", "ih_", "ih_sh", "ih_sh_trend", "ih_trend", "ik_",
  "ik_sh", "ik_sh_trend", "ik_trend", "iwgdp_pt", "iwgdp_z", "lhpwa_z",
  "llisai", "lmig_z", "ln_", "m_", "m_sh", "m_sh_trend", "m_trend", "ncg_z",
  "ncp_z", "ngdpp_z", "ngdpz", "ngdpz_smooth", "nik_z", "nitd_z", "nm_z",
  "nx_z", "p_trend", "pcpis", "ph_p_", "ph_rel", "ph_rel_trend", "pm_ps",
  "pm_ps_trend", "pmstar_", "pmstar_trend", "pn_p_", "pn_rel", "pn_rel_trend",
  "pnt", "pqhpiz", "pstar_", "pstar_trend", "px_ps", "px_ps_trend", "pxstar_",
  "pxstar_trend", "r90d", "r_", "r_anchor", "r_constructed", "r_q", "r_trend",
  "rh_", "rh_anchor", "rh_q", "rh_trend", "rs_", "rs_real", "rs_trend",
  "rshortw_ocr", "rstar_", "rstar_q", "rstar_trend", "rtwi", "tiin", "wcpi",
  "w_real", "w_real_trend", "wrlci_", "x_", "x_sh", "x_sh_trend", "x_trend",
  "y_", "y_trend", "ystar_", "ystar_trend",
  ## legacy exists()-guarded compatibility fallback in .build_log_posterior():
  "log_posterior_old",
  ## diagnostic-plot / summary-table NSE columns, data.table helpers, mirai
  ## worker <<- bindings, and transform-monolith intermediates:
  "..obs_cols", ".worker_Y", ".worker_cm",
  "Band", "Binding", "Calibrated", "Contribution", "Diff", "Difference",
  "Distribution", "ElemLabel", "Equation", "Estimated", "Frequency",
  "Identified", "KL", "Lag_coef", "MeanAbsStrength", "Measure", "Median",
  "Observable", "Parameter", "Period", "Perturbation", "Posterior_Prob",
  "Sample", "Scenario", "Share", "Shock", "Strength", "Welfare", "abs_s",
  "abs_s_plot", "colour", "delta_i_list",
  "dominant_shock",
  "fragile_label", "globally_identified", "high_prec", "ho_gain",
  "ho_important", "identified", "informative", "kl_divergence", "ll_vec",
  "lockdown_ratio", "lp_vec", "min_kl_plot", "moment", "mu_star_plot",
  "ngdpp_z_smooth", "param_group", "param_label", "param_type", "perturbation",
  "pmstar_constructed", "precision", "pxstar_constructed", "rh_constructed",
  "rstar_anchor", "setnames", "singular_value", "value_lo_out", "value_sd",
  "xmax", "xmin"
))
