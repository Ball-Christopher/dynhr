## R/diag-post-d41-innovation.R
## --------------------------------------------------------------------------
## D41 Kalman-filter innovation whiteness (per-observable standardised-
## innovation variance and lag-1 autocorrelation z-statistics).
##
## Moved out of R/diag-orchestrate.R in the 0.9.4 diagnostics refresh. The
## numerical work is kf_innovation_diagnostics() (R/kf-innovation-
## diagnostics.R); this file owns the evaluation point, the likelihood
## settings, the skip rules and the plots.
## --------------------------------------------------------------------------

#' D41. Kalman-filter innovation whiteness
#'
#' Runs \code{\link{kf_innovation_diagnostics}} on \code{data} and reports,
#' per observable \eqn{i}, three z-statistics of the standardised
#' one-step-ahead innovations \eqn{z_{i,t} = v_{i,t}/\sqrt{F_{ii,t}}}: their
#' mean (null 0, sd \eqn{1/\sqrt{T_i}}; added by D41 -- the variance below is
#' demeaned and so blind to a level offset), their sample variance (null 1,
#' \strong{kurtosis-robust} sd \eqn{\sqrt{(\hat\kappa_i - 1)/T_i}}) and their
#' lag-1 autocorrelation (null 0, sd \eqn{1/\sqrt{T_i}}). \eqn{F_t} is the
#' filter's own per-period
#' innovation covariance, so under a correctly specified model each
#' \eqn{z_{i,\cdot}} is iid N(0,1) from the first period (the stationary
#' initialisation is exact, so no burn-in is dropped). Any \eqn{|z| > 4}
#' gives a FAIL badge: at most \eqn{3 n_{obs}} tests at a per-test Gaussian
#' size of \eqn{6.3 \times 10^{-5}}, so the family-wise false-alarm rate stays
#' well under 1\% for any realistic observable count.
#'
#' \strong{Variance test calibration (0.9.4).} The general large-sample
#' variance of a sample variance is \eqn{Var(s^2) = (\mu_4 - \sigma^4)/T =
#' \sigma^4(\kappa - 1)/T}; the familiar \eqn{\sqrt{2/T}} is the special
#' case \eqn{\kappa = 3}. DSGE innovations are routinely fat-tailed even
#' under a correct model (non-Gaussian structural shocks, occasionally-binding
#' constraints, an imperfect linear approximation), and the Gaussian sd then
#' over-rejects badly -- at \eqn{\kappa = 9} (Student \eqn{t_5}) the z-statistic
#' is inflated by a factor of 2. D41 therefore plugs in the sample kurtosis
#' \eqn{\hat\kappa_i} (reported as \code{by_obs$kurt}). The robust sd is used
#' from 20 usable periods per observable; below that the estimated kurtosis is
#' too noisy and the Gaussian sd is kept. The 20-period switch is a PACKAGE
#' CHOICE with no literature source.
#'
#' \strong{Not implemented: the H(h) heteroskedasticity test.} Harvey (1989,
#' ch. 5) and Durbin & Koopman (2012, ch. 2.12/7) pair the normality and
#' Box-Ljung tests on standardised prediction errors with a heteroskedasticity
#' test H(h): the ratio of the sum of squares of the LAST third of the
#' residuals to that of the FIRST third, against \eqn{F(h, h)} with
#' \eqn{h \approx T/3}. That tests a DIFFERENT null -- variance DRIFT across
#' the sample -- from D41's "variance equals 1" against a fixed target. D41
#' does not implement H(h), so a model whose innovation variance is correct on
#' average but drifts over the sample can pass D41's variance check.
#'
#' \strong{Evaluation point.} When \code{compiled} is supplied the model is
#' re-solved at \code{params} (the same steady-state / perturbation /
#' stationarity pipeline as the estimation likelihood), so the decision rule
#' and \eqn{\Sigma_e} describe the same parameter vector. Without
#' \code{compiled}, \code{dr} must already be solved at \code{params}.
#'
#' \strong{Not evaluated (INFO)} when the thin filter cannot reproduce the
#' estimation likelihood: a non-trivial \code{me_extra} or
#' \code{shock_scale}, mixed-frequency \code{obs_aggregation}, a
#' non-Gaussian \code{likelihood}, a model that does not solve at
#' \code{params}, or a transition with a (near-)unit root whose stationary
#' covariance does not exist (the same rule \code{kalman_filter(lik_init =
#' "auto")} uses to fall back to the diffuse initialisation).
#'
#' @param data      Observations: \code{T x n_obs} (columns named by
#'   \code{obs_vars}; extra columns are ignored) or \code{n_obs x T}.
#' @param dr        Decision rule (output of \code{\link{solve_perturbation}}).
#' @param model     Parsed model object.
#' @param params    Named numeric parameter vector (a full vector, or the
#'   estimated subset -- it is merged into \code{model$param_values}).
#' @param obs_vars  Character vector of observed variable names.
#' @param compiled  Optional compiled model; when given, the model is
#'   re-solved at \code{params}.
#' @param me_variance Scalar measurement-error variance used by the
#'   estimation (added to the diagonal of \eqn{F_t}).
#' @param me_extra,shock_scale Per-period tunes used by the estimation. The
#'   thin filter does not implement them; any non-trivial value makes D41
#'   INFO rather than report innovations standardised by the wrong
#'   \eqn{F_t}.
#' @param lik_init  The estimation's \code{lik_init}. D41 always uses the
#'   stationary initialisation (exact for a stationary model); a model that
#'   needs the diffuse one is reported INFO.
#' @param likelihood The estimation's likelihood type; only
#'   \code{"gaussian"} is checked.
#' @param z_crit    Flag threshold on \eqn{|z|} (default 4).
#' @param meta      Optional \code{\link{diag_meta}} provenance descriptor.
#' @references Harvey, A. C. (1989). \emph{Forecasting, Structural Time Series
#'   Models and the Kalman Filter}, ch. 5. Cambridge University Press.
#'   Durbin, J. & Koopman, S. J. (2012). \emph{Time Series Analysis by State
#'   Space Methods}, 2nd ed., ch. 2.12 and 7.
#' @return A \code{dynhr_diagnostic}; \code{$result} is the
#'   \code{kf_innovation_diagnostics} object (\code{NULL} when not evaluated).
#' @noRd
d41_innovation_whiteness <- function(data, dr, model, params, obs_vars,
                                     compiled = NULL,
                                     me_variance = 0,
                                     me_extra = NULL,
                                     shock_scale = NULL,
                                     lik_init = "auto",
                                     likelihood = "gaussian",
                                     z_crit = 4,
                                     meta = NULL) {
  .info <- function(why)
    .make_result(pass = NA,
                 summary = paste0("D41 KF innovation whiteness: not evaluated -- ",
                                  why))

  if (!identical(likelihood %||% "gaussian", "gaussian"))
    return(.info(sprintf(paste0("the estimation likelihood is \"%s\"; the ",
                                "linear-Gaussian innovations do not describe it."),
                         likelihood)))
  if (!is.null(me_extra) && any(is.finite(me_extra) & me_extra != 0))
    return(.info(paste0("the estimation used a per-period me_extra, which the ",
                        "innovation filter does not implement.")))
  if (!is.null(shock_scale) && any(is.finite(shock_scale) & shock_scale != 1))
    return(.info(paste0("the estimation used a per-period shock_scale, which the ",
                        "innovation filter does not implement.")))
  if (!is.null(model$obs_aggregation))
    return(.info(paste0("mixed-frequency obs_aggregation is not implemented by ",
                        "the innovation filter.")))

  Y <- .d41_obs_matrix(data, obs_vars)

  ## Evaluation point: decision rule and params must describe the same theta.
  if (!is.null(compiled)) {
    sol <- .solve_dr_for_theta(model, compiled, cache_system_structure(compiled),
                               params, new.env(parent = emptyenv()),
                               lik_init = lik_init %||% "auto",
                               shock_scale = NULL)
    if (is.null(sol))
      return(.info(paste0("the model does not solve (steady state / BK / ",
                          "stationarity) at the evaluation parameters.")))
    dr <- sol$dr
    params <- sol$params
  } else {
    params <- .apply_theta_to_params(model, params)
  }

  ## Initialisation: mirror kalman_filter(lik_init = "auto") -- stationary
  ## unless a root sits at/over the unit circle so that the Lyapunov P0 is
  ## not a finite PSD covariance.
  TT <- dr$ghx[dr$state_idx, , drop = FALSE]
  if (nrow(TT) > 0L) {
    if (!all(is.finite(TT)))
      return(.info("the state transition matrix is not finite."))
    ev <- eigen(TT, symmetric = FALSE, only.values = TRUE)$values
    if (any(Mod(ev) > 1 - 1e-6)) {
      RR <- dr$ghu[dr$state_idx, , drop = FALSE]
      QQ <- tcrossprod(RR %*% .get_shock_cov(model, dr$exo_names, params), RR)
      P0 <- solve_lyapunov(TT, QQ)
      ## kalman_filter()'s own rule, RELATIVE to P0's scale: the
      ## absolute min eig > -1e-8 flipped with the units of the model.
      if (!.kf_stationary_P0_ok(P0))
        return(.info(sprintf(paste0(
          "the state transition has a unit/explosive root (max |lambda| = ",
          "%.6f); the stationary initialisation does not exist and the thin ",
          "filter has no diffuse phase."), max(Mod(ev)))))
    }
  }

  diag_obj <- kf_innovation_diagnostics(Y, dr = dr, model = model,
                                        params = params, obs_vars = obs_vars,
                                        lik_init = "stationary",
                                        me_variance = me_variance)

  ## Mean test. kf_innovation_diagnostics' variance is DEMEANED, so a level
  ## offset (stale steady state, wrong DD) is invisible to it and reaches its
  ## acf1 only as m^2/(1 + m^2). E[z] = 0 with sd 1/sqrt(T) under the model.
  by_obs <- diag_obj$by_obs
  zm <- diag_obj$z
  n_ok <- rowSums(is.finite(zm))
  by_obs$mean_z <- ifelse(n_ok >= 2L, rowSums(zm, na.rm = TRUE) / pmax(n_ok, 1L),
                          NA_real_)
  by_obs$z_mean <- by_obs$mean_z * sqrt(n_ok)

  ## Variance test, kurtosis-robust (0.9.4). kf_innovation_diagnostics()
  ## calibrates (s^2 - 1) with the GAUSSIAN sd sqrt(2/T). The general
  ## large-sample variance of a sample variance is
  ##   Var(s^2) = (mu_4 - sigma^4)/T = sigma^4 (kappa - 1)/T,
  ## which is 2/T only when kappa = 3. DSGE innovations are routinely
  ## fat-tailed (non-Gaussian structural shocks, occasionally-binding
  ## constraints, an imperfect linear approximation), and then the Gaussian sd
  ## is far too small and the test over-rejects: at kappa = 9 (Student t_5) the
  ## Gaussian z is inflated by sqrt((9-1)/2) = 2x. Recomputed here rather than
  ## in kf_innovation_diagnostics(), whose z_var keeps its documented
  ## Gaussian calibration.
  by_obs <- .d41_robust_var_z(diag_obj$z, by_obs)
  diag_obj$by_obs <- by_obs

  zs <- cbind(by_obs$z_mean, by_obs$z_var, by_obs$z_acf1)
  zmax_i <- apply(abs(zs), 1L, function(r)
    if (any(is.finite(r))) max(r, na.rm = TRUE) else NA_real_)
  flagged <- by_obs$obs_var[which(zmax_i > z_crit)]
  n_z <- sum(is.finite(zs))
  n_flag_z <- sum(abs(zs) > z_crit, na.rm = TRUE)
  max_abs_z <- if (n_z > 0L) max(abs(zs), na.rm = TRUE) else NA_real_
  diag_obj$d41 <- list(z_crit = z_crit, n_z = n_z, n_flagged = n_flag_z,
                       max_abs_z = max_abs_z, flagged = flagged)
  ## No finite statistic (every observable has < 2 usable periods) is not a PASS.
  pass <- if (n_z == 0L) NA else length(flagged) == 0L

  init_note <- if (!identical(lik_init %||% "auto", "auto") &&
                   !identical(lik_init, "stationary"))
    sprintf(" (estimation lik_init = \"%s\"; D41 uses the exact stationary init)",
            lik_init) else ""
  summary_txt <- sprintf(paste0(
    "D41 KF innovation whiteness: %d observable(s), %d z-stat(s), ",
    "max |z| = %.2f, %d flagged (|z| > %g)%s%s."),
    nrow(by_obs), n_z, max_abs_z, n_flag_z, z_crit,
    if (length(flagged)) paste0(" [", paste(flagged, collapse = ", "), "]") else "",
    init_note)

  badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
  llm <- paste(c(
    sprintf("D41 | KF Innovation Whiteness | %s", badge),
    sprintf("  n_obs=%d max_abs_z=%.3f n_flagged=%d z_crit=%g me_variance=%g",
            nrow(by_obs), max_abs_z, n_flag_z, z_crit, me_variance),
    sprintf("  per_obs: %s",
            paste(sprintf("%s T=%d mean=%.3f(z=%.2f) var=%.3f(z=%.2f) acf1=%.3f(z=%.2f)",
                          by_obs$obs_var, by_obs$n_used, by_obs$mean_z,
                          by_obs$z_mean, by_obs$var_z, by_obs$z_var,
                          by_obs$acf1, by_obs$z_acf1),
                  collapse = "; ")),
    if (length(flagged))
      sprintf("  flagged(|z|>%g): %s", z_crit, paste(flagged, collapse = ", ")),
    sprintf("  action: %s",
            if (isTRUE(pass))
              "Standardised innovations are consistent with unit-variance white noise at these parameters."
            else if (is.na(pass))
              "Too few non-missing periods to test."
            else paste0("Innovations deviate from zero-mean unit-variance white noise ",
                        "(mean != 0: level offset -- steady state / measurement ",
                        "constant; var_z > 1: the model under-predicts the forecast-error ",
                        "variance; acf1 != 0: unmodelled persistence). Check the ",
                        "evaluation point, Sigma_e, measurement equations, steady ",
                        "state and measurement error. z_var uses the kurtosis-robust sd ",
                        "sqrt((kappa - 1)/T), so fat tails alone should not flag it; ",
                        "variance DRIFT over the sample is not tested (no H(h) test)."))
  ), collapse = "\n")

  .make_result(
    result = diag_obj, pass = pass,
    plots = .d41_plots(diag_obj, z_crit, meta),
    summary = summary_txt, llm_summary = llm)
}


## Kurtosis-robust z-statistic for "sample variance of the standardised
## innovations equals 1", per observable (row of `z`).
##
##   z_var = (s_i^2 - 1) / sqrt((kappa_i - 1) / T_i),   kappa_i = m4 / m2^2
##
## with m2, m4 the 2nd and 4th central sample moments. kappa = 3 recovers the
## Gaussian sqrt(2/T). The estimated kurtosis is itself noisy in short
## samples, so the robust sd is used only from `n_min` usable periods
## (default 20 -- a PACKAGE CHOICE, no literature source) and kappa is floored
## just above its theoretical minimum of 1 so the sd stays positive. Below
## `n_min` the Gaussian sd is kept and `kurt` is still reported.
.d41_robust_var_z <- function(z, by_obs, n_min = 20L) {
  by_obs$kurt <- NA_real_
  for (i in seq_len(nrow(by_obs))) {
    zi <- z[i, ]
    zi <- zi[is.finite(zi)]
    n_i <- length(zi)
    if (n_i < 4L) next
    d  <- zi - mean(zi)
    m2 <- mean(d^2)
    m4 <- mean(d^4)
    if (!(is.finite(m2) && m2 > 0 && is.finite(m4))) next
    kappa <- m4 / m2^2
    by_obs$kurt[i] <- kappa
    if (n_i < n_min) next
    ## s^2 is the unbiased sample variance, matching kf_innovation_diagnostics.
    s2 <- stats::var(zi)
    by_obs$z_var[i] <- (s2 - 1) / sqrt(max(kappa - 1, 1e-8) / n_i)
  }
  by_obs
}

## Observation matrix in kf_innovation_diagnostics' n_obs x T orientation.
## Named columns (the orchestrator's T x n_obs convention) are selected by
## name, so extra columns and T == n_obs are both handled; unnamed input is
## resolved by dimension, T x n_obs first (the documented orchestrator form).
.d41_obs_matrix <- function(data, obs_vars) {
  Y <- as.matrix(data)
  n_obs <- length(obs_vars)
  if (!is.null(colnames(Y)) && all(obs_vars %in% colnames(Y)))
    return(t(Y[, obs_vars, drop = FALSE]))
  if (!is.null(rownames(Y)) && all(obs_vars %in% rownames(Y)))
    return(Y[obs_vars, , drop = FALSE])
  if (ncol(Y) == n_obs) return(t(Y))
  if (nrow(Y) == n_obs) return(Y)
  .dynhr_abort(sprintf(paste0(
    "D41: data (%d x %d) matches neither T x %d nor %d x T, and its dimnames ",
    "do not contain all of obs_vars."), nrow(Y), ncol(Y), n_obs, n_obs))
}

## Plots: (1) the z-statistics per observable against the flag threshold,
## (2) the standardised innovation series.
.d41_plots <- function(diag_obj, z_crit, meta) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  by_obs <- diag_obj$by_obs
  obs_f <- function(x) factor(x, levels = rev(by_obs$obs_var))
  stat_lv <- c("Mean (null 0)", "Variance (null 1)",
               "Lag-1 autocorrelation (null 0)")
  zdf <- data.frame(
    obs  = obs_f(rep(by_obs$obs_var, 3L)),
    stat = factor(rep(stat_lv, each = nrow(by_obs)), levels = stat_lv),
    z    = c(by_obs$z_mean, by_obs$z_var, by_obs$z_acf1))
  zdf <- zdf[is.finite(zdf$z), , drop = FALSE]
  lim <- max(z_crit + 1, abs(zdf$z), na.rm = TRUE)
  p_z <- ggplot2::ggplot(zdf, ggplot2::aes(x = z, y = obs, colour = stat,
                                           shape = stat)) +
    ggplot2::annotate("rect", xmin = -z_crit, xmax = z_crit,
                      ymin = -Inf, ymax = Inf, fill = dynhr_na_fill,
                      alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, colour = dynhr_na_colour,
                        linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-z_crit, z_crit), linetype = "dashed",
                        colour = dynhr_colours$red, linewidth = 0.4) +
    ggplot2::geom_point(size = 2.6,
                        position = ggplot2::position_dodge(width = 0.45)) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim)) +
    ggplot2::scale_colour_manual(values = c(dynhr_colours$teal,
                                            dynhr_colours$mid_blue,
                                            dynhr_colours$orange),
                                 drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(15, 16, 17), drop = FALSE) +
    theme_dynhr() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "D41: Kalman-filter innovation whiteness",
      subtitle = sprintf(paste0(
        "z-statistics of v/sqrt(F_ii): sqrt(T) x mean, ",
        "(var - 1)/sqrt((kappa - 1)/T), sqrt(T) x acf1; each N(0,1) under the model.\n",
        "Shaded/dashed = +/-%g flag threshold."), z_crit),
      x = "z-statistic", y = NULL, colour = NULL, shape = NULL)

  z <- diag_obj$z
  n_obs <- nrow(z); n_t <- ncol(z)
  long <- data.frame(Period = rep(seq_len(n_t), each = n_obs),
                     obs = factor(rep(by_obs$obs_var, n_t),
                                  levels = by_obs$obs_var),
                     z = as.numeric(z))
  p_s <- ggplot2::ggplot(long, ggplot2::aes(x = Period, y = z)) +
    geom_dynhr_zero() +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted",
                        colour = dynhr_na_colour, linewidth = 0.35) +
    ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.35,
                       na.rm = TRUE) +
    ggplot2::facet_wrap(~ obs, ncol = min(2L, n_obs)) +
    theme_dynhr_compact() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.75)),
                   axis.title.y = ggplot2::element_text(angle = 90, size = ggplot2::rel(0.85))) +
    ggplot2::labs(
      title = "D41: Standardised Kalman-filter innovations",
      subtitle = "v[i,t] / sqrt(F[ii,t]); iid N(0,1) under the model. Dotted = +/-2.",
      x = "Period", y = "Standardised innovation")

  list(zstats = .apply_meta(p_z, meta), innovations = .apply_meta(p_s, meta))
}
