## R/diag-post-d12-shocks.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D12 smoothed shock diagnostics: standardisation (auxiliary residuals),
## mean / variance / normality / outlier checks, informational Ljung-Box.
## --------------------------------------------------------------------------

#' Per-period standard deviation of the smoothed shocks under the model
#'
#' Returns \eqn{sd(\hat\varepsilon_{t|T})}, the square root of the diagonal of
#' \eqn{Var(\hat\varepsilon_{t|T}) = \Sigma_t - Var(\varepsilon_t \mid Y_{1:T})}
#' (law of total variance; Durbin & Koopman 2012, sec. 7.5). Dividing the
#' smoothed shocks by it gives the auxiliary residuals, which are marginally
#' N(0, 1) in every period when the model is correct. Dividing by the shock
#' std \eqn{\sigma_j} instead does NOT: the smoother shrinks, so
#' \eqn{Var(\hat\varepsilon_{t|T}) < \sigma_j^2}, badly so when shocks
#' outnumber observables.
#'
#' \eqn{Var(\varepsilon_t \mid Y)} is read off the ordinary smoother run on a
#' state space augmented with \eqn{\varepsilon_t} as extra states
#' (\eqn{\tilde s_t = (s_t, \varepsilon_t)}, \eqn{\tilde T = diag(T, 0)},
#' \eqn{\tilde R = (R', I)'}, \eqn{\tilde Z = (Z, 0)}), so no new recursion is
#' needed.
#'
#' @param data  T x n_obs data in LEVELS (as for \code{.kalman_smoother_ss}).
#' @param ss    A \code{dsge_ss} (lagged timing).
#' @param d,me_variance,me_extra,shock_scale,lik_init  As for
#'   \code{.kalman_smoother_ss}; pass the SAME values as the smoother run whose
#'   shocks are being standardised.
#' @return T x n_shock matrix (column names = shock names).
#' @noRd
.d12_smoothed_shock_sd <- function(data, ss, d = NULL, me_variance = 0,
                                   me_extra = NULL, shock_scale = NULL,
                                   lik_init = "auto") {
  ns <- ss$n_state; ne <- ss$n_shock; no <- ss$n_obs
  if (is.null(ss$Sigma_e))
    .dynhr_abort("D12: the state space has no Sigma_e; cannot standardise ",
                 "the smoothed shocks.")
  Ta <- matrix(0, ns + ne, ns + ne)
  Ta[seq_len(ns), seq_len(ns)] <- ss$T_mat
  sa <- list(
    T_mat = Ta,
    R_mat = rbind(ss$R_mat, diag(ne)),
    Z_mat = cbind(ss$Z_mat, matrix(0, no, ne)),
    D_mat = ss$D_mat, Sigma_e = ss$Sigma_e,
    n_state = ns + ne, n_shock = ne, n_obs = no,
    state_names = c(ss$state_names, paste0(".d12_eps_", ss$shock_names)),
    shock_names = ss$shock_names)
  ## A hand-built list carries no intercept, so pass the original one.
  d_use <- if (!is.null(d)) d else if (!is.null(ss$d)) ss$d else 0
  sm <- .kalman_smoother_ss(data, sa, d = d_use, me_variance = me_variance,
                            me_extra = me_extra, shock_scale = shock_scale,
                            lik_init = lik_init)
  idx <- ns + seq_len(ne)
  TT  <- dim(sm$smoothed_cov)[3]
  V   <- matrix(vapply(seq_len(TT), function(t)
                  diag(as.matrix(sm$smoothed_cov[idx, idx, t])), numeric(ne)),
                nrow = ne)                                   # ne x T
  prior_var <- matrix(diag(as.matrix(ss$Sigma_e)), ne, TT)   # ne x T
  if (!is.null(shock_scale)) prior_var <- prior_var * as.matrix(shock_scale)^2
  out <- t(sqrt(pmax(prior_var - V, 0)))                     # T x ne
  colnames(out) <- ss$shock_names
  out
}

## Newey-West (Bartlett) long-run variance of x (NA-free).
.d12_lrv <- function(x, lag) {
  n <- length(x); x <- x - mean(x)
  g0 <- sum(x^2) / n
  L <- min(lag, n - 1L)
  s <- g0
  for (k in seq_len(L)) s <- s + 2 * (1 - k / (L + 1)) * sum(x[-seq_len(k)] * x[seq_len(n - k)]) / n
  if (!is.finite(s) || s <= 0) g0 else s
}

## Jarque-Bera statistic and chi2(2) p-value (Jarque & Bera 1987).
.d12_jarque_bera <- function(x) {
  n <- length(x); x <- x - mean(x)
  m2 <- mean(x^2)
  if (n < 3L || m2 <= 0) return(c(statistic = NA_real_, p_value = NA_real_))
  S <- mean(x^3) / m2^1.5
  K <- mean(x^4) / m2^2
  jb <- n / 6 * (S^2 + (K - 3)^2 / 4)
  c(statistic = jb, p_value = stats::pchisq(jb, df = 2, lower.tail = FALSE))
}

#' D12. Smoothed shock diagnostics
#'
#' Checks the smoothed structural shocks \eqn{\hat\varepsilon_{t|T}} after
#' standardising them. The standardisation decides what the checks mean:
#' \describe{
#'   \item{\code{"auxiliary"}}{\code{shock_sd} is a T x n matrix of
#'     \eqn{sd(\hat\varepsilon_{t|T})} (see \code{.d12_smoothed_shock_sd};
#'     the posterior auto-dispatch path supplies it). The standardised series
#'     \eqn{z_t} is marginally N(0, 1) under a correct model, so mean, variance,
#'     normality and outliers are all testable.}
#'   \item{\code{"sigma"}}{\code{shock_sd} is a length-n vector of shock
#'     stds. Smoothed shocks are shrunk, so \eqn{Var(z) \le 1} under a correct
#'     model and only an EXCESS variance is tested.}
#'   \item{\code{"sample"}}{no \code{shock_sd}: each shock is divided by its
#'     own sample sd, so the variance cannot be tested (reported NA).}
#' }
#' Tests per shock (on \eqn{z}): mean = 0 (Newey-West t-test, Bartlett lag
#' \code{max_lag}); variance = 1 (Newey-West t-test on \eqn{z^2 - 1});
#' normality (Jarque-Bera, chi-squared with 2 df); outliers (count of
#' \eqn{|z| >} \code{outlier_threshold} against Binomial(T, 2 Phi(-thr))).
#' The badge FAILs when any p-value falls below \code{alpha / m} (Bonferroni
#' over the m tests actually run).
#'
#' Ljung-Box (df = \code{max_lag}) is reported as INFO but NOT gated.
#' Standardised smoothed disturbances are Harvey & Koopman's (1992)
#' \emph{auxiliary residuals}, and that paper's central caveat is that they
#' are SERIALLY CORRELATED even under a correctly specified model with known
#' parameters -- the smoother's law-of-total-variance shrinkage
#' (\eqn{Var(\hat\epsilon_{t|T}) = \Sigma_e - Var(\epsilon_t|Y)}, Durbin &
#' Koopman 2012, sec. 4.5/7.5) induces the correlation. The chi-squared
#' calibration therefore does not hold: in a 3-shock/1-observable toy model
#' the test rejects 64-98 percent of the time at a nominal 5 percent. This is
#' a structural property, not a finite-sample size distortion, so a Monte
#' Carlo recalibration of the null would not fix it. Harvey & Koopman
#' recommend auxiliary residuals for NORMALITY and OUTLIER/break detection
#' (what the badge gates) and moving serial-correlation testing to the
#' one-step prediction errors -- which is exactly what D41
#' (\code{kf_innovation_diagnostics}) does.
#'
#' @param shocks  T x n matrix of smoothed shocks, or a list with
#'   \code{smoothed_shocks} and optionally \code{shock_sd}.
#' @param shock_names  Character vector (optional; default column names).
#' @param max_lag  Ljung-Box / ACF lag and Newey-West bandwidth (default 8).
#' @param meta  Plot metadata.
#' @param shock_sd  NULL, a length-n vector (shock stds) or a T x n matrix
#'   (per-period sd of the smoothed shock).
#' @param dates  Optional length-T vector for the time axis.
#' @param outlier_threshold  |z| cut-off for outliers (default 3).
#' @param alpha  Family-wise level for the badge (default 0.05).
#' @return dynhr_diagnostic list
#'
#' @references Durbin, J. & Koopman, S. J. (2012). \emph{Time Series Analysis
#'   by State Space Methods}, 2nd ed., sec. 4.5 and 7.5.
#'   Harvey, A. C. & Koopman, S. J. (1992). Diagnostic checking of unobserved-
#'   components time series models. \emph{Journal of Business & Economic
#'   Statistics}, 10(4), 377-389.
#'   Jarque, C. M. & Bera, A. K. (1987). A test for normality of observations
#'   and regression residuals. \emph{Int. Stat. Rev.} 55(2), 163-172.
#'   Ljung, G. M., & Box, G. E. P. (1978). On a measure of lack of fit in
#'   time series models. \emph{Biometrika}, 65(2), 297-303.
#' @noRd
d12_smoothed_shocks <- function(shocks,
                                shock_names = NULL,
                                max_lag = 8,
                                meta = NULL,
                                shock_sd = NULL,
                                dates = NULL,
                                outlier_threshold = 3,
                                alpha = 0.05) {

  if (is.list(shocks) && !is.data.frame(shocks)) {
    if (is.null(shock_sd)) shock_sd <- shocks$shock_sd
    shocks <- shocks$smoothed_shocks
    if (is.null(shocks))
      .dynhr_abort("D12: `shocks` list has no `smoothed_shocks` element.")
  }
  shocks <- as.matrix(shocks)
  n_t <- nrow(shocks)
  n_s <- ncol(shocks)
  if (is.null(shock_names)) {
    shock_names <- if (!is.null(colnames(shocks))) colnames(shocks)
                   else paste0("eps_", seq_len(n_s))
  }
  if (length(shock_names) != n_s)
    .dynhr_abort(sprintf("D12: %d shock names for %d shock columns.",
                         length(shock_names), n_s))

  ## ---- Standardisation --------------------------------------------------
  if (is.null(shock_sd)) {
    std_mode <- "sample"
    sd_mat <- matrix(apply(shocks, 2, stats::sd, na.rm = TRUE), n_t, n_s,
                     byrow = TRUE)
  } else if (is.matrix(shock_sd) && nrow(shock_sd) > 1L) {
    std_mode <- "auxiliary"
    if (!identical(dim(shock_sd), dim(shocks)))
      .dynhr_abort(sprintf("D12: `shock_sd` is %d x %d but `shocks` is %d x %d.",
                           nrow(shock_sd), ncol(shock_sd), n_t, n_s))
    if (!is.null(colnames(shock_sd)) && !is.null(colnames(shocks)))
      shock_sd <- shock_sd[, colnames(shocks), drop = FALSE]
    sd_mat <- shock_sd
  } else {
    std_mode <- "sigma"
    shock_sd <- stats::setNames(as.numeric(shock_sd), names(shock_sd))
    if (!is.null(names(shock_sd)) && !is.null(colnames(shocks)) &&
        all(colnames(shocks) %in% names(shock_sd)))
      shock_sd <- shock_sd[colnames(shocks)]
    if (length(shock_sd) != n_s)
      .dynhr_abort(sprintf("D12: `shock_sd` has %d entries for %d shocks.",
                           length(shock_sd), n_s))
    sd_mat <- matrix(shock_sd, n_t, n_s, byrow = TRUE)   # row t = shock sds
  }
  ## A (near-)zero sd means the shock is determined in that period; there is
  ## nothing to standardise, so the period is dropped for that shock.
  scale_ref <- max(abs(sd_mat), na.rm = TRUE)
  tiny <- !is.finite(sd_mat) | sd_mat <= 1e-10 * max(scale_ref, 1e-300)
  z <- shocks / sd_mat
  z[tiny] <- NA_real_
  colnames(z) <- shock_names

  ## ---- Per-shock statistics ----------------------------------------------
  p_out <- 2 * stats::pnorm(-outlier_threshold)
  per <- lapply(seq_len(n_s), function(j) {
    zj <- z[, j]; ok <- is.finite(zj); x <- zj[ok]; n <- length(x)
    res <- list(n = n, mean = NA_real_, var = NA_real_, mean_p = NA_real_,
                var_p = NA_real_, jb = NA_real_, jb_p = NA_real_,
                n_out = NA_integer_, out_p = NA_real_, out_idx = integer(0),
                lb = list(statistic = NA_real_, p_value = NA_real_, pass = NA))
    if (n < max(10L, max_lag + 2L)) return(res)
    res$mean <- mean(x)
    res$var  <- mean(x^2)          # E[z^2] = 1 is the claim (mean tested separately)
    res$mean_p <- 2 * stats::pnorm(-abs(res$mean / sqrt(.d12_lrv(x, max_lag) / n)))
    if (std_mode != "sample") {
      tv <- (res$var - 1) / sqrt(.d12_lrv(x^2 - 1, max_lag) / n)
      res$var_p <- if (std_mode == "auxiliary") 2 * stats::pnorm(-abs(tv))
                   else stats::pnorm(tv, lower.tail = FALSE)   # excess only
    }
    jb <- .d12_jarque_bera(x)
    res$jb <- unname(jb["statistic"]); res$jb_p <- unname(jb["p_value"])
    res$out_idx <- which(ok)[abs(x) > outlier_threshold]
    res$n_out <- length(res$out_idx)
    res$out_p <- stats::pbinom(res$n_out - 1L, n, p_out, lower.tail = FALSE)
    res$lb <- .ljung_box(x, max_lag = max_lag)
    res
  })
  names(per) <- shock_names
  lb_results <- lapply(per, function(r) r$lb)

  num <- function(f) vapply(per, function(r) as.numeric(r[[f]]), numeric(1))
  shock_stats <- data.frame(
    Shock      = shock_names,
    N          = num("n"),
    Mean       = num("mean"),
    Variance   = num("var"),
    Mean_pval  = num("mean_p"),
    Var_pval   = num("var_p"),
    JB_Stat    = num("jb"),
    JB_pval    = num("jb_p"),
    N_outliers = num("n_out"),
    Outlier_pval = num("out_p"),
    LB_Stat    = vapply(lb_results, function(r) as.numeric(r$statistic), numeric(1)),
    LB_pval    = vapply(lb_results, function(r) as.numeric(r$p_value), numeric(1)),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  pcols <- c("Mean_pval", "Var_pval", "JB_pval", "Outlier_pval")
  pm <- as.matrix(shock_stats[, pcols])
  m_tests <- sum(is.finite(pm))
  crit <- if (m_tests > 0) alpha / m_tests else NA_real_
  flag <- is.finite(pm) & pm < crit
  colnames(flag) <- c("mean", "variance", "normality", "outliers")
  shock_stats$Flags <- apply(flag, 1, function(f) paste(names(f)[f], collapse = ","))
  all_pass <- if (m_tests == 0) NA else !any(flag)

  outliers <- do.call(rbind, lapply(seq_len(n_s), function(j) {
    ix <- per[[j]]$out_idx
    if (!length(ix)) return(NULL)
    data.frame(Shock = shock_names[j], Period = ix,
               Date = if (length(dates) == n_t) as.character(dates[ix]) else as.character(ix),
               z = z[ix, j], row.names = NULL, stringsAsFactors = FALSE)
  }))

  std_label <- switch(std_mode,
    auxiliary = "divided by the model sd of the smoothed shock (auxiliary residuals; N(0,1) under the model)",
    sigma     = "divided by the shock std (smoothed shocks are shrunk: variance <= 1 under the model)",
    sample    = "divided by the sample sd (variance not testable)")

  # --- Plots ---
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    time_x <- if (length(dates) == n_t) dates else seq_len(n_t)
    shock_f <- function(x) factor(x, levels = shock_names)
    ## The compact theme hides the y axis; these panels need the z scale.
    th <- theme_dynhr_compact() +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_text(size = ggplot2::rel(0.75)),
        axis.title.y = ggplot2::element_text(angle = 90, size = ggplot2::rel(0.85)),
        panel.spacing.x = ggplot2::unit(1, "lines"))

    acf_df <- do.call(rbind, lapply(seq_len(n_s), function(j) {
      zj <- z[, j]
      if (sum(is.finite(zj)) < max_lag + 2L) return(NULL)
      a <- stats::acf(zj, lag.max = max_lag, plot = FALSE, na.action = stats::na.pass)
      data.frame(Lag = as.numeric(a$lag[-1]), ACF = as.numeric(a$acf[-1]),
                 Shock = shock_names[j])
    }))
    if (!is.null(acf_df)) {
      acf_df$Shock <- shock_f(acf_df$Shock)
      ci <- stats::qnorm(0.975) / sqrt(n_t)
      p_acf <- ggplot2::ggplot(acf_df, ggplot2::aes(x = Lag, y = ACF)) +
        geom_dynhr_zero() +
        ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
                            colour = dynhr_colours$red, linewidth = 0.4) +
        ggplot2::geom_col(fill = dynhr_colours$mid_blue, width = 0.5) +
        ggplot2::scale_x_continuous(breaks = seq_len(max_lag)) +
        ggplot2::facet_wrap(~ Shock, ncol = min(3L, n_s)) +
        th +
        ggplot2::labs(
          title = "D12: ACF of standardised smoothed shocks (informational)",
          subtitle = sprintf(paste0("Dashed = +/-%.3f, the band for an iid series. ",
                                    "Smoothed shocks are autocorrelated even under a ",
                                    "correct model;\ntest whiteness on the innovations."), ci),
          x = "Lag (periods)", y = "Autocorrelation")
      plots$acf <- .apply_meta(p_acf, meta)
    }

    long <- data.frame(Time = rep(time_x, n_s),
                       z = as.numeric(z),
                       Shock = shock_f(rep(shock_names, each = n_t)))
    long$Outlier <- is.finite(long$z) & abs(long$z) > outlier_threshold
    p_ss <- ggplot2::ggplot(long, ggplot2::aes(x = Time, y = z)) +
      geom_dynhr_zero() +
      ggplot2::geom_hline(yintercept = c(-outlier_threshold, outlier_threshold),
                          linetype = "dashed", colour = dynhr_colours$red,
                          linewidth = 0.4) +
      ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.35,
                         na.rm = TRUE) +
      ggplot2::geom_point(data = long[long$Outlier, , drop = FALSE],
                          colour = dynhr_colours$red, size = 1.4) +
      ggplot2::facet_wrap(~ Shock, ncol = min(3L, n_s)) +
      th +
      ggplot2::labs(title = "D12: Standardised smoothed shocks",
                    subtitle = sprintf("Each shock %s.\nDashed = +/-%g outlier threshold; red points = outliers.",
                                       std_label, outlier_threshold),
                    x = if (length(dates) == n_t) "Date" else "Period",
                    y = "Standardised shock (z)")
    plots$shock_series <- .apply_meta(p_ss, meta)

    qq <- do.call(rbind, lapply(seq_len(n_s), function(j) {
      x <- sort(z[is.finite(z[, j]), j]); n <- length(x)
      if (n < 3L) return(NULL)
      data.frame(Theoretical = stats::qnorm(stats::ppoints(n)), Sample = x,
                 Shock = shock_names[j])
    }))
    if (!is.null(qq)) {
      qq$Shock <- shock_f(qq$Shock)
      p_qq <- ggplot2::ggplot(qq, ggplot2::aes(x = Theoretical, y = Sample)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = dynhr_colours$red, linewidth = 0.4) +
        ggplot2::geom_point(colour = dynhr_colours$mid_blue, size = 0.9) +
        ggplot2::facet_wrap(~ Shock, ncol = min(3L, n_s)) +
        th +
        ggplot2::labs(title = "D12: Normal QQ plot of standardised smoothed shocks",
                      subtitle = "Dashed = N(0,1) reference (45-degree line).",
                      x = "N(0,1) quantile", y = "Standardised shock quantile")
      plots$qq <- .apply_meta(p_qq, meta)
    }
  }

  flagged <- shock_stats$Shock[nzchar(shock_stats$Flags)]
  flag_txt <- paste(sprintf("%s(%s)", flagged,
                            shock_stats$Flags[nzchar(shock_stats$Flags)]),
                    collapse = ", ")
  summary_txt <- sprintf(
    "D12 Smoothed shocks: %d shocks, %d periods, standardisation=%s. %s",
    n_s, n_t, std_mode,
    if (is.na(all_pass)) "Too few periods to test."
    else if (all_pass) sprintf("No shock rejects mean/variance/normality/outlier checks (Bonferroni, %d tests, alpha=%g).",
                               m_tests, alpha)
    else sprintf("Flagged: %s.", flag_txt))

  .make_result(
    result  = list(shock_stats = shock_stats, lb_results = lb_results,
                   standardised = z, standardisation = std_mode,
                   outliers = outliers, alpha_per_test = crit),
    pass    = all_pass,
    plots   = plots,
    summary = summary_txt,
    llm_summary = paste(c(
      sprintf("D12 | Smoothed Shocks | %s",
              if (is.na(all_pass)) "N/A" else if (all_pass) "PASS" else "FAIL"),
      sprintf("  shocks=%d periods=%d standardisation=%s tests=%d alpha_per_test=%.4g",
              n_s, n_t, std_mode, m_tests, crit),
      sprintf("  %s: mean=%.3f var=%.3f jb_p=%.3f outliers=%s lb_p=%.3f(info)",
              shock_stats$Shock, shock_stats$Mean, shock_stats$Variance,
              shock_stats$JB_pval, format(shock_stats$N_outliers),
              shock_stats$LB_pval),
      if (length(flagged)) sprintf("  flagged: %s", flag_txt),
      sprintf("  action: %s",
              if (isTRUE(all_pass))
                "Standardised smoothed shocks are consistent with the model's N(0,1) claim."
              else if (is.na(all_pass)) "Sample too short for the checks."
              else paste0("Inspect flagged shocks: a non-zero mean suggests a missing trend/level ",
                          "term; excess variance or outliers suggest breaks or mis-sized Sigma_e; ",
                          "non-normality suggests fat-tailed shocks."))
    ), collapse = "\n")
  )
}
