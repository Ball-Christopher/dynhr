## R/diag-post-d9-moments.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D9 moment matching; .get_moments_list(), .compute_empirical_moments() helpers
## --------------------------------------------------------------------------

## Normalise any supported moment object to
##   list(std_dev = named sd vector, acf_y = n x K matrix (rows named, cols
##        "lag1".."lagK") or NULL, correlation = n x n matrix or NULL)
##
## Accepted inputs (the first that applies wins):
##   * an object carrying `$moments` (a StochSimulResult / decision rules);
##   * compute_moments() output: `$std_dev`, `$autocorr` (n x n x K array of
##     Corr(y_{t+k}, y_t); only the diagonal is used), `$correlation`;
##   * an already-normalised list: `$std_dev` or `$sigma_y`, `$acf_y`,
##     `$correlation` / `$xcorr_y`.
## `sigma_y` is accepted either as a named SD VECTOR or as a COVARIANCE MATRIX
## (the posterior auto-dispatch path builds data_moments as
## list(sigma_y = cov(data), mean_y = ...)); a matrix is converted with
## sqrt(diag()) -- its diagonal is a variance, not an SD.
.get_moments_list <- function(moments) {
  if (!is.list(moments)) return(NULL)
  if (is.list(moments$moments)) moments <- moments$moments

  sd_src <- if (!is.null(moments$std_dev)) moments$std_dev else moments$sigma_y
  cov_mat <- moments$var_cov
  if (is.matrix(sd_src)) {
    cov_mat <- sd_src
    nm <- rownames(sd_src)
    if (is.null(nm)) nm <- colnames(sd_src)
    sd_src <- stats::setNames(sqrt(pmax(diag(sd_src), 0)), nm)
  }
  if (is.null(sd_src) || is.null(names(sd_src))) return(NULL)
  sigma_y <- stats::setNames(as.numeric(sd_src), names(sd_src))

  acf_mat <- NULL
  ac <- moments$autocorr
  if (is.array(ac) && length(dim(ac)) == 3L) {
    n_lag   <- dim(ac)[3]
    acf_mat <- matrix(NA_real_, dim(ac)[1], n_lag)
    for (k in seq_len(n_lag)) acf_mat[, k] <- diag(ac[, , k])
    rownames(acf_mat) <- if (!is.null(dimnames(ac)[[1]])) dimnames(ac)[[1]] else names(sigma_y)
    colnames(acf_mat) <- paste0("lag", seq_len(n_lag))
  } else if (is.matrix(moments$acf_y)) {
    acf_mat <- moments$acf_y
  } else if (is.matrix(ac)) {
    acf_mat <- ac
  }
  if (!is.null(acf_mat) && (is.null(rownames(acf_mat)) || ncol(acf_mat) == 0L))
    acf_mat <- NULL
  if (!is.null(acf_mat)) colnames(acf_mat) <- paste0("lag", seq_len(ncol(acf_mat)))

  xcorr <- moments$correlation
  if (is.null(xcorr)) xcorr <- moments$xcorr_y
  if (is.null(xcorr) && is.matrix(cov_mat) && !is.null(rownames(cov_mat))) {
    d <- sqrt(pmax(diag(cov_mat), 0))
    xcorr <- cov_mat / outer(d, d)
  }
  if (!is.matrix(xcorr) || is.null(rownames(xcorr))) xcorr <- NULL

  list(std_dev = sigma_y, sigma_y = sigma_y,
       acf_y = acf_mat, autocorr = acf_mat,
       correlation = xcorr, xcorr_y = xcorr,
       var_cov = cov_mat)
}


## Sample moments of a T x n data matrix.
##   std_dev     : sd() (divisor T-1), NA-removed per column.
##   correlation : pairwise-complete Pearson correlation.
##   acf_y       : stats::acf (divisor T, overall mean) with na.pass, so an
##                 internal gap does NOT splice the series together (the old
##                 code dropped NAs first, creating spurious adjacent pairs).
##   n_eff       : number of non-missing observations per column.
.compute_empirical_moments <- function(data, acf_lags = 5) {
  if (!is.matrix(data)) data <- as.matrix(data)
  vnames <- colnames(data)
  n_var  <- ncol(data)
  sds    <- apply(data, 2, stats::sd, na.rm = TRUE)
  names(sds) <- vnames
  cor_mat <- if (n_var > 1L) {
    suppressWarnings(stats::cor(data, use = "pairwise.complete.obs"))
  } else {
    matrix(1, 1, 1, dimnames = list(vnames, vnames))
  }
  acf_mat <- matrix(NA_real_, n_var, acf_lags,
                    dimnames = list(vnames, paste0("lag", seq_len(acf_lags))))
  n_ok <- colSums(!is.na(data))
  for (i in seq_len(n_var)) {
    x <- data[, i]
    if (n_ok[i] > acf_lags + 1 && is.finite(sds[i]) && sds[i] > 0) {
      ac <- stats::acf(x, lag.max = acf_lags, plot = FALSE,
                       na.action = stats::na.pass, demean = TRUE)$acf
      acf_mat[i, ] <- ac[2:(acf_lags + 1)]
    }
  }
  list(std_dev = sds, correlation = cor_mat, acf_y = acf_mat,
       n_eff = stats::setNames(as.integer(n_ok), vnames))
}


#' D9. Moment matching
#'
#' Compares model-implied second moments of the observables with their sample
#' counterparts: standard deviations (the badge), first-order-to-K
#' autocorrelations and contemporaneous cross-correlations (informational).
#'
#' Badge: at least 70\% of the observables with a finite model/data SD ratio
#' have the ratio in \eqn{[0.5, 2]} (within a factor of 2, bounds inclusive).
#' Ratios in \eqn{[0.5, 0.6)} or \eqn{(1.7, 2]} are flagged marginal.  Ratios
#' that are not finite (e.g. an observable loading on a unit root, whose model
#' SD is \code{NaN}) are listed and excluded from the denominator.
#' \code{pass = NA} when no data are supplied (self-comparison), or when no
#' observable has a finite ratio.  The factor-of-2 rule is a coarse
#' plausibility screen, not a test: nothing in the badge accounts for sampling
#' uncertainty.  It is the same naive ratio comparison Dynare users make by
#' hand off \code{moments_varendo}; no standard DSGE package ships a tolerance
#' band for it.
#'
#' \strong{Measurement error (\code{me_var}).} The data SD is the SD of what
#' was actually observed, \eqn{y_t = y^*_t + u_t}, so the like-for-like model
#' quantity is \eqn{\sqrt{Var(y^*) + Var(u)}}.  \code{compute_moments()}
#' returns the state-space SD only and knows nothing about the estimation's
#' measurement error, so until 0.9.4 every ratio was biased DOWN -- most for
#' the observables carrying the most ME.  Supply the estimation's ME variances
#' and they are added to the model variances before the ratio is formed; the
#' posterior path does this automatically.  \code{sd_table} reports both
#' \code{model_sd} (with ME) and \code{model_sd_state} (without).
#'
#' \strong{Posterior-predictive band (INFO).} When posterior draws and a
#' simulation closure are available, D9 additionally reports, per observable, a
#' \code{ppc_level} interval of the SD of a T-length path simulated from the
#' model at draws from the posterior, and whether the observed data SD falls
#' inside it.  That interval carries BOTH parameter uncertainty (across draws)
#' and finite-sample sampling variation (within a draw), which a point
#' comparison against a fixed threshold ignores.  This follows Faust & Gupta
#' (2012) and Herbst & Schorfheide (2016, ch. 6): draw parameters, simulate,
#' compute the statistic, compare its DISTRIBUTION to the single realised
#' value.  It is deliberately NOT an analytic Bartlett/delta-method band --
#' macro series are short and highly autocorrelated, exactly where asymptotic
#' SEs for SDs and ACFs are least reliable, and the literature does not use
#' them for this purpose.  The band is INFO: it is reported and plotted but
#' never moves the badge.
#'
#' @param model_moments  Model-implied moments: \code{compute_moments()} output,
#'   an object with \code{$moments}, or a list with \code{$std_dev} /
#'   \code{$sigma_y} (SD vector), optional \code{$acf_y} (n x K matrix) and
#'   \code{$correlation}.
#' @param data           T x n_obs data matrix (columns named).  When supplied
#'   and \code{data_moments} is NULL, sample moments are computed from it.
#' @param data_moments   Pre-computed sample moments: a list with
#'   \code{$std_dev}, or \code{$sigma_y} as an SD vector or a covariance
#'   matrix.  Takes precedence over \code{data}.
#' @param obs_names      Observables to compare (default: the names common to
#'   model and data).
#' @param metadata       Unused (kept for call compatibility).
#' @param acf_lags       Autocorrelation lags computed from \code{data}.
#' @param key_pairs      Optional character vector \code{"a:b"} restricting the
#'   cross-correlation table to these pairs.
#' @param me_var         Measurement-error VARIANCES to add to the model
#'   variances: NULL (none), a scalar (same variance on every observable, the
#'   \code{me_variance} the estimation uses), a named vector (matched by name;
#'   unnamed observables get 0), an unnamed vector of length
#'   \code{length(obs_names)} (positional), or a square ME covariance matrix
#'   (its diagonal is used).
#' @param ppc_sd_fn      Optional function \code{theta -> named numeric} of
#'   SAMPLE standard deviations computed from a T-length path simulated at
#'   \code{theta}.  Enables the posterior-predictive band.  Follows the
#'   package callback convention: return a non-finite / zero-length value for
#'   an infeasible draw (that draw is dropped), never \code{stop()}.
#' @param ppc_draws      Posterior draws (\code{n_draws x n_params}, named
#'   columns) fed to \code{ppc_sd_fn}.  At least 10 rows.
#' @param ppc_n          Number of draws used for the band (default 200, taken
#'   as an evenly spaced thinning of \code{ppc_draws}).  A PACKAGE CHOICE
#'   trading Monte Carlo error against one model solve + simulation per draw.
#' @param ppc_level      Band coverage (default 0.90).
#' @param meta           Optional \code{diag_meta()} object for plot captions.
#'
#' @return A \code{dynhr_diagnostic} list with \code{result$sd_table}
#'   (variable, model_sd, model_sd_state, me_var, data_sd, ratio, status, and
#'   ppc_lo/ppc_hi/ppc_in when the band ran), \code{result$acf_table}
#'   (variable, lag, model, data, diff), \code{result$corr_table} (pair,
#'   model, data, diff), \code{self_compare}, \code{marginal}, \code{outside},
#'   \code{nonfinite}; \code{pass}; \code{plots} (\code{moments_sd} ratio dot
#'   plot, \code{moments_acf}, \code{moments_corr} when available).
#'
#' @references Schorfheide, F. (2000). Loss function-based evaluation of DSGE models.
#'   \emph{Journal of Applied Econometrics}, 15(6), 645-670.
#'   Faust, J. & Gupta, A. (2012). Posterior predictive analysis for evaluating
#'   DSGE models. NBER Working Paper 17906.
#'   Herbst, E. & Schorfheide, F. (2016). \emph{Bayesian Estimation of DSGE
#'   Models}, ch. 6. Princeton University Press.
#' @noRd
d9_moment_matching <- function(model_moments, data = NULL, data_moments = NULL,
                               obs_names = NULL, metadata = NULL,
                               acf_lags = 5, key_pairs = NULL,
                               me_var = NULL,
                               ppc_sd_fn = NULL, ppc_draws = NULL,
                               ppc_n = 200L, ppc_level = 0.9,
                               meta = NULL) {
  lo <- 0.5; hi <- 2.0; m_lo <- 0.6; m_hi <- 1.7

  mm <- .get_moments_list(model_moments)
  if (is.null(mm)) {
    .dynhr_abort(paste0("d9_moment_matching(): `model_moments` has no named ",
                        "standard deviations ($std_dev, $sigma_y or $moments)."))
  }

  self_compare <- FALSE
  dm <- NULL
  if (!is.null(data_moments)) {
    dm <- .get_moments_list(data_moments)
    if (is.null(dm)) {
      .dynhr_abort(paste0("d9_moment_matching(): `data_moments` has no named ",
                          "standard deviations ($std_dev or $sigma_y)."))
    }
  } else if (!is.null(data)) {
    if (!is.matrix(data)) data <- as.matrix(data)
    dm <- .compute_empirical_moments(data, acf_lags)
  }
  if (is.null(dm)) { dm <- mm; self_compare <- TRUE }

  common <- intersect(names(mm$std_dev), names(dm$std_dev))
  obs_names <- if (!is.null(obs_names)) intersect(obs_names, common) else common
  n_obs <- length(obs_names)

  ## ---- standard deviations ------------------------------------------------
  m_sd_model <- as.numeric(mm$std_dev[obs_names])
  ## Measurement error belongs in the MODEL SD (0.9.4). The data SD is the SD
  ## of what was actually observed, y_t = y*_t + u_t, so the like-for-like
  ## model quantity is sqrt(Var(y*) + Var(u)) -- compute_moments() returns the
  ## STATE-SPACE SD only and knows nothing about the estimation's measurement
  ## error. Omitting it biased every ratio DOWN, the more so for the
  ## observables carrying the most ME.
  me_v <- .d9_align_me_var(me_var, obs_names)
  m_sd <- sqrt(m_sd_model^2 + me_v)
  d_sd <- as.numeric(dm$std_dev[obs_names])
  ratio <- m_sd / pmax(d_sd, 1e-12)
  ok <- is.finite(ratio)
  in_band  <- ok & ratio >= lo & ratio <= hi
  marginal <- in_band & (ratio < m_lo | ratio > m_hi)
  outside  <- ok & !in_band
  status <- ifelse(!ok, "non-finite",
            ifelse(outside, "outside", ifelse(marginal, "marginal", "ok")))
  sd_table <- data.frame(variable = obs_names, model_sd = m_sd,
                         model_sd_state = m_sd_model, me_var = me_v,
                         data_sd = d_sd,
                         ratio = ratio, status = status,
                         stringsAsFactors = FALSE)

  ## ---- posterior-predictive band (INFO only) ------------------------------
  ppc <- .d9_ppc_band(ppc_sd_fn, ppc_draws, obs_names, d_sd, ppc_n, ppc_level)
  if (!is.null(ppc)) {
    sd_table$ppc_lo <- ppc$lo[obs_names]
    sd_table$ppc_hi <- ppc$hi[obs_names]
    sd_table$ppc_in <- ppc$inside[obs_names]
  }

  sd_table <- sd_table[order(-sd_table$data_sd, na.last = TRUE), , drop = FALSE]
  rownames(sd_table) <- NULL

  n_fin <- sum(ok)
  pass <- if (self_compare || n_fin == 0L) NA else sum(in_band) / n_fin >= 0.7

  marginal_vars  <- obs_names[marginal]
  outside_vars   <- obs_names[outside]
  nonfinite_vars <- obs_names[!ok]
  outside_map  <- stats::setNames(outside,  obs_names)
  marginal_map <- stats::setNames(marginal, obs_names)

  ## ---- autocorrelations ---------------------------------------------------
  acf_table <- NULL
  if (!self_compare && !is.null(mm$acf_y) && !is.null(dm$acf_y) && n_obs > 0L) {
    K <- min(ncol(mm$acf_y), ncol(dm$acf_y))
    av <- intersect(obs_names, intersect(rownames(mm$acf_y), rownames(dm$acf_y)))
    if (K > 0L && length(av) > 0L) {
      mod <- as.numeric(mm$acf_y[av, seq_len(K), drop = FALSE])
      dat <- as.numeric(dm$acf_y[av, seq_len(K), drop = FALSE])
      acf_table <- data.frame(variable = rep(av, times = K),
                              lag = rep(seq_len(K), each = length(av)),
                              model = mod, data = dat, diff = mod - dat,
                              stringsAsFactors = FALSE)
    }
  }

  ## ---- contemporaneous cross-correlations ---------------------------------
  corr_table <- NULL
  if (!self_compare && !is.null(mm$correlation) && !is.null(dm$correlation)) {
    cv <- intersect(obs_names, intersect(rownames(mm$correlation),
                                         rownames(dm$correlation)))
    if (length(cv) >= 2L) {
      pr <- utils::combn(cv, 2L)
      a <- pr[1, ]; b <- pr[2, ]
      mod <- mm$correlation[cbind(a, b)]
      dat <- dm$correlation[cbind(a, b)]
      corr_table <- data.frame(pair = paste(a, b, sep = ":"), var1 = a, var2 = b,
                               model = as.numeric(mod), data = as.numeric(dat),
                               diff = as.numeric(mod - dat),
                               stringsAsFactors = FALSE)
      if (!is.null(key_pairs)) {
        rev_pairs <- paste(b, a, sep = ":")
        corr_table <- corr_table[corr_table$pair %in% key_pairs |
                                   rev_pairs %in% key_pairs, , drop = FALSE]
      }
      rownames(corr_table) <- NULL
      if (nrow(corr_table) == 0L) corr_table <- NULL
    }
  }

  ## ---- text ---------------------------------------------------------------
  note <- if (self_compare) " (self-comparison)" else ""
  badge <- if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL"
  summary_parts <- sprintf("D9 Moment matching%s: %d observables", note, n_obs)
  if (!self_compare && n_obs == 0L) {
    summary_parts <- paste0(summary_parts,
      "\n  No observable is common to the model and data moments -- nothing compared.")
  } else if (!self_compare) {
    rng <- if (n_fin > 0L) sprintf("[%.2f, %.2f]", min(ratio[ok]), max(ratio[ok])) else "[n/a]"
    summary_parts <- paste0(summary_parts, sprintf(
      "\n  SD ratio model/data %s -- %d/%d within [0.5, 2] (need >= 70%%)  %s",
      rng, sum(in_band), n_fin, badge))
    if (length(marginal_vars) > 0L)
      summary_parts <- paste0(summary_parts, sprintf("\n  marginal (<0.6 or >1.7): %s",
        paste(sprintf("%s=%.2f", marginal_vars, ratio[marginal]), collapse = ", ")))
    if (length(outside_vars) > 0L)
      summary_parts <- paste0(summary_parts, sprintf("\n  outside [0.5, 2]: %s",
        paste(sprintf("%s=%.2f", outside_vars, ratio[outside]), collapse = ", ")))
    if (length(nonfinite_vars) > 0L)
      summary_parts <- paste0(summary_parts, sprintf(
        "\n  non-finite ratio (excluded): %s", paste(nonfinite_vars, collapse = ", ")))
    if (any(me_v > 0))
      summary_parts <- paste0(summary_parts, sprintf(
        "\n  model SDs include measurement-error variance: %s",
        paste(sprintf("%s=%.4g", obs_names[me_v > 0], me_v[me_v > 0]),
              collapse = ", ")))
    if (!is.null(ppc))
      summary_parts <- paste0(summary_parts, sprintf(paste0(
        "\n  posterior-predictive band (INFO, %d draws, %g%%): %d/%d data SDs ",
        "inside%s"), ppc$n_used, 100 * ppc_level,
        sum(ppc$inside, na.rm = TRUE), sum(!is.na(ppc$inside)),
        if (any(!isTRUE(ppc$inside) & !is.na(ppc$inside)))
          paste0("; outside: ", paste(names(ppc$inside)[!ppc$inside &
                 !is.na(ppc$inside)], collapse = ", ")) else ""))
    if (!is.null(acf_table)) {
      a1 <- acf_table[acf_table$lag == 1L, , drop = FALSE]
      summary_parts <- paste0(summary_parts, sprintf("\n  lag-1 autocorr model/data: %s",
        paste(sprintf("%s %.2f/%.2f", a1$variable, a1$model, a1$data), collapse = ", ")))
    }
    if (!is.null(corr_table)) {
      j <- which.max(abs(corr_table$diff))
      if (length(j) == 1L)
        summary_parts <- paste0(summary_parts, sprintf(
          "\n  largest cross-corr gap: %s model %.2f vs data %.2f",
          corr_table$pair[j], corr_table$model[j], corr_table$data[j]))
    }
  }

  action <- if (self_compare) {
    "Self-comparison mode -- no data provided for benchmarking."
  } else if (is.na(pass)) {
    "No finite model/data SD ratio -- check observable names and model stationarity."
  } else if (isTRUE(pass)) {
    if (length(outside_vars) > 0L)
      sprintf("%d/%d model SDs within [0.5, 2] x data (>= 70%%); still review %s.",
              sum(in_band), n_fin, paste(outside_vars, collapse = ", "))
    else "All model SDs within [0.5, 2] x data SDs. Moment fit is reasonable."
  } else {
    over  <- obs_names[outside & ratio > hi]
    under <- obs_names[outside & ratio < lo]
    parts <- character(0)
    if (length(over) > 0L) parts <- c(parts, paste(paste(sprintf("%s=%.2f", over, ratio[match(over, obs_names)]), collapse = ", "), "over-predicted"))
    if (length(under) > 0L) parts <- c(parts, paste(paste(sprintf("%s=%.2f", under, ratio[match(under, obs_names)]), collapse = ", "), "under-predicted"))
    paste(paste(parts, collapse = "; "),
          "-- check shock variances or measurement equation.")
  }
  flag <- ifelse(status == "outside", " [OUTSIDE]",
          ifelse(status == "marginal", " [MARGINAL]",
          ifelse(status == "non-finite", " [NON-FINITE]", "")))
  llm_summary <- paste(c(
    sprintf("D9 | Moment Matching | %s", badge),
    sprintf("  observables=%d finite=%d within_0.5_2.0=%d/%d marginal=%d outside=%d",
            n_obs, n_fin, sum(in_band), n_fin, length(marginal_vars), length(outside_vars)),
    paste("  sd_ratios (model/data):",
          paste(sprintf("%s: ratio=%.2f%s", obs_names, ratio, flag), collapse = ", ")),
    if (any(me_v > 0))
      paste("  me_var (added to model SD):",
            paste(sprintf("%s=%.4g", obs_names, me_v), collapse = ", ")),
    if (!is.null(ppc))
      paste(sprintf("  ppc_band(%d draws, %g%%, INFO):", ppc$n_used,
                    100 * ppc_level),
            paste(sprintf("%s data_sd=%.3f in [%.3f, %.3f] %s", obs_names,
                          d_sd, ppc$lo[obs_names], ppc$hi[obs_names],
                          ifelse(ppc$inside[obs_names], "IN", "OUT")),
                  collapse = "; ")),
    sprintf("  action: %s", action)
  ), collapse = "\n")

  ## ---- plots --------------------------------------------------------------
  plots_d9 <- list()
  if (!self_compare && n_fin > 0L && requireNamespace("ggplot2", quietly = TRUE)) {
    plots_d9 <- .d9_plots(sd_table, acf_table, corr_table, lo, hi, m_lo, m_hi, meta)
  }

  structure(
    list(
      result  = list(sd_table = sd_table, acf_table = acf_table,
                     corr_table = corr_table, self_compare = self_compare,
                     marginal = marginal_vars, outside = outside_vars,
                     nonfinite = nonfinite_vars, me_var = me_v, ppc = ppc,
                     marginal_map = marginal_map, outside_map = outside_map),
      pass    = pass,
      plots   = plots_d9,
      summary = summary_parts,
      llm_summary = llm_summary
    ),
    class = "dynhr_diagnostic"
  )
}


## Align a measurement-error VARIANCE specification to `obs_names`.
##
## Accepts NULL (no ME), a scalar (the same variance on every observable, the
## `me_variance` scalar the estimation uses), a NAMED vector (matched by name;
## observables it does not name get 0), or an unnamed vector of length
## length(obs_names) (matched by position). A square matrix is read as a
## measurement-error covariance and its diagonal is used. Negative or
## non-finite entries are an error, not a silent zero.
.d9_align_me_var <- function(me_var, obs_names) {
  n <- length(obs_names)
  out <- stats::setNames(rep(0, n), obs_names)
  if (is.null(me_var) || n == 0L) return(out)
  if (is.matrix(me_var) && nrow(me_var) == ncol(me_var)) {
    d <- diag(me_var)
    if (is.null(names(d))) names(d) <- rownames(me_var) %||% colnames(me_var)
    me_var <- d
  }
  v <- as.numeric(me_var)
  names(v) <- names(me_var)
  if (any(!is.finite(v)) || any(v < 0))
    .dynhr_abort("d9_moment_matching(): `me_var` must be finite and non-negative.")
  if (length(v) == 1L && is.null(names(v))) {
    out[] <- v
  } else if (!is.null(names(v))) {
    hit <- intersect(obs_names, names(v))
    out[hit] <- v[hit]
  } else if (length(v) == n) {
    out[] <- v
  } else {
    .dynhr_abort(sprintf(paste0(
      "d9_moment_matching(): `me_var` has length %d with no names; expected 1 ",
      "or %d (one per observable), or a named vector."), length(v), n))
  }
  out
}


## Posterior-predictive band for the DATA standard deviations (Faust & Gupta
## 2012; Herbst & Schorfheide 2016 ch. 6). INFO only -- it never moves the
## badge.
##
## `ppc_sd_fn(theta)` must return a named vector of sample SDs computed from a
## T-length path SIMULATED at theta (T = the estimation sample length), so the
## spread of the returned SDs carries BOTH parameter uncertainty (across
## draws) and finite-sample sampling variation (within a draw). It follows the
## package callback convention: a non-finite / short return for an infeasible
## draw is dropped, it must not stop().
##
## Deliberately NOT an analytic Bartlett/delta-method band: macro series are
## short and highly autocorrelated, exactly where asymptotic SEs for SDs and
## ACFs are least reliable, and no standard DSGE package ships one.
.d9_ppc_band <- function(ppc_sd_fn, ppc_draws, obs_names, d_sd, ppc_n,
                         ppc_level) {
  if (!is.function(ppc_sd_fn) || is.null(ppc_draws)) return(NULL)
  D <- as.matrix(ppc_draws)
  if (nrow(D) < 10L || length(obs_names) == 0L) return(NULL)
  if (!is.finite(ppc_level) || ppc_level <= 0 || ppc_level >= 1)
    .dynhr_abort("d9_moment_matching(): `ppc_level` must be in (0, 1).")
  ppc_n <- as.integer(ppc_n)
  idx <- if (nrow(D) <= ppc_n) seq_len(nrow(D)) else
    round(seq(1, nrow(D), length.out = ppc_n))
  sims <- matrix(NA_real_, length(idx), length(obs_names),
                 dimnames = list(NULL, obs_names))
  for (k in seq_along(idx)) {
    sd_k <- ppc_sd_fn(D[idx[k], ])
    if (is.null(sd_k) || is.null(names(sd_k))) next
    hit <- intersect(obs_names, names(sd_k))
    if (length(hit) == 0L) next
    sims[k, hit] <- as.numeric(sd_k[hit])
  }
  n_used <- sum(apply(sims, 1L, function(r) any(is.finite(r))))
  if (n_used < 10L) return(NULL)
  a <- (1 - ppc_level) / 2
  q <- apply(sims, 2L, function(cl) {
    cl <- cl[is.finite(cl)]
    if (length(cl) < 10L) return(c(NA_real_, NA_real_))
    unname(stats::quantile(cl, c(a, 1 - a)))
  })
  lo <- stats::setNames(q[1, ], obs_names)
  hi <- stats::setNames(q[2, ], obs_names)
  inside <- stats::setNames(is.finite(lo) & is.finite(hi) & is.finite(d_sd) &
                              d_sd >= lo & d_sd <= hi, obs_names)
  inside[!is.finite(lo) | !is.finite(hi)] <- NA
  list(lo = lo, hi = hi, inside = inside, sims = sims, n_used = n_used,
       level = ppc_level)
}


## D9 plots: SD-ratio dot plot (badge geometry), ACF model-vs-data, and
## cross-correlation model-vs-data dumbbells.
.d9_plots <- function(sd_table, acf_table, corr_table, lo, hi, m_lo, m_hi, meta) {
  out <- list()
  src_cols <- c(Model = dynhr_palette_vibrant[1], Data = dynhr_palette_vibrant[4])
  st_cols  <- c(ok = dynhr_palette_vibrant[2], marginal = dynhr_palette_vibrant[4],
                outside = dynhr_palette_vibrant[6])

  st <- sd_table[is.finite(sd_table$ratio), , drop = FALSE]
  st$ratio_plot <- pmin(pmax(st$ratio, 1e-3), 1e3)
  st$variable <- factor(st$variable, levels = rev(st$variable[order(st$ratio_plot)]))
  st$status <- factor(st$status, levels = names(st_cols))
  st$label <- sprintf("%.2f", st$ratio)
  lim <- range(c(st$ratio_plot, lo / 1.5, hi * 1.5))
  p_sd <- ggplot2::ggplot(st, ggplot2::aes(x = .data$ratio_plot, y = .data$variable)) +
    ggplot2::annotate("rect", xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf,
                      fill = dynhr_palette_light[2], alpha = 0.45) +
    ggplot2::geom_vline(xintercept = c(lo, hi), linetype = "dashed",
                        colour = dynhr_na_colour) +
    ggplot2::geom_vline(xintercept = c(m_lo, m_hi), linetype = "dotted",
                        colour = dynhr_na_colour) +
    ggplot2::geom_vline(xintercept = 1, colour = dynhr_na_colour) +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = .data$ratio_plot,
                                       yend = .data$variable, colour = status),
                          linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(colour = status), size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.9, size = 3.2) +
    ggplot2::scale_x_continuous(trans = "log2", limits = lim,
                                breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8),
                                labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8")) +
    ggplot2::scale_colour_manual(values = st_cols, drop = FALSE, name = "Status") +
    theme_dynhr() +
    ggplot2::labs(
      title = "D9: Model / data standard-deviation ratio",
      subtitle = "Shaded: pass band [1/2, 2]; dotted: marginal limits 0.6, 1.7; PASS if >= 70% inside",
      x = "Model SD / data SD (log2 scale)", y = NULL)
  out$moments_sd <- .apply_meta(p_sd, meta)

  ## Posterior-predictive band for the DATA SD (INFO; never gates).
  if (all(c("ppc_lo", "ppc_hi") %in% names(sd_table))) {
    pb <- sd_table[is.finite(sd_table$ppc_lo) & is.finite(sd_table$ppc_hi), ,
                   drop = FALSE]
    if (nrow(pb) > 0L) {
      pb$variable <- factor(pb$variable, levels = rev(pb$variable))
      pb$Verdict <- factor(ifelse(pb$ppc_in %in% TRUE,
                                  "inside", "outside"),
                           levels = c("inside", "outside"))
      p_ppc <- ggplot2::ggplot(pb, ggplot2::aes(y = .data$variable)) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = .data$ppc_lo, xmax = .data$ppc_hi),
          height = 0.25, colour = dynhr_palette_vibrant[1], linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(x = .data$data_sd, colour = .data$Verdict),
                            size = 3) +
        ggplot2::scale_colour_manual(
          values = c(inside = dynhr_palette_vibrant[2],
                     outside = dynhr_palette_vibrant[6]),
          drop = FALSE, name = "Data SD") +
        theme_dynhr() +
        ggplot2::labs(
          title = "D9: Posterior-predictive band for the data SD",
          subtitle = paste0("Bar: SD of a T-length path simulated at posterior ",
                            "draws; dot: observed SD. Informational only."),
          x = "Standard deviation", y = NULL)
      out$moments_ppc <- .apply_meta(p_ppc, meta)
    }
  }

  if (!is.null(acf_table) && any(is.finite(acf_table$model) | is.finite(acf_table$data))) {
    al <- rbind(
      data.frame(variable = acf_table$variable, lag = acf_table$lag,
                 value = acf_table$model, Source = "Model", stringsAsFactors = FALSE),
      data.frame(variable = acf_table$variable, lag = acf_table$lag,
                 value = acf_table$data, Source = "Data", stringsAsFactors = FALSE))
    al$Source <- factor(al$Source, levels = c("Model", "Data"))
    al$variable <- factor(al$variable, levels = unique(acf_table$variable))
    p_acf <- ggplot2::ggplot(al, ggplot2::aes(x = lag, y = value, colour = Source)) +
      ggplot2::geom_hline(yintercept = 0, colour = dynhr_na_colour) +
      ggplot2::geom_line(ggplot2::aes(linetype = Source), linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(shape = Source), size = 2.2, na.rm = TRUE) +
      ggplot2::facet_wrap(~ variable) +
      ggplot2::scale_colour_manual(values = src_cols) +
      ggplot2::scale_x_continuous(breaks = sort(unique(al$lag))) +
      ggplot2::coord_cartesian(ylim = c(-1, 1)) +
      theme_dynhr() +
      ggplot2::labs(title = "D9: Autocorrelation, model vs data",
                    subtitle = "Informational (not part of the badge); no sampling bands",
                    x = "Lag (periods)", y = "Autocorrelation",
                    colour = NULL, linetype = NULL, shape = NULL)
    out$moments_acf <- .apply_meta(p_acf, meta)
  }

  if (!is.null(corr_table)) {
    ct <- corr_table
    ct$pair <- factor(ct$pair, levels = ct$pair[order(abs(ct$diff), na.last = FALSE)])
    cl <- rbind(
      data.frame(pair = ct$pair, value = ct$model, Source = "Model", stringsAsFactors = FALSE),
      data.frame(pair = ct$pair, value = ct$data, Source = "Data", stringsAsFactors = FALSE))
    cl$Source <- factor(cl$Source, levels = c("Model", "Data"))
    p_corr <- ggplot2::ggplot(ct, ggplot2::aes(y = .data$pair)) +
      ggplot2::geom_vline(xintercept = 0, colour = dynhr_na_colour) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$model, xend = .data$data, yend = .data$pair),
                            colour = dynhr_na_colour, linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = cl, ggplot2::aes(x = value, colour = Source, shape = Source),
                          size = 3, na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = src_cols) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      theme_dynhr() +
      ggplot2::labs(title = "D9: Contemporaneous cross-correlation, model vs data",
                    subtitle = "Informational (not part of the badge); largest gaps at top",
                    x = "Correlation", y = NULL, colour = NULL, shape = NULL)
    out$moments_corr <- .apply_meta(p_corr, meta)
  }
  out
}
