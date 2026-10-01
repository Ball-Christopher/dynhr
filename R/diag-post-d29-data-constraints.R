## R/diag-post-d29-data-constraints.R
## --------------------------------------------------------------------------
## Phase H: D29 — Data constraints (Stock-Wright S test at theta)
##
## What D29 measures, for the var / autocovariance moments m of the observed
## series (m-hat from the data, m(theta) from the model): the Stock-Wright
## (2000) S statistic
##      S(theta) = (m-hat - m(theta))' V^{-1} (m-hat - m(theta)),
## asymptotically chi-square(n_moments) at the true theta WHATEVER the
## identification strength (it involves no estimate of theta), i.e. a
## weak-identification-robust test that the data are consistent with the
## model at theta. Reported, never gated: the diagnostic is INFO.
##
## V = Var(m-hat) is by preference MODEL-IMPLIED at theta: the Gaussian
## (Bartlett/Hannan) long-run covariance of the sample autocovariances built
## from the model's own autocovariance function, divided by T -- the
## Christiano-Eichenbaum-Trabandt convention of taking the moment estimator's
## covariance from the model at theta instead of from the sample. A sample
## HAC V is used only when no model autocovariance function is supplied, and
## the result says which was used.
##
## 0.9.4 changes: V switched from sample HAC to model-implied (the HAC S test
## over-rejected badly -- ~13-37% at a nominal 5% for T ~ 200-400 in the
## our simulations, vs ~7% with the model-implied V); D29's own
## identification-strength table was DROPPED in favour of D20 with
## weighting = "sampling" (two divergent strength computations with
## inconsistent thresholds, 2 here vs 1 there); the badge became INFO.
##
## Until 0.9.4: V was the long-run covariance NOT divided by T (every
## strength was sqrt(T) too small); the orchestrator's moment function reused
## the fixed decision rule AND dr$Sigma_e, so the whole Jacobian was zero and
## D29 flagged every parameter; sd/acf1 moments were weighted by the
## covariance of var/acv moments (no delta method); equal-length named
## moment vectors were never aligned by name; and a near-singular V was
## silently replaced by its diagonal.
##
## References:
##   Stock, J. H., & Wright, J. H. (2000). GMM with weak identification.
##     Econometrica, 68(5), 1055-1096.
##   Christiano, L. J., Eichenbaum, M. S., & Trabandt, M. (2016).
##     Unemployment and business cycles. Econometrica, 84(4), 1523-1569.
##     (Moment-estimator covariance taken from the model at theta.)
##   Ruge-Murcia, F. J. (2007). Methods to estimate DSGE models. JEDC, 31(8),
##     2599-2636. (The sample Newey-West/Bartlett V this replaces.)
##   Hannan, E. J. (1970). Multiple Time Series. Wiley, ch. IV (the Gaussian
##     long-run covariance of sample autocovariances, "Bartlett's formula").
##   Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation
##     consistent covariance matrix estimation. Econometrica, 59(3), 817-858.
## --------------------------------------------------------------------------

#' D29. Data Constraints Diagnostic
#'
#' Tests whether the data are consistent with the model at \code{theta} with
#' the identification-robust Stock-Wright S statistic, and reports the
#' standardised gap between each data moment and its model counterpart.
#'
#' \strong{The badge is always INFO.} Stock & Wright (2000) derive the
#' asymptotic \eqn{\chi^2} null of S under a \emph{known or consistently
#' estimated} weight and make no finite-sample size claim; the statistic is
#' well documented to over-reject in short, persistent samples -- in this
#' package's own simulations the sample-HAC version rejected
#' 13-37\% of the time at a nominal 5\% for \eqn{T \approx 200}-400, and the
#' model-implied weight below brings that to roughly 7\%. Even at 7\% a
#' rejection is not strong enough evidence to gate a badge on, so S is
#' reported and never gated.
#'
#' \strong{Weight matrix.} \eqn{V = Var(\hat m)} is model-implied whenever
#' \code{model_acov_fn} is supplied: the Gaussian (Bartlett/Hannan) long-run
#' covariance of the sample autocovariances evaluated at the \emph{model's}
#' autocovariance function at \code{theta}, divided by \eqn{T}. This is the
#' Christiano-Eichenbaum-Trabandt (2016) SMM convention -- take the moment
#' estimator's covariance from the model at \code{theta}, not from the sample
#' -- and it removes the estimation noise in \eqn{\hat\Omega} that drives the
#' over-rejection. Without \code{model_acov_fn} the sample Bartlett HAC
#' (Andrews 1991 plug-in bandwidth) is used instead and
#' \code{result$moment_cov_source} records which.
#'
#' \strong{Identification strength moved to D20.} Until 0.9.4 D29 also
#' reported its own sampling-weighted \eqn{|\theta_i|/SE_i} table with a
#' \code{strength_threshold = 2} gate, duplicating (with an inconsistent
#' threshold) D20's \code{weighting = "sampling"} path. That table is gone:
#' call \code{d20_fisher_identification_strength(..., weighting = "sampling",
#' moment_cov = <this result's moment_cov>)} for data-weighted strengths.
#'
#' @param data            T x n_obs data matrix (observable variables, named
#'   columns).
#' @param model_solve_fn  Function: theta -> named numeric vector of
#'   model-implied POPULATION moments. Names must follow
#'   \code{.compute_data_moments} (\code{var_<obs>}, \code{acvK_<obs>}) or
#'   the \code{sd.<obs>}/\code{sd_<obs>} + \code{acf1.<obs>}/\code{acf1_<obs>}
#'   convention. Moments are aligned by name; an unnamed vector is taken in
#'   \code{.compute_data_moments} order. It must re-solve the model at theta
#'   (a fixed decision rule makes structural parameters look unidentified).
#' @param theta           Named numeric vector of parameter values.
#' @param param_names     Character vector of parameter names (default
#'   \code{names(theta)}).
#' @param moment_names    Ignored when moments are named; otherwise labels.
#' @param max_lag         Maximum autocovariance lag (default 4).
#' @param model_acov_fn   Optional function \code{function(K)} returning the
#'   model-implied autocovariances of the data columns at \code{theta} as an
#'   \code{n_obs x n_obs x (K+1)} array with
#'   \code{g[a, b, k + 1] = Cov(y[a, t + k], y[b, t])} and dimnames matching
#'   \code{colnames(data)}. When supplied, \eqn{Var(\hat m)} is built from it
#'   (preferred; see Details). \code{NULL} falls back to the sample estimator.
#' @param acov_trunc      Truncation \eqn{K} of the Bartlett sum over the
#'   model autocovariance function (default 200). The autocovariances decay
#'   geometrically, so the truncation error is negligible unless the model is
#'   very close to a unit root; \code{result$acov_tail} reports the size of
#'   the retained tail relative to the variance as a check. 200 is a package
#'   choice, not a literature-specified number.
#' @param use_hac         Logical: when \code{model_acov_fn} is NULL, use the
#'   sample Bartlett HAC with the Andrews (1991) plug-in bandwidth (default
#'   TRUE); FALSE uses the i.i.d. estimator. Ignored when
#'   \code{model_acov_fn} is supplied.
#' @param s_level         Size of the reported S test (default 0.05).
#' @param verbose         Print progress messages.
#' @param meta            Plot caption metadata.
#'
#' @return A \code{dynhr_diagnostic} with \code{pass = NA} (always INFO);
#'   \code{result} holds \code{data_moments}, \code{model_moments},
#'   \code{moment_cov} (= Var(m-hat)), \code{moment_cov_source}
#'   (\code{"model-implied"} or \code{"sample-HAC"}/\code{"sample-iid"}),
#'   \code{acov_tail}, \code{S_stat}, \code{S_df}, \code{S_pvalue},
#'   \code{S_reject}, \code{moment_z} (per-moment (m-hat - m)/sd),
#'   \code{bandwidth}, \code{T_obs} and \code{n_moments}.
#'
#' @references
#'   Stock, J. H., & Wright, J. H. (2000). GMM with weak identification.
#'     \emph{Econometrica}, 68(5), 1055-1096.
#'   Christiano, L. J., Eichenbaum, M. S., & Trabandt, M. (2016).
#'     Unemployment and business cycles. \emph{Econometrica}, 84(4),
#'     1523-1569.
#'   Ruge-Murcia, F. J. (2007). Methods to estimate dynamic stochastic general
#'     equilibrium models. \emph{Journal of Economic Dynamics and Control},
#'     31(8), 2599-2636.
#'   Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation
#'     consistent covariance matrix estimation. \emph{Econometrica}, 59(3),
#'     817-858.
#'
#' @noRd
d29_data_driven_constraints <- function(data,
                                         model_solve_fn,
                                         theta,
                                         param_names = NULL,
                                         moment_names = NULL,
                                         max_lag = 4L,
                                         model_acov_fn = NULL,
                                         acov_trunc = 200L,
                                         use_hac = TRUE,
                                         s_level = 0.05,
                                         verbose = FALSE,
                                         meta = NULL) {
  na_result <- function(msg, result = NULL, errored = FALSE)
    .make_result(result = result, pass = NA,
                 summary = paste0("D29 Data constraints: ", msg), errored = errored)

  # ---- 1. Validate ----
  if (is.null(data) || is.null(model_solve_fn) || is.null(theta))
    return(na_result("data, model_solve_fn, and theta are required."))
  if (!is.function(model_solve_fn))
    return(na_result("model_solve_fn must be a function.", errored = TRUE))
  data <- as.matrix(data)
  if (!is.numeric(data) || anyNA(data))
    return(na_result("data must be a numeric matrix without missing values.", errored = TRUE))
  if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_along(theta))
  if (length(param_names) != length(theta))
    return(na_result(sprintf("param_names has %d entries but theta has %d.",
                             length(param_names), length(theta)), errored = TRUE))
  names(theta) <- param_names
  n_par <- length(theta)
  T_obs <- nrow(data)
  n_obs <- ncol(data)
  if (is.null(colnames(data))) colnames(data) <- as.character(seq_len(n_obs))
  max_lag <- as.integer(max_lag)
  if (T_obs <= max_lag + 1L)
    return(na_result(sprintf("T=%d is too short for max_lag=%d.", T_obs, max_lag)))

  # ---- 2. Data moments m-hat and their sampling covariance V = Var(m-hat) ----
  # Preferred: MODEL-IMPLIED V at theta (Bartlett/Hannan long-run covariance
  # of the sample autocovariances, evaluated at the model's own
  # autocovariance function). Falls back to the sample estimator only when no
  # model autocovariance function is available; the result records which.
  m_raw <- .compute_data_moments(data, max_lag = max_lag)
  acov_tail <- NA_real_
  V_raw <- NULL
  cov_source <- NULL
  cov_note <- NULL
  bandwidth <- NA_integer_
  if (is.function(model_acov_fn)) {
    gam <- model_acov_fn(as.integer(acov_trunc))
    if (is.null(gam) || !all(is.finite(gam))) {
      cov_note <- paste0(
        "the model autocovariance function is NULL or non-finite at theta ",
        "(the model does not solve / is not stationary); fell back to the ",
        "sample moment covariance")
    } else {
      mv <- .d29_model_moment_cov(gam, T_obs = T_obs, max_lag = max_lag,
                                  obs_names = colnames(data))
      V_raw <- mv$V
      acov_tail <- mv$tail
      cov_source <- "model-implied"
    }
  }
  if (is.null(V_raw)) {
    if (isTRUE(use_hac)) {
      V_raw <- .d20_moment_sampling_cov(data, max_lag = max_lag)
      bandwidth <- attr(V_raw, "bandwidth")
      attr(V_raw, "bandwidth") <- NULL
      cov_source <- "sample-HAC"
    } else {
      V_raw <- .compute_moment_covariance(data, max_lag = max_lag, use_hac = FALSE) / T_obs
      dimnames(V_raw) <- list(names(m_raw), names(m_raw))
      cov_source <- "sample-iid"
    }
  }

  # ---- 3. Align model moments with data moments (by name) ----
  f0 <- model_solve_fn(theta)
  if (is.null(f0) || length(f0) == 0L)
    return(na_result("model_solve_fn returned nothing at theta.", errored = TRUE))
  al <- .d29_align_moments(f0, m_raw, V_raw, colnames(data), moment_names)
  if (!is.null(al$error))
    return(na_result(al$error, errored = TRUE))
  if (!is.null(al$warning)) .dynhr_warn(al$warning)
  data_moments <- al$m
  moment_cov   <- al$V
  moment_names <- names(data_moments)
  n_mom <- length(data_moments)
  mfun <- function(th) {
    out <- model_solve_fn(th)
    if (length(out) < max(al$idx)) return(rep(NA_real_, n_mom))
    stats::setNames(as.numeric(out[al$idx]), moment_names)
  }
  model_moments <- mfun(theta)
  if (verbose) .dynhr_cat(sprintf("[d29] %d moments (%s) from %d obs x %d variables.\n",
                                  n_mom, al$format, T_obs, n_obs))

  # ---- 4. Model moments must be finite ----
  base_result <- list(data_moments = data_moments, model_moments = model_moments,
                      moment_cov = moment_cov, moment_cov_source = cov_source,
                      acov_tail = acov_tail,
                      bandwidth = bandwidth, T_obs = T_obs, n_moments = n_mom)
  if (!all(is.finite(model_moments)))
    return(na_result(paste0(
      "the model moments are non-finite at theta (the model does not solve / ",
      "is not stationary). The S statistic is undetermined -- a numerical ",
      "failure, not a finding."), result = base_result))

  # ---- 5. V must be a usable (positive definite) weight: no silent fallback ----
  chk <- .d29_omega_check(moment_cov)
  if (!chk$ok)
    return(na_result(sprintf(paste0(
      "the moment covariance Var(m-hat) is not positive definite (%s; ",
      "%d moments, T=%d). Data strength and the S test are undefined -- use ",
      "fewer moments (smaller max_lag / fewer observables) or a longer sample."),
      chk$reason, n_mom, T_obs), result = base_result))
  W <- .robust_Omega_inv_sqrt(moment_cov)
  dimnames(W) <- list(moment_names, moment_names)

  # ---- 6. Stock-Wright S statistic at theta ----
  resid <- data_moments - model_moments
  S_stat <- sum((W %*% resid)^2)
  S_df <- n_mom
  S_pvalue <- stats::pchisq(S_stat, df = S_df, lower.tail = FALSE)
  moment_z <- resid / sqrt(diag(moment_cov))

  # ---- 7. Plots ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$moment_fit <- .apply_meta(.d29_plot_moment_fit(
      moment_z, S_stat, S_df, S_pvalue), meta)
  }

  # ---- 8. Result (always INFO: the S test never gates -- see Details) ----
  result <- c(base_result, list(
    S_stat = S_stat, S_df = S_df, S_pvalue = S_pvalue, S_level = s_level,
    S_reject = S_pvalue < s_level, moment_z = moment_z,
    moment_format = al$format))

  cov_str <- if (identical(cov_source, "model-implied"))
    sprintf("Var(m-hat) is model-implied at theta (Bartlett/Hannan sum truncated at K=%d; retained tail %.1e of the variance).",
            as.integer(acov_trunc), acov_tail)
  else sprintf("Var(m-hat) is the %s SAMPLE estimator%s -- the S test then over-rejects in short or persistent samples; supply model_acov_fn for the model-implied weight.",
               cov_source,
               if (is.na(bandwidth)) "" else sprintf(" (bandwidth %d)", bandwidth))
  if (!is.null(cov_note)) cov_str <- paste0(cov_str, " (", cov_note, ")")
  worst <- names(moment_z)[order(abs(moment_z), decreasing = TRUE)][seq_len(min(3L, n_mom))]
  s_str <- sprintf(
    "S test (Stock-Wright, weak-identification robust; asymptotic) at theta: S=%.2f, df=%d, p=%.3g -- %s at %g%%. INFO ONLY: this statistic over-rejects in finite samples and never gates the badge.",
    S_stat, S_df, S_pvalue,
    if (S_pvalue < s_level) "data REJECT the model moments" else "not rejected",
    100 * s_level)

  .make_result(
    result  = result,
    pass    = NA,
    plots   = plots,
    summary = sprintf("D29 Data constraints: %d params, %d moments, T=%d. %s %s Worst-fitting moments: %s. For parameter identification strength use D20 with weighting=\"sampling\" and this result's moment_cov.",
                      n_par, n_mom, T_obs, s_str, cov_str,
                      paste(sprintf("%s (z=%.2f)", worst, moment_z[worst]), collapse = ", ")),
    llm_summary = paste(c(
      "D29 | Data Constraints (Stock-Wright S) | INFO",
      sprintf("  n_params=%d n_moments=%d T=%d moment_cov=%s S=%.2f df=%d p=%.3g reject_at_%g=%s",
              n_par, n_mom, T_obs, cov_source, S_stat, S_df, S_pvalue,
              s_level, S_pvalue < s_level),
      sprintf("  worst_moments: %s",
              paste(sprintf("%s=%.2f", worst, moment_z[worst]), collapse = ", ")),
      sprintf("  action: %s",
              if (S_pvalue < s_level)
                sprintf(paste0("The data moments %s sit far from the model's at theta. Treat as a ",
                               "HINT, not a verdict: the S test over-rejects in short/persistent ",
                               "samples (info only, never gates)."),
                        paste(head(worst, 2), collapse = " and "))
              else "Data and model moments are consistent at theta (info only; the S test never gates)."),
      sprintf("  note: identification strength is no longer reported here -- use D20 with weighting=\"sampling\" and this result's moment_cov.")
    ), collapse = "\n")
  )
}


# ==========================================================================
# Internal helpers for D29
# ==========================================================================

#' Align model moments with the data moments (D29)
#'
#' Returns \code{idx} (positions in the model output), the matching data
#' moments \code{m} and their covariance \code{V} (delta method when the model
#' reports sd/acf1), the \code{format}, and \code{error}/\code{warning}.
#' @noRd
.d29_align_moments <- function(f0, m_raw, V_raw, obs, moment_names = NULL) {
  fn <- names(f0)
  out <- list(idx = NULL, m = NULL, V = NULL, format = "var/acv",
              error = NULL, warning = NULL)
  if (is.null(fn)) {
    n <- min(length(f0), length(m_raw))
    if (length(f0) != length(m_raw))
      out$warning <- sprintf(
        "d29: model moments are unnamed and have length %d (data: %d); using the first %d in data-moment order (%s, ...).",
        length(f0), length(m_raw), n, paste(head(names(m_raw), 3), collapse = ", "))
    out$idx <- seq_len(n)
    out$m <- m_raw[seq_len(n)]
    out$V <- V_raw[seq_len(n), seq_len(n), drop = FALSE]
    return(out)
  }
  common <- intersect(fn, names(m_raw))
  if (length(common) > 0L) {
    out$idx <- match(common, fn)
    out$m <- m_raw[common]
    out$V <- V_raw[common, common, drop = FALSE]
    return(out)
  }
  # sd/acf1 format: g(m) = (sqrt(var), acv1/var), Var(g) = G V G'.
  keys <- character(0); rows <- list(); gm <- numeric(0)
  for (o in obs) {
    v <- paste0("var_", o); a <- paste0("acv1_", o)
    if (!all(c(v, a) %in% names(m_raw))) next
    sd_nm <- intersect(paste0(c("sd.", "sd_"), o), fn)
    ac_nm <- intersect(paste0(c("acf1.", "acf1_"), o), fn)
    if (length(sd_nm)) {
      g <- stats::setNames(numeric(length(m_raw)), names(m_raw))
      g[v] <- 1 / (2 * sqrt(m_raw[[v]]))
      keys <- c(keys, sd_nm[1]); rows[[length(rows) + 1L]] <- g
      gm <- c(gm, sqrt(m_raw[[v]]))
    }
    if (length(ac_nm)) {
      g <- stats::setNames(numeric(length(m_raw)), names(m_raw))
      g[v] <- -m_raw[[a]] / m_raw[[v]]^2
      g[a] <- 1 / m_raw[[v]]
      keys <- c(keys, ac_nm[1]); rows[[length(rows) + 1L]] <- g
      gm <- c(gm, m_raw[[a]] / m_raw[[v]])
    }
  }
  if (!length(keys)) {
    out$error <- sprintf(paste0(
      "no model moment matches the data moments. Model names: %s. Data names: %s. ",
      "Return var_<obs>/acvK_<obs> (or sd/acf1) moments for the data columns."),
      paste(head(fn, 4), collapse = ", "), paste(head(names(m_raw), 4), collapse = ", "))
    return(out)
  }
  G <- do.call(rbind, rows)
  names(gm) <- keys
  out$idx <- match(keys, fn)
  out$m <- gm
  out$V <- G %*% V_raw[colnames(G), colnames(G)] %*% t(G)
  dimnames(out$V) <- list(keys, keys)
  out$format <- "sd/acf1 (delta method)"
  out
}

#' Model-implied sampling covariance of the var/autocovariance moments (D29)
#'
#' Bartlett's (Hannan 1970, ch. IV) Gaussian long-run covariance of sample
#' autocovariances, evaluated at the MODEL's autocovariance function rather
#' than at sample autocovariances. With
#' \eqn{G_{ab}(h) = Cov(y_{a,t+h}, y_{b,t})} and moments
#' \eqn{\hat m_{a,i} = \hat\gamma_{aa}(i)},
#' \deqn{T\,Cov(\hat m_{a,i}, \hat m_{b,j}) = \sum_d
#'   G_{ab}(d+i-j) G_{ab}(d) + G_{ab}(d+i) G_{ab}(d-j),}
#' truncated at \eqn{|d| \le K}. This is the Christiano-Eichenbaum-Trabandt
#' convention (moment-estimator covariance from the model at theta) and
#' carries no \eqn{\hat\Omega} estimation noise.
#'
#' @param gam Array \code{n_obs x n_obs x (K+1)} with
#'   \code{gam[a, b, k + 1] = G_ab(k)}; dimnames 1-2 are the observable names.
#' @param T_obs Sample size (the long-run covariance is divided by it).
#' @param max_lag Highest autocovariance lag among the moments.
#' @param obs_names Observable names, in data-column order.
#' @return list(V = named n_mom x n_mom covariance in
#'   \code{.compute_data_moments} order, tail = largest retained
#'   \eqn{|G_{ab}(K)| / \sqrt{G_{aa}(0) G_{bb}(0)}}, a truncation check).
#' @noRd
.d29_model_moment_cov <- function(gam, T_obs, max_lag, obs_names) {
  K <- dim(gam)[3L] - 1L
  n_obs <- length(obs_names)
  stopifnot(dim(gam)[1L] == n_obs, dim(gam)[2L] == n_obs, K >= max_lag)

  # Two-sided index: gfull[a, b, h + K + 1] = G_ab(h), G_ab(-h) = G_ba(h).
  gfull <- array(0, c(n_obs, n_obs, 2L * K + 1L))
  for (k in 0:K) {
    gfull[, , K + 1L + k] <- gam[, , k + 1L]
    if (k > 0L) gfull[, , K + 1L - k] <- t(gam[, , k + 1L])
  }

  sd0 <- sqrt(pmax(diag(matrix(gam[, , 1L], n_obs, n_obs)), 0))
  denom <- outer(sd0, sd0)
  denom[denom <= 0] <- Inf
  tail_rel <- max(abs(gam[, , K + 1L]) / denom)

  lags <- 0:max_lag
  nm <- unlist(lapply(lags, function(l)
    paste0(if (l == 0L) "var" else paste0("acv", l), "_", obs_names)))
  n_mom <- length(nm)
  idx_obs <- rep(seq_len(n_obs), times = length(lags))
  idx_lag <- rep(lags, each = n_obs)

  # Sum over |d| <= K - max_lag so every shifted index stays inside gfull.
  D <- K - max_lag
  d <- seq.int(-D, D)
  V <- matrix(0, n_mom, n_mom, dimnames = list(nm, nm))
  for (p in seq_len(n_mom)) {
    a <- idx_obs[p]; i <- idx_lag[p]
    for (q in p:n_mom) {
      b <- idx_obs[q]; j <- idx_lag[q]
      g <- gfull[a, b, ]
      val <- sum(g[d + i - j + K + 1L] * g[d + K + 1L]) +
             sum(g[d + i + K + 1L] * g[d - j + K + 1L])
      V[p, q] <- V[q, p] <- val / T_obs
    }
  }
  list(V = V, tail = tail_rel)
}


#' Model-implied observable autocovariance function at theta (D29)
#'
#' Re-solves the model at \code{theta} and returns
#' \code{g[a, b, k + 1] = Cov(y[a, t + k], y[b, t])} for the observables,
#' \code{k = 0..K}; NULL when the model does not solve. Used by the
#' orchestrator to give D29 a model-implied \eqn{Var(\hat m)}.
#' @noRd
.d29_model_acov_fn <- function(model, compiled, params, obs_names) {
  sys_cache <- cache_system_structure(compiled)
  state <- new.env(parent = emptyenv())
  function(theta, K) {
    pp  <- .apply_theta_to_params(model, theta, params)
    sol <- .solve_dr_for_theta(model, compiled, sys_cache, pp, state)
    if (is.null(sol)) return(NULL)
    K <- as.integer(K)
    mm <- compute_moments(sol$dr, model, n_ar = K, params = sol$params)
    V0 <- mm$var_cov[obs_names, obs_names, drop = FALSE]
    sdo <- outer(sqrt(pmax(diag(V0), 0)), sqrt(pmax(diag(V0), 0)))
    g <- array(NA_real_, c(length(obs_names), length(obs_names), K + 1L),
               dimnames = list(obs_names, obs_names, NULL))
    g[, , 1L] <- V0
    n_o <- length(obs_names)
    for (k in seq_len(K))
      g[, , k + 1L] <- matrix(mm$autocorr[obs_names, obs_names, k], n_o, n_o) * sdo
    g
  }
}


#' Is a moment covariance usable as a GMM weight? (D29)
#'
#' Symmetric, finite, and positive definite with reciprocal condition number
#' (eigenvalue ratio) above \code{n * eps}. A sample covariance of more
#' moments than (effective) observations fails this.
#' @return list(ok, rcond, reason)
#' @noRd
.d29_omega_check <- function(Omega) {
  if (!is.matrix(Omega) || nrow(Omega) != ncol(Omega) || nrow(Omega) == 0L)
    return(list(ok = FALSE, rcond = NA_real_, reason = "not a square matrix"))
  if (!all(is.finite(Omega)))
    return(list(ok = FALSE, rcond = NA_real_, reason = "non-finite entries"))
  if (max(abs(Omega - t(Omega))) > 1e-10 * max(abs(Omega), 1e-300))
    return(list(ok = FALSE, rcond = NA_real_, reason = "not symmetric"))
  ev <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  rc <- if (max(ev) > 0) min(ev) / max(ev) else -Inf
  if (!(rc > nrow(Omega) * .Machine$double.eps))
    return(list(ok = FALSE, rcond = rc,
                reason = sprintf("eigenvalue ratio %.2e", rc)))
  list(ok = TRUE, rcond = rc, reason = "")
}

#' Omega^{-1/2} (inverse transposed Cholesky factor) of a PD covariance
#'
#' \code{W = t(R^{-1})} with \code{Omega = R'R}, so \code{crossprod(W) =
#' Omega^{-1}}. Refuses (error) a matrix that fails
#' \code{.d29_omega_check}: a diagonal or ridge substitute would silently
#' change the weighting (it used to fall back to \code{diag(Omega)}).
#' @noRd
.robust_Omega_inv_sqrt <- function(Omega) {
  chk <- .d29_omega_check(Omega)
  if (!chk$ok)
    .dynhr_abort("Moment covariance is not positive definite (", chk$reason,
                 "); refusing to weight by it (no diagonal fallback).",
                 class = "dynhr_error_singular_moment_cov")
  R <- chol(Omega)
  t(backsolve(R, diag(nrow(R))))
}

#' D29 moment-fit plot: standardised moment gaps and the S statistic
#' @noRd
.d29_plot_moment_fit <- function(z, S, df, p) {
  # group by observable, lag order within (var, acv1, ...)
  obs <- sub("^(var|acv[0-9]+|sd|acf1)[._]", "", names(z))
  z <- z[order(match(obs, unique(obs)), seq_along(z))]
  d <- data.frame(moment = factor(names(z), levels = rev(names(z))),
                  z = as.numeric(z),
                  outside = abs(as.numeric(z)) > stats::qnorm(0.975))
  d$status <- ifelse(d$outside, "|z| > 1.96", "|z| <= 1.96")
  ggplot2::ggplot(d, ggplot2::aes(x = moment, y = z, fill = status)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = c(-1, 1) * stats::qnorm(0.975),
                        linetype = "dashed", colour = dynhr_colours$grey,
                        linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, colour = dynhr_colours$grey, linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("|z| > 1.96" = dynhr_colours$orange,
                                          "|z| <= 1.96" = dynhr_colours$mid_blue),
                               name = NULL) +
    theme_dynhr_diagnostic() +
    ggplot2::labs(
      title = "D29: Data moments vs model moments at theta",
      subtitle = sprintf(paste0("Joint S = %.2f, asymptotically chi2(%d) under the model: p = %.3g\n",
                                "(weak-ID robust; HAC weighting over-rejects in short/persistent samples)"),
                         S, df, p),
      x = NULL, y = "(data - model) / sd(data moment)")
}

#' Compute data moments: variances and autocovariances
#'
#' @param data    T x n_obs matrix
#' @param max_lag Maximum autocovariance lag
#' @return Named numeric vector of moments
#' @noRd
.compute_data_moments <- function(data, max_lag = 4L) {
  T_obs <- nrow(data)
  n_obs <- ncol(data)
  data <- scale(data, scale = FALSE)  # demean

  moments <- c()

  # Variances
  vars <- diag(crossprod(data) / (T_obs - 1))
  names(vars) <- paste0("var_", colnames(data) %||% seq_len(n_obs))
  moments <- c(moments, vars)

  # Autocovariances
  for (lag in seq_len(max_lag)) {
    if (lag >= T_obs) next
    acv <- diag(crossprod(data[(lag + 1):T_obs, , drop = FALSE],
                           data[1:(T_obs - lag), , drop = FALSE])) / (T_obs - 1)
    names(acv) <- paste0("acv", lag, "_", colnames(data) %||% seq_len(n_obs))
    moments <- c(moments, acv)
  }

  moments
}


#' Compute covariance matrix of data moments using HAC or i.i.d. estimator
#'
#' @param data    T x n_obs matrix
#' @param max_lag Maximum lag for HAC truncation
#' @param use_hac Logical: use Newey-West HAC (TRUE) or i.i.d. (FALSE)
#' @return n_mom x n_mom covariance matrix
#' @noRd
.compute_moment_covariance <- function(data, max_lag = 4L, use_hac = TRUE) {
  T_obs <- nrow(data)
  n_obs <- ncol(data)
  data <- scale(data, scale = FALSE)  # demean

  # Build moment vector time series
  # For each t, the "moment observation" is the contribution to the moment vector
  n_mom <- n_obs * (1 + max_lag)  # variances + autocovariances
  moment_ts <- matrix(0, nrow = T_obs, ncol = n_mom)

  col_idx <- 1
  # Variances
  for (j in seq_len(n_obs)) {
    moment_ts[, col_idx] <- data[, j]^2
    col_idx <- col_idx + 1
  }
  # Autocovariances
  for (lag in seq_len(max_lag)) {
    for (j in seq_len(n_obs)) {
      # Contribution at time t: y_{j,t} * y_{j,t-lag}
      # For t <= lag, set to 0
      if (lag < T_obs) {
        moment_ts[(lag + 1):T_obs, col_idx] <- data[(lag + 1):T_obs, j] * data[1:(T_obs - lag), j]
      }
      col_idx <- col_idx + 1
    }
  }

  if (use_hac) {
    # Newey-West HAC estimator
    .newey_west(moment_ts, max_lag = max_lag)
  } else {
    # i.i.d. estimator
    cov(moment_ts, use = "complete.obs")
  }
}


#' Newey-West HAC covariance estimator
#'
#' @param x       T x n matrix
#' @param max_lag Maximum lag for Bartlett kernel truncation
#' @return n x n HAC covariance matrix
#' @noRd
.newey_west <- function(x, max_lag = 4L) {
  T_obs <- nrow(x)
  x_centered <- scale(x, scale = FALSE)
  Gamma0 <- crossprod(x_centered) / T_obs

  hac <- Gamma0
  for (lag in seq_len(max_lag)) {
    if (lag >= T_obs) next
    Gamma_lag <- crossprod(x_centered[(lag + 1):T_obs, , drop = FALSE],
                           x_centered[1:(T_obs - lag), , drop = FALSE]) / T_obs
    weight <- 1 - lag / (max_lag + 1)  # Bartlett kernel
    hac <- hac + weight * (Gamma_lag + t(Gamma_lag))
  }

  hac
}
