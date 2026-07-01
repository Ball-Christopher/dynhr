## R/diag-post-d9-moments.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D9 moment matching; .get_moments_list(), .compute_empirical_moments() helpers
## --------------------------------------------------------------------------

.get_moments_list <- function(moments) {
  if (!is.null(moments$moments)) moments <- moments$moments
  sigma_y <- moments$std_dev
  acf_mat <- NULL
  if (!is.null(moments$autocorr)) {
    ac    <- moments$autocorr
    n_var <- dim(ac)[1]; n_lag <- dim(ac)[3]
    acf_mat <- matrix(NA, n_var, n_lag)
    for (k in seq_len(n_lag)) acf_mat[, k] <- diag(ac[, , k])
    rownames(acf_mat) <- names(sigma_y)
    colnames(acf_mat) <- paste0("lag", seq_len(n_lag))
  }
  xcorr <- moments$correlation
  list(
    std_dev     = sigma_y,  sigma_y    = sigma_y,
    acf_y       = acf_mat,  autocorr   = acf_mat,
    xcorr_y     = xcorr,    correlation = xcorr,
    var_cov     = moments$var_cov
  )
}


.compute_empirical_moments <- function(data, acf_lags = 5) {
  if (!is.matrix(data)) data <- as.matrix(data)
  vnames <- colnames(data)
  n_var  <- ncol(data)
  sds    <- apply(data, 2, sd, na.rm = TRUE)
  names(sds) <- vnames
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  acf_mat <- matrix(NA, n_var, acf_lags,
                    dimnames = list(vnames, paste0("lag", seq_len(acf_lags))))
  for (i in seq_len(n_var)) {
    x <- data[, i][!is.na(data[, i])]
    if (length(x) > acf_lags + 1) {
      ac <- acf(x, lag.max = acf_lags, plot = FALSE)$acf
      acf_mat[i, ] <- ac[2:(acf_lags + 1)]
    }
  }
  list(std_dev = sds, correlation = cor_mat, acf_y = acf_mat)
}


#' D9. Moment matching
#'
#' Compares model-implied standard deviations of observables against the
#' corresponding empirical standard deviations computed from data.  Pass
#' criterion: at least 70\% of observables have a model-to-data SD ratio
#' in [0.5, 2.0] (i.e., within a factor of 2).  Ratios in [0.5, 0.6] or
#' [1.7, 2.0] are flagged as marginal.  \code{pass = NA} when no data or
#' data moments are supplied (self-comparison mode).
#'
#' @param model_moments  Model-implied moment object.  Either a list with
#'   \code{$std_dev} (named numeric vector of standard deviations), or an
#'   object with \code{$moments$std_dev}, or the output of a model solver
#'   that includes standard deviations.
#' @param data           T x n_obs data matrix of observable variables.
#'   When provided and \code{data_moments} is NULL, empirical moments are
#'   computed from \code{data} using \code{.compute_empirical_moments()}.
#' @param data_moments   Pre-computed empirical moment list with \code{$std_dev}.
#'   Takes precedence over \code{data} if both are supplied.
#' @param obs_names      Character vector of observable names to compare.
#'   Defaults to the intersection of model and data observable names.
#' @param metadata       Optional metadata list (deprecated; use \code{meta}).
#' @param acf_lags       Integer: number of autocorrelation lags for empirical
#'   moment computation from \code{data} (default 5).
#' @param key_pairs      Not currently used; reserved for cross-moment checks.
#' @param meta           Optional metadata list for plot annotation.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing \code{sd_table} (model vs data SDs and
#'     ratios), \code{self_compare} (logical), \code{marginal} (variable names
#'     with ratio near the boundary), and \code{outside} (variable names
#'     outside [0.5, 2.0]).}
#'   \item{pass}{Logical -- >= 70\% of observables within [0.5, 2.0] ratio;
#'     \code{NA} in self-comparison mode.}
#'   \item{plots}{ggplot2 bar chart of model vs data standard deviations on
#'     a log scale.}
#'   \item{summary}{Human-readable per-observable summary.}
#'
#' @references Schorfheide, F. (2000). Loss function-based evaluation of DSGE models.
#'   \emph{Journal of Applied Econometrics}, 15(6), 645-670.
#' @noRd
d9_moment_matching <- function(model_moments, data = NULL, data_moments = NULL,
                               obs_names = NULL, metadata = NULL,
                               acf_lags = 5, key_pairs = NULL,
                               meta = NULL) {
    if (!is.null(model_moments$moments)) {
      model_moments <- .get_moments_list(model_moments$moments)
    } else if (!is.null(model_moments$std_dev)) {
      model_moments <- .get_moments_list(model_moments)
    }

    self_compare <- FALSE
    if (!is.null(data) && is.null(data_moments)) {
      if (!is.matrix(data)) data <- as.matrix(data)
      data_moments <- .compute_empirical_moments(data, acf_lags)
    } else if (!is.null(data_moments) && !is.null(data_moments$std_dev)) {
      data_moments <- .get_moments_list(data_moments)
    }
    if (is.null(data_moments)) { data_moments <- model_moments; self_compare <- TRUE }

    all_vars <- names(model_moments$std_dev)
    if (!is.null(obs_names)) {
      obs_names <- intersect(obs_names, intersect(all_vars, names(data_moments$std_dev)))
    } else {
      obs_names <- intersect(all_vars, names(data_moments$std_dev))
    }

    m_sd <- model_moments$std_dev[obs_names]
    d_sd <- data_moments$std_dev[obs_names]
    sd_table <- data.frame(
      variable = obs_names,
      model_sd = as.numeric(m_sd),
      data_sd  = as.numeric(d_sd),
      ratio    = as.numeric(m_sd) / pmax(as.numeric(d_sd), 1e-12),
      stringsAsFactors = FALSE
    )
    sd_table <- sd_table[order(-sd_table$data_sd), ]
    rownames(sd_table) <- NULL

    pass <- if (!self_compare) mean(sd_table$ratio > 0.5 & sd_table$ratio < 2.0) >= 0.7 else NA

    # --- ENHANCED: Add marginal flag for near-boundary observables ---
    in_band     <- sd_table$ratio > 0.5 & sd_table$ratio < 2.0
    marginal    <- (sd_table$ratio > 0.5 & sd_table$ratio < 0.6) |
                   (sd_table$ratio > 1.7 & sd_table$ratio < 2.0)
    outside     <- sd_table$ratio <= 0.5 | sd_table$ratio >= 2.0
    n_marginal  <- sum(marginal & in_band)
    n_outside   <- sum(outside)
    marginal_vars <- sd_table$variable[marginal & in_band]
    outside_vars  <- sd_table$variable[outside]

    note <- if (self_compare) " (self-comparison)" else ""
    summary_parts <- sprintf("D9 Moment matching%s: %d observables", note, length(obs_names))
    if (!self_compare) {
      summary_parts <- paste0(summary_parts, sprintf(
        "\n  Std dev ratios [%.2f, %.2f] -- %d/%d within [0.5,2]  %s%s%s",
        min(sd_table$ratio), max(sd_table$ratio),
        sum(in_band), nrow(sd_table),
        ifelse(pass, "PASS", "FAIL"),
        if (n_marginal > 0) {
          marginal_low  <- sd_table$variable[sd_table$ratio > 0.5 & sd_table$ratio < 0.6]
          marginal_high <- sd_table$variable[sd_table$ratio > 1.7 & sd_table$ratio < 2.0]
          marginal_detail <- c(
            if (length(marginal_low) > 0)
              sprintf("%s [0.5-0.6]", paste(marginal_low, collapse = ", ")),
            if (length(marginal_high) > 0)
              sprintf("%s [1.7-2.0]", paste(marginal_high, collapse = ", "))
          )
          sprintf("  %d marginal: %s", n_marginal,
                  paste(marginal_detail, collapse = "; "))
        } else "",
        if (n_outside > 0)
          sprintf("  %d outside [0.5,2]: %s", n_outside,
                  paste(outside_vars, collapse = ", ")) else ""
      ))
    }

    ratios   <- setNames(sd_table$ratio, sd_table$variable)
    n_obs    <- length(obs_names)

    # Build model-vs-data SD comparison plot (only when real data is available)
    plots_d9 <- list()
    if (!self_compare && requireNamespace("ggplot2", quietly = TRUE) &&
        length(obs_names) > 0 && !all(is.na(d_sd))) {
      sd_long <- data.frame(
        Variable = rep(obs_names, 2L),
        SD       = c(as.numeric(m_sd), as.numeric(d_sd)),
        Source   = rep(c("Model", "Data"), each = length(obs_names)),
        stringsAsFactors = FALSE
      )
      # Use log10 y scale so an outlier (e.g. q_obs) does not compress all
      # other bars. Bars with SD < 1e-6 are clamped to 1e-6 to avoid log(0).
      sd_long$SD <- pmax(sd_long$SD, 1e-6)
      p_sd <- ggplot2::ggplot(sd_long,
                              ggplot2::aes(x = Variable, y = SD, fill = Source)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::scale_y_log10(
          labels = function(x) formatC(x, format = "g", digits = 2)
        ) +
        scale_fill_dynhr_light() +
        theme_dynhr() +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "D9: Model vs data standard deviations (log10 scale)",
          x = NULL, y = "Standard deviation (log scale)"
        )
      plots_d9$moments_sd <- .apply_meta(p_sd, meta)
    }

    # --- ENHANCED: Build per-observable flag maps indexed by variable name ---
    outside_map  <- structure(logical(n_obs), names = obs_names)
    marginal_map <- structure(logical(n_obs), names = obs_names)
    in_band_map  <- structure(logical(n_obs), names = obs_names)
    for (v in obs_names) {
      if (v %in% names(ratios)) {
        r <- ratios[v]
        outside_map[v]  <- r <= 0.5 | r >= 2.0
        marginal_map[v] <- (r > 0.5 & r < 0.6) | (r > 1.7 & r < 2.0)
        in_band_map[v]  <- r > 0.5 & r < 2.0
      }
    }

    structure(
      list(
        result  = list(sd_table = sd_table, self_compare = self_compare,
                        marginal = marginal_vars, outside = outside_vars,
                        marginal_map = marginal_map, outside_map = outside_map),
        pass    = pass,
        plots   = plots_d9,
        summary = summary_parts,
        llm_summary = {
          badge <- if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL"
          ratio_strs <- sprintf("%s: ratio=%.2f%s",
                                obs_names,
                                ratios[obs_names],
                                ifelse(outside_map, " [OUTSIDE]",
                                       ifelse(marginal_map, " [MARGINAL]", "")))
          paste(c(
            sprintf("D9 | Moment Matching | %s", badge),
            sprintf("  observables=%d within_0.5_2.0=%d/%d marginal=%d outside=%d",
                    n_obs, sum(in_band), n_obs, n_marginal, n_outside),
            paste("  sd_ratios (model/data):", paste(ratio_strs, collapse = ", ")),
            sprintf("  action: %s",
                    if (isTRUE(pass))
                      "All model SDs within [0.5, 2.0] x data SDs. Moment fit is reasonable."
                    else if (is.na(pass))
                      "Self-comparison mode -- no data provided for benchmarking."
                    else {
                      over  <- obs_names[outside_map & ratios[obs_names] >= 2.0]
                      under <- obs_names[outside_map & ratios[obs_names] <= 0.5]
                      parts <- character(0)
                      if (length(over)  > 0) parts <- c(parts, paste(sprintf("%s=%.2f", over,  ratios[over]), collapse=", "), "over-predicted")
                      if (length(under) > 0) parts <- c(parts, paste(sprintf("%s=%.2f", under, ratios[under]), collapse=", "), "under-predicted")
                      paste(paste(parts, collapse = "; "),
                            "-- check shock variances or measurement equation.")
                    })
          ), collapse = "\n")
        }
      ),
      class = "dynhr_diagnostic"
    )
}
