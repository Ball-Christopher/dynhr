## R/diag-post-d29-data-constraints.R
## --------------------------------------------------------------------------
## Phase H: D29 — Data-Driven Constraints (Lanne & Luoto 2017)
##
## Implements the Lanne & Luoto (2017) approach for testing identification
## using data moments.  Rather than relying solely on model-implied
## moments, this diagnostic checks whether the observed data actually
## contain enough information to identify the structural parameters.
##
## Algorithm:
##   1. Compute empirical moments from the data (variances, autocovariances,
##      possibly higher-order moments).
##   2. Compute model-implied moments at the calibrated/estimated parameter
##      vector.
##   3. For each parameter direction, compute the "data-identified strength"
##      as the sensitivity of the data moments to that parameter, scaled by
##      the sampling uncertainty of the data moments.
##   4. Flag parameters whose data-identified strength is below threshold.
##   5. Compare data-identified vs model-implied identification (D1/D19) to
##      see which parameters rely on model structure vs. data information.
##
## References:
##   Lanne, M., & Luoto, J. (2017). Data-driven identification of
##     DSGE models. Economics Letters, 158, 35-39.
##   Iskrev, N. (2010). Local identification in DSGE models.
##   Komunjer, I., & Ng, S. (2011). Dynamic identification of DSGE models.
## --------------------------------------------------------------------------

#' D29. Data-Driven Constraints Diagnostic (Lanne & Luoto 2017)
#'
#' Assesses parameter identification using actual data moments rather than
#' purely model-implied quantities.  This diagnostic computes the
#' sensitivity of empirical data moments to each parameter, scaled by the
#' sampling uncertainty of those moments, to determine which parameters
#' are "data-identified."
#'
#' Parameters may be well-identified by the model structure (as in D1/D19)
#' but poorly identified by the data if the data contain little information
#' about them, or vice versa if the model imposes strong cross-equation
#' restrictions that the data help sharpen.
#'
#' @param data            T x n_obs data matrix (observable variables).
#' @param model_solve_fn  Function: theta -> named numeric vector of
#'   model-implied moments.  Should return moments matching those computed
#'   from \code{data}.
#' @param theta           Named numeric vector of parameter values at the
#'   calibration/estimation point.
#' @param param_names     Character vector of parameter names.  If NULL,
#'   derived from \code{theta}.
#' @param moment_names    Character vector of moment names.  If NULL,
#'   derived from the model_solve_fn output.
#' @param max_lag         Maximum autocovariance lag for moment computation
#'   (default 4).
#' @param use_hac         Logical: use HAC (Newey-West) estimator for the
#'   data moment covariance matrix (default TRUE).
#' @param eps             Step size for finite-difference Jacobian (default 1e-5).
#' @param strength_threshold Threshold for data-identified |t|-ratio below
#'   which a parameter is flagged as weakly data-identified (default 2.0).
#' @param verbose         Print progress messages.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing:
#'     \itemize{
#'       \item \code{data_moments} — named numeric vector of empirical moments.
#'       \item \code{model_moments} — model-implied moments at \code{theta}.
#'       \item \code{moment_jacobian} — Jacobian of moments w.r.t. parameters.
#'       \item \code{data_jacobian} — data-scaled Jacobian: J_moment
#'         pre-multiplied by Omega^{-1/2}, the upper-triangular inverse Cholesky
#'         factor of the moment covariance matrix
#'         (\code{t(backsolve(chol(Omega), I))}), with diagonal scaling fallback
#'         when Omega is near-singular.
#'       \item \code{data_ident_strength} — data-identified |t|-ratios.
#'       \item \code{model_ident_strength} — model-implied |t|-ratios
#'         (comparable to D20).
#'       \item \code{moment_cov} — covariance matrix of data moments.
#'       \item \code{weak_data_params} — parameters below threshold.
#'       \item \code{model_vs_data} — comparison table.
#'     }}
#'   \item{pass}{Logical — all parameters have data-identified strength
#'     above threshold.}
#'   \item{plots}{List of ggplot2 objects.}
#'   \item{summary}{Human-readable summary.}
#'
#' @references
#'   Lanne, M., & Luoto, J. (2021). GMM estimation of non-Gaussian structural
#'     vector autoregression. \emph{Journal of Econometrics}, 226(1), 248-270.
#'     (Earlier working paper: 2017.)
#'   Iskrev, N. (2010). Local identification in DSGE models.
#'     \emph{Journal of Monetary Economics}, 57(2), 189-202.
#'   Komunjer, I., & Ng, S. (2011). Dynamic identification of DSGE models.
#'     \emph{Econometrica}, 79(6), 1995-2032.
#'
#' @noRd
d29_data_driven_constraints <- function(data,
                                         model_solve_fn,
                                         theta,
                                         param_names = NULL,
                                         moment_names = NULL,
                                         max_lag = 4L,
                                         use_hac = TRUE,
                                         eps = 1e-5,
                                         strength_threshold = 2.0,
                                         verbose = FALSE,
                                         meta = NULL) {
  # ---- 1. Validate ----
  if (is.null(data) || is.null(model_solve_fn) || is.null(theta)) {
    return(.make_result(
      pass    = NA,
      summary = "D29 Data-Driven Constraints: data, model_solve_fn, and theta are required."
    ))
  }

  if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_along(theta))
  n_par <- length(param_names)
  T_obs <- nrow(data)
  n_obs <- ncol(data)

  # ---- 2. Compute empirical data moments ----
  data_moments <- .compute_data_moments(data, max_lag = max_lag)
  if (is.null(data_moments) || length(data_moments) == 0) {
    return(.make_result(
      pass    = NA,
      summary = "D29 Data-Driven Constraints: Failed to compute data moments."
    ))
  }
  n_mom <- length(data_moments)
  if (is.null(moment_names)) moment_names <- names(data_moments)

  if (verbose) cat(sprintf("[d29] Computed %d moments from %d obs x %d variables.\n",
                           n_mom, T_obs, n_obs))

  # ---- 2b. Quick peek at model moments to detect format mismatch early ----
  # The model_solve_fn may return SDs and ACF1s, while data moments are
  # variances and autocovariances. If so, we'll align them later.
  model_moments_preview <- model_solve_fn(theta)
  model_uses_sd_format <- !is.null(model_moments_preview) &&
    (any(grepl("^sd\\.", names(model_moments_preview))) ||
     any(grepl("^sd_", names(model_moments_preview))))
  model_uses_acf1_format <- !is.null(model_moments_preview) &&
    (any(grepl("^acf1\\.", names(model_moments_preview))) ||
     any(grepl("^acf1_", names(model_moments_preview))))

  # ---- 3. Compute model-implied moments at theta ----
  model_moments <- model_solve_fn(theta)
  if (is.null(model_moments)) {
    return(.make_result(
      pass    = NA,
      summary = "D29 Data-Driven Constraints: model_solve_fn returned NULL at theta."
    ))
  }

  # Ensure same length
  if (length(model_moments) != n_mom) {
    # Try to subset or align by name
    model_named <- !is.null(names(model_moments))
    data_named  <- !is.null(names(data_moments))
    if (model_named && data_named) {
      common <- intersect(names(model_moments), names(data_moments))
      if (length(common) > 0) {
        model_moments <- model_moments[common]
        data_moments <- data_moments[common]
        n_mom <- length(common)
        moment_names <- common
      } else if (model_uses_sd_format && model_uses_acf1_format) {
        # Model uses sd/acf1 format; convert data moments to match
        n_vars <- ncol(data)
        # Data moments are: var_*, acv1_*, acv2_*, acv3_*, acv4_*
        # Model moments are: sd.* (or sd_), acf1.* (or acf1_)
        # Convert: sd = sqrt(var), acf1 = acv1 / var
        data_var <- data_moments[grep("^var_", names(data_moments))]
        data_acv1 <- data_moments[grep("^acv1_", names(data_moments))]
        if (length(data_var) == n_vars && length(data_acv1) == n_vars) {
          data_sd <- sqrt(pmax(data_var, 0))
          data_acf1 <- data_acv1 / pmax(data_var, 1e-16)
          # Rename to match model moment names
          obs_names_data <- sub("^var_", "", names(data_var))
          names(data_sd) <- paste0("sd.", obs_names_data)
          names(data_acf1) <- paste0("acf1.", obs_names_data)
          # Build model-format data moments
          data_moments_converted <- c(data_sd, data_acf1)
          # Align model moments to match
          common <- intersect(names(model_moments), names(data_moments_converted))
          if (length(common) > 0) {
            model_moments <- model_moments[common]
            data_moments <- data_moments_converted[common]
            n_mom <- length(common)
            moment_names <- common
            if (verbose) cat(sprintf("[d29] Converted data moments to sd/acf1 format: %d common moments.\n", n_mom))
          } else {
            # Try name format with underscore instead of dot
            names(data_sd) <- paste0("sd_", obs_names_data)
            names(data_acf1) <- paste0("acf1_", obs_names_data)
            data_moments_converted <- c(data_sd, data_acf1)
            common <- intersect(names(model_moments), names(data_moments_converted))
            if (length(common) > 0) {
              model_moments <- model_moments[common]
              data_moments <- data_moments_converted[common]
              n_mom <- length(common)
              moment_names <- common
              if (verbose) cat(sprintf("[d29] Converted data moments to sd_/acf1_ format: %d common moments.\n", n_mom))
            } else {
              return(.make_result(
                pass = NA,
                summary = sprintf(
                  "D29 Data-Driven Constraints: Moment mismatch. model returns %d moments but data has %d. Could not align sd/acf1 format with data moment names. Model names: %s. Data names: %s.",
                  length(model_moments), n_mom,
                  paste(head(names(model_moments), 4), collapse = ", "),
                  paste(head(names(data_moments), 4), collapse = ", ")
                ),
                errored = TRUE
              ))
            }
          }
        } else {
          return(.make_result(
            pass = NA,
            summary = sprintf(
              "D29 Data-Driven Constraints: Moment mismatch. model returns %d moments but data has %d. Expected %d variance and %d acv1 entries for sd/acf1 conversion.",
              length(model_moments), n_mom, n_vars, n_vars
            ),
            errored = TRUE
          ))
        }
      } else {
        return(.make_result(
          pass = NA,
          summary = sprintf(
            "D29 Data-Driven Constraints: Moment mismatch. model returns %d moments but data has %d. No common moment names. Ensure model_solve_fn returns moments matching the data moment structure (variances + autocovariances at lags 1..%d).",
            length(model_moments), n_mom, max_lag
          ),
          errored = TRUE
        ))
      }
    } else if (!model_named && data_named) {
      # Model moments are unnamed: assign data moment names up to model length
      n_common <- min(length(model_moments), length(data_moments))
      warning(sprintf(
        "d29: model moments are unnamed. Using first %d data moment names as moment_names.",
        n_common))
      model_moments <- model_moments[seq_len(n_common)]
      data_moments  <- data_moments[seq_len(n_common)]
      n_mom <- n_common
      moment_names <- names(data_moments)
      names(model_moments) <- moment_names
    } else {
      # Neither has names: align by position
      n_common <- min(length(model_moments), length(data_moments))
      warning(sprintf(
        "d29: moments are unnamed. Using first %d moments by position.",
        n_common))
      model_moments <- model_moments[seq_len(n_common)]
      data_moments  <- data_moments[seq_len(n_common)]
      n_mom <- n_common
      moment_names <- paste0("m_", seq_len(n_mom))
      names(model_moments) <- moment_names
      names(data_moments)  <- moment_names
    }
  }

  # ---- 4. Numerical Jacobian of moments w.r.t. parameters ----
  J_moment <- .numerical_jacobian(model_solve_fn, theta, eps = eps)
  if (is.null(J_moment)) {
    return(.make_result(
      pass    = NA,
      summary = "D29 Data-Driven Constraints: Failed to compute numerical Jacobian."
    ))
  }
  colnames(J_moment) <- param_names
  rownames(J_moment) <- moment_names %||% paste0("m_", seq_len(nrow(J_moment)))

  # ---- 5. Covariance matrix of data moments (HAC or i.i.d.) ----
  # If moments were converted to sd/acf1 format, compute covariance
  # directly from the aligned moment set.
  if (exists("data_moments_converted", inherits = FALSE) &&
      n_mom < n_obs * (1 + max_lag)) {
    # Build moment time series for the aligned moment types:
    # data SDs and ACF1s (first autocorrelations)
    if (verbose) cat(sprintf("[d29] Recomputing moment covariance for %d aligned moments (sd/acf1).\n", n_mom))
    obs_names_data <- sub("^(sd\\.|sd_|acf1\\.|acf1_)", "", moment_names)
    obs_names_data <- unique(obs_names_data)
    data_demeaned <- scale(data, scale = FALSE)
    T_obs <- nrow(data_demeaned)
    n_aligned <- n_mom
    moment_ts <- matrix(0, nrow = T_obs, ncol = n_aligned)
    colnames(moment_ts) <- moment_names
    for (i in seq_len(n_aligned)) {
      mn <- moment_names[i]
      if (grepl("^sd[\\._]", mn)) {
        obs <- sub("^sd[\\._]", "", mn)
        j <- which(colnames(data_demeaned) == obs)
        if (length(j) == 1) moment_ts[, i] <- data_demeaned[, j]^2
      } else if (grepl("^acf1[\\._]", mn)) {
        obs <- sub("^acf1[\\._]", "", mn)
        j <- which(colnames(data_demeaned) == obs)
        if (length(j) == 1) {
          moment_ts[2:T_obs, i] <- data_demeaned[2:T_obs, j] * data_demeaned[1:(T_obs - 1), j]
        }
      }
    }
    moment_cov <- if (use_hac) .newey_west(moment_ts, max_lag = max_lag) else cov(moment_ts, use = "complete.obs")
  } else {
    moment_cov <- .compute_moment_covariance(data, max_lag = max_lag, use_hac = use_hac)
  }

  # ---- 6. Data-scaled Jacobian ----
  # Standardise: J_data = Omega^{-1/2} * J_moment
  # where Omega is the moment covariance matrix
  # Compute Omega^{-1/2} via Cholesky, fall back to diagonal if singular
  Omega_inv_sqrt <- .robust_Omega_inv_sqrt(moment_cov)

  J_data <- Omega_inv_sqrt %*% J_moment

  # ---- 7. Compute identification strengths ----
  # Adaptive ridge: scale regularisation to the matrix's own diagonal so that
  # near-singular moment matrices (rcond ~ 1e-25) don't crash solve().
  ## Renamed from `.safe_inv`: that name is the PACKAGE-level helper in
  ## R/solve-helpers.R, and a local definition here silently shadowed it for
  ## the rest of this function.
  .d29_ridge_inv <- function(M) {
    ridge <- max(1e-6 * max(abs(diag(M))), 1e-10)
    M_reg <- M + diag(ridge, nrow(M))
    tryCatch(solve(M_reg), error = function(e) {
      sv <- svd(M_reg, nu = 0L, nv = 0L)$d
      sv[sv < .Machine$double.eps * max(sv) * nrow(M)] <- .Machine$double.eps * max(sv) * nrow(M)
      MASS::ginv(M_reg)
    })
  }

  # Model-implied strength (like D20)
  I_model <- crossprod(J_moment)
  I_model_inv <- .d29_ridge_inv(I_model)
  se_model <- sqrt(pmax(diag(I_model_inv), 0))
  model_strength <- abs(theta[param_names]) / pmax(se_model, 1e-16)
  names(model_strength) <- param_names

  # Data-driven strength
  I_data <- crossprod(J_data)
  I_data_inv <- .d29_ridge_inv(I_data)
  se_data <- sqrt(pmax(diag(I_data_inv), 0))
  data_strength <- abs(theta[param_names]) / pmax(se_data, 1e-16)
  names(data_strength) <- param_names

  # ---- 8. Identify weakly data-identified parameters ----
  # Handle NAs in data_strength (from singular/ill-conditioned moment covariance).
  # Parameters with NA strength or strength below threshold are flagged as weak.
  weak_data_params <- param_names[
    is.na(data_strength) | data_strength < strength_threshold
  ]
  # Replace any NA names with actual parameter names
  weak_data_params <- weak_data_params[!is.na(weak_data_params)]

  # Parameters where data identifies better or worse than model
  model_vs_data <- data.frame(
    parameter = param_names,
    model_strength = round(model_strength, 4),
    data_strength = round(data_strength, 4),
    ratio = round(data_strength / pmax(model_strength, 1e-16), 4),
    stringsAsFactors = FALSE
  )

  # ---- 9. Plots ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {

    # Detect whether the model-implied series is available and meaningful.
    # It is absent/uninformative when the Fisher information matrix was singular
    # (se_model ~ 1/1e-16 => model_strength ~ 0). Compare on the DATA scale:
    # if the largest model strength is negligible relative to the data strengths
    # the two-series chart would show invisible model bars, so drop them and use
    # the single-series threshold-coloured chart instead.
    .max_data <- max(data_strength, na.rm = TRUE)
    model_series_ok <- !all(is.na(model_strength)) &&
      is.finite(.max_data) && .max_data > 0 &&
      max(model_strength, na.rm = TRUE) > 0.05 * .max_data

    if (model_series_ok) {
      # Both series available: grouped bar chart coloured by Source
      comp_df <- data.frame(
        Parameter = rep(param_names, 2),
        Source    = rep(c("Model-implied", "Data-driven"), each = n_par),
        Strength  = c(model_strength, data_strength),
        stringsAsFactors = FALSE
      )
      p_sc <- ggplot2::ggplot(
        comp_df, ggplot2::aes(x = Parameter, y = Strength, fill = Source)
      ) +
        ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = strength_threshold,
                           linetype = "dashed", colour = dynhr_colours$red,
                           linewidth = 0.5) +
        ggplot2::scale_fill_manual(
          values = c("Model-implied" = dynhr_colours$mid_blue,
                     "Data-driven"   = dynhr_colours$orange)
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
          title    = "D29: Model-implied vs data-driven identification strength",
          subtitle = sprintf("%d parameters, %d moments, T=%d", n_par, n_mom, T_obs),
          x = NULL, y = "|t|-ratio"
        )
    } else {
      # Model-implied series is unavailable (singular Fisher info / no Jacobian).
      # Drop it entirely and colour each bar by pass / fail instead so the
      # chart is not misleading with an empty or ghost series.
      data_df <- data.frame(
        Parameter  = param_names,
        Strength   = data_strength,
        Identified = data_strength >= strength_threshold,
        stringsAsFactors = FALSE
      )
      # Replace NA strength with 0 for display purposes
      data_df$Strength[is.na(data_df$Strength)] <- 0
      data_df$Identified[is.na(data_df$Identified)] <- FALSE

      p_sc <- ggplot2::ggplot(
        data_df,
        ggplot2::aes(x = Parameter, y = Strength,
                     fill = Identified)
      ) +
        ggplot2::geom_col(colour = "white", linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = strength_threshold,
                           linetype = "dashed", colour = dynhr_colours$red,
                           linewidth = 0.5) +
        ggplot2::scale_fill_manual(
          values = c(`TRUE`  = dynhr_colours$teal,
                     `FALSE` = dynhr_colours$orange),
          labels = c(`TRUE`  = "Identified",
                     `FALSE` = "Weak"),
          name   = NULL
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
          title    = "D29: Data-driven identification strength",
          subtitle = sprintf(
            "%d parameters, %d moments, T=%d  [model-implied series unavailable]",
            n_par, n_mom, T_obs),
          x = NULL, y = "|t|-ratio"
        )
    }
    plots$strength_comparison <- .apply_meta(p_sc, meta)
  }

  # ---- 10. Result ----
  pass <- length(weak_data_params) == 0

  result <- list(
    data_moments           = data_moments,
    model_moments          = model_moments,
    moment_jacobian        = J_moment,
    data_jacobian          = J_data,
    data_ident_strength    = data_strength,
    model_ident_strength   = model_strength,
    moment_cov             = moment_cov,
    weak_data_params       = weak_data_params,
    model_vs_data          = model_vs_data,
    T_obs                  = T_obs,
    n_moments              = n_mom
  )

  weak_str <- if (length(weak_data_params) > 0)
    sprintf("Weak data-identified: %s", paste(weak_data_params, collapse = ", "))
  else "All parameters are data-identified."

  .make_result(
    result  = result,
    pass    = pass,
    plots   = plots,
    summary = sprintf(
      "D29 Data-Driven Constraints: %d params, %d moments, T=%d. %s",
      n_par, n_mom, T_obs, weak_str
    ),
    llm_summary = sprintf(
      "[%s] D29 Data-Driven Constraints n_params=%d n_moments=%d T=%d n_weak=%d",
      if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL",
      n_par, n_mom, T_obs, length(weak_data_params)
    )
  )

}


# ==========================================================================
# Internal helpers for D29
# ==========================================================================

#' Robust Omega^{-1/2} via Cholesky with diagonal fallback
#'
#' @param Omega Covariance matrix
#' @return Omega^{-1/2} matrix
#' @noRd
.robust_Omega_inv_sqrt <- function(Omega) {
  # Check condition number before Cholesky
  if (is.matrix(Omega) && all(is.finite(Omega)) &&
      is.finite(rcond(Omega)) && rcond(Omega) > .Machine$double.eps) {
    R <- chol(Omega)
    # Return R^{-T}: the inverse Cholesky factor, a square-root of the precision
    # J_data = Omega^{-1/2} * J_moment  means we need  t(R^{-1})  = R^{-T}
    return(t(backsolve(R, diag(nrow(R)))))
  }
  # Diagonal fallback
  d <- diag(Omega)
  if (all(d > 0)) {
    return(sqrt(diag(1 / d)))
  }
  # Last resort: regularized identity
  diag(1 / sqrt(pmax(d, 1e-16)))
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
