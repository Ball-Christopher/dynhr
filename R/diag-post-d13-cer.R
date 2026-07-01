## R/diag-post-d13-cer.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
## Updated: D13 cross-equation restrictions â€” full implementation.
##
## D13 tests whether the cross-equation restrictions (CERs) implied by the
## DSGE model's rational-expectations structure are supported by the data.
##
## Approach:
##   1. Build the state-space representation from the decision rules (dr).
##   2. Compute model-implied autocovariances at lags 0..p via Lyapunov.
##   3. Solve the Yule-Walker equations to obtain model-implied VAR(p)
##      coefficients (restricted VAR).
##   4. Estimate an unrestricted VAR(p) from the data via OLS (equation by
##      equation).
##   5. Compare restricted vs. unrestricted coefficients:
##      - Frobenius norm of the difference (scaled by Frobenius norm of
##        the unrestricted coefficients).
##      - Wald-type test using the OLS covariance of the unrestricted VAR.
##   6. Report a pass/fail based on the normalised discrepancy.
##
## References:
##   Del Negro & Schorfheide (2004). Priors from general equilibrium models
##     for VARs. International Economic Review.
##   Fernandez-Villaverde, Rubio-Ramirez & Schorfheide (2016). Solution and
##     estimation methods for DSGE models. Handbook of Macroeconomics.
## --------------------------------------------------------------------------

#' D13. Cross-equation restrictions
#'
#' Tests whether the cross-equation restrictions implied by the DSGE model's
#' rational-expectations structure are consistent with the data. Compares
#' model-implied VAR coefficients (via Yule-Walker from the state-space) to
#' unrestricted OLS-estimated VAR coefficients.
#'
#' @param dr          DecisionRules object from \code{solve_perturbation()}
#'   or \code{stoch_simul()}. Must contain \code{ghx}, \code{ghu},
#'   \code{state_idx}, \code{endo_names}, \code{exo_names}.
#' @param data        Matrix (T x n_obs) of observable time series.
#' @param obs_names   Character vector of observable variable names (must
#'   match column names of \code{data} and entries in \code{dr$endo_names}).
#' @param sigma_e     Optional shock covariance matrix (n_exo x n_exo).
#'   If NULL, computed from \code{model} or assumed identity.
#' @param model       dynhr_mod object (used for shock covariance if
#'   \code{sigma_e} is NULL and compatible with estimation context).
#' @param var_lag     VAR lag order (default 2).
#' @param meta        Optional \code{dynhr_diag_meta} from \code{diag_meta()}.
#' @param wald_p_threshold  p-value threshold below which the null hypothesis
#'   (restrictions hold) is rejected  (default 0.05).
#' @param norm_threshold     Relative Frobenius-norm threshold above which
#'   restrictions are flagged (default 0.5, i.e. 50% discrepancy).
#' @return \code{dynhr_diagnostic} list.
#' @noRd
d13_cross_equation_restrictions <- function(dr            = NULL,
                                             data          = NULL,
                                             obs_names     = NULL,
                                             sigma_e       = NULL,
                                             model         = NULL,
                                             var_lag       = 2L,
                                             meta          = NULL,
                                             wald_p_threshold  = 0.05,
                                             norm_threshold    = 0.5) {

    # ---- 1. Check prerequisites ----
    if (is.null(dr) || is.null(data)) {
      return(.make_result(
        result  = NULL,
        pass    = NA,
        plots   = list(),
        summary = paste(
          "D13 Cross-equation restrictions: provide dr and data.",
          "Build the state-space from decision rules, compute model-implied",
          "autocovariances, solve Yule-Walker for restricted VAR coefficients,",
          "and compare to unrestricted OLS estimates."
        )
      ))
    }

    # ---- 2. Prepare data ----
    if (!is.matrix(data)) data <- as.matrix(data)
    T_obs <- nrow(data)
    if (T_obs <= var_lag + 2L) {
      return(.make_result(
        pass    = NA,
        summary = sprintf(
          "D13: Too few observations (%d) for VAR(%d). Need at least %d.",
          T_obs, var_lag, var_lag + 3L
        )
      ))
    }

    # Determine observable names from data if not provided
    if (is.null(obs_names)) obs_names <- colnames(data)
    if (is.null(obs_names)) {
      n_endo <- nrow(dr$ghx)
      obs_names <- dr$endo_names[seq_len(ncol(data))]
    }

    # Subset data to matching observables
    obs_in_data <- intersect(obs_names, colnames(data))
    if (length(obs_in_data) == 0L) {
      obs_in_data <- colnames(data)[seq_len(min(ncol(data), length(obs_names)))]
    }
    Y <- data[, obs_in_data, drop = FALSE]
    n_obs <- ncol(Y)
    if (n_obs < 1L) {
      return(.make_result(pass = NA, summary = "D13: No valid observables found in data."))
    }

    # ---- 3. Extract state-space from dr ----
    state_idx <- dr$state_idx
    T_mat <- dr$ghx[state_idx, , drop = FALSE]       # n_state x n_state
    n_state <- nrow(T_mat)

    obs_idx <- match(obs_in_data, dr$endo_names)
    if (anyNA(obs_idx)) {
      return(.make_result(pass = NA, summary = "D13: Observables not found in dr$endo_names."))
    }

    Z_mat <- dr$ghx[obs_idx, , drop = FALSE]          # n_obs x n_state
    R_mat <- dr$ghu[state_idx, , drop = FALSE]        # n_state x n_exo
    D_mat <- dr$ghu[obs_idx, , drop = FALSE]          # n_obs x n_exo
    n_exo <- ncol(R_mat)

    # Shock covariance
    if (is.null(sigma_e)) {
      if (!is.null(model)) {
        Sigma_e <- .get_shock_cov(model, dr$exo_names, model$param_values)
      } else {
        Sigma_e <- diag(n_exo)
      }
    } else {
      Sigma_e <- as.matrix(sigma_e)
    }

    # ---- 4. Compute model-implied autocovariances ----
    QQ <- R_mat %*% Sigma_e %*% t(R_mat)
    Sigma_s <- solve_lyapunov(T_mat, QQ)

    # Observable covariance at lag 0
    Gamma_0 <- Z_mat %*% Sigma_s %*% t(Z_mat) + D_mat %*% Sigma_e %*% t(D_mat)

    # Observable cross-covariances at lags 1..var_lag
    Gamma_list <- vector("list", var_lag + 1L)
    Gamma_list[[1L]] <- Gamma_0
    T_power <- diag(n_state)
    for (k in seq_len(var_lag)) {
      T_power <- T_power %*% T_mat
      Gamma_list[[k + 1L]] <- Z_mat %*% T_power %*% Sigma_s %*% t(Z_mat)
    }

    # Build the block-Toeplitz matrix G for Yule-Walker
    G <- matrix(0, nrow = var_lag * n_obs, ncol = var_lag * n_obs)
    for (i in seq_len(var_lag)) {
      for (j in seq_len(var_lag)) {
        lag_idx <- i - j
        if (lag_idx == 0L) {
          block <- Gamma_0
        } else if (lag_idx > 0L) {
          block <- Gamma_list[[lag_idx + 1L]]
        } else {
          block <- t(Gamma_list[[-lag_idx + 1L]])
        }
        row_idx <- ((i - 1L) * n_obs + 1L):(i * n_obs)
        col_idx <- ((j - 1L) * n_obs + 1L):(j * n_obs)
        G[row_idx, col_idx] <- block
      }
    }

    # Right-hand side: vec(Gamma_1, ..., Gamma_p)
    rhs <- matrix(0, nrow = var_lag * n_obs, ncol = n_obs)
    for (k in seq_len(var_lag)) {
      row_idx <- ((k - 1L) * n_obs + 1L):(k * n_obs)
      rhs[row_idx, ] <- t(Gamma_list[[k + 1L]])
    }

    # Solve for model-implied VAR coefficients: A = G^{-1} * rhs
    G_inv <- tryCatch(solve(G), error = function(e) NULL)
    if (is.null(G_inv)) {
      return(.make_result(
        pass    = NA,
        summary = "D13: Cannot invert block-Toeplitz matrix (near-singular)."
      ))
    }
    A_model <- G_inv %*% rhs   # (p*n_obs) x n_obs

    # Reshape: A_model_arr is n_obs x (n_obs * p) for comparison
    A_model_arr <- matrix(0, nrow = n_obs, ncol = n_obs * var_lag)
    for (k in seq_len(var_lag)) {
      row_idx <- ((k - 1L) * n_obs + 1L):(k * n_obs)
      A_model_arr[, row_idx] <- t(A_model[row_idx, , drop = FALSE])
    }

    # Compute model-implied residual covariance
    resid_cov_model <- Gamma_0
    for (k in seq_len(var_lag)) {
      Ak <- A_model_arr[, ((k - 1L) * n_obs + 1L):(k * n_obs)]
      resid_cov_model <- resid_cov_model - Ak %*% t(Gamma_list[[k + 1L]])
    }

    # ---- 5. Estimate unrestricted VAR via OLS ----
    T_eff <- T_obs - var_lag
    Y_lagged <- matrix(NA, nrow = T_eff, ncol = n_obs * var_lag)
    for (k in seq_len(var_lag)) {
      Y_lagged[, ((k - 1L) * n_obs + 1L):(k * n_obs)] <-
        Y[(var_lag - k + 1L):(T_obs - k), , drop = FALSE]
    }
    Y_dep <- Y[(var_lag + 1L):T_obs, , drop = FALSE]

    X <- cbind(1, Y_lagged)
    colnames(X) <- c("const", paste0(rep(obs_in_data, var_lag), "_lag",
                                     rep(seq_len(var_lag), each = n_obs)))

    XtX_inv <- solve(crossprod(X))
    if (is.null(XtX_inv)) {
      return(.make_result(
        pass    = NA,
        summary = "D13: Cannot invert X'X for unrestricted VAR (near-singular)."
      ))
    }
    A_unrest <- XtX_inv %*% crossprod(X, Y_dep)

    resid_unrest <- Y_dep - X %*% A_unrest
    Sigma_u_unrest <- crossprod(resid_unrest) / (T_eff - ncol(X))

    A_unrest_slopes <- A_unrest[-1L, , drop = FALSE]

    # ---- 6. Compare restricted vs unrestricted ----
    diff_mat <- A_unrest_slopes - A_model
    frob_diff <- sqrt(sum(diff_mat^2))
    frob_unrest <- sqrt(sum(A_unrest_slopes^2))
    rel_diff <- if (frob_unrest > 1e-12) frob_diff / frob_unrest else Inf

    # Wald-type test
    n_slopes <- n_obs * var_lag
    XtX_slopes <- XtX_inv[-1L, -1L, drop = FALSE]
    V_ols <- kronecker(Sigma_u_unrest, XtX_slopes)

    diff_vec <- as.vector(diff_mat)
    Wald_stat <- as.numeric(t(diff_vec) %*% solve(V_ols, diff_vec))

    wald_p <- if (!is.na(Wald_stat) && Wald_stat > 0) {
      1 - pchisq(Wald_stat, df = n_slopes * n_obs)
    } else {
      NA_real_
    }

    pass <- if (!is.na(wald_p)) {
      rel_diff < norm_threshold && wald_p > wald_p_threshold
    } else {
      rel_diff < norm_threshold
    }

    # --- ENHANCED: Per-equation Frobenius norm and residual variance ratio ---
    # For each observable equation (column of diff_mat), compute Frobenius norm
    eq_frob <- sqrt(colSums(diff_mat^2))
    names(eq_frob) <- obs_in_data

    # Residual variance ratio: model_var / unrestricted_var per observable
    resid_var_model <- diag(resid_cov_model)
    resid_var_unrest <- diag(Sigma_u_unrest)
    resid_var_ratio <- resid_var_model / pmax(resid_var_unrest, 1e-16)
    names(resid_var_ratio) <- obs_in_data

    # Per-observable R-squared from the unrestricted VAR (how well does data predict each var)
    unrest_r2 <- 1 - resid_var_unrest / diag(stats::var(Y_dep))

    result <- list(
      A_model   = A_model_arr,
      A_unrest  = matrix(A_unrest_slopes, nrow = n_slopes, ncol = n_obs),
      diff_norm_frob = frob_diff,
      rel_diff   = rel_diff,
      Wald_stat  = Wald_stat,
      Wald_p     = wald_p,
      resid_cov_model  = resid_cov_model,
      resid_cov_unrest = Sigma_u_unrest,
      n_obs      = n_obs,
      var_lag    = var_lag,
      T_eff      = T_eff,
      eq_frob    = eq_frob,
      resid_var_ratio = resid_var_ratio,
      unrest_r2  = unrest_r2
    )

    # ---- 7. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .gg <- getNamespace("ggplot2")

      coef_df <- data.frame(
        Lag_coef = rep(rownames(A_unrest_slopes), n_obs),
        Equation = rep(obs_in_data, each = n_slopes),
        Model    = as.vector(A_model),
        Unrest   = as.vector(A_unrest_slopes),
        stringsAsFactors = FALSE
      )
      coef_df$Difference <- coef_df$Unrest - coef_df$Model

      p_coef <- .gg$ggplot(coef_df, .gg$aes(x = Equation, y = Lag_coef)) +
        .gg$geom_tile(.gg$aes(fill = Difference), colour = "white", linewidth = 0.3) +
        .gg$scale_fill_gradient2(low = dynhr_colours$dark_blue,
                                 mid = "white",
                                 high = dynhr_colours$red,
                                 midpoint = 0,
                                 name = "Diff") +
        theme_dynhr(base_size = 10) +
        .gg$theme(axis.text.x = .gg$element_text(angle = 45, hjust = 1)) +
        .gg$labs(title = "D13: VAR coefficient comparison (model vs unrestricted)",
                 subtitle = sprintf("Rel. diff = %.2f | Wald p = %.3f",
                                    rel_diff, if (is.na(wald_p)) -1 else wald_p),
                 x = NULL, y = "Lag coefficient")
      plots$coef_heatmap <- .apply_meta(p_coef, meta)

      rc_df <- data.frame(
        Equation = rep(obs_in_data, 2L),
        Type     = rep(c("Model", "Unrestricted"), each = n_obs),
        Value    = c(diag(resid_cov_model), diag(Sigma_u_unrest)),
        stringsAsFactors = FALSE
      )

      # Identify equations where the model-implied residual variance is negative
      # (mathematically invalid; indicates a near-singular or mis-specified model).
      neg_var_eqs <- obs_in_data[diag(resid_cov_model) < 0]
      p_resid <- .gg$ggplot(rc_df, .gg$aes(x = Equation, y = Value, fill = Type)) +
        .gg$geom_col(position = "dodge", width = 0.7) +
        # Zero reference line makes negative-variance bars immediately visible.
        .gg$geom_hline(yintercept = 0,
                       colour = "grey20", linewidth = 0.5) +
        .gg$scale_fill_manual(values = c("Model" = dynhr_colours$mid_blue,
                                         "Unrestricted" = dynhr_colours$orange)) +
        theme_dynhr_diagnostic() +
        .gg$labs(
          title = "D13: Residual variance comparison (diagonal)",
          subtitle = if (length(neg_var_eqs) > 0)
            sprintf("WARNING: negative model residual variance in: %s (model misspecification)",
                    paste(neg_var_eqs, collapse = ", "))
          else
            "Model vs unrestricted VAR residual variances",
          x = NULL, y = "Residual variance")
      plots$resid_cov <- .apply_meta(p_resid, meta)
    }

    # Identify worst-fit observables (safe indexing for n_obs < 3)
    worst_eq <- names(sort(eq_frob, decreasing = TRUE))
    worst_resid <- names(sort(resid_var_ratio, decreasing = TRUE))
    n_worst <- min(3, n_obs)

    # Flag negative model residual variances — these are mathematically invalid.
    neg_resid_vars <- obs_in_data[diag(resid_cov_model) < 0]
    neg_var_note <- if (length(neg_resid_vars) > 0)
      sprintf(" WARNING: negative model residual variance in %s (misspecification).",
              paste(neg_resid_vars, collapse = ", "))
    else ""

    summary_str <- sprintf(
      "D13 Cross-equation restrictions: VAR(%d) with %d obs x %d vars. %s%s%s%s",
      var_lag, T_eff, n_obs,
      if (!is.na(wald_p)) {
        sprintf("Wald p = %.4f, rel. diff = %.2f. %s",
                wald_p, rel_diff,
                if (isTRUE(pass)) "Restrictions supported by data (PASS)."
                else "Restrictions rejected (FAIL).")
      } else {
        sprintf("Frobenius rel. diff = %.2f. %s",
                rel_diff,
                if (rel_diff < norm_threshold) "Restrictions broadly consistent (PASS)."
                else "Large discrepancy (FAIL).")
      },
      if (!isTRUE(pass) && n_obs > 0) sprintf(
        " Worst-fit eq: %s (eq. frob: %s).",
        paste(worst_eq[1:n_worst], collapse = ", "),
        paste(sprintf("%.2f", eq_frob[worst_eq[1:n_worst]]), collapse = ", ")
      ) else "",
      if (!isTRUE(pass) && n_obs > 0) sprintf(
        " Highest resid var ratio: %s (model/OLS: %s).",
        paste(worst_resid[1:n_worst], collapse = ", "),
        paste(sprintf("%.1f", resid_var_ratio[worst_resid[1:n_worst]]),
              collapse = ", ")
      ) else "",
      neg_var_note
    )

    llm_summary <- paste(c(
      sprintf("D13 | Cross-Equation Restrictions | %s",
              if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "FAIL" else "INFO"),
      sprintf("  VAR(%d) n_obs=%d T=%d", var_lag, n_obs, T_eff),
      sprintf("  rel_frob_diff=%.4f wald_p=%s", rel_diff,
              if (is.na(wald_p)) "NA" else sprintf("%.4f", wald_p)),
      sprintf("  per_eq_frob: %s",
              paste(sprintf("%s=%.2f", names(eq_frob), eq_frob), collapse = ", ")),
      sprintf("  resid_var_ratio(model/OLS): %s",
              paste(sprintf("%s=%.1f", names(resid_var_ratio), resid_var_ratio),
                    collapse = ", ")),
      sprintf("  unrest_r2: %s",
              paste(sprintf("%s=%.2f", names(unrest_r2), unrest_r2), collapse = ", ")),
      sprintf("  worst_equations: %s", paste(worst_eq, collapse = ", ")),
      sprintf("  action: %s",
              if (isTRUE(pass))
                "DSGE cross-equation restrictions are consistent with unrestricted VAR."
              else
                paste0("DSGE restrictions rejected. Key discrepancies in equation(s): ",
                       paste(worst_eq, collapse = ", "),
                       ". These observables have the largest VAR coefficient gap. ",
                       "Consider respecifying price/wage setting, adding shocks, ",
                       "or allowing measurement error for these series."))
    ), collapse = "\n")

    .make_result(
      result      = result,
      pass        = pass,
      plots       = plots,
      summary     = summary_str,
      llm_summary = llm_summary
    )
}
