## R/diag-post-d13-cer.R
## --------------------------------------------------------------------------
## D13 cross-equation restrictions (CER).
##
## D13 tests whether the VAR(p) implied by the solved DSGE model is consistent
## with an unrestricted VAR(p) fitted to the data.
##
## Model (first-order decision rules, deviations from steady state):
##   s_t = T s_{t-1} + R e_t          (T = ghx[state,], R = ghu[state,])
##   y_t = Z s_{t-1} + D e_t          (Z = ghx[obs,],   D = ghu[obs,])
##   e_t ~ (0, Sigma_e)
##
## Approach:
##   1. Sigma_s = Var(s_t) from the Lyapunov equation.
##   2. Observable autocovariances Gamma_k = E[y_t y_{t-k}']:
##        Gamma_0 = Z Sigma_s Z' + D Sigma_e D'
##        Gamma_k = Z T^k Sigma_s Z' + Z T^{k-1} R Sigma_e D'   (k >= 1)
##      (the second term is the covariance between y_t's state part and the
##      contemporaneous shock inside y_{t-k}; it vanishes only when D = 0).
##   3. Population projection of y_t on (y_{t-1}, ..., y_{t-p}) (Yule-Walker):
##        [A_1 ... A_p] = [Gamma_1 ... Gamma_p] G^{-1},
##        G[i, j] = E[y_{t-i} y_{t-j}'] = Gamma_{j-i} (i <= j), Gamma_{i-j}' else
##      Sigma_u = Gamma_0 - sum_k A_k Gamma_k'.
##      This is the pseudo-true value OLS converges to if the model is true,
##      whether or not the model has a finite-order VAR representation.
##   4. Unrestricted VAR(p) with intercept by OLS.
##   5. Compare: relative Frobenius distance and a Wald statistic
##        W = vec(B_ols - B_model)' [Sigma_u_ols (x) (X'X)^{-1}_slopes]^{-1} (...)
##      ~ chi^2(n^2 p) under H0 (model coefficients held fixed).
##
## References:
##   Del Negro & Schorfheide (2004). Priors from general equilibrium models
##     for VARs. International Economic Review.
##   Fernandez-Villaverde, Rubio-Ramirez & Schorfheide (2016). Solution and
##     estimation methods for DSGE models. Handbook of Macroeconomics.
## --------------------------------------------------------------------------

#' Measurement-error variances for D13 from a mode/estimation result
#'
#' Pulls the constant per-observable measurement-error variance the
#' estimation used (\code{mode_result$me_variance}, scalar or one entry per
#' observable) into the named form \code{d13_cross_equation_restrictions()}
#' wants. Returns \code{NULL} when there is no ME, when it is all zero, or
#' when it cannot be matched to \code{obs_names} (D13 then behaves exactly as
#' before). The time-varying \code{me_extra} is deliberately ignored: D13
#' compares stationary second moments, which a per-period ME path does not
#' have.
#' @noRd
.d13_me_var_from_mode <- function(mode_result, obs_names) {
  if (is.null(mode_result) || is.null(obs_names)) return(NULL)
  mv <- mode_result$me_variance
  if (is.null(mv) || !is.numeric(mv) || length(mv) == 0L) return(NULL)
  if (anyNA(mv) || any(!is.finite(mv)) || any(mv < 0)) return(NULL)
  if (all(mv == 0)) return(NULL)
  n_obs <- length(obs_names)
  if (!is.null(names(mv))) {
    keep <- intersect(obs_names, names(mv))   # obs_names order, not mv's
    if (length(keep) == 0L) return(NULL)
    return(mv[keep])
  }
  if (length(mv) == 1L) return(stats::setNames(rep(as.numeric(mv), n_obs), obs_names))
  if (length(mv) == n_obs) return(stats::setNames(as.numeric(mv), obs_names))
  NULL
}


#' D13. Cross-equation restrictions
#'
#' Tests whether the cross-equation restrictions implied by the DSGE model's
#' rational-expectations structure are consistent with the data. The
#' model-implied VAR(p) coefficients are the population Yule-Walker projection
#' of the observables on their own p lags under the first-order solution
#' (including the covariance between the states and the contemporaneous
#' shock loading \code{ghu[obs, ]}); they are compared with an unrestricted
#' OLS VAR(p) (with intercept) by a Wald test with \eqn{n^2 p} degrees of
#' freedom (which sets the PASS/FAIL badge) and a relative Frobenius distance
#' (reported as an effect size).
#'
#' \strong{Measurement error.} Supply \code{me_var} to add classical
#' (white-noise, mutually independent) measurement-error variances to the
#' model-implied \emph{contemporaneous} observable covariance
#' \eqn{\Gamma_0}; the higher-order \eqn{\Gamma_k} are untouched, because
#' white measurement error is uncorrelated across periods. The Yule-Walker
#' projection is then computed from the inflated \eqn{\Gamma_0}, so the
#' model-implied VAR coefficients are attenuated exactly as the data's would
#' be. Without it a model whose only misfit is unmodelled measurement error
#' is penalised: measurement error is standard in DSGE-VAR and Bayesian DSGE
#' estimation (An & Schorfheide 2007), so a "perfect fit, zero ME" comparison
#' is the unrealistic case. \code{me_var = NULL} (the default) reproduces the
#' pre-0.9.4 zero-ME behaviour.
#'
#' The Wald test remains the only statistic that sets the badge. It is a
#' weighted quadratic form in the gap between data and model reduced-form
#' dynamics -- the same family as the GMM/SMM moment-matching criteria used
#' across the DSGE literature (Ruge-Murcia 2007; Christiano, Eichenbaum &
#' Trabandt 2016) -- although the specific OLS-VAR-coefficient covariance
#' \eqn{\Sigma_u \otimes (X'X)^{-1}} is this package's construction rather
#' than a named paper's statistic. No surveyed toolbox (Dynare, IRIS,
#' MacroModelling.jl, RISE) automates a model-vs-data VAR test at all.
#'
#' @param dr          DecisionRules object from \code{solve_perturbation()}
#'   or \code{stoch_simul()}. Must contain \code{ghx}, \code{ghu},
#'   \code{state_idx}, \code{endo_names}.
#' @param data        Matrix (T x n_obs) of observable time series (finite).
#' @param obs_names   Character vector of observable names (entries of
#'   \code{dr$endo_names}). Selects/orders the columns of \code{data} when it
#'   has column names; otherwise it names the columns positionally and must
#'   have length \code{ncol(data)}. Defaults to \code{colnames(data)}.
#' @param sigma_e     Shock covariance matrix (n_exo x n_exo). If NULL it is
#'   taken from \code{model} (\code{shock_cov()}), then from
#'   \code{dr$Sigma_e}; there is no identity fallback.
#' @param model       dynhr_mod object (used for the shock covariance if
#'   \code{sigma_e} is NULL).
#' @param me_var      Optional measurement-error variances, one per
#'   observable: a named numeric vector (names must be entries of
#'   \code{obs_names}; observables not named get 0), or a single unnamed value
#'   applied to every observable. Non-negative and finite. Added to the
#'   diagonal of the model-implied \eqn{\Gamma_0}; see Details. \code{NULL}
#'   (default) means no measurement error.
#' @param var_lag     VAR lag order (default 2).
#' @param meta        Optional \code{dynhr_diag_meta} from \code{diag_meta()}.
#' @param wald_p_threshold  p-value threshold below which the null hypothesis
#'   (restrictions hold) is rejected (default 0.05).
#' @param norm_threshold     Relative Frobenius-norm threshold used for the
#'   badge ONLY when the Wald statistic cannot be computed (default 0.5).
#'   The relative distance is reported always, but it is not a calibrated
#'   test: under a correctly specified model it is dominated by OLS sampling
#'   noise in typical macro samples.
#' @return \code{dynhr_diagnostic} list. \code{result$A_model} and
#'   \code{result$A_unrest} are both \code{n_obs x (n_obs * var_lag)}
#'   (row = equation, columns = lag-1 block, lag-2 block, ...).
#' @noRd
d13_cross_equation_restrictions <- function(dr            = NULL,
                                             data          = NULL,
                                             obs_names     = NULL,
                                             sigma_e       = NULL,
                                             model         = NULL,
                                             me_var        = NULL,
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
    for (fld in c("ghx", "ghu", "state_idx", "endo_names")) {
      if (is.null(dr[[fld]])) .dynhr_abort(sprintf("D13: `dr$%s` is missing.", fld))
    }
    if (length(var_lag) != 1L || !is.finite(var_lag) || var_lag < 1 ||
        var_lag != round(var_lag)) {
      .dynhr_abort("D13: `var_lag` must be a positive integer.")
    }
    var_lag <- as.integer(var_lag)

    # ---- 2. Prepare data ----
    if (!is.matrix(data)) data <- as.matrix(data)
    if (!is.numeric(data)) .dynhr_abort("D13: `data` must be numeric.")
    if (is.null(obs_names)) obs_names <- colnames(data)
    if (is.null(obs_names)) {
      .dynhr_abort("D13: `data` has no column names; supply `obs_names`.")
    }
    if (is.null(colnames(data))) {
      if (length(obs_names) != ncol(data)) {
        .dynhr_abort(sprintf("D13: `data` has no column names and length(obs_names) = %d != ncol(data) = %d.",
                     length(obs_names), ncol(data)))
      }
      colnames(data) <- obs_names
    }
    missing_obs <- setdiff(obs_names, colnames(data))
    if (length(missing_obs) > 0L) {
      .dynhr_abort(sprintf("D13: observables not in colnames(data): %s.",
                   paste(missing_obs, collapse = ", ")))
    }
    obs_in_data <- obs_names
    Y <- data[, obs_in_data, drop = FALSE]
    if (any(!is.finite(Y))) {
      .dynhr_abort("D13: `data` contains non-finite values; the VAR needs a balanced sample.")
    }
    n_obs <- ncol(Y)
    T_obs <- nrow(Y)
    n_reg <- 1L + n_obs * var_lag
    if (T_obs - var_lag <= n_reg) {
      return(.make_result(
        pass    = NA,
        summary = sprintf(
          "D13: Too few observations (%d) for a %d-variable VAR(%d). Need more than %d.",
          T_obs, n_obs, var_lag, var_lag + n_reg
        )
      ))
    }

    # ---- 3. Extract state-space from dr ----
    state_idx <- dr$state_idx
    T_mat <- dr$ghx[state_idx, , drop = FALSE]       # n_state x n_state
    n_state <- nrow(T_mat)

    obs_idx <- match(obs_in_data, dr$endo_names)
    if (anyNA(obs_idx)) {
      .dynhr_abort(sprintf("D13: observables not found in dr$endo_names: %s.",
                   paste(obs_in_data[is.na(obs_idx)], collapse = ", ")))
    }

    Z_mat <- dr$ghx[obs_idx, , drop = FALSE]          # n_obs x n_state
    R_mat <- dr$ghu[state_idx, , drop = FALSE]        # n_state x n_exo
    D_mat <- dr$ghu[obs_idx, , drop = FALSE]          # n_obs x n_exo
    n_exo <- ncol(R_mat)

    # Shock covariance: explicit > model > dr$Sigma_e. No identity default.
    Sigma_e <- if (!is.null(sigma_e)) {
      as.matrix(sigma_e)
    } else if (!is.null(model)) {
      .get_shock_cov(model, dr$exo_names, model$param_values)
    } else if (!is.null(dr$Sigma_e)) {
      as.matrix(dr$Sigma_e)
    } else {
      .dynhr_abort("D13: no shock covariance: supply `sigma_e`, `model`, or a `dr` carrying `Sigma_e`.")
    }
    if (!identical(dim(Sigma_e), c(n_exo, n_exo))) {
      .dynhr_abort(sprintf("D13: shock covariance is %s but dr has %d shocks.",
                   paste(dim(Sigma_e), collapse = "x"), n_exo))
    }
    if (!is.null(rownames(Sigma_e)) && !is.null(dr$exo_names)) {
      if (!setequal(rownames(Sigma_e), dr$exo_names)) {
        .dynhr_abort("D13: rownames(sigma_e) do not match dr$exo_names.")
      }
      Sigma_e <- Sigma_e[dr$exo_names, dr$exo_names, drop = FALSE]
    }
    if (any(!is.finite(Sigma_e))) .dynhr_abort("D13: shock covariance is not finite.")

    # ---- 4. Compute model-implied autocovariances ----
    QQ <- R_mat %*% Sigma_e %*% t(R_mat)
    Sigma_s <- solve_lyapunov(T_mat, QQ)
    RSD <- R_mat %*% Sigma_e %*% t(D_mat)             # Cov(s_t, D e_t)

    Gamma_0 <- Z_mat %*% Sigma_s %*% t(Z_mat) + D_mat %*% Sigma_e %*% t(D_mat)

    # Classical measurement error: white, so it inflates Gamma_0 only. The
    # Gamma_k (k >= 1) are unchanged, which is what attenuates the
    # model-implied VAR coefficients through the Yule-Walker solve below.
    me_vec <- stats::setNames(numeric(n_obs), obs_in_data)
    if (!is.null(me_var)) {
      if (!is.numeric(me_var) || anyNA(me_var) || any(!is.finite(me_var)))
        .dynhr_abort("D13: `me_var` must be finite and numeric.")
      if (any(me_var < 0))
        .dynhr_abort("D13: `me_var` must be non-negative (it is a variance).")
      if (is.null(names(me_var))) {
        if (length(me_var) == 1L) {
          me_vec[] <- me_var
        } else if (length(me_var) == n_obs) {
          me_vec[] <- me_var
        } else {
          .dynhr_abort(sprintf(
            "D13: unnamed `me_var` has length %d; supply 1, %d (in obs_names order), or a named vector.",
            length(me_var), n_obs))
        }
      } else {
        unknown <- setdiff(names(me_var), obs_in_data)
        if (length(unknown) > 0L)
          .dynhr_abort(sprintf("D13: `me_var` names not among the observables: %s.",
                               paste(unknown, collapse = ", ")))
        me_vec[names(me_var)] <- me_var
      }
    }
    Gamma_0 <- Gamma_0 + diag(me_vec, nrow = n_obs)

    # Gamma_list[[k + 1]] = E[y_t y_{t-k}']
    Gamma_list <- vector("list", var_lag + 1L)
    Gamma_list[[1L]] <- Gamma_0
    T_prev <- diag(n_state)                           # T^{k-1}
    for (k in seq_len(var_lag)) {
      T_k <- T_prev %*% T_mat
      Gamma_list[[k + 1L]] <- Z_mat %*% (T_k %*% Sigma_s %*% t(Z_mat) + T_prev %*% RSD)
      T_prev <- T_k
    }

    # G[i, j] = E[y_{t-i} y_{t-j}']
    blk_idx <- function(i) ((i - 1L) * n_obs + 1L):(i * n_obs)
    G <- matrix(0, nrow = var_lag * n_obs, ncol = var_lag * n_obs)
    for (i in seq_len(var_lag)) {
      for (j in seq_len(var_lag)) {
        G[blk_idx(i), blk_idx(j)] <- if (i <= j) Gamma_list[[j - i + 1L]]
                                     else t(Gamma_list[[i - j + 1L]])
      }
    }
    # C[i-block] = E[y_t y_{t-i}'] = Gamma_i
    C_mat <- do.call(cbind, Gamma_list[-1L])          # n_obs x (n_obs p)

    if (any(!is.finite(G)) || rcond(G) < 1e-12) {
      return(.make_result(
        pass    = NA,
        summary = sprintf(
          "D13: model-implied lag covariance matrix is singular (rcond = %.2e); the observables are (near-)collinear under the model (stochastic singularity).",
          if (all(is.finite(G))) rcond(G) else NaN)
      ))
    }
    A_model_arr <- t(solve(G, t(C_mat)))              # n_obs x (n_obs p), G symmetric
    resid_cov_model <- Gamma_0 - A_model_arr %*% t(C_mat)
    resid_cov_model <- (resid_cov_model + t(resid_cov_model)) / 2

    coef_names <- paste0(rep(obs_in_data, var_lag), "_lag",
                         rep(seq_len(var_lag), each = n_obs))
    dimnames(A_model_arr) <- list(obs_in_data, coef_names)
    dimnames(resid_cov_model) <- list(obs_in_data, obs_in_data)

    # ---- 5. Estimate unrestricted VAR via OLS ----
    T_eff <- T_obs - var_lag
    Y_lagged <- matrix(NA_real_, nrow = T_eff, ncol = n_obs * var_lag)
    for (k in seq_len(var_lag)) {
      Y_lagged[, blk_idx(k)] <- Y[(var_lag - k + 1L):(T_obs - k), , drop = FALSE]
    }
    Y_dep <- Y[(var_lag + 1L):T_obs, , drop = FALSE]

    X <- cbind(1, Y_lagged)
    colnames(X) <- c("const", coef_names)
    XtX <- crossprod(X)
    if (rcond(XtX) < 1e-12) {
      return(.make_result(
        pass    = NA,
        summary = sprintf("D13: X'X of the unrestricted VAR is singular (rcond = %.2e); collinear or constant data.",
                          rcond(XtX))
      ))
    }
    XtX_inv <- solve(XtX)
    B_ols <- XtX_inv %*% crossprod(X, Y_dep)          # n_reg x n_obs

    resid_unrest <- Y_dep - X %*% B_ols
    Sigma_u_unrest <- crossprod(resid_unrest) / (T_eff - n_reg)
    dimnames(Sigma_u_unrest) <- list(obs_in_data, obs_in_data)

    A_unrest_arr <- t(B_ols[-1L, , drop = FALSE])     # n_obs x (n_obs p)
    dimnames(A_unrest_arr) <- list(obs_in_data, coef_names)
    se_arr <- sqrt(outer(diag(Sigma_u_unrest), diag(XtX_inv)[-1L]))
    dimnames(se_arr) <- dimnames(A_unrest_arr)

    # ---- 6. Compare restricted vs unrestricted ----
    diff_mat <- A_unrest_arr - A_model_arr            # n_obs x (n_obs p)
    frob_diff <- sqrt(sum(diff_mat^2))
    frob_unrest <- sqrt(sum(A_unrest_arr^2))
    rel_diff <- if (frob_unrest > 1e-12) frob_diff / frob_unrest else Inf

    # Wald: Cov(vec(B_slopes)) = Sigma_u (x) (X'X)^{-1}_slopes, vec column-major
    # over equations, i.e. vec(t(diff_mat)).
    df_wald <- n_obs * n_obs * var_lag
    V_ols <- kronecker(Sigma_u_unrest, XtX_inv[-1L, -1L, drop = FALSE])
    diff_vec <- as.vector(t(diff_mat))
    Wald_stat <- if (rcond(V_ols) > 1e-14) {
      as.numeric(crossprod(diff_vec, solve(V_ols, diff_vec)))
    } else {
      NA_real_
    }
    wald_p <- if (is.finite(Wald_stat)) {
      stats::pchisq(Wald_stat, df = df_wald, lower.tail = FALSE)
    } else {
      NA_real_
    }

    # The badge is the calibrated Wald test. rel_diff is an effect size whose
    # null distribution depends on T and on the persistence of the data (at
    # T = 250 it exceeds 0.5 in most samples drawn from the model itself), so
    # it only decides the badge when the Wald statistic is unavailable.
    pass <- if (!is.na(wald_p)) {
      wald_p > wald_p_threshold
    } else if (is.finite(rel_diff)) {
      rel_diff < norm_threshold
    } else {
      NA
    }

    eq_frob <- sqrt(rowSums(diff_mat^2))
    names(eq_frob) <- obs_in_data
    # Per-equation Wald (Sigma_u[m, m] * (X'X)^{-1}_slopes)
    eq_wald <- vapply(seq_len(n_obs), function(m) {
      d <- diff_mat[m, ]
      as.numeric(crossprod(d, solve(XtX_inv[-1L, -1L, drop = FALSE], d))) /
        Sigma_u_unrest[m, m]
    }, numeric(1))
    eq_wald_p <- stats::pchisq(eq_wald, df = n_obs * var_lag, lower.tail = FALSE)
    names(eq_wald) <- names(eq_wald_p) <- obs_in_data

    resid_var_model <- diag(resid_cov_model)
    resid_var_unrest <- diag(Sigma_u_unrest)
    resid_var_ratio <- resid_var_model / pmax(resid_var_unrest, 1e-16)
    names(resid_var_ratio) <- obs_in_data

    unrest_r2 <- 1 - colSums(resid_unrest^2) /
      colSums(sweep(Y_dep, 2L, colMeans(Y_dep))^2)
    names(unrest_r2) <- obs_in_data

    result <- list(
      A_model   = A_model_arr,
      A_unrest  = A_unrest_arr,
      A_unrest_se = se_arr,
      diff_norm_frob = frob_diff,
      rel_diff   = rel_diff,
      Wald_stat  = Wald_stat,
      Wald_df    = df_wald,
      Wald_p     = wald_p,
      eq_wald    = eq_wald,
      eq_wald_p  = eq_wald_p,
      resid_cov_model  = resid_cov_model,
      resid_cov_unrest = Sigma_u_unrest,
      Sigma_e    = Sigma_e,
      me_var     = me_vec,
      n_obs      = n_obs,
      var_lag    = var_lag,
      T_eff      = T_eff,
      eq_frob    = eq_frob,
      resid_var_ratio = resid_var_ratio,
      unrest_r2  = unrest_r2
    )

    wald_txt <- if (is.na(wald_p)) "Wald n/a" else sprintf("Wald = %.1f (df %d), p = %.3g",
                                                           Wald_stat, df_wald, wald_p)

    # ---- 7. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .gg <- getNamespace("ggplot2")

      coef_df <- data.frame(
        Equation = factor(rep(obs_in_data, times = n_obs * var_lag),
                          levels = obs_in_data),
        Coef     = factor(rep(coef_names, each = n_obs), levels = rev(coef_names)),
        Model    = as.vector(A_model_arr),
        OLS      = as.vector(A_unrest_arr),
        SE       = as.vector(se_arr),
        stringsAsFactors = FALSE
      )
      coef_df$lo <- coef_df$OLS - 1.96 * coef_df$SE
      coef_df$hi <- coef_df$OLS + 1.96 * coef_df$SE
      pt_df <- rbind(
        data.frame(Equation = coef_df$Equation, Coef = coef_df$Coef,
                   Value = coef_df$OLS, Source = "Unrestricted OLS (95% CI)"),
        data.frame(Equation = coef_df$Equation, Coef = coef_df$Coef,
                   Value = coef_df$Model, Source = "Model-implied")
      )
      pal <- c("Unrestricted OLS (95% CI)" = dynhr_colours$orange,
               "Model-implied"             = dynhr_colours$mid_blue)
      shp <- c("Unrestricted OLS (95% CI)" = 16, "Model-implied" = 4)

      p_coef <- .gg$ggplot(coef_df, .gg$aes(y = Coef)) +
        .gg$geom_vline(xintercept = 0, colour = dynhr_colours$grey, linewidth = 0.4) +
        .gg$geom_errorbar(.gg$aes(xmin = lo, xmax = hi), width = 0.3,
                          orientation = "y",
                          colour = dynhr_colours$orange, linewidth = 0.5) +
        .gg$geom_point(data = pt_df,
                       .gg$aes(x = Value, colour = Source, shape = Source),
                       size = 2.4, stroke = 1.1) +
        .gg$scale_colour_manual(values = pal, name = NULL) +
        .gg$scale_shape_manual(values = shp, name = NULL) +
        .gg$facet_wrap(~Equation, scales = "free_x",
                       labeller = .gg$labeller(Equation = function(x) paste("Equation:", x))) +
        theme_dynhr(base_size = 11) +
        .gg$labs(title = sprintf("D13: VAR(%d) coefficients, model-implied vs unrestricted", var_lag),
                 subtitle = sprintf("Rel. Frobenius diff = %.2f | %s (reject if p < %.2f) | bars: OLS +/- 1.96 s.e.",
                                    rel_diff, wald_txt, wald_p_threshold),
                 x = "Coefficient on regressor", y = "Regressor (variable_lag)")
      plots$coef_compare <- .apply_meta(p_coef, meta)

      rc_df <- data.frame(
        Equation = factor(rep(obs_in_data, 2L), levels = obs_in_data),
        Type     = rep(c("Model-implied", "Unrestricted OLS"), each = n_obs),
        Value    = c(resid_var_model, resid_var_unrest),
        stringsAsFactors = FALSE
      )

      neg_var_eqs <- obs_in_data[resid_var_model < 0]
      p_resid <- .gg$ggplot(rc_df, .gg$aes(x = Type, y = Value, fill = Type)) +
        .gg$geom_col(width = 0.7) +
        .gg$facet_wrap(~Equation, scales = "free_y") +
        .gg$geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.5) +
        .gg$scale_fill_manual(values = c("Model-implied" = dynhr_colours$mid_blue,
                                         "Unrestricted OLS" = dynhr_colours$orange),
                              name = NULL) +
        theme_dynhr(base_size = 11) +
        .gg$theme(axis.text.x = .gg$element_blank(),
                  axis.ticks.x = .gg$element_blank()) +
        .gg$labs(
          title = sprintf("D13: VAR(%d) one-step residual variance by equation", var_lag),
          subtitle = if (length(neg_var_eqs) > 0)
            sprintf("WARNING: negative model residual variance in: %s",
                    paste(neg_var_eqs, collapse = ", "))
          else
            "Model-implied innovation variance vs OLS residual variance (squared data units; per-panel scale)",
          x = NULL, y = "Residual variance")
      plots$resid_cov <- .apply_meta(p_resid, meta)
    }

    ord <- order(eq_wald_p)
    worst_eq <- obs_in_data[ord]
    worst_resid <- names(sort(resid_var_ratio, decreasing = TRUE))
    n_worst <- min(3L, n_obs)

    neg_resid_vars <- obs_in_data[resid_var_model < 0]
    neg_var_note <- if (length(neg_resid_vars) > 0)
      sprintf(" WARNING: negative model residual variance in %s.",
              paste(neg_resid_vars, collapse = ", "))
    else ""

    verdict <- if (isTRUE(pass)) "Restrictions supported by data (PASS)."
               else "Restrictions rejected (FAIL)."
    summary_str <- sprintf(
      "D13 Cross-equation restrictions: VAR(%d), %d obs x %d vars. %s, rel. diff = %.2f. %s%s%s%s",
      var_lag, T_eff, n_obs, wald_txt, rel_diff, verdict,
      if (!isTRUE(pass)) sprintf(
        " Worst-fit eq (per-eq Wald p): %s.",
        paste(sprintf("%s (%.3g)", worst_eq[seq_len(n_worst)],
                      eq_wald_p[worst_eq[seq_len(n_worst)]]), collapse = ", ")
      ) else "",
      if (!isTRUE(pass)) sprintf(
        " Highest resid var ratio (model/OLS): %s.",
        paste(sprintf("%s (%.2f)", worst_resid[seq_len(n_worst)],
                      resid_var_ratio[worst_resid[seq_len(n_worst)]]),
              collapse = ", ")
      ) else "",
      neg_var_note
    )

    llm_summary <- paste(c(
      sprintf("D13 | Cross-Equation Restrictions | %s",
              if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "FAIL" else "INFO"),
      sprintf("  VAR(%d) n_obs=%d T=%d me_var=%s", var_lag, n_obs, T_eff,
              if (all(me_vec == 0)) "none"
              else paste(sprintf("%s=%.3g", names(me_vec), me_vec), collapse = " ")),
      sprintf("  rel_frob_diff=%.4f wald=%s df=%d wald_p=%s", rel_diff,
              if (is.na(Wald_stat)) "NA" else sprintf("%.2f", Wald_stat), df_wald,
              if (is.na(wald_p)) "NA" else sprintf("%.4g", wald_p)),
      sprintf("  per_eq_wald_p: %s",
              paste(sprintf("%s=%.3g", names(eq_wald_p), eq_wald_p), collapse = ", ")),
      sprintf("  per_eq_frob: %s",
              paste(sprintf("%s=%.2f", names(eq_frob), eq_frob), collapse = ", ")),
      sprintf("  resid_var_ratio(model/OLS): %s",
              paste(sprintf("%s=%.2f", names(resid_var_ratio), resid_var_ratio),
                    collapse = ", ")),
      sprintf("  unrest_r2: %s",
              paste(sprintf("%s=%.2f", names(unrest_r2), unrest_r2), collapse = ", ")),
      sprintf("  worst_equations: %s", paste(worst_eq, collapse = ", ")),
      sprintf("  action: %s",
              if (isTRUE(pass))
                "DSGE cross-equation restrictions are consistent with unrestricted VAR."
              else
                paste0("DSGE restrictions rejected. Largest per-equation discrepancies: ",
                       paste(worst_eq[seq_len(n_worst)], collapse = ", "),
                       ". Consider respecifying the dynamics feeding these series, ",
                       "adding shocks, or",
                       if (all(me_vec == 0))
                         " passing the estimated measurement-error variances as me_var (none were supplied here)."
                       else " raising the measurement-error variances in me_var."))
    ), collapse = "\n")

    .make_result(
      result      = result,
      pass        = pass,
      plots       = plots,
      summary     = summary_str,
      llm_summary = llm_summary
    )
}
