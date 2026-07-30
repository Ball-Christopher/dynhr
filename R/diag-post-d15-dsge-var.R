## R/diag-post-d15-dsge-var.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
## Updated: D15 DSGE-VAR -- full implementation.
##
## D15 implements the DSGE-VAR approach of Del Negro & Schorfheide (2004),
## which uses the DSGE model to construct a prior for a VAR and then
## compares the restricted (lambda -> Inf) and unrestricted (lambda -> 0) estimates by
## optimising over the tightness parameter lambda.
##
## Approach:
##   1. Build the state-space from the decision rules.
##   2. Compute model-implied autocovariances -> Yule-Walker VAR coefficients
##      and residual covariance (the "DSGE prior").
##   3. For each lambda in a grid {lambda_1, ..., lambda_K}:
##      a. Construct lambda-scaled dummy observations from the DSGE prior.
##      b. Augment actual data with dummy observations.
##      c. Estimate VAR coefficients via OLS on augmented data.
##      d. Compute log marginal likelihood (Laplace / BIC approximation).
##   4. Report the optimal lambda* that maximises the marginal likelihood.
##   5. Interpret lambda*: lambda* ~= Inf -> DSGE restrictions supported;
##      lambda* ~= 0 -> DSGE restrictions rejected.
##
## References:
##   Del Negro & Schorfheide (2004). Priors from general equilibrium models
##     for VARs. International Economic Review, 45(2), 643-673.
##   Del Negro, Schorfheide, Smets & Wouters (2007). On the fit of
##     New Keynesian models. J. Business & Economic Statistics, 25(2).
## --------------------------------------------------------------------------

#' D15. DSGE-VAR tightness (lambda)
#'
#' Estimates the DSGE-VAR model of Del Negro & Schorfheide (2004) and
#' reports the optimal tightness parameter lambda.  Lambda = Inf recovers
#' the DSGE restrictions exactly; lambda near zero favours the unrestricted
#' VAR.  The procedure searches over a grid of lambda values and selects
#' the one maximising the log marginal likelihood.
#'
#' @param dr          DecisionRules object from \code{solve_perturbation()}
#'   or \code{stoch_simul()}. Must contain \code{ghx}, \code{ghu},
#'   \code{state_idx}, \code{endo_names}, \code{exo_names}.
#' @param data        Matrix (T x n_obs) of observable time series.
#' @param obs_names   Character vector of observable variable names.
#' @param sigma_e     Optional shock covariance matrix. If NULL, computed
#'   from \code{model} or assumed identity.
#' @param model       dynhr_mod object (for shock covariance if needed).
#' @param var_lag     VAR lag order (default 2).
#' @param lambda_grid Numeric vector of lambda values to search over.
#'   Default: \code{c(0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10)}.
#' @param meta        Optional \code{dynhr_diag_meta} from \code{diag_meta()}.
#' @return \code{dynhr_diagnostic} list with optimal lambda and comparison.
#' @noRd
d15_dsge_var <- function(dr            = NULL,
                          data          = NULL,
                          obs_names     = NULL,
                          sigma_e       = NULL,
                          model         = NULL,
                          var_lag       = 2L,
                          lambda_grid   = c(0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10),
                          meta          = NULL) {
  
    # ---- 1. Check prerequisites ----
    if (is.null(dr) || is.null(data)) {
      return(.make_result(
        result  = NULL,
        pass    = NA,
        plots   = list(),
        summary = paste(
          "D15 DSGE-VAR: provide dr and data.",
          "Computes the optimal DSGE-VAR tightness lambda by comparing",
          "model-implied prior moments to unrestricted VAR estimates."
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
          "D15: Too few observations (%d) for VAR(%d).", T_obs, var_lag
        )
      ))
    }

    if (is.null(obs_names)) obs_names <- colnames(data)
    obs_in_data <- intersect(obs_names, colnames(data))
    if (length(obs_in_data) < 2L) {
      obs_in_data <- colnames(data)[seq_len(min(ncol(data), length(obs_names)))]
    }
    Y <- data[, obs_in_data, drop = FALSE]
    n_obs <- ncol(Y)

    # ---- 3. Build DSGE prior (model-implied VAR coefficients + residual cov) ----
    state_idx <- dr$state_idx
    T_mat <- dr$ghx[state_idx, , drop = FALSE]
    n_state <- nrow(T_mat)

    obs_idx <- match(obs_in_data, dr$endo_names)
    if (anyNA(obs_idx)) {
      return(.make_result(pass = NA, summary = "D15: Observables not found in dr$endo_names."))
    }

    Z_mat <- dr$ghx[obs_idx, , drop = FALSE]
    R_mat <- dr$ghu[state_idx, , drop = FALSE]
    D_mat <- dr$ghu[obs_idx, , drop = FALSE]
    n_exo <- ncol(R_mat)

    if (is.null(sigma_e)) {
      if (!is.null(model)) {
        Sigma_e <- .get_shock_cov(model, dr$exo_names, model$param_values)
      } else {
        Sigma_e <- diag(n_exo)
      }
    } else {
      Sigma_e <- as.matrix(sigma_e)
    }

    # Compute model-implied autocovariances (as in D13)
    QQ <- R_mat %*% Sigma_e %*% t(R_mat)
    Sigma_s <- solve_lyapunov(T_mat, QQ)

    Gamma_0 <- Z_mat %*% Sigma_s %*% t(Z_mat) + D_mat %*% Sigma_e %*% t(D_mat)

    Gamma_list <- vector("list", var_lag + 1L)
    Gamma_list[[1L]] <- Gamma_0
    T_power <- diag(n_state)
    for (k in seq_len(var_lag)) {
      T_power <- T_power %*% T_mat
      Gamma_list[[k + 1L]] <- Z_mat %*% T_power %*% Sigma_s %*% t(Z_mat)
    }

    # Yule-Walker for DSGE prior
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

    rhs <- matrix(0, nrow = var_lag * n_obs, ncol = n_obs)
    for (k in seq_len(var_lag)) {
      row_idx <- ((k - 1L) * n_obs + 1L):(k * n_obs)
      rhs[row_idx, ] <- t(Gamma_list[[k + 1L]])
    }

    G_rcond <- tryCatch(rcond(G), error = function(e) 0)
    if (!is.finite(G_rcond) || G_rcond < 1e-14)
      return(.make_result(pass = NA,
        summary = sprintf("D15: Block-Toeplitz near-singular (rcond=%.2e). With %d obs > %d shocks the DSGE autocovariance matrix is rank-deficient; increase measurement-error variance or reduce observables.", G_rcond, n_obs, n_obs)))
    G_inv <- tryCatch(solve(G), error = function(e) NULL)
    if (is.null(G_inv)) return(.make_result(
      pass = NA, summary = "D15: Cannot invert block-Toeplitz (near-singular)."
    ))
    A_prior <- G_inv %*% rhs

    # Residual covariance from DSGE prior
    Sigma_u_prior <- Gamma_0
    for (k in seq_len(var_lag)) {
      Ak <- A_prior[((k - 1L) * n_obs + 1L):(k * n_obs), , drop = FALSE]
      Sigma_u_prior <- Sigma_u_prior - t(Ak) %*% t(Gamma_list[[k + 1L]])
    }
    Sigma_u_prior <- (Sigma_u_prior + t(Sigma_u_prior)) / 2  # symmetrise

    # Cholesky of prior residual covariance
    chol_Sigma <- tryCatch(chol(Sigma_u_prior), error = function(e) NULL)
    if (is.null(chol_Sigma)) {
      # Identify which observables contribute most to non-PD by inspecting
      # the minimum eigenvalue of each observable's diagonal block.
      non_pd_hint <- tryCatch({
        ev <- eigen(Sigma_u_prior, symmetric = TRUE, only.values = TRUE)$values
        neg_idx <- which(ev < 0)
        if (length(neg_idx) > 0) {
          sprintf(" Min eigenvalue = %.3e. Check measurement error specification for: %s.",
                  min(ev),
                  paste(obs_in_data[order(diag(Sigma_u_prior))][seq_len(min(3, n_obs))],
                        collapse = ", "))
        } else ""
      }, error = function(e) "")
      return(.make_result(pass = NA, summary = paste0(
        "D15: DSGE prior residual covariance is not positive definite.",
        non_pd_hint,
        " Action: add measurement error (me_var) for the named observables,",
        " or check that shock covariance Sigma_e is correctly specified."
      )))
    }

    # ---- 4. Build actual data matrices for VAR estimation ----
    T_eff <- T_obs - var_lag
    Y_lagged <- matrix(NA, nrow = T_eff, ncol = n_obs * var_lag)
    for (k in seq_len(var_lag)) {
      Y_lagged[, ((k - 1L) * n_obs + 1L):(k * n_obs)] <-
        Y[(var_lag - k + 1L):(T_obs - k), , drop = FALSE]
    }
    Y_dep <- Y[(var_lag + 1L):T_obs, , drop = FALSE]
    X_data <- cbind(1, Y_lagged)  # with intercept
    k_coef <- 1L + n_obs * var_lag  # number of coefficients per equation

    # ---- 5. Search over lambda grid ----
    lambda_grid <- sort(lambda_grid)
    n_lambda <- length(lambda_grid)
    log_ml <- rep(NA_real_, n_lambda)
    A_lambda_list <- vector("list", n_lambda)
    Sigma_u_lambda_list <- vector("list", n_lambda)

    # Pre-compute data moments needed for all lambdas
    XtX_data <- crossprod(X_data)
    XtY_data <- crossprod(X_data, Y_dep)
    YtY_data <- crossprod(Y_dep)

    for (li in seq_len(n_lambda)) {
      lam <- lambda_grid[li]

      # Construct dummy observations: T_dummy = lam * n_obs (scaled by lambda)
      # Following Del Negro & Schorfheide: dummy obs are constructed from the
      # prior moments such that the posterior = (1-lambda)*likelihood + lambda*prior
      # We use the formulation where dummy counts = lambda * T_eff
      T_dummy <- max(1L, round(lam * n_obs))

      # Dummy Y: Y_dummy has covariance Sigma_u_prior (scaled)
      # Dummy X: X_dummy such that beta_prior = (X_dummy'X_dummy)^{-1} X_dummy'Y_dummy
      # We construct: Y_dummy = T_dummy * chol_Sigma / sqrt(T_dummy)... 
      # Actually, the proper construction:
      # Y_dummy = matrix with T_dummy rows drawn so that:
      #   Y_dummy'Y_dummy = T_dummy * Sigma_u_prior  (prior sum of squares)
      #   X_dummy'X_dummy = T_dummy * Gamma_xx       (prior cross of regressors)
      #   X_dummy'Y_dummy = T_dummy * Gamma_xy       (prior cross)
      # where Gamma_xx = [1 0; 0 G] and Gamma_xy = [0; rhs]
      # 
      # Simpler: just use the prior moments directly as artificial observations.
      # The dummy observations are sqrt(lambda) * T of artificial data that match the
      # DSGE prior moments.

      # Build dummy observations from the prior VAR structure
      # X_dummy'X_dummy = lambda * T * diag(1, G)  (intercept + lag coefficients)
      # But since G might have off-diagonals, use the full structure
      
      # Following DSGE-Var literature more closely:
      # The DSGE prior is: beta|Sigma ~ N(beta_DSGE(theta), Sigma (x) (lambdaT * Gamma_xx)^{-1})  
      # and Sigma ~ IW(lambdaT * Sigma_u(theta), lambdaT - k)
      # where Gamma_xx = E[x_t x_t'] and Gamma_xy = E[x_t y_t']
      #
      # We set up dummy observations:
      # Y_dummy = zeros(T_dummy, n_obs)  -- with the right second-moment structure
      # X_dummy = zeros(T_dummy, k_coef)
      
      # Actually, construct minimal sufficient dummy obs via Cholesky factors.
      # We want: 
      #   X_dummy'X_dummy = T_dummy * diag(1, G)  -- block diag for const and lags
      #   X_dummy'Y_dummy = T_dummy * [0; rhs]    -- cross moments
      #   Y_dummy'Y_dummy = T_dummy * Gamma_0     -- dependent variable second moments
      #
      # The approach: create dummy observations so that the dummy-data OLS gives A_prior.
      # This is achieved by:
      #   Y_dummy_i = A_prior' * x_dummy_i + u_i, where u_i ~ N(0, Sigma_u_prior)
      # For the sufficient-statistic approach:
      # We can directly construct the dummy cross-products without creating T_dummy rows.

      # Construct the dummy sufficient statistics directly (more numerically stable):
      # Dummy cross-products (scaled by lambda)
      G_full <- matrix(0, nrow = k_coef, ncol = k_coef)
      G_full[1L, 1L] <- T_dummy  # intercept component
      G_full[-1L, -1L] <- T_dummy * G  # lag coefficient component

      # Add a small ridge for numerical stability
      G_full <- G_full + diag(k_coef) * 1e-10

      XY_dummy <- matrix(0, nrow = k_coef, ncol = n_obs)
      XY_dummy[-1L, ] <- T_dummy * rhs
      # Note: the intercept rows of XY_dummy stay zero because prior has zero mean for intercept

      YY_dummy <- T_dummy * Gamma_0

      # Augmented sufficient statistics
      XtX_aug <- XtX_data + G_full
      XtY_aug <- XtY_data + XY_dummy
      YtY_aug <- YtY_data + YY_dummy

      # Solve for augmented VAR coefficients
      XtX_aug_inv <- tryCatch(solve(XtX_aug), error = function(e) NULL)
      if (is.null(XtX_aug_inv)) {
        log_ml[li] <- -Inf
        next
      }
      A_aug <- XtX_aug_inv %*% XtY_aug
      A_lambda_list[[li]] <- A_aug

      # Residual covariance from augmented data
      resid_aug_cov <- (YtY_aug - t(XtY_aug) %*% A_aug) / (T_eff + T_dummy)
      Sigma_u_lambda_list[[li]] <- resid_aug_cov

      # ---- Log marginal likelihood approximation ----
      # Use the Laplace / BIC-type approximation:
      # ln p(Y|lambda) ~= -(T*M/2)*ln(2pi) - (T/2)*ln|Sigma_u| - (1/2)*tr(Sigma_u^{-1} * SS_res)
      # where SS_res = Y'Y - Y'X (X'X)^{-1} X'Y
      #
      # For BIC approximation:
      # log ML ~= -(T/2) * (M*ln(2pi) + ln|Sigma_u| + M)
      # where M = n_obs, and Sigma_u is the ML residual covariance

      # Actually compute the proper marginal likelihood for the VAR:
      # Using the normal-diffuse prior, the log marginal likelihood is:
      # ln p(Y|lambda) = const + (M/2)*ln|lambdaT*G_xx| - (M/2)*ln|X'X + lambdaT*G_xx| 
      #             + (nu/2)*ln|lambdaT*Sigma_u| - (nu'/2)*ln|S|
      # where nu = lambdaT - k + M + 1, nu' = T + lambdaT - k + M + 1
      # and S = Y'Y + lambdaT*Gamma_0 - (Y'X + lambdaT*Gamma_xy) * (X'X + lambdaT*Gamma_xx)^{-1} * (X'Y + lambdaT*Gamma_xy')

      T_total <- T_eff + T_dummy
      k_coef_per_eq <- k_coef  # coef per equation including intercept
      nu_prior <- T_dummy - k_coef_per_eq + n_obs + 1L
      nu_post  <- T_total - k_coef_per_eq + n_obs + 1L

      ln_det_XtX_prior <- {
        log(det(G_full))
      }
      
      ln_det_XtX_post <- {
        log(det(XtX_aug))
      }

      ln_det_Sp <- {
        log(det(Sigma_u_prior * T_dummy))
      }

      # S = Y'Y + lambdaT*Gamma_0 - (Y'X + lambdaT*Gamma_xy) * (X'X + lambdaT*Gamma_xx)^{-1} * (X'Y + lambdaT*Gamma_xy')
      S_post <- YtY_aug - t(XtY_aug) %*% A_aug
      ln_det_S_post <- {
        log(det(S_post))
      }

      # Log marginal likelihood -- EXACT conjugate matrix-variate
      # Normal-Inverse-Wishart evidence (Del Negro & Schorfheide 2004;
      # equivalently the natural-conjugate BVAR marginal likelihood of
      # Kadiyala & Karlsson 1997 / Banbura, Giannone & Reichlin 2010):
      #
      #   ln p(Y|lambda) = -(n*T/2) ln(pi) + [ln Gamma_n(nu_post/2) - ln Gamma_n(nu_prior/2)]
      #     + (n/2) ln|XtX_prior| - (n/2) ln|XtX_post|
      #     + (nu_prior/2) ln|S_prior| - (nu_post/2) ln|S_post|
      #
      # where Gamma_n(.) is the multivariate gamma function. The
      # log-multivariate-gamma RATIO term is NOT a constant across lambda
      # (nu_prior depends on T_dummy = round(lambda*n_obs)) and was
      # previously omitted here (a BIC/Laplace-style approximation) --
      # ground-truth validated against an independent closed-form oracle
      # in test-mdd-calibration.R (mdd_calibration(case = "dsge_var")).
      mvgamma_ratio <- .d15_lmvgamma_ratio(nu_post, nu_prior, n_obs)

      log_ml_val <- mvgamma_ratio +
                    0.5 * n_obs * ln_det_XtX_prior +
                    0.5 * nu_prior * ln_det_Sp -
                    0.5 * n_obs * ln_det_XtX_post -
                    0.5 * nu_post * ln_det_S_post

      # Add constant term: - (T*n_obs/2) * log(pi)
      log_ml_val <- log_ml_val - 0.5 * T_eff * n_obs * log(pi)

      log_ml[li] <- log_ml_val
    }

    # ---- 6. Find optimal lambda ----
    best_idx <- which.max(log_ml)
    lambda_opt <- lambda_grid[best_idx]
    log_ml_opt <- log_ml[best_idx]

    # Check if lambda_opt is at the boundary
    at_lower <- best_idx == 1L
    at_upper <- best_idx == n_lambda

    # --- ENHANCED: Extend grid when lambda* at boundary ---
    extended_grid <- FALSE
    if (at_lower && n_lambda >= 3L) {
      ext_low <- c(0.01, 0.03, 0.05, 0.10, 0.15, 0.20)
      ext_low <- ext_low[ext_low < lambda_grid[1L]]
      if (length(ext_low) > 0) {
        message("[D15] extending lambda grid lower...")
        ext_log_ml <- .extend_lambda_grid(
          XtX_data, XtY_data, YtY_data, G, rhs, Gamma_0,
          A_prior, Sigma_u_prior, T_eff, n_obs, k_coef,
          ext_low, obs_in_data
        )
        all_lambda <- c(ext_low, lambda_grid)
        all_log_ml <- c(ext_log_ml$log_ml, log_ml)
        all_A <- c(ext_log_ml$A_list, A_lambda_list)
        all_Sigma <- c(ext_log_ml$Sigma_list, Sigma_u_lambda_list)
        new_best_idx <- which.max(all_log_ml)
        lambda_opt <- all_lambda[new_best_idx]
        log_ml_opt <- all_log_ml[new_best_idx]
        lambda_grid <- all_lambda
        log_ml <- all_log_ml
        A_lambda_list <- all_A
        Sigma_u_lambda_list <- all_Sigma
        n_lambda <- length(lambda_grid)
        at_lower <- new_best_idx == 1L
        extended_grid <- TRUE
      }
    }
    if (at_upper) {
      ext_high <- c(20, 50, 100)
      ext_high <- ext_high[ext_high > lambda_grid[n_lambda]]
      if (length(ext_high) > 0) {
        message("[D15] extending lambda grid upper...")
        ext_log_ml <- .extend_lambda_grid(
          XtX_data, XtY_data, YtY_data, G, rhs, Gamma_0,
          A_prior, Sigma_u_prior, T_eff, n_obs, k_coef,
          ext_high, obs_in_data
        )
        all_lambda <- c(lambda_grid, ext_high)
        all_log_ml <- c(log_ml, ext_log_ml$log_ml)
        all_A <- c(A_lambda_list, ext_log_ml$A_list)
        all_Sigma <- c(Sigma_u_lambda_list, ext_log_ml$Sigma_list)
        new_best_idx <- which.max(all_log_ml)
        lambda_opt <- all_lambda[new_best_idx]
        log_ml_opt <- all_log_ml[new_best_idx]
        lambda_grid <- all_lambda
        log_ml <- all_log_ml
        A_lambda_list <- all_A
        Sigma_u_lambda_list <- all_Sigma
        n_lambda <- length(lambda_grid)
        at_upper <- new_best_idx == n_lambda
        extended_grid <- TRUE
      }
    }

    # Residual variance ratio per observable
    best_idx_final <- which(lambda_grid == lambda_opt)
    best_idx_final <- if (length(best_idx_final) == 0) which.min(abs(lambda_grid - lambda_opt)) else best_idx_final[1L]
    A_opt_use <- A_lambda_list[[best_idx_final]]
    Sigma_opt_use <- Sigma_u_lambda_list[[best_idx_final]]
    resid_var_prior <- diag(Sigma_u_prior)
    resid_var_opt <- diag(Sigma_opt_use)
    # Unrestricted = smallest lambda entry (closest to lambda->0)
    unrestricted_idx <- which.min(lambda_grid)
    Sigma_u_unrest <- Sigma_u_lambda_list[[unrestricted_idx]]
    resid_var_unrest <- diag(Sigma_u_unrest)
    resid_ratio_opt_unrest <- resid_var_opt / pmax(resid_var_unrest, 1e-16)
    resid_ratio_prior_unrest <- resid_var_prior / pmax(resid_var_unrest, 1e-16)
    names(resid_ratio_opt_unrest) <- obs_in_data
    names(resid_ratio_prior_unrest) <- obs_in_data
    worst_resid <- names(sort(resid_ratio_opt_unrest, decreasing = TRUE))[1:min(3, n_obs)]

    # Interpret lambda
    interpretation <- if (lambda_opt >= 5) {
      "DSGE restrictions strongly supported (lambda* large)."
    } else if (lambda_opt >= 2) {
      "DSGE restrictions moderately supported (lambda* >= 2)."
    } else if (lambda_opt >= 0.75) {
      "DSGE restrictions weakly supported (lambda* moderate)."
    } else {
      "DSGE restrictions receive little support (lambda* small). Unrestricted VAR preferred."
    }

    if (extended_grid) {
      interpretation <- paste(interpretation,
        sprintf(" Grid extended to [%.2f, %.0f].", min(lambda_grid), max(lambda_grid)))
    }
    if (at_lower) {
      interpretation <- paste(interpretation, "Optimal lambda at lower bound -- DSGE restrictions strongly rejected.")
    }
    if (at_upper) {
      interpretation <- paste(interpretation, "Optimal lambda at upper bound -- DSGE restrictions strongly supported.")
    }

    pass <- if (lambda_opt >= 2) TRUE else if (lambda_opt <= 0.5) FALSE else NA

    # ---- 7. Build results ----
    lambda_results <- data.frame(
      lambda     = lambda_grid,
      log_ml     = log_ml,
      log_ml_rel = log_ml - max(log_ml, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    lambda_results <- lambda_results[is.finite(lambda_results$log_ml), , drop = FALSE]

    result <- list(
      lambda_opt      = lambda_opt,
      log_ml_opt      = log_ml_opt,
      lambda_grid     = lambda_grid,
      log_ml_values   = log_ml,
      lambda_results  = lambda_results,
      A_prior         = A_prior,
      Sigma_u_prior   = Sigma_u_prior,
      A_optimal       = A_opt_use,
      Sigma_u_optimal = Sigma_opt_use,
      n_obs           = n_obs,
      var_lag         = var_lag,
      T_eff           = T_eff,
      k_coef          = k_coef,
      at_lower_bound  = at_lower,
      at_upper_bound  = at_upper,
      extended_grid   = extended_grid,
      resid_ratio_opt_unrest  = resid_ratio_opt_unrest,
      resid_ratio_prior_unrest = resid_ratio_prior_unrest,
      worst_resid_obs = worst_resid,
      # ---- DSGE-implied prior moments + data sufficient statistics.
      # Exposed (beyond what the diagnostic itself needs) so an independent
      # ground-truth oracle can reconstruct the exact lambda-scaled MNIW
      # sufficient statistics (XtX_prior(lambda) = blockdiag(T_dummy,
      # T_dummy*G), XtY_prior(lambda) = T_dummy*[0; rhs], YtY_prior(lambda)
      # = T_dummy*Gamma_0, with T_dummy = round(lambda*n_obs)) without
      # re-deriving the Yule-Walker/Lyapunov machinery -- see
      # mdd_calibration(case = "dsge_var") / test-mdd-calibration.R.
      G_prior         = G,
      rhs_prior       = rhs,
      Gamma0_prior    = Gamma_0,
      XtX_data        = XtX_data,
      XtY_data        = XtY_data,
      YtY_data        = YtY_data
    )

    # ---- 8. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .gg <- getNamespace("ggplot2")

      # (a) Log marginal likelihood vs lambda
      lr_df <- lambda_results
      if (nrow(lr_df) >= 3L) {
        p_ml <- .gg$ggplot(lr_df, .gg$aes(x = lambda, y = log_ml)) +
          .gg$geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.8) +
          .gg$geom_point(colour = dynhr_colours$dark_blue, size = 2) +
          .gg$geom_vline(xintercept = lambda_opt, linetype = "dashed",
                         colour = dynhr_colours$red, linewidth = 0.5) +
          .gg$annotate("text", x = lambda_opt, y = min(log_ml, na.rm = TRUE),
                       label = sprintf("lambda* = %.2f", lambda_opt),
                       hjust = -0.1, vjust = -0.5, size = 3.5,
                       colour = dynhr_colours$red) +
          .gg$scale_x_log10() +
          theme_dynhr_diagnostic() +
          .gg$labs(title = "D15: DSGE-VAR log marginal likelihood",
                   subtitle = interpretation,
                   x = expression(lambda ~ (tightness)),
                   y = "Log marginal likelihood")
        plots$log_ml <- .apply_meta(p_ml, meta)
      }

      # (b) Coefficient comparison: prior, optimal, unrestricted
      best_idx_plot <- if (exists("best_idx_final") && length(best_idx_final) > 0 && best_idx_final <= length(A_lambda_list))
                         best_idx_final else min(best_idx, length(A_lambda_list))
      A_unrest_slopes <- A_lambda_list[[which.min(lambda_grid)]][-1L, , drop = FALSE]
      A_opt_slopes    <- A_lambda_list[[best_idx_plot]][-1L, , drop = FALSE]
      A_prior_slopes  <- A_prior

      coef_compare <- data.frame(
        Coef = rep(paste0("lag", rep(seq_len(var_lag), each = n_obs), "_",
                          rep(obs_in_data, var_lag)), n_obs),
        Eq   = rep(obs_in_data, each = n_obs * var_lag),
        Prior = as.vector(A_prior_slopes),
        Optimal = as.vector(A_opt_slopes),
        Unrestricted = as.vector(A_unrest_slopes),
        stringsAsFactors = FALSE
      )

      coef_long <- stats::reshape(
        coef_compare,
        direction = "long",
        varying = c("Prior", "Optimal", "Unrestricted"),
        v.names = "Value",
        timevar = "Type",
        times = c("Prior", "Optimal", "Unrestricted"),
        idvar = c("Coef", "Eq")
      )
      rownames(coef_long) <- NULL

      p_coef <- .gg$ggplot(coef_long, .gg$aes(x = Coef, y = Value, fill = Type)) +
        .gg$geom_col(position = "dodge", width = 0.7) +
        .gg$facet_wrap(~ Eq, scales = "free_y", ncol = min(n_obs, 4)) +
        .gg$scale_fill_manual(values = c("Prior" = dynhr_colours$grey,
                                         "Optimal" = dynhr_colours$mid_blue,
                                         "Unrestricted" = dynhr_colours$orange)) +
        theme_dynhr_diagnostic() +
        .gg$theme(axis.text.x = .gg$element_text(angle = 45, hjust = 1, size = 7)) +
        .gg$labs(title = sprintf("D15: VAR coefficients (lambda* = %.2f)", lambda_opt),
                 x = NULL, y = "Coefficient")
      plots$coef_comparison <- .apply_meta(p_coef, meta)
    }

    # ---- 9. Summary ----
    # Residual variance interpretation
    resid_detail <- if (exists("resid_ratio_opt_unrest") && n_obs > 0) {
      high <- names(sort(resid_ratio_opt_unrest, decreasing = TRUE))[1:min(3, n_obs)]
      low  <- names(sort(resid_ratio_opt_unrest))[1:min(2, n_obs)]
      sprintf(" Resid var ratio (opt/OLS): %s=%.1f (high), %s=%.1f (low).",
              high[1], resid_ratio_opt_unrest[high[1]],
              low[1], resid_ratio_opt_unrest[low[1]])
    } else ""

    summary_str <- sprintf(
      "D15 DSGE-VAR: optimal lambda* = %.2f (log ML = %.1f).%s %s",
      lambda_opt, log_ml_opt,
      resid_detail,
      interpretation
    )

    llm_summary <- paste(c(
      sprintf("D15 | DSGE-VAR | %s",
              if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "FAIL" else "INFO"),
      sprintf("  lambda_opt=%.4f log_ml_opt=%.1f n_lambda_grid=%d",
              lambda_opt, log_ml_opt, nrow(lambda_results)),
      sprintf("  lambda_range=[%.2f, %.2f] extended=%s",
              min(lambda_grid), max(lambda_grid), if (extended_grid) "yes" else "no"),
      sprintf("  resid_var_ratio(opt/unrest): %s",
              paste(sprintf("%s=%.1f", names(resid_ratio_opt_unrest),
                            resid_ratio_opt_unrest), collapse = ", ")),
      sprintf("  resid_var_ratio(prior/unrest): %s",
              paste(sprintf("%s=%.1f", names(resid_ratio_prior_unrest),
                            resid_ratio_prior_unrest), collapse = ", ")),
      sprintf("  %s", interpretation),
      sprintf("  action: %s",
              if (lambda_opt >= 2)
                "DSGE restrictions well-supported. Consider using DSGE prior in VAR analysis."
              else if (lambda_opt >= 0.75)
                paste0("DSGE restrictions moderate support (lambda*=", round(lambda_opt, 2),
                       "). Review model specification for: ",
                       paste(head(worst_resid, 2), collapse = ", "), ".")
              else
                paste0("DSGE restrictions not supported (lambda*=", round(lambda_opt, 2),
                       "). Largest residual gap in: ",
                       paste(head(worst_resid, 2), collapse = ", "),
                       ". Consider respecifying dynamics for these series."))
    ), collapse = "\n")

    .make_result(
      result      = result,
      pass        = pass,
      plots       = plots,
      summary     = summary_str,
      llm_summary = llm_summary
    )
}


#' Helper: Evaluate lambda grid for D15 grid extension
#' @param XtX_data, XtY_data, YtY_data Pre-computed data sufficient statistics
#' @param G, rhs, Gamma_0 DSGE prior moments
#' @param A_prior, Sigma_u_prior DSGE prior VAR coefficients and residual cov
#' @param T_eff, n_obs, k_coef Dimensions
#' @param lambda_values Numeric vector of lambda values to evaluate
#' @return List with log_ml, A_list, Sigma_list
#' @noRd
.extend_lambda_grid <- function(XtX_data, XtY_data, YtY_data,
                                 G, rhs, Gamma_0,
                                 A_prior, Sigma_u_prior,
                                 T_eff, n_obs, k_coef,
                                 lambda_values, obs_in_data) {
  chol_Sigma_prior <- chol(Sigma_u_prior)
  if (is.null(chol_Sigma_prior)) {
    return(list(log_ml = rep(-Inf, length(lambda_values)),
                A_list = replicate(length(lambda_values), matrix(0, k_coef, n_obs), simplify = FALSE),
                Sigma_list = replicate(length(lambda_values), diag(n_obs), simplify = FALSE)))
  }

  n_lambda <- length(lambda_values)
  log_ml <- rep(NA_real_, n_lambda)
  A_list <- vector("list", n_lambda)
  Sigma_list <- vector("list", n_lambda)

  for (li in seq_len(n_lambda)) {
    lam <- lambda_values[li]
    T_dummy <- max(1L, round(lam * n_obs))

    G_full <- matrix(0, nrow = k_coef, ncol = k_coef)
    G_full[1L, 1L] <- T_dummy
    G_full[-1L, -1L] <- T_dummy * G
    G_full <- G_full + diag(k_coef) * 1e-10

    XY_dummy <- matrix(0, nrow = k_coef, ncol = n_obs)
    XY_dummy[-1L, ] <- T_dummy * rhs
    YY_dummy <- T_dummy * Gamma_0

    XtX_aug <- XtX_data + G_full
    XtY_aug <- XtY_data + XY_dummy
    YtY_aug <- YtY_data + YY_dummy

    XtX_aug_inv <- tryCatch(solve(XtX_aug), error = function(e) NULL)
    if (is.null(XtX_aug_inv)) {
      log_ml[li] <- -Inf
      A_list[[li]] <- matrix(0, k_coef, n_obs)
      Sigma_list[[li]] <- diag(n_obs)
      next
    }
    A_aug <- XtX_aug_inv %*% XtY_aug
    A_list[[li]] <- A_aug

    resid_aug_cov <- (YtY_aug - t(XtY_aug) %*% A_aug) / (T_eff + T_dummy)
    Sigma_list[[li]] <- resid_aug_cov

    T_total <- T_eff + T_dummy
    k_coef_per_eq <- k_coef
    nu_prior <- T_dummy - k_coef_per_eq + n_obs + 1L
    nu_post  <- T_total - k_coef_per_eq + n_obs + 1L

    ln_det_XtX_prior <- log(det(G_full))
    ln_det_XtX_post <- log(det(XtX_aug))
    ln_det_Sp <- log(det(Sigma_u_prior * T_dummy))

    S_post <- YtY_aug - t(XtY_aug) %*% A_aug
    ln_det_S_post <- log(det(S_post))

    mvgamma_ratio <- .d15_lmvgamma_ratio(nu_post, nu_prior, n_obs)

    log_ml_val <- mvgamma_ratio +
                  0.5 * n_obs * ln_det_XtX_prior +
                  0.5 * nu_prior * ln_det_Sp -
                  0.5 * n_obs * ln_det_XtX_post -
                  0.5 * nu_post * ln_det_S_post
    log_ml_val <- log_ml_val - 0.5 * T_eff * n_obs * log(pi)

    log_ml[li] <- log_ml_val
  }

  list(log_ml = log_ml, A_list = A_list, Sigma_list = Sigma_list)
}


#' Helper: log ratio of multivariate gamma functions Gamma_n(v1/2)/Gamma_n(v0/2)
#'
#' The multivariate gamma function is
#' \code{Gamma_n(a) = pi^{n(n-1)/4} * prod_{i=1}^n Gamma(a - (i-1)/2)}.
#' The \code{pi^{n(n-1)/4}} prefactor is identical for \code{Gamma_n(v1/2)}
#' and \code{Gamma_n(v0/2)} (same \code{n}) and cancels in the ratio, so only
#' the sum of \code{lgamma} differences is returned.
#'
#' @param v1,v0 Degrees of freedom (posterior, prior).
#' @param n     Matrix dimension (number of observables).
#' @return Numeric scalar: \code{log(Gamma_n(v1/2) / Gamma_n(v0/2))}.
#' @noRd
.d15_lmvgamma_ratio <- function(v1, v0, n) {
  i <- seq_len(n)
  sum(lgamma((v1 - i + 1) / 2) - lgamma((v0 - i + 1) / 2))
}
