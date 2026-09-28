## R/diag-post-d15-dsge-var.R
## --------------------------------------------------------------------------
## D15 DSGE-VAR(lambda) misspecification check, Del Negro & Schorfheide
## (2004, "DS04").
##
## The DSGE model's population moments Gamma(theta) define a conjugate
## Normal-Inverse-Wishart prior for a VAR(p) with intercept:
##
##   Sigma_u | theta        ~ IW( lambda*T * Sigma_u*(theta), lambda*T - k )
##   Phi | Sigma_u, theta   ~ N( Phi*(theta), Sigma_u (x) (lambda*T*Gamma_XX)^{-1} )
##
## with Phi* = Gamma_XX^{-1} Gamma_XY, Sigma_u* = Gamma_YY - Gamma_YX Phi*,
## T = effective sample (after the p initial lags), k = 1 + n*p, and
## x_t = (1, y_{t-1}', ..., y_{t-p}')'. The prior is proper only for
## lambda*T >= k + n, so lambda_min = (k + n) / T. The marginal likelihood
## (DS04 appendix; Del Negro, Schorfheide, Smets & Wouters 2007 eq. 19) is
##
##   p(Y|lambda) = |lambda T Gxx + X'X|^{-n/2} |S_post|^{-((1+lambda)T-k)/2}
##               / ( |lambda T Gxx|^{-n/2} |lambda T Sigma_u*|^{-(lambda T-k)/2} )
##               * pi^{-nT/2}
##               * prod_i Gamma(((1+lambda)T-k+1-i)/2) / prod_i Gamma((lambda T-k+1-i)/2)
##
## (the 2^{...} factors of DS04 combine with (2 pi)^{-nT/2} into pi^{-nT/2}),
## where S_post = (1+lambda) T Sigma_hat_u(lambda) is the posterior residual
## sum of squares. lambda -> lambda_min: data dominate; lambda large: the DSGE
## cross-equation restrictions dominate. lambda* = argmax p(Y|lambda).
##
## Model-implied moments (state-space s_t = T s_{t-1} + R e_t,
## y_t = Z s_{t-1} + D e_t, Var(e) = Sigma_e):
##   Gamma_h := E[y_{t+h} y_t'] = Z T^h Sigma_s Z' + Z T^{h-1} R Sigma_e D'  (h >= 1)
##   Gamma_0 = Z Sigma_s Z' + D Sigma_e D'
##   Gxx block (i,j) = E[y_{t-i} y_{t-j}'] = Gamma_{j-i} (j >= i), Gamma_{i-j}' (i > j)
##   Gxy block k     = E[y_{t-k} y_t']     = Gamma_k'
## The data are taken as deviations from the model's steady state (the prior
## puts the intercept at zero).
## --------------------------------------------------------------------------

#' D15. DSGE-VAR tightness (lambda)
#'
#' Estimates the DSGE-VAR(lambda) of Del Negro & Schorfheide (2004) and
#' reports the tightness lambda* that maximises the exact marginal
#' likelihood over a grid. Large lambda* supports the DSGE restrictions;
#' lambda* at the minimum admissible value \code{(k + n) / T} favours the
#' unrestricted VAR.
#'
#' \strong{Demeaning (\code{demean}).} DS04's derivation works with data in
#' deviations from the DSGE's own steady state: the dummy-observation prior
#' pins the VAR intercept at zero, so the VAR's implied unconditional mean IS
#' the model steady state. That is also what makes \code{lambda = Inf} a
#' well-defined limit. D15 therefore demeans by \code{dr}'s steady state by
#' default. \strong{This differs from Dynare}, whose \code{prefilter} option
#' demeans with the SAMPLE mean; both conventions are defensible (the sample
#' mean is agnostic about whether the model's steady state matches the data's
#' unconditional mean, which matters when the upstream detrending differs from
#' the model's), but they give different numbers, so a Dynare DSGE-VAR
#' workflow ported here should set \code{demean = "sample"} to compare.
#'
#' \strong{\code{lambda = Inf}.} At the infinite-tightness endpoint the
#' cross-equation restrictions bind exactly and the DSGE-VAR collapses to the
#' VAR(p) implied by the DSGE state space: coefficients
#' \eqn{\Phi^* = \Gamma_{XX}^{-1}\Gamma_{XY}} with a zero intercept and
#' residual covariance \eqn{\Sigma_u^*}. There is no prior left to integrate
#' out, so \code{log p(Y | Inf)} is that fixed Gaussian VAR's own likelihood.
#' (It is the DSGE's VAR(p) APPROXIMATION, not the DSGE itself -- a linear
#' DSGE state space is a VARMA in general.) Del Negro, Schorfheide, Smets &
#' Wouters (2007) use this as the "DSGE" reference endpoint of the lambda
#' sweep, against \code{lambda_min} as the unrestricted-VAR endpoint.
#'
#' @param dr          DecisionRules object from \code{solve_perturbation()}.
#'   Must contain \code{ghx}, \code{ghu}, \code{state_idx}, \code{endo_names},
#'   \code{exo_names}.
#' @param data        Matrix (T x n_obs) of observables, in LEVELS. They are
#'   demeaned here according to \code{demean} (see Details); pass
#'   \code{demean = "none"} if they are already in deviations.
#' @param obs_names   Observable names (columns of \code{data} and entries of
#'   \code{dr$endo_names}). Default \code{colnames(data)}.
#' @param sigma_e     Optional shock covariance. If NULL: from \code{model}
#'   (at \code{params}), else \code{dr$Sigma_e}; never an identity default.
#' @param model       dynhr_mod object (for the shock covariance).
#' @param params      Parameter vector at which \code{dr} was solved (default
#'   \code{model$param_values}).
#' @param var_lag     VAR lag order p (default 2).
#' @param lambda_grid Lambda values (DS04 units: dummy observations =
#'   lambda * T). Values below \code{lambda_min = (k + n) / T} are dropped.
#'   \code{Inf} is a valid, first-class grid point (the DSGE-implied VAR
#'   limit) and is in the default grid.
#' @param demean      How to centre \code{data} before fitting:
#'   \code{"steady_state"} (default) subtracts the model steady state of each
#'   observable, taken from \code{dr$ys}; \code{"sample"} subtracts the
#'   column means (Dynare's \code{prefilter} convention); \code{"none"}
#'   leaves the data as given. \code{"steady_state"} silently degrades to
#'   \code{"sample"} -- and says so in the summary -- when \code{dr} carries
#'   no usable steady state for the observables.
#' @param meta        Optional \code{dynhr_diag_meta} from \code{diag_meta()}.
#' @param verbose     Verbosity for progress messages (NULL = package option).
#' @return \code{dynhr_diagnostic} list with optimal lambda and comparison.
#' @references
#'   Del Negro, M. & Schorfheide, F. (2004). Priors from general equilibrium
#'   models for VARs. \emph{International Economic Review}, 45(2), 643-673.
#'   Del Negro, M., Schorfheide, F., Smets, F. & Wouters, R. (2007). On the fit
#'   of New Keynesian models. \emph{Journal of Business & Economic
#'   Statistics}, 25(2), 123-143.
#' @noRd
d15_dsge_var <- function(dr            = NULL,
                         data          = NULL,
                         obs_names     = NULL,
                         sigma_e       = NULL,
                         model         = NULL,
                         params        = NULL,
                         var_lag       = 2L,
                         lambda_grid   = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10,
                                           Inf),
                         demean        = c("steady_state", "sample", "none"),
                         meta          = NULL,
                         verbose       = NULL) {
  demean <- match.arg(demean)

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
  var_lag <- as.integer(var_lag)
  if (length(var_lag) != 1L || is.na(var_lag) || var_lag < 1L)
    .dynhr_abort("d15_dsge_var(): `var_lag` must be a positive integer.")
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  ## lambda = Inf is a first-class grid point, not an edge case: DNSS (2007)
  ## plot the log marginal likelihood over a lambda sweep with the
  ## DSGE-implied VAR limit as the reference endpoint (see .d15_eval_grid()).
  if (length(lambda_grid) == 0L || anyNA(lambda_grid) || any(lambda_grid <= 0))
    .dynhr_abort("d15_dsge_var(): `lambda_grid` must be positive values (Inf allowed).")

  # ---- 2. Prepare data ----
  if (!is.matrix(data)) data <- as.matrix(data)
  if (is.null(obs_names)) obs_names <- colnames(data)
  if (is.null(obs_names))
    return(.make_result(pass = NA,
      summary = "D15: data has no column names and obs_names was not given."))
  if (is.null(colnames(data))) {
    if (length(obs_names) != ncol(data))
      return(.make_result(pass = NA, summary = sprintf(
        "D15: data has no column names and %d columns, but %d obs_names.",
        ncol(data), length(obs_names))))
    colnames(data) <- obs_names
  }
  miss_data <- setdiff(obs_names, colnames(data))
  if (length(miss_data) > 0L)
    return(.make_result(pass = NA, summary = paste0(
      "D15: observables not found in data columns: ",
      paste(miss_data, collapse = ", "), ".")))
  obs_in_data <- obs_names
  Y <- data[, obs_in_data, drop = FALSE]
  storage.mode(Y) <- "double"
  if (any(!is.finite(Y)))
    return(.make_result(pass = NA,
      summary = "D15: data contain missing or non-finite values; DSGE-VAR needs a balanced panel."))

  ## ---- Demeaning (0.9.4) ------------------------------------------------
  ## DS04's DSGE-VAR derivation works with data in deviations from the DSGE's
  ## own steady state: the dummy-observation prior pins the VAR intercept at
  ## zero, so the implied unconditional mean IS the model steady state, which
  ## is what makes the lambda -> Inf endpoint (the DSGE-implied VAR) a
  ## well-defined limit. Until 0.9.4 the orchestrator handed this function the
  ## raw LEVEL series while the prior still centred at zero, so on any model
  ## with non-zero observable steady states (nk_demo: 0.5, 2, 4) the DSGE
  ## prior was compared against a VAR fitted to un-demeaned data.
  ##
  ## This differs from Dynare's `prefilter`, which demeans with the SAMPLE
  ## mean. Both are defensible -- the sample mean is agnostic about whether
  ## the model's steady state matches the data's unconditional mean, while the
  ## steady state is faithful to DS04's derivation -- but they give different
  ## numbers, so a Dynare DSGE-VAR workflow ported here will not reproduce
  ## exactly unless `demean = "sample"` is requested.
  demean_note <- ""
  ss_used <- NULL
  if (demean == "steady_state") {
    ss_vec <- .d15_obs_steady_state(dr, obs_in_data)
    if (is.null(ss_vec)) {
      demean <- "sample"
      demean_note <- paste0(
        " Data demeaned by the SAMPLE mean: the decision rule carries no ",
        "steady state for these observables (dr$ys), so the model steady ",
        "state was unavailable.")
    } else {
      ss_used <- ss_vec
      Y <- sweep(Y, 2L, ss_vec, "-")
      demean_note <- sprintf(
        " Data demeaned by the model steady state (%s).",
        paste(sprintf("%s=%.4g", obs_in_data, ss_vec), collapse = ", "))
    }
  }
  if (demean == "sample") {
    ss_used <- colMeans(Y)
    Y <- sweep(Y, 2L, ss_used, "-")
    if (!nzchar(demean_note))
      demean_note <- " Data demeaned by the SAMPLE mean (demean = \"sample\")."
  }

  n_obs  <- ncol(Y)
  T_obs  <- nrow(Y)
  k_coef <- 1L + n_obs * var_lag
  T_eff  <- T_obs - var_lag
  if (T_eff <= k_coef) {
    return(.make_result(pass = NA, summary = sprintf(
      "D15: Too few observations (%d) for a VAR(%d) in %d variables (need T - p > %d).",
      T_obs, var_lag, n_obs, k_coef)))
  }

  # ---- 3. DSGE prior moments ----
  obs_idx <- match(obs_in_data, dr$endo_names)
  if (anyNA(obs_idx)) {
    return(.make_result(pass = NA, summary = paste0(
      "D15: observables not found in dr$endo_names: ",
      paste(obs_in_data[is.na(obs_idx)], collapse = ", "), ".")))
  }
  state_idx <- dr$state_idx
  T_mat <- dr$ghx[state_idx, , drop = FALSE]
  R_mat <- dr$ghu[state_idx, , drop = FALSE]
  Z_mat <- dr$ghx[obs_idx, , drop = FALSE]
  D_mat <- dr$ghu[obs_idx, , drop = FALSE]
  n_exo <- ncol(R_mat)

  Sigma_e <- if (!is.null(sigma_e)) {
    as.matrix(sigma_e)
  } else if (!is.null(model)) {
    p_use <- model$param_values
    if (!is.null(params)) p_use[names(params)] <- params   # partial vectors allowed
    .get_shock_cov(model, dr$exo_names, p_use)
  } else if (!is.null(dr$Sigma_e)) {
    as.matrix(dr$Sigma_e)
  } else {
    .dynhr_abort("d15_dsge_var(): no shock covariance -- pass `sigma_e` or `model` ",
                 "(dr$Sigma_e is also absent).")
  }
  if (!identical(dim(Sigma_e), c(n_exo, n_exo)))
    .dynhr_abort(sprintf("d15_dsge_var(): Sigma_e must be %d x %d (got %s).",
                         n_exo, n_exo, paste(dim(Sigma_e), collapse = " x ")))

  mom <- .d15_dsge_moments(T_mat, R_mat, Z_mat, D_mat, Sigma_e, var_lag)
  G       <- mom$Gxx
  rhs     <- mom$Gxy
  Gamma_0 <- mom$Gyy

  G_rcond <- rcond(G)
  if (!is.finite(G_rcond) || G_rcond < 1e-14)
    return(.make_result(pass = NA, summary = sprintf(paste0(
      "D15: model-implied regressor covariance is near-singular (rcond=%.2e). ",
      "With %d observables and %d shocks the DSGE autocovariances are rank-deficient; ",
      "add measurement error or reduce observables."), G_rcond, n_obs, n_exo)))
  A_prior <- solve(G, rhs)
  Sigma_u_prior <- Gamma_0 - t(rhs) %*% A_prior
  Sigma_u_prior <- (Sigma_u_prior + t(Sigma_u_prior)) / 2
  dimnames(Sigma_u_prior) <- list(obs_in_data, obs_in_data)
  ev_prior <- eigen(Sigma_u_prior, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev_prior) <= 1e-12 * max(abs(ev_prior), 1e-300)) {
    return(.make_result(pass = NA, summary = paste0(
      sprintf("D15: DSGE prior residual covariance is not positive definite (min eigenvalue %.3e). ",
              min(ev_prior)),
      "Smallest innovation variances: ",
      paste(obs_in_data[order(diag(Sigma_u_prior))][seq_len(min(3L, n_obs))], collapse = ", "),
      ". Action: add measurement error for these observables (the model has ",
      "fewer shocks than observables), or check Sigma_e.")))
  }

  # ---- 4. Data sufficient statistics ----
  Y_lagged <- matrix(NA_real_, nrow = T_eff, ncol = n_obs * var_lag)
  for (k in seq_len(var_lag)) {
    Y_lagged[, ((k - 1L) * n_obs + 1L):(k * n_obs)] <-
      Y[(var_lag - k + 1L):(T_obs - k), , drop = FALSE]
  }
  Y_dep  <- Y[(var_lag + 1L):T_obs, , drop = FALSE]
  X_data <- cbind(1, Y_lagged)
  XtX_data <- crossprod(X_data)
  XtY_data <- crossprod(X_data, Y_dep)
  YtY_data <- crossprod(Y_dep)

  # Unrestricted (lambda -> 0) benchmark: OLS on the data alone.
  A_ols     <- solve(XtX_data, XtY_data)
  Sigma_ols <- (YtY_data - t(XtY_data) %*% A_ols) / T_eff
  dimnames(Sigma_ols) <- list(obs_in_data, obs_in_data)

  # ---- 5. Lambda grid (DS04 units) ----
  lambda_min <- (k_coef + n_obs) / T_eff
  dropped    <- lambda_grid[lambda_grid < lambda_min]
  lambda_grid <- lambda_grid[lambda_grid >= lambda_min]
  if (length(lambda_grid) == 0L)
    lambda_grid <- lambda_min
  stats_args <- list(T_eff = T_eff, n_obs = n_obs, k_coef = k_coef,
                     G = G, rhs = rhs, Gamma_0 = Gamma_0,
                     XtX_data = XtX_data, XtY_data = XtY_data, YtY_data = YtY_data)
  ev <- .d15_eval_grid(lambda_grid, stats_args)

  extended_grid <- FALSE
  best_idx <- which.max(ev$log_ml)
  if (length(best_idx) == 1L && best_idx == 1L && lambda_grid[1L] > lambda_min * (1 + 1e-8)) {
    ext_low <- lambda_min * (lambda_grid[1L] / lambda_min)^c(0, 1/3, 2/3)
    .dynhr_inform(sprintf("[D15] lambda* at grid lower end; extending grid down to lambda_min = %.3f",
                          lambda_min), verbose = verbose)
    ev_low <- .d15_eval_grid(ext_low, stats_args)
    lambda_grid <- c(ext_low, lambda_grid)
    ev <- .d15_bind_eval(ev_low, ev)
    extended_grid <- TRUE
  }
  best_idx <- which.max(ev$log_ml)
  if (length(best_idx) == 1L && best_idx == length(lambda_grid)) {
    ext_high <- c(20, 50, 100)
    ext_high <- ext_high[ext_high > lambda_grid[length(lambda_grid)]]
    if (length(ext_high) > 0L) {
      .dynhr_inform("[D15] lambda* at grid upper end; extending grid to lambda = ",
                    paste(ext_high, collapse = ", "), verbose = verbose)
      ev_high <- .d15_eval_grid(ext_high, stats_args)
      lambda_grid <- c(lambda_grid, ext_high)
      ev <- .d15_bind_eval(ev, ev_high)
      extended_grid <- TRUE
    }
  }
  log_ml <- ev$log_ml
  n_lambda <- length(lambda_grid)

  # ---- 6. Optimal lambda ----
  best_idx <- which.max(log_ml)
  if (length(best_idx) == 0L)
    return(.make_result(pass = NA,
      summary = "D15: log marginal likelihood is non-finite at every lambda."))
  lambda_opt <- lambda_grid[best_idx]
  log_ml_opt <- log_ml[best_idx]
  at_lower   <- best_idx == 1L
  at_upper   <- best_idx == n_lambda
  at_min     <- at_lower && lambda_opt <= lambda_min * (1 + 1e-8)

  A_opt_use     <- ev$A_list[[best_idx]]
  Sigma_opt_use <- ev$Sigma_list[[best_idx]]
  dimnames(Sigma_opt_use) <- list(obs_in_data, obs_in_data)
  resid_var_prior  <- diag(Sigma_u_prior)
  resid_var_opt    <- diag(Sigma_opt_use)
  resid_var_unrest <- diag(Sigma_ols)
  resid_ratio_opt_unrest   <- resid_var_opt / pmax(resid_var_unrest, 1e-300)
  resid_ratio_prior_unrest <- resid_var_prior / pmax(resid_var_unrest, 1e-300)
  names(resid_ratio_opt_unrest)   <- obs_in_data
  names(resid_ratio_prior_unrest) <- obs_in_data
  worst_resid <- names(sort(resid_ratio_prior_unrest, decreasing = TRUE))[seq_len(min(3L, n_obs))]

  interpretation <- if (lambda_opt >= 5) {
    "DSGE restrictions strongly supported (lambda* large)."
  } else if (lambda_opt >= 2) {
    "DSGE restrictions moderately supported (lambda* >= 2)."
  } else if (lambda_opt >= 0.75) {
    "DSGE restrictions weakly supported (lambda* moderate)."
  } else {
    "DSGE restrictions receive little support (lambda* small); the unrestricted VAR fits better."
  }
  if (at_min) {
    interpretation <- paste(interpretation,
      "lambda* is the minimum admissible value -- DSGE restrictions strongly rejected.")
  } else if (at_upper && is.finite(lambda_opt)) {
    interpretation <- paste(interpretation,
      "lambda* is the largest grid value -- evidence may still rise beyond it.")
  } else if (at_upper) {
    interpretation <- paste(interpretation,
      paste0("lambda* = Inf: the DSGE-implied VAR (cross-equation restrictions ",
             "imposed exactly) has the highest marginal likelihood on the grid."))
  }

  pass <- if (lambda_opt >= 2) TRUE else if (lambda_opt <= 0.5) FALSE else NA

  # ---- 7. Results ----
  lambda_results <- data.frame(
    lambda     = lambda_grid,
    log_ml     = log_ml,
    log_ml_rel = log_ml - log_ml_opt,
    stringsAsFactors = FALSE
  )
  lambda_results <- lambda_results[is.finite(lambda_results$log_ml), , drop = FALSE]

  result <- list(
    lambda_opt      = lambda_opt,
    log_ml_opt      = log_ml_opt,
    lambda_min      = lambda_min,
    lambda_dropped  = dropped,
    lambda_grid     = lambda_grid,
    log_ml_values   = log_ml,
    lambda_results  = lambda_results,
    A_prior         = A_prior,
    Sigma_u_prior   = Sigma_u_prior,
    A_optimal       = A_opt_use,
    Sigma_u_optimal = Sigma_opt_use,
    A_ols           = A_ols,
    Sigma_u_ols     = Sigma_ols,
    Sigma_e         = Sigma_e,
    n_obs           = n_obs,
    var_lag         = var_lag,
    T_eff           = T_eff,
    k_coef          = k_coef,
    demean          = demean,
    demean_offset   = ss_used,
    at_lower_bound  = at_lower,
    at_upper_bound  = at_upper,
    at_lambda_min   = at_min,
    extended_grid   = extended_grid,
    resid_ratio_opt_unrest   = resid_ratio_opt_unrest,
    resid_ratio_prior_unrest = resid_ratio_prior_unrest,
    worst_resid_obs = worst_resid,
    # DSGE-implied population moments + data sufficient statistics, exposed so
    # an independent oracle can rebuild the lambda-scaled NIW prior:
    # XtX_prior = lambda*T_eff*blockdiag(1, G_prior),
    # XtY_prior = lambda*T_eff*[0; rhs_prior], YtY_prior = lambda*T_eff*Gamma0_prior.
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
    lr_df <- lambda_results
    ## lambda = Inf cannot sit on a log10 axis; it is drawn as a horizontal
    ## reference line (the DSGE-implied VAR endpoint) instead.
    lml_inf <- lr_df$log_ml_rel[!is.finite(lr_df$lambda)]
    lr_df   <- lr_df[is.finite(lr_df$lambda), , drop = FALSE]
    if (nrow(lr_df) >= 2L) {
      p_ml <- .gg$ggplot(lr_df, .gg$aes(x = lambda, y = log_ml_rel)) +
        .gg$geom_vline(xintercept = lambda_min, linetype = "dotted",
                       colour = dynhr_colours$grey, linewidth = 0.5) +
        .gg$geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.8) +
        .gg$geom_point(colour = dynhr_colours$dark_blue, size = 2) +
        .gg$geom_vline(xintercept = lambda_opt, linetype = "dashed",
                       colour = dynhr_colours$red, linewidth = 0.5) +
        .gg$annotate("text", x = lambda_opt, y = min(lr_df$log_ml_rel),
                     label = sprintf("lambda* = %s", .d15_fmt_lambda(lambda_opt)),
                     hjust = if (at_upper) 1.05 else -0.05, vjust = 0, size = 3.5,
                     colour = dynhr_colours$red) +
        { if (length(lml_inf) == 1L && is.finite(lml_inf))
            .gg$geom_hline(yintercept = lml_inf, linetype = "longdash",
                           colour = dynhr_colours$teal, linewidth = 0.5) } +
        { if (length(lml_inf) == 1L && is.finite(lml_inf))
            .gg$annotate("text", x = max(lr_df$lambda), y = lml_inf,
                         label = "lambda = Inf (DSGE-implied VAR)",
                         hjust = 1, vjust = -0.4, size = 3,
                         colour = dynhr_colours$teal) } +
        .gg$scale_x_log10() +
        theme_dynhr_diagnostic() +
        .gg$labs(title = "D15: DSGE-VAR marginal likelihood vs tightness",
                 subtitle = sprintf(paste0(
                   "VAR(%d), n = %d, T = %d, prior = lambda x T dummy obs; ",
                   "log ML at lambda* = %.1f\n",
                   "dashed: lambda* = %s; dotted: minimum admissible lambda = (k+n)/T = %.3f"),
                   var_lag, n_obs, T_eff, log_ml_opt,
                   .d15_fmt_lambda(lambda_opt), lambda_min),
                 x = expression(lambda ~ "(DSGE prior weight, log scale)"),
                 y = "log p(Y | lambda) - max (nats)")
      plots$log_ml <- .apply_meta(p_ml, meta)
    }

    # (b) Slope coefficients: DSGE prior, posterior mean at lambda*, OLS.
    reg_lab <- paste0(rep(obs_in_data, var_lag), "(-", rep(seq_len(var_lag), each = n_obs), ")")
    n_reg <- length(reg_lab)
    type_lab <- c("DSGE prior",
                  sprintf("DSGE-VAR (lambda* = %s)", .d15_fmt_lambda(lambda_opt)),
                  "OLS VAR")
    coef_long <- data.frame(
      Regressor = factor(rep(reg_lab, 3L * n_obs), levels = reg_lab),
      Eq        = factor(rep(rep(paste("equation:", obs_in_data), each = n_reg), 3L),
                         levels = paste("equation:", obs_in_data)),
      Type      = factor(rep(type_lab, each = n_reg * n_obs), levels = type_lab),
      Value     = c(as.vector(A_prior), as.vector(A_opt_use[-1L, , drop = FALSE]),
                    as.vector(A_ols[-1L, , drop = FALSE])),
      stringsAsFactors = FALSE
    )
    fill_vals <- stats::setNames(c(dynhr_colours$grey, dynhr_colours$mid_blue,
                                   dynhr_colours$orange), type_lab)
    p_coef <- .gg$ggplot(coef_long, .gg$aes(x = Regressor, y = Value, fill = Type)) +
      .gg$geom_hline(yintercept = 0, colour = dynhr_colours$grey, linewidth = 0.3) +
      .gg$geom_col(position = .gg$position_dodge(width = 0.8), width = 0.75) +
      .gg$facet_wrap(~ Eq, scales = "free_y", ncol = min(n_obs, 3L)) +
      .gg$scale_fill_manual(values = fill_vals, name = NULL) +
      theme_dynhr_diagnostic() +
      .gg$theme(axis.text.x = .gg$element_text(angle = 45, hjust = 1, size = 8),
                legend.position = "bottom") +
      .gg$labs(title = "D15: VAR slope coefficients",
               subtitle = "Each panel is one VAR equation; bars are the lagged regressors",
               x = "Regressor (lag)", y = "Coefficient")
    plots$coef_comparison <- .apply_meta(p_coef, meta)
  }

  # ---- 9. Summary ----
  hi <- names(which.max(resid_ratio_opt_unrest))
  lo <- names(which.min(resid_ratio_opt_unrest))
  resid_detail <- sprintf(
    " Innovation variance ratio (lambda*/OLS): %s=%.2f (max), %s=%.2f (min).",
    hi, resid_ratio_opt_unrest[[hi]], lo, resid_ratio_opt_unrest[[lo]])
  grid_note <- sprintf(" Grid [%.2f, %s] (lambda_min = %.2f%s%s).",
                       min(lambda_grid), .d15_fmt_lambda(max(lambda_grid)),
                       lambda_min,
                       if (extended_grid) "; extended" else "",
                       if (length(dropped)) sprintf("; %d value(s) below lambda_min dropped",
                                                    length(dropped)) else "")

  summary_str <- sprintf(
    "D15 DSGE-VAR: optimal lambda* = %s (log ML = %.1f).%s%s%s %s",
    .d15_fmt_lambda(lambda_opt), log_ml_opt, resid_detail, grid_note,
    demean_note, interpretation)

  llm_summary <- paste(c(
    sprintf("D15 | DSGE-VAR | %s",
            if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "FAIL" else "INFO"),
    sprintf("  lambda_opt=%s log_ml_opt=%.3f lambda_min=%.4f n_lambda_grid=%d",
            .d15_fmt_lambda(lambda_opt), log_ml_opt, lambda_min,
            nrow(lambda_results)),
    sprintf("  lambda_range=[%.3f, %s] extended=%s",
            min(lambda_grid), .d15_fmt_lambda(max(lambda_grid)),
            if (extended_grid) "yes" else "no"),
    sprintf("  demean=%s%s", demean, demean_note),
    sprintf("  resid_var_ratio(lambda*/OLS): %s",
            paste(sprintf("%s=%.3f", names(resid_ratio_opt_unrest),
                          resid_ratio_opt_unrest), collapse = ", ")),
    sprintf("  resid_var_ratio(DSGE prior/OLS): %s",
            paste(sprintf("%s=%.3f", names(resid_ratio_prior_unrest),
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
                     "). Largest DSGE-vs-OLS innovation-variance gap in: ",
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


#' Helper: format a lambda value that may be Inf
#' @noRd
.d15_fmt_lambda <- function(x) {
  if (is.finite(x)) sprintf("%.2f", x) else "Inf"
}


#' Helper: model steady state of the D15 observables
#'
#' \code{dr$ys} is either a named numeric vector or a \code{dynhr_steady}
#' object with \code{$values}. Returns NULL (not an error, and not a silent
#' zero) when it is absent, unnamed, non-finite or missing an observable --
#' the caller then falls back to the sample mean and SAYS so.
#' @noRd
.d15_obs_steady_state <- function(dr, obs_names) {
  ys <- dr$ys
  if (inherits(ys, "dynhr_steady")) ys <- ys$values
  if (is.null(ys) || is.null(names(ys))) return(NULL)
  if (!all(obs_names %in% names(ys))) return(NULL)
  out <- as.numeric(ys[obs_names])
  if (!all(is.finite(out))) return(NULL)
  names(out) <- obs_names
  out
}


#' Helper: DSGE-implied VAR population moments for D15
#'
#' For s_t = T s_{t-1} + R e_t, y_t = Z s_{t-1} + D e_t, Var(e) = Sigma_e,
#' returns Gyy = E[y_t y_t'], Gxx = E[x x'] and Gxy = E[x y_t'] with
#' x = (y_{t-1}', ..., y_{t-p}')' (no intercept).
#' @noRd
.d15_dsge_moments <- function(T_mat, R_mat, Z_mat, D_mat, Sigma_e, var_lag) {
  n <- nrow(Z_mat)
  RS <- R_mat %*% Sigma_e
  Sigma_s <- solve_lyapunov(T_mat, RS %*% t(R_mat))
  Gamma_0 <- Z_mat %*% Sigma_s %*% t(Z_mat) + D_mat %*% Sigma_e %*% t(D_mat)
  Gamma_0 <- (Gamma_0 + t(Gamma_0)) / 2
  # Gam[[h + 1]] = E[y_{t+h} y_t']
  Gam <- vector("list", var_lag + 1L)
  Gam[[1L]] <- Gamma_0
  T_pow_prev <- diag(nrow(T_mat))              # T^{h-1}
  for (h in seq_len(var_lag)) {
    T_pow <- T_pow_prev %*% T_mat              # T^h
    Gam[[h + 1L]] <- Z_mat %*% T_pow %*% Sigma_s %*% t(Z_mat) +
      Z_mat %*% T_pow_prev %*% RS %*% t(D_mat)
    T_pow_prev <- T_pow
  }
  Gxx <- matrix(0, var_lag * n, var_lag * n)
  Gxy <- matrix(0, var_lag * n, n)
  blk <- function(i) ((i - 1L) * n + 1L):(i * n)
  for (i in seq_len(var_lag)) {
    for (j in seq_len(var_lag)) {
      Gxx[blk(i), blk(j)] <- if (j >= i) Gam[[j - i + 1L]] else t(Gam[[i - j + 1L]])
    }
    Gxy[blk(i), ] <- t(Gam[[i + 1L]])
  }
  Gxx <- (Gxx + t(Gxx)) / 2
  list(Gyy = Gamma_0, Gxx = Gxx, Gxy = Gxy)
}


#' Helper: log|M| for a symmetric matrix, NA if not positive definite
#' @noRd
.d15_logdet_pd <- function(M) {
  d <- determinant(M, logarithm = TRUE)
  if (d$sign <= 0 || !is.finite(d$modulus)) return(NA_real_)
  as.numeric(d$modulus)
}


#' Helper: exact DS04 log marginal likelihood at each lambda
#'
#' @param lambda_values Lambda values (must satisfy lambda*T_eff >= k + n).
#' @param s List with T_eff, n_obs, k_coef, G, rhs, Gamma_0, XtX_data,
#'   XtY_data, YtY_data.
#' @return List with log_ml, A_list (posterior mean coefficients, k x n) and
#'   Sigma_list (posterior residual covariance S_post / ((1+lambda) T)).
#' @noRd
.d15_eval_grid <- function(lambda_values, s) {
  n <- s$n_obs; k <- s$k_coef; Te <- s$T_eff
  nl <- length(lambda_values)
  log_ml <- rep(NA_real_, nl)
  A_list <- vector("list", nl)
  Sigma_list <- vector("list", nl)
  A_star     <- solve(s$G, s$rhs)
  Sigma_star <- s$Gamma_0 - t(s$rhs) %*% A_star
  Sigma_star <- (Sigma_star + t(Sigma_star)) / 2
  for (li in seq_len(nl)) {
    ## lambda = Inf: the DSGE cross-equation restrictions bind exactly and the
    ## DSGE-VAR collapses to the VAR(p) implied by the DSGE state space --
    ## coefficients Phi* = Gxx^{-1} Gxy with a ZERO intercept (the data are in
    ## deviations from the steady state) and residual covariance Sigma_u*.
    ## The marginal likelihood is then just that fixed Gaussian VAR's
    ## likelihood, with no prior left to integrate out. DNSS (2007) treat this
    ## as the "DSGE" endpoint of the lambda sweep.
    if (!is.finite(lambda_values[li])) {
      A_inf <- rbind(matrix(0, 1L, n), A_star)
      S_inf <- s$YtY_data - t(s$XtY_data) %*% A_inf -
        t(A_inf) %*% s$XtY_data + t(A_inf) %*% s$XtX_data %*% A_inf
      S_inf <- (S_inf + t(S_inf)) / 2
      ld_Su <- .d15_logdet_pd(Sigma_star)
      A_list[[li]] <- A_inf
      Sigma_list[[li]] <- Sigma_star
      if (!is.na(ld_Su))
        log_ml[li] <- -0.5 * Te * n * log(2 * pi) - 0.5 * Te * ld_Su -
          0.5 * sum(diag(solve(Sigma_star, S_inf)))
      next
    }
    Ts <- lambda_values[li] * Te
    XtX0 <- matrix(0, k, k)
    XtX0[1L, 1L]   <- Ts
    XtX0[-1L, -1L] <- Ts * s$G
    XtY0 <- matrix(0, k, n)
    XtY0[-1L, ] <- Ts * s$rhs
    XtX1 <- XtX0 + s$XtX_data
    XtY1 <- XtY0 + s$XtY_data
    A1   <- solve(XtX1, XtY1)
    S1   <- Ts * s$Gamma_0 + s$YtY_data - t(XtY1) %*% A1
    S1   <- (S1 + t(S1)) / 2
    A_list[[li]] <- A1
    Sigma_list[[li]] <- S1 / (Ts + Te)
    nu0 <- Ts - k
    nu1 <- Ts + Te - k
    ld_X0 <- .d15_logdet_pd(XtX0)
    ld_X1 <- .d15_logdet_pd(XtX1)
    ld_S0 <- .d15_logdet_pd(Ts * Sigma_star)
    ld_S1 <- .d15_logdet_pd(S1)
    if (nu0 <= n - 1 || anyNA(c(ld_X0, ld_X1, ld_S0, ld_S1))) next
    log_ml[li] <- .d15_lmvgamma_ratio(nu1, nu0, n) +
      0.5 * n * (ld_X0 - ld_X1) +
      0.5 * nu0 * ld_S0 - 0.5 * nu1 * ld_S1 -
      0.5 * Te * n * log(pi)
  }
  list(log_ml = log_ml, A_list = A_list, Sigma_list = Sigma_list)
}


#' Helper: concatenate two .d15_eval_grid() results
#' @noRd
.d15_bind_eval <- function(a, b) {
  list(log_ml = c(a$log_ml, b$log_ml),
       A_list = c(a$A_list, b$A_list),
       Sigma_list = c(a$Sigma_list, b$Sigma_list))
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
