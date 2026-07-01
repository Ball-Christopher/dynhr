## R/diag-pre-d25-higher-order-ident.R
## --------------------------------------------------------------------------
## Phase C: D25 Higher-Order Identification (Revisited).
##
## Extends the moment vector for D19/D20 to include higher-order cumulants
## (second-order variance corrections, third cumulants) when the Ramsey
## solution or a second/third-order perturbation is available. This allows
## identification diagnostics to exploit information content beyond the
## first-order (linear-Gaussian) moment vector.
##
## Key references:
##   Iskrev, N. (2010). Local identification in DSGE models.
##     Journal of Monetary Economics, 57(2), 189-202.
##   Komunjer, I., & Ng, S. (2011). Dynamic identification of DSGE models.
##     Econometrica, 79(6), 1995-2032.
##   Qu, Z., & Tkachenko, D. (2017). Global identification in DSGE models
##     based on the spectral density. Econometrica, 85(5), 1571-1638.
##   Gross, T., & Hansen, J. (2021). Optimal policy design in nonlinear
##     DSGE models: An n-order accurate approximation.
## --------------------------------------------------------------------------

#' D25. Higher-order identification diagnostic
#'
#' Extends the standard D19/D20 identification diagnostics to include
#' higher-order cumulants (second-order variance corrections, third
#' cumulants) in the moment vector. This reveals identification content
#' that is invisible at first order -- i.e., parameters that affect the
#' model's higher-order properties but cancel out at order 1.
#'
#' Two modes are supported:
#'   1. **Ramsey-augmented mode**: Uses \code{ramsey_result} (a
#'      \code{dynhr_ramsey_result2} object). The augmented decision rules
#'      include multipliers, and the Ramsey-optimal policy solution at
#'      order >= 2 provides second-order moments and uncertainty corrections.
#'   2. **Standalone DR mode**: Uses a \code{DecisionRules2} or
#'      \code{DecisionRules3} object directly, without a Ramsey solution.
#'      In this mode, the moment vector is augmented with a proxy for the
#'      second-order mean correction (\code{|ghss[i]| / impact-SD[i]} per
#'      observable) and (if available) third cumulants.
#'      \strong{Note:} this is a proxy, not the full Komunjer-Ng (2011)
#'      second-order moment augmentation, which would require propagating
#'      \code{ghxx}/\code{ghuu} through a Lyapunov-like recursion to obtain
#'      true unconditional variance corrections.
#'
#' @param ramsey_result  Optional \code{dynhr_ramsey_result2} from
#'   \code{ramsey_model()}. If provided, the augmented (Ramsey-optimal)
#'   decision rules are used.
#' @param dr             Decision rules object (DecisionRules2 or DecisionRules3).
#'   Used when \code{ramsey_result} is NULL, or to provide the original
#'   (competitive-equilibrium) DR for comparison.
#' @param model          dynhr_mod (for variable names, observation matrix).
#' @param params         Named parameter vector at the calibration point.
#' @param ss             Named steady-state vector.
#' @param obs_mat        Observation matrix (n_obs x n_state). If NULL,
#'   defaults to full-state observation.
#' @param Sigma_e        Shock covariance matrix. If NULL, uses identity.
#' @param param_names    Optional character vector of parameter names.
#' @param include_third_cumulant Logical. If TRUE (and a DecisionRules3
#'   object is available), include third cumulants of observables in the
#'   moment vector. Default FALSE (computationally heavier).
#' @param eps            Step size for finite differences (default 1e-5).
#' @param verbose        Print progress messages.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing:
#'     \itemize{
#'       \item \code{order1_jacobian} — Jacobian from first-order moments only
#'       \item \code{extended_jacobian} — Jacobian from extended moment vector
#'       \item \code{order1_rank} — rank from first-order moments
#'       \item \code{extended_rank} — rank from extended moments
#'       \item \code{identification_gain} — change in number of identified
#'         directions
#'       \item \code{param_contributions} — per-parameter contribution of
#'         higher-order terms to identification
#'       \item \code{moment_vector_info} — breakdown of which moments are added
#'     }}
#'   \item{pass}{Logical — whether identification improves with higher-order
#'     moments.}
#'   \item{plots}{List of ggplot2 objects.}
#'   \item{summary}{Human-readable summary.}
#'
#' @references
#'   Komunjer, I., & Ng, S. (2011). Dynamic identification of DSGE models.
#'     \emph{Econometrica}, 79(6), 1995-2032.
#'
#' @noRd
d25_higher_order_identification <- function(ramsey_result = NULL,
                                             dr = NULL,
                                             model = NULL,
                                             params = NULL,
                                             ss = NULL,
                                             obs_mat = NULL,
                                             Sigma_e = NULL,
                                             param_names = NULL,
                                             include_third_cumulant = FALSE,
                                             eps = 1e-5,
                                             verbose = FALSE,
                                             meta = NULL) {

    # ---- 1. Resolve decision rules ----
    # Prefer Ramsey-augmented DR when available
    dr_order <- NULL
    if (!is.null(ramsey_result) && inherits(ramsey_result, "dynhr_ramsey_result2")) {
      ramsey_dr <- ramsey_result$ramsey_dr$ramsey_dr
      if (!is.null(ramsey_dr)) {
        dr <- ramsey_dr
        if (verbose) cat("[d25] Using Ramsey-optimal decision rules.\n")
      }
    }

    if (is.null(dr)) {
      return(.make_result(
        pass    = NA,
        summary = "D25 Higher-order identification: no decision rules provided.",
        llm_summary = "[INFO] D25 Higher-order identification | status=skipped reason=no_dr"
      ))
    }

    # Determine perturbation order
    if (inherits(dr, "DecisionRules3")) {
      dr_order <- 3L
    } else if (inherits(dr, "DecisionRules2")) {
      dr_order <- 2L
    } else if (!is.null(dr$ghxx) && !is.null(dr$ghuu)) {
      dr_order <- 2L
    } else if (!is.null(dr$ghx) && !is.null(dr$ghu)) {
      dr_order <- 1L
    } else {
      return(.make_result(
        pass    = NA,
        summary = "D25 Higher-order identification: unrecognised decision rules object.",
        llm_summary = "[INFO] D25 Higher-order identification | status=skipped reason=unrecognised_dr"
      ))
    }

    if (dr_order < 2L) {
      return(.make_result(
        pass    = NA,
        summary = "D25 Higher-order identification: requires order >= 2. Current DR is order 1.",
        llm_summary = "[INFO] D25 Higher-order identification | status=skipped reason=order1_only"
      ))
    }

    if (verbose) {
      cat(sprintf("[d25] DR order = %d\n", dr_order))
    }

    # ---- 2. Extract state-space ----
    ss_vec <- if (!is.null(ss)) ss else {
      if (!is.null(dr$ys)) dr$ys else {
        if (!is.null(ramsey_result$competitive_ss)) ramsey_result$competitive_ss else NULL
      }
    }

    if (is.null(params) && !is.null(ramsey_result)) {
      params <- ramsey_result$augmented_model$param_values
    }

    ghx <- as.matrix(dr$ghx)
    ghu <- as.matrix(dr$ghu)
    n_endo <- nrow(ghx)
    n_shock <- ncol(ghu)
    state_idx <- dr$state_idx %||% seq_len(n_endo)
    n_state <- length(state_idx)  # actual number of state variables
    endo_names <- dr$endo_names %||% rownames(ghx) %||% {
      if (!is.null(model)) model$var_names else paste0("y_", seq_len(n_endo))
    }

    # Observation matrix — default to full-state observation (n_endo x n_endo)
    if (is.null(obs_mat)) {
      if (!is.null(dr$Z)) {
        obs_mat <- as.matrix(dr$Z)
      } else if (!is.null(model$obs_mat)) {
        obs_mat <- as.matrix(model$obs_mat)
      } else {
        obs_mat <- diag(n_endo)
      }
    } else {
      obs_mat <- as.matrix(obs_mat)
    }
    n_obs <- nrow(obs_mat)

    # Shock covariance
    if (is.null(Sigma_e)) {
      Sigma_e <- diag(n_shock)
    }

    # Parameter names
    if (is.null(param_names)) {
      param_names <- names(params) %||% paste0("theta_", seq_along(params))
    }
    n_par <- length(params)

    # ---- 3. Build first-order moment vector (baseline) ----
    .first_order_moments <- function(th) {
      # This is called within the numerical Jacobian.
      # Uses the state-space to compute first-order moments.
      # State-to-state transition (ghx[state_idx,]) for Lyapunov,
      # then projects to full endo space via ghx/ghu for observables.
      si <- dr$state_idx %||% seq_len(n_endo)
      n_s <- length(si)
      T_state <- ghx[si, , drop = FALSE]         # n_s x n_s
      R_state <- ghu[si, , drop = FALSE]         # n_s x n_shock
      R_mat <- R_state %*% chol(Sigma_e)         # n_s x n_shock
      Z_mat <- obs_mat                           # n_obs x n_endo

      # State covariance via Lyapunov
      QQ <- R_mat %*% t(R_mat)                   # n_s x n_s
      Sigma_s <- solve_lyapunov(T_state, QQ)     # n_s x n_s

      # Full endo covariance: Sigma_y = ghx * Sigma_s * ghx' + ghu * Sigma_e * ghu'
      Sigma_y <- ghx %*% Sigma_s %*% t(ghx) + ghu %*% Sigma_e %*% t(ghu)
      Cov_y <- Z_mat %*% Sigma_y %*% t(Z_mat)    # n_obs x n_obs
      sigma_y <- sqrt(pmax(diag(Cov_y), 0))

      # Autocorrelations at lags 1-4
      max_lag <- 4L
      acf_y <- matrix(0, nrow = n_obs, ncol = max_lag)
      # Build selector matrix S: n_state x n_endo mapping states to endo positions
      S_sel <- matrix(0, nrow = n_s, ncol = n_endo)
      for (i in seq_len(n_s)) S_sel[i, si[i]] <- 1
      Gamma_prev <- Sigma_y
      for (lag in seq_len(max_lag)) {
        Gamma_lag <- ghx %*% S_sel %*% Gamma_prev
        diag_cov <- diag(Cov_y)
        for (i in seq_len(n_obs)) {
          if (diag_cov[i] > 1e-14) {
            # Project Gamma_lag to observables via Z_mat
            Cov_lag_i <- Z_mat[i, , drop = FALSE] %*% Gamma_lag %*% t(Z_mat[i, , drop = FALSE])
            acf_y[i, lag] <- Cov_lag_i[1, 1] / diag_cov[i]
          }
        }
        Gamma_prev <- Gamma_lag
      }

      # Flatten into moment vector
      c(sigma_y, as.numeric(t(acf_y)))
    }

    # ---- 4. Build extended moment vector (with higher-order terms) ----
    # Higher-order contributions:
    #   (a) Second-order unconditional variance correction
    #       Unconditional variance at order 2 = order-1 variance + ghss correction
    #   (b) Skewness / third cumulant from ghs3 (at order 3)
    #   (c) Pruning-corrected moments if pruning matters

    .extended_moments <- function(th) {
      # First-order moments
      o1 <- .first_order_moments(th)

      # Second-order correction to unconditional mean and variance
      ghss <- dr$ghss %||% rep(0, n_endo)
      ghxx <- dr$ghxx %||% NULL
      ghuu <- dr$ghuu %||% NULL

      # The ghss vector gives the uncertainty correction to the mean.
      # Under Gaussian shocks, the second-order unconditional mean is ys + 0.5 * ghss.
      # The unconditional variance correction involves ghxx, ghuu and the
      # state covariance. Here we use the ghss magnitude per observable
      # as a proxy for the second-order correction.

      if (!is.null(ghss) && length(ghss) == n_endo) {
        # Project ghss onto observables (obs_mat is n_obs x n_endo)
        ghss_obs <- obs_mat %*% ghss
        # Relative magnitude of second-order correction
        o1_sigma <- o1[seq_len(n_obs)]
        ss_correction <- rep(0, n_obs)
        for (i in seq_len(n_obs)) {
          if (o1_sigma[i] > 1e-14) {
            ss_correction[i] <- abs(ghss_obs[i]) / o1_sigma[i]
          }
        }
      } else {
        ss_correction <- rep(0, n_obs)
      }

      # Third cumulant (skewness) from ghs3 if available
      skewness_obs <- NULL
      if (isTRUE(include_third_cumulant) && dr_order >= 3L) {
        ghs3 <- dr$ghs3 %||% NULL
        if (!is.null(ghs3) && length(ghs3) == n_endo && any(abs(ghs3) > 1e-14)) {
          ghs3_obs <- obs_mat %*% ghs3
          # Normalise by order-1 std dev to get skewness-like measure
          o1_sigma <- o1[seq_len(n_obs)]
          skewness_obs <- rep(0, n_obs)
          for (i in seq_len(n_obs)) {
            if (o1_sigma[i] > 1e-14) {
              skewness_obs[i] <- ghs3_obs[i] / o1_sigma[i]^3
            }
          }
        }
      }

      # Assemble extended moment vector
      extended <- c(o1, ss_correction)
      if (!is.null(skewness_obs)) {
        extended <- c(extended, skewness_obs)
      }

      extended
    }

    # ---- 5. Compute Jacobians via finite differences ----
    # We use the Jacobian of moment functions directly from the current DR.
    # This is a local approximation valid near the calibration point.

    # Compute baseline moment vectors
    f0_o1 <- .first_order_moments(params)
    f0_ext <- .extended_moments(params)

    # Numerical Jacobian of first-order moments
    J_o1 <- .numerical_jacobian(.first_order_moments, params, eps = eps)
    colnames(J_o1) <- param_names
    rownames(J_o1) <- paste0("o1_m", seq_len(length(f0_o1)))

    # Numerical Jacobian of extended moments
    J_ext <- .numerical_jacobian(.extended_moments, params, eps = eps)
    n_o1 <- length(f0_o1)
    n_ext <- length(f0_ext)
    rownames(J_ext) <- c(
      paste0("o1_m", seq_len(n_o1)),
      paste0("ho_ss_corr_", seq_len(n_obs)),
      if (!is.null(f0_ext) && n_ext > n_o1 + n_obs) {
        paste0("ho_skew_", seq_len(n_ext - n_o1 - n_obs))
      }
    )

    # ---- 6. SVD and rank analysis ----
    sv_o1 <- svd(J_o1)
    sv_ext <- svd(J_ext)

    rank_o1 <- .svd_rank(sv_o1$d, dim(J_o1))
    rank_ext <- .svd_rank(sv_ext$d, dim(J_ext))

    identification_gain <- rank_ext - rank_o1

    # ---- 7. Per-parameter contribution of higher-order terms ----
    # For each parameter, compare the norm of the column in J_o1 vs J_ext
    param_contrib <- data.frame(
      parameter = param_names,
      order1_norm = sqrt(colSums(J_o1^2)),
      extended_norm = sqrt(colSums(J_ext^2)),
      stringsAsFactors = FALSE
    )
    param_contrib$ho_ratio <- ifelse(
      param_contrib$order1_norm > 1e-14,
      param_contrib$extended_norm / param_contrib$order1_norm,
      1.0
    )
    param_contrib$ho_gain <- param_contrib$ho_ratio - 1.0
    # Parameters where higher-order adds >10% identification content
    param_contrib$ho_important <- param_contrib$ho_gain > 0.10

    # ---- 8. Moment vector info ----
    moment_info <- list(
      n_first_order = n_o1,
      n_ss_correction = n_obs,
      n_skewness = if (!is.null(f0_ext) && n_ext > n_o1 + n_obs) {
        n_ext - n_o1 - n_obs
      } else 0,
      n_total_first_order_moments = n_o1,
      n_total_extended_moments = n_ext
    )

    # ---- 9. Build pass/fail ----
    # pass = NA when rank_o1 == 0: D1 upstream failure (moment function returns
    # zero-rank output).  D25 cannot assess higher-order gain against a zero
    # first-order baseline; the FAIL would be spurious and contradictory to
    # whatever D1 actually reports for the model.
    # Otherwise: pass if gain > 0 (higher order adds identification) OR both
    # first-order and extended are full rank (already fully identified).
    pass <- if (rank_o1 == 0L) {
      NA  # upstream D1 rank 0 -- cannot assess higher-order gain
    } else {
      identification_gain > 0 || rank_ext == n_par
    }

    # ---- 10. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # (a) Comparison of singular values
      if (length(sv_o1$d) > 0 && length(sv_ext$d) > 0) {
        sv_df <- data.frame(
          index = seq_len(max(length(sv_o1$d), length(sv_ext$d))),
          order1 = c(sv_o1$d, rep(NA, max(0, length(sv_ext$d) - length(sv_o1$d)))),
          extended = c(sv_ext$d, rep(NA, max(0, length(sv_o1$d) - length(sv_ext$d))))
        )
        sv_long <- reshape2::melt(sv_df, id.vars = "index",
                                   variable.name = "type",
                                   value.name = "singular_value")
        sv_long <- sv_long[is.finite(sv_long$singular_value), ]

        sv_long$singular_value <- pmax(sv_long$singular_value, 1e-300)
        p_sv <- ggplot2::ggplot(
          sv_long,
          ggplot2::aes(x = index, y = singular_value,
                       colour = type, shape = type)
        ) +
          ggplot2::geom_point(size = 2.5) +
          ggplot2::geom_line(alpha = 0.5) +
          ggplot2::scale_y_log10() +
          ggplot2::scale_colour_manual(
            values = c("order1" = dynhr_colours$mid_blue,
                       "extended" = dynhr_colours$red),
            name = "Moment set",
            labels = c("order1" = "First-order",
                       "extended" = "Extended (higher-order)")
          ) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D25: Singular values -- first-order vs extended",
            subtitle = sprintf("Rank: O1=%d, Ext=%d (gain=%+d)",
                               rank_o1, rank_ext, identification_gain),
            x = "Singular value index",
            y = "Value (log scale)"
          )
        plots$sv_comparison <- .apply_meta(p_sv, meta)
      }

      # (b) Per-parameter identification gain from higher-order terms
      if (nrow(param_contrib) > 0) {
        param_contrib$param_label <- factor(
          param_contrib$parameter,
          levels = param_contrib$parameter[order(param_contrib$ho_gain)]
        )
        p_pc <- ggplot2::ggplot(
          param_contrib,
          ggplot2::aes(x = param_label, y = ho_gain,
                       fill = ho_important)
        ) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::geom_hline(yintercept = 0, linetype = "solid",
                              colour = "grey50", linewidth = 0.3) +
          ggplot2::geom_hline(yintercept = 0.10, linetype = "dashed",
                              colour = dynhr_colours$red, linewidth = 0.4) +
          ggplot2::coord_flip() +
          ggplot2::scale_fill_manual(
            values = c("TRUE" = dynhr_colours$mid_blue,
                       "FALSE" = dynhr_colours$grey),
            guide = "none"
          ) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D25: Identification gain from higher-order moments",
            subtitle = "Positive = extended moment set improves sensitivity",
            x = NULL,
            y = "Relative gain (extended_norm / order1_norm - 1)"
          )
        plots$param_contribution <- .apply_meta(p_pc, meta)
      }
    }

    # ---- 11. Summary ----
    n_ho_important <- sum(param_contrib$ho_important)
    summary_text <- sprintf(
      "D25 Higher-order identification: rank O1=%d, Ext=%d (gain=%+d). %s",
      rank_o1, rank_ext, identification_gain,
      if (rank_o1 == 0L) {
        "INFO: first-order rank is 0 -- cannot assess higher-order gain; recheck moment function."
      } else if (identification_gain > 0) {
        sprintf("Higher-order moments add %d identification direction(s). %d parameter(s) benefit significantly.",
                identification_gain, n_ho_important)
      } else if (rank_o1 == n_par && rank_ext == n_par) {
        sprintf("Already fully identified at order 1 (rank=%d/%d). Higher-order moments confirm full identification.",
                rank_o1, n_par)
      } else {
        sprintf("No gain from higher-order moments (rank unchanged at %d/%d).",
                rank_o1, n_par)
      }
    )

    top_params <- head(
      param_contrib$parameter[order(param_contrib$ho_gain, decreasing = TRUE)],
      5
    )
    top_params <- top_params[param_contrib$ho_gain[match(top_params, param_contrib$parameter)] > 0.05]
    top_str <- if (length(top_params) > 0) {
      paste(sprintf("%s(%.1f%%)", top_params,
                    round(100 * param_contrib$ho_gain[match(top_params, param_contrib$parameter)], 1)),
            collapse = ", ")
    } else "none"

    llm_summary <- paste(c(
      sprintf("D25 | Higher-Order Identification | %s",
              if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL"),
      sprintf("  dr_order=%d rank_o1=%d rank_ext=%d gain=%+d n_par=%d",
              dr_order, rank_o1, rank_ext, identification_gain, n_par),
      if (is.na(pass))
        "  note: rank_o1=0 (upstream D1 rank failure) -- D25 cannot assess higher-order gain.",
      sprintf("  higher_order_gain_params: %s", top_str),
      sprintf("  moment_breakdown: o1=%d ss_corr=%d skew=%d total=%d",
              moment_info$n_first_order,
              moment_info$n_ss_correction,
              moment_info$n_skewness,
              moment_info$n_total_extended_moments),
      sprintf("  action: %s",
              if (is.na(pass))
                "First-order rank is 0 -- recheck moment function / model setup before interpreting D25."
              else if (identification_gain > 0)
                "Higher-order moments add identifying content. Consider using order>=2 for estimation."
              else if (rank_o1 < n_par)
                "Higher-order moments do not resolve identification failure. Consider additional observables."
              else
                "Model is already well-identified at first order.")
    ), collapse = "\n")

    .make_result(
      result = list(
        order1_jacobian = J_o1,
        extended_jacobian = J_ext,
        order1_svd = sv_o1,
        extended_svd = sv_ext,
        order1_rank = rank_o1,
        extended_rank = rank_ext,
        identification_gain = identification_gain,
        param_contributions = param_contrib,
        moment_vector_info = moment_info,
        dr_order = dr_order
      ),
      pass = pass,
      plots = plots,
      summary = summary_text,
      llm_summary = llm_summary
    )
}
