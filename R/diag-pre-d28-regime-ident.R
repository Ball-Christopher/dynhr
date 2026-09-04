## R/diag-pre-d28-regime-ident.R
## --------------------------------------------------------------------------
## Phase H: D28 — Regime-Switching Identification Diagnostic
##
## Assesses parameter identification in Markov-switching DSGE models.
## In a regime-switching setting, identification of structural parameters
## may differ across regimes, and the transition probabilities themselves
## become additional parameters to be identified.
##
## Algorithm:
##   1. Take a set of regime definitions (each with potentially different
##      parameter values or model structures) and a transition matrix.
##   2. For each regime, compute the identification Jacobian (as in D1/D19)
##      using regime-specific moments.
##   3. For the full Markov-switching system, compute the ergodic
##      identification Jacobian by averaging across regimes weighted by
##      the ergodic distribution.
##   4. Report per-regime vs. ergodic identification, and flag parameters
##      whose identification depends on the current regime.
##
## References:
##   Davig & Leeper (2007). Generalizing the Taylor Principle. AER.
##   Foerster et al. (2016). Markov-switching DSGE models: a solution
##     algorithm. J. Econ. Dyn. Control.
##   Iskrev, N. (2010). Local identification in DSGE models.
##   Qu, Z., & Tkachenko, D. (2012). Identification and frequency domain
##     analysis of DSGE models.
## --------------------------------------------------------------------------

#' D28. Regime-Switching Identification Diagnostic
#'
#' Assesses how parameter identification varies across regimes in a
#' Markov-switching DSGE model.  Parameters may be well-identified in
#' one regime but weakly identified in another, and the transition
#' probabilities themselves are also subject to identification analysis.
#'
#' The diagnostic:
#' \enumerate{
#'   \item For each regime state, solves the model at regime-specific
#'     parameter values and computes the identification Jacobian
#'     (moments w.r.t. parameters), following the D1/D19 approach.
#'   \item Computes an ergodic identification Jacobian by averaging
#'     regime-specific moment Jacobians weighted by the ergodic
#'     distribution of the Markov chain.
#'   \item Computes the D20-style identification strength (|t|-ratios)
#'     for both the per-regime and ergodic cases.
#'   \item Flags parameters whose identification is materially
#'     regime-dependent.
#' }
#'
#' @param regime_defs       A list of regime definitions.  Each element must
#'   be a list with:
#'   \describe{
#'     \item{\code{name}}{Character: regime name.}
#'     \item{\code{params}}{Named numeric vector: parameter values in this regime.}
#'     \item{\code{model_solve_fn}}{Function: theta -> moments for this regime.
#'       If omitted, a generic solver is used from \code{model} and \code{dr}.}
#'   }
#' @param transition_matrix Square Markov transition matrix P, where
#'   P[s, s'] = P(regime_{t+1} = s' | regime_t = s). Rows must sum to 1.
#'   rownames/colnames must match \code{names(regime_defs)}.
#' @param model             A dynhr_mod object (optional, used for generic
#'   solver when \code{model_solve_fn} is not provided per regime).
#' @param dr                Baseline DecisionRules object (optional).
#' @param params            Named baseline parameter vector. Used when
#'   \code{regime_defs} entries specify only deviations from baseline.
#' @param param_names       Character vector of parameter names to analyse.
#'   If NULL, derived from \code{params}.
#' @param eps               Step size for finite-difference Jacobian (default 1e-5).
#' @param strength_threshold Threshold for D20-style |t|-ratio below which a
#'   parameter is flagged weakly identified (default 1.0).
#' @param verbose           Print progress messages.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing:
#'     \itemize{
#'       \item \code{regime_jacobians} — named list of Jacobian matrices,
#'         one per regime.
#'       \item \code{ergodic_jacobian} — ergodic-weighted Jacobian matrix.
#'       \item \code{per_regime_strength} — matrix (n_param x n_regimes)
#'         of identification strength.
#'       \item \code{ergodic_strength} — numeric vector of ergodic
#'         identification strength.
#'       \item \code{regime_ranks} — integer vector of rank per regime.
#'       \item \code{ergodic_rank} — rank of the ergodic Jacobian.
#'       \item \code{regime_dependent_params} — character vector of
#'         regime-dependent parameters.
#'       \item \code{transition_matrix} — the transition matrix used.
#'       \item \code{ergodic_dist} — the ergodic distribution.
#'     }}
#'   \item{pass}{Logical — all parameters identified in ergodic sense AND
#'     no parameter both regime-dependent AND weakly identified.}
#'   \item{plots}{List of ggplot2 objects.}
#'   \item{summary}{Human-readable summary.}
#'
#' @noRd
d28_regime_switching_identification <- function(regime_defs,
                                                 transition_matrix,
                                                 model = NULL,
                                                 dr = NULL,
                                                 params = NULL,
                                                 param_names = NULL,
                                                 eps = 1e-5,
                                                 strength_threshold = 1.0,
                                                 verbose = FALSE,
                                                 meta = NULL) {
  # ---- 1. Validate ----
    n_regimes <- length(regime_defs)
    if (n_regimes < 1L) {
      return(.make_result(
        pass    = NA,
        summary = "D28 Regime-Switching Identification: No regime definitions provided."
      ))
    }

    regime_names <- names(regime_defs) %||% paste0("regime", seq_len(n_regimes))

    if (!is.null(transition_matrix)) {
      if (nrow(transition_matrix) != n_regimes || ncol(transition_matrix) != n_regimes) {
        stop(sprintf("transition_matrix must be %d x %d.", n_regimes, n_regimes))
      }
      stopifnot(all(abs(rowSums(transition_matrix) - 1) < 1e-6))
    }

    # ---- 2. Resolve parameter names ----
    if (is.null(param_names)) {
      param_names <- names(params)
      if (is.null(param_names)) {
        # Try from first regime
        param_names <- names(regime_defs[[1]]$params %||% params)
      }
    }
    if (is.null(param_names) || length(param_names) == 0L) {
      return(.make_result(
        pass    = NA,
        summary = "D28 Regime-Switching Identification: No parameter names available."
      ))
    }
    n_par <- length(param_names)

    # ---- 3. Compute ergodic distribution (power iteration) ----
    ergodic_dist <- NULL
    if (!is.null(transition_matrix)) {
      ergodic_dist <- .compute_ergodic_dist_ms(transition_matrix)
      if (verbose) {
        cat(sprintf("[d28] Ergodic distribution: %s\n",
                    paste(sprintf("%s=%.4f", regime_names, ergodic_dist), collapse = ", ")))
      }
    }

    # ---- 4. Per-regime identification ----
    regime_jacobians <- vector("list", n_regimes)
    names(regime_jacobians) <- regime_names
    regime_strength <- matrix(NA_real_, nrow = n_par, ncol = n_regimes,
                              dimnames = list(param_names, regime_names))
    regime_ranks <- integer(n_regimes)
    names(regime_ranks) <- regime_names

    for (r in seq_len(n_regimes)) {
      rd <- regime_defs[[r]]
      r_name <- regime_names[r]

      if (verbose) cat(sprintf("[d28]   Regime %s...", r_name))

      # Get model_solve_fn for this regime
      solve_fn <- rd$model_solve_fn
      if (is.null(solve_fn)) {
        solve_fn <- .make_regime_solve_fn(model, dr, params, rd$params, r_name)
      }

      if (is.null(solve_fn)) {
        if (verbose) cat(" no solve_fn available, skipping.\n")
        regime_ranks[r] <- NA_integer_
        next
      }

      # Get regime-specific theta
      theta_r <- params
      if (!is.null(rd$params)) {
        theta_r[names(rd$params)] <- rd$params
      }

      # Numerical Jacobian
      J <- .numerical_jacobian(function(th) solve_fn(th),
                                theta_r[param_names], eps = eps)

      if (is.null(J)) {
        if (verbose) cat(" Jacobian failed.\n")
        regime_ranks[r] <- NA_integer_
        next
      }

      colnames(J) <- param_names
      regime_jacobians[[r]] <- J

      # Rank via SVD
      sv_J <- svd(J)
      regime_ranks[r] <- .svd_rank(sv_J$d, dim(J))

      # Strength (D20-style); .safe_sym_inv handles the singular case.
      I_raw <- crossprod(J)
      I_reg <- I_raw + diag(1e-10, n_par)
      I_inv <- .safe_sym_inv(I_reg)
      se <- sqrt(pmax(diag(I_inv), 0))
      regime_strength[, r] <- abs(theta_r[param_names]) / pmax(se, 1e-16)

      if (verbose) cat(sprintf(" rank=%d\n", regime_ranks[r]))
    }

    # ---- 5. Ergodic identification ----
    ergodic_jacobian <- NULL
    ergodic_strength <- NULL
    ergodic_rank <- NA_integer_

    if (!is.null(ergodic_dist) && n_regimes > 0) {
      # Weighted average of Jacobians
      valid_regimes <- which(!vapply(regime_jacobians, is.null, logical(1)))
      if (length(valid_regimes) > 0) {
        n_moments <- min(vapply(regime_jacobians[valid_regimes], nrow, integer(1)))
        if (n_moments > 0) {
          ergodic_J <- matrix(0, nrow = n_moments, ncol = n_par)
          weights_total <- sum(ergodic_dist[valid_regimes])
          for (r in valid_regimes) {
            w <- ergodic_dist[r] / weights_total
            ergodic_J <- ergodic_J + w * regime_jacobians[[r]][seq_len(n_moments), , drop = FALSE]
          }
          colnames(ergodic_J) <- param_names
          ergodic_jacobian <- ergodic_J

          # Rank
          sv_e <- svd(ergodic_J)
          ergodic_rank <- .svd_rank(sv_e$d, dim(ergodic_J))

          # Strength (the ergodic Fisher info is often singular: shock-std
          # parameters are collinear; .safe_sym_inv falls back to a pseudo-inverse).
          I_e <- crossprod(ergodic_J)
          I_e_reg <- I_e + diag(1e-10, n_par)
          I_e_inv <- .safe_sym_inv(I_e_reg)
          se_e <- sqrt(pmax(diag(I_e_inv), 0))
          ergodic_strength <- abs(params[param_names]) / pmax(se_e, 1e-16)
          names(ergodic_strength) <- param_names
        }
      }
    }

    # ---- 6. Cross-regime dependence analysis ----
    regime_dependent <- character(0)
    weak_by_regime <- list()

    for (i in seq_len(n_par)) {
      pname <- param_names[i]
      s_i <- regime_strength[i, ]
      s_valid <- s_i[is.finite(s_i) & s_i > 0]

      weak_regimes <- which(s_i < strength_threshold)
      if (length(weak_regimes) > 0) {
        weak_by_regime[[pname]] <- regime_names[weak_regimes]
      }

      # Regime dependence: high variance across regimes
      if (length(s_valid) >= 2) {
        cv <- sd(s_valid) / max(mean(s_valid), 1e-16)
        if (cv > 1.0) {  # Coefficient of variation > 100%
          regime_dependent <- c(regime_dependent, pname)
        }
      }
    }

    # Gate:
    #   - pass = NA when n_regime_dependent == 0 (nothing tested; vacuous)
    #   - FAIL if any param is weak in ALL regimes (model is unidentified
    #     regardless of regime, so regime-switching provides no help)
    #   - FAIL if any param is both regime-dependent AND weakly identified
    #     in at least one regime
    n_regime_dependent <- length(regime_dependent)
    # Params weak in ALL regimes (every regime has s_i < threshold)
    n_valid_regimes <- sum(!is.na(regime_ranks))
    params_weak_in_all <- character(0)
    for (i in seq_len(n_par)) {
      s_i <- regime_strength[i, ]
      n_weak_regimes <- sum(is.finite(s_i) & s_i < strength_threshold, na.rm = TRUE)
      n_finite <- sum(is.finite(s_i))
      if (n_finite > 0 && n_weak_regimes == n_finite) {
        params_weak_in_all <- c(params_weak_in_all, param_names[i])
      }
    }

    pass <- if (n_regime_dependent == 0L) {
      NA  # Nothing to test; D28 is vacuous for this model
    } else if (length(params_weak_in_all) > 0) {
      FALSE  # Some params are globally unidentified regardless of regime
    } else {
      all_ok <- TRUE
      for (pname in regime_dependent) {
        if (pname %in% names(weak_by_regime)) { all_ok <- FALSE; break }
      }
      all_ok
    }

    # ---- 7. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE) && n_regimes > 1) {
      strength_df <- as.data.frame.table(regime_strength,
                                          responseName = "Strength",
                                          stringsAsFactors = FALSE)
      colnames(strength_df) <- c("Parameter", "Regime", "Strength")

      p_rs <- ggplot2::ggplot(
        strength_df, ggplot2::aes(x = Regime, y = Parameter, fill = Strength)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        scale_fill_dynhr_cividis(name = "|t|-ratio", na.value = "grey50") +
        theme_dynhr_diagnostic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
          title = "D28: Regime-dependent identification strength",
          subtitle = sprintf("%d regimes, ergodic rank = %s",
                             n_regimes,
                             if (is.na(ergodic_rank)) "N/A" else as.character(ergodic_rank)),
          x = NULL, y = NULL
        )
      plots$regime_strength <- .apply_meta(p_rs, meta)

      # Ergodic vs per-regime comparison
      if (!is.null(ergodic_strength)) {
        comp_df <- data.frame(
          Parameter = rep(param_names, n_regimes + 1L),
          Regime = c(rep(regime_names, each = n_par), rep("ergodic", n_par)),
          Strength = c(as.vector(regime_strength), ergodic_strength),
          stringsAsFactors = FALSE
        )
        p_ec <- ggplot2::ggplot(
          comp_df, ggplot2::aes(x = Regime, y = Parameter, fill = Strength)
        ) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
          scale_fill_dynhr_cividis(name = "|t|-ratio", na.value = "grey50") +
          theme_dynhr_diagnostic() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(
            title = "D28: Per-regime vs ergodic identification",
            x = NULL, y = NULL
          )
        plots$ergodic_comparison <- .apply_meta(p_ec, meta)
      }
    }

    # ---- 8. Build result ----
    result <- list(
      regime_jacobians          = regime_jacobians,
      ergodic_jacobian          = ergodic_jacobian,
      per_regime_strength       = regime_strength,
      ergodic_strength          = ergodic_strength,
      regime_ranks              = regime_ranks,
      ergodic_rank              = ergodic_rank,
      regime_dependent_params   = regime_dependent,
      weak_params_by_regime     = weak_by_regime,
      params_weak_in_all_regimes = params_weak_in_all,
      transition_matrix         = transition_matrix,
      ergodic_dist              = ergodic_dist,
      n_regimes                 = n_regimes,
      regime_names              = regime_names
    )

    rank_str <- if (is.na(ergodic_rank))
      sprintf("per-regime ranks: [%s]", paste(regime_ranks, collapse = ", "))
    else sprintf("ergodic rank = %d", ergodic_rank)

    dep_str <- if (n_regime_dependent > 0L)
      sprintf("%d regime-dependent param(s): %s", n_regime_dependent,
              paste(regime_dependent, collapse = ", "))
    else "No regime-dependent parameters (D28 is vacuous for this model)."

    .make_result(
      result  = result,
      pass    = pass,
      plots   = plots,
      summary = sprintf("D28 Regime-Switching Identification: %d regimes. %s. %s%s",
                        n_regimes, rank_str, dep_str,
                        if (length(params_weak_in_all) > 0)
                          sprintf(" GLOBALLY WEAK (all regimes): %s.",
                                  paste(params_weak_in_all, collapse = ", "))
                        else ""),
      llm_summary = {
        badge <- if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL"
        paste(c(
          sprintf("D28 | Regime-Switching Identification | %s", badge),
          sprintf("  n_regimes=%d ergodic_rank=%s n_regime_dependent=%d",
                  n_regimes,
                  if (is.na(ergodic_rank)) "NA" else as.character(ergodic_rank),
                  n_regime_dependent),
          if (is.na(pass))
            "  note: No regime-dependent parameters detected; D28 is vacuous for this model.",
          if (length(params_weak_in_all) > 0)
            sprintf("  globally_weak (all regimes): %s",
                    paste(params_weak_in_all, collapse = ", ")),
          sprintf("  action: %s",
                  if (is.na(pass))
                    "No regime-dependent variation detected. Check transition matrix and regime definitions."
                  else if (length(params_weak_in_all) > 0)
                    sprintf("Parameters %s are weakly identified in ALL regimes. Regime-switching does not resolve identification -- add moments or calibrate.",
                            paste(head(params_weak_in_all, 3), collapse = ", "))
                  else if (isTRUE(pass))
                    "Regime-switching identification adequate."
                  else
                    sprintf("Regime-dependent parameters %s are weakly identified in some regime. Consider richer measurement or fewer regime-specific parameters.",
                            paste(head(regime_dependent[regime_dependent %in% names(weak_by_regime)], 3), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}


# ==========================================================================
# Internal helpers for D28
# ==========================================================================

#' Build a regime-specific solve function
#'
#' @param model      dynhr_mod
#' @param dr         Baseline DecisionRules
#' @param params     Baseline params
#' @param regime_params  Regime-specific param overrides
#' @param regime_name    Regime name (for messages)
#' @return Function: theta -> moments, or NULL if not constructable
#' @noRd
.make_regime_solve_fn <- function(model, dr, params, regime_params, regime_name) {
  if (is.null(model) && is.null(dr)) return(NULL)

  # Merge params
  theta_base <- params
  if (!is.null(regime_params)) {
    theta_base[names(regime_params)] <- regime_params
  }

  function(theta) {
    theta_full <- theta_base
    theta_full[names(theta)] <- theta
    {
      dr_new <- .stoch_simul_internal_diag(model, theta_full, dr_order = 1L)
      if (is.null(dr_new) || is.null(dr_new$ghx)) return(NULL)
      moments <- .moments_from_dr(dr_new, model = model, params = theta_full)
      if (is.null(moments)) return(NULL)
      acf_y_diag <- if (length(dim(moments$acf_y)) == 3L) {
        diag(moments$acf_y[, , 1])
      } else {
        diag(as.matrix(moments$acf_y))
      }
      m <- c(
        log_var  = log(pmax(diag(moments$sigma_y), 1e-16)),
        log_acf1 = log(pmax(acf_y_diag, 1e-16))
      )
      names(m) <- paste0("m_", seq_along(m))
      m
    }
  }
}


#' Compute ergodic distribution via power iteration
#'
#' Duplicate-safe wrapper: prefers the existing .compute_ergodic_dist from
#' ramsey-regime.R if available, otherwise uses a local implementation.
#'
#' @param P Transition matrix (rows sum to 1)
#' @return Numeric vector of ergodic probabilities
#' @noRd
.compute_ergodic_dist_ms <- function(P) {
  if (exists(".compute_ergodic_dist", mode = "function")) {
    return(.compute_ergodic_dist(P))
  }
  n <- nrow(P)
  pi <- rep(1/n, n)
  for (iter in seq_len(10000)) {
    pi_new <- pi %*% P
    pi_new <- as.numeric(pi_new)
    if (max(abs(pi_new - pi)) < 1e-14) break
    pi <- pi_new
  }
  pi / sum(pi)
}
